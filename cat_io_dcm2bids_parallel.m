function varargout = cat_io_dcm2bids_parallel(action,varargin)
%cat_io_dcm2bids_parallel. Subject-wise (parallel) processing of cat_io_dcm2bids.
%  The scans of a run are processed in subject packages, i.e. all scans of
%  a BIDS subject (PatientID, for subject labels with site per site, see
%  packages) together, as the sessions of a subject depend on each other
%  (e.g. the T1w of another session for sessions without anatomical scan
%  and the sessions file of the subject). A package converts its images
%  (see convertSources in cat_io_dcm2bids_import), processes its scans (see
%  cat_io_dcm2bids_process), renders the subject and session reports of
%  its subjects (with the slice cache for the rendering of all subjects,
%  see cat_io_dcm2bids_render), and writes the participant reports of its
%  sessions (see cat_io_dcm2bids_participantreport), with one line of the
%  processing times of the package (conversion, QC and export, 
%  preprocessing, and reports) and the subject number of the run.
%  For nproc > 1, nproc background MATLAB processes (workers, see
%  cat_parallelize2) take the next package of a queue until all packages
%  are processed (i.e. long and short packages are balanced), whereas the
%  main process shows the log of each finished package (otherwise the
%  packages are processed one after another by the main process).
%  Outputs of all subjects (tables of the database, participants, DICOM
%  export, report, rendering of all subjects) are written by cat_io_dcm2bids
%  with the merged results of all packages (in the order of the packages,
%  i.e. the same order as without parallel processing).
%  Packages of workers that did not start (e.g. without MATLAB license) are
%  processed by the main process, whereas packages of failed workers (e.g.
%  out of memory) are reported as errors (processed again by the next run).
%  The queue directory catDCM2BIDSpar_* of the database (context, results,
%  and logs of the packages) is removed at the end, but kept after errors,
%  whereas the logs of the workers are kept in the logs directory of the
%  tables of the database.
%
%  varargout = cat_io_dcm2bids_parallel(action,varargin)
%
%  Actions (see the help of the local functions):
%    R  = cat_io_dcm2bids_parallel('run',C,Psnr,nproc)
%    pk = cat_io_dcm2bids_parallel('packages',Psnr[,C])
%    R  = cat_io_dcm2bids_parallel(job)   (package of cat_parallelize2)
%
%  See also cat_io_dcm2bids, cat_io_dcm2bids_process, cat_parallelize2.

  if isstruct(action) % (package with job.data, job.C, job.render, and job.report)
    if ~isfield(action,'report'), action.report = 0; end
    varargout{1} = processPackage(action.C,action.data,action.render,action.report);
    return
  end
  switch action
    case {'run','packages'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_parallel:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function R = run(C,Psnr,nproc)
%run. Process the scan directories Psnr (sorted) with the context C of the
%  run (see cat_io_dcm2bids) in subject packages by nproc background
%  processes (one process for nproc = 0 or one package). R contains the
%  merged results of all packages (see processScans in
%  cat_io_dcm2bids_process) with R.rendered (subject reports rendered by
%  the packages), R.reported (participant reports written by the packages),
%  and R.failed (number of failed packages).
  pk = packages(Psnr,C);
  nw = min(nproc,numel(pk));
  % subject number of the run of each subject directory of the database
  % (the package number, for the subject lines of the packages)
  C.subnum = containers.Map('KeyType','char','ValueType','double'); C.nsub = numel(pk); 
  for k = 1:numel(pk)
    for P = unique(cellfun(@(p) fileparts(fileparts(p)),Psnr(pk{k}),'UniformOutput',false))', C.subnum(P{1}) = k; end
  end
  % subject reports (otherwise the BIDS images are zipped again, see
  % gzipunzipOutputdata) and participant reports (at the end of the run
  % with the CAT preprocessing, see cat_io_dcm2bids) of the packages
  opts   = C.job.opts; 
  render = opts.render.run && opts.output > 1 && opts.gzipi == opts.gzipe; 
  report = isfield(opts,'participantreport') && opts.participantreport.run && opts.output > 1 && ...
    ~(isfield(opts,'cat') && opts.cat.run); 
  if nw < 2 % (one process: the packages one after another)
    res = cell(1,numel(pk)); 
    for k = 1:numel(pk), res{k} = processPackage(C,Psnr(pk{k}),render,report); end
    R = mergeResults(res); R.rendered = render; R.reported = report; R.failed = 0;
    return
  end

  % packages by the workers of cat_parallelize2 with the queue directory in
  % the database and the logs of the packages with colors
  job    = struct('nproc',nw,'data',{Psnr},'C',C,'render',render,'report',report);
  popts  = struct('packages',{pk},'replay',1,'mergefun',@mergeResults,'qroot',C.Pmdbdir, ...
    'qprefix','catDCM2BIDSpar_','logdir',fullfile(C.Ptbldir,'logs'),'logname','log_DCM2BIDS','name','DCM2BIDS');
  [R,status] = cat_parallelize2(job,mfilename,'data',popts);
  failed = nnz(~cellfun('isempty',{status.err}));
  if isempty(R)
    error('cat_io_dcm2bids_parallel:failed','All %d packages failed.',numel(pk));
  end
  R.rendered = render; R.reported = report; R.failed = failed;
  if failed > 0
    cat_io_cprintf('err','  %d package(s) failed, i.e. their scans are processed again by the next run.\n',failed);
  end
end
% =========================================================================
function pk = packages(Psnr,C)
%packages. Subject packages of the scan directories Psnr (sorted,
%  catDCM2BIDSdb/sub-<serial>-<PatientID>/ses-*/snr-*) as cell of indices,
%  i.e. all scan directories of the same BIDS subject in the order of their
%  first scan directory: the same PatientID (BIDS label), i.e. also of 
%  several scanners for subject labels without site (opts.subIDform 1, 
%  e.g. for the T1w of other sessions), and otherwise the same site label
%  and PatientID (opts.subIDform 2, i.e. sub-SITE-PID, with
%  the site label of the scanner, see setupSites in cat_io_dcm2bids_protocols,
%  i.e. scanners with the same label of the site dictionary stay together), 
%  e.g. the scanners of a traveling subject are processed in parallel. 
%  Without context C of the run (job options and sites), the packages are
%  given by the PatientID.
  persite = exist('C','var') && isfield(C,'job') && isfield(C.job.opts,'subIDform') && ...
    C.job.opts.subIDform == 2;
  Psub = cellfun(@(p) fileparts(fileparts(p)),Psnr,'UniformOutput',false);
  [usub,~,si] = unique(Psub,'stable');
  key = usub;
  for ui = 1:numel(usub)
    Pjson = cat_io_dcm2bids_db('getDBscanJSONs',Psnr{find(si == ui,1)});
    for ji = 1:numel(Pjson)
      try
        V = cat_io_json(Pjson{ji});
        if isfield(V,'PatientID') && ~isempty(V.PatientID)
          key{ui} = lower(strrep(cat_io_dcm2bids_bids('bidsLabel',char(V.PatientID)),'_','')); 
          if persite
            V = cat_io_dcm2bids_db('assurePatientDCMfields',V); 
            key{ui} = [cat_io_dcm2bids_protocols('setupSites',V,C.sites,C.job.opts) '|' key{ui}]; 
          end
          break
        end
      catch % (unreadable sidecar: the subject directory as package)
      end
    end
  end
  [~,~,ki] = unique(key(si),'stable');
  pk = arrayfun(@(k) find(ki == k)',1:max([0; ki(:)]),'UniformOutput',false);
end
% =========================================================================
function R = processPackage(C,Psnr,render,report)
%processPackage. Convert the images of the scan directories Psnr (only
%  into their subject directories, see convertSources in
%  cat_io_dcm2bids_import), process the scans (see cat_io_dcm2bids_process),
%  render the subject and session reports of its subjects (render), and
%  write the participant reports of its sessions (report), with the line
%  of the processing times: conversion (DCM2NII), protocol test, export,
%  and QC of the scans (QC/export), preprocessing of the sessions (PP, see
%  R.tpp of processScans), and reports.
  t0      = tic; 
  Psub    = unique(cellfun(@(p) fileparts(fileparts(p)),Psnr,'UniformOutput',false));
  convmsg = cat_io_dcm2bids_import('convertSources',Psnr,C.Pmdbdir,C.job.opts,Psub);
  tconv   = toc(t0); 
  R       = cat_io_dcm2bids_process('processScans',C,Psnr,convmsg);
  t1      = tic; 
  sub     = unique(R.sub(~cellfun('isempty',R.sub) & ~cellfun('isempty',R.BIDSsub)));
  if render && ~isempty(sub)
    ropts = C.job.opts.render; ropts.subfilter = sub; ropts.dbdir = C.Pmdbdir;
    cat_io_dcm2bids_render(fullfile(C.Poutdir,C.job.BIDSdir,C.job.BIDSsubdir), ...
      fullfile(C.Poutdir,[C.job.BIDSdir '-report'],'render',C.job.BIDSsubdir), ropts);
  end
  if report && ~isempty(sub)
    try
      cat_io_dcm2bids_participantreport('writeReports',C.Poutdir,C.job,Psub);
    catch e
      cat_io_cprintf('err','  Participant reports failed: %s (%s, line %d)\n',e.message,e.stack(1).name,e.stack(1).line);
    end
  end
  t   = round([tconv 0 R.tpp toc(t1) toc(t0)]); t(2) = max(0,t(5) - t(1) - t(3) - t(4)); % (QC/export as rest, i.e. the sum fits)
  hms = @(t) char(duration(0,0,t,'Format','hh:mm:ss')); 
  fprintf('%s\n',repmat('-',1,154)); 
  fprintf('%154s\n',sprintf('DCM2NII %s + QC/export %s + PP %s + reports %s = %s',hms(t(1)),hms(t(2)),hms(t(3)),hms(t(4)),hms(t(5)))); 
end
% =========================================================================
function R = mergeResults(res)
%mergeResults. Results of all packages (see processScans in
%  cat_io_dcm2bids_process) in the order of the packages, where the
%  participant data of a subject in several packages keeps the data of
%  the first scan (minimum age, see processScans).
  R = res{1};
  FN = {'rec','varg','dbV','Pdbjson','datatype','pstatus','site','sub','ses','BIDSfile','BIDSsub','intended','bidsroots','convver','quarantine'};
  for ri = 2:numel(res)
    for fni = 1:numel(FN), R.(FN{fni}) = [R.(FN{fni}) res{ri}.(FN{fni})]; end
    for key = res{ri}.srcmap.keys, R.srcmap(key{1}) = res{ri}.srcmap(key{1}); end
    for key = res{ri}.partdata.keys
      v = res{ri}.partdata(key{1});
      if ~isKey(R.partdata,key{1}) || v.V.PatientAge < R.partdata(key{1}).V.PatientAge || ...
          (isnan(R.partdata(key{1}).V.PatientAge) && ~isnan(v.V.PatientAge))
        R.partdata(key{1}) = v;
      end
    end
  end
end

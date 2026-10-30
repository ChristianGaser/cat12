function varargout = cat_io_dcm2bids_parallel(action,varargin)
%cat_io_dcm2bids_parallel. Subject-wise (parallel) processing of cat_io_dcm2bids.
%  The scans of a run are processed in subject packages, i.e. all scans of
%  a BIDS subject (PatientID, for subject labels with site per site, see
%  packages) together, as the sessions of a subject depend on each other
%  (e.g. the T1w of another session for sessions without anatomical scan
%  and the sessions file of the subject). A package converts its images
%  (see convertSources in cat_io_dcm2bids_import), processes its scans (see
%  cat_io_dcm2bids_process), and renders the subject and session reports
%  of its subjects (with the slice cache for the rendering of all subjects,
%  see cat_io_dcm2bids_render).
%  For nproc > 0, nproc background MATLAB processes (workers, see
%  cat_parallelize) take the next package of a queue until all packages
%  are processed (i.e. long and short packages are balanced), whereas the
%  main process shows the log of each finished package.
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
%    R  = cat_io_dcm2bids_parallel(job)   (package of cat_parallelize)
%
%  See also cat_io_dcm2bids, cat_io_dcm2bids_process, cat_parallelize.

  if isstruct(action) % (package with job.data, job.C, and job.render)
    varargout{1} = processPackage(action.C,action.data,action.render);
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
%  the packages) and R.failed (number of failed packages).
  pk = packages(Psnr,C);
  nw = min(nproc,numel(pk));
  % no parallel processes without new packages (e.g. a new export of the
  % processed database), as the start of the workers takes longer
  if nw > 1 && ~any(newPackages(Psnr,pk,C))
    fprintf('DCM2BIDS - no new subjects, i.e. no parallel processes.\n'); nw = 0;
  end
  if nw < 2 % (one process)
    R = processPackage(C,Psnr,0); R.rendered = 0; R.failed = 0;
    return
  end

  % packages by the workers of cat_parallelize with the queue directory in
  % the database and the logs of the packages with colors
  render = C.job.opts.render.run && C.job.opts.output > 1 && C.job.opts.gzipi == C.job.opts.gzipe; % (otherwise the BIDS images are zipped again, see gzipunzipOutputdata)
  job    = struct('nproc',nw,'data',{Psnr},'C',C,'render',render);
  popts  = struct('packages',{pk},'replay',1,'mergefun',@mergeResults,'qroot',C.Pmdbdir, ...
    'qprefix','catDCM2BIDSpar_','logdir',fullfile(C.Ptbldir,'logs'),'logname','log_DCM2BIDS','name','DCM2BIDS');
  [R,status] = cat_parallelize(job,mfilename,'data',popts);
  failed = nnz(~cellfun('isempty',{status.err}));
  if isempty(R)
    error('cat_io_dcm2bids_parallel:failed','All %d packages failed.',numel(pk));
  end
  R.rendered = render; R.failed = failed;
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
%  several scanners for subject labels without site (opts.subIDform 1 and 
%  3, e.g. for the T1w of other sessions), and otherwise the same site
%  label and PatientID (opts.subIDform 2, 4, and 5, e.g. sub-SITE-PID, with
%  the site label of the scanner, see setupSites in cat_io_dcm2bids_protocols,
%  i.e. scanners with the same label of the site dictionary stay together), 
%  e.g. the scanners of a traveling subject are processed in parallel. 
%  Without context C of the run (job options and sites), the packages are
%  given by the PatientID.
  persite = exist('C','var') && isfield(C,'job') && isfield(C.job.opts,'subIDform') && ...
    ismember(C.job.opts.subIDform,[2 4 5]);
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
function new = newPackages(Psnr,pk,C)
%newPackages. Packages pk of the scan directories Psnr with new scans, i.e.
%  sessions that are not in the cumulative scan list of the database
%  (scans.mat, see writeTables in cat_io_dcm2bids_tables) or scans without
%  image (conversion, for the outputs 2 and 3, besides localizer and scouts
%  that are ignored, see processScans in cat_io_dcm2bids_process), and all
%  packages for opts.rerun (the test ends with the first new package, i.e.
%  the following packages stay new).
  new = true(1,numel(pk));
  if C.job.opts.rerun, return; end
  Pmat = fullfile(C.Ptbldir,'scans.mat');
  if exist(Pmat,'file'), S = load(Pmat,'R'); known = unique({S.R.Session}); else, return; end
  if C.job.opts.gzipi, niiext = '.nii.gz'; else, niiext = '.nii'; end
  for k = 1:numel(pk)
    new(k) = false;
    for di = pk{k}
      if ~ismember(fileparts(Psnr{di}(numel(C.Pmdbdir)+2:end)),known), new(k) = true; break; end
      if C.job.opts.output < 2, continue; end
      Pjson = cat_io_dcm2bids_db('getDBscanJSONs',Psnr{di});
      for ji = 1:numel(Pjson)
        parts = cat_io_dcm2bids_db('nameParts',spm_file(Pjson{ji},'basename'));
        if ~any(cat_io_contains(lower(parts),{'localizer','scout'})) && ~exist(spm_file(Pjson{ji},'ext',niiext),'file')
          new(k) = true; break
        end
      end
      if new(k), break; end
    end
    if new(k), return; end % (one new package is enough)
  end
end
% =========================================================================
function R = processPackage(C,Psnr,render)
%processPackage. Convert the images of the scan directories Psnr (only
%  into their subject directories, see convertSources in
%  cat_io_dcm2bids_import), process the scans (see cat_io_dcm2bids_process),
%  and render the subject and session reports of its subjects (render).
  Psub    = unique(cellfun(@(p) fileparts(fileparts(p)),Psnr,'UniformOutput',false));
  convmsg = cat_io_dcm2bids_import('convertSources',Psnr,C.Pmdbdir,C.job.opts,Psub);
  R       = cat_io_dcm2bids_process('processScans',C,Psnr,convmsg);
  sub     = unique(R.sub(~cellfun('isempty',R.sub) & ~cellfun('isempty',R.BIDSsub)));
  if render && ~isempty(sub)
    ropts = C.job.opts.render; ropts.subfilter = sub; ropts.dbdir = C.Pmdbdir;
    cat_io_dcm2bids_render(fullfile(C.Poutdir,C.job.BIDSdir,C.job.BIDSsubdir), ...
      fullfile(C.Poutdir,[C.job.BIDSdir '-report'],'render',C.job.BIDSsubdir), ropts);
  end
end
% =========================================================================
function R = mergeResults(res)
%mergeResults. Results of all packages (see processScans in
%  cat_io_dcm2bids_process) in the order of the packages, where the
%  participant data of a subject in several packages keeps the data of
%  the first scan (minimum age, see processScans).
  R = res{1};
  FN = {'rec','varg','dbV','Pdbjson','datatype','pstatus','site','sub','ses','BIDSfile','BIDSsub','intended','bidsroots','convver'};
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

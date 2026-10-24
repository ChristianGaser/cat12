function varargout = cat_io_dcm2bids_parallel(action,varargin)
%cat_io_dcm2bids_parallel. Subject-wise (parallel) processing of cat_io_dcm2bids.
%  The scans of a run are processed in subject packages, i.e. all scans of
%  a subject (PatientID, also of several scanners and sessions, see
%  packages) together, as the sessions of a subject depend on each other
%  (e.g. the T1w of another session for sessions without anatomical scan
%  and the sessions file of the subject). A package converts its images
%  (see convertSources in cat_io_dcm2bids_import), processes its scans (see
%  cat_io_dcm2bids_process), and renders the subject and session reports
%  of its subjects (with the slice cache for the rendering of all subjects,
%  see cat_io_dcm2bids_render).
%  For nproc > 0, nproc background MATLAB processes (workers, without the
%  Parallel Computing Toolbox) take the next package of a queue until all
%  packages are processed (i.e. long and short packages are balanced),
%  whereas the main process shows the output of each finished package.
%  Outputs of all subjects (tables of the database, participants, DICOM
%  export, report, rendering of all subjects) are written by cat_io_dcm2bids
%  with the merged results of all packages (in the order of the packages,
%  i.e. the same order as without parallel processing).
%  Packages of workers that did not start (e.g. without MATLAB license) are
%  processed by the main process, whereas packages of failed workers (e.g.
%  out of memory) are reported as errors (processed again by the next run).
%  The queue directory catDCM2BIDSpar_* of the database (context, packages,
%  results, and logs) is removed at the end, but kept after errors.
%
%  varargout = cat_io_dcm2bids_parallel(action,varargin)
%
%  Actions (see the help of the local functions):
%    R  = cat_io_dcm2bids_parallel('run',C,Psnr,nproc)
%    pk = cat_io_dcm2bids_parallel('packages',Psnr)
%         cat_io_dcm2bids_parallel('worker',Pq,wi)   (background process)
%
%  See also cat_io_dcm2bids, cat_io_dcm2bids_process.

  switch action
    case {'run','packages','worker'}
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
  pk = packages(Psnr);
  nw = min(nproc,numel(pk));
  if nw > 0 && ispc && contains(matlabroot,' ')
    cat_io_cprintf('warn','  No parallel processes, as the MATLAB directory contains spaces.\n'); nw = 0;
  end
  if nw < 2 % (one process)
    R = processPackage(C,Psnr,0); R.rendered = 0; R.failed = 0;
    return
  end

  % queue directory with the context of the run and the packages
  [~,tmpname] = fileparts(tempname);
  Pq = fullfile(C.Pmdbdir,['catDCM2BIDSpar_' tmpname]);
  D  = dir(fullfile(C.Pmdbdir,'catDCM2BIDSpar_*')); % (older than one day, e.g. after errors)
  for di = 1:numel(D)
    if D(di).isdir && datetime(D(di).datenum,'ConvertFrom','datenum') < datetime('now') - days(1)
      rmdir(fullfile(C.Pmdbdir,D(di).name),'s');
    end
  end
  mkdir(Pq);
  render   = C.job.opts.render.run && C.job.opts.output > 1 && C.job.opts.gzipi == C.job.opts.gzipe; % (otherwise the BIDS images are zipped again, see gzipunzipOutputdata)
  P        = cellfun(@(k) Psnr(k),pk,'UniformOutput',false);
  nthreads = max(1,floor(feature('numcores') / nw)); % computational threads per worker
  [spmdef,catdef] = getDefaults; % SPM and CAT defaults of the main process
  save(fullfile(Pq,'context.mat'),'C','P','render','nthreads','spmdef','catdef');

  % start the workers (background MATLAB processes)
  fprintf('DCM2BIDS - process %d subject packages by %d parallel processes (logs: %s)\n',numel(P),nw,Pq);
  for wi = 1:nw
    cmd = sprintf('addpath(''%s'',''%s''); cat_io_dcm2bids_parallel(''worker'',''%s'',%d);', ...
      spm('dir'),fileparts(mfilename('fullpath')),Pq,wi);
    Plog = fullfile(Pq,sprintf('w%02d.log',wi));
    if ispc
      system(sprintf('start "" /b "%s" -batch "%s" -logfile "%s"',fullfile(matlabroot,'bin','matlab.exe'),cmd,Plog));
    else
      system(sprintf('"%s" -batch "%s" > "%s" 2>&1 &',fullfile(matlabroot,'bin','matlab'),cmd,Plog));
    end
  end
  wclean = onCleanup(@() stopWorkers(Pq,nw)); % (e.g. after an interruption of the main process)

  % output of the finished packages (logs) until all workers are finished
  t0 = datetime('now'); res = cell(1,numel(P)); done = false(1,numel(P)); failed = 0;
  while ~all(done)
    pause(2);
    for k = find(~done)
      Pres = fullfile(Pq,sprintf('p%04d.mat',k)); Perr = fullfile(Pq,sprintf('p%04d_err.mat',k));
      if exist(Pres,'file') || exist(Perr,'file')
        printLog(fullfile(Pq,sprintf('p%04d.log',k)));
        if exist(Pres,'file')
          S = load(Pres,'R'); res{k} = S.R;
        else
          S = load(Perr,'E'); failed = failed + 1;
          cat_io_cprintf('err','  Package %d failed: %s\n',k,S.E.message);
        end
        done(k) = true;
        fprintf('DCM2BIDS - %d of %d packages processed (%s)\n',nnz(done),numel(P),char(duration(datetime('now') - t0)));
      end
    end
    if ~any(workersActive(Pq,nw,t0)), break; end
  end

  % packages that were not processed by the workers: not started packages
  % by the main process (e.g. workers without license), packages of failed
  % workers as errors
  for k = find(~done)
    if exist(fullfile(Pq,sprintf('p%04d.claim',k)),'dir')
      printLog(fullfile(Pq,sprintf('p%04d.log',k))); failed = failed + 1;
      cat_io_cprintf('err','  Package %d failed (the parallel process ended, e.g. out of memory, see %s).\n',k,Pq);
    else
      cat_io_cprintf('warn','  Package %d by the main process.\n',k);
      res{k} = processPackage(C,P{k},render);
    end
  end
  R = mergeResults(res(~cellfun('isempty',res)));
  R.rendered = render; R.failed = failed;
  clear wclean;
  if failed == 0
    rmdir(Pq,'s');
  else
    cat_io_cprintf('err','  %d package(s) failed (logs: %s), i.e. their scans are processed again by the next run.\n',failed,Pq);
  end
end
% =========================================================================
function pk = packages(Psnr)
%packages. Subject packages of the scan directories Psnr (sorted,
%  catDCM2BIDSdb/sub-<serial>-<PatientID>/ses-*/snr-*) as cell of indices,
%  i.e. all scan directories of the same PatientID (BIDS label, i.e. also
%  of several scanners, see the subject labels of opts.subIDform) in the
%  order of their first scan directory.
  Psub = cellfun(@(p) fileparts(fileparts(p)),Psnr,'UniformOutput',false);
  [usub,~,si] = unique(Psub,'stable');
  key = usub;
  for ui = 1:numel(usub)
    Pjson = cat_io_dcm2bids_db('getDBscanJSONs',Psnr{find(si == ui,1)});
    for ji = 1:numel(Pjson)
      try
        V = cat_io_json(Pjson{ji});
        if isfield(V,'PatientID') && ~isempty(V.PatientID)
          key{ui} = lower(strrep(cat_io_dcm2bids_bids('bidsLabel',char(V.PatientID)),'_','')); break
        end
      catch % (unreadable sidecar: the subject directory as package)
      end
    end
  end
  [~,~,ki] = unique(key(si),'stable');
  pk = arrayfun(@(k) find(ki == k)',1:max([0; ki(:)]),'UniformOutput',false);
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
    ropts = C.job.opts.render; ropts.subfilter = sub;
    cat_io_dcm2bids_render(fullfile(C.Poutdir,C.job.BIDSdir,C.job.BIDSsubdir), ...
      fullfile(C.Poutdir,[C.job.BIDSdir '-report'],'render',C.job.BIDSsubdir), ropts);
  end
end
% =========================================================================
function worker(Pq,wi)
%worker. Background process wi of the queue directory Pq: processes the
%  next package that is not processed by another worker (claimed by its
%  directory p####.claim) until all packages are processed, with the result
%  p####.mat (or the error p####_err.mat) and the log p####.log.
  if exist('matlabProcessID','file'), pid = matlabProcessID; else, fn = 'getpid'; pid = feature(fn); end % (older releases)
  fid = fopen(fullfile(Pq,sprintf('w%02d.pid',wi)),'w'); fprintf(fid,'%d',pid); fclose(fid);
  wdone = onCleanup(@() touch(fullfile(Pq,sprintf('w%02d.done',wi))));
  setPrint;
  S = load(fullfile(Pq,'context.mat'));
  evalc('spm(''defaults'',''fmri''); spm_jobman(''initcfg'');'); % (paths and batch system)
  setDefaults(S.spmdef,S.catdef);
  maxNumCompThreads(S.nthreads);
  for k = 1:numel(S.P)
    if ~claim(fullfile(Pq,sprintf('p%04d.claim',k))), continue; end
    diary(fullfile(Pq,sprintf('p%04d.log',k)));
    try
      R = processPackage(S.C,S.P{k},S.render);
      Ptmp = fullfile(Pq,sprintf('p%04d_tmp.mat',k));
      save(Ptmp,'R'); movefile(Ptmp,fullfile(Pq,sprintf('p%04d.mat',k))); % (complete file for the main process)
    catch e
      fprintf('\nERROR: %s\n',e.message);
      for si = 1:numel(e.stack), fprintf('  %s (line %d)\n',e.stack(si).name,e.stack(si).line); end
      E = struct('message',e.message,'stack',e.stack);
      save(fullfile(Pq,sprintf('p%04d_err.mat',k)),'E');
    end
    diary off
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
% =========================================================================
function act = workersActive(Pq,nw,t0)
%workersActive. Workers that are not finished: without done file and with
%  a running process (or without process ID within 5 minutes after the
%  start, e.g. the MATLAB startup with the license check).
  act = false(1,nw);
  for wi = 1:nw
    if exist(fullfile(Pq,sprintf('w%02d.done',wi)),'file'), continue; end
    pid = workerPID(Pq,wi);
    if isnan(pid)
      act(wi) = datetime('now') - t0 < minutes(5);
    else
      act(wi) = processRuns(pid);
    end
  end
end
% =========================================================================
function stopWorkers(Pq,nw)
%stopWorkers. Stop the workers that are still running (e.g. after an
%  interruption of the main process).
  if ~exist(Pq,'dir'), return; end
  for wi = 1:nw
    pid = workerPID(Pq,wi);
    if ~isnan(pid) && ~exist(fullfile(Pq,sprintf('w%02d.done',wi)),'file') && processRuns(pid)
      if ispc, system(sprintf('taskkill /F /PID %d',pid)); else, system(sprintf('kill %d',pid)); end
      cat_io_cprintf('warn','  Stopped parallel process %d (PID %d).\n',wi,pid);
    end
  end
end
% =========================================================================
function pid = workerPID(Pq,wi)
%workerPID. Process ID of the worker wi (NaN if not yet available).
  pid = nan; P = fullfile(Pq,sprintf('w%02d.pid',wi));
  if exist(P,'file'), pid = str2double(fileread(P)); end
end
% =========================================================================
function ok = processRuns(pid)
%processRuns. Test if the process pid runs.
  if ispc
    [~,out] = system(sprintf('tasklist /FI "PID eq %d" /NH',pid)); ok = contains(out,sprintf(' %d ',pid));
  else
    ok = system(sprintf('kill -0 %d 2>/dev/null',pid)) == 0;
  end
end
% =========================================================================
function ok = claim(P)
%claim. Claim a package by its directory (atomic creation, false if
%  another worker created it before).
  try
    ok = logical(java.io.File(P).mkdir());
  catch
    [ok,~,id] = mkdir(P); ok = ok && ~strcmp(id,'MATLAB:MKDIR:DirectoryExists');
  end
end
% =========================================================================
function printLog(P)
%printLog. Print the log of a package (the command line output of its
%  processing).
  if exist(P,'file'), fprintf('%s',fileread(P)); end
end
% =========================================================================
function touch(P)
%touch. Create an empty file.
  fid = fopen(P,'w'); if fid > 0, fclose(fid); end
end
% =========================================================================
function setPrint
%setPrint. Simple printing of cat_io_cprintf (without colors) in the log.
  global cprintferror %#ok<GVMIS>
  cprintferror = 1;
end
% =========================================================================
function [spmdef,catdef] = getDefaults
%getDefaults. SPM and CAT defaults of this process.
  global defaults cat %#ok<GVMIS>
  spmdef = defaults; catdef = cat;
end
% =========================================================================
function setDefaults(spmdef,catdef)
%setDefaults. SPM and CAT defaults of the main process.
  global defaults cat %#ok<GVMIS>
  defaults = spmdef; cat = catdef;
end

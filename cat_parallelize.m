function varargout = cat_parallelize(job,func,datafield,popts)
%cat_parallelize. Parallel processing of a function with job structure by
%  background MATLAB processes (workers, without the Parallel Computing
%  Toolbox). The items of the data field of the job (e.g. the files of
%  job.data) are split into packages that the workers take from a queue
%  until all packages are processed (i.e. long and short packages are
%  balanced). Each package is processed by func(job) with its items in
%  job.(datafield) and job.process_index (the worker mode of the function,
%  see the calling scheme below), where errors are limited to the package.
%  The outputs of the packages are merged in the order of the packages,
%  i.e. in the order of the items (e.g. for SPM dependencies).
%  Packages of workers that did not start (e.g. without MATLAB license) are
%  processed by the main process, whereas packages of failed workers (e.g.
%  out of memory) are reported as errors. An interruption of the main
%  process (Ctrl-C) stops the workers. The queue directory (context,
%  results, and logs of the packages) is removed at the end, but kept after
%  errors. The logs of the workers are written to the log directory (with
%  their links in the command window).
%
%  [out,status] = cat_parallelize(job,func[,datafield,popts])
%
%  job        .. job structure with the number of processes job.nproc and
%                the items in job.(datafield), i.e. a cellstr, a struct
%                array (e.g. job.subj), a cell of lists of the same length
%                (e.g. job.cdata), or a cell of groups of different length
%                (e.g. samples of job.data_surf), and further fields with
%                the same items in job.datafields (cellstr of fieldnames,
%                e.g. in cat_vol_savg)
%  func       .. function name
%  datafield  .. fieldname of the items (default 'data')
%  popts      .. options of the parallel processing (default):
%   .packages .. packages as cell of item indices ({}, see psize)
%   .psize    .. items per package (0 - automatic, i.e. one item per
%                package for up to 20 items per process and otherwise
%                10 packages per process, e.g. for many small surfaces)
%   .openlog  .. open the logs of the workers in the editor (0)
%   .logdir   .. directory of the logs of the workers ('' - the CAT BIDS
%                log directory of the items or the working directory)
%   .logname  .. name of the logs of the workers ('' - log_<func>, i.e.
%                <logname>_<date>_<worker>.txt)
%   .replay   .. print the log of each finished package with the colors of
%                cat_io_cprintf (0), i.e. the workers write markup of the
%                styles into their logs in the queue directory that are 
%                copied without markup to the logs of the workers
%   .printfun .. function handle that prints the line of a finished
%                package (instead of the default line) and returns a value
%                that is collected in status.print, called with a struct
%                with the fields k, n (packages), ndone (finished
%                packages), items, files, log, err, worker, and time (only
%                in the main process)
%   .mergefun .. function handle to merge the outputs of all packages
%                (cell in package order, default: concatenation of the
%                fields, see mergeOutput)
%   .qroot    .. directory of the queue directory ('' - tempdir)
%   .qprefix  .. name prefix of the queue directory ('catPar_')
%   .workdir  .. work directory of the packages (''), i.e. the input files
%                of a package (and their files with the same name, e.g.
%                BIDS sidecars) are copied into the work directory (with
%                their full path) and all new files are copied back after
%                the processing of the package (with the paths of the work
%                directory replaced in text and MAT files and the output),
%                e.g. for functions with many read and write operations of
%                data on a slow disk (e.g. an external hard disk):
%                ''     - none, i.e. the processing in the input directories
%                'auto' - the temporary directory of the system (tempdir) 
%                         if its disk is at least two times faster than the
%                         disk of the input (and has enough free space)
%                <dir>  - directory
%   .workspace.. disk space of a package in the work directory in GB (2)
%   .output   .. outputs of func (-1 - if requested, 0 - none, e.g. for
%                functions with outputs that are not required, 1 - yes)
%   .nthreads .. computational threads per worker (0 - computational
%                threads of this process / workers)
%   .name     .. name in the messages (func)
%  out        .. merged output of the packages (if requested)
%  status     .. struct array of the packages with the fields items, files,
%                err (error message), log, worker (0 - main process), and
%                print (output of printfun)
%
%  Calling scheme in func (the worker mode is given by job.process_index):
%
%   % split job and data into separate processes to save computation time
%   if isfield(job,'nproc') && job.nproc>0 && ~isfield(job,'process_index')
%     if nargout==1
%       varargout{1} = cat_parallelize(job,mfilename,'data_surf');
%     else
%       cat_parallelize(job,mfilename,'data_surf');
%     end
%     return
%   end
%
%  The workers are started as cat_parallelize('worker',Pq,wi) with the
%  queue directory Pq and the worker number wi.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  if ischar(job)
    if strcmp(job,'worker'), worker(func,datafield); return; end
    error('cat_parallelize:unknownAction','Unknown action "%s".',job);
  end
  if ~exist('datafield','var') || isempty(datafield), datafield = 'data'; end
  if ~exist('popts','var') || isempty(popts), popts = struct(); end
  popts = checkOptions(popts,func);

  % fields of the job that were given by the old version
  if ~isfield(job,'CATDir'),   job.CATDir   = fullfile(fileparts(mfilename('fullpath')),'CAT'); end
  if ~isfield(job,'fsavgDir'), job.fsavgDir = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces'); end
  if ~isfield(job,'lazy'),     job.lazy     = 0; end
  if ~isfield(job,'debug'),    job.debug    = cat_get_defaults('extopts.verb') > 2; end

  % items and packages
  if ischar(job.(datafield)), job.(datafield) = cellstr(job.(datafield)); end
  [n,mode] = itemMode(job.(datafield));
  files    = arrayfun(@(i) itemName(job.(datafield),mode,i),1:n,'UniformOutput',false);
  pk       = popts.packages;
  if isempty(pk)
    nw = max(1,min(job.nproc,n));
    if popts.psize > 0
      ps = popts.psize;
    elseif n <= 20 * nw
      ps = 1;
    else
      ps = ceil(n / (10 * nw));
    end
    pk = arrayfun(@(a) a:min(n,a + ps - 1),1:ps:n,'UniformOutput',false);
  end
  np = numel(pk);
  nw = min(job.nproc,np);
  if nw > 0 && ispc && (contains(matlabroot,' ') || contains(spm('dir'),' '))
    cat_io_cprintf('warn',['  No parallel processes, as the MATLAB or SPM directory contains spaces ' ...
      '(processing in this process).\n']);
    nw = 0;
  end
  if popts.output < 0, nout = nargout > 0 && hasOutput(func); else, nout = popts.output > 0; end
  status = struct('items',pk,'files',cellfun(@(k) files(k),pk,'UniformOutput',false), ...
    'err','','log','','worker',0,'print',[]);

  spm('FnBanner',popts.name,cat_version);

  % context of the workers in the queue directory
  ctx.job       = job;
  ctx.func      = func;
  ctx.datafield = datafield;
  ctx.mode      = mode;
  ctx.pk        = pk;
  ctx.nout      = nout;
  ctx.markup    = popts.replay;
  ctx.nthreads  = popts.nthreads;
  if ctx.nthreads == 0, ctx.nthreads = max(1,floor(maxNumCompThreads / max(1,nw))); end % (limit of this process)
  [ctx.spmdef,ctx.catdef] = getDefaults;
  ctx.mpath     = path;
  Pq = queueDir(popts);
  ctx.workroot  = workRoot(popts,files,max(1,nw));
  save(fullfile(Pq,'context.mat'),'-struct','ctx');
  res = cell(1,np); done = false(1,np);

  % start the workers (background MATLAB processes) with their logs in the
  % log directory (for replay with markup in the queue directory, see 
  % syncLogs)
  if nw > 0
    logdir = popts.logdir;
    if isempty(logdir), logdir = bidsLogDir(files); end
    if isempty(logdir), logdir = pwd; end
    if ~exist(logdir,'dir'), mkdir(logdir); end
    logdate = char(datetime('now','Format','yyyyMMdd-HHmmss'));
    Plog = arrayfun(@(wi) fullfile(logdir,sprintf('%s_%s_%02d.txt',popts.logname,logdate,wi)),1:nw,'UniformOutput',false);
    loff = zeros(1,nw); % (copied bytes of the logs with markup)
    fprintf('%s - process %d package(s) by %d parallel processes (queue: %s)\n',popts.name,np,nw,Pq);
    for wi = 1:nw
      if ctx.markup
        touch(Plog{wi}); startWorker(Pq,wi,fullfile(Pq,sprintf('w%02d.log',wi)));
      else
        startWorker(Pq,wi,Plog{wi});
      end
      fprintf('  Process %d: %s\n',wi,spm_file(Plog{wi},'link','edit(''%s'')'));
      pause(0.5 + rand(1)/10); % (starting many processes at once can cause MATLAB errors)
    end
    wclean = onCleanup(@() stopWorkers(Pq,nw,ctx.workroot)); % (e.g. after an interruption of the main process)
    opened = ~popts.openlog | ~canEdit; opened = repmat(opened,1,nw);

    % output of the finished packages until all workers are finished
    try, cat_progress_bar('Init',np,popts.name,'Packages Completed'); end %#ok<TRYNC>
    t0 = datetime('now');
    while ~all(done)
      pause(2);
      for wi = find(~opened)
        if exist(Plog{wi},'file'), edit(Plog{wi}); opened(wi) = true; end
      end
      [res,done,status] = collect(Pq,res,done,status,popts,t0);
      if ctx.markup, loff = syncLogs(Pq,Plog,loff); end
      if ~any(workersActive(Pq,nw,t0)), break; end
    end
    [res,done,status] = collect(Pq,res,done,status,popts,t0); % (packages that were finished at the end)
    for ti = 1:10 % (end of the workers after the last package)
      if all(done) && any(workersActive(Pq,nw,t0)), pause(1); end
    end
    if ctx.markup, syncLogs(Pq,Plog,loff); end
  end

  % packages that were not processed by the workers: not started packages
  % by the main process (e.g. workers without license or nw = 0), packages
  % of failed workers as errors
  for k = find(~done)
    Pclaim = fullfile(Pq,sprintf('p%04d.claim',k));
    if nw > 0 && exist(Pclaim,'dir')
      status(k).worker = claimWorker(Pclaim);
      status(k).log    = fullfile(Pq,sprintf('p%04d.log',k));
      status(k).err    = 'the parallel process ended (e.g. out of memory)';
      if popts.replay, printLog(status(k).log); end
      status(k).print = printPackage(k,status,popts,nan,nnz(done) + 1);
    elseif nw == 0 || claim(Pclaim)
      if nw > 0, cat_io_cprintf('warn','  Package %d by the main process.\n',k); end
      t1 = datetime('now');
      d0 = {get(0,'Diary'),get(0,'DiaryFile')}; % (diary of the user)
      status(k).log = fullfile(Pq,sprintf('p%04d.log',k)); diary(status(k).log);
      try
        res{k} = runPackage(ctx,k,Pq);
      catch e
        status(k).err = e.message;
        cat_io_cprintf('err','\nERROR: %s\n',e.message);
        for si = 1:numel(e.stack), fprintf('  %s (line %d)\n',e.stack(si).name,e.stack(si).line); end
      end
      diary off; if strcmp(d0{1},'on'), diary(d0{2}); end
      status(k).print = printPackage(k,status,popts,datetime('now') - t1,nnz(done) + 1);
    end
    done(k) = true;
  end
  if exist('wclean','var'), clear wclean; end

  % merged output in the order of the packages
  if nout
    res = res(~cellfun('isempty',res));
    if isempty(res)
      varargout{1} = [];
    elseif ~isempty(popts.mergefun)
      varargout{1} = popts.mergefun(res);
    else
      varargout{1} = res{1};
      for ri = 2:numel(res), varargout{1} = mergeOutput(varargout{1},res{ri}); end
    end
  end
  varargout{2} = status;
  try, cat_progress_bar('Clear'); end %#ok<TRYNC>

  failed = ~cellfun('isempty',{status.err});
  if any(failed)
    cat_io_cprintf('err','%s - %d of %d package(s) failed (logs: %s).\n',popts.name,nnz(failed),np,Pq);
  elseif exist(Pq,'dir')
    rmdir(Pq,'s');
  end
  fprintf('%-40s: %30s\n','Completed',spm('time'));
end
% =========================================================================
function popts = checkOptions(popts,func)
%checkOptions. Default options of the parallel processing (see help).
  def = struct('packages',{{}},'psize',0,'openlog',0,'logdir','','logname','', ...
    'replay',0,'printfun',[],'mergefun',[],'qroot','','qprefix','catPar_','workdir','','workspace',2, ...
    'output',-1,'nthreads',0,'name',func);
  FN = fieldnames(def);
  for fi = 1:numel(FN)
    if ~isfield(popts,FN{fi}), popts.(FN{fi}) = def.(FN{fi}); end
  end
  if isempty(popts.logname), popts.logname = ['log_' func]; end
  if isempty(popts.qroot),   popts.qroot   = tempdir; end
end
% =========================================================================
function Pq = queueDir(popts)
%queueDir. New queue directory in popts.qroot, where queue directories
%  that are older than one day are removed (e.g. after errors).
  D = dir(fullfile(popts.qroot,[popts.qprefix '*']));
  for di = 1:numel(D)
    if D(di).isdir && datetime(D(di).datenum,'ConvertFrom','datenum') < datetime('now') - days(1)
      try, rmdir(fullfile(popts.qroot,D(di).name),'s'); end %#ok<TRYNC>
    end
  end
  [~,tmpname] = fileparts(tempname);
  Pq = fullfile(popts.qroot,[popts.qprefix tmpname]);
  mkdir(Pq);
end
% =========================================================================
function [n,mode] = itemMode(x)
%itemMode. Number of items n of the data field x and its mode: 'lists'
%  for a cell of cellstr of the same length (e.g. job.cdata with the item
%  i in each list), 'groups' for a cell of cellstr of different length
%  (e.g. samples, with the items of all groups), and otherwise 'items'
%  (cellstr, struct array, or cell of other items, e.g. subjects).
  mode = 'items'; n = numel(x);
  if iscell(x) && ~isempty(x) && all(cellfun(@iscell,x(:))) && all(cellfun(@iscellstr,x(:)))
    L = cellfun(@numel,x);
    if all(L == L(1))
      mode = 'lists'; n = L(1);
    else
      mode = 'groups'; n = sum(L);
    end
  end
end
% =========================================================================
function x = subsetItems(x,mode,idx)
%subsetItems. Items idx of the data field x of the mode (see itemMode).
  switch mode
    case 'items'
      x = x(idx);
    case 'lists'
      for li = 1:numel(x), x{li} = x{li}(idx); end
    case 'groups'
      c = [0 cumsum(cellfun(@numel,x(:)'))];
      for gi = 1:numel(x), x{gi} = x{gi}(idx(idx > c(gi) & idx <= c(gi+1)) - c(gi)); end
  end
end
% =========================================================================
function name = itemName(x,mode,i)
%itemName. First filename (string) of the item i (for the messages).
  switch mode
    case 'items'
      if iscell(x), name = firstString(x{i}); else, name = firstString(x(i)); end
    case 'lists'
      name = firstString(x{1}{i});
    case 'groups'
      c  = [0 cumsum(cellfun(@numel,x(:)'))];
      gi = find(i > c,1,'last');
      name = firstString(x{gi}{i - c(gi)});
  end
end
% =========================================================================
function s = firstString(v)
%firstString. First string of a (nested) cell or struct.
  s = '';
  if ischar(v)
    if ~isempty(v), s = v(1,:); end
  elseif iscell(v)
    for i = 1:numel(v), s = firstString(v{i}); if ~isempty(s), return; end; end
  elseif isstruct(v)
    FN = fieldnames(v);
    for i = 1:numel(v)
      for fi = 1:numel(FN), s = firstString(v(i).(FN{fi})); if ~isempty(s), return; end; end
    end
  end
end
% =========================================================================
function job = packageJob(S,k)
%packageJob. Job of the package k of the context S with its items in the
%  data field (and the further fields job.datafields) and the worker mode
%  (job.process_index with the indices of the items).
  job = S.job; idx = S.pk{k};
  job.(S.datafield) = subsetItems(job.(S.datafield),S.mode,idx);
  if isfield(job,'datafields')
    for fi = 1:numel(job.datafields)
      if isfield(job,job.datafields{fi})
        [~,m] = itemMode(job.(job.datafields{fi}));
        job.(job.datafields{fi}) = subsetItems(job.(job.datafields{fi}),m,idx);
      end
    end
  end
  job.process_index = idx;
  job.printPID      = 1;
  job.verb          = 1;
end
% =========================================================================
function output = runPackage(S,k,Pq)
%runPackage. Process the package k of the context S (in the work directory
%  S.workroot, see workIn and workOut) with the output of func (or []).
  job = packageJob(S,k); output = []; W = [];
  if ~isempty(S.workroot)
    [~,qname] = fileparts(Pq);
    [job,W]   = workIn(job,S,fullfile(S.workroot,sprintf('catParWork_%s_p%04d',qname,k)));
  end
  try
    if S.nout, output = feval(S.func,job); else, feval(S.func,job); end
  catch e
    if ~isempty(W), workOut(W); end % (e.g. the logs of the error)
    rethrow(e);
  end
  if ~isempty(W), output = workOut(W,output); end
end
% =========================================================================
function root = workRoot(popts,files,nw)
%workRoot. Work directory of the packages (see popts.workdir), where work
%  directories of packages that are older than one day are removed (e.g.
%  of interrupted runs).
  root = '';
  if isempty(popts.workdir) || strcmp(popts.workdir,'native'), return; end
  if strcmp(popts.workdir,'auto')
    Pin = '';
    for fi = 1:numel(files)
      if exist(spm_file(files{fi},'number',''),'file'), Pin = fileparts(spm_file(files{fi},'number','')); break; end
    end
    if isempty(Pin), return; end
    Sin  = cat_io_dcm2bids_helper('diskSpeed',Pin);
    Sloc = cat_io_dcm2bids_helper('diskSpeed',tempdir);
    need = nw * popts.workspace + 5;
    if Sloc.MBps > 2 * Sin.MBps && Sloc.GB > need
      root = tempdir;
      fprintf('  Work directory: %s (%0.0f MB/s, input disk %0.0f MB/s)\n',root,Sloc.MBps,Sin.MBps);
    else
      fprintf('  Work directory: none (input disk %0.0f MB/s, system disk %0.0f MB/s with %0.0f GB free)\n', ...
        Sin.MBps,Sloc.MBps,Sloc.GB);
      return
    end
  else
    root = popts.workdir;
    if ~exist(root,'dir'), mkdir(root); end
    fprintf('  Work directory: %s\n',root);
  end
  D = dir(fullfile(root,'catParWork_*'));
  for di = 1:numel(D)
    if D(di).isdir && datetime(D(di).datenum,'ConvertFrom','datenum') < datetime('now') - days(1)
      try, rmdir(fullfile(root,D(di).name),'s'); end %#ok<TRYNC>
    end
  end
end
% =========================================================================
function [job,W] = workIn(job,S,root)
%workIn. Copy the input files of the package (existing files of the data
%  fields of the job, and their files with the same name, e.g. BIDS
%  sidecars) into the work directory root (with their full path), where
%  the data fields of the job get the copies. W contains the root and the
%  copied files (with date and size) to copy back only new files (see
%  workOut).
  W = struct('root',root,'in',{{}},'date',[],'bytes',[]);
  F = {S.datafield};
  if isfield(job,'datafields'), F = [F; job.datafields(:)]; end
  P = {};
  for fi = 1:numel(F)
    if isfield(job,F{fi}), P = [P; allStrings(job.(F{fi}))]; end %#ok<AGROW>
  end
  P = unique(spm_file(P,'number',''));
  P = P(cellfun(@(p) exist(p,'file') == 2,P));
  for pi = 1:numel(P)
    [pp,ff,ee] = fileparts(P{pi});
    if strcmp(ee,'.gz'), [~,ff] = fileparts(ff); end % (name without .nii.gz)
    D = dir(fullfile(pp,[ff '.*']));
    for di = 1:numel(D)
      if D(di).isdir, continue; end
      Psrc = fullfile(pp,D(di).name); Pdst = workPath(W.root,Psrc);
      if ~exist(fileparts(Pdst),'dir'), mkdir(fileparts(Pdst)); end
      copyfile(Psrc,Pdst);
      Dd = dir(Pdst); W.in{end+1} = Pdst; W.date(end+1) = Dd.datenum; W.bytes(end+1) = Dd.bytes;
    end
  end
  for fi = 1:numel(F)
    if isfield(job,F{fi})
      for pi = 1:numel(P), job.(F{fi}) = replaceStrings(job.(F{fi}),P{pi},workPath(W.root,P{pi})); end
    end
  end
end
% =========================================================================
function output = workOut(W,output)
%workOut. Copy the new and changed files of the work directory W.root back
%  (to their original path), where the paths of the work directory are
%  replaced in text files (e.g. XML reports and logs), MAT files, and the
%  output, and remove the work directory.
  if ~exist(W.root,'dir'), return; end
  roots = {W.root};
  if startsWith(W.root,'/var/'), roots{end+1} = ['/private' W.root]; end % (macOS)
  if startsWith(W.root,'/private/var/'), roots{end+1} = W.root(9:end); end
  D = dir(fullfile(W.root,'**','*')); D = D(~[D.isdir]);
  for di = 1:numel(D)
    Pw = fullfile(D(di).folder,D(di).name);
    ii = find(strcmp(W.in,Pw),1);
    if ~isempty(ii) && W.date(ii) == D(di).datenum && W.bytes(ii) == D(di).bytes, continue; end % (input)
    Po = Pw(numel(W.root)+1:end);
    if ispc, Po = [Po(2) ':' Po(3:end)]; end % (drive, see workPath)
    if ~exist(fileparts(Po),'dir'), mkdir(fileparts(Po)); end
    copyfile(Pw,Po);
    [~,~,ee] = fileparts(Po);
    try
      if any(strcmpi(ee,{'.xml','.txt','.json','.csv','.tsv','.html','.log','.gii'}))
        str = fileread(Po);
        if any(cellfun(@(r) contains(str,r),roots))
          for ri = 1:numel(roots), str = strrep(str,roots{ri},''); end
          fid = fopen(Po,'w'); fwrite(fid,str); fclose(fid);
        end
      elseif strcmpi(ee,'.mat')
        M = load(Po); M2 = M;
        for ri = 1:numel(roots), M2 = replaceStrings(M2,roots{ri},''); end
        if ~isequal(M,M2), save(Po,'-struct','M2'); end
      end
    catch e
      cat_io_cprintf('warn','  Cannot replace the work directory in "%s": %s\n',Po,e.message);
    end
  end
  if nargin > 1
    for ri = 1:numel(roots), output = replaceStrings(output,roots{ri},''); end
  end
  rmdir(W.root,'s');
end
% =========================================================================
function Pw = workPath(root,P)
%workPath. Path of the file P in the work directory root (full path of P
%  below root, without the colon of the drive on Windows).
  if ispc, P = strrep(P,':',''); end
  Pw = fullfile(root,P);
end
% =========================================================================
function S = allStrings(v)
%allStrings. All strings of a (nested) cell or struct.
  S = {};
  if ischar(v)
    S = cellstr(v);
  elseif iscell(v)
    for i = 1:numel(v), S = [S; allStrings(v{i})]; end %#ok<AGROW>
  elseif isstruct(v)
    FN = fieldnames(v);
    for i = 1:numel(v)
      for fi = 1:numel(FN), S = [S; allStrings(v(i).(FN{fi}))]; end %#ok<AGROW>
    end
  end
end
% =========================================================================
function v = replaceStrings(v,a,b)
%replaceStrings. Replace a by b in all strings of a (nested) cell or
%  struct.
  if ischar(v)
    if size(v,1) > 1, v = char(strrep(cellstr(v),a,b)); else, v = strrep(v,a,b); end
  elseif iscell(v)
    for i = 1:numel(v), v{i} = replaceStrings(v{i},a,b); end
  elseif isstruct(v)
    FN = fieldnames(v);
    for i = 1:numel(v)
      for fi = 1:numel(FN), v(i).(FN{fi}) = replaceStrings(v(i).(FN{fi}),a,b); end
    end
  end
end
% =========================================================================
function ok = hasOutput(func)
%hasOutput. Test if the function func has an output.
  try
    ok = nargout(func) ~= 0;
  catch
    ok = true;
  end
end
% =========================================================================
function startWorker(Pq,wi,Plog)
%startWorker. Start the background MATLAB process wi of the queue Pq with
%  the command line output in Plog.
  cmd = sprintf('addpath(''%s''); cat_parallelize(''worker'',''%s'',%d);',fileparts(mfilename('fullpath')),Pq,wi);
  if cat_io_matlabversion < 20191 % (without -batch)
    opt = '-nodesktop -nosplash -r'; cmd = ['try, ' cmd ' catch e, disp(e.message); end; exit;'];
  else
    opt = '-batch';
  end
  if ispc
    system(sprintf('start "" /b "%s" %s "%s" -logfile "%s"',fullfile(matlabroot,'bin','matlab.exe'),opt,cmd,Plog));
  else
    system(sprintf('"%s" %s "%s" > "%s" 2>&1 &',fullfile(matlabroot,'bin','matlab'),opt,cmd,Plog));
  end
end
% =========================================================================
function worker(Pq,wi)
%worker. Background process wi of the queue directory Pq: processes the
%  next package that is not processed by another worker (claimed by its
%  directory p####.claim) until all packages are processed, with the output
%  p####.mat (or the error p####_err.mat) and the log p####.log.
  if ischar(wi), wi = str2double(wi); end
  if exist('matlabProcessID','file'), pid = matlabProcessID; else, fn = 'getpid'; pid = feature(fn); end % (older releases)
  fid = fopen(fullfile(Pq,sprintf('w%02d.pid',wi)),'w'); fprintf(fid,'%d',pid); fclose(fid);
  wdone = onCleanup(@() touch(fullfile(Pq,sprintf('w%02d.done',wi))));
  S = load(fullfile(Pq,'context.mat'));
  path(S.mpath);
  global cprintferror %#ok<GVMIS>
  cprintferror = 1 + S.markup; % (simple printing, or markup of the styles for printLog)
  evalc('spm(''defaults'',''fmri''); spm_jobman(''initcfg'');'); % (batch system)
  setDefaults(S.spmdef,S.catdef);
  warning('off','MATLAB:maxNumCompThreads:Deprecated'); maxNumCompThreads(S.nthreads);
  fprintf('cat_parallelize: process %d (PID %d) of %s\n',wi,pid,S.func);
  for k = 1:numel(S.pk)
    Pclaim = fullfile(Pq,sprintf('p%04d.claim',k));
    if ~claim(Pclaim), continue; end
    touch(fullfile(Pclaim,sprintf('w%02d',wi)));
    diary(fullfile(Pq,sprintf('p%04d.log',k)));
    try
      output = runPackage(S,k,Pq);
      Ptmp = fullfile(Pq,sprintf('p%04d_tmp.mat',k));
      save(Ptmp,'output'); movefile(Ptmp,fullfile(Pq,sprintf('p%04d.mat',k))); % (complete file for the main process)
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
function [res,done,status] = collect(Pq,res,done,status,popts,t0)
%collect. Outputs of the packages that were finished by the workers and
%  their output in the command window.
  for k = find(~done)
    Pres = fullfile(Pq,sprintf('p%04d.mat',k)); Perr = fullfile(Pq,sprintf('p%04d_err.mat',k));
    if ~exist(Pres,'file') && ~exist(Perr,'file'), continue; end
    Pclaim = fullfile(Pq,sprintf('p%04d.claim',k));
    status(k).log    = fullfile(Pq,sprintf('p%04d.log',k));
    status(k).worker = claimWorker(Pclaim);
    if exist(Pres,'file')
      S = load(Pres,'output'); res{k} = S.output; Pend = Pres;
    else
      S = load(Perr,'E'); status(k).err = S.E.message; Pend = Perr;
    end
    D1 = dir(Pclaim); D2 = dir(Pend);
    if isempty(D1) || isempty(D2), dt = nan; else, dt = days(D2.datenum - min([D1.datenum])); end
    if popts.replay, printLog(status(k).log); end
    done(k) = true;
    status(k).print = printPackage(k,status,popts,dt,nnz(done));
    if popts.replay
      fprintf('%s - %d of %d packages processed (%s)\n',popts.name,nnz(done),numel(done),char(duration(datetime('now') - t0)));
    end
    try, cat_progress_bar('Set',nnz(done)); end %#ok<TRYNC>
  end
end
% =========================================================================
function off = syncLogs(Pq,Plog,off)
%syncLogs. Copy the new output of the workers with style markup (logs in
%  the queue directory, see printLog) without markup to their logs Plog,
%  where off contains the copied bytes (without incomplete markup at the
%  end).
  for wi = 1:numel(Plog)
    fid = fopen(fullfile(Pq,sprintf('w%02d.log',wi)),'r'); if fid < 0, continue; end
    fseek(fid,off(wi),'bof'); b = fread(fid,[1 inf],'*uint8'); fclose(fid);
    i17 = find(b == 17,1,'last');
    if ~isempty(i17) && ~any(b(i17:end) == 19), b = b(1:i17-1); end % (incomplete markup)
    if isempty(b), continue; end
    off(wi) = off(wi) + numel(b);
    str = regexprep(char(b),[char(17) '[^' char(18) ']*' char(18) '([^' char(19) ']*)' char(19)],'$1');
    fid = fopen(Plog{wi},'a'); if fid > 0, fwrite(fid,uint8(str)); fclose(fid); end
  end
end
% =========================================================================
function out = printPackage(k,status,popts,dt,ndone)
%printPackage. Output of the finished package k by popts.printfun or the
%  default line (and the error), where dt is the processing time and ndone
%  the number of finished packages.
  out = [];
  if ~isempty(popts.printfun)
    info = status(k); info.k = k; info.n = numel(status); info.ndone = ndone; info.time = dt;
    try
      out = popts.printfun(info); return
    catch e
      cat_io_cprintf('warn','  Output of package %d failed: %s\n',k,e.message);
    end
  end
  name = spm_str_manip(status(k).files{1},'k60');
  if numel(status(k).files) > 1, name = sprintf('%s (+%d)',name,numel(status(k).files) - 1); end
  if isempty(status(k).err)
    if popts.replay, return; end % (the log was printed)
    if isnan(dt), tstr = ''; else, dt.Format = 'hh:mm:ss'; tstr = sprintf(' (%s)',char(dt)); end
    cat_io_cprintf([0 0.5 0],'  %d/%d finished%s: %s\n',k,numel(status),tstr,name);
  else
    cat_io_cprintf('err','  %d/%d failed: %s\n    %s\n',k,numel(status),name,status(k).err);
    if ~isempty(status(k).log) && exist(status(k).log,'file')
      fprintf('    Log: %s\n',spm_file(status(k).log,'link','edit(''%s'')'));
    end
  end
end
% =========================================================================
function a = mergeOutput(a,b)
%mergeOutput. Merge the output b of a package to the outputs a of the
%  previous packages: fields of structures (also of struct arrays, e.g.
%  out.sample(i).files), concatenation of arrays and cells (as row if one
%  of them is a row and otherwise as column, e.g. for packages with one
%  file, where empty entries of cells are removed, e.g. of failed items),
%  and cellstr of different strings.
  if isempty(a), a = b; return; end
  if isempty(b), return; end
  if isstruct(a) && isstruct(b)
    FN = fieldnames(b);
    for j = 1:numel(b)
      for fi = 1:numel(FN)
        if j > numel(a) || ~isfield(a,FN{fi})
          a(j).(FN{fi}) = b(j).(FN{fi});
        else
          a(j).(FN{fi}) = mergeOutput(a(j).(FN{fi}),b(j).(FN{fi}));
        end
      end
    end
  elseif ischar(a) && ischar(b)
    if ~strcmp(a,b), a = [cellstr(a); cellstr(b)]; end
  elseif (iscell(a) && iscell(b)) || ((isnumeric(a) || islogical(a)) && (isnumeric(b) || islogical(b)))
    if isvector(a) && isvector(b)
      if size(a,2) > 1 || size(b,2) > 1, a = [a(:)' b(:)']; else, a = [a(:); b(:)]; end
    else
      try, a = [a; b]; catch, a = [a b]; end
    end
    if iscell(a), a(cellfun('isempty',a)) = []; end
  else
    try, a = [a; b]; end %#ok<TRYNC>
  end
end
% =========================================================================
function logdir = bidsLogDir(files)
%bidsLogDir. CAT BIDS log directory of the files (derivatives/CAT*/log)
%  or '' (the working directory).
  logdir = '';
  try
    BIDS = cat_io_BIDS(files,struct());
    if any([BIDS.isBIDS]), logdir = cat_io_BIDS(BIDS(1),'logdir'); end
  catch
  end
end
% =========================================================================
function ok = canEdit
%canEdit. Test if logs can be opened in the editor.
  ok = ~strcmpi(spm_check_version,'octave') && usejava('jvm') && feature('ShowFigureWindows') && usejava('awt');
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
function stopWorkers(Pq,nw,workroot)
%stopWorkers. Stop the workers that are still running (e.g. after an
%  interruption of the main process) and remove their work directories.
  if ~exist(Pq,'dir'), return; end
  if ~isempty(workroot)
    [~,qname] = fileparts(Pq);
  end
  for wi = 1:nw
    pid = workerPID(Pq,wi);
    if ~isnan(pid) && ~exist(fullfile(Pq,sprintf('w%02d.done',wi)),'file') && processRuns(pid)
      if ispc, system(sprintf('taskkill /F /PID %d',pid)); else, system(sprintf('kill %d',pid)); end
      cat_io_cprintf('warn','  Stopped parallel process %d (PID %d).\n',wi,pid);
    end
  end
  if ~isempty(workroot)
    D = dir(fullfile(workroot,['catParWork_' qname '_p*'])); % (packages of stopped or failed workers)
    for di = 1:numel(D), try, rmdir(fullfile(D(di).folder,D(di).name),'s'); end; end %#ok<TRYNC>
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
function wi = claimWorker(Pclaim)
%claimWorker. Worker of a claimed package (0 - unknown).
  wi = 0; D = dir(fullfile(Pclaim,'w*'));
  if ~isempty(D), wi = str2double(D(1).name(2:end)); end
end
% =========================================================================
function printLog(P)
%printLog. Print the log of a package (the command line output of its
%  processing) with the colors of cat_io_cprintf (style markup of the
%  workers, see cat_io_cprintf with cprintferror = 2).
  if ~exist(P,'file'), return; end
  str = fileread(P);
  [tok,txt] = regexp(str,[char(17) '([^' char(18) ']*)' char(18) '([^' char(19) ']*)' char(19)],'tokens','split');
  for ti = 1:numel(txt)
    fprintf('%s',txt{ti});
    if ti <= numel(tok)
      style = tok{ti}{1};
      if ~isempty(style) && style(1) == '[', style = str2num(style); end %#ok<ST2NM> % RGB, e.g. [0 0 0.5]
      cat_io_cprintf(style,'%s',tok{ti}{2});
    end
  end
end
% =========================================================================
function touch(P)
%touch. Create an empty file.
  fid = fopen(P,'w'); if fid > 0, fclose(fid); end
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

function varargout = cat_run(job)
% Segment a bunch of images
% ______________________________________________________________________
%
%   FORMAT cat_run(job)
%
%   job.channel(n).vols{m}
%   job.channel(n).biasreg
%   job.channel(n).biasfwhm
%   job.channel(n).write
%   job.tissue(k).tpm
%   job.tissue(k).ngaus
%   job.tissue(k).native
%   job.tissue(k).warped
%
% See the user interface for a description of the fields.
%
% based on John Ashburners version of
% spm_preproc8_run.m 2281 2008-10-01 12:52:50Z john $
% ______________________________________________________________________
% Christian Gaser
% $Id$

%#ok<*AGROW>

%rev = '$Rev$';


% split job and data into separate processes to save computation time
if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index'))  
  cat_io_cprintf('warn',...
    ['\nWARNING: Please note that no additional modules in the batch can be run \n' ...
     '         except CAT12 segmentation. Any dependencies will be broken for \n' ...
     '         subsequent modules if you split the job into separate processes.\n\n']);
    
  % rescue original subjects
  job_data = job.data;
  n_subjects = numel(job.data);
  if job.nproc > n_subjects
    job.nproc = n_subjects;
  end
  job.process_index = cell(job.nproc,1);

  % initial splitting of data
  for i=1:job.nproc
    job.process_index{i} = (1:job.nproc:(n_subjects-job.nproc+1))+(i-1);
  end

  % check if all data are covered
  for i=1:rem(n_subjects,job.nproc)
    job.process_index{i} = [job.process_index{i} n_subjects-i+1];
  end

  tmp_array = cell(job.nproc,1); job.printPID = 1; 
    
  logdate   = datestr(now,'YYYYmmdd_HHMMSS');
  for i=1:job.nproc
    fprintf('Running job %d:\n',i);
    for fi=1:numel(job_data(job.process_index{i}))
      fprintf('  %s\n',spm_str_manip(char(job_data(job.process_index{i}(fi))),'a78')); 
    end
    job.data = job_data(job.process_index{i});
         
    % temporary name for saving job information
    tmp_name = [tempname '.mat'];
    tmp_array{i} = tmp_name; 
    def = cat_get_defaults; job = cat_io_checkinopt(job,def); % further job update required here to get the latest cat defaults
    global defaults cat12; %#ok<NUSED,TLEV>
    save(tmp_name,'job','defaults','cat12');
    clear defaults cat12;
    
    % matlab command, cprintferror=1 for simple printing         
    matlab_cmd = sprintf('"global cprintferror; cprintferror=1; addpath %s %s %s %s; load %s; cat_run(job); "',...
      spm('dir'),fullfile(spm('dir'),'toolbox','cat12'),...
        fullfile(spm('dir'),'toolbox','OldNorm'),fullfile(spm('dir'),'toolbox','DARTEL'), tmp_name);

    % log-file for output
    log_name = ['catlog_main_' logdate '_log' sprintf('%02d',i) '.txt'];

    % call matlab with command in the background
    if ispc
      % check for spaces in filenames that will not work with windows systems and background jobs
      if strfind(spm('dir'),' ')
        cat_io_cprintf('warn',...
            ['\nWARNING: No background processes possible because your SPM installation is located in \n' ...
             '         a folder that contains spaces. Please set the number of processes in the GUI \n'...
             '         to ''0''. In order to split your job into different processes,\n' ...
             '         please do not use any spaces in folder names!\n\n']);
         job.nproc = 0; 
         job = update_job(job);
         varargout{1} = run_job(job);
         return; 
      end
      % prepare system specific path for matlab
      export_cmd = ['set PATH=' fullfile(matlabroot,'bin')];
      [status,result] = system(export_cmd);
      system_cmd = ['start matlab -nodesktop -nosplash -r ' matlab_cmd ' -logfile ' log_name];
    else
      % -nodisplay .. nodisplay is without figure output > problem with CAT report ... was there a server problem with -nodesktop?
      system_cmd = [fullfile(matlabroot,'bin') '/matlab -nodesktop -nosplash -r ' matlab_cmd ' -logfile ' log_name ' 2>&1 & '];
    end
    [status,result] = system(system_cmd); 
    cat_check_system_output(status,result);
    
    test = 0; lim = 10; ptime = 0.5;
    while test<lim
      if ~exist(log_name,'file')
        pause(ptime); 
        test = test + ptime; 
        if test>=lim
          cat_io_cprintf('warn','"%s" not exist after %d seconds! Proceed! \n',log_name,lim)
        end
      else 
        test = inf; 
        edit(log_name);
      end
    end

    edit(log_name);
    fprintf('\nCheck %s for logging information.\n',spm_file(log_name,'link','edit(''%s'')'));
    fprintf('_______________________________________________________________\n');

  end

  job = update_job(job);
  varargout{1} = vout_job(job);
  return
end

if isfield(job,'printPID') && job.printPID 
  cat_display_matlab_PID
end
job = update_job(job);

varargout{1} = run_job(job);

return
%_______________________________________________________________________
function job = update_job(job)

  % get defaults
  def = cat_get_defaults;
  if isfield(job.extopts,'restypes')
    def.extopts.restype = (char(fieldnames(job.extopts.restypes))); 
    def.extopts.resval  = job.extopts.restypes.(def.extopts.restype);
  end
  def.opts.fwhm = 1;
  def.nproc     = 0; 
  job = cat_io_checkinopt(job,def);
  if ~isfield(job.extopts,'restypes')
    job.extopts.restypes.(def.extopts.restype) = job.extopts.resval;  
  end
  
  % check range of str variables
  FN = {'NCstr','WMHCstr','LASstr','BVCstr','gcutstr','cleanupstr','mrf'};
  for fni=1:numel(FN)
    if ~isfield(job.extopts,FN{fni})  
      job.extopts.(FN{fni}) = max(0,min(1,job.extopts.(FN{fni})));
    end
  end

  % deselect ROI output and print warning if dartel template was changed
  [pth,nam] = spm_fileparts(job.extopts.darteltpm{1});
  if ~strcmp(nam,'Template_1_IXI555_MNI152')
    warning('DARTEL:template:change',...
      'Dartel template was changed: Please be aware that ROI analysis and other template-specific options cannot be used.');
    job.output.ROI = 0;
  end
  
  % set cat12.bb and vb.vox by Dartel template properties
  Vd       = spm_vol([job.extopts.darteltpm{1} ',1']);
  [bb,vox] = spm_get_bbox(Vd, 'old');  
  if job.extopts.bb(1)>job.extopts.bb(2), bbt=job.extopts.bb(1); job.extopts.bb(1)=job.extopts.bb(2); job.extopts.bb(2)=bbt; clear bbt; end
  if bb(1)>bb(2), bbt=bb(1); bb(1)=bb(2); bb(2)=bbt; clear bbt; end
  job.extopts.bb  = [ max(bb(1,1:3) , bb(1,1:3) ./ ((isinf(bb(1,1:3)) | isnan(bb(1,1:3)))+eps))
                      min(bb(2,1:3) , bb(2,1:3) ./ ((isinf(bb(2,1:3)) | isnan(bb(2,1:3)))+eps)) ];
          
  if isinf(job.extopts.vox) || isnan(job.extopts.vox)
    job.extopts.vox = abs(vox);
  end
  
  
  % prepare tissue priors and number of gaussians for all 6 classes
  [pth,nam,ext] = spm_fileparts(job.opts.tpm{1});
  clsn = numel(spm_vol(fullfile(pth,[nam ext]))); 
  tissue = struct();
  for i=1:clsn;
    tissue(i).ngaus = job.opts.ngaus(i);
    tissue(i).tpm = [fullfile(pth,[nam ext]) ',' num2str(i)];
  end
  
  tissue(1).warped = [job.output.GM.warped  (job.output.GM.mod==1)        (job.output.GM.mod==2)       ];
  tissue(1).native = [job.output.GM.native  (job.output.GM.dartel==1)     (job.output.GM.dartel==2)    ];
  tissue(2).warped = [job.output.WM.warped  (job.output.WM.mod==1)        (job.output.WM.mod==2)       ];
  tissue(2).native = [job.output.WM.native  (job.output.WM.dartel==1)     (job.output.WM.dartel==2)    ];
  tissue(3).warped = [job.output.CSF.warped (job.output.CSF.mod==1)       (job.output.CSF.mod==2)      ];
  tissue(3).native = [job.output.CSF.native (job.output.CSF.dartel==1)    (job.output.CSF.dartel==2)   ];

  % never write class 4-6
  for i=4:6;
    tissue(i).warped = [0 0 0];
    tissue(i).native = [0 0 0];
  end

  job.channel  = struct('vols',{job.data});
  job.tissue   = tissue;
return;

%_______________________________________________________________________
function vout = run_job(job)
  vout   = vout_job(job);

  % load tpm priors 
  tpm = char(cat(1,job.tissue(:).tpm));
  tpm = spm_load_priors8(tpm);

  for subj=1:numel(job.channel(1).vols),
    % __________________________________________________________________
    % Separation for old and new try-catch blocks of matlab. The new
    % try-catch block has to be in a separate file to avoid an error.
    % Both functions finally call cat_run_job.
    % See also cat_run_newcatch and cat_run_newcatch.
    % __________________________________________________________________
    %if job.extopts.ignoreErrors
      if cat_io_matlabversion>20072 
        cat_run_newcatch(job,tpm,subj); 
      else
        % inactive because of unclear error messages
        %cat_run_oldcatch(job,tpm,subj);
        cat_run_job(job,tpm,subj);
      end
    %else
    %  cat_run_job(job,tpm,subj);
    %end
  end

  colormap(gray)
  
  if isfield(job,'nproc') && job.nproc>0 
    fprintf('\n%s',repmat('_',1,72));
    fprintf('\nCAT12 Segmentation job finished.\n');
  end
return
%_______________________________________________________________________

function vout = vout_job(job)
% ----------------------------------------------------------------------
% create output structure for SPM batch mode
% ----------------------------------------------------------------------

n     = numel(job.channel(1).vols);
parts = cell(n,4);

biascorr  = {};
wbiascorr = {};
label     = {};
wlabel    = {};
rlabel    = {};
alabel    = {};

if job.extopts.subfolders
  mrifolder = 'mri';
else
  mrifolder = '';
end

for j=1:n,
    [parts{j,:}] = spm_fileparts(job.channel(1).vols{j});
end


% bias
% ----------------------------------------------------------------------
if job.output.bias.native,
    biascorr = cell(n,1);
    for j=1:n
        biascorr{j} = fullfile(parts{j,1},mrifolder,['m',parts{j,2},'.nii']);
    end
end

if job.output.bias.warped,
    wbiascorr = cell(n,1);
    for j=1:n
        wbiascorr{j} = fullfile(parts{j,1},mrifolder,['wm',parts{j,2},'.nii']);
    end
end

if job.output.bias.dartel==1,
    rbiascorr = cell(n,1);
    for j=1:n
        rbiascorr{j} = fullfile(parts{j,1},mrifolder,['rm',parts{j,2},'.nii']);
    end
end

if job.output.bias.dartel==2,
    abiascorr = cell(n,1);
    for j=1:n
        abiascorr{j} = fullfile(parts{j,1},mrifolder,['rm',parts{j,2},'_affine.nii']);
    end
end


% label
% ----------------------------------------------------------------------
if job.output.label.native,
    label = cell(n,1);
    for j=1:n
        label{j} = fullfile(parts{j,1},mrifolder,['p0',parts{j,2},'.nii']);
    end
end

if job.output.label.warped,
    wlabel = cell(n,1);
    for j=1:n
        wlabel{j} = fullfile(parts{j,1},mrifolder,['wp0',parts{j,2},'.nii']);
    end
end

if job.output.label.dartel==1,
    rlabel = cell(n,1);
    for j=1:n
        rlabel{j} = fullfile(parts{j,1},mrifolder,['rp0',parts{j,2},'.nii']);
    end
end

if job.output.label.dartel==2,
    alabel = cell(n,1);
    for j=1:n
        alabel{j} = fullfile(parts{j,1},mrifolder,['rp0',parts{j,2},'_affine.nii']);
    end
end


% ----------------------------------------------------------------------
param = cell(n,1);
for j=1:n
    param{j} = fullfile(parts{j,1},['cat12_',parts{j,2},'.mat']);
end


% tissues
% ----------------------------------------------------------------------
tiss = struct('p',{},'rp',{},'rpa',{},'wp',{},'mwp',{},'m0wp',{});
for i=1:numel(job.tissue),
    if job.tissue(i).native(1),
        tiss(i).p = cell(n,1);
        for j=1:n
            tiss(i).p{j} = fullfile(parts{j,1},mrifolder,['p',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).native(2),
        tiss(i).rp = cell(n,1);
        for j=1:n
            tiss(i).rp{j} = fullfile(parts{j,1},mrifolder,['rp',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).native(3),
        tiss(i).rpa = cell(n,1);
        for j=1:n
            tiss(i).rpa{j} = fullfile(parts{j,1},mrifolder,['rp',num2str(i),parts{j,2},'_affine.nii']);
        end
    end
    if job.tissue(i).warped(1),
        tiss(i).wp = cell(n,1);
        for j=1:n
            tiss(i).wp{j} = fullfile(parts{j,1},mrifolder,['wp',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).warped(2),
        tiss(i).mwp = cell(n,1);
        for j=1:n
            tiss(i).mwp{j} = fullfile(parts{j,1},mrifolder,['mwp',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).warped(3),
        tiss(i).m0wp = cell(n,1);
        for j=1:n
            tiss(i).m0wp{j} = fullfile(parts{j,1},mrifolder,['m0wp',num2str(i),parts{j,2},'.nii']);
        end
    end
end


% warping fields
% ----------------------------------------------------------------------
if job.output.warps(1),
    fordef = cell(n,1);
    for j=1:n
        fordef{j} = fullfile(parts{j,1},mrifolder,['y_',parts{j,2},'.nii']);
    end
else
    fordef = {};
end

if job.output.warps(2),
    invdef = cell(n,1);
    for j=1:n
        invdef{j} = fullfile(parts{j,1},mrifolder,['iy_',parts{j,2},'.nii']);
    end
else
    invdef = {};
end


% jacobian
% ----------------------------------------------------------------------
if job.output.jacobian.warped,
    jacobian = cell(n,1);
    for j=1:n
        jacobian{j} = '';
    end
else
    jacobian = {};
end


% ----------------------------------------------------------------------
vout  = struct('tiss',tiss,'label',{label},'wlabel',{wlabel},'rlabel',{rlabel},'alabel',{alabel},...
               'biascorr',{biascorr},'wbiascorr',{wbiascorr},'param',{param},...
               'invdef',{invdef},'fordef',{fordef},'jacobian',{jacobian});
%_______________________________________________________________________

%=======================================================================

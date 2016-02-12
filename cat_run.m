function varargout = cat_run(job,arg)
% Segment a bunch of images
% FORMAT cat_run(job)
% job.channel(n).vols{m}
% job.channel(n).biasreg
% job.channel(n).biasfwhm
% job.channel(n).write
% job.tissue(k).tpm
% job.tissue(k).ngaus
% job.tissue(k).native
% job.tissue(k).warped
% job.cat.affreg
% job.cat.reg
% job.cat.samp
% job.cat.warps
% job.cat.darteltpm
% job.cat.print
%
% See the user interface for a description of the fields.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% based on John Ashburners version of
% spm_preproc8_run.m 2281 2008-10-01 12:52:50Z john $
%
% Christian Gaser
% $Id$
%
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

  tmp_array = cell(job.nproc,1);
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
    global defaults cat12; %#ok<NUSED,TLEV>
    save(tmp_name,'job','defaults','cat12');
    clear defaults cat12;
    
    % matlab command          
    matlab_cmd = sprintf('"addpath %s %s %s %s;load %s; cat_run(job); "',spm('dir'),fullfile(spm('dir'),'toolbox','cat12'),...
        fullfile(spm('dir'),'toolbox','OldNorm'),fullfile(spm('dir'),'toolbox','DARTEL'), tmp_name);

    % log-file for output
    log_name = ['catlog_main_' logdate '_log' sprintf('%02d',i) '.txt'];

    % call matlab with command in the background
    if ispc
      % prepare system specific path for matlab
      export_cmd = ['set PATH=' fullfile(matlabroot,'bin')];
      system_cmd = [export_cmd ' & start matlab.bat -nodesktop -nosplash -r ' matlab_cmd ' -logfile ' log_name];
    else
      % -nodisplay .. nodisplay is without figure output > problem with CAT report ... was there a server problem with -nodesktop?
      system_cmd = [fullfile(matlabroot,'bin') '/matlab -nodesktop -nosplash -r ' matlab_cmd ' -logfile ' log_name ' 2>&1 & '];
    end

    [status,result] = system(system_cmd);
    
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

job = update_job(job);

varargout{1} = run_job(job);

return
%_______________________________________________________________________
function job = update_job(job)

  cat12 = struct('species', cat_get_defaults('extopts.species'), ... job.extopts.species,...
             'cat12atlas',cat_get_defaults('extopts.cat12atlas'), ... 
             'darteltpm', job.extopts.darteltpm{1}, ...
             'brainmask', cat_get_defaults('extopts.brainmask'), ...
             'affreg',    job.opts.affreg,...
             'samp',      cat_get_defaults('opts.samp'),...
             'warps',     job.output.warps,...
             'sanlm',     cat_get_defaults('extopts.sanlm'),... % job.extopts.sanlm,...
             'print',     job.extopts.print,...
             'subfolders',cat_get_defaults('extopts.subfolders'),...
             'ngaus',     cat_get_defaults('opts.ngaus'),...
             'reg',       cat_get_defaults('opts.warpreg'),...
             'bb',        cat_get_defaults('extopts.bb'));

  if isfield(job.extopts,'restype')
    cat12.restype = char(fieldnames(job.extopts.restype));
    cat12.resval  = job.extopts.restype.(cat12.restype); 
  else
    cat12.restype = cat_get_defaults('extopts.restype');
    cat12.resval  = cat_get_defaults('extopts.resval');
  end
  if isfield(job.extopts,'sanlm')
    cat12.sanlm = job.extopts.sanlm;
  end
  if ~isfield(cat12,'vox')
    cat12.vox = cat_get_defaults('extopts.vox');
  end
  if ~isfield(job.extopts,'verb')
    job.extopts.verb =  cat_get_defaults('extopts.verb');
  end
  if ~isfield(job.extopts,'APP')
    job.extopts.APP =  cat_get_defaults('extopts.APP');
  end
  if ~isfield(job.extopts,'pbtres')
    job.extopts.pbtres = cat_get_defaults('extopts.pbtres');
  end
  if ~isfield(job.output,'ROI')
    job.output.ROI =  cat_get_defaults('output.ROI');
  end
           
  if ~isfield(job.output,'CSF')
    job.output.CSF =  struct('modulated',cat_get_defaults('output.CSF.mod'),'dartel',cat_get_defaults('output.CSF.dartel'),...
                             'warped',cat_get_defaults('output.CSF.warped'),'native',cat_get_defaults('output.CSF.native'));
  end

  if ~isfield(job.output,'label')
    job.output.label =  cat_get_defaults('output.label');
  end

  % deselect ROI output and print warning if dartel template was changed
  [pth,nam,ext] = spm_fileparts(cat12.darteltpm);
  if ~strcmp(nam,'Template_1_IXI555_MNI152')
    warning('Dartel template was changed: Please be aware that ROI analysis and other template-specific options cannot be used.');
    job.output.ROI = 0;
  end
  
  % set cat12.bb and vb.vox by Dartel template properties
  Vd       = spm_vol([cat12.darteltpm ',1']);
  [bb,vox] = spm_get_bbox(Vd, 'old');  
  if cat12.bb(1)>cat12.bb(2), bbt=cat12.bb(1); cat12.bb(1)=cat12.bb(2); cat12.bb(2)=bbt; clear bbt; end
  if bb(1)>bb(2), bbt=bb(1); bb(1)=bb(2); bb(2)=bbt; clear bbt; end
  cat12.bb  = [ max(bb(1,1:3) , bb(1,1:3) ./ ((isinf(bb(1,1:3)) | isnan(bb(1,1:3)))+eps))
                min(bb(2,1:3) , bb(2,1:3) ./ ((isinf(bb(2,1:3)) | isnan(bb(2,1:3)))+eps)) ];
          
  if isinf(cat12.vox) || isnan(cat12.vox)
    cat12.vox = abs(vox);
  end

  % prepare tissue priors and number of gaussians for all 6 classes
  [pth,nam,ext] = spm_fileparts(job.opts.tpm{1});
  clsn = numel(spm_vol(fullfile(pth,[nam ext]))); 
  tissue = struct();
  for i=1:clsn;
    tissue(i).ngaus = cat12.ngaus(i);
    tissue(i).tpm = [fullfile(pth,[nam ext]) ',' num2str(i)];
  end

  % check whether native field is defined (only defined for expert mode)              
  if ~isfield(job.output.GM,'warped')
    job.output.bias.dartel = 0;
    job.output.bias.native = 0;
    job.output.GM.warped   = 0;
    job.output.WM.warped   = 0; 
  end
  
  tissue(1).warped = [job.output.GM.warped  (job.output.GM.modulated==1)  (job.output.GM.modulated==2) ];
  tissue(1).native = [job.output.GM.native  (job.output.GM.dartel==1)     (job.output.GM.dartel==2)    ];
  tissue(2).warped = [job.output.WM.warped  (job.output.WM.modulated==1)  (job.output.WM.modulated==2) ];
  tissue(2).native = [job.output.WM.native  (job.output.WM.dartel==1)     (job.output.WM.dartel==2)    ];
  tissue(3).warped = [job.output.CSF.warped (job.output.CSF.modulated==1) (job.output.CSF.modulated==2)];
  tissue(3).native = [job.output.CSF.native (job.output.CSF.dartel==1)    (job.output.CSF.dartel==2)   ];

  % never write class 4-6
  for i=4:6;
    tissue(i).warped = [0 0 0];
    tissue(i).native = [0 0 0];
  end

  job.bias     = [cat_get_defaults('output.bias.native')  cat_get_defaults('output.bias.warped') cat_get_defaults('output.bias.dartel')];
  job.label    = [job.output.label.native job.output.label.warped (job.output.label.dartel==1) (job.output.label.dartel==2)];
  job.jacobian = job.output.jacobian.warped;
  job.biasreg  = cat_get_defaults('opts.biasreg');
  job.biasfwhm = cat_get_defaults('opts.biasfwhm');
  job.channel  = struct('vols',{job.data});
  job.cat      = cat12;
  job.warps    = job.output.warps;
  job.tissue   = tissue;
  job.ignoreErrors = cat_get_defaults('extopts.ignoreErrors');
return;

%_______________________________________________________________________
function vout = run_job(job)
  vout   = vout_job(job);

  if ~isfield(job.cat,'fwhm'),    job.cat.fwhm    =  1; end

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
    if job.ignoreErrors
      if cat_io_matlabversion>20072 
        cat_run_newcatch(job,tpm,subj); 
      else
        % inactive because of unclear error messages
        %cat_run_oldcatch(job,tpm,subj);
        cat_run_job(job,tpm,subj);
      end
    else
      cat_run_job(job,tpm,subj);
    end
  end

  colormap(gray)
  
  if isfield(job,'nproc') && job.nproc>0 
    fprintf('\n%s',repmat('_',1,72));
    fprintf('\nCAT12 Segmentation job finished.\n');
  end
return
%_______________________________________________________________________

function vout = vout_job(job)

n     = numel(job.channel(1).vols);
parts = cell(n,4);

biascorr  = {};
wbiascorr = {};
label  = {};
wlabel = {};
rlabel = {};
alabel = {};

if cat_get_defaults('extopts.subfolders')
  mrifolder = 'mri';
else
  mrifolder = '';
end

for j=1:n,
    [parts{j,:}] = spm_fileparts(job.channel(1).vols{j});
end

if job.bias(1),
    biascorr = cell(n,1);
    for j=1:n
        biascorr{j} = fullfile(parts{j,1},mrifolder,['m',parts{j,2},'.nii']);
    end
end

if job.bias(2),
    wbiascorr = cell(n,1);
    for j=1:n
        wbiascorr{j} = fullfile(parts{j,1},mrifolder,['wm',parts{j,2},'.nii']);
    end
end

if job.label(1),
    label = cell(n,1);
    for j=1:n
        label{j} = fullfile(parts{j,1},mrifolder,['p0',parts{j,2},'.nii']);
    end
end

if job.label(2),
    wlabel = cell(n,1);
    for j=1:n
        wlabel{j} = fullfile(parts{j,1},mrifolder,['wp0',parts{j,2},'.nii']);
    end
end

if job.label(3),
    rlabel = cell(n,1);
    for j=1:n
        rlabel{j} = fullfile(parts{j,1},mrifolder,['rp0',parts{j,2},'.nii']);
    end
end

if job.label(4),
    alabel = cell(n,1);
    for j=1:n
        alabel{j} = fullfile(parts{j,1},mrifolder,['rp0',parts{j,2},'_affine.nii']);
    end
end

param = cell(n,1);
for j=1:n
    param{j} = fullfile(parts{j,1},['cat12_',parts{j,2},'.mat']);
end

tiss = struct('p',{},'rp',{},'rpa',{},'wp',{},'mwp',{},'m0wp',{});
for i=1:numel(job.tissue),
    if job.tissue(i).native(1),
        tiss(i).c = cell(n,1);
        for j=1:n
            tiss(i).c{j} = fullfile(parts{j,1},mrifolder,['p',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).native(2),
        tiss(i).rc = cell(n,1);
        for j=1:n
            tiss(i).rc{j} = fullfile(parts{j,1},mrifolder,['rp',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).native(3),
        tiss(i).rca = cell(n,1);
        for j=1:n
            tiss(i).rca{j} = fullfile(parts{j,1},mrifolder,['rp',num2str(i),parts{j,2},'_affine.nii']);
        end
    end
    if job.tissue(i).warped(1),
        tiss(i).wc = cell(n,1);
        for j=1:n
            tiss(i).wc{j} = fullfile(parts{j,1},mrifolder,['wp',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).warped(2),
        tiss(i).mwc = cell(n,1);
        for j=1:n
            tiss(i).mwc{j} = fullfile(parts{j,1},mrifolder,['mwp',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).warped(3),
        tiss(i).m0wc = cell(n,1);
        for j=1:n
            tiss(i).m0wc{j} = fullfile(parts{j,1},mrifolder,['m0wp',num2str(i),parts{j,2},'.nii']);
        end
    end
end

if job.cat.warps(1),
    fordef = cell(n,1);
    for j=1:n
        fordef{j} = fullfile(parts{j,1},mrifolder,['y_',parts{j,2},'.nii']);
    end
else
    fordef = {};
end

if job.cat.warps(2),
    invdef = cell(n,1);
    for j=1:n
        invdef{j} = fullfile(parts{j,1},mrifolder,['iy_',parts{j,2},'.nii']);
    end
else
    invdef = {};
end

if job.jacobian,
    jacobian = cell(n,1);
    for j=1:n
        jacobian{j} = '';
    end
else
    jacobian = {};
end

vout  = struct('tiss',tiss,'label',{label},'wlabel',{wlabel},'rlabel',{rlabel},'alabel',{alabel},...
               'biascorr',{biascorr},'wbiascorr',{wbiascorr},'param',{param},...
               'invdef',{invdef},'fordef',{fordef},'jacobian',{jacobian});
%_______________________________________________________________________

%=======================================================================

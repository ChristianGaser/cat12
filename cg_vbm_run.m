function varargout = cg_vbm_run(job,arg)
% Segment a bunch of images
% FORMAT cg_vbm_run(job)
% job.channel(n).vols{m}
% job.channel(n).biasreg
% job.channel(n).biasfwhm
% job.channel(n).write
% job.tissue(k).tpm
% job.tissue(k).ngaus
% job.tissue(k).native
% job.tissue(k).warped
% job.vbm.affreg
% job.vbm.reg
% job.vbm.samp
% job.vbm.write
% job.vbm.darteltpm
% job.vbm.print
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

% check whether estimation & write should be done
estwrite = isfield(job,'opts');

% set some defaults if segmentations are not estimated
if ~estwrite
    job.opts = struct('biasreg',0.001,'biasfwhm',60,'affreg','mni',...
                      'reg',[0 0.001 0.5 0.025 0.1],'samp',3,'ngaus',[3 3 2 3 4 2]);
end

channel = struct('vols',{job.data});

vbm = struct('species',   cg_vbm_get_defaults('extopts.species'), ... job.extopts.species,...
             'vbm12atlas',cg_vbm_get_defaults('extopts.vbm12atlas'), ... 
             'darteltpm', job.extopts.darteltpm{1}, ...
             'brainmask', cg_vbm_get_defaults('extopts.brainmask'), ...
             'affreg',    job.opts.affreg,...
             'samp',      cg_vbm_get_defaults('opts.samp'),...
             'write',     job.output.warps,...
             'sanlm',     cg_vbm_get_defaults('extopts.sanlm'),... % job.extopts.sanlm,...
             'print',     job.extopts.print,...
             'ngaus',     cg_vbm_get_defaults('opts.ngaus'),...
             'reg',       cg_vbm_get_defaults('opts.warpreg'),...
             'bb',        cg_vbm_get_defaults('extopts.bb'),...
             'vox',       cg_vbm_get_defaults('extopts.vox'));

if isfield(job.extopts,'restype')
  vbm.restype = char(fieldnames(job.extopts.restype));
  vbm.resval  = job.extopts.restype.(vbm.restype); 
else
  vbm.restype = cg_vbm_get_defaults('extopts.restype');
  vbm.resval  = cg_vbm_get_defaults('extopts.resval');
end
           
% set vbm.bb and vb.vox by Dartel template properties
Vd       = spm_vol([vbm.darteltpm ',1']);
[bb,vox] = spm_get_bbox(Vd, 'old');  
if vbm.bb(1)>vbm.bb(2), bbt=vbm.bb(1); vbm.bb(1)=vbm.bb(2); vbm.bb(2)=bbt; clear bbt; end
if bb(1)>bb(2), bbt=bb(1); bb(1)=bb(2); bb(2)=bbt; clear bbt; end
vbm.bb  = [ max(bb(1,1:3) , bb(1,1:3) ./ (isinf(bb(1,1:3)) | isnan(bb(1,1:3))))
            min(bb(2,1:3) , bb(2,1:3) ./ (isinf(bb(2,1:3)) | isnan(bb(2,1:3)))) ];
          
if isinf(vbm.vox) || isnan(vbm.vox)
  vbm.vox = abs(vox);
end



% prepare tissue priors and number of gaussians for all 6 classes
if estwrite
    [pth,nam,ext] = spm_fileparts(job.opts.tpm{1});
    clsn = numel(spm_vol(fullfile(pth,[nam ext]))); 
    tissue = struct();
    for i=1:clsn;
        tissue(i).ngaus = vbm.ngaus(i);
        tissue(i).tpm = [fullfile(pth,[nam ext]) ',' num2str(i)];
    end
end

% write tissue class 1-3              
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

job.bias     = [job.output.bias.native  job.output.bias.warped job.output.bias.affine];
job.label    = [job.output.label.native job.output.label.warped (job.output.label.dartel==1) (job.output.label.dartel==2)];
job.jacobian = job.output.jacobian.warped;

job.biasreg  = job.opts.biasreg;
job.biasfwhm = job.opts.biasfwhm;
job.channel  = channel;
job.vbm      = vbm;
job.warps    = job.output.warps;
job.tissue   = tissue;

if nargin == 1, arg = 'run'; end

switch lower(arg)
    case 'run'
       varargout{1} = run_job(job,estwrite);
    case 'check'
        varargout{1} = check_job(job);
    case 'vfiles'
        varargout{1} = vfiles_job(job);
    case 'vout'
        varargout{1} = vout_job(job);
    otherwise
        error('Unknown argument ("%s").', arg);
end

return
%_______________________________________________________________________

%_______________________________________________________________________
function vout = run_job(job,estwrite)
  vout   = vout_job(job);

  if ~isfield(job.vbm,'fwhm'),    job.vbm.fwhm    =  1; end

  % load tpm priors only for estimate and write
  if estwrite
      tpm = char(cat(1,job.tissue(:).tpm));
      tpm = spm_load_priors8(tpm);
  else
      tpm = '';
  end

  for subj=1:numel(job.channel(1).vols),
    % __________________________________________________________________
    % Separation for old and new try-catch blocks of matlab. The new
    % try-catch block has to be in a separate file to avoid an error.
    % Both functions finally call cg_vbm_run_job.
    % See also cg_vbm_run_newcatch and cg_vbm_run_newcatch.
    % __________________________________________________________________
    matlabversion = version; 
    points = strfind(matlabversion,'.');
    if str2double(matlabversion(1:points(1)-1))<=7 && ...
       str2double(matlabversion(points(1)+1:points(2)-1))<=5
      cg_vbm_run_oldcatch(job,estwrite,tpm,subj);
    else
      cg_vbm_run_newcatch(job,estwrite,tpm,subj);
    end
  end

  colormap(gray)

return
%_______________________________________________________________________


%_______________________________________________________________________
function msg = check_job(job)
msg = {};
if numel(job.channel) >1,
    k = numel(job.channel(1).vols);
    for i=2:numel(job.channel),
        if numel(job.channel(i).vols)~=k,
            msg = {['Incompatible number of images in channel ' num2str(i)]};
            break
        end
    end
elseif numel(job.channel)==0,
    msg = {'No data'};
end
return
%_______________________________________________________________________

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
%jacobian = {};

for j=1:n,
    [parts{j,:}] = spm_fileparts(job.channel(1).vols{j});
end

if job.bias(1),
    biascorr = cell(n,1);
    for j=1:n
        biascorr{j} = fullfile(parts{j,1},['m',parts{j,2},'.nii']);
    end
end

if job.bias(2),
    wbiascorr = cell(n,1);
    for j=1:n
        wbiascorr{j} = fullfile(parts{j,1},['wm',parts{j,2},'.nii']);
    end
end

if job.label(1),
    label = cell(n,1);
    for j=1:n
        label{j} = fullfile(parts{j,1},['p0',parts{j,2},'.nii']);
    end
end

if job.label(2),
    wlabel = cell(n,1);
    for j=1:n
        wlabel{j} = fullfile(parts{j,1},['wp0',parts{j,2},'.nii']);
    end
end

if job.label(3),
    rlabel = cell(n,1);
    for j=1:n
        rlabel{j} = fullfile(parts{j,1},['rp0',parts{j,2},'.nii']);
    end
end

if job.label(4),
    alabel = cell(n,1);
    for j=1:n
        alabel{j} = fullfile(parts{j,1},['rp0',parts{j,2},'_affine.nii']);
    end
end

param = cell(n,1);
for j=1:n
    param{j} = fullfile(parts{j,1},['vbm12_',parts{j,2},'.mat']);
end

tiss = struct('c',{},'rc',{},'rca',{},'wc',{},'mwc',{},'m0wc',{});
for i=1:numel(job.tissue),
    if job.tissue(i).native(1),
        tiss(i).c = cell(n,1);
        for j=1:n
            tiss(i).c{j} = fullfile(parts{j,1},['p',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).native(2),
        tiss(i).rc = cell(n,1);
        for j=1:n
            tiss(i).rc{j} = fullfile(parts{j,1},['rp',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).native(3),
        tiss(i).rca = cell(n,1);
        for j=1:n
            tiss(i).rca{j} = fullfile(parts{j,1},['rp',num2str(i),parts{j,2},'_affine.nii']);
        end
    end
    if job.tissue(i).warped(1),
        tiss(i).wc = cell(n,1);
        for j=1:n
            tiss(i).wc{j} = fullfile(parts{j,1},['wp',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).warped(2),
        tiss(i).mwc = cell(n,1);
        for j=1:n
            tiss(i).mwc{j} = fullfile(parts{j,1},['mwp',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).warped(3),
        tiss(i).m0wc = cell(n,1);
        for j=1:n
            tiss(i).m0wc{j} = fullfile(parts{j,1},['m0wp',num2str(i),parts{j,2},'.nii']);
        end
    end
end

if job.vbm.write(1),
    fordef = cell(n,1);
    for j=1:n
        fordef{j} = fullfile(parts{j,1},['y_',parts{j,2},'.nii']);
    end
else
    fordef = {};
end

if job.vbm.write(2),
    invdef = cell(n,1);
    for j=1:n
        invdef{j} = fullfile(parts{j,1},['iy_',parts{j,2},'.nii']);
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

%_______________________________________________________________________
function vf = vfiles_job(job)
vout = vout_job(job);
vf   = vout.param;
if ~isempty(vout.invdef),     vf = [vf vout.invdef]; end
if ~isempty(vout.fordef),     vf = [vf, vout.fordef]; end
if ~isempty(vout.jacobian),   vf = [vf, vout.jacobian]; end

if ~isempty(vout.biascorr),   vf = [vf, vout.biascorr]; end
if ~isempty(vout.wbiascorr),  vf = [vf, vout.wbiascorr]; end
if ~isempty(vout.label),      vf = [vf, vout.label]; end
if ~isempty(vout.wlabel),     vf = [vf, vout.wlabel]; end
if ~isempty(vout.rlabel),     vf = [vf, vout.rlabel]; end
if ~isempty(vout.alabel),     vf = [vf, vout.alabel]; end

for i=1:numel(vout.tiss)
    if ~isempty(vout.tiss(i).c),   vf = [vf vout.tiss(i).c];   end 
    if ~isempty(vout.tiss(i).rc),  vf = [vf vout.tiss(i).rc];  end 
    if ~isempty(vout.tiss(i).rca), vf = [vf vout.tiss(i).rca]; end
    if ~isempty(vout.tiss(i).wc),  vf = [vf vout.tiss(i).wc];  end
    if ~isempty(vout.tiss(i).mwc), vf = [vf vout.tiss(i).mwc]; end
    if ~isempty(vout.tiss(i).m0wc),vf = [vf vout.tiss(i).m0wc];end
end
vf = reshape(vf,numel(vf),1);
%_______________________________________________________________________

%=======================================================================

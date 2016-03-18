function cat_vol_sanlm(varargin)
% Spatial Adaptive Non Local Means (SANLM) Denoising Filter
%_______________________________________________________________________
% Filter a set of images and add the prefix 'sanlm_'.
% Missing input will call GUI or/and use defaults. 
%
% Input:
% job    - harvested job data structure (see matlabbatch help)
% 
% Output:
% out    - computation results, usually a struct variable.
%
% cat_vol_sanlm(job)
%   job.data   = set of images 
%   job.prefix = prefix for filtered images (default = 'sanlm_') 
%   job.rician = noise distribution
%
% Example:
%   cat_vol_sanlm(struct('data','','prefix','n','rician',0));
%_______________________________________________________________________
% Christian Gaser
% $Id$

if nargin == 0 
    job.data = cellstr(spm_select([1 Inf],'image','select images to filter'));
else
    job = varargin{1};
end

if ~isfield(job,'data') || isempty(job.data)
   job.data = cellstr(spm_select([1 Inf],'image','select images to filter'));
else
   job.data = cellstr(job.data);
end
if isempty(char(job.data)); return; end

def.verb   = 1;         % be verbose
def.prefix = 'sanlm_';  % prefix
def.NCstr  = inf;       % 0 - no denoising, eps - light denoising, 1 - maximum denoising, inf = auto; 
def.rician = 0;         % use inf for GUI
def.local  = 0;         % local weighing (only auto NCstr); 
job = cat_io_checkinopt(job,def);

job.NCstr = max(0,min(1,job.NCstr)) + isinf(job.NCstr)*job.NCstr;           % garanty values from 0 to 1 or inf
if isinf(job.rician), spm_input('Rician noise?',1,'yes|no',[1,0],2); end  % GUI

V = spm_vol(char(job.data));

spm_clf('Interactive'); 
spm_progress_bar('Init',numel(job.data),'SANLM-Filtering','Volumes Complete');
for i = 1:numel(job.data)
    [pth,nm,xt,vr] = spm_fileparts(deblank(V(i).fname));

    src = single(spm_read_vols(V(i)));
    % prevent NaN
    src(isnan(src)) = 0;

    if job.NCstr>0, srco = src + 0; end

    cat_sanlm(src,3,1,job.rician);

    % adaptive global denoising 
    if isinf(job.NCstr) || sign(job.NCstr)==-1;
      Yh     = src>mean(src(:)); % object
      Tth    = mean(src(Yh(:)));
      NCstr  = min(1,max(0, mean( abs(src(Yh(:)) - srco(Yh(:))) ./ Tth ) * 8 * min(2,max(0,abs(job.NCstr))) ));
      NCs    = abs(src - srco) / Tth * 8 * min(2,max(0,abs(job.NCstr))); spm_smooth(NCs,NCs,2);
      NCs    = max(0,min(1,NCs));
    else 
      NCstr  = job.NCstr;
    end
    
    % mix original and noise corrected image
    if job.local
      src = srco.*(1-NCs) + src.*NCs; 
    else
      src = srco*(1-NCstr) + src*NCstr; 
    end
    
    V(i).fname = fullfile(pth,[job.prefix nm '.nii' vr]);
    V(i).descrip = sprintf('%s SANLM filtered (NCstr=%d)',V(i).descrip,job.NCstr);
    
    if job.verb
      fprintf('NCstr = %0.2f: Output %s\n',NCstr,spm_file(V(i).fname,'link','spm_display(''%s'')'));
    end
    
    % use at least float precision
    if  V(i).dt(1)<16, V(i).dt(1) = 16; end 
    spm_write_vol(V(i), src);
    spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');


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

if ~isfield(job,'prefix')
    job.prefix = 'sanlm_';
end

if ~isfield(job,'rician') 
    job.rician = spm_input('Rician noise?',1,'yes|no',[1,0],2);
end

if ~isfield(job,'NCstr')
    job.NCstr = inf;
else
    job.NCstr = max(0,min(1,job.NCstr)) + isinf(job.NCstr);
end

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
    if isinf(job.NCstr);
      job.NCstr = min(1,max(0,mean(abs(src(:) - srco(:)))*10)); 
    end
    
    % mix original and noise corrected image and go back to original resolution
    src = srco*(1-job.NCstr) + src*job.NCstr; 
    
    V(i).fname = fullfile(pth,[job.prefix nm '.nii' vr]);
    V(i).descrip = sprintf('%s SANLM filtered (NCstr=%d)',V(i).descrip,job.NCstr);
    
    % use at least float precision
    if  V(i).dt(1)<16, V(i).dt(1) = 16; end 
    spm_write_vol(V(i), src);
    spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');


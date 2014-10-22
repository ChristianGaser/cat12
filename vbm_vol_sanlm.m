function vbm_vol_sanlm(varargin)
% Spatial Adaptive Non Local Means Denoising Filter
% Filter a set of images and add the prefix 'sanlm_'.
% 
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%
%_______________________________________________________________________
% Christian Gaser
% $Id$

if nargin == 0
    P = spm_select([1 Inf],'image','select images to filter');
    for i=1:size(P,1)
        job.data{i} = deblank(P(i,:));
    end
else
    job = varargin{1};
end

if ~isfield(job,'prefix')
    job.prefix = 'sanlm_';
end

if ~isfield(job,'rician') 
    job.rician = spm_input('Rician noise?',1,'yes|no',[1,0],2);
end

V = spm_vol(char(job.data));

spm_progress_bar('Init',numel(job.data),'SANLM-Filtering','Volumes Complete');
for i = 1:numel(job.data)
    [pth,nm,xt,vr] = spm_fileparts(deblank(V(i).fname));

    src = single(spm_read_vols(V(i)));
    % prevent NaN
    src(isnan(src)) = 0;
    try
        sanlmMex(src,3,1,job.rician);
    catch
        sanlmMex_noopenmp(src,3,1,job.rician);
    end

    V(i).fname = fullfile(pth,[job.prefix nm '.nii' vr]);
    V(i).descrip = sprintf('%s SANLM filtered',V(i).descrip);
    
    % use at least float precision
    if  V(i).dt(1)<16 V(i).dt(1) = 16; end 
    spm_write_vol(V(i), src);
    spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');


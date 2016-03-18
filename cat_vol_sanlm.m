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

def.verb    = 1;         % be verbose
def.prefix  = 'sanlm_';  % prefix
def.postfix = ''; 
def.NCstr   = inf;       % 0 - no denoising, eps - light denoising, 1 - maximum denoising, inf = auto; 
def.rician  = 0;         % use inf for GUI
def.local   = 1;         % local weighing (only auto NCstr); 
job = cat_io_checkinopt(job,def);

if job.NCstr==0, fprintf('HAHA. Nothing to do.\n'); end
if strcmp(job.postfix,'NCstr'), job.postfix = sprintf('_NCstr%0.2f',job.NCstr); end
job.NCstr = max(-2,min(1,job.NCstr)) + isinf(job.NCstr)*job.NCstr;           % garanty values from 0 to 1 or inf
if isinf(job.rician), spm_input('Rician noise?',1,'yes|no',[1,0],2); end  % GUI

V = spm_vol(char(job.data));

spm_clf('Interactive'); 
spm_progress_bar('Init',numel(job.data),'SANLM-Filtering','Volumes Complete');
for i = 1:numel(job.data)
    [pth,nm,xt,vr] = spm_fileparts(deblank(V(i).fname));

    src = single(spm_read_vols(V(i)));
    % prevent NaN
    src(isnan(src)) = 0;

    stime = clock; 
    if isinf(job.NCstr) || job.NCstr~=0, srco = src + 0; end
    
    % filter
    cat_sanlm(src,3,1,job.rician);

    % adaptive global denoising 
    %   vx_vol  = sqrt(sum(V(i).mat(1:3,1:3).^2)); 
    %   ds('d2','',vx_vol,srco/Tth,(srco.*(1-NCstr) + src.*NCstr) /Tth,src/Tth,(srco.*(1-NCs) + src.*NCs) /Tth,50);
    if isinf(job.NCstr) || sign(job.NCstr)==-1
      Yh     = src>mean(src(:)); % object
      Tth    = mean(src(Yh(:)));
      NCstr  = min(1,max(0, mean( abs(src(Yh(:)) - srco(Yh(:))) ./ Tth ) * 10 * min(1,max(0,abs(job.NCstr))) ));
      NC     = min(2,abs(src - srco) ./ max(eps,src) * 10 * min(2,max(0,abs(job.NCstr))));
      NCs    = NC+0; spm_smooth(NCs,NCs,2); NCs = NCs .* cat_stat_nanmean(NCs(Yh(:))) / cat_stat_nanmean(NC(Yh(:)));
      NCs    = max(0,min(1,NCs));
      
    else 
      NCstr  = job.NCstr;
    end
    
    %% mix original and noise corrected image
    if job.local
      src = srco.*(1-NCs) + src.*NCs; 
    else
      src = srco*(1-NCstr) + src*NCstr; 
    end
    
    V(i).fname = fullfile(pth,[job.prefix nm job.postfix '.nii' vr]);
    V(i).descrip = sprintf('%s SANLM filtered (NCstr=%-0.2f)',V(i).descrip,job.NCstr);
    
    if job.verb
      fprintf('NCstr = %-04.2f,%5.0fs: Output %s\n',sign(job.NCstr) * NCstr,etime(clock,stime),spm_file(V(i).fname,'link','spm_display(''%s'')'));
    end
    
    % use at least float precision
    if  V(i).dt(1)<16, V(i).dt(1) = 16; end 
    spm_write_vol(V(i), src);
    spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');


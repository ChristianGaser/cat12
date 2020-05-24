function varagout = cat_vol_createTPM(job)

% opt. 
%   class_number  = [6];      % 0 = auto; 
%   brain_classes = [1,2,3];  % add 7 if extra WMH class
%   brain_closing = 0;        % only if brain_classes without CSF 
%
% Case 1: Shooting template (simple)
%         Input:  Shooting templates
%                 [brain class ids to create brainmask] 
%         Output  TPM (6 classes), T1-simulate[, brainmask]
% 
% Case 2: Shooting template (full) 
%         Input: Shooting templates 
%                u-files 
%                   >> affine/rigid normalized classes
%                   >> mT1 files
%         Output: TPM, T1-simulate, (masked) T1-real, brainmask
%
% Case 3: normalized input >> create template ...
%         Input: wp1 files
%                wm  files
%         Output: TPM, T1-simulate, (masked) T1-real, brainmask

if ~exist('job','var'), job = struct(); end
  
def.template.files    = {}; % eg. Shooting with 2 to 6 clases
%def.template.ufiles   = {}; 
%def.template.mfiles   = {};
def.template.mixing   = [0.1 0.15 0.2 0.25 0.3]; % weighting for template number, e.g., Shooting with 5 templates [0.1 0.15 0.2 0.25 0.3] 

def.segments.wpfiles  = {}; % 1-2, 1-3, 1-5/6, 1-7
def.segments.wmfiles  = {};

def.opts.ssize        = [ 0   1   2   4   ];  % we use multipe smoothing levels to softly include also larger changes ...
def.opts.sweight      = [ 0.4 0.3 0.2 0.1 ];  % ... but with lower weighting 
def.opts.scsize       = [ 1 1 1 2 2 4 ];      % moreover we use a class specific filter size factor to suport smoother head classes


def.opts.ncls         = 6;            % public  - number of classes with background as last layer  
def.opts.verb         = 2;            % public  - be verbose (2=display result)
def.opts.name         = 'myTPM.nii';  % public  - template name
%def.opts.esthdcls     = 1;            % ?       - estimate head classes based on the T1 map if required and possible (if T1 is available)  
def.opts.bcls         = 1:3;          % public  - brain classes to define brainmask - there is maybe some subcortical classes
def.opts.bclosing     = 0.5;          % private - just remove wholes 
def.opts.bsmoohting   = 1;            % private - final smoothing of the brain mask in voxels 
def.opts.minprob      = 0.01;         % private - minimum prob in the whole image 
%def.opts.outdir       = ''; 

job = cat_io_checkinopt(job,def); 

if 1% job.opts.test
  %%
  job.template.files = {
    '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/chimpanzee_Template_0_GS.nii,1';
    '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/chimpanzee_Template_1_GS.nii,1';
    '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/chimpanzee_Template_2_GS.nii,1';
    '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/chimpanzee_Template_3_GS.nii,1';
    '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/chimpanzee_Template_4_GS.nii,1';
    };
  job.template.files = {
    '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/chimpanzee_Template_0_GS.nii,1';
    '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/chimpanzee_Template_1_GS.nii,1';
    '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/chimpanzee_Template_2_GS.nii,1';
    '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/chimpanzee_Template_3_GS.nii,1';
    '/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_animals/chimpanzee_Template_4_GS.nii,1';
    };
  
end

%% helping functions
cell2num = @(x) cell2mat( shiftdim(x(:), -ndims(x{1})) ) ;
clsnorm  = @(x)  shiftdim( num2cell( cell2num(x) ./ repmat( sum( cell2num(x) , ndims(x{1}) + 1) , ...
  [ones(1,ndims(x{1})),numel(x)])  , 1:ndims(x{1})) , -ndims(x{1})) ; 


%% load main data 
if isfield(job,'segments')
  % create average template
elseif  isfield(job,'template')
  %%
  Vtemp0 = spm_vol(job.template.files{1}); 
  vx_vol = [1 1 1]; 
  
  if job.opts.ncls>0
    ncls = job.opts.ncls; 
  else
    ncls = Vtemp0.mat;  
  end
  
  if numel( job.template.mixing ) == 1
    job.template.mixing = repmat( job.template.mixing , 1 , ncls );
  end
  
  % create empty TPM
  Ytpm = cell(1,ncls);
  for ci = 1:ncls
    Ytpm{ci} = zeros(Vtemp0.dim,'single'); 
  end
  
  % load template
  tempdims = 1;
  for ti = 1:numel(job.template.files)
    [pp,ff,ee] = spm_fileparts(job.template.files{ti});
    for ci = 1:ncls
      Vtemp    = spm_vol( fullfile(pp,sprintf('%s%s,%d',ff,ee,ci)) ); 
      if ~isempty(Vtemp)
        Ytemp    = single( spm_read_vols( Vtemp ) ); 
        spm_smooth(Ytemp,Ytemp,job.template.mixing(ci) ./ vx_vol);  
        Ytpm{ci} = Ytpm{ci} + Ytemp .* job.template.mixing(ti);  
        tempdims = ci;
      end
    end
  end
end


%% find background
%  We cannot be sure if we got a background layer and which one it is. 
%  Hence, we use a box at the image boundary to count the asigned voxel to 
%  detect if there is a background in any given layer. If not we have to 
%  create one by the missed voxels, else we have to guaranty that the last 
%  class is the background. 
Ybgb = true(size(Ytpm{1})); Ybgb(3:end-2,3:end-2,3:end-2) = false;
tpmbg = zeros(1,ncls); for ci = 1:ncls, tpmbg(ci) = sum(Ytpm{ci}(Ybgb)) / sum(Ybgb(:)); end; clear Ybgb; 
[t,bg] = max(tpmbg .* (tpmbg>0.5)); clear t;  %#ok<ASGLU> % bg is the old background class, we will need it


%  create brainmask
%  Our background maybe also include CSF or other low intensity boxels. 
%  So we have to create a brain mask and correct it. 
Yb = sum( cell2num( Ytpm( unique( min( setdiff( job.opts.bcls , bg) , numel(Ytpm) )) ) ) , 4); 
Yb = max(Yb, cat_vol_smooth3X(single(cat_vol_morph(Yb>0.1,'lc',job.opts.bclosing)),job.opts.bsmoothing)); 


%  create head classes



%% If all classes are empty then we can create a new background. 
if isempty(bg)
  Ybg = sum( cell2num( Ytpm ) , 4);
  Ybg = max(Ybg(:)) - Ybg;
  Ybg = Ybg .* (1-Yb); % without brain voxel 
  Ytpm( job.opts.ncls ) = Ybg; 
else
  if bg < numel(Ytpm)
    Ybg       = Ytpm{ bg };
    Ybg       = Ybg .* (1-Yb); % without brain voxel 
    Ytpm{end} = Ybg;
    Ytpm{bg}  = Ytpm{bg} .* Yb; 
  end
end
clear tpmbg 



%% normalize probabilites
Ytpm = clsnorm(Ytpm);
 

% - create cat atlas

% - normalize all classes

% - write output 

% - display output



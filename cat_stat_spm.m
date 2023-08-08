function varargout = cat_stat_spm(SPM)
% Workaround to use fsaverage surface as SurfaceID (for displaying results)
% spm_spm is used to estimate the model and the mesh of the 1st file in the model 
% is replaced by the fsaverage brain because this mesh is used for overlaying
% results.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

if nargin == 0
  P = spm_select([1 Inf],'^SPM\.mat$','Select SPM.mat file(s)');
elseif exist('SPM','var') && isfield(SPM,'spmmat')
  P = char(SPM.spmmat); spmmat = SPM.spmmat;
elseif ischar(SPM)
  P = SPM;
end

if exist('P','var')
  for i=1:size(P,1)
    swd = spm_file(P(i,:),'fpath');
    load(fullfile(swd,'SPM.mat'));
    SPM.swd  = swd; 
    cat_stat_spm(SPM);
  end
  if nargout && exist('spmmat','var') 
    varargout{1}.spmmat = spmmat;
  end
  return
end

if ~isfield(SPM,'xY')
  error(sprintf('SPM.mat was not correctly saved. Please check that you have set the following flag in spm_defaults if your files are > 2GB:\ndefaults.mat.format = ''-v7.3'''));
end

fmt = spm_get_defaults('mat.format');
s = whos('SPM');
if s.bytes > 2147483647, fmt = '-v7.3'; end

% older formats don't support large files
spm_get_defaults('mat.format',fmt);

% check for 32k meshes
if SPM.xY.VY(1).dim(1) == 32492 || SPM.xY.VY(1).dim(1) == 64984
  fsavgDir = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces_32k');
else
  fsavgDir = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces');
end

% select underlying surface and prefer shooting template
job.surftype = 1; 
surftype = {'freesurfer',cat_get_defaults('extopts.shootingsurf')};
if ~exist(fullfile(fsavgDir, ['lh.central.' surftype{job.surftype} '.gii']))
  job.surftype = 2; 
end

% check that folder exist and number of vertices fits
if exist(fsavgDir,'dir') == 7 && (SPM.xY.VY(1).dim(1) == 163842 || SPM.xY.VY(1).dim(1) == 327684 || ...
    SPM.xY.VY(1).dim(1) == 655368) || SPM.xY.VY(1).dim(1) == 32492 || SPM.xY.VY(1).dim(1) == 64984
  
  [pp,ff]   = spm_fileparts(SPM.xY.VY(1).fname);

  % find mesh string      
  hemi_ind = [];
  hemi_ind = [hemi_ind strfind(ff,'mesh.')];
  if ~isempty(hemi_ind)
    
    SPM.xY.VY(1).private.private.metadata = struct('name','SurfaceID','value',fullfile(fsavgDir, ['mesh.central.' surftype{job.surftype} '.gii']));
    M0 = gifti({fullfile(fsavgDir, ['lh.central.' surftype{job.surftype} '.gii']), fullfile(fsavgDir, ['rh.central.' surftype{job.surftype} '.gii'])});
    G.faces = [M0(1).faces; M0(2).faces+size(M0(1).vertices,1)];
    G.vertices = [M0(1).vertices; M0(2).vertices];

    % cerebellar lobes?
    if SPM.xY.VY(1).dim(1) == 655368
      M0 = gifti({fullfile(fsavgDir, 'cb.central.freesurfer.gii')});  %, fullfile(fsavgDir, 'rc.central.freesurfer.gii')});
      G.faces = [G.faces; M0(1).faces+2*size(M0(1).vertices,1)];      % ; M0(2).faces+3*size(M0(1).vertices,1)];
      G.vertices = [G.vertices; M0(1).vertices];                      % ; M0(2).vertices];
    end
    clear M0;
    
    SPM.xVol.G = G;
    
    % remove memory demanding faces and vertices which are not necessary
    for i=1:length(SPM.xY.VY)
      SPM.xY.VY(i).private.faces = [];
      SPM.xY.VY(i).private.vertices = [];
    end
    
  else

    % find lh|rh string
    hemi_ind = [];
    hemi_ind = [hemi_ind strfind(ff,'lh.')];
    hemi_ind = [hemi_ind strfind(ff,'rh.')];
    hemi = ff(hemi_ind:hemi_ind+1);
    if ~isempty(hemi)
      SPM.xY.VY(1).private.private.metadata = struct('name','SurfaceID','value',fullfile(fsavgDir,[hemi '.central.' surftype{job.surftype} '.gii']));
      G = fullfile(fsavgDir,[hemi '.central.' surftype{job.surftype} '.gii']);
      SPM.xVol.G = gifti(G);
      
      % remove memory demanding faces and vertices which are not necessary
      for i=1:length(SPM.xY.VY)
        SPM.xY.VY(i).private.faces = [];
        SPM.xY.VY(i).private.vertices = [];
      end
      
    end
  end 
end

if nargout>0
  varargout{1} = spm_spm(SPM);
else
  spm_spm(SPM);
end  
end
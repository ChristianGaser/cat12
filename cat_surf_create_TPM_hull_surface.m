function Phull = cat_surf_create_TPM_hull_surface(tpm,human,skull)
% _________________________________________________________________________
% Creates a surface of the brain and headmask and save the data in one file
% named as "bh.headbrain.$TPM-filename$.gii" in the cat surface directory.
% If the file allready exist only the filename is returned otherwise a new 
% brain/head surface file is created. 
%
%   Phull =  cat_surf_create_TPM_hull_surface(tpm[,human,skull])
%   Shull =  cat_surf_create_TPM_hull_surface(tpm[,human,skull])
% 
%   tpm   .. TPM-filename, TPM-header, SPM-TPM-structure 
%   human .. flag to write animal templates into another directory
%   skull .. create outline without skull
%   Phull .. filename of the surface file 
%   Shull .. export surface in case of reading/writing errors
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$ 


  % check input
  if ~exist('tpm','var')
    Ptpm = fullfile(spm('dir'),'tpm','TPM.nii');
  else
    if ischar(tpm)
      Ptpm = tpm;
      clear tpm; 
    elseif iscell(tpm)
      Ptpm = char(tpm);
      clear tpm;
    else
      if numel(tpm)==6
        Ptpm = tpm(1).fname; 
        clear tpm; 
      else
        % SPM-TPM structure
        try
          Ptpm = tpm.V(1).fname; 
        catch
          Ptpm = tpm(1).fname; 
        end
      end
    end
  end
  if ~exist('human','var')
    species = '';
    human   = 1;
  else
    if isnumeric(human) || islogical(human)
      species = '';
      human   = human > 0;
    else % it is a string
      species = ['.' human];
      human   = strcmp(human,'human');
    end
  end
  if ~exist('skull','var')
    skull = 1;
  end
  
  % define filename
  [pp,Pname,ee] = spm_fileparts(Ptpm); Pname = strrep([Pname ee],'.nii',''); 
  if human
    if skull > 0
      Phull = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces',sprintf('bh.headbrain%s.%s.gii',species,Pname));
    else
      Phull = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces',sprintf('bh.brain%s.%s.gii',species,Pname));
    end    
  else
    if skull > 0
      Phull = fullfile(spm('dir'),'toolbox','cat12','templates_animals_surfaces',sprintf('bh.headbrain%s.%s.gii',species,Pname));
    else
      Phull = fullfile(spm('dir'),'toolbox','cat12','templates_animals_surfaces',sprintf('bh.brain%s.%s.gii',species,Pname));
    end
  end
  
  % nothing to do - just return filename 
  if ~cat_io_rerun(Ptpm,Phull,0,1), return; end
    
  
  % load SPM-TPM-structure
  if ~exist('tpm','var')
    % here we load the full TPM structure 
    istpm = 1; 
    tpm = spm_load_priors8(Ptpm);
  
    % create brainmask surface
    if skull==-2
      Yb = exp(tpm.dat{1}) + exp(tpm.dat{2});
    else
      Yb = exp(tpm.dat{1}) + exp(tpm.dat{2}) + exp(tpm.dat{3});
      % remove SPM CSF eye
      Yb = Yb .* smooth3( cat_vol_morph( cat_vol_morph( Yb , 'lo' , 2), 'd') ); 
    end
  elseif isstruct(tpm)
    % if we got the TPM as file header structure
    istpm = 2; 
    try 
      Yb = spm_read_vol( V(1) ) + spm_read_vol( V(2) ); 
    catch 
      % create brainmask surface
      if skull==-2
        Yb = exp(tpm.dat{1}) + exp(tpm.dat{2});
      else
        Yb = exp(tpm.dat{1}) + exp(tpm.dat{2}) + exp(tpm.dat{3});
        % remove SPM CSF eye
        Yb = Yb .* smooth3( cat_vol_morph( cat_vol_morph( Yb , 'lo' , 2), 'd') ); 
      end
    end
  else
    % otherwise avoid errors
    Yb = zeros(5,5,5);  
  end
  Sh   = isosurface(Yb,0.5);
  
  
  % create head surface based on 5 classes with average probability threshold
  if skull && istpm
    % use all classes but nat the background 
    if istpm == 1
      Yhd = exp(tpm.dat{1});
      for i = 2:numel(tpm.dat)-1
         Yhd = Yhd + exp(tpm.dat{i});
      end
    elseif istpm == 2
      Yhd = spm_read_vol( tpm(1) );
      for i = 2:numel(tpm)-1
         Yhd = Yhd + spm_read_vol( tpm(i) );
      end
    end
    Yhd = cat_vol_morph(Yhd>0.5,'l',[1 0.5]);
    Yhd = ~cat_vol_morph(~Yhd,'l',[1 0.5]);
    if strcmpi(spm_check_version,'octave')
      Shd = spm_mesh_isosurface(smooth3(Yhd),0.5 - ( 0.4 * (skull<0)),0.5);
    else
      Shd = isosurface(Yhd,0.5 - ( 0.4 * (skull<0)) );
    end
    Sh.faces    = [Sh.faces; Shd.faces + size(Sh.vertices,1)];
    Sh.vertices = [Sh.vertices; Shd.vertices];
  end
  if strcmpi(spm_check_version,'octave')
    Sh  = spm_mesh_reduce(Sh,0.2);
  else
    Sh  = reducepatch(Sh,0.2);
  end    
    
  % save surface
  if istpm == 1
    vmat  = tpm.V(1).mat(1:3,:)*[0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];
  elseif istpm == 2
    vmat  = tpm(1).mat(1:3,:)*[0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];
  end  
  if istpm
    Sh.vertices = (vmat*[Sh.vertices' ; ones(1,size(Sh.vertices,1))])'; 
    mati = spm_imatrix(tpm.V(1).mat); 
    if mati(7)<0, Sh.faces = [Sh.faces(:,1) Sh.faces(:,3) Sh.faces(:,2)]; end
    try
      save(gifti(struct('faces',Sh.faces,'vertices',Sh.vertices)),Phull);   
    catch
      if exist(Phull,'file')
        fprintf('Warning: Could not update %s with newer version. Please change write permissions.\n',Phull);
      else
        error('Could not save %s. Please change write permissions.\n',Phull);
      end
    end
  end
end
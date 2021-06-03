function Phull = cat_surf_create_TPM_hull_surface(tpm,human)
% _________________________________________________________________________
% Creates a surface of the brain and headmask and save the data in one file
% named as "bh.headbrain.$TPM-filename$.gii" in the cat surface directory.
% If the file allready exist only the filename is returned otherwise a new 
% brain/head surface file is created. 
%
%   Phull =  cat_surf_create_TPM_hull_surface(tpm)
% 
%   tpm   .. TPM-filename, TPM-header, SPM-TPM-structure 
%   Phull .. filename of the surface file 
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
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
    human = 1;
  end

  
  % define filename
  [pp,Pname,ee] = spm_fileparts(Ptpm); Pname = strrep([Pname ee],'.nii',''); 
  if human
    Phull = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces',sprintf('bh.headbrain.%s.gii',Pname));
  else
    Phull = fullfile(spm('dir'),'toolbox','cat12','templates_animals_surfaces',sprintf('bh.headbrain.%s.gii',Pname));
  end
  
  % nothing to do - just return filename 
  if ~cat_io_rerun(Phull,Ptpm), return; end
    
  
  % load SPM-TPM-structure
  if ~exist('tpm','var')
    tpm  = spm_load_priors8(Ptpm);
  end
  
  
  % create brainmask surface
  Yb   = exp(tpm.dat{1}) + exp(tpm.dat{2}) + exp(tpm.dat{3});
  Sh   = isosurface(Yb,0.5);
  % create head surface based on 5 classes with average probability threshold
  Yhd = exp(tpm.dat{1}) + exp(tpm.dat{2}) + exp(tpm.dat{3}) + exp(tpm.dat{4}) + exp(tpm.dat{5});
  Yhd = cat_vol_morph(Yhd>0.5,'l',[1 0.5]);
  Yhd = ~cat_vol_morph(~Yhd,'l',[1 0.5]);
  Shd = isosurface(Yhd,0.5);
  Sh.faces    = [Sh.faces; Shd.faces + size(Sh.vertices,1)];
  Sh.vertices = [Sh.vertices; Shd.vertices];
  Sh  = reducepatch(Sh,0.2);
      
    
  % save surface
  vmat  = tpm.V(1).mat(1:3,:)*[0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];
  Sh.vertices = (vmat*[Sh.vertices' ; ones(1,size(Sh.vertices,1))])'; 
  mati = spm_imatrix(tpm.V(1).mat); 
  if mati(7)<0, Sh.faces = [Sh.faces(:,1) Sh.faces(:,3) Sh.faces(:,2)]; end
  try
    save(gifti(struct('faces',Sh.faces,'vertices',Sh.vertices)),Phull);   
  catch
    if exist(Phull)
      fprintf('Warning: Could not update %s with newer version. Please change write permissions.\n',Phull);
    else
      error('Could not save %s. Please change write permissions.\n',Phull);
    end
  end
end
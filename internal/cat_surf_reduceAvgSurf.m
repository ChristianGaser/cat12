function cat_surf_reduceAvgSurf(P)
% cat_surf_reduceAvgSurf(P)
% ______________________________________________________________________
% Function to reduce a set of average surface to allow faster
% processing. Although this function worked well, surface processing 
% speed was not improved :/
%
% WARNING: 
%   The reduction based on the geometry of the average is will trend to 
%   have large face. 
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
% ______________________________________________________________________
% $Id$ 

  
  if ~isfield('job','P') || isempty(P)
    job.P = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces','lh.central.freesurfer.gii');
  end
  
  def.reduce   = 1;       % 1 = Matlab; 2=Caret
  def.debug    = 1;       % display debug information
  def.reduceCS = 20000;   % number of faces of the reduce surfaces
  def.CATDir   = fullfile(spm('dir'),'toolbox','cat12','CAT');   
  
  job = cat_io_checkinopt(job,def);
  
  % add system dependent extension to CAT folder
  if ispc
    job.CATDir = [job.CATDir '.w32'];
  elseif ismac
    job.CATDir = [job.CATDir '.maci64'];
  elseif isunix
    job.CATDir = [job.CATDir '.glnx86'];
  end  
  

  
  side = {'lh.','rh.'};
  for si = 1:numel(side)
    %%
    fprintf('Reduce Average Surface %s\n',side{si})
  
    % surface names
    [pp,ff,ee]   = spm_fileparts(job.P); 
    
    ff = strrep(ff,'lh.',side{si});
    job.Po       = fullfile(pp,sprintf('%s%s',ff,ee));
    job.Pr       = fullfile(pp,sprintf('%s.%dk%s',ff,job.reduceCS/1000,ee));

    
    % reduce the original surface
    if job.reduce == 1 % matlab
      CS  = gifti(job.Po); clear CSr;
      CSr.vertices = double(CS.vertices); CSr.faces = double(CS.faces); 
      CSr = reducepatch(CSr,job.reduceCS); 
      save(gifti(struct('faces',int32(CSr.faces),'vertices',single(CSr.vertices))),job.Pr);

      if 0
        % refine to avoid to large faces ... 
        % the problem is that this will introduce further vertices that do
        % not exist in the original surface 
        % the alternative will be to use caret surface reduction that do
        % this more equaly 
        meshres = 2 * 100000/job.reduceCS; 
        cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f',job.Pr,job.Pr,meshres); 
        [ST, RS] = system(fullfile(job.CATDir,cmd)); cat_check_system_output(ST,RS,job.debug);
      end
    elseif job.reduce == 2 % Caret
      % Caret Dir
      % hmmm ... no caret_command :/

      
    else
      error('job.reduce has to be 1 (Matlab) or 2 (Caret). No reduction make no sense.\n');
    end
    
    fprintf('  Display reduced average %s\n',spm_file(job.Pr,'link','cat_surf_display(''%s'')'));

    
    % create mapping between old and new surace
    CSr = gifti(job.Pr);
    vusedo = true(size(CS.vertices,1),1); 
    for vi=1:size(CSr.vertices,1)
      vusedo(find(all([CS.vertices(:,1)==CSr.vertices(vi,1), ...
                       CS.vertices(:,2)==CSr.vertices(vi,2), ...
                       CS.vertices(:,3)==CSr.vertices(vi,3)],2),1,'first'))=0;
    end
        
    % use mapping for other surface
    other = {'sphere','inflated'};
    for oi = 1:numel(other)
      job.Pox{oi} = strrep(job.Po,'central',other{oi}); 
      job.Prx{oi} = strrep(job.Pr,'central',other{oi}); 
      
      CSo = gifti(job.Pox{oi});
      
      CSon.vertices = CSo.vertices(vusedo,:); 
      CSon.faces    = CSr.faces;

      save(gifti(struct('faces',int32(CSon.faces),'vertices',single(CSon.vertices))),job.Prx{oi});

      fprintf('  Display reduced %s %s\n',other{oi},spm_file(job.Pr,'link','cat_surf_display(''%s'')'));
    end
  end

end

function vbm_surf_display(varargin)
% ______________________________________________________________________
% Function to display surfaces. Based on spm_mesh_render.
%
% [Psdata] = vbm_surf_display(job)
% 
% job.data_resample
% job.fwhm
% ______________________________________________________________________
% Robert Dahnke
% $Id$

  if ~exist('varargin','var') 
    if isstruct(varargin)
      P = varargin{1}.data;
    else
      P = varargin{1};
    end
  else
    P = spm_select([1 24],'any','Select surface','','','[lr]h.*');
  end
  P = cellstr(P);
  
  sinfo = vbm_surf_info(P); 
  for i=1:numel(P)
    %% 
    
    % old ... % h = spm_mesh_render(sinfo(i).Pmesh);
    
    if ~strcmp(sinfo(i).Pmesh,sinfo(i).Pdata) && ~isempty(sinfo(i).Pdata)
      % only gifti surface without texture
      h = vbm_surf_render(sinfo(i).Pmesh,'Pcdata',sinfo(i).Pdata);
    else
      % only gifti surface without texture
      h = vbm_surf_render(sinfo(i).Pmesh);
      set(h.patch,'AmbientStrength',0.2,'DiffuseStrength',0.8,'SpecularStrength',0.1)
    end
    
    % textur handling
    set(h.figure,'MenuBar','none','Toolbar','none','Name',spm_file(P{i},'short60'),'NumberTitle','off');
    vbm_surf_render('ColourBar',h.axis,'on');
    switch sinfo(i).texture
      case 'thickness'
        vbm_surf_render('ColourMap',h.axis,jet); 
        vbm_surf_render('clim',h.axis,[0 5]);
      case 'gyrification'
        vbm_surf_render('ColourMap',h.axis,vbm_io_colormaps('curvature',128)); 
        %vbm_surf_render('clim',h.axis,[0 90]);
      case 'logsulc'
        vbm_surf_render('ColourMap',h.axis,vbm_io_colormaps('hotinv',128));
        vbm_surf_render('clim',h.axis,[0 1.5]);
    end

  end
end
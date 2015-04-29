function vbm_surf_display(varargin)
% ______________________________________________________________________
% Function to display surfaces. Warapper  to vbm_surf_render.
%
% [Psdata] = vbm_surf_display(job)
% 
% job.data_resample
% job.fwhm
% ______________________________________________________________________
% Robert Dahnke
% $Id$

  if nargin>0
    if isstruct(varargin)
      P = varargin{1}.data;
    else
      P = varargin{1};
    end
  else
    P = spm_select([1 24],'any','Select surface','','','[lr]h.*');
  end
  if isempty(P), return; end
  P = cellstr(P);
  
  
  %%
  sinfo = vbm_surf_info(P); 
  for i=1:numel(P)
    try

      % old ... % h = spm_mesh_render(sinfo(i).Pmesh);

      if ~strcmp(sinfo(i).Pmesh,sinfo(i).Pdata) && ~isempty(sinfo(i).Pdata)
        % only gifti surface without texture
        h = vbm_surf_render(sinfo(i).Pmesh,'Pcdata',sinfo(i).Pdata);
      else
        % only gifti surface without texture
        h = vbm_surf_render(sinfo(i).Pmesh);
      end

      % textur handling
      set(h.figure,'MenuBar','none','Toolbar','none','Name',spm_file(P{i},'short60'),'NumberTitle','off');
      vbm_surf_render('ColourBar',h.axis,'on');
      if strcmp(sinfo(i).side,'rh'), view(h.axis,[90 0]); end
      switch sinfo(i).texture
        case 'thickness'
          vbm_surf_render('ColourMap',h.axis,jet); 
          clim = iscaling(h.cdata);
          if clim(1)>3 || clim(1)<0  || clim(2)<2 || clim(2)>8  
            vbm_surf_render('clim',h.axis,clim);
          else % default range
            vbm_surf_render('clim',h.axis,[0 5]);
          end
        case {'curvature','gyrification'}
          vbm_surf_render('ColourMap',h.axis,vbm_io_colormaps('curvature',128)); 
          vbm_surf_render('clim',h.axis,iscaling(h.cdata));
        case 'logsulc'
          vbm_surf_render('ColourMap',h.axis,vbm_io_colormaps('hotinv',128));
          vbm_surf_render('clim',h.axis,iscaling(h.cdata));
          %vbm_surf_render('clim',h.axis,[0 1.5]);
        case 'defects'
          %set(h.patch,'DiffuseStrength',0,
        case 'sphere'
          % exist curvature or depth???
        case 'central'
          set(h.patch,'AmbientStrength',0.2,'DiffuseStrength',0.8,'SpecularStrength',0.1)
        otherwise
          if ~isempty(h.cdata)
            clim = iscaling(h.cdata);
            if clim(1)<0
              clim = [-max(abs(clim)) max(abs(clim))];
              vbm_surf_render('ColourMap',h.axis,vbm_io_colormaps('BWR',128)); 
            else
              vbm_surf_render('ColourMap',h.axis,vbm_io_colormaps('hotinv',128)); 
            end
            vbm_surf_render('clim',h.axis,clim);
          end
      end
    catch %#ok<CTCH>
      try
        h = vbm_surf_render(P{i});
      catch
        vbm_io_cprintf('err',sprintf('ERROR: Can''t display surface %s.\n',P{i})); 
      end
    end
  end
end
function clim = iscaling(cdata,plim)
%%
  ASD = min(0.02,max(eps,0.05*std(cdata))/max(abs(cdata))); 
  if ~exist('plim','var'), plim = [ASD 1-ASD]; end 

  bcdata  = [min(cdata) max(cdata)]; 
  range   = bcdata(1):diff(bcdata)/1000:bcdata(2);
  hst     = hist(cdata,range);
  clim(1) = range(max(1,find(cumsum(hst)/sum(hst)>plim(1),1,'first')));
  clim(2) = range(min(numel(range),find(cumsum(hst)/sum(hst)>plim(2),1,'first')));
end




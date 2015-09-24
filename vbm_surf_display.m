function varargout = vbm_surf_display(varargin)
% ______________________________________________________________________
% Function to display surfaces. Wrapper to vbm_surf_render.
%
% [Psdata] = vbm_surf_display(job)
% 
% job.data      .. [rl]h.* surfaces 
% job.colormap  .. colormap
% job.caxis     .. range of the colormap
% job.multisurf .. load both sides, if possible  
% job.view      .. view 
%                   l=left, r=right
%                   a=anterior, p=posterior
%                   s=superior, i=inferior
%
% job.imgprint.do   .. print image (default = 0)
% job.imgprint.type .. render image type (default = png)
% job.dpi           .. print resolution of the image (default = 600 dpi)
%
% ______________________________________________________________________
% Robert Dahnke
% $Id$

  SVNid = '$Rev$';

  if nargin>0
    if isstruct(varargin{1})
      job = varargin{1};
      if ~isfield(job,'data')
        job.data = spm_select([1 24],'any','Select surface','','','[lr]h.*');
        job.imgprint.do    = 0;
        job.imgprint.close = 0;  
      end
    else
      job.data = varargin{1};
    end
  else
    job.data = spm_select([1 24],'any','Select surface','','','[lr]h.*');
    job.imgprint.do    = 0;
    job.imgprint.close = 0;  
  end
  if isempty(job.data), return; end
  job.data = cellstr(job.data);
  
  % scaling options for textures
  def.colormap = '';
  def.caxis    = []; % default/auto, range
  
  % print options ... just a quick output > vbm_surf_print as final function 
  def.imgprint.type  = '-dpng';
  def.imgprint.dpi   = 600;
  def.imgprint.fdpi  = @(x) ['-r' num2str(x)];
  def.imgprint.do    = 1;
  def.imgprint.close = 0;

  % multi-surface output for one subject 
  def.multisurf = 0; % 0 - no; 1 - both hemispheres;
  
  job = checkinopt(job,def);
  
  %%
  sinfo = vbm_surf_info(job.data); 
  spm('FnBanner',mfilename,SVNid); 
  for i=1:numel(job.data)
    
    % load multiple surfaces
    if job.multisurf
      if strcmp(sinfo(i).side,'rh'), oside = 'lh'; else oside = 'rh'; end
      Pmesh = [sinfo(i).Pmesh vbm_surf_rename(sinfo(i).Pmesh,'side',oside)];
      Pdata = [sinfo(i).Pdata vbm_surf_rename(sinfo(i).Pdata,'side',oside)]; 
    else
      Pmesh = sinfo(i).Pmesh;
      Pdata = sinfo(i).Pdata; 
    end
    
    
    
    try
      fprintf('  %s\n',job.data{i});

      if ~all(strcmp(Pmesh,Pdata)) && ~isempty(Pdata) && (~job.multisurf || ~all(cellfun('isempty',Pdata)))
        % only gifti surface without texture
        if isfield(job,'parent')
          h = vbm_surf_render('disp',Pmesh,'Pcdata',Pdata,'parent',job.parent);
        else
          h = vbm_surf_render('disp',Pmesh,'Pcdata',Pdata);
        end  
      else
        % only gifti surface without texture
        if isfield(job,'ah')
          h = vbm_surf_render(Pmesh,'parent',job.parent);
        else
          h = vbm_surf_render(Pmesh);
        end
      end
      
      
      %% textur handling
      set(h.figure,'MenuBar','none','Toolbar','none','Name',spm_file(job.data{i},'short60'),'NumberTitle','off');
      vbm_surf_render('ColourBar',h.axis,'on');
      if ~job.multisurf && strcmp(sinfo(i).side,'rh'), view(h.axis,[90 0]); end
      
      
      % colormap
      if isempty(job.colormap)
        vbm_surf_render('ColourMap',h.axis,jet(256)); 
      else
        vbm_surf_render('ColourMap',h.axis,eval(job.colormap));
      end
      
      % scaling
      if isempty(job.caxis)
        switch sinfo(i).texture
          case {'defects','sphere'}
            % no texture
          case {'central'}
            % default curvature
            set(h.patch,'AmbientStrength',0.2,'DiffuseStrength',0.8,'SpecularStrength',0.1)
          case ''
            % no texture name
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
          otherwise
            %%
            ranges = {
              ... name single group
              'thickness'     [0.5  4.0]  [0.5  4.0]
              'gyruswidthWM'  [0.5  8.0]  [1.0  7.0]
              'gyruswidth'    [1.0 12.0]  [1.5 11.0]
              'sulcuswidth'   [0.0  3.0]  [0.0  3.0]
              'gyrification'  [0.0  1.0]  [0.0  0.5]
              'logsulc'       [0.0  1.5]  [0.0  1.5]
            };

            texturei = find(cellfun('isempty',strfind(ranges(:,1),sinfo(i).texture))==0,1,'first');

            if ~isempty(texturei)
              vbm_surf_render('clim',h.axis,ranges{texturei,3});
            else
              clim = iscaling(h.cdata);  
              vbm_surf_render('clim',h.axis,round(clim));
            end
        end     
      else
        vbm_surf_render('clim',h.axis,job.caxis);
      end
    catch %#ok<CTCH>
      try
        h = vbm_surf_render(job.data{i});
      catch %#ok<CTCH>
        vbm_io_cprintf('err',sprintf('ERROR: Can''t display surface %s.\n',job.data{i})); 
      end
    end
    
    
    
    
    %% view
    viewname = '';
    if isfield(job,'view')
      switch lower(job.view)
        case {'r','right'},                 view([  90   0]); viewname = '.r';
        case {'l','left'},                  view([ -90   0]); viewname = '.l';
        case {'t','s','top','superior'},    view([   0  90]); viewname = '.s';
        case {'b','i','bottom','inferior'}, view([-180 -90]); viewname = '.i'; 
        case {'f','a','front','anterior'},  view([-180   0]); viewname = '.a';
        case {'p','back','posterior'},      view([   0   0]); viewname = '.p';
        otherwise
          if isnumeric(job.view) && size(job.view)==2
            view(job.view); viewname = sprintf('.%04dx%04d',mod(job.view,360));
          else
            error('Unknown view.\n')
          end
      end
    end    
    
    
    
    
    %% print
    if job.imgprint.do 
      %%
      pfname = fullfile(sinfo(i).pp,sprintf('%s.%s%s.%s',[sinfo(i).ff,viewname,job.imgprint.type(3:end)]));
      print(h.figure , job.imgprint.type , job.imgprint.fdpi(job.imgprint.dpi) , pfname ); 
      
      if job.imgprint.close
        close(h.figure);
      end
    end
    
    if nargout>0
      varargout{1}{i} = h;
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




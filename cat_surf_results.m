function cat_surf_results(action,varargin)

%cat_surf_results to visualize results based on log P-maps
%
%_______________________________________________________________________
% Christian Gaser
% $Id: cat_surf_results.m 938 2016-05-19 08:35:43Z gaser $

global H

%-Input parameters
%--------------------------------------------------------------------------
if ~nargin, action = 'Disp'; end

pos = cell(2,1);

if ~ischar(action)
    varargin = {action varargin{:}};
    action   = 'Disp';
end

varargout = {[]};
H.clip = [];

%-Action
%--------------------------------------------------------------------------
switch lower(action)
    
    %-Display
    %======================================================================
    case 'disp'

        % positions & font size
        ws = spm('Winsize','Graphics');
        FS = spm('FontSizes');
        
        % figure and 5 views
        H.pos{1} = struct(...
            'fig',   [10  10  2*ws(3) ws(3)],...   % figure
            'cbar',  [0.400 0.550 0.200 0.300],... % colorbar for correlation matrix
            'view1', [0.075 0.450 0.325 0.325],... % surface view
            'view2', [0.075 0.050 0.325 0.325],... % surface view
            'view3', [0.600 0.450 0.325 0.325],... % surface view
            'view4', [0.600 0.050 0.325 0.325],... % surface view
            'view5', [0.300 0.200 0.400 0.400]);   % surface view   

        H.pos{2} = struct(...
            'fig',   [2*ws(3)+10  10  0.3*ws(3) ws(3)],...   % figure
            'left',  [0.100 0.925 0.800 0.050],... % select left hemisphere
            'right', [0.100 0.875 0.800 0.050],... % select right hemisphere
            'surf',  [0.100 0.805 0.800 0.050],... % 
            'thresh',[0.100 0.750 0.800 0.050],... % 
            'tview', [0.100 0.700 0.800 0.050],... % 
            'nocbar',[0.100 0.650 0.800 0.050],... % 
            'inv',   [0.100 0.600 0.800 0.050],... % 
            'ovmin', [0.100 0.400 0.800 0.150],... % 
            'ovmax', [0.100 0.250 0.800 0.150],... % 
            'save',  [0.100 0.100 0.800 0.050],... % 
            'close', [0.100 0.050 0.800 0.050],... % close button
            'text',  [0.100 0.150 0.800 0.200]);   % textbox   

        % 4 views
        H.pos{3} = struct(...
            'cbar',  [0.400 0.550 0.200 0.300],... % colorbar for correlation matrix
            'view1', [0.150 0.450 0.325 0.325],... % surface view
            'view2', [0.150 0.050 0.325 0.325],... % surface view
            'view3', [0.525 0.450 0.325 0.325],... % surface view
            'view4', [0.525 0.050 0.325 0.325]);   % surface view   

        % create figures
        for i=1:2
          H.figure(i) = figure(i+1);
          clf(H.figure(i));
        
          set(H.figure(i),'MenuBar','none','Position',H.pos{i}.fig,...
            'Name','Results','NumberTitle','off');
        end
        
        % define S
        H.S{1}.name  = ''; H.S{1}.side = 'lh';
        H.S{2}.name = '';  H.S{2}.side = 'rh';
        
        % add button for closing all windows
        H.close = uicontrol(H.figure(2),...
                'string','Close','Units','normalized',...
                'position',H.pos{2}.close,...
                'style','Pushbutton','HorizontalAlignment','center',...
                'callback','for i=2:26, try close(i); end; end;',...
                'ToolTipString','Close windows',...
                'Interruptible','on','Enable','on');
        
        % select results
        H.left = uicontrol(H.figure(2),...
                'string','Select left hemisphere data','Units','normalized',...
                'position',H.pos{2}.left,...
                'style','Pushbutton','HorizontalAlignment','center',...
                'callback',{@select_data,1},...
                'ToolTipString','Select resulst for left hemisphere',...
                'Interruptible','on','Enable','on');
        
        H.right = uicontrol(H.figure(2),...
                'string','Select right hemisphere data','Units','normalized',...
                'position',H.pos{2}.right,...
                'style','Pushbutton','HorizontalAlignment','center',...
                'callback',{@select_data,2},...
                'ToolTipString','Select resulst for right hemisphere',...
                'Interruptible','on','Enable','on');
        
        str  = { 'Underlying surface...','central','inflated','Dartel'};
        tmp  = { {@select_surf, 1},...
                 {@select_surf, 2},...
                 {@select_surf, 3}};
        
        H.surf = uicontrol(H.figure(2),...
                'string',str,'Units','normalized',...
                'position',H.pos{2}.surf,'UserData',tmp,...
                'style','PopUp','HorizontalAlignment','center',...
                'callback','spm(''PopUpCB'',gcbo)',...
                'ToolTipString','Underlying surface',...
                'Interruptible','on','Enable','off');

        str  = { 'Threshold...','P<0.05','P<0.01','P<0.001'};
        tmp  = { {@select_thresh, 1.3},...
                 {@select_thresh, 2},...
                 {@select_thresh, 3}};
        
        H.thresh = uicontrol(H.figure(2),...
                'string',str,'Units','normalized',...
                'position',H.pos{2}.thresh,'UserData',tmp,...
                'style','PopUp','HorizontalAlignment','center',...
                'callback','spm(''PopUpCB'',gcbo)',...
                'ToolTipString','Threshold',...
                'Interruptible','on','Visible','off');
        
        H.tview = uicontrol(H.figure(2),...
                'string','Disable top view','Units','normalized',...
                'position',H.pos{2}.tview,...
                'style','CheckBox','HorizontalAlignment','center',...
                'callback',{@checkbox_tview},...
                'ToolTipString','Disable top view in the image center',...
                'Interruptible','on','Visible','off');

        H.inv = uicontrol(H.figure(2),...
                'string','Invert results','Units','normalized',...
                'position',H.pos{2}.inv,...
                'style','CheckBox','HorizontalAlignment','center',...
                'callback',{@checkbox_inv},...
                'ToolTipString','Invert results',...
                'Interruptible','on','Visible','off');

        H.nocbar = uicontrol(H.figure(2),...
                'string','Disable colorbar','Units','normalized',...
                'position',H.pos{2}.nocbar,...
                'style','CheckBox','HorizontalAlignment','center',...
                'callback',{@checkbox_nocbar},...
                'ToolTipString','Disable colorbar',...
                'Interruptible','on','Visible','off');

        H.save = uicontrol(H.figure(2),...
                'string','Save','Units','normalized',...
                'position',H.pos{2}.save,...
                'style','Pushbutton','HorizontalAlignment','center',...
                'callback',{@save_image},...
                'ToolTipString','Save png image',...
                'Interruptible','on','Enable','off');

        if nargin == 3
          H.S{1}.name = varargin{1};
          H.S{2}.name = varargin{2};
          H.logP = 1;
          
          for ind=1:2
            H.S{ind}.Y    = spm_data_read(spm_data_hdr_read(H.S{ind}.name));
            H.S{ind}.info = cat_surf_info(H.S{ind}.name,1); 
            
            % check whether name contains 'log' tha indicates a logP file
            for i=1:size(H.S{ind}.name,1)
              if isempty(strfind(H.S{ind}.info(i).ff,'log'))
                H.logP = 0;
              end
            end
          end
          H.disable_tview = 0;
          H.show_inv = 0;
          H.inverted = 0;
          H.disable_cbar = 0;

          display_results_all;

          set(H.surf,'Enable','on');
          set(H.save,'Enable','on');
          set(H.tview,'Visible','on');
          set(H.nocbar,'Visible','on');
          if min(min(H.S{1}.Y(:),H.S{2}.Y(:))) < 0
            set(H.inv,'Visible','on');
          end
        end

    %-ColourBar
    %======================================================================
    case {'colourbar', 'colorbar'}
        if isempty(varargin), varargin{1} = gca; end
        if length(varargin) == 1, varargin{2} = 'on'; end
        H = getHandles(varargin{1});
        d   = getappdata(H.patch(1),'data');
        col = getappdata(H.patch(1),'colourmap');
        if strcmpi(varargin{2},'off')
            if isfield(H,'colourbar') && ishandle(H.colourbar)
                delete(H.colourbar);
                H = rmfield(H,'colourbar');
                setappdata(H.axis,'handles',H);
            end
            return;
        end
        if isempty(d) || ~any(d(:)), varargout = {H}; return; end
        if isempty(col), col = jet(256); end
        if ~isfield(H,'colourbar') || ~ishandle(H.colourbar)
%            H.colourbar = colorbar('peer',gca,'NorthOutside');
            H.colourbar = colorbar('NorthOutside');
            set(H.colourbar,'Tag','');
            set(get(H.colourbar,'Children'),'Tag','');
        end
        c(1:size(col,1),1,1:size(col,2)) = col;
        ic = findobj(H.colourbar,'Type','image');
        clim = getappdata(H.patch(1), 'clim');
        if isempty(clim), clim = [false NaN NaN]; end
        
        if size(d,1) > size(d,2), d = d'; end

        % Update colorbar colors if clipping is used
        H.clip = getappdata(H.patch(1), 'clip')
        if ~isempty(H.clip)
            if ~isnan(H.clip(2)) && ~isnan(H.clip(3))
                ncol = length(col);
                col_step = (clim(3) - clim(2))/ncol;
                cmin = max([1,ceil((H.clip(2)-clim(2))/col_step)]);
                cmax = min([ncol,floor((H.clip(3)-clim(2))/col_step)]);
                col(cmin:cmax,:) = repmat([0.5 0.5 0.5],(cmax-cmin+1),1);
                c(1:size(col,1),1,1:size(col,2)) = col;
            end
        end
        if size(d,1) > 1
            set(ic,'CData',c(1:size(d,1),:,:));
            set(ic,'YData',[1 size(d,1)]);
            set(H.colourbar,'YLim',[1 size(d,1)]);
            set(H.colourbar,'YTickLabel',[]);
        else
            set(ic,'CData',c);
            clim = getappdata(H.patch(1),'clim');
            if isempty(clim), clim = [false min(d) max(d)]; end
            set(ic,'YData',clim(2:3));
            set(H.colourbar,'YLim',clim(2:3));
        end
        setappdata(H.axis,'handles',H);
        
    %-ColourMap
    %======================================================================
    case {'colourmap', 'colormap'}
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
        if length(varargin) == 1
            varargout = { getappdata(H.patch(1),'colourmap') };
            return;
        else
            setappdata(H.patch(1),'colourmap',varargin{2});
            d = getappdata(H.patch(1),'data');
            H = updateTexture(H,d);
        end
        if nargin>1
            colormap(varargin{2});
        end
        
    %-CLim
    %======================================================================
    case 'clim'
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
        if length(varargin) == 1
            c = getappdata(H.patch,'clim');
            if ~isempty(c), c = c(2:3); end
            varargout = { c };
            return;
        else
            if strcmp(varargin{2},'on') || isempty(varargin{2}) || any(~isfinite(varargin{2}))
                setappdata(H.patch,'clim',[false NaN NaN]);
            else
                setappdata(H.patch,'clim',[true varargin{2}]);
            end
            d = getappdata(H.patch,'data');
            H = updateTexture(H,d);
          
        end
        
        if nargin>1 && isnumeric(varargin{2}) && numel(varargin{2})==2
            caxis(H.axis,varargin{2});
        else
            caxis(H.axis,[min(d),max(d)])
        end
        
    %-CLip
    %======================================================================
    case 'clip'
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
        if length(varargin) == 1
            c = getappdata(H.patch,'clip');
            if ~isempty(c), c = c(2:3); end
            varargout = { c };
            return;
        else
            if isempty(varargin{2}) || any(~isfinite(varargin{2}))
                for ind = 1:5
                  setappdata(H.patch(ind),'clip',[false NaN NaN]);
                end
            else
                for ind = 1:5
                  setappdata(H.patch(ind),'clip',[true varargin{2}]);
                end
            end
            for ind = 1:5
              d = getappdata(H.patch,'data');
              H = updateTexture(H,ind,d);
            end
        end

    end       
        
%-----------------------------------------------------------------------
function H = select_thresh(thresh)
%-----------------------------------------------------------------------
global H

H.clip = [true -thresh thresh];
for ind=1:(5 - H.disable_tview)
  setappdata(H.patch(ind),'clip',H.clip);
  d = getappdata(H.patch(ind),'data');
  H = updateTexture(H,ind,d);
end

%-----------------------------------------------------------------------
function H = select_surf(surf)
%-----------------------------------------------------------------------
global H

for ind=1:2
  switch surf
  case 1
    H.S{ind}.info(1).Pmesh = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces',[H.S{ind}.info(1).side '.central.freesurfer.gii']);
  case 2
    H.S{ind}.info(1).Pmesh = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces',[H.S{ind}.info(1).side '.inflated.freesurfer.gii']);
  case 3
    H.S{ind}.info(1).Pmesh = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces',[H.S{ind}.info(1).side '.central.Template_T1_IXI555_MNI152.gii']);
  end
end

display_results_all;
for ind = 1:(5 - H.disable_tview)
  setappdata(H.patch(ind),'clip',[true NaN NaN]);
  d = getappdata(H.patch(ind),'data');
  H = updateTexture(H,ind,d);
end

%-----------------------------------------------------------------------
function display_results_all(obj, event_obj)
%-----------------------------------------------------------------------
global H

if (size(H.S{1}.Y) > 1 | size(H.S{2}.Y) > 1) & min(min(H.S{1}.Y(:),H.S{2}.Y(:))) < 0
  disp('Warning: Only results with positive values are displayed!');
end

% clear larger area and set background color to update labels and title
H.axis = axes('Parent',H.figure(1),'Position',[-.1 -.1 1.1 1.1],'Color',[1 1 1]);
cla(H.axis);

H.renderer = get(H.figure(1),'Renderer');
set(H.figure(1),'Renderer','OpenGL');

%-Compute mesh curvature
%------------------------------------------------------------------
g = gifti(fullfile(spm('dir'),'toolbox','cat12','templates_surfaces',[H.S{1}.info(1).side '.mc.central.freesurfer.gii']));
H.S{1}.curv = g.cdata;
g = gifti(fullfile(spm('dir'),'toolbox','cat12','templates_surfaces',[H.S{2}.info(1).side '.mc.central.freesurfer.gii']));
H.S{2}.curv = g.cdata;

if ~H.disable_tview
  ind = 1;
  display_results(5, H.pos{1}.view5, [ 0 90]);
else
  ind = 3;
end
display_results(1, H.pos{ind}.view1, [ 90 0]);
display_results(2, H.pos{ind}.view2, [-90 0]);
display_results(3, H.pos{ind}.view3, [-90 0]);
display_results(4, H.pos{ind}.view4, [ 90 0]);

% add colorbar
%H.cbar = axes('Position',H.pos{1}.cbar,'Parent',H.figure);
%cat_surf_results('Colourbar', H); 

H.S{1}.thresh = min(H.S{1}.Y(H.S{1}.Y(:)>0));
H.S{1}.thresh = min(H.S{1}.thresh,min(H.S{2}.Y(H.S{2}.Y(:)>0)));

H.S{1}.min = min(min(H.S{1}.Y(:),H.S{2}.Y(:)));
H.S{1}.max = max(max(H.S{1}.Y(:),H.S{2}.Y(:)));

if H.S{1}.min < 0
  mnx = max(abs([H.S{1}.min,H.S{1}.max]));
  H.S{1}.min = -mnx;
  H.S{1}.max =  mnx;
end

H.S{1}.max = round(1.1*H.S{1}.max);
if H.S{1}.min < 0
  H.S{1}.min = round(1.1*H.S{1}.min);
else
  H.S{1}.min = round(0.9*H.S{1}.min);
end

if H.logP 
  set(H.thresh,'Visible','on');
end

for ind=1:(5 - H.disable_tview)
  setappdata(H.patch(ind), 'clim', [true H.S{1}.min H.S{1}.max]);
  d = getappdata(H.patch(ind),'data');
  H = updateTexture(H,ind,d);
end

if ~H.disable_cbar
  H = show_colorbar(H);
end

% show slider for range of results
if size(d,1)==1

  % allow slider a more extended range
  mnx = 1.5*max(abs([H.S{1}.min H.S{1}.max]));

  sliderPanel(...
        'Parent'  , H.figure(2), ...
        'Title'   , 'Result min', ...
        'Position', H.pos{2}.ovmin, ...
        'Backgroundcolor', [0.8 0.8 0.8],...
        'Min'     , -mnx, ...
        'Max'     , mnx, ...
        'Value'   , H.S{1}.min, ...
        'FontName', 'Verdana', ...
        'FontSize', 8, ...
        'NumFormat', '%f', ...
        'Callback', @slider_clim_min);

  sliderPanel(...
        'Parent'  , H.figure(2), ...
        'Title'   , 'Result max', ...
        'Position', H.pos{2}.ovmax, ...
        'Backgroundcolor', [0.8 0.8 0.8],...
        'Min'     , -mnx, ...
        'Max'     , mnx, ...
        'Value'   , H.S{1}.max, ...
        'FontName', 'Verdana', ...
        'FontSize', 8, ...
        'NumFormat', '%f', ...
        'Callback', @slider_clim_max);
end

%-----------------------------------------------------------------------
function H = show_colorbar(H)
%-----------------------------------------------------------------------

% show colorbar
figure(H.figure(1))
if numel(H.S{1}.info) == 1
  H.cbar = axes('Parent',H.figure(1),'Position',H.pos{1}.cbar,'Color',[0.5 0.5 0.5],'Visible','off');
  clim = getappdata(H.patch(1), 'clim');
  axis(H.cbar,'off'); caxis([clim(2) clim(3)]);
  if H.logP, title('p-value');end
  H.colourbar = colorbar('peer',H.cbar,'Northoutside');
  colormap(getappdata(H.patch(1),'col'));
  
  if H.logP
    XTick = get(H.colourbar,'XTick');

    XTickLabel = [];
    for i=1:length(XTick)
      if XTick(i) > 0
        XTickLabel = char(XTickLabel,remove_zeros(sprintf('%.g',10^(-XTick(i)))));
      elseif XTick(i) < 0
        XTickLabel = char(XTickLabel,remove_zeros(sprintf('-%.g',10^(XTick(i)))));
      else
        XTickLabel = char(XTickLabel,'');
      end
    end
    set(H.colourbar,'XTickLabel',XTickLabel(2:end,:),'XTick',XTick);
  end
else
  H.cbar = axes('Parent',H.figure(1),'Position',H.pos{1}.cbar+[0 0.10 0 0],'Color',[0.5 0.5 0.5],'Visible','off');
  if numel(H.S{1}.info) ==3
    cb = [7 1 1 4 2 2 7;...
          7 1 6 7 5 2 7;...
          7 7 3 3 3 7 7];
  else
    cb = [7 1 1 4 2 2 7;...
          7 1 1 4 2 2 7];
  end
  imagesc(cb);
  colormap([1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 1 1 1])
  axis(H.cbar,'off'); axis('image');  
end

%-----------------------------------------------------------------------
function display_results(ind, win, vw)
%-----------------------------------------------------------------------
global H

if ind < 5
  M  = gifti(H.S{round(ind/2)}.info(1).Pmesh);
  Mc.cdata = H.S{round(ind/2)}.Y;
else
  Ml = gifti(H.S{1}.info(1).Pmesh);
  Mr = gifti(H.S{2}.info(1).Pmesh);
  Mcl.cdata = H.S{1}.Y;
  Mcr.cdata = H.S{2}.Y;
  M.faces = [Ml.faces; Mr.faces + size(Ml.vertices,1)];
  M.vertices = [Ml.vertices; Mr.vertices];
  M.mat = Ml.mat;
  Mc.cdata = [Mcl.cdata; Mcr.cdata];
end

if isfield(Mc,'cdata')
  M.cdata = Mc.cdata;
else
  M.cdata = []; 
end

H.axis = axes('Position',win,'Parent',H.figure(1),'Visible','off');
H.figure(1) = ancestor(H.axis,'figure');
figure(H.figure(1)); axes(H.axis);

if isfield(M,'facevertexcdata')
  H.cdata = M.facevertexcdata;
else
  H.cdata = []; 
end

if ~isfield(M,'vertices') || ~isfield(M,'faces')
  error('cat_surf_results:nomesh','ERROR:cat_surf_render: No input mesh.');
end

%% -Patch
%------------------------------------------------------------------
P = struct('vertices',M.vertices, 'faces',double(M.faces));
H.patch(ind) = patch(P,...
            'FaceColor',        [0.6 0.6 0.6],...
            'EdgeColor',        'none',...
            'FaceLighting',     'gouraud',...
            'SpecularStrength', 0.7,...
            'AmbientStrength',  0.4,...
            'DiffuseStrength',  0.6,...
            'SpecularExponent', 10,...
            'Clipping',         'off',...
            'DeleteFcn',        {@myDeleteFcn, H.renderer},...
            'Visible',          'off',...
            'Tag',              'CATSurfRender',...
            'Parent',           H.axis);
setappdata(H.patch(ind),'patch',P);

%-Compute mesh curvature
%------------------------------------------------------------------
if ind < 5
  curv = H.S{round(ind/2)}.curv;
else
  curv = [H.S{1}.curv; H.S{2}.curv];
end

setappdata(H.patch(ind),'curvature',curv);

%-Apply texture to mesh
%------------------------------------------------------------------
if isfield(M,'facevertexcdata')
  T = M.facevertexcdata;
elseif isfield(M,'cdata')
  T = M.cdata;
else
  T = [];
end
H = updateTexture(H,ind,T);

axis(H.axis,'image');
axis(H.axis,'off');
view(H.axis,vw);
material(H.figure(1),'dull');

% default lighting
H.light(1) = camlight; set(H.light(1),'Parent',H.axis); 
if ismac
  % switch off local light (camlight)
  caml = findall(gcf,'Type','light','Style','local');     
  set(caml,'visible','off');
            
  % set inner light
  H.light(2) = light('Position',[0 0 0]); 
  set(H.patch(ind),'BackFaceLighting','unlit');
end
        
setappdata(H.axis,'handles',H);
set(H.patch(ind),'Visible','on');
camlight(H.light(1))

%==========================================================================
function [H, C] = updateTexture(H,ind,v,col)

%-Project data onto surface mesh
%--------------------------------------------------------------------------
if size(v,2) < size(v,1)
  v = v';
end
v(isinf(v)) = NaN;

setappdata(H.patch(ind),'data',v);

%-Get colourmap
%--------------------------------------------------------------------------
if ~exist('col','var')
  if size(v,1) == 1
    col = jet(256);
  else
    % use RGB colormap
    col = zeros(256,3,size(v,1));
    for i=1:3
      col(:,i,i) = 1;
    end
  end
end

setappdata(H.patch(ind),'col',col);

if ~exist('FaceColor','var') || isempty(FaceColor), FaceColor = 'interp'; end
%setappdata(H.colourmap,'colourmap',col);

%-Get curvature
%--------------------------------------------------------------------------
curv = getappdata(H.patch(ind),'curvature');

if size(curv,2) == 1
    th = 0.15;
    curv((curv<-th)) = -1.5*th;
    curv((curv>th))  =  0.1*th;
    curv = 0.5*(curv + th)/(2*th);
    curv = 0.5 + repmat(curv,1,3);
end

%-Create RGB representation of data according to colourmap
%--------------------------------------------------------------------------
C = zeros(size(v,2),3);
clim = getappdata(H.patch(ind), 'clim');
if isempty(clim), clim = [false NaN NaN]; end
mi = clim(2); ma = clim(3);
if any(v(:))
    if ~clim(1), mi = min(v(:)); ma = max(v(:)); end
    % don't allow negative values for multiple maps
    if size(v,1) > 1 & mi < 0
      if ~isempty(H.clip)
        H.clip(2) = -Inf;
      else
        H.clip = [true -Inf 0];
      end
    end
    for i=1:size(v,1)
        C = C + squeeze(ind2rgb(floor(((v(i,:)-mi)/(ma-mi))*size(col,1)),col(:,:,i)));
    end
end

%H.clip = getappdata(H.patch(ind), 'clip');
if ~isempty(H.clip)
    v(v>H.clip(2) & v<H.clip(3)) = NaN;
    setappdata(H.patch(ind), 'clip', [true H.clip(2) H.clip(3)]);
end

setappdata(H.patch(ind), 'clim', [true mi ma]);

%-Build texture by merging curvature and data
%--------------------------------------------------------------------------
if size(v,1) > 1 % RGB
  for i=1:size(v,1)
    C(:,i) = any(v(i,:),1)' .* C(:,i);
  end
else
  C = repmat(any(v,1),3,1)' .* C;
end

% replace regions below threshold by curvature
ind0 = repmat(~any(v,1),3,1)';
C(ind0) = curv(ind0);

set(H.patch(ind), 'FaceVertexCData',C, 'FaceColor',FaceColor);

%-Update the colourbar
%--------------------------------------------------------------------------
if isfield(H,'colourbar')
%    H = cat_surf_results('Colourbar',H);
end

%-----------------------------------------------------------------------
function select_data(obj, event_obj, ind)
%-----------------------------------------------------------------------
global H

if isempty(H.S{2}.name)
  str = 'log';
else
  str = 'log';
end

side = '';
str_side = {'left','right'};

H.logP = 1;

while ~strcmp(side,H.S{ind}.side)
  H.S{ind}.name = spm_select([1 3],'mesh',['Select log P map for ' str_side{ind} ' hemisphere'],'','',str);
  H.S{ind}.Y    = spm_data_read(spm_data_hdr_read(H.S{ind}.name));
  H.S{ind}.info = cat_surf_info(H.S{ind}.name,1); 
  side = H.S{ind}.side;
  for i=1:size(H.S{ind}.name,1)
    % check whether name contains 'log' tha indicates a logP file
    if isempty(strfind(H.S{ind}.info(i).ff,'log'))
      H.logP = 0;
    end
    if ~strcmp(H.S{ind}.info(i).side,H.S{ind}.side)
      fprintf('%s does not contain %s hemisphere data.\n',H.S{ind}.name(i,:),str_side{ind});
      side = '';
    end
  end
end

H.show_inv = 0;
H.disable_tview = 0;
H.inverted = 0;
H.disable_cbar = 0;
H.clip = [false NaN NaN];

% enable display button if both sides are defined
if ~isempty(H.S{1}.name) & ~isempty(H.S{2}.name)
  display_results_all;
  set(H.surf,'Enable','on');
  set(H.save,'Enable','on');
  set(H.tview,'Visible','on');
  set(H.nocbar,'Visible','on');
  if min(min(H.S{1}.Y(:),H.S{2}.Y(:))) < 0
    set(H.inv,'Visible','on');
  end
end

%==========================================================================
function save_image(obj,event_obj,filename)

global H
  %%
  
  if ~exist('filename','var')

    nm = H.S{1}.info(1).ff;
    filename = [nm '.png'];
    
    % end with _0???.ext?
    if length(nm) > 4
      if strcmp(nm(length(nm)-4:length(nm)-3),'_0') 
    
        SPM_name = fullfile(H.S{1}.info(1).pp, 'SPM.mat');
    
        % SPM.mat exist?
        if exist(SPM_name,'file')
          load(SPM_name);
          xCon = SPM.xCon;
          Ic = str2double(nm(length(nm)-3:length(nm)));
          str_num = deblank(xCon(Ic).name);

          % replace spaces with "_" and characters like "<" or ">" with "gt" or "lt"
          str_num(strfind(str_num,' ')) = '_';
          strpos = strfind(str_num,' > ');
          if ~isempty(strpos), str_num = [str_num(1:strpos-1) '_gt_' str_num(strpos+1:end)]; end
          strpos = strfind(str_num,' < ');
          if ~isempty(strpos), str_num = [str_num(1:strpos-1) '_lt_' str_num(strpos+1:end)]; end
          strpos = strfind(str_num,'>');
          if ~isempty(strpos), str_num = [str_num(1:strpos-1) 'gt' str_num(strpos+1:end)]; end
          strpos = strfind(str_num,'<');
          if ~isempty(strpos), str_num = [str_num(1:strpos-1) 'lt' str_num(strpos+1:end)]; end
          str_num = spm_str_manip(str_num,'v');
        
          if ~isempty(H.clip)
            if isnan(H.clip(3))
              str_thresh = '_';
            else
              str_thresh = sprintf('P%g_',round(1000*10^(-H.clip(3)))/10);
            end
          else
            str_thresh = '_';
          end
          filename = ['logP_' str_thresh str_num '.png'];
        end
      end
    end

    filename = uiputfile({...
      '*.png' 'PNG files (*.png)'}, 'Save as', filename);
  else
    [pth,nam,ext] = fileparts(filename);
    if isempty(pth), pth = cd; end
    if ~strcmp({'.gii','.png'},ext), nam = [nam ext]; end
    if isempty(nam)
      filename = uiputfile({...
        '*.png' 'PNG files (*.png)'}, 'Save as',nam);
    else
      filename = fullfile(pth,[nam '.png']);
    end
  end
  
  try
    FS = get(H.colourbar,'FontSize');
    FU = get(H.colourbar,'FontUnits');
 
    set(H.colourbar,'FontUnits','Normalized','FontSize',0.1);
    changed_font = 1;
  catch
    changed_font = 0;
  end
  
  if isdeployed
      deployprint(H.figure(1), '-dpng', '-opengl', filename);
  else
      print(H.figure(1), '-dpng', '-r300', '-opengl',filename);
  end

  if changed_font
    set(H.colourbar,'FontUnits',FU,'FontSize',FS);
  end

%==========================================================================
function slider_clim_min(hObject, evt)
global H

figure(H.figure(1))
val = get(hObject, 'Value');
c = getappdata(H.patch(1),'clim');
for ind = 1:(5 - H.disable_tview)
  setappdata(H.patch(ind),'clim',[true val c(3)]);
  d = getappdata(H.patch(ind),'data');
  H = updateTexture(H,ind,d);
end

% delete colorbar and display again
set(get(H.cbar,'Title'),'Visible','off')
if size(d,1) == 1
  set(H.colourbar,'Visible','off');
  if ~H.disable_cbar
    H = show_colorbar(H);
  end
end

%==========================================================================
function slider_clim_max(hObject, evt)
global H

figure(H.figure(1))
val = get(hObject, 'Value');
c = getappdata(H.patch(1),'clim');
for ind = 1:(5 - H.disable_tview)
  setappdata(H.patch(ind),'clim',[true c(2) val]);
  d = getappdata(H.patch(ind),'data');
  H = updateTexture(H,ind,d);
end

% delete colorbar and display again
set(get(H.cbar,'Title'),'Visible','off')
if size(d,1) == 1
  set(H.colourbar,'Visible','off');
  if ~H.disable_cbar
    H = show_colorbar(H);
  end
end

%==========================================================================
function checkbox_tview(obj, event_obj)
global H
  
H.disable_tview = get(H.tview,'Value');
display_results_all;

%==========================================================================
function checkbox_inv(obj, event_obj)
global H
  
H.show_inv = get(H.inv,'Value');

if H.show_inv & ~H.inverted
  H.S{1}.Y = -H.S{1}.Y;
  H.S{2}.Y = -H.S{2}.Y;
  H.inverted = 1;
end

if ~H.show_inv & H.inverted
  H.S{1}.Y = -H.S{1}.Y;
  H.S{2}.Y = -H.S{2}.Y;
  H.inverted = 0;
end
display_results_all;

%==========================================================================
function checkbox_nocbar(obj, event_obj)
global H
  
H.disable_cbar = get(H.nocbar,'Value');

if H.disable_cbar
  d = getappdata(H.patch(1),'data');

  % delete colorbar and title
  if size(d,1) == 1
    set(get(H.cbar,'Title'),'Visible','off')
    set(H.colourbar,'Visible','off')
  else % delete only axis
    cla(H.cbar);
  end
else
  H = show_colorbar(H);
end

%==========================================================================
function H = getHandles(H)
if ~nargin || isempty(H), H = gca; end
if ishandle(H) && ~isappdata(H,'handles')
    a = H; clear H;
    H.axis     = a;
    H.figure(1)   = ancestor(H.axis,'figure');
    H.patch    = findobj(H.axis,'type','patch');
    H.light    = findobj(H.axis,'type','light');
    H.rotate3d = rotate3d(H.figure(1));
    setappdata(H.axis,'handles',H);
elseif ishandle(H)
    H = getappdata(H,'handles');
else
    H = getappdata(H.axis,'handles');
end

%==========================================================================
function myDeleteFcn(obj,evt,renderer)
try rotate3d(get(obj,'parent'),'off'); end
set(ancestor(obj,'figure'),'Renderer',renderer);

%==========================================================================
function s=remove_zeros(s)

pos = length(s);
while pos>1
  if strcmp(s(pos),'0')
    s(pos)='';
    pos = pos-1;
  else break
  end
end

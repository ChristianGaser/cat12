function cat_surf_results(action,varargin)

%cat_surf_results to visualize results based on log P-maps
%
%_______________________________________________________________________
% Christian Gaser
% $Id: cat_surf_results.m 938 2016-05-19 08:35:43Z gaser $

global pos H S clip

%-Input parameters
%--------------------------------------------------------------------------
if ~nargin, action = 'Disp'; end

pos = cell(2,1);

if ~ischar(action)
    varargin = {action varargin{:}};
    action   = 'Disp';
end

varargout = {[]};
clip = [];

%-Action
%--------------------------------------------------------------------------
switch lower(action)
    
    %-Display
    %======================================================================
    case 'disp'

        % positions & font size
        ws = spm('Winsize','Graphics');
        FS = spm('FontSizes');
        
        pos{1} = struct(...
            'fig',   [10  10  2*ws(3) ws(3)],...   % figure
            'cbar',  [0.400 0.550 0.200 0.300],... % colorbar for correlation matrix
            'view1', [0.075 0.450 0.325 0.325],... % surface view
            'view2', [0.075 0.050 0.325 0.325],... % surface view
            'view3', [0.600 0.450 0.325 0.325],... % surface view
            'view4', [0.600 0.050 0.325 0.325],... % surface view
            'view5', [0.300 0.200 0.400 0.400]);   % surface view   


        pos{2} = struct(...
            'fig',   [2*ws(3)+10  10  0.3*ws(3) ws(3)],...   % figure
            'left',  [0.100 0.925 0.800 0.050],... % select left hemisphere
            'right', [0.100 0.875 0.800 0.050],... % select right hemisphere
            'surf',  [0.100 0.805 0.800 0.050],... % 
            'thresh',[0.100 0.765 0.800 0.050],... % 
            'save',  [0.100 0.725 0.800 0.050],... % 
            'close', [0.100 0.600 0.800 0.050],... % close button
            'text',  [0.100 0.550 0.800 0.200]);   % textbox   

        % create figures
        for i=1:2
          H.figure(i) = figure(i+1);
          clf(H.figure(i));
        
          set(H.figure(i),'MenuBar','none','Position',pos{i}.fig,...
            'Name','Results','NumberTitle','off');
        end
        
        % define S
        S{1}.name  = ''; S{1}.side = 'lh';
        S{2}.name = '';  S{2}.side = 'rh';
        
        % add button for closing all windows
        H.close = uicontrol(H.figure(2),...
                'string','Close','Units','normalized',...
                'position',pos{2}.close,...
                'style','Pushbutton','HorizontalAlignment','center',...
                'callback','for i=2:26, try close(i); end; end;',...
                'ToolTipString','Close windows',...
                'Interruptible','on','Enable','on');
        
        % select results
        H.left = uicontrol(H.figure(2),...
                'string','Select left hemisphere data','Units','normalized',...
                'position',pos{2}.left,...
                'style','Pushbutton','HorizontalAlignment','center',...
                'callback',{@select_data,1},...
                'ToolTipString','Select resulst for left hemisphere',...
                'Interruptible','on','Enable','on');
        
        H.right = uicontrol(H.figure(2),...
                'string','Select right hemisphere data','Units','normalized',...
                'position',pos{2}.right,...
                'style','Pushbutton','HorizontalAlignment','center',...
                'callback',{@select_data,2},...
                'ToolTipString','Select resulst for right hemisphere',...
                'Interruptible','on','Enable','on');
        
        str  = { 'Threshold...','P<0.05','P<0.01','P<0.001'};
        tmp  = { {@select_thresh, 1.3},...
                 {@select_thresh, 2},...
                 {@select_thresh, 3}};
        
        H.thresh = uicontrol(H.figure(2),...
                'string',str,'Units','normalized',...
                'position',pos{2}.thresh,'UserData',tmp,...
                'style','PopUp','HorizontalAlignment','center',...
                'callback','spm(''PopUpCB'',gcbo)',...
                'ToolTipString','Threshold',...
                'Interruptible','on','Visible','off');
        
        str  = { 'Underlying surface...','central','inflated','Dartel'};
        tmp  = { {@select_surf, 1},...
                 {@select_surf, 2},...
                 {@select_surf, 3}};
        
        H.surf = uicontrol(H.figure(2),...
                'string',str,'Units','normalized',...
                'position',pos{2}.surf,'UserData',tmp,...
                'style','PopUp','HorizontalAlignment','center',...
                'callback','spm(''PopUpCB'',gcbo)',...
                'ToolTipString','Underlying surface',...
                'Interruptible','on','Enable','off');

        H.save = uicontrol(H.figure(2),...
                'string','Save','Units','normalized',...
                'position',pos{2}.save,...
                'style','Pushbutton','HorizontalAlignment','center',...
                'callback',{@save_image},...
                'ToolTipString','Save png image',...
                'Interruptible','on','Enable','off');

        if nargin == 3
          S{1}.name = varargin{1};
          S{2}.name = varargin{2};
          for ind=1:2
            S{ind}.Y    = spm_data_read(spm_data_hdr_read(S{ind}.name));
            S{ind}.info = cat_surf_info(S{ind}.name,1); 
          end
          set(H.surf,'Enable','on');
          set(H.save,'Enable','on');
          display_results_all;
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
        clip = getappdata(H.patch(1), 'clip')
        if ~isempty(clip)
            if ~isnan(clip(2)) && ~isnan(clip(3))
                ncol = length(col);
                col_step = (clim(3) - clim(2))/ncol;
                cmin = max([1,ceil((clip(2)-clim(2))/col_step)]);
                cmax = min([ncol,floor((clip(3)-clim(2))/col_step)]);
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
global H clip

clip = [true -thresh thresh];
for ind=1:5
  try
    setappdata(H.patch(ind),'clip',clip);
    d = getappdata(H.patch(ind),'data');
    H = updateTexture(H,ind,d);
  end
end

%-----------------------------------------------------------------------
function H = select_surf(surf)
%-----------------------------------------------------------------------
global H S

for ind=1:2
  switch surf
  case 1
    S{ind}.info(1).Pmesh = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces',[S{ind}.info(1).side '.central.freesurfer.gii']);
  case 2
    S{ind}.info(1).Pmesh = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces',[S{ind}.info(1).side '.inflated.freesurfer.gii']);
  case 3
    S{ind}.info(1).Pmesh = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces',[S{ind}.info(1).side '.central.Template_T1_IXI555_MNI152.gii']);
  end
end

display_results_all;
for ind = 1:5
  setappdata(H.patch(ind),'clip',[true NaN NaN]);
  d = getappdata(H.patch(ind),'data');
  H = updateTexture(H,ind,d);
end

%-----------------------------------------------------------------------
function display_results_all(obj, event_obj)
%-----------------------------------------------------------------------
global pos H S

% clear larger area and set background color to update labels and title
H.axis = axes('Parent',H.figure(1),'Position',[-.1 -.1 1.1 1.1],'Color',[1 1 1]);
cla(H.axis);

H.renderer = get(H.figure(1),'Renderer');
set(H.figure(1),'Renderer','OpenGL');

%-Compute mesh curvature
%------------------------------------------------------------------
g = gifti(fullfile(spm('dir'),'toolbox','cat12','templates_surfaces',[S{1}.info(1).side '.mc.central.freesurfer.gii']));
S{1}.curv = g.cdata;
g = gifti(fullfile(spm('dir'),'toolbox','cat12','templates_surfaces',[S{2}.info(1).side '.mc.central.freesurfer.gii']));
S{2}.curv = g.cdata;

display_results(1, pos{1}.view1, [ 90 0]);
display_results(2, pos{1}.view2, [-90 0]);
display_results(3, pos{1}.view3, [-90 0]);
display_results(4, pos{1}.view4, [ 90 0]);
display_results(5, pos{1}.view5, [ 0 90]);

% add colorbar
%H.cbar = axes('Position',pos{1}.cbar,'Parent',H.figure);
%cat_surf_results('Colourbar', H); 

S{1}.thresh = min(S{1}.Y(S{1}.Y(:)>0));
S{1}.thresh = min(S{1}.thresh,min(S{2}.Y(S{2}.Y(:)>0)));

S{1}.min = min(min(S{1}.Y(:),S{2}.Y(:)));
S{1}.max = max(max(S{1}.Y(:),S{2}.Y(:)));

if S{1}.min < 0
  mnx = max(abs([S{1}.min,S{1}.max]));
  S{1}.min = -mnx;
  S{1}.max =  mnx;
end

if S{1}.thresh < 0.1 
  set(H.thresh,'Visible','on');
end

for ind=1:5
  setappdata(H.patch(ind), 'clim', [true S{1}.min S{1}.max]);
  d = getappdata(H.patch(ind),'data');
  H = updateTexture(H,ind,d);
end

if numel(S{1}.info) == 1
  H.cbar = axes('Parent',H.figure(1),'Position',pos{1}.cbar,'Color',[0.5 0.5 0.5],'Visible','off');
  axis(H.cbar,'off'); caxis([S{1}.min,S{1}.max]); title('p-value');
  H.colourbar = colorbar('peer',H.cbar,'Northoutside');
  colormap(getappdata(H.patch(ind),'col'));
  
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

%-----------------------------------------------------------------------
function display_results(ind, win, vw)
%-----------------------------------------------------------------------
global H S clip

if ind < 5
  M  = gifti(S{round(ind/2)}.info(1).Pmesh);
  Mc.cdata = S{round(ind/2)}.Y;
else
  Ml = gifti(S{1}.info(1).Pmesh);
  Mr = gifti(S{2}.info(1).Pmesh);
  Mcl.cdata = S{1}.Y;
  Mcr.cdata = S{2}.Y;
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
  curv = S{round(ind/2)}.curv;
else
  curv = [S{1}.curv; S{2}.curv];
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

global clip

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
    for i=1:size(v,1)
      col(:,i,i) = (1:256)/256;
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
    curv((curv<-th)) = -th;
    curv((curv>th))  =  th;
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
    for i=1:size(v,1)
        C = C + squeeze(ind2rgb(floor(((v(i,:)-mi)/(ma-mi))*size(col,1)),col(:,:,i)));
    end
end

%clip = getappdata(H.patch(ind), 'clip');
if ~isempty(clip)
    v(v>clip(2) & v<clip(3)) = NaN;
    setappdata(H.patch(ind), 'clip', [true clip(2) clip(3)]);
end

setappdata(H.patch(ind), 'clim', [true mi ma]);

%-Build texture by merging curvature and data
%--------------------------------------------------------------------------
C = repmat(~any(v,1),3,1)' .* curv + repmat(any(v,1),3,1)' .* C;

set(H.patch(ind), 'FaceVertexCData',C, 'FaceColor',FaceColor);

%-Update the colourbar
%--------------------------------------------------------------------------
if isfield(H,'colourbar')
%    H = cat_surf_results('Colourbar',H);
end

%-----------------------------------------------------------------------
function select_data(obj, event_obj, ind)
%-----------------------------------------------------------------------
global H S

if isempty(S{2}.name)
  str = 'log';
else
  str = 'log';
end

side = '';
str_side = {'left','right'};

while ~strcmp(side,S{ind}.side)
  S{ind}.name = spm_select([1 3],'mesh',['Select log P map for ' str_side{ind} ' hemisphere'],'','',str);
  S{ind}.Y    = spm_data_read(spm_data_hdr_read(S{ind}.name));
  S{ind}.info = cat_surf_info(S{ind}.name,1); 
  side = S{ind}.side;
  for i=1:size(S{ind}.name,1)
    if ~strcmp(S{ind}.info(i).side,S{ind}.side)
      fprintf('%s does not contain %s hemisphere data.\n',S{ind}.name(i,:),str_side{ind});
      side = '';
    end
  end
end

% enable display button if both sides are defined
if ~isempty(S{1}.name) & ~isempty(S{2}.name)
  set(H.surf,'Enable','on');
  set(H.save,'Enable','on');
  display_results_all;
end

%==========================================================================
function save_image(obj,event_obj,filename)
global H
  %%
  
  if ~exist('filename','var')
    filename = uiputfile({...
      '*.png' 'PNG files (*.png)'}, 'Save as');
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

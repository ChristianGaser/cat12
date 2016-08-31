function cat_surf_results(action,varargin)

%cat_surf_results to visualize results based on log P-maps
%
%_______________________________________________________________________
% Christian Gaser
% $Id$

global H

%-Input parameters
%--------------------------------------------------------------------------
if ~nargin, action = 'Disp'; end

if ~ischar(action)
    varargin = {action varargin{:}};
    action   = 'Disp';
end

H.clip = [];
H.clim = [];
H.bkg_col = [1 1 1];
H.transp = 0;
H.data_sel = 0;

%-Action
%--------------------------------------------------------------------------
switch lower(action)
    
    %-Display
    %======================================================================
    case 'disp'

        % positions & font size
        ws = spm('Winsize','Graphics');
        FS = spm('FontSizes');
        
        % different positions for views with 4 and 5 images
        H.viewpos = {[0.075 0.450 0.325 0.325;  0.150 0.450 0.325 0.325],...
                     [0.075 0.050 0.325 0.325;  0.150 0.050 0.325 0.325],...
                     [0.600 0.450 0.325 0.325;  0.525 0.450 0.325 0.325],...
                     [0.600 0.050 0.325 0.325;  0.525 0.050 0.325 0.325],...
                     [0.300 0.200 0.400 0.400;  0.300 2.000 0.400 0.400]};

        % figure 1
        H.pos{1} = struct(...
            'fig',   [10  10  2*ws(3) ws(3)],...   % figure
            'cbar',  [0.400 0.550 0.200 0.300; 0.440 0.700 0.120 0.120]);   % colorbar   

        % figure 2
        H.pos{2} = struct(...
            'fig',   [2*ws(3)+10 10 0.6*ws(3) ws(3)],...   % figure
            'left',  [0.050 0.925 0.425 0.050],... % select left hemisphere
            'right', [0.525 0.925 0.425 0.050],... % select right hemisphere
            'surf',  [0.050 0.855 0.425 0.050],... % 
            'atlas', [0.525 0.855 0.425 0.050],... % 
            'thresh',[0.525 0.800 0.425 0.050],... % 
            'cmap',  [0.050 0.800 0.425 0.050],... % 
            'tview', [0.050 0.750 0.425 0.050],... % 
            'nocbar',[0.050 0.700 0.425 0.050],... % 
            'bkg',   [0.525 0.750 0.425 0.050],... % 
            'transp',[0.525 0.700 0.425 0.050],... % 
            'inv',   [0.525 0.650 0.425 0.050],... % 
            'info',  [0.050 0.650 0.425 0.050],... % 
            'ovmin', [0.050 0.450 0.425 0.150],... % 
            'ovmax', [0.525 0.450 0.425 0.150],... % 
            'save',  [0.050 0.050 0.425 0.050],... % 
            'close', [0.525 0.050 0.425 0.050],... % close button
            'text',  [0.050 0.150 0.425 0.200]);   % textbox   

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
                'ToolTipString','Select results for left hemisphere (log-p maps)',...
                'Interruptible','on','Enable','on');
        
        H.right = uicontrol(H.figure(2),...
                'string','Select right hemisphere data','Units','normalized',...
                'position',H.pos{2}.right,...
                'style','Pushbutton','HorizontalAlignment','center',...
                'callback',{@select_data,2},...
                'ToolTipString','Select results for right hemisphere (log-p maps)',...
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

        str  = { 'Threshold...','No threshold','P<0.05','P<0.01','P<0.001'};
        tmp  = { {@select_thresh, 0},...
                 {@select_thresh, 1.3},...
                 {@select_thresh, 2},...
                 {@select_thresh, 3}};
        
        H.thresh = uicontrol(H.figure(2),...
                'string',str,'Units','normalized',...
                'position',H.pos{2}.thresh,'UserData',tmp,...
                'style','PopUp','HorizontalAlignment','center',...
                'callback','spm(''PopUpCB'',gcbo)',...
                'ToolTipString','Threshold',...
                'Interruptible','on','Visible','off');
        
        str  = { 'Colormap...','jet','hot','hsv','cold-hot'};
        tmp  = { {@select_cmap, 1},...
                 {@select_cmap, 2},...
                 {@select_cmap, 3},...
                 {@select_cmap, 4}};
        
        H.cmap = uicontrol(H.figure(2),...
                'string',str,'Units','normalized',...
                'position',H.pos{2}.cmap,'UserData',tmp,...
                'style','PopUp','HorizontalAlignment','center',...
                'callback','spm(''PopUpCB'',gcbo)',...
                'ToolTipString','Threshold',...
                'Interruptible','on','Visible','off');

        str  = { 'Atlas labeling...','Desikan-Killiany DKT40','Destrieux 2009'};
        tmp  = { {@select_atlas, 1},...
                 {@select_atlas, 2}};
        
        H.atlas = uicontrol(H.figure(2),...
                'string',str,'Units','normalized',...
                'position',H.pos{2}.atlas,'UserData',tmp,...
                'style','PopUp','HorizontalAlignment','center',...
                'callback','spm(''PopUpCB'',gcbo)',...
                'ToolTipString','Atlas Labeling',...
                'Interruptible','on','Visible','off');

        H.tview = uicontrol(H.figure(2),...
                'string','Hide top view','Units','normalized',...
                'position',H.pos{2}.tview,...
                'style','CheckBox','HorizontalAlignment','center',...
                'callback',{@checkbox_tview},...
                'ToolTipString','Hide top view in the image center',...
                'Interruptible','on','Visible','off');

        H.inv = uicontrol(H.figure(2),...
                'string','Invert results','Units','normalized',...
                'position',H.pos{2}.inv,...
                'style','CheckBox','HorizontalAlignment','center',...
                'callback',{@checkbox_inv},...
                'ToolTipString','Invert results',...
                'Interruptible','on','Visible','off');

        H.bkg = uicontrol(H.figure(2),...
                'string','Black background','Units','normalized',...
                'position',H.pos{2}.bkg,...
                'style','CheckBox','HorizontalAlignment','center',...
                'callback',{@checkbox_bkg},...
                'ToolTipString','Black background',...
                'Interruptible','on','Visible','off');

        H.transp = uicontrol(H.figure(2),...
                'string','Transparent overlay','Units','normalized',...
                'position',H.pos{2}.transp,...
                'style','CheckBox','HorizontalAlignment','center',...
                'callback',{@checkbox_transp},...
                'ToolTipString','Transparent overlay',...
                'Interruptible','on','Visible','off');

        H.info = uicontrol(H.figure(2),...
                'string','Show filename','Units','normalized',...
                'position',H.pos{2}.info,...
                'style','CheckBox','HorizontalAlignment','center',...
                'callback',{@checkbox_info},...
                'ToolTipString','Show file information in image',...
                'Interruptible','on','Visible','off');

        H.nocbar = uicontrol(H.figure(2),...
                'string','Hide colorbar','Units','normalized',...
                'position',H.pos{2}.nocbar,...
                'style','CheckBox','HorizontalAlignment','center',...
                'callback',{@checkbox_nocbar},...
                'ToolTipString','Hide colorbar',...
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
          
          [pth{1},nm1,ext1] = spm_fileparts(varargin{1});
          [pth{2},nm2,ext2] = spm_fileparts(varargin{2});
          
          % SPM.mat found for both hemispheres
          if strcmp([nm1 ext1],'SPM.mat') && strcmp([nm2 ext2],'SPM.mat')
            H.logP = 0;
            
            for ind=1:2
              swd1 = pwd;
              spm_figure('GetWin','Interactive');
              cd(pth{ind})
              xSPM.swd = pwd;
              [xSPM,v] = spm_getSPM(xSPM);
              cd(swd1);
              
              dat = struct('XYZ', v.XYZ,...
                        't',   v.Z',...
                        'mat', v.M,...
                        'dim', v.DIM,...
                        'dat', v.Z');
              
              H.S{ind}.info = cat_surf_info(H.S{ind}.name,0); 
              g = gifti(H.S{ind}.info.Pmesh);

        mat    = v.M;
        V = g.vertices;
    XYZ        = double(inv(mat)*[V';ones(1,size(V,1))]);
%
    H.S{ind}.Y     = spm_sample_vol(Y,XYZ(1,:),XYZ(2,:),XYZ(3,:),0)';

              H.S{ind}.Y = spm_mesh_project(g.vertices,dat)';
            end
          else
          
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
          end
          
          H.disable_tview = 0;
          H.show_neg = 1;
          H.disable_cbar = 0;
          H.show_transp = 0;
          H.black_bkg = 0;
          H.show_info = 0;
          
          display_results_all;
          
          set(H.surf,'Enable','on');
          set(H.save,'Enable','on');
          set(H.tview,'Visible','on');
          set(H.nocbar,'Visible','on');
          set(H.bkg,'Visible','on');
          set(H.transp,'Visible','on');
          set(H.info,'Visible','on');
        
          if min(min(H.S{1}.Y(:),H.S{2}.Y(:))) < 0
            set(H.inv,'Visible','on');
          end
          
          if (size(H.S{1}.name,1) == 1) && (size(H.S{2}.name,1) == 1)
            set(H.cmap,'Visible','on');
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
        H.clip = getappdata(H.patch(1), 'clip');
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
        if numel(H.S{1}.info) > 1
            set(ic,'CData',c(1:numel(H.S{1}.info),:,:));
            set(ic,'YData',[1 numel(H.S{1}.info)]);
            set(H.colourbar,'YLim',[1 numel(H.S{1}.info)]);
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

H.thresh_value = thresh;

if H.show_neg
  H.clip = [true -thresh thresh];
else
  H.clip = [true 0 thresh];
end

% rather use NaN values for zero threshold
if thresh == 0
  H.clip = [false NaN NaN];
end

for ind=1:5
  setappdata(H.patch(ind),'clip',H.clip);
  col = getappdata(H.patch(ind),'col');
  d = getappdata(H.patch(ind),'data');
  H = updateTexture(H,ind,d,col,H.show_transp);
end

set(H.atlas,'Visible','on');

if ~H.disable_cbar
  H = show_colorbar(H);
end

%-----------------------------------------------------------------------
function H = select_cmap(cmap)
%-----------------------------------------------------------------------
global H

switch cmap
  case 1
    col = jet(256);
  case 2
    col = hot(256);
  case 3
    col = hsv(256);
  case 4
    col = [1-hot(128);(hot(128))];
end

for ind=1:5
  setappdata(H.patch(ind),'col',col);
  d = getappdata(H.patch(ind),'data');
  H = updateTexture(H,ind,d,col,H.show_transp);
end

if ~H.disable_cbar
  H = show_colorbar(H);
end

%-----------------------------------------------------------------------
function H = select_atlas(atlas)
%-----------------------------------------------------------------------
global H

% get threshold from clipping
thresh = [0 0];
if ~isempty(H.clip)
  if ~isnan(H.clip(2)) && ~isnan(H.clip(3))
    thresh = [H.clip(2:3)];
  end
end

for ind = [1 3]
  if atlas == 1 % DKT40 atlas
    atlas_name = fullfile(spm('dir'),'toolbox','cat12','atlases_surfaces',[H.S{round(ind/2)}.info(1).side ....
      '.aparc_DKT40JT.freesurfer.annot']);
  else % Destrieux
    atlas_name = fullfile(spm('dir'),'toolbox','cat12','atlases_surfaces',[H.S{round(ind/2)}.info(1).side ...
      '.aparc_a2009s.freesurfer.annot']);
  end
  [vertices, rdata, colortable, rcsv] = cat_io_FreeSurfer('read_annotation',atlas_name);

  M = getappdata(H.patch(ind),'patch');
  A       = spm_mesh_adjacency(M.faces);
  A       = A + speye(size(A));
  d = getappdata(H.patch(ind),'data');

  % apply thresholds
  dp = d > thresh(2); indp = find(dp);
  dn = d < thresh(1); indn = find(dn);
  
  % atlas name
  if atlas == 1
    atlas_name = 'Desikan-Killiany DKT40 Atlas';
  elseif atlas == 2
    atlas_name = 'Destrieux 2009 Atlas';
  end

  % go through pos. effects
  if ~isempty(indp)
  
    C = find_connected_component(A, dp);
    C = C(indp);
    rdata2 = rdata(indp);
  
    fprintf('\n\n______________________________________________________\n');
    fprintf('%s: Positive effects in %s',atlas_name,H.S{round(ind/2)}.info(1).side);
    fprintf('\n______________________________________________________\n\n');
  
    if H.logP, fprintf('%7s\t%8s\t%s\n','P-value','Size','Overlap of atlas region');
    else,      fprintf('%7s\t%8s\t%s\n','Value  ','Size','Overlap of atlas region'); end

    for i = 1:max(C)
      N = find(C == i);
      k = length(N);
    
      dmax = d(indp); dmax = max(dmax(N));
      
      if H.logP, fprintf('\n%1.5f\t%8d',10^(-dmax),k);
      else,      fprintf('\n%6.1f\t%8d',dmax,k); end
      
      Nrdata = rdata2(N);
      roi_size = zeros(size(rcsv,1)-1,1);
      
      for j=2:size(rcsv,1)
        ind3 = find(Nrdata == rcsv{j,1});
        roi_size(j-1) = 100*length(ind3)/k;
      end

      % sort wrt size
      [ii, jj] = sort(roi_size,'descend');
      jj(ii==0) = [];
      
      for j=1:length(jj)
        if roi_size(jj(j)) > 1
          if j==1, fprintf('\t%3.1f%s\t%s\n',roi_size(jj(j)),'%',rcsv{jj(j)+1,2});
          else,    fprintf('%7s\t%8s\t%3.1f%s\t%s\n','       ','        ',...
                roi_size(jj(j)),'%',rcsv{jj(j)+1,2}); 
          end
        end
      end

    end
  end
      
  % go through neg. effects
  if ~isempty(indn)

    C = find_connected_component(A, dn);
    C = C(indn);
    rdata2 = rdata(indn);

    fprintf('\n\n______________________________________________________\n');
    fprintf('%s: Negative effects in %s',atlas_name,H.S{round(ind/2)}.info(1).side);
    fprintf('\n______________________________________________________\n\n');
  
    if H.logP, fprintf('%7s\t%8s\t%s\n','P-value','Size','Overlap of atlas region');
    else,      fprintf('%7s\t%8s\t%s\n','Value  ','Size','Overlap of atlas region'); end

    for i = 1:max(C)
      N = find(C == i);
      k = length(N);
    
      dmin = d(indn); dmin = min(dmin(N));
      if H.logP, fprintf('\n%1.5f\t%8d',10^(dmin),k);
      else,      fprintf('\n%6.1f\t%8d',-dmin,k); end

      Nrdata = rdata2(N);
      roi_size = zeros(size(rcsv,1)-1,1);
      for j=2:size(rcsv,1)
        ind3 = find(Nrdata == rcsv{j,1});
        roi_size(j-1) = 100*length(ind3)/k;
      end

      % sort wrt size
      [ii, jj] = sort(roi_size,'descend');
      jj(ii==0) = [];
      
      for j=1:length(jj)
        if roi_size(jj(j)) > 1
          if j==1, fprintf('\t%3.1f%s\t%s\n',roi_size(jj(j)),'%',rcsv{jj(j)+1,2});
          else,    fprintf('%7s\t%8s\t%3.1f%s\t%s\n','       ','        ',...
                roi_size(jj(j)),'%',rcsv{jj(j)+1,2}); 
          end
        end
      end
      
    end
  end
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

g{1} = gifti(H.S{1}.info(1).Pmesh);
g{2} = gifti(H.S{2}.info(1).Pmesh);

for ind = 1:5
  if ind < 5
    M  = g{round(ind/2)};
  else
    M.faces = [g{1}.faces; g{2}.faces + size(g{1}.vertices,1)];
    M.vertices = [g{1}.vertices; g{2}.vertices];
    M.mat = g{1}.mat;
  end

  set(H.patch(ind),'Vertices',M.vertices);
  set(H.patch(ind),'Faces',M.faces);
end

%-----------------------------------------------------------------------
function display_results_all(obj, event_obj)
%-----------------------------------------------------------------------
global H

if (size(H.S{1}.Y) > 1 | size(H.S{2}.Y) > 1) & min(min(H.S{1}.Y(:),H.S{2}.Y(:))) < 0
  disp('Warning: Only results with positive values are displayed!');
end

% clear larger area and set background color to update labels and title
H.Ha = axes('Parent',H.figure(1),'Position',[-.1 -.1 1.1 1.1],'Color',H.bkg_col);
cla(H.Ha);

H.renderer = get(H.figure(1),'Renderer');
set(H.figure(1),'Renderer','OpenGL');

%-Compute mesh curvature
%------------------------------------------------------------------
g = gifti(fullfile(spm('dir'),'toolbox','cat12','templates_surfaces',[H.S{1}.info(1).side '.mc.central.freesurfer.gii']));
H.S{1}.curv = g.cdata;
g = gifti(fullfile(spm('dir'),'toolbox','cat12','templates_surfaces',[H.S{2}.info(1).side '.mc.central.freesurfer.gii']));
H.S{2}.curv = g.cdata;

vv = [90 0; -90 0; -90 0; 90 0; 0 90];
for ind = 1:5
  display_results(ind, H.viewpos{ind}(H.disable_tview+1,:), vv(ind,:));
end

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

H.clim = [true H.S{1}.min H.S{1}.max];
for ind=1:5
  setappdata(H.patch(ind), 'clim', [true H.S{1}.min H.S{1}.max]);
  col = getappdata(H.patch(ind), 'col');
  d = getappdata(H.patch(ind),'data');
  H = updateTexture(H,ind,d,col,H.show_transp);
end

% only show threshold popup if log-name was found and minimal value > 0 is lt 1.3
if H.logP && (H.S{1}.thresh < 1.3)
  set(H.thresh,'Visible','on');
end

if numel(H.S{1}.info)==1
  % get sure that image is thresholded and there are at least 20% zero/NaN areas
  if (sum(d~=0)/numel(d) < 0.8)         
    set(H.atlas,'Visible','on');
  end
end

if ~H.disable_cbar
  H = show_colorbar(H);
end

% show slider for range of results
if numel(H.S{1}.info)==1

  % allow slider a more extended range
  mnx = 1.5*max(abs([H.S{1}.min H.S{1}.max]));

  sliderPanel(...
        'Parent'  , H.figure(2), ...
        'Title'   , 'Overlay min', ...
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
        'Title'   , 'Overlay max', ...
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
  if ~isfield(H,'cbar') || ~ishandle(H.cbar)
    H.cbar = axes('Parent',H.figure(1),'Position',H.pos{1}.cbar(1,:),'Color',[0.5 0.5 0.5],'Visible','off');
    H.colourbar = colorbar('peer',H.cbar,'Northoutside');
  end
  if H.logP, title('p-value','Color',1-H.bkg_col);end
  clim = getappdata(H.patch(1), 'clim');
  axis(H.cbar,'off'); caxis([clim(2) clim(3)]);
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
  set(H.colourbar,'XColor',1-H.bkg_col,'YColor',1-H.bkg_col);
  
  % Update colorbar colors if clipping is used
  clip = getappdata(H.patch(1), 'clip');
  col = getappdata(H.patch(1), 'col');
  if ~isempty(clip)
    if ~isnan(clip(2)) && ~isnan(clip(3))
      ncol = length(col);
      col_step = (clim(3) - clim(2))/ncol;
      cmin = max([1,ceil((clip(2)-clim(2))/col_step)]);
      cmax = min([ncol,floor((clip(3)-clim(2))/col_step)]);
      col(cmin:cmax,:) = repmat([0.5 0.5 0.5],(cmax-cmin+1),1);
      colormap(col);
    end
  end

else

  if ~isfield(H,'cbar') || ~ishandle(H.cbar)
    H.cbar = axes('Parent',H.figure(1),'Position',H.pos{1}.cbar(2,:),'Color',[0.5 0.5 0.5],'Visible','off');
  end
  
  % RGB colorbar
  if numel(H.S{1}.info) ==3
    cb = [8 1 1 4 2 2 8;...
          8 1 6 7 5 2 8;...
          8 8 3 3 3 8 8];
  else %RG colorbar
    cb = [8 1 1 4 2 2 8;...
          8 1 1 4 2 2 8];
  end
  imagesc(cb);
  colormap([1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 1 1 1; H.bkg_col])
  axis(H.cbar,'off'); axis('image');  
end

%-----------------------------------------------------------------------
function display_results(ind, win, vw)
%-----------------------------------------------------------------------
global H

% rescue old color before a new H.patch is created
try
  col = getappdata(H.patch(ind), 'col');
catch
  col = [];
end

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
setappdata(H.patch(ind),'axis',H.axis);

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

if isempty(col)
  H = updateTexture(H,ind,T);
else
  H = updateTexture(H,ind,T,col);
end

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
function [H, C] = updateTexture(H,ind,v,col,transp)

%-Project data onto surface mesh
%--------------------------------------------------------------------------
if size(v,2) < size(v,1)
  v = v';
end
v(isinf(v)) = NaN;

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

setappdata(H.patch(ind),'data',v);
setappdata(H.patch(ind),'col',col);

if ~exist('FaceColor','var') || isempty(FaceColor), FaceColor = 'interp'; end

%-Get curvature
%--------------------------------------------------------------------------
curv = getappdata(H.patch(ind),'curvature');

if size(curv,2) == 1
    th = 0.15;
    curv((curv<-th)) = -2*th;
    curv((curv>th))  =  0.1*th;
    curv = 0.5*(curv + th)/(2*th);
    curv = 0.5 + repmat(curv,1,3);
    curv = curv/max(curv(:));
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
    if size(v,1) > 1 && mi < 0
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

if ~isempty(H.clip)
    v(v>H.clip(2) & v<H.clip(3)) = NaN;
    setappdata(H.patch(ind), 'clip', [true H.clip(2) H.clip(3)]);
end

setappdata(H.patch(ind), 'clim', [true mi ma]);
H.clim = [true mi ma];

%-Build texture by merging curvature and data
%--------------------------------------------------------------------------
if size(v,1) > 1 % RGB
  for i=1:size(v,1)
    C(:,i) = any(v(i,:),1)' .* C(:,i);
  end
else
  C = repmat(any(v,1),3,1)' .* C;
end

% add curvature pattern if transparency is defined
if nargin > 4
  if transp
    C = (0.5+0.5*curv) .* C;
  end
end

% replace regions below threshold by curvature
ind0 = repmat(~any(v,1),3,1)';
C(ind0) = curv(ind0);

set(H.patch(ind), 'FaceVertexCData',C, 'FaceColor',FaceColor);

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
  H.S{ind}.info = cat_surf_info(H.S{ind}.name,1); 
  try
    H.S{ind}.Y    = spm_data_read(spm_data_hdr_read(H.S{ind}.name));
  catch
    error('No cdata found.');
  end
  
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

% increase counter for selected data
H.data_sel = H.data_sel + 1;

H.disable_tview = 0;
H.disable_cbar = 0;
H.show_neg = 1;
H.clip = [false NaN NaN];
H.show_transp = 0;
H.black_bkg = 0;
H.show_info = 0;

% display if both sides are defined
if ~isempty(H.S{1}.name)  &&  ~isempty(H.S{2}.name) && H.data_sel == 2

  display_results_all;
  set(H.surf,'Enable','on');
  set(H.save,'Enable','on');
  set(H.tview,'Visible','on');
  set(H.nocbar,'Visible','on');
  set(H.bkg,'Visible','on');
  set(H.transp,'Visible','on');
  set(H.info,'Visible','on');

  if min(min(H.S{1}.Y(:),H.S{2}.Y(:))) < 0
    set(H.inv,'Visible','on');
  end

  if (size(H.S{1}.name,1) == 1) && (size(H.S{2}.name,1) == 1)
    set(H.cmap,'Visible','on');
  end

  % reset counter for selected data
  H.data_sel = 0;

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
 
    set(H.colourbar,'FontUnits','Normalized','FontSize',0.15);
    changed_font = 1;
  catch
    changed_font = 0;
  end
  
  % keep background color
  set(H.figure(1),'InvertHardcopy','off');
  
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
for ind = 1:5
  setappdata(H.patch(ind),'clim',[true val c(3)]);
  col = getappdata(H.patch(ind),'col');
  d = getappdata(H.patch(ind),'data');
  H = updateTexture(H,ind,d,col,H.show_transp);
end

% update colorbar 
if numel(H.S{1}.info) == 1 && ~H.disable_cbar
  H = show_colorbar(H);
end

H.clim = [true val c(3)];

%==========================================================================
function slider_clim_max(hObject, evt)
global H

figure(H.figure(1))
val = get(hObject, 'Value');
c = getappdata(H.patch(1),'clim');
for ind = 1:5
  setappdata(H.patch(ind),'clim',[true c(2) val]);
  col = getappdata(H.patch(ind),'col');
  d = getappdata(H.patch(ind),'data');
  H = updateTexture(H,ind,d,col,H.show_transp);
end

% update colorbar 
if numel(H.S{1}.info) == 1 && ~H.disable_cbar
  H = show_colorbar(H);
end

H.clim = [true c(2) val];

%==========================================================================
function checkbox_tview(obj, event_obj)
global H
  
H.disable_tview = get(H.tview,'Value');

for ind = 1:5
  Ha = getappdata(H.patch(ind),'axis');
  set(Ha,'position',H.viewpos{ind}(H.disable_tview+1,:));
end

%==========================================================================
function checkbox_inv(obj, event_obj)
global H
  
H.show_inv = get(H.inv,'Value');

for ind=1:5
  setappdata(H.patch(ind),'clip',H.clip);
  col = getappdata(H.patch(ind),'col');
  d = getappdata(H.patch(ind),'data');
  H = updateTexture(H,ind,-d,col,H.show_transp);
end

%==========================================================================
function checkbox_transp(obj, event_obj)
global H
  
H.show_transp = get(H.transp,'Value');

for ind=1:5
  col = getappdata(H.patch(ind),'col');
  d = getappdata(H.patch(ind),'data');
  H = updateTexture(H,ind,d,col,H.show_transp);
end

% update colorbar 
if numel(H.S{1}.info) == 1 && ~H.disable_cbar
  H = show_colorbar(H);
end

%==========================================================================
function checkbox_bkg(obj, event_obj)
global H
  
H.black_bkg = get(H.bkg,'Value');

if H.black_bkg
  H.bkg_col = [0 0 0];
else
  H.bkg_col = [1 1 1];
end

set(H.Ha,'Color',H.bkg_col);
set(get(H.cbar,'Title'),'Color',1-H.bkg_col);

if H.show_info
  set(get(getappdata(H.patch(1),'axis'),'Title'),'Color',1-H.bkg_col);
  set(get(getappdata(H.patch(3),'axis'),'Title'),'Color',1-H.bkg_col);
end

if numel(H.S{1}.info) == 1
  set(H.colourbar,'XColor',1-H.bkg_col,'YColor',1-H.bkg_col);
end

if ~H.disable_cbar
  H = show_colorbar(H);
end

%==========================================================================
function checkbox_info(obj, event_obj)
global H
  
H.show_info = get(H.info,'Value');

if H.show_info
  set(get(getappdata(H.patch(1),'axis'),'Title'),'String',...
      spm_str_manip(H.S{1}.name,'k50d'),'Interpreter', 'none','Color',1-H.bkg_col)
  set(get(getappdata(H.patch(3),'axis'),'Title'),'String',...
      spm_str_manip(H.S{2}.name,'k50d'),'Interpreter', 'none','Color',1-H.bkg_col)
else
  set(get(getappdata(H.patch(1),'axis'),'Title'),'String','')
  set(get(getappdata(H.patch(3),'axis'),'Title'),'String','')
end

%==========================================================================
function checkbox_noneg(obj, event_obj)
global H
  
H.show_neg = 1 - get(H.noneg,'Value');

clim = getappdata(H.patch(1), 'clim');

if H.show_neg
  if isfield(H,'thresh_value')
  tmp = H.thresh_value
    H.clip = [true -H.thresh_value H.thresh_value];
  else
    H.clip = [true -Inf -Inf];
  end
else
  if ~isempty(H.clip)
    if ~isnan(H.clip(2)) && ~isnan(H.clip(3))
      H.clip(3) = 0;
    else
      H.clip(3) = 0;
    end
  else
    H.clip = [true -Inf 0];
  end
end

for ind=1:5
  setappdata(H.patch(ind),'clim',clim);
  setappdata(H.patch(ind),'clip',H.clip);
  col = getappdata(H.patch(ind),'col');
  d = getappdata(H.patch(ind),'data');
  H = updateTexture(H,ind,d,col,H.show_transp);
end

if ~H.disable_cbar
  H = show_colorbar(H);
end

%==========================================================================
function checkbox_nocbar(obj, event_obj)
global H
  
H.disable_cbar = get(H.nocbar,'Value');

if H.disable_cbar
  % delete colorbar and title
  if numel(H.S{1}.info) == 1
    set(H.colourbar,'Visible','off')  
    set(get(H.cbar,'Title'),'Visible','off')
  else % delete only axis
    cla(H.cbar);
  end
else
  if numel(H.S{1}.info) == 1
    set(get(H.cbar,'Title'),'Visible','on')
    set(H.colourbar,'Visible','on')  
  else
    H = show_colorbar(H);
  end
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

%==========================================================================
function C = find_connected_component(A, T);
% find connected components 
% FORMAT C = find_connected_component(A,T)
% A        - a [nxn[ (reduced) adjacency matrix
% T        - a [nx1] data vector (using NaNs or logicals), n = #vertices
%
% C        - a [nx1] vector of cluster indices
%
% modified version from spm_mesh_clusters.m 5065 2012-11-16 20:00:21Z guillaume
%


%-Input parameters
%--------------------------------------------------------------------------
if ~islogical(T)
  T   = ~isnan(T);
end
  
A1 = A;
A1(~T,:) = [];
A1(:,~T) = [];

%-And perform Dulmage-Mendelsohn decomposition to find connected components
%--------------------------------------------------------------------------
[p,q,r] = dmperm(A1);
N       = diff(r);
CC      = zeros(size(A1,1),1);
for i = 1:length(r)-1
  CC(p(r(i):r(i+1)-1)) = i;
end
C       = NaN(numel(T),1);
C(T)    = CC;

%-Sort connected component labels according to their size
%--------------------------------------------------------------------------
[N,ni]  = sort(N(:), 1, 'descend');
[ni,ni] = sort(ni);
C(T)    = ni(C(T));

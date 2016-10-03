function y = cat_surf_results(action,varargin)
% Visualise results for both hemispheres of surface-based analysis (preferable on log P-maps)
%
% FORMAT y = cat_surf_results('Disp',leftSurface,rightSurface,[show_warning])
% leftSurface  - a GIfTI filename/object or patch structure
% rightSurface - a GIfTI filename/object or patch structure
% show_warning - (default 1) show warning if filenames differ (it is very likely that even if they are 
%                in different folders for each hemisphere and have the same name)
%
% y            - adjusted or predicted response
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

% set start values
H.clip     = [];
H.clim     = [];
H.bkg_col  = [0 0 0];
H.transp   = 1;
H.data_sel = [0 0];
H.data_n   = [1 3];
H.cursor_mode = 1;
H.XTick    = [];
y          = [];

%-Action
%--------------------------------------------------------------------------
switch lower(action)
    
    %-Display
    %======================================================================
    case 'disp'

        % positions & font size
        ws = spm('Winsize','Graphics');
        FS = spm('FontSizes');
        
        % result window with 5 surface views and alternative positions without top view
        H.viewpos = {[0.075 0.450 0.325 0.325;  0.150 0.550 0.325 0.325],... % lh medial
                     [0.075 0.050 0.325 0.325;  0.150 0.150 0.325 0.325],... % lh lateral
                     [0.600 0.450 0.325 0.325;  0.525 0.550 0.325 0.325],... % rh medial
                     [0.600 0.050 0.325 0.325;  0.525 0.150 0.325 0.325],... % rh lateral
                     [0.300 0.180 0.400 0.400;  0.300 2.000 0.400 0.400]};   % lh+rh top
                     
        % figure 1 with result window
        H.pos{1} = struct(...
            'fig',   [10  10  2*ws(3) ws(3)],...   % figure
            'cbar',  [0.400 -0.180 0.200 0.300; 0.440 0.050 0.120 0.120]);   % colorbar   

        % figure 2 with GUI
        H.pos{2} = struct(...
          'fig',   [2*ws(3)+10 10 0.6*ws(3) ws(3)],... 
          'sel',   [0.290 0.930 0.425 0.060],...
          'surf',  [0.050 0.855 0.425 0.050],'atlas', [0.525 0.855 0.425 0.050],... 
          'cursor',[0.050 0.800 0.425 0.050],'thresh',[0.525 0.800 0.425 0.050],... 
          'cmap',  [0.050 0.750 0.425 0.050],...
          'tview', [0.050 0.700 0.425 0.050],'bkg',   [0.525 0.700 0.425 0.050],... 
          'nocbar',[0.050 0.650 0.425 0.050],'transp',[0.525 0.650 0.425 0.050],... 
          'info',  [0.050 0.600 0.425 0.050],'inv',   [0.525 0.600 0.425 0.050],... 
          'ovmin', [0.050 0.400 0.425 0.150],'ovmax', [0.525 0.400 0.425 0.150],... 
          'save',  [0.050 0.050 0.425 0.050],'close', [0.525 0.050 0.425 0.050]);   

        % figure 3 with data plot
        H.pos{3} = struct(...
            'fig',   [10  30+ws(3)  ws(3) 0.5*ws(3)]);

        % create figures 1+2
        for i=1:2
          H.figure(i) = figure(i+11);
          clf(H.figure(i));
        
          set(H.figure(i),'MenuBar','none','Position',H.pos{i}.fig,...
            'Name','Results','NumberTitle','off');
        end
        
        % define S structure that contains information for lh and rh
        H.S{1}.name = ''; H.S{1}.side = 'lh';
        H.S{2}.name = ''; H.S{2}.side = 'rh';
        
        % closing all windows
        H.close = uicontrol(H.figure(2),...
                'string','Close','Units','normalized',...
                'position',H.pos{2}.close,...
                'style','Pushbutton','HorizontalAlignment','center',...
                'callback','for i=2:26, try close(i); end; end;',...
                'ToolTipString','Close windows',...
                'Interruptible','on','Enable','on');

        % select results for lh and rh
        H.sel = uicontrol(H.figure(2),...
                'string','Select surface data','Units','normalized',...
                'position',H.pos{2}.sel,...
                'style','Pushbutton','HorizontalAlignment','center',...
                'callback',@select_data,...
                'ToolTipString','Select results (up to 3) for both hemispheres (e.g. log-p maps)',...
                'Interruptible','on','Enable','on');
                
        str  = { 'Underlying surface...','central','inflated','Dartel'};
        tmp  = { {@select_surf, 1},...
                 {@select_surf, 2},...
                 {@select_surf, 3}};
        
        % underlying surface
        H.surf = uicontrol(H.figure(2),...
                'string',str,'Units','normalized',...
                'position',H.pos{2}.surf,'UserData',tmp,...
                'style','PopUp','HorizontalAlignment','center',...
                'callback','spm(''PopUpCB'',gcbo)',...
                'ToolTipString','Underlying surface',...
                'Interruptible','on','Visible','off');

        str  = { 'Threshold...','No threshold','P<0.05','P<0.01','P<0.001'};
        tmp  = { {@select_thresh, 0},...
                 {@select_thresh, 1.3},...
                 {@select_thresh, 2},...
                 {@select_thresh, 3}};
        
        % threshold
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
        
        % colormap
        H.cmap = uicontrol(H.figure(2),...
                'string',str,'Units','normalized',...
                'position',H.pos{2}.cmap,'UserData',tmp,...
                'style','PopUp','HorizontalAlignment','center',...
                'callback','spm(''PopUpCB'',gcbo)',...
                'ToolTipString','Threshold',...
                'Interruptible','on','Visible','off');

        str  = { 'Atlas labeling...','Desikan-Killiany DK40','Destrieux 2009'};
        tmp  = { {@select_atlas, 1},...
                 {@select_atlas, 2}};
        
        % atla for labeling
        H.atlas = uicontrol(H.figure(2),...
                'string',str,'Units','normalized',...
                'position',H.pos{2}.atlas,'UserData',tmp,...
                'style','PopUp','HorizontalAlignment','center',...
                'callback','spm(''PopUpCB'',gcbo)',...
                'ToolTipString','Atlas Labeling',...
                'Interruptible','on','Visible','off');

        str  = { 'Plot Data...','Nothing','Desikan-Killiany DK40',...
                 'Destrieux 2009','Plot data at vertex','Plot mean data inside cluster'};
        tmp  = { {@select_cursor, 0},...
                 {@select_cursor, 1},...
                 {@select_cursor, 2},...
                 {@select_cursor, 3},...
                 {@select_cursor, 4}};
        
        % data cursor for data plotting and atlas names
        H.cursor = uicontrol(H.figure(2),...
                'string',str,'Units','normalized',...
                'position',H.pos{2}.cursor,'UserData',tmp,...
                'style','PopUp','HorizontalAlignment','center',...
                'callback','spm(''PopUpCB'',gcbo)',...
                'ToolTipString','Data Cursor Mode',...
                'Interruptible','on','Visible','off');

        % disable top view
        H.tview = uicontrol(H.figure(2),...
                'string','Hide top view','Units','normalized',...
                'position',H.pos{2}.tview,...
                'style','CheckBox','HorizontalAlignment','center',...
                'callback',{@checkbox_tview},...
                'ToolTipString','Hide top view in the image center',...
                'Interruptible','on','Visible','off');

        % invert results
        H.inv = uicontrol(H.figure(2),...
                'string','Invert results','Units','normalized',...
                'position',H.pos{2}.inv,...
                'style','CheckBox','HorizontalAlignment','center',...
                'callback',{@checkbox_inv},...
                'ToolTipString','Invert results',...
                'Interruptible','on','Visible','off');

        % white background
        H.bkg = uicontrol(H.figure(2),...
                'string','White background','Units','normalized',...
                'position',H.pos{2}.bkg,...
                'style','CheckBox','HorizontalAlignment','center',...
                'callback',{@checkbox_bkg},...
                'ToolTipString','White background',...
                'Interruptible','on','Visible','off');

        % transparent view
        H.transp = uicontrol(H.figure(2),...
                'string','Disable transparency','Units','normalized',...
                'position',H.pos{2}.transp,...
                'style','CheckBox','HorizontalAlignment','center',...
                'callback',{@checkbox_transp},...
                'ToolTipString','Disable transparent overlay',...
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
                'Interruptible','on','Visible','off');

        if nargin >= 3
        
          H.S{1}.name = varargin{1};
          H.S{2}.name = varargin{2};
          
          [pth{1},nm1,ext1] = spm_fileparts(H.S{1}.name(1,:));
          [pth{2},nm2,ext2] = spm_fileparts(H.S{2}.name(1,:));
          
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
                             't', v.Z',...
                           'mat', v.M,...
                           'dim', v.DIM,...
                           'dat', v.Z');
              
              H.S{ind}.info = cat_surf_info(H.S{ind}.name,0); 
              g = gifti(H.S{ind}.info.Pmesh);

              mat    = v.M;
              V = g.vertices;
              XYZ        = double(inv(mat)*[V';ones(1,size(V,1))]);
              H.S{ind}.Y = spm_sample_vol(Y,XYZ(1,:),XYZ(2,:),XYZ(3,:),0)';
              H.S{ind}.Y = spm_mesh_project(g.vertices,dat)';
            end
          else
          
            H.logP = 1;
          
            for ind=1:2

              H.S{ind}.info = cat_surf_info(H.S{ind}.name,1); 

              % cdata found?
              try
                H.S{ind}.Y = spm_data_read(spm_data_hdr_read(H.S{ind}.name));
              catch
                if ind == 1
                  try
                    H.S{ind}.Y = zeros(size(H.S{2}.Y));
                  catch
                    error('No data in surfaces found.');
                  end
                else
                  H.S{ind}.Y = zeros(size(H.S{1}.Y));
                end
                
                % reset name
                H.S{ind}.name = '';
              end
                
              % check whether name contains 'log' that indicates a logP file
              for i=1:size(H.S{ind}.name,1)
                if isempty(strfind(H.S{ind}.info(i).ff,'log'))
                  H.logP = 0;
                end
              end
            end
          end
          
          H.disable_tview = 0;
          H.show_neg      = 1;
          H.disable_cbar  = 0;
          H.show_transp   = 1;
          H.white_bgk     = 0;
          H.show_info     = 0;
          
          display_results_all;
          
          set(H.surf,'Visible','on');
          set(H.save,'Visible','on');
          set(H.tview,'Visible','on');
          set(H.nocbar,'Visible','on');
          set(H.bkg,'Visible','on');
          set(H.transp,'Visible','on');
          set(H.info,'Visible','on');
          set(H.cursor,'Visible','on');
        
          if min(min(H.S{1}.Y(:)),min(H.S{2}.Y(:))) < 0
            set(H.inv,'Visible','on');
          end
          
          if (size(H.S{1}.name,1) == 1) && (size(H.S{2}.name,1) == 1)
            set(H.cmap,'Visible','on');
          end
        
          H.rdata{1} = [];
          H.rdata{2} = [];
          for ind = 1:2
            atlas_name = fullfile(spm('dir'),'toolbox','cat12','atlases_surfaces',[H.S{ind}.info(1).side ....
              '.aparc_DK40.freesurfer.annot']);
            [vertices, rdata0, colortable, rcsv1] = cat_io_FreeSurfer('read_annotation',atlas_name);
            H.rdata{1} = [H.rdata{1} rdata0];
            atlas_name = fullfile(spm('dir'),'toolbox','cat12','atlases_surfaces',[H.S{ind}.info(1).side ....
              '.aparc_a2009s.freesurfer.annot']);
            [vertices, rdata0, colortable, rcsv2] = cat_io_FreeSurfer('read_annotation',atlas_name);
            H.rdata{2} = [H.rdata{2} rdata0];
          end
          H.rcsv{1} = rcsv1;
          H.rcsv{2} = rcsv2;

          H.dcm_obj = datacursormode(H.figure(1));
          set(H.dcm_obj, 'Enable','on', 'SnapToDataVertex','on', ...
            'DisplayStyle','datatip', 'Updatefcn',{@myDataCursorAtlas, H});
        
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

% atlas name
if atlas == 1
  atlas_name = 'Desikan-Killiany DK40 Atlas';
elseif atlas == 2
  atlas_name = 'Destrieux 2009 Atlas';
end

for ind = [1 3]

  % atlas data
  rcsv = H.rcsv{atlas};
  rdata = H.rdata{atlas}(:,round(ind/2));

  M = getappdata(H.patch(ind),'patch');
  A       = spm_mesh_adjacency(M.faces);
  A       = A + speye(size(A));
  d = getappdata(H.patch(ind),'data');

  % apply thresholds
  dp = d > thresh(2); indp = find(dp);
  dn = d < thresh(1); indn = find(dn);
  
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

if (size(H.S{1}.Y) > 1 | size(H.S{2}.Y) > 1) & min(min(H.S{1}.Y(:)),min(H.S{2}.Y(:))) < 0
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

H.S{1}.min = min(min(H.S{1}.Y(:)),min(H.S{2}.Y(:)));
H.S{1}.max = max(max(H.S{1}.Y(:)),max(H.S{2}.Y(:)));

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

% only show threshold popup if log-name was found and minimal value > 0 is lt 2
if H.logP & (H.S{1}.thresh < 2)
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
  mnx = 2*max(abs([H.S{1}.min H.S{1}.max]));

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
  if H.logP, title(H.cbar,'p-value','Color',1-H.bkg_col);end
  clim = getappdata(H.patch(1), 'clim');
  axis(H.cbar,'off'); caxis([clim(2) clim(3)]);
  colormap(getappdata(H.patch(1),'col'));
  
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

  if H.logP
    XTick = get(H.colourbar,'XTick');
    
    % save original XTick values
    if isempty(H.XTick), H.XTick = XTick; end

    % if threshold is 1.3 (p<0.05) change XTick accordingly and correct by 0.3
    if ~isempty(clip)
      if clip(3) == 1.3
        XTick_step = round((clim(3)-clim(2))/5);
        if clip(2) == -1.3
          XTick = [(round(clim(2))-0.3):XTick_step:-1.3 0 1.3:XTick_step:(round(clim(3))+0.3)];
        else
          XTick = [0 1.3:XTick_step:(round(clim(3))+0.3)];
        end
      end
    else
      % rescue original XThick values if clipping is changed
      if ~isempty(H.XTick), XTick = H.XTick; end
    end
    
    % change XTickLabel
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
    
  end % end H.logP
  
  set(H.colourbar,'XColor',1-H.bkg_col,'YColor',1-H.bkg_col);
  
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
  
  % check whether number of data for lh/rh differ and fill with zeros
  diff_size_Y = size(H.S{1}.Y,2) - size(H.S{2}.Y,2);
  if diff_size_Y > 0
    Mcr.cdata = [Mcr.cdata zeros(size(H.S{2}.Y,1),1)];
  end
  if diff_size_Y < 0
    Mcl.cdata = [Mcl.cdata; zeros(size(H.S{1}.Y,1),1)];
  end
  
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
function select_data(obj, event_obj)
%-----------------------------------------------------------------------
global H

H.logP = 1;
lh = []; rh = [];

P = spm_select([1 6],'mesh',['Select up to 3 (log P) maps for left and right hemisphere']);
info = cat_surf_info(P);

n = size(P,1);

for i=1:n

  % check whether name contains 'log' that indicates a logP file
  if isempty(strfind(info(i).ff,'log'))
    H.logP = 0;
  end

  % check where left and right hemisphere data were found
  if strcmp(info(i).side, 'lh')
    lh = [lh i];
  elseif strcmp(info(i).side, 'rh')
    rh = [rh i];
  end
end

H.S{1}.name = P(lh,:);
H.S{2}.name = P(rh,:);
H.S{1}.info = cat_surf_info(H.S{1}.name,1); 
H.S{2}.info = cat_surf_info(H.S{2}.name,1); 

if isempty(H.S{2}.name)
  H.S{2}.name = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces','rh.central.freesurfer.gii');
elseif isempty(H.S{1}.name)
  H.S{1}.name = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces','lh.central.freesurfer.gii');
end

cat_surf_results('disp',H.S{1}.name, H.S{2}.name);

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
  
H.show_transp = ~get(H.transp,'Value');

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
  
H.white_bgk = get(H.bkg,'Value');

if H.white_bgk
  H.bkg_col = [1 1 1];
else
  H.bkg_col = [0 0 0];
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
if ishandle(H) & ~isappdata(H,'handles')
    a = H; clear H;
    H.axis     = a;
    H.figure(1)   = ancestor(H.axis,'figure');
    H.patch    = findobj(H.axis,'type','patch');
    H.light    = findobj(H.axis,'type','light');
    H.rotate3d = rotate3d(H.figure(1));
    setappdata(H.axis,'handles',H);
else
    H = getappdata(H,'handles');
end

%==========================================================================
function select_cursor(mode)

global H y

H.cursor_mode = mode;

dcm_obj = datacursormode(H.figure(1));
switch H.cursor_mode

  case 0 % disable and delete datatip
    set(dcm_obj, 'Enable','off', 'SnapToDataVertex','on', ...
        'DisplayStyle','datatip', 'Updatefcn',{@myDataCursorAtlas,H});
    figure(H.figure(1))
    delete(findall(gca,'Type','hggroup','HandleVisibility','off'));
    try, close(H.figure(3)); end
    
  case {1,2}
    set(dcm_obj, 'Enable','on', 'SnapToDataVertex','on', ...
        'DisplayStyle','datatip', 'Updatefcn',{@myDataCursorAtlas,H});
    try, close(H.figure(3)); end
        
  case {3,4}
    figure(H.figure(1))
    try
      delete(findall(gca,'Type','hggroup','HandleVisibility','off'));
    end
    
    H.figure(3) = figure(3+11);
    set(H.figure(3),'MenuBar','none','Position',H.pos{3}.fig,...
        'Name','Plot','NumberTitle','off');
    
    SPM_found = 1;
    for i=1:2
      SPM_name = fullfile(H.S{i}.info(1).pp, 'SPM.mat');
    
      % SPM.mat exist?
      if exist(SPM_name,'file')
        load(SPM_name);
        H.SPM{i} = SPM;
        if i==1
          H.Ic = spm_input('Which contrast?',1,'m',{SPM.xCon.name});
          str   = 'predicted or adjusted values?';
          H.predicted = spm_input(str,2,'b',{'predicted','adjusted'},[1 0]);
        end
      else
        SPM_found = 0;
        spm('alert!','No SPM.mat file found.\nPlease check that you have not moved your files or your result file was moved from the folder where the SPM.mat is stored.',1);
      end
    end
    if SPM_found
      set(dcm_obj, 'Enable','on', 'SnapToDataVertex','on', ...
        'DisplayStyle','datatip', 'Updatefcn',{@myDataCursorCluster,H.cursor_mode-3});
    end


end

%==========================================================================
function txt = myDataCursorCluster(obj,evt,plot_mean)
global H y

pos = get(evt,'Position');

i = ismember(get(H.patch(1),'vertices'),pos,'rows');
node = find(i);
ind = 1;
node_list = 1:numel(get(H.patch(1),'vertices'));

if isempty(node)
  i = ismember(get(H.patch(3),'vertices'),pos,'rows');
  node = find(i);
  ind = 3;
  node_list = 1:numel(get(H.patch(3),'vertices'));
end

% get threshold from clipping
thresh = [0 0];
if ~isempty(H.clip)
  if ~isnan(H.clip(2)) && ~isnan(H.clip(3))
    thresh = [H.clip(2:3)];
  end
end

if plot_mean
  
  found_node = [];
  cluster_number = 0;
  cluster_side = 0;

if ~isfield(H,'A')
    M = getappdata(H.patch(ind),'patch');
    A       = spm_mesh_adjacency(M.faces);
    A       = A + speye(size(A));
    H.A = A;
else
    A = H.A;
    d = getappdata(H.patch(ind),'data');

    % apply thresholds
    dp = d > thresh(2); indp = find(dp);
    dn = d < thresh(1); indn = find(dn);
    
    % go through pos. effects
    if ~isempty(indp)
  
      C = find_connected_component(A, dp);
      C = C(indp);
      node_list2 = node_list(indp);
      
      for i = 1:max(C)
        N = find(C == i);          
        XYZ = node_list2(N);
        found_node = find(XYZ == node);
        if ~isempty(found_node)
          cluster_number = i;
          cluster_side = ind;
          break;
        end
      end
    end
  
    % go through neg. effects if no node was found
    if ~isempty(indn) && ~isempty(found_node)

      C = find_connected_component(A, dn);
      C = C(indn);
      node_list2 = node_list(indp);
  
      for i = 1:max(C)
        N = find(C == i);          
        XYZ = node_list2(N);
        found_node = find(XYZ == node);
        if ~isempty(found_node)
          cluster_number = i;
          cluster_side = ind;
          break;
        end
      end
    end
  end

  if isempty(found_node)
    txt = {'Cursor outside of cluster'};
  else
    if cluster_side == 1
      txt = {sprintf('lh: Cluster %d',cluster_number)};
    else
      txt = {sprintf('rh: Cluster %d',cluster_number)};
    end
  end
else
  % use single node as region 
  XYZ = node;
  txt = {sprintf('Node %d',node)};
end

y = get_cluster_data(H,XYZ,ind);
figure(H.figure(3))
if plot_mean
  y = mean(y,2);
  if isempty(found_node)
    y(:) = 0;
  end
end
plot(y)

%==========================================================================
function y = get_cluster_data(H,XYZ, ind)

SPM = H.SPM{round(ind/2)};
Ic = H.Ic;
predicted = H.predicted;

% get raw data and whiten
y = spm_data_read(SPM.xY.VY,'xyz',XYZ);
y = spm_filter(SPM.xX.K,SPM.xX.W*y);
R   = spm_sp('r',SPM.xX.xKXs,y);

beta   = spm_data_read(SPM.Vbeta,'xyz',XYZ);
ResMS = spm_data_read(SPM.VResMS,'xyz',XYZ);

% compute contrast of parameter estimates and 90% C.I.
%------------------------------------------------------------------
cbeta = SPM.xCon(Ic).c'*beta;

% predicted or adjusted response
%------------------------------------------------------------------
if predicted

  % fitted (predicted) data (Y = X1*beta)
  %--------------------------------------------------------------
  % this should be SPM.xX.xKXs.X instead of SPM.xX.X below
  Y = SPM.xX.X*SPM.xCon(Ic).c*pinv(SPM.xCon(Ic).c)*beta;
else

  % fitted (corrected)  data (Y = X1o*beta)
  %--------------------------------------------------------------
  Y = spm_FcUtil('Yc',SPM.xCon(Ic),SPM.xX.xKXs,beta);

end

% adjusted data
%------------------------------------------------------------------
y     = Y + R;

if 0
H.y{i} = H.y{i} - SPM.xX.xKXs.X * beta;
X0     = SPM.xX.xKXs.X(:,[SPM.xX.iB SPM.xX.iG]);
X0     = X0(:,any(X0));
        %-Compute regional response in terms of first eigenvariate
        %--------------------------------------------------------------------------
        [m,n]   = size(H.y{i});
        if m > n
          [v,s,v] = svd(H.y{i}'*H.y{i});
          s       = diag(s);
          v       = v(:,1);
          u       = y*v/sqrt(s(1));
        else
          [u,s,u] = svd(H.y{i}*H.y{i}');
          s       = diag(s);
          u       = u(:,1);
          v       = y'*u/sqrt(s(1));
        end
        d       = sign(sum(v));
        u       = u*d;
        v       = v*d;
        Y       = u*sqrt(s(1)/n); 
          
        %-Set in structure
        %--------------------------------------------------------------------------
        xY.y    = y;
        xY.u    = Y;
        xY.v    = v;
        xY.s    = s;
 
        display_VOI(xY,i);
        
end

%==========================================================================
function txt = myDataCursorAtlas(obj,evt,H)

pos = get(evt,'Position');

if H.cursor_mode == 1
  txt = {'Desikan DK40'};
else
  txt = {'Destrieux 2009'};
end

i = ismember(get(H.patch(1),'vertices'),pos,'rows');
node = find(i);
ind = 1;

if isempty(node)
  i = ismember(get(H.patch(3),'vertices'),pos,'rows');
  node = find(i);
  ind = 2;
end

rdata_pos = H.rdata{H.cursor_mode}(node,ind);

rcsv = H.rcsv{H.cursor_mode};

for j=2:size(rcsv,1)
  if rdata_pos == rcsv{j,1}
    txt = {txt{:} [H.S{ind}.side ' ' rcsv{j,2}]};
    j = size(rcsv,1);
  end
end

if 0
hMe = findobj(H.axis,'Tag','CrossBar');
if ~isempty(hMe)
    ws = warning('off');
    spm_XYZreg('SetCoords',pos,get(hMe,'UserData'));
    warning(ws);
end
end

%==========================================================================
function myDeleteFcn(obj,evt,renderer)
try rotate3d(get(obj,'parent'),'off'); end
set(ancestor(obj,'figure'),'Renderer',renderer);

%==========================================================================
function s=remove_zeros(s)

pos = length(s);
while pos > 1
  if strcmp(s(pos),'0')
    s(pos)=''
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

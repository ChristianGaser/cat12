function varargout = cat_surf_results(action, varargin)
% Visualise results for both hemispheres of surface-based analysis 
% (preferable on log P-maps). 
%
% FORMAT y = cat_surf_results('Disp',Surface)
% Surface  - a GIfTI filename/object or patch structure
%
% y        - adjusted, predicted or raw response
%
%
% Further actions for batch mode: 
%  * cat_surf_results('disp',filename(s)) 
%  Init display for selected file(s)
%
%  * cat_surf_results('rgb',filename(s)) 
%  Init display for RGB overlay of up to 3 selected file(s)
%
%  * cat_surf_results('batch',job) 
%  See cat_conf_stools.
%
%  * cat_surf_results('surface',1..4) 
%  Select surface type.
%  1 - FSaverage, 2 - Inflated, 3 - Dartel, 4 - Flatmap
%
%  * cat_surf_results('texture',0..3)
%  Select surface underlay. 
%  0 - no transparancy, 1 - mean curvature, 2 - sulcal depth, 3 - nothing
%
%  * cat_surf_results('border',0..3)
%  Select border overlay. 
%  0 - no border, 1 - Desikan-Killiany DK40, 2- Destrieux 2009, 3 - HCP Multi-Modal Parcellation
%
%  * cat_surf_results('atlas',1..3,min_extent,min_overlap)
%  Select atlas labeling and apply minimum cluster extent and atlas overlap. 
%  1 - Desikan-Killiany DK40, 2- Destrieux 2009, 3 - HCP Multi-Modal Parcellation
%
%  * cat_surf_results('view',1..3)
%  Select render view.
%  1 - topview, 2 - bottomview, 3 - sideview
%
%  * cat_surf_results('ov',structure)
%  Define OV fields for cat_vol_slice_overlay.m.
%  Example:
%  struct('atlas','cat12_neuromorphometrics','slices_str','-55:5:60','xy',[5 5],'transform','axial')
%
%  * cat_surf_results('colorbar')
%  Disable colorbar.
%
%  * cat_surf_results('colormap',1..4)
%  Select overlay colormap.
%  1 - jet, 2 - hot, 3 - hsv, 4 - cold-hot
%
%  * cat_surf_results('ctitle',string)
%  Set title for colorbar.
%
%  * cat_surf_results('colormap','customized',cmap)
%  use customized overlay colormap.
%  cmap can be a nx3 variable with the colormap values
%  or a function that creates a colormap (e.g. jet)
%
%  * cat_surf_results('invcolormap',0..1)
%  Default (0) or inverts colormap (1). Toggles without input.
%  
%  * cat_surf_results('background',0..2)
%  White (1) or black (0|2) background. Toggles without input.
%
%  * cat_surf_results('showfilename',0..1); 
%  Show (1) or not show (0) surface name in figure. Toggles without input.
%
%  * cat_surf_results('clim',[mn mx]);
%  Define clim to define view range. 
%
%  * cat_surf_results('clip',[mn mx]);
%  Define clip to limit view range. This is the range that is not displayed
%  with colors and thresholded.
%
%  * cat_surf_results('threshold',val);
%  Define statistical threshold. 
%  0 - none, -log10(0.05) - 0.05, 2 - 0.01, 3 - 0.001
%
%  * cat_surf_results('forward');
%  * cat_surf_results('next');
%  Select next file
%
%  * cat_surf_results('backward');
%  * cat_surf_results('prev');
%  Select previous file
%
%  * cat_surf_results('transparency');
%  Disable transparency. 
%
%  * cat_surf_results('hide_neg',0..1); 
%  Hide negative results (1) or show everything (0). Toggles without input.
%
%  * cat_surf_results('print',fpart); 
%  Save render as png image file. The parameter fpart has to be a structure
%  with the fields 'outdir','prefix', and 'suffix' but can be optionally
%  skipped to estimate the filename automatically.
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

global H OV y x

% ignore this warning writing gifti with int32 (eg. cat_surf_createCS:580 > gifti/subsref:45)
warning off MATLAB:subscripting:noSubscriptsSpecified

%-Input parameters
%--------------------------------------------------------------------------
if ~nargin, clearvars -global H OV x y; action = 'Disp'; end

if ~ischar(action)
  varargin = {action varargin{:}};
  action = 'batch';
end


%-Action
%--------------------------------------------------------------------------
switch lower(action)
  
  %-Display or display as RGB overlay
  %======================================================================
  case {'disp','rgb'}
    
    % if nifti file is given use select_data function to prepare temporary
    % mesh
    if nargin
      [pth,nam,ext] = fileparts(varargin{1}(1,:));
      if strcmp(ext,'.nii')
        if strcmp(lower(action),'rgb')
          fprintf('RGB overlay is not supported for volumes.\n');
          return
        end
        select_data([],[],varargin{1});
        return
      end
    end
    
    % close old figures
    if ishandle(12),close(12);end;if ishandle(22),close(22);end;

    % remove any existing data
    if exist('H','var')
      if isfield(H,'S'), H = rmfield(H,'S'); end
      if isfield(H,'pos'), H = rmfield(H,'pos'); end
    end


    % set start values
    y              = [];
    H.clip         = [];
    H.clim         = [];
    H.XTick        = [];
    H.bkg_col      = [0 0 0];
    H.show_inv     = 0;
    H.no_neg       = 0;
    H.show_transp  = 1; 
    H.n_surf       = 1;
    H.thresh_value = 0;
    H.cursor_mode  = 1;
    H.texture_mode = 1;
    H.border_mode  = 0;
    H.str32k       = '';
    H.SPM_found    = 1;
    H.surf_sel     = 1;
    H.results_sel  = 1;
    H.isfsavg      = 1;
    H.fixscl       = 0;
    H.col          = [.8 .8 .8; 1 .5 .5];
    H.cmap_col     = jet(256);
    H.FS           = cat_get_defaults('extopts.fontsize');
        
    clearvars -global OV
    % positions
    WS = spm('Winsize', 'Graphics');
    H.SS = get(0, 'Screensize');
    if 2.6 * WS(3) > H.SS(3)
      WS(3) = WS(3) / (2.6 * WS(3) / H.SS(3));
    end
    
    % result window with 5 surface views and alternative positions without top view and  only with lateral views
    eh = 0.03;
    H.viewpos = {[0.025 0.450+eh 0.375 0.375;  0.025 0.450+eh 0.375 0.375;  0.025 2.000+eh 0.375 0.375],... % lh medial
                 [0.025 0.025+eh 0.375 0.375;  0.025 0.025+eh 0.375 0.375;  0.175 0.350+eh 0.175 0.350],... % lh lateral
                 [0.600 0.450+eh 0.375 0.375;  0.600 0.450+eh 0.375 0.375;  0.600 2.000+eh 0.375 0.375],... % rh medial
                 [0.600 0.025+eh 0.375 0.375;  0.600 0.025+eh 0.375 0.375;  0.675 0.350+eh 0.175 0.350],... % rh lateral
                 [0.300 0.150+eh 0.400 0.500;  0.300 2.000+eh 0.400 0.500;  0.300 2.000+eh 0.400 0.500],... % lh+rh top
                 [0.400 0.725+eh 0.200 0.225;  0.400 0.300+eh 0.200 0.225;  0.400 0.725+eh 0.200 0.225]};   % data plot
    
    % change size and position of flatmaps for >= R20014b
    if spm_check_version('matlab', '8.4') >= 0
      H.viewpos{2}(3, :) = [-0.075 0.150 0.650 0.650]; % lh lateral
      H.viewpos{4}(3, :) = [0.425 0.150 0.650 0.650];  % rh lateral
    end
    
    % figure 1 with result window
    H.pos{1} = struct( ...
      'fig',  [10 10 round(2.6*WS(3)) WS(3)], ... % figure
      'cbar', [0.400 -0.125 0.200 0.300; 0.440 0.025 0.120 0.120]);% colorbar
    
    % figure 2 with GUI
    H.pos{2} = struct(...
      'prev',   [0.050 0.965 0.050 0.025],...
      'sel',    [0.100 0.935 0.800 0.050],...
      'next',   [0.900 0.965 0.050 0.025],...
      'nam',    [0.050 0.875 0.900 0.050],...
      'surf',   [0.050 0.800 0.425 0.050],'mview',   [0.525 0.800 0.425 0.050],... 
      'text',   [0.050 0.750 0.425 0.050],'thresh',  [0.525 0.750 0.425 0.050],... 
      'cmap',   [0.050 0.700 0.425 0.050],'atlas',   [0.525 0.700 0.425 0.050],...
      'cursor', [0.050 0.650 0.425 0.050],'border',  [0.525 0.650 0.425 0.050],...
      'nocbar', [0.050 0.600 0.425 0.050],'transp',  [0.525 0.600 0.425 0.050],... 
      'info',   [0.050 0.550 0.425 0.050],'bkg',     [0.525 0.550 0.425 0.050],... 
      'inv',    [0.050 0.500 0.425 0.050],'hide_neg',[0.525 0.500 0.425 0.050],...
      'xy',     [0.050 0.325 0.425 0.100],'trans',   [0.525 0.325 0.425 0.100],...
      'labels', [0.050 0.345 0.425 0.030],'str3',    [0.525 0.345 0.100 0.030],'slice',   [0.625 0.345 0.310 0.030],...
      'fixscl', [0.050 0.220 0.425 0.050],'scaling', [0.050 0.220 0.425 0.050],...
      'ovmin',  [0.050 0.125 0.425 0.100],'ovmax',   [0.525 0.125 0.425 0.100],... 
      'save',   [0.050 0.050 0.425 0.050],'close',   [0.525 0.050 0.425 0.050]);
        
    H.figure = figure(23);
    clf(H.figure);
  
    set(H.figure, 'MenuBar', 'none', 'Position', H.pos{1}.fig, ...
      'Name', 'CAT Results', 'NumberTitle', 'off', 'Renderer', 'OpenGL');
      
    H.panel(1) = uipanel('Position',[0 0 2/2.6 1],'units','normalized','BackgroundColor',...
      H.bkg_col,'BorderType','none'); 
    H.panel(2) = uipanel('Position',[2/2.6 0 0.6/2.6 1],'units','normalized','BorderType','none','BackgroundColor',H.col(1,:)); 
    
    % define S structure that contains information for lh and rh
    H.S{1}.name = ''; H.S{1}.side = 'lh';
    H.S{2}.name = ''; H.S{2}.side = 'rh';
     
    
    % Extra button
    if cat_get_defaults('extopts.expertgui')>1
      % Scaling (definition from cat_conf_stools)
      labels = {'SD2','SD4','SD8','%100','%99.99','min-max','0-max'};
      str = [{['Datarange ' char(133)]},labels];
      tmp = {}; 
      for il=1:numel(labels)
        tmp = [tmp {{ @(x)cat_surf_results('clims',x),labels{il} }}]; %#ok<AGROW>
      end

      H.scaling = uicontrol(H.panel(2), ...
        'String', str, 'Units', 'normalized', ...
        'Position', H.pos{2}.scaling, 'Userdata', tmp, ...
        'Style', 'PopUp', 'HorizontalAlignment', 'center', ...
        'Callback', 'spm(''PopUpCB'',gcbo)', ...
        'FontSize',H.FS,...
        'ToolTipString', 'Data range limits', ...
        'Interruptible', 'on', 'Enable', 'off');
    end
      
    H.prev = uicontrol(H.panel(2), ...
      'String', '<', 'Units', 'normalized', ...
      'Position', H.pos{2}.prev, ...
      'Style', 'Pushbutton', 'HorizontalAlignment', 'center', ...
      'Callback', @prev_file, ...
      'FontSize',H.FS,'ForegroundColor','black',...
      'ToolTipString', 'Next file', ...
      'Interruptible', 'on', 'Visible', 'off');
    
    H.next = uicontrol(H.panel(2), ...
      'String', '>', 'Units', 'normalized', ...
      'Position', H.pos{2}.next, ...
      'Style', 'Pushbutton', 'HorizontalAlignment', 'center', ...
      'Callback', @next_file, ...
      'FontSize',H.FS,'ForegroundColor','black',...
      'ToolTipString', 'Next file', ...
      'Interruptible', 'on', 'Visible', 'off');

    % closing all windows
    H.close = uicontrol(H.panel(2), ...
      'String', 'Close', 'Units', 'normalized', ...
      'Position', H.pos{2}.close, ...
      'Style', 'Pushbutton', 'HorizontalAlignment', 'center', ...
      'Callback', 'clear -globalvar H OV;for i=10:30, if ishandle(i); close(i); end; end', ...
      'FontSize',H.FS,'ForegroundColor','red',...
      'ToolTipString', 'Close windows', ...
      'Interruptible', 'on', 'Enable', 'on');

    % select results for lh and rh
    H.sel = uicontrol(H.panel(2), ...
      'String', 'Select Data for Surface Overlay', 'Units', 'normalized', ...
      'Position', H.pos{2}.sel, ...
      'Style', 'Pushbutton', 'HorizontalAlignment', 'center', ...
      'Callback', @select_data, ...
      'FontSize',H.FS,...
      'ToolTipString', 'Select results (up to 24) for both hemispheres (e.g. log-p volume or surface maps)', ...
      'Interruptible', 'on', 'Enable', 'on');
    
    H.save = uicontrol(H.panel(2), ...
      'String', 'Save', 'Units', 'normalized', ...
      'Position', H.pos{2}.save, ...
      'Style', 'Pushbutton', 'HorizontalAlignment', 'center', ...
      'Callback', {@save_image}, ...
      'FontSize',H.FS,...
      'ToolTipString', 'Save png image', ...
      'Interruptible', 'on', 'Enable', 'off');

    str = {'Surface', 'FSaverage', 'Inflated', 'Dartel', 'Flatmap'};
    tmp = {{@select_surf, 1}, ...
           {@select_surf, 2}, ...
           {@select_surf, 3}, ...
           {@select_surf, 4}};
    
    % underlying surface
    H.surf = uicontrol(H.panel(2), ...
      'String', str, 'Units', 'normalized', ...
      'Position', H.pos{2}.surf, 'UserData', tmp, ...
      'Style', 'PopUpMenu', 'HorizontalAlignment', 'center', ...
      'Callback', 'spm(''PopUpCB'',gcbo)', ...
      'FontSize',H.FS,...
      'ToolTipString', 'Underlying Surface', ...
      'Interruptible', 'on', 'Enable', 'off');
    
    str = {'Threshold', 'No threshold', 'P<0.05', 'P<0.01', 'P<0.001'};
    tmp = {{@select_thresh, 0}, ...
           {@select_thresh, -log10(0.05)}, ...
           {@select_thresh, 2}, ...
           {@select_thresh, 3}};
    
    % threshold
    H.thresh = uicontrol(H.panel(2), ...
      'String', str, 'Units', 'normalized', ...
      'Position', H.pos{2}.thresh, 'UserData', tmp, ...
      'Style', 'PopUpMenu', 'HorizontalAlignment', 'center', ...
      'Callback', 'spm(''PopUpCB'',gcbo)', ...
      'FontSize',H.FS,...
      'ToolTipString', 'Threshold', ...
      'Interruptible', 'on', 'Enable', 'off');
    
    str = {'Colormap', 'jet', 'hot', 'hsv', 'cold-hot', 'customized'};
    tmp = {{@select_cmap, 1}, ...
           {@select_cmap, 2}, ...
           {@select_cmap, 3}, ...
           {@select_cmap, 4}, ...
           {@select_cmap, 5}};
    
    % colormap
    H.cmap = uicontrol(H.panel(2), ...
      'String', str, 'Units', 'normalized', ...
      'Position', H.pos{2}.cmap, 'UserData', tmp, ...
      'Style', 'PopUpMenu', 'HorizontalAlignment', 'center', ...
      'Callback', 'spm(''PopUpCB'',gcbo)', ...
      'FontSize',H.FS,...
      'ToolTipString', 'Threshold', ...
      'Interruptible', 'on', 'Enable', 'off');
    
    str = {'Atlas Labeling', 'Desikan-Killiany DK40', 'Destrieux 2009', 'HCP Multi-Modal Parcellation'};
    tmp = {{@select_atlas, 1}, ...
           {@select_atlas, 2}, ...
           {@select_atlas, 3}};
    
    % atlas for labeling
    H.atlas = uicontrol(H.panel(2), ...
      'String', str, 'Units', 'normalized', ...
      'Position', H.pos{2}.atlas, 'UserData', tmp, ...
      'Style', 'PopUpMenu', 'HorizontalAlignment', 'center', ...
      'Callback', 'spm(''PopUpCB'',gcbo)', ...
      'FontSize',H.FS,...
      'ToolTipString', 'Atlas Labeling', ...
      'Interruptible', 'on', 'Enable', 'off');
    
    str = {'Data Cursor', 'Disable data cursor', 'Atlas regions: All Atlases', ...
      'Atlas regions: Desikan-Killiany DK40', 'Atlas regions: Destrieux 2009', ...
      'Atlas region: HCP Multi-Modal Parcellation', 'Plot data at vertex', ...
      'Plot mean data inside cluster', 'Enable/Disable rotate3d'};
    tmp = {{@select_cursor, 0}, ...
           {@select_cursor, 1}, ...
           {@select_cursor, 2}, ...
           {@select_cursor, 3}, ...
           {@select_cursor, 4}, ...
           {@select_cursor, 5}, ...
           {@select_cursor, 6}, ...
           {@select_cursor, 7}};
    
    % data cursor for data plotting and atlas names
    H.cursor = uicontrol(H.panel(2), ...
      'String', str, 'Units', 'normalized', ...
      'Position', H.pos{2}.cursor, 'UserData', tmp, ...
      'Style', 'PopUpMenu', 'HorizontalAlignment', 'center', ...
      'Callback', 'spm(''PopUpCB'',gcbo)', ...
      'FontSize',H.FS,...
      'ToolTipString', 'Data Cursor Mode', ...
      'Interruptible', 'on', 'Enable', 'off');
    
    str = {'View', 'Show top view', 'Show bottom view', 'Show only lateral and medial views'};
    tmp = {{@select_view, 1}, ...
           {@select_view, -1}, ...
           {@select_view, 2}};
    
    % view
    H.mview = uicontrol(H.panel(2), ...
      'String', str, 'Units', 'normalized', ...
      'Position', H.pos{2}.mview, 'UserData', tmp, ...
      'Style', 'PopUpMenu', 'HorizontalAlignment', 'center', ...
      'Callback', 'spm(''PopUpCB'',gcbo)', ...
      'FontSize',H.FS,...
      'ToolTipString', 'Select View', ...
      'Interruptible', 'on', 'Enable', 'off');
    
    str = {'Underlay', 'Mean curvature', 'Sulcal depth', 'Nothing'};
    tmp = {{@select_texture, 1}, ...
           {@select_texture, 2}, ...
           {@select_texture, 3}};
   
    % underlying texture
    H.text = uicontrol(H.panel(2), ...
      'String', str, 'Units', 'normalized', ...
      'Position', H.pos{2}.text, 'UserData', tmp, ...
      'Style', 'PopUpMenu', 'HorizontalAlignment', 'center', ...
      'Callback', 'spm(''PopUpCB'',gcbo)', ...
      'FontSize',H.FS,...
      'ToolTipString', 'Select Underlying Texture', ...
      'Interruptible', 'on', 'Enable', 'off');
    
    str = {'Atlas Border', 'No Overlay', 'Overlay Desikan-Killiany DK40', 'Overlay Destrieux 2009', ...
         'Overlay HCP Multi-Modal Parcellation'};
    tmp = {{@select_border, 0}, ...
           {@select_border, 1}, ...
           {@select_border, 2}, ...
           {@select_border, 3}};
    
    % atlas for border overlay
    H.border = uicontrol(H.panel(2), ...
      'String', str, 'Units', 'normalized', ...
      'Position', H.pos{2}.border, 'UserData', tmp, ...
      'Style', 'PopUpMenu', 'HorizontalAlignment', 'center', ...
      'Callback', 'spm(''PopUpCB'',gcbo)', ...
      'FontSize',H.FS,...
      'ToolTipString', 'Atlas Border Overlay', ...
      'Interruptible', 'on', 'Enable', 'off');
    
    % invert results
    H.inv = uicontrol(H.panel(2), ...
      'String', 'Invert colormap', 'Units', 'normalized', ...
      'BackgroundColor',H.col(1,:),...
      'Position', H.pos{2}.inv, ...
      'Style', 'CheckBox', 'HorizontalAlignment', 'center', ...
      'Callback', {@checkbox_inv}, ...
      'FontSize',H.FS,...
      'ToolTipString', 'Invert results', ...
      'Interruptible', 'on', 'Enable', 'off');
    
    % show only results for pos. contrast
    H.hide_neg = uicontrol(H.panel(2), ...
      'String', 'Hide neg. results', 'Units', 'normalized', ...
      'BackgroundColor',H.col(1,:),...
      'Position', H.pos{2}.hide_neg, ...
      'Style', 'CheckBox', 'HorizontalAlignment', 'center', ...
      'Callback', {@checkbox_hide_neg}, ...
      'FontSize',H.FS,...
      'ToolTipString', 'Hide neg. results', ...
      'Interruptible', 'on', 'Enable', 'off');
    
    % white background
    H.bkg = uicontrol(H.panel(2), ...
      'String', 'White background', 'Units', 'normalized', ...
      'BackgroundColor',H.col(1,:),...
      'Position', H.pos{2}.bkg, ...
      'Style', 'CheckBox', 'HorizontalAlignment', 'center', ...
      'Callback', {@checkbox_bkg}, ...
      'FontSize',H.FS,...
      'ToolTipString', 'White background', ...
      'Interruptible', 'on', 'Enable', 'off');
    
    % transparent view
    H.transp = uicontrol(H.panel(2), ...
      'String', 'Disable transparency', 'Units', 'normalized', ...
      'BackgroundColor',H.col(1,:),...
      'Position', H.pos{2}.transp, ...
      'Style', 'CheckBox', 'HorizontalAlignment', 'center', ...
      'Callback', {@checkbox_transp}, ...
      'FontSize',H.FS,...
      'ToolTipString', 'Disable transparent overlay', ...
      'Interruptible', 'on', 'Enable', 'off');
    
    H.info = uicontrol(H.panel(2), ...
      'String', 'Show filename', 'Units', 'normalized', ...
      'BackgroundColor',H.col(1,:),...
      'Position', H.pos{2}.info, ...
      'Style', 'CheckBox', 'HorizontalAlignment', 'center', ...
      'Callback', {@checkbox_info}, ...
      'FontSize',H.FS,...
      'ToolTipString', 'Show file information in image', ...
      'Interruptible', 'on', 'Enable', 'off');
    
    if cat_get_defaults('extopts.expertgui')<3
      H.nocbar = uicontrol(H.panel(2), ...
        'String', 'Hide colorbar', 'Units', 'normalized', ...
        'BackgroundColor',H.col(1,:),...
        'Position', H.pos{2}.nocbar, ...
        'Style', 'CheckBox', 'HorizontalAlignment', 'center', ...
        'Callback', {@checkbox_nocbar}, ...
        'FontSize',H.FS,...
        'ToolTipString', 'Hide colorbar', ...
        'Interruptible', 'on', 'Enable', 'off');
    else
      str = {'Colorbar', 'none', 'default', 'histogram'};
      tmp = {{@(x) cat_surf_results('colorbar',x), 0}, ...
         {@(x) cat_surf_results('colorbar',x), 1}, ...
         {@(x) cat_surf_results('colorbar',x), 2}};

      % colormap
      H.nocbar = uicontrol(H.panel(2), ...
        'String', str, 'Units', 'normalized', ...
        'Position', H.pos{2}.nocbar, 'UserData', tmp, ...
        'Style', 'PopUp', 'HorizontalAlignment', 'center', ...
        'Callback', 'spm(''PopUpCB'',gcbo)', ...
        'FontSize',H.FS,...
        'ToolTipString', 'Threshold', ...
        'Interruptible', 'on', 'Enable', 'off');
    end
    
    H.fix = uicontrol(H.panel(2), ...
      'String', 'Fix scaling', 'Units', 'normalized', ...
      'BackgroundColor',H.col(1,:),...
      'Position', H.pos{2}.fixscl, ...
      'Style', 'CheckBox', 'HorizontalAlignment', 'center', ...
      'Callback', {@checkbox_fixscl}, ...
      'FontSize',H.FS,...
      'ToolTipString', 'Fix scaling', ...
      'Interruptible', 'on', 'Visible', 'off');    
        
    if nargin >= 2
      
      if isempty(varargin{1}), return; end
      
      H.S{1}.name = varargin{1};
      H.S{2}.name = varargin{1};
      
      try
        Y = spm_data_read(spm_data_hdr_read(H.S{1}.name));
        all_zero = 0;
      catch
        cat_io_cprintf('err','No data in surfaces found or surfaces have different mesh structure (32k vs. 164k).\n');
        if nargout, varargout{1} = []; end
        G = gifti(H.S{1}.name);
        Y = zeros(size(G.vertices,1),1);
        all_zero = 1;
      end
                  
      [pth1, nm1, ext1] = spm_fileparts(H.S{1}.name(1,:));
                    
      % read meshes
      H.S{1}.info = cat_surf_info(H.S{1}.name, 1);
      H.S{2}.info = H.S{1}.info;            
      H.S{1}.info(1).side = 'lh';
      H.S{2}.info(1).side = 'rh';
      H.n_surf = numel(H.S{1}.info);
      if H.S{1}.info(1).nvertices == 64984
        H.str32k = '_32k';
      else
        H.str32k = '';
      end

      if H.n_surf > 3 && strcmp(lower(action),'rgb')
        fprintf('For RGB overlay a maximum of 3 surfaces is allowed.\n');
        return
      end
      
      H.S{1}.info(1).side = 'lh';
      H.S{1}.info(1).Pmesh = fullfile(fileparts(mfilename('fullpath')), ...
          ['templates_surfaces' H.str32k], 'lh.central.freesurfer.gii');
      H.S{2}.info(1).side = 'rh';
      H.S{2}.info(1).Pmesh = fullfile(fileparts(mfilename('fullpath')), ....
          ['templates_surfaces' H.str32k], 'rh.central.freesurfer.gii');

      for ind = 1:2
        H.S{ind}.M = gifti(H.S{ind}.info(1).Pmesh);            
        % get adjacency information
        H.S{ind}.A = spm_mesh_adjacency(H.S{ind}.M);
      end

      H.logP = ones(H.n_surf,1);

      for i=1:H.n_surf
        if isempty(strfind(H.S{1}.info(i).ff, 'log'))
          H.logP(i) = 0;
        end
      end

      H.nY2 = size(Y,1)/2;
      H.S{1}.Y = Y(1:H.nY2, :);
      H.S{2}.Y = Y((H.nY2+1):end, :);

      % correct H.isvol if more than one surface was given
      if isfield(H,'isvol') && H.n_surf ~= numel(H.isvol)
        for i=2:H.n_surf
          H.isvol(i) = H.isvol(1);
        end
      end
      
      % delete temporary files that were created for volume mapping
      for i=1:H.n_surf
        if ~isfield(H,'isvol')
          H.isvol(i) = 0;
        elseif H.isvol(i)
          delete(deblank(H.S{1}.name(i,:)))
        end
      end
      
      % if size of cdata does not fit to mesh size load the underlying
      % mesh instead of fsaverage mesh and divide hemispheres into lh/rh
      if size(Y,1) ~= (size(H.S{1}.M.faces,1)+4)
        Mg = gifti(deblank(H.S{1}.name(1,:)));
        sz_faces2 = size(Mg.faces,1)/2;
        sz_vertices2 = size(Mg.vertices,1)/2;
        H.S{1}.M.faces = Mg.faces(1:sz_faces2,:);
        H.S{1}.M.vertices = Mg.vertices(1:sz_vertices2,:);
        H.isfsavg = 0;
      end

      if ~H.isfsavg
        H.S{2}.M.faces = Mg.faces(((sz_faces2+1):end),:) - sz_vertices2;
        H.S{2}.M.vertices = Mg.vertices((sz_vertices2+1):end,:);
      end
                          
      % rescue original name for later result selection
      H.S1 = H.S{1};
      H.S2 = H.S{2};
      
      H.view = 1;
      H.show_transp = 1;
      H.white_bgk = 0;
      H.show_info = 0;
      if all_zero
        H.disable_cbar = 1;
      else
        H.disable_cbar = 0;
      end
      
      % result selection or RGB overlay if more than one result was loaded
      if H.n_surf > 1
        % pre-select 1st mesh if we cannot use RGB overlay
        if ~strcmp(lower(action),'rgb')
          sel = 1;
          if isempty(H.S1.name)
            spm('alert*','Do not mix meshes with different resolutions (i.e. 164k vs. 32k)');
            if nargout, varargout{1} = []; end
            return
          end
          H.S{1}.name = H.S1.name(sel, :);
          H.S{2}.name = H.S2.name(sel, :);
          H.S{1}.Y = H.S1.Y(:, sel);
          H.S{2}.Y = H.S2.Y(:, sel);
        end
        
        % delete old selection ui
        delete(H.sel);
        H = rmfield(H, 'sel');
        
        str = cell(1, H.n_surf + 2);
        tmp = cell(1, H.n_surf + 1);
        str{1} = 'Select Result ';
        [C,C2] = spm_str_manip( H.S1.name , 'C');
        for s = 1:H.n_surf
          str{s + 1} = spm_str_manip(H.S1.name(s,:), 'k60d');
          tmp{s} = {@select_results, s};
        end
        
        % print selected filename
        H.nam = axes('Parent', H.panel(2), 'Position', H.pos{2}.nam);
        cla(H.nam);
        axis(H.nam, 'off')
        text(0.5, 0.5, spm_str_manip(H.S{1}.name, 'k60d'), 'Parent', H.nam, 'Interpreter', 'none', ...
          'FontSize', H.FS, 'HorizontalAlignment', 'center');
        
        % set # of surfaces back to "1" if we cannot use RGB overlay
        if ~strcmp(lower(action),'rgb'), H.n_surf = 1; end
 
        % new selection ui
        str{s + 2} = 'Select new data';
        tmp{s + 1} = {@select_data};
        H.sel = uicontrol(H.panel(2), ...
          'String', str, 'Units', 'normalized', ...
          'Position', H.pos{2}.sel, 'UserData', tmp, ...
          'Style', 'Popup', 'HorizontalAlignment', 'center', ...
          'Callback', 'spm(''PopUpCB'',gcbo)', ...
           'FontSize',H.FS,...
          'ToolTipString', 'Select results', ...
          'Interruptible', 'on', 'Enable', 'on');
          
        % enable fixing of scale
        set(H.fix, 'Visible', 'on');
      end
      
      display_results_all;
      
      H.SPM_found = 1;
      SPM_name = fullfile(H.S{1}.info(1).pp, 'SPM.mat');
      
      % SPM.mat exist?
      if ~isempty(H.S{1}.name)
        H.SPM_found = 0;
      end

      % Don't allow plot functions for RGB or volume maps or if SPM.mat was not found
      if (H.n_surf > 1 && H.SPM_found) || H.isvol(1) || strcmp(lower(action),'rgb')
        str = {'Data Cursor', 'Disable data cursor', 'Atlas regions: All atlases',...
          'Atlas regions: Desikan-Killiany DK40', 'Atlas regions: Destrieux 2009',...
          'Atlas region: HCP Multi-Modal Parcellation', 'Enable/Disable rotate3d'};
        tmp = {{@select_cursor, 0}, ...
               {@select_cursor, 1}, ...
               {@select_cursor, 2}, ...
               {@select_cursor, 3}, ...
               {@select_cursor, 4}, ...
               {@select_cursor, 7}};
        
        set(H.cursor,'String', str, 'UserData', tmp);
      end
      
      % enable some menus only if mesh data can be assumed to be resampled
      if (length(H.S{1}.Y) == 32492 || length(H.S{1}.Y) == 163842 || length(H.S{1}.Y) == 40962) && H.isfsavg
        set(H.surf,   'Enable', 'on');
        set(H.text,   'Enable', 'on');
        set(H.cursor, 'Enable', 'on');
        set(H.border, 'Enable', 'on');
      end
      
      set(H.save,   'Enable', 'on');
      set(H.mview,  'Enable', 'on');
      if ~all_zero
        set(H.nocbar, 'Enable', 'on');
        set(H.cmap,   'Enable', 'on');
        set(H.inv,    'Enable', 'on');
        set(H.transp, 'Enable', 'on');
      else
        set(H.slider_min, 'Visible','off');
        set(H.slider_max, 'Visible','off');
      end
      set(H.bkg,    'Enable', 'on');
      set(H.info,   'Enable', 'on');
      if isfield(H,'scaling') &&  isvalid(H.scaling)
        set(H.scaling, 'Enable', 'on');
      end
            
      H.rdata{1} = [];
      H.rdata{2} = [];
      H.rdata{3} = [];
      for ind = 1:2
        atlas_name = fullfile(fileparts(mfilename('fullpath')), ['atlases_surfaces' H.str32k], ...
        [H.S{ind}.info(1).side '.aparc_DK40.freesurfer.annot']);
        [vertices, rdata0, colortable, rcsv1] = cat_io_FreeSurfer('read_annotation', atlas_name);
        H.rdata{1} = [H.rdata{1} rdata0];
        atlas_name = fullfile(fileparts(mfilename('fullpath')), ['atlases_surfaces' H.str32k], ...
        [H.S{ind}.info(1).side '.aparc_a2009s.freesurfer.annot']);
        [vertices, rdata0, colortable, rcsv2] = cat_io_FreeSurfer('read_annotation', atlas_name);
        H.rdata{2} = [H.rdata{2} rdata0];
        atlas_name = fullfile(fileparts(mfilename('fullpath')), ['atlases_surfaces' H.str32k], ...
        [H.S{ind}.info(1).side '.aparc_HCP_MMP1.freesurfer.annot']);
        [vertices, rdata0, colortable, rcsv3] = cat_io_FreeSurfer('read_annotation', atlas_name);
        H.rdata{3} = [H.rdata{3} rdata0];
      end
      H.rcsv{1} = rcsv1;
      H.rcsv{2} = rcsv2;
      H.rcsv{3} = rcsv3;
      
      H.dcm_obj = datacursormode(H.figure);
      set(H.dcm_obj, 'Enable', 'on', 'SnapToDataVertex', 'on',...
        'DisplayStyle', 'datatip', 'Updatefcn', {@myDataCursorAtlas, H});
      % sometimes interpreter cannot be set for some Matlab versions
      try set(dcm_obj, 'Interpreter', 'none'); end
      
    end
    
    if nargout, varargout{1} = y; end
    
  %-ColourBar
  %======================================================================
  case {'colourbar', 'colorbar'}
    % RD202003 Colorbars with histogram does not work stable
    if 0 %nargin>1
      if varargin{1} == 2
        cat_surf_results('hist',1); 
      else
        cat_surf_results('hist',0); 
      end
      if varargin{1} == get(H.nocbar, 'Value')  
        cat_surf_results('colorbar');
      end
    else
      set(H.nocbar, 'Value', ~get(H.nocbar, 'Value') );
      checkbox_nocbar;
      %if ~get(H.nocbar, 'Value')
      %  cat_surf_results('hist',0); 
      %end
    end
     %{  
    %if isempty(varargin), varargin{1} = gca; end
    if length(varargin) == 1, varargin{1} = 'on'; end
    switch varargin{1}
      case 0, varargin{1} = 'off';
      case 1, varargin{1} = 'on'; 
      case 2, varargin{1} = 'on'; cat_surf_results('hist'); 
    end
    
    %H = getHandles(varargin{1});
    d = getappdata(H.patch(1), 'data');
    col = getappdata(H.patch(1), 'colourmap');
    if strcmpi(varargin{1}, 'off')
      if isfield(H, 'colourbar') && ishandle(H.colourbar)
        delete(H.colourbar);
        H = rmfield(H, 'colourbar');
        setappdata(H.axis, 'handles', H);
      end
      return;
    end
    if isempty(d) || ~any(d(:)), varargout = {H}; return; end
    if isempty(col), col = jet(256); end
    if ~isfield(H, 'colourbar') || ~ishandle(H.colourbar)
      %      H.colourbar = colorbar('peer',gca,'NorthOutside');
      H.colourbar = colorbar('NorthOutside');
      set(H.colourbar, 'Tag', '');
      set(get(H.colourbar, 'Children'), 'Tag', '');
    end
    c(1:size(col, 1), 1, 1:size(col, 2)) = col;
    ic = findobj(H.colourbar, 'Type', 'image');
    clim = getappdata(H.patch(1), 'clim');
    if isempty(clim), clim = [false NaN NaN]; end
    
    if size(d, 1) > size(d, 2), d = d'; end
    
    % Update colorbar colors if clipping is used
    H.clip = getappdata(H.patch(1), 'clip');
    if ~isempty(H.clip)
      if ~isnan(H.clip(2)) && ~isnan(H.clip(3))
        ncol = length(col);
        col_step = (clim(3) - clim(2)) / ncol;
        cmin = max([1, ceil((H.clip(2) - clim(2)) / col_step)]);
        cmax = min([ncol, floor((H.clip(3) - clim(2)) / col_step)]);
        col(cmin:cmax, :) = repmat([0.5 0.5 0.5], (cmax - cmin + 1), 1);
        c(1:size(col, 1), 1, 1:size(col, 2)) = col;
      end
    end
    if H.n_surf > 1
      set(ic, 'CData', c(1:H.n_surf, :, :));
      set(ic, 'YData', [1 H.n_surf]);
      set(H.colourbar, 'YLim', [1 H.n_surf]);
      set(H.colourbar, 'YTickLabel', []);
    else
      set(ic, 'CData', c);
      clim = getappdata(H.patch(1), 'clim');
      if isempty(clim), clim = [false min(d) max(d)]; end
      set(ic, 'YData', clim(2:3));
      set(H.colourbar, 'YLim', clim(2:3));
    end
    setappdata(H.axis, 'handles', H);
    
    if nargout, varargout{1} = y; end
     %}
    
  %-ColourMapTitle
  %======================================================================
  case {'ctitle'}

    if ~H.disable_cbar
      
      title(H.cbar, varargin{1}, 'Color', 1 - H.bkg_col);
    end

  %-ColourMap
  %======================================================================
  case {'colourmap', 'colormap'}
    if isempty(varargin), varargin{1} = gca; end
    if isobject(varargin{1})
      H = getHandles(varargin{1});
      if length(varargin) == 1
        varargout = {getappdata(H.patch(1), 'colourmap')};
        return;
      else
        setappdata(H.patch(1), 'colourmap', varargin{2});
        d = getappdata(H.patch(1), 'data');
        H = updateTexture(H, d);
      end
      if nargin > 1
        colormap(varargin{2});
      end
    else
      cm = varargin{1}; 
      switch cm
      case {1,2,3,4,5}, cmap = cm; 
      case 'jet',       cmap = 1; 
      case 'hot',       cmap = 2; 
      case 'hsv',       cmap = 3; 
      case 'cold-hot',  cmap = 4; 
      case 'customized',cmap = 5; 
      otherwise
        spm('alert*','Unknown colormap');
        if nargout, varargout{1} = []; end
        return
      end      
      if nargin==3
        select_cmap(cmap,varargin{2});
      else
        select_cmap(cmap);
      end
    end
    
    
  %-FileSelection forward
  %======================================================================
  case {'forward','next'}

    next_file([],[]);

  %-FileSelection backward
  %======================================================================
  case {'backward','previous'}

    prev_file([],[]);

  %-OV
  %======================================================================
  case 'ov'      
    
    OV = varargin{1};

    if isfield(H,'Pvol_sel')
      H.OV = OV;
      H = update_slice_overlay(H);
    else
      disp('This option can only be used for volume rendering.')
    end

  %-CLim
  %======================================================================
  case 'clim'      
    
    caxis(H.axis, varargin{1});
    for ind = 1:5
      setappdata(H.patch(ind), 'clim', [true varargin{1}]);
      d = getappdata(H.patch(ind), 'data');
      col = getappdata(H.patch(ind), 'col');
      H = updateTexture(H, ind, d, col, H.show_transp);
    end    

    if ~H.isvol(H.results_sel)
      set(H.atlas, 'Enable', 'on');
    end

    if ~H.disable_cbar
      H = show_colorbar(H);
    end

    if isfield(H,'Pvol_sel')
      H = update_slice_overlay(H);
    end

%-CLip
  %======================================================================
  case 'clip'

    if nargin < 2
      c = getappdata(H.patch(5), 'clip');
      if ~isempty(c), c = c(2:3); end
      varargout = {c};
      return;
    else
      if isempty(varargin{1}) || any(~isfinite(varargin{1}))
        H.clip = [false NaN NaN];
      else
        H.clip = [true varargin{1}];
      end
      for ind = 1:5
        setappdata(H.patch(ind), 'clip', H.clip);
      end
      for ind = 1:5
        col = getappdata(H.patch(ind), 'col');
        H = updateTexture(H, ind);
      end
    end

    if ~H.disable_cbar
      H = show_colorbar(H);
    end

    if isfield(H,'Pvol_sel')
      H = update_slice_overlay(H);
    end

  %-CLims
  %======================================================================
  case 'clims'
    if nargin>1, H.datascale    = varargin{1}; end
    if nargin>2, H.datascaleval = varargin{2}; end
    
    c = getappdata(H.patch(5), 'data');
    %%
    if isfield(H,'datascale')
    switch H.datascale
      case 'default'
      H.clim(2:3) = [max(min(c(:)),0) max(c(:))];
      case {'minmax','min-max'}
      H.clim(2:3) = [min(c(:)) max(c(:))];
      case {'0max','0-max'}
      H.clim(2:3) = [0 max(c(:))];
      otherwise
      switch H.datascale(1)
        case '%'
        if str2double(H.datascale(2:end))>0 && str2double(H.datascale(2:end))<100 
          H.datascaleval = str2double(H.datascale(2:end))/100;
        end
        [cc,cx] = cat_stat_histth(c(c(:)~=0),H.datascaleval); %#ok<ASGLU>
        H.clim(2:3) = cx;
        case 'M'
        if str2double(H.datascale(4:end))>0  
          H.datascaleval = str2double(H.datascale(4:end));
        end
        mv = cat_stat_nanmean(c(:)); sv = cat_stat_nanstd(c(:)); 
        H.clim(2:3) = [ max(min(c(:)),mv - H.datascaleval * sv) , min(max(c(:)), mv + H.datascaleval * sv)];
        case 'S'
        if str2double(H.datascale(3:end))>0  
          H.datascaleval = str2double(H.datascale(3:end));
        end
        mv = cat_stat_nanmean(c(:)); sv = cat_stat_nanstd(c(:)); 
        H.clim(2:3) = [ mv - H.datascaleval * sv ,  mv + H.datascaleval * sv];
        case 'C'
          nclim = str2num(H.datascale(2:end)); %#ok<ST2NM>
          if numel(nclim)~=2, error('incorrect clim input'); end
          H.clim(2:3) = nclim;
        otherwise
          error('unkown H.datascale %s.\n',H.datascale);
      end
    end
    end

    %% update textures of each patch
    for j=1:numel(H.patch) 
      setappdata(H.patch(j), 'clim',H.clim);
      H = updateTexture(H, j); 
    end
    if isvalid(H.str_min)
      set(H.str_min, 'String', sprintf('%g',H.clim(2)));
      set(H.str_max, 'String', sprintf('%g',H.clim(3)));
    end
    show_colorbar(H); 
  
    
  %- show histogram
  %======================================================================
  case 'hist'
    if nargin>1
      if varargin{1}~=any(isempty(findobj('tag','cat_surf_results_hist'))) 
        cat_surf_results('hist')
      end
    else
    
      if numel(H.patch)>=5 && H.patch(1).isvalid &&  H.patch(3).isvalid &&  H.patch(5).isvalid

        if nargin>1, draw = varargin{1}; else, draw = ~any(isempty(findobj('tag','cat_surf_results_hist'))); end

        % move elements if histogram is added or removed
        if draw==0 || draw~=2
        top = H.patch(5).Parent; 
        pos = get(top,'Position'); 
        set(top,'Position',pos + sign(any(isempty(findobj('tag','cat_surf_results_hist')))-0.5) * [0 0.13 0 0]);
        end

        if any(isnan(H.clim)), cat_surf_results('clims','default'); end

        % draw/update histgram / or remove it
        mode = 0; 
        if any(isempty(findobj('tag','cat_surf_results_hist'))) || draw==2 
        if draw==2
          delete(findobj('tag','cat_surf_results_hist'));
          if isfield(H,'hist'); H = rmfield(H,'hist'); end
        end

        % print colors (red, green/dark-green
        color  = {[1 0 0],[0 1 0]/(1+H.bkg_col(1))};  
        linet  = {'-','--'};
        % position of the right and the left text box
        tabpos = 0.010; 
        if mode
          tpos = {[0.49 tabpos 0.065 0.03],[0.565 tabpos 0.065 0.03],[0.4 tabpos 0.12 0.03]};
        else
          tpos = {[0.4 tabpos 0.12 0.06],[0.445 tabpos 0.065 0.06],[0.497 tabpos 0.065 0.06],[0.55 tabpos 0.065 0.06],[0.60 tabpos 0.065 0.06]};
        end
        % histogram axis 
        H.histax = axes('Parent', H.panel(1), 'Position', [0.4 0.102 0.20 0.15],'Visible', 'off', 'tag','cat_surf_results_hist'); 
        try
          xlim(H.histax,H.clim(2:3).*[1 1+eps]);
        catch
          disp(1);
        end
        hold on;
        % standard text
        if mode 
          H.dtxt(3).ax  = axes('Parent', H.panel(1), 'Position',tpos{3}, 'Visible', 'off','tag','cat_surf_results_text');
          H.dtxt(3).txt = text(0,1,sprintf('%s %s %s: \n%s - %s:\n%s:', 'mean',char(177),'std','min','max','median'),'color',[0.5 0.5 0.5],'Parent',H.dtxt(3).ax);
        else
          H.dtxt(3).ax(1) = axes('Parent', H.panel(1), 'Position',tpos{1}, 'Visible', 'off','tag','cat_surf_results_text');
          H.dtxt(3).ax(2) = axes('Parent', H.panel(1), 'Position',tpos{2}, 'Visible', 'off','tag','cat_surf_results_text');
          H.dtxt(3).ax(3) = axes('Parent', H.panel(1), 'Position',tpos{3}, 'Visible', 'off','tag','cat_surf_results_text');
          H.dtxt(3).ax(4) = axes('Parent', H.panel(1), 'Position',tpos{4}, 'Visible', 'off','tag','cat_surf_results_text');
          H.dtxt(3).ax(5) = axes('Parent', H.panel(1), 'Position',tpos{5}, 'Visible', 'off','tag','cat_surf_results_text');
          text(0,1,sprintf('side'),'color',[0.5 0.5 0.5],'Parent',H.dtxt(3).ax(1));
          text(0,1,sprintf('\n\nleft'),'color',color{1},'Parent',H.dtxt(3).ax(1));
          text(0,1,sprintf('\n\n\n\nright'),'color',color{2},'Parent',H.dtxt(3).ax(1));
          text(0,1,'min','color',[0.5 0.5 0.5],'HorizontalAlignment','center','Parent',H.dtxt(3).ax(2));
          text(0,1,['mean ' char(177) ' std'],'color',[0.5 0.5 0.5],'HorizontalAlignment','center','Parent',H.dtxt(3).ax(3));
          text(0,1,'median','color',[0.5 0.5 0.5],'HorizontalAlignment','center','Parent',H.dtxt(3).ax(4));
          text(0,1,'max','color',[0.5 0.5 0.5],'HorizontalAlignment','right','Parent',H.dtxt(3).ax(5));
        end
        for i=1:2
          side = getappdata(H.patch( i*2 - 1 ), 'data');

          if ~all(isnan(side))
            % histogram plot may fail due to NAN or whatever ...
            [d,h] = hist( side(~isinf(side(:)) & ~isnan(side(:)) & side(:)<3.4027e+38 & side(:)>-3.4027e+38 & side(:)<H.clim(3) & side(:)>H.clim(2) ), ...
              H.clim(2) : diff(H.clim(2:3))/100 : H.clim(3) );
            d = d./numel(side);
            % plot histogram line and its median
            med = cat_stat_nanmedian(side(:));
            quantile = [h(find(cumsum(d)/sum(d)>0.25,1,'first')),h(find(cumsum(d)/sum(d)>0.75,1,'first'))]; 
            % print histogram
            line(H.histax,h,d,'color',color{i},'LineWidth',1);
            % print median
            line(H.histax,[med med],[0 d(find(h>=med,1,'first'))],'color',color{i},'linestyle',linet{i});
            % print quantile 
            if numel(quantile)>1
              fill(H.histax,[quantile(1)   quantile(2)   quantile(2)      quantile(1)],...
                    max(d)*(0.08*[i i (i+1) (i+1)] + 0.16),color{i});
              %
              if mode
              H.dtxt(i).ax  = axes('Parent', H.panel(1), 'Position',tpos{i}, 'Visible', 'off','tag','cat_surf_results_text');
              H.dtxt(i).txt = text(0,1,sprintf('%10.3f %s %0.3f\n%10.3f - %0.3f\n',...
                cat_stat_nanmean(side(:)),char(177),cat_stat_nanstd(side(:)),...
                min(side(:)),max(side(:)),med),...
                'color',color{i},'HorizontalAlignment','center','Parent',H.dtxt(i).ax);
              else
              text(0,1,sprintf('%s%0.3f',sprintf(repmat('\n',1,i*2)),min(side(:))),'color',color{i},'HorizontalAlignment','center','Parent',H.dtxt(3).ax(2));
              text(0,1,sprintf('%s%0.3f%s%0.3f',sprintf(repmat('\n',1,i*2)),cat_stat_nanmean(side(:)),char(177),...
                cat_stat_nanstd(side(:))),'color',color{i},'HorizontalAlignment','center','Parent',H.dtxt(3).ax(3));
              text(0,1,sprintf('%s%0.3f',sprintf(repmat('\n',1,i*2)),cat_stat_nanmedian(side(:))),'color',color{i},'HorizontalAlignment','center','Parent',H.dtxt(3).ax(4));
              text(0,1,sprintf('%s%0.3f',sprintf(repmat('\n',1,i*2)),max(side(:))),'color',color{i},'HorizontalAlignment','right','Parent',H.dtxt(3).ax(5));
              end
            end
          end
        end

        else
        delete(findobj('tag','cat_surf_results_text'));
        delete(findobj('tag','cat_surf_results_hist'));
        if isfield(H,'hist'); H = rmfield(H,'hist'); end
        end 

      end
    end
    
    
    
  %- set surface 
  %======================================================================
  case 'surface'
    surface = varargin{1};
    if any(surface == 1:4)
      select_surf(surface);
    end
     
    
  %- set texture 
  %======================================================================
  case 'texture'
    texture = varargin{1};
    if any(texture == 1:3)
      select_texture(texture);
    elseif texture == 0
      cat_surf_results('transparency',0);
    end
  
    
  %- set border 
  %======================================================================
  case 'border'
    border_mode = varargin{1};
    if any(border_mode == 0:3)
      select_border(border_mode);
    end
  
  %- set atlas
  %======================================================================
  case 'atlas'
    atlas_mode = varargin{1};
    if nargin > 2
      min_extent = varargin{2};
    else
      min_extent = 1;
    end
    if nargin > 3
      min_overlap = varargin{3};
    else
      min_overlap = 1;
    end
    
    if any(atlas_mode == 1:3)
      select_atlas(atlas_mode, min_extent, min_overlap);
    end

  %- set view
  %======================================================================
  case 'view'
    view = varargin{1};
    if any(view == [1,-1,2])
      select_view(view);
    end
  
    
  %- set background
  %======================================================================
  case 'background'
    if nargin>1
      switch varargin{1}
        case {1,'white'}
          if get(H.bkg, 'Value')==0
            H.bkg_val = [1 1 1];
            cat_surf_results('background');
          end
        case {0,2,'black'}
          if get(H.bkg, 'Value')==1
            H.bkg_val = [0 0 0];
            cat_surf_results('background');
          end
        otherwise
           error('Unknown background option'); 
      end
    else
      set(H.bkg, 'Value', ~get(H.bkg, 'Value') );
      checkbox_bkg;
    end
    
    
  %- set showfilename
  %======================================================================
  case 'showfilename'
    if nargin>1
      if varargin{1} ~= get(H.bkg, 'Value')==0  
        cat_surf_results('showfilename');
      end
    else
      set(H.info, 'Value', ~get(H.info, 'Value') );
      checkbox_info;
    end
    
    
  %- set transparency
  %======================================================================
  case 'transparency'
    if nargin>1
      if varargin{1} ~= get(H.transp, 'Value')==0  
        cat_surf_results('transparency');
      end
    else
      set(H.transp, 'Value', ~get(H.transp, 'Value') );
      checkbox_transp;
    end 
    
    
  %- set inverse colormap
  %======================================================================
  case 'invcolormap'
    if nargin>1
      if varargin{1} ~= get(H.inv, 'Value')==0  
        cat_surf_results('invcolormap');
      end
    else
      set(H.inv, 'Value', ~get(H.inv, 'Value') );
      checkbox_inv;
    end
    
    
  %- set threshold
  %======================================================================
  case 'threshold'
    if nargin>1
      select_thresh(varargin{1});
    else
      select_thresh(0);
    end 
    
    
  %- set hide negative values
  %======================================================================
  case 'hide_neg'
    if nargin>1
      if varargin{1} ~= get(H.hide_neg, 'Value')==0  
        cat_surf_results('hide_neg');
      end
    else
      set(H.hide_neg, 'Value', ~get(H.hide_neg, 'Value') );
      checkbox_hide_neg;
    end
   
    
  %- SPM/CAT batch mode
  %======================================================================
  case 'batch'
    job = varargin{1};
    
    if any(~cellfun(@exist,job.cdata)) 
      if size(job.cdata,1)==1,  cat_io_cprintf('err', 'Input file does not exist!\n'); 
      else,                     cat_io_cprintf('err', 'Input files do not exist!\n');
      end
      if nargout, varargout{1}.png = {''}; end
      return
    end
    
    % create window
    select_data([],[],char(job.cdata));
    
    %% set parameter
    % RD202003: not working ... ,'colorbar'
    %if isfield(job,'render') && isfield(job.render,'clims'), job.render.clim = job.render.clims; end
    FN = {'surface','view','texture','transparency','invcolormap','colormap','clims','background','showfilename'}; 
    for fni=1:numel(FN)
      %%
      if isfield(job,'render') && isfield(job.render,FN{fni}) 
        try
          cat_surf_results(FN{fni},job.render.(FN{fni})); 
        end
      end
    end

    FN = {'threshold','hide_neg'}; 
    for fni=1:numel(FN)
      if isfield(job,'stat') && isfield(job.stat,FN{fni})
      cat_surf_results(FN{fni},job.stat.(FN{fni})); 
      end
    end
    
    %% save result
    if isfield(job,'fparts')
      files = cat_surf_results('print',job.fparts);
    else
      files = cat_surf_results('print');
    end
    varargout{1}.png = files; 
    
    % close figure after export
    close(H.figure); 
    clear -globalvar H; 
    
    
  %- save image
  %======================================================================
  case 'print'
    
    if nargin>1
      fparts = varargin{1};
    else
      fparts.outdir = '';
      fparts.prefix = '';
      fparts.suffix = '';
    end
    
    if iscell(fparts.outdir), fparts.outdir = fparts.outdir{1}; end
    [pp,ff] = spm_fileparts( H.S1.name );
    if isempty(fparts.outdir)
      fparts.outdir = pp; 
    end
    [tmp, pathname, ext] = spm_fileparts(pp);
    filename = fullfile(fparts.outdir,[fparts.prefix pathname ext '_' ff fparts.suffix '.png']); %#ok<AGROW>
    if ~exist(fparts.outdir,'dir')
      mkdir(fparts.outdir)
    end
    save_image(1,1,filename);
    
    varargout{1} = filename; 
    
  otherwise   
    error('Unknown action "%s"!\n',action); 
end 

     

%==========================================================================
function Ho = select_thresh(thresh)
global H OV

H.thresh_value = thresh;
H.clip = [true -thresh thresh];

if ~H.isvol
  mn = min(min(H.S{1}.Y(:)), min(H.S{2}.Y(:)));
  mx = max(max(H.S{1}.Y(:)), max(H.S{2}.Y(:)));
else
  mn = H.Pvol_mn{H.results_sel};
  mx = H.Pvol_mx{H.results_sel};
end

if H.logP(H.results_sel) && (H.S{1}.thresh < -log10(0.05))
  set(H.thresh, 'Enable', 'on');
  if mn < 0 && H.n_surf == 1
    set(H.hide_neg, 'Enable', 'on');
%    set(H.hide_neg, 'Value', 0);
  end
else
  set(H.thresh,   'Enable', 'off');
  set(H.hide_neg, 'Enable', 'off');
%  set(H.hide_neg, 'Value', 0);
end

if mn > H.clip(2) && mx < H.clip(3)
  set(H.slider{1}, 'Visible', 'off');
  set(H.slider{2}, 'Visible', 'off');
else
  set(H.slider{1}, 'Visible', 'on');
  set(H.slider{2}, 'Visible', 'on');
end

H.no_neg = get(H.hide_neg, 'Value');

clim = getappdata(H.patch(1), 'clim');
if H.logP(H.results_sel)
  clim(2) = floor(clim(2));
  clim(3) = ceil(clim(3));
end

% rather use NaN values for zero threshold
if thresh == 0
  H.clip = [false NaN NaN];
end

if H.no_neg
  H.clip = [true -Inf thresh];
  clim = [true 0 clim(3)];
  set(H.slider_min, 'Min', 0);
  set(H.slider_max, 'Min', 0);
end

if (thresh > 0)
  set(H.slider_min, 'Value', thresh)
  set(H.str_min, 'String', sprintf('%g',thresh));
end

if mn > -thresh
  H.clim = [true thresh clim(3)];
  set(H.slider_min, 'Value', thresh);
  set(H.str_min, 'String', sprintf('%g',thresh));
else
  H.clim = [true -clim(3) clim(3)];
  set(H.slider_min, 'Min', -clim(3));
  set(H.slider_min, 'Value', -clim(3));
  set(H.str_min, 'String', sprintf('%g',-clim(3)));
end
  
for ind = 1:5
  if mn > -thresh
    setappdata(H.patch(ind), 'clim', [true thresh clim(3)]);
  elseif thresh == 0
    setappdata(H.patch(ind), 'clim', [true -clim(3) clim(3)]);
  end
  
  setappdata(H.patch(ind), 'clip', H.clip);
  setappdata(H.patch(ind), 'clim', H.clim);
  col = getappdata(H.patch(ind), 'col');
  d = getappdata(H.patch(ind), 'data');
  H = updateTexture(H, ind, d, col, H.show_transp);
end

if ~H.isvol(H.results_sel)
  set(H.atlas, 'Enable', 'on');
end

if ~H.disable_cbar
  H = show_colorbar(H);
end

if isfield(H,'Pvol_sel')
  if (mn < 0 && mn > -thresh) && (mx >= 0 && mx < thresh)
    try close(OV.fig); end
    try close(OV.fig_mip); end
  else
    H = update_slice_overlay(H);
  end
end

if nargout, Ho = H; end

%==========================================================================
function disphist(type)
global H; 

i = get(H.sel,'value');
if isfield(H,'patch')
if i==0
  d = getappdata(H.patch(1),'data');
elseif i<=numel(H.patch)
  d = getappdata(H.patch(i),'data');
end
cat_stat_histth(d(d(:)~=0),1,type); 
end


%==========================================================================
function Ho = select_cmap(cmap,cust_cmap)
global H

switch cmap
  case 1
    col = jet(256);
  case 2
    col = hot(256);
  case 3
    col = hsv(256);
  case 4
    col = [1 - hot(128); (hot(128))];
  case 5
    col = [];
    while size(col,2)~=3
      if nargin == 1
        tmp = inputdlg('Colormap','Input');
      else
        tmp{1} = cust_cmap;
      end
      % if a variable is given set is as global
      try eval(['global ' char(tmp{1})]); end
      
      try
        % is it a variable?
        col = eval(tmp{1});
        %  or a function?
        if isempty(col)
          col = feval(tmp{1});
        end
        % get sure that min/max is between 0..1
        mn = min(col(:));
        mx = max(col(:));
        if mx > 1 || mn<0
          col = (col-mn)/(mx-mn);
        end
      catch
        fprintf('Wrong input. You have to define a colormap as variable with size nx3.\n')
      end
    end
end

H.cmap_col = col;

for ind = 1:5
  setappdata(H.patch(ind), 'col', col);
  d = getappdata(H.patch(ind), 'data');
  H = updateTexture(H, ind, d, col, H.show_transp);
end

if ~H.disable_cbar
  H = show_colorbar(H);
end

if isfield(H,'Pvol_sel')
  H = update_slice_overlay(H);
end

H.cmap_col = col;
if nargout, Ho = H; end

%==========================================================================
function Ho = select_atlas(atlas, min_extent, min_overlap)
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
elseif atlas == 3
  atlas_name = 'HCP Multi-Modal Parcellation';
end

% go through left and right hemisphere
for ind = [1 3]
  
  % atlas data
  rcsv = H.rcsv{atlas};
  rdata = H.rdata{atlas}(:, round(ind / 2));
  
  A = H.S{round(ind / 2)}.A;
  A = A + speye(size(A));
  d0 = getappdata(H.patch(ind), 'data');
  
  % go through all surfaces
  for indsurf = 1:H.n_surf
    d = d0(indsurf, :);
    
    % apply thresholds
    dp = d > thresh(2); indp = find(dp);
    dn = d < thresh(1); indn = find(dn);
    
    % go through pos. effects
    if ~isempty(indp)
      
      C = find_connected_component(A, dp);
      C = C(indp);
      rdata2 = rdata(indp);
      
      fprintf('\n\n______________________________________________________\n');
      fprintf('%s: Positive effects in %s', atlas_name, H.S{round(ind / 2)}.info(1).side);
      fprintf('\n%s', spm_str_manip(H.S{round(ind / 2)}.info(indsurf).fname, 'k50d'));
      fprintf('\n______________________________________________________\n\n');
      
      if H.logP(H.results_sel), fprintf('%7s\t%16s\t%s\n', 'P-value', 'Cluster-Size', 'Overlap of atlas region');
      else, fprintf('%7s\t%16s\t%s\n', 'Value  ', 'Cluster-Size', 'Overlap of atlas region'); end
      
      for i = 1:max(C)
        N = find(C == i);
        k = length(N);
        
        dmax = d(indp); dmax = max(dmax(N));
        
        if k >= min_extent
          if H.logP(H.results_sel), fprintf('\n%1.5f\t%16d', 10^(-dmax), k);
          else, fprintf('\n%.5f\t%16d', dmax, k); end
          
          Nrdata = rdata2(N);
          roi_size = zeros(size(rcsv, 1) - 1, 1);
          
          for j = 2:size(rcsv, 1)
            ind3 = find(Nrdata == rcsv{j, 1});
            roi_size(j - 1) = 100 * length(ind3) / k;
          end
          
          % sort wrt size
          [ii, jj] = sort(roi_size, 'descend');
          jj(ii == 0) = [];
          
          for j = 1:length(jj)
            if roi_size(jj(j)) >= min_overlap
              if j == 1, fprintf('\t%3.0f%s\t%s\n', roi_size(jj(j)), '%', rcsv{jj(j) + 1, 2});
              else, fprintf('%7s\t%16s\t%3.0f%s\t%s\n', '       ', '        ', ...
                  roi_size(jj(j)), '%', rcsv{jj(j) + 1, 2});
              end
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
      fprintf('%s: Negative effects in %s', atlas_name, H.S{round(ind / 2)}.info(1).side);
      fprintf('\n%s', spm_str_manip(H.S{round(ind / 2)}.info(indsurf).fname, 'k50d'));
      fprintf('\n______________________________________________________\n\n');
      
      if H.logP(H.results_sel), fprintf('%7s\t%16s\t%s\n', 'P-value', 'Cluster-Size', 'Overlap of atlas region');
      else, fprintf('%7s\t%16s\t%s\n', 'Value  ', 'Cluster-Size', 'Overlap of atlas region'); end
      
      for i = 1:max(C)
        N = find(C == i);
        k = length(N);
        
        dmin = d(indn); dmin = min(dmin(N));

        if k >= min_extent
          if H.logP(H.results_sel), fprintf('\n%1.5f\t%16d', 10^(dmin), k);
          else, fprintf('\n%.5f\t%16d', -dmin, k); end
          
          Nrdata = rdata2(N);
          roi_size = zeros(size(rcsv, 1) - 1, 1);
          for j = 2:size(rcsv, 1)
            ind3 = find(Nrdata == rcsv{j, 1});
            roi_size(j - 1) = 100 * length(ind3) / k;
          end
          
          % sort wrt size
          [ii, jj] = sort(roi_size, 'descend');
          jj(ii == 0) = [];
          
          for j = 1:length(jj)
            if roi_size(jj(j)) >= min_overlap
              if j == 1, fprintf('\t%3.0f%s\t%s\n', roi_size(jj(j)), '%', rcsv{jj(j) + 1, 2});
              else, fprintf('%7s\t%16s\t%3.0f%s\t%s\n', '       ', '        ', ...
                  roi_size(jj(j)), '%', rcsv{jj(j) + 1, 2});
              end
            end
          end
        end
      end
    end
  end
end
if nargout, Ho = H; end

%==========================================================================
function Ho = select_results(sel)
global H

H.results_sel = sel;

clearDataCursorPlot(H);

H.S{1}.name = H.S1.name(sel, :);
H.S{2}.name = H.S2.name(sel, :);
H.S{1}.Y = H.S1.Y(:, sel);
H.S{2}.Y = H.S2.Y(:, sel);

H.S{1}.info = cat_surf_info(H.S{1}.name, 1);          
H.S{2}.info = H.S{1}.info;            
H.S{1}.info(1).side = 'lh';
H.S{2}.info(1).side = 'rh';

% check whether data for left or right hemipshere are all non-zero
ind1 = find(H.S{1}.Y(:) ~= 0);
ind2 = find(H.S{2}.Y(:) ~= 0);

% estimate min value > 0 and min/max values
if any(ind1) && any(ind2)
  H.S{1}.thresh = min(abs(H.S{1}.Y(abs(H.S{1}.Y(:)) > 0)));
  H.S{1}.thresh = min(H.S{1}.thresh, min(abs(H.S{2}.Y(abs(H.S{2}.Y(:)) > 0))));
  H.S{1}.min = min(min(H.S{1}.Y(~isinf(H.S{1}.Y))), min(H.S{2}.Y(~isinf(H.S{2}.Y))));
  H.S{1}.max = max(max(H.S{1}.Y(~isinf(H.S{1}.Y))), max(H.S{2}.Y(~isinf(H.S{2}.Y))));
elseif isempty(ind1)
  H.S{1}.thresh = min(abs(H.S{2}.Y(abs(H.S{2}.Y(:)) > 0)));
  H.S{1}.min = min(H.S{2}.Y(~isinf(H.S{2}.Y)));
  H.S{1}.max = max(H.S{2}.Y(~isinf(H.S{2}.Y)));
elseif isempty(ind2)
  H.S{1}.thresh = min(abs(H.S{1}.Y(abs(H.S{1}.Y(:)) > 0)));
  H.S{1}.min = min(H.S{1}.Y(~isinf(H.S{1}.Y)));
  H.S{1}.max = max(H.S{1}.Y(~isinf(H.S{1}.Y)));
end

mn = H.S{1}.min;

% deal with neg. values
if H.S{1}.min < 0
  mnx = max(abs([H.S{1}.min, H.S{1}.max]));
  H.S{1}.min = - mnx;
  H.S{1}.max = mnx;
end

% add 10% to min/max values
H.S{1}.max = round(1100 * H.S{1}.max) / 1000;
if H.S{1}.min < 0
  H.S{1}.min = round(1100 * H.S{1}.min) / 1000;
else
  H.S{1}.min = round(900 * H.S{1}.min) / 1000;
end

% correct lower clim to "0" if no values are exceeding threshold
if mn > -H.thresh_value
  H.clim = [true H.thresh_value H.S{1}.max];
else
  H.clim = [true H.S{1}.min H.S{1}.max];
end

if H.logP(H.results_sel)
  H.clim(2) = floor(H.clim(2));
  H.clim(3) = ceil(H.clim(3));
end

% only apply thresholds that are slightly larger than zero
if H.S{1}.thresh > 0.00015 && H.thresh_value == 0
  H.clip = [true -H.S{1}.thresh H.S{1}.thresh];
else
  H.clip = [true -H.thresh_value H.thresh_value];
end

% rather use NaN values for zero threshold
if H.thresh == 0
  H.clip = [false NaN NaN];
end

if H.no_neg
  H.clip = [true -Inf H.clip(3)];
  set(H.slider_min, 'Value', 0);
  set(H.str_min, 'String', 0);
  set(H.slider_min, 'Min', 0);
  set(H.slider_max, 'Min', 0);
end

H.n_surf = 1;

for ind = 1:5
  if H.S{1}.thresh > 0.00015
    setappdata(H.patch(ind), 'clip', H.clip);
  end
  
  % update clim only for non-fixed scaling
  if ~H.fixscl
    setappdata(H.patch(ind), 'clim', H.clim);
  end
  col = getappdata(H.patch(ind), 'col');
  
  if ind > 4
    d = [H.S{1}.Y; H.S{2}.Y];
  else
    d = H.S{round(ind / 2)}.Y;
  end
  
  H = updateTexture(H, ind, d, col, H.show_transp);
end

% correct value of slider if no values are exceeding threshold
if H.S{1}.min > - H.thresh_value
  set(H.slider_min, 'Value', 0);
  set(H.str_min, 'String', 0);
  set(H.slider_min, 'Min', 0);
  set(H.slider_max, 'Min', 0);
end

% update sliders for non-fixed scaling
if ~H.fixscl
  if H.logP(H.results_sel)
    H.clim(2) = floor(H.clim(2));
    H.clim(3) = ceil(H.clim(3));
  end
  
  set(H.slider_min, 'Value', H.clim(2));
  set(H.slider_max, 'Value', H.clim(3));
  if H.clim(2) > 0
    set(H.slider_min, 'Min', ceil(0.5*H.clim(2)), 'Max', ceil(2*H.clim(3)));
    set(H.slider_max, 'Min', ceil(0.5*H.clim(2)), 'Max', ceil(2*H.clim(3)));
  else
    set(H.slider_min, 'Min', ceil(2*H.clim(2)), 'Max', ceil(2*H.clim(3)));
    set(H.slider_max, 'Min', ceil(2*H.clim(2)), 'Max', ceil(2*H.clim(3)));
  end
  set(H.str_min, 'String', sprintf('%g',H.clim(2)));
  set(H.str_max, 'String', sprintf('%g',H.clim(3)));
end

% update file information and colorbar
checkbox_info;

% only show threshold popup if log-name was found and minimal value > 0 is
% < -log10(0.05)
if H.logP(H.results_sel) && (H.S{1}.thresh < -log10(0.05))
  set(H.thresh, 'Enable', 'on');
else
  set(H.thresh, 'Enable', 'off');
end

if min(min(H.S{1}.Y(:)), min(H.S{2}.Y(:))) < 0 && H.n_surf == 1
  set(H.hide_neg, 'Enable', 'on');
  set(H.hide_neg, 'Value', 0);
else
  set(H.hide_neg, 'Enable', 'off');
  set(H.hide_neg, 'Value', 0);
end

% enable/disable atlas widget
if H.isvol(sel)
  H.Pvol_sel = H.Pvol{sel};
  if isfield(H,'Pvol_sel')
    H = update_slice_overlay(H);
  end
  set(H.atlas, 'Enable', 'off');
else
  set(H.atlas, 'Enable', 'on');
end

if ~H.disable_cbar
  H = show_colorbar(H);
end

% Don't allow plot functions for volume data
if H.isvol(sel)
  str = {'Data Cursor', 'Disable data cursor', 'Atlas regions: All Atlases', ...
    'Atlas regions: Desikan-Killiany DK40', 'Atlas regions: Destrieux 2009', ...
    'Atlas region: HCP Multi-Modal Parcellation', 'Enable/Disable rotate3d'};
  tmp = {{@select_cursor, 0}, ...
         {@select_cursor, 1}, ...
         {@select_cursor, 2}, ...
         {@select_cursor, 3}, ...
         {@select_cursor, 4}, ...
         {@select_cursor, 7}};
else
  str = {'Data Cursor', 'Disable data cursor', 'Atlas regions: All Atlases', ...
    'Atlas regions: Desikan-Killiany DK40', 'Atlas regions: Destrieux 2009', ...
    'Atlas region: HCP Multi-Modal Parcellation', 'Plot data at vertex', ...
    'Plot mean data inside cluster', 'Enable/Disable rotate3d'};
  tmp = {{@select_cursor, 0}, ...
         {@select_cursor, 1}, ...
         {@select_cursor, 2}, ...
         {@select_cursor, 3}, ...
         {@select_cursor, 4}, ...
         {@select_cursor, 5}, ...
         {@select_cursor, 6}, ...
         {@select_cursor, 7}};
end

set(H.cursor,'String', str, 'UserData', tmp);

if H.results_sel < numel(H.Pvol)
  set(H.next,'Visible','on')
else
  set(H.next,'Visible','off')
end

if H.results_sel > 1
  set(H.prev,'Visible','on')
else
  set(H.prev,'Visible','off')
end

% print selected filename
cla(H.nam);
axis(H.nam, 'off')
text(0.5, 0.5, spm_str_manip(H.S{1}.name, 'k60d'), 'Parent', H.nam, 'Interpreter', 'none', ...
  'FontSize', H.FS, 'HorizontalAlignment', 'center');
if nargout, Ho = H; end

%==========================================================================
function H = update_slice_overlay(H,file_save)
global OV

if ~isfield(H,'Pvol_sel') || isempty(H.Pvol{1}), return; end

if ~H.isvol
  mn = min(min(H.S{1}.Y(:)), min(H.S{2}.Y(:)));
  mx = max(max(H.S{1}.Y(:)), max(H.S{2}.Y(:)));
else
  mn = H.Pvol_mn{H.results_sel};
  mx = H.Pvol_mx{H.results_sel};
end

% check whether values are above threshold and close windows and return if not
if (mn < 0 && mn > -H.thresh_value) && (mx >= 0 && mx < H.thresh_value)
  try close(OV.fig); end
  try close(OV.fig_mip); end
  return
end

% correct position of overlay window
OV.pos = [10 H.SS(4)];

% display image as overlay
OV.reference_image = char(cat_get_defaults('extopts.shootingT1'));

OV.FS = 0.12;
OV.bkg_col = H.bkg_col;

if H.show_transp
  OV.opacity = 0.6;                                      
  OV.reference_range = [0.2 1.6];                        
else
  OV.opacity = Inf;                                      
  OV.reference_range = [0.2 1.0];                        
end

% clipping if (meaningful) defined
if ~isempty(H.clip) && H.clip(2) ~= H.clip(3)
  if ~isnan(H.clip(2)) && ~isnan(H.clip(3))
    OV.clip = H.clip(2:3);
  else
    OV.clip = [Inf -Inf];
  end
  OV.func = sprintf('i1(i1>%f & i1<%f)=NaN;',OV.clip(1),OV.clip(2));
else
  OV.func = 'i1(i1==0)=NaN;';
end

if isfield(H,'cmap_col')
  if size(H.cmap_col,2) == 3
    OV.cmap = H.cmap_col;                                      
  else
    OV.cmap = jet;
  end
else
  OV.cmap = jet;
end
OV.overview = []; % don't show slice overviev

OV.range = H.clim(2:3);

OV.name = char(H.Pvol{H.results_sel});
if nargin < 2
  OV.save = 'none';
else
  OV.save = file_save;
  [pth,name,ext] = fileparts(file_save);
  OV.save_image = fullfile(pth,[strrep(name,'slices_','mip_') ext]);
end

if H.disable_cbar
  OV.cbar = [];
else
  OV.cbar = 2;
end

% show MIP and keep position if window exists
OV.fig_mip = 12;
if ishandle(OV.fig_mip)
  cat_vol_img2mip(OV);
else
  cat_vol_img2mip(OV);
  pos = get(gcf,'Position');
  set(gcf,'Position',[H.SS(3) - pos(3) 0 pos(3:4)]);
end

% default OV fields
OV.atlas = 'none';
OV.slices_str = '-55:5:60';
OV.xy = [5 5];
OV.transform = char('axial');

if isfield(H,'OV')
  ov_field = {'atlas','slices_str','xy','transform'};
  for i=1:numel(ov_field)
    if isfield(H.OV,ov_field{i})
      OV.(ov_field{i}) = H.OV.(ov_field{i});
    end
  end
end

% show MIP and keep position if window exists
OV.fig = 22;

OV.func = [];
if ishandle(OV.fig)
  try, cat_vol_slice_overlay(OV); end
else
  OV.pos(1:2) = [10 1200];
  try, cat_vol_slice_overlay(OV); end
end


% disable file saving now
OV.save = 'none';
OV.save_image = '';

% select numbers of columns/rows
OV.xy_sel = get_xy(OV);
str = 'Slice Columns/Rows';
H.OV_xy = uicontrol(H.panel(2), ...
  'String', cellstr(char(str,num2str(OV.xy_sel))), 'Units', 'normalized', ...
  'Position', H.pos{2}.xy, 'Callback', @select_OV_xy, ...
  'Style', 'PopUpMenu', 'HorizontalAlignment', 'center', ...
  'FontSize',H.FS,...
  'ToolTipString', 'Select columns and rows', ...
  'Interruptible', 'on', 'Enable', 'on');

% select orientation
H.OV_trans = uicontrol(H.panel(2), ...
  'String', {'Slice Orientation','axial','sagittal','coronal'}, 'Units', 'normalized', ...
  'Position', H.pos{2}.trans, 'Callback', @select_OV_trans, ...
  'Style', 'PopUpMenu', 'HorizontalAlignment', 'center', ...
  'FontSize',H.FS,...
  'ToolTipString', 'Select orientation', ...
  'Interruptible', 'on', 'Enable', 'on');

% define slices
uicontrol(H.panel(2),'Style', 'Text', 'Units', 'normalized','FontSize',H.FS,...
  'BackgroundColor',H.col(1,:),...
  'String', 'Slices', 'Position', H.pos{2}.str3);
H.OV_slice = uicontrol(H.panel(2),'Style', 'edit', 'Units', 'normalized',...
  'ToolTipString', 'Define slices. Use Inf, NaN or empty string for automatically estimating slices with local maxima', ...
  'String', OV.slices_str, 'Position', H.pos{2}.slice, 'Callback', @OV_slice);

H.OV_labels = uicontrol(H.panel(2), ...
  'String', 'Hide slice labels', 'Units', 'normalized', ...
  'BackgroundColor',H.col(1,:),...
  'Position', H.pos{2}.labels, ...
  'Style', 'CheckBox', 'HorizontalAlignment', 'center', ...
  'Callback', {@checkbox_OV_labels}, ...
  'FontSize',H.FS,...
  'ToolTipString', 'Hide slice labels', ...
  'Interruptible', 'on', 'Enable', 'on');
    
H.OV_save = uicontrol(H.panel(2), ...
  'String', 'Save png', 'Units', 'normalized', ...
  'Position', H.pos{2}.save, 'Callback', @save_image, ...
  'Style', 'PushButton', 'HorizontalAlignment', 'center', ...
  'FontSize',H.FS,...
  'ToolTipString', 'Save overlay as png files', ...
  'Interruptible', 'on', 'Enable', 'on');

%==========================================================================
function checkbox_OV_labels(obj, event_obj)
global H OV

hide_labels = get(H.OV_labels, 'Value');
if hide_labels
  OV.labels = [];
else
  OV.labels.size = H.FS;
end
cat_vol_slice_overlay(OV);

%==========================================================================
function xy = get_xy(OV)

% get number of images with optional colorbar
n = numel(str2num(OV.slices_str));
if ~isempty(OV.cbar), n = n + 1; end

nn = round(n^0.4);
if n > 8, x = nn:round(n / nn); else x = 1:n; end
xy = [];
for i = 1:length(x)
    y = round(n / x(i));
  % check whether y is to small
  while y * x(i) < n, y = y + 1; end
  if i > 2
    if y * x(i - 1) < n, xy = [xy; [x(i) y]]; end
  else xy = [xy; [x(i) y]]; end
end

% change order of x and y
for i = 1:size(xy, 2)
  yx = [xy(i, 2) xy(i, 1)];
  xy = [xy; yx];
end

% remove duplicates
xy = [[n 1];xy];
xy = unique(xy, 'rows');

%==========================================================================
function select_OV_xy(hObject, event)
global OV H

value = get(hObject, 'value');
if value > 1
  OV.xy_sel = get_xy(OV);
  set(H.OV_xy,'String', cellstr(char('Slice Columns/Rows',num2str(OV.xy_sel))));
  OV.xy = OV.xy_sel(value-1,:);
  cat_vol_slice_overlay(OV);
end

%==========================================================================
function select_OV_trans(hObject, event)
global OV

value = get(hObject, 'value');
if value > 1
  str = get(hObject,'String');
  OV.transform = str{value};
  cat_vol_slice_overlay(OV);
end

%==========================================================================
function OV_slice(hObject, event)
global OV H

slices_str = get(hObject, 'String');
if ~isempty(str2num(slices_str)) || strcmp(slices_str,'') || isempty(slices_str)
  OV.slices_str = slices_str;
  if ~isfinite(str2num(slices_str))
    OV.slices_str = '';
  end
  OV = cat_vol_slice_overlay(OV);
  OV.xy_sel = get_xy(OV);
  set(H.OV_xy,'String', cellstr(char('Slice Columns/Rows',num2str(OV.xy_sel))));
else
  fprintf('You have to define numbers\n');
end

%==========================================================================
function Ho = select_surf(surf)
global H

H.surf_sel = surf;

for ind = 1:2
  switch surf
    case 1
      H.S{ind}.info(1).Pmesh = fullfile(fileparts(mfilename('fullpath')), ['templates_surfaces' H.str32k], ...
      [H.S{ind}.info(1).side '.central.freesurfer.gii']);
    case 2
      H.S{ind}.info(1).Pmesh = fullfile(fileparts(mfilename('fullpath')), ['templates_surfaces' H.str32k], ...
      [H.S{ind}.info(1).side '.inflated.freesurfer.gii']);
    case 3
      H.S{ind}.info(1).Pmesh = fullfile(fileparts(mfilename('fullpath')), ['templates_surfaces' H.str32k], ...
      [H.S{ind}.info(1).side '.central.' cat_get_defaults('extopts.shootingsurf') '.gii']);
    case 4
      H.S{ind}.info(1).Pmesh = fullfile(fileparts(mfilename('fullpath')), ['templates_surfaces' H.str32k], ...
      [H.S{ind}.info(1).side '.patch.freesurfer.gii']);
  end
  H.S{ind}.M = gifti(H.S{ind}.info(1).Pmesh);
end

for ind = 1:5
  if ind < 5 % single hemisphere views
    M = H.S{round(ind / 2)}.M;
  else
    M.faces = [H.S{1}.M.faces; H.S{2}.M.faces + size(H.S{1}.M.vertices, 1)];
    M.vertices = [H.S{1}.M.vertices; H.S{2}.M.vertices];
    M.mat = H.S{1}.M.mat;
  end
  
  set(H.patch(ind), 'Vertices', M.vertices);
  set(H.patch(ind), 'Faces', M.faces);
  
  % rescale axis except for flatmaps
  Ha = getappdata(H.patch(ind), 'axis');
  axes(Ha);
  
  if surf < 4
    axis(Ha, 'image');
    axis(Ha, 'off');
  end
end

% only show lateral views for flatmaps
if surf == 4
  select_view(3)
elseif H.view == 3
  select_view(1)
end

% don't show data cursor, view functions and data plot that will not work for flatmaps
if surf == 4
  set(H.cursor, 'Enable', 'off');
  set(H.mview, 'Enable', 'off');
  clearDataCursorPlot(H);
else
  set(H.cursor, 'Enable', 'on');
  set(H.mview, 'Enable', 'on');
end

if nargout, Ho = H; end

%==========================================================================
function display_results_all(obj, event_obj)
global H

if (size(H.S{1}.Y) > 1 | size(H.S{2}.Y) > 1) & min(min(H.S{1}.Y(:)), min(H.S{2}.Y(:))) < 0
  disp('Warning: Only results with positive values are displayed!');
end

% clear larger area and set background color to update labels and title
H.Ha = axes('Parent', H.panel(1), 'Position', [-.1 -.1 1.1 1.1], 'Color', H.bkg_col);
cla(H.Ha);

H.renderer = get(H.figure, 'Renderer');
set(H.figure, 'Renderer', 'OpenGL');

%-Get mesh curvature and sulcal depth
%------------------------------------------------------------------
for i = 1:2
  g1 = gifti(fullfile(fileparts(mfilename('fullpath')), ['templates_surfaces' H.str32k], [H.S{i}.info(1).side '.mc.freesurfer.gii']));
  g2 = gifti(fullfile(fileparts(mfilename('fullpath')), ['templates_surfaces' H.str32k], [H.S{i}.info(1).side '.sqrtsulc.freesurfer.gii']));
  H.S{i}.curv = cell(3, 1);
  H.S{i}.curv{1} = g1.cdata;
  H.S{i}.curv{2} = g2.cdata;
  H.S{i}.curv{3} = zeros(size(g2.cdata));
end

if H.view == 1 % top view
  vv = [90 0; -90 0; -90 0; 90 0; 0 90];
else % bottom view
  vv = [90 0; -90 0; -90 0; 90 0; 0 -90];
end

for ind = 1:5
  display_results(ind, H.viewpos{ind}(abs(H.view), :), vv(ind, :));
end

% prepare dataplot axes
H.dataplot = axes('Position', H.viewpos{6}(abs(H.view), :), 'Parent', H.panel(1), 'Color', H.bkg_col);
H.figure = ancestor(H.dataplot, 'figure');
try axes(H.dataplot); end
axis off

% check whether data for left or right hemipshere are all non-zero
ind1 = H.S{1}.Y(:) ~= 0;
ind2 = H.S{2}.Y(:) ~= 0;

% estimate min value > 0 and min/max values
if any(ind1) && any(ind2)
  H.S{1}.thresh = min(abs(H.S{1}.Y(abs(H.S{1}.Y(:)) > 0)));
  H.S{1}.thresh = min(H.S{1}.thresh, min(abs(H.S{2}.Y(abs(H.S{2}.Y(:)) > 0))));
  H.S{1}.min = min(min(H.S{1}.Y(~isinf(H.S{1}.Y(:)) & ind1)), min(H.S{2}.Y(~isinf(H.S{2}.Y(:)) & ind2)));
  H.S{1}.max = max(max(H.S{1}.Y(~isinf(H.S{1}.Y(:)) & ind1)), max(H.S{2}.Y(~isinf(H.S{2}.Y(:)) & ind2)));
elseif all(ind1==0)
  H.S{1}.thresh = min(abs(H.S{2}.Y(abs(H.S{2}.Y(:)) > 0)));
  H.S{1}.min = min(H.S{2}.Y(~isinf(H.S{2}.Y(:)) & ind2));
  H.S{1}.max = max(H.S{2}.Y(~isinf(H.S{2}.Y(:)) & ind2));
elseif all(ind2==0)
  H.S{1}.thresh = min(abs(H.S{1}.Y(abs(H.S{1}.Y(:)) > 0)));
  H.S{1}.min = min(H.S{1}.Y(~isinf(H.S{1}.Y(:)) & ind1));
  H.S{1}.max = max(H.S{1}.Y(~isinf(H.S{1}.Y(:)) & ind1));
end

% deal with neg. values
if H.S{1}.min < 0
  mnx = max(abs([H.S{1}.min, H.S{1}.max]));
  H.S{1}.min = - mnx;
  H.S{1}.max = mnx;
end

if H.logP(H.results_sel)
  mnS = floor(H.S{1}.min);
  mxS = ceil(H.S{1}.max);
else
  mnS = H.S{1}.min;
  mxS = H.S{1}.max;
end

H.clim = [true mnS mxS];

% only apply thresholds that are slightly larger than zero
if H.S{1}.thresh > 0.00015
  H.clip = [true -H.S{1}.thresh H.S{1}.thresh];
end

for ind = 1:5
  if H.S{1}.thresh > 0.00015
    setappdata(H.patch(ind), 'clip', H.clip);
  end
  setappdata(H.patch(ind), 'clim', [true mnS mxS]);
  col = getappdata(H.patch(ind), 'col');
  d = getappdata(H.patch(ind), 'data');
  H = updateTexture(H, ind, d, col, H.show_transp);
end

if isfield(H,'Pvol_sel')
  H = update_slice_overlay(H);
end

% only show threshold popup if log-name was found and minimal value > 0 is < -log10(0.05)
if H.logP(H.results_sel) && (H.S{1}.thresh < -log10(0.05))
  set(H.thresh, 'Enable', 'on');
end

if min(min(H.S{1}.Y(:)), min(H.S{2}.Y(:))) < 0 && H.n_surf == 1
  set(H.hide_neg, 'Enable', 'on');
  set(H.hide_neg, 'Value', 0);
end

if H.n_surf == 1 && ~H.isvol(H.results_sel)
  % get sure that image is thresholded and there are at least 20% zero/NaN areas
  if (sum(d ~= 0) / numel(d) < 0.8)
    set(H.atlas, 'Enable', 'on');
  end
end

if ~H.disable_cbar
  H = show_colorbar(H);
end

% show slider for range of results
if H.n_surf == 1
  
  % allow slider a more extended range
  max_abs = ceil(2 * max(abs([H.S{1}.min H.S{1}.max])));
  if H.S{1}.min < 0
    mnx = [-max_abs max_abs];
  else
    mnx = [0 max_abs];
  end

  try
    [H.slider_min, H.slider{1}, H.str_min] = sliderPanel( ...
      'Parent', H.panel(2), ...
      'Title', 'Overlay min', ...
      'Position', H.pos{2}.ovmin, ...
      'Backgroundcolor', H.col(1,:), ...
      'Min', mnx(1), ...
      'Max', mnx(2), ...
      'Value', mnS, ...
      'FontName', 'Verdana', ...
      'FontSize', H.FS-1, ...
      'NumFormat', '%g', ...
      'Callback', @slider_clim_min);

    [H.slider_max, H.slider{2}, H.str_max] = sliderPanel( ...
      'Parent', H.panel(2), ...
      'Title', 'Overlay max', ...
      'Position', H.pos{2}.ovmax, ...
      'Backgroundcolor', H.col(1,:), ...
      'Min', mnx(1), ...
      'Max', mnx(2), ...
      'Value', mxS, ...
      'FontName', 'Verdana', ...
      'FontSize', H.FS-1, ...
      'NumFormat', '%g', ...
      'Callback', @slider_clim_max);
  end
end

%==========================================================================
function H = show_colorbar(H)

% show colorbar
if H.n_surf == 1
   
  if isfield(H, 'cbar')
    delete(findobj('tag','cat_surf_results_colorbar'));
    H = rmfield(H, 'cbar');
  end
  
  H.cbar = axes('Parent', H.panel(1), 'Position', H.pos{1}.cbar(1, :), 'Color', H.bkg_col, 'Visible', 'off','tag','cat_surf_results_colorbar');
  H.colourbar = colorbar('peer', H.cbar, 'Northoutside');
  
  if H.logP(H.results_sel), title(H.cbar, 'p-value', 'Color', 1 - H.bkg_col); end
  clim = getappdata(H.patch(1), 'clim');
  axis(H.cbar, 'off');
  
  col = H.cmap_col;
  figure(H.figure)
  colormap(col);

  if clim(3) > clim(2)
    caxis([clim(2) clim(3)]);
  end
  
  % Update colorbar colors if clipping is used
  clip = getappdata(H.patch(1), 'clip');
  if ~isempty(clip)
    if ~isnan(clip(2)) && ~isnan(clip(3))
      ncol = length(col);
      col_step = (clim(3) - clim(2)) / ncol;
      cmin = max([1, ceil((clip(2) - clim(2)) / col_step)]);
      cmax = min([ncol, floor((clip(3) - clim(2)) / col_step)]);
      col(cmin:cmax, :) = repmat([0.5 0.5 0.5], (cmax - cmin + 1), 1);
      colormap(col);
    end
  end
  
  if H.logP(H.results_sel)
    
    XTick = get(H.colourbar, 'XTick');
    
    % save original XTick values
    if isempty(H.XTick), H.XTick = XTick; end

    % if threshold is between 1.3..1.4 (p<0.05) change XTick at threshold
    % to (-)log10(0.05)
    if ~isempty(clip) && abs(clip(3)) >= 1.3 && abs(clip(3)) <= 1.4 && clim(3) > clim(2)
      XTick_step = ceil((clim(3) - clim(2)) / numel(XTick));
      if clip(2) <= - 1.3 && clip(2) >= - 1.4 && clim(2) < 0
        XTick = [round(clim(2)):XTick_step:round(clim(3))];
        mid = find(XTick==0);
        if ~isempty(mid)
          XTick(mid-1:mid+1) = [log10(0.05) 0 -log10(0.05)];
        else
          XTick(XTick == -1) =  log10(0.05);
          XTick(XTick == 1)  = -log10(0.05);
        end
      else
        XTick = [0:XTick_step:round(clim(3))];
        XTick(2) = -log10(0.05);
      end
      
    else
      mn = floor(min(XTick));
      mx = ceil(max(XTick));

      % only allow integer values
      XTick = floor(mn:mx);
    end
    
    % change XTickLabel
    XTickLabel = [];
    for i = 1:length(XTick)
      if XTick(i) > 0
        if XTick(i) > 7
          % use 1E-x notation
          XTickLabel = char(XTickLabel, remove_zeros(sprintf('%g', 10^(-XTick(i)))));
        else
          % use 0.000x notation
          XTickLabel = char(XTickLabel, remove_zeros(sprintf('%3.7f', 10^(-XTick(i)))));
        end
      elseif XTick(i) < 0
        if XTick(i) < -7
          % use 1E-x notation
          XTickLabel = char(XTickLabel, remove_zeros(sprintf('-%g', 10^(XTick(i)))));
        else
          % use 0.000x notation
          XTickLabel = char(XTickLabel, remove_zeros(sprintf('-%3.7f', 10^(XTick(i)))));
        end
      else
        XTickLabel = char(XTickLabel, '');
      end
    end
    
    set(H.colourbar, 'XTickLabel', XTickLabel(2:end, :), 'XTick', XTick);
  end % end H.logP
  
  set(H.colourbar, 'XColor', 1-H.bkg_col, 'YColor', 1-H.bkg_col);

  try set(H.colourbar, 'TickDirection','out'); end
  try
    set(H.colourbar, 'TickLength', 0.01);
  catch
    set(H.colourbar, 'TickLength', [0 0]);
  end
  
  %{
  if isfield(H,'hist')
    cat_surf_results('hist')
    cat_surf_results('hist')
  end
  %}
else
  delete(findobj('tag','cat_surf_results_hist'));
  
  if ~isfield(H, 'cbar') || ~ishandle(H.cbar)
    H.cbar = axes('Parent', H.panel(1), 'Position', H.pos{1}.cbar(2, :), 'Color', H.bkg_col);
    try
      set(H.cbar, 'Enable', 'off');
    end
  end
  
  % RGB colorbar
  if H.n_surf == 3
    cb = [8 1 1 4 2 2 8; ...
        8 1 6 7 5 2 8; ...
        8 8 3 3 3 8 8];
  else %RG colorbar
    cb = [8 1 1 4 2 2 8; ...
        8 1 1 4 2 2 8];
  end
  imagesc(cb);
  colormap([1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 1 1 1; H.bkg_col]);
  axis(H.cbar, 'off'); axis('image');
end

%==========================================================================
function display_results(ind, win, vw)
global H

% rescue old color before a new H.patch is created
try
  col = getappdata(H.patch(ind), 'col');
catch
  col = [];
end

col = H.cmap_col;

if ind < 5 % single hemisphere views
  M = H.S{round(ind / 2)}.M;
  Mc.cdata = H.S{round(ind / 2)}.Y;
else
  Ml = H.S{1}.M;
  Mr = H.S{2}.M;
  Mcl.cdata = H.S{1}.Y;
  Mcr.cdata = H.S{2}.Y;
  
  % check whether number of data for lh/rh differ and fill with zeros
  diff_size_Y = size(H.S{1}.Y, 2) - size(H.S{2}.Y, 2);
  if diff_size_Y > 0
    Mcr.cdata = [Mcr.cdata zeros(size(H.S{2}.Y, 1), 1)];
  end
  if diff_size_Y < 0
    Mcl.cdata = [Mcl.cdata; zeros(size(H.S{1}.Y, 1), 1)];
  end
  
  M.faces = [Ml.faces; Mr.faces + size(Ml.vertices, 1)];
  M.vertices = [Ml.vertices; Mr.vertices];
  M.mat = Ml.mat;
  Mc.cdata = [Mcl.cdata; Mcr.cdata];
end

if isfield(Mc, 'cdata')
  M.cdata = Mc.cdata;
else
  M.cdata = [];
end

H.axis = axes('Position', win, 'Parent', H.panel(1), 'Visible', 'off');
H.figure = ancestor(H.axis, 'figure');
%axes(H.axis);

if isfield(M, 'facevertexcdata')
  H.cdata = M.facevertexcdata;
else
  H.cdata = [];
end

if ~isfield(M, 'vertices') || ~isfield(M, 'faces')
  error('cat_surf_results:nomesh', 'ERROR:cat_surf_results: No input mesh.');
end

%% -Patch
%------------------------------------------------------------------
P = struct('vertices', M.vertices, 'faces', double(M.faces));
H.patch(ind) = patch(P, ...
  'FaceColor', [0.6 0.6 0.6], ...
  'EdgeColor', 'none', ...
  'FaceLighting', 'gouraud', ...
  'SpecularStrength', 0.1, ...
  'AmbientStrength', 1.0, ...
  'DiffuseStrength', 0.6, ...
  'SpecularExponent', 15, ...
  'Clipping', 'off', ...
  'DeleteFcn', {@myDeleteFcn, H.renderer}, ...
  'Visible', 'off', ...
  'Tag', 'CATSurfRender', ...
  'Parent', H.axis);
setappdata(H.patch(ind), 'patch', P);
setappdata(H.patch(ind), 'axis', H.axis);

%-Apply texture to mesh
%------------------------------------------------------------------
if isfield(M, 'facevertexcdata')
  T = M.facevertexcdata;
elseif isfield(M, 'cdata')
  T = M.cdata;
else
  T = [];
end

if isempty(col)
  H = updateTexture(H, ind, T);
else
  H = updateTexture(H, ind, T, col);
end

axis(H.axis, 'image');
axis(H.axis, 'off');
view(H.axis, vw);
material(H.figure, 'dull');

% default lighting
H.light(1) = camlight('headlight'); set(H.light(1), 'Parent', H.axis);
setappdata(H.axis, 'handles', H);
set(H.patch(ind), 'Visible', 'on');
camlight(H.light(1),'headlight')

%==========================================================================
function [H, C] = updateTexture(H, ind, v, col, transp)

%-Project data onto surface mesh
%--------------------------------------------------------------------------
if nargin<3
  v = getappdata(H.patch(ind), 'data');
end
if nargin<5
  transp = H.transp;
end
if size(v, 2) < size(v, 1)
  v = v';
end
v(isinf(v)) = NaN;


%-Get colourmap
%--------------------------------------------------------------------------
if ~exist('col', 'var')
  if size(v, 1) == 1
    col = jet(256);
  else
    % use RGB colormap
    col = zeros(256, 3, size(v, 1));
    for i = 1:3
      col(:, i, i) = 1;
    end
  end
end

setappdata(H.patch(ind), 'data', v);
setappdata(H.patch(ind), 'col', col);

if ~exist('FaceColor', 'var') || isempty(FaceColor), FaceColor = 'interp'; end

%-Get curvature
%--------------------------------------------------------------------------
if ind < 5 % single hemisphere views
  curv = H.S{round(ind / 2)}.curv{H.texture_mode};
else
  curv = [H.S{1}.curv{H.texture_mode}; H.S{2}.curv{H.texture_mode}];
end

if size(curv,2) == 1
  
  % emphasize mean curvature values by using sqrt
  indneg = find(curv < 0);
  curv(indneg) = - ((-curv(indneg)).^0.5);
  indpos = find(curv > 0);
  curv(indpos) = (curv(indpos).^0.5);
  curv = curv - min(curv);
  
  curv = 0.5 + repmat(curv, 1, 3);
  curv = curv / max(curv(:));
  
  % for sulcal depth (with no neg. values) use inverted values
  if H.texture_mode == 2
    curv = 1 - curv;
  end
end

%-Create RGB representation of data according to colourmap
%--------------------------------------------------------------------------
C = zeros(size(v, 2), 3);
clim = getappdata(H.patch(ind), 'clim');
if isempty(clim) || numel(clim)==1, clim = [false NaN NaN]; end
mi = clim(2); ma = clim(3);

if any(v(:))
  if size(v, 1) > 1 && size(col,3) == 1
    % use RGB colormap
    col = zeros(256, 3, size(v, 1));
    for i = 1:3
      col(:, i, i) = 1;
    end
  end
  
  if ~clim(1), mi = min(v(:)); ma = max(v(:)); end
  % don't allow negative values for multiple maps
  if size(v, 1) > 1 && mi < 0
    if ~isempty(H.clip)
      H.clip(2) = - Inf;
    else
      H.clip = [true -Inf 0];
    end
  end
  for i = 1:size(v, 1)
    C = C + squeeze(ind2rgb(floor(((v(i, :) - mi) / (ma - mi)) * (size(col, 1)-1)), col(:, :, i)));
  end
end

if ~isempty(H.clip)
  v(v > H.clip(2) & v < H.clip(3)) = NaN;
  setappdata(H.patch(ind), 'clip', [true H.clip(2) H.clip(3)]);
end

setappdata(H.patch(ind), 'clim', [true mi ma]);
H.clim = [true mi ma];

  
%-Build texture by merging curvature and data
%--------------------------------------------------------------------------
if size(v, 1) > 1 % RGB
  for i = 1:size(v, 1)
    C(:, i) = any(v(i, :), 1)' .* C(:, i);
  end
else
  C = repmat(any(v, 1), 3, 1)' .* C;
end

% add curvature pattern if transparency is defined
if nargin > 4
  if transp && size(C, 1) == size(curv, 1)
    C = (0.5 + 0.5 * curv) .* C;
  end
end

% replace regions below threshold by curvature
ind0 = repmat(~any(v, 1), 3, 1)';
if size(C, 1) == size(curv, 1)
  C(ind0) = curv(ind0);
else
  C(ind0) = 0;
end

if H.texture_mode == 3
    C(ind0) = 0.5;
end

%-Add/delete atlas border
%--------------------------------------------------------------------------

if ~H.border_mode
  for k=1:5
  try
    h3 = getappdata(H.patch(k), 'h3');
    if ~isempty(h3)
      for i=1:size(h3,1)
        delete(h3(i))
      end
    end
  end
  end
end


set(H.patch(ind), 'FaceVertexCData', C, 'FaceColor', FaceColor, 'FaceLighting','none');

if H.texture_mode == 3
  set(H.patch(ind), 'FaceLighting','none');
end

if H.border_mode
  if ind == 1
  for k=1:2
    rdata = H.rdata{H.border_mode}(:, k);
    datarange = unique(rdata);
    Hi = hist(rdata(:),datarange);
    datarange(Hi<2) = [];
  
    t = datarange;
    M = H.S{k}.M;
    H.S{k}.Cm = cell(numel(t),1);
    for i=1:numel(t)
      T = zeros(size(rdata));
      T(rdata == datarange(i)) = 1;
      try
        Cm = spm_mesh_contour(M,struct('T',T,'t',0.5));
        H.S{k}.Cm{i} = Cm;
      catch
        fprintf('Please update SPM12 for using that function.\n');
        break
      end
    end
  end
  for k=1:5
    h3 = getappdata(H.patch(k), 'h3');
    if ~isempty(h3)
    for i=1:size(h3,1)
      delete(h3(i))
    end
    end
    setappdata(H.patch(k), 'h3', []);
  end
  end
  
  Ha = getappdata(H.patch(ind), 'axis');
  hold on

  if ind < 5 % single hemisphere views
    Cm = H.S{round(ind / 2)}.Cm;
    col = jet(size(Cm,1));
    for j=1:size(Cm,1)
%      col_sel = col(j,:);
      col_sel = 'k';
      if ~isempty(Cm{j})
        for i=1:size(Cm,2)
          h3 = plot3(Cm{j}(i).xdata,Cm{j}(i).ydata,Cm{j}(i).zdata,'Color',col_sel,'LineWidth',2);
          setappdata(H.patch(ind), 'h3', [getappdata(H.patch(ind), 'h3'); h3]);
          set(h3,'Parent',Ha);
        end
      end
    end
  else
    for k = 1:2
      Cm = H.S{k}.Cm;
      col = jet(size(Cm,1));
      for j=1:size(Cm,1)
%        col_sel = col(j,:);
        col_sel = 'k';
        if ~isempty(Cm{j})
          for i=1:size(Cm,2)
            h3 = plot3(Cm{j}(i).xdata,Cm{j}(i).ydata,Cm{j}(i).zdata,'Color',col_sel,'LineWidth',2);
            setappdata(H.patch(ind), 'h3', [getappdata(H.patch(ind), 'h3'); h3]);
            set(h3,'Parent',Ha);
          end
        end
      end
    end
  end
  
  hold off

end

if 0 %isfield(H,'histax')
  cat_surf_results('hist')
  cat_surf_results('hist')
end
  
%==========================================================================
function next_file(obj, event_obj)
global H

if  H.results_sel + 1 <= numel(H.Pvol)
  H.results_sel = H.results_sel + 1;
else
  fprintf('Last file already selected.\n');
  return
end

H.Pvol_sel = deblank(H.Pvol{H.results_sel});
fprintf('New selected file: %s\n',H.S1.name(H.results_sel,:));
select_results(H.results_sel);

%==========================================================================
function prev_file(obj, event_obj)
global H

if  H.results_sel - 1 >= 1
  H.results_sel = H.results_sel - 1;
else
  fprintf('First file already selected.\n');
  return
end

H.Pvol_sel = deblank(H.Pvol{H.results_sel});
fprintf('New selected file: %s\n',H.S1.name(H.results_sel,:));
select_results(H.results_sel);
    
%==========================================================================
function select_data(obj, event_obj, P)
global H

if ~exist('P','var')
  P = spm_select([1 50], {'mesh','image'}, 'Select up to 50 volume or surface maps');
end

n = size(P, 1);
if n == 0; return; end

% correct filename and extension for volumes
for i = 1:n
  % volume found?
  if ~isempty(strfind(P(i,:),'.nii,'))
    P(i,:) = strrep(P(i,:),'.nii,1','.nii  ');
  end
end

H.isvol = zeros(n,1);
H.logP  = zeros(n,1);
H.Pvol  = cell(n,1);
P0      = cell(n,1);

% check for volumes
for i = 1:n
  if ~isempty(strfind(P(i,:),'.nii'))
    % names for template and thickness file and output 
    Pvol = deblank(P(i,:));
    [pp,ff,ee] = spm_fileparts(Pvol);
    Pmesh_lh  = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces_32k',['lh.central.' cat_get_defaults('extopts.shootingsurf') '.gii']);
    Pthick_lh = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces_32k',['lh.thickness.' cat_get_defaults('extopts.shootingsurf')]);
    Pout_lh   = fullfile(pp,['lh.',ff '.gii']);
    Pmesh_rh  = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces_32k',['rh.central.' cat_get_defaults('extopts.shootingsurf') '.gii']);
    Pthick_rh = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces_32k',['rh.thickness.' cat_get_defaults('extopts.shootingsurf')]);
    Pout_rh   = fullfile(pp,['rh.',ff '.gii']);
    Pout      = fullfile(pp,['mesh.',ff '.resampled_32k.gii']);
        
    % map 3D volume to template surface inside cortical band with maxabs mapping function
    cmd = sprintf('CAT_3dVol2Surf -maxabs -steps 7 -start -0.6 -end 0.6 -thickness "%s" "%s" "%s" "%s"',Pthick_lh, Pmesh_lh, Pvol, Pout_lh);
    cat_system(cmd,0);
    cmd = sprintf('CAT_3dVol2Surf -maxabs -steps 7 -start -0.6 -end 0.6 -thickness "%s" "%s" "%s" "%s"',Pthick_rh, Pmesh_rh, Pvol, Pout_rh);
    cat_system(cmd,0);

    % combine left and right hemipshere data
    M_lh = gifti(Pout_lh);
    M_rh = gifti(Pout_rh);
    delete(Pout_lh); delete(Pout_rh); % delete temporary files
    M.cdata = [M_lh.cdata; M_rh.cdata];
    M.private.metadata = struct('name','SurfaceID','value',Pout);
    save(gifti(M), Pout, 'Base64Binary');
    P0{i} = Pout; 
    H.isvol(i) = 1;
    H.Pvol{i} = Pvol;
    vol = spm_read_vols(spm_vol(Pvol));
    H.Pvol_mn{i} = min(vol(vol ~= 0));
    H.Pvol_mx{i} = max(vol(vol ~= 0));
    
    % print warning if no SPM.mat file was found
    if ~exist(fullfile(pp, 'SPM.mat'), 'file') && ~exist(fullfile(pp, 'vSPM.mat'), 'file')
      fprintf('Warning: Please ensure that %s is spatially registered to the new MNI152NLin2009cAsym template space.\n',[ff ee]);
    end
  else
    P0{i} = deblank(P(i,:));
  end
end
H.Pvol_sel = H.Pvol{1};

P0   = char(P0);
info = cat_surf_info(P0,0);

for i = 1:n
  if info(i).nvertices == 64984
    H.str32k = '_32k';
  else
    H.str32k = '';
  end
  
  % check whether name contains 'log' that indicates a logP file
  if isempty(strfind(info(i).ff, 'log'))
    H.logP(i) = 0;
  end
  
  if strcmp(info(i).side, 'lh') || strcmp(info(i).side, 'rh')
    spm('alert*','Display of separate hemispheres is not supported anymore');
    if nargout, varargout{1} = []; end
    return
  end

end

H.S{1}.name = P0;
H.S{2}.name = P0;

cat_surf_results('disp', P0);

if numel(H.Pvol) > 1
  set(H.next,'Visible','on')
end

%==========================================================================
function save_image(obj, event_obj, filename)

global H
%% 

try
  dcm_obj = datacursormode(H.figure);
  set(dcm_obj, 'Enable', 'off');
  delete(findall(gcf,'Type','hggroup'));
end

if ~exist('filename', 'var')
  
  nm = H.S{1}.info(1).ff;
  [pp, nm] = spm_fileparts(H.S{1}.name);
  pathname = spm_fileparts(pp);
  filename = [pathname '_' nm '.png'];
  
  % end with _0???.ext?
  if length(nm) > 4
    if strcmp(nm(length(nm) - 4:length(nm) - 3), '_0')
      
      SPM_name = fullfile(pp, 'SPM.mat');
      
      % SPM.mat exist?
      if exist(SPM_name, 'file')
        load(SPM_name);
        xCon = SPM.xCon;
        Ic = str2double(nm(length(nm) - 3:length(nm)));
        str_num = deblank(xCon(Ic).name);
        
        % replace spaces with "_" and characters like "<" or ">" with "gt" or "lt"
        str_num(strfind(str_num, ' ')) = '_';
        strpos = strfind(str_num, ' > ');
        if ~isempty(strpos), str_num = [str_num(1:strpos - 1) '_gt_' str_num(strpos + 1:end)]; end
        strpos = strfind(str_num, ' < ');
        if ~isempty(strpos), str_num = [str_num(1:strpos - 1) '_lt_' str_num(strpos + 1:end)]; end
        strpos = strfind(str_num, '>');
        if ~isempty(strpos), str_num = [str_num(1:strpos - 1) 'gt' str_num(strpos + 1:end)]; end
        strpos = strfind(str_num, '<');
        if ~isempty(strpos), str_num = [str_num(1:strpos - 1) 'lt' str_num(strpos + 1:end)]; end
        str_num = spm_str_manip(str_num, 'v');
        
        if ~isempty(H.clip)
          if isnan(H.clip(3))
            str_thresh = '_';
          else
            str_thresh = sprintf('P%g_', round(1000 * 10^(-H.clip(3))) / 10);
          end
        else
          str_thresh = '_';
        end
        filename = ['logP_' str_thresh str_num '.png'];
      end
    end
  end
  
  [filename, newpth] = uiputfile({ ...
    '*.png' 'PNG files (*.png)'}, 'Save as', filename);
else
  [pth, nam, ext] = fileparts(filename);
  if isempty(pth), pth = cd; end
  if ~strcmp({'.gii', '.png'}, ext), nam = [nam ext]; end
  if isempty(nam)
    [filename, newpth] = uiputfile({ ...
      '*.png' 'PNG files (*.png)'}, 'Save as', nam);
  else
    filename = [nam '.png'];
    newpth = pth;
  end
end

% keep background color
set(H.figure, 'InvertHardcopy', 'off', 'PaperPositionMode', 'auto');

% posc to correct margin of the figure [r t l b] if no dataplot is present
hh0 = getframe(H.figure);
sz0 = [0 0 size(hh0.cdata,[2 1])];
pos = round(sz0.*H.panel(1).Position);

if ~(isfield(H, 'dataplot') && strcmpi(get(H.dataplot,'Visible'),'on'))
  posc = [70 70 70 70]; pos = [pos(1) + posc(1), pos(2) + posc(2) + posc(4), pos(3) - posc(1) - posc(3), pos(4) - posc(2) - posc(4)];
end

img = frame2im(hh0);
img = img(pos(2)+1:pos(2)+pos(4),pos(1)+1:pos(1)+pos(3),:);
if H.results_sel ~= 4 && ~isfield(H, 'dataplot')
  % crop image if it's not a flatmap
  sz = size(img);
  img = img(round(0.1*sz(1):sz(1)),round(0.05*sz(2):0.95*sz(2)),:);
end

col = colormap;
imwrite(img,col,fullfile(newpth,filename));
fprintf('Image %s saved.\n',fullfile(newpth,filename));

H = update_slice_overlay(H,fullfile(newpth,['slices_' filename]));

% write dataplot window
if isfield(H, 'dataplot') && strcmpi(get(H.dataplot,'Visible'),'on')
  filename = ['plot_' filename];
  pos = round(getpixelposition(H.dataplot)); 
  
  % increase position to also include labels
  pos = round(pos.*[0.9 0.925 1.25 1.3]);
  hh = getframe(H.figure,pos);
  img_plot = frame2im(hh);
  imwrite(img_plot,col,fullfile(newpth,filename));
  fprintf('Dataplot %s saved.\n',fullfile(newpth,filename));
end

dcm_obj = datacursormode(H.figure);
set(dcm_obj, 'Enable', 'on');

%==========================================================================
function slider_clim_min(hObject, evt)
global H

val = get(hObject, 'Value');
c = getappdata(H.patch(1), 'clim');

% prevent range exceeding
if val > c(3)
  return
end

for ind = 1:5
  setappdata(H.patch(ind), 'clim', [true val c(3)]);
  col = getappdata(H.patch(ind), 'col');
  d = getappdata(H.patch(ind), 'data');
  H = updateTexture(H, ind, d, col, H.show_transp);
end

% update colorbar
if H.n_surf == 1 && ~H.disable_cbar
  H = show_colorbar(H);
end

H.clim = [true val c(3)];

if isfield(H,'Pvol_sel')
  H = update_slice_overlay(H);
end

%==========================================================================
function slider_clim_max(hObject, evt)
global H

val = get(hObject, 'Value');
c = getappdata(H.patch(1), 'clim');

% prevent range exceeding
if val < c(2)
  return
end

for ind = 1:5
  setappdata(H.patch(ind), 'clim', [true c(2) val]);
  col = getappdata(H.patch(ind), 'col');
  d = getappdata(H.patch(ind), 'data');
  H = updateTexture(H, ind, d, col, H.show_transp);
end

% update colorbar
if H.n_surf == 1 && ~H.disable_cbar
  H = show_colorbar(H);
end

H.clim = [true c(2) val];

if isfield(H,'Pvol_sel')
  H = update_slice_overlay(H);
end

%==========================================================================
function checkbox_inv(obj, event_obj)
global H

H.show_inv = get(H.inv, 'Value');

for ind = 1:5
  setappdata(H.patch(ind), 'clip', H.clip);
  col = getappdata(H.patch(ind), 'col');
  setappdata(H.patch(ind), 'col', flipud(col));
  d = getappdata(H.patch(ind), 'data');
  H = updateTexture(H, ind, d, flipud(col), H.show_transp);
end

H.cmap_col = flipud(col);

if isfield(H,'Pvol_sel')
  H = update_slice_overlay(H);
end

if ~H.disable_cbar
  H = show_colorbar(H);
end

%==========================================================================
function checkbox_fixscl(obj, event_obj)
global H

H.fixscl = get(H.fix, 'Value');

%==========================================================================
function Ho = checkbox_hide_neg(obj, event_obj)
global H

H.no_neg = get(H.hide_neg, 'Value');

thresh = H.thresh_value;
clip = getappdata(H.patch(1), 'clip');
clim = getappdata(H.patch(1), 'clim');

% get min value for both hemispheres
mn = min(min(min(getappdata(H.patch(1), 'data'))), min(min(getappdata(H.patch(3), 'data'))));

if H.no_neg
  H.clip = [true -Inf thresh];
  H.clim = [true thresh clim(3)];
  set(H.slider_min, 'Value', 0);
  set(H.str_min, 'String', 0);
  set(H.slider_min, 'Min', 0);
  set(H.slider_max, 'Min', 0);
else
  H.clip = [true -thresh thresh];
  if mn < -thresh
    H.clim = [true -clim(3) clim(3)];
    set(H.slider_min, 'Value', -clim(3));
    set(H.str_min, 'String', sprintf('%g',-clim(3)));
    set(H.slider_min, 'Min', ceil(-2*clim(3)));
    set(H.slider_max, 'Min', ceil(-2*clim(3)));
  end
end

for ind = 1:5
  setappdata(H.patch(ind), 'clip', H.clip);
  setappdata(H.patch(ind), 'clim', H.clim);
  col = getappdata(H.patch(ind), 'col');
  d = getappdata(H.patch(ind), 'data');
  mn = min(mn, min(d(:)));
  H = updateTexture(H, ind, d, col, H.show_transp);
end

% correct value of slider if no values are exceeding threshold
if mn > -thresh && H.n_surf == 1
  set(H.slider_min, 'Value', 0);
  set(H.str_min, 'String', 0);
  set(H.slider_min, 'Min', 0);
  set(H.slider_max, 'Min', 0);
end

if ~H.isvol(H.results_sel)
  set(H.atlas, 'Enable', 'on');
end

if ~H.disable_cbar
  H = show_colorbar(H);
end

if isfield(H,'Pvol_sel')
  H = update_slice_overlay(H);
end

if nargout, Ho = H; end

%==========================================================================
function checkbox_transp(obj, event_obj)
global H

H.show_transp = ~get(H.transp, 'Value');

for ind = 1:5
  col = getappdata(H.patch(ind), 'col');
  d = getappdata(H.patch(ind), 'data');
  H = updateTexture(H, ind, d, col, H.show_transp);
end

if isfield(H,'Pvol_sel')
  H = update_slice_overlay(H);
end

% update colorbar
if H.n_surf == 1 && ~H.disable_cbar
  H = show_colorbar(H);
end

%==========================================================================
function checkbox_bkg(obj, event_obj)
global H

H.white_bgk = get(H.bkg, 'Value');

if H.white_bgk
  H.bkg_col = [1 1 1];
else
  H.bkg_col = [0 0 0];
end

set(H.Ha, 'Color', H.bkg_col);
%set(get(H.cbar, 'Title'), 'Color', 1 - H.bkg_col);

title = findobj('tag','cat_surf_result_title'); 
if ~isempty(title)
  set( get( title ,'children'),'Color',  1 - H.bkg_col);
  %set(get(getappdata(H.patch(1), 'axis'), 'Title'), 'Color', 1 - H.bkg_col);
  %set(get(getappdata(H.patch(3), 'axis'), 'Title'), 'Color', 1 - H.bkg_col);
end

if H.n_surf == 1
  set(H.colourbar, 'XColor', 1 - H.bkg_col, 'YColor', 1 - H.bkg_col);
end

if ~H.disable_cbar
  H = show_colorbar(H);
end

if isfield(H, 'dataplot')
  try set(H.dataplot, 'XColor', 1 - H.bkg_col, 'YColor', 1 - H.bkg_col, 'Color', H.bkg_col); end
end

if isfield(H,'Pvol_sel')
  H = update_slice_overlay(H);
end

%==========================================================================
function checkbox_info(obj, event_obj)
global H

H.show_info = get(H.info, 'Value');

if H.show_info
  delete(findobj('tag','cat_surf_result_title'));

  ax = axes('Parent',H.panel(1),'Position',[0.5 0.0 0.9 0.03],'visible','off','tag','cat_surf_result_title','Color',H.bkg_col);  
  text(0,1,spm_str_manip(H.S{1}.name, 'k150d'),'HorizontalAlignment','center','interpreter','none','Color', 1 - H.bkg_col,'Parent',ax);
          
else
  delete(findobj('tag','cat_surf_result_title'));
  set(get(getappdata(H.patch(1), 'axis'), 'Title'), 'String', '')
  set(get(getappdata(H.patch(3), 'axis'), 'Title'), 'String', '')
end

%==========================================================================
function checkbox_nocbar(obj, event_obj)
global H OV

H.disable_cbar = get(H.nocbar, 'Value');

if H.disable_cbar
  % delete colorbar and title
  if H.n_surf == 1
    set(H.colourbar, 'Visible', 'off')
    set(get(H.cbar, 'Title'), 'Visible', 'off')
  else % delete only axis
    cla(H.cbar);
  end
else
  if H.n_surf == 1
    set(get(H.cbar, 'Title'), 'Visible', 'on')
    set(H.colourbar, 'Visible', 'on')
    H = show_colorbar(H);
  else
    H = show_colorbar(H);
  end
end

if isfield(H,'Pvol_sel')
  if H.disable_cbar
    OV.cbar = [];
  end
  H = update_slice_overlay(H);
end

%==========================================================================
function H = getHandles(H)

if ~nargin || isempty(H), H = gca; end
if ishandle(H) && ~isappdata(H, 'handles')
  a = H; clear H;
  H.axis = a;
  H.figure = ancestor(H.axis, 'figure');
  H.patch = findobj(H.axis, 'type', 'patch');
  H.light = findobj(H.axis, 'type', 'light');
  H.rotate3d = rotate3d(H.panel(1));
  setappdata(H.axis, 'handles', H);
else
  H = getappdata(H, 'handles');
end

%==========================================================================
function select_view(view)
global H

% check that view changed
if view ~= H.view
  
  if view > 0 % top view
    vv = [90 0; -90 0; -90 0; 90 0; 0 90];
  else % bottom view
    vv = [90 0; -90 0; -90 0; 90 0; 0 -90];
  end
  
  for ind = 1:5
    Ha = getappdata(H.patch(ind), 'axis');
    set(Ha, 'Position', H.viewpos{ind}(abs(view), :), 'View', vv(ind, :));
  end
  
  axes(Ha);
  camlight(H.light(1),'headlight')
  
  if isfield(H, 'dataplot')
    set(H.dataplot, 'Position', H.viewpos{6}(abs(view), :), 'Parent', H.panel(1), 'Color', H.bkg_col);
  end
  
  % save view
  H.view = view;
end

%==========================================================================
function select_texture(texture_mode)
global H

% check that view changed
if texture_mode ~= H.texture_mode
  
  H.texture_mode = texture_mode;
  
  for ind = 1:5
    col = getappdata(H.patch(ind), 'col');
    d = getappdata(H.patch(ind), 'data');
    H = updateTexture(H, ind, d, col, H.show_transp);
  end
  
end

%==========================================================================
function select_border(border_mode)
global H

H.border_mode = border_mode;

if H.border_mode
  set(H.surf, 'Enable', 'off');
  if (length(H.S{1}.Y) == 32492 || length(H.S{1}.Y) == 163842 || length(H.S{1}.Y) == 40962) && H.isfsavg
    fprintf('To change underlying surface, disable Atlas Border Overlay.\n');
  end
else
  if (length(H.S{1}.Y) == 32492 || length(H.S{1}.Y) == 163842 || length(H.S{1}.Y) == 40962) && H.isfsavg
    set(H.surf,   'Enable', 'on');
  end
end

for ind = 1:5
  col = getappdata(H.patch(ind), 'col');
  d = getappdata(H.patch(ind), 'data');
  H = updateTexture(H, ind, d, col, H.show_transp);
end

%==========================================================================
function select_cursor(cursor_mode)

global H

dcm_obj = datacursormode(H.figure);
H.cursor_mode = cursor_mode;

switch H.cursor_mode
  
  case 0 % disable and delete datatip
    rotate3d off;
    
    clearDataCursorPlot(H);
  case {1, 2, 3, 4}
    clearDataCursorPlot(H);
    set(dcm_obj, 'Enable', 'on', 'SnapToDataVertex', 'on', ...
      'DisplayStyle', 'datatip', 'Updatefcn', {@myDataCursorAtlas, H});
    % sometimes interpreter cannot be set for some Matlab versions
    try set(dcm_obj, 'Interpreter', 'none'); end
  case {5, 6}
    
    try
      delete(findall(gcf,'Type','hggroup'));
    end
    
    SPM_found = 1;
    SPM_name = fullfile(H.S{1}.info(1).pp, 'SPM.mat');
    
    % SPM.mat exist?
    if exist(SPM_name, 'file')
    
      load(SPM_name);

      % if analysis was moved we have to correct header structure
      SPM.VResMS = spm_data_hdr_read(fullfile(H.S{1}.info(1).pp,SPM.VResMS.fname));
      Vbeta = spm_data_hdr_read(fullfile(H.S{1}.info(1).pp,SPM.Vbeta(1).fname));
      for j=2:numel(SPM.Vbeta)
        Vbeta(j) = spm_data_hdr_read(fullfile(H.S{1}.info(1).pp,SPM.Vbeta(j).fname));
      end
      SPM.Vbeta = Vbeta;
      H.SPM{1} = SPM;

      str = 'predicted, adjusted or raw values?';
      H.predicted = spm_input(str, 1, 'm', {'predicted', 'adjusted', 'raw', 'adjusted raw'}, [1 0 -1 -2]);
      
      % ask for contrast for predicted or adjusted data
      if H.predicted >= 0
        H.Ic = spm_input('Which contrast?', 2, 'm', {SPM.xCon.name});
      end
    elseif ~isempty(H.S{1}.name)
      SPM_found = 0;
      spm('alert!', 'No SPM.mat file found. Please check that you have not moved your files or your result file was moved from the folder where the SPM.mat is stored.', 1);
    end
    
    if SPM_found && ~exist(SPM.xY.VY(1).fname,'file')
      SPM_found = 0;
      spm('alert!', 'No data found. Please check that you have not moved your files.', 1);
    end

    if SPM_found
      set(dcm_obj, 'Enable', 'on', 'SnapToDataVertex', 'on', ...
        'DisplayStyle', 'datatip', 'Updatefcn', {@myDataCursorCluster});
      % sometimes interpreter cannot be set for some Matlab versions
      try set(dcm_obj, 'Interpreter', 'none'); end
      if H.predicted > -2
        fprintf('The values are available at the MATLAB command line as variable ''y''\n');
      else
        fprintf('The adjusted raw values are available at the MATLAB command line as variable ''y''\n');
      end
    end
  case 7 % enable/disable rotate3d
    clearDataCursorPlot(H);
    rotate3d;
    fprintf('Use mouse to rotate views and right mouse click for more options.\n');
end

%==========================================================================
function clearDataCursorPlot(H)

if isfield(H, 'dataplot')
  cla(H.dataplot);
  axis(H.dataplot, 'off');
end

try
  dcm_obj = datacursormode(H.figure);
  set(dcm_obj, 'Enable', 'off');
  delete(findall(gcf,'Type','hggroup'));
end

%==========================================================================
function txt = myDataCursorCluster(obj, evt)
global H y x

% first entries are atlases
plot_mean = H.cursor_mode - 5;
pos = get(evt, 'Position');

i = ismember(get(H.patch(1), 'vertices'), pos, 'rows');
node = find(i);
ind = 1;
node_list = 1:numel(get(H.patch(1), 'vertices'));

if isempty(node)
  i = ismember(get(H.patch(3), 'vertices'), pos, 'rows');
  node = find(i);
  ind = 3;
  node_list = 1:numel(get(H.patch(3), 'vertices'));
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
  
  A = H.S{round(ind / 2)}.A;
  A = A + speye(size(A));
  d = getappdata(H.patch(ind), 'data');
  
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
  if ~isempty(indn) && isempty(found_node)
    
    C = find_connected_component(A, dn);
    C = C(indn);
    node_list2 = node_list(indn);
    
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
    cstr = 'neg. ';
  else
    cstr = '';
  end
  
  if isempty(found_node)
    txt = {'Cursor outside of cluster'};
  else
    if cluster_side == 1
      txt = {sprintf('lh: Cluster %s%d', cstr, cluster_number)};
    else
      txt = {sprintf('rh: Cluster %s%d', cstr, cluster_number)};
    end
  end
else
  % print value at cursor
  XYZ = node;
  value = H.S{round(ind / 2)}.Y(node);
  if H.logP
    if value < 0
      txt = {sprintf('p=-%g', 10^(value))};
    else
      txt = {sprintf('p=%g', 10^(-value))};
    end
  else
    txt = {sprintf('Value %g', value)};
  end
end

% for merged meshes we only have one SPM.mat with data from both hemispheres
% add offset for right hemisphere
if round(ind / 2) == 2
  XYZ = XYZ + H.nY2;
end

% always one mesh
ind = 1;
[y, cbeta, CI] = get_cluster_data(H, XYZ, ind);

% if no cluster was selected set data to zero
if plot_mean && isempty(found_node)
  y(:) = 0;
  cbeta(:) = 0;
  CI(:) = 0;
end

cla(H.dataplot)
hold(H.dataplot, 'on')

set(H.dataplot, 'XColor', 1 - H.bkg_col, 'YColor', 1 - H.bkg_col,...
    'YGrid','on','Visible','on');

xX = H.SPM{1}.xX;
iH = xX.iH;

Ic = [];
nm = H.S{1}.info(1).ff;
% end with _0???.ext?
if length(nm) > 4
  if strcmp(nm(length(nm) - 4:length(nm) - 3), '_0')
    Ic = str2double(nm(length(nm) - 3:length(nm)));
  elseif isfield(H,'Ic')
    Ic = H.Ic;
  end
else
  if isfield(H,'Ic')
    Ic = H.Ic;
  end
end

if H.predicted >=0
  ystr = 'contrast estimate';
  groupcolor = jet(size(cbeta,1));
  for i=1:size(cbeta,1)
    h = bar(H.dataplot, i, cbeta(i,:));
    set(h, 'FaceColor', groupcolor(i,:), 'FaceAlpha', 0.75);
  end
  
  % standard error
  %--------------------------------------------------------------
  CI = CI / 2;
  for j = 1:length(cbeta)
    line([j j], ([CI(j) -CI(j)] + cbeta(j)), 'LineWidth', 6, 'Color', H.col(2, :), 'Parent', H.dataplot)
  end
  set(H.dataplot, 'XLim', [0.4 (length(cbeta) + 0.6)], 'XTicklabel', '', 'XTick', [], 'YGrid','off')
  ylim(H.dataplot,'auto')
else
  ystr = 'raw data';
  

  c0 = H.SPM{1}.xCon(Ic).c;
  
  [indi, indj] = find(c0~=0);
  ind_X = unique(indi)';
  X = H.SPM{1}.xX.X;
  X = X(:,ind_X);

  covariate = 0;
  
  % check for covariates
  if ~isempty(H.SPM{1}.xX.iC) && numel(ind_X) <= 2
    for i=1:numel(ind_X)
      % contrast is defined at entries of iC
      if ~isempty(find(ind_X(i) == H.SPM{1}.xX.iC))
        covariate = 1;
      else
        covariate = 0;
      end
    end
  end
  c0 = c0(ind_X,:);

  % show scatter plot and linear fit
  if (covariate && numel(ind_X)<=2) || (exist('x','var') && numel(x)==size(X,1))
    axes(H.dataplot);
    cla
    
    xx = cell(numel(ind_X),1);

    % get column where groups are coded and use last found column
    if ~isempty(H.SPM{1}.xX.iH)
      n_groups = length(H.SPM{1}.xX.iH);
    elseif ~isempty(H.SPM{1}.xX.iC)
      n_groups = length(H.SPM{1}.xX.iC);
    end
    [rw,cl] = find(H.SPM{1}.xX.I == n_groups); 
    group_col = max(cl);
    
    % use existing x-variable if available
    if exist('x','var') && numel(x)==size(X,1)
      xx_array = [min(x) max(x)]; 
      for i=1:numel(ind_X)
        ind_data{i} = find(any(xX.X(:,xX.iH(i)),2));
        xx{i} = x(ind_data{i});
      end
    else
      xx_array = [min(X(X~=0)) max(X(X~=0))]; 
      for i=1:numel(ind_X)
        xx{i} = X(H.SPM{1}.xX.I(:,group_col)==i,i);
      end
    end

    col = [1 0 0;0 0 1; 0 1 0];
    yy = cell(numel(iH),1);
    
    for i=1:numel(iH)
      ind_data{i} = find(any(xX.X(:,xX.iH(i)),2));
      yy{i} = y(ind_data{i},:);
    end

    hold on
    for i=1:numel(ind_X)
      x2 = xx{i};
      y2 = mean(yy{i},2);
      ind_finite = isfinite(x2);
      x2 = x2(ind_finite); y2 = y2(ind_finite);
      if ~isempty(x2)
        plot(x2,y2,'.','MarkerSize',10,'Color',col(i,:));
        P = polyfit(x2,y2,1);     
        % plot trend line
        plot(xx_array,polyval(P,xx_array),'Color',col(i,:),'LineWidth',2);
      end
    end
    hold off
    
    if numel(ind_X)==2
      fprintf('Red: 1st group\tBlue: 2nd group\n'); 
    end
    
    xlim(H.dataplot,'auto')
    ylim(H.dataplot,'auto')
    set(H.dataplot,'XTickMode','auto');
    
  elseif ~isempty(iH) && numel(iH)>1 && (H.predicted < 0)
    if ~isempty(xX.iB)
      yy = y;
      vstruct = struct('style',4,'darkmode',all(1 - H.bkg_col),'datasymbol','.','I',xX.I);
    else
      yy = cell(numel(iH),1);

      for i=1:numel(iH)
        ind_data{i} = find(any(xX.X(:,xX.iH(i)),2));
        yy{i} = y(ind_data{i},:);
      end
      vstruct = struct('style',4,'darkmode',all(1 - H.bkg_col),'datasymbol','.');
    end
    
    axes(H.dataplot);
    cat_plot_boxplot(yy,vstruct);
    set(H.dataplot, 'XColor', 1 - H.bkg_col, 'YColor', 1 - H.bkg_col,...
      'Color', H.bkg_col, 'YGrid','on','Visible','on');    
  else
    h = plot(H.dataplot, y);
    set(h, 'Color', H.col(2, :))
    set(H.dataplot, 'XLim', [0 (length(y))], 'XTicklabel', '', 'XTick', [])
  end

end

% plot group coding for Anovas with more than 1 group and native data
if 0 %~isempty(iH) && numel(iH)>1 && (H.predicted < 0) && isempty(xX.iB)
  yl = get(H.dataplot,'YLim');
  pcol = gray(numel(iH)+2);
  for i=1:numel(iH)
    ind_data = find(any(xX.X(:,xX.iH(i)),2));
    
    % plot only if ind_data is a continuous row
    if all(diff(ind_data)==1)
      line([min(ind_data)-0.5 max(ind_data)+0.5], [yl(1) yl(1)], 'LineWidth', 6, 'Color', pcol(i+1,:), 'Parent', H.dataplot)
    end
  end
  
end

if ~isempty(Ic)
  xlabel(H.dataplot, H.SPM{round(ind / 2)}.xCon(Ic).name, 'FontSize', H.FS+3, 'Color', 1 - H.bkg_col, 'Interpreter', 'none')
end

if plot_mean
  ylabel(H.dataplot, sprintf('mean %s\ninside cluster',ystr), 'FontSize', H.FS+3, 'Color', 1 - H.bkg_col)
else
  ylabel(H.dataplot, sprintf('%s',ystr), 'FontSize', H.FS+3, 'Color', 1 - H.bkg_col)
end

hold(H.dataplot, 'off')

assignin('base', 'y', y);

%==========================================================================
function [y, cbeta, CI] = get_cluster_data(H, XYZ, ind)

SPM = H.SPM{round(ind / 2)};

Ic = [];
if isfield(H,'Ic')
  Ic = H.Ic;
else
  nm = H.S{1}.info(1).ff;
  
  % end with _0???.ext?
  if length(nm) > 4
    if strcmp(nm(length(nm) - 4:length(nm) - 3), '_0')
      Ic = str2double(nm(length(nm) - 3:length(nm)));
    end
  end
end

% get raw data and whiten
y = spm_data_read(SPM.xY.VY, 'xyz', XYZ);

% for adjusted or predicted data use model estimations
%------------------------------------------------------------------
if H.predicted >= 0

  y = spm_filter(SPM.xX.K, SPM.xX.W * y);
  R = spm_sp('r', SPM.xX.xKXs, y);
  
  beta  = spm_data_read(SPM.Vbeta,'xyz',XYZ);
  ResMS = spm_data_read(SPM.VResMS, 'xyz', XYZ);
  
  ResMS = mean(ResMS, 2);
  Bcov  = ResMS * SPM.xX.Bcov;
  
  % compute contrast of parameter estimates and 90% C.I.
  %------------------------------------------------------------------
  cbeta = SPM.xCon(Ic).c' * beta;
  cbeta = mean(cbeta, 2);
  
  CI = 1.6449; % = spm_invNcdf(1 - 0.05);
  SE = sqrt(diag(SPM.xCon(Ic).c' * Bcov * SPM.xCon(Ic).c));
  CI = CI * SE;
  
  % predicted or adjusted response
  %------------------------------------------------------------------
  if H.predicted == 1
    
    % fitted (predicted) data (Y = X1*beta)
    %--------------------------------------------------------------
    % this should be SPM.xX.xKXs.X instead of SPM.xX.X below
    Y = SPM.xX.X * SPM.xCon(Ic).c * pinv(SPM.xCon(Ic).c) * beta;
  else
    
    % fitted (corrected)  data (Y = X1o*beta)
    %--------------------------------------------------------------
    Y = spm_FcUtil('Yc', SPM.xCon(Ic), SPM.xX.xKXs, beta);
    
  end

  y = Y + R;
  y = cat_stat_nanmean(y, 2);
else
  y = cat_stat_nanmean(y, 2);

  % also adjust raw data
  if H.predicted==-2
    if ~isempty(Ic)     
      
      c0 = SPM.xCon(Ic).c;
      
      [indi, indj] = find(c0~=0);
      ind_X = unique(indi)';
    
      covariate = 0;
      
      % check for covariates
      if ~isempty(SPM.xX.iC) && numel(ind_X) <= 2
        for i=1:numel(ind_X)
          % contrast is defined at entries of iC
          if ~isempty(find(ind_X(i) == SPM.xX.iC))
            covariate = 1;
          else
            covariate = 0;
          end
        end
      end
      c0 = c0(ind_X,:);
      repeated_anova = ~isempty(SPM.xX.iB);

      % define subject effects and potential covariates
      G_columns = [SPM.xX.iB SPM.xX.iC];
      
      % only consider nuisance parameters and parameters where
      % contrast is defined
      for i=1:numel(ind_X)
        G_columns(find(G_columns==ind_X(i))) = [];
      end
      
      if ~isempty(G_columns)
        % remove nuisance effects from data
        G = SPM.xX.X(:,G_columns);
        G = G - mean(G);
        y = y - G*(pinv(G)*y);
      end

      % get column where groups are coded and use last found column
      if repeated_anova
        [rw,cl] = find(SPM.xX.I == length(SPM.xX.iB)); % find column which codes subject factor (length(xX.iB) -> n_subj)
        subj_col = cl(1);
        if subj_col == 3
          group_col = 2;
        else
          group_col = 3;
        end           
        n_groups = max(SPM.xX.I(:,group_col));
      elseif ~isempty(SPM.xX.iH)
        n_groups = length(SPM.xX.iH);
      elseif ~isempty(SPM.xX.iC)
        n_groups = length(SPM.xX.iC);
      end
      [rw,cl] = find(SPM.xX.I == n_groups); 
      
      if ~repeated_anova
        group_col = max(cl);
      end

      if repeated_anova || ((n_groups > 1) && covariate)
        beta0  = spm_data_read(SPM.Vbeta,'xyz',XYZ);
        beta   = mean(beta0,2);
        mean_group = zeros(n_groups,1);
        count_times = 1;
        for i=1:n_groups
          ind_group = find(SPM.xX.I(:,group_col) == i);
          if repeated_anova
            % find subjects effects in that group
            ind_subj = unique(SPM.xX.I(ind_group,subj_col));
            n_subj_group = numel(ind_subj);
            n_times = max(SPM.xX.I(ind_group,4));
            mean_group(i) = sum(beta(SPM.xX.iH(count_times:(count_times+n_times-1))))/n_times + ...
              sum(beta(SPM.xX.iB(ind_subj)))/n_subj_group;
            count_times = count_times + n_times;
          else
            mean_group(i) = beta(SPM.xX.iH(i));
          end
          y(ind_group,:) = y(ind_group,:) - mean(y(ind_group,:)) + mean_group(i);
        end
      end

    else
      fprintf('No adjustement for raw data possible!\n');
    end
  end
  cbeta = [];
  CI = [];
end

%==========================================================================
function txt = myDataCursorAtlas(obj, evt, H)

pos = get(evt, 'Position');

i = ismember(get(H.patch(1), 'vertices'), pos, 'rows');
node = find(i);
ind = 1;

if isempty(node)
  i = ismember(get(H.patch(3), 'vertices'), pos, 'rows');
  node = find(i);
  ind = 2;
end

if H.cursor_mode > 1

  sel_atlas = H.cursor_mode - 1;
  switch sel_atlas
    case 1
      txt = {'Desikan DK40'};
    case 2
      txt = {'Destrieux 2009'};
    case 3
      txt = {'HCP_MMP1'};
  end
  
  rdata_pos = H.rdata{sel_atlas}(node, ind);
  
  if rdata_pos == 0
    txt = {''};
    return
  end
  
  rcsv = H.rcsv{sel_atlas};
  
  for j = 2:size(rcsv, 1)
    if rdata_pos == rcsv{j, 1}
      txt = {txt{:} [H.S{ind}.side ' ' rcsv{j, 2}]};
      continue
    end
  end
else

  txt = {''};
  for sel_atlas=1:3
    switch sel_atlas
      case 1
        atlas_str = 'Desikan DK40: ';
      case 2
        atlas_str = 'Destrieux 2009: ';
      case 3
        atlas_str = 'HCP_MMP1: ';
    end
    
    rdata_pos = H.rdata{sel_atlas}(node, ind);
    if rdata_pos == 0
      txt = {''};
      return
    end

    rcsv = H.rcsv{sel_atlas};
    
    for j = 2:size(rcsv, 1)
      if rdata_pos == rcsv{j, 1}
        txt = {txt{:} [atlas_str H.S{ind}.side ' ' rcsv{j, 2}]};
        continue
      end
    end
  end
  txt = txt(2:4);
end

%==========================================================================
function myDeleteFcn(obj, evt, renderer)
try rotate3d(get(obj, 'parent'), 'off'); end
set(ancestor(obj, 'figure'), 'Renderer', renderer);

%==========================================================================
function s = remove_zeros(s)

pos = length(s);
while pos > 1
  if strcmp(s(pos), '0')
    s(pos) = '';
    pos = pos - 1;
  else break
  end
end

%==========================================================================
function C = find_connected_component(A, T);
% find connected components
% FORMAT C = find_connected_component(A,T)
% A      - a [nxn[ (reduced) adjacency matrix
% T      - a [nx1] data vector (using NaNs or logicals), n = #vertices
%
% C      - a [nx1] vector of cluster indices
%
% modified version from spm_mesh_clusters.m 5065 2012-11-16 20:00:21Z guillaume
%


%-Input parameters
%--------------------------------------------------------------------------
if ~islogical(T)
  T = ~isnan(T);
end

A1 = A;
A1(~T, :) = [];
A1(:, ~T) = [];

%-And perform Dulmage-Mendelsohn decomposition to find connected components
%--------------------------------------------------------------------------
[p, q, r] = dmperm(A1);
N = diff(r);
CC = zeros(size(A1, 1), 1);
for i = 1:length(r) - 1
  CC(p(r(i):r(i + 1) - 1)) = i;
end
C = NaN(numel(T), 1);
C(T) = CC;

%-Sort connected component labels according to their size
%--------------------------------------------------------------------------
[N, ni] = sort(N(:), 1, 'descend');
[ni, ni] = sort(ni);
C(T) = ni(C(T));

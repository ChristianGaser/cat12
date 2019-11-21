function varargout = cat_surf_results(action, varargin)
% Visualise results for both hemispheres of surface-based analysis 
% (preferable on log P-maps). 
%
% FORMAT y = cat_surf_results('Disp',leftSurface,rightSurface)
% leftSurface  - a GIfTI filename/object or patch structure
% rightSurface - a GIfTI filename/object or patch structure
%
% y            - adjusted, predicted or raw response
%
%
% Futher actions for batch mode: 
%  * cat_surf_results('batch',job) 
%    See cat_conf_stools. 
%
%  * cat_surf_results('surface',1..4) 
%    Select surface type.
%
%  * cat_surf_results('texture',0..2)
%    Select surface underlay. 
%    0 - none, 1 - mean curvature, 2 - sulcal depth
%
%  * cat_surf_results('view',1..3)
%    Select render view.
%    1 - topview, 2 - bottomview, 3 - sideview
%
%  * cat_surf_results('colormap',1..4)
%    Select overlay colormap.
%    1 - jet, 2 - hot, 3 - hsv, 4 - cold-hot
%
%  * cat_surf_results('invcolormap'[,0..1])
%    Default (0) or inverts colormap (1). Toggles without input.
%    
%  * cat_surf_results('background',[,0..2])
%    White (1) or black (0|2) background. Toggles without input.
%
%  * cat_surf_results('showfilename'[,0..1]); 
%    Show (1) or not show (0) surface name in figure. Toggles without input.
%
%  * cat_surf_results('threshold'[,0..4]);
%    Define statistical threshold. 
%    0 - none, 1.3 - 0.05, 2 - 0.01, 3 - 0.001
%
%  * cat_surf_results('hide_neg',[,0..1]); 
%    Hide negative results (1) or show everything (0). Toggles without input.
%
%_______________________________________________________________________
% Christian Gaser & Robert Dahnke (batch-mode)
% $Id$

global H y

%-Input parameters
%--------------------------------------------------------------------------
if ~nargin, action = 'Disp'; end

if ~ischar(action)
    varargin = {action varargin{:}};
    action = 'Disp';
end


%-Action
%--------------------------------------------------------------------------
switch lower(action)
    
    %-Display
    %======================================================================
    case 'disp'
      
        % remove any existing data
        if isfield(H,'S')
           H = rmfield(H,'S');
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
        H.col          = [.8 .8 .8; 1 .5 .5];
        H.FS           = cat_get_defaults('extopts.fontsize');
        H.n_surf       = 1;
        H.thresh_value = 0;
        H.cursor_mode  = 1;
        H.text_mode    = 1;
        H.border_mode  = 0;
        H.is32k        = 0;
        H.str32k       = '';
        H.SPM_found    = 1;
        H.surf_sel     = 1;
% colorbar addon histogram & values
% 
      
        % positions
        WS = spm('Winsize', 'Graphics');
        SS = get(0, 'Screensize');
        if 2.6 * WS(3) > SS(3)
            WS(3) = WS(3) / (2.6 * WS(3) / SS(3));
        end
        
        % result window with 5 surface views and alternative positions without top view and  only with lateral views
        H.viewpos = {[0.025 0.450 0.375 0.375;  0.025 0.450 0.375 0.375;  0.025 2.000 0.375 0.375],... % lh medial
                     [0.025 0.025 0.375 0.375;  0.025 0.025 0.375 0.375;  0.175 0.350 0.175 0.350],... % lh lateral
                     [0.600 0.450 0.375 0.375;  0.600 0.450 0.375 0.375;  0.600 2.000 0.375 0.375],... % rh medial
                     [0.600 0.025 0.375 0.375;  0.600 0.025 0.375 0.375;  0.675 0.350 0.175 0.350],... % rh lateral
                     [0.300 0.150 0.400 0.500;  0.300 2.000 0.400 0.500;  0.300 2.000 0.400 0.500],... % lh+rh top
                     [0.400 0.750 0.200 0.225;  0.400 0.300 0.200 0.225;  0.400 0.750 0.200 0.225]};   % data plot
        
        % change size and position of flatmaps for >= R20014b
        if spm_check_version('matlab', '8.4') >= 0
            H.viewpos{2}(3, :) = [-0.075 0.150 0.650 0.650]; % lh lateral
            H.viewpos{4}(3, :) = [0.425 0.150 0.650 0.650]; % rh lateral
        end
        
        % figure 1 with result window
        H.pos{1} = struct( ...
            'fig', [10 10 round(2.6*WS(3)) WS(3)], ... % figure
            'cbar', [0.400 -0.150 0.200 0.300; 0.440 0.025 0.120 0.120]);% colorbar
        
        % figure 2 with GUI
        H.pos{2} = struct(...
          'fig',    [2*WS(3)+10 10 0.6*WS(3) WS(3)],... 
          'sel',    [0.050 0.935 0.900 0.050],...[0.290 0.930 0.425 0.050],...
          'nam',    [0.050 0.875 0.900 0.050],...
          'scaling',[0.050 0.260 0.425 0.050],... % developer
          'surf',   [0.050 0.800 0.425 0.050],'mview',   [0.525 0.800 0.425 0.050],... 
          'text',   [0.050 0.725 0.425 0.050],'thresh',  [0.525 0.725 0.425 0.050],... 
          'cmap',   [0.050 0.650 0.425 0.050],'atlas',   [0.525 0.650 0.425 0.050],...
          'cursor', [0.050 0.575 0.425 0.050],'border',  [0.525 0.575 0.425 0.050],...
          'nocbar', [0.050 0.500 0.425 0.050],'transp',  [0.525 0.500 0.425 0.050],... 
          'info',   [0.050 0.425 0.425 0.050],'bkg',     [0.525 0.425 0.425 0.050],... 
          'inv',    [0.050 0.350 0.425 0.050],'hide_neg',[0.525 0.350 0.425 0.050],...
          'ovmin',  [0.050 0.125 0.425 0.100],'ovmax',   [0.525 0.125 0.425 0.100],... 
          'save',   [0.050 0.050 0.425 0.050],'close',   [0.525 0.050 0.425 0.050]);
        
  		  H.figure = figure(22);
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
            
        % closing all windows
        H.close = uicontrol(H.panel(2), ...
            'String', 'Close', 'Units', 'normalized', ...
            'Position', H.pos{2}.close, ...
            'Style', 'Pushbutton', 'HorizontalAlignment', 'center', ...
            'Callback', 'close(22);', ...
            'FontSize',H.FS,'ForegroundColor','red',...
            'ToolTipString', 'Close windows', ...
            'Interruptible', 'on', 'Enable', 'on');

        % select results for lh and rh
        H.sel = uicontrol(H.panel(2), ...
            'String', 'Select Surface Data', 'Units', 'normalized', ...
            'Position', H.pos{2}.sel, ...
            'Style', 'Pushbutton', 'HorizontalAlignment', 'center', ...
            'Callback', @select_data, ...
            'FontSize',H.FS,...
            'ToolTipString', 'Select results (up to 24) for both hemispheres (e.g. log-p maps)', ...
            'Interruptible', 'on', 'Enable', 'on');
        
        H.save = uicontrol(H.panel(2), ...
            'String', 'Save', 'Units', 'normalized', ...
            'Position', H.pos{2}.save, ...
            'Style', 'Pushbutton', 'HorizontalAlignment', 'center', ...
            'Callback', {@save_image}, ...
            'FontSize',H.FS,...
            'ToolTipString', 'Save png image', ...
            'Interruptible', 'on', 'Enable', 'off');

        str = {'Surface', 'Central', 'Inflated', 'Dartel', 'Flatmap'};
        tmp = {{@select_surf, 1}, ...
               {@select_surf, 2}, ...
               {@select_surf, 3}, ...
               {@select_surf, 4}};
        
        % underlying surface
        H.surf = uicontrol(H.panel(2), ...
            'String', str, 'Units', 'normalized', ...
            'Position', H.pos{2}.surf, 'UserData', tmp, ...
            'Style', 'PopUp', 'HorizontalAlignment', 'center', ...
            'Callback', 'spm(''PopUpCB'',gcbo)', ...
            'FontSize',H.FS,...
            'ToolTipString', 'Underlying Surface', ...
            'Interruptible', 'on', 'Enable', 'off');
        
        str = {'Threshold', 'No threshold', 'P<0.05', 'P<0.01', 'P<0.001'};
        tmp = {{@select_thresh, 0}, ...
               {@select_thresh, 1.3}, ...
               {@select_thresh, 2}, ...
               {@select_thresh, 3}};
        
        % threshold
        H.thresh = uicontrol(H.panel(2), ...
            'String', str, 'Units', 'normalized', ...
            'Position', H.pos{2}.thresh, 'UserData', tmp, ...
            'Style', 'PopUp', 'HorizontalAlignment', 'center', ...
            'Callback', 'spm(''PopUpCB'',gcbo)', ...
            'FontSize',H.FS,...
            'ToolTipString', 'Threshold', ...
            'Interruptible', 'on', 'Enable', 'off');
        
        str = {'Colormap', 'jet', 'hot', 'hsv', 'cold-hot'};
        tmp = {{@select_cmap, 1}, ...
               {@select_cmap, 2}, ...
               {@select_cmap, 3}, ...
               {@select_cmap, 4}};
        
        % colormap
        H.cmap = uicontrol(H.panel(2), ...
            'String', str, 'Units', 'normalized', ...
            'Position', H.pos{2}.cmap, 'UserData', tmp, ...
            'Style', 'PopUp', 'HorizontalAlignment', 'center', ...
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
            'Style', 'PopUp', 'HorizontalAlignment', 'center', ...
            'Callback', 'spm(''PopUpCB'',gcbo)', ...
            'FontSize',H.FS,...
            'ToolTipString', 'Atlas Labeling', ...
            'Interruptible', 'on', 'Enable', 'off');
        
        str = {'Data Cursor', 'Disable data cursor', 'Atlas regions: Desikan-Killiany DK40', ...
            'Atlas regions: Destrieux 2009', 'Atlas region: HCP Multi-Modal Parcellation', ...
            'Plot data at vertex', ...
            'Plot mean data inside cluster', 'Enable/Disable rotate3d'};
        tmp = {{@select_cursor, 0}, ...
               {@select_cursor, 1}, ...
               {@select_cursor, 2}, ...
               {@select_cursor, 3}, ...
               {@select_cursor, 4}, ...
               {@select_cursor, 5}, ...
               {@select_cursor, 6}};
        
        % data cursor for data plotting and atlas names
        H.cursor = uicontrol(H.panel(2), ...
            'String', str, 'Units', 'normalized', ...
            'Position', H.pos{2}.cursor, 'UserData', tmp, ...
            'Style', 'PopUp', 'HorizontalAlignment', 'center', ...
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
            'Style', 'PopUp', 'HorizontalAlignment', 'center', ...
            'Callback', 'spm(''PopUpCB'',gcbo)', ...
            'FontSize',H.FS,...
            'ToolTipString', 'Select View', ...
            'Interruptible', 'on', 'Enable', 'off');
        
        str = {'Underlay', 'Mean curvature', 'Sulcal depth'};
        tmp = {{@select_texture, 1}, ...
               {@select_texture, 2}};
        
        % underlying texture
        H.text = uicontrol(H.panel(2), ...
            'String', str, 'Units', 'normalized', ...
            'Position', H.pos{2}.text, 'UserData', tmp, ...
            'Style', 'PopUp', 'HorizontalAlignment', 'center', ...
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
            'Style', 'PopUp', 'HorizontalAlignment', 'center', ...
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
        
        if cat_get_defaults('extopts.expertgui')<2
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
                
        if nargin >= 3
            
            H.S{1}.name = varargin{1};
            H.S{2}.name = varargin{2};
            
            [pth{1}, nm1, ext1] = spm_fileparts(H.S{1}.name(1, :));
            [pth{2}, nm2, ext2] = spm_fileparts(H.S{2}.name(1, :));
            
            % SPM.mat found for both hemispheres (not working yet)
            if strcmp([nm1 ext1], 'SPM.mat') || strcmp([nm2 ext2], 'SPM.mat')
                H.logP = 0;
                
                if strcmp([nm1 ext1], 'SPM.mat')
                    ind = 1;
                else ind = 2; end
                
                swd1 = pwd;
                spm_figure('GetWin', 'Interactive');
                cd(pth{ind})
                xSPM.swd = pwd;
                [xSPM, v] = spm_getSPM(xSPM);
                cd(swd1);
                
                dat = struct('XYZ', v.XYZ, ...
                    't', v.Z', ...
                    'mat', v.M, ...
                    'dim', v.DIM, ...
                    'dat', v.Z');
                
                H.S{ind}.info = cat_surf_info(H.S{ind}.name, 0);
                g = gifti(H.S{ind}.info.Pmesh);
                
                mat = v.M;
                V = g.vertices;
                XYZ = double(inv(mat) * [V'; ones(1, size(V, 1))]);
                H.S{ind}.Y = spm_sample_vol(Y, XYZ(1, :), XYZ(2, :), XYZ(3, :), 0)';
                H.S{ind}.Y = spm_mesh_project(g.vertices, dat)';
            else
                
                H.logP = 1;
                meshes_merged = H.merged;
                
                for ind = 1:2
                    
                    % read meshes
                    H.S{ind}.info = cat_surf_info(H.S{ind}.name, 1);
                    
                    if H.S{ind}.info(1).nvertices == 64984
                        H.str32k = '_32k';
                        H.is32k = 1;
                    else
                        H.str32k = '';
                        H.is32k = 0;
                    end

                    if ~strcmp(H.S{ind}.info(1).side, 'lh') && ~strcmp(H.S{ind}.info(1).side, 'rh')%strcmp(H.S{ind}.info(1).side, 'mesh')
                        meshes_merged = 1;
                        if ind == 1
                            H.S{ind}.info(1).side = 'lh';
                            H.S{ind}.info(1).Pmesh = fullfile(spm('dir'), 'toolbox', 'cat12', ...
                                ['templates_surfaces' H.str32k], 'lh.central.freesurfer.gii');
                        else
                            H.S{ind}.info(1).side = 'rh';
                            H.S{ind}.info(1).Pmesh = fullfile(spm('dir'), 'toolbox', 'cat12', ....
                                ['templates_surfaces' H.str32k], 'rh.central.freesurfer.gii');
                        end
                    end
                    H.S{ind}.M = gifti(H.S{ind}.info(1).Pmesh);
                    
                    % get adjacency information
                    H.S{ind}.A = spm_mesh_adjacency(H.S{ind}.M);
                    
                    % cdata found?
                    if meshes_merged
                        if ind == 1
                            try
                                Y = spm_data_read(spm_data_hdr_read(H.S{ind}.name));
                            catch
                                error('No data in surfaces found.');
                            end
                            if H.is32k
                                H.S{1}.Y = Y(1:32492, :);
                                H.S{2}.Y = Y(32493:end, :);
                            else
                                H.S{1}.Y = Y(1:163842, :);
                                H.S{2}.Y = Y(163843:end, :);
                            end
                        end
                    else
                        gind = gifti(H.S{ind}.name);
                        if isfield(gind,'cdata')
                            H.S{ind}.Y = spm_data_read(spm_data_hdr_read(H.S{ind}.name));
                        else
                            if ind == 1
                                gind = gifti(H.S{2}.name);
                                if isfield(gind,'cdata')
                                    if isnumeric(gind.cdata)
                                        H.S{1}.Y = zeros(size(gind.cdata));
                                    else
                                        H.S{1}.Y = zeros(gind.cdata.dim);
                                    end
                                    H.S{1}.name = '';
                                else
                                    error('No data in surfaces found.');
                                end
                            else
                                gind = gifti(H.S{1}.name);
                                if isfield(gind,'cdata')
                                    if isnumeric(gind.cdata)
                                        H.S{2}.Y = zeros(size(gind.cdata));
                                    else
                                        H.S{2}.Y = zeros(gind.cdata.dim);
                                    end
                                    H.S{2}.name = '';
                                else
                                    error('No data in surfaces found.');
                                end
                            end
                        end
                    end
                                      
                    % check whether name contains 'log' that indicates a logP file
                    for i = 1:size(H.S{ind}.name, 1)
                        if isempty(strfind(H.S{ind}.info(i).ff, 'log'))
                            H.logP = 0;
                        end
                    end
                end                    
            end
            
            % rescue original name for later result selection
            H.S1 = H.S{1};
            H.S2 = H.S{2};
            
            H.n_surf = max(numel(H.S{1}.info), numel(H.S{2}.info));
            H.view = 1;
            H.show_transp = 1;
            H.disable_cbar = 0;
            H.white_bgk = 0;
            H.show_info = 0;
            
            % result selection or RGB overlay if more than one result was loaded
            if H.n_surf > 1
                % pre-select 1st mesh if we cannot use RGB overlay
                %if H.n_surf > 3
                if 1
                    sel = 1;
                    if isempty(H.S1.name)
                        error('Do not mix meshes with different resolutions (i.e. 164k vs. 32k)');
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
                if 1, H.n_surf = 1; end
                %if H.n_surf > 3, H.n_surf = 1; end
                
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
            end
            
            display_results_all;
            
            H.SPM_found = 1;
            for i = 1:2
                SPM_name = fullfile(H.S{i}.info(1).pp, 'SPM.mat');
                
                % SPM.mat exist?
                if ~isempty(H.S{i}.name)
                    H.SPM_found = 0;
                end
            end

            % Don't allow plot functions for RGB maps or if SPM.mat was not found
            if H.n_surf > 1 && H.SPM_found
                str = {'Data Cursor', 'Disable data cursor', 'Atlas regions: Desikan-Killiany DK40', ...
                    'Atlas regions: Destrieux 2009', 'Atlas region: HCP Multi-Modal Parcellation', ...
                    'Enable/Disable rotate3d'};
                tmp = {{@select_cursor, 0}, ...
                       {@select_cursor, 1}, ...
                       {@select_cursor, 2}, ...
                       {@select_cursor, 3}, ...
                       {@select_cursor, 5}};
                
                H.cursor = uicontrol(H.panel(2), ...
                     'String', str, 'Units', 'normalized', ...
                    'Position', H.pos{2}.cursor, 'UserData', tmp, ...
                    'Style', 'PopUp', 'HorizontalAlignment', 'center', ...
                    'Callback', 'spm(''PopUpCB'',gcbo)', ...
                    'FontSize',H.FS,...
                    'ToolTipString', 'Data Cursor Mode', ...
                    'Interruptible', 'on', 'Enable', 'off');
            end
            
            % enable some menus only if mesh data can be assumed to be resampled
            if (length(H.S{1}.Y) == 32492 || length(H.S{1}.Y) == 163842)
                set(H.surf, 'Enable', 'on');
                set(H.text, 'Enable', 'on');
                set(H.cursor, 'Enable', 'on');
                set(H.border, 'Enable', 'on');
            end
            
            set(H.save, 'Enable', 'on');
            set(H.mview, 'Enable', 'on');
            set(H.nocbar, 'Enable', 'on');
            set(H.bkg, 'Enable', 'on');
            set(H.transp, 'Enable', 'on');
            set(H.info, 'Enable', 'on');
            if isfield(H,'scaling')
              set(H.scaling, 'Enable', 'on');
            end
            
            if min(min(H.S{1}.Y(:)), min(H.S{2}.Y(:))) < 0 & H.n_surf == 1
                set(H.inv, 'Enable', 'on');
                set(H.hide_neg, 'Enable', 'on');
                set(H.hide_neg, 'Value', 0);
            end
            
            if H.n_surf == 1
                set(H.cmap, 'Enable', 'on');
            end
            
            H.rdata{1} = [];
            H.rdata{2} = [];
            H.rdata{3} = [];
            for ind = 1:2
                atlas_name = fullfile(spm('dir'), 'toolbox', 'cat12', ['atlases_surfaces' H.str32k], ...
                [H.S{ind}.info(1).side '.aparc_DK40.freesurfer.annot']);
                [vertices, rdata0, colortable, rcsv1] = cat_io_FreeSurfer('read_annotation', atlas_name);
                H.rdata{1} = [H.rdata{1} rdata0];
                atlas_name = fullfile(spm('dir'), 'toolbox', 'cat12', ['atlases_surfaces' H.str32k], ...
                [H.S{ind}.info(1).side '.aparc_a2009s.freesurfer.annot']);
                [vertices, rdata0, colortable, rcsv2] = cat_io_FreeSurfer('read_annotation', atlas_name);
                H.rdata{2} = [H.rdata{2} rdata0];
                atlas_name = fullfile(spm('dir'), 'toolbox', 'cat12', ['atlases_surfaces' H.str32k], ...
                [H.S{ind}.info(1).side '.aparc_HCP_MMP1.freesurfer.annot']);
                [vertices, rdata0, colortable, rcsv3] = cat_io_FreeSurfer('read_annotation', atlas_name);
                H.rdata{3} = [H.rdata{3} rdata0];
            end
            H.rcsv{1} = rcsv1;
            H.rcsv{2} = rcsv2;
            H.rcsv{3} = rcsv3;
            
            H.dcm_obj = datacursormode(H.figure);
            set(H.dcm_obj, 'Enable', 'on', 'SnapToDataVertex', 'on', ...
                'DisplayStyle', 'datatip', 'Updatefcn', {@myDataCursorAtlas, H});
            
        end
        
        if nargout, varargout{1} = y; end
        
    %-ColourBar
    %======================================================================
    case {'colourbar', 'colorbar'}
        if nargin>1
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
            if ~get(H.nocbar, 'Value')
              cat_surf_results('hist',0); 
            end
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
            %            H.colourbar = colorbar('peer',gca,'NorthOutside');
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
            case {1,2,3,4},  cmap = cm; 
            case 'jet',      cmap = 1; 
            case 'hot',      cmap = 2; 
            case 'hsv',      cmap = 3; 
            case 'cold-hot', cmap = 4; 
            otherwise
              error('Unknown colormap\n');
          end          
          select_cmap(cmap);
        end
      
        
    %-CLim
    %======================================================================
    case 'clim'
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
        if length(varargin) == 1
            c = getappdata(H.patch, 'clim');
            if ~isempty(c), c = c(2:3); end
            varargout = {c};
            return;
        else
            if strcmp(varargin{2}, 'on') || isempty(varargin{2}) || any(~isfinite(varargin{2}))
                setappdata(H.patch, 'clim', [false NaN NaN]);
            else
                setappdata(H.patch, 'clim', [true varargin{2}]);
            end
            d = getappdata(H.patch, 'data');
            H = updateTexture(H, d);
            
        end
        
        if nargin > 1 && isnumeric(varargin{2}) && numel(varargin{2}) == 2
            caxis(H.axis, varargin{2});
        else
            caxis(H.axis, [min(d), max(d)])
        end
        
        
    %-CLip
    %======================================================================
    case 'clip'
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
        if length(varargin) == 1
            c = getappdata(H.patch, 'clip');
            if ~isempty(c), c = c(2:3); end
            varargout = {c};
            return;
        else
            if isempty(varargin{2}) || any(~isfinite(varargin{2}))
                for ind = 1:5
                    setappdata(H.patch(ind), 'clip', [false NaN NaN]);
                end
            else
                for ind = 1:5
                    setappdata(H.patch(ind), 'clip', [true varargin{2}]);
                end
            end
            for ind = 1:5
                d = getappdata(H.patch, 'data');
                H = updateTexture(H, ind, d);
            end
        end
    
        
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
              otherwise
                error('unkown H.datascale %s.\n',H.datascale);
            end
        end
      end

      %% update textures of each patch
      for i=1:numel(H.patch) 
        setappdata(H.patch(i), 'clim',H.clim);
        H = updateTexture(H, i ); 
      end
      set(H.str_min, 'String', sprintf('%g',H.clim(2)));
      set(H.str_max, 'String', sprintf('%g',H.clim(3)));
      show_colorbar(H); 
    
      
    %- print histogram
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
            if mode
              tpos = {[0.49 0.015 0.065 0.03],[0.565 0.015 0.065 0.03],[0.4 0.015 0.12 0.03]};
            else
              tpos = {[0.4 0.015 0.12 0.06],[0.445 0.015 0.065 0.06],[0.497 0.015 0.065 0.06],[0.55 0.015 0.065 0.06],[0.60 0.015 0.065 0.06]};
            end
            % histrogram axis 
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
                  fill(H.histax,[quantile(1)   quantile(2)   quantile(2)          quantile(1)],...
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
        texure = varargin{1};
        if any(texure == 1:2)
          select_texture(texure);
        elseif texure == 0
          cat_surf_results('transparency',0);
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
                        cat_surf_results('background');
                    end
                case {0,2,'black'}
                    if get(H.bkg, 'Value')==1
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
        
        % create window
        select_data([],[],char(job.cdata));
        
        % set parameter
        FN = {'surface','view','texture','transparency','colorbar','colormap','invcolormap','background','showfilename','clims'}; 
        for fni=1:numel(FN)
          %%
          if isfield(job,'render') && isfield(job.render,FN{fni})
            cat_surf_results(FN{fni},job.render.(FN{fni})); 
          end
        end
        FN = {'threshold','hide_neg'}; 
        for fni=1:numel(FN)
          if isfield(job,'stat') && isfield(job.stat,FN{fni})
            cat_surf_results(FN{fni},job.stat.(FN{fni})); 
          end
        end
        
        % save result
        if isfield(job,'fparts')
          fparts = job.fparts; 
          files = cat_surf_results('print',fparts);
        else
          files = cat_surf_results('print');
        end
        varargout{1}.png = files; 
        
        % close figure after export
        clear -globalvar H; 
        close(22); 
        
        
    %- save image
    %======================================================================
    case 'print'
        
        if nargin>1
          fparts = varargin{1};
        else
          fparts.outdir = {''};
          fparts.prefix = '';
          fparts.suffix = '';
        end
        
        if nargin>2
          imgs = varargin{2};
        else
          imgs = inf; 
        end
        maximg = numel(H.S1.info);
        if isinf(imgs), imgs = 1:maximg; end 
        imgs(imgs<0 | imgs>maximg) = []; 
        
        %% print images
        for fi=1:numel(imgs)
          if fi>1, select_results(imgs(fi)); end
          
          [pp,ff] = spm_fileparts( H.S1.name(imgs(fi),:) );
          if isempty(fparts.outdir{1})
            fparts.outdir{1} = pp; 
          end
          filenames{imgs(fi)} = fullfile(fparts.outdir{1},[fparts.prefix ff fparts.suffix '.png']); %#ok<AGROW>
          
          save_image(1,1,filenames{imgs(fi)});
          
          % display image
          fprintf('  Display %s\n',...
            spm_file(filenames{imgs(fi)},'link',[...
              'try, delete(314); end; fh=figure(314); img = imread(''%s''); pos = get(fh,''Position'');'...
              'pos(4) = pos(3) * size(img,1)./size(img,2);' ...
              'set(fh,''name'',''cat_surf_result_png'', '...
              '  ''menubar'',''none'',''toolbar'',''none'',''Position'',pos); '...
              'image(img); set(gca,''Position'',[0 0 1 1],''visible'',''off''); '])); 
        end
        
        varargout{1} = filenames; 
        
        
  otherwise   
        error('Unknown action "%s"!\n',action); 
end 

       

%-----------------------------------------------------------------------
function Ho = select_thresh(thresh)
%-----------------------------------------------------------------------
global H

H.thresh_value = thresh;
H.clip = [true -thresh thresh];

H.no_neg = get(H.hide_neg, 'Value');

% get min value for both hemispheres
min_d = min(min(min(getappdata(H.patch(1), 'data'))), min(min(getappdata(H.patch(3), 'data'))));
clim = getappdata(H.patch(1), 'clim');

% rather use NaN values for zero threshold
if thresh == 0
    H.clip = [false NaN NaN];
end

if H.no_neg
    H.clip = [true -Inf thresh];
    clim = [true 0 clim(3)];
    set(H.slider_min, 'Value', 0);
end

for ind = 1:5
    if min_d > -thresh
        setappdata(H.patch(ind), 'clim', [true thresh clim(3)]);
    elseif thresh == 0
        setappdata(H.patch(ind), 'clim', [true -clim(3) clim(3)]);
    end
    
    setappdata(H.patch(ind), 'clip', H.clip);
    col = getappdata(H.patch(ind), 'col');
    d = getappdata(H.patch(ind), 'data');
    min_d = min(min_d, min(d(:)));
    H = updateTexture(H, ind, d, col, H.show_transp);
end

set(H.slider_min, 'Value', H.clim(2))
set(H.str_min, 'String', sprintf('%g',H.clim(2)));

set(H.atlas, 'Enable', 'on');

if ~H.disable_cbar
    H = show_colorbar(H);
end
if nargout, Ho = H; end

%-----------------------------------------------------------------------
function disphist(type)
%-----------------------------------------------------------------------
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


%-----------------------------------------------------------------------
function Ho = select_cmap(cmap)
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
        col = [1 - hot(128); (hot(128))];
end

for ind = 1:5
    setappdata(H.patch(ind), 'col', col);
    d = getappdata(H.patch(ind), 'data');
    H = updateTexture(H, ind, d, col, H.show_transp);
end

if ~H.disable_cbar
    H = show_colorbar(H);
end
if nargout, Ho = H; end

%-----------------------------------------------------------------------
function Ho = select_atlas(atlas)
%-----------------------------------------------------------------------
global H

% get threshold from clipping
thresh = [0 0];
if ~isempty(H.clip)
    if ~isnan(H.clip(2)) & ~isnan(H.clip(3))
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
        dp = d >= thresh(2); indp = find(dp);
        dn = d <= thresh(1); indn = find(dn);
        
        % go through pos. effects
        if ~isempty(indp)
            
            C = find_connected_component(A, dp);
            C = C(indp);
            rdata2 = rdata(indp);
            
            fprintf('\n\n______________________________________________________\n');
            fprintf('%s: Positive effects in %s', atlas_name, H.S{round(ind / 2)}.info(1).side);
            fprintf('\n%s', spm_str_manip(H.S{round(ind / 2)}.info(indsurf).fname, 'k50d'));
            fprintf('\n______________________________________________________\n\n');
            
            if H.logP, fprintf('%7s\t%8s\t%s\n', 'P-value', 'Size', 'Overlap of atlas region');
            else, fprintf('%7s\t%8s\t%s\n', 'Value  ', 'Size', 'Overlap of atlas region'); end
            
            for i = 1:max(C)
                N = find(C == i);
                k = length(N);
                
                dmax = d(indp); dmax = max(dmax(N));
                
                if H.logP, fprintf('\n%1.5f\t%8d', 10^(-dmax), k);
                else, fprintf('\n%6.1f\t%8d', dmax, k); end
                
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
                    if roi_size(jj(j)) > 1
                        if j == 1, fprintf('\t%3.0f%s\t%s\n', roi_size(jj(j)), '%', rcsv{jj(j) + 1, 2});
                        else, fprintf('%7s\t%8s\t%3.0f%s\t%s\n', '       ', '        ', ...
                                roi_size(jj(j)), '%', rcsv{jj(j) + 1, 2});
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
            
            if H.logP, fprintf('%7s\t%8s\t%s\n', 'P-value', 'Size', 'Overlap of atlas region');
            else, fprintf('%7s\t%8s\t%s\n', 'Value  ', 'Size', 'Overlap of atlas region'); end
            
            for i = 1:max(C)
                N = find(C == i);
                k = length(N);
                
                dmin = d(indn); dmin = min(dmin(N));
                if H.logP, fprintf('\n%1.5f\t%8d', 10^(dmin), k);
                else, fprintf('\n%6.1f\t%8d', -dmin, k); end
                
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
                    if roi_size(jj(j)) > 1
                        if j == 1, fprintf('\t%3.0f%s\t%s\n', roi_size(jj(j)), '%', rcsv{jj(j) + 1, 2});
                        else, fprintf('%7s\t%8s\t%3.0f%s\t%s\n', '       ', '        ', ...
                                roi_size(jj(j)), '%', rcsv{jj(j) + 1, 2});
                        end
                    end
                end
            end
        end
    end
end
if nargout, Ho = H; end

%-----------------------------------------------------------------------
function Ho = select_results(sel)
%-----------------------------------------------------------------------
global H

clearDataCursorPlot(H);

H.S{1}.name = H.S1.name(sel, :);
H.S{2}.name = H.S2.name(sel, :);
H.S{1}.Y = H.S1.Y(:, sel);
H.S{2}.Y = H.S2.Y(:, sel);

H.S{1}.info = cat_surf_info(H.S{1}.name, 0);
H.S{2}.info = cat_surf_info(H.S{2}.name, 0);

% check whether data for left or right hemipshere are all non-zero
ind1 = find(H.S{1}.Y(:) ~= 0);
ind2 = find(H.S{2}.Y(:) ~= 0);

% estimate min value > 0 and min/max values
if ~isempty(ind1) & ~isempty(ind2)
    H.S{1}.thresh = min(H.S{1}.Y(H.S{1}.Y(:) > 0));
    H.S{1}.thresh = min(H.S{1}.thresh, min(H.S{2}.Y(H.S{2}.Y(:) > 0)));
    H.S{1}.min = min(min(H.S{1}.Y(~isinf(H.S{1}.Y))), min(H.S{2}.Y(~isinf(H.S{2}.Y))));
    H.S{1}.max = max(max(H.S{1}.Y(~isinf(H.S{1}.Y))), max(H.S{2}.Y(~isinf(H.S{2}.Y))));
elseif isempty(ind1)
    H.S{1}.thresh = min(H.S{2}.Y(H.S{2}.Y(:) > 0));
    H.S{1}.min = min(H.S{2}.Y(~isinf(H.S{2}.Y)));
    H.S{1}.max = max(H.S{2}.Y(~isinf(H.S{2}.Y)));
elseif isempty(ind2)
    H.S{1}.thresh = min(H.S{1}.Y(H.S{1}.Y(:) > 0));
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

% only apply thresholds that are slightly larger than zero
if H.S{1}.thresh > 0.00015 & H.thresh_value == 0
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
end

H.n_surf = 1;

for ind = 1:5
    if H.S{1}.thresh > 0.00015
        setappdata(H.patch(ind), 'clip', H.clip);
    end
    setappdata(H.patch(ind), 'clim', H.clim);
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
end

set(H.slider_min, 'Value', H.clim(2));
set(H.slider_max, 'Value', H.clim(3));
set(H.str_min, 'String', sprintf('%g',H.clim(2)));
set(H.str_max, 'String', sprintf('%g',H.clim(3)));

% update file information and colorbar
checkbox_info;

if ~H.disable_cbar
    H = show_colorbar(H);
end

% print selected filename
cla(H.nam);
axis(H.nam, 'off')
text(0.5, 0.5, spm_str_manip(H.S{1}.name, 'k60d'), 'Parent', H.nam, 'Interpreter', 'none', ...
    'FontSize', H.FS, 'HorizontalAlignment', 'center');
if nargout, Ho = H; end

%-----------------------------------------------------------------------
function Ho = select_surf(surf)
%-----------------------------------------------------------------------
global H

H.surf_sel = surf;

for ind = 1:2
    switch surf
        case 1
            H.S{ind}.info(1).Pmesh = fullfile(spm('dir'), 'toolbox', 'cat12', ['templates_surfaces' H.str32k], ...
            [H.S{ind}.info(1).side '.central.freesurfer.gii']);
        case 2
            H.S{ind}.info(1).Pmesh = fullfile(spm('dir'), 'toolbox', 'cat12', ['templates_surfaces' H.str32k], ...
            [H.S{ind}.info(1).side '.inflated.freesurfer.gii']);
        case 3
            H.S{ind}.info(1).Pmesh = fullfile(spm('dir'), 'toolbox', 'cat12', ['templates_surfaces' H.str32k], ...
            [H.S{ind}.info(1).side '.central.Template_T1_IXI555_MNI152_GS.gii']);
        case 4
            H.S{ind}.info(1).Pmesh = fullfile(spm('dir'), 'toolbox', 'cat12', ['templates_surfaces' H.str32k], [H.S{ind}.info(1).side '.patch.freesurfer.gii']);
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

%-----------------------------------------------------------------------
function display_results_all(obj, event_obj)
%-----------------------------------------------------------------------
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
    g1 = gifti(fullfile(spm('dir'), 'toolbox', 'cat12', ['templates_surfaces' H.str32k], [H.S{i}.info(1).side '.mc.freesurfer.gii']));
    g2 = gifti(fullfile(spm('dir'), 'toolbox', 'cat12', ['templates_surfaces' H.str32k], [H.S{i}.info(1).side '.sqrtsulc.freesurfer.gii']));
    H.S{i}.curv = cell(2, 1);
    H.S{i}.curv{1} = g1.cdata;
    H.S{i}.curv{2} = g2.cdata;
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
ind1 = find(H.S{1}.Y(:) ~= 0);
ind2 = find(H.S{2}.Y(:) ~= 0);

% estimate min value > 0 and min/max values
if ~isempty(ind1) && ~isempty(ind2)
    H.S{1}.thresh = min(H.S{1}.Y(H.S{1}.Y(:) > 0));
    tmp = min(H.S{2}.Y(H.S{2}.Y(:) > 0));
    if ~isempty(tmp)
      H.S{1}.thresh = min(H.S{1}.thresh, tmp);
    end
    H.S{1}.min = min(min(H.S{1}.Y(~isinf(H.S{1}.Y))), min(H.S{2}.Y(~isinf(H.S{2}.Y))));
    H.S{1}.max = max(max(H.S{1}.Y(~isinf(H.S{1}.Y))), max(H.S{2}.Y(~isinf(H.S{2}.Y))));
elseif isempty(ind1)
    H.S{1}.thresh = min(H.S{2}.Y(H.S{2}.Y(:) > 0));
    H.S{1}.min = min(H.S{2}.Y(~isinf(H.S{2}.Y)));
    H.S{1}.max = max(H.S{2}.Y(~isinf(H.S{2}.Y)));
elseif isempty(ind2)
    H.S{1}.thresh = min(H.S{1}.Y(H.S{1}.Y(:) > 0));
    H.S{1}.min = min(H.S{1}.Y(~isinf(H.S{1}.Y)));
    H.S{1}.max = max(H.S{1}.Y(~isinf(H.S{1}.Y)));
end

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

H.clim = [true H.S{1}.min H.S{1}.max];

% only apply thresholds that are slightly larger than zero
if H.S{1}.thresh > 0.00015
    H.clip = [true -H.S{1}.thresh H.S{1}.thresh];
end

for ind = 1:5
    if H.S{1}.thresh > 0.00015
        setappdata(H.patch(ind), 'clip', H.clip);
    end
    setappdata(H.patch(ind), 'clim', [true H.S{1}.min H.S{1}.max]);
    col = getappdata(H.patch(ind), 'col');
    d = getappdata(H.patch(ind), 'data');
    H = updateTexture(H, ind, d, col, H.show_transp);
end

% only show threshold popup if log-name was found and minimal value > 0 is < 1
if H.logP & (H.S{1}.thresh < 1)
    set(H.thresh, 'Enable', 'on');
end

if H.n_surf == 1
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
    mnx = ceil(2 * max(abs([H.S{1}.min H.S{1}.max])));
    
    [H.slider_min, tmp, H.str_min] = sliderPanel( ...
        'Parent', H.panel(2), ...
        'Title', 'Overlay min', ...
        'Position', H.pos{2}.ovmin, ...
        'Backgroundcolor', H.col(1,:), ...
        'Min', -mnx, ...
        'Max', mnx, ...
        'Value', H.S{1}.min, ...
        'FontName', 'Verdana', ...
        'FontSize', H.FS-1, ...
        'NumFormat', '%g', ...
        'Callback', @slider_clim_min);
    
    [H.slider_max, tmp, H.str_max] = sliderPanel( ...
        'Parent', H.panel(2), ...
        'Title', 'Overlay max', ...
        'Position', H.pos{2}.ovmax, ...
        'Backgroundcolor', H.col(1,:), ...
        'Min', -mnx, ...
        'Max', mnx, ...
        'Value', H.S{1}.max, ...
        'FontName', 'Verdana', ...
        'FontSize', H.FS-1, ...
        'NumFormat', '%g', ...
        'Callback', @slider_clim_max);
end

%-----------------------------------------------------------------------
function H = show_colorbar(H)
%-----------------------------------------------------------------------

% show colorbar
if H.n_surf == 1
   
    if isfield(H, 'cbar')
        delete(findobj('tag','cat_surf_results_colorbar'));
        H = rmfield(H, 'cbar');
    end
    
    H.cbar = axes('Parent', H.panel(1), 'Position', H.pos{1}.cbar(1, :), 'Color', H.bkg_col, 'Visible', 'off','tag','cat_surf_results_colorbar');
    H.colourbar = colorbar('peer', H.cbar, 'Northoutside');
    
    if H.logP, title(H.cbar, 'p-value', 'Color', 1 - H.bkg_col); end
    clim = getappdata(H.patch(1), 'clim');
    axis(H.cbar, 'off');
    
    if clim(3) > clim(2)
        caxis([clim(2) clim(3)]);
    end
    
    col = getappdata(H.patch(1), 'col');
    colormap(col);
    
    % Update colorbar colors if clipping is used
    clip = getappdata(H.patch(1), 'clip');
    if ~isempty(clip)
        if ~isnan(clip(2)) & ~isnan(clip(3))
            ncol = length(col);
            col_step = (clim(3) - clim(2)) / ncol;
            cmin = max([1, ceil((clip(2) - clim(2)) / col_step)]);
            cmax = min([ncol, floor((clip(3) - clim(2)) / col_step)]);
            col(cmin:cmax, :) = repmat([0.5 0.5 0.5], (cmax - cmin + 1), 1);
            colormap(col);
        end
    end
    
    if H.logP
        
        XTick = get(H.colourbar, 'XTick');
        
        % save original XTick values
        if isempty(H.XTick), H.XTick = XTick; end

        % if threshold is between 1.3..1.4 (p<0.05) change XTick accordingly and correct by 0.3
        if ~isempty(clip)
            if clip(3) >= 1.3 & clip(3) <= 1.4
                XTick_step = ceil((clim(3) - clim(2)) / 5);
                if clip(2) <= - 1.3 & clip(2) >= - 1.4
                    XTick = [(round(clim(2)) - 0.3):XTick_step: - 1.3 0 1.3:XTick_step:(round(clim(3)) + 0.3)];
                else
                    XTick = [0 1.3:XTick_step:(round(clim(3)) + 0.3)];
                end
            else
                mn = floor(min(XTick));
                mx = ceil(max(XTick));
    
                % only allow integer values
                XTick = floor(mn:mx);
%                if ~isempty(H.XTick), XTick = H.XTick; end
            end
        else
            % rescue original XThick values if clipping is changed
            if ~isempty(H.XTick), XTick = H.XTick; end
        end
        
        % change XTickLabel
        XTickLabel = [];
        for i = 1:length(XTick)
            if XTick(i) > 0
                XTickLabel = char(XTickLabel, remove_zeros(sprintf('%.g', 10^(-XTick(i)))));
            elseif XTick(i) < 0
                XTickLabel = char(XTickLabel, remove_zeros(sprintf('-%.g', 10^(XTick(i)))));
            else
                XTickLabel = char(XTickLabel, '');
            end
        end

        set(H.colourbar, 'XTickLabel', XTickLabel(2:end, :), 'XTick', XTick);
    end % end H.logP
    
    set(H.colourbar, 'XColor', 1-H.bkg_col, 'YColor', 1-H.bkg_col, 'TickDirection','out');
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
        H.cbar = axes('Parent', H.panel(1), 'Position', H.pos{1}.cbar(2, :), 'Color', H.bkg_col, 'Enable', 'off');
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
axes(H.axis);

if isfield(M, 'facevertexcdata')
    H.cdata = M.facevertexcdata;
else
    H.cdata = [];
end

if ~isfield(M, 'vertices') || ~isfield(M, 'faces')
    error('cat_surf_results:nomesh', 'ERROR:cat_surf_render: No input mesh.');
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
    curv = H.S{round(ind / 2)}.curv{H.text_mode};
else
    curv = [H.S{1}.curv{H.text_mode}; H.S{2}.curv{H.text_mode}];
end

if size(curv, 2) == 1
    
    % emphasize mean curvature values by using sqrt
    %  if H.text_mode==1
    if 1
        indneg = find(curv < 0);
        curv(indneg) = - ((-curv(indneg)).^0.5);
        indpos = find(curv > 0);
        curv(indpos) = (curv(indpos).^0.5);
        curv = curv - min(curv);
    end
    
    curv = 0.5 + repmat(curv, 1, 3);
    curv = curv / max(curv(:));
    
    % for sulcal depth (with no neg. values) use inverted values
    if H.text_mode == 2
        curv = 1 - curv;
    end
end

%-Create RGB representation of data according to colourmap
%--------------------------------------------------------------------------
C = zeros(size(v, 2), 3);
clim = getappdata(H.patch(ind), 'clim');
if isempty(clim), clim = [false NaN NaN]; end
mi = clim(2); ma = clim(3);

if any(v(:))
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
        C = C + squeeze(ind2rgb(floor(((v(i, :) - mi) / (ma - mi)) * size(col, 1)), col(:, :, i)));
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

%-Add atlas border
%--------------------------------------------------------------------------
if H.border_mode
    if ind < 5 % single hemisphere views
        A = H.S{round(ind / 2)}.A;
        A = sparse(1:size(H.S{round(ind / 2)}.M.vertices, 1), 1:size(H.S{round(ind / 2)}.M.vertices, 1), 1 ./ sum(A, 2)) * A;
        rdata = H.rdata{H.border_mode}(:, round(ind / 2));
        C0 = (A - speye(size(A))) * double(rdata);
        C(round(C0) ~= 0, :) = 0;
    else
        C0 = [];
        for i = 1:2
            A = H.S{i}.A;
            A = sparse(1:size(H.S{i}.M.vertices, 1), 1:size(H.S{i}.M.vertices, 1), 1 ./ sum(A, 2)) * A;
            rdata = H.rdata{H.border_mode}(:, i);
            C0 = [C0; (A - speye(size(A))) * double(rdata)];
        end
        C(round(C0) ~= 0, :) = 0;
    end
end

set(H.patch(ind), 'FaceVertexCData', C, 'FaceColor', FaceColor);
if isfield(H,'histax')
  cat_surf_results('hist')
  cat_surf_results('hist')
end
  
%-----------------------------------------------------------------------
function select_data(obj, event_obj, P)
%-----------------------------------------------------------------------
global H

H.logP = 1;
lh = []; rh = []; lh_rh = [];

if ~exist('P','var')
  P = spm_select([1 24], 'mesh', 'Select up to 24 maps for left and right hemisphere');
end
info = cat_surf_info(P);

n = size(P, 1);

for i = 1:n
    
    if info(i).nvertices == 64984
        H.str32k = '_32k';
        H.is32k = 1;
    else
        H.str32k = '';
        H.is32k = 0;
    end
    
    % check whether name contains 'log' that indicates a logP file
    if isempty(strfind(info(i).ff, 'log'))
        H.logP = 0;
    end
    
    % check where left and right hemisphere data were found
    if strcmp(info(i).side, 'lh')
        lh = [lh i];
    elseif strcmp(info(i).side, 'rh')
        rh = [rh i];
    elseif strcmp(info(i).side, 'mesh')
        lh_rh = [lh_rh i];
    else % asume mesh
        lh_rh = [lh_rh i];
    end
end

% merged meshes?
if ~isempty(lh_rh)
    % don't mix mesh and lh+rh data
    if ~isempty(lh) || ~isempty(rh)
        error('Mixing of left and right with merged surface data is not supported.');
    end
    
    H.merged = 1;
    H.S{1}.name = P(lh_rh, :);
    H.S{2}.name = P(lh_rh, :);
else % lh or rh meshes
    H.merged = 0;
    H.S{1}.name = P(lh, :);
    H.S{2}.name = P(rh, :);
    if numel(lh) > 1 || numel(rh) > 1
        msg = 'Warning: Display of multiple left and right surface data will not work correctly. Please only use merged hemisphere data for multiple selected surface data or just select one left and right data set.';
        h = spm('alert*', msg, '', spm('CmdLine'), 1);
        select_data;
    end
end

if isempty(H.S{2}.name)
    H.S{2}.name = fullfile(spm('dir'), 'toolbox', 'cat12', ['templates_surfaces' H.str32k], 'rh.central.freesurfer.gii');
elseif isempty(H.S{1}.name)
    H.S{1}.name = fullfile(spm('dir'), 'toolbox', 'cat12', ['templates_surfaces' H.str32k], 'lh.central.freesurfer.gii');
end

cat_surf_results('disp', H.S{1}.name, H.S{2}.name);

%==========================================================================
function save_image(obj, event_obj, filename)

global H
%%

dcm_obj = datacursormode(H.figure);

set(dcm_obj, 'Enable', 'off');

try
    delete(findall(gca, 'Type', 'hggroup', 'HandleVisibility', 'off'));
end

if ~exist('filename', 'var')
    
    nm = H.S{1}.info(1).ff;
    [pp, nm] = spm_fileparts(H.S{1}.name);
    filename = [nm '.png'];
    
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
        filename = fullfile(pth, [nam '.png']);
        newpth = pth;
    end
end

% keep background color
set(H.figure, 'InvertHardcopy', 'off', 'PaperPositionMode', 'auto');

pos = getpixelposition(H.panel(1));
hh = getframe(H.figure,pos);

img = frame2im(hh);
if H.surf_sel ~= 4
    % crop image if it's not a flatmap
    sz = size(img);
    img = img(round(0.1*sz(1):sz(1)),round(0.05*sz(2):0.95*sz(2)),:);
end

col = colormap;
imwrite(img,col,fullfile(newpth,filename));

%==========================================================================
function slider_clim_min(hObject, evt)
global H

val = get(hObject, 'Value');
c = getappdata(H.patch(1), 'clim');

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

%==========================================================================
function slider_clim_max(hObject, evt)
global H

val = get(hObject, 'Value');
c = getappdata(H.patch(1), 'clim');

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

if ~H.disable_cbar
    H = show_colorbar(H);
end

%-----------------------------------------------------------------------
function Ho = checkbox_hide_neg(obj, event_obj)
%-----------------------------------------------------------------------
global H

H.no_neg = get(H.hide_neg, 'Value');

thresh = H.thresh_value;
clip = getappdata(H.patch(1), 'clip');
clim = getappdata(H.patch(1), 'clim');

% get min value for both hemispheres
min_d = min(min(min(getappdata(H.patch(1), 'data'))), min(min(getappdata(H.patch(3), 'data'))));

if H.no_neg
    H.clip = [true -Inf thresh];
    H.clim = [true thresh clim(3)];
    set(H.slider_min, 'Value', 0);
else
    H.clip = [true -thresh thresh];
    if min_d < -thresh
        H.clim = [true -clim(3) clim(3)];
        set(H.slider_min, 'Value', -clim(3));
    end
end

for ind = 1:5
    setappdata(H.patch(ind), 'clip', H.clip);
    setappdata(H.patch(ind), 'clim', H.clim);
    col = getappdata(H.patch(ind), 'col');
    d = getappdata(H.patch(ind), 'data');
    min_d = min(min_d, min(d(:)));
    H = updateTexture(H, ind, d, col, H.show_transp);
end

% correct value of slider if no values are exceeding threshold
if min_d > -thresh & H.n_surf == 1
    set(H.slider_min, 'Value', 0);
end

set(H.atlas, 'Enable', 'on');

if ~H.disable_cbar
    H = show_colorbar(H);
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

% update colorbar
if H.n_surf == 1 & ~H.disable_cbar
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
    try, set(H.dataplot, 'XColor', 1 - H.bkg_col, 'YColor', 1 - H.bkg_col, 'Color', H.bkg_col); end
end

%==========================================================================
function checkbox_info(obj, event_obj)
global H

H.show_info = get(H.info, 'Value');

if H.show_info
    delete(findobj('tag','cat_surf_result_title'));
   %   axes('Parent', H.panel(1), 'Position',tpos{i}, 'Visible', 'off','tag','cat_surf_results_text');
    ax = axes('Parent',H.panel(1),'Position',[0.5 0.82 0.9 0.05],'visible','off','tag','cat_surf_result_title');  
    text(0,1,spm_str_manip(H.S{1}.name, 'k150d'),'HorizontalAlignment','center','interpreter','none','Color', 1 - H.bkg_col,'Parent',ax);
                    
    %set(get(getappdata(H.patch(1), 'axis'), 'Title'), 'String', ...
    %    spm_str_manip(H.S{1}.name, 'k50d'), 'Interpreter', 'none', 'Color', 1 - H.bkg_col)
    %set(get(getappdata(H.patch(3), 'axis'), 'Title'), 'String', ...
    %    spm_str_manip(H.S{2}.name, 'k50d'), 'Interpreter', 'none', 'Color', 1 - H.bkg_col)
else
    delete(findobj('tag','cat_surf_result_title'));
    set(get(getappdata(H.patch(1), 'axis'), 'Title'), 'String', '')
    set(get(getappdata(H.patch(3), 'axis'), 'Title'), 'String', '')
end

%==========================================================================
function checkbox_nocbar(obj, event_obj)
global H

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

%==========================================================================
function H = getHandles(H)
if ~nargin || isempty(H), H = gca; end
if ishandle(H) & ~isappdata(H, 'handles')
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
function select_texture(text_mode)
global H

% check that view changed
if text_mode ~= H.text_mode
    
    if text_mode == 1 % mean curvature
        H.text_mode = 1;
    else % sulcal depth
        H.text_mode = 2;
    end
    
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
    case {1, 2, 3}
        clearDataCursorPlot(H);
        set(dcm_obj, 'Enable', 'on', 'SnapToDataVertex', 'on', ...
            'DisplayStyle', 'datatip', 'Updatefcn', {@myDataCursorAtlas, H});
    case {4, 5}
        
        try
            delete(findall(gca, 'Type', 'hggroup', 'HandleVisibility', 'off'));
        end
        
        SPM_found = 1;
        for i = 1:(2-H.merged)
            SPM_name = fullfile(H.S{i}.info(1).pp, 'SPM.mat');
            
            % SPM.mat exist?
            if exist(SPM_name, 'file')
                load(SPM_name);

                % if analysis was moved we have to correct header structure
                SPM.VResMS = spm_data_hdr_read(fullfile(H.S{i}.info(1).pp,SPM.VResMS.fname));
                Vbeta = spm_data_hdr_read(fullfile(H.S{i}.info(1).pp,SPM.Vbeta(1).fname));
                for j=2:numel(SPM.Vbeta)
                  Vbeta(j) = spm_data_hdr_read(fullfile(H.S{i}.info(1).pp,SPM.Vbeta(j).fname));
                end
                SPM.Vbeta = Vbeta;
                H.SPM{i} = SPM;

                if i == 1
                    str = 'predicted, adjusted or raw values?';
                    H.predicted = spm_input(str, 1, 'b', {'predicted', 'adjusted', 'raw'}, [1 0 -1]);
                    
                    % ask for contrast for predicted or adjusted data
                    if H.predicted >= 0
                        H.Ic = spm_input('Which contrast?', 2, 'm', {SPM.xCon.name});
                    end
                end
            elseif ~isempty(H.S{i}.name)
                SPM_found = 0;
                spm('alert!', 'No SPM.mat file found.\nPlease check that you have not moved your files or your result file was moved from the folder where the SPM.mat is stored.', 1);
            end
        end
        
        if SPM_found
            set(dcm_obj, 'Enable', 'on', 'SnapToDataVertex', 'on', ...
                'DisplayStyle', 'datatip', 'Updatefcn', {@myDataCursorCluster});
            fprintf('The values are available at the MATLAB command line as variable ''y''\n');
        end
    case 6 % enable/disable rotate3d
        clearDataCursorPlot(H);
        rotate3d;
        disp('Use mouse to rotate views.');
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
    delete(findall(gca, 'Type', 'hggroup', 'HandleVisibility', 'off'));
end

%==========================================================================
function txt = myDataCursorCluster(obj, evt)
global H y

% first entries are atlases
plot_mean = H.cursor_mode - 4;
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
    if ~isnan(H.clip(2)) & ~isnan(H.clip(3))
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
    dp = d >= thresh(2); indp = find(dp);
    dn = d <= thresh(1); indn = find(dn);
    
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
    if ~isempty(indn) & isempty(found_node)
        
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
    % use single node as region
    XYZ = node;
    txt = {sprintf('Node %d', node)};
end

% for merged meshe we only have one SPM.mat with data from both hemispheres
if H.merged
    % add offset for right hemisphere
    if round(ind / 2) == 2
        if H.is32k
            XYZ = XYZ + 32492;
        else
            XYZ = XYZ + 163842;
        end
    end
    
    % always one mesh
    ind = 1;
end

[y, cbeta, CI] = get_cluster_data(H, XYZ, ind);

% if no cluster was selected set data to zero
if plot_mean & isempty(found_node)
    y(:) = 0;
    cbeta(:) = 0;
    CI(:) = 0;
end

cla(H.dataplot)
hold(H.dataplot, 'on')

set(H.dataplot, 'XColor', 1 - H.bkg_col, 'YColor', 1 - H.bkg_col,...
      'YGrid','on','Visible','on');

if H.predicted >=0
  ystr = 'contrast estimate';
  h = bar(H.dataplot, cbeta);
  set(h, 'FaceColor', H.col(1, :))
  
  % standard error
  %--------------------------------------------------------------
  CI = CI / 2;
  for j = 1:length(cbeta)
      line([j j], ([CI(j) -CI(j)] + cbeta(j)), 'LineWidth', 6, 'Color', H.col(2, :), 'Parent', H.dataplot)
  end
  set(H.dataplot, 'XLim', [0.4 (length(cbeta) + 0.6)], 'XTicklabel', '', 'XTick', [], 'YGrid','off')

else
  ystr = 'raw data';
  h = plot(H.dataplot, y);
  set(h, 'Color', H.col(2, :))
  set(H.dataplot, 'XLim', [0 (length(y))], 'XTicklabel', '', 'XTick', [])
end

xX = H.SPM{1}.xX;
iH = xX.iH;

% plot group coding for Anovas with more than 1 group and native data
if ~isempty(iH) & numel(iH)>1 & (H.predicted < 0) & isempty(xX.iB)
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

nm = H.S{1}.info(1).ff;

Ic = [];
% end with _0???.ext?
if length(nm) > 4
    if strcmp(nm(length(nm) - 4:length(nm) - 3), '_0')
        Ic = str2double(nm(length(nm) - 3:length(nm)));
    end
else
  if isfield(H,'Ic')
      Ic = H.Ic;
  end
end

if ~isempty(Ic)
    xlabel(H.dataplot, H.SPM{round(ind / 2)}.xCon(Ic).name, 'FontSize', H.FS+3, 'Color', 1 - H.bkg_col)
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
end

predicted = H.predicted;

% get raw data and whiten
y = spm_data_read(SPM.xY.VY, 'xyz', XYZ);

% for adjusted or predicted data use model estimations
%------------------------------------------------------------------
if predicted >= 0

    y = spm_filter(SPM.xX.K, SPM.xX.W * y);
    R = spm_sp('r', SPM.xX.xKXs, y);
    
    beta   = spm_data_read(SPM.Vbeta,'xyz',XYZ);
    ResMS = spm_data_read(SPM.VResMS, 'xyz', XYZ);
    
    ResMS = mean(ResMS, 2);
    Bcov = ResMS * SPM.xX.Bcov;
    
    % compute contrast of parameter estimates and 90% C.I.
    %------------------------------------------------------------------
    cbeta = SPM.xCon(Ic).c' * beta;
    cbeta = mean(cbeta, 2);
    
    CI = 1.6449; % = spm_invNcdf(1 - 0.05);
    SE = sqrt(diag(SPM.xCon(Ic).c' * Bcov * SPM.xCon(Ic).c));
    CI = CI * SE;
    
    % predicted or adjusted response
    %------------------------------------------------------------------
    if predicted == 1
        
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
else
    cbeta = [];
    CI = [];
end

y = cat_stat_nanmean(y, 2);

%==========================================================================
function txt = myDataCursorAtlas(obj, evt, H)

pos = get(evt, 'Position');

if H.cursor_mode == 1
    txt = {'Desikan DK40'};
elseif H.cursor_mode == 2
    txt = {'Destrieux 2009'};
elseif H.cursor_mode == 3
    txt = {'HCP_MMP1'};
end

i = ismember(get(H.patch(1), 'vertices'), pos, 'rows');
node = find(i);
ind = 1;

if isempty(node)
    i = ismember(get(H.patch(3), 'vertices'), pos, 'rows');
    node = find(i);
    ind = 2;
end

rdata_pos = H.rdata{H.cursor_mode}(node, ind);

rcsv = H.rcsv{H.cursor_mode};

for j = 2:size(rcsv, 1)
    if rdata_pos == rcsv{j, 1}
        txt = {txt{:} [H.S{ind}.side ' ' rcsv{j, 2}]};
        j = size(rcsv, 1);
    end
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

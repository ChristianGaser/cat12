function varargout = cat_surf_results(action, varargin)
% Visualise results for both hemispheres of surface-based analysis 
% (preferable on log P-maps). 
%
% FORMAT y = cat_surf_results('Disp',Surface)
% Surface  - a GIfTI filename/object or patch structure
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

global Hr y

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
        y               = [];
        Hr.clip         = [];
        Hr.clim         = [];
        Hr.XTick        = [];
        Hr.bkg_col      = [0 0 0];
        Hr.show_inv     = 0;
        Hr.no_neg       = 0;
        Hr.show_transp  = 1; 
        Hr.col          = [.8 .8 .8; 1 .5 .5];
        Hr.FS           = cat_get_defaults('extopts.fontsize');
        Hr.n_surf       = 1;
        Hr.thresh_value = 0;
        Hr.cursor_mode  = 1;
        Hr.text_mode    = 1;
        Hr.border_mode  = 0;
        Hr.str32k       = '';
        Hr.SPM_found    = 1;
        Hr.surf_sel     = 1;
        Hr.isfsavg      = 1;
        Hr.fixscl       = 0;
        
% colorbar addon histogram & values
% 
      
        % positions
        WS = spm('Winsize', 'Graphics');
        SS = get(0, 'Screensize');
        if 2.6 * WS(3) > SS(3)
            WS(3) = WS(3) / (2.6 * WS(3) / SS(3));
        end
        
        % result window with 5 surface views and alternative positions without top view and  only with lateral views
        Hr.viewpos ={[0.025 0.450 0.375 0.375;  0.025 0.450 0.375 0.375;  0.025 2.000 0.375 0.375],... % lh medial
                     [0.025 0.025 0.375 0.375;  0.025 0.025 0.375 0.375;  0.175 0.350 0.175 0.350],... % lh lateral
                     [0.600 0.450 0.375 0.375;  0.600 0.450 0.375 0.375;  0.600 2.000 0.375 0.375],... % rh medial
                     [0.600 0.025 0.375 0.375;  0.600 0.025 0.375 0.375;  0.675 0.350 0.175 0.350],... % rh lateral
                     [0.300 0.150 0.400 0.500;  0.300 2.000 0.400 0.500;  0.300 2.000 0.400 0.500],... % lh+rh top
                     [0.400 0.750 0.200 0.225;  0.400 0.300 0.200 0.225;  0.400 0.750 0.200 0.225]};   % data plot
        
        % change size and position of flatmaps for >= R20014b
        if spm_check_version('matlab', '8.4') >= 0
            Hr.viewpos{2}(3, :) = [-0.075 0.150 0.650 0.650]; % lh lateral
            Hr.viewpos{4}(3, :) = [0.425 0.150 0.650 0.650];  % rh lateral
        end
        
        % figure 1 with result window
        Hr.pos{1} = struct( ...
            'fig', [10 10 round(2.6*WS(3)) WS(3)], ... % figure
            'cbar', [0.400 -0.150 0.200 0.300; 0.440 0.025 0.120 0.120]);% colorbar
        
        % figure 2 with GUI
        Hr.pos{2} = struct(...
          'fig',    [2*WS(3)+10 10 0.6*WS(3) WS(3)],... 
          'sel',    [0.050 0.935 0.900 0.050],...[0.290 0.930 0.425 0.050],...
          'nam',    [0.050 0.875 0.900 0.050],...
          'surf',   [0.050 0.800 0.425 0.050],'mview',   [0.525 0.800 0.425 0.050],... 
          'text',   [0.050 0.725 0.425 0.050],'thresh',  [0.525 0.725 0.425 0.050],... 
          'cmap',   [0.050 0.650 0.425 0.050],'atlas',   [0.525 0.650 0.425 0.050],...
          'cursor', [0.050 0.575 0.425 0.050],'border',  [0.525 0.575 0.425 0.050],...
          'nocbar', [0.050 0.500 0.425 0.050],'transp',  [0.525 0.500 0.425 0.050],... 
          'info',   [0.050 0.425 0.425 0.050],'bkg',     [0.525 0.425 0.425 0.050],... 
          'inv',    [0.050 0.350 0.425 0.050],'hide_neg',[0.525 0.350 0.425 0.050],...
          'fixscl', [0.050 0.260 0.425 0.050],'scaling', [0.050 0.260 0.425 0.050],...
          'ovmin',  [0.050 0.125 0.425 0.100],'ovmax',   [0.525 0.125 0.425 0.100],... 
          'save',   [0.050 0.050 0.425 0.050],'close',   [0.525 0.050 0.425 0.050]);
        
        Hr.figure = figure(22);
        clf(Hr.figure);
    
        set(Hr.figure, 'MenuBar', 'none', 'Position', Hr.pos{1}.fig, ...
            'Name', 'CAT Results', 'NumberTitle', 'off', 'Renderer', 'OpenGL');
          
        Hr.panel(1) = uipanel('Position',[0 0 2/2.6 1],'units','normalized','BackgroundColor',...
            Hr.bkg_col,'BorderType','none'); 
        Hr.panel(2) = uipanel('Position',[2/2.6 0 0.6/2.6 1],'units','normalized','BorderType','none','BackgroundColor',Hr.col(1,:)); 
        
        % define S structure that contains information for lh and rh
        Hr.S{1}.name = ''; Hr.S{1}.side = 'lh';
        Hr.S{2}.name = ''; Hr.S{2}.side = 'rh';
       
        
        % Extra button
        if cat_get_defaults('extopts.expertgui')>1
          % Scaling (definition from cat_conf_stools)
          labels = {'SD2','SD4','SD8','%100','%99.99','min-max','0-max'};
          str = [{['Datarange ' char(133)]},labels];
          tmp = {}; 
          for il=1:numel(labels)
            tmp = [tmp {{ @(x)cat_surf_results('clims',x),labels{il} }}]; %#ok<AGROW>
          end

          Hr.scaling = uicontrol(Hr.panel(2), ...
              'String', str, 'Units', 'normalized', ...
              'Position', Hr.pos{2}.scaling, 'Userdata', tmp, ...
              'Style', 'PopUp', 'HorizontalAlignment', 'center', ...
              'Callback', 'spm(''PopUpCB'',gcbo)', ...
              'FontSize',Hr.FS,...
              'ToolTipString', 'Data range limits', ...
              'Interruptible', 'on', 'Enable', 'off');
        end
            
        % closing all windows
        Hr.close = uicontrol(Hr.panel(2), ...
            'String', 'Close', 'Units', 'normalized', ...
            'Position', Hr.pos{2}.close, ...
            'Style', 'Pushbutton', 'HorizontalAlignment', 'center', ...
            'Callback', 'close(22);', ...
            'FontSize',Hr.FS,'ForegroundColor','red',...
            'ToolTipString', 'Close windows', ...
            'Interruptible', 'on', 'Enable', 'on');

        % select results for lh and rh
        Hr.sel = uicontrol(Hr.panel(2), ...
            'String', 'Select Surface Data', 'Units', 'normalized', ...
            'Position', Hr.pos{2}.sel, ...
            'Style', 'Pushbutton', 'HorizontalAlignment', 'center', ...
            'Callback', @select_data, ...
            'FontSize',Hr.FS,...
            'ToolTipString', 'Select results (up to 24) for both hemispheres (e.g. log-p maps)', ...
            'Interruptible', 'on', 'Enable', 'on');
        
        Hr.save = uicontrol(Hr.panel(2), ...
            'String', 'Save', 'Units', 'normalized', ...
            'Position', Hr.pos{2}.save, ...
            'Style', 'Pushbutton', 'HorizontalAlignment', 'center', ...
            'Callback', {@save_image}, ...
            'FontSize',Hr.FS,...
            'ToolTipString', 'Save png image', ...
            'Interruptible', 'on', 'Enable', 'off');

        str = {'Surface', 'Central', 'Inflated', 'Dartel', 'Flatmap'};
        tmp = {{@select_surf, 1}, ...
               {@select_surf, 2}, ...
               {@select_surf, 3}, ...
               {@select_surf, 4}};
        
        % underlying surface
        Hr.surf = uicontrol(Hr.panel(2), ...
            'String', str, 'Units', 'normalized', ...
            'Position', Hr.pos{2}.surf, 'UserData', tmp, ...
            'Style', 'PopUp', 'HorizontalAlignment', 'center', ...
            'Callback', 'spm(''PopUpCB'',gcbo)', ...
            'FontSize',Hr.FS,...
            'ToolTipString', 'Underlying Surface', ...
            'Interruptible', 'on', 'Enable', 'off');
        
        str = {'Threshold', 'No threshold', 'P<0.05', 'P<0.01', 'P<0.001'};
        tmp = {{@select_thresh, 0}, ...
               {@select_thresh, 1.3}, ...
               {@select_thresh, 2}, ...
               {@select_thresh, 3}};
        
        % threshold
        Hr.thresh = uicontrol(Hr.panel(2), ...
            'String', str, 'Units', 'normalized', ...
            'Position', Hr.pos{2}.thresh, 'UserData', tmp, ...
            'Style', 'PopUp', 'HorizontalAlignment', 'center', ...
            'Callback', 'spm(''PopUpCB'',gcbo)', ...
            'FontSize',Hr.FS,...
            'ToolTipString', 'Threshold', ...
            'Interruptible', 'on', 'Enable', 'off');
        
        str = {'Colormap', 'jet', 'hot', 'hsv', 'cold-hot'};
        tmp = {{@select_cmap, 1}, ...
               {@select_cmap, 2}, ...
               {@select_cmap, 3}, ...
               {@select_cmap, 4}};
        
        % colormap
        Hr.cmap = uicontrol(Hr.panel(2), ...
            'String', str, 'Units', 'normalized', ...
            'Position', Hr.pos{2}.cmap, 'UserData', tmp, ...
            'Style', 'PopUp', 'HorizontalAlignment', 'center', ...
            'Callback', 'spm(''PopUpCB'',gcbo)', ...
            'FontSize',Hr.FS,...
            'ToolTipString', 'Threshold', ...
            'Interruptible', 'on', 'Enable', 'off');
        
        str = {'Atlas Labeling', 'Desikan-Killiany DK40', 'Destrieux 2009', 'HCP Multi-Modal Parcellation'};
        tmp = {{@select_atlas, 1}, ...
               {@select_atlas, 2}, ...
               {@select_atlas, 3}};
        
        % atlas for labeling
        Hr.atlas = uicontrol(Hr.panel(2), ...
            'String', str, 'Units', 'normalized', ...
            'Position', Hr.pos{2}.atlas, 'UserData', tmp, ...
            'Style', 'PopUp', 'HorizontalAlignment', 'center', ...
            'Callback', 'spm(''PopUpCB'',gcbo)', ...
            'FontSize',Hr.FS,...
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
        Hr.cursor = uicontrol(Hr.panel(2), ...
            'String', str, 'Units', 'normalized', ...
            'Position', Hr.pos{2}.cursor, 'UserData', tmp, ...
            'Style', 'PopUp', 'HorizontalAlignment', 'center', ...
            'Callback', 'spm(''PopUpCB'',gcbo)', ...
            'FontSize',Hr.FS,...
            'ToolTipString', 'Data Cursor Mode', ...
            'Interruptible', 'on', 'Enable', 'off');
        
        str = {'View', 'Show top view', 'Show bottom view', 'Show only lateral and medial views'};
        tmp = {{@select_view, 1}, ...
               {@select_view, -1}, ...
               {@select_view, 2}};
        
        % view
        Hr.mview = uicontrol(Hr.panel(2), ...
            'String', str, 'Units', 'normalized', ...
            'Position', Hr.pos{2}.mview, 'UserData', tmp, ...
            'Style', 'PopUp', 'HorizontalAlignment', 'center', ...
            'Callback', 'spm(''PopUpCB'',gcbo)', ...
            'FontSize',Hr.FS,...
            'ToolTipString', 'Select View', ...
            'Interruptible', 'on', 'Enable', 'off');
        
        str = {'Underlay', 'Mean curvature', 'Sulcal depth'};
        tmp = {{@select_texture, 1}, ...
               {@select_texture, 2}};
        
        % underlying texture
        Hr.text = uicontrol(Hr.panel(2), ...
            'String', str, 'Units', 'normalized', ...
            'Position', Hr.pos{2}.text, 'UserData', tmp, ...
            'Style', 'PopUp', 'HorizontalAlignment', 'center', ...
            'Callback', 'spm(''PopUpCB'',gcbo)', ...
            'FontSize',Hr.FS,...
            'ToolTipString', 'Select Underlying Texture', ...
            'Interruptible', 'on', 'Enable', 'off');
        
        str = {'Atlas Border', 'No Overlay', 'Overlay Desikan-Killiany DK40', 'Overlay Destrieux 2009', ...
               'Overlay HCP Multi-Modal Parcellation'};
        tmp = {{@select_border, 0}, ...
               {@select_border, 1}, ...
               {@select_border, 2}, ...
               {@select_border, 3}};
        
        % atlas for border overlay
        Hr.border = uicontrol(Hr.panel(2), ...
            'String', str, 'Units', 'normalized', ...
            'Position', Hr.pos{2}.border, 'UserData', tmp, ...
            'Style', 'PopUp', 'HorizontalAlignment', 'center', ...
            'Callback', 'spm(''PopUpCB'',gcbo)', ...
            'FontSize',Hr.FS,...
            'ToolTipString', 'Atlas Border Overlay', ...
            'Interruptible', 'on', 'Enable', 'off');
        
        % invert results
        Hr.inv = uicontrol(Hr.panel(2), ...
            'String', 'Invert colormap', 'Units', 'normalized', ...
            'BackgroundColor',Hr.col(1,:),...
            'Position', Hr.pos{2}.inv, ...
            'Style', 'CheckBox', 'HorizontalAlignment', 'center', ...
            'Callback', {@checkbox_inv}, ...
            'FontSize',Hr.FS,...
            'ToolTipString', 'Invert results', ...
            'Interruptible', 'on', 'Enable', 'off');
        
        % show only results for pos. contrast
        Hr.hide_neg = uicontrol(Hr.panel(2), ...
            'String', 'Hide neg. results', 'Units', 'normalized', ...
            'BackgroundColor',Hr.col(1,:),...
            'Position', Hr.pos{2}.hide_neg, ...
            'Style', 'CheckBox', 'HorizontalAlignment', 'center', ...
            'Callback', {@checkbox_hide_neg}, ...
            'FontSize',Hr.FS,...
            'ToolTipString', 'Hide neg. results', ...
            'Interruptible', 'on', 'Enable', 'off');
        
        % white background
        Hr.bkg = uicontrol(Hr.panel(2), ...
            'String', 'White background', 'Units', 'normalized', ...
            'BackgroundColor',Hr.col(1,:),...
            'Position', Hr.pos{2}.bkg, ...
            'Style', 'CheckBox', 'HorizontalAlignment', 'center', ...
            'Callback', {@checkbox_bkg}, ...
            'FontSize',Hr.FS,...
            'ToolTipString', 'White background', ...
            'Interruptible', 'on', 'Enable', 'off');
        
        % transparent view
        Hr.transp = uicontrol(Hr.panel(2), ...
            'String', 'Disable transparency', 'Units', 'normalized', ...
            'BackgroundColor',Hr.col(1,:),...
            'Position', Hr.pos{2}.transp, ...
            'Style', 'CheckBox', 'HorizontalAlignment', 'center', ...
            'Callback', {@checkbox_transp}, ...
            'FontSize',Hr.FS,...
            'ToolTipString', 'Disable transparent overlay', ...
            'Interruptible', 'on', 'Enable', 'off');
        
        Hr.info = uicontrol(Hr.panel(2), ...
            'String', 'Show filename', 'Units', 'normalized', ...
            'BackgroundColor',Hr.col(1,:),...
            'Position', Hr.pos{2}.info, ...
            'Style', 'CheckBox', 'HorizontalAlignment', 'center', ...
            'Callback', {@checkbox_info}, ...
            'FontSize',Hr.FS,...
            'ToolTipString', 'Show file information in image', ...
            'Interruptible', 'on', 'Enable', 'off');
        
        if cat_get_defaults('extopts.expertgui')<2
          Hr.nocbar = uicontrol(Hr.panel(2), ...
              'String', 'Hide colorbar', 'Units', 'normalized', ...
              'BackgroundColor',Hr.col(1,:),...
              'Position', Hr.pos{2}.nocbar, ...
              'Style', 'CheckBox', 'HorizontalAlignment', 'center', ...
              'Callback', {@checkbox_nocbar}, ...
              'FontSize',Hr.FS,...
              'ToolTipString', 'Hide colorbar', ...
              'Interruptible', 'on', 'Enable', 'off');
        else
          str = {'Colorbar', 'none', 'default', 'histogram'};
          tmp = {{@(x) cat_surf_results('colorbar',x), 0}, ...
                 {@(x) cat_surf_results('colorbar',x), 1}, ...
                 {@(x) cat_surf_results('colorbar',x), 2}};

          % colormap
          Hr.nocbar = uicontrol(Hr.panel(2), ...
              'String', str, 'Units', 'normalized', ...
              'Position', Hr.pos{2}.nocbar, 'UserData', tmp, ...
              'Style', 'PopUp', 'HorizontalAlignment', 'center', ...
              'Callback', 'spm(''PopUpCB'',gcbo)', ...
              'FontSize',Hr.FS,...
              'ToolTipString', 'Threshold', ...
              'Interruptible', 'on', 'Enable', 'off');
        end
        
        Hr.fix = uicontrol(Hr.panel(2), ...
            'String', 'Fix scaling', 'Units', 'normalized', ...
            'BackgroundColor',Hr.col(1,:),...
            'Position', Hr.pos{2}.fixscl, ...
            'Style', 'CheckBox', 'HorizontalAlignment', 'center', ...
            'Callback', {@checkbox_fixscl}, ...
            'FontSize',Hr.FS,...
            'ToolTipString', 'Fix scaling', ...
            'Interruptible', 'on', 'Visible', 'off');        
                
        if nargin >= 2
            
            if isempty(varargin{1}), return; end
            
            Hr.S{1}.name = varargin{1};
            Hr.S{2}.name = varargin{1};
            
            [pth{1}, nm1, ext1] = spm_fileparts(Hr.S{1}.name(1, :));
            
            % SPM.mat found for both hemispheres (not working yet)
            if strcmp([nm1 ext1], 'SPM.mat')
                Hr.logP = 0;
                
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
                
                Hr.S{ind}.info = cat_surf_info(Hr.S{ind}.name, 0);
                g = gifti(Hr.S{ind}.info.Pmesh);
                
                mat = v.M;
                V = g.vertices;
                XYZ = double(inv(mat) * [V'; ones(1, size(V, 1))]);
                Hr.S{ind}.Y = spm_sample_vol(Y, XYZ(1, :), XYZ(2, :), XYZ(3, :), 0)';
                Hr.S{ind}.Y = spm_mesh_project(g.vertices, dat)';
            else
                
                Hr.logP = 1;
                
                % read meshes
                Hr.S{1}.info = cat_surf_info(Hr.S{1}.name, 1);                    
                Hr.S{2}.info = Hr.S{1}.info;                    
                if Hr.S{1}.info(1).nvertices == 64984
                    Hr.str32k = '_32k';
                else
                    Hr.str32k = '';
                end

                Hr.S{1}.info(1).side = 'lh';
                Hr.S{1}.info(1).Pmesh = fullfile(spm('dir'), 'toolbox', 'cat12', ...
                        ['templates_surfaces' Hr.str32k], 'lh.central.freesurfer.gii');
                Hr.S{2}.info(1).side = 'rh';
                Hr.S{2}.info(1).Pmesh = fullfile(spm('dir'), 'toolbox', 'cat12', ....
                        ['templates_surfaces' Hr.str32k], 'rh.central.freesurfer.gii');

                for ind = 1:2
                    Hr.S{ind}.M = gifti(Hr.S{ind}.info(1).Pmesh);                    
                    % get adjacency information
                    Hr.S{ind}.A = spm_mesh_adjacency(Hr.S{ind}.M);
                end

                if isempty(strfind(Hr.S{1}.info(1).ff, 'log'))
                    Hr.logP = 0;
                end

                try
                    Y = spm_data_read(spm_data_hdr_read(Hr.S{1}.name));
                catch
                    error('No data in surfaces found or surfaces have different mesh structure (32k vs. 164k).');
                end
                Hr.nY2 = size(Y,1)/2;
                Hr.S{1}.Y = Y(1:Hr.nY2, :);
                Hr.S{2}.Y = Y((Hr.nY2+1):end, :);
                
                % if size of cdata does not fit to mesh size load the underlying
                % mesh instead of fsaverage mesh and divide hemispheres into lh/rh
                if size(Y,1) ~= (size(Hr.S{1}.M.faces,1)+4)
                  Mg = gifti(deblank(Hr.S{1}.name(1,:)));
                  sz_faces2 = size(Mg.faces,1)/2;
                  sz_vertices2 = size(Mg.vertices,1)/2;
                  Hr.S{1}.M.faces = Mg.faces(1:sz_faces2,:);
                  Hr.S{1}.M.vertices = Mg.vertices(1:sz_vertices2,:);
                  Hr.isfsavg = 0;
                end

                if ~Hr.isfsavg
                  Hr.S{2}.M.faces = Mg.faces(((sz_faces2+1):end),:) - sz_vertices2;
                  Hr.S{2}.M.vertices = Mg.vertices((sz_vertices2+1):end,:);
                end
                                      
            end
            
            % rescue original name for later result selection
            Hr.S1 = Hr.S{1};
            Hr.S2 = Hr.S{2};
            
            Hr.n_surf = numel(Hr.S{1}.info);
            Hr.view = 1;
            Hr.show_transp = 1;
            Hr.disable_cbar = 0;
            Hr.white_bgk = 0;
            Hr.show_info = 0;
            
            % result selection or RGB overlay if more than one result was loaded
            if Hr.n_surf > 1
                % pre-select 1st mesh if we cannot use RGB overlay
                %if Hr.n_surf > 3
                if 1
                    sel = 1;
                    if isempty(Hr.S1.name)
                        error('Do not mix meshes with different resolutions (i.e. 164k vs. 32k)');
                    end
                    Hr.S{1}.name = Hr.S1.name(sel, :);
                    Hr.S{2}.name = Hr.S2.name(sel, :);
                    Hr.S{1}.Y = Hr.S1.Y(:, sel);
                    Hr.S{2}.Y = Hr.S2.Y(:, sel);
                end
                
                % delete old selection ui
                delete(Hr.sel);
                H = rmfield(H, 'sel');
                
                str = cell(1, Hr.n_surf + 2);
                tmp = cell(1, Hr.n_surf + 1);
                str{1} = 'Select Result ';
                [C,C2] = spm_str_manip( Hr.S1.name , 'C');
                for s = 1:Hr.n_surf
                    str{s + 1} = spm_str_manip(Hr.S1.name(s,:), 'k60d');
                    tmp{s} = {@select_results, s};
                end
                
                % print selected filename
                Hr.nam = axes('Parent', Hr.panel(2), 'Position', Hr.pos{2}.nam);
                cla(Hr.nam);
                axis(Hr.nam, 'off')
                text(0.5, 0.5, spm_str_manip(Hr.S{1}.name, 'k60d'), 'Parent', Hr.nam, 'Interpreter', 'none', ...
                    'FontSize', Hr.FS, 'HorizontalAlignment', 'center');
                
                % set # of surfaces back to "1" if we cannot use RGB overlay
                if 1, Hr.n_surf = 1; end
                %if Hr.n_surf > 3, Hr.n_surf = 1; end
                
                % new selection ui
                str{s + 2} = 'Select new data';
                tmp{s + 1} = {@select_data};
                Hr.sel = uicontrol(Hr.panel(2), ...
                    'String', str, 'Units', 'normalized', ...
                    'Position', Hr.pos{2}.sel, 'UserData', tmp, ...
                    'Style', 'Popup', 'HorizontalAlignment', 'center', ...
                    'Callback', 'spm(''PopUpCB'',gcbo)', ...
                   'FontSize',Hr.FS,...
                    'ToolTipString', 'Select results', ...
                    'Interruptible', 'on', 'Enable', 'on');
                    
              % enable fixing of scale
              set(Hr.fix, 'Visible', 'on');
            end
            
            display_results_all;
            
            Hr.SPM_found = 1;
            SPM_name = fullfile(Hr.S{1}.info(1).pp, 'SPM.mat');
            
            % SPM.mat exist?
            if ~isempty(Hr.S{1}.name)
                Hr.SPM_found = 0;
            end

            % Don't allow plot functions for RGB maps or if SPM.mat was not found
            if Hr.n_surf > 1 && Hr.SPM_found
                str = {'Data Cursor', 'Disable data cursor', 'Atlas regions: Desikan-Killiany DK40', ...
                    'Atlas regions: Destrieux 2009', 'Atlas region: HCP Multi-Modal Parcellation', ...
                    'Enable/Disable rotate3d'};
                tmp = {{@select_cursor, 0}, ...
                       {@select_cursor, 1}, ...
                       {@select_cursor, 2}, ...
                       {@select_cursor, 3}, ...
                       {@select_cursor, 5}};
                
                Hr.cursor = uicontrol(Hr.panel(2), ...
                     'String', str, 'Units', 'normalized', ...
                    'Position', Hr.pos{2}.cursor, 'UserData', tmp, ...
                    'Style', 'PopUp', 'HorizontalAlignment', 'center', ...
                    'Callback', 'spm(''PopUpCB'',gcbo)', ...
                    'FontSize',Hr.FS,...
                    'ToolTipString', 'Data Cursor Mode', ...
                    'Interruptible', 'on', 'Enable', 'off');
            end
            
            % enable some menus only if mesh data can be assumed to be resampled
            if (length(Hr.S{1}.Y) == 32492 || length(Hr.S{1}.Y) == 163842 || length(Hr.S{1}.Y) == 40962) && Hr.isfsavg
                set(Hr.surf,   'Enable', 'on');
                set(Hr.text,   'Enable', 'on');
                set(Hr.cursor, 'Enable', 'on');
                set(Hr.border, 'Enable', 'on');
            end
            
            set(Hr.save,   'Enable', 'on');
            set(Hr.mview,  'Enable', 'on');
            set(Hr.nocbar, 'Enable', 'on');
            set(Hr.bkg,    'Enable', 'on');
            set(Hr.transp, 'Enable', 'on');
            set(Hr.info,   'Enable', 'on');
            if isfield(H,'scaling')
              set(Hr.scaling, 'Enable', 'on');
            end
            
            if min(min(Hr.S{1}.Y(:)), min(Hr.S{2}.Y(:))) < 0 & Hr.n_surf == 1
                set(Hr.inv, 'Enable', 'on');
                set(Hr.hide_neg, 'Enable', 'on');
                set(Hr.hide_neg, 'Value', 0);
            end
            
            Hr.rdata{1} = [];
            Hr.rdata{2} = [];
            Hr.rdata{3} = [];
            for ind = 1:2
                atlas_name = fullfile(spm('dir'), 'toolbox', 'cat12', ['atlases_surfaces' Hr.str32k], ...
                [Hr.S{ind}.info(1).side '.aparc_DK40.freesurfer.annot']);
                [vertices, rdata0, colortable, rcsv1] = cat_io_FreeSurfer('read_annotation', atlas_name);
                Hr.rdata{1} = [Hr.rdata{1} rdata0];
                atlas_name = fullfile(spm('dir'), 'toolbox', 'cat12', ['atlases_surfaces' Hr.str32k], ...
                [Hr.S{ind}.info(1).side '.aparc_a2009s.freesurfer.annot']);
                [vertices, rdata0, colortable, rcsv2] = cat_io_FreeSurfer('read_annotation', atlas_name);
                Hr.rdata{2} = [Hr.rdata{2} rdata0];
                atlas_name = fullfile(spm('dir'), 'toolbox', 'cat12', ['atlases_surfaces' Hr.str32k], ...
                [Hr.S{ind}.info(1).side '.aparc_HCP_MMP1.freesurfer.annot']);
                [vertices, rdata0, colortable, rcsv3] = cat_io_FreeSurfer('read_annotation', atlas_name);
                Hr.rdata{3} = [Hr.rdata{3} rdata0];
            end
            Hr.rcsv{1} = rcsv1;
            Hr.rcsv{2} = rcsv2;
            Hr.rcsv{3} = rcsv3;
            
            Hr.dcm_obj = datacursormode(Hr.figure);
            set(Hr.dcm_obj, 'Enable', 'on', 'SnapToDataVertex', 'on', ...
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
            if varargin{1} == get(Hr.nocbar, 'Value')  
                cat_surf_results('colorbar');
            end
        else
            set(Hr.nocbar, 'Value', ~get(Hr.nocbar, 'Value') );
            checkbox_nocbar;
            if ~get(Hr.nocbar, 'Value')
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
        d = getappdata(Hr.patch(1), 'data');
        col = getappdata(Hr.patch(1), 'colourmap');
        if strcmpi(varargin{1}, 'off')
            if isfield(H, 'colourbar') && ishandle(Hr.colourbar)
                delete(Hr.colourbar);
                H = rmfield(H, 'colourbar');
                setappdata(Hr.axis, 'handles', H);
            end
            return;
        end
        if isempty(d) || ~any(d(:)), varargout = {H}; return; end
        if isempty(col), col = jet(256); end
        if ~isfield(H, 'colourbar') || ~ishandle(Hr.colourbar)
            %            Hr.colourbar = colorbar('peer',gca,'NorthOutside');
            Hr.colourbar = colorbar('NorthOutside');
            set(Hr.colourbar, 'Tag', '');
            set(get(Hr.colourbar, 'Children'), 'Tag', '');
        end
        c(1:size(col, 1), 1, 1:size(col, 2)) = col;
        ic = findobj(Hr.colourbar, 'Type', 'image');
        clim = getappdata(Hr.patch(1), 'clim');
        if isempty(clim), clim = [false NaN NaN]; end
        
        if size(d, 1) > size(d, 2), d = d'; end
        
        % Update colorbar colors if clipping is used
        Hr.clip = getappdata(Hr.patch(1), 'clip');
        if ~isempty(Hr.clip)
            if ~isnan(Hr.clip(2)) && ~isnan(Hr.clip(3))
                ncol = length(col);
                col_step = (clim(3) - clim(2)) / ncol;
                cmin = max([1, ceil((Hr.clip(2) - clim(2)) / col_step)]);
                cmax = min([ncol, floor((Hr.clip(3) - clim(2)) / col_step)]);
                col(cmin:cmax, :) = repmat([0.5 0.5 0.5], (cmax - cmin + 1), 1);
                c(1:size(col, 1), 1, 1:size(col, 2)) = col;
            end
        end
        if Hr.n_surf > 1
            set(ic, 'CData', c(1:Hr.n_surf, :, :));
            set(ic, 'YData', [1 Hr.n_surf]);
            set(Hr.colourbar, 'YLim', [1 Hr.n_surf]);
            set(Hr.colourbar, 'YTickLabel', []);
        else
            set(ic, 'CData', c);
            clim = getappdata(Hr.patch(1), 'clim');
            if isempty(clim), clim = [false min(d) max(d)]; end
            set(ic, 'YData', clim(2:3));
            set(Hr.colourbar, 'YLim', clim(2:3));
        end
        setappdata(Hr.axis, 'handles', H);
        
        if nargout, varargout{1} = y; end
       %}
        
        
    %-ColourMap
    %======================================================================
    case {'colourmap', 'colormap'}
        if isempty(varargin), varargin{1} = gca; end
        if isobject(varargin{1})
            H = getHandles(varargin{1});
            if length(varargin) == 1
                varargout = {getappdata(Hr.patch(1), 'colourmap')};
                return;
            else
                setappdata(Hr.patch(1), 'colourmap', varargin{2});
                d = getappdata(Hr.patch(1), 'data');
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
            c = getappdata(Hr.patch, 'clim');
            if ~isempty(c), c = c(2:3); end
            varargout = {c};
            return;
        else
            if strcmp(varargin{2}, 'on') || isempty(varargin{2}) || any(~isfinite(varargin{2}))
                setappdata(Hr.patch, 'clim', [false NaN NaN]);
            else
                setappdata(Hr.patch, 'clim', [true varargin{2}]);
            end
            d = getappdata(Hr.patch, 'data');
            H = updateTexture(H, d);
            
        end
        
        if nargin > 1 && isnumeric(varargin{2}) && numel(varargin{2}) == 2
            caxis(Hr.axis, varargin{2});
        else
            caxis(Hr.axis, [min(d), max(d)])
        end
        
        
    %-CLip
    %======================================================================
    case 'clip'
        if isempty(varargin), varargin{1} = gca; end
        H = getHandles(varargin{1});
        if length(varargin) == 1
            c = getappdata(Hr.patch, 'clip');
            if ~isempty(c), c = c(2:3); end
            varargout = {c};
            return;
        else
            if isempty(varargin{2}) || any(~isfinite(varargin{2}))
                for ind = 1:5
                    setappdata(Hr.patch(ind), 'clip', [false NaN NaN]);
                end
            else
                for ind = 1:5
                    setappdata(Hr.patch(ind), 'clip', [true varargin{2}]);
                end
            end
            for ind = 1:5
                d = getappdata(Hr.patch, 'data');
                H = updateTexture(H, ind, d);
            end
        end
    
        
    case 'clims'
      if nargin>1, Hr.datascale    = varargin{1}; end
      if nargin>2, Hr.datascaleval = varargin{2}; end
      
      c = getappdata(Hr.patch(5), 'data');
      %%
      if isfield(H,'datascale')
        switch Hr.datascale
          case 'default'
            Hr.clim(2:3) = [max(min(c(:)),0) max(c(:))];
          case {'minmax','min-max'}
            Hr.clim(2:3) = [min(c(:)) max(c(:))];
          case {'0max','0-max'}
            Hr.clim(2:3) = [0 max(c(:))];
          otherwise
            switch Hr.datascale(1)
              case '%'
                if str2double(Hr.datascale(2:end))>0 && str2double(Hr.datascale(2:end))<100 
                  Hr.datascaleval = str2double(Hr.datascale(2:end))/100;
                end
                [cc,cx] = cat_stat_histth(c(c(:)~=0),Hr.datascaleval); %#ok<ASGLU>
                Hr.clim(2:3) = cx;
              case 'M'
                if str2double(Hr.datascale(4:end))>0  
                  Hr.datascaleval = str2double(Hr.datascale(4:end));
                end
                mv = cat_stat_nanmean(c(:)); sv = cat_stat_nanstd(c(:)); 
                Hr.clim(2:3) = [ max(min(c(:)),mv - Hr.datascaleval * sv) , min(max(c(:)), mv + Hr.datascaleval * sv)];
              case 'S'
                if str2double(Hr.datascale(3:end))>0  
                  Hr.datascaleval = str2double(Hr.datascale(3:end));
                end
                mv = cat_stat_nanmean(c(:)); sv = cat_stat_nanstd(c(:)); 
                Hr.clim(2:3) = [ mv - Hr.datascaleval * sv ,  mv + Hr.datascaleval * sv];
              otherwise
                error('unkown Hr.datascale %s.\n',Hr.datascale);
            end
        end
      end

      %% update textures of each patch
      for j=1:numel(Hr.patch) 
        setappdata(Hr.patch(j), 'clim',Hr.clim);
        H = updateTexture(H, j); 
      end
      set(Hr.str_min, 'String', sprintf('%g',Hr.clim(2)));
      set(Hr.str_max, 'String', sprintf('%g',Hr.clim(3)));
      show_colorbar(H); 
    
      
    %- print histogram
    %======================================================================
    case 'hist'
      if nargin>1
        if varargin{1}~=any(isempty(findobj('tag','cat_surf_results_hist'))) 
          cat_surf_results('hist')
        end
      else
      
        if numel(Hr.patch)>=5 && Hr.patch(1).isvalid &&  Hr.patch(3).isvalid &&  Hr.patch(5).isvalid

          if nargin>1, draw = varargin{1}; else, draw = ~any(isempty(findobj('tag','cat_surf_results_hist'))); end

          % move elements if histogram is added or removed
          if draw==0 || draw~=2
            top = Hr.patch(5).Parent; 
            pos = get(top,'Position'); 
            set(top,'Position',pos + sign(any(isempty(findobj('tag','cat_surf_results_hist')))-0.5) * [0 0.13 0 0]);
          end

          if any(isnan(Hr.clim)), cat_surf_results('clims','default'); end

          % draw/update histgram / or remove it
          mode = 0; 
          if any(isempty(findobj('tag','cat_surf_results_hist'))) || draw==2 
            if draw==2
              delete(findobj('tag','cat_surf_results_hist'));
              if isfield(H,'hist'); H = rmfield(H,'hist'); end
            end

            % print colors (red, green/dark-green
            color  = {[1 0 0],[0 1 0]/(1+Hr.bkg_col(1))};  
            linet  = {'-','--'};
            % position of the right and the left text box
            if mode
              tpos = {[0.49 0.015 0.065 0.03],[0.565 0.015 0.065 0.03],[0.4 0.015 0.12 0.03]};
            else
              tpos = {[0.4 0.015 0.12 0.06],[0.445 0.015 0.065 0.06],[0.497 0.015 0.065 0.06],[0.55 0.015 0.065 0.06],[0.60 0.015 0.065 0.06]};
            end
            % histogram axis 
            Hr.histax = axes('Parent', Hr.panel(1), 'Position', [0.4 0.102 0.20 0.15],'Visible', 'off', 'tag','cat_surf_results_hist'); 
            try
              xlim(Hr.histax,Hr.clim(2:3).*[1 1+eps]);
            catch
              disp(1);
            end
            hold on;
            % standard text
            if mode 
              Hr.dtxt(3).ax  = axes('Parent', Hr.panel(1), 'Position',tpos{3}, 'Visible', 'off','tag','cat_surf_results_text');
              Hr.dtxt(3).txt = text(0,1,sprintf('%s %s %s: \n%s - %s:\n%s:', 'mean',char(177),'std','min','max','median'),'color',[0.5 0.5 0.5],'Parent',Hr.dtxt(3).ax);
            else
              Hr.dtxt(3).ax(1) = axes('Parent', Hr.panel(1), 'Position',tpos{1}, 'Visible', 'off','tag','cat_surf_results_text');
              Hr.dtxt(3).ax(2) = axes('Parent', Hr.panel(1), 'Position',tpos{2}, 'Visible', 'off','tag','cat_surf_results_text');
              Hr.dtxt(3).ax(3) = axes('Parent', Hr.panel(1), 'Position',tpos{3}, 'Visible', 'off','tag','cat_surf_results_text');
              Hr.dtxt(3).ax(4) = axes('Parent', Hr.panel(1), 'Position',tpos{4}, 'Visible', 'off','tag','cat_surf_results_text');
              Hr.dtxt(3).ax(5) = axes('Parent', Hr.panel(1), 'Position',tpos{5}, 'Visible', 'off','tag','cat_surf_results_text');
              text(0,1,sprintf('side'),'color',[0.5 0.5 0.5],'Parent',Hr.dtxt(3).ax(1));
              text(0,1,sprintf('\n\nleft'),'color',color{1},'Parent',Hr.dtxt(3).ax(1));
              text(0,1,sprintf('\n\n\n\nright'),'color',color{2},'Parent',Hr.dtxt(3).ax(1));
              text(0,1,'min','color',[0.5 0.5 0.5],'HorizontalAlignment','center','Parent',Hr.dtxt(3).ax(2));
              text(0,1,['mean ' char(177) ' std'],'color',[0.5 0.5 0.5],'HorizontalAlignment','center','Parent',Hr.dtxt(3).ax(3));
              text(0,1,'median','color',[0.5 0.5 0.5],'HorizontalAlignment','center','Parent',Hr.dtxt(3).ax(4));
              text(0,1,'max','color',[0.5 0.5 0.5],'HorizontalAlignment','right','Parent',Hr.dtxt(3).ax(5));
            end
            for i=1:2
              side = getappdata(Hr.patch( i*2 - 1 ), 'data');
              
              if ~all(isnan(side))
                % histogram plot may fail due to NAN or whatever ...
                [d,h] = hist( side(~isinf(side(:)) & ~isnan(side(:)) & side(:)<3.4027e+38 & side(:)>-3.4027e+38 & side(:)<Hr.clim(3) & side(:)>Hr.clim(2) ), ...
                  Hr.clim(2) : diff(Hr.clim(2:3))/100 : Hr.clim(3) );
                d = d./numel(side);
                % plot histogram line and its median
                med = cat_stat_nanmedian(side(:));
                quantile = [h(find(cumsum(d)/sum(d)>0.25,1,'first')),h(find(cumsum(d)/sum(d)>0.75,1,'first'))]; 
                % print histogram
                line(Hr.histax,h,d,'color',color{i},'LineWidth',1);
                % print median
                line(Hr.histax,[med med],[0 d(find(h>=med,1,'first'))],'color',color{i},'linestyle',linet{i});
                % print quantile 
                if numel(quantile)>1
                  fill(Hr.histax,[quantile(1)   quantile(2)   quantile(2)          quantile(1)],...
                              max(d)*(0.08*[i i (i+1) (i+1)] + 0.16),color{i});
                  %
                  if mode
                    Hr.dtxt(i).ax  = axes('Parent', Hr.panel(1), 'Position',tpos{i}, 'Visible', 'off','tag','cat_surf_results_text');
                    Hr.dtxt(i).txt = text(0,1,sprintf('%10.3f %s %0.3f\n%10.3f - %0.3f\n',...
                      cat_stat_nanmean(side(:)),char(177),cat_stat_nanstd(side(:)),...
                      min(side(:)),max(side(:)),med),...
                      'color',color{i},'HorizontalAlignment','center','Parent',Hr.dtxt(i).ax);
                  else
                    text(0,1,sprintf('%s%0.3f',sprintf(repmat('\n',1,i*2)),min(side(:))),'color',color{i},'HorizontalAlignment','center','Parent',Hr.dtxt(3).ax(2));
                    text(0,1,sprintf('%s%0.3f%s%0.3f',sprintf(repmat('\n',1,i*2)),cat_stat_nanmean(side(:)),char(177),...
                      cat_stat_nanstd(side(:))),'color',color{i},'HorizontalAlignment','center','Parent',Hr.dtxt(3).ax(3));
                    text(0,1,sprintf('%s%0.3f',sprintf(repmat('\n',1,i*2)),cat_stat_nanmedian(side(:))),'color',color{i},'HorizontalAlignment','center','Parent',Hr.dtxt(3).ax(4));
                    text(0,1,sprintf('%s%0.3f',sprintf(repmat('\n',1,i*2)),max(side(:))),'color',color{i},'HorizontalAlignment','right','Parent',Hr.dtxt(3).ax(5));
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
                    if get(Hr.bkg, 'Value')==0
                        cat_surf_results('background');
                    end
                case {0,2,'black'}
                    if get(Hr.bkg, 'Value')==1
                        cat_surf_results('background');
                    end
              otherwise
                   error('Unknown background option'); 
            end
        else
            set(Hr.bkg, 'Value', ~get(Hr.bkg, 'Value') );
            checkbox_bkg;
        end
        
        
    %- set showfilename
    %======================================================================
    case 'showfilename'
        if nargin>1
            if varargin{1} ~= get(Hr.bkg, 'Value')==0  
                cat_surf_results('showfilename');
            end
        else
            set(Hr.info, 'Value', ~get(Hr.info, 'Value') );
            checkbox_info;
        end
        
        
    %- set transparency
    %======================================================================
    case 'transparency'
        if nargin>1
            if varargin{1} ~= get(Hr.transp, 'Value')==0  
                cat_surf_results('transparency');
            end
        else
            set(Hr.transp, 'Value', ~get(Hr.transp, 'Value') );
            checkbox_transp;
        end 
        
        
    %- set inverse colormap
    %======================================================================
    case 'invcolormap'
        if nargin>1
            if varargin{1} ~= get(Hr.inv, 'Value')==0  
                cat_surf_results('invcolormap');
            end
        else
            set(Hr.inv, 'Value', ~get(Hr.inv, 'Value') );
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
            if varargin{1} ~= get(Hr.hide_neg, 'Value')==0  
                cat_surf_results('hide_neg');
            end
        else
            set(Hr.hide_neg, 'Value', ~get(Hr.hide_neg, 'Value') );
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
        clear -globalvar Hr; 
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
        maximg = numel(Hr.S1.info);
        if isinf(imgs), imgs = 1:maximg; end 
        imgs(imgs<0 | imgs>maximg) = []; 
        
        %% print images
        for fi=1:numel(imgs)
          if fi>1, select_results(imgs(fi)); end
          
          [pp,ff] = spm_fileparts( Hr.S1.name(imgs(fi),:) );
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
global Hr

Hr.thresh_value = thresh;
Hr.clip = [true -thresh thresh];

Hr.no_neg = get(Hr.hide_neg, 'Value');

% get min value for both hemispheres
min_d = min(min(min(getappdata(Hr.patch(1), 'data'))), min(min(getappdata(Hr.patch(3), 'data'))));
clim = getappdata(Hr.patch(1), 'clim');

% rather use NaN values for zero threshold
if thresh == 0
    Hr.clip = [false NaN NaN];
end

if Hr.no_neg
    Hr.clip = [true -Inf thresh];
    clim = [true 0 clim(3)];
    set(Hr.slider_min, 'Value', 0);
end

for ind = 1:5
    if min_d > -thresh
        setappdata(Hr.patch(ind), 'clim', [true thresh clim(3)]);
    elseif thresh == 0
        setappdata(Hr.patch(ind), 'clim', [true -clim(3) clim(3)]);
    end
    
    setappdata(Hr.patch(ind), 'clip', Hr.clip);
    col = getappdata(Hr.patch(ind), 'col');
    d = getappdata(Hr.patch(ind), 'data');
    min_d = min(min_d, min(d(:)));
    H = updateTexture(H, ind, d, col, Hr.show_transp);
end

set(Hr.slider_min, 'Value', Hr.clim(2))
set(Hr.str_min, 'String', sprintf('%g',Hr.clim(2)));

set(Hr.atlas, 'Enable', 'on');

if ~Hr.disable_cbar
    H = show_colorbar(H);
end
if nargout, Ho = H; end

%-----------------------------------------------------------------------
function disphist(type)
%-----------------------------------------------------------------------
  global Hr; 
  
  i = get(Hr.sel,'value');
  if isfield(H,'patch')
    if i==0
      d = getappdata(Hr.patch(1),'data');
    elseif i<=numel(Hr.patch)
      d = getappdata(Hr.patch(i),'data');
    end
    cat_stat_histth(d(d(:)~=0),1,type); 
  end


%-----------------------------------------------------------------------
function Ho = select_cmap(cmap)
%-----------------------------------------------------------------------
global Hr

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
    setappdata(Hr.patch(ind), 'col', col);
    d = getappdata(Hr.patch(ind), 'data');
    H = updateTexture(H, ind, d, col, Hr.show_transp);
end

if ~Hr.disable_cbar
    H = show_colorbar(H);
end
if nargout, Ho = H; end

%-----------------------------------------------------------------------
function Ho = select_atlas(atlas)
%-----------------------------------------------------------------------
global Hr

% get threshold from clipping
thresh = [0 0];
if ~isempty(Hr.clip)
    if ~isnan(Hr.clip(2)) & ~isnan(Hr.clip(3))
        thresh = [Hr.clip(2:3)];
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
    rcsv = Hr.rcsv{atlas};
    rdata = Hr.rdata{atlas}(:, round(ind / 2));
    
    A = Hr.S{round(ind / 2)}.A;
    A = A + speye(size(A));
    d0 = getappdata(Hr.patch(ind), 'data');
    
    % go through all surfaces
    for indsurf = 1:Hr.n_surf
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
            fprintf('%s: Positive effects in %s', atlas_name, Hr.S{round(ind / 2)}.info(1).side);
            fprintf('\n%s', spm_str_manip(Hr.S{round(ind / 2)}.info(indsurf).fname, 'k50d'));
            fprintf('\n______________________________________________________\n\n');
            
            if Hr.logP, fprintf('%7s\t%8s\t%s\n', 'P-value', 'Size', 'Overlap of atlas region');
            else, fprintf('%7s\t%8s\t%s\n', 'Value  ', 'Size', 'Overlap of atlas region'); end
            
            for i = 1:max(C)
                N = find(C == i);
                k = length(N);
                
                dmax = d(indp); dmax = max(dmax(N));
                
                if Hr.logP, fprintf('\n%1.5f\t%8d', 10^(-dmax), k);
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
            fprintf('%s: Negative effects in %s', atlas_name, Hr.S{round(ind / 2)}.info(1).side);
            fprintf('\n%s', spm_str_manip(Hr.S{round(ind / 2)}.info(indsurf).fname, 'k50d'));
            fprintf('\n______________________________________________________\n\n');
            
            if Hr.logP, fprintf('%7s\t%8s\t%s\n', 'P-value', 'Size', 'Overlap of atlas region');
            else, fprintf('%7s\t%8s\t%s\n', 'Value  ', 'Size', 'Overlap of atlas region'); end
            
            for i = 1:max(C)
                N = find(C == i);
                k = length(N);
                
                dmin = d(indn); dmin = min(dmin(N));
                if Hr.logP, fprintf('\n%1.5f\t%8d', 10^(dmin), k);
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
global Hr

clearDataCursorPlot(H);

Hr.S{1}.name = Hr.S1.name(sel, :);
Hr.S{2}.name = Hr.S2.name(sel, :);
Hr.S{1}.Y = Hr.S1.Y(:, sel);
Hr.S{2}.Y = Hr.S2.Y(:, sel);

% check whether data for left or right hemipshere are all non-zero
ind1 = find(Hr.S{1}.Y(:) ~= 0);
ind2 = find(Hr.S{2}.Y(:) ~= 0);

% estimate min value > 0 and min/max values
if ~isempty(ind1) & ~isempty(ind2)
    Hr.S{1}.thresh = min(Hr.S{1}.Y(Hr.S{1}.Y(:) > 0));
    Hr.S{1}.thresh = min(Hr.S{1}.thresh, min(Hr.S{2}.Y(Hr.S{2}.Y(:) > 0)));
    Hr.S{1}.min = min(min(Hr.S{1}.Y(~isinf(Hr.S{1}.Y))), min(Hr.S{2}.Y(~isinf(Hr.S{2}.Y))));
    Hr.S{1}.max = max(max(Hr.S{1}.Y(~isinf(Hr.S{1}.Y))), max(Hr.S{2}.Y(~isinf(Hr.S{2}.Y))));
elseif isempty(ind1)
    Hr.S{1}.thresh = min(Hr.S{2}.Y(Hr.S{2}.Y(:) > 0));
    Hr.S{1}.min = min(Hr.S{2}.Y(~isinf(Hr.S{2}.Y)));
    Hr.S{1}.max = max(Hr.S{2}.Y(~isinf(Hr.S{2}.Y)));
elseif isempty(ind2)
    Hr.S{1}.thresh = min(Hr.S{1}.Y(Hr.S{1}.Y(:) > 0));
    Hr.S{1}.min = min(Hr.S{1}.Y(~isinf(Hr.S{1}.Y)));
    Hr.S{1}.max = max(Hr.S{1}.Y(~isinf(Hr.S{1}.Y)));
end

mn = Hr.S{1}.min;

% deal with neg. values
if Hr.S{1}.min < 0
    mnx = max(abs([Hr.S{1}.min, Hr.S{1}.max]));
    Hr.S{1}.min = - mnx;
    Hr.S{1}.max = mnx;
end

% add 10% to min/max values
Hr.S{1}.max = round(1100 * Hr.S{1}.max) / 1000;
if Hr.S{1}.min < 0
    Hr.S{1}.min = round(1100 * Hr.S{1}.min) / 1000;
else
    Hr.S{1}.min = round(900 * Hr.S{1}.min) / 1000;
end

% correct lower clim to "0" if no values are exceeding threshold
if mn > -Hr.thresh_value
    Hr.clim = [true Hr.thresh_value Hr.S{1}.max];
else
    Hr.clim = [true Hr.S{1}.min Hr.S{1}.max];
end

% only apply thresholds that are slightly larger than zero
if Hr.S{1}.thresh > 0.00015 & Hr.thresh_value == 0
    Hr.clip = [true -Hr.S{1}.thresh Hr.S{1}.thresh];
else
    Hr.clip = [true -Hr.thresh_value Hr.thresh_value];
end

% rather use NaN values for zero threshold
if Hr.thresh == 0
    Hr.clip = [false NaN NaN];
end

if Hr.no_neg
    Hr.clip = [true -Inf Hr.clip(3)];
    set(Hr.slider_min, 'Value', 0);
end

Hr.n_surf = 1;

for ind = 1:5
    if Hr.S{1}.thresh > 0.00015
        setappdata(Hr.patch(ind), 'clip', Hr.clip);
    end
    
    % update clim only for non-fixed scaling
    if ~Hr.fixscl
      setappdata(Hr.patch(ind), 'clim', Hr.clim);
    end
    col = getappdata(Hr.patch(ind), 'col');
    
    if ind > 4
        d = [Hr.S{1}.Y; Hr.S{2}.Y];
    else
        d = Hr.S{round(ind / 2)}.Y;
    end
    
    H = updateTexture(H, ind, d, col, Hr.show_transp);
end

% correct value of slider if no values are exceeding threshold
if Hr.S{1}.min > - Hr.thresh_value
    set(Hr.slider_min, 'Value', 0);
end

% update sliders for non-fixed scaling
if ~Hr.fixscl
  set(Hr.slider_min, 'Value', Hr.clim(2));
  set(Hr.slider_max, 'Value', Hr.clim(3));
  set(Hr.str_min, 'String', sprintf('%g',Hr.clim(2)));
  set(Hr.str_max, 'String', sprintf('%g',Hr.clim(3)));
end

% update file information and colorbar
checkbox_info;

if ~Hr.disable_cbar
    H = show_colorbar(H);
end

% print selected filename
cla(Hr.nam);
axis(Hr.nam, 'off')
text(0.5, 0.5, spm_str_manip(Hr.S{1}.name, 'k60d'), 'Parent', Hr.nam, 'Interpreter', 'none', ...
    'FontSize', Hr.FS, 'HorizontalAlignment', 'center');
if nargout, Ho = H; end

%-----------------------------------------------------------------------
function Ho = select_surf(surf)
%-----------------------------------------------------------------------
global Hr

Hr.surf_sel = surf;

for ind = 1:2
    switch surf
        case 1
            Hr.S{ind}.info(1).Pmesh = fullfile(spm('dir'), 'toolbox', 'cat12', ['templates_surfaces' Hr.str32k], ...
            [Hr.S{ind}.info(1).side '.central.freesurfer.gii']);
        case 2
            Hr.S{ind}.info(1).Pmesh = fullfile(spm('dir'), 'toolbox', 'cat12', ['templates_surfaces' Hr.str32k], ...
            [Hr.S{ind}.info(1).side '.inflated.freesurfer.gii']);
        case 3
            Hr.S{ind}.info(1).Pmesh = fullfile(spm('dir'), 'toolbox', 'cat12', ['templates_surfaces' Hr.str32k], ...
            [Hr.S{ind}.info(1).side '.central.Template_T1_IXI555_MNI152_GS.gii']);
        case 4
            Hr.S{ind}.info(1).Pmesh = fullfile(spm('dir'), 'toolbox', 'cat12', ['templates_surfaces' Hr.str32k], [Hr.S{ind}.info(1).side '.patch.freesurfer.gii']);
    end
    Hr.S{ind}.M = gifti(Hr.S{ind}.info(1).Pmesh);
end

for ind = 1:5
    if ind < 5 % single hemisphere views
        M = Hr.S{round(ind / 2)}.M;
    else
        M.faces = [Hr.S{1}.M.faces; Hr.S{2}.M.faces + size(Hr.S{1}.M.vertices, 1)];
        M.vertices = [Hr.S{1}.M.vertices; Hr.S{2}.M.vertices];
        M.mat = Hr.S{1}.M.mat;
    end
    
    set(Hr.patch(ind), 'Vertices', M.vertices);
    set(Hr.patch(ind), 'Faces', M.faces);
    
    % rescale axis except for flatmaps
    Ha = getappdata(Hr.patch(ind), 'axis');
    axes(Ha);
    
    if surf < 4
        axis(Ha, 'image');
        axis(Ha, 'off');
    end
end

% only show lateral views for flatmaps
if surf == 4
    select_view(3)
elseif Hr.view == 3
    select_view(1)
end

% don't show data cursor, view functions and data plot that will not work for flatmaps
if surf == 4
    set(Hr.cursor, 'Enable', 'off');
    set(Hr.mview, 'Enable', 'off');
    clearDataCursorPlot(H);
else
    set(Hr.cursor, 'Enable', 'on');
    set(Hr.mview, 'Enable', 'on');
end
if nargout, Ho = H; end

%-----------------------------------------------------------------------
function display_results_all(obj, event_obj)
%-----------------------------------------------------------------------
global Hr

if (size(Hr.S{1}.Y) > 1 | size(Hr.S{2}.Y) > 1) & min(min(Hr.S{1}.Y(:)), min(Hr.S{2}.Y(:))) < 0
    disp('Warning: Only results with positive values are displayed!');
end

% clear larger area and set background color to update labels and title
Hr.Ha = axes('Parent', Hr.panel(1), 'Position', [-.1 -.1 1.1 1.1], 'Color', Hr.bkg_col);
cla(Hr.Ha);

Hr.renderer = get(Hr.figure, 'Renderer');
set(Hr.figure, 'Renderer', 'OpenGL');

%-Get mesh curvature and sulcal depth
%------------------------------------------------------------------
for i = 1:2
    g1 = gifti(fullfile(spm('dir'), 'toolbox', 'cat12', ['templates_surfaces' Hr.str32k], [Hr.S{i}.info(1).side '.mc.freesurfer.gii']));
    g2 = gifti(fullfile(spm('dir'), 'toolbox', 'cat12', ['templates_surfaces' Hr.str32k], [Hr.S{i}.info(1).side '.sqrtsulc.freesurfer.gii']));
    Hr.S{i}.curv = cell(2, 1);
    Hr.S{i}.curv{1} = g1.cdata;
    Hr.S{i}.curv{2} = g2.cdata;
end

if Hr.view == 1 % top view
    vv = [90 0; -90 0; -90 0; 90 0; 0 90];
else % bottom view
    vv = [90 0; -90 0; -90 0; 90 0; 0 -90];
end

for ind = 1:5
    display_results(ind, Hr.viewpos{ind}(abs(Hr.view), :), vv(ind, :));
end

% prepare dataplot axes
Hr.dataplot = axes('Position', Hr.viewpos{6}(abs(Hr.view), :), 'Parent', Hr.panel(1), 'Color', Hr.bkg_col);
Hr.figure = ancestor(Hr.dataplot, 'figure');
try axes(Hr.dataplot); end
axis off

% check whether data for left or right hemipshere are all non-zero
ind1 = find(Hr.S{1}.Y(:) ~= 0);
ind2 = find(Hr.S{2}.Y(:) ~= 0);

% estimate min value > 0 and min/max values
if ~isempty(ind1) && ~isempty(ind2)
    Hr.S{1}.thresh = min(Hr.S{1}.Y(Hr.S{1}.Y(:) > 0));
    tmp = min(Hr.S{2}.Y(Hr.S{2}.Y(:) > 0));
    if ~isempty(tmp)
      Hr.S{1}.thresh = min(Hr.S{1}.thresh, tmp);
    end
    Hr.S{1}.min = min(min(Hr.S{1}.Y(~isinf(Hr.S{1}.Y))), min(Hr.S{2}.Y(~isinf(Hr.S{2}.Y))));
    Hr.S{1}.max = max(max(Hr.S{1}.Y(~isinf(Hr.S{1}.Y))), max(Hr.S{2}.Y(~isinf(Hr.S{2}.Y))));
elseif isempty(ind1)
    Hr.S{1}.thresh = min(Hr.S{2}.Y(Hr.S{2}.Y(:) > 0));
    Hr.S{1}.min = min(Hr.S{2}.Y(~isinf(Hr.S{2}.Y)));
    Hr.S{1}.max = max(Hr.S{2}.Y(~isinf(Hr.S{2}.Y)));
elseif isempty(ind2)
    Hr.S{1}.thresh = min(Hr.S{1}.Y(Hr.S{1}.Y(:) > 0));
    Hr.S{1}.min = min(Hr.S{1}.Y(~isinf(Hr.S{1}.Y)));
    Hr.S{1}.max = max(Hr.S{1}.Y(~isinf(Hr.S{1}.Y)));
end

% deal with neg. values
if Hr.S{1}.min < 0
    mnx = max(abs([Hr.S{1}.min, Hr.S{1}.max]));
    Hr.S{1}.min = - mnx;
    Hr.S{1}.max = mnx;
end

% add 10% to min/max values
Hr.S{1}.max = round(1100 * Hr.S{1}.max) / 1000;
if Hr.S{1}.min < 0
    Hr.S{1}.min = round(1100 * Hr.S{1}.min) / 1000;
else
    Hr.S{1}.min = round(900 * Hr.S{1}.min) / 1000;
end

Hr.clim = [true Hr.S{1}.min Hr.S{1}.max];

% only apply thresholds that are slightly larger than zero
if Hr.S{1}.thresh > 0.00015
    Hr.clip = [true -Hr.S{1}.thresh Hr.S{1}.thresh];
end

for ind = 1:5
    if Hr.S{1}.thresh > 0.00015
        setappdata(Hr.patch(ind), 'clip', Hr.clip);
    end
    setappdata(Hr.patch(ind), 'clim', [true Hr.S{1}.min Hr.S{1}.max]);
    col = getappdata(Hr.patch(ind), 'col');
    d = getappdata(Hr.patch(ind), 'data');
    H = updateTexture(H, ind, d, col, Hr.show_transp);
end

% only show threshold popup if log-name was found and minimal value > 0 is < 1
if Hr.logP & (Hr.S{1}.thresh < 1)
    set(Hr.thresh, 'Enable', 'on');
end

if Hr.n_surf == 1
    % get sure that image is thresholded and there are at least 20% zero/NaN areas
    if (sum(d ~= 0) / numel(d) < 0.8)
        set(Hr.atlas, 'Enable', 'on');
    end
end

if ~Hr.disable_cbar
    H = show_colorbar(H);
end

% show slider for range of results
if Hr.n_surf == 1
    
    % allow slider a more extended range
    mnx = ceil(2 * max(abs([Hr.S{1}.min Hr.S{1}.max])));
    
    [Hr.slider_min, tmp, Hr.str_min] = sliderPanel( ...
        'Parent', Hr.panel(2), ...
        'Title', 'Overlay min', ...
        'Position', Hr.pos{2}.ovmin, ...
        'Backgroundcolor', Hr.col(1,:), ...
        'Min', -mnx, ...
        'Max', mnx, ...
        'Value', Hr.S{1}.min, ...
        'FontName', 'Verdana', ...
        'FontSize', Hr.FS-1, ...
        'NumFormat', '%g', ...
        'Callback', @slider_clim_min);
    
    [Hr.slider_max, tmp, Hr.str_max] = sliderPanel( ...
        'Parent', Hr.panel(2), ...
        'Title', 'Overlay max', ...
        'Position', Hr.pos{2}.ovmax, ...
        'Backgroundcolor', Hr.col(1,:), ...
        'Min', -mnx, ...
        'Max', mnx, ...
        'Value', Hr.S{1}.max, ...
        'FontName', 'Verdana', ...
        'FontSize', Hr.FS-1, ...
        'NumFormat', '%g', ...
        'Callback', @slider_clim_max);
end

%-----------------------------------------------------------------------
function H = show_colorbar(H)
%-----------------------------------------------------------------------

% show colorbar
if Hr.n_surf == 1
   
    if isfield(H, 'cbar')
        delete(findobj('tag','cat_surf_results_colorbar'));
        H = rmfield(H, 'cbar');
    end
    
    Hr.cbar = axes('Parent', Hr.panel(1), 'Position', Hr.pos{1}.cbar(1, :), 'Color', Hr.bkg_col, 'Visible', 'off','tag','cat_surf_results_colorbar');
    Hr.colourbar = colorbar('peer', Hr.cbar, 'Northoutside');
    
    if Hr.logP, title(Hr.cbar, 'p-value', 'Color', 1 - Hr.bkg_col); end
    clim = getappdata(Hr.patch(1), 'clim');
    axis(Hr.cbar, 'off');
    
    if clim(3) > clim(2)
        caxis([clim(2) clim(3)]);
    end
    
    col = getappdata(Hr.patch(1), 'col');
    colormap(col);
    
    % Update colorbar colors if clipping is used
    clip = getappdata(Hr.patch(1), 'clip');
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
    
    if Hr.logP
        
        XTick = get(Hr.colourbar, 'XTick');
        
        % save original XTick values
        if isempty(Hr.XTick), Hr.XTick = XTick; end

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
%                if ~isempty(Hr.XTick), XTick = Hr.XTick; end
            end
        else
            % rescue original XThick values if clipping is changed
            if ~isempty(Hr.XTick), XTick = Hr.XTick; end
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

        set(Hr.colourbar, 'XTickLabel', XTickLabel(2:end, :), 'XTick', XTick);
    end % end Hr.logP
    
    set(Hr.colourbar, 'XColor', 1-Hr.bkg_col, 'YColor', 1-Hr.bkg_col, 'TickDirection','out');
    try
      set(Hr.colourbar, 'TickLength', 0.01);
    catch
      set(Hr.colourbar, 'TickLength', [0 0]);
    end
    
    %{
    if isfield(H,'hist')
      cat_surf_results('hist')
      cat_surf_results('hist')
    end
    %}
else
    delete(findobj('tag','cat_surf_results_hist'));
    
    if ~isfield(H, 'cbar') || ~ishandle(Hr.cbar)
        Hr.cbar = axes('Parent', Hr.panel(1), 'Position', Hr.pos{1}.cbar(2, :), 'Color', Hr.bkg_col, 'Enable', 'off');
    end
    
    % RGB colorbar
    if Hr.n_surf == 3
        cb = [8 1 1 4 2 2 8; ...
              8 1 6 7 5 2 8; ...
              8 8 3 3 3 8 8];
    else %RG colorbar
        cb = [8 1 1 4 2 2 8; ...
              8 1 1 4 2 2 8];
    end
    imagesc(cb);
    colormap([1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 1 1 1; Hr.bkg_col]);
    axis(Hr.cbar, 'off'); axis('image');
end

%-----------------------------------------------------------------------
function display_results(ind, win, vw)
%-----------------------------------------------------------------------
global Hr

% rescue old color before a new Hr.patch is created
try
    col = getappdata(Hr.patch(ind), 'col');
catch
    col = [];
end

if ind < 5 % single hemisphere views
    M = Hr.S{round(ind / 2)}.M;
    Mc.cdata = Hr.S{round(ind / 2)}.Y;
else
    Ml = Hr.S{1}.M;
    Mr = Hr.S{2}.M;
    Mcl.cdata = Hr.S{1}.Y;
    Mcr.cdata = Hr.S{2}.Y;
    
    % check whether number of data for lh/rh differ and fill with zeros
    diff_size_Y = size(Hr.S{1}.Y, 2) - size(Hr.S{2}.Y, 2);
    if diff_size_Y > 0
        Mcr.cdata = [Mcr.cdata zeros(size(Hr.S{2}.Y, 1), 1)];
    end
    if diff_size_Y < 0
        Mcl.cdata = [Mcl.cdata; zeros(size(Hr.S{1}.Y, 1), 1)];
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

Hr.axis = axes('Position', win, 'Parent', Hr.panel(1), 'Visible', 'off');
Hr.figure = ancestor(Hr.axis, 'figure');
%axes(Hr.axis);

if isfield(M, 'facevertexcdata')
    Hr.cdata = M.facevertexcdata;
else
    Hr.cdata = [];
end

if ~isfield(M, 'vertices') || ~isfield(M, 'faces')
    error('cat_surf_results:nomesh', 'ERROR:cat_surf_render: No input mesh.');
end

%% -Patch
%------------------------------------------------------------------
P = struct('vertices', M.vertices, 'faces', double(M.faces));
Hr.patch(ind) = patch(P, ...
    'FaceColor', [0.6 0.6 0.6], ...
    'EdgeColor', 'none', ...
    'FaceLighting', 'gouraud', ...
    'SpecularStrength', 0.1, ...
    'AmbientStrength', 1.0, ...
    'DiffuseStrength', 0.6, ...
    'SpecularExponent', 15, ...
    'Clipping', 'off', ...
    'DeleteFcn', {@myDeleteFcn, Hr.renderer}, ...
    'Visible', 'off', ...
    'Tag', 'CATSurfRender', ...
    'Parent', Hr.axis);
setappdata(Hr.patch(ind), 'patch', P);
setappdata(Hr.patch(ind), 'axis', Hr.axis);

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

axis(Hr.axis, 'image');
axis(Hr.axis, 'off');
view(Hr.axis, vw);
material(Hr.figure, 'dull');

% default lighting
Hr.light(1) = camlight('headlight'); set(Hr.light(1), 'Parent', Hr.axis);
setappdata(Hr.axis, 'handles', H);
set(Hr.patch(ind), 'Visible', 'on');
camlight(Hr.light(1),'headlight')

%==========================================================================
function [H, C] = updateTexture(H, ind, v, col, transp)

%-Project data onto surface mesh
%--------------------------------------------------------------------------
if nargin<3
  v = getappdata(Hr.patch(ind), 'data');
end
if nargin<5
  transp = Hr.transp;
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

setappdata(Hr.patch(ind), 'data', v);
setappdata(Hr.patch(ind), 'col', col);

if ~exist('FaceColor', 'var') || isempty(FaceColor), FaceColor = 'interp'; end

%-Get curvature
%--------------------------------------------------------------------------
if ind < 5 % single hemisphere views
    curv = Hr.S{round(ind / 2)}.curv{Hr.text_mode};
else
    curv = [Hr.S{1}.curv{Hr.text_mode}; Hr.S{2}.curv{Hr.text_mode}];
end

if size(curv, 2) == 1
    
    % emphasize mean curvature values by using sqrt
    %  if Hr.text_mode==1
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
    if Hr.text_mode == 2
        curv = 1 - curv;
    end
end

%-Create RGB representation of data according to colourmap
%--------------------------------------------------------------------------
C = zeros(size(v, 2), 3);
clim = getappdata(Hr.patch(ind), 'clim');
if isempty(clim), clim = [false NaN NaN]; end
mi = clim(2); ma = clim(3);

if any(v(:))
    if ~clim(1), mi = min(v(:)); ma = max(v(:)); end
    % don't allow negative values for multiple maps
    if size(v, 1) > 1 && mi < 0
        if ~isempty(Hr.clip)
            Hr.clip(2) = - Inf;
        else
            Hr.clip = [true -Inf 0];
        end
    end
    for i = 1:size(v, 1)
        C = C + squeeze(ind2rgb(floor(((v(i, :) - mi) / (ma - mi)) * size(col, 1)), col(:, :, i)));
    end
end

if ~isempty(Hr.clip)
    v(v > Hr.clip(2) & v < Hr.clip(3)) = NaN;
    setappdata(Hr.patch(ind), 'clip', [true Hr.clip(2) Hr.clip(3)]);
end

setappdata(Hr.patch(ind), 'clim', [true mi ma]);
Hr.clim = [true mi ma];

  
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

%-Add/delete atlas border
%--------------------------------------------------------------------------

if ~Hr.border_mode
  for k=1:5
    try
      h3 = getappdata(Hr.patch(k), 'h3');
      if ~isempty(h3)
        for i=1:size(h3,1)
          delete(h3(i))
        end
      end
    end
  end
end

set(Hr.patch(ind), 'FaceVertexCData', C, 'FaceColor', FaceColor);

if Hr.border_mode
  if ind == 1
    for k=1:2
      rdata = Hr.rdata{Hr.border_mode}(:, k);
      datarange = 0:max(rdata(:));
      Hi = hist(rdata(:),datarange);
      indo = [1 find(Hi==0)];
      datarange(indo) = [];
  
      t = datarange;
      M = Hr.S{k}.M;
      Hr.S{k}.Cm = cell(numel(t),1);
      for i=1:numel(t)
        T = zeros(size(rdata));
        T(rdata == datarange(i)) = 1;
        Cm = spm_mesh_contour(M,struct('T',T,'t',0.5));
        Hr.S{k}.Cm{i} = Cm;
      end
    end
    for k=1:5
      h3 = getappdata(Hr.patch(k), 'h3');
      if ~isempty(h3)
        for i=1:size(h3,1)
          delete(h3(i))
        end
      end
      setappdata(Hr.patch(k), 'h3', []);
    end
  end
  
  Ha = getappdata(Hr.patch(ind), 'axis');
  hold on

  if ind < 5 % single hemisphere views
    Cm = Hr.S{round(ind / 2)}.Cm;
    for j=1:size(Cm,1)
      if ~isempty(Cm{j})
				for i=1:size(Cm,2)
					h3 = plot3(Cm{j}(i).xdata,Cm{j}(i).ydata,Cm{j}(i).zdata,'k-','LineWidth',2);
					setappdata(Hr.patch(ind), 'h3', [getappdata(Hr.patch(ind), 'h3'); h3]);
					set(h3,'Parent',Ha);
				end
		  end
    end
  else
    for k = 1:2
      Cm = Hr.S{k}.Cm;
      for j=1:size(Cm,1)
        if ~isempty(Cm{j})
					for i=1:size(Cm,2)
						h3 = plot3(Cm{j}(i).xdata,Cm{j}(i).ydata,Cm{j}(i).zdata,'k-','LineWidth',2);
						setappdata(Hr.patch(ind), 'h3', [getappdata(Hr.patch(ind), 'h3'); h3]);
						set(h3,'Parent',Ha);
					end
				end
      end
    end
  end
  
  hold off

end

if isfield(H,'histax')
  cat_surf_results('hist')
  cat_surf_results('hist')
end
  
%-----------------------------------------------------------------------
function select_data(obj, event_obj, P)
%-----------------------------------------------------------------------
global Hr

Hr.logP = 1;

if ~exist('P','var')
  P = spm_select([1 24], 'mesh', 'Select up to 24 maps for left and right hemisphere');
end
info = cat_surf_info(P,0);

n = size(P, 1);

for i = 1:n
    
    if info(i).nvertices == 64984
        Hr.str32k = '_32k';
    else
        Hr.str32k = '';
    end
    
    % check whether name contains 'log' that indicates a logP file
    if isempty(strfind(info(i).ff, 'log'))
        Hr.logP = 0;
    end
    
    if strcmp(info(i).side, 'lh') | strcmp(info(i).side, 'rh')
        error('Display of separate hemispheres is not supported anymore');
    end

end

Hr.S{1}.name = P;
Hr.S{2}.name = P;

cat_surf_results('disp', P);

%==========================================================================
function save_image(obj, event_obj, filename)

global Hr
%%

dcm_obj = datacursormode(Hr.figure);

set(dcm_obj, 'Enable', 'off');

try
    delete(findall(gca, 'Type', 'hggroup', 'HandleVisibility', 'off'));
end

if ~exist('filename', 'var')
    
    nm = Hr.S{1}.info(1).ff;
    [pp, nm] = spm_fileparts(Hr.S{1}.name);
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
                
                if ~isempty(Hr.clip)
                    if isnan(Hr.clip(3))
                        str_thresh = '_';
                    else
                        str_thresh = sprintf('P%g_', round(1000 * 10^(-Hr.clip(3))) / 10);
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
set(Hr.figure, 'InvertHardcopy', 'off', 'PaperPositionMode', 'auto');

pos = getpixelposition(Hr.panel(1));
hh = getframe(Hr.figure,pos);

img = frame2im(hh);
if Hr.surf_sel ~= 4
    % crop image if it's not a flatmap
    sz = size(img);
    img = img(round(0.1*sz(1):sz(1)),round(0.05*sz(2):0.95*sz(2)),:);
end

col = colormap;
imwrite(img,col,fullfile(newpth,filename));

%==========================================================================
function slider_clim_min(hObject, evt)
global Hr

val = get(hObject, 'Value');
c = getappdata(Hr.patch(1), 'clim');

for ind = 1:5
    setappdata(Hr.patch(ind), 'clim', [true val c(3)]);
    col = getappdata(Hr.patch(ind), 'col');
    d = getappdata(Hr.patch(ind), 'data');
    H = updateTexture(H, ind, d, col, Hr.show_transp);
end

% update colorbar
if Hr.n_surf == 1 && ~Hr.disable_cbar
    H = show_colorbar(H);
end

Hr.clim = [true val c(3)];

%==========================================================================
function slider_clim_max(hObject, evt)
global Hr

val = get(hObject, 'Value');
c = getappdata(Hr.patch(1), 'clim');

for ind = 1:5
    setappdata(Hr.patch(ind), 'clim', [true c(2) val]);
    col = getappdata(Hr.patch(ind), 'col');
    d = getappdata(Hr.patch(ind), 'data');
    H = updateTexture(H, ind, d, col, Hr.show_transp);
end

% update colorbar
if Hr.n_surf == 1 && ~Hr.disable_cbar
    H = show_colorbar(H);
end

Hr.clim = [true c(2) val];

%==========================================================================
function checkbox_inv(obj, event_obj)
global Hr

Hr.show_inv = get(Hr.inv, 'Value');

for ind = 1:5
    setappdata(Hr.patch(ind), 'clip', Hr.clip);
    col = getappdata(Hr.patch(ind), 'col');
    setappdata(Hr.patch(ind), 'col', flipud(col));
    d = getappdata(Hr.patch(ind), 'data');
    H = updateTexture(H, ind, d, flipud(col), Hr.show_transp);
end

if ~Hr.disable_cbar
    H = show_colorbar(H);
end

%==========================================================================
function checkbox_fixscl(obj, event_obj)
global Hr

Hr.fixscl = get(Hr.fix, 'Value');

%==========================================================================
function Ho = checkbox_hide_neg(obj, event_obj)
global Hr

Hr.no_neg = get(Hr.hide_neg, 'Value');

thresh = Hr.thresh_value;
clip = getappdata(Hr.patch(1), 'clip');
clim = getappdata(Hr.patch(1), 'clim');

% get min value for both hemispheres
min_d = min(min(min(getappdata(Hr.patch(1), 'data'))), min(min(getappdata(Hr.patch(3), 'data'))));

if Hr.no_neg
    Hr.clip = [true -Inf thresh];
    Hr.clim = [true thresh clim(3)];
    set(Hr.slider_min, 'Value', 0);
else
    Hr.clip = [true -thresh thresh];
    if min_d < -thresh
        Hr.clim = [true -clim(3) clim(3)];
        set(Hr.slider_min, 'Value', -clim(3));
    end
end

for ind = 1:5
    setappdata(Hr.patch(ind), 'clip', Hr.clip);
    setappdata(Hr.patch(ind), 'clim', Hr.clim);
    col = getappdata(Hr.patch(ind), 'col');
    d = getappdata(Hr.patch(ind), 'data');
    min_d = min(min_d, min(d(:)));
    H = updateTexture(H, ind, d, col, Hr.show_transp);
end

% correct value of slider if no values are exceeding threshold
if min_d > -thresh & Hr.n_surf == 1
    set(Hr.slider_min, 'Value', 0);
end

set(Hr.atlas, 'Enable', 'on');

if ~Hr.disable_cbar
    H = show_colorbar(H);
end
if nargout, Ho = H; end

%==========================================================================
function checkbox_transp(obj, event_obj)
global Hr

Hr.show_transp = ~get(Hr.transp, 'Value');

for ind = 1:5
    col = getappdata(Hr.patch(ind), 'col');
    d = getappdata(Hr.patch(ind), 'data');
    H = updateTexture(H, ind, d, col, Hr.show_transp);
end

% update colorbar
if Hr.n_surf == 1 & ~Hr.disable_cbar
    H = show_colorbar(H);
end

%==========================================================================
function checkbox_bkg(obj, event_obj)
global Hr

Hr.white_bgk = get(Hr.bkg, 'Value');

if Hr.white_bgk
    Hr.bkg_col = [1 1 1];
else
    Hr.bkg_col = [0 0 0];
end

set(Hr.Ha, 'Color', Hr.bkg_col);
%set(get(Hr.cbar, 'Title'), 'Color', 1 - Hr.bkg_col);

title = findobj('tag','cat_surf_result_title'); 
if ~isempty(title)
  set( get( title ,'children'),'Color',  1 - Hr.bkg_col);
  %set(get(getappdata(Hr.patch(1), 'axis'), 'Title'), 'Color', 1 - Hr.bkg_col);
    %set(get(getappdata(Hr.patch(3), 'axis'), 'Title'), 'Color', 1 - Hr.bkg_col);
end

if Hr.n_surf == 1
    set(Hr.colourbar, 'XColor', 1 - Hr.bkg_col, 'YColor', 1 - Hr.bkg_col);
end

if ~Hr.disable_cbar
    H = show_colorbar(H);
end

if isfield(H, 'dataplot')
    try, set(Hr.dataplot, 'XColor', 1 - Hr.bkg_col, 'YColor', 1 - Hr.bkg_col, 'Color', Hr.bkg_col); end
end

%==========================================================================
function checkbox_info(obj, event_obj)
global Hr

Hr.show_info = get(Hr.info, 'Value');

if Hr.show_info
    delete(findobj('tag','cat_surf_result_title'));
   %   axes('Parent', Hr.panel(1), 'Position',tpos{i}, 'Visible', 'off','tag','cat_surf_results_text');
    ax = axes('Parent',Hr.panel(1),'Position',[0.5 0.82 0.9 0.05],'visible','off','tag','cat_surf_result_title');  
    text(0,1,spm_str_manip(Hr.S{1}.name, 'k150d'),'HorizontalAlignment','center','interpreter','none','Color', 1 - Hr.bkg_col,'Parent',ax);
                    
    %set(get(getappdata(Hr.patch(1), 'axis'), 'Title'), 'String', ...
    %    spm_str_manip(Hr.S{1}.name, 'k50d'), 'Interpreter', 'none', 'Color', 1 - Hr.bkg_col)
    %set(get(getappdata(Hr.patch(3), 'axis'), 'Title'), 'String', ...
    %    spm_str_manip(Hr.S{2}.name, 'k50d'), 'Interpreter', 'none', 'Color', 1 - Hr.bkg_col)
else
    delete(findobj('tag','cat_surf_result_title'));
    set(get(getappdata(Hr.patch(1), 'axis'), 'Title'), 'String', '')
    set(get(getappdata(Hr.patch(3), 'axis'), 'Title'), 'String', '')
end

%==========================================================================
function checkbox_nocbar(obj, event_obj)
global Hr

Hr.disable_cbar = get(Hr.nocbar, 'Value');

if Hr.disable_cbar
    % delete colorbar and title
    if Hr.n_surf == 1
        set(Hr.colourbar, 'Visible', 'off')
        set(get(Hr.cbar, 'Title'), 'Visible', 'off')
    else % delete only axis
        cla(Hr.cbar);
    end
else
    if Hr.n_surf == 1
        set(get(Hr.cbar, 'Title'), 'Visible', 'on')
        set(Hr.colourbar, 'Visible', 'on')
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
    Hr.axis = a;
    Hr.figure = ancestor(Hr.axis, 'figure');
    Hr.patch = findobj(Hr.axis, 'type', 'patch');
    Hr.light = findobj(Hr.axis, 'type', 'light');
    Hr.rotate3d = rotate3d(Hr.panel(1));
    setappdata(Hr.axis, 'handles', H);
else
    H = getappdata(H, 'handles');
end

%==========================================================================
function select_view(view)
global Hr

% check that view changed
if view ~= Hr.view
    
    if view > 0 % top view
        vv = [90 0; -90 0; -90 0; 90 0; 0 90];
    else % bottom view
        vv = [90 0; -90 0; -90 0; 90 0; 0 -90];
    end
    
    for ind = 1:5
        Ha = getappdata(Hr.patch(ind), 'axis');
        set(Ha, 'Position', Hr.viewpos{ind}(abs(view), :), 'View', vv(ind, :));
    end
    
    axes(Ha);
    camlight(Hr.light(1),'headlight')
    
    if isfield(H, 'dataplot')
        set(Hr.dataplot, 'Position', Hr.viewpos{6}(abs(view), :), 'Parent', Hr.panel(1), 'Color', Hr.bkg_col);
    end
    
    % save view
    Hr.view = view;
end

%==========================================================================
function select_texture(text_mode)
global Hr

% check that view changed
if text_mode ~= Hr.text_mode
    
    if text_mode == 1 % mean curvature
        Hr.text_mode = 1;
    else % sulcal depth
        Hr.text_mode = 2;
    end
    
    for ind = 1:5
        col = getappdata(Hr.patch(ind), 'col');
        d = getappdata(Hr.patch(ind), 'data');
        H = updateTexture(H, ind, d, col, Hr.show_transp);
    end
    
end

%==========================================================================
function select_border(border_mode)
global Hr

Hr.border_mode = border_mode;

if Hr.border_mode
  set(Hr.surf,   'Enable', 'off');
  if (length(Hr.S{1}.Y) == 32492 || length(Hr.S{1}.Y) == 163842 || length(Hr.S{1}.Y) == 40962) && Hr.isfsavg
    fprintf('To change underlying surface, disable Atlas Border Overlay.\n');
  end
else
  if (length(Hr.S{1}.Y) == 32492 || length(Hr.S{1}.Y) == 163842 || length(Hr.S{1}.Y) == 40962) && Hr.isfsavg
    set(Hr.surf,   'Enable', 'on');
  end
end

for ind = 1:5
    col = getappdata(Hr.patch(ind), 'col');
    d = getappdata(Hr.patch(ind), 'data');
    H = updateTexture(H, ind, d, col, Hr.show_transp);
end


%==========================================================================
function select_cursor(cursor_mode)

global Hr

dcm_obj = datacursormode(Hr.figure);
Hr.cursor_mode = cursor_mode;

switch Hr.cursor_mode
    
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
        SPM_name = fullfile(Hr.S{i}.info(1).pp, 'SPM.mat');
        
        % SPM.mat exist?
        if exist(SPM_name, 'file')
            load(SPM_name);

            % if analysis was moved we have to correct header structure
            SPM.VResMS = spm_data_hdr_read(fullfile(Hr.S{1}.info(1).pp,SPM.VResMS.fname));
            Vbeta = spm_data_hdr_read(fullfile(Hr.S{1}.info(1).pp,SPM.Vbeta(1).fname));
            for j=2:numel(SPM.Vbeta)
              Vbeta(j) = spm_data_hdr_read(fullfile(Hr.S{1}.info(1).pp,SPM.Vbeta(j).fname));
            end
            SPM.Vbeta = Vbeta;
            Hr.SPM{1} = SPM;

            str = 'predicted, adjusted or raw values?';
            Hr.predicted = spm_input(str, 1, 'b', {'predicted', 'adjusted', 'raw'}, [1 0 -1]);
            
            % ask for contrast for predicted or adjusted data
            if Hr.predicted >= 0
                Hr.Ic = spm_input('Which contrast?', 2, 'm', {SPM.xCon.name});
            end
        elseif ~isempty(Hr.S{1}.name)
            SPM_found = 0;
            spm('alert!', 'No SPM.mat file found.\nPlease check that you have not moved your files or your result file was moved from the folder where the SPM.mat is stored.', 1);
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
    cla(Hr.dataplot);
    axis(Hr.dataplot, 'off');
end

try
    dcm_obj = datacursormode(Hr.figure);
    set(dcm_obj, 'Enable', 'off');
    delete(findall(gca, 'Type', 'hggroup', 'HandleVisibility', 'off'));
end

%==========================================================================
function txt = myDataCursorCluster(obj, evt)
global Hr y

% first entries are atlases
plot_mean = Hr.cursor_mode - 4;
pos = get(evt, 'Position');

i = ismember(get(Hr.patch(1), 'vertices'), pos, 'rows');
node = find(i);
ind = 1;
node_list = 1:numel(get(Hr.patch(1), 'vertices'));

if isempty(node)
    i = ismember(get(Hr.patch(3), 'vertices'), pos, 'rows');
    node = find(i);
    ind = 3;
    node_list = 1:numel(get(Hr.patch(3), 'vertices'));
end

% get threshold from clipping
thresh = [0 0];
if ~isempty(Hr.clip)
    if ~isnan(Hr.clip(2)) & ~isnan(Hr.clip(3))
        thresh = [Hr.clip(2:3)];
    end
end

if plot_mean
    
    found_node = [];
    cluster_number = 0;
    cluster_side = 0;
    
    A = Hr.S{round(ind / 2)}.A;
    A = A + speye(size(A));
    d = getappdata(Hr.patch(ind), 'data');
    
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

% for merged meshes we only have one SPM.mat with data from both hemispheres
% add offset for right hemisphere
if round(ind / 2) == 2
    XYZ = XYZ + Hr.nY2;
end

% always one mesh
ind = 1;

[y, cbeta, CI] = get_cluster_data(H, XYZ, ind);

% if no cluster was selected set data to zero
if plot_mean & isempty(found_node)
    y(:) = 0;
    cbeta(:) = 0;
    CI(:) = 0;
end

cla(Hr.dataplot)
hold(Hr.dataplot, 'on')

set(Hr.dataplot, 'XColor', 1 - Hr.bkg_col, 'YColor', 1 - Hr.bkg_col,...
      'YGrid','on','Visible','on');

if Hr.predicted >=0
  ystr = 'contrast estimate';
  h = bar(Hr.dataplot, cbeta);
  set(h, 'FaceColor', Hr.col(1, :))
  
  % standard error
  %--------------------------------------------------------------
  CI = CI / 2;
  for j = 1:length(cbeta)
      line([j j], ([CI(j) -CI(j)] + cbeta(j)), 'LineWidth', 6, 'Color', Hr.col(2, :), 'Parent', Hr.dataplot)
  end
  set(Hr.dataplot, 'XLim', [0.4 (length(cbeta) + 0.6)], 'XTicklabel', '', 'XTick', [], 'YGrid','off')

else
  ystr = 'raw data';
  h = plot(Hr.dataplot, y);
  set(h, 'Color', Hr.col(2, :))
  set(Hr.dataplot, 'XLim', [0 (length(y))], 'XTicklabel', '', 'XTick', [])
end

xX = Hr.SPM{1}.xX;
iH = xX.iH;

% plot group coding for Anovas with more than 1 group and native data
if ~isempty(iH) & numel(iH)>1 & (Hr.predicted < 0) & isempty(xX.iB)
    yl = get(Hr.dataplot,'YLim');
    pcol = gray(numel(iH)+2);
    for i=1:numel(iH)
        ind_data = find(any(xX.X(:,xX.iH(i)),2));
        
        % plot only if ind_data is a continuous row
        if all(diff(ind_data)==1)
            line([min(ind_data)-0.5 max(ind_data)+0.5], [yl(1) yl(1)], 'LineWidth', 6, 'Color', pcol(i+1,:), 'Parent', Hr.dataplot)
        end
    end
    
end

nm = Hr.S{1}.info(1).ff;

Ic = [];
% end with _0???.ext?
if length(nm) > 4
    if strcmp(nm(length(nm) - 4:length(nm) - 3), '_0')
        Ic = str2double(nm(length(nm) - 3:length(nm)));
    end
else
  if isfield(H,'Ic')
      Ic = Hr.Ic;
  end
end

if ~isempty(Ic)
    xlabel(Hr.dataplot, Hr.SPM{round(ind / 2)}.xCon(Ic).name, 'FontSize', Hr.FS+3, 'Color', 1 - Hr.bkg_col)
end

if plot_mean
    ylabel(Hr.dataplot, sprintf('mean %s\ninside cluster',ystr), 'FontSize', Hr.FS+3, 'Color', 1 - Hr.bkg_col)
else
    ylabel(Hr.dataplot, sprintf('%s',ystr), 'FontSize', Hr.FS+3, 'Color', 1 - Hr.bkg_col)
end

hold(Hr.dataplot, 'off')

assignin('base', 'y', y);

%==========================================================================
function [y, cbeta, CI] = get_cluster_data(H, XYZ, ind)

SPM = Hr.SPM{round(ind / 2)};

Ic = [];
if isfield(H,'Ic')
    Ic = Hr.Ic;
end

predicted = Hr.predicted;

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

if Hr.cursor_mode == 1
    txt = {'Desikan DK40'};
elseif Hr.cursor_mode == 2
    txt = {'Destrieux 2009'};
elseif Hr.cursor_mode == 3
    txt = {'HCP_MMP1'};
end

i = ismember(get(Hr.patch(1), 'vertices'), pos, 'rows');
node = find(i);
ind = 1;

if isempty(node)
    i = ismember(get(Hr.patch(3), 'vertices'), pos, 'rows');
    node = find(i);
    ind = 2;
end

rdata_pos = Hr.rdata{Hr.cursor_mode}(node, ind);

rcsv = Hr.rcsv{Hr.cursor_mode};

for j = 2:size(rcsv, 1)
    if rdata_pos == rcsv{j, 1}
        txt = {txt{:} [Hr.S{ind}.side ' ' rcsv{j, 2}]};
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

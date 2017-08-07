%==========================================================================
% PALE - Point And Label setting Extension
%==========================================================================
% Plugin for spm_image to add points to labels to loadet image. Just start
% it with the context-menu entry: Start PALE.
% ADD:
% To add a Point click on screen while holding one of the keys 1-9, which 
% represents the boundarys.
% Alternatively you can choose the boundary with the GUI. But then you have
% to deselect the boundary manually with the deselect button. You can change
% the current slice with the keys 'e' and 'd'. 
% DELETE: 
% You can delete in two ways.
% When you want to remove the last points of the current boundary, you can
% use 'r'. Also you can click on the point you want to remove. For better
% Visibility of them, pressing 'f' makes all visible points bigger. Moving
% the crosshair undo this effect.
%==========================================================================

function ret = spm_ov_pale(varargin)
    global st;
    switch varargin{1}
        case 'context_menu'
            uimenu(varargin{3}, 'Label', 'Start PALE', 'Callback',@start);
        case 'redraw'
            % precalculate matrix 
            is = inv(st.Space);
            a = is(1:3,1:3);
            b = is(1:3,4);
            % color declaration - change on wish
            colors = [ 1 0 0; 0 0 1; 0 1 0; 0 1 1; 1 1 0; 1 0 1; 0 1 0.5; 1 0.5 1; 0.5 0 1];
            fields = fieldnames(st.vols{1}.pale.labelData);
            if(st.vols{1}.pale.currentBoundary == 0)
                % remove all old points
                if(~isempty(st.vols{1}.pale.lines))
                    % delete object
                    delete(st.vols{1}.pale.lines(:));
                    % delete addresses
                    st.vols{1}.pale.lines = [];
                    st.vols{1}.pale.lineC = {};
                end
                % screen all points in bound
                for j=1:9
                    color = colors(j,:);
                    for i=1:size(st.vols{1}.pale.labelData.(fields{st.vols{1}.pale.currentLabel}).(['Boundary' num2str(j)]),2)
                        for k=1:3
                            if ~inSlice(st.centre(k), st.vols{1}.pale.labelData.(fields{st.vols{1}.pale.currentLabel}).(['Boundary' num2str(j)])(:,i), k)
                                continue;
                            end
                            display_point(st.vols{1}.pale.labelData.(fields{st.vols{1}.pale.currentLabel}).(['Boundary' num2str(j)])(:,i), st.bb, k, color, a, b, j, i);
                        end
                    end
                end
            else
                XYZmm = spm_orthviews('Pos');
                 % when any label is selected, write in it
                if isempty(st.vols{1}.pale.labelData.(fields{st.vols{1}.pale.currentLabel}).(['Boundary' num2str(st.vols{1}.pale.currentBoundary)]))
                    st.vols{1}.pale.labelData.(fields{st.vols{1}.pale.currentLabel}).(['Boundary' num2str(st.vols{1}.pale.currentBoundary)])(:,end+1) = XYZmm(:);
                elseif  ~isequal(round(st.vols{1}.pale.labelData.(fields{st.vols{1}.pale.currentLabel}).(['Boundary' num2str(st.vols{1}.pale.currentBoundary)])(:,end),0), round(XYZmm(:),0))
                    st.vols{1}.pale.labelData.(fields{st.vols{1}.pale.currentLabel}).(['Boundary' num2str(st.vols{1}.pale.currentBoundary)])(:,end+1) = XYZmm(:);
                    for k=1:3
                        display_point(XYZmm(:,end), st.bb,k, colors(st.vols{1}.pale.currentBoundary,:),a,b, st.vols{1}.pale.currentBoundary, size(st.vols{1}.pale.labelData.(fields{st.vols{1}.pale.currentLabel}).(['Boundary' num2str(st.vols{1}.pale.currentBoundary)]),2));
                    end
                end
            end
        otherwise
    end
end


function res = inSlice(slice, pos, axe, bb)
    if round(slice) == round(pos(axe))
        res = 1;
        return;
    end
    res = 0;
    return;
end

function display_point(point, bb, axis, color, a,b, boundary, id)
    global st;
    % transform point into screen coordinates with precalculated matrices
    pos = a * point(:) + b;
    % depending on axis draw point
    switch axis
    case 1
        % this view can be flipped, so st.mode contains this data
        if st.mode == 0
            st.vols{1}.pale.lines(end+1) = line(st.vols{1}.ax{3}.ax, pos(3)-bb(1,3)+1, pos(2)-bb(1,2)+1, 'Marker', 's', 'MarkerSize', st.vols{1}.pale.markerSize, 'Color', color, 'MarkerFaceColor', color, 'ButtonDownFcn',@lineCallback);
        else
            st.vols{1}.pale.lines(end+1) = line(st.vols{1}.ax{3}.ax, bb(2,2)+1-pos(2), pos(3)-bb(1,3)+1, 'Marker', 's', 'MarkerSize', st.vols{1}.pale.markerSize, 'Color', color, 'MarkerFaceColor', color, 'ButtonDownFcn',@lineCallback);
        end
    case 2
        st.vols{1}.pale.lines(end+1) = line(st.vols{1}.ax{2}.ax, pos(1)-bb(1,1)+1, pos(3)-bb(1,3)+1, 'Marker', 's', 'MarkerSize', st.vols{1}.pale.markerSize, 'Color', color, 'MarkerFaceColor', color, 'ButtonDownFcn',@lineCallback);
    case 3
        st.vols{1}.pale.lines(end+1) = line(st.vols{1}.ax{1}.ax, pos(1)-bb(1,1)+1, pos(2)-bb(1,2)+1, 'Marker', 's', 'MarkerSize', st.vols{1}.pale.markerSize, 'Color', color, 'MarkerFaceColor', color, 'ButtonDownFcn',@lineCallback);
    end
    % save the point object in pale structure to enable removing them later
    s = struct();
    s.boundary = boundary;
    s.id = id;
    st.vols{1}.pale.lineC{end+1} = s;
end

function start(varargin)
    global st;
    % activate pale
    if isfield(st.vols{1}, 'pale') == 0
        st.paleUI_handle = pale_ui;
        set(st.fig, 'CloseRequestFcn', @close);
        set(st.fig, 'KeyPressFcn', @key_pressed);
        set(st.fig, 'KeyReleaseFcn', @key_released);
    else 
        disp('PALE already running!');
    end
end

function lineCallback(varargin)
global st;
line = get(gcf, 'CurrentObject');
try
    for a=1:size(st.vols{1}.pale.lines,2)
        anyLine = get(st.vols{1}.pale.lines(a));
        if(anyLine.XData == line.XData && anyLine.YData == line.YData)
            st.vols{1}.pale.labelData.(['Region' num2str(st.vols{1}.pale.currentLabel)]).(['Boundary' num2str(st.vols{1}.pale.lineC{a}.boundary)])(:,st.vols{1}.pale.lineC{a}.id) = [];
            delete(st.vols{1}.pale.lines(a))
            st.vols{1}.pale.lines(a) = [];
            st.vols{1}.pale.lineC{a} = [];
            st.vols{1}.pale.lineC = st.vols{1}.pale.lineC(~cellfun('isempty',st.vols{1}.pale.lineC));
            return;
        end
end
catch ME
    redraw();
end
end

function key_pressed(hObject, eventdata, handles)
    global st;
    if isfield(st.vols{1}, 'pale') ~= 0
        % when key is in possible boundaryID limits
        if(str2double(eventdata.Character) <= 9 && str2double(eventdata.Character) > 0)
            st.vols{1}.pale.currentBoundary = str2double(eventdata.Character);
            st.vols{1}.pale.lastBoundary = str2double(eventdata.Character);
        end
        switch eventdata.Character
            case 'e'
                pos = spm_orthviews('pos');
                % get axis
                currentAxis = get(gcf,'CurrentObject');
                for i=1:3
                    if currentAxis.Position == st.vols{1}.ax{i}.ax.Position
                        if(i == 1)
                            axis = 3;
                        elseif(i == 2)
                            axis = 2;
                        else
                            axis = 1;
                        end
                        break;
                    end
                end
                pos(axis) = pos(axis)+1;
                spm_orthviews('reposition',pos);
            case 'd'
                pos = spm_orthviews('pos');
                % get axis
                currentAxis = get(gcf,'CurrentObject');
                for i=1:3
                    if currentAxis.Position == st.vols{1}.ax{i}.ax.Position
                        if(i == 1)
                            axis = 3;
                        elseif(i == 2)
                            axis = 2;
                        else
                            axis = 1;
                        end
                        break;
                    end
                end
                pos(axis) = pos(axis)-1;
                spm_orthviews('reposition',pos);
            case 'r'
                % remove last
                % prevent errors during empty label data
                if isempty(st.vols{1}.pale.lineC)
                    return;
                end
                % continue when boundary is the same as last selected
                if st.vols{1}.pale.lineC{end}.boundary ~= st.vols{1}.pale.lastBoundary
                    return;
                end
                refObj = st.vols{1}.pale.lineC{end};
                refEndObj = st.vols{1}.pale.lineC{end};
                st.vols{1}.pale.labelData.(['Region' num2str(st.vols{1}.pale.currentLabel)]).(['Boundary' num2str(st.vols{1}.pale.lastBoundary)])(:,end) = [];
                % remove all point duplicities from other orthviews
                while refObj.id == refEndObj.id && refObj.boundary == refEndObj.boundary
                    delete(st.vols{1}.pale.lines(end))
                    st.vols{1}.pale.lines(end) = [];
                    st.vols{1}.pale.lineC{end} = [];
                    st.vols{1}.pale.lineC = st.vols{1}.pale.lineC(~cellfun('isempty',st.vols{1}.pale.lineC));
                    if(isempty(st.vols{1}.pale.lineC))
                       break;
                    end
                    refEndObj = st.vols{1}.pale.lineC{end};
                end
            case 'f'
                % delete mode -> draw marker bigger
                % draw 
                if st.vols{1}.pale.markerSize == 3
                    st.vols{1}.pale.markerSize = 7;
                    redraw();
                end
            otherwise
        end
    end
end

function key_released(hObject, eventdata, handles)
    global st;
    if isfield(st.vols{1}, 'pale') ~= 0
        if(str2double(eventdata.Key) > 0 && str2double(eventdata.Key) < 10)
            st.vols{1}.pale.currentBoundary = 0;
        end
        st.vols{1}.pale.markerSize = 3;
    end
end

function close(hObject, eventdata)
global st;
if(isfield(st, 'paleUI_handle'))
    delete(st.paleUI_handle);
    st = rmfield(st, 'paleUI_handle');
end
delete(hObject);
end

function redraw()
    global st;
    % precalculate matrix 
    is = inv(st.Space);
    a = is(1:3,1:3);
    b = is(1:3,4);
    % color declaration - change on wish
    colors = [ 1 0 0; 0 0 1; 0 1 0; 0 1 1; 1 1 0; 1 0 1; 0 1 0.5; 1 0.5 1; 0.5 0 1];
    fields = fieldnames(st.vols{1}.pale.labelData);
    if(st.vols{1}.pale.currentBoundary == 0)
        % remove all old points
        if(~isempty(st.vols{1}.pale.lines))
            % delete object
            delete(st.vols{1}.pale.lines(:));
            % delete addresses
            st.vols{1}.pale.lines = [];
            st.vols{1}.pale.lineC = {};
        end
        % screen all points in bound
        for j=1:9
            color = colors(j,:);
            for i=1:size(st.vols{1}.pale.labelData.(fields{st.vols{1}.pale.currentLabel}).(['Boundary' num2str(j)]),2)
                for k=1:3
                    if ~inSlice(st.centre(k), st.vols{1}.pale.labelData.(fields{st.vols{1}.pale.currentLabel}).(['Boundary' num2str(j)])(:,i), k)
                        continue;
                    end
                    display_point(st.vols{1}.pale.labelData.(fields{st.vols{1}.pale.currentLabel}).(['Boundary' num2str(j)])(:,i), st.bb, k, color, a, b, j, i);
                end
            end
        end
    end
end
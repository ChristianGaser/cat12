function varargout = pale_ui(varargin)
% PALE_UI MATLAB code for pale_ui.fig
%      PALE_UI, by itself, creates a new PALE_UI or raises the existing
%      singleton*.
%
%      H = PALE_UI returns the handle to a new PALE_UI or the handle to
%      the existing singleton*.
%
%      PALE_UI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PALE_UI.M with the given input arguments.
%
%      PALE_UI('Property','Value',...) creates a new PALE_UI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pale_ui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pale_ui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pale_ui

% Last Modified by GUIDE v2.5 10-May-2017 00:40:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pale_ui_OpeningFcn, ...
                   'gui_OutputFcn',  @pale_ui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before pale_ui is made visible.
function pale_ui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pale_ui (see VARARGIN)

% Choose default command line output for pale_ui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pale_ui wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global st;
splittedPath = strsplit(st.vols{1}.fname,'/');
% initalize pale main struct
st.vols{1}.pale = struct();
st.vols{1}.pale.image               = {cell2mat(splittedPath(end))};
st.vols{1}.pale.path                = {st.vols{1}.fname};
st.vols{1}.pale.labelData           = struct();
st.vols{1}.pale.currentLabel        = 1;
st.vols{1}.pale.lines               = [];
st.vols{1}.pale.removeMode          = 0;
st.vols{1}.pale.currentBoundary     = 0;
st.vols{1}.pale.lastBoundary        = 0;
st.vols{1}.pale.lineC               = {};
st.vols{1}.pale.markerSize          = 3;

% initalize first label
newBoundary = struct();
for i=1:9
    newBoundary.(['Boundary' num2str(i)]) = [];
end
st.vols{1}.pale.labelData.('Region1') = newBoundary;
handles.popupmenu3.String = {'Region 1'};


% --- Outputs from this function are returned to the command line.
function varargout = pale_ui_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1. : new
function pushbutton1_Callback(hObject, eventdata, handles)
global st;
if isfield(st.vols{1}, 'pale') ~= 0
    if ~isempty('st.vols{1}.pale.lines')
        delete(st.vols{1}.pale.lines(:));
        st.vols{1}.pale.lines = [];
        st.vols{1}.pale.lineC = {};
    end
end
splittedPath = strsplit(st.vols{1}.fname,'/');

% initalize pale main struct
st.vols{1}.pale = struct();
st.vols{1}.pale.image               = {cell2mat(splittedPath(end))};
st.vols{1}.pale.path                = {st.vols{1}.fname};
st.vols{1}.pale.labelData           = struct();
st.vols{1}.pale.currentLabel        = 1;
st.vols{1}.pale.lines               = [];
st.vols{1}.pale.removeMode          = 0;
st.vols{1}.pale.currentBoundary     = 0;
st.vols{1}.pale.lastBoundary        = 0;
st.vols{1}.pale.lineC               = {};
st.vols{1}.pale.markerSize          = 3;

% initalize first label
newBoundary = struct();
for i=1:9
    newBoundary.(['Boundary' num2str(i)]) = [];
end
st.vols{1}.pale.labelData.('Region1') = newBoundary;
handles.popupmenu3.String = {'Region 1'};

% --- Executes on button press in pushbutton2. : load file
function pushbutton2_Callback(hObject, eventdata, handles)
global st;
try
    [file, path] = uigetfile();
    load([path, file]);
    st.vols{1}.pale = pale;
    % update dropdowns : rid
    RIDs = fieldnames(st.vols{1}.pale.labelData);
    for i=1:size(RIDs,1)
        tmp=get(handles.popupmenu3,'string');
        tmp{end+1}=RIDs{i};
        set(handles.popupmenu3,'string',tmp);
    end
catch ME 
    disp('Error: Invalid Selection');
end

% --- Executes on button press in pushbutton3. : close
function pushbutton3_Callback(hObject, eventdata, handles)
global st;
if isfield(st.vols{1}, 'pale') ~= 0
    handles.popupmenu3.String = {' - '};
    handles.popupmenu3.Value = 1;
    if ~isempty(st.vols{1}.pale.lines)
        delete(st.vols{1}.pale.lines(:));
        st.vols{1}.pale.lines = [];
        st.vols{1}.pale.lineC = {};
    end
    st.vols{1} = rmfield(st.vols{1}, 'pale');
else
    disp('you have to press: "new" or "load" a file');
end

% --- Executes on button press in pushbutton4.: save
function pushbutton4_Callback(hObject, eventdata, handles)
global st; 
if isfield(st.vols{1}, 'pale') ~= 0
    name = strsplit(st.vols{1}.pale.image{1}, '.');
    [file,path] = uiputfile([name{1} '_PALE' '.mat']);
    pale = st.vols{1}.pale;
    % delete lines
    pale.lines = [];
    pale.lineC = {};
    % save in file
    try
        save([path file], 'pale');
    catch ME 
        disp('Error: Invalid Selection');
    end
else
    disp('you have to press: "new" or "load" a file');
end

% --- Executes on selection change in popupmenu2. : Boundary
function popupmenu2_Callback(hObject, eventdata, handles)
global st;
if isfield(st.vols{1}, 'pale') ~= 0
    st.vols{1}.pale.currentBoundary = hObject.Value;
else
    disp('you have to press: "new" or "load" a file');
end

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
hObject.String = {'Boundary 1'};
for i=1:8
    tmp=get(hObject,'string');
    elements = size(tmp,1);
    tmp{end+1}=['Boundary ' num2str(elements+1)];
    set(hObject,'string',tmp);
end

% --- Executes on selection change in popupmenu3. : Region
function popupmenu3_Callback(hObject, eventdata, handles)
global st;
if isfield(st.vols{1}, 'pale') ~= 0
    st.vols{1}.pale.currentLabel = hObject.Value;
else
    disp('you have to press: "new" or "load" a file');
end

% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
hObject.String = {'Region 1'};


% --- Executes on button press in togglebutton1. : remove Mode
function togglebutton1_Callback(hObject, eventdata, handles)
global st;
if isfield(st.vols{1}, 'pale') ~= 0
    st.vols{1}.pale.removeMode = hObject.Value;
else
    disp('you have to press: "new" or "load" a file');
end

% --- Executes on button press in pushbutton7. : add Label
function pushbutton7_Callback(hObject, eventdata, handles)
global st;
if isfield(st.vols{1}, 'pale') ~= 0
    % create empty Boundary
    newBoundary = struct();
    for j=1:9
        newBoundary.(['Boundary' num2str(j)]) = [];
    end
    % add new Label
    allLabels = fieldnames(st.vols{1}.pale.labelData);
    st.vols{1}.pale.labelData.(['Region' num2str(size(allLabels,1)+1)]) = newBoundary;
    % update popupmenue
    tmp=get(handles.popupmenu3,'string');
    elements = size(tmp,1);
    tmp{end+1}=['Region ' num2str(elements+1)];
    set(handles.popupmenu3,'string',tmp);
else
    disp('you have to press: "new" or "load" a file');
end

% --- Executes on button press in pushbutton8. : remove specific Label
function pushbutton8_Callback(hObject, eventdata, handles)
global st;
if isfield(st.vols{1}, 'pale') ~= 0
    try
        input = (inputdlg('Enter Labels ID:'));
        labels = fieldnames(st.vols{1}.pale.labelData);
        if(input{1} <= size(labels,1) && input{1} > 0)
            % delete points on screen
            delete(st.vols{1}.pale.lines(:));
            st.vols{1}.pale.lines = [];
            st.vols{1}.pale.lineC = {};
            % override label with empty boundaries
            newBoundary = struct();
            for i=1:9
                newBoundary.(['Boundary' num2str(i)]) = [];
            end
            st.vols{1}.pale.labelData.(labels{i}) = newBoundary;
        end
    catch ME
        disp('Error: Invalid Selection');
    end
else
    disp('you have to press: "new" or "load" a file');
end
    
function pushbutton9_Callback(hObject, eventdata, handles)
global st;
if isfield(st.vols{1}, 'pale') ~= 0
    st.vols{1}.pale.currentBoundary = 0;
end

% --- Executes on button press in pushbutton10. : remove last
function pushbutton10_Callback(hObject, eventdata, handles)
global st;
if isfield(st.vols{1}, 'pale') ~= 0
    % prevent errors during empty label data
    if isempty(st.vols{1}.pale.lineC)
        return;
    end
    % continue when boundary is the same as last selected
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
else
    disp('you have to press: "new" or "load" a file');
end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
global st;
if isfield(st.vols{1}, 'pale') ~= 0
    if ~isempty('st.vols{1}.pale.lines')
        delete(st.vols{1}.pale.lines(:));
        st.vols{1}.pale.lines = [];
        st.vols{1}.pale.lineC = {};
    end
    st.vols{1} = rmfield(st.vols{1}, 'pale');
end
delete(hObject);


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pale_ui_help;

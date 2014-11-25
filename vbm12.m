function varargout = vbm12(varargin)
% VBM12 M-file for vbm12.fig
%      VBM12, by itself, creates a new VBM12 or raises the existing
%      singleton*.
%
%      H = VBM12 returns the handle to a new VBM12 or the handle to
%      the existing singleton*.
%
%      VBM12('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VBM12.M with the given input arguments.
%
%      VBM12('Property','Value',...) creates a new VBM12 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DEM_demo_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vbm12_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above vbm12title to modify the response to help vbm12

% Last Modified by GUIDE v2.5 25-Nov-2014 10:10:09


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vbm12_OpeningFcn, ...
                   'gui_OutputFcn',  @vbm12_OutputFcn, ...
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


% This creates the 'background' image
if ~nargin
  ha = axes('units','normalized','position',[0 0.87 1 0.13]);
  uistack(ha,'bottom');
  I=imread(fullfile(spm('dir'),'toolbox','vbm12','html','contact.jpg'));
  hi = imagesc(I);
  text(80,100,'VBM12 Toolbox','Color',[1 1 1],'Fontsize',22,'Fontweight','bold');
  set(ha,'handlevisibility','off','visible','off');
end


% --- Executes just before vbm12 is made visible.
function vbm12_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vbm12 (see VARARGIN)

% Choose default command line output for vbm12
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = vbm12_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
spm_jobman('interactive','','spm.tools.vbm.estwrite');

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
spm_jobman('interactive','','spm.tools.vbm.tools.showslice');

% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
spm('PopUpCB',gcbo);

% --- Executes on button press in pushbutton84.
function pushbutton84_Callback(hObject, eventdata, handles)
P = spm_select([1 Inf],'^SPM\.mat$','Select SPM.mat file(s)');
for i=1:size(P,1)
    swd      = spm_file(P(i,:),'fpath');
    load(fullfile(swd,'SPM.mat'));
    SPM.swd  = swd;
    vbm_spm(SPM);
end

% --- Executes on button press in pushbutton155.
function pushbutton155_Callback(hObject, eventdata, handles)
spm_jobman('interactive','','spm.tools.vbm.tools.surfextract');

% --- Executes on button press in pushbutton154.
function pushbutton154_Callback(hObject, eventdata, handles)
spm_jobman('interactive','','spm.tools.vbm.tools.surfresamp');

% --- Executes on button press in pushbutton156.
function pushbutton156_Callback(hObject, eventdata, handles)
P=spm_select([1 24],'gifti','Select surface');
for i=1:size(P,1)
    h = spm_mesh_render(deblank(P(i,:)));
    set(h.figure,'MenuBar','none','Toolbar','none','Name',spm_file(P(i,:),'short40'),'NumberTitle','off');
    spm_mesh_render('ColourMap',h.axis,jet); spm_mesh_render('ColourBar',h.axis,'on');
end

% --- Executes on button press in pushbutton157.
function pushbutton157_Callback(hObject, eventdata, handles)
spm_jobman('interactive','','spm.tools.vbm.tools.T2x');

% --- Executes on button press in pushbutton158.
function pushbutton158_Callback(hObject, eventdata, handles)
cg_slice_overlay;

% --- Executes on button press in pushbutton159.
function pushbutton159_Callback(hObject, eventdata, handles)
spm_jobman('interactive','','spm.tools.vbm.tools.F2x');

% --- Executes on button press in pushbutton160.
function pushbutton160_Callback(hObject, eventdata, handles)
spm('alert',evalc('cg_vbm_update(1)'),'VBM Update');


% --- Executes on button press in pushbutton161.
function pushbutton161_Callback(hObject, eventdata, handles)
try,open(fullfile(spm('dir'),'toolbox','vbm12','VBM12-Manual.pdf'));end


% --- Executes on button press in pushbutton162.
function pushbutton162_Callback(hObject, eventdata, handles)
try,web('http://dbm.neuro.uni-jena.de/vbm','-browser');end

% --- Executes on button press in pushbutton164.
function pushbutton164_Callback(hObject, eventdata, handles)
close(gcf);

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
spm_jobman('interactive','','spm.tools.vbm.tools.long');



% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton165.
function pushbutton165_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton165 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on pushbutton1 and none of its controls.
function pushbutton1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function pushbutton161_CreateFcn(hObject, eventdata, handles)


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
if exist(fullfile(spm('dir'),'toolbox','TFCE'))
    set(hObject,'enable','on');
else
    set(hObject,'enable','off');
end

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

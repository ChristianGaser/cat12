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

% Last Modified by GUIDE v2.5 03-Dec-2014 09:15:35


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

% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


%-------------------------------------------------------------------

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
spm_jobman('interactive','','spm.tools.vbm.estwrite');

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
spm_jobman('interactive','','spm.tools.vbm.tools.long');

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
spm_jobman('interactive','','spm.tools.vbm.tools.showslice');

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
spm_jobman('interactive','','spm.stats.factorial_design');

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
P = spm_select([1 Inf],'^SPM\.mat$','Select SPM.mat file(s)');
for i=1:size(P,1)
    swd      = spm_file(P(i,:),'fpath');
    load(fullfile(swd,'SPM.mat'));
    SPM.swd  = swd;
    vbm_spm(SPM);
end

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
spm_jobman('interactive','','spm.tools.vbm.tools.surfextract');

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
spm_jobman('interactive','','spm.tools.vbm.tools.surfresamp');

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
P=spm_select([1 24],'gifti','Select surface');
for i=1:size(P,1)
  if 0 % use spm
    h = spm_mesh_render(deblank(P(i,:)));
    set(h.figure,'MenuBar','none','Toolbar','none','Name',spm_file(P(i,:),'short40'),'NumberTitle','off');
    spm_mesh_render('ColourMap',h.axis,jet); spm_mesh_render('ColourBar',h.axis,'on');
  else 
    %% use vbm
    h = vbm_mesh_render(deblank(P(i,:)));
    set(h.figure,'MenuBar','none','Toolbar','none','Name',spm_file(P(i,:),'short40'),'NumberTitle','off');
    vbm_mesh_render('ColourMap',h.axis,jet); vbm_mesh_render('ColourBar',h.axis,'on');
  end
end

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
cg_slice_overlay;

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
close(gcf);

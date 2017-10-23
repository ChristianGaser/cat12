function varargout = cat12(varargin)
% ______________________________________________________________________
% CAT12 Toolbox wrapper to start CAT with different user modes or 
% default files.  Changing the user mode requires restarting of CAT and
% SPM.  The expert user mode allows to control further parameters and  
% semi-evaluated functions, whereas the developer mode contain parameter
% for internal tests and unsafe functions.
% 
%   cat12(action)
%   
%   CAT user modes:
%     action = ['default','expert','developer'] 
%
%   CAT default files for other species (in development):
%     action = ['oldwoldmonkeys'|'greaterapes']
%
%   CAT start with own default files:
%     action = 'select' 
%     action = 'mypath/cat_defaults_mydefaults'
%
% ______________________________________________________________________
% Christian Gaser, Robert Dahnke
% $Id$

% CAT12 M-file for cat12.fig
%      CAT12, by itself, creates a new CAT12 or raises the existing
%      singleton*.
%
%      H = CAT12 returns the handle to a new CAT12 or the handle to
%      the existing singleton*.
%
%      CAT12('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CAT12.M with the given input arguments.
%
%      CAT12('Property','Value',...) creates a new CAT12 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DEM_demo_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cat12_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above title to modify the response to help cat12

% Last Modified by GUIDE v2.5 08-Sep-2017 11:56:03

if nargin==0 
  spm_cat12;
  return;
elseif nargin==1 && ~strcmp(varargin{1},'fig')
  spm_cat12(varargin{1});
  return;
elseif nargin==2 && ~strcmp(varargin{1},'fig')
  spm_cat12(varargin{1},varargin{2});
  return;
end

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1; 
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cat12_OpeningFcn, ...
                   'gui_OutputFcn',  @cat12_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ~strcmp(varargin{1},'fig') && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



% --- Executes just before cat12 is made visible.
function cat12_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cat12 (see VARARGIN)

% Choose default command line output for cat12
handles.output = hObject; 

% Update handles structure
guidata(hObject, handles);

% enable/disable different menus if TFCE is installed or not
if exist(fullfile(spm('dir'),'toolbox','TFCE'))
    set(handles.popupmenu2,'String',{'Treshold-Free Cluster Enhancement...','Call TFCE Toolbox'});
else
    set(handles.popupmenu2,'String',{'Treshold-Free Cluster Enhancement...','Install TFCE Toolbox'});
end

% --- Outputs from this function are returned to the command line.
function varargout = cat12_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

spm_clf('Interactive'); 

% Get default command line output from handles structure
varargout{1} = handles.output;

expert  = cat_get_defaults('extopts.expertgui'); 
species = cat_get_defaults('extopts.species'); 
switch expert
  case 1, set(handles.CAT,'color', [0.85 0.85 0.85]);
  case 2, set(handles.CAT,'color', [0.93 0.93 0.93]); 
end

FS = spm('FontSizes');

% This creates the 'background' image
handles.ha = axes('units','normalized','position',[0 0.87 1 0.13]);
I = imread(fullfile(spm('dir'),'toolbox','cat12','html','images','contact.jpg'));
imagesc(I);
axis off; 
text(80,140,'Computational Anatomy Toolbox','Color',[1 1 1],'Fontsize',FS(14),'Fontweight','bold');
switch species
  case 'human',           speciesdisp = ''; 
  case 'ape_greater',     speciesdisp = ' (greater apes)';
  case 'ape_lesser',      speciesdisp = ' (lesser apes)';
  case 'monkey_oldworld', speciesdisp = ' (oldworld monkeys)'; 
  case 'monkey_newworld', speciesdisp = ' (newworld monkeys)'; 
  case 'dog',             speciesdisp = ' (dogs)'; 
  otherwise               speciesdisp = ''; 
end
switch expert
  case 1, text(80,90,['Expert Mode'    speciesdisp],'Color',[0.1 0.7 1.0],'Fontsize',FS(10),'Fontweight','bold'); 
  case 2, text(80,90,['Developer Mode' speciesdisp],'Color',[1.0 0.0 0.0],'Fontsize',FS(10),'Fontweight','bold');
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function CAT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CAT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


%-------------------------------------------------------------------

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
spm_jobman('interactive','','spm.tools.cat.estwrite');

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
spm_jobman('interactive','','spm.tools.cat.tools.long');

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
spm_jobman('interactive','','spm.tools.cat.tools.showslice');

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
spm_jobman('interactive','','spm.stats.factorial_design');

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
P = spm_select([1 Inf],'^SPM\.mat$','Select SPM.mat file(s)');

% workaround to use fsaverage surface as SurfaceID (for displaying results)
for i=1:size(P,1)
    swd      = spm_file(P(i,:),'fpath');
    load(fullfile(swd,'SPM.mat'));
    SPM.swd  = swd;

    fsavgDir = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces');
    
    % check that folder exist and number of vertices fits
    if exist(fsavgDir,'dir') == 7 && (SPM.xY.VY(1).dim(1) == 163842 || SPM.xY.VY(1).dim(1) == 327684 || SPM.xY.VY(1).dim(1) == 655368)
      [pp,ff]   = spm_fileparts(SPM.xY.VY(1).fname);

      % find mesh string      
      hemi_ind = [];
      hemi_ind = [hemi_ind strfind(ff,'mesh.')];
      if ~isempty(hemi_ind)
        
        SPM.xY.VY(1).private.private.metadata = struct('name','SurfaceID','value',fullfile(fsavgDir, 'mesh.central.freesurfer.gii'));
        M0 = gifti({fullfile(fsavgDir, 'lh.central.freesurfer.gii'), fullfile(fsavgDir, 'rh.central.freesurfer.gii')});
        G.faces = [M0(1).faces; M0(2).faces+size(M0(1).vertices,1)];
        G.vertices = [M0(1).vertices; M0(2).vertices];

        % cerebellar lobes?
        if SPM.xY.VY(1).dim(1) == 655368
          M0 = gifti({fullfile(fsavgDir, 'lc.central.freesurfer.gii'), fullfile(fsavgDir, 'rc.central.freesurfer.gii')});
          G.faces = [G.faces; M0(1).faces+2*size(M0(1).vertices,1); M0(2).faces+3*size(M0(1).vertices,1)];
          G.vertices = [G.vertices; M0(1).vertices; M0(2).vertices];
        end
        
        SPM.xVol.G = G;
        
        % remove memory demanding faces and vertices which are not necessary
        for i=1:length(SPM.xY.VY)
          SPM.xY.VY(i).private.faces = [];
          SPM.xY.VY(i).private.vertices = [];
        end
        
        save(fullfile(swd,'SPM.mat'),'SPM', '-v7.3');
      else

        % find lh|rh string
        hemi_ind = [];
        hemi_ind = [hemi_ind strfind(ff,'lh.')];
        hemi_ind = [hemi_ind strfind(ff,'rh.')];
        hemi = ff(hemi_ind:hemi_ind+1);
        if ~isempty(hemi)
          SPM.xY.VY(1).private.private.metadata = struct('name','SurfaceID','value',fullfile(fsavgDir,[hemi '.central.freesurfer.gii']));
          G = fullfile(fsavgDir,[hemi '.central.freesurfer.gii']);
          SPM.xVol.G = gifti(G);
          
          % remove memory demanding faces and vertices which are not necessary
          for i=1:length(SPM.xY.VY)
            SPM.xY.VY(i).private.faces = [];
            SPM.xY.VY(i).private.vertices = [];
          end
          
          save(fullfile(swd,'SPM.mat'),'SPM', '-v7.3');
        end
      end
    end

    spm_spm(SPM);
    
end

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
spm_jobman('interactive','','spm.tools.cat.stools.surfcalc');

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
spm_jobman('interactive','','spm.tools.cat.stools.surfresamp');

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
cat_surf_display;

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
F = spm_figure('FindWin','Menu');
% close SPM windows, if no Menu window exist
if isempty(F)
  spm('Quit')
end
close(gcf);

% --- Executes on selection change in popupmenu11.
function popupmenu11_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu11 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu11


% --- Executes during object creation, after setting all properties.
function popupmenu11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu12.
function popupmenu12_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu12 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu12


% --- Executes during object creation, after setting all properties.
function popupmenu12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton169.
function pushbutton169_Callback(hObject, eventdata, handles)
spm_jobman('interactive','','spm.tools.cat.tools.calcvol');


% --- Executes on button press in pushbutton170.
function pushbutton170_Callback(hObject, eventdata, handles)
cat_io_senderrormail;


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
% --- Executes during object creation, after setting all properties.

% Determine the selected data set.
if get(hObject,'Value') == 2
    if exist(fullfile(spm('dir'),'toolbox','TFCE'))
        % call TFCE toolbox 
        spm_TFCE;
    else % install TFCE toolbox
        d0 = spm('Dir');
        d = fullfile(spm('Dir'),'toolbox'); 
        s = unzip('http://www.neuro.uni-jena.de/tfce/tfce_latest.zip', d);
        addpath(d0);
        rehash
        rehash toolboxcache;
        toolbox_path_cache
        eval(['spm fmri;clear cat_version;spm_cat12']);
    end
end


% --- Executes when CAT is resized.
function CAT_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to CAT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when uipanel27 is resized.
function uipanel27_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to uipanel27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu15.
function popupmenu15_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu15 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu15


% --- Executes during object creation, after setting all properties.
function popupmenu15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu16.
function popupmenu16_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu16 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu16


% --- Executes during object creation, after setting all properties.
function popupmenu16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton8.
function pushbutton8_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu17.
function popupmenu17_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu17 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu17


% --- Executes during object creation, after setting all properties.
function popupmenu17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton174.
function pushbutton174_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton174 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cat_stat_analyze_ROIs;

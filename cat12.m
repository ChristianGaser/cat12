function varargout = cat12(varargin)
% ______________________________________________________________________
% CAT12 Toolbox wrapper to call CAT functions.
% 
%   cat12 
%     .. start with CAT default parameter file
%   cat12('gui')
%     .. start with default file of another species (in development)
%   cat12(species) 
%     .. start with default file of another species (in development)
%        species = ['oldwoldmonkey'|'newwoldmonkey'|'greaterape'|'lesserape']
%   cat12('mypath/cat_defaults_mydefaults') 
%     .. start CAT with another default parameter file
% ______________________________________________________________________
% Christian Gaser
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

% Last Modified by GUIDE v2.5 21-Nov-2015 23:32:29

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


% --- Outputs from this function are returned to the command line.
function varargout = cat12_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% This creates the 'background' image
ha = axes('units','normalized','position',[0 0.87 1 0.13]);
uistack(ha,'bottom');
I = imread(fullfile(spm('dir'),'toolbox','cat12','html','images','contact.jpg'));
hi = imagesc(I);
text(80,140,'Computational Anatomy Toolbox','Color',[1 1 1],'Fontsize',20,'Fontweight','bold');
set(ha,'handlevisibility','off','visible','off');


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
for i=1:size(P,1)
    swd      = spm_file(P(i,:),'fpath');
    load(fullfile(swd,'SPM.mat'));
    SPM.swd  = swd;

    fsavgDir = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces');
    
    % check that folder exist and number of vertices fits
    if exist(fsavgDir) == 7 & SPM.xY.VY(1).dim(1) == 163842
      [pp,ff]   = spm_fileparts(SPM.xY.VY(1).fname);
      
      % find lh|rh string
      hemi_ind = [];
      hemi_ind = [hemi_ind strfind(ff,'lh')];
      hemi_ind = [hemi_ind strfind(ff,'rh')];
      hemi = ff(hemi_ind:hemi_ind+1);
      if ~isempty(hemi)
        SPM.xY.VY(1).private.private.metadata = struct('name','Name','value',fullfile(fsavgDir,[hemi '.central.freesurfer.gii']));
        
        % remove memory demanding faces and vertices which are not necessary
        for i=1:length(SPM.xY.VY)
          SPM.xY.VY(i).private.faces = [];
          SPM.xY.VY(i).private.vertices = [];
        end
        
        save(fullfile(swd,'SPM.mat'),'SPM');
      end
    end

    spm_spm(SPM);
    
    % workaround to use fsaverage surface for display
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

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
cat_vol_slice_overlay;

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
F = spm_figure('FindWin','Menu');
% close SPM windows, if no Menu window exist
if isempty(F)
  spm('Quit')
end
close(gcf);


% --- Executes on selection change in popupmenu9.
function popupmenu9_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu9 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu9


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

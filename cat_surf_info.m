function [varargout] = cat_surf_info(P,read)
% ______________________________________________________________________
% Extact surface information from filename.
%
% sinfo = cat_surf_info(P,readsurf)
%
% sinfo(i). 
%   pp        .. filepath
%   ff        .. filename
%   ee        .. filetype
%   exist     .. exist file?
%   ftype     .. filetype [0=no surface,1=gifti,2=freesurfer]
% 
%   statready .. ready for statistik (^s#mm.*.gii) [0|1]
%   side      .. hemishphere [lh|rh] 
%   datatype  .. [0=nosurf/file|1=mesh|2=data|3=surf]
%                only with readsurf==1 and with surf=mesh+data
%   dataname  .. datafieldname [central|thickness|s3thickness|myclalc...]
%   texture   .. textureclass [central|sphere|thickness|...]
%   resampled .. meshspace [0|1] 
%   template  .. template or individual mesh [0|1] 
%   name      .. name of the dataset
% ______________________________________________________________________
% Robert Dahnke
% $Id$

%#ok<*RGXP1>

  if isempty(P) && nargout>0, varargout{1} = {}; return; end
  
  if nargin<2, read = 0; end

  P = cellstr(P);
  
  sinfo = struct(...
    'fname','',...      % full filename
    'pp','',...         % filepath
    'ff','',...         % filename
    'ee','',...         % filetype
    'exist','',...      % exist
    'fdata','',...      % datainfo (filesize)
    'ftype','',...      % filetype [0=no surface,1=gifti,2=freesurfer]
    ...
    'statready',0,...   % ready for statistik (^s#mm.*.gii)
    'side','',...       % hemishphere
    'name','',...       % subject/template name
    'datatype','',...   % datatype [0=nosurf/file|1=mesh|2=data|3=surf] with surf=mesh+data
    'dataname','',...   % datafieldname [central|thickness|s3thickness...]
    'texture','',...    % textureclass [central|sphere|thickness|...]
    'resampled','',...  % dataspace
    'template','',...   % individual surface or tempalte
    'nvertices',[],...  % number vertices
    'nfaces',[],...     % number faces
    ...
    'Pmesh','',...      % meshfile
    'Psphere','',...    % meshfile
    'Pdefects','',...   % meshfile
    'Pdata',''...       % datafile
  );

  if isempty(P), return; end
  
  for i=1:numel(P)
    [pp,ff,ee] = spm_fileparts(P{i});
    sinfo(i).fdata = dir(P{i});
    
    sinfo(i).fname = P{i};
    sinfo(i).exist = exist(P{i},'file'); 
    sinfo(i).pp = pp;
    switch ee
      case {'.xml','.txt','.html','.csv'}
        sinfo(i).ff = ff;
        sinfo(i).ee = ee;
        sinfo(i).ftype = 0;
        continue
      case '.gii'
        sinfo(i).ff = ff;
        sinfo(i).ee = ee;
        sinfo(i).ftype = 1;
        if sinfo(i).exist && read
          S = gifti(P{i});
        end
      otherwise
        sinfo(i).ff = [ff ee];
        sinfo(i).ee = '';
        sinfo(i).ftype = 0;
        if sinfo(i).exist && read
          clear S; 
          try
            S = cat_io_FreeSurfer('read_surf',P{1}); 
          end
          try
            S.cdata = cat_io_FreeSurfer('read_surf_data',P{1}); 
          end
        end
        if exist('S','var')
          sinfo(i).ftype = 2;
        end
    end
    
    % name
    [tmp,noname,name] = spm_fileparts(sinfo(i).ff); 
    if isempty(name), name=''; else name = name(2:end); end
    sinfo(i).name = name;
   
    sinfo(i).statready = ~isempty(regexp(noname,'^s(?<smooth>\d+)mm\..*')); 
    
    % side
    if     strfind(noname,'lh'), sinfo(i).side='lh'; sidei = strfind(noname,'lh.');
    elseif strfind(noname,'rh'), sinfo(i).side='rh'; sidei = strfind(noname,'rh.');
    else
      sinfo(i).side = spm_input('Hemisphere',1,'lh|rh');
      sidei = strfind(noname,[sinfo(i).side '.']);
    end
    if sidei>0
      sinfo(i).preside = noname(1:sidei-1);
      sinfo(i).posside = noname(sidei+numel(sinfo(i).side)+1:end);
    else
      sinfo(i).preside = '';
      sinfo(i).posside = noname;
    end
    
    % datatype
    if sinfo(i).exist && read
      switch num2str([isfield(S,'vertices'),isfield(S,'cdata')],'%d%d')
        case '00',  sinfo(i).datatype  = 0;
        case '01',  sinfo(i).datatype  = 1;
        case '10',  sinfo(i).datatype  = 2;
        case '11',  sinfo(i).datatype  = 3;
      end
    else
      sinfo(i).datatype = -1;
    end
    
    % dataname
    sinfo(i).dataname  = strrep(sinfo(i).posside,'.resampled','');
    
    % special datatypes
    FN = {'thickness','central','sphere','defects','gyrification','sqrtsulc','frac',...
          'gyruswidth','gyruswidthWM','sulcuswidth','WMdepth','CSFdepth','GWMdepth',...
          'depthWM','depthGWM','depthCSF','depthWMg',...
          'hulldist'};
    sinfo(i).texture = '';
    for fi=1:numel(FN)
      if strfind(sinfo(i).dataname,FN{fi}), sinfo(i).texture = FN{fi}; end
    end   
        
    % template
    sinfo(i).template  = ~isempty(strfind(lower(sinfo(i).ff),'.template')); 

    % resampled
    sinfo(i).resampled = ~isempty(strfind(sinfo(i).posside,'.resampled'));
    if sinfo(i).template,  sinfo(i).resampled = 1; end
    
    
    
    
    % find Mesh and Data Files
    %  -----------------------------------------------------------------
    sinfo(i).Pmesh = '';
    sinfo(i).Pdata = '';
    % here we know that the gifti is a surf
    if sinfo(i).statready 
      sinfo(i).Pmesh = sinfo(i).fname;
      sinfo(i).Pdata = sinfo(i).fname;
    end
    % if we have read the gifti than we can check for the fields
    if isempty(sinfo(i).Pmesh) && sinfo(i).exist && read && isfield(S,'vertices')
      sinfo(i).Pmesh = sinfo(i).fname; 
    end
    if isempty(sinfo(i).Pdata) && sinfo(i).exist && read && isfield(S,'cdata')
      sinfo(i).Pdata = sinfo(i).fname;
    end
    % if the dataname is central we got a mesh or surf datafile
    if isempty(sinfo(i).Pdata) 
      switch sinfo(i).texture
        case {'defects'} % surf
          sinfo(i).Pmesh = sinfo(i).fname;
          sinfo(i).Pdata = sinfo(i).fname;
        case {'central','sphere'} % only mesh
          sinfo(i).Pmesh = sinfo(i).fname;
          sinfo(i).Pdata = '';
        case {'thickness','gyrification','frac','logsulc','GWMdepth','WMdepth','CSFdepth',...
             'depthWM','depthGWM','depthCSF','depthWMg',...
             'gyruswidth','gyruswidthWM','sulcuswidth'} % only thickness
          sinfo(i).Pdata = sinfo(i).fname;
      end
    end
    % if we still dont know what kind of datafile, we can try to find a
    % mesh surface
    if isempty(sinfo(i).Pmesh) 
      % template mesh handling !!!
      Pmesh = char(cat_surf_rename(sinfo(i),'dataname','central','ee','.gii'));
      if exist(Pmesh,'file')
        sinfo(i).Pmesh = Pmesh;
        sinfo(i).Pdata = sinfo(i).fname;
      end
    end
    % if we got still no mesh than we can find an average mesh
    % ...
    if isempty(sinfo(i).Pmesh) %&& sinfo(i).ftype==1
      sinfo(i).Pmesh = ...
        fullfile(spm('dir'),'toolbox','cat12','templates_surfaces',[sinfo(i).side '.central.freesurfer.gii']);
      sinfo(i).Pdata = sinfo(i).fname;
    end
    
    [ppm,ffm,eem] = fileparts(sinfo(i).Pmesh);
    sinfo(i).Psphere  = fullfile(ppm,strrep([ffm eem],'.central.','.sphere.'));
    sinfo(i).Pdefects = fullfile(ppm,strrep([ffm eem],'.central.','.defects.'));
    if ~exist(sinfo(i).Psphere ,'file'), sinfo(i).Psphere  = ''; end
    if ~exist(sinfo(i).Pdefects,'file'), sinfo(i).Pdefects = ''; end

    
    if sinfo(i).exist && read
      if isfield(S,'vertices'), 
        sinfo(i).nvertices = size(S.vertices,1);
      else
        if ~isempty(sinfo(i).Pmesh)
          S2 = gifti(sinfo(i).Pmesh);
          if isfield(S2,'vertices'), S.vertices = S2.vertices; else S.vertices = []; end
          if isfield(S2,'faces'),    S.faces    = S2.faces;    else S.faces = []; end
        end
        sinfo(i).nvertices = size(S.vertices,1);
      end
      if isfield(S,'faces'),    sinfo(i).nfaces    = size(S.faces,1); end
      if isfield(S,'cdata'),    sinfo(i).ncdata    = size(S.cdata,1); end
    end
    
    sinfo(i).catxml = fullfile(pp,['cat_' sinfo(i).name '*.xml']);
    if ~exist(sinfo(i).catxml,'file'), sinfo(i).catxml = ''; end 
    
    if nargout>1
      varargout{2}{i} = S; 
    else
      clear S
    end
  end
  varargout{1} = sinfo; 
end
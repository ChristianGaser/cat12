function sinfo = vbm_surf_info(P,read)
% ______________________________________________________________________
% Extact surface information from filename.
%
% sinfo = vbm_surf_info(P,readsurf)
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

  if nargin<2, read = 0; end

  P = cellstr(P);
  
  sinfo = struct(...
    'fname','',...      % full filename
    'pp','',...         % filepath
    'ff','',...         % filename
    'ee','',...         % filetype
    'exist','',...      % exist
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
    'Pmesh','',...      % meshfile file
    'Pvalue','',...     % texture file
    'nvertices',[],...  % number vertices
    'nfaces',[]...      % number faces
  );

  if isempty(P), return; end
    
  for i=1:numel(P)
    [pp,ff,ee] = spm_fileparts(P{i});

    sinfo(i).fname = P{i};
    sinfo(i).exist = exist(P{i},'file'); 
    sinfo(i).pp = pp;
    switch ee
      case {'.xml','.txt','.html','.csv'}
        sinfo(i).ff = ff;
        sinfo(i).ee = ee;
        sinfo(i).surform = 0;
        continue
      case '.gii'
        sinfo(i).ff = ff;
        sinfo(i).ee = ee;
        sinfo(i).surform = 1;
        if sinfo(i).exist && read
          S = gifti(P{i});
        end
      otherwise
        sinfo(i).ff = [ff ee];
        sinfo(i).ee = '';
        sinfo(i).surform = 2;
        if sinfo(i).exist && read
          clear S; 
          try
            S = vbm_io_FreeSurfer('read_surf',P{1}); 
          end
          try
            S.cdata = vbm_io_FreeSurfer('read_surf_data',P{1}); 
          end
        end
    end
    
    % name
    [tmp,noname,name] = spm_fileparts(sinfo(i).ff); 
    if isempty(name), name=''; else name = name(2:end); end
    sinfo(i).name = name;
   
    sinfo(i).statready = ~isempty(regexp(noname,'^s(?<smooth>\d+)mm\..*')); 
    
    % side
    if     strfind(noname,'lh.'), sinfo(i).side='lh'; sidei = strfind(noname,'lh.');
    elseif strfind(noname,'rh.'), sinfo(i).side='rh'; sidei = strfind(noname,'rh.');
    else                          sinfo(i).side='';   sidei = 0;
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
      switch char([isfield('vertices',S),isfield('cdata',S)])
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
    FN = {'thickness','central','sphere','defects','curv','frac'};
    for fi=1:numel(FN)
      if strfind(sinfo(i).dataname,FN{fi}), sinfo(i).texture = FN{fi}; end
    end   
        
    % template
    sinfo(i).template  = ~isempty(regexp(noname,'.*\.template\..*')); 

    % resampled
    sinfo(i).resampled = ~isempty(regexp(noname,'.*\.resampled\..*'));
    if sinfo(i).template,  sinfo(i).resampled = 1; end
    
    if sinfo(i).template
    % template mesh
      sinfo(i).Pmesh = char(vbm_surf_rename(sinfo(i),'dataname','central'));
    else
      sinfo(i).Pmesh = char(vbm_surf_rename(sinfo(i),'dataname','central'));
    end
    
    if ~strcmp(sinfo(i).Pmesh,sinfo(i).fname);
      sinfo(i).Pdata = sinfo(i).fname;
    else
      sinfo(i).Pdata = '';
    end
    % if resampled > matching avg-meshs?
  end
  
end
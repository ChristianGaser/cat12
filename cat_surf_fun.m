function varargout = cat_surf_fun(action,S,varargin)
% Function to collect surface functions.
% 
% varargout = cat_surf_fun(action,S,varargin)
%
% * Distance estimation:
%   D   = cat_surf_fun('dist',S);     % Estimates the distance between the 
%                                       local vertices of the faces of S.
%   D   = cat_surf_fun('vdist',S);    % Estimates the distance of each 
%                                       voxel to the surface
%
%   A   = cat_surf_fun('area',S);     % estimate surface area (faces)
%
% * Data mapping:
%   V   = cat_surf_fun(S,F);          % map facedata to vertices
%   C   = cat_surf_fun('cdatamapping',S1,S2,cdata[,opt]); 
%                                     % map texture cdata from S2 to S1.
%
% * Surface (data) rendering:
%   V   = cat_surf_fun('surf2vol',S,varargin{1}); 
%                                     % render surface (data) in a volume
%
% * Surface modification:
%   HS  = cat_surf_fun('hull',S);     % estimate hull surface
%   HS  = cat_surf_fun('core',S);     % estimate core surface
%   IS  = cat_surf_fun('inner',S,T);  % estimate inner surface
%   OS  = cat_surf_fun('outer',S,T);  % estimate outer surface
%
% * Surface modification 2:
%   cat_surf_saveICO(S,T,Pcentral,subdir,writeTfs,C)
%   Creates and save the white and pial surfaces based on the displacement 
%   by the half thickness along the surface normals and use the white and 
%   pial surfaces to create the layer4 surface.
%   Saves also the thickness file "pbt" and "thickness".
%
% * Surface area mapping
%   - Estimate nearest neighbor mapping between two surfaces S1 and S2
%      edgemap = cat_surf_fun('createEdgemap',S1,S2);
%   - Applay mapping 
%      cdata2 = cat_surf_fun('useEdgemap',cdata,edgemap);
%
% * Other helping functions:
%   E   = cat_surf_fun('graph2edge',T); % get edges of a triangulation T
%
% * Test functions
%   'cdatamappingtst'
%
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
% ______________________________________________________________________
% $Id$ 

%#ok<*ASGLU>

  switch lower(action)
    case {'dist','distance'}
      varargout{1} = cat_surf_dist(S);
    case 'surf2surf'
      % create mapping between similar surfaces for area projection
      varargout{1} = cat_surf_surf2surf(S,varargin{1});
    case 'area'
       % area estiamtion 
       varargout{1} = cat_surf_area(S);
    case {'smoothcdata','smoothtexture'}
      switch nargin
        case 2
          if isfield(S,'cdata'), varargout{1} = cat_surf_smoothtexture(S,S.cdata,1); end
        case 3, varargout{1} = cat_surf_smoothtexture(S,varargin{1});
        case 4, varargout{1} = cat_surf_smoothtexture(S,varargin{1},varargin{2});
        case 5, varargout{1} = cat_surf_smoothtexture(S,varargin{1},varargin{2},varargin{3});
      end
    case 'maparea'    
      % do the final area projection
      if nargout==1, varargout{1} = cat_surf_maparea(S,varargin{1}); end
      if nargout==2, [varargout{1},varargout{2}] = cat_surf_maparea(S,varargin{1}); end
    case 'hull'
      if nargout==1, varargout{1} = cat_surf_hull(S); end
      if nargout==2, [varargout{1},varargout{2}] = cat_surf_hull(S); end
    case 'core'
      if nargout==1, varargout{1} = cat_surf_core(S,varargin{1}); end
      if nargout==2, [varargout{1},varargout{2}] = cat_surf_core(S,varargin{1}); end
    case {'tfs','tmin','tmax'}
      if numel(varargin)==1
        if nargout==1
          varargout{1} = cat_surf_thickness(action,S,varargin{1}); 
        else
          cat_surf_thickness(action,S,varargin{1});
        end
      else
        if nargout==1
          varargout{1} = cat_surf_thickness(action,S); 
        else
          cat_surf_thickness(action,S);
        end
      end
    case {'inner','outer','white','pial','innervar','outervar','whitevar','pialvar'}
      if numel(varargin)==1
        switch nargout % surface & texture input
          case 0, cat_surf_GMboundarySurface(action,S,varargin{1});
          case 1, varargout{1} = cat_surf_GMboundarySurface(action,S,varargin{1}); 
          case 2, [varargout{1},varargout{2}] = cat_surf_GMboundarySurface(action,S,varargin{1}); 
        end
      else % file input
        switch nargout
          case 0, cat_surf_GMboundarySurface(action,S);
          case 1, varargout{1} = cat_surf_GMboundarySurface(action,S); 
          case 2, [varargout{1},varargout{2}] = cat_surf_GMboundarySurface(action,S); 
        end
      end
    case 'saveico'
      switch nargin
        case 5, cat_surf_saveICO(S,varargin{1},varargin{2},varargin{3});
        case 6, cat_surf_saveICO(S,varargin{1},varargin{2},varargin{3},varargin{4});
        case 7, cat_surf_saveICO(S,varargin{1},varargin{2},varargin{3},varargin{4},varargin{5});
        case 8, cat_surf_saveICO(S,varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},varargin{6});
      end
    case 'collisioncorrection'
      switch nargin
        case 4
          [varargout{1},varargout{2}] = cat_surf_collision_correction(S,varargin{1},varargin{2});
        case 5
          [varargout{1},varargout{2}] = cat_surf_collision_correction(S,varargin{1},varargin{2},varargin{3});
        case 6
          [varargout{1},varargout{2}] = cat_surf_collision_correction(S,varargin{1},varargin{2},varargin{3},varargin{4});
      end
    case 'vdist'
      [varargout{1},varargout{2}] = cat_surf_vdist(S,varargin);
    case 'surf2vol'
      if nargin>2
        [varargout{1},varargout{2},varargout{3}] = cat_surf_surf2vol(S,varargin{1});
      else
        [varargout{1},varargout{2},varargout{3}] = cat_surf_surf2vol(S);
      end
    case 'graph2edge'
      switch nargout
        case 0, cat_surf_edges(S); 
        case 1, varargout{1} = cat_surf_edges(S); 
        case 2, [varargout{1},varargout{2}] = cat_surf_edges(S); 
      end
    case 'cdatamappingtst'
      cat_surf_cdatamappingtst;
    case 'createedgemap'
      varargout{1} = cat_surf_surf2surf(S,varargin{1});
    case 'useedgemap'
      varargout{1} = cat_surf_maparea(S,varargin{1});
    case 'gmv'
      varargout{1} = cat_surf_gmv(S,varargin{1});
    case 'cdatamapping' 
      if nargin<3, varargin{3} = ''; end
      if nargin<4, varargin{4} = struct(); end
      if nargout>1
        [varargout{1},varargout{2}] = cat_surf_cdatamapping(S,varargin{1},varargin{2},varargin{3});
      else
        varargout{1} = cat_surf_cdatamapping(S,varargin{1},varargin{2},varargin{3});
      end 
    otherwise
      error('Unknow action "%s"!\n',action); 
  end
    
end

%% area mapping concept
%  ------------------------------------------------------------------------
%  We need two functions, one that create the mapping between two very 
%  close surfaces (cat_surf_surf2surf) and another one (cat_surf_maparea) 
%  that finally performs the mapping. 
%  ------------------------------------------------------------------------
%  Robert Dahnke 2019/04
function edgemap = cat_surf_surf2surf(S1,S2,normalize)
%  ------------------------------------------------------------------------
%  Create mapping by nearest neighbor.
%  ------------------------------------------------------------------------
%  Robert Dahnke 2019/04
  if ~exist('normalize','var'), normalize=1; end

  % 0. Normalize input
  if normalize
    S1.vertices = S1.vertices/max(S1.vertices(:));
    S2.vertices = S2.vertices/max(S2.vertices(:)) * 0.98; % need slight difference for Delaunay! 
  end
  
  % 1. Delaunay triangulation for fast search
  D1 = delaunayn(double(S1.vertices)); 
  D2 = delaunayn(double(S2.vertices)); 
 
  % 2. find minimal relation between the vertices of S1 and S2 as nearest neighbor
  E1(:,1) = 1:size(S1.vertices,1);
  E2(:,2) = 1:size(S2.vertices,1);
  E1(:,2) = dsearchn(double(S2.vertices),D2,double(S1.vertices)); 
  E2(:,1) = dsearchn(double(S1.vertices),D1,double(S2.vertices)); 
  E = unique([E1;E2],'rows'); 
  
  % 3. estimate edge length as weighting function 
  EL = sum( ( S1.vertices(E(:,1),:) - S2.vertices(E(:,2),:) ).^2 , 2) .^ 0.5; 
  
  % 4. Estimate the number of all by c-function
  %    - outgoing edges of S1 (connections from v element of S1) and
  %    - incoming edges of S2 (connections from v element of S2) 
  nE1 = zeros(max(E(:,1)),1,'single'); EL1 = zeros(max(E(:,1)),1,'single');
  nE2 = zeros(max(E(:,2)),1,'single'); EL2 = zeros(max(E(:,2)),1,'single');
  for i=1:size(E,1)
    nE1( E(i,1) ) =  nE1( E(i,1) ) + 1;
    EL1( E(i,1) ) =  EL1( E(i,1) ) + EL(i);
    nE2( E(i,2) ) =  nE2( E(i,2) ) + 1; 
    EL2( E(i,2) ) =  EL2( E(i,2) ) + EL(i);
  end
  
  % 5. Create a weighting function to project data from Si2St and St2Si.
  edgemap.edges     = E;
  edgemap.dist      = EL;
  edgemap.nvertices = [size(S1.vertices,1),size(S2.vertices,1)];
  edgemap.nforward  = 1  ./ nE1(E(:,1)); 
  edgemap.nbackward = 1  ./ nE2(E(:,2));
  edgemap.dforward  = EL ./ EL1(E(:,1)); 
  edgemap.dbackward = EL ./ EL2(E(:,2)); 
end
function varargout = cat_surf_maparea(varargin)
%  Apply graph-based mapping
%  ------------------------------------------------------------------------
%  use a c-function to process cat_surf_surf2surf mapping function 
%
%   cdata = cat_surf_maparea(cdatain,edgemap[,weighting])
%  
%   cdata     .. texture values at the output surface
%   edgemap   .. mapping structure between two surfaces
%   direction .. direction of cdata mapping if both surfaces have cdata 
%                with direction: 'forward' == '', 'backward' == 'invers'
%   weighting .. type of weighting: 
%                 'num'  .. by number of vertices
%                 'dist' .. by distance to the vertices (default)
%
%  ------------------------------------------------------------------------
%  Robert Dahnke 2019/04

  cdata     = varargin{1}; 
  edgemap = varargin{2};
  if nargin>2 
    dir = varargin{3}; 
  else
    dir = '';
  end
  switch dir
    case {'','forward'},        idir = 0; 
    case {'invers','backward'}, idir = 1;
    otherwise, error('Unkown mapping direction %s.\n',dir);
  end

  varargout{1} = cat_surf_edgemap(edgemap,cdata,idir);
end
function cdata2 = cat_surf_edgemap(edgemap,cdata,idir)
  if idir==0
    cdata2 = zeros(edgemap.nvertices(2),1,'single');
    for i=1:size(edgemap.edges,1)
      cdata2(edgemap.edges(i,2)) = cdata2(edgemap.edges(i,2)) + ...
        cdata(edgemap.edges(i,1)) * edgemap.dforward(i);
    end
  else
    cdata2 = zeros(edgemap.nvertices(1),1,'single');
    for i=1:size(edgemap.edges,1)
      cdata2(edgemap.edges(i,1)) =  cdata2(edgemap.edges(i,1)) + ...
        cdata2(edgemap.edges(i,2)) * edgemap.dbackward(i);
    end    
  end
end

function VV = cat_surf_gmv(IS,OS)
%%
  V = double([IS.vertices;OS.vertices]);

  % create Delaunay triangulation 
  D  = delaunayn(V); 
  
  % classify and remove non GM tetraeder
  DS = D>size(IS.vertices,1);
  D( sum(DS,2)==0 | sum(DS,2)==4 ,:)  = [];  
  clear DS;
  
  % estimate tetraeder volume
  DV = tetraedervolume(D,V);
  
  %% map volume to faces
  VV = zeros(size(IS.vertices,1),1,'single'); DS = VV;
  DF = D; DF(DF>size(IS.vertices,1)) = DF(DF>size(IS.vertices,1)) - size(IS.vertices,1);
  for i=1:numel(DF)
    VV(DF(i)) = VV(DF(i)) + DV( D(i) );  
    DS(DF(i)) = DS(DF(i)) + 1; 
  end
  VV = VV./max(1,DS); 
 
end
function DV = tetraedervolume(D,V)
% estimate tetraeder volume by the Cayley-Menger determinant

  % edgelength
  r = sum( ( V(D(:,1),:) - V(D(:,2),:) ).^2 , 2).^0.5; 
  p = sum( ( V(D(:,2),:) - V(D(:,3),:) ).^2 , 2).^0.5; 
  q = sum( ( V(D(:,3),:) - V(D(:,1),:) ).^2 , 2).^0.5; 
  a = sum( ( V(D(:,1),:) - V(D(:,4),:) ).^2 , 2).^0.5; 
  b = sum( ( V(D(:,2),:) - V(D(:,4),:) ).^2 , 2).^0.5; 
  c = sum( ( V(D(:,3),:) - V(D(:,4),:) ).^2 , 2).^0.5; 

  % volume
  DV = zeros(size(D,1),1);
  for i=1:size(D,1)
    DM = [ 0    r(i) q(i) a(i) 1;...
           r(i) 0    p(i) b(i) 1;...
           q(i) p(i) 0    c(i) 1;...
           a(i) b(i) c(i) 0    1;...
           1    1    1    1    0];
    DV(i) = sqrt(det( DM .^2 ) / 288);
  end
end
function varargout = cat_surf_GMboundarySurface(type,varargin)
  if strfind(type,'var')
    varout=1; type = strrep(type,'var',''); 
  else
    varout=0;
  end
  switch type
    case {'white','inner'}, direction = -0.5;
    case {'pial' ,'outer'}, direction =  0.5;
  end
  
  if nargin==2
    %% use filenames
    [pp,ff,ee] = spm_fileparts(varargin{1}); 
    
    if strcmp(ee,'')
      Praw = cat_io_FreeSurfer('fs2gii',varargin{1}); 
      Praw = Praw{1};
    else
      Praw   = varargin{1};
    end
    Pthick = cat_io_strrep(Praw,{'central','.gii'},{'pbt',''});
    Ptype  = cat_io_strrep(Praw,'central',type);
    
    cmd = sprintf('CAT_Central2Pial "%s" "%s" "%s" %0.2f',Praw,Pthick,Ptype,direction); 
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,1);

    if strcmp(ee,'')
      Ptype = cat_io_FreeSurfer('gii2fs',Ptype); 
    end
    
    % filename
    if varout
      % load surface 
      varargout{1} = gifti(Ptype); 

      % delete temp files
      delete(Ptype);
    else
      varargout{1} = Ptype; 
    end
  else
    % write temp files ...
    Praw   = 'central.';
    Pthick = strrep(Praw,'central','pbt');
    Ptype  = strrep(Praw,'central',type);
   
    cmd = sprintf('CAT_Central2Pial "%s" "%s" %0.2f',Praw,Pthick,direction); 
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,1);
    
    % load surface 
    varargout{1} = gifti(Ptype); 
    
    % delete temp files
    delete(Praw,Pthick,Ptype);
  end
end

function cat_surf_cdatamappingtst

%% Testdata
   Psubcentral  = ['/Volumes/vbmDB/MRData/vbm12tst/results/deffiles/cg_vbm_defaults_template/template_NKI/'...
     'surf/lh.central.NKI_HC_NKI_1013090_T1_SD000000-RS00.gii'];
   PsubsphereA  = strrep(Psubcentral,'central','sphere.reg');              
   %Psubthick    = strrep(strrep(Psubcentral,'central','pbt'),'.gii','');               
   Psubthickres = strrep(strrep(Psubcentral,'central','thickness.resampled'),'lh.','s15mm.lh.'); 
   Psubtmp      = strrep(Psubcentral,'central','tmp'); 
   Pavgtmp      = strrep(strrep(Psubcentral,'central','tmp.resampled'),'lh.','s15mm.lh.'); 
 
   %Pavgcentral  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces','lh.central.freesurfer.gii'));
   PavgsphereA  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces','lh.sphere.freesurfer.gii'); 
   PavgDKT40    = fullfile(spm('dir'),'toolbox','cat12','atlases_surfaces','lh.aparc_DKT40JT.freesurfer.annot');
   
%% Test 1 - avg2sub - ok
   Ssub = gifti(PsubsphereA);
   Savg = gifti(PavgsphereA); 
   [vertices, label, colortable]  = cat_io_FreeSurfer('read_annotation',PavgDKT40); 
   Savg.cdata = label; 
   
   S3 = gifti(Psubcentral); 
   S3.cdata = cat_surf_fun('cdatamapping',Ssub,Savg,'nearest');
   save(gifti(S3),Psubtmp);
   
%% Test2 - sub2avg - ok
   Savg = gifti(PavgsphereA); 
   Ssub = gifti(PsubsphereA);
   %Ssub.cdata = cat_io_FreeSurfer('read_surf_data',Psubthick); 
   Ssub.cdata = cat_surf_fun('area',gifti(Psubcentral));
   
   S3 = gifti(Psubthickres); 
   mapping = {'directed'}; %,'undirected'}; %'nearest',
   for mi = 1:numel(mapping)
     S3.cdata  = cat_surf_fun('cdatamapping',Savg,Ssub,mapping{mi},1);
     S3.cdata  = spm_mesh_smooth(struct('vertices',S3.vertices,'faces',S3.faces),double(S3.cdata'),5);
     fprintf('mapping = %10s: A(sub) = %0.2f vs. A(avg) = %0.2f\n',mapping{mi},sum(Ssub.cdata(:)),sum(S3.cdata(:))); 
     save(gifti(S3),Pavgtmp); cat_surf_display(Pavgtmp)
   end
   
end

% nearest connection between to surfaces
function varargout = cat_surf_cdatamapping(S1,S2,cdata,opt) 
  if ischar(S1), S1 = gifti(S1); end
  if ischar(S2), S2 = gifti(S2); end
  if ischar(cdata)
    Pcdata = cdata;
    [pp,ff,ee] = spm_fileparts(cdata); 
    switch ee
      case '.annot'
        [vertices, cdata]  = cat_io_FreeSurfer('read_annotation',Pcdata); 
        clear vertices
      case '.gii'
        Scdata = gifti(S2); 
        if isfield(Scdata,'cdata')
          cdata = SX.cdata;
        else
          error('cat_surf_fun:cdatamapping:noTexture','No texture found in "%s"!\n',Pcdata);
        end
      otherwise
        cdata =  cat_io_FreeSurfer('read_surf_data',Pcdata);   
    end
  end
  
  if ~exist('cdata','var') || isempty(cdata)
    if isfield(S2,'cdata'), cdata = S2.cdata; end
  end
  
  if ~exist('opt','var'), opt = struct(); end
  def.method = 'nearest';
  def.verb   = 0; 
  def.smooth = 0; 
  opt        = cat_io_checkinopt(opt,def);
  
  if opt.verb, stime1 = cat_io_cmd(sprintf('Data-mapping (%s)',method)); fprintf('\n'); end
  
  % prepare vertices
  S1.vertices = S1.vertices ./ repmat(max(S1.vertices),size(S1.vertices,1),1)*1.1; % *100 
  S2.vertices = S2.vertices ./ repmat(max(S2.vertices),size(S2.vertices,1),1); 
  verticesS1  = double(S1.vertices - repmat(mean(S1.vertices),size(S1.vertices,1),1)); 
  verticesS2  = double(S2.vertices - repmat(mean(S2.vertices),size(S2.vertices,1),1)); 
  
  
  % estimate mapping
  switch opt.method
    case {'nearest'}
      [varargout{2},varargout{3}] = dsearchn([verticesS2;inf(1,3)],double([S2.faces ones(size(S2.faces,1),1)*(size(S2.vertices,1)+1)]),verticesS1);
      varargout{1} = cdata(varargout{2}); 
    case {'undirected','directed'}
      %% use the surface as delauny graph
      switch opt.method 
        case 'directed'
          if opt.verb,  stime = cat_io_cmd('  Edge-Estimation (Nearest)','g5',''); end
          nextS2fromS1 = dsearchn([verticesS2;inf(1,3)],double([S2.faces ones(size(S2.faces,1),1)*(size(S2.vertices,1)+1)]),verticesS1);
          nextS1fromS2 = dsearchn([verticesS1;inf(1,3)],double([S1.faces ones(size(S1.faces,1),1)*(size(S1.vertices,1)+1)]),verticesS2);
          tmp = nextS1fromS2; nextS1fromS2 = nextS2fromS1; nextS2fromS1 = tmp;
          nearestedges = [ (1:numel(nextS2fromS1))', nextS2fromS1; nextS1fromS2 , (1:numel(nextS1fromS2))' ]; 
          nearestedges = unique(nearestedges,'rows');
        case 'undirected'
          if opt.verb,  stime = cat_io_cmd('  Edge-Estimation (Delaunay','g5',''); end
          % nearest is required too
          nextS2fromS1 = dsearchn([verticesS2;inf(1,3)],double([S2.faces ones(size(S2.faces,1),1)*(size(S2.vertices,1)+1)]),verticesS1);
          nextS1fromS2 = dsearchn([verticesS1;inf(1,3)],double([S1.faces ones(size(S1.faces,1),1)*(size(S1.vertices,1)+1)]),verticesS2);
          tmp = nextS1fromS2; nextS1fromS2 = nextS2fromS1; nextS2fromS1 = tmp;
          nearestedges  = [ (1:numel(nextS2fromS1))', nextS2fromS1; nextS1fromS2 , (1:numel(nextS1fromS2))' ]; 
          nearestedges1 = unique(nearestedges,'rows');
          % delauany
          triangulation = delaunayn([verticesS2;verticesS1]);              % delaunay triangulation
          nearestedges  = cat_surf_fun('graph2edge',triangulation);        % get edges 
          nearestedges(sum(nearestedges<=size(verticesS2,1),2)~=1,:)=[];   % only edges between S1 and S2
          nearestedges(:,2) = nearestedges(:,2) - size(verticesS2,1); 
          nearestedges = unique([nearestedges;nearestedges1],'rows');
      end
      if opt.verb, stime = cat_io_cmd('  Weighting','g5','',1,stime); end
      
      if 0
        %% my little testset
        nextS1fromS2 = [1; 1; 3; 4; 4; 4; 5; 5]; 
        nextS2fromS1 = [1; 3; 3; 5; 8; 8];
        cdata        = [1 1 1 1 1 1]';
        nearestedges = [ (1:numel(nextS2fromS1))', nextS2fromS1; nextS1fromS2 , (1:numel(nextS1fromS2))' ]; 
        nearestedges = unique(nearestedges,'rows');
      end
      
      
      %% simplify edges 1
      if 0
        % simpler, but much slower 
        nearestedges = [nearestedges, ones(size(nearestedges,1),1)]; % default weight
        [NeighborsS1,NidS1]  = hist(nearestedges(:,1),1:1:max(nearestedges(:,1)));
        for ni=NidS1(NeighborsS1>1)
          NumNi = nearestedges(:,1)==ni; 
          nearestedges(NumNi,3) =  nearestedges(NumNi,3) ./ sum(NumNi);
        end
      else
        % faster 
        %nearestedges = [nearestedges, ones(size(nearestedges,1),1)]; % default weight
        dist = sum( (S2.vertices(nearestedges(:,1),:) - S1.vertices(nearestedges(:,2),:)).^2 , 2) .^ 0.5; 
        nearestedges = [nearestedges, dist]; % default weight
        list = [1; find(nearestedges(1:end-1,1)~=nearestedges(2:end,1))+1; size(nearestedges,1)]; 
        for ni=1:numel(list)-1
          %nearestedges(list(ni):list(ni+1)-1,3) = nearestedges(list(ni):list(ni+1)-1,3) ./ (list(ni+1) - list(ni)); 
          nearestedges(list(ni):list(ni+1)-1,3) = nearestedges(list(ni):list(ni+1)-1,3) ./ sum(nearestedges(list(ni):list(ni+1)-1,3)); 
        end
      end
      if opt.verb, stime = cat_io_cmd('  Mapping','g5','',1,stime); end

      %%
      if 0
        % correct & simple, but very slow
        varargout{1} = zeros(1,max(nearestedges(:,2)));
        for ni=1:size(nearestedges,1)
          varargout{1}(nearestedges(ni,2)) = varargout{1}(nearestedges(ni,2)) + ...
            cdata(nearestedges(ni,1))' .* nearestedges(ni,3)';
        end
      else
        varargout{1} = zeros(1,max(nearestedges(:,2)));
        if 0
          list = [1; find(nearestedges(1:end-1,2)~=nearestedges(2:end,2))+1; size(nearestedges,1)+1]; 
          for ni=1:numel(list)-1
            varargout{1}(nearestedges(list(ni),2)) = varargout{1}(nearestedges(list(ni),2)) + ...
              sum(cdata(nearestedges(list(ni):list(ni+1)-1,1)) .*  nearestedges(list(ni):list(ni+1)-1,3));
          end
        else
          nearestedges2 = sortrows([nearestedges(:,2) nearestedges(:,1) nearestedges(:,3)]);  
          list = [1; find(nearestedges2(1:end-1,1)~=nearestedges2(2:end,1))+1; size(nearestedges2,1)+1]; 
          for ni=1:numel(list)-1
            varargout{1}(nearestedges2(list(ni),1)) = varargout{1}(nearestedges2(list(ni),1)) + ...
              sum(cdata(nearestedges2(list(ni):list(ni+1)-1,2)) .*  nearestedges2(list(ni):list(ni+1)-1,3));
          end
        end
      end
      if numel(varargout{1})<20, disp(varargout{1}); end
      if opt.verb, cat_io_cmd(' ','g5','',1,stime); end
  end
  
  % default smoothing???
  if opt.smooth
    varargout{1}  = spm_mesh_smooth(struct('vertices',S3.vertices,'faces',S3.faces),double(varargout{1}'),opt.smooth);
  end
  
  if isfield(opt,'fname')
    save(gifti(struct('vertices',S1.vertices,'faces',S1.faces,'cdata',varargout{1})),opt.fname); 
  end
  
  if opt.verb, cat_io_cmd('','','',1,stime1); end
end

function [E,uE] = cat_surf_edges(T)
  if isstruct(T) && isfield(T,'faces')
    T = T.faces;
  end

  T = sort(T,2); E = []; 
  for i=1:size(T,2)-1
    E = [E; T(:,[i i+1])]; %#ok<AGROW>
  end
  [E,uE] = unique(E,'rows');
end

function D = cat_surf_dist(S)
% Estimate the distance between the vertices of the faces of S.
% D = [c,a,b] = [d(AB),d(BC),d(CA)]

  D = [sum( (S.vertices(S.faces(:,1),:) - S.vertices(S.faces(:,2),:)).^2 , 2) .^ 0.5, ...
       sum( (S.vertices(S.faces(:,2),:) - S.vertices(S.faces(:,3),:)).^2 , 2) .^ 0.5, ...
       sum( (S.vertices(S.faces(:,3),:) - S.vertices(S.faces(:,1),:)).^2 , 2) .^ 0.5]; 
     
end

function [AV,AF] = cat_surf_area(S)
% Calculate surface area of the faces AF (Horonsche Form) and map it to the
% vertices AV.

  % facearea (Horonsche Form)
  D = cat_surf_dist(S);
  facesp = sum(D,2) / 2;  % s = (a + b + c) / 2;
  AF = (facesp .* (facesp - D(:,1)) .* (facesp - D(:,2)) .* (facesp - D(:,3))).^0.5; % area=sqrt(s*(s-a)*(s-b)*(s-c));
  
  % numerical (to small point differences) and mapping problems (crossing of streamlines)
  % -> correction because this is theoretical not possible (Laplace field theory)
  AF(AF==0) = eps; % to small values
  AF = abs(AF);    % streamline-crossing
    
  AV = cat_surf_F2V(S,AF);
end

function data = cat_surf_F2V(S,odata)
%% mapping of facedata to vertices

  data   = zeros(size(S.vertices,1),1);
  [v,f]  = sort(S.faces(:)); 
  [f,fj] = ind2sub(size(S.faces),f);  
  far    = odata(f);
  for i=1:numel(v)
    data(v(i)) = data(v(i)) + far(i)/3; 
  end

  %  data = data ./ vcount; %size(S.vertices,2); % Schwerpunkt... besser Voronoi, aber wie bei ner Oberflaeche im Raum???
  
end

function A = cat_surf_smoothtexture(S,A,smooth,Amax)
%  create smooth area texture files
%  ---------------------------------------------------------------------
  debug = 0;
  
  if ~exist('smooth','var'), smooth=1; end

  % temporary file names
  Pname  = tempname; 
  Pmesh  = [Pname 'mesh'];
  Parea  = [Pname 'area'];

  if exist('Amax','var');  A = min(A,Amax); end 
  
  % write surface and textures
  cat_io_FreeSurfer('write_surf',Pmesh,S);
  cat_io_FreeSurfer('write_surf_data',Parea,A);
  
  % smooth textures
  cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s"',Pmesh,Parea,smooth,Parea);
  [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,debug);
  
  % load smoothed textures
  A  = cat_io_FreeSurfer('read_surf_data',Parea);
  
  % delete temporary file
  delete(Parea);
end

function [SH,V] = cat_surf_hull(S)
%% hull creation

  % render surface points
  Vi = cat_surf_fun('surf2vol',S);
  
  % fill mesh
  V  = cat_vol_morph(Vi,'ldc',mean(size(Vi))/6); clear Vi; % closing 
  V  = cat_vol_smooth3X(V,2);    % smoothing
  SH = isosurface(V,0.4);        % create hull 
  V  = V>0.4;
  
  % final mesh operations
  SH.vertices = [SH.vertices(:,2) SH.vertices(:,1) SH.vertices(:,3)]; % matlab flip
  SH.vertices = SH.vertices + repmat(min(S.vertices),size(SH.vertices,1),1) - 5;
end

function PTN = cat_surf_thickness(action,PS,PT)
  if ~exist('T','var')
    % create inner and outer surfaces
    PIS  = cat_surf_fun('inner',PS);  % estimate inner surface
    POS  = cat_surf_fun('outer',PS);  % estimate outer surface
  else
    % create inner and outer surfaces
    PIS  = cat_surf_fun('inner',PS,PT);  % estimate inner surface
    POS  = cat_surf_fun('outer',PS,PT);  % estimate outer surface
  end
  
  % load surfaces
  IS = gifti(PIS); 
  OS = gifti(POS);
  
  % edgemap
  % create mapping between 
  Pedgemap = cat_io_strrep(PS,{'.central.';'.gii'},{'.edgemapnative.';'.mat'});
  if 0%exist(Pedgemap,'file') % ... you have to test if central is older than the edgemap to use this 
    load(Pedgemap,'edgemap'); 
  else
    %%
    stime2  = clock;
    fprintf('  Estimate mapping for native surface');
    edgemap = cat_surf_surf2surf(IS,OS,0); 
    %edgemap.dist = sum ( (IS.vertices(edgemap.edges(:,1),:) - OS.vertices(edgemap.edges(:,2),:)).^2 , 2).^0.5;  
    %save(Pedgemap,'edgemap'); 
    fprintf(' takes %ds\n',round(etime(clock,stime2))); 
  end

  % create thickness metric mapping matrix
  switch lower(action)
    case {'tfs','tmin'}
      Tnear = inf(edgemap.nvertices(1),2,'single'); 
      for i=1:size(edgemap.edges,1)
        Tnear(edgemap.edges(i,1),1) = min( [ Tnear(edgemap.edges(i,1),1) edgemap.dist(i)  ] ) ;
        Tnear(edgemap.edges(i,2),2) = min( [ Tnear(edgemap.edges(i,2),2) edgemap.dist(i)  ] ) ;
      end
    case 'tmax'
      Tfar  = zeros(edgemap.nvertices(1),2,'single'); 
      for i=1:size(edgemap.edges,1)
        Tfar(edgemap.edges(i,1),1)  = max( [ Tfar(edgemap.edges(i,1),1) edgemap.dist(i)  ] ) ;
        Tfar(edgemap.edges(i,2),2)  = max( [ Tfar(edgemap.edges(i,2),2) edgemap.dist(i)  ] ) ;
      end
  end
    
  switch lower(action)
    case 'tfs'
      TN  = mean(Tnear,2);
      PTN = cat_io_strrep(PS,{'.central.';'.gii'},{'.thicknessfs.';''});
    case 'tmin'
      TN  = min(Tnear,[],2);
      PTN = cat_io_strrep(PS,{'.central.';'.gii'},{'.thicknessmin.';''});
    case 'tmax'
      TN  = max(Tfar,[],2);
      PTN = cat_io_strrep(PS,{'.central.';'.gii'},{'.thicknessmax.';''});
  end
  
  % save smoothed textures
  cat_io_FreeSurfer('write_surf_data',PTN,TN);
end

function [SH,V] = cat_surf_core(S,opt)
%% core creation
%  This is much more complicated that the hull definition. So I will need 
%  different types of core definitions. However, I first have to find one
%  (or multiple) anatomical definitions.
%
%  For estimation I can use different techniques: 
%   * morphological operations 
%     > very inaccurate and error-prone 
%     > use of distance & smoothing functions 
%   * smoothing with/without boundaries
%   * anatomical information from volume or better surface atlas maps 
%   * use of other measures such as 
%     - thickness (no)
%     - sulcal depth or outward folding GI (maybe)
%     - curvature (not really)
%
%   * use of percentual scalings
%   * use of multiple threshold levels and averaging to avoid using only 
%     one threshold (=multiband) 
%   * definition as fractal dimension measure?
%

  def.type = 1; 
  def.th   = 0.15;
  opt = cat_io_checkinopt(opt,def); 
  
  
  % render surface points
  Vi = cat_surf_fun('surf2vol',S);
   
  %% break gyri
  if opt.type == 1
    %%
    Vd  = cat_vbdist(single(Vi<0.5));
    
    %%
    Vdn = Vd ./ max(Vd(:));
    Vdn = cat_vol_laplace3R(Vdn,Vdn>0 & Vdn<0.8 ,0.001);
    Vdn = min(opt.th + cat_vol_morph(Vdn>opt.th,'l'),Vdn); 
    V   = Vdn > opt.th; 
    
    %%
    SH  = isosurface(Vdn,opt.th);        % create hull 
  elseif opt.type == 2
    SiGI = S; SiGI.cdata = opt.iGI;
    V    = cat_surf_fun('surf2vol',SiGI,struct('pve',3));
    V    = V ./ max(V(:));
    SH   = isosurface(V,opt.th);        % create hull 
  else
    Vs = cat_vol_smooth3X(Vi,8);    % smoothing
    V  = ~cat_vol_morph(Vs<0.5,'ldc',min(size(Vi))/(6*1.5));  % opening 
    V  = cat_vol_smooth3X(V,6);    % smoothing
    SH = isosurface(V .* smooth3(Vi),0.6);        % create hull 
    V  = min(V>0.6,Vi==1);
    V  = cat_vol_smooth3X(V,4);    % smoothing
    V  = min(V,Vi);
    V  = cat_vol_laplace3R(V,Vi>0 & V<0.9,0.4);
  end   %clear Vi;

  % final mesh operations
  SH.vertices = [SH.vertices(:,2) SH.vertices(:,1) SH.vertices(:,3)]; % matlab flip
  SH.vertices = SH.vertices + repmat(min(S.vertices),size(SH.vertices,1),1) - 5;
end
function cat_surf_saveICO(S,Tpbt,Pcs,subdir,writeTfs,C)
% Creates and save the white and pial surfaces based on the displacement by
% the half thickness along the surface normals and use the inner and outer
% surfaces to create the layer4 surface.
% Saves also the thickness file.

  [pp,ff,ee] = spm_fileparts(Pcs);
  if ~exist('subdir','var')
    subdir = '';
  else
    if ~exist(fullfile(pp,subdir),'dir')
      mkdir(fullfile(pp,subdir)); 
    end
  end
 
  if isfield(S,'vmat') && isfield(S,'mati') 
    S.vertices = (S.vmat * [S.vertices' ; ones(1,size(S.vertices,1))])';
    if S.mati(7)<0, S.faces = [S.faces(:,1) S.faces(:,3) S.faces(:,2)]; end
  end
  
  % normalized surface normals
  N = spm_mesh_normals(S);                             
 
  % inner and outer surface
  VI  = S.vertices + N .* repmat(Tpbt/2,1,3); 
  VO  = S.vertices - N .* repmat(Tpbt/2,1,3); 
 
  % filenames
  Pcentral = fullfile(pp,subdir,[ff ee]);   
  Pwhite   = fullfile(pp,subdir,strrep([ff ee],'central','white'));   
  Ppial    = fullfile(pp,subdir,strrep([ff ee],'central','pial'));   
  Pthick   = fullfile(pp,subdir,strrep(ff,'central','thickness'));   
  Ppbt     = fullfile(pp,subdir,strrep(ff,'central','pbt'));   
  Pcol     = fullfile(pp,subdir,strrep(ff,'central','collision'));   
  Player4  = fullfile(pp,subdir,strrep([ff ee],'central','layer4'));   
  
  % save surfaces
  save(gifti(struct('faces',S.faces,'vertices',S.vertices)),Pcentral,'Base64Binary');
  save(gifti(struct('faces',S.faces,'vertices',VI)),Pwhite,'Base64Binary');
  save(gifti(struct('faces',S.faces,'vertices',VO)),Ppial,'Base64Binary');

  % save thickness
  cat_io_FreeSurfer('write_surf_data',Ppbt,Tpbt);
  if exist('writeTfs','var') && writeTfs
    cmd = sprintf('CAT_SurfDistance -mean "%s" "%s" "%s"',Pwhite,Ppial,Pthick);
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
    fprintf('Display thickness: %s\n',spm_file(Pthick ,'link','cat_surf_display(''%s'')'));
  end
  if exist('C','var')
    cat_io_FreeSurfer('write_surf_data',Pcol,C);
  end
  
  % final correction of central surface in highly folded areas with high mean curvature
  cmd = sprintf('CAT_Central2Pial -equivolume -weight 1 "%s" "%s" "%s" 0', ...
                   Pcentral,Ppbt,Player4);
  [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);
  
  % display something to click
  fprintf('\n    Display surface: %s\n',spm_file(Ppbt  ,'link','cat_surf_display(''%s'')'));
end
function [SN,TN] = cat_surf_collision_correction(S,T,Y,debug,Pcs)
% Delaunay based collision dectection:
% 1) Correction for local curvature
% 2) Creation of a Delaunay graph 
% 3) Differentiation of intra (edges between surface points) and inter
%    surface edges (e.g., edges between to opposite gyri or sulci) and 
%    removal of intra-surface edges. 
% 4) Use of the inter-surfaces edges to detect collision by normal 
%    transformations of the half thickness to obtain the inner and outer 
%    surfaces.
% 5) Further correction of possible flips by normal transformation
%
% This is a prototype that allows correction of the worst things but not 
% all collisions. 
%
%   [SN,TN] = cat_surf_collision_correction(S,T,Y,debug,Pcs)
%
%   SN    .. new surface
%   TN    .. new thickness
%   S     .. original surface
%   T     .. original thickness
%   Y     .. segmentation map or intensity normalized images 
%            for intra/inter surface edge definition
%   debug .. option to write the un- and corrected cortical surfaces in a
%            subdirectory
%   Pcs   .. central surface file name
%
% Robert Dahnke 201909


  % correction for flipped faces to have always the same normal direction
  if S.mati(7)<0, S.faces = [S.faces(:,1) S.faces(:,3) S.faces(:,2)]; S.mati(7) = - S.mati(7); flipped = 1; else flipped = 0; end

  if ~exist('debug','var'), debug = 1; end
  write = debug & exist('Pcs','var');
  
  % larger surface need more smoothing to avoid triangulation problems 
  sf = round( (size(S.faces,1) / 50000).^0.5 );  % empirical value 
  if max(Y(:))<1.5, Y = Y.*2+1; end               
  
  if debug, fprintf('\n'); end
  if write, cat_surf_saveICO(S,T,Pcs,sprintf('pre_collcorr_%0.0fk',round( size(S.faces,1)/1000 / 10) * 10 ),0); end
  stime = cat_io_cmd(sprintf('    Delaunay triangulation of %d vertices:',size(S.vertices,1)),'g5','',debug); 
  
  % smothing fucntions for data and surfaces
  M   = spm_mesh_smooth(S);   % for spm_smoothing matrix
  rms = @(x) sum(x.^2).^0.5;  % for error handling of mad vertices
  smoothsurf = @(V,s) [ ...   % simple surface smoothing 
    spm_mesh_smooth(M,double(V(:,1)),s) , ...
    spm_mesh_smooth(M,double(V(:,2)),s) , ...
    spm_mesh_smooth(M,double(V(:,3)),s) ];
  
  
  
  %% Curvature based correction
  %  The local mean curvature describes simplified the inverse radius of a
  %  fitting sphere that describes one limit of a normal transformation. 
  %  If the thickness is larger that the radius we have to correct it and 
  %  update the central surface point. 
  %  The empirical values are not further validated and only tested for the
  %  "single_subj" and "HR075" datasets (RD20190912).
  %
  V  = S.vertices; 
  SS = S; SS.vertices = smoothsurf(V,2);        % empirical value 
  C  = spm_mesh_curvature( SS ); clear SS; 
  C  = spm_mesh_smooth(M,C,4);                  % empirical value
  N  = spm_mesh_normals(S);     
  
  % cortical thickness is limit by the local curvature because C = 1/R 
  TI = T/2;  TI(C>0) = min( TI(C>0) , 0.5./abs( C(C>0) ) ); % inner thickness between central and white surface
  TO = T/2;  TO(C<0) = min( TO(C<0) , 0.5./abs( C(C<0) ) ); % outer thickness between central and outer surface
  clear C;
  
  % smooth result
  TI = single(spm_mesh_smooth(M,double(TI),5)); % empirical value
  TO = single(spm_mesh_smooth(M,double(TO),5)); % empirical value
  
  
  % estimate the new cortical thickness and correct the inner and outer 
  % surface to update the central point
  VO = V - N .* repmat(TO,1,3);
  VI = V + N .* repmat(TI,1,3);
  
  S.vertices = mean(cat(3,VI,VO),3); 
  if 0
  % debugging function that directly returns this first correction
    TN = TO + TI; 
    SN = S; 
    if write, cat_surf_saveICO(SN,TN,Pcs,sprintf('post_collcorr_%0.0fk',round( size(S.faces,1)/1000 / 10) * 10 ),0); else, fprintf('\n'); end
    return;
  else
    T  = TO + TI; 
  end
  clear VO VI TI TO N; 
  
  
  
  
  %% edge creation based on Delaunay graph 
  %  smooth central surface to create robust Delaunay
  
  % surface smoothing (as loop to correct for outlier due to incorrect surfaces)
  VS  = double(S.vertices); 
  for i = 1:100*sf
    VSi = smoothsurf(VS,2); 
    VM  = rms(VSi-VS)<2; 
    VS(VM,:) = VSi(VM,:); 
  end
  VS = smoothsurf(VS,1);
  if ~debug, clear VSi VM; end
  
  % delaunay graph
  D  = single(delaunayn( VS )); 
  if ~debug, clear VS; end            
  
  % decompose delaunay graph into its edges
  E  = uint32(cat_surf_edges(D));   
  nE = size(E,1); 
  if ~debug, clear D; end
  if debug, cat_io_cprintf('g5',sprintf('%5.0fs\n',etime(clock,stime))); end
  
  
  
  
  %% Remove intra-surface edges 
  %  if you remove to much then the correction will not work
  %  if you do not remove enough then it will add sulci in regions without sulci
  V  = S.vertices;
  
  
  % remove edge that we know from the surface - super save
  stime = clock; 
  EF = uint32(cat_surf_edges(S.faces));           
  E  = setdiff(E,EF,'rows'); clear EF; 
  % remove edges between neigbors of each point - relative save
  nlevel   = 1 + ( size(S.faces,1) > 80000 ) + ( size(S.faces,1) > 50000 ); % use higher levels only for large surfaces
  [NE,MED] = spm_mesh_neighbours(M); 
  if nlevel>1, [NE2,MED2] = spm_mesh_neighbours(M); NE = [NE NE2]; MED = [MED MED2]; end
  if nlevel>2, [NE2,MED2] = spm_mesh_neighbours(M); NE = [NE NE2]; MED = [MED MED2]; end
  clear NE2 MED2
  for i=2:size(NE,2)
    E = setdiff(E,[NE(:,1) NE(:,i)],'rows'); 
    E = setdiff(E,[NE(:,i) NE(:,1)],'rows'); 
  end
  clear NE
  if debug, cat_io_cprintf('g5',sprintf('    remove edges by surface:  %10d > %10d (%0.2f%%%%) %9.0fs',nE,size(E,1),size(E,1)./nE,etime(clock,stime))); end
  
  
  % remove edge by distance - this is not clear but it helps
  stime = clock; 
  LE  = sum( (V(E(:,1),:) - V(E(:,2),:)).^2 , 2) .^ 0.5; % length of edge
  DE  = min( LE > max(0.5,min(1,max(MED(:)))) , min( LE - T(E(:,1))*0.33 , LE - T(E(:,2))*0.33 ));
  NEd = abs(DE); clear LE DE MED;
  
  % remove by angle .. sum(NEa)./numel(NEa), figure, hist( S1alpha, -180:1:180)
  %   N(S)alpha  .. angle between the (smoothed) normals of each edge
  %                 (~0? = surface edge; ~180? between surface edge)       
  %   S[12]alpha .. angle between the edge and the first normal    
  %                 (~0?/~180? = surface edge; ~90? = between surface edge)
  % 
  N  = spm_mesh_normals(S);                 
  NS = N; for i=1:80*sf, NSS = smoothsurf(NS,1); NM = rms(NS - NSS)<0.5; NS(NM,:) = NSS(NM,:); end 
  Nalpha  = [angle(NS(E(:,1),:), NS(E(:,2),:)), ...
             angle(NS(E(:,2),:), NS(E(:,1),:))]; 
  SNalpha = [angle(N(E(:,1),:),  V(E(:,1),:) - V(E(:,2),:)), ...
             angle(N(E(:,2),:),  V(E(:,2),:) - V(E(:,1),:))]; 
  NEna    = mean(Nalpha/180,2);                                    % figure, hist( NEna , 0:0.01:1);
  NEsa    = (abs(90  - SNalpha)/90  + abs(90  - SNalpha)/90)/2;    % figure, hist( NEsa , 0:0.01:1);
  
  % remove by intensity given by the centroids of the edges
  VC  = cat_surf_centroid(V,E); 
  IC  = isocolors(Y,VC); clear VC;
  % outer surface intensity
  VO  = V - N .* repmat(T/2,1,3); 
  VOC = cat_surf_centroid(VO,E); 
  IO  = isocolors(Y,VOC); clear VOC VO; 
  % inner surface intensity
  VI  = V + N .* repmat(T/2,1,3) + 0.1; % GM/WM  
  VIC = cat_surf_centroid(VI,E); 
  II  = isocolors(Y,VIC); clear VIC VI; 
  VI  = V + N .* repmat(T/2,1,3) + 0.5; % save WM  
  VIC = cat_surf_centroid(VI,E); 
  II  = max(II,isocolors(Y,VIC)); clear VIC VI; % use max to get WM value 
  VI  = V + N .* repmat(T/2,1,3) + 1.0; % supersave WM  
  VIC = cat_surf_centroid(VI,E); 
  II  = max(II,isocolors(Y,VIC)); clear VIC VI; % use max to get WM value 
  % combine all intensities 
  NEi = 1 - min(1,max(abs(diff([II IC IO],1,2)),[],2)); 
  %ET  = mean([II IC IO],2)>2.25; % edge classification 
  if ~debug, clear II IC IO; end

  % combine all measures by product to remove many things
  % I also though about an adaptive threshold but it is not so easy ...
  NE = prod( [NEd NEi NEna NEsa] ,2); 
  NE = NE < .05; %max(eps,mean(NE) - 1*std(NE)); 
  E (NE,:) = []; if exist('ET','var'), ET(NE) = []; end
  if debug, cat_io_cprintf('g5',sprintf('\n    remove edges by surface:  %10d > %10d (%0.2f%%%%) %9.0fs',...
      nE,size(E,1),size(E,1)./nE,etime(clock,stime))); end

  
 %fprintf('\nsf = %0.2f',sf);
  SNalpha = [angle(N(E(:,1),:),  V(E(:,1),:) - V(E(:,2),:)), ...
             angle(N(E(:,2),:),  V(E(:,2),:) - V(E(:,1),:))]; 
           
           
  %% Iterative correction routine
  TN = T; SN = S; 
  maxiter = 10; TCsumo = inf; 
  for j=1:maxiter*10
    V   = single(SN.vertices);

    % update surface normales  
    NS  = spm_mesh_normals(SN); 
    
    % inner and outer surface
    VO  = V - N .* repmat(TN/2,1,3);
    VI  = V + N .* repmat(TN/2,1,3);

    
    % first correction step that works but also could be improved
    if 1 % general debugging setting to ignore this block
    
      % edgelength of the central, inner, and outer surface
      LE  = sum( (V(E(:,1),:)  - V(E(:,2),:)).^2  , 2) .^ 0.5;
      LEO = sum( (VO(E(:,1),:) - VO(E(:,2),:)).^2 , 2) .^ 0.5;
      LEI = sum( (VI(E(:,1),:) - VI(E(:,2),:)).^2 , 2) .^ 0.5;

      
      % correct surface 
      % - there is the problem that not all points are connected, 
      %   so it is necessary to smooth to include unconnected neighbors  
      % - in general, only inter-surface-edges are expected here, those
      %   distance describes the maximal local thickness
      % - we also add some sulcus-width to avoid collisions but the effect
      %   will be small due to the smoothing
      % - LEC, LEOC, and LEIC represent the distance error by collisions
      %   of each edge
      sulcuswidth = 0.2;
      LEC  = max(0, (TN(E(:,1))/2 + TN(E(:,2))/2) - max(0,LE - sulcuswidth )); % * (maxiter/10 - j) / (maxiter/10)) );
      LEOC = LEO*0; %max(0, sulcuswidth - LEO);
      LEIC = LEI*0; %max(0, sulcuswidth - LEI);

      % map the collision edge error to the surface by using the maximum
      simple = 0; 
      if simple 
        % Simple correction without side differentiation:
        % This corrects towards the CS without further differentiation of 
        % the reason of the correction. I.e. and outer surface collision 
        % correction will also move the inner surface towards the central 
        % surfaces. 
        %
        % This was my first step and can be removed later (RD201909). 

        TC  = TN*0; 
        for ni=1:size(E,1)
          TC(E(ni,1)) = max([TC(E(ni,1)), LEC(ni), LEOC(ni), LEIC(ni)]); 
          TC(E(ni,2)) = max([TC(E(ni,2)), LEC(ni), LEOC(ni), LEIC(ni)]); 
        end 
        TCsum = max(TC(TC>0)); 

        % smooth the correction map and reestimate the inner and outer surface
        % to estimate a corrected Tlink thickness 
        TC  = max(0,min(TC,single(spm_mesh_smooth(M,double(TC),max(0,2*sf)))*0.75 + single(spm_mesh_smooth(M,double(TC),max(0,0.5*sf)))*0.75 )); 
        TNC = max(0,TN - TC);                       % estimate corrected thickness  
        VOC = V - NS .* repmat(TNC/2,1,3);          % outer surface
        VIC = V + NS .* repmat(TNC/2,1,3);          % inner surface

        % thickness smoothness is a problem by the number of iterations so I have to limit it! 
        % however we apply this correction only in affected areas
        TNC = sum( (VIC - VOC).^2 , 2) .^ 0.5;      % Tlink 
        TNS = single(spm_mesh_smooth(M,double(TNC),max(1,min(10,sf)) )); % * ((maxiter*10) - j)/(maxiter*10) ));clear TNC; 
        TC  = single(spm_mesh_smooth(M,double(TC ),max(1,min(10,sf)) )); % * ((maxiter*10) - j)/(maxiter*10) ));
        TN(TC>0) = TNS(TC>0);                      % only in corrected areas; 
        VOC = V - NS .* repmat(TN/2,1,3);          % outer surface
        VIC = V + NS .* repmat(TN/2,1,3);          % inner surface
      
        % filter surfaces
        VOS = smoothsurf(VOC,1); VOC(TC>0) = VOC(TC>0)*0.5 + 0.5*VOS(TC>0); clear VOS;
        VIS = smoothsurf(VIC,1); VIC(TC>0) = VIC(TC>0)*0.5 + 0.5*VIS(TC>0); clear VIS;
        
      else
        % Complex correction with side differentiation: 
        % Here, the angle between of the inter-surface edges is used to
        % differentiate between the edges within a suclus (outer) or within
        % a gyrus (inner). 
        TOC = TN*0; TIC = TN*0;
        for ni=1:size(E,1)
          TOC(E(ni,1)) = max([TOC(E(ni,1)), LEC(ni) .* (SNalpha(ni,1)<=90), LEOC(ni)]); 
          TOC(E(ni,2)) = max([TOC(E(ni,2)), LEC(ni) .* (SNalpha(ni,2)<=90), LEOC(ni)]); 
          TIC(E(ni,1)) = max([TIC(E(ni,1)), LEC(ni) .* (SNalpha(ni,1)>90), LEIC(ni)]); 
          TIC(E(ni,2)) = max([TIC(E(ni,2)), LEC(ni) .* (SNalpha(ni,2)>90), LEIC(ni)]); 
        end 

        % estimate first corrected inner and outer thickness 
        TOC = max(0,min(TOC,single(spm_mesh_smooth(M,double(TOC),max(0,4.0*sf))*1.50 ) + ...
                            single(spm_mesh_smooth(M,double(TOC),max(0,2.0*sf))*0.75 ) + ...
                            single(spm_mesh_smooth(M,double(TOC),max(0,0.5*sf))*0.25 ) ));
        TIC = max(0,min(TIC,single(spm_mesh_smooth(M,double(TIC),max(0,4.0*sf))*1.50 ) + ...
                            single(spm_mesh_smooth(M,double(TIC),max(0,2.0*sf))*0.75 ) + ...
                            single(spm_mesh_smooth(M,double(TIC),max(0,0.5*sf))*0.25 ) ));
        TC  = TOC + TIC; TCsum = mean(TC(TC>0)); 
        
        % estiamte new inner and outer surfaces
        VOC = V - NS .* repmat( max(0,TN/2 - TOC/2 - TIC/2),1,3);    % outer surface 
        VIC = V + NS .* repmat( max(0,TN/2 - TIC/2 - TOC/2),1,3);    % inner surface
        clear TIC TOC;
  
        % filter surfaces - that is not so easy and will "deform" some
        % cortical structures .. so you need a mask .. maybe by curvature
        % information?
     %   VOS = smoothsurf(VOC,1); VOC(TC>0) = VOS(TC>0); clear VOS;
     %   VIS = smoothsurf(VIC,1); VIC(TC>0) = VIS(TC>0); clear VIS;

        % estimate new smooth thickness function - here smoothing is OK!
        TNC = sum( (VIC - VOC).^2 , 2) .^ 0.5;      % Tlink 
        TNS = single(spm_mesh_smooth(M,double(TNC),max(1,min(10,sf)) ));
        TC  = single(spm_mesh_smooth(M,double(TC ),max(1,min(10,sf)) )); 
% nutze diese information erst um ausreißer zu erfassen ... in einer 4 mm dicken regione wird kein 1 mm sulcus sein!        
        
        TN(TC>0) = TNS(TC>0);                      % only in corrected areas; 
        clear TNS TNC
        
        % update central surface
        VN  = mean(cat(3,VIC,VOC),3); 
        NS  = spm_mesh_normals(SN); 
       % VNS = smoothsurf(VN,1); VN(TC>0) = VNS(TC>0); clear VNS;
         
        % update inner and outer surfaces
        VOC = VN - NS .* repmat(TN/2,1,3);          % outer surface
        VIC = VN + NS .* repmat(TN/2,1,3);          % inner surface

% The idea was to improve the description here but the smoothing did now 
% work and was overall too slow. 
%{        
        % filter surfaces
      %  VOS = smoothsurf(VOC,1); VOC(TC>0) = VOS(TC>0); clear VOS;
      %  VIS = smoothsurf(VIC,1); VIC(TC>0) = VIS(TC>0); clear VIS;
        
        % thickness smoothness is a problem by the number of iterations so I have to limit it! 
        % however we apply this correction only in affected areas
        TNC = sum( (VIC - VOC).^2 , 2) .^ 0.5; 
        TNS = single(spm_mesh_smooth(M,double(TNC),max(1,min(10,sf)) )); % * ((maxiter*10) - j)/(maxiter*10) ));clear TNC; 
        TOC = single(spm_mesh_smooth(M,double(TOC),max(1,min(10,sf)) )); % * ((maxiter*10) - j)/(maxiter*10) ));
        TIC = single(spm_mesh_smooth(M,double(TIC),max(1,min(10,sf)) )); % * ((maxiter*10) - j)/(maxiter*10) ));
        TN(TC>0) = TNS(TC>0);

        VN  = mean(cat(3,VIC,VOC),3); 
        VOC = VN - NS .* repmat(TN/2 - TOC,1,3);          % outer surface
        VIC = VN + NS .* repmat(TN/2 - TIC,1,3);          % inner surface
      
        % filter surfaces
        %VOS = smoothsurf(VOC,1); VOC(TC>0) = VOC(TC>0)*0.0 + 1.0*VOS(TC>0); clear VOS;
        %VIS = smoothsurf(VIC,1); VIC(TC>0) = VIC(TC>0)*0.0 + 1.0*VIS(TC>0); clear VIS;
%}

      end
    else
      VOC = VO;
      VIC = VI;
    end
    % upate thickness and surface
 %   TN  = sum( (VIC - VOC).^2 , 2) .^ 0.5;
 %   SN.vertices = mean(cat(3,VIC,VOC),3); 
    
 
 
 
    %% detect reminding small overlaps in the following way
    %  - the CS is expected to be correct and overlaps happens afters only for the IS and OS. 
    %  - this means that the angle between the edge and the normal will flip alpha(CS) 
    %  - I am not sure if this is really working (RD20190912)
    if 1
      NCalpha = angle(N(E(:,1),:),  V(E(:,1),:)   - V(E(:,2),:)); 
      NOalpha = angle(N(E(:,1),:),  VOC(E(:,1),:) - VOC(E(:,2),:)); 
      NIalpha = angle(N(E(:,1),:),  VIC(E(:,1),:) - VIC(E(:,2),:)); 
      % flipped angles
      NOcorr  = (NCalpha - NOalpha) > 90; 
      NIcorr  = (NCalpha - NIalpha) > 90;
     % Ncorr   = TNC * 0; Ncorr(E(NOcorr | NIcorr,1)) = 1; 
      
      %% correction vector and correction for outer surface
      if 1
        sulcuswidth = 0.2; 
        NNcorr  = VOC(E(NOcorr,2),:) - VOC(E(NOcorr,1),:); 
        VOC(E(NOcorr,1),:) = VOC(E(NOcorr,1),:) + NNcorr*(0.5 + sulcuswidth); 
        VOC(E(NOcorr,2),:) = VOC(E(NOcorr,2),:) - NNcorr*(0.5 + sulcuswidth); 
        % correction vector and correction for inner surface
        NNcorr  = VIC(E(NIcorr,2),:) - VIC(E(NIcorr,1),:); 
        VIC(E(NIcorr,1),:) = VIC(E(NIcorr,1),:) + NNcorr*(0.5 + sulcuswidth); 
        VIC(E(NIcorr,2),:) = VIC(E(NIcorr,2),:) - NNcorr*(0.5 + sulcuswidth);
      else
        % neu - geht wohl nicht
        NNcorr  = VOC(E(NOcorr,2),:) - VOC(E(NOcorr,1),:); 
        NNcorrd = sum( NNcorr.^2 , 2) .^ 0.5;
        VOC(E(NOcorr,1),:) = V(E(NOcorr,1),:) - NS(E(NOcorr,1),:) .* repmat(max(0,TN(E(NOcorr,1),:)/2 - NNcorrd/2),1,3); 
        VOC(E(NOcorr,2),:) = V(E(NOcorr,2),:) - NS(E(NOcorr,2),:) .* repmat(max(0,TN(E(NOcorr,2),:)/2 - NNcorrd/2),1,3); 
        % correction vector and correction for inner surface
        NNcorr  = VIC(E(NIcorr,2),:) - VIC(E(NIcorr,1),:); 
        NNcorrd = sum( NNcorr.^2 , 2) .^ 0.5;
        VIC(E(NIcorr,1),:) = V(E(NIcorr,1),:) + NS(E(NIcorr,1),:) .* repmat(max(0,TN(E(NIcorr,1),:)/2 - NNcorrd/2),1,3); 
        VIC(E(NIcorr,2),:) = V(E(NIcorr,2),:) + NS(E(NIcorr,2),:) .* repmat(max(0,TN(E(NIcorr,2),:)/2 - NNcorrd/2),1,3); 
      end
      %NC  = TN - sum( (VIC - VOC).^2 , 2) .^ 0.5; TCsum = max(NC(NC>0)); 
    end
    
    
    
    
    %% upate thickness and surface
    TN  = sum( (VIC - VOC).^2 , 2) .^ 0.5;
    SN.vertices = mean(cat(3,VIC,VOC),3); 
    
    
    %% this is just for display and loop settings
    %  SX.vertices = VOC; SX.faces = S.faces; SX.facevertexcdata = TC; cat_surf_render2(SX);
    stopiterth = 0.0002; 
    if debug && ( j==1 || mod(j,10)==1 || TCsum<0.01 || abs(TCsumo - TCsum)<stopiterth ) 
      TNM = TN>(mean(TN(:)) - 2*std(TN(:))) & TN<(mean(TN(:)) + 2*std(TN(:)));
      cat_io_cprintf('g5',sprintf('\n    mean reminding overlap: %8.2f mm (Tlink: %4.2f%s%4.2f mm) %9.0fs',...
        TCsum,mean(TN(TNM)),char(177),std(TN(TNM)),etime(clock,stime) )); stime = clock;
    end
    if TCsum<0.01 || abs(TCsumo - TCsum)<stopiterth, break; end
    TCsumo = TCsum;
  end
  
  % export cortical surfaces
  if write, cat_surf_saveICO(SN,TN,Pcs,sprintf('post_collcorr_%0.0fk',round( size(S.faces,1)/1000 / 10) * 10 ),0); else fprintf('\n'); end
  
  % flip back
  if flipped, SN.faces = [SN.faces(:,1) SN.faces(:,3) SN.faces(:,2)]; SN.mati(7) = - SN.mati(7); end
  

end
function C = cat_surf_centroid(V,F,n)
% _________________________________________________________________________
% calculates the centroid of a region
% _________________________________________________________________________

  if ~exist('n','var'), n=1; end
  
  ndim = size(F,2);
  
  switch ndim
    case 2, ET = [1,2];
    case 3, ET = [1,2;1,3];
    case 4, ET = [1,2;1,3;1,4]; 
  end
    
  C = repmat( V(F(:,1),:) , 1, 1, n); 
  for e = 1:size(ET,1)
    ed = diff( cat( 3 , V(F(:,ET(e,1)),:) , V(F(:,ET(e,2)), :) ) , 1 , 3 ); 
    for ni = 1:numel(n)
      C(:,:,ni) = C(:,:,ni) + ni/(n+1) * ed;
    end
  end
end  
function [Yd,Yv] = cat_surf_vdist(S,V,M,opt)
% CAT surface rendering with PVE by distance approximation.
%
% [Yd,Yv] = cat_surf_render(S,V,opt)
%
%   Yd    .. distance map
%   Yv    .. surface to PVE map rendering
%   S     .. surface with verices and faces
%   V     .. given volume or SPM volume structure
%   opt   .. option structure
%    .res .. higher surface resolution (0-default,1-interp)
% 

% Improve speed by voxel-based pp of distance parts, if only Yv is relevant? 
%

  def.res  = 0;  
  def.fast = 0; 
  opt = cat_io_checkinopt(opt,def);

  % improve surface resolution?
  if opt.res
    % ...
    S = cat_surf_fun('interp',S);
  end
  
  %% setup volume and transform vertices 
  if ~exist('V','var')
  % if not given create any volume
    Y  = false( round(max(S.vertices,[],1) - min(S.vertices)) + 10 );     
    Sv = S.vertices - repmat( min(S.vertices,[],1) - 5 , size(S.vertices,1)  , 1 );
  elseif isstruct(V)
    % modify coordinates by orientation matrix
    Y  = false( V.dims );    
    Sv = [S.vertices' ones(1,size(S.vertices,1))] .* V.mat;  
  else
    % simply center the surface in the given volume
    os = round( (size(V) - (max(S.vertices,[],1) - min(S.vertices,[],1))) / 2 ); 
    Sv = S.vertices - repmat( min(S.vertices,[],1) + os , size(S.vertices,1)  , 1 );
  end

  
  %% estimate surface normals to have negative distances inside the surface
  Sn = spm_mesh_normals(S); 
  
  %% distance estimation 
  [VB,Svia] = unique( Sv , 'rows' );      % required for delaunay
  VN = Sn(Svia,:);                          % clear Sn
  VB = double(VB);                        % needed for delaunayn
  T  = delaunayn( VB );                   % delaunayn graph for faster dsearchn processing
  [VR(:,1),VR(:,2),VR(:,3)] = ind2sub(size(Y),1:numel(Y)); % x-y may be flipped!
  [VID,VDD] = dsearchn(VB,T,VR); clear T; % search nearest point with its distance
  VB = single(VB); 
  VM = VR - VB( VID ,:);                  % vector from surface point to voxel
  
  %% estimate if voxel is inside S 
  %  ... this is much to slow ... the convertation to cell should not be possible
  %  ... and this is unused
  %{
  VRc   = mat2cell(VR,ones(size(VR,1),1));
  VNc   = mat2cell(VN(VID,:),ones(size(VR,1),1));
  VRNa  = cellfun( @(u,v) acosd( (u*v') / (sum(u'.^2)^.5 * sum(v'.^2)^.5)) ,VRc,VNc,'UniformOutput',false);
  VSD   = cell2mat(VRNa);
  %}
  
  %% estimate surface normals to use a weighted
  VMVR = mat2cell( cat( VM,VN,ones(size(VM)) , 3 ) , ones(1,size(VM,1)) ); 
  VMVR = cellfun(@shiftdim,VMVR); 
  VSD  = cellfun(@det,VMVR);
  
  Yd = reshape(VDD,size(Y)) .* sign(90-reshape(VSD,size(Y))); 
  Yv = min(1,max(0,Yd + 0.5)); 
end

function [V,vmat,vmati] = cat_surf_surf2vol(S,opt)
%% render inner surface area 
%  Render the volume V with V==1 within the surface. 
%  Use type=1 to render also the surface area with 0.5.
%  The transformation imat to create 
%  SH.vertices = [SH.vertices(:,2) SH.vertices(:,1) SH.vertices(:,3)];     % matlab flip
%  SH.vertices = SH.vertices + imat;

  if ~exist('opt','var'), opt = struct(); end
  def.debug  = 0;
  def.pve    = 1;     % 0 - no PVE; 1 - PVE;
                      % 2 - fill with surface texture values without interpolation and masking (==4)    
                      % 3 - fill with surface texture values with    interpolation and masking (==5)
  def.refine = 0.6;
  def.bdist  = 5; 
  def.res    = 1; % not yet ...
  
  opt = cat_io_checkinopt(opt,def);
  
  % save a temporary version of S and refine it
  Praw = [tempname '.gii'];
  save(gifti(struct('vertices',S.vertices,'faces',S.faces)),Praw); 
  
  cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f',Praw,Praw,opt.refine); 
  [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.debug);

  So = S; S = gifti(Praw);
  delete(Praw);
  
  %% render surface points
  if opt.pve > 1
    % get surface data or give error
    if isfield(So,'cdata')
      cdata = So.cdata; 
    elseif  isfield(So,'facevertexcdata')
      cdata = So.facevertexcdata; 
    else
      error('cat_surf_fun:cat_surf_surf2vol:No datafield for filling');
    end
    
    % render data 
    V    = false( round(max(S.vertices,[],1) - min(S.vertices)) + opt.bdist*2 );     
    vmat = -[min(S.vertices(:,1)) min(S.vertices(:,2)) min(S.vertices(:,3))] + opt.bdist; 
    I    = sub2ind(size(V),...
          max(1,min(size(V,1),round(S.vertices(:,1) + vmat(1)))),...
          max(1,min(size(V,2),round(S.vertices(:,2) + vmat(2)))),...
          max(1,min(size(V,3),round(S.vertices(:,3) + vmat(3)))));
    V(I) = 1; 
    V    = cat_vol_morph(V,'lc',1);  % closeing 
    
    % data filling
    Vv   = zeros( round(max(S.vertices,[],1) - min(S.vertices)) + opt.bdist*2 ,'single');  % same size so S and not So   
    I    = sub2ind(size(V),...
          max(1,min(size(V,1),round(So.vertices(:,1) + vmat(1)))),...
          max(1,min(size(V,2),round(So.vertices(:,2) + vmat(2)))),...
          max(1,min(size(V,3),round(So.vertices(:,3) + vmat(3)))));
    Vv(I) = cdata;
    if opt.pve == 2 || opt.pve == 4
      [D,I] = vbdist(single(V )); Vv = Vv(I); clear D; 
      if opt.pve<4
        [D,I] = vbdist(single(~V)); Vv = Vv(I); clear D,
      end
    else
      Vv    = cat_vol_approx(Vv); 
    end
    
    % final masking
    if opt.pve > 3
      V    = Vv .* V; 
    else
      V    = Vv; 
    end 
    clear Vv; 
    
  elseif opt.pve == 0
    %%
    V    = false( round(max(S.vertices,[],1) - min(S.vertices)) + opt.bdist*2 );     
    vmat = -[min(S.vertices(:,1)) min(S.vertices(:,2)) min(S.vertices(:,3))] + opt.bdist; 
    I    = sub2ind(size(V),...
          max(1,min(size(V,1),round(S.vertices(:,1) + vmat(1)))),...
          max(1,min(size(V,2),round(S.vertices(:,2) + vmat(2)))),...
          max(1,min(size(V,3),round(S.vertices(:,3) + vmat(3)))));
    V(I) = 1; 

    V    = cat_vol_morph(V,'lc',1);  % closeing 
    V(I) = 0;                        % remove points of the surface
    V    = cat_vol_morph(V,'lab');   % final closing
  else %if opt.pve == 1
    %% fast PVE estimation by rendering multiple layer 
    
    Sn = spm_mesh_normals(S); 
    
    V    = zeros( round(max(S.vertices,[],1) - min(S.vertices)) + opt.bdist*2 ,'single');     
    vmat = -[min(S.vertices(:,1)) min(S.vertices(:,2)) min(S.vertices(:,3))] + opt.bdist; 
    
    offset = -0.25:0.25:1.0;
    for oi = 1:numel(offset)
      I = sub2ind(size(V),...
          max(1,min(size(V,1),round(S.vertices(:,1) + Sn(:,1)*offset(oi) + vmat(1)))),...
          max(1,min(size(V,2),round(S.vertices(:,2) + Sn(:,2)*offset(oi) + vmat(2)))),...
          max(1,min(size(V,3),round(S.vertices(:,3) + Sn(:,3)*offset(oi) + vmat(3)))));
      V(I) = min(1,max( V(I) , oi./numel(offset))); 
    end
    V(cat_vol_morph(V>=0.99,'lc',1) & V==0)=1;  % closeing 
  end
  vmati = repmat(min(S.vertices),size(S.vertices,1),1) - 5; 
  %%
  %SH = isosurface(V,0.6);        % create hull 
  %SH.vertices = [SH.vertices(:,2) SH.vertices(:,1) SH.vertices(:,3)]; % matlab flip
  %SH.vertices = SH.vertices + repmat(min(S.vertices),size(SH.vertices,1),1) - 5;
end
function alpha = angle(N1,N2)
  if 1 % fast version
    alpha = acosd( dot(N1,N2,2) ./ ( sum(N1.^2,2).^.5 .* sum(N2.^2,2).^.5 ));
  else
    %%
    alpha = zeros(size(N1,1),1);
    for i=1:size(N1,1)
      a = N1(i,:); b = N2(i,:); 
      alpha(i) = acosd( dot(a,b) / (norm(a) * norm(b)) ); 
    end
  end
end

function N = patchnormals(FV) % 
% Vertex normals of a triangulated mesh, area weighted, left-hand-rule 
% N = patchnormals(FV) - struct with fields, faces Nx3 and vertices Mx3 
% N: vertex normals as Mx3
%
% https://de.mathworks.com/matlabcentral/fileexchange/24330-patch-normals
% by Dirk-Jan Kroon

  %face corners index 
  A = FV.faces(:,1); 
  B = FV.faces(:,2); 
  C = FV.faces(:,3);

  %face normals 
  n = cross(FV.vertices(A,:)-FV.vertices(B,:),FV.vertices(C,:)-FV.vertices(A,:)); %area weighted

  %vertice normals 
  N = zeros(size(FV.vertices)); %init vertex normals 
  for i = 1:size(FV.faces,1) %step through faces (a vertex can be reference any number of times) 
    N(A(i),:) = N(A(i),:) + n(i,:); %sum face normals 
    N(B(i),:) = N(B(i),:) + n(i,:); 
    N(C(i),:) = N(C(i),:) + n(i,:); 
  end
end

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
        otherwise, varargout{1} = cat_surf_smoothtexture(S,varargin{:});
      end
    case 'maparea'    
      % do the final area projection
      if nargout==1, varargout{1} = cat_surf_maparea(S,varargin{:}); end
      if nargout==2, [varargout{1},varargout{2}] = cat_surf_maparea(S,varargin{:}); end
    case 'hull'
      if nargout==1, varargout{1} = cat_surf_hull(S); end
      if nargout==2, [varargout{1},varargout{2}] = cat_surf_hull(S); end
    case 'core'
      if nargout==1, varargout{1} = cat_surf_core(S,varargin{:}); end
      if nargout==2, [varargout{1},varargout{2}] = cat_surf_core(S,varargin{:}); end
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
          case 0, cat_surf_GMboundarySurface(action,S,varargin{:});
          case 1, varargout{1} = cat_surf_GMboundarySurface(action,S,varargin{:}); 
          case 2, [varargout{1},varargout{2}] = cat_surf_GMboundarySurface(action,S,varargin{:}); 
        end
      else % file input
        switch nargout
          case 0, cat_surf_GMboundarySurface(action,S);
          case 1, varargout{1} = cat_surf_GMboundarySurface(action,S); 
          case 2, [varargout{1},varargout{2}] = cat_surf_GMboundarySurface(action,S); 
        end
      end
    case 'evalcs'
      varargout{1} = cat_surf_evalCS(S,varargin{:});
    case 'createinneroutersurface'
      cat_surf_createinneroutersurface(S,varargin{:});
    case 'show_orthview'
      cat_surf_show_orthview(S,varargin{:});
    case 'saveico'
      cat_surf_saveICO(S,varargin{:});
    case 'collisioncorrection'
      [varargout{1},varargout{2},varargout{3}] = cat_surf_collision_correction(S,varargin{:});
    case 'vdist'
      [varargout{1},varargout{2}] = cat_surf_vdist(S,varargin);
    case 'surf2vol'
      if nargin>2
        [varargout{1},varargout{2},varargout{3}] = cat_surf_surf2vol(S,varargin{:});
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
      varargout{1} = cat_surf_surf2surf(S,varargin{:});
    case 'useedgemap'
      varargout{1} = cat_surf_maparea(S,varargin{:});
    case 'gmv'
      varargout{1} = cat_surf_gmv(S,varargin{:});
    case 'cdatamapping' 
      if nargin<3, varargin{3} = ''; end
      if nargin<4, varargin{4} = struct(); end
      if nargout>1
        [varargout{1},varargout{2}] = cat_surf_cdatamapping(S,varargin{:});
      else
        varargout{1} = cat_surf_cdatamapping(S,varargin{:});
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
function [IS,OS] = cat_surf_createinneroutersurface(S,T,Yp0)
  if ~exist('Yp0','var')
    % render surface
    
    
    Yp0 = 0; 
  end
  
  % call laplace 
  L = cat_surf_laplace(Yp0);
  
  % create streamlines
  IS.vertices = cat_surf_steams(L  ,T/2);
  OS.vertices = cat_surf_steams(1-L,T/2);
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
  
  if nargin>=2
    %% use filenames
    [pp,ff,ee] = spm_fileparts(varargin{1}); 
    
    if strcmp(ee,'')
      Praw = cat_io_FreeSurfer('fs2gii',varargin{1}); 
      Praw = Praw{1};
    else
      Praw   = varargin{1};
    end
    if nargin==3
      Pthick = varargin{2};
    else
      Pthick = cat_io_strrep(Praw,{'central','.gii'},{'pbt',''});
    end
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
function res = cat_surf_evalCS(CS,T,Ym,Ypp,Pcentral,verb)
% cat_surf_evalCS in cat_surf_fun
% _________________________________________________________________________
% print out some values
%   res = cat_surf_fun('evalCS',CS[,T,Ym,Yppt])
%   res = cat_surf_evalCS(CS[,T,Ym,Yppt])
%    
%   CS  .. central surface
%   T   .. cortical thickness
%   Ym  .. intensity normalized file with BG=0, CSF=1/3, GM=2/3, and WM=1
%   Ypp .. percent position map 
%   Pcentral .. number of classes for further thickness evaluation or a
%               given filename to detect specific thickness phantom rules
%   verb .. .print results
%
%   res .. structure with data fields of the printed values
% _________________________________________________________________________
% Robert Dahnke 201909

% - maybe also save and link (surface + hist) some files in future
% - the Layer4 handling with the global variables is horrible 

  global vmati mati vmat
 
  QMC       = cat_io_colormaps('marks+',17);
  color     = @(m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
  rms       = @(x) mean( x.^2 ).^0.5;
  rate      = @(x,best,worst) min(6,max(1, max(0,x-best) ./ (worst-best) * 5 + 1));
  
  if ~exist('verb','var'), verb = 1; end

  if 0 %isfield(CS,'vmat') && isfield(CS,'mati') 
    CS.vertices = (CS.vmat * [CS.vertices' ; ones(1,size(CS.vertices,1))])';
    if CS.mati(7)<0, CS.faces = [CS.faces(:,1) CS.faces(:,3) CS.faces(:,2)]; end
  end
  
  N  = spm_mesh_normals(CS);   % normalized surface normals                           
  M  = spm_mesh_smooth(CS);    % smoothing matrix
  if exist('T','var')
    VI = CS.vertices - N .* repmat(T/2,1,3); % white surface
    VO = CS.vertices + N .* repmat(T/2,1,3); % pial surface
  end
  
  if exist('Pcentral','var') && ischar(Pcentral)
    Player4   = strrep(Pcentral,'.central.','.layer4.'); 
    Ppbt      = strrep(Pcentral(1:end-4),'.central.','.pbt.'); 
    Pcentralx = strrep(Pcentral,'.central.','.centralx.'); 
    Player4x  = strrep(Pcentral,'.central.','.layer4x.'); 
    Ppbtx     = strrep(Pcentral(1:end-4),'.central.','.pbtx.'); 
  
    if ~isempty(vmati) && ~isempty(mati)
      if exist(Player4,'file')
        uL4 = 1; 
        L4 = gifti(Player4);
      else
   %     if ~exist(Pcentral,'file') && ~exist(Ppbt,'file') && ~isempty(vmat)
          CS1 = CS; CS1.vertices = (vmat*[CS.vertices' ; ones(1,size(CS.vertices,1))])'; 
          if mati(7)<0, CS1.faces = [CS1.faces(:,1) CS1.faces(:,3) CS1.faces(:,2)]; end
          save(gifti(struct('faces',CS1.faces,'vertices',CS1.vertices)),Pcentralx,'Base64Binary');
          cat_io_FreeSurfer('write_surf_data',Ppbtx,T);  
          cmd = sprintf('CAT_Central2Pial -equivolume -weight 1 "%s" "%s" "%s" 0', ...
                         Pcentralx,Ppbtx,Player4x);
        %else
        %  cmd = sprintf('CAT_Central2Pial -equivolume -weight 1 "%s" "%s" "%s" 0', ...
        %                 Pcentral,Ppbt,Player4x);
        %end
        try
          [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);
          L4 = gifti(Player4x);
          uL4 = 1; 
          delete(Player4x);
        catch
          uL4 = 1; 
        end
        if exist(Pcentralx,'file'), delete(Pcentralx); end
        if exist(Ppbtx,'file'), delete(Ppbtx); end
      end
    else
      uL4 = 0;
    end
  else
    uL4 = 0;
  end
  if uL4
    warning off MATLAB:subscripting:noSubscriptsSpecified
    VC = (vmati*[L4.vertices' ; ones(1,size(L4.vertices,1))])';
  end
  
  % local intensities
  if exist('Ym','var')
    II = isocolors2(Ym,VI); 
    IO = isocolors2(Ym,VO); 
    % local adaption for GM intensity changes by myelination 
    IIs = single(spm_mesh_smooth(M,double(II),round(100 * sqrt(size(CS.faces,1)/180000)))); 
    IOs = single(spm_mesh_smooth(M,double(II),round(100 * sqrt(size(CS.faces,1)/180000)))); 
    % normalization
    II  = II./(IIs/mean(IIs)) - 2.5/3; clear IIs;
    IO  = IO./(IOs/mean(IOs)) - 1.5/3; clear IOs;
    if uL4
      ML  = spm_mesh_smooth(L4);    % smoothing matrix
      IC  = isocolors2(Ym,VC); 
      ICs = single(spm_mesh_smooth(ML,double(IC),round(100 * sqrt(size(CS.faces,1)/180000)))); 
      IC  = IC./(ICs/mean(ICs)) - mean(ICs); clear ICs;
    end
    % output
    if verb
      fprintf('\n    Local intensity RMSE (lower=better): ')
      if uL4
        cat_io_cprintf( color( rate( mean( [rms(II),rms(IC),rms(IO)] ) , 0.05 , 0.20 )) , sprintf('%0.4f ',mean( [rms(II),rms(IC),rms(IO)] )) ); 
      else
        cat_io_cprintf( color( rate( mean( [rms(II),rms(IO)] ) , 0.05 , 0.20 )) , sprintf('%0.4f ',mean( [rms(II),rms(IO)] )) ); 
      end
      cat_io_cprintf( color( rate( rms(II) , 0.05 , 0.20 )) , sprintf('(IS=%0.4f,',rms(II)) ); 
      if uL4, cat_io_cprintf( color( rate( rms(IC) , 0.05 , 0.20 )) , sprintf('L4=%0.4f,',rms(IC)) ); end
      cat_io_cprintf( color( rate( rms(IO) , 0.05 , 0.20 )) , sprintf('OS=%0.4f)\n',rms(IO)) ); 
    end
    res.RMSE_Ym_white  = rms(II);
    if uL4, res.RMSE_Ym_layer4 = rms(IC); end
    res.RMSE_Ym_pial   = rms(IO);
    clear II IO; 
  end
  
  if exist('Ypp','var')
    % Yppi analysis
    II = isocolors2(Ypp,VI);          II = II - 1.0;
    IC = isocolors2(Ypp,CS.vertices); IC = IC - 0.5; 
    IO = isocolors2(Ypp,VO);          IO = IO - 0.0;
    % output
    if verb
      fprintf('    Local position  RMSE (lower=better): '); 
      cat_io_cprintf( color( rate( mean( [rms(IC),rms(II),rms(IO)]) , 0.05 , 0.30 )) ,sprintf('%0.4f ',mean( [rms(IC),rms(II),rms(IO)] )) ); 
      cat_io_cprintf( color( rate( rms(II) , 0.05 , 0.30 )) , sprintf('(IS=%0.4f,' ,rms(II)) ); 
      cat_io_cprintf( color( rate( rms(IC) , 0.05 , 0.30 )) , sprintf('CS=%0.4f,'  ,rms(IC)) ); 
      cat_io_cprintf( color( rate( rms(IO) , 0.05 , 0.30 )) , sprintf('OS=%0.4f)\n',rms(IO)) ); 
    end
    res.RMSE_Ypp_white   = rms(II);
    res.RMSE_Ypp_pial    = rms(IO);
    res.RMSE_Ypp_central = rms(IC);
  end 
  
  % CAT_SelfIntersect  surface_file output_values_file
  if exist('Pcentral','var') && ischar(Pcentral)
    [pp,ff,ee] = spm_fileparts(Pcentral);
    
    %Pselfw = fullfile(pp,strrep(ff,'central','whiteselfintersect'));
    %Pselfp = fullfile(pp,strrep(ff,'central','pialselfintersect'));
    
    %if ~exist(Pselfw,'file')  
      Pwhite = fullfile(pp,strrep([ff ee],'central','whitex'));   
      Ppial  = fullfile(pp,strrep([ff ee],'central','pialx'));   
      Pselfw = fullfile(pp,strrep(ff,'central','whiteselfintersect'));
      Pselfp = fullfile(pp,strrep(ff,'central','pialselfintersect'));
    
      % save surfaces
      save(gifti(struct('faces',CS.faces,'vertices',VI)),Pwhite,'Base64Binary');
      save(gifti(struct('faces',CS.faces,'vertices',VO)),Ppial,'Base64Binary');

      % write self intersection maps
      %tic
      cmd = sprintf('CAT_SelfIntersect "%s" "%s"',Pwhite,Pselfw); 
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);
      cmd = sprintf('CAT_SelfIntersect "%s" "%s"',Ppial,Pselfp); 
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);
      %toc
    
    %end
    
      selfw = cat_io_FreeSurfer('read_surf_data',Pselfw);
      selfp = cat_io_FreeSurfer('read_surf_data',Pselfp);

      area = cat_surf_area(CS);
      
      res.white_self_interection_area = sum((selfw(:)>0) .* area(:)) / 100;
      res.pial_self_interection_area  = sum((selfp(:)>0) .* area(:)) / 100;
      res.white_self_interections     = res.white_self_interection_area / sum(area(:)/100) * 100;
      res.pial_self_interections      = res.pial_self_interection_area  / sum(area(:)/100) * 100;
      
      if verb
        fprintf('    Self intersections (white,pial):     '); 
        cat_io_cprintf( color( rate( res.white_self_interections , 0 , 50 )) , ...
          sprintf('%0.2f%%%% (%0.2f cm%s) ',res.white_self_interections,res.white_self_interection_area,char(178))); 
        cat_io_cprintf( color( rate( res.pial_self_interections , 0 , 50 )) , ...
          sprintf('%0.2f%%%% (%0.2f cm%s)\n',res.white_self_interections,res.pial_self_interection_area,char(178))); 
      end

      delete(Pwhite);
      delete(Ppial); 
      delete(Pselfw);
      delete(Pselfp); 

  end
  
  % thickness analysis
  if exist('T','var')
    if exist('Tclasses','var') && ~isempty(Pcentral)
      if ischar(Pcentral)
        if strfind(Pcentral,'dilated1.5-2.5mm')
          T = cat_stat_histth(T,0.95);
          Pcentral = 3;
        else
          Pcentral = 0;
        end
      end
      
      if Pcentral>0
        if Pcentral>7 || Pcentral<1
          warning('Tclasses has to be between 2 and 7.');
          Pcentral = min(7,max(3,Pcentral)); 
        end
        [mn,sd] = kmeans3D(T,Pcentral);
        if verb
          fprintf('    Thickness mean (%d class(es)):       ',Pcentral)
          fprintf('%7.4f',mn); fprintf('\n'); 
          fprintf('    Thickness std  (%d class(es)):       ',Pcentral)
          fprintf('%7.4f',sd); fprintf('\n');
        end
        res.thickness_mean_nclasses = mn;
        res.thickness_std_nclasses  = sd;
      end
    end
      
    res.thickness_mn_sd_md_mx = [mean(T),std(T),median(T),max(T)];
    if verb
      fprintf('    Thickness values:                    %0.4f%s%0.4f (md=%0.4f,mx=%0.4f)\n',...
        res.thickness_mn_sd_md_mx(1),char(177),res.thickness_mn_sd_md_mx(2:end)); 
    end
  end
  
  % curvature analyis
  if 0
    C = abs(spm_mesh_curvature(CS));
    res.abscurv_mn_sd_md_mx = [mean(C),std(C),median(C),max(C)];
    if verb
      fprintf('    Abs mean curvature values:           %0.4f%s%0.4f (md=%0.4f,mx=%0.4f)\n',...
        res.abscurv_mn_sd_md_mx(1),char(177),res.abscurv_mn_sd_md_mx(2:end)); 
    end
  end
  
  % surface values
  EC  = size(CS.vertices,1) + size(CS.faces,1) - size(spm_mesh_edges(CS),1);
  res.euler_characteristic = EC; 
  if verb
    fprintf('    Faces / Final Euler number:          '); 
    cat_io_cprintf( color( rate( 1 - max(0,size(CS.faces,1)/300000) , 0 , 0.9 )) , sprintf('%d / ',size(CS.faces,1)));
    cat_io_cprintf( color( rate( abs(EC-2) , 0 , 30 )) , sprintf('%d',EC));
    fprintf('\n'); 
  end
end
function cat_surf_saveICO(S,Tpbt,Pcs,subdir,writeTfs,Pm,C)
% _________________________________________________________________________
% Save surface data for debuging:
% Creates and save the white and pial surfaces based on the displacement by
% the half thickness along the surface normals and use the inner and outer
% surfaces to create the layer4 surface.
% Saves also the thickness file.
%
%   cat_surf_saveICO(S,Tpbt,Pcs,subdir,writeTfs[,C,Pm])
%
%   S         .. central surface
%   T         .. cortical thickness
%   Pcs       .. central surface file name (with full path)
%   subdir    .. addition subdirecty in the standard surface directory
%   writeTfs  .. estimate FreeSurfer thickness metric
% [in development]
%   C         .. write further surface data
%   Pm        .. intensity file/volume to map data to the surfaces
% _________________________________________________________________________
% Robert Dahnke 201908

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
 
  % surface filenames
  Pcentral = fullfile(pp,subdir,[ff ee]);   
  Pwhite   = fullfile(pp,subdir,strrep([ff ee],'central','white'));   
  Ppial    = fullfile(pp,subdir,strrep([ff ee],'central','pial'));   
  Pthick   = fullfile(pp,subdir,strrep(ff,'central','thickness'));   
  Ppbt     = fullfile(pp,subdir,strrep(ff,'central','pbt'));   
  PintIS   = fullfile(pp,subdir,strrep(ff,'central','Ym-white'));
  PintOS   = fullfile(pp,subdir,strrep(ff,'central','Ym-pial'));
  PintL4   = fullfile(pp,subdir,strrep(ff,'central','Ym-L4'));
  Pcol     = fullfile(pp,subdir,strrep(ff,'central','collision'));   
  Player4  = fullfile(pp,subdir,strrep([ff ee],'central','layer4'));   
  Pselfw   = fullfile(pp,subdir,strrep([ff ee],'central','whiteselfintersect'));
  Pselfp   = fullfile(pp,subdir,strrep([ff ee],'central','pialselfintersect'));

  % save surfaces
  save(gifti(struct('faces',S.faces,'vertices',S.vertices)),Pcentral,'Base64Binary');
  save(gifti(struct('faces',S.faces,'vertices',VI)),Pwhite,'Base64Binary');
  save(gifti(struct('faces',S.faces,'vertices',VO)),Ppial,'Base64Binary');

  % write self intersection maps
  cmd = sprintf('CAT_SelfIntersect "%s" "%s"',Pwhite,Pselfw); 
  [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);
  cmd = sprintf('CAT_SelfIntersect "%s" "%s"',Ppial,Pselfp); 
  [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);
  
  % save thickness
  cat_io_FreeSurfer('write_surf_data',Ppbt,Tpbt);
  if exist('writeTfs','var') && ~isempty(writeTfs) && writeTfs
    cmd = sprintf('CAT_SurfDistance -mean "%s" "%s" "%s"',Pwhite,Ppial,Pthick);
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-2);
    fprintf('Display thickness: %s\n',spm_file(Pthick ,'link','cat_surf_display(''%s'')'));
  end
  if exist('C','var') && ~isempty(C)
    cat_io_FreeSurfer('write_surf_data',Pcol,C);
  end
  
  % final correction of central surface in highly folded areas with high mean curvature
  cmd = sprintf('CAT_Central2Pial -equivolume -weight 1 "%s" "%s" "%s" 0', ...
                   Pcentral,Ppbt,Player4);
  [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);
  
  % save intensities
  if ~exist('Pm','var')
    % volume filenames for spm_orthview
    sinfo = cat_surf_info(Pcentral);
    
    if cat_get_defaults('extopts.subfolders')
      Pm = fullfile(spm_str_manip(pp,'h'),'mri',['m' sinfo.name '.nii']); 
    else
      Pm = fullfile(pp,['m' sinfo.name '.nii']); 
    end
    if ~exist(Pm,'file')
      Pm = fullfile(spm_str_manip(pp,'hh'),'mri',['m' sinfo.name '.nii']); 
    end
    if ~exist(Pm,'file')
      Pm = fullfile(spm_str_manip(pp,'h'),[sinfo.name '.nii']); 
    end
    if ~exist(Pm,'file')
      Pm = ''; 
    end
  end
  
  if ~isnumeric( Pm ) && exist(Pm,'file')
    % use the file data
    cmd = sprintf('CAT_3dVol2Surf -cubic -steps 1 -start 0 -end 0 "%s" "%s" "%s"',Pwhite , Pm, PintIS);
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);
    cmd = sprintf('CAT_3dVol2Surf -cubic -steps 1 -start 0 -end 0 "%s" "%s" "%s"',Ppial  , Pm, PintOS);
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);
    cmd = sprintf('CAT_3dVol2Surf -cubic -steps 1 -start 0 -end 0 "%s" "%s" "%s"',Player4, Pm, PintL4);
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);
  elseif ndims(Pm)==3
    % use a given volume 
    SL = gifti(Player4);
    warning off MATLAB:subscripting:noSubscriptsSpecified
    SL.vertices = (S.vmati*[SL.vertices' ; ones(1,size(SL.vertices,1))])';
    if S.mati(7)<0,  SL.faces = [SL.faces(:,1) SL.faces(:,3) SL.faces(:,2)]; end
    VL = SL.vertices;

    int = isocolors(Pm,VO); cat_io_FreeSurfer('write_surf_data',PintOS,int);
    int = isocolors(Pm,VI); cat_io_FreeSurfer('write_surf_data',PintIS,int);
    int = isocolors(Pm,VL); cat_io_FreeSurfer('write_surf_data',PintL4,int);
  end
  
  % display something to click
  fprintf('\n    Display surface:  %s\n',spm_file(Ppbt  ,'link','cat_surf_display(''%s'')'));
  if ~isempty(Pm)
    fprintf('    Show in orthview: %s\n',spm_file(Pm ,'link',...
      [ sprintf('cat_surf_fun(''show_orthview'',{''%s'';''%s'';''%s''},',Pcentral,Ppial,Pwhite) '''%s'')']));
  end
  
  if 0
    if ndims(Pm)==3, Ym=Pm; else, Ym=spm_read_vols(spm_vol(Pm)); end
    res = cat_surf_evalCS(S,Tpbt,Ym,Ypp,Tclasses)
  end
end
function [SN,TN,E] = cat_surf_collision_correction(S,T,Y,Ypp,Yl4,opt) 
% _________________________________________________________________________
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
%   [SN,TN,E] = cat_surf_collision_correction(S,T,Y[,debug,E,Pcs])
%
%   SN      .. new surface
%   TN      .. new thickness
%   S       .. original surface
%   T       .. original thickness
%   Y       .. segmentation map or intensity normalized images 
%              for intra/inter surface edge definition
%   Yl4     .. layer4 intensity surface map
%   opt     .. parameter structure
%    .debug .. option to write the un- and corrected cortical surfaces in a
%              subdirectory
%    .verb  ..
%    .E     .. Delaunay edge map (from previous run or empty matrix as input)
%    .Pcs   .. central surface file name to write debugging files
%    ... experimental settings
%    .smoothCSinput  ..
%    .model          ..
%    .PVEcorr        ..
%    .slowdown       ..
%    .
% _________________________________________________________________________
% Robert Dahnke 201909


% -------------------------------------------------------------------------
% Todo:
% - full support of parameters & recall 
% - support of different correction models 
%   (e.g. only collision vs. intensity correction)
% - full documentation and detailed comments
% - helping boundary surfaces? 
%   > partially implemented 
%   > very slow +20-60s for each just for creation  
% - stable subset/list (also as internal/external error measure
% - use of gradients and divergence rather simple intensity information
% - use of Ypp
% - optimization of the layer 4 (layer concept)
% - face-flipping correction that can not be handled by Delaunay because of
%   its neighbor limits
% - improved evaluation concept
% - improved validation concept
% - fast mapping c-function for edge to surf that combine multiple values
%   given by an index map by different functions (mean,min,max,std)
% - surface filter sub-function to remove outlier
% - triangle height rather than edge distance (or combination)
% -------------------------------------------------------------------------

  if ~exist('opt','var'), opt = struct(); end 

  % default variables 
  def.Pcs               = '';     % filename to write debugging output data
  def.debug             = 1;      % debugging output vs. memory optimization 
  def.verb              = 1;      % display debugging information 
  def.E                 = [];     % inter-surface Delaunay edges (of a previous run) 
  def.boundarySurfaces  = 0;      % use inner and outer boundary surface to improve the Delaunay graph 
  def.smoothCSinput     = 0;      % smooth the input CS for more stable Delaunay triangulation in case of locally oversampled surfaces 
  def.PVEcorr           = 1;      % correction of PVE values for 2 boundaries in on voxel (experimental)
  def.slowdown          = 1;      % slowdown may stabilize the process over the iterations   
  def.model             = 2;      % 0 - only collcorr, 1 - only intopt, 2 - both  
  def.vx_vol            = 1; 
  opt                   = cat_io_checkinopt(opt,def); clear def;
  def.write             = opt.debug & ~isempty(opt.Pcs);
  opt                   = cat_io_checkinopt(opt,def); 
  
  
  % helping smoothing functions for data and surfaces
  M   = spm_mesh_smooth(S);         % for spm_smoothing matrix
  rms = @(x) mean( x.^2 ) .^ 0.5;   % for error handling of mad vertices
  smoothsurf = @(V,s) [ ...         % simple surface smoothing 
    spm_mesh_smooth(M,double(V(:,1)),s) , ...
    spm_mesh_smooth(M,double(V(:,2)),s) , ...
    spm_mesh_smooth(M,double(V(:,3)),s) ];
  
  
  % detection and correction for flipped faces to have always the same normal direction
  lim = 1:round(size(S.vertices,1)/1000):size(S.vertices,1); 
  N   = spm_mesh_normals(S);   
  VOl = S.vertices(lim,:) - N(lim,:) .* repmat(T(lim)/2,1,3); 
  VIl = S.vertices(lim,:) + N(lim,:) .* repmat(T(lim)/2,1,3); 
  YOl = isocolors2(Y,VOl); 
  YIl = isocolors2(Y,VIl); 
  flipped = mean(YOl) > mean(YIl); 
  clear N VOl VIl YOl YIl lim; 
  if flipped, S.faces = [S.faces(:,1) S.faces(:,3) S.faces(:,2)]; S.mati(7) = - S.mati(7); end
  
  
  % larger surface need more smoothing to avoid triangulation problems 
  sf = round( sqrt( size(S.faces,1) / 50000) );  % ### empirical value 
  if max(Y(:))<1.5, Y = Y.*2+1; else, Y = max(1,Y); end              
  if opt.debug, fprintf('\n'); end
  if opt.write, cat_surf_saveICO(S,T,Pcs,sprintf('pre_collcorr_%0.0fk',round( size(S.faces,1)/1000 / 10) * 10 ),0); end
  stime = cat_io_cmd(sprintf('    Delaunay triangulation of %d vertices (sf=%d):',size(S.vertices,1),sf),'g5','',opt.debug); 

  
  

  
  
  %% Creation of the inter surface edges based on a Delaunay graph 
  %  ----------------------------------------------------------------------
  %  There is a short cut to apply further iterations without processing
  %  the graph again. 
  %  ----------------------------------------------------------------------
  if isfield(opt,'E') && isempty(opt.E)
    
    % Early versions used a smoothed surface to reduce problems due to
    % artifacts. However, the improved input meshed (createCS2 pipeline)
    % do not need smoothing and surface smoothing can not be combined 
    % with helping surfaces (inner and outer surface points). 
    VS = double(S.vertices); 
    if opt.smoothCSinput
      % Surface smoothing as loop to correct for outlier due to incorrect surfaces.
      % Using the smoothing directly create some extrem large spikes - don't know why (RD20190912).
      % The smoothing is not required in newer version (RD20190922)
      for i = 1:opt.smoothCSinput                           
        VSi = smoothsurf(VS,2); 
        VM  = rms(VSi - VS)<2; 
        VS(VM,:) = VSi(VM,:); 
      end
      VS = smoothsurf(VS,1);
    end
    if ~opt.debug, clear VSi VM; end

    
    % helping boundary surface - (uncorrected) inner or/and outer surface 
    % this was slow (20-60 seconds) and did not work so simple/fast ... need further work
    % maybe just use low resolution surfaces (1 mm)
    if opt.boundarySurfaces
      if opt.smoothCSinput, error('Initial surface smoothing "smoothCSinput" can not be combined with "boundarySurfaces".\n'), end
      VB = S.vertices;                             
      if opt.boundarySurfaces == 1 || opt.boundarySurfaces == 3
        UOS = isosurface(Y,1.5); 
        VS  = [ VS ; UOS.vertices ]; 
        VB  = [ VB ; UOS.vertices ]; 
        if ~opt.debug, clear UOS; end 
      end
      if opt.boundarySurfaces == 2 || opt.boundarySurfaces == 3
        UIS = isosurface(Y,1.5); 
        VS  = [ VS ; UIS.vertices ]; 
        VB  = [ VB ; UIS.vertices ]; 
        if ~opt.debug, clear UIS; end 
      end
    end
    
    
    % Delaunay graph
    D  = single(delaunayn( VS )); 
    if ~opt.debug, clear VS; end            

    % decompose delaunay graph into its edges
    E  = uint32(cat_surf_edges(D));   
    nE = size(E,1); 
    if ~opt.debug, clear D; end
    if opt.debug, cat_io_cprintf('g5',sprintf('%5.0fs\n',etime(clock,stime))); end


    % separate helping boundary surfaces
    if opt.boundarySurfaces == 2 || opt.boundarySurfaces == 3
      Etmep = sum( E>( size(S.vertices,1) + any(opt.boundarySurfaces==[1,3]).*size(UOS.vertices,1) ) , 2 )>0; 
      EUIS  = E( Etmep , : ); 
      EUIS( sum( EUIS>numel(S.vertices) , 2 )~=1, :) = []; % remove all edges that are not between the CS and the UIS
      EUIS  = sort(EUIS,2); 
      E( Etmep , : ) = [];
    end
    if opt.boundarySurfaces == 1 || opt.boundarySurfaces == 3
      Etmep = sum( E>( size(S.vertices,1) ) , 2 )>0; 
      EUOS  = E( Etmep , : ); 
      EUOS( sum( EUOS>numel(S.vertices) , 2 )~=1, :) = []; % remove all edges that are not between the CS and the UOS
      EUOS  = sort(EUOS,2); 
      E( Etmep , : ) = []; 
    end
    

    %% Remove intra-surface edges 
    %  --------------------------------------------------------------------
    %  If we remove too much then the correction will not work.
    %  If we do not remove enough then it will add sulci in regions without sulci

    V  = S.vertices;

    % remove edge that we know from the surface - super save
    stime = clock; 
    EF = uint32(cat_surf_edges(S.faces));           
    E  = setdiff(E,EF,'rows'); clear EF; 
    
    % remove edges between neigbors of each point - relative save
    % get neighbor matrix
    [NE,MED] = spm_mesh_neighbours(M); 
    nNE = size(NE,1);

    % extra element that link on itself and replace the 0 in NE that does not allow matrix indexing 
    NIL = nNE + 1; 
    NE(NIL,:)  = NIL * ones(1,size(NE,2)); NE(NE==0) = nNE+1;
    
    % further levels
    % use higher levels only for large surfaces (use sqrt to compensate area grow factor)
    nlevel = 2; % max(2,round( 1 + sqrt( ceil( size(S.faces,1) / 300000 ) ))); 
    for nli = 1:nlevel % nice idear but not working yet
      for ni = 1:size(NE,2)
        NEN = sum( NE == repmat( NE(NE(:,ni),ni) , 1 , size(NE,2) ),2)>0; 
        NE  = [ NE  min( NIL , NE(:,ni)  + NIL*NEN) ];  %#ok<AGROW>
        MED = [ MED MED(:,ni) ];  %#ok<AGROW>
      end
      
      % sort entries
      [NE,NEsi] = sort(NE,2); MED = MED(NEsi); clear NEsi
      NILi = min([ size(NE,2) , find( sum(NE == NIL,1) >= size(NE,1)*0.5 , 1, 'first') - 1]);
      NE  = NE(:,1:NILi); 
      MED = MED(:,1:NILi);  
    end
    NE(NE==NIL) = 0; 

    % remove edges from the neigbor list
    for i=2:size(NE,2)
      E = setdiff(E,[NE(:,1) NE(:,i)],'rows'); 
      E = setdiff(E,[NE(:,i) NE(:,1)],'rows'); 
    end
    clear NE
    if opt.debug
      cat_io_cprintf('g5',sprintf('    remove edges by surface (l%d):%8d > %9d (%0.2f%%%%) %9.0fs',...
        nlevel,nE,size(E,1),size(E,1)./nE,etime(clock,stime)));
    else
      clear nE
    end


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
               angle(NS(E(:,2),:), NS(E(:,1),:))]; clear NS
    SNalpha = [angle(N(E(:,1),:),  V(E(:,1),:) - V(E(:,2),:)), ...
               angle(N(E(:,2),:),  V(E(:,2),:) - V(E(:,1),:))]; 
    NEna    = mean(Nalpha/180,2); clear Nalpha                       % figure, hist( NEna , 0:0.01:1);
    NEsa    = (abs(90  - SNalpha)/90  + abs(90  - SNalpha)/90)/2;    % figure, hist( NEsa , 0:0.01:1);
    clear SNalpha; 

    % remove by intensity given by the centroids of the edges
    VC  = cat_surf_centroid(V,E); 
    IC  = isocolors2(Y,VC); clear VC;
    % outer surface intensity
    VO  = V - N .* repmat(T/2,1,3); 
    VOC = cat_surf_centroid(VO,E); 
    IO  = isocolors2(Y,VOC); clear VOC VO; 
    % inner surface intensity
    VI  = V + N .* repmat(T/2,1,3) + 0.1; % GM/WM  
    VIC = cat_surf_centroid(VI,E); 
    II  = isocolors2(Y,VIC); clear VIC VI; 
    VI  = V + N .* repmat(T/2,1,3) + 0.5; % save WM  
    VIC = cat_surf_centroid(VI,E); 
    II  = max(II,isocolors2(Y,VIC)); clear VIC VI; % use max to get WM value 
    VI  = V + N .* repmat(T/2,1,3) + 1.0; % supersave WM  
    VIC = cat_surf_centroid(VI,E); 
    II  = max(II,isocolors2(Y,VIC)); clear VIC VI; % use max to get WM value 
    % combine all intensities 
    NEi = 1 - min(1,max(abs(diff([II IC IO],1,2)),[],2)); 
    %ET  = mean([II IC IO],2)>2.25; % edge classification 
    if ~opt.debug, clear II IC IO; end

    % combine all measures by product to remove many things
    % I also though about an adaptive threshold but it is not so easy ...
    NE = prod( [NEd NEi NEna*2 NEsa*2] ,2); % 1.75 % larger values > remove less
    NE = NE < .05; %05; %max(eps,mean(NE) - 1*std(NE)); % smaller values > remove less
    E (NE,:) = []; %if exist('ET','var'), ET(NE) = []; end
    if opt.debug
      cat_io_cprintf('g5',sprintf('\n    remove edges by intensity:            > %9d (%0.2f%%%%) %9.0fs',...
        size(E,1),size(E,1)./nE,etime(clock,stime))); stime = clock;
    else
      clear NE NEd NEi NEna NEsa
    end

   %fprintf('\nsf = %0.2f',sf);
  else
    N  = spm_mesh_normals(S);   
  end
  
  
  
  if opt.debug
    cat_io_cprintf('g5','\n    Prepare Optimization:'); stime = clock; 
  end
  
  %% updated measures
  SNalpha = [angle(N(E(:,1),:),  V(E(:,1),:) - V(E(:,2),:)), ...
             angle(N(E(:,2),:),  V(E(:,2),:) - V(E(:,1),:))]; 

  VC  = cat_surf_centroid(V,E); 
  IC  = isocolors(Y,VC); clear VC;

  OE  = min(1,(min(SNalpha,[],2)<90) + IC<2.15);
  IE  = min(1,(max(SNalpha,[],2)>90) + IC>2.15);
      
  if 1
    % avoid PBT overestimation in gyri (well thickness is correct but
    % measures non-linear/non-orthogonal)
    TN = single(spm_mesh_smooth(M,double(T), sf * 20 )); 
    T  = min(T,TN); 
  end
  
  TN = T; SN = S; TCsum = 0; %#ok<NASGU>
  TCsumo = inf; TNold = inf;
  maxiter  = 20;  % main number of iterations 
  maxiter2 = 20;  % limit of adapting the mixing model
  
  Yl4 = single(spm_mesh_smooth(M,double(Yl4),sf/4 * 100));
  
  % I did not manage to use curvate ...
  %C   = spm_mesh_curvature(S); 
  %C   = spm_mesh_smooth(M,C,1);

  % PVE doubleside correction: 
  % If a voxel contain 38% GM and 62% CSF and has one boundary, it is approximately at the 38% position of the voxel.  
  % If the same voxel contain two boundaries, the each boundary is approximately at the 19% position of that voxel.    
  % Hence, I try to measure the filling effect in regions of two boundaries by the local minimum/maximum to estimate 
  % double the PVE effect (like a sharpening).
  % However, this is relative slow ...
  if opt.PVEcorr 
    Ypvec = cat_vol_localstat(max(1,cat_vol_localstat(min(2,max(1,Y)),Y>1,2,3)),Y>1,2,2);
    Ypvew = cat_vol_localstat(max(2,cat_vol_localstat(min(3,max(2,Y)),Y>1,2,2)),Y>1,2,3);
    Y = max(1,min(3,Y - ((max(Ypvec,Y))-Y) + (Y-(min(Ypvew,Y)))));
  end
  
  if opt.debug
    stime = cat_io_cmd('  Optimize surface:','g5','',opt.verb,stime); fprintf('\n');
  end
      
  %% Iterative correction routine
  for j=1:maxiter+1
    V   = single(SN.vertices);

    % update surface normales  
    N   = spm_mesh_normals(SN); 
    
    % inner and outer surface
    VO  = V - N .* repmat(TN/2,1,3);
    VI  = V + N .* repmat(TN/2,1,3);

    
    % First correction step that works but also could be improved.
    % ---------------------------------------------------------------------
    % Complex side specific correction by the inter-surface edges, that 
    % used the angle between the edges and normals to define edges within 
    % a suclus (outer) or within a gyrus (inner). 
    % In general, only inter-surface edges are expected here, those
    % distance describes the maximal local thickness.  We also add some 
    % sulcus-width to avoid collisions but the effect will be small due 
    % to the smoothing.
    % There are problems that not all points have a inter-surface edge, 
    % so it is necessary to smooth to include unconnected neighbors  
    % LEC, LEOC, and LEIC represent the distance error by collisions
    % of each edge.
    % ---------------------------------------------------------------------
    
    % edgelength of the central, inner, and outer surface
    % ###
    %   The edgelength is just the simples measure - the high of the
    %   tetraeder would be more excact.
    % ###
    LE  = sum( (V(E(:,1),:)  - V(E(:,2),:)).^2  , 2) .^ 0.5;          % distance between the central surface (~ thickness/2 + thickness/2)
    LEO = sum( (VO(E(:,1),:) - VO(E(:,2),:)).^2 , 2) .^ 0.5;          % distance between the outer surface (~?minimal/maximal distance) 
    LEI = sum( (VI(E(:,1),:) - VI(E(:,2),:)).^2 , 2) .^ 0.5;          % distance between the inner surface (~ minimal/maximal distance)
    if opt.boundarySurfaces == 1 || opt.boundarySurfaces == 3
      LEUOS = sum( (VB(EUOS(:,1),:)  - VB(EUOS(:,2),:)).^2  , 2) .^ 0.5;  % distance b
    end
    if opt.boundarySurfaces == 2 || opt.boundarySurfaces == 3
      LEUIS = sum( (VB(EUIS(:,1),:)  - VB(EUIS(:,2),:)).^2  , 2) .^ 0.5;
    end

    % estimate error for each Delaunay edge 
    % (sum local thickness and sulcuswidth vs. the length of the edge) 
    %sulcuswidth = 0.0; % worse results with additional width
    LECP = max(0, LE - ( TN(E(:,1))/2 + TN(E(:,2))/2 + 0.02 ) ) / 2; % - sulcuswidth )); 
    LEC  = max(-inf, (TN(E(:,1))/2 + TN(E(:,2))/2) - (LE - 0.02) ) / 2; % - sulcuswidth )); 
    LEOC = LEO .* max(-inf, 0.02 - LEO) / 2; %clear LEO; % minimum distance between points (rare spaecial case) 
    LEIC = LEI .* max(-inf, 0.02 - LEI) / 2; %clear LEI; % minimum distance between points (rare spaecial case)
    TNP  = repmat(TN,3,1); %2 + any(opt.boundarySurfaces == [1,3]) + any(opt.boundarySurfaces == [2,3]) ,1);
    if opt.boundarySurfaces == 1 || opt.boundarySurfaces == 3
      LEUOC = max(0, LEUOS - ( TNP(EUOS(:,1))/2 ) );
    end
    if opt.boundarySurfaces == 2 || opt.boundarySurfaces == 3
      LEUIC = max(0, LEUIS - ( TNP(EUIS(:,1))/2 ) );
    end
    
    
    
    % map the Delaunay edge correction to the vertices (simple maximum)
    % ###
    %   You may (also) use some intensity information here! ... added
    %   Moreover, a loop is very slow and the estimation of a mapping
    %   would be better! But how? ... partially implemented
    % ###
    %{
            EOid = TN*0; EIid = TN*0; 
            for ni=1:size(E,1)
              EOid(E(ni,1)) = EOid(E(ni,1)) .* , LEC(ni) .* OE(ni), LEOC(ni)]); 
              EOid(E(ni,2)) = max([EOid(E(ni,2)), LEC(ni) .* OE(ni), LEOC(ni)]); 
              EIid(E(ni,1)) = max([TIC(E(ni,1)), LEC(ni) .* IE(ni), LEIC(ni)]); 
              EIid(E(ni,2)) = max([TIC(E(ni,2)), LEC(ni) .* IE(ni), LEIC(ni)]); 
            end
    %}
    OE = min(1,(min(SNalpha,[],2)<60) + IC<2.15);
    IE = min(1,(max(SNalpha,[],2)>60) + IC>2.15); 
        
    TOC = TN*0; TIC = TN*0; TOCP = TN*0; TICP = TN*0; %PVE_LEOC = TN*0; PVE_LEIC = TN*0;
    app = 1; 
    if app == 1
      for ni=1:size(E,1)
       % OE = min(1,(min(SNalpha(ni,:))<60) + IC(ni)<2.15);
       % IE = min(1,(max(SNalpha(ni,:))>60) + IC(ni)>2.15); 
        
        %{
        PVE_LEOC(E(ni,1)) = max(0,opt.vx_vol - LEOC(ni)); 
        PVE_LEOC(E(ni,2)) = max(0,opt.vx_vol - LEOC(ni)); 
        PVE_LEIC(E(ni,1)) = max(0,opt.vx_vol - LEIC(ni)); 
        PVE_LEIC(E(ni,2)) = max(0,opt.vx_vol - LEIC(ni)); 
        %}
        
        TOC(E(ni,1)) = max([TOC(E(ni,1)), LEC(ni) .* OE(ni)]); %, -inf*LEOC(ni)]); 
        TOC(E(ni,2)) = max([TOC(E(ni,2)), LEC(ni) .* OE(ni)]); %, -inf*LEOC(ni)]); 
        TIC(E(ni,1)) = max([TIC(E(ni,1)), LEC(ni) .* IE(ni), LEIC(ni)]); 
        TIC(E(ni,2)) = max([TIC(E(ni,2)), LEC(ni) .* IE(ni), LEIC(ni)]); 
        
        % with angle weighting ...
        TOCP(E(ni,1)) = max([TOCP(E(ni,1)),LECP(ni) .* OE(ni)]); 
        TOCP(E(ni,2)) = max([TOCP(E(ni,1)),LECP(ni) .* OE(ni)]); 
        TICP(E(ni,1)) = max([TOCP(E(ni,1)),LECP(ni) .* IE(ni)]); 
        TICP(E(ni,2)) = max([TOCP(E(ni,1)),LECP(ni) .* IE(ni)]); 
      end
    else
      for ni=1:size(E,1)
        TOC(E(ni,1)) = max([TOC(E(ni,1)), LEC(ni) .* OE(ni), LEOC(ni)]); 
        TOC(E(ni,2)) = max([TOC(E(ni,2)), LEC(ni) .* OE(ni), LEOC(ni)]); 
        TIC(E(ni,1)) = max([TIC(E(ni,1)), LEC(ni) .* IE(ni), LEIC(ni)]); 
        TIC(E(ni,2)) = max([TIC(E(ni,2)), LEC(ni) .* IE(ni), LEIC(ni)]); 
      end
    end
    if opt.boundarySurfaces == 1 || opt.boundarySurfaces == 3
      for ni=1:size(EUOS,1)
         TOC( mod( EUOS(ni,1)-1 , size(SN,1) )+1) = LEUOC(ni); 
      end
    end
    if opt.boundarySurfaces == 2 || opt.boundarySurfaces == 3
      for ni=1:size(EUOS,1)
         TOC( mod(EUIS(ni,1)-1, size(SN,1) )+1) = LEUIC(ni); 
      end
    end
    clear LEC LEOC LIOC;
    
    if opt.slowdown
      slowdown  = max(1,2/j); 
    else
      slowdown  = 1; 
    end
      
    if opt.model == 0 || opt.model == 2
      TOC  = single( spm_mesh_smooth(M,double(TOC), 1 ))*1.4;   TOC = TOC / (slowdown/2); 
      TIC  = single( spm_mesh_smooth(M,double(TIC), 1 ))*1.2;   TIC = TIC / (slowdown/2); 
      if opt.verb, fprintf('\n  TIC: %0.2f%s%0.2f, TOC: %0.2f%s%0.2f',mean(TIC),char(177),std(TIC),mean(TOC),char(177),std(TOC)); end
    end
    %%
    if opt.model
      % filter limits
      TOCP  = single( spm_mesh_smooth(M,double(TOCP), 1 ))*1.0;% 1.5  
      TICP  = single( spm_mesh_smooth(M,double(TICP), 1 ))*0.8;

      % correction for intensities ...
      YI    = isocolors2(Y,VI); 
      YO    = isocolors2(Y,VO);  
      YppO  = isocolors2(Ypp,VO);  
     
      if opt.model == 1, fprintf('\n'); end
      if opt.verb, fprintf('  YIC: %0.2f%s%0.2f, YOC: %0.2f%s%0.2f',mean(YI),char(177),std(YI),mean(YO),char(177),std(YO)); end 

      WMth  = 3; YI   = max( -TICP , max(-1, min(0.5, YI - ((WMth/2 + Yl4/2) )  ))  ) / (slowdown);
      CSFth = 1; YO   = max( -TOCP , max(-1, min(0.5, ((CSFth/2 + Yl4/2) ) - YO ))  ) / (slowdown);% + 2*C
      CSFth = 0; Yppc = max( -TOCP , max(-0.05, min(0.05, 0.01 - YppO ))  ) / (slowdown);% + 2*C
      
      if 1
        YC = isocolors2(Y,V ); 
        YC = max( -0.5, min( 0.5, YC - Yl4 )) / (slowdown);
        YI = YI * 0.8 + 0.2 * YC; 
        YO = (YO * 0.8 - 0.2 * YC) .* min(1,YppO*20); 
      end
      YO = YO * 0.8 - 0.2 * Yppc; 
        
      if opt.verb, fprintf(', YIC: %0.2f%s%0.2f, YOC: %0.2f%s%0.2f',mean(YI),char(177),std(YI),mean(YO),char(177),std(YO)); end

      VOC = V - N .* repmat( TN/2 - YO ,1,3);    % outer surface 
      VIC = V + N .* repmat( TN/2 - YI ,1,3);    % inner surface

      YIC   = isocolors2(Y,VIC); 
      YOC   = isocolors2(Y,VOC);  
      YppOC = isocolors2(Ypp,VOC);  

      YI = YI .* ( abs(YI - (WMth/2 + Yl4/2))  > abs(YIC - (WMth/2 + Yl4/2)));
      YO = YO .* ( abs((CSFth/2 + Yl4/2) - YO) > abs((CSFth/2 + Yl4/2) - YOC) & YppOC>0);

      % filter correction 
      YO = single(spm_mesh_smooth(M,double(YO), sf )); 
      YI = single(spm_mesh_smooth(M,double(YI), sf ));

% if the point is/was perfect then do not change      
% if the new point is perfect / better than the old then simply use it?
      
      % combine 
      if opt.model == 1 % only intensity 
        TIC  = YI; 
        TOC  = YO; 
      elseif opt.model % combine
        mixing = max(0,-0.2 + 1.0*(j/maxiter2));%0.8;
        TIC  = mean( cat( 4, TIC*(1-mixing) , YI*mixing), 4); 
        TOC  = mean( cat( 4, TOC*(1-mixing) , YO*mixing), 4); 
      end  
    end
    
    
    % estimate first corrected inner and outer thickness 
    % Different levels of smoothing were use to have more effect on neighbors. 
    TOC  = single( spm_mesh_smooth(M,double( TOC  ), sf*2 ));  
    TIC  = single( spm_mesh_smooth(M,double( TIC  ), sf*2 ));
    TC   = TOC + TIC; TCsum = rms(TC(TC>0));


    % correction in specified areas that also include a general
    % smoothness constrains of the cortical thickness
    TNC = TN - TOC/2 - TIC/2; 
    TNC = single(spm_mesh_smooth(M,double(TNC), sf*max(0.5,2 - 3*(j/maxiter2))) ); 
    TN(TC>0) = TNC(TC>0);   
    clear TC TNC flim;

    
    % estimate new inner and outer surfaces
    VOC = V - N .* repmat( TN/2 - TOC ,1,3);    % outer surface 
    VIC = V + N .* repmat( TN/2 - TIC ,1,3);    % inner surface
    clear TOC TIC;
    
%    svf = 0.05; 
%    %VOC = VOC*(1-svf) + svf*smoothsurf(VOC,1); 
%    VIC = VIC*(1-svf) + svf*smoothsurf(VIC,1); 
    
    % update thickness and surface
    TN  = sum( (VIC - VOC).^2 , 2) .^ 0.5;
    SN.vertices = mean(cat(3,VIC,VOC),3); 
    
    
    %% this is just for display and loop settings
    %  SX.vertices = VOC; SX.faces = S.faces; SX.facevertexcdata = TC; cat_surf_render2(SX);
    stopiterth = 0.00005; 
    if opt.debug && ( j==1 || mod(j,5)==1 || abs(TCsum)<0.01 || abs(TCsumo - TCsum)<stopiterth ) 
      TNM = TN>(mean(TN(:)) - 2*std(TN(:))) & TN<(mean(TN(:)) + 2*std(TN(:)));
      if ~opt.verb, fprintf('\n'); end
      cat_io_cprintf('g5',sprintf('    reminding overlap:      %8.4f mm (Tlink: %4.2f%s%4.2f mm) %9.0fs',...
        TCsum,mean(TN(TNM)),char(177),std(TN(TNM)),etime(clock,stime) )); stime = clock;
    end
    if ( TCsum<0.005 || abs(TCsumo - TCsum)<stopiterth) && abs( mean(TN(NM)) - TNold )<0.001, break; end
    TCsumo = TCsum; TNold = mean(TN(NM)); 
  end
  
  % export cortical surfaces
  if opt.write, cat_surf_saveICO(SN,TN,Pcs,sprintf('post_collcorr_%0.0fk',round( size(S.faces,1)/1000 / 10) * 10 ),0); else fprintf('\n'); end
  
  
  %% flip back
  if flipped, SN.faces = [SN.faces(:,1) SN.faces(:,3) SN.faces(:,2)]; SN.mati(7) = - SN.mati(7); end
  

end
function cat_surf_show_orthview(Psurf,Pm)
  fg = spm_figure('GetWin','Graphics');
  %fg = spm_figure('Create','SurfaceOverlay');%,Psurf);
  spm_figure('clear')

  id = 1;
  global st

  [pp,ff,ee] = spm_fileparts(Pm);
  hhm = spm_orthviews('Image',spm_vol(Pm));
  spm_orthviews('Caption',hhm,{'m*.nii (Intensity Normalized T1)'},'FontWeight','Bold');
  if ff(1)=='m', spm_orthviews('window',hhm,[0.3 1.03]); caxis([0.3,1.03]); end
  spm_orthviews('AddContext'); % need the context menu for mesh handling

  ov_mesh = 1; 
  for ix=1:numel(Psurf) 
    
    if ov_mesh && exist(Psurf{ix},'file')
      try
        spm_ov_mesh('display',id,Psurf{ix});
      catch
        fprintf('Please update to a newer version of spm12 for using this contour overlay\n');
        ov_mesh = 0;
        continue;
      end
    end
  end

  %% change line style
  if ov_mesh
    styles = {'b-','g-','r-','c-','m-','y-','w-','b.-','g.-','r.-','c.-','m.-','y.-','w.-'}; % need more if meshes were added
    names  = {'central';'pial';'white';'';'';''};
    hM = findobj(st.vols{1}.ax{1}.cm,'Label','Mesh');
    UD = get(hM,'UserData');
    UD.width = [repmat(0.75,1,numel(UD.width) - numel(Psurf))  repmat(0.5,1,numel(Psurf))]; 
    UD.style = styles; %(1:numel(Psurf)); % need more if meshes were added
    set(hM,'UserData',UD);
    spm_ov_mesh('redraw',id);

    % TPM legend
    cc = axes('Position',[0.55 0.4 0.02 0.01],'Parent',fg);
    text(cc,0,1,[spm_str_manip(pp,'t') '/' ff ':']); 
    axis(cc,'off')
    for ix=1:numel(Psurf) 
      cc = axes('Position',[0.55 0.4 - 0.02*ix 0.02 0.01],'Parent',fg);
      plot(cc,[0 1],[1 1],styles{ix}); 
      text(cc,1.2,1,names{ix}); 
      axis(cc,'off')
    end
  end
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
  def.debug  = 0;     % debugging output vs. memory optimization 
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
function V = isocolors2(R,V,opt)
% ______________________________________________________________________
% calculates an interpolated value of a vertex in R  
% We have to calculate everything with double, thus larger images will 
% cause memory issues.
% ______________________________________________________________________
  
  if isempty(V), return; end
  if ndims(R)~=3,  error('MATLAB:isocolor2:dimsR','Only 2 or 3 dimensional input of R.'); end
  if ~exist('opt','var'), opt=struct(); end
  
  def.interp = 'linear';
  opt = cat_io_checkinopt(opt,def);
  
  if  isa(R,'double'), R = single(R); end
  if ~isa(V,'double'), V = double(V); VD=0; else VD=1; end
  
  nV   = size(V,1);
  ndim = size(V,2);
  
  switch opt.interp
    case 'nearest'
      V = max(1,min(round(V),repmat(size(R),nV,1))); 
      V = R(sub2ind(size(R),V(:,2),V(:,1),V(:,3)));
    case 'linear'
      nb  = repmat(shiftdim(double([0 0 0;0 0 1;0 1 0;0 1 1;1 0 0;1 0 1;1 1 0;1 1 1]'),-1),nV,1);  
      enb = repmat(shiftdim((ones(8,1,'double')*[size(R,2),size(R,1),size(R,3)])',-1),nV,1);  

      % calculate the weight of a neigbor (volume of the other corner) and
      w8b = reshape(repmat(V,1,2^ndim),[nV,ndim,2^ndim]); clear V;
      % if the streamline is near the boundary of the image you could be out of range if you add 1 
      n8b = min(floor(w8b) + nb,enb); clear enb
      n8b = max(n8b,1);
      w8b = flipdim(prod(abs(n8b - w8b),2),3);        

      % multiply this with the intensity value of R
      V = sum(R(sub2ind(size(R),n8b(:,2,:),n8b(:,1,:),n8b(:,3,:))) .* w8b,3);
  end  
  if ~VD, V = single(V); end
end
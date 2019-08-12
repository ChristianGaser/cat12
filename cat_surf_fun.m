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
    Pthick = cat_io_strrep(Praw,{'central','.gii'},{'thickness',''});
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
    Pthick = strrep(Praw,'central','thickness');
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
   %Psubthick    = strrep(strrep(Psubcentral,'central','thickness'),'.gii','');               
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
  Sn = patchnormals(S); 
  
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
    
    Sn = patchnormals(S); 
    Sn = Sn ./ repmat(sum(Sn.^2,2).^0.5,1,3); % normalize
    
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

function N = patchnormals(FV) 
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

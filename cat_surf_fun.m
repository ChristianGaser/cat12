function varargout = cat_surf_fun(action,S,varargin)
% Function to collect small surface functions.
% 
% varargout = cat_surf_fun(action,S)
%
%   A   = cat_surf_fun('area',S);     % estimate surface area (faces)
%   V   = cat_surf_F2V(S,F);          % map facedata to vertices
%   HS  = cat_surf_fun('hull',S);     % estimate (optimized) hull surface
%   V   = cat_surf_surf2vol(S,varargin{1}); % render surface in a volume
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
% ______________________________________________________________________
% $Id: cat_surf_createCS.m 937 2016-05-12 08:56:23Z gaser $ 

  switch action
    case 'distance'
      varargout{1} = cat_surf_dist(S);
    case 'area'
      varargout{1} = cat_surf_area(S);
    case 'hull'
      if nargout==1, varargout{1} = cat_surf_hull(S); end
      if nargout==2, [varargout{1},varargout{2}] = cat_surf_hull(S); end
    case 'surf2vol'
      if nargin>2
        [varargout{1},varargout{2},varargout{3}] = cat_surf_surf2vol(S,varargin);
      else
        [varargout{1},varargout{2},varargout{3}] = cat_surf_surf2vol(S);
      end
  end
    
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
  
  % numberical (to small point diffences) and mapping broblems (crossing of streamlines)
  % -> correction because this is theoretical not possible (laplace field theorie)
  AF(AF==0) = eps; % to small values
  AF = abs(AF);    % streamline-crossing
    
  AV = cat_surf_F2V(S,AF);
end
function data = cat_surf_F2V(S,odata)
%% mapping of facedata to vertices

  data   = zeros(size(S.vertices,1),1);
  [v,f]  = sort(S.faces(:)); 
  [f,fj] = ind2sub(size(S.faces),f);  %#ok<ASGLU>
  far = odata(f);
  for i=1:numel(S.faces), data(v(i)) = data(v(i)) + far(i); end

  data = data / size(S.vertices,2); % Schwerpunkt... besser Voronoi, aber wie bei ner Oberfl?che im Raum???
end
function [SH,V] = cat_surf_hull(S)
%% hull creation

  % render surface points
  V = false( round(max(S.vertices,[],1) - min(S.vertices))+10 );     
  I = sub2ind(size(V),round(S.vertices(:,1) - min(S.vertices(:,1)) + 5),...
                      round(S.vertices(:,2) - min(S.vertices(:,2)) + 5),...
                      round(S.vertices(:,3) - min(S.vertices(:,3)) + 5));
  V(I) = 1; clear I; 
  
  % 
  V  = cat_vol_morph(V,'lc',8);  % closeing 
  V  = cat_vol_smooth3X(V,2);    % smoothing
  SH = isosurface(V,0.4);        % create hull 
  V  = V>0.4;
  
  SH.vertices = [SH.vertices(:,2) SH.vertices(:,1) SH.vertices(:,3)]; % matlab flip
  SH.vertices = SH.vertices + repmat(min(S.vertices),size(SH.vertices,1),1) - 5;
end
function [V,vmat,vmati] = cat_surf_surf2vol(S,opt)
%% render inner surface area 
%  Render the volume V with V==1 within the surface. 
%  Use type=1 to render also the surface area with 0.5.
%  The transformation imat to create 
%  SH.vertices = [SH.vertices(:,2) SH.vertices(:,1) SH.vertices(:,3)];     % matlab flip
%  SH.vertices = SH.vertices + imat;

  if ~exist('opt','var'), opt = struct(); end
  def.debug  = 1;
  def.type   = 0; 
  def.refine = 0.8;
  def.bdist  = 5; 
  def.res    = 1; % not yet ...
  
  opt = cat_io_checkinopt(opt,def);
  
  % save a temporary version of S and refine it
  Praw = [tempname '.gii'];
  save(gifti(S),Praw);
  
  cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f',Praw,Praw,opt.refine); 
  [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.debug);

  S = gifti(Praw);
  delete(Praw);
  
  %% render surface points
  V    = false( round(max(S.vertices,[],1) - min(S.vertices))+10 );     
  vmat = -[min(S.vertices(:,1)) min(S.vertices(:,2)) min(S.vertices(:,3))] + opt.bdist; 
  I    = sub2ind(size(V),...
        max(1,min(size(V,1),round(S.vertices(:,1) + vmat(1)))),...
        max(1,min(size(V,2),round(S.vertices(:,2) + vmat(2)))),...
        max(1,min(size(V,3),round(S.vertices(:,3) + vmat(3)))));
  V(I) = 1; 
  
  V    = cat_vol_morph(V,'lc',1);  % closeing 
  V(I) = 0;
  V    = cat_vol_morph(V,'lab');
  if opt.type==1
    Vd = cat_vol_morph(V,'d',1); 
    V  = single(V);
    V(intersect(I,find(Vd>0))) = 0.5;
    V  = cat_vol_smooth3X(V,0.6);    % smoothing
  end
  
  vmati = repmat(min(S.vertices),size(S.vertices,1),1) - 5; 
  %%
  %SH = isosurface(V,0.6);        % create hull 
  %SH.vertices = [SH.vertices(:,2) SH.vertices(:,1) SH.vertices(:,3)]; % matlab flip
  %SH.vertices = SH.vertices + repmat(min(S.vertices),size(SH.vertices,1),1) - 5;
end

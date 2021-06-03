function varargout = cat_vol_gbdist(I,M,bth,side,interp,BP,IVR,IMinf)
%cat_vol_gbdist. Voxel-based distance estimation.
% ______________________________________________________________________
% Calculates the distance of all points p with I(p)>lth and I(p)<hth to 
% an isosurface generated on I with the threshold bth. 
%
%  [DD,ID] = cat_vol_gbdist(I,M[,bth,side,interp])
%     I      .. intenity image for borders
%     M      .. mask for calculation
%     bth    .. border threshold
%     side   .. side with positive distance
%     interp .. use intern interpolation to # voxel
%
% See also cat_vbdist, cat_vol_eidist.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id: cat_surf_epivolsurf.m 1842 2021-06-01 14:41:58Z gaser $ 


  def.bth     = 0.5; 
  def.side    = 1; 
  def.interp  = 1; 
  def.BP      = []; 
 
  if nargin>2 && isstruct(bth)
  else
  % check input
    if ~exist('bth','var'), bth=eps; end
    if ~exist('side','var'), side=1; else side=sign(side); end 
    if side==0, error('side must be a positive or negativ number'); end
    if ~exist('interp','var'), interp=0; end
    if ~exist('BP','var'), BP=[]; end
    if ~isempty(BP) && size(BP,2)~=3, error('only 3d'); end 
  %  if ~exist('IVR','var'); IVR=1; end
    if ~exist('IMinf','var'); IMinf=0; end
  end
  if isa(I,'double'), I=single(I); end
  if ~exist('M','var'), M = I<bth; end
  I(isnan(I))=0; 
  sizeI = single(size(I));
  signI = I<bth; 
  
% find voxels   
  VRi=find(M); clear M;                           % index of the voxels in R 
  [VR(:,2),VR(:,1),VR(:,3)]=ind2sub(sizeI,VRi);   % xyz   of the voxels in R

% surface reconstruction, interpolationen and downsampling if there are to many points  
  %I=smooth3(I,'gaussian',3,0.3*IVR); 
  VB=isosurface(I,bth);
  if isempty(VB.vertices), varargout{1}=zeros(size(I),'single'); varargout{2}=zeros(size(I),'int32'); return; end % no surface...
  if     size(VB.vertices,1)<=150000, VBI=cat_surf_meshinterp(VB,interp,'dist',1); 
  elseif size(VB.vertices,1)<=450000, VBI=cat_surf_meshinterp(VB,interp,'dist',sqrt(2)); 
  end
  
  faceth=600000;
  if exist('VBI','var')
    if ~isempty(VBI) && size(VBI.faces,1)>faceth, VBI=reducepatch(VBI,faceth); end;
    VB=VB.vertices; VBI=VBI.vertices;
    VB=unique([VB;VBI],'rows'); clear VBI; 
  else
    if ~isempty(VB) && size(VB.faces,1)>faceth, VB=reducepatch(VB,faceth); end;
    VB=VB.vertices;
  end
  clear I;
  
  % add other boundary points
  if ~isempty(BP), VB=unique([VB;BP],'rows'); end
  
  try %#ok<TRYNC> error 
    [x,xi]=intersect(VR,VB,'rows'); 
    VR(xi,:)=[]; VRi(xi)=[]; 
    clear x xi;
  end
  
  % create delaunayn triangularisation - error if failed
  save('gbdisttmp.mat','VR','VRi','signI'); clear VR Ri signI;
  T=delaunayn(VB); if ~exist('T','var'), error('Delaunayn triangulation failed.'); end
  load('gbdisttmp.mat','VR'); 
  
  try
    [VID,VDD] = dsearchn(VB,T,VR);
  catch  %#ok<CTCH>
    npoints     = size(VB,1);
    maxpartsize = 1000;                                   
    parts       = max(1,ceil(npoints / maxpartsize)); 
    partsize    = floor(npoints/(parts)); 
    VDD=zeros(size(VR,1),1,'single'); 
    if nargout==2, VID=zeros(size(VR,1),1,'single'); end
    for p = 1:parts
      [l,h] = getrange(p,partsize,parts,npoints);
      [VIDt,VDDt] = dsearchn(VB,T,VR(l:h,:));  
      VDD(l:h) = single(VDDt); 
      if nargout==2, VID(l:h) = single(VIDt); end
    end
  end
  clear VR T;
  load('gbdisttmp.mat','VRi','signI'); delete('gbdisttmp.mat');
  
  if IMinf
    varargout{1}=inf(sizeI,'single');
  else
    varargout{1}=zeros(sizeI,'single');  
  end;
  varargout{1}(VRi) = VDD;
  varargout{1} = side * varargout{1} .* (1 - 2*single(signI)); % setup side 
  varargout{1}(isnan(varargout{1}))=0;                       % remove nans
  
  if nargout==2
    varargout{2} = reshape(1:prod(sizeI),sizeI); 
    varargout{2}(VRi) = sub2ind(sizeI,max(1,min(sizeI(1),round(VB(VID,2)))),max(1,min(sizeI(2),round(VB(VID,1)))),max(1,min(sizeI(3),round(VB(VID,3))))); 
  end
end
function S=cat_surf_meshinterp(S,interp,method,distth)  
  if ~exist('interp','var'), interp = 1; else interp=single(interp); end
  if interp==0, return, end
  if ~exist('method','var'), method = 'linear'; end

  if ~isfield(S,'vertices')        || size(S.vertices,1)==0,        warning('Meshinterp:NoVertices','S has no vertices'); return; end
  if ~isfield(S,'faces')           || size(S.faces,1)==0,           warning('Meshinterp:NoFaces','S has no faces');       return; end
  if ~isfield(S,'facevertexcdata') || size(S.facevertexcdata,1)==0; else C=S.facevertexcdata; end
  if exist('C','var'); CT=(size(C,1)==size(S.vertices,1))+1; else CT=0; end
  
  V=S.vertices; F=single(S.faces); clear S; 
  
  for i=1:interp
    nV=single(size(V,1)); nF=single(size(F,1)); 
    
    NF=(1:nF)';
    
    switch method
      case 'linear'
       
        % addition vertices (middle of the edge)
        V1 = V(F(:,1),:) + 0.5*diff(cat(3,V(F(:,1),:),V(F(:,2),:)),1,3);
        V2 = V(F(:,2),:) + 0.5*diff(cat(3,V(F(:,2),:),V(F(:,3),:)),1,3);
        V3 = V(F(:,3),:) + 0.5*diff(cat(3,V(F(:,3),:),V(F(:,1),:)),1,3);

        % new faces which replace the old one
        F1 = [F(:,1),  nV + 2*nF + NF, nV +        NF];
        F2 = [F(:,2),  nV +        NF, nV +   nF + NF];
        F3 = [F(:,3),  nV +   nF + NF, nV + 2*nF + NF];
        F4 = [nV + NF, nV + 2*nF + NF, nV +   nF + NF];

        % colors
        if     CT==2, C=[C;nanmean(C(F(:,1),:),C(F(:,2),:));nanmean(C(F(:,2),:),C(F(:,3),:));nanmean(C(F(:,3),:),C(F(:,1),:))]; %#ok<AGROW>
        elseif CT==1, C=repmat(C,4,1);
        end
        
        V = [V;V1;V2;V3];  clear V1 V2 V3;    %#ok<AGROW>
        F = [F1;F2;F3;F4]; clear F1 F2 F3 F4; 

        % remove double vertices
        if CT==0, [V,F]   = reduce_points(V,F);
        else      [V,F,C] = reduce_points(V,F,C);
        end
        
      case 'dist'
        if ~exist('distth','var'), distth=sqrt(2); end
        
        % addition vertices (middle of the edge)
        E1 = diff(cat(3,V(F(:,1),:),V(F(:,2),:)),1,3);
        E2 = diff(cat(3,V(F(:,2),:),V(F(:,3),:)),1,3);
        E3 = diff(cat(3,V(F(:,3),:),V(F(:,1),:)),1,3);
        
        V1 = V(F(:,1),:) + repmat((sum(E1.^2,2).^0.5)>=distth,1,3) .* (0.5*E1);
        V2 = V(F(:,2),:) + repmat((sum(E2.^2,2).^0.5)>=distth,1,3) .* (0.5*E2);
        V3 = V(F(:,3),:) + repmat((sum(E3.^2,2).^0.5)>=distth,1,3) .* (0.5*E3);

        % new faces which replace the old one
        F1 = [F(:,1),  nV + 2*nF + NF, nV +        NF];
        F2 = [F(:,2),  nV +        NF, nV +   nF + NF];
        F3 = [F(:,3),  nV +   nF + NF, nV + 2*nF + NF];
        F4 = [nV + NF, nV + 2*nF + NF, nV +   nF + NF];
        
        % colors
        if     CT==2, C=[C;nanmean(C(F(:,1),:),C(F(:,2),:));nanmean(C(F(:,2),:),C(F(:,3),:));nanmean(C(F(:,3),:),C(F(:,1),:))]; %#ok<AGROW>
        elseif CT==1, C=repmat(C,4,1);
        end
        
        V = [V;V1;V2;V3];  clear V1 V2 V3;    %#ok<AGROW>
        F = [F1;F2;F3;F4]; clear F1 F2 F3 F4; 

        
        % remove double vertices
        if CT==0, [V,F]   = reduce_points(V,F);
        else      [V,F,C] = reduce_points(V,F,C);
        end

        % remove degnerated faces
        F((F(:,1)==F(:,2)) | (F(:,1)==F(:,3)) | (F(:,2)==F(:,3)),:)=[];
        
      otherwise 
        error('ERROR: Unknown method "%s"',method); 
    end
  end
  S.vertices = V; S.faces = double(F); if exist('C','var'), S.facevertexcdata = C; end
end
function [V,F,C]=reduce_points(V,F,C)
  try
    [V,~,j]  = unique(V, 'rows'); 
  catch %#ok<CTCH>
    V=single(V);
    [V,~,j]  = unique(V, 'rows'); 
  end
  j(end+1) = nan;
  F(isnan(F)) = length(j);
  if size(F,1)==1, F = j(F)'; if exist('C','var'), C=j(C)'; end    
  else             F = j(F);  if exist('C','var'), C=j(C);  end    
  end
end

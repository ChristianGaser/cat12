function [Yth1,S]=vbm_surf_createCS(V,Ym,Ya,YMF,opt)
% ______________________________________________________________________
% Surface creation and thickness estimation.
%
% [Yth1,S]=vbm_surf_createCS(V,Ym,Ya,YMF,opt)
%
% Yth1 = thickness map
% S    = structure with surfaces, like the left hemishere, that contain
%        vertices, faces, GM thickness (th1), and the transformation to
%        map for nifti space (vmat) and back (vmati).
% V    = spm_vol-structure 
% Ym   = the (local) intensity, noise, and bias corrected T1 image
% Ya   = the atlas map with the ROIs for left and right hemispheres
%        (this is generated with vbm_vol_partvol)
% YMF  = a logical map with the area that has to be filed
%        (this is generated with vbm_vol_partvol)
%   
% opt.surf       = {'left','right'[,'cerebellum','brain']} - side
%    .reduceCS   = 100000 - number of faces
%    .isosmooth  = 0.3    - initial surface smoothing
%    .FSoutput   = 1      - do not remove FreeSurfer surface and 
%                           thickness files
%
% Options set by cg_vbm_defaults.m
%    .interpV    = 0.5    - mm-resolution for estimation
% 
% Here we used the intensity normalized image Ym, rather that the Yp0
% image, because it has more information about sulci that we need 
% especialy for asymetrical sulci.
% Furthermore, all non-cortical regions and blood vessels were removed 
% (for left and right surface). Blood vessels (with high contrast) can 
% lead to strong error in the topology correction. Higher resolution 
% also helps to reduce artifacts.
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
% ______________________________________________________________________
% $Id$ 
 
  if ~exist('opt','var'), opt=struct(); end
  vx_vol = sqrt(sum(V.mat(1:3,1:3).^2));
  
  def.surf      = {'left','right'}; % {'left','right','cerebellum','brain'}
  def.reduceCS  = 100000;  
  def.isosmooth = 0.3;
  def.FSoutput  = 1; 
  opt           = vbm_io_updateStruct(opt,def);
  opt.interpV   = cg_vbm_get_defaults('extopts.pbtres');
  opt.interpV   = max(0.5,min([min(vx_vol),opt.interpV,1]));
  opt.CATDir    = fullfile(spm('dir'),'toolbox','vbm12','CAT'); 
  opt.CATpath   = sprintf('PATH=$PATH:%s',opt.CATDir);
  [pp,ff]       = spm_fileparts(V.fname);
  
  % get both sides in the atlas map
  NS = @(Ys,s) Ys==s | Ys==s+1; 
  
  % filling
  Ymf  = max(Ym,min(1,YMF)); 
  Ymfs = vbm_vol_smooth3X(Ymf,1); 
  Ytmp = vbm_vol_morph(YMF,'d',3) & Ymfs>2.3/3;
  Ymf(Ytmp) = max(min(Ym(Ytmp),1),Ymfs(Ytmp)); clear Ytmp Ymfs YMF Ym; 
  Ymf = Ymf*3;
    
  % removing blood vessels, and other regions
  Yth1=zeros(size(Ymf),'single'); 
  
  for si=1:numel(opt.surf)
    % surface filenames
    Pcsr = fullfile(pp,sprintf('CSr.%s.%s',opt.surf{si},ff)); % raw
    Pcss = fullfile(pp,sprintf('CSs.%s.%s',opt.surf{si},ff)); % sphere
    Pcsf = fullfile(pp,sprintf('CSf.%s.%s',opt.surf{si},ff)); % fiducial
    Pcst = fullfile(pp,sprintf('CSt.%s.%s',opt.surf{si},ff)); % thickness

    % reduce for object area
    switch opt.surf{si}
      case {'L','left'},       Ymfs = max(1,Ymf .* ~(NS(Ya,3) | NS(Ya,7) | NS(Ya,11) | NS(Ya,13)) .* (mod(Ya,2)==1)); 
      case {'R','right'},      Ymfs = max(1,Ymf .* ~(NS(Ya,3) | NS(Ya,7) | NS(Ya,11) | NS(Ya,13)) .* (mod(Ya,2)==0));      
      case {'C','cerebellum'}, Ymfs = max(1,Ymf .* NS(Ya,3));
      case {'B','brain'},      Ymfs = max(1,Ymf);
    end 
    
    
    %% thickness estimation
    [Ymfs,BB]   = vbm_vol_resize(Ymfs,'reduceBrain',vx_vol,2,Ymf>1);    % removing of background
    [Ymfs,resI] = vbm_vol_resize(Ymfs,'interp',V,opt.interpV);          % interpolate volume

    % pbt calculation
    [Yth1i,Yppi] = vbm_vol_pbt(Ymfs,struct('resV',opt.interpV)); clear Ymfs;       
    Yth1i(Yth1i>10)=0; Yppi(isnan(Yppi))=0; 

    % removing non brain objects
    Ymm = Yppi>=0.5 & ~vbm_vol_morph(Yppi>=0.5,'labopen',1); 
    Yppi(Ymm) = 0.5-eps; Yth1i(Ymm) = 0;
    
    % smoothing only in save areas...
    if opt.isosmooth>0
      Yppsi  = Yppi+0;  spm_smooth(Yppi ,Yppsi ,opt.isosmooth/opt.interpV); M=Yppi<0.8; Yppi(M) =Yppsi(M);  clear Yppsi  M;
      Yth1si = Yth1i+0; spm_smooth(Yth1i,Yth1si,opt.isosmooth/opt.interpV); M=Yppi>0.2; Yth1i(M)=Yth1si(M); clear Yth1si M;
    end      
    
    Yth1t = vbm_vol_resize(Yth1i,'deinterp',resI);                    % back to original resolution
    Yth1t = vbm_vol_resize(Yth1t,'dereduceBrain',BB);                 % adding of background
    Yth1  = max(Yth1,Yth1t);                                          % save on main image
    clear Yth1t;
    
    %% surface coordinate transformations
    vmatBB = spm_imatrix(resI.hdrN.mat); vmatBBV = spm_imatrix(V.mat);
    vmatBB(1:3) = vmatBB(1:3) + vmatBBV(7:9).*BB.BB([1,3,5]);
    vmatBB = spm_matrix(vmatBB);

    vmat  = vmatBB(1:3,:)*[0 1 0 0;1 0 0 0;0 0 1 0; 0 0 0 1];
    vmati = inv([vmat; 0 0 0 1]);vmati(4,:)=[];    

    % surface generation 
    CS  = isosurface(Yppi,0.5); CSO=CS;

    % only one major object
    vbm_io_FreeSurfer('write_surf',Pcsr,CS); 
    [SS,SR] = system(sprintf('%s;CAT_SeparatePolygon ''%s'' ''%s'' -1;',opt.CATpath,Pcsr,Pcsr)); RS(SS,SR)
    CS = vbm_io_FreeSurfer('read_surf',Pcsr); 

    % correct the number of vertices depending on the number of major objects
    if opt.reduceCS>0, 
      switch opt.surf{si}
        case {'B','brain'}, CS = reducepatch(CS,opt.reduceCS*2); 
        otherwise,          CS = reducepatch(CS,opt.reduceCS);
      end
    end

    % only one major object
    CS.vertices = (vmat*[CS.vertices' ; ones(1,size(CS.vertices,1))])'; 
    vbm_io_FreeSurfer('write_surf',Pcsr,CS); 

   % write Yppi file with 1 mm resolution for the final deformation
    Yppt = vbm_vol_resize(Yppi,'deinterp',resI);                         % back to original resolution
    Yppt = vbm_vol_resize(Yppt,'dereduceBrain',BB);                     % adding of background
    Vpp  = vbm_io_writenii(V,Yppt,'pp','percentage position map','uint8',[0,1/255],[1 0 0 0],0,[]);
    clear Yppt;

    % spherical surface mapping
    sphsmooth=max(1,round(size(CS.vertices,1)/100000));
    [SS,SR] = system(sprintf('%s;CAT_SeparatePolygon ''%s'' ''%s'' -1;CAT_Surf2Sphere ''%s'' ''%s'' 5 %d',...
      opt.CATpath,Pcsr,Pcsr,Pcsr,Pcss,sphsmooth)); RS(SS,SR)

    % correction and refinement skript only if the spherical projection works
    % otherwise we have to use the original surface what takes longer
    CSX = vbm_io_FreeSurfer('read_surf',Pcss);
    if isnan(CSX.vertices(1)) 
      CS=CSO; sphsmooth=max(2,round(size(CS.vertices,1)/100000));
      CS.vertices = (vmat*[CS.vertices' ; ones(1,size(CS.vertices,1))])'; 
      vbm_io_FreeSurfer('write_surf',Pcsr,CS); 
      [SS,SR] = system(sprintf('%s;CAT_SeparatePolygon ''%s'' ''%s'' -1;CAT_Surf2Sphere ''%s'' ''%s'' 5 %d',...
        opt.CATpath,Pcsr,Pcsr,Pcsr,Pcss,sphsmooth)); RS(SS,SR)
    end

    % topology correction and surface optimation 
    clear CSX;
    [SS,SR] = system(sprintf('%s;CAT_FixTopology -n 81920 -refine_length 1.5 "%s" "%s" "%s"',...
      opt.CATpath,Pcsr,Pcss,Pcsf)); RS(SS,SR);
    [SS,SR] = system(sprintf('%s;CAT_DeformSurf_ui "%s" "%s" "%s" 128 0.05',...
      opt.CATpath,Vpp.fname,Pcsf,Pcsf)); RS(SS,SR);

    % read final surface and map thickness data
    CS = vbm_io_FreeSurfer('read_surf',Pcsf); 
    CS.vertices = (vmati*[CS.vertices' ; ones(1,size(CS.vertices,1))])'; 
    CS.facevertexcdata = isocolors2(Yth1i,CS.vertices); 
    if opt.FSoutput, vbm_io_FreeSurfer('write_surf_data',Pcst,CS.facevertexcdata); end
    
    % visualize a side
    % csp=patch(CS); view(3), camlight, lighting phong, axis equal off; set(csp,'facecolor','interp','edgecolor','none')

    % create output structure
    S.(opt.surf{si}).vertices = CS.vertices;
    S.(opt.surf{si}).faces    = CS.faces;
    S.(opt.surf{si}).th1      = CS.facevertexcdata;
    S.(opt.surf{si}).vmat     = vmat;
    S.(opt.surf{si}).vmati    = vmati;
    clear Yth1i Yppi;

    % we have to delete the origal faces, because they have a different number of vertices after
    % CAT_FixTopology because of the refinement!
    % and we have to delete maybe the freesurfer surface
    delete(Pcsr);  
    delete(Pcss);
    delete(Vpp.fname);
    if ~opt.FSoutput, delete(Pcsf); end 
  end
end
 
function V=isocolors2(R,V,opt)
% ______________________________________________________________________
% calculates a linear interpolatet value of a vertex in R  
% We have to calculate everything with double, else bigger images will 
% lead to problems.
% ______________________________________________________________________
  
  if isempty(V), return; end
  if ~ismatrix(V), error('MATLAB:isocolor2:dimsV','Vertices dimension error. Check V.'); end
  if ndims(R)~=3, error('MATLAB:isocolor2:dimsR','Only 2 or 3 dimensional input of R.'); end
  if ~exist('opt','var'), opt=struct(); end
  
  def.interp = 'linear';
  opt = checkinopt(opt,def);
  
  if  isa(R,'double'), R = single(R); end
  if ~isa(V,'double'), V = double(V); VD=0; else VD=1; end
  
  nV   = size(V,1);
  ndim = size(V,2);
  
  switch opt.interp
    case 'nearest'
      V = max(1,min(round(V),repmat(ndim,nV,1))); 
      V = R(sub2ind(size(R),V(:,2),V(:,1),V(:,3)));
    case 'linear'
      nb  = repmat(shiftdim(double([0 0 0;0 0 1;0 1 0;0 1 1;1 0 0;1 0 1;1 1 0;1 1 1]'),-1),nV,1);  
      enb = repmat(shiftdim((ones(8,1,'double')*[size(R,2),size(R,1),size(R,3)])',-1),nV,1);  

      % calculate the weight of a neigbor (volume of the other corner) and
      w8b = reshape(repmat(V,1,2^ndim),[nV,ndim,2^ndim]); clear V;
      % if the streamline ist near the boundery of the image you could be out of range if you add 1 
      n8b = min(floor(w8b) + nb,enb); clear enb
      n8b = max(n8b,1);
      w8b = flipdim(prod(abs(n8b - w8b),2),3);        

      % multiply this with the intensity-value of R
      V = sum(R(sub2ind(size(R),n8b(:,2,:),n8b(:,1,:),n8b(:,3,:))) .* w8b,3);
  end  
  if ~VD, V = single(V); end
end                   
  

function RS(status,result)
  if status==1 || ~isempty(strfind(result,'ERROR')), error('VBM:system_error',result); end
end
%{
function CS=SeparatePolygon(CS)
  fac = spm_mesh_label(CS.faces); faci=hist(fac(:),1:max(fac(:))); [tmp,faci]=max(faci); 
  V1 = CS.faces(fac==faci,:); V1=unique(V1(:)); % points of the first object
  Vn = CS.faces(fac~=faci,:); Vn=unique(Vn(:)); % points of the other objects
  [Vb,iV1,iVn] = intersect(V1,Vn); Vn(iVn)=[];  % exclusive points of the other objects
  CS.faces(fac~=faci,:)=[];                     % remove faces of the other objects
  CS.vertices(Vn,:)=[];                         % remove vertices of the other objects   
  for i=numel(Vn):-1:1, CS.faces(CS.faces>Vn(i))=CS.faces(CS.faces>Vn(i))-1; end
end
%}
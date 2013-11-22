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
% opt.surf       = {'lh','rh'[,'cerebellum','brain']} - side
%    .reduceCS   = 100000 - number of faces
%    .isosmooth  = 0.3    - initial surface smoothing
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
  
  opt.debug     = 1;
  def.surf      = {'lh','rh'}; % {'lh','rh','cerebellum','brain'}
  def.reduceCS  = 100000;  
  def.isosmooth = 0.3;
  opt           = vbm_io_updateStruct(opt,def);
  opt.interpV   = cg_vbm_get_defaults('extopts.pbtres');
  opt.interpV   = max(0.5,min([min(vx_vol),opt.interpV,1]));
  opt.fsavgDir  = fullfile(spm('dir'),'toolbox','vbm12','fsaverage'); 
  opt.CATDir    = fullfile(spm('dir'),'toolbox','vbm12','CAT');   

  % add system dependent extension to CAT folder
  if ispc
    opt.CATDir = [opt.CATDir '.w32'];
  elseif ismac
    opt.CATDir = [opt.CATDir '.maci64'];
  elseif isunix
    opt.CATDir = [opt.CATDir '.glx'];
  end  

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
  Yth1 = zeros(size(Ymf),'single'); 
  
  for si=1:numel(opt.surf)
    % surface filenames
    Praw       = fullfile(pp,sprintf('%s.central.nofix.%s',opt.surf{si},ff));    % raw
    Psphere0   = fullfile(pp,sprintf('%s.sphere.nofix.%s',opt.surf{si},ff));     % sphere.nofix
    Pcentral   = fullfile(pp,sprintf('%s.central.%s',opt.surf{si},ff));          % fiducial
    Pthick     = fullfile(pp,sprintf('%s.thickness.%s',opt.surf{si},ff));        % thickness
    Psphere    = fullfile(pp,sprintf('%s.sphere.%s',opt.surf{si},ff));           % sphere
    Pspherereg = fullfile(pp,sprintf('%s.sphere.reg.%s',opt.surf{si},ff));       % sphere.reg
    Pfsavg     = fullfile(opt.fsavgDir,sprintf('%s.smoothwm',opt.surf{si}));     % fsaverage smoothwm
    Pfsavgsph  = fullfile(opt.fsavgDir,sprintf('%s.sphere',opt.surf{si}));       % fsaverage sphere

    % reduce for object area
    switch opt.surf{si}
      case {'L','lh'},         Ymfs = max(1,Ymf .* ~(NS(Ya,3) | NS(Ya,7) | NS(Ya,11) | NS(Ya,13)) .* (mod(Ya,2)==1)); 
      case {'R','rh'},         Ymfs = max(1,Ymf .* ~(NS(Ya,3) | NS(Ya,7) | NS(Ya,11) | NS(Ya,13)) .* (mod(Ya,2)==0));      
      case {'C','cerebellum'}, Ymfs = max(1,Ymf .* NS(Ya,3));
      case {'B','brain'},      Ymfs = max(1,Ymf);
    end 
    
    
    %% thickness estimation
    [Ymfs,BB]   = vbm_vol_resize(Ymfs,'reduceBrain',vx_vol,2,Ymf>1);    % removing of background
    [Ymfs,resI] = vbm_vol_resize(Ymfs,'interp',V,opt.interpV);          % interpolate volume

    % pbt calculation
    fprintf('\n%s:\nThickness estimation\n',opt.surf{si});
    [Yth1i,Yppi] = vbm_vol_pbt(Ymfs,struct('resV',opt.interpV)); clear Ymfs;       
    Yth1i(Yth1i>10)=0; Yppi(isnan(Yppi))=0; 

    % removing non brain objects
    Ymm = Yppi>=0.5 & ~vbm_vol_morph(Yppi>=0.5,'labopen',1); 
    Yppi(Ymm) = 0.5-eps; Yth1i(Ymm) = 0;
    
    % smoothing only in safe areas...
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

    % use only one major object
    vbm_io_FreeSurfer('write_surf',Praw,CS);
    cmd = sprintf('CAT_SeparatePolygon "%s" "%s" -1',Praw,Praw);
    [ST, RS] = system(fullfile(opt.CATDir,cmd)); check_system_output(ST,RS,opt.debug);
    CS = vbm_io_FreeSurfer('read_surf',Praw); 

    % correct the number of vertices depending on the number of major objects
    if opt.reduceCS>0, 
      switch opt.surf{si}
        case {'B','brain'}, CS = reducepatch(CS,opt.reduceCS*2); 
        otherwise,          CS = reducepatch(CS,opt.reduceCS);
      end
    end

    % only one major object
    CS.vertices = (vmat*[CS.vertices' ; ones(1,size(CS.vertices,1))])'; 
    vbm_io_FreeSurfer('write_surf',Praw,CS); 

    % write Yppi file with 1 mm resolution for the final deformation
    Yppt = vbm_vol_resize(Yppi,'deinterp',resI);                         % back to original resolution
    Yppt = vbm_vol_resize(Yppt,'dereduceBrain',BB);                      % adding of background
    Vpp  = vbm_io_writenii(V,Yppt,'pp','percentage position map','uint8',[0,1/255],[1 0 0 0],0,[]);
    clear Yppt;

    % spherical surface mapping
    disp('Initial spherical mapping');
    cmd = sprintf('CAT_SeparatePolygon "%s" "%s" -1',Praw,Praw);
    [ST, RS] = system(fullfile(opt.CATDir,cmd)); check_system_output(ST,RS,opt.debug);
    cmd = sprintf('CAT_Surf2Sphere "%s" "%s" 5',Praw,Psphere0);
    [ST, RS] = system(fullfile(opt.CATDir,cmd)); check_system_output(ST,RS,opt.debug);

    % use correction and refinement script only if the spherical projection works
    % otherwise we have to use the original surface which takes much longer
    CSX = vbm_io_FreeSurfer('read_surf',Psphere0);
    if isnan(CSX.vertices(1)) 
      disp('Spherical mapping with original high-resoluted surface');
      CS = CSO;
      CS.vertices = (vmat*[CS.vertices' ; ones(1,size(CS.vertices,1))])'; 
      vbm_io_FreeSurfer('write_surf',Praw,CS); 
      cmd = sprintf('CAT_SeparatePolygon "%s" "%s" -1',Praw,Praw);
      [ST, RS] = system(fullfile(opt.CATDir,cmd));  check_system_output(ST,RS,opt.debug);
      cmd = sprintf('CAT_Surf2Sphere "%s" "%s" 5',Praw,Psphere0);
      [ST, RS] = system(fullfile(opt.CATDir,cmd)); check_system_output(ST,RS,opt.debug);
    end

    % topology correction and surface refinement 
    clear CSX;
    disp('Topology correction');
    cmd = sprintf('CAT_FixTopology -n 81920 -refine_length 1.5 "%s" "%s" "%s"',Praw,Psphere0,Pcentral);
    [ST, RS] = system(fullfile(opt.CATDir,cmd)); check_system_output(ST,RS,opt.debug);
    disp('Surface refinement');
    th = 128;
    cmd = sprintf('CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .5 avg -0.15 0.15 .1 .1 15 0 "%g" "%g" n 0 0 0 250 0.01 0.0',Vpp.fname,Pcentral,Pcentral,th,th);
    [ST, RS] = system(fullfile(opt.CATDir,cmd)); check_system_output(ST,RS,opt.debug);
    cmd = sprintf('CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .5 avg -0.05 0.05 .1 .1 15 0 "%g" "%g" n 0 0 0 250 0.01 0.0',Vpp.fname,Pcentral,Pcentral,th,th);
    [ST, RS] = system(fullfile(opt.CATDir,cmd)); check_system_output(ST,RS,opt.debug);


    % spherical surface mapping of corrected surface
    disp('Final spherical mapping');
    cmd = sprintf('CAT_Surf2Sphere "%s" "%s" 5',Pcentral,Psphere);
    [ST, RS] = system(fullfile(opt.CATDir,cmd)); 
    check_system_output(ST,RS,opt.debug);

    % spherical registration to fsaverage
    disp('Spherical registration');
    cmd = sprintf('CAT_WarpSurf -i "%s" -is "%s" -t "%s" -ts "%s" -ws "%s"',Pcentral,Psphere,Pfsavg,Pfsavgsph,Pspherereg);
    [ST, RS] = system(fullfile(opt.CATDir,cmd)); check_system_output(ST,RS,opt.debug);

    % read final surface and map thickness data
    disp('Project thickness map to surface');
    CS = vbm_io_FreeSurfer('read_surf',Pcentral); 
    CS.vertices = (vmati*[CS.vertices' ; ones(1,size(CS.vertices,1))])'; 
    CS.facevertexcdata = isocolors2(Yth1i,CS.vertices); 
    vbm_io_FreeSurfer('write_surf_data',Pthick,CS.facevertexcdata);
    
    % visualize a side
    % csp=patch(CS); view(3), camlight, lighting phong, axis equal off; set(csp,'facecolor','interp','edgecolor','none')

    % create output structure
    S.(opt.surf{si}).vertices = CS.vertices;
    S.(opt.surf{si}).faces    = CS.faces;
    S.(opt.surf{si}).th1      = CS.facevertexcdata;
    S.(opt.surf{si}).vmat     = vmat;
    S.(opt.surf{si}).vmati    = vmati;
    clear Yth1i Yppi;

    % we have to delete the original faces, because they have a different number of vertices after
    % CAT_FixTopology!
    delete(Praw);  
    delete(Psphere0);
    delete(Vpp.fname);
  end
end
 
function V = isocolors2(R,V,opt)
% ______________________________________________________________________
% calculates a linear interpolated value of a vertex in R  
% We have to calculate everything with double, thus larger images will 
% cause memory issues.
% ______________________________________________________________________
  
  if isempty(V), return; end
  if ~ismatrix(V), error('MATLAB:isocolor2:dimsV','Vertices dimension error. Check V.'); end
  if ndims(R)~=3,  error('MATLAB:isocolor2:dimsR','Only 2 or 3 dimensional input of R.'); end
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
  

function check_system_output(status,result,debugON)
  if status==1 || ~isempty(strfind(result,'ERROR')), error('VBM:system_error',result); end
  if nargin > 2
    if debugON, disp(result); end
  end
end

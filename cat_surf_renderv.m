function varargout = cat_surf_renderv(S,facevertexcdata,opt)
% cat_surf_renderv. Simple push rendering for brain surfaces.
% Surface rendering with openGL is currently not working on surfaces. 
% This function was designe as simple replacement of the matlab surface 
% render in the CAT report functions. It renders only the points and do not
% fill faces. It is just a simple solution and still in development. 
%
%  h   = cat_surf_renderv(S[,facevertexcdata,opt])
%  img = cat_surf_renderv(S[,facevertexcdata,opt])
%  
%  S                .. file or surface mesh with vertices, faces, and 
%                      facevertexdata/cdata field
%  facevertexcdata  .. file or surface texture (vertexdata only)
%  opt              .. render parameter
%   .mat            .. orientation matrix (e.g. for rigid registration)
%   .view           .. render view = ['l'|'r'|'t'|'d'|'f'|b'] (default 'l')
%   .interp         .. resolution paramter (default 1.4) 
%                      (higer = more pixel but slower)
%   .bd             .. image boundary (default 0.5)
%   .h              .. figure/axis handle
% 
% See cat_main_reportfig for example use.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$ 


  dbs = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

  if ~exist('opt','var'), opt = struct(); end

  % default parameter
  def.mat        = eye(4); 
  def.view       = 'l'; 
  def.clim       = [0 6]; 
  def.interp     = 1.4;
  def.h          = []; 
  def.bd         = 0.5;
  def.CATrefMesh = 0;     % not working at all - texture resmaple problem
  def.optimize   = 0;     % not really working - in development

  opt = cat_io_checkinopt(opt,def); 

  opt.bd = opt.bd * 2^def.interp;

  
  % load/set data
  if ischar(S) || iscell(S)
    Ps = S; 
    S  = gifti(char(Ps));
  end
  if ischar(facevertexcdata) || iscell(facevertexcdata)
    %[pp,ff,ee] = cat_
    %facevertexcdata = gifti
  end
  
  if ~isempty(facevertexcdata)
    S.facevertexcdata = facevertexcdata;
  elseif isfield(S,'facevertexcdata')
    S.facevertexcdata = S.facevertexcdata; 
    facevertexcdata   = S.facevertexcdata;
  elseif isfield(S,'cdata')
    S.facevertexcdata = S.cdata; 
    facevertexcdata   = S.cdata;
  end
  
  
  % rotation settings
  if ischar( opt.view )
    switch opt.view
      case 't', view = [  0   0  90]; 
      case 'd', view = [  0   0 -90]; 
      case 'r', view = [  0  90   0];
      case 'l', view = [180 -90   0];
      case 'f', view = [  0   0   0]; 
      case 'b', view = [  0   0   0]; 
    end
  else
    view = opt.view; 
  end
  
  
  % main
  if debug, tic; end
  Sx = S;

  % rotation  
  Sx.vertices = [Sx.vertices, ones(size(Sx.vertices,1),1)] * ...
    (spm_matrix([0 0 0 deg2rad( view ) 1 1 1 0 0 0]) ); 
  Sx.vertices = Sx.vertices(:,1:3); 

  
  if opt.CATrefMesh
    % surface refinement with CAT_RefineMesh .. not yet working
    Sxo = Sx; 
    Pcentral = sprintf('%s.gii',tempname);        
    save(gifti(struct('faces',Sx.faces,'vertices',Sx.vertices)),Pcentral);    
    cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f 0',Pcentral,Pcentral,1); 
    cat_system(cmd,0);
    Sx = gifti(Pcentral);
    delete(Pcentral); 
    Sx.facevertexcdata = cat_surf_fun('surf2surf',Sxo,Sx,Sxo.facevertexcdata); % ... not working 
  end
  
  
  %  Optimize:
  %  The general idea is to remove faces that are not vissible and those
  %  will not need further refinement but the problem is that they can 
  %  effect the normals ...
  %  in the best case this could half runtime
  if 1 %opt.optimize
    cat_mesh_smooth = @(M,C,s) single( spm_mesh_smooth(M,double(C),s) );
  
    for ii = 1:round(opt.interp)
      
      if 0
      % remove backface faces .. the test is slower than the improvement
        M    = spm_mesh_smooth(Sx);  
        Nv   = spm_mesh_normals(struct('vertices',Sx.vertices,'faces',Sx.faces));
        NvA  = abs( cat_surf_fun('angle',Nv,repmat([0 0 1],size(Nv,1),1))); clear Nv; 
        NvA  = cat_mesh_smooth(M,NvA,10); clear M; 
        back = sum( NvA(Sx.faces) < 45 , 2); % the angle depend on surface orientatation definition
        Sx.faces(back > 2,:) = [];
        clear NvA; 
      end
      
      
      if 0
      % the idea was to remove small faces does not need further refinement
      % but is also not working
        [x,A] = cat_surf_fun('area',Sx); clear x; 
        Sx.faces(A < 0.25 / 2^opt.interp,:) = [];
      end

      if 1
      % removel faces those points are not visible 
        M    = spm_mesh_smooth(Sx);  
        Vb   = round(Sx.vertices * 2^(opt.interp)); %max(2,min(8,8 / 2^ii)));
        Vb(:,3) = -Vb(:,3);
        Vb   = Vb - repmat( min( Vb ) , size(Vb,1) , 1 ) + opt.bd*2; 
        imgz = inf( round(max( Vb(:,1:2) ) + opt.bd*4),'single'); 
        for vi = 1:size(Vb,1)
          i = round(Vb(vi,1));
          j = round(Vb(vi,2));
          if Vb(vi,3) < imgz(i,j) 
            imgz(i,j) = Vb(vi,3); 
          end
        end
        hidden = zeros(size(Sx.facevertexcdata),'single');
        for vi = 1:size(Vb,1)
          i = round(Vb(vi,1));
          j = round(Vb(vi,2));
          if Vb(vi,3) > imgz(i,j) + 1
            hidden(vi) = 1; 
          end
        end
        hidden = cat_mesh_smooth(M,hidden,10); 
        Sx.faces( sum( hidden(Sx.faces),2) > 0.999,:) = [];
      end

      % mesh interpolation 
      Sx = cat_surf_fun('meshinterp',Sx,1);
    end
  else
    Sx = cat_surf_fun('meshinterp',Sx,opt.interp);
  end
 
  
  % estimate normals
  % we here use a direct camlight 
  N   = abs( spm_mesh_normals(struct('vertices',Sx.vertices,'faces',Sx.faces)));
  NA  = cat_surf_fun('angle',N,repmat([0 0 1],size(N,1),1));
  NA  = abs( NA - 90 ); 
  
  
  % increase sampling by adding displaced points to avoid holes
  % ... limit this to voxel that are orientated to the cam
  Sx.vertices(:,3) = -Sx.vertices(:,3);
  Sx.vertices = Sx.vertices - repmat( min( Sx.vertices ) , size(Sx.vertices,1) , 1 ) + opt.bd; 
  Sx.vertices = Sx.vertices * 2^opt.interp; 
  method = 21; 
  if method == 1
    Sx.vertices = [Sx.vertices; Sx.vertices - 0.25.*N; Sx.vertices - 0.5.*N; Sx.vertices - 1.0.*N]; 
    cdatax = repmat(Sx.facevertexcdata,4,1); 
    NA     = repmat(NA,4,1); 
  elseif method == 2
    Sx.vertices = [floor(Sx.vertices(:,1)*1)/1 , floor(Sx.vertices(:,2)*1)/1 , Sx.vertices(:,3) ; 
                    floor(Sx.vertices(:,1)*1)/1 ,  ceil(Sx.vertices(:,2)*1)/1 , Sx.vertices(:,3) ; 
                     ceil(Sx.vertices(:,1)*1)/1 , floor(Sx.vertices(:,2)*1)/1 , Sx.vertices(:,3) ; 
                     ceil(Sx.vertices(:,1)*1)/1 ,  ceil(Sx.vertices(:,2)*1)/1 , Sx.vertices(:,3) ; 
                     Sx.vertices - 0.25.*N; Sx.vertices - 0.5.*N; Sx.vertices - 1.0.*N
                     ]; 
    cdatax = repmat(Sx.facevertexcdata,7,1); 
    NA     = repmat(NA,7,1); 
  elseif method == 21
    Sx.vertices = [floor(Sx.vertices(:,1)*1)/1 , floor(Sx.vertices(:,2)*1)/1 , Sx.vertices(:,3) ; 
                   floor(Sx.vertices(:,1)*1)/1 ,  ceil(Sx.vertices(:,2)*1)/1 , Sx.vertices(:,3) ; 
                    ceil(Sx.vertices(:,1)*1)/1 , floor(Sx.vertices(:,2)*1)/1 , Sx.vertices(:,3) ; 
                    ceil(Sx.vertices(:,1)*1)/1 ,  ceil(Sx.vertices(:,2)*1)/1 , Sx.vertices(:,3) ; 
                    ]; 
    cdatax = repmat(Sx.facevertexcdata,4,1); 
    NA     = repmat(NA,4,1); 
  elseif method == 3
    Sx.vertices = [Sx.vertices; Sx.vertices - 0.5.*N]; 
    cdatax       = repmat(Sx.facevertexcdata,2,1); 
    NA           = repmat(NA,2,1); 
  end
    

  % render image
  imgc = zeros( round(max( Sx.vertices(:,1:2) ) + opt.bd*2 ),'single'); % texture map
  imgz =   inf( round(max( Sx.vertices(:,1:2) ) + opt.bd*2 ),'single'); % z-depth map (visibility)
  imgn = zeros( round(max( Sx.vertices(:,1:2) ) + opt.bd*2 ),'single'); % normal map (lightning)
  for vi = 1:size(Sx.vertices,1)
    i = round( Sx.vertices(vi,1) ); 
    j = round( Sx.vertices(vi,2) );
    if Sx.vertices(vi,3) < imgz(i,j) 
      imgc(i,j) = cdatax(vi); 
      imgn(i,j) = NA(vi); 
      imgz(i,j) = Sx.vertices(vi,3); 
    end
  end
  
  bgm = repmat( isinf( imgz ) , 1,1,3);
  bgm = cat_vol_morph( bgm ,'o');
  
  % normalize range to 255 for RGB convertation
  imgtn   = min(255,max(0,round( ( imgc - opt.clim(1)  ) / opt.clim(2)                * 255 ))); clear imgc;
  imgzn   = min(255,max(0,round( ( imgz - min(imgz(:)) ) / max( imgz(imgz(:)<inf ) )  * 255 ))); %clear imgz;
  imgnn   = min(255,max(0,round( ( imgn - min(imgn(:)) ) / max( imgn(imgn(:)<inf ) )  * 255 ))); clear imgn;

  % median
  if 0 % slow and inoptimal
    imgtn = median2d(imgtn,10);
    imgzn = median2d(imgzn,10);
    imgnn = median2d(imgnn,10);
  end
  
  % smooth
  if 1
    imgtn = smooth2d(imgtn,0.3);
    imgnn = smooth2d(imgnn,0.3);
    imgzn = smooth2d(imgzn,0.3);
    clear imgs; 
  end
  
  % convert images to RGB
  imgRGB  = ind2rgb( imgtn , jet(256) );  clear imgtn; 
  imgzRGB = ind2rgb( imgzn , gray(256) ); clear imgzn; 
  imgnRGB = ind2rgb( imgnn , gray(256) ); clear imgnn; 

  % set background
  bgc  = get(opt.h,'color'); 
  bg   = cat(3,bgc(1) * ones( size(imgz)), bgc(2) * ones( size(imgz)), bgc(3) * ones( size(imgz)) ); 
  img  = min(255,max(0,imgRGB .* (0.2 + max(0,imgnRGB*1.05 - 0.05).^0.5 * 1.0) - (imgzRGB/8))) ; 
  img  = img.*(1-bgm) + bg.*bgm; 
  
  if debug 
    %figure(4534); imagesc( img ); axis equal off; % subplot(2,1,1); image( imgRGB ); subplot(2,1,2); 
    toc
  end

  % set output
  if ~isempty( opt.h )
    image( opt.h , img ); 
    axis(opt.h,'equal','off'); 
     
    varargout{1}.h     = opt.h; 
    varargout{1}.cdata = facevertexcdata; 
  else
    varargout{1} = img; 
  end
end
function imgc = median2d(imgc,th)
  imgs = repmat(single(imgc),1,1,3); 
  imgs = cat_vol_median3(imgs); %,true(size(imgs)),true(size(imgs)),th);
  imgc = imgs(:,:,2);  
end
function img = smooth2d(img,th)
  imgs = smooth3(repmat(img,1,1,3)); 
  img  = img*(1-th) + th*imgs(:,:,2); 
  img  = cat_vol_ctype(img);
end

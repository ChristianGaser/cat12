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

 %#ok<*AGROW>

  if ~exist('opt','var'), opt=struct(); end
  vx_vol = sqrt(sum(V.mat(1:3,1:3).^2));
  
  def.debug     = cg_vbm_get_defaults('extopts.debug');
  def.surf      = {'lh','rh'}; % {'lh','rh','cerebellum','brain'}
  def.reduceCS  = 100000;  
  def.isosmooth = 0.3;
  opt           = vbm_io_updateStruct(def,opt);
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
    opt.CATDir = [opt.CATDir '.glnx86'];
  end  

  [pp,ff]   = spm_fileparts(V.fname);

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
    try

      % surface filenames
      Praw       = fullfile(pp,sprintf('%s.central.nofix.%s',opt.surf{si},ff));    % raw
      Psphere0   = fullfile(pp,sprintf('%s.sphere.nofix.%s',opt.surf{si},ff));     % sphere.nofix
      Pcentral   = fullfile(pp,sprintf('%s.central.%s',opt.surf{si},ff));          % fiducial
      Pthick     = fullfile(pp,sprintf('%s.thickness.%s',opt.surf{si},ff));        % thickness
      Pdefects   = fullfile(pp,sprintf('%s.defects.%s',opt.surf{si},ff));          % defects
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
      if si==1, fprintf('\n'); end
      fprintf('%s:\n',opt.surf{si});
      stime = vbm_io_cmd('  Thickness estimation');

      [Ymfs,BB]   = vbm_vol_resize(Ymfs,'reduceBrain',vx_vol,2,Ymfs>1);   % removing of background
      [Ymfs,resI] = vbm_vol_resize(Ymfs,'interp',V,opt.interpV);          % interpolate volume

      % pbt calculation
      [Yth1i,Yppi] = vbm_vol_pbt(Ymfs,struct('resV',opt.interpV)); clear Ymfs;       
      Yth1i(Yth1i>10)=0; Yppi(isnan(Yppi))=0; 

      Yth1t = vbm_vol_resize(Yth1i,'deinterp',resI);                      % back to original resolution
      Yth1t = vbm_vol_resize(Yth1t,'dereduceBrain',BB);                   % adding of background
      Yth1  = max(Yth1,Yth1t);                                            % save on main image
      clear Yth1t;
      fprintf('%4.0fs\n',etime(clock,stime)); 


      %% Write Ypp for final deformation
      %  Write Yppi file with 1 mm resolution for the final deformation, 
      %  because CAT_DeformSurf_ui can not handle higher resolutions and
      %  will create cauliflower surfaces. 
      %  
      %  #################################################################
      %  Here we still have a problem. The image with 0 mm looks good, and
      %  its coordinatats are correct. Check Reg shows pefect matching, but
      %  something still went wrong... but if the error is not in the picture
      %  this meas that the error is in the surface and that vmat and vmati
      %  are not correct!
      %  #################################################################
      usePPmap = 0;
      if usePPmap
        Yppt = vbm_vol_resize(Yppi,'deinterp',resI);                        % back to original resolution
        Yppt = vbm_vol_resize(Yppt,'dereduceBrain',BB);                     % adding of background
        Vpp  = vbm_io_writenii(V,Yppt,'pp','percentage position map','uint8',[0,1/255],[1 0 0 0],0,[]);
        clear Yppt;

        Vpp1 = Vpp; 
        Vpp1.fname    = fullfile(pp,['pp1' ff ee]);
        vmat2         = spm_imatrix(Vpp1.mat);
        Vpp1.dim(1:3) = round(Vpp1.dim .* abs(vmat2(7:9)));
        vmat2(7:9)    = sign(vmat2(7:9)).*[1 1 1];
        Vpp1.mat      = spm_matrix(vmat2);

        Vpp1 = spm_create_vol(Vpp1); 
        for x3 = 1:Vpp1.dim(3),
          M    = inv(spm_matrix([0 0 -x3 0 0 0 1 1 1])*inv(Vpp1.mat)*Vpp.mat);
          v    = spm_slice_vol(Vpp,M,Vpp1.dim(1:2),1);       
          Vpp1 = spm_write_plane(Vpp1,v,x3);
        end;
        clear M v x3; 
      end


      %% surface coordinate transformations
      stime = vbm_io_cmd('  Create initial surface');
      vmatBBV = spm_imatrix(V.mat);

      vmat  = V.mat(1:3,:)*[0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];
      vmati = inv([vmat; 0 0 0 1]); vmati(4,:)=[];    

      % surface generation and coordinate adaption for the orinal image
      % this works well...
      CS  = isosurface(Yppi,0.5); %clear Yppi;
      CS.vertices = CS.vertices .* repmat(abs(opt.interpV ./ vmatBBV([8,7,9])),size(CS.vertices,1),1);
      CS.vertices = CS.vertices + repmat( BB.BB([3,1,5]) - 1,size(CS.vertices,1),1); 
      %CS.vertices(:,1:2) = CS.vertices(:,2:-1:1);
      CSO = CS; % save old surface, for later correction of reduction error

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
      CS.vertices = (vmat*[CS.vertices' ; ones(1,size(CS.vertices,1))])'; 
      vbm_io_FreeSurfer('write_surf',Praw,CS); 
      fprintf('%4.0fs\n',etime(clock,stime)); 


      %% spherical surface mapping 1 of the uncorrected surface for topology correction
      stime = vbm_io_cmd('  Initial spherical mapping');
      cmd = sprintf('CAT_SeparatePolygon "%s" "%s" -1',Praw,Praw);
      [ST, RS] = system(fullfile(opt.CATDir,cmd)); check_system_output(ST,RS,opt.debug);
      cmd = sprintf('CAT_Surf2Sphere "%s" "%s" 5',Praw,Psphere0);
      [ST, RS] = system(fullfile(opt.CATDir,cmd)); check_system_output(ST,RS,opt.debug);
      fprintf('%4.0fs\n',etime(clock,stime)); 

      %% mark defects and save as gifti
      cmd = sprintf('CAT_MarkDefects "%s" "%s" "%s"',Praw,Psphere0,Pdefects);
      [ST, RS] = system(fullfile(opt.CATDir,cmd)); check_system_output(ST,RS,opt.debug);
      cmd = sprintf('CAT_AddValuesToSurf "%s" "%s" "%s"',Praw,Pdefects,[Pdefects '.gii']);
      [ST, RS] = system(fullfile(opt.CATDir,cmd)); check_system_output(ST,RS,opt.debug);
      
      %% topology correction and surface refinement 
      str = '  Topology correction and surface refinement'; fprintf('%s:%s',str,repmat(' ',1,67-length(str))); stime = clock;
      cmd = sprintf('CAT_FixTopology -n 81920 -refine_length 1.5 "%s" "%s" "%s"',Praw,Psphere0,Pcentral);
      [ST, RS] = system(fullfile(opt.CATDir,cmd)); check_system_output(ST,RS,opt.debug);
      if usePPmap
        % surface refinement by surface deformation based on the PP map
        th = 128;
        cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .5 ' ...
                       'avg -0.15 0.15 .1 .1 15 0 "%g" "%g" n 0 0 0 250 0.01 0.0'], ...
                       Vpp1.fname,Pcentral,Pcentral,th,th);
        [ST, RS] = system(fullfile(opt.CATDir,cmd)); check_system_output(ST,RS,opt.debug);
        cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .5 ' ...
                       'avg -0.05 0.05 .1 .1 15 0 "%g" "%g" n 0 0 0 250 0.01 0.0'], ...
                       Vpp1.fname,Pcentral,Pcentral,th,th);
        [ST, RS] = system(fullfile(opt.CATDir,cmd)); check_system_output(ST,RS,opt.debug);
      else
        % surface refinement by simple smoothing
        cmd = sprintf('CAT_BlurSurfHK "%s" "%s" 2',Pcentral,Pcentral);
        [ST, RS] = system(fullfile(opt.CATDir,cmd)); check_system_output(ST,RS,opt.debug);
      end
      fprintf('%4.0fs\n',etime(clock,stime)); 


      %% spherical surface mapping 2 of corrected surface
      stime = vbm_io_cmd('  Spherical mapping');
      cmd = sprintf('CAT_Surf2Sphere "%s" "%s" 5',Pcentral,Psphere);
      [ST, RS] = system(fullfile(opt.CATDir,cmd)); check_system_output(ST,RS,opt.debug);
      fprintf('%4.0fs\n',etime(clock,stime)); 

      % spherical registration to fsaverage
      stime = vbm_io_cmd('  Spherical registration');
      cmd = sprintf('CAT_WarpSurf -i "%s" -is "%s" -t "%s" -ts "%s" -ws "%s"',Pcentral,Psphere,Pfsavg,Pfsavgsph,Pspherereg);
      [ST, RS] = system(fullfile(opt.CATDir,cmd)); check_system_output(ST,RS,opt.debug);
      fprintf('%4.0fs\n',etime(clock,stime)); 

      % read final surface and map thickness data
      CS  = vbm_io_FreeSurfer('read_surf',Pcentral); 
      CS.vertices = (vmati*[CS.vertices' ; ones(1,size(CS.vertices,1))])'; 
      CS.facevertexcdata = isocolors2(Yth1,CS.vertices); 
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
      delete(Pdefects);  
      delete(Psphere0);
      if usePPmap
        delete(Vpp.fname);
        delete(Vpp1.fname);
      end
    catch err
      switch err.identifier
        case 'VBM:surf_createCS:system_error'
          vbm_io_cprintf('error','ERR\n')
        case 'VBM:surf_createCS:segmenationfault'
          vbm_io_cprintf('error','ERR\n')
        otherwise
          rethrow(err);
      end
    end
  end  
end
 
function V = isocolors2(R,V,opt)
% ______________________________________________________________________
% calculates a linear interpolated value of a vertex in R  
% We have to calculate everything with double, thus larger images will 
% cause memory issues.
% ______________________________________________________________________
  
  if isempty(V), return; end
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
  if status==1
    if ~isempty(strfind(result,'ERROR'))
      error('VBM:surf_createCS:system_error',result); 
    elseif ~isempty(strfind(result,'Segmentation fault'))
      error('VBM:surf_createCS:segmenationfault',result); 
    end
  end
  if nargin > 2
    if debugON, disp(result); end
  end
end

function [Yth1,S,Psurf] = cat_surf_createCS(V,Ym,Ya,YMF,opt)
% ______________________________________________________________________
% Surface creation and thickness estimation.
%
% [Yth1,S]=cat_surf_createCS(V,Ym,Ya,YMF,opt)
%
% Yth1 = thickness map
% S    = structure with surfaces, like the left hemishere, that contains
%        vertices, faces, GM thickness (th1), and the transformation to
%        map to nifti space (vmat) and back (vmati).
% V    = spm_vol-structure 
% Ym   = the (local) intensity, noise, and bias corrected T1 image
% Ya   = the atlas map with the ROIs for left and right hemispheres
%        (this is generated with cat_vol_partvol)
% YMF  = a logical map with the area that has to be filled
%        (this is generated with cat_vol_partvol)
%   
% opt.surf       = {'lh','rh'[,'cerebellum','brain']} - side
%    .reduceCS   = 100000 - number of faces
%
% Options set by cat_defaults.m
%    .interpV    = 0.5    - mm-resolution for thickness estimation
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

  % set defaults
  vx_vol = sqrt(sum(V.mat(1:3,1:3).^2));
  if ~exist('opt','var'), opt=struct(); end
  def.verb      = 2; 
  def.surf      = {'lh','rh'}; % {'lh','rh','cerebellum','brain'}
  def.interpV   = max(0.25,min([min(vx_vol),opt.interpV,1]));
  def.reduceCS  = 100000;  
  def.usePPmap  = 1; 
  def.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 
  def.CATDir    = fullfile(spm('dir'),'toolbox','cat12','CAT');   
  opt           = cat_io_updateStruct(def,opt);

  Psurf = struct(); 
  
  % add system dependent extension to CAT folder
  if ispc
    opt.CATDir = [opt.CATDir '.w32'];
  elseif ismac
    opt.CATDir = [opt.CATDir '.maci64'];
  elseif isunix
    opt.CATDir = [opt.CATDir '.glnx86'];
  end  

  % correction for 'n' prefix for noise corrected and/or interpolated files
  [pp,ff]   = spm_fileparts(V.fname);

  if cat_get_defaults('extopts.subfolders')
    surffolder = 'surf';
    mrifolder = 'mri';
    pp = spm_str_manip(pp,'h'); % remove 'mri' in pathname that already exists
  else
    surffolder = '';
    mrifolder = '';
  end

  if ff(1)=='n'
    if (exist(fullfile(pp,[ff(2:end) '.nii']), 'file')) || (exist(fullfile(pp,[ff(2:end) '.img']), 'file'))
      ff = ff(2:end);
    end
  end

  % get both sides in the atlas map
  NS = @(Ys,s) Ys==s | Ys==s+1; 
    
  % filling
  Ymf  = max(Ym,min(0.95,YMF)); 
  Ymfs = cat_vol_smooth3X(Ymf,1); 
  Ytmp = cat_vol_morph(YMF,'d',3) & Ymfs>2.3/3;
  Ymf(Ytmp) = max(min(Ym(Ytmp),0),Ymfs(Ytmp)); clear Ytmp Ymfs YMF; 
  Ymf = Ymf*3;
    
  % removing blood vessels, and other regions
  Yth1 = zeros(size(Ymf),'single'); 
  if opt.expertgui > 1
    Ywd  = zeros(size(Ymf),'single'); 
    Ycd  = zeros(size(Ymf),'single'); 
  end
  
  [D,I] = cat_vbdist(single(Ya>0)); Ya = Ya(I); % for sides
  
  for si=1:numel(opt.surf)
   
    % surface filenames
    Praw       = fullfile(pp,surffolder,sprintf('%s.central.nofix.%s.gii',opt.surf{si},ff));    % raw
    Psphere0   = fullfile(pp,surffolder,sprintf('%s.sphere.nofix.%s.gii',opt.surf{si},ff));     % sphere.nofix
    Pcentral   = fullfile(pp,surffolder,sprintf('%s.central.%s.gii',opt.surf{si},ff));          % fiducial
    Pthick     = fullfile(pp,surffolder,sprintf('%s.thickness.%s',opt.surf{si},ff));            % thickness / GM depth
    Pgwo       = fullfile(pp,surffolder,sprintf('%s.depthWMo.%s',opt.surf{si},ff));             % gyrus width / GWM depth / gyral span
    Pgw        = fullfile(pp,surffolder,sprintf('%s.depthGWM.%s',opt.surf{si},ff));             % gyrus width / GWM depth / gyral span
    Pgww       = fullfile(pp,surffolder,sprintf('%s.depthWM.%s',opt.surf{si},ff));              % gyrus witdh of the WM / WM depth
    Pgwwg      = fullfile(pp,surffolder,sprintf('%s.depthWMg.%s',opt.surf{si},ff));             % gyrus witdh of the WM / WM depth
    Psw        = fullfile(pp,surffolder,sprintf('%s.depthCSF.%s',opt.surf{si},ff));             % sulcus width / CSF depth / sulcal span
    Pdefects0  = fullfile(pp,surffolder,sprintf('%s.defects.%s',opt.surf{si},ff));              % defects temporary file
    Pdefects   = fullfile(pp,surffolder,sprintf('%s.defects.%s.gii',opt.surf{si},ff));          % defects
    Psphere    = fullfile(pp,surffolder,sprintf('%s.sphere.%s.gii',opt.surf{si},ff));           % sphere
    Pspherereg = fullfile(pp,surffolder,sprintf('%s.sphere.reg.%s.gii',opt.surf{si},ff));       % sphere.reg
    Pfsavg     = fullfile(opt.fsavgDir,sprintf('%s.central.freesurfer.gii',opt.surf{si}));      % fsaverage central
    Pfsavgsph  = fullfile(opt.fsavgDir,sprintf('%s.sphere.freesurfer.gii',opt.surf{si}));       % fsaverage sphere    

    surffile = {'Praw','Psphere0','Pcentral','Pthick','Pgw','Pgww','Psw',...
      'Pdefects0','Pdefects','Psphere','Pspherereg','Pfsavg','Pfsavgsph'};
    for sfi=1:numel(surffile)
      eval(sprintf('Psurf(si).%s = %s;',surffile{sfi},surffile{sfi})); 
    end
    
    % reduce for object area
    switch opt.surf{si}
      case {'L','lh'},         Ymfs = Ymf .* (Ya>0) .* ~(NS(Ya,3) | NS(Ya,7) | NS(Ya,11) | NS(Ya,13)) .* (mod(Ya,2)==1); Yside = mod(Ya,2)==1;
      case {'R','rh'},         Ymfs = Ymf .* (Ya>0) .* ~(NS(Ya,3) | NS(Ya,7) | NS(Ya,11) | NS(Ya,13)) .* (mod(Ya,2)==0); Yside = mod(Ya,2)==0;      
      case {'C','cerebellum'}, Ymfs = Ymf .* (Ya>0) .* NS(Ya,3); Yside = NS(Ya,3)>0;
      case {'B','brain'},      Ymfs = Ymf .* (Ya>0); Yside = true(size(Ya));
    end 
    
    %% thickness estimation
    if si==1, fprintf('\n'); end
    fprintf('%s:\n',opt.surf{si});
    stime = cat_io_cmd(sprintf('  Thickness estimation (%0.2f mm%s)',opt.interpV,char(179)));
    
    [Ymfs,Yside,BB] = cat_vol_resize({Ymfs,Yside},'reduceBrain',vx_vol,4,Ymfs>1.1); % removing background
    [Ymfs,resI]     = cat_vol_resize(Ymfs,'interp',V,opt.interpV);                  % interpolate volume
    Yside           = cat_vol_resize(Yside,'interp',V,opt.interpV)>0.5;             % interpolate volume
   
    % pbt calculation
    [Yth1i,Yppi] = cat_vol_pbt(max(1,Ymfs),struct('resV',opt.interpV)); % avoid underestimated thickness in gyris
    if ~opt.expertgui, clear Ymfs; end
    Yth1i(Yth1i>10)=0; Yppi(isnan(Yppi))=0;  
    [D,I] = cat_vbdist(Yth1i,Yside); Yth1i = Yth1i(I); clear D I;       % add further values around the cortex
    Yth1t = cat_vol_resize(Yth1i,'deinterp',resI); clear Yth1i;         % back to original resolution
    Yth1t = cat_vol_resize(Yth1t,'dereduceBrain',BB);                   % adding background
    Yth1  = max(Yth1,Yth1t);                                            % save on main image
    clear Yth1t;
    fprintf('%4.0fs\n',etime(clock,stime)); 
    
    %% PBT estimation of the gyrus and sulcus width 
    if opt.expertgui > 1
      %% gyrus width / WM depth
      %  For the WM depth estimation it is better to use the L4 boundary
      %  and correct later for thickness, because the WM is very thin in
      %  gyral regions and will cause bad values. 
      %  On the other side we do not want the whole filled block of the 
      %  Yppi map and so we have to mix both the original WM map and the
      %  Yppi map. 
      %  As far as there is no thickness in pure WM regions there will
      %  be no correction. 
      %
      %    figure, isosurface(smooth3(Yppi),0.5,Yth1i), axis equal off
      stime = cat_io_cmd('  WM depth estimation');
      [Yar,Ymr,BB] = cat_vol_resize({Ya,Ym},'reduceBrain',vx_vol,BB.BB);    % removing background
      Yar   = uint8(cat_vol_resize(Yar,'interp',V,opt.interpV,'nearest'));  % interpolate volume
      Ymr   = cat_vol_resize(Ymr,'interp',V,opt.interpV);                   % interpolate volume
      switch opt.surf{si}
        case {'L','lh'}, 
          Ymr = Ymr .* (Yar>0) .* ~(NS(Yar,3) | NS(Yar,7) | NS(Yar,11) | NS(Yar,13)) .* (mod(Yar,2)==1);
          Ynw = smooth3(cat_vol_morph(NS(Yar,5) | NS(Yar,9) | NS(Yar,15) | NS(Yar,23),'d',2) | ...
                 (cat_vol_morph(Yppi==1,'e',2) & Ymr>1.7/3 & Ymr<2.5/3) & (mod(Yar,2)==1)); 
        case {'R','rh'},
          Ymr = Ymr .* (Yar>0) .* ~(NS(Yar,3) | NS(Yar,7) | NS(Yar,11) | NS(Yar,13)) .* (mod(Yar,2)==0);    
          Ynw = smooth3(cat_vol_morph(NS(Yar,5) | NS(Yar,9) | NS(Yar,15) | NS(Yar,23),'d',2) | ...
                 (cat_vol_morph(Yppi==1,'e',2) & Ymr>1.7/3 & Ymr<2.5/3) & (mod(Yar,2)==0)); 
        case {'C','cerebellum'}, Ymr = Ymr .* (Yar>0) .* NS(Yar,3);
        case {'B','brain'},      Ymr = Ymr .* (Yar>0);
      end 
     % clear Yar; 
      %%
      Yppis = Yppi .* (1-Ynw) + max(0,min(1,Ymr*3-2)) .* Ynw;           % adding real WM map 
      Ywdt  = cat_vbdist(1-Yppis);                                      % estimate distance map to central/WM surface
      Ywdt  = cat_vol_pbtp(max(2,4-Ymfs),Ywdt,inf(size(Ywdt),'single'))*opt.interpV;
      [D,I] = cat_vbdist(single(Ywdt>0),Yside); Ywdt = Ywdt(I); clear D I;    % add further values around the cortex
      %%
      Ywdt  = cat_vol_median3(Ywdt); Ywdt = smooth3(Ywdt);              % smoothing
      Ywdt  = cat_vol_resize(Ywdt,'deinterp',resI); 
      Ywdt  = cat_vol_resize(Ywdt,'dereduceBrain',BB);                  % adding background
      Ywd   = max(Ywd,Ywdt); 
      clear Ywdt;
      
      %% sulcus width / CSF depth
      %  for the CSF depth we cannot use the origal data, because of
      %  sulcal blurring, but we got the PP map at half distance and
      %  correct later for half thickness
      fprintf('%4.0fs\n',etime(clock,stime)); 
      stime = cat_io_cmd('  CSF depth estimation');
      YM    = smooth3(cat_vol_morph(Ymfs>0.5,'o',4))>0.5;               % smooth CSF/background-skull boundary 
      Yppis = min(Ymr,Yppi); Yppis(isnan(Yppis))=0;                     % we want also CSF within the ventricle (for tests)
      Ycdt  = cat_vbdist(Yppis,YM); clear Yppis                         % distance to the cental/CSF-GM boundary
      Ycdt  = cat_vol_pbtp(max(2,Ymfs),Ycdt,inf(size(Ycdt),'single'))*opt.interpV; 
      Ycdt(~YM)=0;
      [D,I] = cat_vbdist(single(Ycdt>0),Yside); Ycdt = Ycdt(I); clear D I;    % add further values around the cortex
      Ycdt  = cat_vol_median3(Ycdt); Ycdt = smooth3(Ycdt);              % smoothing
      Ycdt  = cat_vol_resize(Ycdt,'deinterp',resI); 
      Ycdt  = cat_vol_resize(Ycdt,'dereduceBrain',BB); 
      Ycd   = max(Ycd,Ycdt); 
      clear Ycdt;
      fprintf('%4.0fs\n',etime(clock,stime));
      clear Ymr;
    end
    clear Ymfs;
    
    
    %% Write Ypp for final deformation
    %  Write Yppi file with 1 mm resolution for the final deformation, 
    %  because CAT_DeformSurf_ui cannot handle higher resolutions.
    if opt.usePPmap
      Yppt = cat_vol_resize(Yppi,'deinterp',resI);                        % back to original resolution
      Yppt = cat_vol_resize(Yppt,'dereduceBrain',BB);                     % adding of background
      Vpp  = cat_io_writenii(V,Yppt,'','pp','percentage position map','uint8',[0,1/255],[1 0 0 0]);
      clear Yppt;

      Vpp1 = Vpp; 
      Vpp1.fname    = fullfile(pp,mrifolder,['pp1' ff '.nii']);
      vmat2         = spm_imatrix(Vpp1.mat);
      Vpp1.dim(1:3) = round(Vpp1.dim .* abs(vmat2(7:9)));
      vmat2(7:9)    = sign(vmat2(7:9)).*[1 1 1];
      Vpp1.mat      = spm_matrix(vmat2);

      Vpp1 = spm_create_vol(Vpp1); 
      for x3 = 1:Vpp1.dim(3),
        M    = inv(spm_matrix([0 0 -x3 0 0 0 1 1 1])*inv(Vpp1.mat)*Vpp.mat); %#ok<MINV>
        v    = spm_slice_vol(Vpp,M,Vpp1.dim(1:2),1);       
        Vpp1 = spm_write_plane(Vpp1,v,x3);
      end;
      clear M v x3; 
    end

    %% surface coordinate transformations
    stime = cat_io_cmd('  Create initial surface','g5','',opt.verb); fprintf('\n');
    vmatBBV = spm_imatrix(V.mat);

    vmat  = V.mat(1:3,:)*[0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];
    vmati = inv([vmat; 0 0 0 1]); vmati(4,:) = [];    

    % if we can use the PP map we can start with a surface that is close to WM surface because this might minimize severe
    % topology defects. Otherwise we use a threshold of 0.5 which is the central surface.
    % However, this approach did not really improved topology correction, thus we again use a value of 0.5
    if opt.usePPmap, th_initial = 0.5; else th_initial = 0.5; end
    [tmp,CS.faces,CS.vertices] = cat_vol_genus0(Yppi,th_initial);
    clear tmp Yppi;

    % correction for the boundary box used within the surface creation process 
    CS.vertices = CS.vertices .* repmat(abs(opt.interpV ./ vmatBBV([8,7,9])),size(CS.vertices,1),1);
    CS.vertices = CS.vertices +  repmat( BB.BB([3,1,5]) - 1,size(CS.vertices,1),1);

    fprintf('%s %4.0fs\n',repmat(' ',1,66),etime(clock,stime)); 
    
    % correct the number of vertices depending on the number of major objects
    if opt.reduceCS>0, 
      switch opt.surf{si}
        case {'B','brain'}, CS = reducepatch(CS,opt.reduceCS*2); 
        otherwise,          CS = reducepatch(CS,opt.reduceCS);
      end
      stime = cat_io_cmd(sprintf('  Reduce surface to %d faces:',size(CS.faces,1)),'g5','',opt.verb); 
    end
    
    % transform coordinates
    CS.vertices = (vmat*[CS.vertices' ; ones(1,size(CS.vertices,1))])'; 
    save(gifti(struct('faces',CS.faces,'vertices',CS.vertices)),Praw);

    olddir = pwd;
    cd(opt.CATDir);
    
    % after reducepatch many triangles have very large area which causes isses for resampling
    % RefineMesh addds triangles in those areas
    cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f',Praw,Praw,2); 
    [ST, RS] = system(cmd); cat_check_system_output(ST,RS,opt.debug);

    % remove some unconnected meshes
    cmd = sprintf('CAT_SeparatePolygon "%s" "%s" -1',Praw,Praw); % CAT_SeparatePolygon works here
    [ST, RS] = system(cmd); cat_check_system_output(ST,RS,opt.debug);

    % spherical surface mapping 1 of the uncorrected surface for topology correction
    cmd = sprintf('CAT_Surf2Sphere "%s" "%s" 5',Praw,Psphere0);
    [ST, RS] = system(cmd); cat_check_system_output(ST,RS,opt.debug);

    % mark defects and save as gifti 
    if opt.debug == 2 
      cmd = sprintf('CAT_MarkDefects -binary "%s" "%s" "%s"',Praw,Psphere0,Pdefects0); 
      [ST, RS] = system(cmd); cat_check_system_output(ST,RS,opt.debug);
      cmd = sprintf('CAT_AddValuesToSurf "%s" "%s" "%s"',Praw,Pdefects0,Pdefects);
      [ST, RS] = system(cmd); cat_check_system_output(ST,RS,opt.debug);
    end
   
    %% topology correction and surface refinement 
    stime = cat_io_cmd('  Topology correction and surface refinement','g5','',opt.verb,stime);
    cmd = sprintf('CAT_FixTopology -deform -n 81920 -refine_length 2 "%s" "%s" "%s"',Praw,Psphere0,Pcentral);
    [ST, RS] = system(cmd); cat_check_system_output(ST,RS,opt.debug);
    
    if opt.usePPmap
      % we use thickness values to get from the initial (white matter) surface to the central surface
      % the extent depends on the inital threshold of the surface creation
      extent = th_initial - 0.5;
      if extent ~= 0
        CS = gifti(Pcentral);
        CS.vertices = (vmati*[CS.vertices' ; ones(1,size(CS.vertices,1))])';
        facevertexcdata = isocolors2(Yth1,CS.vertices); 
        cat_io_FreeSurfer('write_surf_data',Pthick,facevertexcdata);

        cmd = sprintf('CAT_Central2Pial "%s" "%s" "%s" "%g"',Pcentral,Pthick,Pcentral,extent);
        [ST, RS] = system(cmd); cat_check_system_output(ST,RS,opt.debug);
      end
      
      % surface refinement by surface deformation based on the PP map
      th = 0.5;
      cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .1 ' ...
                     'avg -0.01 0.01 .1 .1 5 0 "%g" "%g" n 0 0 0 100 0.01 0.0'], ...
                     Vpp1.fname,Pcentral,Pcentral,th,th);
      [ST, RS] = system(cmd); cat_check_system_output(ST,RS,opt.debug);
      
      % need some more refinement because some vertices are distorted after CAT_DeformSurf
      cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f',Pcentral,Pcentral,1.5);
      [ST, RS] = system(cmd); cat_check_system_output(ST,RS,opt.debug);
      
      cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .5 ' ...
                     'avg -0.1 0.1 .1 .1 5 0 "%g" "%g" n 0 0 0 100 0.01 0.0'], ...
                     Vpp1.fname,Pcentral,Pcentral,th,th);
      [ST, RS] = system(cmd); cat_check_system_output(ST,RS,opt.debug);
    else
      % surface refinement by simple smoothing
      cmd = sprintf('CAT_BlurSurfHK "%s" "%s" %0.2f',Pcentral,Pcentral,2);
      [ST, RS] = system(cmd); cat_check_system_output(ST,RS,opt.debug);
    end
    
    %% spherical surface mapping 2 of corrected surface
    stime = cat_io_cmd('  Spherical mapping with areal smoothing','g5','',opt.verb,stime); 
    cmd = sprintf('CAT_Surf2Sphere "%s" "%s" 10',Pcentral,Psphere);
    [ST, RS] = system(cmd); cat_check_system_output(ST,RS,opt.debug);
    
    % spherical registration to fsaverage template
    stime = cat_io_cmd('  Spherical registration','g5','',opt.verb,stime);
    cmd = sprintf('CAT_WarpSurf -type 0 -i "%s" -is "%s" -t "%s" -ts "%s" -ws "%s"',Pcentral,Psphere,Pfsavg,Pfsavgsph,Pspherereg);
    [ST, RS] = system(cmd); cat_check_system_output(ST,RS,opt.debug);
    
    % read final surface and map thickness data
    stime = cat_io_cmd('  Thickness / Depth mapping','g5','',opt.verb,stime);
    if ~opt.usePPmap || ((th_initial - 0.5) == 0)
      CS = gifti(Pcentral);
      CS.vertices = (vmati*[CS.vertices' ; ones(1,size(CS.vertices,1))])';
      facevertexcdata = isocolors2(Yth1,CS.vertices); 
      cat_io_FreeSurfer('write_surf_data',Pthick,facevertexcdata);
    end
    
    % map WM and CSF width data (corrected by thickness)
    if opt.expertgui > 1
      %%
      facevertexcdata2  = isocolors2(Ywd,CS.vertices); 
      facevertexcdata2c = max(eps,facevertexcdata2 - facevertexcdata/2);
      cat_io_FreeSurfer('write_surf_data',Pgwo,facevertexcdata2c); % gyrus width WM only
      facevertexcdata2c = correctWMdepth(CS,facevertexcdata2c,100,0.2);
      cat_io_FreeSurfer('write_surf_data',Pgww,facevertexcdata2c); % gyrus width WM only
      facevertexcdata3c = facevertexcdata2c + facevertexcdata; % );
      cat_io_FreeSurfer('write_surf_data',Pgw,facevertexcdata3c); % gyrus width (WM and GM)
      facevertexcdata4 = estimateWMdepthgradient(CS,facevertexcdata2c);
      cat_io_FreeSurfer('write_surf_data',Pgwwg,facevertexcdata4); % gyrus width WM only > gradient
      % smooth resampled values
      try
        cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s"',Pcentral,Pgwwg,3,Pgwwg);
        [ST, RS] = system(cmd); cat_check_system_output(ST,RS,opt.debug);
      end
      %%
      %clear facevertexcdata2 facevertexcdata2c facevertexcdata3c facevertexcdata4; 
      % just a test ... problem with other species ...
      %norm = sum(Ymf(:)>0.5) / prod(vx_vol) / 1000 / 1400;
      %norm = mean([2 1 1].*diff([min(CS.vertices);max(CS.vertices)])); 
      %norm = mean([2 1 1].*std(CS.vertices)); % maybe the hull surface is better...
 
      facevertexcdata3 = isocolors2(Ycd,CS.vertices); 
      facevertexcdata3 = max(eps,facevertexcdata3 - facevertexcdata/2); 
      cat_io_FreeSurfer('write_surf_data',Psw,facevertexcdata3);
    end
    fprintf('%4.0fs\n',etime(clock,stime)); 
    
    % visualize a side
    % csp=patch(CS); view(3), camlight, lighting phong, axis equal off; set(csp,'facecolor','interp','edgecolor','none')

    % create output structure
    S.(opt.surf{si}).vertices = CS.vertices;
    S.(opt.surf{si}).faces    = CS.faces;
    S.(opt.surf{si}).vmat     = vmat;
    S.(opt.surf{si}).vmati    = vmati;
    S.(opt.surf{si}).th1    = facevertexcdata;
    if opt.expertgui > 1
      S.(opt.surf{si}).th2    = facevertexcdata2;
      S.(opt.surf{si}).th3    = facevertexcdata3;
    end
    clear Yth1i
    
    cd(olddir);

    % we have to delete the original faces, because they have a different number of vertices after
    % CAT_FixTopology!
    delete(Praw);  
    if opt.debug == 2
      delete(Pdefects0);  
    end
    delete(Psphere0);
    if opt.usePPmap
      delete(Vpp.fname);
      delete(Vpp1.fname);
    end
    clear CS
  end  
  
  if opt.debug && opt.verb
    for si=1:numel(Psurf)
      fprintf('Display thickness: %s\n',spm_file(Psurf(si).Pthick,'link','cat_surf_display(''%s'')'));
    end
  end
end

%=======================================================================
function [cdata,i] = correctWMdepth(CS,cdata,iter,lengthfactor)
% ______________________________________________________________________
% Correct deep WM depth values that does not fit to the local thickness 
% of the local gyri.
% 
% lengthfactor should be between 0.2 and 0.4
% ______________________________________________________________________

  if ~exist('lengthfactor','var'), lengthfactor = 1/3; end
  if ~exist('iter','var'), iter = 100; end

  %%
  SV  = CS.vertices;                                                          % Surface Vertices 
  SE  = unique([CS.faces(:,1:2);CS.faces(:,2:3);CS.faces(:,3:-2:1)],'rows');  % Surface Edges
  SEv = single(diff(cat(3,SV(SE(:,1),:),SV(SE(:,2),:)),1,3));                 % Surface Edge Vector
  SEL = sum(SEv.^2,2).^0.5;                                                   % Surface Edge Length  
  clear SEv

  
  %%
  i=0; cdatac = cdata+1; pc = 1; oc = 0; 
  while i<iter && pc~=oc; 
  %%
    pc = sum( abs(cdata - cdatac)>0.05 ); 
    i=i+1; cdatac = cdata;
    
    M  = (cdatac(SE(:,1)) - SEL(SE(:,1))*lengthfactor ) > cdatac(SE(:,2)); 
    cdata(SE(M,1)) = cdatac(SE(M,2)) + SEL(SE(M,1))*lengthfactor; 
    M  = (cdata(SE(:,2)) - SEL(SE(:,2))*lengthfactor ) > cdatac(SE(:,1));
    cdata(SE(M,2)) = cdatac(SE(M,1)) + SEL(SE(M,1))*lengthfactor; 
    oc = sum( abs(cdata - cdatac)>0.05 );
    
    %fprintf('%d - %8.2f - %d\n',i,sum( abs(cdata - cdatac)>0.05 ),pc~=oc)
    
  end
  
end
%=======================================================================
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
  opt = cat_io_checkinopt(opt,def);
  
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
     %=======================================================================
function cdata = estimateWMdepthgradient(CS,cdata)
% ______________________________________________________________________
% Estimates the maximum local gradient of a surface. 
% Major use is the WM depth that grows with increasing sulcal depth. 
% It measures the amount of WM behind the cortex, but more relevant is
% the amout of WM fibers that this reagion will add to the WM depth. 
% The width of the street next to a house gives not the connectivity of
% this house, but the width of the entrance does!
% This measure can be improved by furhter information of sulcal depth.
% ______________________________________________________________________

  %%
  SV  = CS.vertices;                                                          % Surface Vertices 
  SE  = unique([CS.faces(:,1:2);CS.faces(:,2:3);CS.faces(:,3:-2:1)],'rows');  % Surface Edges
  SEv = single(diff(cat(3,SV(SE(:,1),:),SV(SE(:,2),:)),1,3));                 % Surface Edge Vector
  SEL = sum(SEv.^2,2).^0.5;                                                   % Surface Edge Length  
  clear SEv

  
  %%
  cdata_l = inf(size(cdata),'single'); 
  cdata_h = zeros(size(cdata),'single'); 
  for i=1:size(SE,1)
    val = (cdata(SE(i,2)) - cdata(SE(i,1)))*SEL(SE(i,1));
    cdata_l(SE(i,1)) = min([cdata_l(SE(i,1)),val]);
    cdata_h(SE(i,1)) = max([cdata_h(SE(i,2)),val]);
  end
  cdata = cdata_h - cdata_l; 
end
              

function [Yth,S,P,res] = cat_surf_createCS5(V,V0,Ym,Yp0,Ya,YMF,Yb0,opt,job)
% ______________________________________________________________________
% Surface creation and thickness estimation.
%
% [Yth,S,P] = cat_surf_createCS4(V,V0,Ym,Ya,YMF,Ypb0,opt,job)
%
% Yth    .. thickness map
% S      .. structure with surfaces, like the left hemisphere, that contains
%           vertices, faces, GM thickness (th1)
% P  .. name of surface files
% res    .. intermediate and final surface creation information
% V      .. spm_vol-structure of internally interpolated image
% V0     .. spm_vol-structure of original image
% Ym     .. the (local) intensity, noise, and bias corrected T1 image
% Ya     .. the atlas map with the ROIs for left and right hemispheres
%           (this is generated with cat_vol_partvol)
% YMF    .. a logical map with the area that has to be filled
%           (this is generated with cat_vol_partvol)
% Ytemplate .. Shooting template to improve cerebellar surface
%              reconstruction
% Yb    .. modified mask from gcut
%   
% opt.surf       = {'lh','rh'[,'lc','rc']} - side
%
% Options set by cat_defaults.m
%    .interpV    = 0.5    - mm-resolution for thickness estimation
% 
% Here we use the intensity normalized image Ym, rather than the Yp0
% image, because it has more information about sulci that we need 
% especially for asymmetrical sulci.
% Furthermore, all non-cortical regions and blood vessels are removed 
% (for left and right surface). Blood vessels (with high contrast) can 
% lead to strong errors in the topology correction. Higher resolution 
% also helps to reduce artifacts.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$ 

  %#ok<*AGROW,*STREMP,*ASGLU,*SFLD,*STFLD>

  % get both sides in the atlas map
  NS = @(Ys,s) Ys==s | Ys==s+1; 
  if ~exist('opt','var'), opt = struct(); end                 % create variable if not exist
  
  % Turn off gifti data warning in gifti/subsref (line 45)
  %   Warning: A value of class "int32" was indexed with no subscripts specified. 
  %            Currently the result of this operation is the indexed value itself, 
  %            but in a future release, it will be an error. 
  warning('off','MATLAB:subscripting:noSubscriptsSpecified');
  cstime = clock; %#ok<*CLOCK>

  % test-variables that should be (partially) removed later
  skip_registration       = isfield(opt,'surf') && isscalar(opt.surf); % skip spherical registration for quick tests
  create_white_pial       = 2; % uses only the quick WM and Pial surface estimation (0-no,1-yes,2-improve)
  myelinCorrection        = .3; % .25 - sight correction, 1 - maximum correction
  setcut2zero             = 0; % works but result in worse values as could be expected
   
  % surface output and evaluation parameter 
  res   = struct('lh',struct(),'rh',struct()); 
  Yth   = zeros(size(Yp0),'single');  % initialize WM/CSF thickness/width/depth maps
  S     = struct();
  

  % set defaults
  % set debugging variable
  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end
  vx_vol                  = sqrt(sum(V.mat(1:3,1:3).^2));               % further interpolation based on internal resolution 
  def.verb                = cat_get_defaults('extopts.expertgui');      % 0-none, 1-minimal, 2-default, 3-details, 4-debug
  def.surf                = {'lh','rh'};                                % surface reconstruction setting with {'lh','rh','cb'} 
  % There is a new SPM approach spm_mesh_reduce that is maybe more robust. 
  % Higher resolution is at least required for animal preprocessing that is given by cat_main.
  def.LAB                 = cat_get_defaults('extopts.LAB');  % brain regions 
  % RD20250306: Tfs has large issues currently with some corrected defects
  def.useprior            = ''; 
  def.thick_limit         = cat_get_defaults('extopts.thick_limit');       % 6mm upper limit for thickness (similar limit as used in Freesurfer)
  def.thick_measure       = cat_get_defaults('extopts.thick_measure');     % 0-PBT; 1-Tfs (Freesurfer method using mean(TnearIS,TnearOS)) ##########
  def.foldingcorrection   = cat_get_defaults('extopts.foldingcorrection'); % tickness correction that is influence by folding
  def.fsavgDir            = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces'); 
  def.outputpp.native     = 0;  % output of Ypp map for cortical orientation in EEG/MEG 
  def.outputpp.warped     = 0;
  def.outputpp.dartel     = 0;
  def.vdist               = 2; 
  def.wmnoise             = nan; 
  def.eidist              = 0; 

  % options that rely on other options
  opt                     = cat_io_updateStruct(def,opt); clear def; 
  opt.vol                 = any(~cellfun('isempty',strfind(opt.surf,'v')));   % only volume-based thickness estimation  
  opt.surf                = cat_io_strrep(opt.surf,'v','');                   % after definition of the 'vol' varialbe we simplify 'surf'
  opt.interpV             = max(0.1,min([opt.interpV,1.5]));                  % general limitation of the PBT resolution
  opt.reconres            = opt.interpV;
  switch opt.SRP
    case 0 % pbtsimpleC, full resolution
      use_cat_vol_pbtsimple = 0; 
    case 2 % default
      use_cat_vol_pbtsimple = 1;   
      create_white_pial     = 2 - isscalar(opt.surf); % quick for tests
  end


  % apply the modified mask from gcut
  % for non-gcut approaches or inverse weighting Yb0 only contains ones
  Yp0   = Yp0 .* (Yb0>0.5);
  Ym    = Ym  .* (Yb0>0.5);


  Ywmr = cat_vol_resize(Ym .* (Yp0>2.75),'reduceV',1,2,16,'meanm'); 
  Ywmrstd = cat_vol_localstat(Ywmr,Ywmr>0,2,4); 
  opt.wmnoise  = cat_stat_nanmean(Ywmrstd(Ywmrstd>0)); 


  % enlarge atlas definition 
  [~,I] = cat_vbdist(single(Ya>0)); Ya=Ya(I); clear I;  


  % improve PVE (important for SPM25)
  % - correct MRF classification problems, i.e., realign GM voxels with
  %   very high/low raw intensities to WM/CSF (thin WM in BUSS01)
  % - this might include blood vessels but these we have to handle anyway
  % - this might be better after denoising but before segmentation 
  %   is this done in case of AMAP?
  % - in case of SPM 0.75 mm seems to help for MRF? or 
  %   to run SPM without MRF
  if opt.SRP > 0
    % sharpening of raw image
    Ym = Ym + (Ym - smooth3(Ym)) / 2; 
  
    if cat_stat_nanmedian(Ym(round(Yp0(:))==1)) < cat_stat_nanmedian(Ym(round(Yp0(:))==3)) % T1w
      Yp0 = max(Yp0, 3*(Yp0>1.9 & Ym*3>2.5 & Ym*3<3.5) + 2.5*(Yp0>1.9 & Ym*3>2.25 & Ym*3<3.5));
    elseif cat_stat_nanmedian(Ym(round(Yp0(:))==1)) > cat_stat_nanmedian(Ym(round(Yp0(:))==3)) % T2w
      Yp0 = max(Yp0, 3*(Yp0>1.9 & Ym*3>2.5 & Ym*3<3.5) + 2.5*(Yp0>1.9 & Ym*3>2.25 & Ym*3<3.5));
    end
  end


  % simple filling
  [Yp0f,Ymf] = fillVentricle(Yp0,Ym,Ya,YMF,vx_vol); clear Ym YMF; 
  

  % position balancing map 
  % - this map will be used later to optimise the position map to have
  %   equal object/non-objects parts to stabilise the reconstruction
  %   (ie. thickening thin and thinning thick structures)
  Ypb = cat_vol_morph( Yp0f > 1.5 , 'ldc', 8 ) & cat_vol_morph( Yp0f < 2.5 , 'ldc', 8 ); 
  

  % prepare file and directory names
  [P,mridir,surfdir,ff] = cat_surf_createCS_fun('setFileNames',V0,job,opt); 
  

  % main loop for each surface structure 
  for si = 1:numel(opt.surf)
   
    % print something
    if si==1, fprintf('\n'); end; fprintf('%s:\n',opt.surf{si});
    clear Ynocerebrum

    % prepare longitudinal case if required 
    useprior = cat_surf_createCS_fun('setupprior',opt, job.BIDS(1).surfdir,P,si);


    % reduce for object area
    Ynocerebrum = ~(NS(Ya,opt.LAB.CB) | NS(Ya,opt.LAB.BV) | NS(Ya,opt.LAB.ON) | NS(Ya,opt.LAB.MB) | NS(Ya,opt.LAB.BS)); 
    switch opt.surf{si}
      case {'lh'}
        Ymfs   = Ymf  .* (Ya>0) .* Ynocerebrum .* (mod(Ya,2)==1); 
        Yp0fs  = Yp0f .* (Ya>0) .* Ynocerebrum .* (mod(Ya,2)==1); 
        Yside  = single(mod(Ya,2)==1); 
      case {'rh'}  
        Ymfs   = Ymf  .* (Ya>0) .* Ynocerebrum .* (mod(Ya,2)==0); 
        Yp0fs  = Yp0f .* (Ya>0) .* Ynocerebrum .* (mod(Ya,2)==0); 
        Yside  = single(mod(Ya,2)==0);  
      case {'cb'}
        Ymfs   = Ymf  .* (Ya>0) .* NS(Ya,opt.LAB.CB); 
        Yp0fs  = Yp0f .* (Ya>0) .* NS(Ya,opt.LAB.CB); 
        Yside  = zeros(size(Yp0fs),'single'); % new full cerebellum reconstruction
    end 
    
    
    % RD2025: mark cutting regions to avoid cortical modelling (i.e. to avoid/reduce thickness estimates)
    Ycutregion = min(1,max(0,Yp0f - 2 + NS(Ya,opt.LAB.VT)) * 1 .* ...
      smooth3( NS(Ya,opt.LAB.ON) | NS(Ya,opt.LAB.MB) | NS(Ya,opt.LAB.VT) | NS(Ya,opt.LAB.TH) | ~Yside )); 


    % removing background (smoothing to remove artifacts)
    [Yp0fs,Ymfs,Ycutregion,Ypbs,BB] = cat_vol_resize({Yp0fs,Ymfs,Ycutregion,Ypb}, 'reduceBrain', vx_vol, 4, smooth3(Yp0fs)>1.5); 
    
    % interpolation 
    imethod         = 'cubic'; % cubic is be better in general but we have to consider interpolation artifacts
    [Yp0fs,resI]    = cat_vol_resize(max(1,Yp0fs),'interp',V,opt.interpV,imethod);                  % interpolate volume
    Ymfs            = cat_vol_resize(max(1,Ymfs),'interp',V,opt.interpV,imethod);                   % interpolate volume 
    Ypbs            = cat_vol_resize(Ypbs,'interp',V,opt.interpV,imethod) > 0.5;                     
    Ycutregion      = cat_vol_resize(Ycutregion,'interp',V,opt.interpV,imethod);
    Ycutregiond     = cat_vol_smooth3X( cat_vol_morph(Ycutregion,'dd',1.5,opt.interpV), 4);         % a smooth version to mix maps
   
    
    % surface coordinate transformation matrix
    [Vmfs,Smat] = createOutputFileStructures(V,V0,resI,BB,opt,mridir,ff,si); 


    % thickness and position map estimation 
    Vmfs.dt = [16 1];
    spm_write_vol(Vmfs, Yp0fs);


    % Myelination corretion for the CS40 that is also part of cat_vol_pbtsimpleCS4
    % #### need final cleanup ####
    if myelinCorrection>0 && ~use_cat_vol_pbtsimple
      stime = cat_io_cmd(sprintf('  Myelin correction (%0.1f%%)',myelinCorrection*100),'g5'); 
      opt.myelinCorrection = myelinCorrection; 
      Yp0fs = myelincorrection(Yp0fs, vx_vol, opt,P,Vmfs,si); 
      spm_write_vol(Vmfs, Yp0fs);
      fprintf('%5.0fs\n',etime(clock,stime)); 
    end  


    % Cerebellar mid layber segmentation 
    % #### This is not fully working ###
    if cat_io_contains(opt.surf{si},'cb')
      stime = cat_io_cmd('  Cerebellar preparation','g5'); 
      if 1
        %%
        if 0
          Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));
          spm_write_vol(Vmfs, max(0,(Ymfs-1)/2));
          Vgm = Vmfs; Vgm.fname = spm_file(Vgm.fname,'suffix','_seg1'); spm_write_vol(Vgm,Yp0toC(Yp0fs,2));
          Vwm = Vmfs; Vwm.fname = spm_file(Vwm.fname,'suffix','_seg2'); spm_write_vol(Vwm,Yp0toC(Yp0fs,3));
          Vwm = Vmfs; Vwm.fname = spm_file(Vwm.fname,'prefix','c_'); spm_write_vol(Vwm,(Ymfs - 1)/2);
          Vwm = Vmfs; Vwm.fname = spm_file(Vwm.fname,'suffix','_pcereb','prefix','c_'); spm_write_vol(Vwm,Yp0fs>1);
          %%
          Pcb = suit_amap( struct('gray',{{Vmfs.fname}}) );
          Vcba  = spm_vol(Pcb.cba{1}); 
          Yp0fs = spm_read_vols(Vcba);
        else
          %%
          Pcb = suit_amap('run',Yp0fs,max(0,Ymfs - 1)/2,vx_vol);
          Yp0fs = min(3,max(1,Pcb.Yml * 3));
          Yppii = Pcb.Yta;
        end

      else
        Yp0fso   = Yp0fs;
        
  
        %% skull-stripping
        Yp0fs    = Yp0fso;
        [~,Ycbd] = cat_vol_downcut( single(cat_vol_morph( smooth3(Ymfs)>2.1,'l')),(Yp0fs.^.1) - 1,1); 
        cbdth    = cat_stat_nanmedian(Ycbd(Yp0fs>1.1 & Yp0fs<1.25)); 
        Ycbs     = Yp0fs>1.1 & Ycbd < cbdth; 
        
        % enhanced segmentation 
        Ymfscbw  = cat_vol_localstat( Ymfs-1, (Yp0fs>1.25 & Yp0fs<2.75)  & Ycbs, 2/mean(opt.interpV), 8); 
        Ymfscbw  = max(Ymfscbw, cat_vol_localstat( Ymfs-2, (Yp0fs>2.75) & Ycbs, 2/mean(opt.interpV), 8)*.5 + ...
                                cat_vol_localstat( Ymfs-2, (Yp0fs>2.75) & Ycbs, 2/mean(opt.interpV), 3)*.5); 
        Ymfscbw  = cat_vol_approx( Ymfscbw ,'rec',2); 
        
        % composition
        Yp0fs    = max(1, real(0.95+0.05.*(1-Ycbd./cbdth).^.01)  .* smooth3(Ycbs)  .*  max(1,Ymfs) ./ Ymfscbw);
        cat_sanlm( Yp0fs , 1,3); 
      end
      fprintf('%5.0fs\n',etime(clock,stime)); 
    end
    
    %%
    stime = cat_io_cmd(sprintf('  Thickness estimation (%0.2f mm%s)',opt.interpV,native2unicode(179, 'latin1')),'g5');  % fprintf('\n'); 
    if use_cat_vol_pbtsimple
      [Yth1i,Yppi,Ymfsc] = cat_vol_pbtsimpleCS4(Yp0fs, opt.interpV,...
        struct('myelinCorrection',myelinCorrection, 'verb',1, 'gyrusrecon',1, 'wmnoise',opt.wmnoise, 'eidist',opt.eidist));
    else
      cmd = sprintf('CAT_VolThicknessPbt  -correct-voxelsize 0   -median-filter 2   -downsample 0 "%s" "%s" "%s"', ...
        Vmfs.fname, P(si).Pgmt, P(si).Pppm);
      cat_system(cmd,opt.verb-3);
      Vgmt  = spm_vol(P(si).Pgmt); Yth1i = spm_read_vols(Vgmt); 
      Vppi  = spm_vol(P(si).Pppm); Yppi  = spm_read_vols(Vppi); 
      Ymfsc = Yp0fs; 

      % correction of general offset in mm 
      Yth1i = max(0,Yth1i - 0.28); 
    end
    if cat_io_contains(opt.surf{si},'cb')
      % here we use the normalized T1w data directly 
      Yppi  = Yppii; %max(0,min(2,Yp0fs-1)/2); 
      Yp0fs = Yp0fs;
      Ymfsc = Yp0fs;
    end

    %% prepare thickness output
    Yth1t = cat_vol_resize(Yth1i,'deinterp',resI);                         % back to original resolution
    Yth1t = cat_vol_resize(Yth1t+1,'dereduceBrain',BB)-1;                  % adding background
    Yth   = max(Yth,Yth1t .* Yside);                                       % save on main image
    
    % prepare Ypp
    Ymfsc = Yp0fs .* Ycutregion + (1-Ycutregion) .* Ymfsc; % ##################
    Vppm  = Vmfs;  Vppm.fname = P(si).Pppm;  Vppm.dt(1) = 16; Vppm.pinfo(1) = 1; Vpp = spm_write_vol(Vppm, Yppi); 
    Yppi  = max(0,min(1,(cat_vol_smooth3X(Yppi,4)-.5)*4 + .5)) .* Ycutregiond  + (1-Ycutregiond) .* Yppi; 


    %%
    if opt.vol
      S = struct(); P = '';
      if opt.verb<2, fprintf('%5.0fs',etime(clock,stime)); end %#ok<*DETIM>
      continue; % ############# deactive only for tests ################
    end
 
    
    %  surface creation 
    %  --------------------------------------------------------------------
    %  In case of the central surface frequency patter between sulci and 
    %  gyri is more harmonized and even lower surface reconstruction 
    %  resolution (eg. 1.5 and 2.0) are possible. 
    %  The resolution of 1.2 mm result in a similar number of vertices/faces
    %  as the CS2 pipeline. High resolution surfaces are quite large (disk space) 
    %  and 
    %  --------------------------------------------------------------------
    if useprior 
      fprintf('\n');
      stime = cat_io_cmd('  Load and refine subject average surface','g5','',opt.verb); %if opt.verb>2, fprintf('\n'); end
      CS  = loadSurf(P(si).Pcentral); 
    else
      % optimized downsampling of the the Ypp map and
      if isscalar(opt.surf), time_sr = clock; end % temporary for tests 
      [Vppmi,rel] = exportPPmap( Yp0 .* Ypb, Ymfsc, Yppi , Vmfs, Ypbs, 0, si, opt, '', ff);
      Vp0 = Vmfs; Vp0.fname = spm_file(P(si).Pp0,'suffix','_tmp'); 
      spm_write_vol(Vp0, double(Yp0fs));
 
      % static/dynamic surface threshold
      if 0 %job.extopts.SRP == 0
        % default values
        gycon = .5;
      else
        % using a global dynamic threshold to avoid topology defects that
        % is typically resulting in longer time for surface creation and
        % severe closing of sulci.
        % In case of low resolution reconstruction the threshold has to be
        % higher to avoid fatal closing.
        relwm = sum(max(0,Ymf(:)-2)) / sum(min(1,Ymf(:))); 
        gycon = max(0,1 - max(0,relwm * 3 - 1) * 2);
        gycon = min(.85,.65 - .1*gycon  +  max(0,mean(vx_vol-.5)/3)); 
      end
     
      if cat_io_contains(opt.surf{si},'cb')
        %% cerebellar mid layber segmentation 
        %  RD20260320: WARNING: Be careful with edits here because this block is working magicaly.
        Yppi     = (cat_vol_median3(Ymfs,Yp0fs>1,Yp0fs>0,-0.1)); 
        [~,Ycbd] = cat_vol_downcut( single(smooth3(Yp0fs)>2.1),(Yp0fs.^.1) - 1,1); 
        cbdth    = cat_stat_nanmedian(Ycbd(Yppi>1.1 & Yppi<1.5)); 
        Ycbs     = Yppi>1.1 & Ycbd < cbdth; 
        Ymfscbw  = cat_vol_localstat( Yppi-1, (Yp0fs>1.25 & Yp0fs<2.75) & Ycbs, min(10,2 / opt.interpV), 8); 
        Ymfscbw  = cat_vol_approx( Ymfscbw ); 
        Yppi     = (0.95+0.05.*max(eps,1-Ycbd./cbdth).^.01)  .* smooth3(Ycbs) .* (Yppi-1) ./ Ymfscbw / 2; 
        %Yppi     = Yppi + 0.5*(Yppi-smooth3(Yppi)); 
        %Yppi     = cat_vol_median3(Yppi,Yppi>0); 
        if 0
          %%
          [~,~,Yde1] = cat_vol_eidist(single(Yppi>.9),Yppi.^10e-8); 
          [~,~,Yde2] = cat_vol_eidist(single(Yppi<.1),(1-Yppi).^10e-8); 
          Yppi = (Yppi + (Yde2./(Yde1+Yde2)))/2; 
        end

        %%
        for i=1:4, cat_sanlm(Yppii,1,3); Yppi = Yppii; end 
      end

%      cmd = sprintf('CAT_VolMarchingCubes  "%s" "%s"  -thresh .75 -iter 1 -median-filter 0', Vppmi.fname, P(si).Pcentral); %#ok<NASGU>
%      CS  = loadSurf(P(si).Pcentral); 
%      AA  = sum( cat_surf_fun( 'area', CS )); % in mm2 
      
     %%
      vxth     = 200 * (1+2*cat_io_contains(opt.surf{si},'cb'));  % voxel threshold of accepted changes to adjust the topology (but in mm!) 
      ECdiffth = 10;   % accaptable EC of the final surface 
      stime = cat_io_cmd(sprintf('  Create initial surface (%0.2f mm)',opt.interpV),'g5','',opt.verb,stime); %if opt.verb>2, fprintf('\n'); end
      for ECiter = 1 % this loop is to allow futher topology iterations if the first one was not good enough
        for thi = 1; %*(1+1*cat_io_contains(opt.surf{si},'cb')) % with this loop we further increase the threshold of the Ypp map by .1 (eg. more WM like surface)
          % Main initial surface creation with topology correction.
          if exist(P(si).Pcentral,'file'), delete(P(si).Pcentral); end % have to delete it to get useful error messages in case of reprocessing/testing
          if cat_io_contains(opt.surf{si},'cb')
            % RD202603: just as quick placeholder ...
            opt0=opt; opt0.reconres = 0.25; 
            gycon = .55 + 0.025*(thi);  %max(.5,sum(Yppi(:)>.1 & Yppi(:)<.5) ./ sum(Yppi(:)>.1 & Yppi(:)<.9)) + 0.05*(thi); 
            [Vppmi,rel] = exportPPmap( Yp0 .* Ypb, Ymfsc, Pcb.Yta , Vmfs, Ypbs, 0, si, opt0, '', ff); % have to rewrite it
%            spm_write_vol(Vppmi, (Pcb.Yppi*.5 + 0.5*(Pcb.Yta))  ); 
            if 1
              cmd = sprintf('CAT_VolMarchingCubes  "%s" "%s"  -thresh "%0.4f" -iter 1 -median-filter "%d"',  ...
                Vppmi.fname, P(si).Pcentral, gycon, 0); %#ok<NASGU>
            else
              spm_write_vol(Vp0,  Yp0fs ); 
              cmd = sprintf('CAT_VolMarchingCubes  "%s" "%s"  -thresh "%0.4f" -iter %d -label "%s" -median-filter "%d"',  ...
                Vppmi.fname, P(si).Pcentral, gycon, ECiter, Vp0.fname, min(2,round( (opt.wmnoise*20) / opt.interpV ) )); 
            end
          else
            gycon = gycon + .1*(thi>1);
            [Vppmi,rel] = exportPPmap( Yp0 .* Ypb, Ymfsc, Yppi , Vmfs, Ypbs, 0, si, opt, '', ff); % have to rewrite it
            cmd = sprintf('CAT_VolMarchingCubes  "%s" "%s"  -thresh "%0.4f" -iter %d -label "%s" -median-filter "%d"',  ...
              Vppmi.fname, P(si).Pcentral, gycon, ECiter, Vp0.fname, min(2,round( (opt.wmnoise*20) / opt.interpV ) )); %#ok<NASGU>
          end  
          txt = evalc('cat_system(cmd ,3)');

          % save surface characteristics 
          res.ECsurf       = opt.surf;
          res.EC(si)       = str2double( txt(max(strfind(txt,':'))+1 : max(strfind(txt,'('))-2) ); 
          res.ECmodvx(si)  = str2double( txt(max(strfind(txt,'('))+1 : max(strfind(txt,'voxel'))-1) );
          res.ECmodwmp(si) = res.ECmodvx(si) / sum( round(Yp0fs(:))==3 ) * 100;
 
          % evaluate surface characteristics
          if gycon>.85  ||  res.ECmodvx(si) * prod(vx_vol) <= vxth  &&  abs(res.EC(si) - 2) <= ECdiffth 
            if opt.verb>1 && thi>0, cat_io_cprintf([0 .5 0],sprintf( ...
              '\n    Accepted configuration (EC=%d, %5dvx, %0.2f%%) for Ypp>%0.2f   ', ...
              res.EC(si), res.ECmodvx(si), res.ECmodwmp(si), gycon));
            end
            break;      
          else
            if opt.verb>1, cat_io_cprintf([.7 0 0.5],sprintf( ...
              '\n    Sulcal blurring (EC=%3d, %5dvx ~ %0.2f%%) for Ypp>%0.2f        ', ...
              res.EC(si), res.ECmodvx(si), res.ECmodwmp(si), gycon));
            end
          end
        end

        % evaluate surface characteristics
        if gycon>.85  ||  res.ECmodvx(si) * prod(vx_vol) <= vxth  &&  abs(res.EC(si) - 2) <= ECdiffth 
          break; 
        end
      end

      if abs(res.EC(si) - 2) > ECdiffth
        cat_io_addwarning([mfilename ':EC'],sprintf('Incorrect Euler Number i.e. reminding topological defects (EC=%d).',res.EC(si)) ,1,[1 1]);
      end
      if res.ECmodvx(si)  * prod(vx_vol) > vxth
        cat_io_addwarning([mfilename ':MarchingCubes'],sprintf('Increased topological adjustments of %d voxels (~%0.2f%%%% of WM).',...
          res.ECmodvx(si),  res.ECmodwmp(si)), 1,[1 1]);
      end
      if cat_io_contains(opt.surf{si},'cb')
        Vppm = spm_write_vol(Vppm, Pcb.Yppi); % rewrite unweighted position map
      else
        Vppm = spm_write_vol(Vppm, Yppi); % rewrite unweighted position map
      end

      % load and check surface and 
      % define lower surface resolution based on the surface area of the subject
      % use 80000 as lower limit what is similar to vdist=4
      CS  = loadSurf(P(si).Pcentral); 
      A   = sum( cat_surf_fun( 'area', CS )); % in mm2 
      nCS = min( size(CS.faces,1) , max( 80000, min( 360000, round(A * 8 / opt.vdist^2) ))); 
      nCS = nCS * (1+3*cat_io_contains(opt.surf{si},'cb'));
      EC0 = size(CS.vertices,1) + size(CS.faces,1) - size(spm_mesh_edges(CS),1);  
   
      % surface reduction 
      % - low surface resolution are highly important for fast processing and smaller disk usage
      % - less critical for thickness but surface quality (intensity/position values) get worse
      % - it is important to keep the topology 
      % - a more aggressive 
      stimesr = clock; 
      nCS0 = size(CS.faces,1); 
      saveSurf(CS,P(si).Pcentral); 
      cmd  = sprintf('CAT_SurfReduce -aggr "7" -ratio "%0.2f" "%s" "%s"', nCS/nCS0 , P(si).Pcentral, P(si).Pcentral); %5-0.05*opt.reconres);
      cat_system(cmd ,opt.verb-3);  
      CS = loadSurf(P(si).Pcentral); 
      nCS1 = size(CS.faces,1); 
      
      if nCS1 > nCS * 1.5 
        cat_io_addwarning('cat_surf_createCS4:ReductionFailed', ...
          sprintf('Surface reduction incomplete due to topology constrains (%d).',(nCS1/1000-nCS/1000) ./ (nCS0/1000-nCS/1000)),0,[0 1]); 
      end 
      cat_io_cmd(sprintf('  Reduce surface (nF: %0.0fk > ~%0.0fk)', nCS0/1000, nCS1/1000),'g5','',opt.verb,stime + (clock - stimesr));
      stime = stimesr; 
    end


    % remove artifacts 
    % - this can maybe partially avoided by the more aggressive surface reduction
    % - although the filter reduces some artefacts in problematic cases the surface values are worse and the issue are still not fully solved
    % - tried to use it after correction of self-intersections but this made it even worse
    % - minor improvement of values (IE+.001,PE+0.002) but visually strong reduction of local outliers >> useful
    if ~cat_io_contains(opt.surf{si},'cb') % opt.SRP > 0
      stime = cat_io_cmd('  Filter topology artifacts','g5','',opt.verb,stime);
      CS = smoothArt(Yth1i,P,CS,Smat,Vppm,1,si,opt,1,0); % last options - method & refine
    end
    saveSurf(CS,P(si).Pcentral); 
    
    % get PBT thickness and update cutregion
    if setcut2zero
      facevertexcdatanocut = 1 - cat_surf_fun('isocolors',Ycutregion,CS.vertices,Smat.matlabIBB_mm); 
      facevertexcdatanocut = max(eps,min(1,(cat_surf_fun('smoothcdata',CS,facevertexcdatanocut,4) - .5) * 2)); 
    else
      facevertexcdatanocut = 1; 
    end
    facevertexcdata = cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm) .* facevertexcdatanocut; 
    cat_io_FreeSurfer('write_surf_data',P(si).Ppbt,facevertexcdata);
    


    % surface deformation (after smoothArt) 
    % CAT_DeformSurf in combination with smoothing correction CAT_Central2Pial 
    % allows PP values that are closer to the central position (std 0.049 vs. 0.73)
    % what improves also intensity and position values of the inner/outer surfaces!
    stime = cat_io_cmd('  Optimize surface','g5','',opt.verb,stime); 
    if opt.SRP==0
      cmd = sprintf('CAT_SurfDeform -remove_intersect -iter 75 -isovalue 0.5 -w1 0.1 -w2 0.1 -w3 1.0 -sigma 0.2 "%s" "%s" "%s"', Vppm.fname, P(si).Pcentral, P(si).Pcentral);
      cat_system(cmd,opt.verb-3);
    else
      % further iteration allow (only) slight improvements 
      cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none  0  1  -1  .1 ' ...           
                'avg  -0.1  0.1 .2  .1  5  0 "0.5"  "0.5"  n 0  0  0 %d  %g  0.0 0'], ...    
                 Vppm.fname,P(si).Pcentral,P(si).Pcentral,50,0.001);  
      cat_system(cmd,opt.verb-3);
      cmd = sprintf('CAT_Central2Pial -equivolume -weight 0.5 "%s" "%s" "%s" 0',P(si).Pcentral,P(si).Ppbt,P(si).Pcentral);
      cat_system(cmd,opt.verb-3);
    end
    
    % Remove self-intersections and update thickness (as the surface topology changes too)
% ######### need further refinement    
    if ~useprior %&& opt.SRP > 2 % NOT WORKING IN CS40
      cmd = sprintf('CAT_SurfRemoveIntersections "%s" "%s"',P(si).Pcentral,P(si).Pcentral);
      cat_system(cmd,opt.verb-3);
      cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none  0  1  -1  .1 ' ...           
                  'avg  -0.1  0.1 .2  .1  5  0 "0.5"  "0.5"  n 0  0  0 %d  %g  0.0 0'], ...    
                   Vppm.fname,P(si).Pcentral,P(si).Pcentral,25,0.001);  
        cat_system(cmd,opt.verb-3);
        cmd = sprintf('CAT_Central2Pial -equivolume -weight 0.5 "%s" "%s" "%s" 0',P(si).Pcentral,P(si).Ppbt,P(si).Pcentral);
        cat_system(cmd,opt.verb-3);
      CS = loadSurf(P(si).Pcentral);
      facevertexcdatanocut = 1 - cat_surf_fun('isocolors',Ycutregion,CS.vertices,Smat.matlabIBB_mm); 
      facevertexcdata = max(eps,cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm) .* facevertexcdatanocut); 
      cat_io_FreeSurfer('write_surf_data',P(si).Ppbt,facevertexcdata);
      cat_io_FreeSurfer('write_surf_data',P(si).Pthick,facevertexcdata);
    end


    %% Collision correction by Delaunay triangularization
    %  --------------------------------------------------------------------
    %  New self-intersection correction that uses different detections of
    %  self-intersections (SIDs; RY vs. PBT) with/without further optimization. 
    %  It does not fully avoid self-intersections because some are already 
    %  in the CS and some other required strong changes that result in worse
    %  thickness results.
    %  RD202508: Optimization does not fully work on the simple PP map.
    %            The optimization torwards 0.05 and 0.95 on the interpolated maps still result in light thickness overestimation. 
    if opt.SRP > 2
      facevertexcdatac = facevertexcdata; 
      % quick correction with local optimization 
      [CS,facevertexcdatac,SIs] = cat_surf_fun('collisionCorrectionPBT',CS,facevertexcdatac .* facevertexcdatanocut,Ymfsc,Yppi, ... 
        struct('optimize',0,'verb',isscalar(opt.surf)*2,'mat',Smat.matlabIBB_mm,'vx_vol',vx_vol,'CS4',0)); % CS4 is not working
      if isscalar(opt.surf), fprintf('\b\b'); end
      % accurate but slow function (especially for multiple jobs) 
      [CS,facevertexcdatac,SIs] = cat_surf_fun('collisionCorrectionRY' ,CS,facevertexcdatac .* facevertexcdatanocut,Yppi,...
        struct('Pcs',P(si).Pcentral,'verb',isscalar(opt.surf)*2,'mat',Smat.matlabIBB_mm,'accuracy',1/2^2));
      if isscalar(opt.surf), fprintf('\b\b'); end
      facevertexcdata = max(eps,max(facevertexcdata,facevertexcdatac) .* facevertexcdatanocut); 
      cat_io_FreeSurfer('write_surf_data',P(si).Ppbt,facevertexcdata);
      saveSurf(CS,P(si).Pcentral); 
    end
   

    % evaluate and save results
    if isempty(stime), stime = clock; end
    fprintf('%5.0fs',etime(clock,stime)); stime = []; if 1, fprintf('\n'); end %debug


    % create white and central surfaces
    if create_white_pial
      stime = cat_io_cmd('  Create pial and white surface','g5','',opt.verb,stime); 
      %if create_white_pial == 1 
        % quick
        cmd = sprintf('CAT_Central2Pial "%s" "%s" "%s" 0.5',P(si).Pcentral,P(si).Ppbt,P(si).Ppial);
        cat_system(cmd,opt.verb-3);
        cmd = sprintf('CAT_Central2Pial "%s" "%s" "%s" -0.5',P(si).Pcentral,P(si).Ppbt,P(si).Pwhite);
        cat_system(cmd,opt.verb-3); 
      if create_white_pial>1
        % write segmentation if necessary (and delete it later)
        if ~exist(spm_file(P(si).Pp0,'suffix','_tmp'),'file')
          Vp0 = Vmfs; Vp0.fname = spm_file(P(si).Pp0,'suffix','_tmp'); 
          spm_write_vol(Vp0, Yp0fs);
        end

        % estimate pial and white surface with a limited thickness map
        facevertexcdata = min(6,max(eps, cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm) .* facevertexcdatanocut)); 
        cat_io_FreeSurfer('write_surf_data', P(si).Pthick,   facevertexcdata);
        spm_write_vol(Vmfs, Yp0fs);
        cmd = sprintf('CAT_Surf2PialWhite "%s" "%s" "%s" "%s" "%s"', ...
          P(si).Pcentral, P(si).Pthick, Vmfs.fname, P(si).Pp0, P(si).Pwhite);
        cat_system(cmd,opt.verb-3); 
        
        % create new central surface and update thickness values
        cmd = sprintf('CAT_AverageSurfaces "%s" "%s" -avg "%s" ', ...
          P(si).Ppial, P(si).Pwhite, P(si).Pcentral);
        cat_system(cmd,opt.verb-3); 
        CS = loadSurf(P(si).Pcentral);
        facevertexcdata = max(eps, cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm) .* facevertexcdatanocut); 
        cat_io_FreeSurfer('write_surf_data', P(si).Ppbt, facevertexcdata);
      end
    end


    if isscalar(opt.surf)
      % only for test visualization  
      fprintf('%5.0fs',etime(clock,stime)); 
      cat_surf_createCS_fun('quickeval',V0,Vpp,Ymfs,Yppi,CS,P,Smat,res,opt,EC0,si,time_sr,4); 
      return
    end
    if exist(P(si).Pthick,'file'), delete(P(si).Pthick); end
        

    %% skip that part if a prior image is defined
    if ~useprior && ~skip_registration
      % spherical surface mapping of the final corrected surface with minimal areal smoothing
      stime = cat_io_cmd('  Spherical mapping with areal smoothing','g5','',opt.verb,stime); 
      cmd = sprintf('CAT_Surf2Sphere "%s" "%s" %d',P(si).Pcentral,P(si).Psphere,6);
      cat_system(cmd,opt.verb-3);
      
      % spherical registration to fsaverage template
      stime = cat_io_cmd('  Spherical registration','g5','',opt.verb,stime); %#ok<NASGU>
      cmd = sprintf('CAT_SurfWarp -steps 2 -avg -i "%s" -is "%s" -t "%s" -ts "%s" -ws "%s"', ...
        P(si).Pcentral,P(si).Psphere,P(si).Pfsavg,P(si).Pfsavgsph,P(si).Pspherereg);
      cat_system(cmd,opt.verb-3);
    end  


    if debug || cat_get_defaults('extopts.expertgui')>1, fprintf('\n'); end
    if opt.thick_measure == 1
      % use central surface and thickness to estimate Freesurfer thickness
      cmd = sprintf('CAT_SurfDistance -mean -thickness "%s" "%s" "%s"',P(si).Ppbt,P(si).Pcentral,P(si).Pthick);
      cat_system(cmd,opt.verb-3);
      
      % apply upper thickness limit
      % here the 5 mm thickness limit of FreeSurfer might be better rather then our 6 mm 
      facevertexcdata = max(eps, cat_io_FreeSurfer('read_surf_data',P(si).Pthick) .* facevertexcdatanocut);  
      cat_io_FreeSurfer('write_surf_data', P(si).Pthick, min(opt.thick_limit,facevertexcdata) );  
    else
      % otherwise simply use the original values of the PBT map  
      % WARNING: 
      %   The values of the ?h.pbt.* files are used to estimate further 
      %   surfaces and therefore corrected for self-intersections!
      %   We use here the orignal PBT values as they are less depending on 
      %   local highly individual features.  
      facevertexcdata = max(eps, cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm) .* facevertexcdatanocut); 
      cat_io_FreeSurfer('write_surf_data', P(si).Ppbt,   facevertexcdata);
      cat_io_FreeSurfer('write_surf_data', P(si).Pthick, min(opt.thick_limit,facevertexcdata) );  
    end
    fprintf('\n');
    

    % correct thickness based on folding pattern, but smaller thickness values 
    % are corrected less strongly than larger thickness values
    if opt.foldingcorrection
      cmd = sprintf('CAT_SurfCorrectThicknessFolding -slope 1.0 -max "%f" "%s" "%s" "%s"', opt.thick_limit, P(si).Pcentral, P(si).Pthick, P(si).Pthick);
      cat_system(cmd,opt.verb-3);
      % only minimum
      facevertexcdatac = cat_io_FreeSurfer('read_surf_data',P(si).Pthick);
      facevertexcdata  = min(facevertexcdata,facevertexcdatac); 
      cat_io_FreeSurfer('write_surf_data', P(si).Pthick, min(opt.thick_limit,facevertexcdata) );  
    end


    % final surface evaluation 
    res.(opt.surf{si}).createCS_final = cat_surf_fun('evalCS', ...
      loadSurf(P(si).Pcentral), cat_io_FreeSurfer('read_surf_data',P(si).Ppbt), cat_io_FreeSurfer('read_surf_data',P(si).Pthick), ...
      Ymfs, Yppi, P(si).Pcentral, Smat.matlabIBB_mm, debug + (cat_get_defaults('extopts.expertgui')>1), cat_get_defaults('extopts.expertgui')>1);
 

    % average final values
    FNres = fieldnames( res.(opt.surf{si}).createCS_final );
    for fnr = 1:numel(FNres)
      if ~isfield(res,'final') || ~isfield(res.final,FNres{fnr})
        res.final.(FNres{fnr}) = res.(opt.surf{si}).createCS_final.(FNres{fnr}) / numel(opt.surf);
      else
        res.final.(FNres{fnr}) = res.final.(FNres{fnr}) + res.(opt.surf{si}).createCS_final.(FNres{fnr}) / numel(opt.surf);
      end
    end
    
    % create output structure
    S.(opt.surf{si}) = struct('faces',CS.faces,'vertices',CS.vertices,'th1',facevertexcdata); 
    clear facevertexcdata Yth1i CS; 
    
    if exist('Vp0','var') && exist(Vp0.fname ,'file'), delete(Vp0.fname); end
    if exist('Vppm','var') && exist(Vppm.fname ,'file'), delete(Vppm.fname); end
    if ~debug && exist('Vpp','var') && exist(Vpp.fname ,'file') && ~opt.outputpp.native, delete(Vpp.fname); end
    if ~debug && exist('Vgmt','var') && exist(Vgmt.fname ,'file'), delete(Vgmt.fname); end
    if ~debug && exist('Vmfs','var') && exist(Vmfs.fname ,'file'), delete(Vmfs.fname); end
  
    % processing time per side for manual tests
    if si == numel(opt.surf) && si == 1
      cat_io_cmd('  ','g5','',opt.verb);
      fprintf('%5ds\n',round(etime(clock,cstime)));
    end
  end  
  

  % calculate surface quality parameters for all surfaces
  res = cat_surf_createCS_fun('addSurfaceQualityMeasures',res,opt);
 

  % print final stats
  cat_surf_createCS_fun('evalProcessing',res,opt,P,V0); 
  %evalProcessing(res,opt,P,V0)

end
%=======================================================================
function saveSurf(CS,P)
  save(gifti(struct('faces',CS.faces,'vertices',CS.vertices)),P,'Base64Binary');
end
%=======================================================================
function CS1 =loadSurf(P)
  % add 1s because sometimes surface is not yet ready to read...
  if ~exist(P,'file')
    pause(3)
    if ~exist(P,'file')
      pause(1)
      error('Surface file %s could not be found due to previous processing errors.',P);
    end 
  end 
  
  try
    CS = gifti(P);
  catch
    error('Surface file %s could not be read due to previous processing errors.',P);
  end
  
  CS1.vertices = CS.vertices; CS1.faces = CS.faces; 
  if isfield(CS,'cdata'), CS1.cdata = CS.cdata; end
end
%==========================================================================
function [Vmfs,Smat] = createOutputFileStructures(V,V0,resI,BB,opt,mridir,ff,si)
    matI              = spm_imatrix(V.mat); 
    matI(7:9)         = sign( matI(7:9))   .* repmat( opt.interpV , 1 , 3); 
    matiBB            = spm_imatrix(V.mat   * [eye(4,3) [ (BB.BB([1,3,5])' - 1) ; 1]]); 
    matIBB            = matiBB; 
    matIBB(7:9)       = sign( matiBB(7:9)) .* repmat( opt.interpV , 1 , 3); 
    Smat.matlabi_mm   = V.mat * [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];                % CAT internal space
    Smat.matlabI_mm   = spm_matrix(matI) * [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];     % PBT interpolated space
    Smat.matlabIBB_mm = spm_matrix(matIBB) * [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];   % PBT interpolated
    Smat.matlabiBB_mm = spm_matrix(matiBB) * [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];   % PBT interpolated

    Vmfs = resI.hdrN;
    Vmfs.pinfo = V0.pinfo;
    Vmfs.fname = fullfile(mridir, sprintf('%s_seg-%s.nii',ff,opt.surf{si}));
    if isfield(Vmfs,'dat'),     Vmfs = rmfield(Vmfs,'dat'); end
    if isfield(Vmfs,'private'), Vmfs = rmfield(Vmfs,'private'); end
    matiBB = spm_imatrix(V.mat   * [eye(4,3) [ (BB.BB([1,3,5])' - 1) ; 1]]); 
    Vmfs.mat(1:3,4) = matiBB(1:3); 
end
%==========================================================================
function [Vpp,Vppr] = prepareDownsampling(Vmfs,Ypps,pp0_surffolder,ff,opt,si)
%prepareDownsampling. Create nifti structure for low resolution.  
    Vpp        = Vmfs; 
    Vpp.pinfo  = repmat([1;0], 1,size(Ypps,3));
    Vpp.dat(:,:,:) = smooth3(single(Ypps)); 
    Vpp.dt(1)  = 16; 
    Vppr       = Vpp; 
    imat       = spm_imatrix(Vpp.mat); 
    Vppr.dim   = round(Vpp.dim .* opt.interpV./opt.reconres);
    imat(7:9)  = opt.reconres .* sign(imat(7:9));
    Vppr.mat   = spm_matrix(imat); clear imat;  
    Vppr.fname = fullfile(pp0_surffolder, sprintf('%s.pp.%s.nii',opt.surf{si},ff));  
end            
%==========================================================================
function [Yp0f,Ymf] = fillVentricle(Yp0,Ym,Ya,YMF,vx_vol)
%fillVentricle. Simple filling of ventricles by a closing around a mask.

  NS  = @(Ys,s) Ys==s | Ys==s+1; 
  LAB = cat_get_defaults('extopts.LAB');

  % simple filling by the YMF mask  
  Yp0f = max(Yp0  ,min(1,YMF & ~NS(Ya,LAB.HC) & ~( cat_vol_morph( NS(Ya,LAB.HC),'dd',2,vx_vol)))); 
  Ymf  = max(Ym   ,min(1,YMF & ~NS(Ya,LAB.HC) & ~( cat_vol_morph( NS(Ya,LAB.HC),'dd',2,vx_vol)))); 

  % close gaps in Yp0f
  Yp0fs = cat_vol_smooth3X(Yp0f,1);
  Ytmp  = cat_vol_morph(YMF,'dd',3,vx_vol) & Yp0fs>2.1/3;
  Yp0f(Ytmp) = max(min(Yp0(Ytmp),0),Yp0fs(Ytmp)); clear Ytmp Yp0fs; 
  Yp0f  = Yp0f * 3;
  
  % close gaps in Ymfs
  Ymfs  = cat_vol_smooth3X(Ymf,1);
  Ytmp  = cat_vol_morph(YMF,'dd',3,vx_vol) & Ymfs>2.1/3;
  Ymf(Ytmp) = max(min(Ym(Ytmp),0),Ymfs(Ytmp)); clear Ytmp Ymfs; 
  Ymf   = Ymf * 3;
end
%==========================================================================
function CS = smoothArt(Yth1i,P,CS,Smat,Vppm,facevertexcdatanocut,si,opt,method,refine)
%% Topology artifact correction (after first deformation):
%  Although the topology is fine after correction the geometry often
%  suffers by small snake-like objects that were deformed close to the 
%  central layer, resulting in underestimation of the FS thickness.
    
  saveSurf(CS,P(si).Pcentral); 
  cmd = sprintf('CAT_SurfDeform "%s" "%s" "%s"', Vppm.fname, P(si).Pcentral, P(si).Pcentral);
  cat_system(cmd,opt.verb-3);
  cmd = sprintf('CAT_SurfRemoveIntersections "%s" "%s"',P(si).Pcentral,P(si).Pcentral);
  cat_system(cmd,opt.verb-3);
  cmd = sprintf('CAT_SurfDeform "%s" "%s" "%s"', Vppm.fname, P(si).Pcentral, P(si).Pcentral);
  cat_system(cmd,opt.verb-3);
  CS = loadSurf(P(si).Pcentral);
 
  pbtthick = cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm) .* facevertexcdatanocut; 
  cat_io_FreeSurfer('write_surf_data',P(si).Ppbt,pbtthick);
  cmd = sprintf('CAT_SurfDistance -mean -thickness "%s" "%s" "%s"', P(si).Ppbt, P(si).Pcentral, P(si).Pthick);
  cat_system(cmd,opt.verb-3);
  fsthick  = cat_io_FreeSurfer('read_surf_data',P(si).Pthick);  
  
  % define outlier maps
  tart   = (pbtthick - fsthick)>.5 & pbtthick<2;                        % artifacts from topology correction   
  [~,dI] = unique(CS.vertices,'rows'); SM = true(size(CS.vertices,1),1); SM(dI) = false; % vertices with same coordinates
  iarea  = 1./cat_surf_fun('area',CS);                                                   % face area to identify tiny faces (=artifacs)
  
  if method % default function - working
    CS = cat_surf_smoothr(CS,tart | SM | iarea>10000, size(CS.vertices,1)/10, 10); % local
    CS = cat_surf_smoothr(CS,tart>=0, 1, 1); % all
    saveSurf(CS,P(si).Pcentral); 
  else
    % this is not really working
    cat_io_FreeSurfer('write_surf_data',P(si).Pmsk,tart | SM | iarea>1000);  % create a mask for filtering
    cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s"', P(si).Pcentral, P(si).Pmsk, round(size(CS.vertices,1)/100), P(si).Pmsk); % local
    cat_system(cmd,opt.verb-3); delete(P(si).Pmsk);
    cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g"', P(si).Pcentral, P(si).Pmsk, 10); % all 
    cat_system(cmd,opt.verb-3); delete(P(si).Pmsk);
  end 

  if refine
    cmd  = sprintf('CAT_Central2Pial -equivolume -weight 0.55 "%s" "%s" "%s" 0',P(si).Pcentral,P(si).Ppbt,P(si).Pcentral);
    cat_system(cmd,opt.verb-3);
  end

  cmd = sprintf('CAT_SurfDeform -iter 100 "%s" "%s" "%s"', Vppm.fname, P(si).Pcentral, P(si).Pcentral);
  cat_system(cmd,opt.verb-3);
  
  CS = loadSurf(P(si).Pcentral);
end
%==========================================================================
function [Vppm,rel] = exportPPmap( Yp0, Yp0fs, Yppi, Vmfs, Ypbs, optimize, si, opt, surffolder, ff)
%exportPPmap. Prepare image for resolution reduction.    
%
% Local optimization helps to keep thin structures. 
%
% [Vppm,Ypps] = exportPPmap( Yp0, Yp0fs, Yppi, Yth1i, Vmfs, vx_vol, si, opt, surffolder, ff)
%
% I tried here multiple things but most are not working or even contraproductive. 
%
% * Enhancement by sharpening as inverse smoothing: 
%   Was not robust and did not support good control.
%
% * Enhandement by local optimization to balance the sulcal-gyral span:
%   Partially working but only very small weightings are possible.
%   However, the values (distance+thickness) can help to control other 
%   modifications, such as the min-max-based enlargements of gyri/sulci. 
%
% * Divergence-based skeleton: 
%   Caused bridge defects but the information is useful to limit other things.
%
% * Intensity-based morphometry (gray-dilate) as minimum/maximum filter:
%   This works if the 



  % Yp0fs = cat_vol_median3(Yp0fs,Yppi>0 & Yppi<1);

  %% == optimized downsampling ==
  % extrem values:  .1^.3 = .5012, .6^1.3 = 0.5148
  rel  = max(.3,min(1.3, (nnz(Yp0(:)>2.5/3) / nnz(Yp0(:)>.5/3 & Yp0(:)<1.5/3)) ));
  switch optimize
    case 1 % local ... not working!
    %% balance suclus/gyrus thickness - make thin fat and fat thin - needs full resolution !
      Yid = cat_vbdist( single(Yppi < 0.45), Ypbs)/2  +  cat_vbdist( single(Yppi < 0.55), Ypbs)/2; Yid = Yid .* Ypbs; 
      Yod = cat_vbdist( single(Yppi > 0.45), Ypbs)/2  +  cat_vbdist( single(Yppi > 0.55), Ypbs)/2; Yod = Yod .* Ypbs; 
      % remove outliers
      Yit = cat_vol_pbtp( single(3 - (Yppi>.5)), Yid, 0*Yid); Yit(Yit > 10) = 0; 
      Yot = cat_vol_pbtp( single(2 + (Yppi>.5)), Yod, 0*Yod); Yot(Yot > 10) = 0; 
      % approximation 
      Yit = cat_vol_approx(Yit); 
      Yot = cat_vol_approx(Yot); 
      % "simple" correction:
      % - this has to be use very carefully and "Yppi.^max(.95,min(1.05, Yit ./ Yot ))" 
      %   presents already the stronges correction !!!
      Ypps = min(0,max( (Yp0fs)-2, Yppi)).^max(.3,min(1.3, cat_vol_smooth3X( Yit ./ Yot ,2 ) ) );
    case 2 % global ... is less good could be further improved
      Ypps = min(1,max(Yp0fs-2, Yppi.^max(.3,min(1.3,rel)) )); 
    case 3 % ########### need futher evaluation / test ! #########
      [Yppsr,Ypbsr,resYpp] = cat_vol_resize({Yppi * 2,single(Ypbs)},'reduceV',opt.interpV,4,32,'meanm');
      Ypppr = cat_vol_localstat(Yppsr,Ypbsr>.1,1,1,8); 
      Ypppr(Ypbsr>.5) = 1; 
      Ypppr = cat_vol_approx(Ypppr); 
      Yppp  = cat_vol_resize(Ypppr,'dereduceV',resYpp); 
      Yppp  = cat_vol_smooth3X(Yppp,4); 
      Ypps  = Yppi.^max(.5,min(1.0,Yppp)); % max( max(0,Yp0fs-2), Yppi).^max(.3,min(1.3,Yppp)); 
    otherwise
      Ypps = Yppi; 
  end
 
  %% Final downsampling:
  % - We combine here the average position with the GM/WM interface.
  fs = 0; %.35; 
  if opt.interpV == opt.reconres
    Vppm           = Vmfs;
    Ytx            = cat_vol_smooth3X( Ypps , fs );
  else
    [Vpp,Vppr]     = prepareDownsampling(Vmfs,Ypps,surffolder,ff,opt,si);
    Vpp.dat(:,:,:) = cat_vol_smooth3X( Ypps , fs ); 
    [Vppm,Yt]      = cat_vol_imcalc(Vpp,Vppr,'i1',struct('interp',5,'verb',0,'mask',-1));
    Vppm.pinfo     = Vmfs.pinfo;
    Ytx            = Yt; 
    Vppm
  end
 
  if isfield(Vppm,'dat'), Vppm = rmfield(Vppm,{'dat'}); end
  Vppm = spm_write_vol(Vppm, Ytx ); 
end
% ======================================================================
function pvet = estimatePVEsize( Yp0 , fast ) 
%estimatePVEsize. Estimate the voxel-size of partial volume label map.
%
%  pvet = estimatePVEsize( Yp0 [, fast]) 
% 
%  Yp0  .. label map
%  pvet .. partial volume effect size in voxel
%  fast .. use fast approximation (default) or full estimation (=0)

  if ~exist('fast','var'), fast = 1; end

  if fast
  % fast approximation
    Ypvet = cat_vbdist(single( Yp0<1.1 | (Yp0>1.9 & Yp0<2.1) | Yp0>2.9 )) * 2.5; 
    Ypvet(Ypvet>100 | Ypvet<0) = 0;
  else
  % accurate estimation 
    Ycgd  = cat_vbdist(single(Yp0<1.1), Yp0<1.9); Ycgd(Ycgd>100 | Ycgd<0) = 0; 
    Ygcd  = cat_vbdist(single(Yp0>1.9), Yp0>1.1); Ygcd(Ygcd>100 | Ygcd<0) = 0; 
    Ygwd  = cat_vbdist(single(Yp0<2.1), Yp0<2.9); Ygwd(Ygwd>100 | Ygwd<0) = 0; 
    Ywgd  = cat_vbdist(single(Yp0>2.9), Yp0>2.1); Ywgd(Ywgd>100 | Ywgd<0) = 0; 

    Ypvet = max(Ycgd + Ygcd, Ygwd + Ywgd); 
  end

  % final evaluation
  pvet  = max(1,median(Ypvet(Ypvet(:)>0) - 1)); 
end
% ======================================================================
function Yp0 = myelincorrection(Yp0,vx_vol,opt,P,Vmfs,si)
    % quick estimation of the cortical thickness
    
    if 1 
    % quicker but sufficient version (~10s) 
      % parameters 
      opt.verb      = 0; 
      opt.levels    = 2; 
     
      opt.pvet0 = estimatePVEsize( Yp0 , 0); 
      opt.pvet  = max(mean(vx_vol),min(2 / mean(vx_vol),opt.pvet0)); 

      % from cat-vol_pbtsimpleCS4
      [Ycd, Ywd] = cat_vol_cwdist(Yp0, opt);
    
      % projection-based thickness mapping
      Ygmt0 = cat_vol_pbtp( round(Yp0) , Ywd, Ycd);
      Ygmt0 = cat_vol_approx(Ygmt0,'rec',2); 
    
      % reestimation of the CSF distance 
      Ypp   = min(1,min(Ygmt0,Ycd) ./ max(eps,Ygmt0)); Ypp(Yp0>2.5 & Ypp==0) = 1; 
    else
      % alternative save but slow version (45s)
      cmd = sprintf('CAT_VolThicknessPbt  -correct-voxelsize 0   -median-filter 2   -downsample 0 "%s" "%s" "%s"', Vmfs.fname, P(si).Pgmt, P(si).Pppm);
      cat_system(cmd,opt.verb-3);
      Vgmt0  = spm_vol(P(si).Pgmt); Ygmt0 = spm_read_vols(Vgmt0); 
      Vpp    = spm_vol(P(si).Pppm); Ypp  = spm_read_vols(Vpp); 
      Ycd    = Ygmt0 .* Ypp; 
    end
  
    % reestimation of the CSF distance 
    Ycdc2 = cat_vbdist( single( max(Yp0<=1, 1 - Ycd - Ypp) ), true(size(Ycd)) );
    Ycdc2(Ycdc2 > 6 / mean(vx_vol)) = 0; 
   
    % estimate the full tissue thickness (we needed the GM thickness and WM to reconstruct the sulcus)
    Ybmt  = cat_vol_pbtp( min(3,4 - min(2,Yp0)), Ycdc2, Ycdc2*inf); 
    Ybmt  = cat_vol_approx(Ybmt); 

    % estimate correction area
    medgmt = median(Ygmt0(:)); 
    try iqrgmt = iqr(Ygmt0(:)); catch, iqrgmt = std(Ygmt0(:)); end
    YenoughWM  = Ycdc2 < Ybmt - 1.5;
    YthinnerGM = max(0,medgmt - 1.5*iqrgmt - Ygmt0); 
    Yclose2CSF = Ycdc2>0 & Ycdc2<(medgmt - 1.5*iqrgmt); 
    Ygmwmpve   = cat_vol_morph(Yp0>2 & Yp0<2.9,'do',1); % | smooth3(Yp0>2 & Yp0<2.9)>.7);
    Ycor = YenoughWM & YthinnerGM & Yclose2CSF & Ygmwmpve;
    clear YenoughWM YthinnerGM Yclose2CSF Ygmwmpve;

    Yp0  = max(min(Yp0,2),max(Yp0>=2.95,Yp0 - smooth3( Ycor ) * opt.myelinCorrection));
    clear Ycdc2 Ybmt Ycor; 

end
% ======================================================================
function [Ycd, Ywd] = cat_vol_cwdist(Yp0,opt)
%cat_vol_cwdist. Estimation of CSF and WM distance in a label map Yp0.
% 
% [Ycd, Ywd] = cat_vol_cwdist(Yp0,opt)
%
% Ycd, Ywd         .. CSF and WM distance maps
% opt              .. parameter structure
%  .levels         .. number of dual distance measurements
%  .range          .. limitation to avoid bias by interpolation overshoot 
%
    
  def.levels = 4; 
  def.rangeE = 0.05; % extimation limit extention - smaller better
  def.range  = 0.40; % addition range - larger better but limited by interpolation artifacts
  def.vxs    = 1;

  opt = cat_io_checkinopt(opt,def); 
  vxc = .5 - opt.range; % the more PVE the less we need to correct for

  range  = min(.45,opt.range);    % defined by the PVE boundary of the AMAP, i.e range of .25 - .75  
  rangeE = min(.45,opt.rangeE);    % defined by the PVE boundary of the AMAP, i.e range of .25 - .75  
  
  % The idea is that the we use here the full range of of the 1.5 and 2.5 
  % AMAP class to define the full thickness. However, we measure still
  % from the .5er boundary and we have to handle the 
  YMM = Yp0 >= 1.5 - rangeE | isnan(Yp0); % range for WMD
  YMC = Yp0 <= 2.5 + rangeE | isnan(Yp0); % range for CSFD
  
  % multi-level distance estimation
  Ycd = zeros(size(Yp0),'single'); 
  Ywd = zeros(size(Yp0),'single'); 
  hss = opt.levels; % number of opt.levels (as pairs)
  for si = 1:hss
    offset = max(0,min(.4, range * si/(hss+1))); 

    % CSF dist
    [Ycdl,YI]  = cat_vbdist(single(Yp0 < ( 1.5 - offset)), YMC ); 
    Ycdl       = (Ycdl - min(vxc,max(0,- (Yp0(YI) - (1.5-offset))*opt.pvet))) .* (YMM & YMC); 
    [Ycdh,YI]  = cat_vbdist(single(Yp0 < ( 1.5 + offset)), YMC ); 
    Ycdh       = (Ycdh - min(vxc,max(0,- (Yp0(YI) - (1.5+offset))*opt.pvet))) .* (YMM & YMC); 
    Ycd        = Ycd + .5/hss .* Ycdl  +  .5/hss .* Ycdh; 

    % WM distances
    [Ywdl,YI]  = cat_vbdist(single(Yp0 > ( 2.5 - offset)), YMM );
    Ywdl       = (Ywdl - min(vxc,max(0, (Yp0(YI) - (2.5-offset))*opt.pvet))) .* (YMM & YMC); 
    [Ywdh,YI]  = cat_vbdist(single(Yp0 > ( 2.5 + offset)), YMM ); 
    Ywdh       = (Ywdh - min(vxc,max(0, (Yp0(YI) - (2.5+offset))*opt.pvet))) .* (YMM & YMC); 
    Ywd        = Ywd + .5/hss .* Ywdl  +  .5/hss .* Ywdh;
  end

  % endpoint PVE correction 
  Ycd = Ycd - min(.5,max(-.0, (Yp0>2.5-0*opt.rangeE & Yp0<2.5+opt.rangeE) .* ( (Yp0-2.5)*opt.pvet*1 + 0*opt.vxs/2) ));
  Ywd = Ywd + min(.0,max(-.5, (Yp0>1.5-opt.rangeE & Yp0<1.5+0*opt.rangeE) .* ( (Yp0-1.5)*opt.pvet*1 - 0*opt.vxs/2) ));
end
% ======================================================================

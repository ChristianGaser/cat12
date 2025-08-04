function [Yth,S,P,res] = cat_surf_createCS4(V,V0,Ym,Yp0,Ya,YMF,Yb0,opt,job)
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
  use_cat_vol_pbtsimple   = 1;
  skip_registration       = isfield(opt,'surf') && isscalar(opt.surf); % skip spherical registration for quick tests
  create_white_pial       = 1; % uses only the quick WM and Pial surface estimation 

  myelinCorrection        = .2; % .25 - sight correction, 1 - maximum correction
  setcut2zero             = 0; % works but result in worse values as could be expected
  
  % surface output and evaluation parameter 
  res   = struct('lh',struct(),'rh',struct()); 
  Yth   = zeros(size(Yp0),'single');  % initialize WM/CSF thickness/width/depth maps
  S     = struct();
  

  % set defaults
  % set debugging variable
  dbs   = dbstatus; debug = 1; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end
  vx_vol                  = sqrt(sum(V.mat(1:3,1:3).^2));               % further interpolation based on internal resolution 
  def.verb                = cat_get_defaults('extopts.expertgui');      % 0-none, 1-minimal, 2-default, 3-details, 4-debug
  def.surf                = {'lh','rh'};                                % surface reconstruction setting with {'lh','rh','cb'} 
  % There is a new SPM approach spm_mesh_reduce that is maybe more robust. 
  % Higher resolution is at least required for animal preprocessing that is given by cat_main.
  def.LAB                 = cat_get_defaults('extopts.LAB');  % brain regions 
  % RD20250306: Tfs has large issues currently with some corrected defects
  def.useprior            = ''; 
  def.thick_limit         = 6;                                % 6mm upper limit for thickness (similar limit as used in Freesurfer)
  def.foldingcorrection   = 1;                                % tickness correction that is influence by folding
  def.thick_measure       = 0;                                % 0-PBT; 1-Tfs (Freesurfer method using mean(TnearIS,TnearOS)) ##########
  def.fsavgDir            = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces'); 
  def.outputpp.native     = 0;  % output of Ypp map for cortical orientation in EEG/MEG 
  def.outputpp.warped     = 0;
  def.outputpp.dartel     = 0;
  def.vdist               = 2; 
  
  % options that rely on other options
  opt                     = cat_io_updateStruct(def,opt); clear def; 
  opt.vol                 = any(~cellfun('isempty',strfind(opt.surf,'v')));   % only volume-based thickness estimation  
  opt.surf                = cat_io_strrep(opt.surf,'v','');                   % after definition of the 'vol' varialbe we simplify 'surf'
  opt.interpV             = max(0.1,min([opt.interpV,1.5]));                  % general limitation of the PBT resolution
  switch opt.SRP
    case 0 % pbtsimpleC, full resolution
      use_cat_vol_pbtsimple = 0; 
      myelinCorrection      = 0; 
      opt.reconres          = .5; 
    case 1 % myelincorrection 
      use_cat_vol_pbtsimple = 1; 
      myelinCorrection      = 0.3; 
      opt.reconres          = .6; 
    case 2 % default
      use_cat_vol_pbtsimple = 1; 
      myelinCorrection      = 0; 
      opt.reconres          = .6; 
  end


  % apply the modified mask from gcut
  % for non-gcut approaches or inverse weighting Yb0 only contains ones
  Yp0   = Yp0 .* (Yb0>0.5);
  Ym    = Ym  .* (Yb0>0.5);


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
    if 1 
      % sharpening of raw image
      Ym = Ym + (Ym - smooth3(Ym)) / 2; 
      %Ym(Yp0>2 & Yp0<3) = Ymx(Yp0>2 & Yp0<3); 
      %Ym(Yp0>1 & Yp0<2) = Ymx(Yp0>1 & Yp0<2); 
    end
 
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
  [P,pp0,mrifolder,surffolder,surfdir,ff] = setFileNames(V0,job,opt); 
  

  % main loop for each surface structure 
  for si = 1:numel(opt.surf)
   
    % prepare longitudinal case if required 
    useprior = setupprior(opt,surfdir,P,si);


    %% reduce for object area
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
    % print something
    if si==1, fprintf('\n'); end; fprintf('%s:\n',opt.surf{si});
    clear Ynocerebrum

    
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
    [Vmfs,Smat] = createOutputFileStructures(V,V0,resI,BB,opt,pp0,mrifolder,ff,si); 


    %% thickness and position map estimation 
    stime = cat_io_cmd(sprintf('  Thickness estimation (%0.2f mm%s)',opt.interpV,native2unicode(179, 'latin1')),'g5'); stimet = stime; % fprintf('\n'); 
    if use_cat_vol_pbtsimple
      if all(opt.interpV ~= vx_vol) 
        %% reduction of interpolation artifacts 
        %  - remove (i) slight over/underestimation close to boundaries and 
        %    (ii) larger voxel steps by utilizing the median filter but  
        %    try to prevent blurring of thin sulci/gyri
        %  - allow slighly smoother surfaces (better position with similar intensity) 
        Yfr   = cat_vol_morph( round(Yp0fs)==2 ,'d',2);
        Yp0fs = max( 1 , min( 3,  ...              % general limits
          min( Yp0fs .* (3 - 2*(Yp0fs<1.5)) , ...  % protect CSF
           max( Yp0fs .* (Yp0fs>2.5) , ...         % protect WM 
           cat_vol_median3(Yp0fs,Yfr) ))));        % cleanup otherwise
        % sharpening
        Yp0fs = max( 1 , min( 3, Yp0fs + .5 * (round(Yp0fs)~=2) .* (Yp0fs - smooth3(Yp0fs)) )); 
      end

      % WM geometry and topology correction 
      % - this might lower thickness estimates and it might close small sulci
      % - takes about 1 minute so I would try to avoid it
      % - in case of amap it can be run for 2.25 and 2.75
      % ... we will try to avoid it as is can also create new problems and take 30s 
      if 1 
        Yp0fs = (cat_vol_morph(Yp0fs*2 - 5,'wmtc',1,1,0) + 5)/2; % correct at 2.75
        %Yp0fs = (cat_vol_morph(Yp0fs*2 - 4,'wmtc',1,1,0) + 4)/2; % correct at 2.25 ... can close gyri better avoid
      end      
      
      % thickness and position estimation
      [Yth1i,Yppi,Ymfsc] = cat_vol_pbtsimpleCS4(Yp0fs, opt.interpV,struct('myelinCorrection',myelinCorrection,'verb',1,'gyrusrecon',1));
      
    else
      %% Write PP
      Vmfs.dt = [16 1];
      spm_write_vol(Vmfs, Yp0fs);
      cmd = sprintf('CAT_VolThicknessPbt -median-filter 2 -sharpen 0.02 -fwhm 3 -downsample 0 "%s" "%s" "%s"', Vmfs.fname, P(si).Pgmt, P(si).Pppm);
      cat_system(cmd,opt.verb-3);
      Vgmt = spm_vol(P(si).Pgmt); Yth1i = spm_read_vols(Vgmt); 
      Vppi = spm_vol(P(si).Pppm); Yppi  = spm_read_vols(Vppi); 
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

    %% RD20250330: The aim is to smoothing cutting regions to avoid defects there. 
    if 0;; % just for manual tests
      opt.reconres    = .6;  % although 1.5-2.0 mm can work in general, defects correction can result in very severe and unpredictable problems
      opt.vdist       = 2;   % controls surface reduction (linear value >> sqare for area, (8~2 mm), 4~1.2 mm, 2~1 mm,  1~0.7, .5~0.5)
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
      reconattempts = round((1.2 - .7)*10); % for reduction of .1 mm per step with 1.2 mm in the worst case
      for li = 1:reconattempts
        if li > 1, opt.reconres = opt.reconres + 0.1; end

        % optimized downsampling of the the Ypp map and 
        if isscalar(opt.surf), time_sr = clock; end % temporary for tests 
        [Vppmi,rel] = exportPPmap( Yp0 .* Ypb, Ymfsc, Yppi, Vmfs, Ypbs, vx_vol, si, opt, surffolder, ff);
        
        % Main initial surface creation with topology correction.
        stime = cat_io_cmd(sprintf('  Create initial surface (%0.2f mm)',opt.reconres),'g5','',opt.verb,stime); %if opt.verb>2, fprintf('\n'); end
        if exist(P(si).Pcentral,'file'), delete(P(si).Pcentral); end % have to delete it to get useful error messages in case of reprocessing/testing
        cmd  = sprintf('CAT_VolMarchingCubes "%s" "%s" -thresh "%0.4f" ', Vppmi.fname, P(si).Pcentral, .55); %5-0.05*opt.reconres);
        Vpp = spm_write_vol(Vppm, Yppi);
        cat_system(cmd ,opt.verb-3);
  
  
        % load and check surface
        % use 80000 as lower limit what is similar to vdist=4
        CS  = loadSurf(P(si).Pcentral); 
        A   = sum( cat_surf_fun( 'area', CS )); % in mm2 
        nCS = min( size(CS.faces,1) , max( 80000, min( 360000, round(A * 2 * 4 / opt.vdist^2) ))); 
        if opt.SRP >= 2, nCS = min(300000,nCS); end
        EC0 = size(CS.vertices,1) + size(CS.faces,1) - size(spm_mesh_edges(CS),1);  
        if EC0 ~= 2
          if li < reconattempts
            % try again with another resolution
            continue
          else
            cat_io_addwarning('cat_surf_createCS4:TopologyWarning', ...
             sprintf('Extracted surface might have small topology issues (Euler count = %d).',EC0),0,[0 1]); 
          end
        end      
      
        % surface reduction 
        % - low surface resolution are highly important for fast processing and smaller disk usage
        % - less critical for thickness but surface quality (intensity/position values) get worse
        % - it is important to keep the topology 
        % - a more aggressive 
        stimesr = clock; 
        nCS0 = size(CS.faces,1); 
        saveSurf(CS,P(si).Pcentral); 
        cmd  = sprintf('CAT_SurfReduce -aggr "5" -ratio "%0.2f" "%s" "%s"', nCS/nCS0 , P(si).Pcentral, P(si).Pcentral); %5-0.05*opt.reconres);
        cat_system(cmd ,opt.verb-3);
        CS = loadSurf(P(si).Pcentral); 
        nCS1 = size(CS.faces,1); 
        saveSurf(CS,P(si).Pcentral); 
        
        if li == reconattempts  ||  nCS1 < nCS * 1.5 
          % all done
          cat_io_cmd(sprintf('  Reduce surface (nF: %0.0fk > ~%0.0fk)', nCS0/1000, nCS1/1000),'g5','',opt.verb,stime);
          break
        elseif li < reconattempts  
          % if the surface canot be reduced then try again with lower resoution
          continue
        elseif li == reconattempts && nCS1 > nCS * 1.5 
          cat_io_addwarning('cat_surf_createCS4:ReductionFailed', ...
            sprintf('Surface reduction incomplete due to topology constrains (%d).',(nCS1/1000-nCS/1000) ./ (nCS0/1000-nCS/1000)),0,[0 1]); 
        end 
        stime = stimesr; 
      end
    end


    % remove artifacts 
    % - this can maybe partially avoided by the more aggressive surface reduction
    % - although the filter reduces some artefacts in problematic cases the surface values are worse and the issue are still not fully solved
    % - tried to use it after correction of self-intersections but this made it even worse
    % >> better not use now
    if 0 
      stime = cat_io_cmd('  Filter topology artifacts','g5','',opt.verb,stime);
      CS = smoothArt(Yth1i,P,CS,Smat,Vppm,1,si,opt,1,1); % last options - method & refine
    end

    % get PBT thickness and update cutregion
    if setcut2zero
      facevertexcdatanocut = 1 - cat_surf_fun('isocolors',Ycutregion,CS.vertices,Smat.matlabIBB_mm); 
      facevertexcdatanocut = max(0.5,min(1,(cat_surf_fun('smoothcdata',CS,facevertexcdatanocut,4) - .5) * 2)); 
    else
      facevertexcdatanocut = 1; 
    end
    facevertexcdata = cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm) .* facevertexcdatanocut; 
    cat_io_FreeSurfer('write_surf_data',P(si).Ppbt,facevertexcdata);
    cat_io_FreeSurfer('write_surf_data',P(si).Pthick,facevertexcdata);


    % surface deformation  
    stime = cat_io_cmd(sprintf('  Optimize surface (nF: %0.0fk)', size(CS.faces,1)/1000),'g5','',opt.verb,stime); 
    for li = 1:2 % ... interation 1 with additional correction, iteration 2 without ... only small improvements for multiple iterations
      cmd = sprintf('CAT_SurfDeform -iter 100 "%s" "%s" "%s"', Vppm.fname, P(si).Pcentral, P(si).Pcentral);
      cat_system(cmd,opt.verb-3);
      CS = loadSurf(P(si).Pcentral);
 
      if opt.SRP > 0 && li == 1
        for xi = 1:2 % interation 1 with optimization, 2 without .. 
          % quick correction with local optimization 
          [CS,facevertexcdatac,SIs] = cat_surf_fun('collisionCorrectionPBT',CS,facevertexcdata .* facevertexcdatanocut,Ymfsc,Yppi, ... 
            struct('optimize',xi==1,'verb',isscalar(opt.surf)*2,'mat',Smat.matlabIBB_mm,'vx_vol',vx_vol,'CS4',0)); % CS4 is not working
          if isscalar(opt.surf), fprintf('\b\b'); end
          % accurate but slow function (especially for multiple jobs) 
          if xi == 1 %&& ~isscalar(opt.surf) 
            [CS,facevertexcdatac,SIs] = cat_surf_fun('collisionCorrectionRY' ,CS,facevertexcdatac .* facevertexcdatanocut,Yppi,...
              struct('Pcs',P(si).Pcentral,'verb',isscalar(opt.surf)*2,'mat',Smat.matlabIBB_mm,'accuracy',1/2^2));
            if isscalar(opt.surf), fprintf('\b\b'); end
          end
        end
        facevertexcdata = max(facevertexcdata,facevertexcdatac) .* facevertexcdatanocut; 
        cat_io_FreeSurfer('write_surf_data',P(si).Ppbt,facevertexcdata);
        saveSurf(CS,P(si).Pcentral); 
      end
    end
   

    % evaluate and save results
    if isempty(stime), stime = clock; end
    fprintf('%5.0fs',etime(clock,stime)); stime = []; if 1, fprintf('\n'); end %debug


    % create white and central surfaces
    if create_white_pial
      stime = cat_io_cmd('  Create pial and white surface','g5','',opt.verb,stime); 
      cmd = sprintf('CAT_Central2Pial "%s" "%s" "%s" 0.5',P(si).Pcentral,P(si).Ppbt,P(si).Ppial);
      cat_system(cmd,opt.verb-3);
      cmd = sprintf('CAT_Central2Pial "%s" "%s" "%s" -0.5',P(si).Pcentral,P(si).Ppbt,P(si).Pwhite);
      cat_system(cmd,opt.verb-3); 
    end


    if isscalar(opt.surf)
    %% ========= only for test visualization ========= 
      if 1 % opt.thick_measure == 1 % allways here 
        cmd = sprintf('CAT_SurfDistance -mean -thickness "%s" "%s" "%s"',P(si).Ppbt,P(si).Pcentral,P(si).Pthick);
        cat_system(cmd,opt.verb-3);
        % apply upper thickness limit
        facevertexcdata = cat_io_FreeSurfer('read_surf_data',P(si).Pthick) .* facevertexcdatanocut;  
        facevertexcdata(facevertexcdata > opt.thick_limit) = opt.thick_limit;
        cat_io_FreeSurfer('write_surf_data',P(si).Pthick,facevertexcdata);  
      end
 
      fprintf('\n');
      res.(opt.surf{si}).createCS_final = cat_surf_fun('evalCS', ...
        loadSurf(P(si).Pcentral), cat_io_FreeSurfer('read_surf_data',P(si).Ppbt), facevertexcdata, ...
        Ymfs,Yppi,P(si).Pcentral,Smat.matlabIBB_mm,2,0);
      CS2 = CS; CS2.cdata = facevertexcdata; H = cat_surf_render2(CS2);
      cat_surf_render2('clim',H,[0 6]); 
      cat_surf_render2('view',H,cat_io_strrep(opt.surf{si},{'lh','rh','ch'},{'right','left','back'})); 
      cat_surf_render2('ColourBar',H,'on');
      title(sprintf('CS4%d, R=%0.2f, nF=%0.0fk, EC=%d, \n TH=%0.3f±%0.3f, IE=%0.3f, PE=%0.3f, ptime=%0.0fs, time=%s', ...
        opt.SRP, opt.reconres, size(CS.faces,1)/1000, EC0, ...
        mean( facevertexcdata ), std(facevertexcdata), ...
        mean( [ res.(opt.surf{si}).createCS_final.RMSE_Ym_white,  res.(opt.surf{si}).createCS_final.RMSE_Ym_layer4,   res.(opt.surf{si}).createCS_final.RMSE_Ym_pial ] ) , ...
        mean( [ res.(opt.surf{si}).createCS_final.RMSE_Ypp_white, res.(opt.surf{si}).createCS_final.RMSE_Ypp_central, res.(opt.surf{si}).createCS_final.RMSE_Ypp_pial ] ) , ...
        etime(clock,time_sr), datetime))
      subtitle( strrep( spm_str_manip(P(si).Pcentral,'a90') ,'_','\_'))
      fprintf('    Runtime:                             %0.0fs\n',etime(clock,time_sr)); 


      % preparation for SPM orthviews links
      Po = P(si).Pm; if ~exist(Po,'file'); Po = V0.fname; end
      if ~exist(Po,'file')  && exist([V0.fname '.gz'],'file'), Po = [V0.fname '.gz']; end
      Porthfiles = ['{', sprintf('''%s'',''%s''',P(si).Ppial, P(si).Pwhite ) '}'];
      Porthcolor = '{''-g'',''-r''}';  
      Porthnames = '{''white'',''pial''}'; 
      fprintf('  Show surfaces in orthview:   %s | %s | %s | (%s) | %s \n', ...
        spm_file([opt.surf{si} '.pbt'],'link', sprintf('H=cat_surf_display(''%s'');',P(si).Ppbt)), ...
        spm_file([opt.surf{si} '.thick'],'link', sprintf('H=cat_surf_display(''%s'');',P(si).Pthick)), ...
        spm_file('segmentation' ,'link', sprintf('cat_surf_fun(''show_orthview'',%s,''%s'',%s,%s)',Porthfiles,P(si).Pp0, Porthcolor,Porthnames)), ...
        spm_file('ppmap'        ,'link', sprintf('cat_surf_fun(''show_orthview'',%s,''%s'',%s,%s)',Porthfiles,Vpp.fname, Porthcolor,Porthnames)), ...
        spm_file('original'     ,'link', sprintf('cat_surf_fun(''show_orthview'',%s,''%s'',%s,%s)',Porthfiles,Po,        Porthcolor,Porthnames)));
      
      subtitle( strrep( spm_str_manip( P(si).Pcentral,'a90') ,'_','\_'))
      fprintf('    Runtime:                             %0.0fs\n',etime(clock,time_sr)); 
   
      return
    end



%% == TEST STOP ==


    % skip that part if a prior image is defined
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
      % estimate Freesurfer thickness measure Tfs using mean(Tnear1,Tnear2)
      if 0
        % estimate Freesurfer thickness based on WM and Pial surface ... not ready yet
        cmd = sprintf('CAT_SurfDistance -mean "%s" "%s" "%s"',P(si).Pwhite,P(si).Ppial,P(si).Pthick);
        cat_system(cmd,opt.verb-3);
      else 
        % use central surface and thickness to estimate Freesurfer thickness
        cmd = sprintf('CAT_SurfDistance -mean -thickness "%s" "%s" "%s"',P(si).Ppbt,P(si).Pcentral,P(si).Pthick);
        cat_system(cmd,opt.verb-3);
      end
      
      % apply upper thickness limit
      % here the 5 mm thickness limit of FreeSurfer might be better rather then our 6 mm 
      facevertexcdata = cat_io_FreeSurfer('read_surf_data',P(si).Pthick) .* facevertexcdatanocut;  
      facevertexcdata(facevertexcdata > opt.thick_limit) = opt.thick_limit; 
      cat_io_FreeSurfer('write_surf_data',P(si).Pthick,facevertexcdata);  
    else
      % otherwise simply use the original values of the PBT map  
      % WARNING: 
      %   The values of the ?h.pbt.* files are used to estimate further 
      %   surfaces and therefore corrected for self-intersections!
      %   We use here the orignal PBT values as they are less depending on 
      %   local highly individual features.  
      facevertexcdata = cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm) .* facevertexcdatanocut; 
      cat_io_FreeSurfer('write_surf_data',P(si).Ppbt,facevertexcdata);
    end
    fprintf('\n');
    

    % correct thickness based on folding pattern 
    if opt.foldingcorrection
      cmd = sprintf('CAT_SurfCorrectThicknessFolding -max "%f" "%s" "%s" "%s"',opt.thick_limit,P(si).Pcentral,P(si).Ppbt,P(si).Pthick);
      cat_system(cmd,opt.verb-3);
    else
      copyfile(P(si).Ppbt,P(si).Pthick,'f'); 
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
    
    if exist(Vppm.fname ,'file'), delete(Vppm.fname); end
    if debug && exist(Vpp.fname ,'file') && ~opt.outputpp.native, delete(Vpp.fname); end
  
    % processing time per side for manual tests
    if si == numel(opt.surf) && si == 1
      cat_io_cmd('  ','g5','',opt.verb);
      fprintf('%5ds\n',round(etime(clock,cstime)));
    end
  end  
  

  % calculate surface quality parameters for all surfaces
  res = addSurfaceQualityMeasures(res,opt);
 

  % print final stats
  evalProcessing(res,opt,P,V0,si)

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
function [Vmfs,Smat] = createOutputFileStructures(V,V0,resI,BB,opt,pp0,mrifolder,ff,si)
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
    Vmfs.fname = fullfile(pp0,mrifolder, sprintf('%s_seg-%s.nii',ff,opt.surf{si}));
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
function res = addSurfaceQualityMeasures(res,opt)
%addSurfaceQualityMeasures. Measures to describe surface properties. 
  res.mnth = []; res.sdth = []; 
  res.mnRMSE_Ypp = []; res.mnRMSE_Ym = []; 
  res.SIw = []; res.SIp = []; res.SIwa = []; res.SIpa = []; 
  for si=1:numel(opt.surf)
    if any(strcmp(opt.surf{si},{'lh','rh'}))
      if isfield(res.(opt.surf{si}).createCS_final,'fsthickness_mn_sd_md_mx') && ... 
        ~isnan( res.(opt.surf{si}).createCS_final.fsthickness_mn_sd_md_mx(1) )
        res.mnth      = [ res.mnth  res.(opt.surf{si}).createCS_final.fsthickness_mn_sd_md_mx(1) ]; 
        res.sdth      = [ res.sdth  res.(opt.surf{si}).createCS_final.fsthickness_mn_sd_md_mx(2) ]; 
      else
        res.mnth      = [ res.mnth  res.(opt.surf{si}).createCS_final.thickness_mn_sd_md_mx(1) ];
        res.sdth      = [ res.sdth  res.(opt.surf{si}).createCS_final.thickness_mn_sd_md_mx(2) ]; 
      end
      res.mnRMSE_Ym   = [ res.mnRMSE_Ym   mean([...
        res.(opt.surf{si}).createCS_final.RMSE_Ym_layer4 ...
        res.(opt.surf{si}).createCS_final.RMSE_Ym_white ...
        res.(opt.surf{si}).createCS_final.RMSE_Ym_pial ]) ];
      res.mnRMSE_Ypp  = [ res.mnRMSE_Ypp  mean([...
        res.(opt.surf{si}).createCS_final.RMSE_Ypp_central ...
        res.(opt.surf{si}).createCS_final.RMSE_Ypp_white ...
        res.(opt.surf{si}).createCS_final.RMSE_Ypp_pial ]) ];
      if isfield(res.(opt.surf{si}).createCS_final,'white_self_interections')
        res.SIw     = [ res.SIw  res.(opt.surf{si}).createCS_final.white_self_interections ]; 
        res.SIp     = [ res.SIp  res.(opt.surf{si}).createCS_final.pial_self_interections  ]; 
        res.SIwa    = [ res.SIwa res.(opt.surf{si}).createCS_final.white_self_interection_area ]; 
        res.SIpa    = [ res.SIpa res.(opt.surf{si}).createCS_final.pial_self_interection_area  ]; 
      end
    end
  end
  
  % final res structure
  res.EC          = NaN; 
  res.defect_size = NaN;
  res.defect_area = NaN;
  res.defects     = NaN;
  res.mnth        = mean(res.mnth); 
  res.sdth        = mean(res.sdth); 
  res.RMSE_Ym     = mean(res.mnRMSE_Ym);
  res.RMSE_Ypp    = mean(res.mnRMSE_Ypp);
end
%==========================================================================
function evalProcessing(res,opt,P,V0,si)
  if opt.verb && ~opt.vol  
  % display some evaluation 
  % - For normal use we limited the surface measures.  
  % - Surface intensity would be interesting as cortical measure similar to thickness (also age dependent).
  %   Especially the outer surface will describe the sulcal blurring in children. 
  %   But the mixing of surface quality and anatomical features is problematic. 
  % - The position value describes how good the transformation of the PBT map into a surface worked. 
  %   Also the position values depend on age. Children have worse pial values due to sulcal blurring but
  %   the white surface is may effected by aging, e.g., by WMHs.
  % - However, for both intensity and position some (average) maps would be also interesting. 
  %   Especially, some Kappa similar measure that describes the differences to the Ym or Ypp would be nice.
  % - What does the Euler characteristic say?  Wouldn't the defect number more useful for users? 

    if any(~cellfun('isempty',strfind(opt.surf,'cb'))), cbtxt = 'cerebral '; else, cbtxt = ''; end
    fprintf('Final %ssurface processing results: \n', cbtxt);
      
    % function to estimate the number of interactions of the surface deformation: d=distance in mm and a=accuracy 
    QMC    = cat_io_colormaps('marks+',17);
    color  = @(m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
    rate   = @(x,best,worst) min(6,max(1, max(0,x-best) ./ (worst-best) * 5 + 1));
  
    if cat_get_defaults('extopts.expertgui')
    % color output currently only for expert ...
      if isfield(res.(opt.surf{si}).createCS_final,'fsthickness_mn_sd_md_mx')
        fprintf('  Average thickness (FS):                     ');
      else
        fprintf('  Average thickness (PBT):                    ');
      end
      cat_io_cprintf( color( rate( abs( res.mnth - 2.5 ) , 0 , 2.0 )) , sprintf('%0.4f'     , res.mnth ) );  fprintf(' %s ',native2unicode(177, 'latin1'));
      cat_io_cprintf( color( rate( abs( res.sdth - 0.5 ) , 0 , 1.0 )) , sprintf('%0.4f mm\n', res.sdth ) );

      fprintf('  Surface intensity / position RMSE:          ');
      cat_io_cprintf( color( rate( mean(res.mnRMSE_Ym)  , 0.05 , 0.3 ) ) , sprintf('%0.4f / ', mean(res.mnRMSE_Ym) ) );
      cat_io_cprintf( color( rate( mean(res.mnRMSE_Ypp) , 0.05 , 0.3 ) ) , sprintf('%0.4f\n',  mean(res.mnRMSE_Ypp) ) );
        
      if isfield(res.(opt.surf{si}).createCS_final,'white_self_interections')
        fprintf('  Pial/white self-intersections:              ');
        cat_io_cprintf( color( rate(  mean([res.SIw,res.SIp]) , 0 , 20 ) ) , sprintf('%0.2f%%%% (%0.2f mm%s)\n'  , mean([res.SIw,res.SIp]) , mean([res.SIwa,res.SIpa]) , char(178) ) );
      end
    else
      fprintf('  Average thickness:                          %0.4f %s %0.4f mm\n' , res.mnth, native2unicode(177, 'latin1'), res.sdth);
    end
    
    for si=1:numel(P)
      fprintf('  Display thickness:          %s\n',spm_file(P(si).Pthick,'link','cat_surf_display(''%s'')'));
    end
    
    %% surfaces in spm_orthview
    if exist(P(si).Pm,'file'), Po = P(si).Pm; else, Po = V0.fname; end
    if ~exist(Po,'file') && exist([V0.fname '.gz'],'file'), Po = [V0.fname '.gz']; end
    
    Porthfiles = '{'; Porthcolor = '{'; Porthnames = '{';
    for si=1:numel(P)
      Porthfiles = [ Porthfiles , sprintf('''%s'',''%s'',',P(si).Ppial, P(si).Pwhite )]; 
      Porthcolor = [ Porthcolor , '''-g'',''-r'',' ]; 
      Porthnames = [ Porthnames , '''pial'',''white'',' ];
    end
    Porthfiles = [ Porthfiles(1:end-1) '}']; 
    Porthcolor = [ Porthcolor(1:end-1) '}']; 
    Porthnames = [ Porthnames(1:end-1) '}']; 
  
    if 1 %debug 
      fprintf('  Show surfaces in orthview:  %s\n',spm_file(Po ,'link',...
        sprintf('cat_surf_fun(''show_orthview'',%s,''%s'',%s,%s)',Porthfiles,Po,Porthcolor,Porthnames))) ;
    end

  end
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
function [P,pp0,mrifolder,pp0_surffolder,surffolder,ff] = setFileNames(V0,job,opt) 
  
  [mrifolder, ~, surffolder] = cat_io_subfolders(V0.fname,job);
  
  % get original filename without 'n'
  [pp0,ff] = spm_fileparts(V0.fname);
  
  % correct '../' parts in directory for BIDS structure
  [stat, val] = fileattrib(fullfile(pp0,surffolder));
  if stat, pp0_surffolder = val.Name; else, pp0_surffolder = fullfile(pp0,surffolder); end
  if ~exist(fullfile(pp0_surffolder),'dir'), mkdir(fullfile(pp0_surffolder)); end

  % surface filenames
  for si = 1:numel(opt.surf)
    P(si).Pm         = fullfile(pp0,mrifolder, sprintf('m%s.nii',ff));                                 % raw
    P(si).Pp0        = fullfile(pp0,mrifolder, sprintf('p0%s.nii',ff));                                % labelmap
    P(si).Pdefects   = fullfile(pp0,mrifolder, sprintf('defects_%s.nii',ff));                          % defect
    P(si).Pcentral   = fullfile(pp0_surffolder,sprintf('%s.central.%s.gii',opt.surf{si},ff));          % central
    P(si).Pcentralh  = fullfile(pp0_surffolder,sprintf('%s.centralh.%s.gii',opt.surf{si},ff));         % central
    P(si).Pcentralr  = fullfile(pp0_surffolder,sprintf('%s.central.resampled.%s.gii',opt.surf{si},ff));% central .. used in inactive path
    P(si).Ppial      = fullfile(pp0_surffolder,sprintf('%s.pial.%s.gii',opt.surf{si},ff));             % pial (GM/CSF)
    P(si).Pwhite     = fullfile(pp0_surffolder,sprintf('%s.white.%s.gii',opt.surf{si},ff));            % white (WM/GM)
    P(si).Pthick     = fullfile(pp0_surffolder,sprintf('%s.thickness.%s',opt.surf{si},ff));            % FS thickness / GM depth
    P(si).Pmsk       = fullfile(pp0_surffolder,sprintf('%s.msk.%s',opt.surf{si},ff));                  % msk
    P(si).Ppbt       = fullfile(pp0_surffolder,sprintf('%s.pbt.%s',opt.surf{si},ff));                  % PBT thickness / GM depth
    P(si).Psphere    = fullfile(pp0_surffolder,sprintf('%s.sphere.%s.gii',opt.surf{si},ff));           % sphere
    P(si).Pspherereg = fullfile(pp0_surffolder,sprintf('%s.sphere.reg.%s.gii',opt.surf{si},ff));       % sphere.reg
    P(si).Pgmt       = fullfile(pp0,mrifolder, sprintf('%s_thickness-%s.nii',ff,opt.surf{si}));        % temp thickness
    P(si).Pppm       = fullfile(pp0,mrifolder, sprintf('%s_ppm-%s.nii',ff,opt.surf{si}));              % temp position map
    P(si).Pfsavg     = fullfile(opt.fsavgDir,  sprintf('%s.central.freesurfer.gii',opt.surf{si}));     % fsaverage central
    P(si).Pfsavgsph  = fullfile(opt.fsavgDir,  sprintf('%s.sphere.freesurfer.gii',opt.surf{si}));      % fsaverage sphere    
  end
end
%==========================================================================
function useprior = setupprior(opt,surffolder,P,si)
%setupprior. prepare longitidunal files

  % use surface of given (average) data as prior for longitudinal mode
  if isfield(opt,'useprior') && ~isempty(opt.useprior) 
    % RD20200729: delete later ... && exist(char(opt.useprior),'file') 
    % if it not exist than filecopy has to print the error
    [pp1,ff1] = spm_fileparts(opt.useprior);
    % correct '../' parts in directory for BIDS structure
    [stat, val] = fileattrib(fullfile(pp1,surffolder));
    if stat, pp1_surffolder = val.Name; else, pp1_surffolder = fullfile(pp1,surffolder);  end
    
    % try to copy surface files from prior to individual surface data 
    useprior = 1;
    useprior = useprior & copyfile(fullfile(pp1_surffolder,sprintf('%s.central.%s.gii',opt.surf{si},ff1)),P(si).Pcentral,'f');
    useprior = useprior & copyfile(fullfile(pp1_surffolder,sprintf('%s.sphere.%s.gii',opt.surf{si},ff1)),P(si).Psphere,'f');
    useprior = useprior & copyfile(fullfile(pp1_surffolder,sprintf('%s.sphere.reg.%s.gii',opt.surf{si},ff1)),P(si).Pspherereg,'f');
    
    if ~useprior
      warn_str = sprintf('Surface files for %s not found. Move on with individual surface extraction.\n',pp1_surffolder);
      fprintf('\nWARNING: %s',warn_str);
      cat_io_addwarning('cat_surf_createCS4:noPiorSurface', warn_str);
    else
      fprintf('\nUse existing surface as prior and thus skip many processing steps:\n%s\n',pp1_surffolder);
    end      
  else
    useprior = 0;
  end
end
function CS = smoothArt(Yth1i,P,CS,Smat,Vppm,facevertexcdatanocut,si,opt,method,refine)
%% Topology artifact correction (after first deformation):
%  Although the topology is fine after correction the geometry often
%  suffers by small snake-like objects that were deformed close to the 
%  central layer, resulting in underestimation of the FS thickness.
  pbtthick = cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm) .* facevertexcdatanocut; 
  cat_io_FreeSurfer('write_surf_data',P(si).Ppbt,pbtthick);
  cmd = sprintf('CAT_SurfDistance -mean -thickness "%s" "%s" "%s"', P(si).Ppbt, P(si).Pcentral, P(si).Pthick);
  cat_system(cmd,opt.verb-3);
  fsthick  = cat_io_FreeSurfer('read_surf_data',P(si).Pthick);  
  
  % define outlier maps
  tart   = (pbtthick - fsthick)>.5 | (pbtthick - fsthick)>.5;                        % artifacts from topology correction   
  [~,dI] = unique(CS.vertices,'rows'); SM = true(size(CS.vertices,1),1); SM(dI) = false; % vertices with same coordinates
  iarea  = 1./cat_surf_fun('area',CS);                                                   % face area to identify tiny faces (=artifacs)
  
  if method % old function 
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
    %saveSurf(CS,P(si).Pcentral);
    cmd  = sprintf('CAT_Central2Pial -equivolume -weight 0.55 "%s" "%s" "%s" 0',P(si).Pcentral,P(si).Ppbt,P(si).Pcentral);
    cat_system(cmd,opt.verb-3);
    if 0
      CS2  = loadSurf(P(si).Pcentral);
      tart = repmat( tart | SM | iarea>1000 , 1,3);
      CS.vertices = CS.vertices.*(tart) + (~tart).*CS2.vertices; 
      saveSurf(CS,P(si).Pcentral);
    end
  end

  cmd = sprintf('CAT_SurfDeform -iter 100 "%s" "%s" "%s"', Vppm.fname, P(si).Pcentral, P(si).Pcentral);
  cat_system(cmd,opt.verb-3);
  
  CS = loadSurf(P(si).Pcentral);
end
%==========================================================================
function [Vppm,rel] = exportPPmap( Yp0, Yp0fs, Yppi, Vmfs, Ypbs, vx_vol, si, opt, surffolder, ff)
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
  optimize = (opt.SRP>0) .* 3; 
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
      %fprintf('(Yppp=%0.2f)',mean(Yppp(Yppi(:)>0 & Yppi(:)<1))); 
    otherwise
      Ypps = Yppi; 
  end
 
  %% Final downsampling:
  % - We combine here the average position with the GM/WM interface.
  %Ypps = cat_vol_median3(Ypps,Ypps>0 & Ypps<1);
  fs = 0; %.35; 
  if opt.interpV == opt.reconres
    Ytx            = Ypps; % smooth3( Ypps ); % max(Ypps,max(0,min(1,Yp0fs-2)) .^ .25); % smooth3?
    Vppm           = Vmfs;
  else
    [Vpp,Vppr]     = prepareDownsampling(Vmfs,Ypps,surffolder,ff,opt,si);
    Vpp.dat(:,:,:) = cat_vol_smooth3X( Ypps , fs ); 
    [Vppm,Yt]      = cat_vol_imcalc(Vpp,Vppr,'i1',struct('interp',5,'verb',0,'mask',-1));
    Vpp.dat(:,:,:) = cat_vol_smooth3X( max(0,min(1,Yp0fs-2)) .^ rel, fs ); 
    [~,Yw]         = cat_vol_imcalc(Vpp,Vppr,'i1',struct('interp',5,'verb',0,'mask',-1));   
    Vpp.dat(:,:,:) = cat_vol_smooth3X( max(0,min(1,2-Yp0fs)) .^ rel, fs );
    [~,Yc]         = cat_vol_imcalc(Vpp,Vppr,'i1',struct('interp',5,'verb',0,'mask',-1));
    Vppm.pinfo     = Vmfs.pinfo;
    Ytx            = Yt; %min(max(0,1-Yc), max(Yt,Yw)); 
  end
 
  if isfield(Vppm,'dat'), Vppm = rmfield(Vppm,{'dat'}); end
  Vppm = spm_write_vol(Vppm, Ytx ); 

  if 0
    Vppm.fname = fullfile(surffolder,sprintf('%s.pp.%s.r%0.2f.surfrecon.nii',opt.surf{si},ff,opt.reconres));  
    spm_write_vol(Vppm, Ytx ); 
  end
end

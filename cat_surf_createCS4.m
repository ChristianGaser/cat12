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
  skip_registration       = isfield(opt,'surf') && isscalar(opt.surf); % skip spherical registration for quick tests
  create_white_pial       = 1; % uses only the quick WM and Pial surface estimation 

  myelinCorrection        = .5; % .25 - sight correction, 1 - maximum correction
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
  def.reconres            = .5;                                
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
    case 1 % myelincorrection 
      use_cat_vol_pbtsimple = 1; 
      myelinCorrection      = 0; 
    case 2 % default
      use_cat_vol_pbtsimple = 1; 
      myelinCorrection      = 0.5; 
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
  [P,pp0,mrifolder,surffolder,surfdir,ff] = cat_surf_createCS_fun('setFileNames',V0,job,opt); 
  

  % main loop for each surface structure 
  for si = 1:numel(opt.surf)
   
    % print something
    if si==1, fprintf('\n'); end; fprintf('%s:\n',opt.surf{si});
    clear Ynocerebrum

    % prepare longitudinal case if required 
    useprior = cat_surf_createCS_fun('setupprior',opt,surfdir,P,si);


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
      cmd = sprintf('CAT_VolThicknessPbt -median-filter 2 -downsample 0 "%s" "%s" "%s"', Vmfs.fname, P(si).Pgmt, P(si).Pppm);
      cat_system(cmd,opt.verb-3);
      Vgmt = spm_vol(P(si).Pgmt); Yth1i = spm_read_vols(Vgmt); 
      Vppi = spm_vol(P(si).Pppm); Yppi  = spm_read_vols(Vppi); 
      Ymfsc = Yp0fs; 

      % correction of general offset in mm 
      Yth1i = max(0,Yth1i - 0.3); 
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
        cat_system(cmd ,opt.verb-3);
        Vppm = spm_write_vol(Vppm, Yppi);
        if 0 % minor improvement of values (IE-0.003, PE=-0.007) but not visually but 35s (50%) increased processing time >> not useful
          cmd = sprintf('CAT_SurfDeform -iter 100 "%s" "%s" "%s"', Vppm.fname, P(si).Pcentral, P(si).Pcentral);
          cat_system(cmd,opt.verb-3);
        end    
        
  
        % load and check surface
        % use 80000 as lower limit what is similar to vdist=4
        CS  = loadSurf(P(si).Pcentral); 
        A   = sum( cat_surf_fun( 'area', CS )); % in mm2 
        nCS = min( size(CS.faces,1) , max( 80000, min( 360000, round(A * 2 * 4 / opt.vdist^2) ))); 
        if opt.SRP >= 2, nCS = min(300000,nCS); end
        EC0 = size(CS.vertices,1) + size(CS.faces,1) - size(spm_mesh_edges(CS),1);  
        % test for box-error 
        % in some cases (SPM preprocessing) the topology defects are so 
        % severe that the topology correction just fill everything
        IS  = cat_surf_fun('isocolors',Yppi,CS.vertices,Smat.matlabIBB_mm); 
        if EC0 ~= 2 || mean(IS(:)) < .4 
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
        cmd  = sprintf('CAT_SurfReduce -aggr "7" -ratio "%0.2f" "%s" "%s"', nCS/nCS0 , P(si).Pcentral, P(si).Pcentral); %5-0.05*opt.reconres);
        cat_system(cmd ,opt.verb-3);  
        CS = loadSurf(P(si).Pcentral); 
        nCS1 = size(CS.faces,1); 
        
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
    % - minor improvement of values (IE+.001,PE+0.002) but visually strong reduction of local outliers >> useful
    stime = cat_io_cmd('  Filter topology artifacts','g5','',opt.verb,stime);
    CS = smoothArt(Yth1i,P,CS,Smat,Vppm,1,si,opt,1,1); % last options - method & refine
    saveSurf(CS,P(si).Pcentral); 


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



    %% surface deformation (after smoothArt) 
    % CAT_DeformSurf in combination with smoothing correction CAT_Central2Pial 
    % allows PP values that are closer to the central position (std 0.049 vs. 0.73)
    % what improves also intensity and position values of the inner/outer surfaces!
    stime = cat_io_cmd('  Optimize surface','g5','',opt.verb,stime); 
    switch 2
      case 1
        cmd = sprintf('CAT_SurfDeform -iter 25 "%s" "%s" "%s"', Vppm.fname, P(si).Pcentral, P(si).Pcentral);
        cat_system(cmd,opt.verb-3);
      case 2
        % further iteration allow (only) slight improvements 
        cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none  0  1  -1  .1 ' ...           
                  'avg  -0.1  0.1 .2  .1  5  0 "0.5"  "0.5"  n 0  0  0 %d  %g  0.0 0'], ...    
                   Vppm.fname,P(si).Pcentral,P(si).Pcentral,25,0.001);  
        cat_system(cmd,opt.verb-3);
        cmd = sprintf('CAT_Central2Pial -equivolume -weight 0.5 "%s" "%s" "%s" 0',P(si).Pcentral,P(si).Ppbt,P(si).Pcentral);
        cat_system(cmd,opt.verb-3);
    end
    
    % Remove self-intersections
    cmd = sprintf('CAT_SurfRemoveIntersections "%s" "%s"',P(si).Pcentral,P(si).Pcentral);
    cat_system(cmd,opt.verb-3);
    CS = loadSurf(P(si).Pcentral);


    %% Collision correction by Delaunay triangularization
    %  --------------------------------------------------------------------
    %  New self-intersection correction that uses different detections of
    %  self-intersections (SIDs; RY vs. PBT) with/without further optimization. 
    %  It does not fully avoid self-intersections because some are already 
    %  in the CS and some other required strong changes that result in worse
    %  thickness results.
    %  RD202508: Optimization does not fully work on the simple PP map.
    %            The optimization torwards 0.05 and 0.95 on the interpolated maps still result in light thickness overestimation. 
    if opt.SRP > 0 
      facevertexcdatac = facevertexcdata; 
      % quick correction with local optimization 
      [CS,facevertexcdatac,SIs] = cat_surf_fun('collisionCorrectionPBT',CS,facevertexcdatac .* facevertexcdatanocut,Ymfsc,Yppi, ... 
        struct('optimize',0,'verb',isscalar(opt.surf)*2,'mat',Smat.matlabIBB_mm,'vx_vol',vx_vol,'CS4',0)); % CS4 is not working
      if isscalar(opt.surf), fprintf('\b\b'); end
      % accurate but slow function (especially for multiple jobs) 
      [CS,facevertexcdatac,SIs] = cat_surf_fun('collisionCorrectionRY' ,CS,facevertexcdatac .* facevertexcdatanocut,Yppi,...
        struct('Pcs',P(si).Pcentral,'verb',isscalar(opt.surf)*2,'mat',Smat.matlabIBB_mm,'accuracy',1/2^2));
      if isscalar(opt.surf), fprintf('\b\b'); end
      facevertexcdata = max(facevertexcdata,facevertexcdatac) .* facevertexcdatanocut; 
      cat_io_FreeSurfer('write_surf_data',P(si).Ppbt,facevertexcdata);
      saveSurf(CS,P(si).Pcentral); 
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
      % only for test visualization  
      cat_surf_createCS_fun('quickeval',V0,Vpp,Ymfs,Yppi,CS,P,Smat,res,opt,EC0,si,time_sr,4); 
      return
    end

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
      cat_io_FreeSurfer('write_surf_data', P(si).Pthick, min(opt.thick_limit,facevertexcdata) );  
    else
      % otherwise simply use the original values of the PBT map  
      % WARNING: 
      %   The values of the ?h.pbt.* files are used to estimate further 
      %   surfaces and therefore corrected for self-intersections!
      %   We use here the orignal PBT values as they are less depending on 
      %   local highly individual features.  
      facevertexcdata = cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm) .* facevertexcdatanocut; 
      cat_io_FreeSurfer('write_surf_data', P(si).Ppbt,   facevertexcdata);
      cat_io_FreeSurfer('write_surf_data', P(si).Pthick, min(opt.thick_limit,facevertexcdata) );  
    end
    fprintf('\n');
    

    % correct thickness based on folding pattern 
    if opt.foldingcorrection
      cmd = sprintf('CAT_SurfCorrectThicknessFolding -max "%f" "%s" "%s" "%s"', opt.thick_limit, P(si).Pcentral, P(si).Pthick, P(si).Pthick);
      cat_system(cmd,opt.verb-3);
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
    if ~debug && exist(Vgmt.fname ,'file'), delete(Vgmt.fname); end
    if ~debug && exist(Vmfs.fname ,'file'), delete(Vmfs.fname); end
  
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
    Vppm           = Vmfs;
    Ytx            = cat_vol_smooth3X( Ypps , fs );
  else
    [Vpp,Vppr]     = prepareDownsampling(Vmfs,Ypps,surffolder,ff,opt,si);
    Vpp.dat(:,:,:) = cat_vol_smooth3X( Ypps , fs ); 
    [Vppm,Yt]      = cat_vol_imcalc(Vpp,Vppr,'i1',struct('interp',5,'verb',0,'mask',-1));
    %Vpp.dat(:,:,:) = cat_vol_smooth3X( max(0,min(1,Yp0fs-2)) .^ rel, fs ); 
    %[~,Yw]         = cat_vol_imcalc(Vpp,Vppr,'i1',struct('interp',5,'verb',0,'mask',-1));   
    %Vpp.dat(:,:,:) = cat_vol_smooth3X( max(0,min(1,2-Yp0fs)) .^ rel, fs );
    %[~,Yc]         = cat_vol_imcalc(Vpp,Vppr,'i1',struct('interp',5,'verb',0,'mask',-1));
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

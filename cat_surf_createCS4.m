function [Yth,S,P,res] = cat_surf_createCS4(V,V0,Ym,Yp0,Ya,YMF,Yb0,opt,job)
% ______________________________________________________________________
% Surface creation and thickness estimation.
% Here we focus on the segmentation label map as robust standard. 
% All non-cortical regions and blood vessels are removed. Especially, 
% uncorrected high-intensity blood vessels and strong noise can result 
% in severe geometric problems due to inoptimal topology correction. 
% Thin areas with partial volume effects between GM and WM are precorrected
% to avoid underestimations. 
%
% Compared to previous pipelines (cat_surf_createCS[2]), the preparation of 
% maps (e.g. filling of structures and blood-vessle/myelin correction) was 
% reworked. Also the PBT-reconstruction was improved by better multiple 
% distance estimation to better consider partial volume effects and further
% smoothness constrains. Finally, also the surface reconstruction itself
% was strongly updated, by an approach with internal topology correction
% rather than the post-correction by the spherical harmonics. However, this
% cause sulcal blurring in some (low-quality) cases, requiring addition 
% tests. The reconstructed high-resolution surface is then optimized for 
% mesh size by downsampling and local refinement and repositioned to the 
% central position to reconstruct also a general white and pial surface. 
% Mesh resolution is set to about 1 mm to support anatomical details as 
% well as fast, memory/disk-space optimal output without severe self-
% intersections of the additional white/pial surface. 
%
%  [Yth,S,P,res] = cat_surf_createCS4(V,V0,Ym,Ya,YMF,Yb0,opt,job)
% 
%  Output parameters:
%   Yth    .. thickness map
%   S      .. structure with surfaces, like the left hemisphere, that 
%             contains vertices, faces, GM thickness (th1)
%   P      .. structure with name of surface files
%   res    .. intermediate and final surface creation information
%
%  Input parameters:
%   V      .. spm_vol-structure of internally interpolated image
%   V0     .. spm_vol-structure of original image
%   Ym     .. the (local) intensity, noise, and bias corrected T1 image
%             with [0,1/3,2/3,1] for [BG,CSF,GM,WM]
%   Yp0    .. brain tissue label map with background, CSF, GM and WM with 
%             values of [0,1,2,3]
%   Ya     .. the atlas map with the ROIs for left and right hemispheres
%             (this is generated with cat_vol_partvol)
%   YMF    .. a logical map with the area that has to be filled
%             (this is generated with cat_vol_partvol)
%   Yb0    .. modified brain mask from gcut
%   opt    .. structure with parameter settings (see bellow)
%   job    .. main preprocessing job structure 
%
%  opt-structure:
%   (1) main fields
%    .surf               .. surface reconstruction setting with 
%                             {'lh','rh'[,'cb']} 
%    .useprior           .. longitudinal prior file
%    .folding correction .. correction of gyral/sulcal folding bias 
%                           (default = 1)
%    .interpV            .. PBT resolution in mm (default .5)
%    .reconres           .. surface reconstruction resolution to run the 
%                           marching cubes (default: defined by opt.interpV)
%    .vdist              .. Mesh resolution parameter in mm that constrols 
%                           the typical distancebetween mesh vertices by
%                           downsampling and refinement. We aim in 1 mm
%                           mesh resolution (about 2 vertices/mm2) to balance 
%                           anatomical details, mesh properties and 
%                           computational/memory/space requirements demands. 
%    .myelinCorrection   .. correction of thin myelinated (GM-WM PVE) areas 
%                           (0-no, 0.25-light, 0.5-strong, 1.0-maximum)
%                           (default = .3)
%    .create_white_pial  .. 
%    .skip_registration  .. debug
% 
%   (2) further/internal parameters and settings controled by SPM defaults
%    .verb          .. verbose level defined by CAT_defaults expert level
%                        0-none, 1-minimal, 2-default, 3-details, 4-debug
%    .LAB           .. atlas defintion from CAT_defaults
%    .thick_limit   .. default is 6 mm upper limit for thickness (similar limit as used in Freesurfer)
%    .thick_measure .. 
%    .fsavgDir      .. CAT surface template directory 
% 
% See also:
%   cat_surf_ceateCS, cat_surf_ceateCS2, cat_surf_createCS_fun, 
%   cat_vol_pbt2, cat_vol_pbtsimple[CS4], cat_surf_fun
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$ 


  % get both sides in the atlas map
  NS = @(Ys,s) Ys==s | Ys==s+1; 
  if ~exist('opt','var'), opt = struct(); end 
  cstime = clock; %#ok<*CLOCK>
    
  % surface output and evaluation parameter 
  res   = struct('lh',struct(),'rh',struct()); 
  Yth   = zeros(size(Yp0),'single');  % thickness map
  S     = struct();
  

  % set debugging variable
  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end
  vx_vol                  = sqrt(sum(V.mat(1:3,1:3).^2));                  % further interpolation based on internal resolution 

  % set defaults
  def.surf                = {'lh','rh'}; % surface reconstruction setting with {'lh','rh','cb'} 
  def.useprior            = ''; % 
  def.vdist               = 2;  % surface resolution parameter
  def.foldingcorrection   = 1;  % tickness correction that is influence by folding
  def.myelinCorrection    = .3; % .25 - sight correction, 1 - maximum correction
                                % good value between .25 and .5
  def.create_white_pial   = 1;  % uses only the quick WM and Pial surface estimation (0-no,1-yes,2-improve,3-improve+update)
                                % - 
                                % - (isfield(opt,'surf') && isscalar(opt.surf));  
                                % 
  % default set by cat_defaults or quite fix/rarely-used parameters
  def.verb                = cat_get_defaults('extopts.expertgui');         % 0-none, 1-minimal, 2-default, 3-details, 4-debug
  def.LAB                 = cat_get_defaults('extopts.LAB');               % brain regions 
  def.thick_limit         = cat_get_defaults('extopts.thick_limit');       % 6 mm upper limit for thickness (similar limit as used in Freesurfer)
  def.thick_measure       = cat_get_defaults('extopts.thick_measure');     % 0-PBT; 1-Tfs (Freesurfer method using mean(TnearIS,TnearOS)) 
  def.fsavgDir            = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces'); 
  def.outputpp.native     = 0;  % output of native Ypp map for cortical orientation in EEG/MEG, i.e., don't delete intermediate files
  %def.outputpp.warped     = 0; % not implemented
  %def.outputpp.dartel     = 0; % not implemented
  def.skip_registration   = isfield(opt,'surf') && isscalar(opt.surf); % skip spherical registration for quick tests
  
  % options that rely on other options
  opt                     = cat_io_updateStruct(def,opt); clear def; 
  opt.vol                 = any(~cellfun('isempty',strfind(opt.surf,'v')));   % only volume-based thickness and position map estimation  
  opt.surf                = cat_io_strrep(opt.surf,'v','');                   % after definition of the 'vol' varialbe we simplify 'surf'
  opt.interpV             = max(0.1,min([opt.interpV,2]));                    % general limitation of the PBT resolution
  opt.reconres            = opt.interpV;                                      % surface reconstruction resolution for marching cubes
  opt.denoise             = 0;                                                % helps partially but cause problems with local thickness 
                                                                              % - underestimation in simple phantom cases but
                                                                              %   gyral overestimation in real data by breaking WM peaks
                                                                              % - sharpening helps but not sufficient enough
  
  % apply the modified mask from gcut
  % for non-gcut approaches or inverse weighting Yb0 only contains ones
  Yp0 = Yp0 .* (Yb0>0.5);
  Ym  = Ym  .* (Yb0>0.5);
  
  % enlarge atlas definition 
  [~,I] = cat_vbdist(single(Ya>0)); Ya=Ya(I); clear I;  


  %% denoise to recude artifacts and WM hyperintensities 
  %  quite save with clear improvement in surface creation (less defects and smaller corrections) 
  %  thickness overestimation in real low-quality data in case of overfiltering!
  if opt.SRP > 1  &&  opt.denoise
    fast = 1; sharp = 1; 
    Yp0  = denoiseSanlm(Yp0, vx_vol, fast, sharp);
    Ym   = denoiseSanlm(Ym,  vx_vol, fast, sharp);
  end
  
  % estimate reminding noise
  Ywmr = cat_vol_resize(Ym .* (Yp0>2.75/3),'reduceV',vx_vol,2,16,'meanm'); 
  Ywmrstd = cat_vol_localstat(Ywmr,Ywmr>0,2,4); clear Ywmr;
  opt.wmnoise = min(1/6,cat_stat_nanmean(Ywmrstd(Ywmrstd>0))); 
  
  %% The ORNLM filter modifies the overall intensities and we have to rescale the image!
  % This rescaling needs also a slight correction for the noise to avoid underestimation. 
  % However, the filtering migh also removes too many details that cause than overestimations
  % and it is therefore only used in experimental pipeline
  printcount = 0; fprintf('\n'); 
  if opt.SRP > 1  &&  opt.denoise
    fprintf('Denoising (CNR=%0.2f%%) ',opt.wmnoise * 100); 
    printcount = printcount + 1; 
    [Yp0,Ym] = cat_surf_createCS_denoise(Yp0, Ym, opt.wmnoise);
  end


  %% BVC - might be problematic and only relevant to avoid blurring in T2/PD/FLAIR/BVs !
  if job.inv_weighting
    Ygc = cat_vol_localstat(Yp0,Yb0,1,2);
    Ygc = cat_vol_localstat(Ygc,Yb0,1,3);
    Yd  = max(0,(Yp0 - min(Yp0,Ygc)) .* cat_vol_morph(Yb0,'e')); clear Ygc;
    Yd  = cat_ornlm(Yd,1,1,.05); 
    Yd(smooth3(Yd)<.125) = 0; 
    if printcount, mspace = ', '; else, mspace = ''; end
    fprintf('%sBVC (CNR=%0.2f%%)', mspace, median(Yd(Yd(:)>.1))*100); 
    printcount = printcount + 1; 
    Yp0 = min(Yp0, max(2,Yp0 - Yd .* ( cat_vol_smooth3X(Yp0,4) > 2.125/3 & Yp0>2.5/3 & cat_vol_smooth3X(Yp0,4)<2.5))); clear Yd;
    Yp0 = min(Yp0,cat_ornlm(Yp0,1,1,opt.wmnoise));
  end


  %% WMH and PVS correction 
  %  Eg.  NISALS_UTR_SP30T_als3_T1w  dataset with Swiss cheese WM and wmnoise ~ 0.04.
  %  Blood vessels should be already removed here and we should also not apply this in non T1w data! 
  %  This can cause blurring of sulci and is controlled by the wmnoise value!
  if opt.wmnoise < 0.03  &&  ~job.inv_weighting  &&  opt.SRP > 0 
    if printcount, mspace = ', '; else, mspace = ''; end
    fprintf('%sWMHC (%0.2f%%)', mspace, max(0,min(1,opt.wmnoise*100 - 3))); 
    Yp0 = closeWMHandPVS(Yp0, Ya, NS, job.extopts.LAB.CB, opt.wmnoise, vx_vol);
  end
  

  %% simple filling
  %  In contrast to the old fill used in CS2, this pipeline requires stronger 
  %  filling as no spherical correction will remove tiny wormholes. 
  %  However, the filling is critical in closing small sucli!
  [Yp0f,Ymf] = fillVentricle(Yp0,Ym,Ya,YMF,vx_vol,opt); %clear Ym YMF; 


  % prepare file and directory names
  [P,mridir,surfdir,ff] = cat_surf_createCS_fun('setFileNames',V0,job,opt);
  % ensure surface directory exists (may not be created when job.output.surface=0,
  % which causes CAT_VolMarchingCubes to fail silently on Windows)
  if ~exist(surfdir,'dir'), mkdir(surfdir); end  

  % main loop for each surface structure 
  for si = 1:numel(opt.surf)
   
    % print something
    if si==1, fprintf('\n'); end; fprintf('%s:\n',opt.surf{si});
    
    % prepare longitudinal case if required 
    useprior = cat_surf_createCS_fun('setupprior',opt, job.BIDS(1).surfdir,P,si);

    % reduce for object area
    iscerebellum = cat_io_contains(opt.surf{si},'cb') && exist('suit_amap','file'); 
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
        Yside  = zeros(size(Yp0fs),'single'); 
    end 
    

    % removing background (smoothing to remove artifacts)
    [Yp0fs,Ymfs,BB] = cat_vol_resize({Yp0fs,Ymfs}, 'reduceBrain', vx_vol, 4, smooth3(Yp0fs)>1.5); 
    

    % interpolation 
    % cubic is be better than linear although we have to consider interpolation artifacts
    [Yp0fs,resI]    = cat_vol_resize(max(1,Yp0fs),'interp',V,opt.interpV,'cubic'); 
    Ymfs            = cat_vol_resize(max(1,Ymfs),'interp',V,opt.interpV,'cubic');  
    

    if iscerebellum 
    % Run cerebellum specific SUIT segmentation to update the tissue 
    % segmenation designed for the 3 cerebelar layer and their partial 
    % volume.  
      stime = cat_io_cmd('  Cerebellar preparation','g5');
      [~,Ycbd] = cat_vol_downcut( single(cat_vol_morph( smooth3(Ymfs)>2.1,'l')),(Yp0fs.^.1) - 1,1); 
      cbdth    = cat_stat_nanmedian(Ycbd(Yp0fs>1.1 & Yp0fs<1.25)); 
      Ycbs     = Yp0fs>1.1 & Ycbd < cbdth; 
      Pcb      = suit_amap('run',Yp0fs .* Ycbs,max(0,(Ycbs .* Ymfs) - 1)/2,vx_vol);
      Yp0fs    = min(3,max(1,(Pcb.Yml) * 3));
      fprintf('%5.0fs\n',etime(clock,stime)); 
    end
    

    % surface coordinate transformation matrix
    [Vmfs,Smat] = createOutputFileStructures(V,V0,resI,BB,opt,mridir,ff,si); 
    

    % Myelination corretion 
    % **** this is not fully save working in all cases and can cause severe topology issues (phantom) ****
    if opt.myelinCorrection > 0  &&  ~iscerebellum  &&  opt.SRP>0  
      stime = cat_io_cmd(sprintf('  Myelin correction (%0.1f%%)',opt.myelinCorrection*100),'g5'); 
      Yp0fs = myelincorrection(Yp0fs, vx_vol, opt, P, Vmfs, si, opt.SRP>0); 
      fprintf('%5.0fs\n',etime(clock,stime)); 
    end  
    
    
    stime = cat_io_cmd(sprintf('  Thickness estimation (CS4%d,%0.2f mm%s)', opt.SRP, opt.interpV,native2unicode(179, 'latin1')));
    Vmfs.dt = [16 1]; spm_write_vol(Vmfs, Yp0fs );
    if opt.SRP > 0
      [Yth1i,Yppi] = cat_vol_pbtsimpleCS4(Yp0fs, opt.interpV,...
        struct('verb',1, 'gyrusrecon',~iscerebellum, 'eidist',1, ...
        'NVBC',~iscerebellum, 'denoise',~iscerebellum, 'wmnoise', opt.wmnoise)); 
      Vppm = Vmfs; Vppm.fname = P(si).Pppm; spm_write_vol(Vppm, Yppi);
    else 
      cmd = sprintf('CAT_VolThicknessPbt  -range 0.45  -correct-voxelsize 0  "%s" "%s" "%s"', Vmfs.fname, P(si).Pgmt, P(si).Pppm);
      cat_system(cmd,3);
      Vgmt  = spm_vol(P(si).Pgmt); Yth1i = spm_read_vols(Vgmt); 
      % correction of general offset in mm 
      % Estimated for various interpolation resolution of 
      %  sub-ds000113-1.3.0-01_ses-forrestgump_desc-snr30_rf30_2_res-07mm_thickness20mm-30mm_T1w.nii.gz
      Yth1i = max(0,min(10,Yth1i + 0.25 - 0.35*mean(opt.interpV) )); % update
      Vppm = Vmfs; Vppm.fname = P(si).Pppm; 
      Vppm = spm_vol(P(si).Pppm); Yppi = single( spm_read_vols(Vppm) ); 
    end
  

    % denoise position map 
    % Can correct some outliers that might become topology defects.
    % We use the minimum here to avoid sulcal blurring what is more risiky
    % for the current pipeline. As the last tests showed slighly worse 
    % results denoising we avoid in the default pipelines.
    if opt.SRP > 0 
      Yppis = cat_ornlm(Yppi + 1,1,1,min(1/4,opt.wmnoise / 2)) - 1;
      Yppis = (Yppis - min(Yppis(:))) ./ (max(Yppis(:)) - min(Yppis(:))); 
      Yppi  = min(Yppi,Yppis); clear Yppis
    end


    % prepare thickness output
    Yth1t = cat_vol_resize(Yth1i,'deinterp',resI);                         % back to original resolution
    Yth1t = cat_vol_resize(Yth1t+1,'dereduceBrain',BB) - 1;                % adding background
    Yth   = max(Yth,Yth1t .* Yside); clear Yth1t                           % save on main image


    %%
    if opt.vol
      % only voxelbased thickness map
      S = struct(); P = '';
      if opt.verb<2, fprintf('%5.0fs',etime(clock,stime)); end %#ok<*DETIM>
      continue; 
    end
 
    
    %%  surface creation 
    %  --------------------------------------------------------------------
    %  In case of the central surface the frequency patter between sulci and 
    %  gyri is more harmonized and even lower resolution (eg. 1.5 and 2.0) 
    %  for surface reconstruction are  are possible. 
    %  --------------------------------------------------------------------

    %  Surface creation control variables to detect severe geometrical 
    %  changes of the topology correction (e.g. closing of large sulci).
    res.help = [
      'All fields come with values for each hemisphere si and applied position ' ...
      'map thresholds thi. \n' ...
      '  area_gt(si)       = surface area of the uncorrected initial surface \n' ...
      '  thi(si)           = final threshold \n' ...
      '  pe(si,thi)        = position error of the central surface \n' ...
      '  ie(si,thi)        = segment intensity error of the central surface \n' ...
      '  pie(si,thi)       = position and intensity error of the central surface \n' ...
      '  gycon(si,thi)     = surface reconstruction error as product of area error and surface genus \n'...
      '  ECf(si,thi)       = Euler Characteristic\n'...
      '  genus(si,thi)     = Surface genus number abs(EC-2)\n' ...
      '  area_tc(si,thi)   = area of the topology corrected surface\n' ...
      '  area_uc(si,thi)   = area of the topology uncorrected surface\n'...
      '  surferr(si,thi)   = \n'...
      '  surferrgt(si,thi) = \n'...
      '  ECmodvx(si)       = \n'...
      ];
    
    if isscalar(opt.surf), time_sr = clock; end % temporary for tests 
    areaerr  = .005;  % accepted percentage area changes by topology correction ... 0.5% are good
    genuserr = 6;     % accaptable genus error of the final surface 
    final    = 0;     % flag to use specific settings in the last run 

    if useprior 
    % Longitudinal processing block that has to set up all the variables   
    % that control the surface creation. 
      stime = cat_io_cmd('  Load and refine subject average surface','g5','',opt.verb,stime);
      res.EC(si)            = 0; 
      res.ECmodvx(si)       = res.EC(si);
      % use average to define gt area
      CS                    = loadSurf(P(si).Pcentral); 
      thi                   = 1; 
      res.area_gt(si)       = sum(cat_surf_fun('area',CS)) / 100; 
      res.thi(si)           = 1; 
      res.pe(si,thi)        = 0; 
      res.ie(si,thi)        = 0; 
      res.pie(si,thi)       = 0; 
      res.gycon(si,thi)     = .5; 
      res.ECf(si,thi)       = size(CS.vertices,1) + size(CS.faces,1) - size(spm_mesh_edges(CS),1); 
      res.genus(si,thi)     = abs(res.ECf(si,thi) - 2);
      res.area_tc(si,thi)   = res.area_gt(si);
      res.area_uc(si,thi)   = res.area_gt(si);
      res.surferr(si,thi)   = 0;
      res.surferrgt(si,thi) = 0;
      res.ECmodvx(si)       = 0;
      clear CS; 
    else
      % optimized downsampling of the the Ypp map and
      Vp0 = Vmfs; Vp0.fname = spm_file(P(si).Pp0,'suffix','_tmp'); 
      spm_write_vol(Vp0, Yp0fs);

      % static/dynamic surface threshold
      % Using a global dynamic threshold to avoid topology defects that
      % is typically resulting in multiple reconstruction loops and often
      % cause severe closing of sulci. 
      % In case of low resolution reconstruction the threshold has to be
      % higher to avoid fatal closing of even major sulci. The threshold
      % has to be higher because the WM in gyri is often thicker then the
      % CSF in sulci, so that higher threshold result in a more balance
      % representation of object/non-object.
      if iscerebellum
        % for the cerebellum we uses the segmentation directly 
        gth   = .65;
        gycon = .65; 
        %% export map
        opt0  = opt; opt0.reconres = min(0.5,opt.reconres);   % minimum resolution 
        Yppi  = Pcb.Yta + Pcb.Yppi; 
        Yppi  = cat_ornlm(Yppi,3,4,.33);
        Yppi  = max(0,Yppi - prctile(Yppi(:),0.01)) ./ ( prctile(Yppi(:),99.99) - prctile(Yppi(:),0.01) );
        Yppi  = min(1,max(0,Yppi + .5*(Yppi - smooth3(Yppi))));
        Vppmi = exportPPmap( Yppi , Vmfs, opt0); % write temporary file only for surface reconstruction 
      else
        gth   = .5; % default value for the uncorrected surface
        relwm = sum(max(0,Ymf(:)-2)) / sum(min(1,Ymf(:))); 
        gycon = max(0,1 - max(0,relwm * 3 - 1) * 2);
        gycon = max(0.55,min(.85,.65 - .1*gycon)); 
        
        % export map
        if opt.SRP > 0
          if opt.SRP > 0
            % sharpening operation - this improves intensity and position values
            Yppis = smooth3(Yppi);
            Yppii = ( min(1,max(0,Yppi + (Yppis - ...
              cat_vol_smooth3X(Yppis,2/mean(vx_vol))/3 - ...
              cat_vol_smooth3X(Yppis,4/mean(vx_vol))/3 - ...
              cat_vol_smooth3X(Yppis,8/mean(vx_vol))/3  ))));
          else
            Yppii = Yppi; 
          end
          opt0  = opt; opt0.reconres = min(0.5,opt.reconres);   % minimum resolution 
          Vppmi = exportPPmap( Yppii .^ max(1,1 + opt.wmnoise*100 ) , Vmfs, opt0); % write temporary file only for surface reconstruction 
          Vp0   = exportPPmap( Yp0fs , Vp0, opt0); % write temporary file only for surface reconstruction 
        else
          Vppmi = spm_vol(P(si).Pppm); 
        end
      end
     
      % create "gold-standard" surface 
      cmd = sprintf('CAT_VolMarchingCubes  "%s" "%s" -thresh "%0.2f" -fast', Vppmi.fname, P(si).Pcentral,gth);
      cat_system(cmd ,opt.verb-3);  
      CSG05 = loadSurf(P(si).Pcentral);
      res.area_gt(si) = sum(cat_surf_fun('area',CSG05)) / 100; 
      if ~debug, clear CSG05; end 

      %% Iterative runs to create the central surface.
      stime = cat_io_cmd(sprintf('  Create initial surface (%0.2f mm)',opt.interpV),'g5','',opt.verb,stime); 
      for thi = 1:round(8/mean(opt.interpV))  % with this loop we further increase the threshold of the Ypp map by .1 (eg. more WM like surface)
        %% Main initial surface creation with topology correction.
        if exist(P(si).Pcentral,'file'), delete(P(si).Pcentral); end % have to delete it to get useful error messages in case of reprocessing/testing
        gycon = min(.9, gycon + .04 / opt.interpV * (thi>1)); % + (opt.wmnoise*2)); % higher threshold in case of noisy data (may issues with PVEs)
 
        % create basic surface to know the real surface geometry
        cmd = sprintf('CAT_VolMarchingCubes  "%s" "%s" -thresh "%0.4f" -fast', Vppmi.fname, P(si).Pcentral, gycon);
        cat_system(cmd ,opt.verb-3);  
        CSG{thi} = loadSurf(P(si).Pcentral); %#ok<AGROW>

        if iscerebellum
          % RD202603: just as quick placeholder ...
          cmd = sprintf('CAT_VolMarchingCubes  "%s" "%s"  -thresh "%0.4f" -iter 1 -median-filter "%d"',  ...
            Vppmi.fname, P(si).Pcentral, gycon, 0); %#ok<NASGU>
          txt = evalc('cat_system(cmd ,3)');
        else
          cmd = sprintf('CAT_VolMarchingCubes  "%s" "%s"  -thresh "%0.4f" -iter %d -label "%s" -median-filter "%d"',  ...
            Vppmi.fname, P(si).Pcentral, gycon, 1+final, Vp0.fname, min(1,floor( .5 / opt.interpV ) )); %#ok<NASGU>
          txt = evalc('cat_system(cmd ,3)'); 
        end
        CST{thi} = loadSurf(P(si).Pcentral); %#ok<AGROW>

        % position and intensity errors
        pe  = max(eps,cat_surf_fun('isocolors',Yppi,CST{thi}.vertices,Smat.matlabIBB_mm)) - .5;
        ie  = max(eps,cat_surf_fun('isocolors',Ymfs,CST{thi}.vertices,Smat.matlabIBB_mm)) - 2; 
        ppe = sign(pe) .* abs(pe).^2; 
        iie = sign(ie) .* abs(ie).^2; 
        
        % save surface characteristics 
        res.pe(si,thi)        = cat_stat_nanmean(pe.^2).^.5; clear pe; 
        res.ie(si,thi)        = cat_stat_nanmean(ie.^2).^.5; clear ie; 
        res.pie(si,thi)       = (cat_stat_nanmean(ppe < -0.1) + cat_stat_nanmean(iie < -0.1))/2 * 100; clear ppe iie;
        res.gycon(si,thi)     = gycon; 
        res.genus(si,thi)     = abs(str2double( txt(max(strfind(txt,':'))+1 : max(strfind(txt,'('))-2) ) - 2); 
        res.area_tc(si,thi)   = sum(cat_surf_fun('area',CST{thi})) / 100;
        res.area_uc(si,thi)   = sum(cat_surf_fun('area',CSG{thi})) / 100;
        res.surferr(si,thi)   = abs(res.area_tc(si,thi) - res.area_uc(si,thi)) / res.area_uc(si,thi) * (res.genus(si,thi)+1)/genuserr;
        res.surferrgt(si,thi) = abs(res.area_tc(si,thi) - res.area_gt(si))     / res.area_gt(si)     * (res.genus(si,thi)+1)/genuserr - (gycon - gth)/16;
        res.ECf(si,thi)       = str2double( txt(max(strfind(txt,':'))+1 : max(strfind(txt,'('))-2) ); 
        res.ECmodvx(si)       = str2double( txt(max(strfind(txt,'('))+1 : max(strfind(txt,'voxel'))-1) );
        res.SE(si,thi)        = mean( [ res.surferr(si,thi)*100 res.surferrgt(si,thi)*100 res.pie(si,thi) ].^2 ).^.5;
       
        %% evaluate surface characteristics
        if opt.verb>1 && cat_get_defaults('extopts.expertgui')>0
          cat_io_cprintf([min(1,.5 + .5*res.SE(si,thi)) max(0,.5 - .5*res.SE(si,thi)) .5],sprintf( ... 
            '\n    Run Ypp>%0.2f (EC=%3d, SE=%6.3f%%)                             ', ...
             res.gycon(si,thi), res.ECf(si,thi), res.SE(si,thi) ));
        end
        % stop iteration ... value over .8 are critical with dynamic surface thresholding
        if final || gycon > .8 || res.SE(si,thi) < .5, break; end
        % in case we havent stopped for a good reason (low blurring and acceptable euler number),
        % we try to find the best threshold so far and do a final run with further internal interation  
        if gycon > .75 
          err = res.surferr(si,:) + res.surferrgt(si,:); 
          thi = find( err == min(err), 1, 'first'); %#ok<FXSET>
          final = 1; 
        end
      end
      

      % use best surface 
      if thi > 1 
        err = res.SE(si,:); 
        thi = find( err == min(err), 1, 'first');
        saveSurf( CST{thi} , P(si).Pcentral); 
      end
      res.EC(si)  = res.ECf(si,thi); 
      res.thi(si) = thi; 
      if res.genus(si,thi) > genuserr
        cat_io_addwarning([mfilename ':EC'],sprintf('Incorrect Euler Number i.e. reminding topological defects (EC=%d). ', ...
          res.EC(si)),1,[1 2*(res.surferr(si,thi) < areaerr)]);
      end
      if res.surferr(si,thi) > areaerr
        cat_io_addwarning([mfilename ':MarchingCubes:' opt.surf{si}],sprintf('Increased topological adjustments of %0.2f%%%% in %s surface. ',...
          res.surferr(si,thi),opt.surf{si}), 1,[1 2]);
      end
      
      % load and check surface and define lower surface resolution based on the surface area of the subject.
      % Use 80000 as lower limit what is similar to vdist=4.
      % More vertices will increase the number of self-intersections in folded parts.
      % Hence, less vertices are not just faster but also cause less conflicts.
      CS   = loadSurf(P(si).Pcentral); 
      CS.facevertexcdata = max(eps,cat_surf_fun('isocolors',Yppi,CS.vertices,Smat.matlabIBB_mm)); 
      A    = sum( cat_surf_fun( 'area', CS )); % in mm2 
      nCS0 = size(CS.faces,1); 
      nCS  = min( size(CS.faces,1) , max( 80000, min( 360000, round(A * 4 / opt.vdist^2) ))); % ###### A * 8 
      nCS  = nCS * (1 + 3*iscerebellum);
   
      
      % surface reduction 
      % - low surface resolution are highly important for fast processing and smaller disk usage
      % - less critical for thickness but surface quality (intensity/position values) get worse
      % - it is important avoid self-intersetions
      if nCS < nCS0 
        stimer = clock; 
        cmd = sprintf('CAT_SurfReduce -aggr "7" -ratio "%0.2f" "%s" "%s"', nCS/nCS0 , P(si).Pcentral, P(si).Pcentral); 
        cat_system(cmd ,opt.verb-3);  
        cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f', P(si).Pcentral, P(si).Pcentral, opt.vdist * 2);
        cat_system(cmd ,opt.verb-3);

        % load
        CS  = loadSurf(P(si).Pcentral); 
        CS.facevertexcdata = max(eps,cat_surf_fun('isocolors',Yppi,CS.vertices,Smat.matlabIBB_mm)); % #### for tests
        nCS1 = size(CS.faces,1); 
        cat_io_cmd(sprintf('  Reduce surface (nF: %0.0fk > ~%0.0fk)', nCS0/1000, nCS1/1000),'g5','',...
          opt.verb,stime - [0 0 0 0 0 etime(stimer,stime)]);
        stime = stimer;
      
        % evaluation
        if nCS1 > nCS * 1.5 
          cat_io_addwarning('cat_surf_createCS4:ReductionFailed', ...
            sprintf('Surface reduction incomplete due to topology constrains (%d).',(nCS1/1000-nCS/1000) ./ (nCS0/1000-nCS/1000)),0,[0 1]); 
        end 
      end
      % == end surface creation == 


      %% Surface position corretion 
      stime = cat_io_cmd('  Optimize surface','g5','',opt.verb,stime); 
      if opt.SRP > 0
        exportPPmap( Yppi , Vmfs, opt0); % write temporary file only for surface reconstruction 
      end

      % Surface deformation: 
      %  This is the function that I would like to use but is is not working good enough,
      %  ie. the position values are not close enough to the center and deviating localy 
      %  extremly what cause problems in developing the inner and outer surface by using 
      %  half thickness as start value. 
      %  Also the old surface deformation (CAT_DeformSurf) is not ideal as it is too
      %  smooth and loosing anatomical details.
      cmd = sprintf('CAT_SurfDeform -iter 75 -isovalue 0.5 -w1 0.1 -w2 0.1 -w3 1.0 -sigma 0.1 "%s" "%s" "%s"', ...
        P(si).Pppm, P(si).Pcentral, P(si).Pcentral);
      cat_system(cmd,opt.verb-3);
      CS = loadSurf(P(si).Pcentral); 
      CS.facevertexcdata = max(eps,cat_surf_fun('isocolors',Yppi,CS.vertices,Smat.matlabIBB_mm)); % #### for tests
      if opt.SRP > 0
        %% ==== Conceptual solution of the issue that the surface creation is biased ====
        %  To move the vertices to their correct position, we utilize the MATLAB isosurface function. 
        %  We cannot use the CSG05 as this is meandering to strong representing incorrect isovalues.
        %  The alignment is done by a nearest point search utilizing the dealunay triangulation. 
        %  The problem is that matlab isosurface returns grid-aligned vertices that result in bad
        %  face properties. 
        
        %CSM = isosurface(interp3(Yppi,1),0.5); CSM.vertices = CSM.vertices/2; % interpoltion can further improve but with high comp. costs
        CSM = isosurface(Yppi,0.5); 
        % add arificial boundary points to avoid error in delaunayn, eg. in MRART sub-561646_acq-headmotion1_T1w: 
        %  QH6417 qhull precision error (qh_merge_twisted): twisted facet f4322286 does not contain pinched vertices.  
        sz  = size(Yppi); 
        CSM.vertices(end+1:end+8,:) = [sz.*[0 0 0]; sz.*[0 0 1]; sz.*[0 1 0]; sz.*[1 0 0]; 
                                       sz.*[1 1 1]; sz.*[1 1 0]; sz.*[1 0 1]; sz.*[0 1 1]]; 
        CSM = cat_surf_fun('smat',CSM,Smat.matlabIBB_mm);
        
        D1  = delaunayn(double(CSM.vertices)); % CSG05
        NN  = dsearchn(double(CSM.vertices),D1,double(CS.vertices)); 
        % We now mix the intial surface and the MATALB isosurface aligned positions.
        % We move vertices based on their deviation from the ideal popsition, avoiding so loosing of nice face properties of the initial mesh.  
        dev = repmat( min(1,abs(cat_surf_fun('isocolors',Yppi,CS.vertices,Smat.matlabIBB_mm) - .5)).^.5 , 1,3);
        CS.vertices = CS.vertices.*(1-dev) + (dev).*CSM.vertices(NN,:);
        CS.facevertexcdata = max(eps,cat_surf_fun('isocolors',Yppi,CS.vertices,Smat.matlabIBB_mm)); % #### for tests
     

        % Anyway, we have to filter the surface to avoid bad faces but also noisy surface parts and "wormholes". 
        % - this can maybe partially avoided by the more aggressive surface reduction
        % - although the filter reduces some artefacts in problematic cases the surface values are worse and the issue are still not fully solved
        % - tried to use it after correction of self-intersections but this made it even worse
        % - minor improvement of values (IE+.001,PE+0.002) but visually strong reduction of local outliers >> useful
        % ### This is quite slow and could be improved. 
        CS  = smoothArt(Yth1i,P,CS,Smat,Vppm,1,si,opt,1,0); % last options - method & refine

        % After removing problematic parts it is again possible to fit to the closest position (causing some but less issues as before)
        dev = repmat( min(1,abs(cat_surf_fun('isocolors',Yppi,CS.vertices,Smat.matlabIBB_mm) - .5)).^.5 , 1,3);
        NN  = dsearchn(double(CSM.vertices),D1,double(CS.vertices)); 
        CS.vertices = CS.vertices.*(1-dev) + (dev).*CSM.vertices(NN,:);
        CS.facevertexcdata = max(eps,cat_surf_fun('isocolors',Yppi,CS.vertices,Smat.matlabIBB_mm)); % #### for tests

        if ~debug, clear NN D1 dev CSM; end
      end

      
      % save our final surface
      saveSurf(CS,P(si).Pcentral); 
      facevertexcdata = max(eps,cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm)); 
      cat_io_FreeSurfer('write_surf_data',P(si).Ppbt,facevertexcdata);

      if opt.SRP > 0
      % #### not sure if this is usefull ... could be done even before ... helpful for surface deformation
      % write extended position map to optimize for the inner/outer surfaces
        Ymfss = cat_vol_smooth3X(Ymfs,2);
        Yppi(Ymfs>2.5) = max(Yppi(Ymfs>2.5),(Ymfss(Ymfs>2.5)-1.5)*1); 
        Yppi(Ymfs<1.5) = min(Yppi(Ymfs<1.5),(min(1.5,Ymfss(Ymfs<1.5))-1.5)*1); % negative
        Vppm = Vmfs; Vppm.fname = P(si).Pppm; spm_write_vol(Vppm, Yppi);
      end
      
    end


    % get PBT thickness 
    CS = loadSurf(P(si).Pcentral);
    facevertexcdata = max(eps,cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm)); 
    cat_io_FreeSurfer('write_surf_data',P(si).Ppbt,facevertexcdata);


    % create white and central surfaces
    if opt.create_white_pial  
      stime = cat_io_cmd('  Create pial and white surface','g5','',opt.verb,stime); 
      if useprior || opt.create_white_pial == 1 
        if 0
          %% quick interal preview - this cuts partially into the tissues
          cmd = sprintf('CAT_Central2Pial "%s" "%s" "%s" 0.5',P(si).Pcentral,P(si).Ppbt,P(si).Ppial);
          cat_system(cmd,opt.verb-3);
          cmd = sprintf('CAT_Central2Pial "%s" "%s" "%s" -0.5',P(si).Pcentral,P(si).Ppbt,P(si).Pwhite);
          cat_system(cmd,opt.verb-3); 

        else
          %% default reconstruction 
          %  that use further adaptations/deformations to reduce the self-intersections

          % optimize surface for basic self-intersections
          if 1
            CS2 = cat_surf_fun('collisionCorrectionPBT',CS,facevertexcdata,Ymfs,Yppi, ...
              struct('optimize',0,'verb',0,'mat',Smat.matlabIBB_mm,'vx_vol',vx_vol,'CS4',1));
            saveSurf(CS2,P(si).Pcentral);
          end

          % (1) create simple boundary surfaces  
          cmd = sprintf('CAT_Central2Pial "%s" "%s" "%s" 0.5',P(si).Pcentral,P(si).Ppbt,P(si).Ppial);
          cat_system(cmd,opt.verb-3);
          cmd = sprintf('CAT_Central2Pial "%s" "%s" "%s" -0.5',P(si).Pcentral,P(si).Ppbt,P(si).Pwhite);
          cat_system(cmd,opt.verb-3); 
          CSw0 = loadSurf(P(si).Pwhite);
          CSp0 = loadSurf(P(si).Ppial);
          saveSurf(CS,P(si).Pcentral);
        
          cat_io_FreeSurfer('write_surf_data', P(si).Pthick, facevertexcdata);
          % (2) create refined boundary surfaces
          if useprior
            % no refinement
            cmd = sprintf('CAT_SurfDeform -iter 100 -isovalue 1 -w1 0.0001 -w2 0.01 -w3 1.0 -sigma 0.001 "%s" "%s" "%s"', ...
              P(si).Pppm, P(si).Pcentral, P(si).Pwhite);
            cat_system(cmd,opt.verb-3);
            cmd = sprintf('CAT_SurfDeform -iter 100 -isovalue 0 -w1 0.0001 -w2 0.01 -w3 1.0 -sigma 0.001 "%s" "%s" "%s"', ...
              P(si).Pppm, P(si).Pcentral, P(si).Ppial);
            cat_system(cmd,opt.verb-3);
          else
            % this function refines the surface
            Vp0 = Vmfs; Vp0.fname = spm_file(P(si).Pp0,'suffix','_tmp'); 
            spm_write_vol(Vp0, Yp0fs);
            cmd = sprintf('CAT_Surf2PialWhite "%s" "%s" "%s" "%s" "%s"', ...
              P(si).Pcentral, P(si).Pthick, Vp0.fname, P(si).Ppial, P(si).Pwhite);
            cat_system(cmd,opt.verb-3);
            delete(Vp0.fname); 
          end
          CSw1 = loadSurf(P(si).Pwhite);
          CSp1 = loadSurf(P(si).Ppial);
     
          % (3) mix simple and refined boundary surfaces
          %     0 .. 1 - mixing value from simple to enhanced boundary reconstruction
          %              The simple version often has good intensity/position values
          %              but also many (small) self-intersections. 
          mix = .5; %.5; 
          CSw0.vertices = CSw0.vertices.*(1-mix) + mix.*CSw1.vertices;
          CSp0.vertices = CSp0.vertices.*(1-mix) + mix.*CSp1.vertices;
          
          saveSurf(CSw0,P(si).Pwhite); 
          saveSurf(CSp0,P(si).Ppial); 
        end

      else % if 2 (refine but do not update thickness) or 3 (refine and update)
        %% write segmentation if necessary (and delete it later)
        Vp0 = Vmfs; Vp0.fname = spm_file(P(si).Pp0,'suffix','_tmp'); 
        spm_write_vol(Vp0, Yp0fs);
      
        % estimate pial and white surface with a limited thickness map
        facevertexcdatafs = min(6,max(eps, cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm))); 
        cat_io_FreeSurfer('write_surf_data', P(si).Pthick, facevertexcdatafs);
        spm_write_vol(Vmfs, Yp0fs);
        cmd = sprintf('CAT_Surf2PialWhite "%s" "%s" "%s" "%s" "%s"', ...
          P(si).Pcentral, P(si).Pthick, Vp0.fname, P(si).Ppial, P(si).Pwhite);
        cat_system(cmd,opt.verb-3); 
        delete(Vp0.fname); 

        % create new central surface and update thickness values
        cmd = sprintf('CAT_AverageSurfaces "%s" "%s" -avg "%s" ', ...
          P(si).Ppial, P(si).Pwhite, P(si).Pcentral);
        cat_system(cmd,opt.verb-3); 
        CS0 = loadSurf(P(si).Pcentral);

        % If we want to use the white/pial not just as show element than 
        % we have to get the updated thickness and position!
        if opt.create_white_pial > 2
          %%
          cmd = sprintf('CAT_SurfDistance -mean -thickness "%s" "%s" "%s"', P(si).Ppbt, P(si).Pcentral, P(si).Pthick);
          cat_system(cmd,opt.verb-3);
          fsthick = cat_io_FreeSurfer('read_surf_data',P(si).Pthick);  
          fsthick = spm_mesh_smooth(CS,fsthick,4);
          facevertexcdata = min( fsthick + .15 , facevertexcdata ); 
        
          % check for problems (there where issues in case of difficult geometry and thickness)
          facevertexcdata0 = max(eps, cat_surf_fun('isocolors',Yth1i,CS0.vertices,Smat.matlabIBB_mm)); 
          if all(facevertexcdata0(:)<100)
            cat_io_FreeSurfer('write_surf_data', P(si).Ppbt, min( facevertexcdata0 , facevertexcdata ) );
            clear facevertexcdata0;
          else
            % rewrite CS if something went wrong
            % - we keep the pial and white so far 
            cat_io_addwarning('PialWhiteSurfErr','Creation of sophisticated white and pial surface failed. Use simple version.',2);
            % rewrite previous output
            saveSurf( CST{thi} , P(si).Pcentral); 
            cmd = sprintf('CAT_Central2Pial "%s" "%s" "%s" 0.5',P(si).Pcentral,P(si).Ppbt,P(si).Ppial);
            cat_system(cmd,opt.verb-3);
            cmd = sprintf('CAT_Central2Pial "%s" "%s" "%s" -0.5',P(si).Pcentral,P(si).Ppbt,P(si).Pwhite);
            cat_system(cmd,opt.verb-3); 
          end
        end
      end

      if 0
        %% eval
        Vppm = spm_vol(P(si).Pppm); Yppi = spm_read_vols(Vppm); 
        cat_surf_createCS_fun('quickeval',V0,Vppm,Ymfs,Yppi,CS,P,Smat,res,opt,res.EC(si),si,time_sr,4); 
      end

    end


    
    %% final surface error
    % - with topology but ignoring possible increasemnt!
    CS = loadSurf(P(si).Pcentral);
    res.area_final(si)   = sum(cat_surf_fun('area',CS)) / 100;
    res.surferrfinal(si) = max( 0 , res.area_gt(si) - res.area_final(si) ) / res.area_gt(si) * (res.genus(si,thi)+1)/genuserr;


    % create thickness map
    if opt.thick_measure == 1
      % use central surface and thickness to estimate Freesurfer thickness
      cmd = sprintf('CAT_SurfDistance -mean -thickness "%s" "%s" "%s"', P(si).Ppbt, P(si).Pcentral, P(si).Pthick);
      cat_system(cmd,opt.verb-3);
    else
      % just write PBT thickness
      cat_io_FreeSurfer('write_surf_data', P(si).Pthick, facevertexcdata);
    end


    % correct thickness based on folding pattern, but smaller thickness values 
    % are corrected less strongly than larger thickness values
    if opt.foldingcorrection
      cmd = sprintf('CAT_SurfCorrectThicknessFolding -slope 1.0 -max "%f" "%s" "%s" "%s"', ...
        opt.thick_limit, P(si).Pcentral, P(si).Pthick, P(si).Pthick);
      cat_system(cmd,opt.verb-3);
    end


    if isscalar(opt.surf)
      % only for test visualization  
      fprintf('%5.0fs',etime(clock,stime)); 
      Vppm = spm_vol(P(si).Pppm); Yppi = spm_read_vols(Vppm); 
      cat_surf_createCS_fun('quickeval',V0,Vppm,Ymfs,Yppi,CS,P,Smat,res,opt,res.EC(si),si,time_sr,4); 
      return
    end
        

    %% skip that part if a prior image is defined
    if ~useprior && ~opt.skip_registration
      % spherical surface mapping of the final corrected surface with minimal areal smoothing
      stime = cat_io_cmd('  Spherical mapping with areal smoothing','g5','',opt.verb,stime); 
      cmd = sprintf('CAT_Surf2Sphere "%s" "%s" %d',P(si).Pcentral,P(si).Psphere, ... 
        round( (5 + round( sqrt( size(CS.faces,1) / 10000 ) + 1 ))) ); % use dynamic filter interations depending on surface size (300k with value 10)
      cat_system(cmd,opt.verb-3);
      
      % spherical registration to fsaverage template
      stime = cat_io_cmd('  Spherical registration','g5','',opt.verb,stime); 
      cmd = sprintf('CAT_SurfWarp -steps 2 -avg -i "%s" -is "%s" -t "%s" -ts "%s" -ws "%s"', ...
        P(si).Pcentral,P(si).Psphere,P(si).Pfsavg,P(si).Pfsavgsph,P(si).Pspherereg);
      cat_system(cmd,opt.verb-3);
      if opt.verb, fprintf('%5.0fs',etime(clock,stime)); end %#ok<*DETIM>
    end  


    if debug || cat_get_defaults('extopts.expertgui')>1, fprintf('\n'); end
    if opt.thick_measure == 1
      % use central surface and thickness to estimate Freesurfer thickness
      cmd = sprintf('CAT_SurfDistance -mean -thickness "%s" "%s" "%s"',P(si).Ppbt,P(si).Pcentral,P(si).Pthick);
      cat_system(cmd,opt.verb-3);
    
      % apply upper thickness limit
      % here the 5 mm thickness limit of FreeSurfer might be better rather then our 6 mm 
      facevertexcdatafs = max(eps, cat_io_FreeSurfer('read_surf_data',P(si).Pthick));  
    else
      % otherwise simply use the original values of the PBT map  
      % WARNING: 
      %   The values of the ?h.pbt.* files are used to estimate further 
      %   surfaces and therefore corrected for self-intersections!
      %   We use here the orignal PBT values as they are less depending on 
      %   local highly individual features.  
      facevertexcdatafs = facevertexcdata; 
    end
    cat_io_FreeSurfer('write_surf_data', P(si).Pthick, min(opt.thick_limit,facevertexcdatafs) );
    fprintf('\n');
    

    % final surface evaluation 
    Vppm = spm_vol(P(si).Pppm); Yppi = spm_read_vols(Vppm); 
    res.(opt.surf{si}).createCS_final = cat_surf_fun('evalCS', ...
      loadSurf(P(si).Pcentral), cat_io_FreeSurfer('read_surf_data',P(si).Ppbt), cat_io_FreeSurfer('read_surf_data',P(si).Pthick), ...
      Ymfs, Yppi, P(si).Pcentral, Smat.matlabIBB_mm, debug + (cat_get_defaults('extopts.expertgui')>1), 1);
    clear Yppi; 

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
  cat_io_cmd('  ','g5','',opt.verb);
  fprintf('%5ds\n',round(etime(clock,cstime)));


  % calculate surface quality parameters for all surfaces
  res = cat_surf_createCS_fun('addSurfaceQualityMeasures',res,opt);
 
  % print final stats
  cat_surf_createCS_fun('evalProcessing',res,opt,P,V0); 

end
%=======================================================================
function saveSurf(CS,P)
  save(gifti(struct('faces',CS.faces,'vertices',CS.vertices)),P,'Base64Binary'); %#ok<USENS>
end
%=======================================================================
function CS1 = loadSurf(P)
  % add 5s + 10s because sometimes surface is not yet ready to read (i.e. NAS)
  if ~exist(P,'file')
    pause(5)
    if ~exist(P,'file')
      pause(10+rand(1))
      try
        CS = gifti(P);
        CS1.vertices = CS.vertices; CS1.faces = CS.faces; 
        if isfield(CS,'cdata'), CS1.cdata = CS.cdata; end
      catch
        error('Surface file %s could not be read due to previous processing errors.',P);
      end
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
%createOutputFileStructures. Create basic structure V and transformations Smat.
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
function [Yp0f,Ymf] = fillVentricle(Yp0,Ym,Ya,YMF,vx_vol,opt)
%fillVentricle. Simple filling of ventricles by a closing around a mask.

  NS  = @(Ys,s) Ys==s | Ys==s+1; 
  LAB = cat_get_defaults('extopts.LAB');
  
  %% extend filling map 
  %  Relevant in case of larger ventricular GM in some thickness phantoms
  %  that caused severe topology errors those corrected often blurred sulci.
  %  However, this can of course also close thin slulci!
  if opt.SRP>0 
    % find simple holes 
    Yholes = Yp0<2.75/3 & ~cat_vol_morph(Yp0<2.75/3,'lo',1,vx_vol) | Yp0<2.25/3 & ~cat_vol_morph(Yp0<2.25/3,'l'); 
    Yholes = smooth3(Yholes)>.4  &  ~cat_vol_morph( cat_vol_morph( Yp0<2.25/3 ,'do',2),'d'); 
    Ync    = cat_vol_morph( NS(Ya,LAB.HC) | NS(Ya,LAB.PH), 'dd',5,vx_vol);

    % go to lower resolution for laplace bottleneck closing
    res = 1.2; 
    [Yfm,resR] = cat_vol_resize(smooth3(YMF & ~Ync ), 'reduceV', vx_vol, res, 16, 'median');  
    Yp0r       = cat_vol_resize(max(1/3 * cat_vol_morph(Yp0>1/3,'lc'),Yp0), 'reduceV', vx_vol, res, 16, 'median');  
    Yholesr    = cat_vol_resize(Yholes, 'reduceV', vx_vol, res, 16, 'max');  clear Yholes
    
    % define the regions
    Yfm = round( (Yfm>.5 | Yholesr>.25 ) + 1); clear Yholes; 
    Yfm(cat_vol_morph( Yp0r>2.75/3 ,'do',1,resR.vx_volr) ) = 1.9; 
    Yfm(Yp0r<.5/3) = 0; %clear Yp0r;
   
    % applay laplace filter and back to original resolution 
    Yfm = cat_vol_laplace3R(single(Yfm), Yfm==1, .005*mean(resR.vx_volr)); 
    Yfm = cat_vol_resize(Yfm,'dereduceV',resR);

    % final mask
    Yfm = Yfm>1.9 | cat_vol_morph( cat_vol_morph( cat_vol_morph(smooth3(Yfm)>1.899,'do',3,vx_vol) & Yp0>2.25/3,'ldc',1,vx_vol) ,'e',2,vx_vol); 
    Yfm = cat_vol_morph(Yfm,'ldc',1,vx_vol) & ~Ync; 
    Yas = cat_vol_smooth3X( YMF | NS(Ya,LAB.BG) | NS(Ya,LAB.TH) | NS(Ya,LAB.MB)  | NS(Ya,LAB.VT) ,8); % deep structures
    Yfm = Yfm & (cat_vol_smooth3X(Yp0 | Yas*3,8)>2/3) & (Yas>.25) ; % assure that it is generally closer to WM (avoid occipital blurring)
    
    % combine closing maps
    YMF = max(YMF, (0.01/opt.wmnoise).* cat_vol_morph( cat_vol_morph( Yp0>2.75/3 ,'de',1.5,vx_vol) | Yfm,'dc',1,vx_vol) ); 
  end


  %% simple filling by the YMF mask  
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
function Yp0 = closeWMHandPVS(Yp0, Ya, NS, LABCB, wmnoise, vx_vol)
%closeWMHandPVS. Correct strong WMHs and PVSs.
%  Eg.  NISALS_UTR_SP30T_als3_T1w  dataset with Swiss cheese WM and wmnoise ~ 0.04.
%  Blood vessels should be already removed here and we should also not apply this in non T1w data! 
%  This can cause blurring of sulci (eg. INDI_HC_NWK_sub13411). 
    
  Yp0o = Yp0; vxs = mean(vx_vol); 
  
  %% first we close the GM and the GM-WM segment
  Yp0    = Yp0o;
  Yp0hgm = single( 1 + (Yp0>2.25/3) ); Yp0hgm(cat_vol_smooth3X(Yp0)<1.75/3 & Yp0<2.25) = 0;
  Yp0hgm = cat_vol_laplace3R(Yp0hgm,Yp0hgm==1,0.25); % small values blurr small sulci too
  Yp0hgm = Yp0hgm & (cat_vol_smooth3X(Yp0,4/vxs)<2.125/3); % avoid blurres of sulci
  Yp0hgm = max( min(2/3,smooth3(Yp0hgm)*2/3) , (smooth3(Yp0hgm)/2).^2 ); % first part close pure GM holes, second fills GM+tissue
  Yp0    = max(Yp0,  ~NS(Ya, LABCB) .* Yp0hgm); % prepare correction 
  Yp0    = max(Yp0o, min(1,Yp0 + smooth3( max(0,smooth3(Yp0-Yp0o) * 4 -.25) .* (Yp0o<1.9/3)).^.5 )); % corrected CSF was WM 
  
  %% next we fix the WM
  Yp0hgm = single( 1 + (Yp0>2.75/3) ); Yp0hgm(cat_vol_smooth3X(Yp0)<2.25/3 & Yp0<2.75) = 0;
  Yp0hgm = cat_vol_laplace3R(Yp0hgm,Yp0hgm==1,0.1); 
  Yp0hgm = max( Yp0hgm, cat_vol_morph(Yp0hgm,'ldc',1)); 
  Yp0    = max(Yp0,min(1,smooth3(Yp0hgm-1)));

  % final mixing beside the cerebellum
  corr   = min(1,wmnoise*100 - 3) .* (1 - NS(Ya, LABCB));
  Yp0    = Yp0o .* (1-corr) + (corr) .* Yp0; 
end
%==========================================================================
function Vppm = exportPPmap( Yppi, Vmfs, opt)
%exportPPmap. Write position map. 
  
  if opt.interpV == opt.reconres
    Vppm           = Vmfs;
    Ytx            = Yppi;
  else
    Vpp            = Vmfs; 
    Vpp.pinfo      = repmat([1;0], 1,size(Yppi,3));
    Vpp.dat(:,:,:) = single(Yppi); 
    Vpp.dt(1)      = 16; 
    Vppr           = Vpp; 

    imat           = spm_imatrix(Vpp.mat); 
    Vppr.dim       = round(Vpp.dim .* opt.interpV./opt.reconres);
    imat(7:9)      = opt.reconres .* sign(imat(7:9));
    Vppr.mat       = spm_matrix(imat); clear imat;  
    Vppr.fname     = spm_file(Vmfs.fname,'suffix','_pptmp');
    
    spm_smooth(Yppi,Yppi, repmat( max(0.2,min(1,opt.reconres/opt.interpV)), 1,3) );
    
    Vpp.dat(:,:,:) = Yppi; 
    [Vppm,Ytx]     = cat_vol_imcalc(Vpp,Vppr,'i1',struct('interp',5,'verb',0,'mask',-1));
    Vppm.pinfo     = Vmfs.pinfo;
  end
 
  if isfield(Vppm,'dat'), Vppm = rmfield(Vppm,{'dat'}); end
  if exist('Vppr','var') && exist(Vppr.fname,'file'), delete(Vppr.fname); end
  if exist('Vppm','var') && exist(Vppm.fname,'file'), delete(Vppm.fname); end
  spm_write_vol(Vppm, Ytx );
  Vppm = spm_vol(Vppm.fname); 
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
function Yp0 = myelincorrection(Yp0,vx_vol,opt,P,Vmfs,si,quick) 
%myelincorrection. Correct high intensity GM range if suitable. 

  % quick estimation of the cortical thickness
  if ~exist('quick','var'), quick = 1; end
 
  if quick 
  % quicker but sufficient version (~10s) 
    % parameters 
    opt.verb      = 0; 
    opt.levels    = 2; 
   
    opt.pvet0 = estimatePVEsize( Yp0 , 0); 
    opt.pvet  = max(mean(vx_vol),min(2 / mean(vx_vol),opt.pvet0)); 

    % from cat-vol_pbtsimpleCS4
    [Ycd, Ywd] = cat_vol_cwdist(Yp0, opt);
    % set unvisited points as holes/handels
    Yp0(isinf(Ywd)) = 1; Yp0(isinf(Ycd)) = 3; 
    Ywd(isinf(Ywd)) = 0; Ycd(isinf(Ycd)) = 0;
  
    % projection-based thickness mapping
    Ygmt0 = cat_vol_pbtp( round(Yp0) , Ywd, Ycd);  
    Yp0(Ygmt0(:)>10e10) = 1; 
    Ymsk = Ygmt0(:)>10e10; 
    Ywd(Ymsk) = 0; Ycd(Ymsk) = 0; Ygmt0(Ymsk) = 0; clear Ymsk

    Ygmt0 = cat_vol_approx(Ygmt0,'rec',2); 
  
    % reestimation of the CSF distance 
    Ypp   = min(1,min(max(0,Ygmt0-Ywd),Ycd) ./ max(eps,Ygmt0)); Ypp(Yp0>2.5 & Ypp==0) = 1; 
  else
    % alternative save but slow version (45s)
    %[Vmfs,Smat] = createOutputFileStructures(V,V0,resI,BB,opt,mridir,ff,si); 
    Vmfs.dt = [16 1]; spm_write_vol(Vmfs, Yp0 );
    cmd = sprintf('CAT_VolThicknessPbt  -correct-voxelsize 0   -median-filter 2   -downsample 0 "%s" "%s" "%s"', Vmfs.fname, P(si).Pgmt, P(si).Pppm);
    cat_system(cmd,opt.verb-3);
    Vgmt0  = spm_vol(P(si).Pgmt); Ygmt0 = spm_read_vols(Vgmt0); 
    Ygmt0  = max(0,Ygmt0 - 0.56*mean(opt.interpV) ); 
    Vpp    = spm_vol(P(si).Pppm); Ypp  = spm_read_vols(Vpp); 
    Ycd    = Ygmt0 .* Ypp; 
  end

  % reestimation of the CSF distance 
  Ycdc2 = cat_vbdist( single( max(Yp0<=1, 1 - Ycd - Ypp) ), true(size(Ycd)) );
  Ycdc2(Ycdc2 > 6 / mean(vx_vol)) = 0; 
 
  % estimate the full tissue thickness (we needed the GM thickness and WM to reconstruct the sulcus)
  Ybmt  = cat_vol_pbtp( min(3,4 - min(2,Yp0)), Ycdc2, Ycdc2*inf); Ybmt(Ybmt>1000) = 0; 
  Ybmt  = cat_vol_approx(Ybmt); 

  % estimate correction area
  medgmt = median(Ygmt0(:)); 
  try iqrgmt = iqr(Ygmt0(:)); catch, iqrgmt = std(Ygmt0(:)); end
  YenoughWM  = Ycdc2>0 & Ycdc2 < Ybmt - 1.5;
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
    
  def.levels = 2; 
  def.rangeE = 0.4; % extimation limit extention - smaller better
  def.range  = 0.4; % addition range - larger better but limited by interpolation artifacts
  def.vxs    = 1;

  opt  = cat_io_checkinopt(opt,def); 
  vxc  = .5 - opt.range; % the more PVE the less we need to correct for
  
  range  = min(.4,opt.range);    % defined by the PVE boundary of the AMAP, i.e range of .25 - .75  
  rangeE = min(.4,opt.rangeE);   % defined by the PVE boundary of the AMAP, i.e range of .25 - .75  
  
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
    [Ycdl,YI]  = cat_vbdist3(single(Yp0 < ( 1.5 - offset)), YMC ); 
    Ycdl       = (Ycdl - min(vxc,max(0,- (Yp0(YI) - (1.5-offset))*opt.pvet ))) .* (YMM & YMC); 
    [Ycdh,YI]  = cat_vbdist3(single(Yp0 < ( 1.5 + offset)), YMC ); 
    Ycdh       = (Ycdh - min(vxc,max(0,- (Yp0(YI) - (1.5+offset))*opt.pvet ))) .* (YMM & YMC); 
    Ycd        = Ycd + .5/hss .* Ycdl  +  .5/hss .* Ycdh; 

    % WM distances
    [Ywdl,YI]  = cat_vbdist3(single(Yp0 > ( 2.5 - offset)), YMM );
    Ywdl       = (Ywdl - min(vxc,max(0, (Yp0(YI) - (2.5-offset))*opt.pvet ))) .* (YMM & YMC); 
    [Ywdh,YI]  = cat_vbdist3(single(Yp0 > ( 2.5 + offset)), YMM ); 
    Ywdh       = (Ywdh - min(vxc,max(0, (Yp0(YI) - (2.5+offset))*opt.pvet ))) .* (YMM & YMC); 
    Ywd        = Ywd + .5/hss .* Ywdl  +  .5/hss .* Ywdh;
  end

  % endpoint PVE correction 
  Ycd = Ycd - min(.5,max(-.0, (Yp0>2.5-0*opt.rangeE & Yp0<2.5+opt.rangeE) .* ( (Yp0-2.5)*opt.pvet*1 + 0*opt.vxs/2) ));
  Ywd = Ywd + min(.0,max(-.5, (Yp0>1.5-opt.rangeE & Yp0<1.5+0*opt.rangeE) .* ( (Yp0-1.5)*opt.pvet*1 - 0*opt.vxs/2) ));
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
function [Yp0s,Yms] = cat_surf_createCS_denoise(Yp0,Ym, wmnoise)
%cat_surf_createCS_denoise. Apply ORNLM denoising. 
% The ORNLM filter modifies the overall intensities and we have to rescale the image!

  if wmnoise == 0, Yp0s = Yp0; Yms = Ym; return; end

  Yp0x   = smooth3(Yp0);
  Yp0s   = cat_ornlm(Yp0,1,1,min(1/4, wmnoise ));
  Yms    = cat_ornlm(Ym ,1,1,min(1/4, wmnoise ));
  Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));

  % sharpening
  fs   = min(.5,wmnoise ); 
  Yms1 = smooth3(Yms); 
  Yms2 = cat_vol_smooth3X(Yms,2); 
  Yms  = Yms  + fs*smooth3(Ym  - Yms1) + fs*smooth3(Ym  - Yms2);
  Yp0s = Yp0s + fs*smooth3(Yp0 - Yms1) + fs*smooth3(Ym  - Yms2); 
  clear Yms1 Yms2; 

  % correct intensities Yp0
  Tth.T3thx = [ ...
    median(Yp0s( Yp0toC(Yp0x(:)*3,1)>.5))/2, ...
    median(Yp0s( Yp0toC(Yp0x(:)*3,1)>.5)), ...
    median(Yp0s( Yp0toC(Yp0x(:)*3,2)>.5)), ...
    median(Yp0s( Yp0toC(Yp0x(:)*3,3)>.5)) 4/3]; 
  Tth.T3th  = 0:1/3:4/3; 
  Yp0s = max( cat_vol_morph(cat_vol_morph(Yp0x>2.75/3,'de'),'ldc'), cat_main_gintnormi(Yp0s/3,Tth));
  Yp0s = min( 1 - 2/3*cat_vol_morph(cat_vol_morph(Yp0x<1.25/3,'de'),'ldc') , Yp0s);
  
  % correct intensities Ym
  Tth.T3thx = [ ...
    median(Yms( Yp0toC(Yp0x(:)*3,1)>.5)), ...
    median(Yms( Yp0toC(Yp0x(:)*3,1)>.5)), ...
    median(Yms( Yp0toC(Yp0x(:)*3,2)>.5)), ...
    median(Yms( Yp0toC(Yp0x(:)*3,3)>.5)) 4/3]; 
  Tth.T3th  = 0:1/3:4/3; 
  Yms = max( cat_vol_morph(Ym/3,'de'), cat_main_gintnormi(Yms/3,Tth));
  Yms = min( 1 - 2/3*cat_vol_morph(cat_vol_morph(Yp0x<1.25/3,'de'),'ldc') , Yms);

end
function [Yms,noise] = denoiseSanlm(Ym,vx_vol,fast,sharpen)

  if ~exist('vx_vol','var'),  vx_vol  = ones(1,3); end
  if ~exist('fast','var'),    fast    = 1; end
  if ~exist('sharpen','var'), sharpen = 1; end
 
  Yms = Ym; cat_sanlm(Yms,1,3); 
  
  if fast
    [Ymr,red] = cat_vol_resize(Yms,'reduceV',1,2,16,'meanm'); 
    Ymrs = Ymr + 0; cat_sanlm(Ymrs,1,3); Ymr = Ymr - Ymrs;
    Ymr1 = cat_vol_resize(Ymr+1,'dereduceV',red)-1; 

    [Ymr,red] = cat_vol_resize(Yms(2:end,2:end,2:end),'reduceV',1,2,16,'meanm'); 
    Ymrs = Ymr + 0; cat_sanlm(Ymrs,1,3); Ymr = Ymr - Ymrs;
    Ymr2 = Ymr1; Ymr2(2:end,2:end,2:end) = cat_vol_resize(Ymr+1,'dereduceV',red)-1; 
    
    Yms2 = Yms - Ymr1*4 - Ymr2*4;
    clear Ymr1 Ymr2
  else
    Yms2 = Yms;        
    for x = 1:2
      for y = 1:2
        for z = 1:2
          [Ymr,red] = cat_vol_resize(Yms(x:end,y:end,z:end),'reduceV',1,2,16,'meanm'); 
          Ymrs = Ymr + 0; cat_sanlm(Ymrs,1,3); Ymr = Ymr - Ymrs;
          Yms2(x:end,y:end,z:end) = Yms2(x:end,y:end,z:end) - (cat_vol_resize(Ymr + 1,'dereduceV',red) - 1); 
          clear Ymr Ymrs
        end
      end
    end
  end
  Ymsk  = Yms ~= 0; 
  cat_sanlm(Yms2,1,3); 
  Yw    = cat_vol_approx( (Yms>1.5/3 ) .* ( log10(Yms*10) ./ log10(Yms2*10))); 
  Yw    = Yw ./ cat_stat_nanmean(Yw(Yms(:)>1.5/3)); 
  Yms2a = Yms*.5 + .5*Yms2 ./ Yw;

  noise = cat_stat_nanmean( (Ym(Ymsk(:)) - Yms2(Ymsk(:))).^2 )^.5; 
  fs    = min(.5,noise/2); 
  Yms2  = Yms2a + fs * 2 * smooth3(Ym - cat_vol_smooth3X(Yms2a,2)) + fs * 4 * smooth3(Ym - cat_vol_smooth3X(Yms2a,4)); clear Yms2a

  %% use lowres filtering only for highres data with smooth transition
  vxmix = .5; % * max(0,min(1,mean(vx_vol)-.5)); % .5 to 1.
  Yms   = Yms.*(vxmix) + (1-vxmix).*Yms2; clear Yms2 

  cat_sanlm(Yms,1,3); 

  Ymsk  = Yms ~= 0; 
  noise = cat_stat_nanmean( (Ym(Ymsk(:)) - Yms(Ymsk(:))).^2 )^.5; 
  clear Ymsk; 

  %% sharpening
  if sharpen && noise > .001
    fs    = min(.5,noise); 
    Yms1  = smooth3(Yms); 
    Yms2  = cat_vol_smooth3X(Yms,2); 
    Yms   = Yms  + fs*smooth3(Yms  - Yms1) + fs*2*smooth3(Yms  - Yms2); 
    clear Yms1 Yms2
  end
end
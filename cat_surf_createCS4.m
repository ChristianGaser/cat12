function [Yth,S,Psurf,res] = cat_surf_createCS4(V,V0,Ym,Ya,YMF,Ytemplate,Yb0,opt,job)
% ______________________________________________________________________
% Surface creation and thickness estimation.
%
% [Yth,S,Psurf] = cat_surf_createCS3(V,V0,Ym,Ya,YMF,Ytemplate,opt)
%
% Yth    .. thickness map
% S      .. structure with surfaces, like the left hemisphere, that contains
%           vertices, faces, GM thickness (th1)
% Psurf  .. name of surface files
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
  
  % Turn off gifti data warning in gifti/subsref (line 45)
  %   Warning: A value of class "int32" was indexed with no subscripts specified. 
  %            Currently the result of this operation is the indexed value itself, 
  %            but in a future release, it will be an error. 
  warning('off','MATLAB:subscripting:noSubscriptsSpecified');
  cstime = clock;

  % test-variables that should be (partially) removed later
  use_cat_vol_pbtsimple   = 1;
  use_fullpreprocess      = 1;
  use_high_resolution     = 0; % force .5 mm surface creation and processing otherwise use input resolution  
  skip_registration       = 1; % skip spherical registration for quick tests
  create_white_pial       = 1; % uses only the quick WM and Pial surface estimation 

  myelinCorrection        = .4; % .25 - sight correction, 1 - maximum correction
  deepFilling             = 1; % fill deep holes that were not closed by the YMF map or WMHs and PVSs
  sharpenGyri             = 1; % refine thin WM gyri 
  dynamicThresholds       = 1; % estimate athropy dependend smoothing and position creation thresholds
                               % (i.e. closer to pial in case less WM, closer to WM in case of much WM) 
  setcut2zero             = 1; % works but result in worse values as could be expected
  create_surf             = 1; % V1: fancy-strange topology correction, dynamic Ypp threshold  
                               % V2: morph-based topology correction, fixed Ypp threshold (worse values, thick but more broken thin gyri)
  deformsurf              = 1; % V0: no, need stronger surface smoohting ( bad surface values - need .5 Ypp threshold )
                               % V1: old version - better values but self-intersections issues
                               % V2: new version - worse  values but less self-intersection issues (similar values - surface show steps)
  if opt.SRP<2 
    myelinCorrection      = 0; 
    deepFilling           = 0; 
    sharpenGyri           = 0;
    dynamicThresholds     = 0;
    setcut2zero           = 0; 
    create_surf           = 1; 
    deformsurf            = 1; 
  end
  if opt.SRP>2
    create_surf           = 2; 
    deformsurf            = 2; 
  end

  % == dynamic parameters ==
  % The basic ideas is that the surface creation of atrophic braings should 
  % profit by lower position thresholds to reduce risk of breaking thin 
  % structures, whereas less atrophic brains rquire higher thresholds to 
  % avoid blurring of thin sulci. 
  %
  % large CSF areas (ventricle/external CSF) are biasing the estimation of CSF in the sulci
  if dynamicThresholds
    YlargeCSF               = cat_vol_morph(Yb0 & Ym<1.5/3,'o',2); 
    % extended WM vs. extended sulcal CSF volume to quantify thin sulci in relation to the gyri
    createsurface_ppth0     = sum(Yb0(:) & ~YlargeCSF(:) & ~YMF(:) & Ym(:)>2.25/3) / ...
                              sum(Yb0(:) & ~YlargeCSF(:) & ~YMF(:) & Ym(:)<1.75/3);
    % close-opening for surface creation depending on atrophy 
    createsurface_scl0      = sum(Yb0(:) & ~YlargeCSF(:) & ~YMF(:) & Ym(:)<2.75/3) / ...
                              sum(Yb0(:) & ~YlargeCSF(:) & ~YMF(:) & Ym(:)>2.75/3);
    createsurface_ppth      = max(.5,min(1.5, createsurface_ppth0 )); 
    createsurface_scl       = max(.9,min(1.5, createsurface_scl0 ));
    % blood vessels (BV) requires higher (WM-closer) position threshold and stronger opening
    createsurface_BV        = sum(NS(Ya(:),7)) / sum(Yb0(:)); 
    % prefere open if many blood vessels
    createsurface_scl       = max(0.9,createsurface_scl - createsurface_BV*30); 
    createsurface_ppth      = createsurface_ppth + createsurface_BV * 10; 
    % adapt smoothing parameters for deformation and blood vessels
    createsurface_prefwhm   = -4 * (2 -   (deformsurf>0)); 
    createsurface_postfwhm  =  2 * (3 - 2*(deformsurf>0)); % lower resolution / lower PVE need more smoothing (required especially without deformsurf)    
    createsurface_prefwhm   = createsurface_prefwhm  * max( 1,min(2,createsurface_BV / .005 * createsurface_scl0));
    createsurface_postfwhm  = createsurface_postfwhm * max( 1,min(2,createsurface_BV / .005 * createsurface_scl0));
    % isosurface threshold
    th_initial = 0.50 * min(1.333,max( 0.666, createsurface_ppth / createsurface_scl ) );
  
    if opt.verb
      fprintf('\n  bv=%0.02f, scl=%0.2f(%0.2f), ppth=%0.2f(%0.2f), fwhm=%0.1f/%0.1f, thi=%.2f ', ...
        createsurface_BV, createsurface_scl, createsurface_scl0, ...
        createsurface_ppth, createsurface_ppth0, ...
        createsurface_prefwhm, createsurface_postfwhm, ...
        th_initial);
    end
  else
  % fixed parameters
    createsurface_scl       =  .9; 
    createsurface_prefwhm   = -.2; 
    createsurface_postfwhm  =  2; 
    th_initial              =  0.55; % .6 shows more blue gyral breaks 
    fprintf('\n scl=%0.2f, fwhm=%0.1f/%0.1f, thi=%.2f ', ...
      createsurface_scl, createsurface_prefwhm, createsurface_postfwhm, th_initial); 
  end
  %Ywmhs = Ym<2.5/3 & cat_vol_morph(Ym>2.5/3,'lc'); rWMHV = sum(Ywmhs(:))./sum(Ym(:)>2.5/3); clear Ywmhs; 
  

  % set debugging variable
  dbs   = dbstatus; debug = 1; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end
  S = struct();

  % set defaults
  if ~exist('opt','var'), opt = struct(); end                 % create variable if not exist
  vx_vol        = sqrt(sum(V.mat(1:3,1:3).^2));               % further interpolation based on internal resolution 
  def.verb      = cat_get_defaults('extopts.expertgui');      % 0-none, 1-minimal, 2-default, 3-details, 4-debug
  def.surf      = {'lh','rh'};                                % surface reconstruction setting with {'lh','rh','cb'} 
  % There is a new SPM approach spm_mesh_reduce that is maybe more robust. 
  % Higher resolution is at least required for animal preprocessing that is given by cat_main.
  def.LAB                 = cat_get_defaults('extopts.LAB');  % brain regions 
  def.SPM                 = 0;                                % surface-reconstration based on SPM segmentation input (see cat_main)
  def.pbtlas              = 0;                                % myelination correction option (in development - not working correctly in all data, RD201907)  
  def.pbtmethod           = 'pbtsimple';                      % projection-based thickness (PBT) estimation ('pbt2x' (with minimum setting), 'pbt2', or 'pbtsimple')
  def.sharpenCB           = 1;                                % sharpening function for the cerebellum (in development, RD2017-2019)
  % RD20250306: Tfs has large issues currently with some corrected defects
  def.thick_limit         = 5;                                % 5mm upper limit for thickness (same limit as used in Freesurfer)
  def.thick_measure       = 0;                                % 0-PBT; 1-Tfs (Freesurfer method using mean(TnearIS,TnearOS))
  
  def.fsavgDir            = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces'); 
  def.close_parahipp      = 0;
  def.outputpp.native     = 0;  % output of Ypp map for cortical orientation in EEG/MEG 
  def.outputpp.warped     = 0;
  def.outputpp.dartel     = 0;
  def.skip_registration   = skip_registration;

  % options that rely on other options
  opt                     = cat_io_updateStruct(def,opt); clear def; 
  opt.vol                 = any(~cellfun('isempty',strfind(opt.surf,'v')));   % only volume-based thickness estimation  
  opt.surf                = cat_io_strrep(opt.surf,'v','');                   % after definition of the 'vol' varialbe we simplify 'surf'
  opt.interpV             = max(0.1,min([opt.interpV,1.5]));                  % general limitation of the PBT resolution
  
  % for testing only  
  skip_registration = debug | opt.skip_registration;
  
  % Another parameter to control runtime is the accuracy of the surface
  % deformation. As far as we primary adapt the mesh resolution above, it 
  % is useful to use sqrt values rather than linear values to avoid square
  % processing times for higher quality levels. Otherwise, we can simple 
  % avoid changes here
  
  % function to estimate the number of interactions of the surface deformation: d=distance in mm and a=accuracy 
  QMC    = cat_io_colormaps('marks+',17);
  color  = @(m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
  rate   = @(x,best,worst) min(6,max(1, max(0,x-best) ./ (worst-best) * 5 + 1));
  
  if exist('job','var')
    [mrifolder, reportfolder, surffolder, labelfolder] = cat_io_subfolders(V0.fname,job);
  else
    [mrifolder, reportfolder, surffolder, labelfolder] = cat_io_subfolders(V0.fname);
  end
  
  % get original filename without 'n'
  [pp0,ff] = spm_fileparts(V0.fname);
  
  % correct '../' parts in directory for BIDS structure
  [stat, val] = fileattrib(fullfile(pp0,surffolder));
  if stat
    pp0_surffolder = val.Name;
  else
    pp0_surffolder = fullfile(pp0,surffolder);
  end

  if ~exist(fullfile(pp0_surffolder),'dir'), mkdir(fullfile(pp0_surffolder)); end

  % apply the modified mask from gcut
  % for non-gcut approaches or inverse weighting Yb0 only contains ones
  Ym = Ym.*(Yb0>0.5);

  % noise reduction for higher resolutions (>=1 mm full correction, 1.5 mm as lower limit)
  % (added 20160920 ~R1010 due to severe sulcus reconstruction problems with 1.5 Tesla data)
if use_fullpreprocess
  Yms = Ym + 0; cat_sanlm(Yms,3,1);
  mf  = min(1,max(0,3-2*mean(vx_vol,2))); 
  Ym  = mf * Yms  +  (1-mf) * Ym;
  clear Yms;
end
  %% filling
  Ymf  = max(Ym,min(1,YMF & ~NS(Ya,opt.LAB.HC) & ~( cat_vol_morph( NS(Ya,opt.LAB.HC),'dd',2,vx_vol)))); 
  
  if deepFilling || sharpenGyri
    % estimate CC to avoid changes in that range
    Ycc  = cat_vol_morph( ...
             cat_vol_morph(mod(Ya,2)==0 & Ya>0,'d') & ...
             cat_vol_morph(mod(Ya,2)==1,'d') & Ymf>2.5/3 ,'l'); 
    Ycc  = smooth3(cat_vol_morph( Ycc ,'dd', 15) & Ymf<.95)>.5; 
  end

  if deepFilling
    % filling of deep wholes 
    [~,D] = cat_vol_downcut( single(cat_vol_morph( Ymf<2/3,'ldo',1)) , 1 - Ymf.^10, .01); 
    Ymf1  = min(1,Ymf + max(0,smooth3( D==max(D(:)) ) .* (1 - NS(Ya,opt.LAB.CB) - Ycc)));
    Ymf1(Ymf<.95 & ~cat_vol_morph(Ymf<0.95,'ldo',1.5,vx_vol)) = 1; 
  else
    Ymf1 = Ymf; 
  end

  %%
  if sharpenGyri
    %% open to avoid BVs
    
    % define WMH/PVS correction map by high CSF and low WM proportion
    Yp0toC = @(c) 1-min(1,abs(max(1,Ymf1*3)-c));
    %Ycsf = Yp0toC(1); spm_smooth(Ycsf,Ycsf,8./vx_vol); 
    %Ygm  = Yp0toC(2); spm_smooth(Ygm ,Ygm ,2./vx_vol);
    %Ywm  = Yp0toC(3); spm_smooth(Ygm ,Ygm ,4./vx_vol);
    six = [1 2 4]; sd=numel(six)+1; Ywm = Yp0toC(3)/sd;  
    for si = six, Ywms = Yp0toC(3); spm_smooth(Ywms ,Ywms ,4./vx_vol); Ywm = Ywm + Ywms/sd; end
    %Yvt  = single(cat_vol_morph(NS(Ya,opt.LAB.VT),'dd',5,vx_vol)); spm_smooth(Yvt ,Yvt ,16./vx_vol);
    %Ycsf = max(Ycsf,min(2,Yvt*8)); clear Yvt; 
    
    % thickness and position estimate to correct (too) thick regions 
    [Ygmt,Ypp] = cat_vol_pbtsimple(Ymf1*3,vx_vol,struct('supersimple',1));
    Ygmt = cat_vol_approx(Ygmt .* (Ygmt<10)); 

    % potential region of change - partial volume voxels with real GM but not WM neighbor (ie. thin gyri)
    Ypve  = min(3,cat_vol_smooth3X(Ymf1>2/3 & Ymf1<2.5/3 & ~(cat_vol_morph( Ymf1<2.25/3 ,'d') & cat_vol_morph( Ymf1>2.75/3 ,'d')),4) * 4); 
   
    %% potential change
    for xii=1:1
      Ymgc0 = cat_vol_morph(Ymf1 .* (Ymf1>0/3),'gc',1,vx_vol); 
      Ymgo0 = cat_vol_morph(Ymf1 .* (Ymf1>0/3),'go',1,vx_vol); 
      Ymgc  = cat_vol_morph(Ymf1 .* (Ymf1>2/3),'gc',1,vx_vol); 
      Ymgo  = cat_vol_morph(Ymf1 .* (Ymf1>2/3),'go',1,vx_vol); 
      Ymc   = cat_vol_median3( Ymgc - Ymgo + Ymgc0 - Ymgo0 , Ymf>0 & Ymf<1,Ymf>=0,.2); 
      clear Ymgc0 Ymgo0 Ymgc Ymgo
  
      % final correction
      Yecb = 1 * NS(Ya,opt.LAB.CB);
      if xii==1, Ymf2 = Ymf1; end
      Ymf2 = max(Ymf2,min(2.7/3,Ymf2 + ( Ymc .* 3*(1-th_initial) .* smooth3( Ywm .* Ypp.^.1 .* (Ygmt+Yecb).^2 .* Ypve).^2  ).^2));
      Ymf2 = max(Ymf2,min(3.2/3,Ymf2 + ( Ymc .* 1*(1-th_initial) .*        ( Ywm .* Ypp.^.1 .* (Ygmt+Yecb).^2 .* Ypve).^2  )   ));
      
      Ymfc = Ymf2>2.5/3 & smooth3(Ymf2)<2.5/3 & Ypve; Ymfc(smooth3(Ymfc)<.3)=0; 
      Ymf2 = max(Ymf2,Ymfc*3.5/3); % PVE voxel that disappear for light changes
      
      cat_sanlm(Ymf2,1,3);
      clear Yecb
    end

   
    %% final update
    Ymf   = Ymf1 * (th_initial-1/3)*3 + (1 - (th_initial-1/3)*3) * Ymf2;
    clear Xmf1 Ymf2 Ypve Ygmt Ypp;
  end
 
  Ymfs = cat_vol_smooth3X(Ymf,1);
  Ytmp = cat_vol_morph(YMF,'dd',3,vx_vol) & Ymfs>2.1/3;
  Ymf(Ytmp) = max(min(Ym(Ytmp),0),Ymfs(Ytmp)); clear Ytmp Ymfs; 
  Ymf = Ymf * 3;

if use_fullpreprocess
  % removing fine WM structures in the hippocampus area to reduce topological and geometrical defects (added RD20190912)
  % use erode to reduce probability of cutting other gyri
  HCmask = cat_vol_morph( NS(Ya,opt.LAB.HC) , 'de', 1.5, vx_vol) & ~YMF; % RD202501: open only not filled regions
  Ymf( HCmask ) =  min(2,Ymf( HCmask )); clear HCmask; 
end

  % surface output and evaluation parameter 
  Psurf = struct(); 
  res   = struct('lh',struct(),'rh',struct()); 
  
  
  %% reduction of artifact, blood vessel, and meninges next to the cortex in SPM segmentations 
  %  (are often visible as very thin structures that were added to the WM 
  %  or removed from the brain)
  if ~opt.SPM && use_fullpreprocess
    Ydiv  = cat_vol_div(Ymf,vx_vol); 
    Ycsfd = cat_vbdist(single(Ymf<1.8),Ymf>1,vx_vol);
    Yctd  = cat_vbdist(single(Ymf<0.5),Ymf>0,vx_vol); 
    Ysroi = Ymf>2  &  Yctd<10  & Ycsfd>0 & Ycsfd<2.0 & ...
            cat_vol_morph(~NS(Ya,opt.LAB.HC) & ~NS(Ya,opt.LAB.HI) & ...
              ~NS(Ya,opt.LAB.PH) & ~NS(Ya,opt.LAB.VT),'erode',4); 
    if ~debug, clear Ycsfd Yctd; end 
    Ybv   = cat_vol_morph(Ymf+Ydiv./max(1,Ymf)>3.5,'d') & Ymf>2; 
    Ymf(Ybv) = 1.4; 
    Ymfs  = cat_vol_median3(Ymf,Ysroi | Ybv,Ymf>eps & ~Ybv,0.1); % median filter
    Ymf   = min(Ymf, mf * Ymfs  +  (1-mf) * Ymf);
    clear Ysroi Ydiv Ybv

  end
  if ~debug, clear Ysroi Ymfs Yctd Ybv Ymfs; end
  
  % cleanup and smoothing of the hippocampus amygdala to remove high
  % frequency structures that we cannot preprocess yet
  if use_fullpreprocess
    Ymf = hippocampus_amygdala_cleanup(Ymf,Ya,vx_vol,opt.close_parahipp,1); % last var = doit 
  end
  
  % Sharpening of thin structures in the cerebellum (gyri and sulci)
  if ~opt.SPM && opt.sharpenCB && any(~cellfun('isempty',strfind(opt.surf,'cb')))
    Ymf = sharpen_cerebellum(Ym,Ymf,Ytemplate,Ya,vx_vol,opt.verb);
  end
  
  % main loop for each surface structure 
  for si = 1:numel(opt.surf)
   
    % surface filenames
    Pm         = fullfile(pp0,mrifolder, sprintf('m%s',ff));    % raw
    Pdefects   = fullfile(pp0,mrifolder, sprintf('defects_%s.nii',ff));    % defect
    Pcentral   = fullfile(pp0_surffolder,sprintf('%s.central.%s.gii',opt.surf{si},ff));          % central
    Pcentralr  = fullfile(pp0_surffolder,sprintf('%s.central.resampled.%s.gii',opt.surf{si},ff));%#ok<NASGU> % central .. used in inactive path
    Player4    = fullfile(pp0_surffolder,sprintf('%s.layer4.%s.gii',opt.surf{si},ff));           % layer4
    PintL4     = fullfile(pp0_surffolder,sprintf('%s.intlayer4.%s',opt.surf{si},ff));            % layer4 intensity
    Ppial      = fullfile(pp0_surffolder,sprintf('%s.pial.%s.gii',opt.surf{si},ff));             % pial (GM/CSF)
    Pwhite     = fullfile(pp0_surffolder,sprintf('%s.white.%s.gii',opt.surf{si},ff));            % white (WM/GM)
    Pthick     = fullfile(pp0_surffolder,sprintf('%s.thickness.%s',opt.surf{si},ff));            % FS thickness / GM depth
    Ppbt       = fullfile(pp0_surffolder,sprintf('%s.pbt.%s',opt.surf{si},ff));                  % PBT thickness / GM depth
    Psphere    = fullfile(pp0_surffolder,sprintf('%s.sphere.%s.gii',opt.surf{si},ff));           % sphere
    Pspherereg = fullfile(pp0_surffolder,sprintf('%s.sphere.reg.%s.gii',opt.surf{si},ff));       % sphere.reg
    Pfsavg     = fullfile(opt.fsavgDir,  sprintf('%s.central.freesurfer.gii',opt.surf{si}));     % fsaverage central
    Pfsavgsph  = fullfile(opt.fsavgDir,  sprintf('%s.sphere.freesurfer.gii',opt.surf{si}));      % fsaverage sphere    
    
    % use surface of given (average) data as prior for longitudinal mode
    if isfield(opt,'useprior') && ~isempty(opt.useprior) 
      % RD20200729: delete later ... && exist(char(opt.useprior),'file') 
      % if it not exist than filecopy has to print the error
      [pp1,ff1] = spm_fileparts(opt.useprior);
      % correct '../' parts in directory for BIDS structure
      [stat, val] = fileattrib(fullfile(pp1,surffolder));
      if stat
        pp1_surffolder = val.Name;
      else
        pp1_surffolder = fullfile(pp1,surffolder);
      end
      
      % try to copy surface files from prior to individual surface data 
      useprior = 1;
      useprior = useprior & copyfile(fullfile(pp1_surffolder,sprintf('%s.central.%s.gii',opt.surf{si},ff1)),Pcentral,'f');
      useprior = useprior & copyfile(fullfile(pp1_surffolder,sprintf('%s.sphere.%s.gii',opt.surf{si},ff1)),Psphere,'f');
      useprior = useprior & copyfile(fullfile(pp1_surffolder,sprintf('%s.sphere.reg.%s.gii',opt.surf{si},ff1)),Pspherereg,'f');
      
      if ~useprior
        warn_str = sprintf('Surface files for %s not found. Move on with individual surface extraction.\n',pp1_surffolder);
        fprintf('\nWARNING: %s',warn_str);
        cat_io_addwarning('cat_surf_createCS3:noPiorSurface', warn_str);
      else
        fprintf('\nUse existing surface as prior and thus skip many processing steps:\n%s\n',pp1_surffolder);
      end      
    else
      useprior = 0;
    end
    
    % add the variables defined in "surffile" to the "Psurf" output variable
    surffile = {'Pcentral','Pthick','Ppbt','Psphere','Pspherereg','Pfsavg','Pfsavgsph','Pwhite','Ppial'};
    for sfi=1:numel(surffile)
      eval(sprintf('Psurf(si).%s = %s;',surffile{sfi},surffile{sfi})); 
    end
    
    
    %% reduce for object area
    Ynocerebrum = ~(NS(Ya,opt.LAB.CB) | NS(Ya,opt.LAB.BV) | NS(Ya,opt.LAB.ON) | NS(Ya,opt.LAB.MB) | NS(Ya,opt.LAB.BS)); 
    switch opt.surf{si}
      case {'lh'}
        Ymfs  = Ymf .* (Ya>0) .* Ynocerebrum .* (mod(Ya,2)==1); 
        Yside = mod(Ya,2)==1; 
      case {'rh'}  
        Ymfs  = Ymf .* (Ya>0) .* Ynocerebrum .* (mod(Ya,2)==0); 
        Yside = mod(Ya,2)==0;  
      case {'cb'}
        Ymfs  = Ymf .* (Ya>0) .* NS(Ya,opt.LAB.CB); 
        Yside = true(size(Ymfs)); % new full cerebellum reconstruction
    end 
    clear Ynocerebrum
    % mark cutting regions to avoid cortical modelling (i.e. to avoid/reduce thickness estimates)
    Ycutregion = min(1,max(0,Ymf - 2 + NS(Ya,opt.LAB.VT)) * 1 .* ...
      smooth3( NS(Ya,opt.LAB.ON) | NS(Ya,opt.LAB.MB) | NS(Ya,opt.LAB.VT) | NS(Ya,opt.LAB.TH) | ~Yside )); 

    % print something
    if si==1, fprintf('\n'); end
    fprintf('%s:\n',opt.surf{si});
    
    % check for cerebellar processing  
    iscerebellum = strcmp(opt.surf{si},'cb');
    if ~iscerebellum 
      % RD202107:  Use atlas and Shooting template information to close the
      %            hippocampal gyrus but open the hippocampal region.
      mask_parahipp = opt.close_parahipp & NS(Ya,opt.LAB.PH) | NS(Ya,opt.LAB.HC) | NS(Ya,opt.LAB.VT); 
      if isempty(Ytemplate)
        VT = NS(Ya,opt.LAB.PH) | NS(Ya,opt.LAB.HC); % this does not work but also should not create problems
      else
        VT = Ytemplate .* NS(Ya,opt.LAB.PH) | NS(Ya,opt.LAB.HC);
      end
      HC = NS(Ya,opt.LAB.HC); 
    end
   
    % removing background (smoothing to remove artifacts)
    switch opt.surf{si}
      case {'lh','rh'}
        [Ymfs,Ysidei,mask_parahipp,VT,HC,Ycutregion,BB] = ...
          cat_vol_resize({Ymfs,Yside,mask_parahipp,VT,HC,Ycutregion},'reduceBrain',vx_vol,4,smooth3(Ymfs)>1.5); 
      case {'cb'}
        [Ymfs,Ysidei,Ycutregion,BB] = ...
          cat_vol_resize({Ymfs,Yside,Ycutregion},'reduceBrain',vx_vol,4,smooth3(Ymfs)>1.5); 
    end
    
    % interpolation 
    imethod         = 'cubic'; % cubic should be better in general - however, linear is better for small thickness (version?)
    [Ymfs,resI]     = cat_vol_resize(max(1,Ymfs),'interp',V,opt.interpV,imethod);                  % interpolate volume
    Ysidei          = cat_vol_resize(Ysidei,'interp',V,opt.interpV,imethod)>0.5;                   % interpolate volume (small dilatation)
    Ycutregion      = cat_vol_resize(Ycutregion,'interp',V,opt.interpV,imethod);
    Ycutregiond     = cat_vol_smooth3X( cat_vol_morph(Ycutregion,'dd',1.5,opt.interpV), 4);        % a smooth version to mix maps
    
    if ~iscerebellum
      % get dilated mask of gyrus parahippocampalis and hippocampus of both sides
      %mask_parahipp = smooth3( cat_vol_morph(mask_parahipp,'dd',2/opt.interpV) );
      VT = cat_vol_resize(VT,'interp',V,opt.interpV); 
      HC = cat_vol_resize(HC,'interp',V,opt.interpV); 
      mask_parahipp = cat_vol_resize(mask_parahipp,'interp',V,opt.interpV)>0.5;          % interpolate volume
    end 
    
    % PVE with background will lead to a slight underestimation?
    Ymfs = min(3,max(1,Ymfs));
    
    % Write PP
    Vmfs = resI.hdrN;
    Vmfs.fname = fullfile(pp0,mrifolder, sprintf('%s_seg-%s.nii',ff,opt.surf{si}));
    Vmfs.pinfo = V0.pinfo;
    Vmfs = rmfield(Vmfs,'dat');
    Vmfs = rmfield(Vmfs,'private');
    matiBB            = spm_imatrix(V.mat   * [eye(4,3) [ (BB.BB([1,3,5])' - 1) ; 1]]); 
    Vmfs.mat(1:3,4) = matiBB(1:3); 
    spm_write_vol(Vmfs, Ymfs);
    gmt_name = fullfile(pp0,mrifolder, sprintf('%s_thickness-%s.nii',ff,opt.surf{si}));
    ppm_name = fullfile(pp0,mrifolder, sprintf('%s_ppm-%s.nii',ff,opt.surf{si}));

    stime = cat_io_cmd(sprintf('  Thickness estimation (%0.2f mm%s)',opt.interpV,native2unicode(179, 'latin1')),'g5'); stimet = stime;
    if strcmp(opt.pbtmethod,'pbtsimple') 
      if use_cat_vol_pbtsimple
        [Yth1i,Yppi,Ymfsc] = cat_vol_pbtsimple(Ymfs, opt.interpV, struct('classic', opt.SRP<1,'myelinCorrection',myelinCorrection));
        Ymfs = Ymfs .* Ycutregion + (1-Ycutregion) .* Ymfsc; clear Ymfsc
        Vppm = Vmfs;
        Vppm.fname = ppm_name;
        spm_write_vol(Vppm, Yppi);
        Vppm = spm_vol(ppm_name);        
      else
        cmd = sprintf('CAT_VolThicknessPbt -median-filter 2 -sharpen 0.02 -fwhm 3 -downsample 0 "%s" "%s" "%s"', Vmfs.fname, gmt_name, ppm_name);
        cat_system(cmd,opt.verb-3);
        Vgmt = spm_vol(gmt_name); Yth1i = spm_read_vols(Vgmt);      
        Vppm = spm_vol(ppm_name); Yppi = spm_read_vols(Vppm);
      end
    else 
      [Yth1i,Yppi,Ymfs] = cat_vol_pbt(Ymfs,struct('method',opt.pbtmethod,'resV',opt.interpV,'vmat',...
        V.mat(1:3,:)*[0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1],'pbtlas',opt.pbtlas)); % avoid underestimated thickness in gyri
    end   
    % RD20250330: The aim is to smoothing cutting regions to avoid defects there. 
    Yppi = max(0,min(1,(cat_vol_smooth3X(Yppi,4)-.5)*4 + .5)) .* Ycutregiond  + (1-Ycutregiond) .* Yppi; 
  

    Yth = cat_vol_resize(Yth1i,'deinterp',resI);                       % back to original resolution
    Yth = cat_vol_resize(Yth+1,'dereduceBrain',BB)-1;                  % adding background
    
    % use interpolated ppm image instead of 0.5mm version
    if ~use_high_resolution
      %% RD202503: Smoothing before resolution reduction clearly improves the results.
      if create_surf == 2 && deformsurf
        % RD202503: Adaptation for thickness to be closer to CSF/WM in thick/thin regions - need further tests!
        Ypps = max(Yppi>=1., Yppi .* min(2,max(.5, Yth1i / (1.25 / opt.interpV) ) )); 
      else
        Ypps = Yppi; 
      end
      if create_surf == 2 && deformsurf  ||  1  % special filling
        for xi = 1:2
          [~,D] = cat_vol_downcut( single(cat_vol_morph( Ypps<.5,'ldo',4)), 1 - Ypps, .02); 
          D = cat_vol_median3(D,D>1000);
          Yfill = cat_vol_morph( Ymfs>2.2 & Ypps<.7 & smooth3(D)>1000 , 'd' );  clear D; 
          Ypps  = max(Ypps,Yfill); clear Yfill; 
        end
      end
      Ypps  = smooth3(Ypps);
      Ypp   = Yppi .* Ycutregiond  +  (1-Ycutregiond) .* Ypps;   
      Yfill = Ypp < th_initial  &  ~cat_vol_morph(Ypp < th_initial,'ldo',1); Ypp( Yfill ) = max(Ypp(Yfill),th_initial); 
      Ypp  = cat_vol_resize( Ypp ,'deinterp',resI);                       % back to original resolution
      Ypp  = cat_vol_resize(Ypp,'dereduceBrain',BB);
      Vpp  = cat_io_writenii(V,Ypp,'',sprintf('%s.pp',opt.surf{si}) ,...
        'initial percentage position map - V2 - pial=1/3, central=1/2, white=2/3, 0 and 1 to stablize white/pial values','int8',[0,1/100],[1 0 0]); 
      Vppm = Vpp(1);
    end
    clear M v x3 Ycutregiond  

    if opt.vol
      S = struct(); Psurf = '';
      if opt.verb<2, fprintf('%5.0fs',etime(clock,stime)); end
      continue; 
    end
    
    
    %% surface creation 
    %  --------------------------------------------------------------------
    %  Surface create should be at 0.5 mm to support a useful description
    %  in even narrow sulci, i.e., for 1 mm thickness the sulci would be 
    %  2 mm wide were a 1 mm percentage position map is strongly limited.
    %  The surface is reconstructed by marching cubes on a binary version 
    %  of the PBT radial position map Yppi to use a simple voxel-based
    %  topology correction to avoid small defects and incorrect
    %  triangulation (e.g., matlab isosurface function). Next, the
    %  complexity of the surface is reduced by the spm_mesh_reduce function
    %  (more accurate and stable than the matlab reduce function?) to
    %  improve the performance of the following main topology correction. 
    %  However, the reduction can partially lead to very large faces that
    %  need an additional refine. 
    %  Voxelbased topology optimization is important for the hippo. gyrus. 
    %  --------------------------------------------------------------------
    nosurfopt = iscerebellum;  %#ok<NASGU>
    if ~useprior 
      fprintf('%s%4.0fs\n',repmat(' ',1,67*opt.verb>1),etime(clock,stimet)); % add space in case of details
      if ~use_high_resolution
        stime = cat_io_cmd(sprintf('  Create initial surface (%0.2fx%0.2fx%0.2fmm)',vx_vol),'g5','',opt.verb); %if opt.verb>2, fprintf('\n'); end
      else
        stime = cat_io_cmd(sprintf('  Create initial surface (%0.2fmm)',opt.interpV),'g5','',opt.verb); %if opt.verb>2, fprintf('\n'); end
      end
    else
      fprintf('\n');
      stime = cat_io_cmd('  Load and refine subject average surface','g5','',opt.verb); %if opt.verb>2, fprintf('\n'); end
    end
    % surface coordinate transformation matrix
    matI              = spm_imatrix(V.mat); 
    matI(7:9)         = sign( matI(7:9))   .* repmat( opt.interpV , 1 , 3); 
    matiBB            = spm_imatrix(V.mat   * [eye(4,3) [ (BB.BB([1,3,5])' - 1) ; 1]]); 
    matIBB            = matiBB; 
    matIBB(7:9)       = sign( matiBB(7:9)) .* repmat( opt.interpV , 1 , 3); 
    Smat.matlabi_mm   = V.mat * [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];                % CAT internal space
    Smat.matlabI_mm   = spm_matrix(matI) * [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];     % PBT interpolated space
    Smat.matlabIBB_mm = spm_matrix(matIBB) * [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];   % PBT interpolated
    Smat.matlabiBB_mm = spm_matrix(matiBB) * [0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];   % PBT interpolated
    

    
    %% 
    if ~useprior 

      % Main initial surface creation using cat_vol_genus0
      % cat_vol_genus0 uses a "simple" marching cube without use of isovalues
      % that is used in the MATLAB isosurface function. Our test showed that
      % the surface deformation allows the same or better accuracy and also
      % that the meshes of cat_vol_genus0 are more regular and also allow 
      % voxel-based topology optimization.  
      % We can set dist to 0.9 for the cortex because using adaptive
      % thresholds this was found to be the optimal dist threshold

      %cmd = sprintf('CAT_MarchingCubesGenus0 -fwhm "3" -thresh "%g" "%s" "%s"',th_initial,Vppm.fname,Pcentral);
      % -pre-fwhm = -1: masked smoothing with 1mm
      % -post-fwhm = 2: psot smoothing with 2mm (use lower threshold of 0.49)
      % -scl-opening ... 0.9 too less for old subjects with thin WM and WMHs >> 1.5 
      %cmd = sprintf('CAT_VolMarchingCubes -local-smoothing 10 -scl-opening "0.9" -median-filter "2" -pre-fwhm "-2" -post-fwhm "1.5" -thresh "%g" "%s" "%s"',th_initial,Vppm.fname,Pcentral);
      if ~deformsurf, th_initial = .5; end
      if exist(Pcentral,'file'), delete(Pcentral); end % have to delete it to get useful error messages in case of reprocessing/testing
      if create_surf == 1 % old version
        cmd = sprintf('CAT_VolMarchingCubesOld -local-smoothing 10 -scl-opening "%d" -median-filter "2" -pre-fwhm "%0.1f" -post-fwhm "%0.1f" -thresh "%g" "%s" "%s" "%s"', ...
          createsurface_scl, createsurface_prefwhm, createsurface_postfwhm, th_initial, Vppm.fname, Pcentral, Pdefects);
      else 
        cmd = sprintf('CAT_VolMarchingCubes "%s" "%s" -thresh "%0.4f" ', ... -pre-fwhm "%0.4f" -post-fwhm "%0.4f"
          Vppm.fname, Pcentral, .6); % , createsurface_prefwhm, createsurface_postfwhm 
      end
      cat_system(cmd,opt.verb-3);
     
      if ~use_high_resolution  
      % Update Ypp without offsets for surface reconstruction
        Ypp  = cat_vol_resize(Yppi,'deinterp',resI);                      
        Ypp  = cat_vol_resize(Ypp,'dereduceBrain',BB);
        Vpp  = cat_io_writenii(V,Ypp,'',sprintf('%s.pp',opt.surf{si}) ,...
          'percentage position map - V2 - pial=1/3, central=1/2, white=2/3, 0 and 1 to stablize white/pial values','int8',[0,1/100],[1 0 0]); 
      end
    
      
      %%
      if deformsurf==1
        % old surface deformation function 
        stepsize = [0.1 0.01 0.001]; 
        for ss = 1:numel(stepsize)
          cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none  0  1  -1  .1 ' ...           
                      'avg  -0.1  0.1 .2  .1  5  0 "0.5"  "0.5"  n 0  0  0 %d  %g  0.0 0'], ...    
                      Vppm.fname,Pcentral,Pcentral,20,stepsize(ss) );
          cat_system(cmd,opt.verb-3);
        end
      elseif deformsurf==2
        % new surface reconstruction function
        % w1:    0.05 (schärfer) - 0.15 (smoother)
        % w2:    gradientstärke:   0.1 geringer effect
        % w3:    isovalue:         0.8 - 1.2 "schrittweite"                 >> 4 
        % sigma: filterbreite      0.2* - 0.4 (und w3 zusammen reduzieren)  >> 0.05
        deformsurfversion = 3; 
        if deformsurfversion == 1
        % original single call
          cmd = sprintf('CAT_SurfDeform "%s" "%s" "%s" ',Vppm.fname,Pcentral,Pcentral);
          cat_system(cmd,opt.verb-3);
        elseif deformsurfversion == 2 
        % modified single call (better)
          cmd = sprintf('CAT_SurfDeform -iter 150 -w1 0.01 -w3 4.0 -sigma 0.05 "%s" "%s" "%s" ',Vppm.fname,Pcentral,Pcentral);
          cat_system(cmd,opt.verb-3);
        else
        % multi call with some adaptions (best) 
          iters    = [50  40  30   20   10]; 
          stepsize = [0.2 0.1 0.05 0.02 0.01]; 
          for ss = 1:numel(stepsize)
            cmd = sprintf('CAT_SurfDeform -iter %d -w1 "%0.8f" -w3 "%0.8f" -verbose -sigma 0.05 "%s" "%s" "%s" ', ...    
                        iters(ss), stepsize(ss)/10, stepsize(ss)*20, Vppm.fname, Pcentral, Pcentral); %
            [~,txt] = cat_system(cmd,opt.verb - 2);
          end
        end
      end

      CS = loadSurf(Pcentral);
      if setcut2zero
        facevertexcdatanocut = 1 - cat_surf_fun('isocolors',Ycutregion,CS.vertices,Smat.matlabIBB_mm); 
        facevertexcdatanocut = max(0.5,min(1,(cat_surf_fun('smoothcdata',CS,facevertexcdatanocut,4) - .5) * 2)); % * median(vx_vol)
      else
        facevertexcdatanocut = 1; 
      end
      EC0 = size(CS.vertices,1) + size(CS.faces,1) - size(spm_mesh_edges(CS),1);
      if EC0 ~= 2
        warning('cat_surf_createCS3:CAT_MarchingCubesGenus0', ...
           'Extracted surface might have small topology issues (Euler count = %d).\n',EC0); 
      end      
                        
      % evaluate and save results
      if isempty(stime), stime = clock; end
      fprintf('%5.0fs',etime(clock,stime)); stime = []; if 1, fprintf('\n'); end %debug
    else
    % for use of average prior in long. pipeline we have to deform the average mesh to current pp distance map        
      CS = loadSurf(Pcentral);
      correct_mesh = 1;
      cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none  0  1  -1  .1 ' ...   
                     'avg  %0.3f  %0.3f .2  .1  5  0 "0.5"  "0.5"  n 0  0  0 %d %g  0.0 0 %d'], ...          
                      Vppm.fname,Pcentral,Pcentral,-0.1,0.1,100,0.01,correct_mesh); 
      cat_system(cmd,opt.verb-3);
    end

    % get thickness data
    facevertexcdata = cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm) .* facevertexcdatanocut; 
    cat_io_FreeSurfer('write_surf_data',Ppbt,facevertexcdata);

    if opt.SRP>0
      %% call collision correction
      %  RD202108: Use further iterations if self-intersections are still very high.  
      %            (test data was an high resolution ex-vivo chimp PD image that had still strong SIs after first correction) 
      stime  = cat_io_cmd('  Reduction of surface collisions with optimization:','g5','',opt.verb,stime); 
      verblc = 0; %debug;  
      for xi = 1:2 %opt.SRP % interation 1 with optimization, 2 without .. 
        [CS,facevertexcdata,SIs] = cat_surf_fun('collisionCorrectionPBT',CS,facevertexcdata,Ymfs,Yppi,...
          struct('optimize',xi==1,'verb',verblc,'mat',Smat.matlabIBB_mm,'vx_vol',vx_vol,'CS4',1)); 
        if verblc, fprintf('\b\b'); end
        if xi == 1 && opt.SRP>1 % this function is quite slow if multiple jobs are running!
          [CS,facevertexcdata,SIs] = cat_surf_fun('collisionCorrectionRY' ,CS,facevertexcdata,Yppi,...
            struct('Pcs',Pcentral,'verb',verblc,'mat',Smat.matlabIBB_mm,'accuracy',1/2^2));
          if verblc, fprintf('\b\b'); end
        end
      end
      facevertexcdata = facevertexcdata .* facevertexcdatanocut; 
      saveSurf(CS,Pcentral); cat_io_FreeSurfer('write_surf_data',Ppbt,facevertexcdata);
    end

    % skip that part if a prior image is defined
    if ~useprior && ~skip_registration
      %% spherical surface mapping of the final corrected surface
      %  with minimal areal smoothing
      stime = cat_io_cmd('  Spherical mapping with areal smoothing','g5','',opt.verb,stime); 
      cmd = sprintf('CAT_Surf2Sphere "%s" "%s" %d',Pcentral,Psphere,6);
      cat_system(cmd,opt.verb-3);
      
      % spherical registration to fsaverage template
      stime = cat_io_cmd('  Spherical registration','g5','',opt.verb,stime);
      cmd = sprintf('CAT_SurfWarp -steps 2 -avg -i "%s" -is "%s" -t "%s" -ts "%s" -ws "%s"', ...
        Pcentral,Psphere,Pfsavg,Pfsavgsph,Pspherereg);
      cat_system(cmd,opt.verb-3);
    end  

    % create white and central surfaces
    if create_white_pial
      stime = cat_io_cmd('  Create pial and white surface','g5','',opt.verb,stime); 
      cmd = sprintf('CAT_Central2Pial "%s" "%s" "%s" 0.5',Pcentral,Ppbt,Ppial);
      cat_system(cmd,opt.verb-3);
      cmd = sprintf('CAT_Central2Pial "%s" "%s" "%s" -0.5',Pcentral,Ppbt,Pwhite);
      cat_system(cmd,opt.verb-3);
    end

    % estimate Freesurfer thickness measure Tfs using mean(Tnear1,Tnear2)
    if debug || cat_get_defaults('extopts.expertgui')>1, fprintf('\n'); end
    if opt.thick_measure == 1
      % not ready yet
      if 0
        cmd = sprintf('CAT_SurfDistance -mean "%s" "%s" "%s"',Pwhite,Ppial,Pthick);
        cat_system(cmd,opt.verb-3);
      else % use central surface and thickness
        cmd = sprintf('CAT_SurfDistance -mean -thickness "%s" "%s" "%s"',Ppbt,Pcentral,Pthick);
        cat_system(cmd,opt.verb-3);
      end
      
      % apply upper thickness limit
      facevertexcdata = cat_io_FreeSurfer('read_surf_data',Pthick) .* facevertexcdatanocut;  
      facevertexcdata(facevertexcdata > opt.thick_limit) = opt.thick_limit;
      cat_io_FreeSurfer('write_surf_data',Pthick,facevertexcdata);  
      
      % final surface evaluation 
      res.(opt.surf{si}).createCS_final = cat_surf_fun('evalCS', ...
        loadSurf(Pcentral), cat_io_FreeSurfer('read_surf_data',Ppbt), facevertexcdata, ...
        Ymfs,Yppi,Pcentral,Smat.matlabIBB_mm,debug + (cat_get_defaults('extopts.expertgui')>1),cat_get_defaults('extopts.expertgui')>1);
    else % otherwise simply copy ?h.pbt.* to ?h.thickness.*
      copyfile(Ppbt,Pthick,'f'); 
      res.(opt.surf{si}).createCS_final = cat_surf_fun('evalCS', ...
        loadSurf(Pcentral), cat_io_FreeSurfer('read_surf_data',Ppbt), [], ...
        Ymfs, Yppi, Pcentral, Smat.matlabIBB_mm, debug, cat_get_defaults('extopts.expertgui')>1);
    end
    fprintf('\n');
    
    % average final values
    FNres = fieldnames( res.(opt.surf{si}).createCS_final );
    for fnr = 1:numel(FNres)
      if ~isfield(res,'final') || ~isfield(res.final,FNres{fnr})
        res.final.(FNres{fnr}) = res.(opt.surf{si}).createCS_final.(FNres{fnr}) / numel(opt.surf);
      else
        res.final.(FNres{fnr}) = res.final.(FNres{fnr}) + res.(opt.surf{si}).createCS_final.(FNres{fnr}) / numel(opt.surf);
      end
    end
    if isfield(res.(opt.surf{si}),'createCS_resampled') 
      FNres = fieldnames( res.(opt.surf{si}).createCS_resampled );
      for fnr = 1:numel(FNres)
        if isfield(res.(opt.surf{si}),'createCS_resampled') 
          if ~isfield(res,'createCS_resampled') || ~isfield(res.createCS_resampled,FNres{fnr}) 
            res.resampled.(FNres{fnr}) = res.(opt.surf{si}).createCS_resampled.(FNres{fnr}) / numel(opt.surf);
          else
            res.resampled.(FNres{fnr}) = res.resampled.(FNres{fnr}) + res.(opt.surf{si}).createCS_resampled.(FNres{fnr}) / numel(opt.surf);
          end
        end
      end
    end
    
    
    % create output structure
    S.(opt.surf{si}) = struct('faces',CS.faces,'vertices',CS.vertices,'th1',facevertexcdata); clear facevertexcdata; 
    clear Yth1i

    % we have to delete the original faces, because they have a different 
    % number of vertices after CAT_FixTopology!
    %if exist(Vppm.fname ,'file'), delete(Vppm.fname); end
    if debug && ~use_high_resolution && exist(Vpp.fname ,'file') && ~opt.outputpp.native, delete(Vpp.fname); end
    
    clear CS

    % processing time per side for manual tests
    if si == numel(opt.surf) && si == 1
      cat_io_cmd('  ','g5','',opt.verb);
      fprintf('%5ds\n',round(etime(clock,cstime)));
    end
  end  
  
  %% calculate surface quality parameters for all surfaces
  mnth = []; sdth = []; mnRMSE_Ypp = []; mnRMSE_Ym = []; sdRMSE_Ym = []; sdRMSE_Ypp = []; 
  SIw = []; SIp = []; SIwa = []; SIpa = []; 
  for si=1:numel(opt.surf)
    if any(strcmp(opt.surf{si},{'lh','rh'}))
      if isfield(res.(opt.surf{si}).createCS_final,'fsthickness_mn_sd_md_mx') && ... 
        ~isnan( res.(opt.surf{si}).createCS_final.fsthickness_mn_sd_md_mx(1) )
        mnth      = [ mnth  res.(opt.surf{si}).createCS_final.fsthickness_mn_sd_md_mx(1) ]; 
        sdth      = [ sdth  res.(opt.surf{si}).createCS_final.fsthickness_mn_sd_md_mx(2) ]; 
      else
        mnth      = [ mnth  res.(opt.surf{si}).createCS_final.thickness_mn_sd_md_mx(1) ];
        sdth      = [ sdth  res.(opt.surf{si}).createCS_final.thickness_mn_sd_md_mx(2) ]; 
      end
      mnRMSE_Ym   = [ mnRMSE_Ym   mean([...
        res.(opt.surf{si}).createCS_final.RMSE_Ym_layer4 ...
        res.(opt.surf{si}).createCS_final.RMSE_Ym_white ...
        res.(opt.surf{si}).createCS_final.RMSE_Ym_pial ]) ];
      sdRMSE_Ym   = [ sdRMSE_Ym  std([...
        res.(opt.surf{si}).createCS_final.RMSE_Ym_layer4 ...
        res.(opt.surf{si}).createCS_final.RMSE_Ym_white ...
        res.(opt.surf{si}).createCS_final.RMSE_Ym_pial ]) ];
      mnRMSE_Ypp  = [ mnRMSE_Ypp  mean([...
        res.(opt.surf{si}).createCS_final.RMSE_Ypp_central ...
        res.(opt.surf{si}).createCS_final.RMSE_Ypp_white ...
        res.(opt.surf{si}).createCS_final.RMSE_Ypp_pial ]) ];
      sdRMSE_Ypp  = [ sdRMSE_Ypp  std([...
        res.(opt.surf{si}).createCS_final.RMSE_Ypp_central ...
        res.(opt.surf{si}).createCS_final.RMSE_Ypp_white ...
        res.(opt.surf{si}).createCS_final.RMSE_Ypp_pial ]) ];
      if isfield(res.(opt.surf{si}).createCS_final,'white_self_interections')
        SIw     = [ SIw  res.(opt.surf{si}).createCS_final.white_self_interections ]; 
        SIp     = [ SIp  res.(opt.surf{si}).createCS_final.pial_self_interections  ]; 
        SIwa    = [ SIwa res.(opt.surf{si}).createCS_final.white_self_interection_area ]; 
        SIpa    = [ SIpa res.(opt.surf{si}).createCS_final.pial_self_interection_area  ]; 
      end
    end
  end
  
  % final res structure
  res.Smat        = Smat; 
  res.EC          = NaN; 
  res.defect_size = NaN;
  res.defect_area = NaN;
  res.defects     = NaN;
  res.RMSE_Ym     = mean(mnRMSE_Ym);
  res.RMSE_Ypp    = mean(mnRMSE_Ypp);
    
  
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
      
    if cat_get_defaults('extopts.expertgui')
    % color output currently only for expert ...
      if isfield(res.(opt.surf{si}).createCS_final,'fsthickness_mn_sd_md_mx')
        fprintf('  Average thickness (FS):                     ');
      else
        fprintf('  Average thickness (PBT):                    ');
      end
      cat_io_cprintf( color( rate( abs( mean(mnth) - 2.5 ) , 0 , 2.0 )) , sprintf('%0.4f'  , mean(mnth) ) );  fprintf(' %s ',native2unicode(177, 'latin1'));
      cat_io_cprintf( color( rate( abs( mean(sdth) - 0.5 ) , 0 , 1.0 )) , sprintf('%0.4f mm\n', mean(sdth) ) );

      fprintf('  Surface intensity / position RMSE:          ');
      cat_io_cprintf( color( rate( mean(mnRMSE_Ym)  , 0.05 , 0.3 ) ) , sprintf('%0.4f / ', mean(mnRMSE_Ym) ) );
      cat_io_cprintf( color( rate( mean(mnRMSE_Ypp) , 0.05 , 0.3 ) ) , sprintf('%0.4f\n', mean(mnRMSE_Ypp) ) );
        
      if isfield(res.(opt.surf{si}).createCS_final,'white_self_interections')
        fprintf('  Pial/white self-intersections:              ');
        cat_io_cprintf( color( rate(  mean([SIw,SIp]) , 0 , 20 ) ) , sprintf('%0.2f%%%% (%0.2f mm%s)\n'  , mean([SIw,SIp]) , mean([SIwa,SIpa]) , char(178) ) );
      end
    else
      fprintf('  Average thickness:                          %0.4f %s %0.4f mm\n' , mean(mnth), native2unicode(177, 'latin1'), mean(sdth));
    end
    
    for si=1:numel(Psurf)
      fprintf('  Display thickness:          %s\n',spm_file(Psurf(si).Pthick,'link','cat_surf_display(''%s'')'));
    end
    
    %% surfaces in spm_orthview
    if exist(Pm,'file'), Po = Pm; else, Po = V0.fname; end
    
    Porthfiles = '{'; Porthcolor = '{'; Porthnames = '{';
    for si=1:numel(Psurf)
      Porthfiles = [ Porthfiles , sprintf('''%s'',''%s'',',Psurf(si).Ppial, Psurf(si).Pwhite )]; 
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
function Ymf = hippocampus_amygdala_cleanup(Ymf,Ya,vx_vol,close_parahipp,doit)
%% Amygdala hippocampus smoothing. 
%  We use a median filter to remove the nice details of the hippocampus
%  that will cause topology errors and self-intersections. 
%  Currently, I have no CAT ROI for Amygdala - but it would be more 
%  robust to filter (simple smoothing) this region strongly because 
%  the "random" details especially in good data introduce more variance.
%  (RD 201912)
%  RD202107: Added closing of the parahippocampal gyrus. 
%
%    Ymf = hippocampus_amygdala_cleanup(Ymf,Ya[,doit])
% 
%    Ymf  .. intensity normalized (WM=3,GM=2,CSF=1) filled and 
%            skull-stripped image
%    Ya   .. CAT atlas map
%    doit .. do it (default = 1) 
   
  if ~exist('doit','var'), doit = 1; end
  
  if doit
    % remove side definition 
    NS   = @(Ys,s) Ys==s | Ys==s+1; 
    LAB  = cat_get_defaults('extopts.LAB');  
   
     
    %% RD202107: Close CSF region between hippocampus and amygdala
    %  --------------------------------------------------------------------
    %  This could be part of cat_vol_partvol to improve the 
    %  detection of the lower arms of the ventricles.  The region has 
    %  slightly increased variance, but much less than the topology problems
    %  parahippocampal region.  Nevertheless, the cuts in some brains that 
    %  trigger the variance are much deeper than we would generally expect
    %  (e.g. FSaverage) and therefore closing seems appropriate.
    %  --------------------------------------------------------------------
    if 1  
      Ysv  = NS(Ya,LAB.PH) & Ymf<1.9 & Ymf>0.5; 
      Ysv(smooth3(Ysv)<0.5) = 0; 
      Ysv  = cat_vol_morph( Ysv , 'do' , 1.4 , vx_vol);
      Ysvd = cat_vol_morph( Ysv , 'dd' , 5   , vx_vol); 
      Ysvc = cat_vol_morph( (Ysvd & Ymf>2) | Ysv , 'lc', 2 , vx_vol);
      Ysvc = smooth3(Ysvc);
      Ymf  = max( Ymf , Ysvc * 3); 
      clear Ysv Ysvc Ysvd; 
    end


    
    Yphg = false(size(Ymf)); 
    
    
    
    %% strong cleanup by median filter within the hippocampus
    Ymsk = cat_vol_morph( NS(Ya,LAB.PH) | NS(Ya,LAB.ON) | NS(Ya,LAB.BS) , 'dd', 2 );
    Ymsk = Ymf>0 & cat_vol_morph( NS(Ya,LAB.HC) , 'dd' , 3 , vx_vol ) & ~Ymsk & ~Yphg; 
    Ymf  = cat_vol_median3( Ymf , Ymsk ); 
    Ymf  = cat_vol_median3( Ymf , Ymsk ); 

    % further cleanup by smoothing
    Ymsk = NS(Ya,LAB.PH) | NS(Ya,LAB.ON) | NS(Ya,LAB.BS);
    Ymsk = Ymf>0 & NS(Ya,LAB.HC) & ~Ymsk & ~Yphg; 
    Ymsk = smooth3(Ymsk); 
    Ymf  = min(Ymf,3-Ymsk); 
  end
end


%==========================================================================
function Ymf = sharpen_cerebellum(Ym,Ymf,Ytemplate,Ya,vx_vol,verb,doit)
%% Sharpening of thin structures in the cerebellum (gyri and sulci)
%  Processing of the cerebellum needs a complete updated with a simplified
%  reconstruction rather than the maximum level of detail. However, the 
%  main branches are quite similare and their is a lot of variance in 
%  aging and between MR protocols. 
%  (RD 2015-201912)  

  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end
  
  if ~exist('doit','var'), doit = 1; end
  if ~exist('verb','var'), verb = 0; end
  
  if doit
    if verb>2, fprintf('\n'); stime = cat_io_cmd('  Sharpen cerebellum'); end
    %Ytemplate = Ytemplate * 3; 
    
    LAB  = cat_get_defaults('extopts.LAB');  
    NS   = @(Ys,s) Ys==s | Ys==s+1; 
    Ywmd = cat_vbdist(single(Ymf>2.5),Ymf>1,vx_vol);

    %% normalized Ym divergence
    Ydiv = cat_vol_div(max(2,Ymf)); %Ydivl  = cat_vol_div(Ymf,vx_vol); 
    Ydc  = cat_vol_localstat(abs(Ydiv),Ymf>0,1,3);
    Ydc  = cat_vol_localstat(Ydc,Ymf>0,1,1);
    Ydiv = Ydiv ./ Ydc; clear Ydc;
    Ydiv( abs(Ydiv) > 1.5 ) = 0;
    
    %% normalized Ytemplate divergence
    if ~isempty(Ytemplate)
      Ydivt = cat_vol_div(max(2,Ytemplate));
      Ydct  = cat_vol_localstat(abs(Ydivt),Ytemplate>0,1,3);
      Ydct  = cat_vol_localstat(Ydct,Ymf>0,1,1);
      Ydivt = Ydivt ./ Ydct; clear Ydct;
      Ydivt( abs(Ydivt) > 1.5 ) = 0;
    else
      Ytemplate = Ymf; 
      Ydivt     = Ydiv; 
    end
    
    %% bias-correction based
    % WM 
    Ycsfd = cat_vbdist(single(Ymf<1.8),Ymf>1,vx_vol);
    Ymsk = ((cat_vol_morph(NS(Ya,LAB.CB),'e',3) | Ymf) & ( (Ym-Ydiv).*(Ytemplate/3-Ydivt) )>2/3 ) |  ...
           (NS(Ya,LAB.PH) & ( Ymf>2.2 | (Ymf>2 & Ydiv<-0.01) ) ) | ...                  % hippocampal gyri
           (NS(Ya,LAB.CT) & ( Ymf>2.2 | (Ymf>2 & Ydiv<-0.01 & ...
              Ycsfd>cat_stat_nanmean(Ycsfd(Ycsfd(:)>0 & Ycsfd(:)<100)) )*1.0) );            % distant gyri and sulci in the cerebrum
    Yi   = cat_vol_localstat(Ymf,Ymsk,1,3);
    % GM
    Ymsk = (NS(Ya,LAB.CB) & ( Ymf>1.9 & Ymf<2.2 & Ycsfd>0 & Ydiv>-0.05) ) | ...         % sulci and gyri in the cerebellum 
           (NS(Ya,LAB.PH) & ( Ymf>1.3 & Ymf<2.2 & Ycsfd>0 ) ) | ...                     % hippocampal gyri
           (NS(Ya,LAB.CT) & ( Ymf>1.3 & Ymf<2.2 & Ycsfd>0 & ...
              Ywmd>cat_stat_nanmean(Ywmd(Ywmd(:)>0 & Ywmd(:)<100))*0.2 ) );                 % distant gyri and sulci in the cerebrum
    Yi   = Yi + cat_vol_localstat(Ymf,Yi==0 & Ymsk,1,1)/2*3;%& ( cat_vol_morph(Yi==0,'e') & Ymf>2.2)
    Yi   = cat_vol_localstat(Yi,Yi>0,1,3);
    Yi   = cat_vol_localstat(Yi,Yi>0,1,1); 
    if ~debug, clear Ywmd Ymsk Ycsfd; end
    % CSF - instable and not required
    Ywi = cat_vol_approx(Yi,'nn',1,1,struct('lfO',2)); clear Yi
    
    %% only cerebellum
    Ycb = cat_vol_smooth3X( NS(Ya,LAB.CB) , 2 ); 
    Ywi = 3*ones(size(Ywi),'single').*(1-Ycb) + Ycb.*Ywi; clear Ycb 
    Ymf = Ymf./Ywi * 3; 
    
    % denoising result
    Ycb = Ymf .* NS(Ya,LAB.CB); cat_sanlm(Ycb,3,1); Ymf(NS(Ya,LAB.CB)) = Ycb(NS(Ya,LAB.CB)); clear Ycb
    
    
    %% sharpening (RD 201912)
    Ycb = (NS(Ya,LAB.CB)>0.5) .* max(0,min(1,min(2,max(-1,(Ymf/3).^2 - 0.1*Ydiv) .* max(-1,Ytemplate/3 - 0.02*Ydivt - 0.03*Ydiv )*3 - 1)/3 + 2/3)); 
    if ~debug, clear Ydiv; end
    for i=1:3, Ycb = max(0,min(1, Ycb - smooth3(cat_vol_median3(Ycb,Ycb>0,Ycb>0) - Ycb) )); end
    Ycb = min(1,Ycb); 
    cat_sanlm(Ycb,3,1); 
    
    %% final mixing
    Ymsk = cat_vol_smooth3X(NS(Ya,LAB.CB) & Ycb<1.9/3,0.5); 
    Ycb  = Ycb.*(1-Ymsk) + Ymsk.*Ymf/3;
    Ycs  = NS(Ya,LAB.CB) .* cat_vol_smooth3X( NS(Ya,LAB.CB) , 4 ); 
    Ymf  = Ymf.*(1-Ycs) + Ycs.*Ycb*3; clear Ycs 
    
    if verb>2, fprintf('%5.0fs\n',etime(clock,stime)); end
  end
end

              

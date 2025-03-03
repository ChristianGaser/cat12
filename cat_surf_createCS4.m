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

  % Turn off gifti data warning in gifti/subsref (line 45)
  %   Warning: A value of class "int32" was indexed with no subscripts specified. 
  %            Currently the result of this operation is the indexed value itself, 
  %            but in a future release, it will be an error. 
  warning('off','MATLAB:subscripting:noSubscriptsSpecified');
  cstime = clock;

  use_cat_vol_pbtsimple = 1;
  use_fullpreprocess = 1;
  use_high_resolution = 1;
  skip_registration = 1;
  
  % remove mask_parahipp
  % check Yth

  % set debugging variable
  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end
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
  
  % isosurface threshold
  th_initial = 0.50; 
  
  % too slow currently
  create_white_pial = opt.SRP==1;
  
  % Another parameter to control runtime is the accuracy of the surface
  % deformation. As far as we primary adapt the mesh resolution above, it 
  % is useful to use sqrt values rather than linear values to avoid square
  % processing times for higher quality levels. Otherwise, we can simple 
  % avoid changes here
  
  % function to estimate the number of interactions of the surface deformation: d=distance in mm and a=accuracy 
  QMC    = cat_io_colormaps('marks+',17);
  color  = @(m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
  rate   = @(x,best,worst) min(6,max(1, max(0,x-best) ./ (worst-best) * 5 + 1));
  
  % some internal overview for developers
  if opt.verb>2 
    fprintf('\nSurface reconstruction:              %s\n',....
      sprintf('%s',char( cellfun(@(x) [x ' '],opt.surf,'UniformOutput',0) )')); 
    fprintf('  PBT resolution:                    %0.3f\n',opt.interpV);
  end  
  
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

  %% get both sides in the atlas map
  NS = @(Ys,s) Ys==s | Ys==s+1; 
  
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
  
  % filling
  Ymf  = max(Ym,min(1,YMF & ~NS(Ya,opt.LAB.HC) & ~( cat_vol_morph( NS(Ya,opt.LAB.HC),'dd',2,vx_vol)))); 
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
    Pcentral   = fullfile(pp0_surffolder,sprintf('%s.central.%s.gii',opt.surf{si},ff));          % central
    Pcentralr  = fullfile(pp0_surffolder,sprintf('%s.central.resampled.%s.gii',opt.surf{si},ff));%#ok<NASGU> % central .. used in inactive path
    Player4    = fullfile(pp0_surffolder,sprintf('%s.layer4.%s.gii',opt.surf{si},ff));           % layer4
    PintL4     = fullfile(pp0_surffolder,sprintf('%s.intlayer4.%s',opt.surf{si},ff));            % layer4 intensity
    Ppial      = fullfile(pp0_surffolder,sprintf('%s.pial.%s.gii',opt.surf{si},ff));             % pial (GM/CSF)
    Pwhite     = fullfile(pp0_surffolder,sprintf('%s.white.%s.gii',opt.surf{si},ff));            % white (WM/GM)
    Pthick     = fullfile(pp0_surffolder,sprintf('%s.thickness.%s',opt.surf{si},ff));            % FS thickness / GM depth
    Ppbt       = fullfile(pp0_surffolder,sprintf('%s.pbt.%s',opt.surf{si},ff));                  % PBT thickness / GM depth
    Pgwo       = fullfile(pp0_surffolder,sprintf('%s.depthWMo.%s',opt.surf{si},ff));             % gyrus width / GWM depth / gyral span
    Pgw        = fullfile(pp0_surffolder,sprintf('%s.depthGWM.%s',opt.surf{si},ff));             % gyrus width / GWM depth / gyral span
    Pgww       = fullfile(pp0_surffolder,sprintf('%s.depthWM.%s',opt.surf{si},ff));              % gyrus width of the WM / WM depth
    Pgwwg      = fullfile(pp0_surffolder,sprintf('%s.depthWMg.%s',opt.surf{si},ff));             % gyrus width of the WM / WM depth
    Psw        = fullfile(pp0_surffolder,sprintf('%s.depthCSF.%s',opt.surf{si},ff));             % sulcus width / CSF depth / sulcal span
    Psphere    = fullfile(pp0_surffolder,sprintf('%s.sphere.%s.gii',opt.surf{si},ff));           % sphere
    Pspherereg = fullfile(pp0_surffolder,sprintf('%s.sphere.reg.%s.gii',opt.surf{si},ff));       % sphere.reg
    Pfsavg     = fullfile(opt.fsavgDir,  sprintf('%s.central.freesurfer.gii',opt.surf{si}));     % fsaverage central
    Pfsavgsph  = fullfile(opt.fsavgDir,  sprintf('%s.sphere.freesurfer.gii',opt.surf{si}));      % fsaverage sphere    
    
    % use surface of given (average) data as prior for longitudinal mode
    if isfield(opt,'useprior') && ~isempty(opt.useprior) 
      % RD20200729: delete later ... && exist(char(opt.useprior),'file') 
      % if it not exist than filecopy has to print the error
      priorname = opt.useprior;
      [pp1,ff1] = spm_fileparts(priorname);
      % correct '../' parts in directory for BIDS structure
      [stat, val] = fileattrib(fullfile(pp1,surffolder));
      if stat
        pp1_surffolder = val.Name;
      else
        pp1_folder = fullfile(pp1,surffolder);
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
    surffile = {'Pcentral','Pthick','Ppbt','Pgw','Pgww','Psw',...
      'Psphere','Pspherereg','Pfsavg','Pfsavgsph','Pwhite','Ppial'};
    for sfi=1:numel(surffile)
      eval(sprintf('Psurf(si).%s = %s;',surffile{sfi},surffile{sfi})); 
    end
    
    %% reduce for object area
    switch opt.surf{si}
      case {'lh'},  Ymfs = Ymf .* (Ya>0) .* ~(NS(Ya,opt.LAB.CB) | NS(Ya,opt.LAB.BV) | NS(Ya,opt.LAB.ON) | NS(Ya,opt.LAB.MB) | NS(Ya,opt.LAB.BS)) .* (mod(Ya,2)==1); Yside = mod(Ya,2)==1; 
      case {'rh'},  Ymfs = Ymf .* (Ya>0) .* ~(NS(Ya,opt.LAB.CB) | NS(Ya,opt.LAB.BV) | NS(Ya,opt.LAB.ON) | NS(Ya,opt.LAB.MB) | NS(Ya,opt.LAB.BS)) .* (mod(Ya,2)==0); Yside = mod(Ya,2)==0;  
      case {'cb'},  Ymfs = Ymf .* (Ya>0) .*   NS(Ya,opt.LAB.CB); Yside = true(size(Ymfs)); % new full cerebellum reconstruction
    end 
   
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
      case {'lh','rh'},  [Ymfs,Ysidei,mask_parahipp,VT,HC,BB] = cat_vol_resize({Ymfs,Yside,mask_parahipp,VT,HC},'reduceBrain',vx_vol,4,smooth3(Ymfs)>1.5); 
      case {'cb'},       [Ymfs,Ysidei,BB] = cat_vol_resize({Ymfs,Yside},'reduceBrain',vx_vol,4,smooth3(Ymfs)>1.5); 
    end
    
    % interpolation 
    imethod         = 'cubic'; % cubic should be better in general - however, linear is better for small thickness (version?)
    [Ymfs,resI]     = cat_vol_resize(max(1,Ymfs),'interp',V,opt.interpV,imethod);                  % interpolate volume
    Ysidei          = cat_vol_resize(Ysidei,'interp',V,opt.interpV,imethod)>0.5;                   % interpolate volume (small dilatation)
    
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

    stime = cat_io_cmd(sprintf('  Thickness estimation (%0.2f mm%s)',opt.interpV,native2unicode(179, 'latin1'))); stimet = stime;
    if strcmp(opt.pbtmethod,'pbtsimple') 
      if use_cat_vol_pbtsimple
        [Yth1i,Yppi] = cat_vol_pbtsimple(Ymfs,opt.interpV);
        Vppm = Vmfs;
        Vppm.fname = ppm_name;
        spm_write_vol(Vppm, Yppi);
        Vppm = spm_vol(ppm_name);        
      else
        cmd = sprintf('CAT_VolThicknessPbt -median-filter 2 -sharpen 0.02 -fwhm 3 -downsample 0 "%s" "%s" "%s"',Vmfs.fname,gmt_name, ppm_name);
        cat_system(cmd,opt.verb-3);
        Vgmt = spm_vol(gmt_name); Yth1i = spm_read_vols(Vgmt);      
        Vppm = spm_vol(ppm_name); Yppi = spm_read_vols(Vppm);
      end
    else 
      [Yth1i,Yppi,Ymfs] = cat_vol_pbt(Ymfs,struct('method',opt.pbtmethod,'resV',opt.interpV,'vmat',...
        V.mat(1:3,:)*[0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1],'pbtlas',opt.pbtlas)); % avoid underestimated thickness in gyri
    end  

    Yth = cat_vol_resize(Yth1i,'deinterp',resI);                         % back to original resolution
    Yth = cat_vol_resize(Yth+1,'dereduceBrain',BB)-1;                  % adding background
    
    % use interpolated ppm image instead of 0.5mm version
    if ~use_high_resolution
      Ypp = cat_vol_resize(Yppi + 1,'deinterp',resI);                       % back to original resolution
      Ypp = cat_vol_resize(Yppt,'dereduceBrain',BB) - 1;
      Vpp = cat_io_writenii(V,Yppt,'',sprintf('%s.pp',opt.surf{si}) ,...
        'percentage position map - V2 - pial=1/3, central=1/2, white=2/3, 0 and 1 to stablize white/pial values','uint8',[0,1/255],[1 0 0]); %clear Yppt
      Vppm = Vpp(1);
    end
    
    clear M v x3; 

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
    if opt.verb==1, fprintf('%s %4.0fs\n',repmat(' ',1,66),etime(clock,stimet)); end
    nosurfopt = iscerebellum;  %#ok<NASGU>
    if ~useprior 
      stime = cat_io_cmd('  Create initial surface','g5','',opt.verb); %if opt.verb>2, fprintf('\n'); end
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
      cmd = sprintf('CAT_VolMarchingCubes -local-smoothing 10 -scl-opening "0.9" -median-filter "2" -pre-fwhm "-2" -post-fwhm "1.5" -thresh "%g" "%s" "%s"',th_initial,Vppm.fname,Pcentral);
      cat_system(cmd,opt.verb-3);

      % Collins-without: 2.5996 ± 0.6292 mm, 0.0798 / 0.1096, 10.29% (121.38 cm²) 24.19% (285.46 cm²)
      % Collins-with:    2.5713 ± 0.6525 mm, 0.0723 / 0.0934,  8.51% (98.93 cm²)  23.79% (276.42 cm²)
      cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none  0  1  -1  .1 ' ...           
                  'avg  -0.1  0.1 .2  .1  5  0 "0.5"  "0.5"  n 0  0  0 %d  %g  0.0 0'], ...    
                  Vppm.fname,Pcentral,Pcentral,50,0.001);
      cat_system(cmd,opt.verb-3);

      %
            
      CS = loadSurf(Pcentral);
      EC0 = size(CS.vertices,1) + size(CS.faces,1) - size(spm_mesh_edges(CS),1);
      
      if EC0 ~= 2
        warning('cat_surf_createCS3:CAT_MarchingCubesGenus0', ...
           'Extracted surface might have small topology issues (Euler count = %d).\n',EC0); 
      end      
                  
      % evaluate and save results
      if isempty(stime), stime = clock; end
      fprintf('%5.0fs',etime(clock,stime)); stime = []; if 1, fprintf('\n'); end %debug
    end

    % for use of average prior in long. pipeline we have to deform the average mesh to current pp distance map        
    if useprior
      CS = loadSurf(Pcentral);
      correct_mesh = 1;
      cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none  0  1  -1  .1 ' ...   
                     'avg  %0.3f  %0.3f .2  .1  5  0 "0.5"  "0.5"  n 0  0  0 %d %g  0.0 0 %d'], ...          
                      Vppm.fname,Pcentral,Pcentral,-0.1,0.1,100,0.01,correct_mesh); 
      cat_system(cmd,opt.verb-3);
    end
    
    % get thickness data
    facevertexcdata = cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm); 
    cat_io_FreeSurfer('write_surf_data',Ppbt,facevertexcdata);
        
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

    copyfile(Ppbt,Pthick,'f');  
    res.(opt.surf{si}).createCS_final = cat_surf_fun('evalCS',loadSurf(Pcentral),cat_io_FreeSurfer('read_surf_data',Ppbt),[],Ymfs,Yppi,Pcentral,Smat.matlabIBB_mm,debug,cat_get_defaults('extopts.expertgui')>1);
    
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
  
  % calculate surface quality parameters for all surfaces
  mnth = []; sdth = []; mnRMSE_Ypp = []; mnRMSE_Ym = []; sdRMSE_Ym = []; sdRMSE_Ypp = []; 
  SIw = []; SIp = []; SIwa = []; SIpa = []; 
  for si=1:numel(opt.surf)
    if any(strcmp(opt.surf{si},{'lh','rh'}))
      if isfield(res.(opt.surf{si}).createCS_final,'fsthickness_mn_sd_md_mx')
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
    if any(~cellfun('isempty',strfind(opt.surf,'cb'))), cbtxt = 'cerebral '; else cbtxt = ''; end
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

              

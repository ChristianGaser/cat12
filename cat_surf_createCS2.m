function [Yth1,S,Psurf,EC,defect_size,res] = cat_surf_createCS2(V,V0,Ym,Ya,Yp0,YMF,opt)
% ______________________________________________________________________
% Surface creation and thickness estimation.
%
% [Yth1,S,Psurf,EC]=cat_surf_createCS(V,V0,Ym,Ya,YMF,opt)
%
% Yth1  = thickness map
% S     = structure with surfaces, like the left hemishere, that contains
%        vertices, faces, GM thickness (th1), and the transformation to
%        map to nifti space (vmat) and back (vmati).
% Psurf = name of surface files
% EC    = Euler characteristics
% defect_size = size of topology defects
% V     = spm_vol-structure of internally interpolated image
% V0    = spm_vol-structure of original image
% Ym    = the (local) intensity, noise, and bias corrected T1 image
% Ya    = the atlas map with the ROIs for left and right hemispheres
%        (this is generated with cat_vol_partvol)
% Yp0   = label image for surface deformation
% YMF   = a logical map with the area that has to be filled
%        (this is generated with cat_vol_partvol)
%   
% opt.surf       = {'lh','rh'[,'lc','rc']} - side
%    .reduceCS   = 100000 - number of faces
%
% Options set by cat_defaults.m
%    .interpV    = 0.5    - mm-resolution for thickness estimation
% 
% Here we used the intensity normalized image Ym, rather that the Yp0
% image, because it has more information about sulci that we need 
% especially for asymmetrical sulci.
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

%#ok<*AGROW,*STREMP,*ASGLU,*SFLD,*STFLD>

  % Turn off gifti data warning in gifti/subsref (line 45)
  %   Warning: A value of class "int32" was indexed with no subscripts specified. 
  %            Currently the result of this operation is the indexed value itself, 
  %            but in a future release, it will be an error. 
  warning('off','MATLAB:subscripting:noSubscriptsSpecified');
  cstime = clock; 

  % variables to tranfer from MATLAB to image coordinates used by loadSurf and saveSurf subfunctions
  global vmat vmati mati

  % surface evaluation parameter 
  res = struct('euler_characteristic',nan,'defect_size_promile',nan,'lh',struct(),'rh',struct()); 
  
  % set debuging variable
  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end
 
  % set defaults
  if ~exist('opt','var'), opt = struct(); end                 % create variable if not exist
  vx_vol        = sqrt(sum(V.mat(1:3,1:3).^2));               % further interpolation based on internal resolution 
  vx_vol0       = sqrt(sum(V0.mat(1:3,1:3).^2));              % final surface resolution based on original image resolution
  def.verb      = max(2,1 + cat_get_defaults('extopts.expertgui')); % 0-none, 1-minimal, 2-default, 3-details, 4-debug
  def.surf      = {'lh','rh'};                                % surface reconstruction setting with {'lh','rh','rc','lc'} and
                                                              % additional fast option 'fst' (e.g., 'lhfst') without topo.-corr. and sphere-reg. 
  
  % reducepatch has some issues with self intersections and should only be used for "fast" option
  % There is a new SPM approach spm_mesh_reduce that is maybe more robust. 
  % Higher resolutions are at least required for animal preprocessing that is given by cat_main.
  def.LAB                 = cat_get_defaults('extopts.LAB');  % brain regions 
  def.SPM                 = 0;                                % surface-reconstration based on SPM segmentation input (see cat_main)
  def.pbtlas              = 0;                                % myelination correction option (in development - not working correctly in all data, RD201907)  
  def.pbtmethod           = 'pbt2x';                          % projection-based thickness (PBT) estimation ('pbt2x' (with minimum setting) or 'pbt2')
  def.WMT                 = 0;                                % pbt-based WM/CSF width/depth/thickness estimation 
  def.sharpenCB           = 0;                                % sharpening function for the cerebellum (in development, RD2017?)
  def.thick_measure       = 1;                                % 0-PBT; 1-Tfs (Freesurfer method using mean(TnearIS,TnearOS))
  def.thick_limit         = 5;                                % 5mm upper limit for thickness (same limit as used in Freesurfer)
  def.new_release         = 0;                                % developer flag to test new functionality for new release (currently not used)
  def.collcorr            = 1;                                % correction of surface collisions: 0 - none; 1 - approach A; 2 - approach B, 3 - both, 4 - extract_pial_white
  def.fsavgDir            = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 
  def.add_parahipp        = cat_get_defaults('extopts.add_parahipp');
  def.scale_cortex        = cat_get_defaults('extopts.scale_cortex');
  def.close_parahipp      = cat_get_defaults('extopts.close_parahipp');
  def.write_debugsurfs    = cat_get_defaults('extopts.expertgui')>1;          % here we may need a better definition 
   
  % options that rely on other options
  opt                     = cat_io_updateStruct(def,opt); clear def; 
  opt.fast                = any(~cellfun('isempty',strfind(opt.surf,'fst'))); % fast registration without topo.-corr. and sphere-reg.
  opt.vol                 = any(~cellfun('isempty',strfind(opt.surf,'v')));   % only volume-based thickness estimation  
  opt.surf                = cat_io_strrep(opt.surf,{'sfst','fst','v'},'');    % after defintion of the 'fast' and 'vol' varialbe we simplify 'surf'
  opt.interpV             = max(0.1,min([opt.interpV,1.5]));                  % general limition of the PBT resolution
  opt.interpVold          = opt.interpV;                                      % save the default setting because interpV is later modified?
  opt.extract_pial_white          = opt.collcorr==1;                          % estimate pial and white matter surface (in development and very slow!)
  opt.force_no_selfintersections  = opt.extract_pial_white;                   % exact estimation requires this setting

  % distance between vertices that can be set directly by "vdist" or indirectly by "interpV"  
  % - surface should have more than 80k faces to look nice, whereas more than 400k does not improve the visual quality 
  % - controlled by power function to avoid a quatratic grow of the number of faces
  % - vdisto = [4 2 1 0.5] => [100k 200k 400k 800k] faces
  %max( 0.5 , min( 2 , min( opt.interpV , mean(vx_vol0) ))); % use sqare to use sqrt in general  
  if ~isfield(opt,'vdist') || opt.vdist == 0, opt.vdist  = 4/3; end
  opt.vdisto = opt.vdist; 
  opt.vdist  = sqrt(opt.vdist * 2); % here we use the sqrt to support linear mesh resolution increase (default input is [4 2 1 0.5])
  
  % Another parameter to control runtime is the accuracy of the surface
  % deformation. As far as we primary adapt the mesh resolution above, it 
  % is useful to use sqrt values rather than linear values to avoid squart
  % processing times for higher quality levels. Otherwise, we can simple 
  % avoid changes here ... so we can define this parameter utilizing the 
  % vdist parameter by simply divide it by 100.
  %def.surfaccuracy        = 0.01; % no adaption here otherwise processing time will not simply double 
  def.surfaccuracy = opt.vdist / 200; 
  def.reduceCS     = (300000 * sqrt(4/3 * 2) ) ./ opt.vdist; % to test ... fprintf('%g ',300000 * sqrt(1.3*2) ./ (( [4 2 1.3 1 0.5] * 2).^0.5))
  opt              = cat_io_updateStruct(def,opt);
  
  if opt.write_debugsurfs, opt.verb = 3; end
  
  % function to estimate the number of interations of the surface deformation: d=distance in mm and a=accuracy 
  moveth = @(d,a) [ round(d / a) , a ]; 
  
  if opt.verb>2
    fprintf('\nSurface reconstruction:              %s\n',....
      sprintf('%s',char( cellfun(@(x) [x ' '],opt.surf,'UniformOutput',0) )')); 
    fprintf('  PBT resolution:                    %0.3f\n',opt.interpV);
    fprintf('  lower face limit:                  %g\n',opt.reduceCS);
    fprintf('  maximal vertex distance:           %0.3f mm\n',opt.vdist);
    fprintf('  optimization stepsize:             %0.3f mm',opt.surfaccuracy);
  end  
  

  Psurf = struct(); 
  
  % correction for 'n' prefix for noise corrected and/or interpolated files
  [pp,ff]   = spm_fileparts(V.fname);

  if cat_get_defaults('extopts.subfolders')
    if strcmp(opt.pbtmethod,'pbt3')
      surffolder = sprintf('surf_%s_%0.2f',opt.pbtmethod,opt.interpV);
    elseif strcmp(opt.pbtmethod,'pbt2xf')
      opt.pbtmethod = 'pbt2x';
      surffolder = sprintf('surf_%s_%0.2f',opt.pbtmethod,opt.interpV);
    else
      surffolder = 'surf2';
    end
    mrifolder = 'mri';
    pp = spm_str_manip(pp,'h'); % remove 'mri' in pathname that already exists
    if ~exist(fullfile(pp,surffolder),'dir'), mkdir(fullfile(pp,surffolder)); end
  else
    surffolder = '';
    mrifolder = '';
  end

  if ff(1)=='n'
    if (exist(fullfile(pp,[ff(2:end) '.nii']), 'file')) || (exist(fullfile(pp,[ff(2:end) '.img']), 'file'))
      ff = ff(2:end);
    end
  end

  
  
  %% get both sides in the atlas map
  NS = @(Ys,s) Ys==s | Ys==s+1; 
    
  % noise reduction for higher resolutions (>=1 mm full correction, 1.5 mm as lower limit)
  % (added 20160920 ~R1010 due to servere sulcus reconstruction problems with 1.5 Tesla data)
  Yms = Ym + 0; cat_sanlm(Yms,3,1);
  mf  = min(1,max(0,3-2*mean(vx_vol,2))); 
  Ym  = mf * Yms  +  (1-mf) * Ym;
  clear Yms;
   
  % filling
  Ymf  = max(Ym,min(1,YMF)); 
  Ymfs = cat_vol_smooth3X(Ymf,1); 
  Ytmp = cat_vol_morph(YMF,'d',3) & Ymfs>2.3/3;
  Ymf(Ytmp) = max(min(Ym(Ytmp),0),Ymfs(Ytmp)); clear Ytmp Ymfs; 
  Ymf = Ymf*3;
  
  % removing fine WM structures in the hippocampus area to reduce topo and geometrie defects (added RD20190912)
  % use erode to reduce probability of cutting other gyri
  HCmask = cat_vol_morph( NS(Ya,opt.LAB.HC) , 'de', 1.5, vx_vol); 
  Ymf( HCmask ) =  min(2,Ymf( HCmask ));

  
  
  %% reduction of artifact, blood vessel, and meninges next to the cortex in SPM segmentations 
  %  (are often visible as very thin structures that were added to the WM 
  %  or removed from the brain)
  if ~opt.SPM
    Ydiv  = cat_vol_div(Ymf,vx_vol); 
    Ycsfd = cat_vbdist(single(Ymf<1.5),Ymf>1,vx_vol);
    Yctd  = cat_vbdist(single(Ymf<0.5),Ymf>0,vx_vol); 
    Ysroi = Ymf>2  &  Yctd<10  & Ycsfd>0 & Ycsfd<2 & ...
            cat_vol_morph(~NS(Ya,opt.LAB.HC) & ~NS(Ya,opt.LAB.HI) & ...
              ~NS(Ya,opt.LAB.PH) & ~NS(Ya,opt.LAB.VT),'erode',4); 
    Ybv   = cat_vol_morph(Ymf+Ydiv./max(1,Ymf)>3.5,'d') & Ymf>2; 
    Ymf(Ybv) = 1.4; 
    Ymfs  = cat_vol_median3(Ymf,Ysroi | Ybv,Ymf>eps & ~Ybv,0.1); % median filter
    Ymf   = mf * Ymfs  +  (1-mf) * Ymf;

    %% closing of small WMHs and blood vessels
    %vols = [sum(round(Ymf(:))==1 & Ya(:)>0) sum(round(Ymf(:))==2)  sum(round(Ymf(:))==3)] / sum(round(Ymf(:))>0); 
    %volt = min(1,max(0,mean([ (vols(1)-0.20)*5  (1 - max(0,min(0.3,vols(3)-0.2))*10) ]))); 
    %Ywmh = cat_vol_morph(Ymf>max(2.2,2.5 - 0.3*volt),'lc',volt); 
    %Ymf  = max(Ymf,smooth3(Ywmh)*2.9); 
  
    % gaussian filter? ... only in tissue regions
    %Ymfs = cat_vol_smooth3X(max(1,Ymf),0.5*min(1,max(0,1.5-mean(vx_vol)))); 
    %Ymf(Ymf>1) = Ymfs(Ymf>1);
  end
  if ~debug, clear Ysroi Ymfs Yctd Ybv Ymfs; end

  
  %% sharpening of thin structures in the cerebellum (gyri and sulci)
  % WARNING: this will change cortical thickness!
  if ~opt.SPM && opt.sharpenCB
    Ydiv = cat_vol_div(Ymf); %Ydivl  = cat_vol_div(Ymf,vx_vol); 
    Ywmd = cat_vbdist(single(Ymf>2.5),Ymf>1,vx_vol);
    if 0
      %% divergence based
      %  this works in principle but gyral crones and sulcal values are overestimated ... need limit
      Ymsk = (NS(Ya,opt.LAB.CB) & ((Ymf<2.8 & Ymf>2.0          ) | (Ymf<1.9 & Ymf>1.2         )) ) | ... sulci and gyri in the cerebellum 
             (NS(Ya,opt.LAB.CT) & ((Ymf<2.8 & Ymf>2.0 & Ycsfd>3) | (Ymf<1.9 & Ymf>1.2 & Ywmd>3)) ) | ... distant gyri and sulci in the cerebrum
             (NS(Ya,opt.LAB.PH) & ((Ymf<2.8 & Ymf>2.0 & Ycsfd>3) | (Ymf<1.9 & Ymf>1.2 & Ywmd>3)) );
      Ymf  = min(3,max( min(1,Ymf) , Ymf - (abs(Ydivl) .* Ydiv) .* Ymsk));
    end
    
    if 1
      %% biascorrection based
      % WM 
      Ymsk = ((NS(Ya,opt.LAB.CB) | YMF) & ( Ymf>2.2 | (Ymf>2 & Ydiv<-0.01) ) ) | ...          % sulci and gyri in the cerebellum 
             (NS(Ya,opt.LAB.PH) & ( Ymf>2.2 | (Ymf>2 & Ydiv<-0.01) ) ) | ...                  % hippocampal gyri
             (NS(Ya,opt.LAB.CT) & ( Ymf>2.2 | (Ymf>2 & Ydiv<-0.01 & ...
                Ycsfd>cat_stat_nanmean(Ycsfd(Ycsfd(:)>0 & Ycsfd(:)<100)) )*1.0) );            % distant gyri and sulci in the cerebrum
      Yi   = cat_vol_localstat(Ymf,Ymsk,1,3);
      % GM
      Ymsk = (NS(Ya,opt.LAB.CB) & ( Ymf>1.9 & Ymf<2.2 & Ycsfd>0 & Ydiv>-0.05) ) | ...         % sulci and gyri in the cerebellum 
             (NS(Ya,opt.LAB.PH) & ( Ymf>1.3 & Ymf<2.2 & Ycsfd>0 ) ) | ...                     % hippocampal gyri
             (NS(Ya,opt.LAB.CT) & ( Ymf>1.3 & Ymf<2.2 & Ycsfd>0 & ...
                Ywmd>cat_stat_nanmean(Ywmd(Ywmd(:)>0 & Ywmd(:)<100))*0.2 ) );                 % distant gyri and sulci in the cerebrum
      Yi   = Yi + cat_vol_localstat(Ymf,Yi==0 & Ymsk,1,1)/2*3;
      Yi   = cat_vol_localstat(Yi,Yi>0,1,3);
      Yi   = cat_vol_localstat(Yi,Yi>0,1,1); 
      if ~debug, clear Ywmd; end
      % CSF - instable and not required
      %Ymsk = NS(Ya,opt.LAB.VT) & Ymf>=0.5 & Ymf<1.5;                                         % sulci and gyri in the cerebellum 
      %Yi  = Yi + cat_vol_localstat(Ymf,Yi==0 & Ymsk,1,3)*3;
      
      Ywi = cat_vol_approx(Yi,'nn',1,2,struct('lfO',2)); 
      
      Ymf = Ymf./Ywi * 3; 
      if ~debug, clear Ywi Yi; end
    end
    if ~debug, clear Ymsk; end
  end
  if ~debug, clear Ydiv Ycsfd; end
  
  % initialize WM/CSF thickness/weidth/depth maps
  Yth1 = zeros(size(Ymf),'single'); 
  if opt.WMT > 1
    Ywd  = zeros(size(Ymf),'single'); 
    Ycd  = zeros(size(Ymf),'single'); 
  end
  
  % complete atlas map 
  [D,I] = cat_vbdist(single(Ya>0)); Ya = Ya(I); clear D;  
  
  % use sum of EC's and defect sizes for all surfaces, thus set values initially to 0
  EC = 0;
  defect_size = 0;

  
  % main loop for each surface structure 
  for si=1:numel(opt.surf)
   
    % surface filenames
    Praw       = fullfile(pp,surffolder,sprintf('%s.central.nofix.%s.gii',opt.surf{si},ff));    % raw
    Psphere0   = fullfile(pp,surffolder,sprintf('%s.sphere.nofix.%s.gii',opt.surf{si},ff));     % sphere.nofix
    Pcentral   = fullfile(pp,surffolder,sprintf('%s.central.%s.gii',opt.surf{si},ff));          % central
    Pcentralr  = fullfile(pp,surffolder,sprintf('%s.central.resampled.%s.gii',opt.surf{si},ff));% central
    Player4    = fullfile(pp,surffolder,sprintf('%s.layer4.%s.gii',opt.surf{si},ff));           % layer4
    Ppial      = fullfile(pp,surffolder,sprintf('%s.pial.%s.gii',opt.surf{si},ff));             % pial (GM/CSF)
    Pwhite     = fullfile(pp,surffolder,sprintf('%s.white.%s.gii',opt.surf{si},ff));            % white (WM/GM)
    Pthick     = fullfile(pp,surffolder,sprintf('%s.thickness.%s',opt.surf{si},ff));            % FS thickness / GM depth
    Ppbt       = fullfile(pp,surffolder,sprintf('%s.pbt.%s',opt.surf{si},ff));                  % PBT thickness / GM depth
    Pgwo       = fullfile(pp,surffolder,sprintf('%s.depthWMo.%s',opt.surf{si},ff));             % gyrus width / GWM depth / gyral span
    Pgw        = fullfile(pp,surffolder,sprintf('%s.depthGWM.%s',opt.surf{si},ff));             % gyrus width / GWM depth / gyral span
    Pgww       = fullfile(pp,surffolder,sprintf('%s.depthWM.%s',opt.surf{si},ff));              % gyrus witdh of the WM / WM depth
    Pgwwg      = fullfile(pp,surffolder,sprintf('%s.depthWMg.%s',opt.surf{si},ff));             % gyrus witdh of the WM / WM depth
    Psw        = fullfile(pp,surffolder,sprintf('%s.depthCSF.%s',opt.surf{si},ff));             % sulcus width / CSF depth / sulcal span
    Pdefects0  = fullfile(pp,surffolder,sprintf('%s.defects.%s',opt.surf{si},ff));              % defects temporary file
    Pmask      = fullfile(pp,surffolder,sprintf('%s.mask.%s',opt.surf{si},ff));                 % mask
    Ptemp      = fullfile(pp,surffolder,sprintf('%s.temp.%s',opt.surf{si},ff));                 % temporary file
    Pdefects   = fullfile(pp,surffolder,sprintf('%s.defects.%s.gii',opt.surf{si},ff));          % defects
    Psphere    = fullfile(pp,surffolder,sprintf('%s.sphere.%s.gii',opt.surf{si},ff));           % sphere
    Pspherereg = fullfile(pp,surffolder,sprintf('%s.sphere.reg.%s.gii',opt.surf{si},ff));       % sphere.reg
    Pfsavg     = fullfile(opt.fsavgDir, sprintf('%s.central.freesurfer.gii',opt.surf{si}));     % fsaverage central
    Pfsavgsph  = fullfile(opt.fsavgDir, sprintf('%s.sphere.freesurfer.gii',opt.surf{si}));      % fsaverage sphere    
    Pfsavgmask = fullfile(opt.fsavgDir, sprintf('%s.mask',opt.surf{si}));                       % fsaverage mask    
    
    surffile = {'Praw','Psphere0','Pcentral','Pthick','Ppbt','Pgw','Pgww','Psw',...
      'Pdefects0','Pdefects','Psphere','Pspherereg','Pfsavg','Pfsavgsph','Pwhite','Ppial'};
    for sfi=1:numel(surffile)
      eval(sprintf('Psurf(si).%s = %s;',surffile{sfi},surffile{sfi})); 
    end
        
    % reduce for object area
    switch opt.surf{si}
      case {'lh'},  Ymfs = Ymf .* (Ya>0) .* ~(NS(Ya,opt.LAB.CB) | NS(Ya,opt.LAB.BV) | NS(Ya,opt.LAB.ON) | NS(Ya,opt.LAB.MB)) .* (mod(Ya,2)==1); Yside = mod(Ya,2)==1; 
      case {'rh'},  Ymfs = Ymf .* (Ya>0) .* ~(NS(Ya,opt.LAB.CB) | NS(Ya,opt.LAB.BV) | NS(Ya,opt.LAB.ON) | NS(Ya,opt.LAB.MB)) .* (mod(Ya,2)==0); Yside = mod(Ya,2)==0;  
      case {'lc'},  Ymfs = Ymf .* (Ya>0) .*   NS(Ya,opt.LAB.CB).* (mod(Ya,2)==1); Yside = mod(Ya,2)==1; 
      case {'rc'},  Ymfs = Ymf .* (Ya>0) .*   NS(Ya,opt.LAB.CB).* (mod(Ya,2)==0); Yside = mod(Ya,2)==0; 
    end 
    
    switch opt.surf{si}
      case {'lh','rh'}, opt.interpV = opt.interpVold; 
      case {'lc','rc'}, opt.interpV = opt.interpVold / 2 ; 
    end 
    
    % check for cerebellar hemis 
    iscerebellum = strcmp(opt.surf{si},'lc') || strcmp(opt.surf{si},'rc');
    
    % scaling factor for reducing patches and refinement for cerebellar  
    % hemis according to voxel size or 1 for cerebrum
    scale_cerebellum  = 1 + (iscerebellum * max(1,min(3,1/mean(vx_vol,2))));
    
    % get dilated mask of gyrus parahippocampalis and hippocampus of both sides
    if ~iscerebellum
      mask_parahipp = cat_vol_morph(NS(Ya,opt.LAB.PH) | NS(Ya,opt.LAB.HC),'d',6);
    end
    
    %% thickness estimation
    
    % print something
    if si==1, fprintf('\n'); end
    switch opt.fast
      case 1, fprintf('%s - fast without registration:\n',opt.surf{si});
      case 0, fprintf('%s:\n',opt.surf{si});
    end
    stime = cat_io_cmd(sprintf('  Thickness estimation (%0.2f mm%s)',opt.interpV,char(179))); stimet =stime;
    
    % removing background (smoothing to remove artifacts)
    switch opt.surf{si}
      case {'lh','rh'},  [Ymfs,Yside,mask_parahipp,BB] = cat_vol_resize({Ymfs,Yside,mask_parahipp},'reduceBrain',vx_vol,4,smooth3(Ymfs)>1.5); 
      case {'lc','rc'},  [Ymfs,Yside,BB] = cat_vol_resize({Ymfs,Yside},'reduceBrain',vx_vol,4,smooth3(Ymfs)>1.5); 
    end
    
    imethod         = 'cubic'; % cubic should be better in general - however, linear is better for small thickness (version?)
    [Ymfs,resI]     = cat_vol_resize(max(1,Ymfs),'interp',V,opt.interpV,imethod);                  % interpolate volume
    Yside           = cat_vol_resize(Yside,'interp',V,opt.interpV,imethod)>0.5;                    % interpolate volume (small dilatation)
    
    if ~iscerebellum
      mask_parahipp   = cat_vol_resize(mask_parahipp,'interp',V,opt.interpV)>0.5;          % interpolate volume
    end 
    
    % PVE with background will lead to a light underestimation?
    Ymfs = min(3,max(1,Ymfs));

    
    %% pbt calculation
    if strcmp(opt.pbtmethod,'pbt3')
      [Yth1i,Yppi] = cat_vol_pbt3(Ymfs,struct('method',opt.pbtmethod,'cb',iscerebellum,'resV',opt.interpV,...
        'vmat',V.mat(1:3,:)*[0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1])); % avoid underestimated thickness in gyri
    else
      [Yth1i,Yppi,Ymfs] = cat_vol_pbt(Ymfs,struct('method',opt.pbtmethod,'resV',opt.interpV,'vmat',...
        V.mat(1:3,:)*[0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1],'pbtlas',opt.pbtlas)); % avoid underestimated thickness in gyri
    end  
    %%
    if ~opt.WMT && ~debug, clear Ymfs; end
    Yth1i(Yth1i>10)=0; Yppi(isnan(Yppi))=0;  
    [D,I] = cat_vbdist(Yth1i,Yside); Yth1i = Yth1i(I); clear D I;       % add further values around the cortex
    Yth1t = cat_vol_resize(Yth1i,'deinterp',resI); clear Yth1i;         % back to original resolution
    Yth1t = cat_vol_resize(Yth1t,'dereduceBrain',BB);                   % adding background
    Yth1  = max(Yth1,Yth1t);                                            % save on main image
    clear Yth1t;
    
    if opt.vol
      S = struct(); Psurf = '';
      fprintf('%5.0fs',etime(clock,stime)); 
      continue
    end
    
    %% PBT estimation of the gyrus and sulcus width 
    if opt.WMT > 1 
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
        case {'lh'} 
          Ymr = Ymr .* (Yar>0) .* ~(NS(Yar,3) | NS(Yar,7) | NS(Yar,11) | NS(Yar,13)) .* (mod(Yar,2)==1);
          Ynw = smooth3(cat_vol_morph(NS(Yar,5) | NS(Yar,9) | NS(Yar,15) | NS(Yar,23),'d',2) | ...
                 (cat_vol_morph(Yppi==1,'e',2) & Ymr>1.7/3 & Ymr<2.5/3) & (mod(Yar,2)==1)); 
        case {'rh'}
          Ymr = Ymr .* (Yar>0) .* ~(NS(Yar,3) | NS(Yar,7) | NS(Yar,11) | NS(Yar,13)) .* (mod(Yar,2)==0);    
          Ynw = smooth3(cat_vol_morph(NS(Yar,5) | NS(Yar,9) | NS(Yar,15) | NS(Yar,23),'d',2) | ...
                 (cat_vol_morph(Yppi==1,'e',2) & Ymr>1.7/3 & Ymr<2.5/3) & (mod(Yar,2)==0)); 
        case {'lc'}
          Ymr = Ymr .* (Yar>0) .* NS(Yar,3) .* (mod(Yar,2)==1);
        case {'rc'}
          Ymr = Ymr .* (Yar>0) .* NS(Yar,3) .* (mod(Yar,2)==0);
      end 
     % clear Yar; 
      %%
      Yppis = Yppi .* (1-Ynw) + max(0,min(1,Ymr*3-2)) .* Ynw;                         % adding real WM map 
      Ywdt  = cat_vol_eidist(1-Yppis,ones(size(Yppis),'single'));                     % estimate distance map to central/WM surface
      Ywdt  = cat_vol_pbtp(max(2,4-Ymfs),Ywdt,inf(size(Ywdt),'single'))*opt.interpV;
      [D,I] = cat_vbdist(single(Ywdt>0.01),Yppis>0); Ywdt = Ywdt(I); clear D I Yppis; % add further values around the cortex
      Ywdt  = cat_vol_median3(Ywdt,Ywdt>0.01,Ywdt>0.01);                    
      Ywdt = cat_vol_localstat(Ywdt,Ywdt>0.1,1,1);     % smoothing
      Ywdt  = cat_vol_resize(Ywdt,'deinterp',resI);                                   % back to original resolution
      Ywdt  = cat_vol_resize(Ywdt,'dereduceBrain',BB);                                % adding background
      Ywd   = max(Ywd,Ywdt); 
      clear Ywdt;
      
      %% sulcus width / CSF depth
      %  for the CSF depth we cannot use the origal data, because of
      %  sulcal blurring, but we got the PP map at half distance and
      %  correct later for half thickness
      fprintf('%5.0fs',etime(clock,stime)); 
      stime = cat_io_cmd('  CSF depth estimation');
      YM    = single(smooth3(cat_vol_morph(Ymr<0.1,'o',4))<0.5); YM(YM==0)=nan;       % smooth CSF/background-skull boundary 
      Yppis = Yppi .* ((Ymr+0.25)>Yppi) + min(1,Ymr*3-1) .* ((Ymr+0.25)<=Yppi);       % we want also CSF within the ventricle (for tests)
      Ycdt  = cat_vol_eidist(Yppis,YM);                                               % distance to the cental/CSF-GM boundary
      Ycdt  = cat_vol_pbtp(max(2,Ymfs),Ycdt,inf(size(Ycdt),'single'))*opt.interpV; Ycdt(isnan(Ycdt))=0;
      [D,I] = cat_vbdist(single(Ycdt>0),Yppis>0 & Yppis<3); Ycdt = Ycdt(I); clear D I Yppis; % add further values around the cortex
      Ycdt  = cat_vol_median3(Ycdt,Ycdt>0.01,Ycdt>0.01);                              % median filtering
      Ycdt = cat_vol_localstat(Ycdt,Ycdt>0.1,1,1);                                    % smoothing
      Ycdt  = cat_vol_resize(Ycdt,'deinterp',resI);                                   % back to original resolution
      Ycdt  = cat_vol_resize(Ycdt,'dereduceBrain',BB);                                % adding background
      Ycd   = max(Ycd,Ycdt); 
      clear Ycdt;
      %fprintf('%5.0fs',etime(clock,stime));
      clear Ymr;
    end
    
    %if ~debug, clear Ymfs; else Yppio=Yppi; end
    %if debug, Yppio=Yppi; end
    fprintf('%5.0fs\n',etime(clock,stime));
    
    %% Replace isolated voxels and holes in Ypp by its median value
    
    % indicate isolated holes and replace by median of the neighbors
    Yppi(Yppi<0.35 & ~cat_vol_morph(Yppi<1,'l'))=1;  % close major wholes in the WM 
    Ymsk = Yppi==0 & cat_vol_morph(Yppi>0.9,'d',1); % filter small wholes close to the WM
    Yppi = cat_vol_median3(single(Yppi),Ymsk,~Ymsk); 
    
    %% indicate isolated objects and replace by median of the neighbors
    Yppi(Yppi>0.65 & cat_vol_morph(Yppi==0,'l'))=0;
    Ymsk = Yppi>0.95 & cat_vol_morph(Yppi<0.1,'d',1); 
    Yppi = cat_vol_median3(single(Yppi),Ymsk,~Ymsk);
    if ~debug, clear Ymsk; end
    
    %% Write Ypp for final deformation
    %  Write Yppi file with 1 mm resolution for the final deformation, 
    %  because CAT_DeformSurf achieved better results using that resolution
    Yppt = cat_vol_resize(Yppi,'deinterp',resI);                        % back to original resolution
    Yppt = cat_vol_resize(Yppt,'dereduceBrain',BB);                     % adding of background

    % scale Yppt so that backgrounds remains 0 and WM 1, but cortical band is 
    % now in the range of 0.1..0.9
    if opt.extract_pial_white
      % mask hemispheres and regions
      switch opt.surf{si}
        case {'lh'},  Yp0s = Yp0 .* (Ya>0) .* ~(NS(Ya,opt.LAB.CB) | NS(Ya,opt.LAB.ON) | NS(Ya,opt.LAB.MB)) .* (mod(Ya,2)==1);  
        case {'rh'},  Yp0s = Yp0 .* (Ya>0) .* ~(NS(Ya,opt.LAB.CB) | NS(Ya,opt.LAB.ON) | NS(Ya,opt.LAB.MB)) .* (mod(Ya,2)==0);   
        case {'lc'},  Yp0s = Yp0 .* (Ya>0) .*   NS(Ya,opt.LAB.CB).* (mod(Ya,2)==1);  
        case {'rc'},  Yp0s = Yp0 .* (Ya>0) .*   NS(Ya,opt.LAB.CB).* (mod(Ya,2)==0); 
      end
      %Ymfs = cat_vol_resize(Ymfs,'deinterp',resI); clear Yth1i;         % back to original resolution
      %Ymfs = cat_vol_resize(Ymfs,'dereduceBrain',BB);                   % adding background
    
      Vyp0s  = cat_io_writenii(V,Yp0s,'','yp0s','scaled image','uint8',[0,1/255],[1 0 0 0]);

      indi = find((Yppt>0) & (Yppt<0.99999));
      Yppt(indi) = 0.1 + (0.8*Yppt(indi));
    end
    Vpp  = cat_io_writenii(V,Yppt,'','pp' ,'percentage position map','uint8',[0,1/255],[1 0 0 0]);
    Vmcs = cat_io_writenii(V,max(1/3,min(1,Ym)),'','mcs','intensity normalized map for central surface recontruction','uint8',[0,1/255],[1 0 0 0]);
    
    Vpp1 = Vpp; 
    Vpp1.fname    = fullfile(pp,mrifolder,['pp1' ff '.nii']);
    vmat2         = spm_imatrix(Vpp1.mat);
    Vpp1.dim(1:3) = round(Vpp1.dim .* abs(vmat2(7:9)*(1 + iscerebellum)));   % use double resolution in case of cerebellum
    vmat2(7:9)    = sign(vmat2(7:9)).*[1 1 1]/(1 + iscerebellum);            % use double resolution in case of cerebellum
    Vpp1.mat      = spm_matrix(vmat2);
    Vpp1 = spm_create_vol(Vpp1); 
    for x3 = 1:Vpp1.dim(3)
      M    = inv(spm_matrix([0 0 -x3 0 0 0 1 1 1]) * inv(Vpp1.mat) * Vpp.mat); %#ok<MINV>
      v    = spm_slice_vol(Vpp,M,Vpp1.dim(1:2),1);       
      Vpp1 = spm_write_plane(Vpp1,v,x3);
    end
    clear M v x3; 

    
    
    %% surface creation 
    %  --------------------------------------------------------------------
    %  Surface create should be at 0.5 mm to support a useful description
    %  in even narrow sulci, i.e., for 1 mm thickness the sulci would be 
    %  2 mm wide were a 1 mm percental position map is strongly limited.
    %  The surface is reconstructed by marching cubes on a binary version 
    %  of the PBT radial position map Yppi to use a simple voxelbased
    %  topology correction to avoid small defects and incorrected
    %  triangulation (e.g., matlab isosurface function). Next, the
    %  complexity of the surface is reduced by the spm_mesh_reduce function
    %  (more accurate and stable than the matlab reduce function?) to
    %  improve the performance of the following main topology correction. 
    %  However, the reduction can partially lead to very large faces that
    %  need an additianal refine. 
    %  --------------------------------------------------------------------
    if opt.verb==1, fprintf('%s %4.0fs\n',repmat(' ',1,66),etime(clock,stimet)); end
    stime = cat_io_cmd('  Create initial surface','g5','',opt.verb); %if opt.verb>2, fprintf('\n'); end
    
    
    % surface coordinate transformations that are used in the "saveCS" and "loadCS" functions  
    mati  = spm_imatrix(V.mat); 
    vmat  = V.mat(1:3,:)*[0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];
    vmati = inv([vmat; 0 0 0 1]); vmati(4,:) = [];    

    
    % scaling correction to reduce topology errors in the parahippocampus
    % ### not tested for cerebellum yet (RD201909) ###
    if ~iscerebellum
      ind0 = find(Yppi<=0);
      Yppi = opt.scale_cortex * Yppi;
      
      % smooth mask to have smooth border
      mask_parahipp_smoothed = zeros(size(mask_parahipp));
      spm_smooth(double(mask_parahipp),mask_parahipp_smoothed,[4 4 4] / opt.interpV);
      Yppi  = Yppi + opt.add_parahipp/opt.scale_cortex * mask_parahipp_smoothed;
      Yppi(ind0) = 0; clear ind0; %#ok<FNDSB>
      
      % optionally apply closing inside mask for parahippocampal gyrus to get rid 
      % of the holes that lead to large cuts in gyri after topology correction
      if opt.close_parahipp
        tmp = cat_vol_morph(Yppi,'labclose',1);
        Yppi(mask_parahipp) = tmp(mask_parahipp); %#ok<NASGU>
        if ~debug, clear tmp; end 
      end
    end
    
    % Marching cubes surface creation and correction for the boundary box 
    % used within the surface creation process.  It is better to use the 
    % full voxel resolution in combination with surface reduction rather   
    % then using lower voxel resolutions! Moreover, it was NOT necessary to 
    % use surface deformation before the surface reduction!
    evalc('clear CS; [tmp,CS.faces,CS.vertices] = cat_vol_genus0(Yppi,0.5);');
    CS.vertices = CS.vertices .* repmat(abs(opt.interpV ./ mati([8,7,9])),size(CS.vertices,1),1);
    CS.vertices = CS.vertices +  repmat( BB.BB([3,1,5]) - 1,size(CS.vertices,1),1);
    %if ~debug, clear Yppi; end

    % reduce resolution with higher resolution for cerebellum and fast option
    CS = spm_mesh_reduce(CS, 81920 / (1 + (opt.vdist>2)) * scale_cerebellum * ( 1 + opt.fast) ); 
    saveSurf(CS,Praw);

    % remove unconnected meshes
    cmd = sprintf('CAT_SeparatePolygon "%s" "%s" -1',Praw,Praw); 
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
     
    % refine super-large faces with adaption for cerebellum and fast option
    cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f',Praw,Praw,3 / scale_cerebellum / ( 1 + opt.fast) ); 
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

    % Create a smooth surface for the topology correction. 
    % It don't has to be perfect because it will replaced completelly!
    for li = 2.^(0:2) 
      cmds = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none  0  1  -1  .1 ' ...   
                     'avg  %0.3f  %0.3f .2  .1  2  0 "0.5"  "0.5"  n 0  0  0 %d %g  0.0 0'], ...          
                      Vpp1.fname,Praw,Praw,-0.2/li,0.2/li,moveth(1/li,opt.surfaccuracy*4/li));
    end
    [ST, RS] = cat_system(cmds); cat_check_system_output(ST,RS,opt.verb-3);
      
    % load surf and map thickness
    CS = loadSurf(Praw);
    facevertexcdata = isocolors2(Yth1,CS.vertices); 
    

    
    if opt.fast 
    %% Fast processing without topology correction and sperical registration
    %  --------------------------------------------------------------------
    %  The one and only fast option that is equal to the init surface but 
    %  with collision correction. For fast surface and thickness outputs 
    %  for visual analysis. 
    %  --------------------------------------------------------------------
   
        % correction for surface collision of the IS and OS
        if opt.collcorr > 1      
          stime = cat_io_cmd(sprintf('  Correction of surface collisions:'),'g5','',opt.verb,stime); 
          CS1.vertices = CS.vertices; CS1.faces = CS.faces; CS1.vmat = vmat; CS1.vmati = vmati; CS1.mati = mati; 
          [CS1,facevertexcdata] = cat_surf_fun('collisionCorrection',CS1,...
            facevertexcdata,Ymf,Yppt,cat_get_defaults('extopts.expertgui')==2,Pcentral);
          CS.vertices = CS1.vertices; CS.faces = CS1.faces; clear CS1;
          collcorrstr = 'collcorr'; 
        else
          collcorrstr = '';
        end

        % save final data and datastructure
        saveSurf(CS,Pcentral);
        cat_io_FreeSurfer('write_surf_data',Ppbt,facevertexcdata);
        S.(opt.surf{si}) = struct('faces',CS.faces,'vertices',CS.vertices,'vmat',vmat,...
            'vmati',vmati,'mati',mati,'th1',facevertexcdata);

        % create white and central surfaces
        cat_surf_fun('white',Pcentral);
        cat_surf_fun('pial',Pcentral);

        % evaluate and save results
        CS1.vertices = CS.vertices; CS1.faces = CS.faces; CS1.vmat = vmat; CS1.vmati = vmati; CS1.mati = mati; 
        fprintf('%5.0fs',etime(clock,stime)); stime = []; 
        if opt.write_debugsurfs
          cat_surf_fun('saveico',CS1,isocolors2(Yth1,CS1.vertices),Pcentral,...
            sprintf('createCS_1_initfast_%s_pbtres%0.2fmm_vdist%0.2fmm',collcorrstr,opt.interpVold,opt.vdist)); 
        end     
        res.(opt.surf{si}).createCS_1_initfast = cat_surf_fun('evalCS',CS1,isocolors2(Yth1,CS1.vertices),Ym,Yppt,Pcentral,opt.verb-2);
        clear CS1 

          
        %  map WM and CSF width data (corrected by thickness)
        %  cat_surf_parameters and removed here (RD201909)
        if opt.WMT > 1
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
            try %#ok<TRYNC>
              cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s"',Pcentral,Pgwwg,3,Pgwwg);
              [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
            end

            facevertexcdata3 = isocolors2(Ycd,CS.vertices);
            facevertexcdata3 = max(eps,facevertexcdata3 - facevertexcdata/2); 
            cat_io_FreeSurfer('write_surf_data',Psw,facevertexcdata3);

            setfield(S.(opt.surf{si}),'th2',nan(size(facevertexcdata)));  
            setfield(S.(opt.surf{si}),'th3',nan(size(facevertexcdata)));
        end

        if ~debug
          delete(Vpp.fname);
          delete(Vpp1.fname);
        end

        % estimate Euler characteristics: EC = #vertices + #faces - #edges
        EC0 = size(CS.vertices,1) + size(CS.faces,1) - size(spm_mesh_edges(CS),1);
        EC  = EC + abs(EC0);


        % estimate Freesurfer thickness measure Tfs using mean(Tnear1,Tnear2)
        if opt.thick_measure == 1
          stime = cat_io_cmd('  Tfs thickness estimation:','g5','',opt.verb,stime);
          cmd = sprintf('CAT_SurfDistance -mean -thickness "%s" "%s" "%s"',Ppbt,Pcentral,Pthick);
          [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

          % apply upper thickness limit
          facevertexcdata = cat_io_FreeSurfer('read_surf_data',Pthick);  
          facevertexcdata(facevertexcdata > opt.thick_limit) = opt.thick_limit;
          cat_io_FreeSurfer('write_surf_data',Pthick,facevertexcdata);  
        else % otherwise simply copy ?h.pbt.* to ?h.thickness.*
          copyfile(Ppbt,Pthick);
        end

        fprintf('%5.0fs',etime(clock,stime)); 
        cat_surf_fun('evalCS',CS,facevertexcdata,Ym,Yppt,Pcentral,opt.verb-2);
        
        clear CS
        continue
    end
    %% evaluate and save results
    CS1.vertices = CS.vertices; CS1.faces = CS.faces; CS1.vmat = vmat; CS1.vmati = vmati; CS1.mati = mati; 
    fprintf('%5.0fs',etime(clock,stime)); stime = []; 
    if opt.write_debugsurfs
      cat_surf_fun('saveico',CS1,isocolors2(Yth1,CS1.vertices),Pcentral,sprintf('createCS_1_init_pbtres%0.2fmm_vdist%0.2fmm',opt.interpVold,opt.vdist)); 
    end     
    res.(opt.surf{si}).createCS_init = cat_surf_fun('evalCS',CS1,isocolors2(Yth1,CS1.vertices),Ym,Yppt,Pcentral,opt.verb-2);
    clear CS1 
    
    
    
        
    
    %% Topology correction and surface refinement
    %  --------------------------------------------------------------------
    %  This topology correction creates a complettely new surface based on  
    %  spherical hormonic functions resulting in a relative unbalanced
    %  local resolution (i.e., oversampled in the insula) that is corrected
    %  in the next block.  However, this also means the resoltion of the
    %  input surface don't have to be super high (see above). 
    %  --------------------------------------------------------------------
    stime = cat_io_cmd('  Topology correction:','g5','',opt.verb,stime); 
    
    % spherical surface mapping 1 of the uncorrected surface for topology correction
    cmd = sprintf('CAT_Surf2Sphere "%s" "%s" 5',Praw,Psphere0);
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

    % estimate size of topology defects 
    cmd = sprintf('CAT_MarkDefects "%s" "%s" "%s"',Praw,Psphere0,Pdefects0); 
    [ST, RS] = cat_system(cmd);
    defect_sizes = cat_io_FreeSurfer('read_surf_data',Pdefects0);
    defect_size0 = sum(defect_sizes > 0) / length(defect_sizes) * 100; % percent
    defect_area0 = sum(defect_sizes > 0) / length(defect_sizes) .* ...
      sum(cat_surf_fun('area',CS)) / opt.interpV / 100; % cm2
    
    defect_size  = defect_size + defect_size0;
    
    % estimate Euler characteristics: EC = #vertices + #faces - #edges
    EC0 = size(CS.vertices,1) + size(CS.faces,1) - size(spm_mesh_edges(CS),1);
    EC  = EC + abs(EC0);
    
    % topology correction and surface refinement
    % Higher -n will result in larger but still unbalanced meshes and the 
    % refine_lenght parameter is more important to obtain nice meshes.
    if opt.verb>3, fprintf('\n'); end
    cmd = sprintf('CAT_FixTopology -lim 512 -bw 1024 -n %d -refine_length %g "%s" "%s" "%s"',...
      81920, opt.vdist / scale_cerebellum,Praw,Psphere0,Pcentral);
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
    
    
    
 
    %% Optimize topology corrected mesh: 
    %  --------------------------------------------------------------------
    %  Oversampling can lead to problems in the normal transformation to 
    %  obtain the inner and outer surface, especially in the Insula/Amygdala.
    %  The idea is to reduce the sampling of the mesh and then to refine 
    %  the mesh again. 
    %  However it is highly essential that most regions where not corrected
    %  and a specific masking of the Insula (relative small triangles in a
    %  specific area on one/all of the main cortical surfaces or flipping) 
    %  is possible useful.
    %  ... seams that this is working and it takes only a few seconds!
    %  --------------------------------------------------------------------
    stime = cat_io_cmd('  Surface optimization and refinement:','g5','',opt.verb,stime); 

    % refinement - important for sulci .. here we need a lot of details with a similar resolution as the Insula 
    cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f',Pcentral,Pcentral,0.8 / scale_cerebellum); 
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

    % surface refinement (this time even before reduction)
    for li = 2.^(0:1)
      cmds = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none  0  1  -1  .1 ' ...           
                     'avg  -0.1  0.1 .2  .1  %d  0 "0.5"  "0.5"  n 0  0  0 %d  %0.2f  0.0 0'], ...    
                      Vpp1.fname,Pcentral,Pcentral,1/li,moveth(0.4/li,opt.surfaccuracy * 4/li));
      [ST, RS] = cat_system(cmds); cat_check_system_output(ST,RS,opt.verb-3);
    end
    
    % reduce - as far as the Insula/Amygdala is not so heavily folded compared to sulci this region is first reduced 
    CS = loadSurf(Pcentral);
    CS = spm_mesh_reduce(struct('vertices',CS.vertices,'faces',CS.faces),...
      min( max( 81920 , opt.reduceCS/2 ) , 81920 * 4 ) * scale_cerebellum);
    saveSurf(CS,Pcentral); 

     % remove unconnected meshes
    cmd = sprintf('CAT_SeparatePolygon "%s" "%s" -1',Pcentral,Pcentral); 
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
   
    % refinement - guaranty our default resolution
    cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f',Pcentral,Pcentral,opt.vdist / scale_cerebellum);
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
    
    % surface defomation for relaxation after reduction and refinement
    for li = 2.^(1:2)
      cmds = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none  0  1  -1  .1 ' ...           
                     'avg  -0.1  0.1 .2  .1  %d  0 "0.5"  "0.5"  n 0  0  0 %d  %0.2f  0.0 0'], ...    
                      Vpp1.fname,Pcentral,Pcentral,1/li,moveth(0.4/li,opt.surfaccuracy*4/li));
      [ST, RS] = cat_system(cmds); cat_check_system_output(ST,RS,opt.verb-3);
    end 
    
    % read final surface and map thickness data
		CS = loadSurf(Pcentral);
    facevertexcdata = isocolors2(Yth1,CS.vertices); 
    cat_io_FreeSurfer('write_surf_data',Ppbt,facevertexcdata);

    
    
    
    
    %% final correction of central surface in highly folded areas 
    %  with high mean curvature with weight of 0.7 and further refinement
    %  of the mesh and its vertices based on the position map 
    cmd = sprintf('CAT_Central2Pial -equivolume -weight 0.7 "%s" "%s" "%s" 0',Pcentral,Ppbt,Pcentral);
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

    % we need some refinement because some vertices are too large to be deformed with high accuracy
    cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f 1',Pcentral,Pcentral,opt.vdist / scale_cerebellum); % adaption for cerebellum
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

    % surface refinement by surface deformation based on the PP map
    for li = 2.^(0:1)
      cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .1 ' ...
                     'avg -0.1 0.1 .2 .1 %d 0 "0.5" "0.5" n 0 0 0 %d %0.2f 0.0 0'], ...
                     Vpp.fname,Pcentral,Pcentral,1/li,moveth(0.2,opt.surfaccuracy * 2/li)); 
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
    end
    
    % need some more refinement because some vertices are distorted after CAT_DeformSurf
    cmd = sprintf('CAT_RefineMesh "%s" "%s" %0.2f 1',Pcentral,Pcentral,opt.vdist / scale_cerebellum); % adaption for cerebellum
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

    % final surface refinement
    for li = 2.^(0:1)
      cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .1 ' ...
                     'avg -0.05 0.05 .5 .1 %d 0 "0.5" "0.5" n 0 0 0 %d %0.2f 0.0 0'], ...
                     Vpp.fname,Pcentral,Pcentral,1/li,moveth(0.1/li,opt.surfaccuracy * 2/li)); 
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
    end

    % read final surface and map thickness data
    CS = loadSurf(Pcentral);
    facevertexcdata = isocolors2(Yth1,CS.vertices); 
    cat_io_FreeSurfer('write_surf_data',Ppbt,facevertexcdata);
   
    % evaluate and save results
    CS1.vertices = CS.vertices; CS1.faces = CS.faces; CS1.vmat = vmat; CS1.vmati = vmati; CS1.mati = mati; 
    fprintf('%5.0fs',etime(clock,stime)); stime = []; 
    if opt.write_debugsurfs
      cat_surf_fun('saveico',CS1,isocolors2(Yth1,CS1.vertices),Pcentral,sprintf('createCS_2_refined_pbtres%0.2fmm_vdist%0.2fmm',opt.interpVold,opt.vdist)); 
    end
    res.(opt.surf{si}).createCS_2_refine = cat_surf_fun('evalCS' ,CS1,isocolors2(Yth1,CS1.vertices),Ym,Yppt,Pcentral,opt.verb-2);
    res.(opt.surf{si}).createCS_final    = res.(opt.surf{si}).createCS_2_refine; 
    clear CS1
    
    
    
    
    %% Collision correction by Delaunay triangularization
    if opt.collcorr > 1  
      stime = cat_io_cmd('  Correction of surface collisions:','g5','',opt.verb,stime); 

      % final correction of central surface in highly folded areas with high mean curvature
      cmd = sprintf('CAT_Central2Pial -equivolume -weight 1 "%s" "%s" "%s" 0',Pcentral,Ppbt,Player4);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,0);
      Sl4 = loadSurf(Player4);
      Yl4 = isocolors2(Ymf,Sl4.vertices); delete(Player4), clear Sl4; 
      
      %% call collision correction
      if exist('CSO','var'), CS = CSO; facevertexcdata = facevertexcdatao; else, CSO = CS; facevertexcdatao = facevertexcdata; end
      CS1.vertices = CS.vertices; CS1.faces = CS.faces; CS1.vmat = vmat; CS1.vmati = vmati; CS1.mati = mati; clear CS;  
      [CS1,facevertexcdata,E] = cat_surf_fun('collisionCorrection',CS1,facevertexcdata,Ymf,Yppt,Yl4); % do not write
      %[CS1,facevertexcdata,E] = cat_surf_fun('collisionCorrection',CS1,facevertexcdata,max(1,Ym*3),Yppt,Yl4); % do not write
      CS.vertices = CS1.vertices; CS.faces = CS1.faces; clear CS1; 
      if ~debug, clear Yl4; end
      
      %saveSurf(CS,Pcentral); 
      cat_io_FreeSurfer('write_surf_data',Ppbt,facevertexcdata);
      
      % evaluate and save results
      CS1.vertices = CS.vertices; CS1.faces = CS.faces; CS1.vmat = vmat; CS1.mati = mati; 
      fprintf('%5.0fs',etime(clock,stime)); stime = []; 
      if opt.write_debugsurfs
        cat_surf_fun('saveico',CS1,facevertexcdata,Pcentral,sprintf('createCS_3_collcorr_%0.2fmm_vdist%0.2fmm',opt.interpVold,opt.vdist)); 
      end
      res.(opt.surf{si}).createCS_3_collcorr = cat_surf_fun('evalCS' ,CS1,facevertexcdata,Ym,Yppt,Pcentral,opt.verb-2);
      res.(opt.surf{si}).createCS_final    = res.(opt.surf{si}).createCS_3_collcorr; 
      clear CS1;
      
      
      
      
    elseif opt.collcorr == 1
    %% final correction of cortical thickness using pial and WM surface
    %  This does not work right now. Although, there are maybe no self
    %  intersections, the deformation to pial and white are still
    %  incomplete, resulting in high thickness underestimation! 
    
      % estimation of pial surface
      th2 = 0.1; % GM/CSF border
      cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .2 ' ...
                     'avg -0.05 0.05 .1 .1 5 0 "%g" "%g" n 0 0 0 %d %0.04f 0.0 1'], ...
                     Vpp.fname,Pcentral,Ppial,th2,th2,moveth(0.3,opt.surfaccuracy / 10));
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

      % if deformation stopped earlier call CAT_Central2Pial to get closer to the pial surface
      if ~isempty(strfind(RS,'Stopped after'))
        stime = cat_io_cmd('  Estimation of pial surface','g5','',opt.verb,stime);
        cmd = sprintf('CAT_Central2Pial -check_intersect "%s" "%s" "%s" 0.3',Pcentral,Ppbt,Ppial);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

        cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .2 ' ...
                     'avg -0.05 0.05 .1 .1 5 0 "%g" "%g" n 0 0 0 %d %0.04f 0.0 1'], ...
                     Vpp.fname,Ppial,Ppial,th2,th2,moveth(0.3,opt.surfaccuracy / 10));
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
      end

      stime = cat_io_cmd('  Correction of pial surface in highly folded areas','g5','',opt.verb,stime);
      cmd = sprintf('CAT_Central2Pial -equivolume -weight 0.5 "%s" "%s" "%s" 0',Ppial,Ppbt,Ppial);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

      % estimation of white matter surface
      th2 = 0.9; % GM/WM border
      cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .2 ' ...
                     'avg -0.05 0.05 .1 .1 5 0 "%g" "%g" n 0 0 0 %d %0.04f 0.0 1'], ...
                     Vpp.fname,Pcentral,Pwhite,th2,th2,moveth(0.3,opt.surfaccuracy / 10));
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

      % if deformation stopped earlier call CAT_Central2Pial to get closer to the white matter surface
      if ~isempty(strfind(RS,'Stopped after')) 
        stime = cat_io_cmd('  Estimation of white matter surface','g5','',opt.verb,stime);
        cmd = sprintf(['CAT_Central2Pial -check_intersect "%s" "%s" "%s" -0.3'], ...
                       Pcentral,Ppbt,Pwhite);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

        cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none 0 1 -1 .2 ' ...
                     'avg -0.05 0.05 .1 .1 5 0 "%g" "%g" n 0 0 0 %d %0.04f 0.0 1'], ...
                     Vpp.fname,Pwhite,Pwhite,th2,th2,moveth(0.3,opt.surfaccuracy / 10));
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
      end

      stime = cat_io_cmd('  Correction of white matter surface in highly folded areas','g5','',opt.verb,stime);
      cmd = sprintf('CAT_Central2Pial -equivolume -weight 0.3 "%s" "%s" "%s" 0',Pwhite,Ppbt,Pwhite);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

      % update central surface as average between white and pial surface
      cmd = sprintf('CAT_AverageSurfaces -avg "%s" "%s" "%s"',Pcentral,Pwhite,Ppial);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

      updateThickness = 0;
      if updateThickness % more testing necessary to also correct thickness
        % correction of cortical thickness
        cmd = sprintf('CAT_Hausdorff  "%s" "%s" "%s"',Pwhite,Ppial,Ppbt);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
        saveiconame = 'new';
      else
        saveiconame = 'pbt';
      end
      facevertexcdata = cat_io_FreeSurfer('read_surf_data',Ppbt);
     
      % evaluate and save results
      CS = loadSurf(Pcentral); CS1.vertices = CS.vertices; CS1.faces = CS.faces; CS1.vmat = vmat; CS1.mati = mati; 
      fprintf('%5.0fs',etime(clock,stime)); stime = []; 
      if opt.write_debugsurfs
        cat_surf_fun('saveico',CS1,facevertexcdata,Pcentral,...
          sprintf('debug_3_extract_pial_white_pbtres_%s_pbtres%0.2fmm_vdist%0.2fmm',saveiconame,opt.interpVold,opt.vdist));  
      end
      res.(opt.surf{si}).createCS_3_extract_pial_white = cat_surf_fun('evalCS' ,CS1,facevertexcdata,Ym,Yppt,Pcentral,opt.verb-2);
      res.(opt.surf{si}).createCS_final                = res.(opt.surf{si}).extract_pial_white; 
      clear CS1 
    end

    
    % Further error measures: 
    % - Map local noise pattern?  
    %   >> Detection of movement artifacts?
    % - Map local bias field?     
    %   >> This finally codes the heat noise of the original image >> input by cat_main 
    %   >> Relevant for local interpretation of results, e.g. low signal (high noise) ares are less reliable. 
    % ... both aspects are more general QA topics and not so relevant here
    %
    % - Map collision error?      
    %   >> output of collision detection? - not independent
    %   >> only relevant for me? 
    % - Map distance error? 
    %   >> To what? RAW surface? This would include topology defect errors.
    %      Possibly mask by defect map? And then? 
    %   >> estimate distance of each GM voxels to its closes surface point
    %      . describes the number of voxels out of thickness range (missing GM voxels - but also artefacts) 
    %      estimate distance of each ~GM brain voxels within the local thickness
    %      . describes the number of voxels the are in the cortical surface ribbon but should not 
    %      ... the biggest problem are here the artifacts ...
    %      ... use defect map to avoid counting in some region 
% map defects to the final surface 
    if 0 %opt.verb > 2 
% #################### 
% Torbens error ...
% this does not work because the RAW surface will be deleted and the individual central has another structure       
% you have to use the Delaunay mapping to transfer the values ...
      %cmd = sprintf('CAT_MarkDefects -binary "%s" "%s" "%s"',Praw,Psphere0,Pdefects0); 
      %[ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb);
      
      cmd = sprintf('CAT_AddValuesToSurf "%s" "%s" "%s"',Praw,Pdefects0,Pdefects);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
      
      % read data 
      %varargout = cat_surf_cdatamapping(S1,S2,cdata,opt) 
      % map data
    end
    delete(Pdefects0);  

    
    
    % Test without surface registration - just a shortcut for manual tests! 
    if 0
      cat_io_cmd('  ','g5','',opt.verb,stime);  
      S.(opt.surf{si}) = struct('faces',CS.faces,'vertices',CS.vertices,'vmat',vmat,'vmati',vmati,'mati',mati,'th1',facevertexcdata); 
      if si==numel(opt.surf) && si == 1
        cat_io_cmd('  ','g5','',opt.verb,cstime);
        sprintf('%5ds\n',round(etime(clock,cstime)));
      end
      continue
    end
    
    
    
    
    %% spherical surface mapping and registration of the final corrected surface
    %  use more iterations for higher surfaces 
    stime = cat_io_cmd('  Spherical mapping with areal smoothing','g5','',opt.verb,stime); 
    cmd = sprintf('CAT_Surf2Sphere "%s" "%s" %d',Pcentral,Psphere,5 + round( size(CS.faces,1) / 100000 ));
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
    
    % spherical registration to fsaverage template
    stime = cat_io_cmd('  Spherical registration','g5','',opt.verb,stime);
    if opt.vdist>2 % low quality 
      cmd = sprintf(['CAT_WarpSurf -i "%s" -is "%s" -t "%s" -ts "%s" -ws "%s" ' ...
        '-size 256 128 -loop 1 -steps 1 -runs 2'],...
        Pcentral,Psphere,Pfsavg,Pfsavgsph,Pspherereg);
    else
      cmd = sprintf('CAT_WarpSurf -steps 2 -avg -i "%s" -is "%s" -t "%s" -ts "%s" -ws "%s"', ...
        Pcentral,Psphere,Pfsavg,Pfsavgsph,Pspherereg);
    end
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
    
    % setting of thickness values to zero for masked area (by the inverse transformation)
    % does not work properly for all data and is moreover part of cat_surf_resample ...

    
    % display some evaluation 
    fprintf('%5.0fs\n',etime(clock,stime)); 
    fprintf('  Surface Euler number:                  %d\n',EC0);
    fprintf('  Overall size of topology defects:      %0.2f%% (~%0.2f cm%s)',...
      defect_size0,defect_area0,char(178));

    
    % evaluate and save results
    if opt.write_debugsurfs
      % filenames for resmapling
      Presamp   = fullfile(pp,surffolder,sprintf('%s.tmp.resampled.%s'    ,opt.surf{si},ff));  
      Ppbtr     = fullfile(pp,surffolder,sprintf('%s.pbt.resampled.%s'    ,opt.surf{si},ff));  
      Ppbtr_gii = [Ppbtr '.gii'];
      
      % resample values using warped sphere 
      cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s" "%s" "%s"',Pcentral,Pspherereg,Pfsavgsph,Presamp,Ppbt,Ppbtr);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
      
      if 0 
        % resample surface using warped sphere with better surface quality (using Spherical harmonics)
        % ###
        % deactivated because the resampling of the surface alone leads to displacements of the textures (RD20190927)!
        % ###
        cmd = sprintf('CAT_ResampleSphericalSurfSPH -n 327680 "%s" "%s" "%s"',Pcentral,Pspherereg,Presamp);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);

        % resample surface according to freesurfer sphere
        cmd = sprintf('CAT_ResampleSurf "%s" NULL "%s" "%s"',Presamp,Pfsavgsph,Presamp);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3); 
      end
      
      % add values to resampled surf and save as gifti
      cmd = sprintf('CAT_AddValuesToSurf "%s" "%s" "%s"',Presamp,Ppbtr,Ppbtr_gii); 
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3); 
      if exist(Ppbtr,'file'), delete(Ppbtr); end

      % remove path from metadata to allow that files can be moved (pathname is fixed in metadata) 
      [pp2,ff2,ex2] = spm_fileparts(Ppbtr_gii); 
      g = gifti(Ppbtr_gii);
      g.private.metadata = struct('name','SurfaceID','value',[ff2 ex2]);
      save(g, Ppbtr_gii, 'Base64Binary');
      
      % intensity based evaluation
      CSr = loadSurf(Ppbtr_gii);
      CSr = struct('vertices',CSr.vertices,'faces',CSr.faces,'cdata',CSr.cdata,'vmat',vmat,'mati',mati); 
      cat_surf_fun('saveico',CSr,CSr.cdata,Pcentralr,sprintf('createCS_4_resampled_pbtres%0.2fmm_vdist%0.2fmm',opt.interpVold,opt.vdist)); 
      res.(opt.surf{si}).createCS_resampled = cat_surf_fun('evalCS',CSr,CSr.cdata,Ym,Yppt,Pcentralr);
      clear CSr CS1
    end
    if ~isfield( res.(opt.surf{si}),'createCS_final')
      res.(opt.surf{si}).createCS_final = cat_surf_fun('evalCS',loadSurf(Pcentral),cat_io_FreeSurfer('read_surf_data',Ppbt),Ym,Yppt,Pcentral);
    else 
      fprintf('\n'); 
    end
    
    
    %% average final values
    FNres = fieldnames( res.(opt.surf{si}).createCS_final );
    for fnr = 1:numel(FNres)
      if ~isfield(res,'final') || ~isfield(res.final,FNres{fnr})
        res.final.(FNres{fnr}) = res.(opt.surf{si}).createCS_final.(FNres{fnr}) / numel(opt.surf);
      else
        res.final.(FNres{fnr}) = res.final.(FNres{fnr}) + res.(opt.surf{si}).createCS_final.(FNres{fnr}) / numel(opt.surf);
      end
    end
    if opt.write_debugsurfs
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
    
    
    
    % estimate Freesurfer thickness measure Tfs using mean(Tnear1,Tnear2)
    if opt.thick_measure == 1
      % not ready yet
      if 0
%      if opt.extract_pial_white && ~opt.fast % use white and pial surfaces
        cmd = sprintf('CAT_SurfDistance -mean "%s" "%s" "%s"',Pwhite,Ppial,Pthick);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
      else % use central surface and thickness
        cmd = sprintf('CAT_SurfDistance -mean -thickness "%s" "%s" "%s"',Ppbt,Pcentral,Pthick);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
      end
      
      % apply upper thickness limit
			facevertexcdata = cat_io_FreeSurfer('read_surf_data',Pthick);  
			facevertexcdata(facevertexcdata > opt.thick_limit) = opt.thick_limit;
			cat_io_FreeSurfer('write_surf_data',Pthick,facevertexcdata);  
    else % otherwise simply copy ?h.pbt.* to ?h.thickness.*
      copyfile(Ppbt,Pthick);
    end
    
    
    
    
    %% WM and CSF thickness
    %  Will hopefully be improved in future and may become part of 
    %  cat_surf_parameters and removed here (RD201909)
    if opt.WMT > 1
      % map WM and CSF width data (corrected by thickness)
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
      try %#ok<TRYNC>
        cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s"',Pcentral,Pgwwg,3,Pgwwg);
        [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.verb-3);
      end
      
      facevertexcdata3 = isocolors2(Ycd,CS.vertices); 
      facevertexcdata3 = max(eps,facevertexcdata3 - facevertexcdata/2); 
      cat_io_FreeSurfer('write_surf_data',Psw,facevertexcdata3);
    end
    
    
    % create output structure
    S.(opt.surf{si}) = struct('faces',CS.faces,'vertices',CS.vertices,'vmat',vmat,...
        'vmati',vmati,'mati',mati,'th1',facevertexcdata);
    if opt.WMT > 1
      S.(opt.surf{si}) = setfield(S.(opt.surf{si}),'th2',facevertexcdata2);
      S.(opt.surf{si}) = setfield(S.(opt.surf{si}),'th3',facevertexcdata3);
    end
    clear Yth1i

    
    % we have to delete the original faces, because they have a different number of vertices after
    % CAT_FixTopology!
    if exist(Praw      ,'file'), delete(Praw); end
    if exist(Psphere0  ,'file'), delete(Psphere0); end
    if exist(Vpp.fname ,'file'), delete(Vpp.fname); end
    if exist(Vpp1.fname,'file'), delete(Vpp1.fname); end
    if opt.verb > 2 && exist(Pdefects0,'file'), delete(Pdefects0); end
    if opt.extract_pial_white && exist(Vyp0s.fname,'file'), delete(Vyp0s.fname); end
    if exist(Vmcs.fname,'file'), delete(Vmcs.fname); end
    clear CS

    % create white and central surfaces
    cat_surf_fun('white',Pcentral);
    cat_surf_fun('pial',Pcentral);
    
    % processing time per side for manaual tests
    if si==numel(opt.surf) && si == 1
      cat_io_cmd('  ','g5','',opt.verb,cstime);
      sprintf('%5ds\n',round(etime(clock,cstime)));
    end
  end  
  
  % calculate mean EC and defect size for all surfaces
  EC          = round(EC / numel(opt.surf));
  defect_size = round(defect_size / numel(opt.surf));

  if opt.verb
    for si=1:numel(Psurf)
      fprintf('Display thickness: %s\n',spm_file(Psurf(si).Pthick,'link','cat_surf_display(''%s'')'));
    end
  end
end

%=======================================================================
function saveSurf(CS,P)
  global vmat mati
  
  CS.vertices = (vmat*[CS.vertices' ; ones(1,size(CS.vertices,1))])'; 
  if mati(7)<0, CS.faces = [CS.faces(:,1) CS.faces(:,3) CS.faces(:,2)]; end
  save(gifti(struct('faces',CS.faces,'vertices',CS.vertices)),P,'Base64Binary');
end

%=======================================================================
function CS = loadSurf(P)
  global vmati mati
  CS = gifti(P);
  warning off MATLAB:subscripting:noSubscriptsSpecified
  CS.vertices = (vmati*[CS.vertices' ; ones(1,size(CS.vertices,1))])';
  if mati(7)<0, CS.faces = [CS.faces(:,1) CS.faces(:,3) CS.faces(:,2)]; end
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
% calculates an interpolated value of a vertex in R  
% We have to calculate everything with double, thus larger images will 
% cause memory issues.
% ______________________________________________________________________
  
  if isempty(V), return; end
  if ndims(R)~=3,  error('MATLAB:isocolor2:dimsR','Only 2 or 3 dimensional input of R.'); end
  if ~exist('opt','var'), opt=struct(); end
  
  def.interp = 'linear';
  opt = cat_io_checkinopt(opt,def);
  
  if  isa(R,'double'), R = single(R); end
  if ~isa(V,'double'), V = double(V); VD=0; else, VD=1; end
  
  nV   = size(V,1);
  ndim = size(V,2);
  
  switch opt.interp
    case 'nearest'
      V = max(1,min(round(V),repmat(size(R),nV,1))); 
      V = R(sub2ind(size(R),V(:,2),V(:,1),V(:,3)));
    case 'linear'
      nb  = repmat(shiftdim(double([0 0 0;0 0 1;0 1 0;0 1 1;1 0 0;1 0 1;1 1 0;1 1 1]'),-1),nV,1);  
      enb = repmat(shiftdim((ones(8,1,'double')*[size(R,2),size(R,1),size(R,3)])',-1),nV,1);  

      % calculate the weight of a neigbor (volume of the other corner) and
      w8b = reshape(repmat(V,1,2^ndim),[nV,ndim,2^ndim]); clear V;
      % if the streamline is near the boundary of the image you could be out of range if you add 1 
      n8b = min(floor(w8b) + nb,enb); clear enb
      n8b = max(n8b,1);
      w8b = flipdim(prod(abs(n8b - w8b),2),3);        

      % multiply this with the intensity value of R
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
% the amount of WM fibers that this region will add to the WM depth. 
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
              

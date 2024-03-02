function [Yth,S,Psurf,res] = cat_surf_createCS3(V,V0,Ym,Ya,YMF,Ytemplate,Yb0,opt,job)
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
% $Id: 2558 2024-02-28 $ 

  %#ok<*AGROW,*STREMP,*ASGLU,*SFLD,*STFLD>

  % Turn off gifti data warning in gifti/subsref (line 45)
  %   Warning: A value of class "int32" was indexed with no subscripts specified. 
  %            Currently the result of this operation is the indexed value itself, 
  %            but in a future release, it will be an error. 
  warning('off','MATLAB:subscripting:noSubscriptsSpecified');
  cstime = clock; 
  
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
  def.pbtmethod           = 'CAT_VolThicknessPbt';                      % projection-based thickness (PBT) estimation ('pbt2x' (with minimum setting), 'pbt2', or 'pbtsimple')
  def.sharpenCB           = 1;                                % sharpening function for the cerebellum (in development, RD2017-2019)
  def.thick_measure       = 1;                                % 0-PBT; 1-Tfs (Freesurfer method using mean(TnearIS,TnearOS))
  def.thick_limit         = 5;                                % 5mm upper limit for thickness (same limit as used in Freesurfer)
  def.surf_measures       = 1;                                % 0 - none, 1 - only thickness, 2 - expert maps (myelin,defects), 3 - developer (WMT,CSFT, ...),
                                                              % 4 - debug output, 5 - debug extended (more substeps and mex output)

  def.fsavgDir            = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces'); 
  def.add_parahipp        = cat_get_defaults('extopts.add_parahipp');
  def.scale_cortex        = cat_get_defaults('extopts.scale_cortex');
  def.close_parahipp      = cat_get_defaults('extopts.close_parahipp');
  def.outputpp.native     = 0;  % output of Ypp map for cortical orientation in EEG/MEG 
  def.outputpp.warped     = 0;
  def.outputpp.dartel     = 0;
    
  % options that rely on other options
  opt                     = cat_io_updateStruct(def,opt); clear def; 
  opt.vol                 = any(~cellfun('isempty',strfind(opt.surf,'v')));   % only volume-based thickness estimation  
  opt.surf                = cat_io_strrep(opt.surf,'v','');                   % after definition of the 'vol' varialbe we simplify 'surf'
  opt.interpV             = max(0.1,min([opt.interpV,1.5]));                  % general limitation of the PBT resolution
  
  % for testing only  
  skip_registration = debug & 1;
  
  % we should test this more detailed!
%  opt.add_parahipp       = 0;
%  opt.scale_cortex       = 1.0;
%  fprintf('Set add_parahipp to %d and scale_cortex to %g. Not yet fully tested!\n',opt.add_parahipp,opt.scale_cortex)

  % isosurface threshold
  th_initial = 0.49; % we have to slightly lower the threshold if we use post-smoothing of 2mm with CAT_VolMarchingCubes
  
  % too slow currently
  create_white_pial = opt.SRP==1;
  
  % Another parameter to control runtime is the accuracy of the surface
  % deformation. As far as we primary adapt the mesh resolution above, it 
  % is useful to use sqrt values rather than linear values to avoid square
  % processing times for higher quality levels. Otherwise, we can simple 
  % avoid changes here

  if opt.surf_measures > 4, opt.verb = 3; end
  
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
  Yms = Ym + 0; cat_sanlm(Yms,3,1);
  mf  = min(1,max(0,3-2*mean(vx_vol,2))); 
  Ym  = mf * Yms  +  (1-mf) * Ym;
  clear Yms;
   
  % filling
  Ymf  = max(Ym,min(1,YMF & ~NS(Ya,opt.LAB.HC) & ~( cat_vol_morph( NS(Ya,opt.LAB.HC),'d',2) & NS(Ya,opt.LAB.TH) ))); 
  Ymfs = cat_vol_smooth3X(Ymf,1); 
  Ytmp = cat_vol_morph(YMF,'d',3) & Ymfs>2.3/3;
  Ymf(Ytmp) = max(min(Ym(Ytmp),0),Ymfs(Ytmp)); clear Ytmp Ymfs; 
  Ymf = Ymf * 3;
  
  % removing fine WM structures in the hippocampus area to reduce topological and geometrical defects (added RD20190912)
  % use erode to reduce probability of cutting other gyri
  HCmask = cat_vol_morph( NS(Ya,opt.LAB.HC) , 'de', 1.5, vx_vol); 
  Ymf( HCmask ) =  min(2,Ymf( HCmask )); clear HCmask; 

  % surface output and evaluation parameter 
  Psurf = struct(); 
  res   = struct('lh',struct(),'rh',struct()); 
  % initialize WM/CSF thickness/width/depth maps
  Yth   = zeros(size(Ymf),'single'); 
  Ypp   = -ones(size(Ymf),'single');   
  
  
  %% reduction of artifact, blood vessel, and meninges next to the cortex in SPM segmentations 
  %  (are often visible as very thin structures that were added to the WM 
  %  or removed from the brain)
  if ~opt.SPM
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
    
  
  % cleanup and smoothing of the hippocampus amygdala to remove high
  % frequency structures that we cannot preprocess yet
  Ymf = hippocampus_amygdala_cleanup(Ymf,Ya,vx_vol,opt.close_parahipp,1); % last var = doit  
  
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
      case {'lh'},  Ymfs = Ymf .* (Ya>0) .* ~(NS(Ya,opt.LAB.CB) | NS(Ya,opt.LAB.BV) | NS(Ya,opt.LAB.ON) | NS(Ya,opt.LAB.MB)) .* (mod(Ya,2)==1); Yside = mod(Ya,2)==1; 
      case {'rh'},  Ymfs = Ymf .* (Ya>0) .* ~(NS(Ya,opt.LAB.CB) | NS(Ya,opt.LAB.BV) | NS(Ya,opt.LAB.ON) | NS(Ya,opt.LAB.MB)) .* (mod(Ya,2)==0); Yside = mod(Ya,2)==0;  
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
      mask_parahipp = NS(Ya,opt.LAB.PH) | NS(Ya,opt.LAB.HC) | NS(Ya,opt.LAB.VT); 
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
    
    if opt.close_parahipp && ~iscerebellum
      %% RD202107:  Additional close of parahippocampus for the WM.
      %             Dynamic closing size.
      tmp = Ymfs>2.5 | (Ymfs/2 .* VT)>1.0; 
      tmp(smooth3(tmp)<0.3) = 0; 
      tmp = cat_vol_morph(tmp,'lab'); % remove small dots 
      tmp = cat_vol_morph(tmp,'close',round(1/opt.interpV)); % close holes
      Ymfs(mask_parahipp) = max(Ymfs(mask_parahipp),2.6 * tmp(mask_parahipp) .* (1-HC(mask_parahipp))); 
      if ~debug, clear tmp; end 
    end
    


    stime = cat_io_cmd(sprintf('  Thickness estimation (%0.2f mm%s)',opt.interpV,native2unicode(179, 'latin1'))); stimet = stime;
    if strcmp(opt.pbtmethod,'pbtsimple') 
      [Yth1i,Yppi] = cat_vol_pbtsimple(Ymfs,opt.interpV); 
    elseif strcmp(opt.pbtmethod,'CAT_VolThicknessPbt')
      Vtmp = resI.hdrN;
      Vtmp.fname = fullfile(pp0_surffolder,sprintf('%s.Ymfs.%s.nii',opt.surf{si},ff));
      gmt_name = fullfile(pp0_surffolder,sprintf('%s.gmt.%s.nii',opt.surf{si},ff));
      ppm_name = fullfile(pp0_surffolder,sprintf('%s.ppm.%s.nii',opt.surf{si},ff));
      Vtmp.pinfo = Vtmp.pinfo(:,1);
      spm_write_vol(Vtmp, Ymfs);
      cmd = sprintf('CAT_VolThicknessPbt "%s" "%s" "%s"',Vtmp.fname,gmt_name, ppm_name);
      fprintf('\n%s\n\n',cmd);
      cat_system(cmd,opt.verb-3);      
      Yth1i = single(spm_read_vols(spm_vol(gmt_name)));
      Yppi = single(spm_read_vols(spm_vol(ppm_name)));
    else 
      [Yth1i,Yppi,Ymfs] = cat_vol_pbt(Ymfs,struct('method',opt.pbtmethod,'resV',opt.interpV,'vmat',...
        V.mat(1:3,:)*[0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1],'pbtlas',opt.pbtlas)); % avoid underestimated thickness in gyri
    end  

    % back to internal resolution (at least for thickness)
    Yth1i(Yth1i>10) = 0; Yppi(isnan(Yppi)) = -1;                           % general thickness limit  
    [D,I] = cat_vbdist(Yth1i,Ysidei); Yth1i = Yth1i(I); clear D I;         % add further values around the cortex
    Yth1t = cat_vol_resize(Yth1i,'deinterp',resI);                         % back to original resolution
    Yth1t = cat_vol_resize(Yth1t+1,'dereduceBrain',BB)-1;                  % adding background
    Yth  = max(Yth,Yth1t .* Yside);                                        % save on main image
    if ~debug, clear Yth1t Ysidei; end
    %if debug, fprintf('\n%26.22f | %26.22f | %26.22f\n',std(Ymfs(:)),std(Yth1i(:)),std(Yppi(:))); end

    % use median and subtle smoothing to minimize noise in
    Yth1i  = cat_vol_median3(Yth1i);                    
    spm_smooth(Yth1i,Yth1i,[1 1 1] / opt.interpV);
    
    %% PBT estimation of the gyrus and sulcus width 
    if opt.surf_measures > 3 
      Ywd = zeros(size(Ymf),'single'); 
      Ycd = zeros(size(Ymf),'single'); 
      % CG202105 call of subfunction not yet working because of missing variables
%      [Ywdt,Ycdt,stime] = cat_surf_createCS2wdcd(Ya,Ym,Ywd,Ycd,stime); 
      NS = @(Ys,s) Ys==s | Ys==s+1; 
     
      stime = cat_io_cmd('  WM depth estimation','g5','',opt.verb-1,stime); 
      
      [Yar,Ymr] = cat_vol_resize({Ya,Ym},'reduceBrain',vx_vol,BB.BB);       % removing background
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
        case {'cb'}
          Ymr = Ymr .* (Yar>0) .* NS(Yar,3);
          Ynw = true(size(Ymr)); 
      end 
      clear Yar; 
    
      %
      Yppis = Yppi .* (1-Ynw) + max(0,min(1,Ymr*3-2)) .* Ynw; clear Ynw;              % adding real WM map 
      Ywdt  = cat_vol_eidist(1-Yppis,ones(size(Yppis),'single'));                     % estimate distance map to central/WM surface
      Ywdt  = cat_vol_pbtp(max(2,4-Ymfs),Ywdt,inf(size(Ywdt),'single'))*opt.interpV;
      [D,I] = cat_vbdist(single(Ywdt>0.01),Yppis>0); Ywdt = Ywdt(I); clear D I Yppis; % add further values around the cortex
      Ywdt  = cat_vol_median3(Ywdt,Ywdt>0.01,Ywdt>0.01);                    
      Ywdt  = cat_vol_localstat(Ywdt,Ywdt>0.1,1,1);                                   % smoothing
      Ywdt  = cat_vol_resize(Ywdt,'deinterp',resI);                                   % back to original resolution
      Ywdt  = cat_vol_resize(Ywdt,'dereduceBrain',BB);                                % adding background
      Ywdt  = max(Ywd,Ywdt); 
      clear Ywd;
    
      % sulcus width / CSF depth
      %  for the CSF depth we cannot use the original data, because of
      %  sulcal blurring, but we got the PP map at half distance and
      %  correct later for half thickness
      stime = cat_io_cmd('  CSF depth estimation','g5','',opt.verb-1,stime); 
      YM    = single(smooth3(cat_vol_morph(Ymr<0.1,'o',4))<0.5); YM(YM==0)=nan;       % smooth CSF/background-skull boundary 
      Yppis = Yppi .* ((Ymr+0.25)>Yppi) + min(1,Ymr*3-1) .* ((Ymr+0.25)<=Yppi);       % we want also CSF within the ventricle (for tests)
      Ycdt  = cat_vol_eidist(Yppis,YM);                                               % distance to the central/CSF-GM boundary
      Ycdt  = cat_vol_pbtp(max(2,Ymfs),Ycdt,inf(size(Ycdt),'single'))*opt.interpV; Ycdt(isnan(Ycdt))=0;
      [D,I] = cat_vbdist(single(Ycdt>0),Yppis>0 & Yppis<3); Ycdt = Ycdt(I); clear D I Yppis; % add further values around the cortex
      Ycdt  = cat_vol_median3(Ycdt,Ycdt>0.01,Ycdt>0.01);                              % median filtering
      Ycdt  = cat_vol_localstat(Ycdt,Ycdt>0.1,1,1);                                   % smoothing
      Ycdt  = cat_vol_resize(Ycdt,'deinterp',resI);                                   % back to original resolution
      Ycdt  = cat_vol_resize(Ycdt,'dereduceBrain',BB);                                % adding background
      Ycdt  = max(Ycd,Ycdt); 

      clear Ywd Ycd; 
    end
    if ~useprior, fprintf('%5.0fs\n',etime(clock,stime)); end
    
    
    
    %% Replace isolated voxels and holes in Ypp by its median value
    %  RD20210401: this has now only low or no effects
    if 1
      % indicate isolated holes and replace by median of the neighbors
      Yppi(Yppi<0.35 & ~cat_vol_morph(Yppi<1,'l')) = 1;  % close major wholes in the WM 
      Ymsk = Yppi==-1 & cat_vol_morph(Yppi>0.9,'d',1); % filter small wholes close to the WM
      Yppi = cat_vol_median3(single(Yppi)+1,Ymsk,~Ymsk)-1; 

      % indicate isolated objects and replace by median of the neighbors
      Yppi(Yppi>0.65 & cat_vol_morph(Yppi==-1,'l')) = -1;
      Ymsk = Yppi>0.95 & cat_vol_morph(Yppi<-0.9,'d',1); 
      Yppi = cat_vol_median3(single(Yppi)+1,Ymsk,~Ymsk)-1;
      if ~debug, clear Ymsk; end
    end
    
    if 1
      % RD20210401: but there are some new background dots
      Ymsk = Yppi==-1 & cat_vol_morph(Yppi>-1,'lc');  % close major wholes in the WM 
      Yppi = cat_vol_median3(single(Yppi)+1,Ymsk,~Ymsk)-1; 
    end     
    
    %%  Write Ypp for final deformation
    %  Ytt has the internal resolution (e.g. <1 mm) and is only used temparary 
    %  Write Yppi file with 1 mm resolution for the final deformation, 
    %  because CAT_DeformSurf achieved better results using that resolution
    %  RD20210401: The +1 is important to avoid problems with negative
    %              values and boundary effects with zeros.
    %###################
    % RD20210401: Why are we using only the 1 and not the 0.5 mm resolution?
    %             Was there still some resoluion issue or was the interpolation as smoothing helpfull?  
    Yppt = cat_vol_resize(Yppi + 1,'deinterp',resI);                       % back to original resolution
    Yppt = cat_vol_resize(Yppt,'dereduceBrain',BB) - 1;                    % adding of background
    
    % update Ypp with potentially new hemispheric maps
    Ypp(Yppt>-1) = max(Ypp(Yppt>-1),Yppt(Yppt>-1)); 
    
    % This pp-map combines both hemispheres and is saved in the original image 
    % resolution. This map is not used, but just kept for debugging
    % or testing
    if debug
      Vpp  = cat_io_writenii(V0,Ypp,mrifolder,sprintf('%s.pp',opt.surf{si}) ,...
        'percentage position map - pial=0, central=0.5, white=1','uint8',[0,1/255],...
        min([1 1 2],[1 opt.outputpp.warped opt.outputpp.dartel]),opt.trans);
      Vpp  = Vpp(1); 
    end
    
    % keep largest structure
    Yppt1 = Yppt >= th_initial;
    Yppt1 = Yppt1 - cat_vol_morph(Yppt1, 'lab', 1, vx_vol);
    Yppt(Yppt1 > 0) = -1;

    % fill remaining holes
    Yppt1 = Yppt < th_initial;
    Yppt1 = Yppt1 - cat_vol_morph(Yppt1, 'lab', 1, vx_vol);
    Yppt(Yppt1 > 0) = 1;

    % Map with only one hemisphere that is internally used for surface
    % processing
    Vpp_side = cat_io_writenii(V,Yppt,'',sprintf('%s.ppt',opt.surf{si}) ,...
      'percentage position map - pial=0, central=0.5, white=1','uint8',[0,1/255],[1 0 0]);
    Vpp_side = Vpp_side(1);

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
      %% scaling correction to reduce topology errors in the parahippocampus
      % ### not tested for cerebellum yet (RD201909) ###
      if ~iscerebellum
        scale_cortex = opt.scale_cortex; %/0.7 * 0.9;

        ind0   = find(Yppi<=0);
        Yppisc = scale_cortex * Yppi;

        % smooth mask to have smooth border
        mask_parahipp_smoothed = zeros(size(mask_parahipp));
        spm_smooth(double(mask_parahipp),mask_parahipp_smoothed,[4 4 4] / opt.interpV);
        Yppisc  = Yppisc + opt.add_parahipp / scale_cortex * mask_parahipp_smoothed;
        Yppisc(ind0) = 0; clear ind0; %#ok<FNDSB>

        % optionally apply closing inside mask for parahippocampal gyrus to get rid 
        % of the holes that lead to large cuts in gyri after topology correction
        if opt.close_parahipp && ~iscerebellum
          tmp = cat_vol_morph(Yppisc .* (mask_parahipp | (Ymfs/2 .* VT)>1)  ,'close',round(1/opt.interpV)); 
          Yppisc(mask_parahipp) = tmp(mask_parahipp) .* (1-HC(mask_parahipp)); 
          if ~debug, clear tmp; end 
        end
      else
        scale_cortex = 0.7;
        Yppisc = scale_cortex * Yppi;
        
        %% thickness depending cortical scaling - this seems to work but need further tests (RD201911)
        %  RD202103 ... this changes a lot !
        Yth1i  = cat_vol_localstat(Yth1i,Yth1i>0,2,2);
        Yts    = cat_vol_approx(Yth1i,2);  
        Yts    = 1 + max(-0.5,min(0.5,(Yts - mean(Yth1i(:))) / (2 * mean(Yth1i(:)))  )); 
        Yppisc = max(0.55 .* (Yppi>=1),min(1.5, Yppisc .* Yts )); % factor 1 is Ypp at 0.5 > limit 0.55 and 1.5 
        scale_cortex = scale_cortex * median( Yts(:) );
      end

      clear mask_parahipp_smoothed;

      % we use another scaling and can therefore update this map
      Yppisc = Yppisc .* Ymfs/2; 
      Yppisc(smooth3(Yppisc)<0.3) = 0; 
      Yppisc(smooth3(Yppisc)>0.7) = 1; 
      Yppisc( Yppisc<0.5 & ~cat_vol_morph(Yppisc<0.5,'l')) = 1;  % close major wholes in the WM 
      Yppisc( Yppisc>0.5 & ~cat_vol_morph(Yppisc>0.5,'l')) = 0;  % remove small dots


      %%
      % Marching cubes surface creation and correction for the boundary box 
      % used within the surface creation process.  It is better to use the 
      % full voxel resolution in combination with surface reduction rather   
      % then using lower voxel resolutions! Moreover, it was NOT necessary to 
      % use surface deformation before the surface reduction!
      % #####
      %   I am not sure if the topologoy correction is optimal.
      %   Moreover, a correction should also change the Yppi to avoid self-intersections. 
      %   Maybe a smooth adaptation similar to "mask_parahipp_smoothed" can be used here. 
      %   However, this is quite complex and I miss the time go on ... 
      %   RD201911
      % #####
      if iscerebellum
          
        %% region-growing
        % Ylt = single( cat_vol_morph(Yppi<0.1,'l') +  2 * cat_vol_morph(Yppi>0.9,'l') );
        Ylt = single( cat_vol_morph(Ymfs<1.9,'l') +  2 * cat_vol_morph(Ymfs>2.5,'l') );
        [Ylt,D] = cat_vol_simgrow(Ylt,Ymfs+Yppi,0.05); Ylt(D>10) = 0;  % Yppi is to similiar
        [Ylt,D] = cat_vol_simgrow(Ylt,Ymfs+Yppi,0.10); Ylt(D>10) = 0; 
        [Ylt,D] = cat_vol_simgrow(Ylt,Ymfs+Yppi,0.20); Ylt(D>10) = 0; 
        [Ylt,D] = cat_vol_simgrow(Ylt,Ymfs+Yppi,0.50); Ylt(D>10) = 0; 
        Ylts    = smooth3(Ylt);
        

        %% relative position map between core and hull
        %  to control opening and closing operations to reduce topology defects
        %  the topology correction is critical part where we need a quite
        %  robust inner core to avoid superlarge wholes 

        % the hull is quite simple in general 
        Ycbh  = cat_vol_morph(Yppi>0.1 & Ymfs>1.9,'lc',4); 
        % the core is more complicated because we have assure that something
        % is there (at least the highest 20% of the hull distance)
        Ycbhd = cat_vbdist( single( ~Ycbh ) ); 
        Ycbhd = Ycbhd/max(Ycbhd(Ycbhd<1000)); 
        Ycbc  = cat_vol_smooth3X(Yppi + Ymfs/3 + (Ycbhd*0.8) ,2)>2; 
        Ycbc  = cat_vol_smooth3X( cat_vol_morph( cat_vol_morph( cat_vol_morph( Ycbc ,'o',1) ,'l',[2 0.2]), 'lc',2),2) >0.5;
        %% moreover we can use the Laplace filter for subtile openen and closing
        if debug, tic; end
        Yltw  = (Yppi>0.6 | cat_vol_morph(Yppi>0.9,'dc',1) | Ycbc)/2 + Ycbc/2; 
        Yltw  = cat_vol_laplace3R(single(Yltw),Yltw>0 & Yltw<1,0.01);
        Yltw  = max(Yltw,Ycbc | cat_vol_morph( cat_vol_morph( Yltw>0.25 , 'do' , 1.5 ) , 'l' , [2 0.2])); 
        Yltw  = cat_vol_laplace3R(single(Yltw),Yltw>0 & Yltw<1,0.001);
        if debug, toc; end
        %%
        if debug, tic; end
        Yltc  = Ycbh/2 + (Yltw>0.002)/2;
        Yltc  = cat_vol_laplace3R(single(Yltc),Yltc>0 & Yltc<1,0.01);
        Yltw(Yltc>0.99) = 1; 
        if debug, toc; end
        %%
        Ycbc  = cat_vol_smooth3X( cat_vol_smooth3X( Yltw , 4)> 0.6 & Yltw>0.5 ,2)>.5 ;
        Yppisc = max(0.55 .* (Yppi>=1),min(1.5, Yppisc .* Yts )); % factor 1 is Ypp at 0.5 > limit 0.55 and 1.5 
        if debug, clear Yts; end

        %Ycbc  = cat_vol_morph( Ycbc , 'd') & Yppi>0.5;
        %% distance mapping
        Ycbhd = cat_vbdist( single( ~Ycbh ) , ~Ycbc ); 
        Ycbcd = cat_vbdist( single(  Ycbc ) ,  Ycbh ); 

        Ycbpp = min(Ycbh,max(Ycbc,Ycbhd ./ max(eps,Ycbhd+Ycbcd))); 
        Ycbpp = cat_vol_localstat(Ycbpp,Ycbpp>0,1,1); 
        Ycbcd = Ycbcd ./ mean(Ycbcd(Ycbcd(:)>0 & Ycbcd(:)<1000));
        Ycbhd = Ycbhd ./ mean(Ycbhd(Ycbhd(:)>0 & Ycbhd(:)<1000));
        if ~debug, clear Ycbc Ycbh; end

        %% 
        Ycbth  = @(lth,pth,mth,cth) max(Ylt,Ylts)>(1.5*lth) & (Yppi - ( 0.5 - Ycbpp )/2 )>pth & Ymfs>mth & ...
                   ( (cth>=0).*(Ycbpp>cth)  | (cth<0).*(Ycbpp<-cth) ); 
        Ycbth2 = @(lth,pth,mth,cth) max(Ylt,Ylts)>(1.5*lth) & (Yppi - ( 0.5 - Ycbpp )/2 )>pth & Ymfs>mth & ...
                   ( (cth>=0).*(Ycbcd<cth)  | (cth<0).*(Ycbcd>-cth) ); 
        if ~debug, clear Ylt Ycbcd; end

        %%
        Ycbm   = Yppi .* min( 0.4 + Ycbhd/5, 0.5 + Ycbpp*0.60); % in sum larger than one to have save core structure
        Ycbms  = smooth3(Ycbm);
        % opening by smoothing in outer regions
        Ycbppt = max(0,min(1,Ycbpp * 5)); 
        Ycbm   = Ycbms .* (1-Ycbppt) + Ycbm .* Ycbppt; 
        % closing by smoothing in inner regions
        Ycbppt = max(0,min(1,(1 - Ycbpp ) * 5));
        Ycbm   = Ycbms .* (1-Ycbppt) + Ycbm .* Ycbppt; 
        clear Ycbms; 
        clear Ycbhd Ycbpp

        Ycbm = cat_vol_resize(Ycbm,'deinterp',resI);                       % back to original resolution
        Ycbm = cat_vol_resize(Ycbm,'dereduceBrain',BB);                    % adding of background

        Ycbm1 = Ycbm >= th_initial;
        Ycbm1 = Ycbm1 - cat_vol_morph(Ycbm1, 'lab', 1, vx_vol);
        Ycbm(Ycbm1 > 0) = -1;
    
        % fill remaining holes
        Ycbm1 = Ycbm < th_initial;
        Ycbm1 = Ycbm1 - cat_vol_morph(Ycbm1, 'lab', 1, vx_vol);
        Ycbm(Ycbm1 > 0) = 1;

        spm_write_vol(Vpp_side,Ycbm);
         
        %cmd = sprintf('CAT_MarchingCubesGenus0 -fwhm "3" -thresh "%g" "%s" "%s"',th_initial,Vpp_side.fname,Pcentral);
        % -pre-fwhm = -1: masked smoothing with 1mm
        % -post-fwhm = 2: psot smoothing with 2mm (use lower threshold of 0.49)
        cmd = sprintf('CAT_VolMarchingCubes -pre-fwhm "-1" -post-fwhm "2" -thresh "%g" "%s" "%s"',th_initial,Vpp_side.fname,Pcentral);
        fprintf('\n%s\n\n',cmd);
        cat_system(cmd,opt.verb-3);      
      
      else
        % Main initial surface creation using cat_vol_genus0
        % cat_vol_genus0 uses a "simple" marching cube without use of isovalues
        % that is used in the MATLAB isosurface function. Our test showed that
        % the surface deformation allows the same or better accuracy and also
        % that the meshes of cat_vol_genus0 are more regular and also allow 
        % voxel-based topology optimization.  
        % We can set dist to 0.9 for the cortex because using adaptive
        % thresholds this was found to be the optimal dist threshold

        %cmd = sprintf('CAT_MarchingCubesGenus0 -fwhm "3" -thresh "%g" "%s" "%s"',th_initial,Vpp_side.fname,Pcentral);
        % -pre-fwhm = -1: masked smoothing with 1mm
        % -post-fwhm = 2: psot smoothing with 2mm (use lower threshold of 0.49)
        cmd = sprintf('CAT_VolMarchingCubes -pre-fwhm "-1" -post-fwhm "2" -thresh "%g" "%s" "%s"',th_initial,Vpp_side.fname,Pcentral);
        cat_system(cmd,opt.verb-3);      
        
      end
      %
            
      CS = loadSurf(Pcentral);
      EC0 = size(CS.vertices,1) + size(CS.faces,1) - size(spm_mesh_edges(CS),1);
      
      if EC0 ~= 2
        warning('cat_surf_createCS3:CAT_MarchingCubesGenus0', ...
           'Extracted surface might have small topology issues (Euler count = %d).\n',EC0); 
      end      

      if ~debug, clear mask_parahipp_smoothed; end
                  
      % evaluate and save results
      if isempty(stime), stime = clock; end
      fprintf('%5.0fs',etime(clock,stime)); stime = []; if 1, fprintf('\n'); end %debug
      %{
      res.(opt.surf{si}).createCS_init = cat_surf_fun('evalCS',CS,cat_surf_fun('isocolors',CS,Yth1i,Smat.matlabIBB_mm),[],Ymfs,Yppi,Pcentral,Smat.matlabIBB_mm,debug);
      if 0 %debug 
        % save surface for further evaluation 
        cat_surf_fun('saveico',CS,cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm),Pcentral,sprintf('createCS_1_init_pbtres%0.2fmm',opt.interpV),Ymfs,Smat.matlabIBB_mm); 
      else
        fprintf('\n'); 
      end
      %}
    end

    % for use of average prior in long. pipeline we have to deform the average mesh to current pp distance map        
    if useprior
      CS = loadSurf(Pcentral);
      correct_mesh = 1;
      cmd = sprintf(['CAT_DeformSurf "%s" none 0 0 0 "%s" "%s" none  0  1  -1  .1 ' ...   
                     'avg  %0.3f  %0.3f .2  .1  5  0 "0.5"  "0.5"  n 0  0  0 %d %g  0.0 0 %d'], ...          
                      Vpp_side.fname,Pcentral,Pcentral,-0.1,0.1,100,0.01,correct_mesh); 
      cat_system(cmd,opt.verb-3);
    end
    
    % get thickness data
    facevertexcdata = cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm); 
    cat_io_FreeSurfer('write_surf_data',Ppbt,facevertexcdata);
    
    % Test without surface registration - just a shortcut for manual tests! 
    if 0 %skip_registration
      cat_io_cmd('  ','g5','',opt.verb,stime);  
      S.(opt.surf{si}) = struct('faces',CS.faces,'vertices',CS.vertices,'th1',facevertexcdata); clear facevertexcdata; 
      if si==numel(opt.surf) && si == 1
        cat_io_cmd('  ','g5','',opt.verb,cstime);
        sprintf('%5ds\n',round(etime(clock,cstime)));
      end
      res.(opt.surf{si}).createCS_final = cat_surf_fun('evalCS',loadSurf(Pcentral),cat_io_FreeSurfer('read_surf_data',Ppbt),facevertexcdata,Ymfs,Yppi,Pcentral,Smat.matlabIBB_mm,debug,cat_get_defaults('extopts.expertgui')>1);
      continue
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
      cmd = sprintf('CAT_WarpSurf -steps 2 -avg -i "%s" -is "%s" -t "%s" -ts "%s" -ws "%s"', ...
        Pcentral,Psphere,Pfsavg,Pfsavgsph,Pspherereg);
      cat_system(cmd,opt.verb-3);
    end  

    % create white and central surfaces
    if create_white_pial
      stime = cat_io_cmd('  Create pial and white surface','g5','',opt.verb,stime); 
      cmd = sprintf('CAT_Central2Pial -check_intersect "%s" "%s" "%s" 0.5',Pcentral,Ppbt,Ppial);
      cat_system(cmd,opt.verb-3);
      cmd = sprintf('CAT_Central2Pial -check_intersect "%s" "%s" "%s" -0.5',Pcentral,Ppbt,Pwhite);
      cat_system(cmd,opt.verb-3);
    elseif opt.SRP>=2
      %% call collision correction
      %  RD202108: Use further iterations if self-intersections are still very high.  
      %            (test data was an high resolution ex-vivo chimp PD image that had still strong SIs after first correction) 
      if opt.SRP==2
        stime = cat_io_cmd('  Reduction of surface collisions:','g5','',opt.verb,stime); 
      else
        stime = cat_io_cmd('  Reduction of surface collisions with optimization:','g5','',opt.verb,stime); 
      end
      CS = loadSurf(Pcentral); 
      SIOs = 100; SIs = 80; maxiter = 1; iter = 0; verblc = debug;  % maxiter =2
      while SIs>5 && SIs<SIOs*0.9 && iter<maxiter
        SIOs = SIs; iter = iter + 1; 
        [CS,facevertexcdata,SIs] = cat_surf_fun('collisionCorrectionPBT',CS,facevertexcdata,Ymfs,Yppi,...
          struct('optimize',iter<2 && opt.SRP>=2,'verb',verblc,'mat',Smat.matlabIBB_mm,'vx_vol',vx_vol)); 
        if verblc, fprintf('\b\b'); end
        if strcmpi(spm_check_version,'octave') && iter == 1
          cat_io_addwarning('cat_surf_createCS3:nofullSRP','Fine correction of surface collisions is not yet available under Octave.',2)
        else % ############### not working any longer
         %[CS,facevertexcdata,SIs] = cat_surf_fun('collisionCorrectionRY' ,CS,facevertexcdata,Ymfs,struct('Pcs',Pcentral,'verb',verblc,'mat',Smat.matlabIBB_mm,'accuracy',1/2^3)); 
        end
      end
      saveSurf(CS,Pcentral);
      cat_io_FreeSurfer('write_surf_data',Pthick,facevertexcdata);  
    end

    
    % write myelination map (Ypp intensity of layer 4)  
    if opt.surf_measures > 1 
      cmd = sprintf('CAT_Central2Pial -equivolume -weight 1 "%s" "%s" "%s" 0',Pcentral,Ppbt,Player4);
      cat_system(cmd,0);
      L4  = gifti(Player4);
      L4v = cat_surf_fun('isocolors',Ymf,L4,Smat.matlabi_mm); clear L4
      cat_io_FreeSurfer('write_surf_data',PintL4,L4v); clear L4v
    end
    
    
    % estimate Freesurfer thickness measure Tfs using mean(Tnear1,Tnear2)
    if opt.thick_measure == 1
      stime = cat_io_cmd('  Estimate final thickness','g5','',opt.verb,stime); 
      cmd = sprintf('CAT_SurfDistance -mean -thickness "%s" "%s" "%s"',Ppbt,Pcentral,Pthick);
      cat_system(cmd,opt.verb-3);
      
      % apply upper thickness limit
      facevertexcdata = cat_io_FreeSurfer('read_surf_data',Pthick);  
      facevertexcdata(facevertexcdata > opt.thick_limit) = opt.thick_limit;
      cat_io_FreeSurfer('write_surf_data',Pthick,facevertexcdata);  
      
      % final surface evaluation 
      %if debug || cat_get_defaults('extopts.expertgui')>1, fprintf('\n'); end
      if debug
        cat_surf_fun('saveico',CS,cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm),Pcentral,'',Ymfs,Smat.matlabIBB_mm); 
      end
      res.(opt.surf{si}).createCS_final = cat_surf_fun('evalCS',loadSurf(Pcentral),cat_io_FreeSurfer('read_surf_data',Ppbt),facevertexcdata,Ymfs,Yppi,Pcentral,Smat.matlabIBB_mm,debug,cat_get_defaults('extopts.expertgui')>1);
    else % otherwise simply copy ?h.pbt.* to ?h.thickness.*
      copyfile(Ppbt,Pthick,'f');
  
      % final surface evaluation 
      % save surface for further evaluation 
      if debug
        cat_surf_fun('saveico',CS,cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm),Pcentral,'',Ymfs,Smat.matlabIBB_mm); 
      end
      res.(opt.surf{si}).createCS_final = cat_surf_fun('evalCS',loadSurf(Pcentral),cat_io_FreeSurfer('read_surf_data',Ppbt),[],Ymfs,Yppi,Pcentral,Smat.matlabIBB_mm,debug,cat_get_defaults('extopts.expertgui')>1);
    
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
    
    
    %% WM and CSF thickness
    %  Will hopefully be improved in future and may become part of 
    %  cat_surf_parameters and removed here (RD201909)
    if opt.surf_measures > 3 
      % map WM and CSF width data (corrected by thickness)
      facevertexcdata2  = cat_surf_fun('isocolors',Ywdt,CS.vertices,Smat.matlabi_mm); 
      facevertexcdata2c = max(eps,facevertexcdata2 - facevertexcdata/2);
      cat_io_FreeSurfer('write_surf_data',Pgwo,facevertexcdata2c); % gyrus width WM only
      facevertexcdata2c = correctWMdepth(CS,facevertexcdata2c,100,0.2);
      cat_io_FreeSurfer('write_surf_data',Pgww,facevertexcdata2c); % gyrus width WM only
      facevertexcdata3c = facevertexcdata2c + facevertexcdata; % );
      cat_io_FreeSurfer('write_surf_data',Pgw,facevertexcdata3c); clear facevertexcdata3c; % gyrus width (WM and GM)
      facevertexcdata4 = estimateWMdepthgradient(CS,facevertexcdata2c); clear facevertexcdata2c; 
      cat_io_FreeSurfer('write_surf_data',Pgwwg,facevertexcdata4); clear facevertexcdata4; % gyrus width WM only > gradient
      
      % smooth resampled values
      try %#ok<TRYNC>
        cmd = sprintf('CAT_SmoothDiffusion -fwhm "%g" -values "%s" "%s" "%s"',3,Pgwwg,Pcentral,Pgwwg);
        cat_system(cmd,opt.verb-3);
      end
      
      facevertexcdata3 = cat_surf_fun('isocolors',Ycdt,CS.vertices,Smat.matlabi_mm); 
      facevertexcdata3 = max(eps,facevertexcdata3 - facevertexcdata/2); 
      cat_io_FreeSurfer('write_surf_data',Psw,facevertexcdata3);
    end
    
    
    % create output structure
    S.(opt.surf{si}) = struct('faces',CS.faces,'vertices',CS.vertices,'th1',facevertexcdata); clear facevertexcdata; 
    if opt.surf_measures > 3
      S.(opt.surf{si}) = setfield(S.(opt.surf{si}),'th2',facevertexcdata2); clear facevertexcdata2; 
      S.(opt.surf{si}) = setfield(S.(opt.surf{si}),'th3',facevertexcdata3); clear facevertexcdata3; 
    end
    clear Yth1i

    % we have to delete the original faces, because they have a different 
    % number of vertices after CAT_FixTopology!
    if exist(Vpp_side.fname ,'file'), delete(Vpp_side.fname); end
    if debug && exist(Vpp.fname ,'file') && ~opt.outputpp.native, delete(Vpp.fname); end
    
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
  
    if debug 
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
function CS1 = loadSurf(P)

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

%=======================================================================
function [cdata,i] = correctWMdepth(CS,cdata,iter,lengthfactor)
% ______________________________________________________________________
% Correct deep WM depth values that does not fit to the local thickness 
% of the local gyri.
% 
% length factor should be between 0.2 and 0.4
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
  while i<iter && pc~=oc 
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


    
    %% RD202107: closing of parahippocampal gyri
    %  --------------------------------------------------------------------
    %  We are mostly interested in closing of holes in deep parahippocampal 
    %  regions. Hence, we will limit our operations by an enlargements of
    %  other structures (like sulcal depth).
    %  Try close parahippocampal gyrus by increasing the intensity of 
    %  GM-WM structures before we filter in the hippocampus.
    %  --------------------------------------------------------------------
    if close_parahipp
      % Limitation by sulcal depth like map.
      dd   = 1; % this parameter shoulb depend on brain size to work in primates too
      Yphn = cat_vol_morph( NS(Ya,LAB.BS) | NS(Ya,LAB.CB) | NS(Ya,LAB.ON) | Ymf==0,'dd',dd * 8, vx_vol) | ... 
             cat_vol_morph( NS(Ya,LAB.HC),'de',3);  % initial definition with extensiion (~ sulcal depth)
      Yphn = cat_vol_smooth3X(Yphn,2)>0.5;          % soften boundaries
      Yphn = cat_vol_morph( Yphn , 'dc' , 4);       % close some regions

      % mask by "obvious" structures           
      Ymsk = Ymf>2.1 & cat_vol_morph( NS(Ya,LAB.HC) |  NS(Ya,LAB.PH), 'dd' , 3 , vx_vol ) & ~Yphn; 
      Ymsk(smooth3(Ymsk)<0.5) = 0; 
      Ymsk = cat_vol_morph( Ymsk , 'dc' , 1.5 , vx_vol ); 

      % we need a wider region close to the hippocampus and parahippocampal gyrus
      Yg   = cat_vol_grad(Ymf);
      Ydiv = cat_vol_div(Ymf);
      Yphi = cat_vol_morph( NS(Ya,LAB.HC), 'dd' , 5 , vx_vol ) & cat_vol_morph( NS(Ya,LAB.PH), 'dd' , 5 , vx_vol ) & ~Yphn; 
      Yphx =  ( Yphi | cat_vol_morph( NS(Ya,LAB.PH), 'dd' , 2 , vx_vol )) & -Ydiv>(2-Ymf-Yg/3) & ~Yphn; clear Yphi; 
      Yphx = cat_vol_morph( Yphx & mod(Ya,2)==0, 'l' ) | cat_vol_morph( Yphx & mod(Ya,2)==1, 'l' );
      Yphx = cat_vol_morph( Yphx , 'dc' , 1 , vx_vol ); 
      clear Yphs Ydiv Yg

      % only one structure per side and only close to the hippocampus and parahippocampal gyrus and closer to the ventricle 
      % (closing of the deep parahippocampal gyrus)
      Yphg = cat_vol_morph( (Ymsk | Yphx ) & mod(Ya,2)==0, 'l' ) | cat_vol_morph( (Ymsk | Yphx ) & mod(Ya,2)==1, 'l' ); 
      Yphg = Yphg & ~Yphn & cat_vol_morph(NS(Ya,LAB.HC) | NS(Ya,LAB.PH) , 'dd', 3 );

      % to avoid that our hard changes effect surrounding areas we keep only the inner part in this region 
      Yphs = cat_vol_morph( cat_vol_morph(Yphg | NS(Ya,LAB.HC) | NS(Ya,LAB.PH),'dc',10) ,'de',2,vx_vol); % we start with the 5 mm of the Yphi!
      Yphg = Yphg & Yphs; 

      % final corretion
      str  = 0.5; % higher values = stronger correction)
      Ymf  = ( (Ymf./3) .^ (1 - str * smooth3(Yphg))) * 3; 
      Ymf  = max( Ymf , min(3,smooth3(Yphg * 4)) );
    else
      Yphg = false(size(Ymf)); 
    end
    
    
    
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

%==========================================================================
function cdata = estimateWMdepthgradient(CS,cdata)
% _________________________________________________________________________
% Estimates the maximum local gradient of a surface. 
% Major use is the WM depth that grows with increasing sulcal depth. 
% It measures the amount of WM behind the cortex, but more relevant is
% the amount of WM fibers that this region will add to the WM depth. 
% The width of the street next to a house gives not the connectivity of
% this house, but the width of the entrance does!
% This measure can be improved by further information of sulcal depth.
% _________________________________________________________________________

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
              
%=======================================================================
function varargout = cat_vol_genus0opt(Yo,th,limit,debug)
% cat_vol_genus0opt: Voxel-based topology optimization and surface creation 
%   The correction of large defects is often not optimal and this function
%   uses only small corrections. 
% 
%    [Yc,S] = cat_vol_genus0vol(Yo[,limit,debug])
%  
%    Yc    .. corrected volume 
%    Yo    .. original volume 
%    S     .. surface 
%    th    .. threshold for creating surface
%    limit .. maximum number of voxels to correct a defect (default = 30)
%    debug .. print details.  
%

  if nargin < 2, th = 0.5; end
  if nargin < 3, limit = 30; end
  if nargin < 4, debug = 0; end
  
  Yc = Yo; nooptimization = limit==0;  %#ok<NASGU>
  if limit==0
    % use all corrections
    if nargout>1
      txt = evalc(sprintf('[Yc,S.faces,S.vertices] = cat_vol_genus0(Yo,th,nooptimization);'));
    else
      txt = evalc(sprintf('Yc = cat_vol_genus0(Yo,th,nooptimization);'));
    end
    
    if debug, fprintf(txt); end
  else
    % use only some corrections
    txt = evalc(sprintf('Yc = cat_vol_genus0(Yo,th,nooptimization);'));
    
    % remove larger corrections
    Yvxcorr = abs(Yc - (Yo > th))>0; 
    Yvxdef  = spm_bwlabel( double( Yvxcorr ) ); clear Yppiscrc; 
    Yvxdef  = cat_vol_morph(Yvxdef,'l',[inf limit]) > 0; % large corrections that we remove 
    
    if debug
      fprintf(txt); 
      fprintf('  Number of voxels of genus-topocorr: %d\n  Finally used corrections:  %0.2f%%\n', ...
        sum(Yvxcorr(:)) , 100 * sum(Yvxcorr(:) & ~Yvxdef(:)) / sum(Yvxcorr(:)) );
    end
    
    Yc = Yc & ~Yvxdef; 
  
    % final surface creation without correction
    if nargout>1
      evalc(sprintf('[Yt,S.faces,S.vertices] = cat_vol_genus0( single(Yc) ,th,1);')); 
    end
  
  end
  
  varargout{1} = Yc; 
  if nargout>1, varargout{2} = S; end
end

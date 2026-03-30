function [Yth,S,P,res] = cat_surf_createCS4(V,V0,Ym,Yp0,Ya,YMF,Yb0,opt,job)
% ______________________________________________________________________
% Surface creation and thickness estimation.
%
% [Yth,S,P] = cat_surf_createCS4(V,V0,Ym,Ya,YMF,Yb0,opt,job)
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
  def.verb                = cat_get_defaults('extopts.expertgui');         % 0-none, 1-minimal, 2-default, 3-details, 4-debug
  def.surf                = {'lh','rh'};                                   % surface reconstruction setting with {'lh','rh','cb'} 
  def.LAB                 = cat_get_defaults('extopts.LAB');               % brain regions 
  def.useprior            = ''; 
  def.thick_limit         = cat_get_defaults('extopts.thick_limit');       % 6mm upper limit for thickness (similar limit as used in Freesurfer)
  def.thick_measure       = cat_get_defaults('extopts.thick_measure');     % 0-PBT; 1-Tfs (Freesurfer method using mean(TnearIS,TnearOS)) 
  def.foldingcorrection   = cat_get_defaults('extopts.foldingcorrection'); % tickness correction that is influence by folding
  def.fsavgDir            = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces'); 
  def.outputpp.native     = 0;  % output of Ypp map for cortical orientation in EEG/MEG 
  def.outputpp.warped     = 0;
  def.outputpp.dartel     = 0;
  def.vdist               = 2;  % surface resolution parameter
  def.myelinCorrection    = .5; % .25 - sight correction, 1 - maximum correction
  def.create_white_pial   = 2 - (isfield(opt,'surf') && isscalar(opt.surf));   % uses only the quick WM and Pial surface estimation (0-no,1-yes,2-improve)
  def.skip_registration   = isfield(opt,'surf') && isscalar(opt.surf); % skip spherical registration for quick tests
  

  % options that rely on other options
  opt                     = cat_io_updateStruct(def,opt); clear def; 
  opt.vol                 = any(~cellfun('isempty',strfind(opt.surf,'v')));   % only volume-based thickness estimation  
  opt.surf                = cat_io_strrep(opt.surf,'v','');                   % after definition of the 'vol' varialbe we simplify 'surf'
  opt.interpV             = max(0.1,min([opt.interpV,2]));                  % general limitation of the PBT resolution
  opt.reconres            = opt.interpV;
  
  % apply the modified mask from gcut
  % for non-gcut approaches or inverse weighting Yb0 only contains ones
  Yp0 = Yp0 .* (Yb0>0.5);
  Ym  = Ym  .* (Yb0>0.5);


  % enlarge atlas definition 
  [~,I] = cat_vbdist(single(Ya>0)); Ya=Ya(I); clear I;  


  % simple filling
  [Yp0f,Ymf] = fillVentricle(Yp0,Ym,Ya,YMF,vx_vol); clear Ym YMF; 
  

  % prepare file and directory names
  [P,mridir,~,ff] = cat_surf_createCS_fun('setFileNames',V0,job,opt); 
  

  % main loop for each surface structure 
  for si = 1:numel(opt.surf)
   
    % print something
    if si==1, fprintf('\n'); end; fprintf('%s:\n',opt.surf{si});
    
    % prepare longitudinal case if required 
    useprior = cat_surf_createCS_fun('setupprior',opt, job.BIDS(1).surfdir,P,si);

    % reduce for object area
    iscerebellum = cat_io_contains(opt.surf{si},'cb'); 
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
    

    % Myelination corretion for the CS40 that is also part of cat_vol_pbtsimpleCS4
    if opt.myelinCorrection > 0 
      stime = cat_io_cmd(sprintf('  Myelin correction (%0.1f%%)',opt.myelinCorrection*100),'g5'); 
      Vmfs.dt = [16 1]; spm_write_vol(Vmfs, Yp0fs);
      Yp0fs = myelincorrection(Yp0fs, vx_vol, opt,P,Vmfs,si); 
      fprintf('%5.0fs\n',etime(clock,stime)); 
    end  
    
    
    stime = cat_io_cmd(sprintf('  Thickness estimation (%0.2f mm%s)',opt.interpV,native2unicode(179, 'latin1')),'g5'); 
    Vmfs.dt = [16 1]; spm_write_vol(Vmfs, Yp0fs);
    if opt.SRP > 0
      [Yth1i,Yppi] = cat_vol_pbtsimpleCS4(Yp0fs, opt.interpV,...
        struct('myelinCorrection',0, 'verb',1, 'gyrusrecon',~iscerebellum, ...
        'eidist',1,'NVBC',~iscerebellum,'denoise',1 & ~iscerebellum)); 
      Vppm = Vmfs; Vppm.fname = P(si).Pppm; spm_write_vol(Vppm, Yppi);
    else
      cmd = sprintf('CAT_VolThicknessPbt  -correct-voxelsize 0   -median-filter 2   -downsample 0  "%s" "%s" "%s"', ...
        Vmfs.fname, P(si).Pgmt, P(si).Pppm);
      cat_system(cmd,opt.verb-3);
      Vgmt  = spm_vol(P(si).Pgmt); Yth1i = spm_read_vols(Vgmt); 
      % correction of general offset in mm 
      Yth1i = max(0,Yth1i - 0.28); 
    end
    

    % prepare thickness output
    Yth1t = cat_vol_resize(Yth1i,'deinterp',resI);                         % back to original resolution
    Yth1t = cat_vol_resize(Yth1t+1,'dereduceBrain',BB) - 1;                % adding background
    Yth   = max(Yth,Yth1t .* Yside); clear Yth1t                           % save on main image


    if opt.vol
      % only voxelbased thickness map
      S = struct(); P = '';
      if opt.verb<2, fprintf('%5.0fs',etime(clock,stime)); end %#ok<*DETIM>
      continue; 
    end
 
    
    %  surface creation 
    %  --------------------------------------------------------------------
    %  In case of the central surface the frequency patter between sulci and 
    %  gyri is more harmonized and even lower resolution (eg. 1.5 and 2.0) 
    %  for surface reconstruction are  are possible. 
    %  --------------------------------------------------------------------
    if isscalar(opt.surf), time_sr = clock; end % temporary for tests 
    if useprior 
      stime = cat_io_cmd('  Load and refine subject average surface','g5','',opt.verb,stime);
      res.EC(si)       = 0; 
      res.ECmodvx(si)  = res.EC;
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
        gth   = .75;
        gycon = .70; 
        % export map
        opt0  = opt; opt0.reconres = min(0.5,opt.reconres);   % minimum resolution 
        Yppi  = Pcb.Yppi; 
        Yppi  = min(1,max(0,Yppi + .5*(Yppi - smooth3(Yppi))));
        Vppmi = exportPPmap( Yppi , Vmfs, opt0); % write temporary file only for surface reconstruction 
      else
        gth   = .5;
        relwm = sum(max(0,Ymf(:)-2)) / sum(min(1,Ymf(:))); 
        gycon = max(0,1 - max(0,relwm * 3 - 1) * 2);
        gycon = max(0.4,min(.85,.6 - .2*gycon)); 
        % export map
        Vppmi = spm_vol(P(si).Pppm); 
      end
     
      % create "gold-standard" surface 
      cmd = sprintf('CAT_VolMarchingCubes  "%s" "%s" -thresh "%0.2f" -fast', Vppmi.fname, P(si).Pcentral,gth);
      cat_system(cmd ,opt.verb-3);  
      CSG05 = loadSurf(P(si).Pcentral);
      res.area_gt(si) = sum(cat_surf_fun('area',CSG05)) / 100; 
      if ~debug, clear CSG05; end 

      %% Iterative runs to create the central surface.
      areaerr  = .5;    % accepted percentage area changes by topology correction 
      genuserr = 6;     % accaptable EC of the final surface 
      stime = cat_io_cmd(sprintf('  Create initial surface (%0.2f mm)',opt.interpV),'g5','',opt.verb,stime); 
      for ECiter = 1 %:2- 1*iscerebellum % this loop is to allow futher topology iterations if the first one was not good enough
        for thi = 1:round(8/mean(vx_vol*2))  % with this loop we further increase the threshold of the Ypp map by .1 (eg. more WM like surface)
          %% Main initial surface creation with topology correction.
          if exist(P(si).Pcentral,'file'), delete(P(si).Pcentral); end % have to delete it to get useful error messages in case of reprocessing/testing
          gycon = gycon + .05/mean(vx_vol*2).*(thi>1);
   
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
              Vppmi.fname, P(si).Pcentral, gycon, ECiter, Vp0.fname, min(1,floor( .5 / opt.interpV ) )); %#ok<NASGU>
            txt = evalc('cat_system(cmd ,3)'); 
          end
          CST{thi} = loadSurf(P(si).Pcentral); %#ok<AGROW>
        
          % save surface characteristics 
          res.gycon(si,thi)     = gycon; 
          res.genus(si,thi)     = abs(str2double( txt(max(strfind(txt,':'))+1 : max(strfind(txt,'('))-2) ) - 2); 
          res.area_tc(si,thi)   = sum(cat_surf_fun('area',CST{thi})) / 100;
          res.area_uc(si,thi)   = sum(cat_surf_fun('area',CSG{thi})) / 100;
          res.surferr(si,thi)   = abs(res.area_tc(si,thi) - res.area_uc(si,thi)) / res.area_uc(si,thi) * 100 * (res.genus(si,thi)+1)/genuserr;
          res.surferrgt(si,thi) = abs(res.area_tc(si,thi) - res.area_gt(si))     / res.area_gt(si)     * 100 * (res.genus(si,thi)+1)/genuserr - (gycon - .5)/2;
          res.ECf(si,thi)       = str2double( txt(max(strfind(txt,':'))+1 : max(strfind(txt,'('))-2) ); 
          res.ECmodvx(si)       = str2double( txt(max(strfind(txt,'('))+1 : max(strfind(txt,'voxel'))-1) );
         
          %% evaluate surface characteristics
          if gycon > .8  ||  ( res.surferr(si,thi) < areaerr  &&  res.surferrgt(si,thi) < areaerr &&  res.genus(si,thi) < genuserr  )
            break;      
          elseif opt.verb>1
            cat_io_cprintf([.7 0 0.5],sprintf( ...
              '\n    Sulcal blurring for Ypp>%0.2f (EC=%3d,ATE=%5.1f%%,ATE0=%5.1f%%)  ', ...
               res.gycon(si,thi), res.ECf(si,thi), res.surferr(si,thi), res.surferrgt(si,thi) ));
          end
        end

        % evaluate surface characteristics
        if gycon > .85  ||   ( res.surferr(si,thi) < areaerr && res.surferrgt(si,thi) < areaerr  &&  res.genus(si,thi) < genuserr  )
          break; 
        end
      end

      %% use best surface 
      if thi > 1 
        err = res.surferr(si,:) + res.surferrgt(si,:); 
        thi = find( err == min(err), 1, 'first');
        saveSurf( CST{thi} , P(si).Pcentral); 
      end
      res.EC(si)  = res.ECf(si,thi); 
      res.thi(si) = thi; 
      if opt.verb>1
        cat_io_cprintf([0 .5 0],sprintf( ...
                '\n    Final configuration Ypp>%0.2f (EC=%3d,ATE=%5.1f%%,ATE0=%5.1f%%)  ', ...
                res.gycon(si,thi), res.ECf(si,thi), res.surferr(si,thi), res.surferrgt(si,thi) ));
      end
      if res.genus > genuserr
        cat_io_addwarning([mfilename ':EC'],sprintf('Incorrect Euler Number i.e. reminding topological defects (EC=%d). ', ...
          res.EC(si)),1,[1 2*(res.surferr(si,thi) < areaerr)]);
      end
      if res.surferr(si,thi) > areaerr
        cat_io_addwarning([mfilename ':MarchingCubes'],sprintf('Increased topological adjustments of %0.2f%%%%. ',...
          res.surferr(si,thi)), 1,[1 2]);
      end
      
      % load and check surface and 
      % define lower surface resolution based on the surface area of the subject
      % use 80000 as lower limit what is similar to vdist=4
      CS   = loadSurf(P(si).Pcentral); 
      A    = sum( cat_surf_fun( 'area', CS )); % in mm2 
      nCS0 = size(CS.faces,1); 
      nCS  = min( size(CS.faces,1) , max( 80000, min( 360000, round(A * 8 / opt.vdist^2) ))); 
      nCS  = nCS * (1 + 3*iscerebellum);
      
      % surface reduction 
      % - low surface resolution are highly important for fast processing and smaller disk usage
      % - less critical for thickness but surface quality (intensity/position values) get worse
      % - it is important to keep the topology 
      % - a more aggressive 
      stimesr = clock;
      if nCS < nCS0 
        cmd  = sprintf('CAT_SurfReduce -aggr "7" -ratio "%0.2f" "%s" "%s"', nCS/nCS0 , P(si).Pcentral, P(si).Pcentral); 
        cat_system(cmd ,opt.verb-3);  
        CS   = loadSurf(P(si).Pcentral); 
      end
      nCS1 = size(CS.faces,1); 
      clear CS

      if nCS1 > nCS * 1.5 
        cat_io_addwarning('cat_surf_createCS4:ReductionFailed', ...
          sprintf('Surface reduction incomplete due to topology constrains (%d).',(nCS1/1000-nCS/1000) ./ (nCS0/1000-nCS/1000)),0,[0 1]); 
      end 
      cat_io_cmd(sprintf('  Reduce surface (nF: %0.0fk > ~%0.0fk)', nCS0/1000, nCS1/1000),'g5','',opt.verb,stime + (clock - stimesr));
      stime = stimesr; 
      
    end

    
    %% surface deformation 
    % CAT_DeformSurf in combination with smoothing correction CAT_Central2Pial 
    % allows PP values that are closer to the central position (std 0.049 vs. 0.73)
    % what improves also intensity and position values of the inner/outer surfaces!
    if useprior, remove_intersect = ' '; else, remove_intersect = ' -remove_intersect '; end
    stime = cat_io_cmd('  Optimize surface','g5','',opt.verb,stime); 
     cmd = sprintf('CAT_SurfDeform %s -iter 75 -isovalue 0.5 -w1 0.1 -w2 0.1 -w3 1.0 -sigma 1.2 "%s" "%s" "%s"', ...
      remove_intersect, P(si).Pppm, P(si).Pcentral, P(si).Pcentral);
    cat_system(cmd,opt.verb-3);


    % get PBT thickness 
    CS = loadSurf(P(si).Pcentral);
    facevertexcdata = max(eps,cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm)); 
    cat_io_FreeSurfer('write_surf_data',P(si).Ppbt,facevertexcdata);
    

    %% create white and central surfaces
    if opt.create_white_pial  
      stime = cat_io_cmd('  Create pial and white surface','g5','',opt.verb,stime); 
      if useprior || opt.create_white_pial == 1 
        % quick
        cmd = sprintf('CAT_Central2Pial "%s" "%s" "%s" 0.5',P(si).Pcentral,P(si).Ppbt,P(si).Ppial);
        cat_system(cmd,opt.verb-3);
        cmd = sprintf('CAT_Central2Pial "%s" "%s" "%s" -0.5',P(si).Pcentral,P(si).Ppbt,P(si).Pwhite);
        cat_system(cmd,opt.verb-3); 
      else
        % write segmentation if necessary (and delete it later)
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
        CS = loadSurf(P(si).Pcentral);
        facevertexcdata = max(eps, cat_surf_fun('isocolors',Yth1i,CS.vertices,Smat.matlabIBB_mm)); 
        cat_io_FreeSurfer('write_surf_data', P(si).Ppbt, facevertexcdata);
      end
    end


    if isscalar(opt.surf)
      % only for test visualization  
      fprintf('%5.0fs',etime(clock,stime)); 
      Vppm = spm_vol(P(si).Pppm); Yppi  = spm_read_vols(Vppm); 
      if ~exist(P(si).Pthick,'file'), cat_io_FreeSurfer('write_surf_data', P(si).Pthick, facevertexcdata); end
      cat_surf_createCS_fun('quickeval',V0,Vppm,Ymfs,Yppi,CS,P,Smat,res,opt,res.EC(si),si,time_sr,4); 
      if exist(P(si).Pthick,'file'), delete(P(si).Pthick); end
      return
    end
        

    %% skip that part if a prior image is defined
    if ~useprior && ~opt.skip_registration
      % spherical surface mapping of the final corrected surface with minimal areal smoothing
      stime = cat_io_cmd('  Spherical mapping with areal smoothing','g5','',opt.verb,stime); 
      cmd = sprintf('CAT_Surf2Sphere "%s" "%s" %d',P(si).Pcentral,P(si).Psphere,6);
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
    

    % correct thickness based on folding pattern, but smaller thickness values 
    % are corrected less strongly than larger thickness values
    if opt.foldingcorrection
      cmd = sprintf('CAT_SurfCorrectThicknessFolding -slope 1.0 -max "%f" "%s" "%s" "%s"', ...
        opt.thick_limit, P(si).Pcentral, P(si).Pthick, P(si).Pthick);
      cat_system(cmd,opt.verb-3);
      % only minimum thickness to avoid overcorrections
      facevertexcdatafsc = cat_io_FreeSurfer('read_surf_data',P(si).Pthick);
      facevertexcdatafs = min(facevertexcdatafs,facevertexcdatafsc); clear facevertexcdatafsc
      cat_io_FreeSurfer('write_surf_data', P(si).Pthick, min(opt.thick_limit,facevertexcdatafs) ); 
    end


    % final surface evaluation 
    Vppm = spm_vol(P(si).Pppm); Yppi = spm_read_vols(Vppm); 
    res.(opt.surf{si}).createCS_final = cat_surf_fun('evalCS', ...
      loadSurf(P(si).Pcentral), cat_io_FreeSurfer('read_surf_data',P(si).Ppbt), cat_io_FreeSurfer('read_surf_data',P(si).Pthick), ...
      Ymfs, Yppi, P(si).Pcentral, Smat.matlabIBB_mm, debug + (cat_get_defaults('extopts.expertgui')>1), cat_get_defaults('extopts.expertgui')>1);
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
  Vppm = spm_write_vol(Vppm, Ytx );
  spm_smooth(Vppm,Vppm,repmat(.2,1,3));
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
    
  def.levels = 2; 
  def.rangeE = 0.4; % extimation limit extention - smaller better
  def.range  = 0.4; % addition range - larger better but limited by interpolation artifacts
  def.vxs    = 1;

  opt  = cat_io_checkinopt(opt,def); 
  vxc  = .5 - opt.range; % the more PVE the less we need to correct for
  
  range  = min(.4,opt.range);    % defined by the PVE boundary of the AMAP, i.e range of .25 - .75  
  rangeE = min(.4,opt.rangeE);    % defined by the PVE boundary of the AMAP, i.e range of .25 - .75  
  
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
    Ycdl       = (Ycdl - min(vxc,max(0,- (Yp0(YI) - (1.5-offset))*opt.pvet ))) .* (YMM & YMC); 
    [Ycdh,YI]  = cat_vbdist(single(Yp0 < ( 1.5 + offset)), YMC ); 
    Ycdh       = (Ycdh - min(vxc,max(0,- (Yp0(YI) - (1.5+offset))*opt.pvet ))) .* (YMM & YMC); 
    Ycd        = Ycd + .5/hss .* Ycdl  +  .5/hss .* Ycdh; 

    % WM distances
    [Ywdl,YI]  = cat_vbdist(single(Yp0 > ( 2.5 - offset)), YMM );
    Ywdl       = (Ywdl - min(vxc,max(0, (Yp0(YI) - (2.5-offset))*opt.pvet ))) .* (YMM & YMC); 
    [Ywdh,YI]  = cat_vbdist(single(Yp0 > ( 2.5 + offset)), YMM ); 
    Ywdh       = (Ywdh - min(vxc,max(0, (Yp0(YI) - (2.5+offset))*opt.pvet ))) .* (YMM & YMC); 
    Ywd        = Ywd + .5/hss .* Ywdl  +  .5/hss .* Ywdh;
  end

  % endpoint PVE correction 
  Ycd = Ycd - min(.5,max(-.0, (Yp0>2.5-0*opt.rangeE & Yp0<2.5+opt.rangeE) .* ( (Yp0-2.5)*opt.pvet*1 + 0*opt.vxs/2) ));
  Ywd = Ywd + min(.0,max(-.5, (Yp0>1.5-opt.rangeE & Yp0<1.5+0*opt.rangeE) .* ( (Yp0-1.5)*opt.pvet*1 - 0*opt.vxs/2) ));
end
%==========================================================================
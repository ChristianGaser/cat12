function cat_run_job1639(job,tpm,subj)
% run CAT 
% ______________________________________________________________________
%
% Initialization functions for the CAT preprocessing
%  * creation of the subfolder structure (if active)
%  * check of image resolution (avoid scans with very low resolution)
%  * interpolation 
%  * affine preprocessing (APP)
%    >> cat_run_job_APP_init
%    >> cat_run_job_APP_final
%  * affine registration
%  * initial SPM preprocessing
%
%   cat_run_job(job,tpm,subj)
% 
%   job  .. SPM job structure with main parameter
%   tpm  .. tissue probability map (hdr structure)
%   subj .. file name
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

%#ok<*WNOFF,*WNON>

  job.test_warnings = 0; % just for tests

  % if there is a breakpoint in this file set debug=1 and do not clear temporary variables 
  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end
  
  clearvars -global cat_err_res;
  global cat_err_res; % for CAT error report
  cat_err_res.stime        = clock; 
  cat_err_res.cat_warnings = cat_io_addwarning('reset'); % reset warnings     
  stime  = clock; 
  stime0 = stime; % overall processing time
  
  % create subfolders if not exist
  pth = spm_fileparts(job.channel(1).vols{subj}); 
  [mrifolder, reportfolder, surffolder, labelfolder] = cat_io_subfolders(job.channel(1).vols{subj},job);

  if job.extopts.subfolders
  
    folders = {mrifolder,reportfolder};
    warning('off', 'MATLAB:MKDIR:DirectoryExists');
    for i=1:numel(folders)
      if ~exist(fullfile(pth,folders{i}),'dir')
        [stat, msg] = mkdir(fullfile(pth,folders{i}));
        if ~stat
          fprintf('%s: Error while creating directory %s\n\n\n',msg,fullfile(pth,folders{i}));
          return
        end
      end
    end
  
    if ~exist(fullfile(pth,surffolder),'dir') && job.output.surface
      [stat, msg] = mkdir(fullfile(pth,surffolder));
      if ~stat
        fprintf('%s: Error while creating directory %s\n\n\n',msg,fullfile(pth,surffolder));
        return
      end
    end
  
    if ~exist(fullfile(pth,labelfolder),'dir') && job.output.ROI
      [stat, msg] = mkdir(fullfile(pth,labelfolder));
      if ~stat
        fprintf('%s: Error while creating directory %s\n\n\n',msg,fullfile(pth,labelfolder));
        return
      end
    end
    
  end
  
  % create subject-wise diary file with the command-line output
  [pp,ff,ee,ex] = spm_fileparts(job.data{subj});  %#ok<ASGLU>
  catlog = fullfile(pth,reportfolder,['catlog_' ff '.txt']);
  if exist(catlog,'file'), delete(catlog); end % write every time a new file, turn this off to have an additional log file
  
  % check if not another diary is already written that is not the default- or catlog-file. 
  if ~strcmpi(spm_check_version,'octave')
    olddiary = spm_str_manip( get(0,'DiaryFile') , 't');
    usediary = ~isempty(strfind( olddiary , 'diary' )) || ~isempty(strfind( olddiary , 'catlog_' )); 
    if usediary
      diary(catlog); 
      diary on; 
    else  
      cat_io_cprintf('warn',sprintf('External diary log is written to "%s".\n',get(0,'DiaryFile'))); 
    end
  else
    usediary = 0; 
  end
  
  % print current CAT release number and subject file
  [n,r] = cat_version;
  str  = sprintf('%s r%s: %d/%d',n,r,subj,numel(job.channel(1).vols));
  str2 = spm_str_manip(job.channel(1).vols{subj}(1:end-2),['a' num2str(70 - length(str))]);
  cat_io_cprintf([0.2 0.2 0.8],'\n%s\n%s: %s%s\n%s\n',...
        repmat('-',1,72),str,...
        repmat(' ',1,70 - length(str) - length(str2)),str2,...
        repmat('-',1,72));
  clear r str str2
    
  
  %  -----------------------------------------------------------------
  %  separation of full CAT preprocessing and SPM segmentation
  %  preprocessing (running DARTEL and PBT with SPM segmentation)
  %  -----------------------------------------------------------------
  [pp,ff,ee,ex] = spm_fileparts(job.data{subj}); 
  if exist(fullfile(pp,['c1' ff(3:end) ee]),'file') && ...
     exist(fullfile(pp,['c2' ff(3:end) ee]),'file') && ...
     exist(fullfile(pp,['c3' ff(3:end) ee]),'file') && ...
     exist(fullfile(pp,[ff(3:end) '_seg8.mat']),'file') && strcmp(ff(1:2),'c1')
     
      job.data{subj}          = fullfile(pp,[ff ee]); 
      job.channel.vols{subj}  = fullfile(pp,[ff ee]); 

      % prepare SPM preprocessing structure 
      images = job.channel(1).vols{subj};
      for n=2:numel(job.channel)
        images = char(images,job.channel(n).vols{subj});
      end

      obj.image    = spm_vol(images);
      obj.fwhm     = job.opts.fwhm;
      obj.biasreg  = cat(1,job.opts.biasreg);
      obj.biasfwhm = cat(1,job.opts.biasfwhm);
      obj.tol      = job.opts.tol;
      obj.lkp      = [];
      obj.reg      = job.opts.warpreg;
      obj.samp     = job.opts.samp;              
      spm_check_orientations(obj.image);

      if all(isfinite(cat(1,job.tissue.ngaus)))
          for k=1:numel(job.tissue)
              obj.lkp = [obj.lkp ones(1,job.tissue(k).ngaus)*k];
          end
      end

      Pseg8 = fullfile(pp,[ff(3:end) '_seg8.mat']); 
      if ~exist(Pseg8,'file')
        error('cat_run_job1639:SPMpp_MissSeg8mat','Can''t find "%s" file!',Pseg8);
      end
      res = load(Pseg8);
      
      % load tpm priors 
      tpm = spm_load_priors8(res.tpm);
      obj.lkp = res.lkp; 
      obj.tpm = tpm; 
      
      % Special cases with different class numbers in case of SPM input
      if max(obj.lkp)==6
      % default cases
      elseif max(obj.lkp)==3
        cat_io_addwarning('SPMpp_PostMortem','Detected only 3 classes that are interpretated as GM, WM, and CSF/background.',0,[0 1])
      elseif max(obj.lkp)==4 
        cat_io_addwarning('SPMpp_SkullStripped','Detected only 4 classes that are interpretated as GM, WM, CSF, and background',0,[0 1])
      else
        cat_io_addwarning('SPMpp_AtypicalClsNumber',sprintf('Atypical number of input classes (max(lkp)=%d).',max(obj.lkp)),2,[0 1])
      end
      
      
      cfname  = fullfile(pp,[ff ee]);
      ofname  = fullfile(pp,[ff(3:end) ee]); 
      nfname  = fullfile(pp,mrifolder,['n' ff '.nii']); 
      copyfile(ofname,nfname,'f'); 

      Ysrc0    = single(spm_read_vols(obj.image)); 
      Ylesion  = single(isnan(Ysrc0) | isinf(Ysrc0) | Ysrc0==0); clear Ysrc0;
      
     
      job.channel(1).vols{subj}  = [nfname ex];
      job.channel(1).vols0{subj} = [ofname ex];
      res.image  = spm_vol([nfname ex]);
      res.image0 = spm_vol([ofname ex]);
      res.imagec = spm_vol([cfname ex]);
      res.spmpp  = 1; 
      job.spmpp  = 1; 
      
      cat_err_res.obj = obj; 
  else

      %  -----------------------------------------------------------------
      %  check resolution properties
      %  -----------------------------------------------------------------
      %  There were some images that should not be processed. So we have  
      %  to check for large slice thickness and low spatial resolution.
      %  RD201909: I tried 4x4x4 and 1x1x8 mm data with default and NLM 
      %            interpolation. Also NLM shows less edges and more 
      %            correct surfaces, the thickness results are worse and 
      %            the limits are ok.
      %  RD202007: Low resolution data is now allowed if ignoreErrer > 1.
      %            Tested again NLM and Boeseflug interpolation and there
      %            are many artefacts and simple spline interpolation is 
      %            more save. 
      %  RD202107: Print warning for reslimit/2 and alert for reslimit.  
      %  -----------------------------------------------------------------
      for n=1:numel(job.channel) 
        V = spm_vol(job.channel(n).vols{subj});
        vx_vol = sqrt(sum(V.mat(1:3,1:3).^2));

        % maximum [ slice-thickness , volume^3 , anisotropy ]
        reslimits    = [5 4 8];
        
        % too thin slices
        if any( vx_vol > reslimits(1)/2 ) || job.test_warnings
          mid = [mfilename 'cat_run_job:TooLowResolution']; 
          msg = sprintf(['Voxel resolution should be better than %d mm in any dimension for \\\\n' ...
            'reliable preprocessing! This image has a resolution of %0.2fx%0.2fx%0.2f mm%s. '], ... 
            reslimits(1),vx_vol,native2unicode(179, 'latin1'));
          cat_io_addwarning(mid,msg,1 + any( vx_vol > reslimits(1) ) ,[0 1],vx_vol);
        end
        
        % too small voxel volume (smaller than 3x3x3 mm3)
        if prod(vx_vol) > (reslimits(2)/2)^3 || job.test_warnings
          mid = [mfilename 'cat_run_job:TooLargeVoxelVolume']; 
          msg = sprintf(['Voxel volume should be smaller than %d mm%s (around %dx%dx%d mm%s) for \\\\n' ...
                  'reliable preprocessing! This image has a voxel volume of %0.2f mm%s. '], ...
                  reslimits(2)^3,native2unicode(179, 'latin1'),reslimits(2),reslimits(2),reslimits(2),...
                  native2unicode(179, 'latin1'),prod(vx_vol),native2unicode(179, 'latin1'));
          cat_io_addwarning(mid,msg,1 + (prod(vx_vol) > reslimits(2)^3),[0 1],vx_vol);
        end
        
        % anisotropy
        if max(vx_vol) / min(vx_vol) > reslimits(3)/2 || job.test_warnings
          mid = [mfilename 'cat_run_job:TooStrongAnisotropy'];
          msg = sprintf(['Voxel anisotropy (max(vx_size)/min(vx_size)) should be smaller than %d for \\\\n' ...
                  'reliable preprocessing! This image has a resolution %0.2fx%0.2fx%0.2f mm%s \\\\nand a anisotropy of %0.2f. '], ...
                  reslimits(3),vx_vol,native2unicode(179, 'latin1'),max(vx_vol)/min(vx_vol));
          cat_io_addwarning(mid,msg,1 + (max(vx_vol) / min(vx_vol) > reslimits(3)/3),[0 1],vx_vol);
        end
      end
      
      % save original file name 
      for n=1:numel(job.channel) 
        job.channel(n).vols0{subj} = job.channel(n).vols{subj};
      end
      
     
      % always create the n*.nii image because of the real masking of the
      % T1 data for spm_preproc8 that includes rewriting the image!
      for n=1:numel(job.channel) 
        [pp,ff,ee] = spm_fileparts(job.channel(n).vols{subj}); 
        ofname  = fullfile(pp,[ff ee]); 
        nfname  = fullfile(pp,mrifolder,['n' ff '.nii']);
        if strcmp(ee,'.nii')
          if ~copyfile(ofname,nfname,'f')
            spm('alert!',sprintf('ERROR: Check file permissions for folder %s.\n',fullfile(pp,mrifolder)),'',spm('CmdLine'),1);
          end
        elseif strcmp(ee,'.img')
          V = spm_vol(job.channel(n).vols{subj});
          Y = spm_read_vols(V);
          V.fname = nfname; 
          spm_write_vol(V,Y);
          clear Y; 
        end
        job.channel(n).vols{subj} = nfname;

        %% denoising
        if job.extopts.NCstr~=0
          NCstr.labels = {'none','full','light','medium','strong','heavy'};
          NCstr.values = {0 1 2 -inf 4 5}; 
          stime = cat_io_cmd(sprintf('SANLM denoising (%s)',...
            NCstr.labels{find(cell2mat(NCstr.values)==job.extopts.NCstr,1,'first')}));
          cat_vol_sanlm(struct('data',nfname,'verb',0,'prefix','','NCstr',job.extopts.NCstr)); 
          fprintf('%5.0fs\n',etime(clock,stime));   
        end

        %% skull-stripping detection
        %  ------------------------------------------------------------
        %  Detect skull-stripping or defaceing because it strongly 
        %  affects SPM segmentation that expects gaussian distribution! 
        %  If a brain mask was used than we expect 
        %   - many zeros (50% for small background - 80-90% for large backgrounds)
        %   - a smaller volume because of missing skull (below 2500 cm3)
        %   - only one object (the masked regions)
        %   - only one background (not in every case?)
        %   - less variance of tissue intensity (only 3 brain classes)
        %  ------------------------------------------------------------
        VFn   = spm_vol(nfname); 
        YF    = spm_read_vols(VFn); 
        Oth   = cat_stat_nanmean(YF(YF(:)~=0 & YF(:)>cat_stat_nanmean(YF(:)))); 
        F0vol = cat_stat_nansum(YF(:)~=0) * prod(vx_vol) / 1000; 
        F0std = cat_stat_nanstd(YF(YF(:)>0.5*Oth & YF(:)>0)/Oth); 
        YFC = YF~=0; 
        if sum(YFC(:)>0)<numel(YFC)*0.9 && sum(YFC(:)>0)>numel(YFC)*0.1  % if there is a meanful background
          YFC = ~cat_vol_morph(YF~=0,'lc',1);                            % close noisy background
        end
        [YL,numo] = spm_bwlabel(double(YF~=0),26);  clear YL;            %#ok<ASGLU> % number of objects
        [YL,numi] = spm_bwlabel(double(YFC==0),26); clear YL;            %#ok<ASGLU> % number of background regions 
        ppe.affreg.skullstrippedpara = [sum(YF(:)==0)/numel(YF) numo numi F0vol F0std]; 
        ppe.affreg.skullstripped = ...
          ppe.affreg.skullstrippedpara(1)>0.5 && ...                     % many zeros
          ppe.affreg.skullstrippedpara(2)<15  && ...                     % only a few objects
          ppe.affreg.skullstrippedpara(3)<10 && ...                      % only a few background regions 
          F0vol<2500 && F0std<0.5;                                       % many zeros and not too big
        ppe.affreg.skullstripped = ppe.affreg.skullstripped || ...
          sum([ppe.affreg.skullstrippedpara(1)>0.8 F0vol<1500 F0std<0.4])>1; % or 2 extreme values
        if ~debug, clear YFC F0vol F0std numo numi; end 
        % not automatic detection in animals
        ppe.affreg.skullstripped = ppe.affreg.skullstripped && strcmp(job.extopts.species,'human');
        
        %% high intensity background (MP2Rage)
        % RD202008: improved object detection with gradient
        [YF,R]= cat_vol_resize(YF,'reduceV',vx_vol,2,32,'meanm');
        YFm   = cat_stat_histth(YF,[0.95 0.95],struct('scale',[0 1])); 
        Yg    = cat_vol_grad(YFm,R.vx_volr); 
        gth   = max(0.05,min(0.2,median(Yg(Yg(:)>median(Yg(Yg(:)>0.1)))))); % object edge threshold
        YOB   = abs(YFm)>0.1 & Yg>gth;                                   % high intensity object edges 
        YOB   = cat_vol_morph(YOB,'ldc',8/mean(R.vx_volr));                 % full object
        % image pricture frame to test for high intensity background in case of defaced data
        hd    = max(3,round(0.03 * size(YF))); 
        YCO   = true(size(YF)); YCO(hd(1):end-hd(1)+1,hd(2):end-hd(2)+1,hd(3):end-hd(3)+1) = false; 
        % background  
        if sum(YOB(:)>0)<numel(YOB)*0.9 && sum(YOB(:)>0)>numel(YOB)*0.1  % if there is a meanful background
          YBG = ~cat_vol_morph(YOB,'lc',2/mean(R.vx_volr));                 % close noisy background
        else
          YBG = ~cat_vol_morph(YOB,'lc',2/mean(R.vx_volr));  
          msg = [mfilename 'Detection of background failed.']; 
          cat_io_addwarning('cat_run_job:failedBGD',msg,1,[0 1]);
        end
        ppe.affreg.highBGpara = [ ...
          cat_stat_nanmedian( YFm( YBG(:) > 1/3 )) ... normal background
          cat_stat_nanmedian( YFm( YCO(:) > 1/3 )) ... pricture frame background 
          cat_stat_nanstd( YFm(YBG(:)) > 1/3)]; % I am not sure if we should use the std, because inverted images are maybe quite similar
        ppe.affreg.highBG     = ...
          ppe.affreg.highBGpara(1) > 1/5 || ...
          ppe.affreg.highBGpara(2) > 1/5; 
        
        
        %% Interpolation
        %  -----------------------------------------------------------------
        %  The interpolation can help reducing problems for morphological
        %  operations for low resolutions and strong isotropic images. 
        %  Especially for Dartel registration a native resolution larger than the Dartel 
        %  resolution helps to reduce normalization artifacts of the
        %  deformations. Furthermore, even if artifacts can be reduced by the final smoothing
        %  it is much better to avoid them.  

        % prepare header of resampled volume
        Vi        = spm_vol(job.channel(n).vols{subj}); 
        vx_vol    = sqrt(sum(Vi.mat(1:3,1:3).^2));
%        vx_vol    = round(vx_vol*10^2)/10^2; % avoid small differences 

        % we have to look for the name of the field due to the GUI job struct generation! 
        restype   = char(fieldnames(job.extopts.restypes));
        switch restype
          case 'native'
            vx_voli  = vx_vol;
          case 'fixed' 
            vx_voli  = min(vx_vol ,job.extopts.restypes.(restype)(1) ./ ...
                       ((vx_vol > (job.extopts.restypes.(restype)(1)+job.extopts.restypes.(restype)(2)))+eps));
            vx_voli  = max(vx_voli,job.extopts.restypes.(restype)(1) .* ...
                       ( vx_vol < (job.extopts.restypes.(restype)(1)-job.extopts.restypes.(restype)(2))));
          case 'best'
            best_vx  = max( min(vx_vol) ,job.extopts.restypes.(restype)(1)); 
            vx_voli  = min(vx_vol ,best_vx ./ ((vx_vol > (best_vx + job.extopts.restypes.(restype)(2)))+eps));
          case 'optimal'
            %%
            aniso   = @(vx_vol) (max(vx_vol) / min(vx_vol)^(1/3))^(1/3);                                              % penetration factor
            volres  = @(vx_vol) repmat( round( aniso(vx_vol) * prod(vx_vol)^(1/3) * 10)/10 , 1 , 3);                  % volume resolution
            optresi = @(vx_vol) min( job.extopts.restypes.(restype)(1) , max( median(vx_vol) , volres(vx_vol) ) );    % optimal resolution 
            optdiff = @(vx_vol) abs( vx_vol - optresi(vx_vol) ) < job.extopts.restypes.(restype)(2);                  % tolerance limites
            optimal = @(vx_vol) vx_vol .* optdiff(vx_vol) + optresi(vx_vol) .* (1 - optdiff(vx_vol) );                % final optimal resolution 
            vx_voli = optimal(vx_vol); 
          otherwise 
            error('cat_run_job:restype','Unknown resolution type ''%s''. Choose between ''fixed'',''native'',''optimal'', and ''best''.',restype)
        end

        % interpolation 
        if any( (vx_vol ~= vx_voli) )  
          stime = cat_io_cmd(sprintf('Internal resampling (%4.2fx%4.2fx%4.2fmm > %4.2fx%4.2fx%4.2fmm)',vx_vol,vx_voli));
         
          if 1
            imat      = spm_imatrix(Vi.mat); 
            Vi.dim    = round(Vi.dim .* vx_vol./vx_voli);
            imat(7:9) = vx_voli .* sign(imat(7:9));
            Vi.mat    = spm_matrix(imat); clear imat; 
            Vn        = spm_vol(job.channel(n).vols{subj}); 
            cat_vol_imcalc(Vn,Vi,'i1',struct('interp',2,'verb',0,'mask',-1));
          else
            %% Small improvement for CAT12.9 that uses the cat_vol_resize function rather than the simple interpolation. 
            %  However, postive effects only in case of strong reductions >2, ie. it is nearly useless.  
            jobr              = struct(); 
            jobr.data         = {Vi.fname}; 
            jobr.interp       = -3005; % spline with smoothing in case of downsampling;  default without smoothing -5; 
            jobr.verb         = debug; 
            jobr.lazy         = 0; 
            jobr.prefix       = ''; 
            jobr.restype.res  = vx_voli; % use other resolution for test  
            Pr = cat_vol_resize(jobr); 

            if 0
              % test reinterpolation and estimate the RMSE
              jobr.data         = Pr.res;
              jobr              = rmfield(jobr,'restype'); 
              jobr.restype.Pref = {V.fname};  
              jobr.prefix       = 'I'; 
              Pre = cat_vol_resize(jobr); 
              disp('.'); 

              Yo = spm_read_vols(V);
              Yr = spm_read_vols(spm_vol(Pre.res{1}));
              fprintf('%16.8f\n',sqrt( mean((Yo(:) - Yr(:)).^2)));
            end            
          end
          vx_vol = vx_voli;
        
          fprintf('%5.0fs\n',etime(clock,stime));     
        else
          vx_vol = sqrt(sum(Vi.mat(1:3,1:3).^2));
        end

        clear Vi Vn;
        
        
        %% APP bias correction (APPs1 and APPs2)
        %  ------------------------------------------------------------
        %  Bias correction is essential for stable affine registration 
        %  but also the following preprocessing. This approach uses the
        %  SPM Unified segmentation for intial bias correction of the 
        %  input data with different FWHMs (low to high frequency) and 
        %  resolutions (low to high).
        %  As far as SPM finally also gives us a nice initial segmentation 
        %  why not use it for a improved maximum based correction?!
        %  ------------------------------------------------------------
        if (job.extopts.APP==1 || job.extopts.APP==2) 
           job.subj = subj;
           [Ym,Ybg,WMth] = cat_run_job_APP_SPMinit(job,tpm,ppe,n,...
             ofname,nfname,mrifolder,ppe.affreg.skullstripped);
        end
        
        
      end

      
      % prepare SPM preprocessing structure 
      images = job.channel(1).vols{subj};
      for n=2:numel(job.channel)
        images = char(images,job.channel(n).vols{subj});
      end
      obj.image    = spm_vol(images);
      obj.fwhm     = job.opts.fwhm;
      obj.biasreg  = job.opts.biasreg;
      obj.biasfwhm = job.opts.biasfwhm;
      obj.tpm      = tpm;   
      obj.reg      = job.opts.warpreg;
      obj.samp     = job.opts.samp;              
      obj.tol      = job.opts.tol;
      obj.lkp      = [];
      if ~strcmp('human',job.extopts.species) 
        % RD202105: There are multiple problems in primates and increased 
        %           accuracy is maybe better (eg. 0.5 - 0.66) 
        scannorm   = 0.7; %prod(obj.image.dims .* vx_vol).^(1/3) / 20; % variance from typical field fo view to normalized parameter 
        obj.samp   = obj.samp * scannorm; % normalize by voxel size 
        obj.fwhm   = obj.fwhm * scannorm; 
      end
      if all(isfinite(cat(1,job.tissue.ngaus)))
        for k=1:numel(job.tissue)
          obj.lkp = [obj.lkp ones(1,job.tissue(k).ngaus)*k];
        end
      end
      spm_check_orientations(obj.image);
      cat_err_res.obj = obj; 
      
      %% Initial affine registration.
      %  -----------------------------------------------------------------
      [pp,ff] = spm_fileparts(job.channel(1).vols{subj});
      Pbt     = fullfile(pp,mrifolder,['brainmask_' ff '.nii']);
      Pb      = char(job.extopts.brainmask);
      Pt1     = char(job.extopts.T1);
      
      if ~isempty(job.opts.affreg)      
        % first affine registration (with APP)
        
        % load template and remove the skull if the image is skull-stripped
        try 
          VG = spm_vol(Pt1);
        catch
          pause(rand(1))
          VG = spm_vol(Pt1);
        end
        VF = spm_vol(obj.image(1));

        % skull-stripping of the template
        if ppe.affreg.skullstripped || job.extopts.gcutstr<0
          % print a warning for all users that did not turn off skull-stripping 
          % because processing of skull-stripped data is not standard!
          if job.extopts.gcutstr>=0 || job.test_warnings
            msg = [...
              'Detected skull-stripped or strongly masked image. Skip APP. \\n' ...
              'Use skull-stripped initial affine registration template and \\n' ...
              'TPM without head tissues (class 4 and 5)!']; 
            if job.extopts.verb>1 && job.extopts.expertgui
               msg = [msg sprintf(['\\\\n  BG: %0.2f%%%%%%%% zeros; %d object(s); %d background region(s) \\\\n' ...
                '  %4.0f cm%s; normalized SD of all tissues %0.2f'],...
                ppe.affreg.skullstrippedpara(1:4),native2unicode(179, 'latin1'),ppe.affreg.skullstrippedpara(5))]; 
            end  
            cat_io_addwarning([mfilename 'cat_run_job:skullStrippedInputWithSkullStripping'],msg,1,[0 1],ppe.affreg.skullstrippedpara);
            
          elseif job.extopts.gcutstr<0 && ~ppe.affreg.skullstripped || job.test_warnings
            cat_io_addwarning([mfilename 'cat_run_job:noSkullStrippingButSkull'],[...
                'Skull-Stripping is deactivated but skull was detected. \\n' ...
                'Go on without skull-stripping what possibly will fail.'],1,[0 1],ppe.affreg.skullstrippedpara);
          end

          % skull-stripping of the template
          VB = spm_vol(Pb);
          [VB2,YB] = cat_vol_imcalc([VG,VB],Pbt,'i1 .* i2',struct('interp',3,'verb',0,'mask',-1)); 
          VB2.dat(:,:,:) = eval(sprintf('%s(YB/max(YB(:))*255);',spm_type(VB2.dt))); 
          VB2.pinfo      = repmat([1;0],1,size(YB,3));
          VG             = cat_spm_smoothto8bit(VB2,0.5);
          clear VB2 YB; 
        end

        % Rescale images so that globals are better conditioned
        VF.pinfo(1:2,:) = VF.pinfo(1:2,:)/spm_global(VF);
        VG.pinfo(1:2,:) = VG.pinfo(1:2,:)/spm_global(VG);

        % APP step 1 rough bias correction and preparation of the affine 
        % registration
        % --------------------------------------------------------------
        % Already for the rough initial affine registration a simple  
        % bias corrected and intensity scaled image is required, because
        % large head intensities can disturb the whole process.
        % --------------------------------------------------------------
        % ds('l2','',vx_vol,Ym, Yt + 2*Ybg,obj.image.private.dat(:,:,:)/WMth,Ym,60)
        if (job.extopts.APP == 1070 || job.extopts.APP == 1144) && ~ppe.affreg.highBG 
          stime = cat_io_cmd('APP: Rough bias correction'); 
          try
            Ysrc  = single(obj.image.private.dat(:,:,:)); 
            if job.extopts.APP == 1070 
              [Ym,Yt,Ybg,WMth] = cat_run_job_APP_init1070(Ysrc,vx_vol,job.extopts.verb); %#ok<ASGLU>
            else % new version R1144
              [Ym,Yt,Ybg,WMth,bias,Tth,ppe.APPi] = cat_run_job_APP_init(...
                Ysrc,vx_vol,struct('verb',job.extopts.verb,'APPstr',job.opts.biasstr));  %#ok<ASGLU>
            end
          catch apperr
            %% very simple affine preprocessing ... only simple warning
            cat_io_addwarning([mfilename ':APPerror'],'APP failed. Use simple scaling.',1,[0 0],apperr);
            [Ym,Yt,Ybg,WMth] = APPmini(obj,VF); %#ok<ASGLU>
            if cat_get_defaults('extopts.send_info')
              urlinfo = sprintf('%s%s%s%s%s%s%s%s%s%s',cat_version,'%2F',computer,'%2F','errors',...
                 '%2F','cat_run_job:failedAPP','%2F','WARNING: APP failed. Use simple scaling.','cat_run_job');
              cat_io_send_to_server(urlinfo);
            end
          end
          APPRMS = checkAPP(Ym,Ysrc); 
          if APPRMS>1 || job.test_warnings
            if job.extopts.ignoreErrors < 1 
              fprintf('\n'); 
              error('cat_run_job:APPerror','Detect problems in APP preprocessing (APPRMS: %0.4f). Do not use APP results. ',APPRMS);
            else
              cat_io_addwarning([mfilename ':APPerror'],...
                sprintf('Detect problems in APP preprocessing (APPRMS: %0.4f). \\\\nDo not use APP results. ',APPRMS),1,[0 1],APPRMS);
            end 
          end
      
          if ~debug, clear Yt; end

          if job.extopts.setCOM && ~( isfield(job,'useprior') && ~isempty(job.useprior) ) && ~ppe.affreg.highBG
            
          else
            stime = cat_io_cmd('Affine registration','','',1,stime); 
          end
          
          % write data to VF
          VF.dt         = [spm_type('UINT8') spm_platform('bigend')];
          VF.dat(:,:,:) = cat_vol_ctype(Ym * 200,'uint8'); 
          VF.pinfo      = repmat([1;0],1,size(Ym,3));
          clear WI; 

          % smoothing
          resa  = obj.samp*2; % definine smoothing by sample size
          VF1   = spm_smoothto8bit(VF,resa);
          VG1   = spm_smoothto8bit(VG,resa);

        elseif job.extopts.setCOM && ~( isfield(job,'useprior') && ~isempty(job.useprior) ) && ~ppe.affreg.highBG
          % standard approach (no APP) with static resa value and no VG smoothing
          stime = cat_io_cmd('Coarse affine registration');
          resa  = 8;
          VF1   = spm_smoothto8bit(VF,resa);
          VG1   = VG; 
          [Ym,Yt,Ybg,WMth] = APPmini(obj,VF);
        else
          stime = cat_io_cmd('Skip initial affine registration due to high-intensity background','','',1);  
          VF = spm_vol(obj.image(1));
          [Ym,Yt,Ybg,WMth] = APPmini(obj,VF); %#ok<ASGLU>
        end

        %% prepare affine parameter 
        aflags     = struct('sep',obj.samp,'regtype','subj','WG',[],'WF',[],'globnorm',1);
        aflags.sep = max(aflags.sep,max(sqrt(sum(VG(1).mat(1:3,1:3).^2))));
        aflags.sep = max(aflags.sep,max(sqrt(sum(VF(1).mat(1:3,1:3).^2))));

        % use affine transformation of given (average) data for longitudinal mode
        if isfield(job,'useprior') && ~isempty(job.useprior)
          priorname = job.useprior{1};
          [pp,ff,ee,ex] = spm_fileparts(priorname);  %#ok<ASGLU>
          catxml = fullfile(pp,reportfolder,['cat_' ff '.xml']);
          
          % check that file exists and get affine transformation
          if exist(catxml,'file')
            fprintf('\nUse affine transformation from:\n%s\n',priorname);
            stime    = cat_io_cmd(' ',' ','',job.extopts.verb); 
            xml      = cat_io_xml(catxml);
            % sometimes xml file does not contain affine transformation
            if ~isfield(xml,'SPMpreprocessing')
              cat_io_cprintf('warn',sprintf('WARNING: File %s does not contain successful affine transformation. Use individual affine transformation\n',catxml));
              Affine   = eye(4); 
              useprior = 0;
            else
              Affine   = xml.SPMpreprocessing.Affine;
              affscale = 1;
              useprior = 1;
            end
          else
            cat_io_cprintf('warn',sprintf('WARNING: File %s not found. Use individual affine transformation\n',catxml));
            useprior = 0;
          end
        else
          Affine   = eye(4); 
          useprior = 0;
        
          % correct origin using COM and invert translation and use it as starting value
          if job.extopts.setCOM && ~ppe.affreg.highBG 
            fprintf('%5.0fs\n',etime(clock,stime)); stime = clock;  
            Affine_com        = cat_vol_set_com(VF1);
            Affine_com(1:3,4) = -Affine_com(1:3,4);
          else
            Affine_com = eye(4);
          end

          if ~ppe.affreg.highBG ... strcmp('human',job.extopts.species) && 
            % affine registration
            try
              spm_plot_convergence('Init','Coarse affine registration','Mean squared difference','Iteration');
            catch
              spm_chi2_plot('Init','Coarse affine registration','Mean squared difference','Iteration');
            end

            warning off
            try 
              [Affine0, affscale]  = spm_affreg(VG1, VF1, aflags, Affine_com); Affine = Affine0;
            catch
              affscale = 0; 
            end
            if affscale>3 || affscale<0.5
              cat_io_cmd('Coarse affine registration failed. Try fine affine registration.','','',1,stime);
              Affine = Affine_com; 
            end
            warning on
          end
        end
        
        %% APP step 2 - brainmasking and second tissue separated bias correction  
        %  ---------------------------------------------------------
        %  The second part of APP maps a brainmask to native space and 
        %  refines it by morphologic operations and region-growing to
        %  adapt for worse initial affine alignments. It is important
        %  that the mask covers the whole brain, whereas additional
        %  masked head is here less problematic.
        %  ---------------------------------------------------------
        %    ds('l2','',vx_vol,Ym,Yb,Ym,Yp0,90)
        
        
        % fine affine registration 
        if ~useprior && ~ppe.affreg.highBG ... strcmp('human',job.extopts.species) && 
          aflags.sep = obj.samp/2; 
          aflags.sep = max(aflags.sep,max(sqrt(sum(VG(1).mat(1:3,1:3).^2))));
          aflags.sep = max(aflags.sep,max(sqrt(sum(VF(1).mat(1:3,1:3).^2))));
          
          stime = cat_io_cmd('Affine registration','','',1,stime); 
          if job.extopts.APP>0
            VF.dt         = [spm_type('UINT8') spm_platform('bigend')];
            VF.pinfo      = repmat([1;0],1,size(Ym,3));
            VF.dat(:,:,:) = cat_vol_ctype(Ym*200); 
          end
          VF1 = spm_smoothto8bit(VF,aflags.sep);
          VG1 = spm_smoothto8bit(VG,aflags.sep);

          try
            spm_plot_convergence('Init','Affine registration','Mean squared difference','Iteration');
          catch
            spm_chi2_plot('Init','Affine registration','Mean squared difference','Iteration');
          end
          warning off
          [Affine1,affscale1] = spm_affreg(VG1, VF1, aflags, Affine, affscale);
          warning on
          if ~any(any(isnan(Affine1(1:3,:)))) && affscale1>0.5 && affscale1<3, Affine = Affine1; end
        else
          Affine1 = Affine; 
        end
        clear VG1 VF1
       
      else
        VF = spm_vol(obj.image(1));
        [Ym,Yt,Ybg,WMth] = APPmini(obj,VF); %#ok<ASGLU>
        %[Ym,Yt,Ybg,WMth] = cat_run_job_APP_init1070(single(obj.image.private.dat(:,:,:)),vx_vol,job.extopts.verb); %#ok<ASGLU>
        if ~debug, clear Yt; end
      end
          
      
      
      %% Lesion masking as zero values of the orignal image (2018-06):
      %  We do not use NaN and -INF because (i) most images are only (u)int16
      %  and do not allow such values, (ii) NaN can be part of the background
      %  of resliced images, and (iii) multiple options are not required here. 
      %  Zero values can also occure by poor data scaling or processing in the 
      %  background but also by other (large) CSF regions and we have to remove  
      %  these regions later. 
      %  We further discussed to use a separate mask images but finally decided
      %  to keep this as simple as possible using no additional options!
      %  Moreover, we have to test here anyway to create warnings in case
      %  of inoptimal settings (e.g. no SLC but possible large lesions).
      obj.image0 = spm_vol(job.channel(1).vols0{subj});
      Ysrc0      = spm_read_vols(obj.image0); 
      Ylesion    = single(Ysrc0==0 | isnan(Ysrc0) | isinf(Ysrc0)); clear Ysrc0; 
      Ylesion(smooth3(Ylesion)<0.5)=0; % general denoising 
      if any( obj.image0.dim ~= obj.image.dim )
        mat      = obj.image0.mat \ obj.image.mat;
        Ylesion  = smooth3(Ylesion); 
        Ylesionr = zeros(obj.image.dim,'single'); 
        for i=1:obj.image.dim(3)
          Ylesionr(:,:,i) = single(spm_slice_vol(Ylesion,mat*spm_matrix([0 0 i]),obj.image.dim(1:2),[1,NaN]));
        end
        Ylesion = Ylesionr>0.5; clear Ylesionr;
      end
      % use brainmask
      VFa = VF; VFa.mat = Affine * VF.mat; %Fa.mat = res0(2).Affine * VF.mat;
      if isfield(VFa,'dat'), VFa = rmfield(VFa,'dat'); end
      [Vmsk,Yb] = cat_vol_imcalc([VFa,spm_vol(Pb)],Pbt,'i2',struct('interp',3,'verb',0,'mask',-1)); clear Vmsk;  %#ok<ASGLU>
      Ylesion = Ylesion & ~cat_vol_morph(Yb<0.9,'dd',5); clear Yb; 
      % check settings 
      % RD202105: in primates the data, template and affreg is often inoptimal so we skip this test  
      if sum(Ylesion(:))/1000 > 1 && strcmp('human',job.extopts.species)
        fprintf('%5.0fs\n',etime(clock,stime)); stime = []; 
        if ~job.extopts.SLC
          % this could be critical and we use a warning for >1 cm3 and an alert in case of >10 cm3
          cat_io_addwarning([mfilename ':StrokeLesionButNoCorrection'],sprintf( ...
           ['There are %0.2f cm%s of zeros within the brain but Stroke Lesion \\\\n', ...
            'Correction (SLC) inactive (available in the expert mode). '], ...
            sum(Ylesion(:))/1000,native2unicode(179, 'latin1')),1 + (sum(Ylesion(:))/1000 > 10),[0 1]);  
          clear Ylesion; 
        else
          cat_io_cprintf('note',sprintf('SLC: Found masked region of %0.2f cm%s. \n', sum(Ylesion(:))/1000,native2unicode(179, 'latin1'))); 
        end
      end
      
      %% APP for spm_maff8
      %  optimize intensity range
      %  we have to rewrite the image, because SPM reads it again 
      if job.extopts.APP>0
          % WM threshold
          Ysrc = single(obj.image.private.dat(:,:,:)); 
          Ysrc(isnan(Ysrc) | isinf(Ysrc)) = min(Ysrc(:));

          if job.extopts.APP==1070 || job.extopts.APP==1144 
            % APPinit is just a simple bias correction for affreg and should
            % not be used further although it maybe helps in some cases!
            Ymc = Ysrc; 
          else
            bth = min( [ mean(single(Ysrc( Ybg(:)))) - 2*std(single(Ysrc( Ybg(:)))) , ...
                         mean(single(Ysrc(~Ybg(:)))) - 4*std(single(Ysrc(~Ybg(:)))) , ...
                         min(single(Ysrc(~Ybg(:)))) ]); 
            % use bias corrected image with original intensities 
            Ymc = Ym * abs(diff([bth,WMth])) + bth; 
            clear bth
          end
          
          % set variable and write image
          obj.image.dat(:,:,:)         = Ymc;  
          obj.image.pinfo              = repmat([255;0],1,size(Ymc,3));
          obj.image.private.dat(:,:,:) = Ymc; 

          obj.image.dt    = [spm_type('FLOAT32') spm_platform('bigend')];
          obj.image.pinfo = repmat([1;0],1,size(Ymc,3));

          % mask the background
          % RD20211229: Masking of unwanted regions in differnt cases.
          % In standard cross-secitonal processing the background of the 
          % SPM TPM is quite smooth and a hard masking works (pre R1900). 
          % However, in longitudinal TPMs with its hard background setting 
          % the SPM US have enough values. Hence, we have to include some 
          % (random) values close (10-15 mm) to the brain ("corona").  
          % It would be also possible to test the smoothness of the TPM 
          % backgroud class to avoid problems with "own" hard TPMs.
          isSPMtpm = strcmp(job.opts.tpm , fullfile(spm('dir'),'tpm','TPM.nii') ) && strcmp(job.extopts.species,'human'); 
          if exist('Ybg','var') && job.extopts.setCOM ~= 120 % setCOM == 120 - useCOM,useMaffreg,noMask
            if ~isempty(job.useprior) || job.extopts.new_release
              % new minimal masking approach in longitidunal processing to avoid backgound peak erros and for future releases 
              Ymsk        = cat_vol_morph( ~Ybg ,'dd',10,vx_vol) & ...          % remove voxels far from head
                              ~( Ybg & rand(size(Ybg))>0.5) & ...               % have a noisy corona
                              ~( cat_vol_grad( Ysrc , vx_vol)==0  &  Ysrc==0 ); % remove voxel that are 0 and have no gradient
            else
              % RD20220103: old cross-sectional setting with small correction for own TPMs
              if isSPMtpm
                Ymsk      = ~Ybg; % old default - mask background
              else
                cat_io_addwarning([mfilename ':noSPMTPM-noBGmasking'],...
                  'Different TPM detected - deactivated background masking!',1,[1 2]);  
                Ymsk      = []; % new special case for other TPMs
              end
            end
            if ~isempty( Ymsk )
              obj.msk       = VF; 
              obj.msk.pinfo = repmat([255;0],1,size(Ybg,3));
              obj.msk.dt    = [spm_type('uint8') spm_platform('bigend')];
              obj.msk.dat   = uint8( Ymsk ); 
              obj.msk       = spm_smoothto8bit(obj.msk,0.1); 
            end
          end 
          clear Ysrc Ymsk; 
      end

      

      
      %% Fine affine Registration with automatic selection in case of multiple TPMs. 
      %  This may not work for non human data (or very small brains).
      %  This part should be an external (coop?) function?
      if useprior
        stime = cat_io_cmd('SPM preprocessing 1 (estimate 1 - use prior):','','',1,stime); 
      elseif job.extopts.setCOM == 10 % no maffreg
        stime = cat_io_cmd('SPM preprocessing 1 (estimate 1 - use no TPM registration):','','',1,stime); 
      else
        stime = cat_io_cmd('SPM preprocessing 1 (estimate 1 - TPM registration):','','',1,stime); 
      end
      if ~isempty(job.opts.affreg) && ~useprior && job.extopts.setCOM ~= 10 % setcom == 10 - never use ... && strcmp('human',job.extopts.species)
        if numel(job.opts.tpm)>1
          %% merging of multiple TPMs
          obj2 = obj; obj2.image.dat(:,:,:) = max(0.0,Ym);
          [Affine,obj.tpm,res0] = cat_run_job_multiTPM(job,obj2,Affine,ppe.affreg.skullstripped,1); %#ok<ASGLU>
          Affine3 = Affine; 
        elseif 1 %if strcmp(job.extopts.species,'human')
          %% only one TPM (old approach); 
          spm_plot_convergence('Init','Fine affine registration','Mean squared difference','Iteration');
          warning off 
          
          try
            Affine2 = spm_maff8(obj.image(1),obj.samp,(obj.fwhm+1)*16,obj.tpm,Affine ,job.opts.affreg,80);
          catch
            Affine2 = spm_maff8(obj.image(1),obj.samp,(obj.fwhm+1)*16,obj.tpm,Affine ,job.opts.affreg);
          end
          scl1 = abs(det(Affine1(1:3,1:3)));
          scl2 = abs(det(Affine2(1:3,1:3)));

          if any(any(isnan(Affine2(1:3,:))))
            try
              Affine2 = spm_maff8(obj.image(1),obj.samp,(obj.fwhm+1)*4,obj.tpm,Affine ,job.opts.affreg,80);
            catch
              Affine2 = spm_maff8(obj.image(1),obj.samp,(obj.fwhm+1)*4,obj.tpm,Affine ,job.opts.affreg);
            end
            if any(any(isnan(Affine2(1:3,:)))) 
              Affine2 = Affine; 
            end
          else
            % check for > 10% larger scaling 
            if scl1 > 1.1*scl2 && job.extopts.setCOM ~= 11 % setcom == 11 - use always 
              stime = cat_io_cmd('  Use previous affine registration.','warn','',1,stime);
              %fprintf('\n  First fine affine registration failed.\n  Use affine registration from previous step.                ');
              Affine2 = Affine1;
              scl2 = scl1;
            end
          end
          try
            Affine3 = spm_maff8(obj.image(1),obj.samp,obj.fwhm,obj.tpm,Affine2,job.opts.affreg,80);
          catch
            Affine3 = spm_maff8(obj.image(1),obj.samp,obj.fwhm,obj.tpm,Affine2,job.opts.affreg);
          end

          if ~any(any(isnan(Affine3(1:3,:))))
            scl3 = abs(det(Affine3(1:3,1:3)));
            % check for > 5% larger scaling 
            if scl2 > 1.05*scl3 && job.extopts.setCOM ~= 11 % setcom == 11 - use always
              stime = cat_io_cmd('  Use previous affine registration.','warn','',1,stime);
              %fprintf('\n  Final fine affine registration failed.\n  Use fine affine registration from previous step.                ');
              Affine = Affine2;
            else
              Affine = Affine3;
            end
          else % Affine3 failed, use Affine2
            Affine = Affine2;
          end
          warning on  
        else
          Affine2 = Affine1; 
          Affine3 = Affine1; 
        end
        
        %% test for flipping 
        %fliptest = 2; 
        %[ppe.affreg.flipped, ppe.affreg.flippedval,stime] = cat_vol_testflipping(obj,Affine,fliptest,stime);
        
        if 0
          %% visual control for development and debugging
          VFa = VF; VFa.mat = AffineMod * VF.mat; %Fa.mat = res0(2).Affine * VF.mat;
          if isfield(VFa,'dat'), VFa = rmfield(VFa,'dat'); end
          [Vmsk,Yb] = cat_vol_imcalc([VFa,spm_vol(Pb)],Pbt,'i2',struct('interp',3,'verb',0,'mask',-1));  
          %[Vmsk,Yb] = cat_vol_imcalc([VFa;obj.tpm.V(1:3)],Pbt,'i2 + i3 + i4',struct('interp',3,'verb',0));  
          %[Vmsk,Yb] = cat_vol_imcalc([VFa;obj.tpm.V(5)],Pbt,'i2',struct('interp',3,'verb',0));  
          ds('d2sm','',1,Ym,Ym*0.5 + 0.5*Ym.*(Yb>0.5),round(size(Yb,3)*0.6))
        end
        
       
        if isfield(ppe.affreg,'skullstripped') && ~ppe.affreg.skullstripped 
          %% affreg with brainmask
          if debug 
            [Affine,Ybi,Ymi,Ym0] = cat_run_job_APRGs(Ym,Ybg,VF,Pb,Pbt,Affine,vx_vol,obj,job); %#ok<ASGLU>
          else
            [Affine,Ybi] = cat_run_job_APRGs(Ym,Ybg,VF,Pb,Pbt,Affine,vx_vol,obj,job);
          end
        end
    
        if ppe.affreg.skullstripped || job.extopts.gcutstr<0
          %% update number of SPM gaussian classes 
          Ybg = 1 - spm_read_vols(obj.tpm.V(1)) - spm_read_vols(obj.tpm.V(2)) - spm_read_vols(obj.tpm.V(3));
          noCSF = job.extopts.gcutstr == -2; 
          if 1
            for k=1:3 - noCSF
              obj.tpm.dat{k}     = spm_read_vols(obj.tpm.V(k));
              obj.tpm.V(k).dt(1) = 64;
              obj.tpm.V(k).dat   = double(obj.tpm.dat{k});
              obj.tpm.V(k).pinfo = repmat([1;0],1,size(Ybg,3));
            end
          end

          obj.tpm.V(4 - noCSF).dat = Ybg;
          obj.tpm.dat{4 - noCSF}   = Ybg; 
          obj.tpm.V(4 - noCSF).pinfo = repmat([1;0],1,size(Ybg,3));
          obj.tpm.V(4 - noCSF).dt(1) = 64;
          obj.tpm.dat(5 - noCSF:6) = []; 
          obj.tpm.V(5 - noCSF:6)   = []; 
          obj.tpm.bg1(4 - noCSF)   = obj.tpm.bg1(6);
          obj.tpm.bg2(4 - noCSF)   = obj.tpm.bg1(6);
          obj.tpm.bg1(5 - noCSF:6) = [];
          obj.tpm.bg2(5 - noCSF:6) = [];
          %obj.tpm.V = rmfield(obj.tpm.V,'private');
          
          % tryed 3 peaks per class, but BG detection error require manual 
          % correction (set 0) that is simple with only one class  
          if noCSF
            job.opts.ngaus = ([job.tissue(1:3).ngaus])'; 
          else
            job.opts.ngaus = [([job.tissue(1:3).ngaus])';1]; % 3*ones(4,1);1; 
          end
          obj.lkp        = [];
          for k=1:numel(job.opts.ngaus)
            job.tissue(k).ngaus = job.opts.ngaus(k);
            obj.lkp = [obj.lkp ones(1,job.tissue(k).ngaus)*k];
          end
        end
      end
      
      % adpation parameter for affine registration? 0.98 and 1.02?
      if isfield(job.extopts,'affmod') && any(job.extopts.affmod)
        AffineUnmod = Affine; 
        if numel(job.extopts.affmod)>6, job.extopts.affmod = job.extopts.affmod(1:6); end % remove too many
        if numel(job.extopts.affmod)<3, job.extopts.affmod(end+1:3) = job.extopts.affmod(1); end % isotropic
        if numel(job.extopts.affmod)<6, job.extopts.affmod(end+1:6) = 0; end % add translation
        fprintf('\n  Modify affine regitration (S=[%+3d%+3d%+3d], T=[%+3d%+3d%+3d])',job.extopts.affmod);
        sf   = (100 - job.extopts.affmod(1:3)) / 100;  
        imat = spm_imatrix(Affine); 
        COMc = [eye(4,3), [ 0; -24 / mean(imat(7:9)); -12 / mean(imat(7:9)); 1]  ]; 
        imat = spm_imatrix(Affine * COMc); 
        imat(1:3) = imat(1:3) - job.extopts.affmod(4:6); 
        imat(7:9) = imat(7:9) .* sf;  
        AffineMod = spm_matrix(imat) / COMc; 
        
        res.AffineUnmod = AffineUnmod; 
        res.AffineMod   = AffineMod;
      else
        AffineMod = Affine;
      end 
      obj.Affine  = AffineMod;
      cat_err_res.obj = obj; 
    
      
      %% SPM preprocessing 1
      %  ds('l2','a',0.5,Ym,Ybg,Ym,Ym,140);
      %  ds('l2','a',0.5,Ysrc/WMth,Yb,Ysrc/WMth,Yb,140);
      warning off 
      try 
        %% inital estimate
        stime = cat_io_cmd('SPM preprocessing 1 (estimate 2):','','',job.extopts.verb-1,stime);
        obj.tol = job.opts.tol; 
        
        % RD202012:  Missclassification of GM as CSF and BG as tissue:
        %  We observed problems with high-quality data (e.g. AVGs) and
        %  interpolated low resolution data (single_subT1=Collins), 
        %  where (low-intensity) GM was missclassified as CSF but also 
        %  miss-classification of background. The problems where caused
        %  by the US (or better the way we use it here) and higher
        %  accuracy (increased number of minimum iterations in 
        %  cat_spm_preproc8) was essential. Nevertheless, some
        %  cases still cause severe errors at 3 mm sample size but 
        %  not for other resolutions (eg. 4, 6, 2 mm). In addition, the
        %  log-likelihood became NaN in such cases. Hence, I added a 
        %  little loop her to test other resolutions for samp. We keep
        %  the output here quit simple to avoid confusion. samp is a
        %  rarely used expert parameter and other resolutions are only 
        %  used as backup and the effects should be not too strong for 
        %  normal data without strong bias. 

        % sampling resolution definition
        if      round(obj.samp) == 3, samp = [obj.samp 4 2]; 
        elseif  round(obj.samp) == 2, samp = [obj.samp 3 4];
        elseif  ~strcmp(job.extopts.species,'human')
                                      samp = [obj.samp obj.samp*2 obj.samp/2];
        else,                         samp = [obj.samp 3 2]; 
        end 

        if job.opts.redspmres 
          image1 = obj.image; 
          [obj.image,redspmres]  = cat_vol_resize(obj.image,'interpv',1);
        end
        
        % run loop until you get a non NaN
        % #### additional threshold is maybe also helpful ####
        warning off; % turn off "Warning: Using 'state' to set RANDN's internal state causes RAND ..."
        for sampi = 1:numel(samp)
          obj.samp = samp(sampi); 
          try 
            res = cat_spm_preproc8(obj);
            if any(~isnan(res.ll))
              break
            else
              stime = cat_io_cmd(sprintf('SPM preprocessing 1 (estimate %d):',...
                2 + sampi),'caution','',job.extopts.verb-1,stime);
            end
          catch  
            % RD202110: Catch real errors of cat_spm_preproc8 and try a 
            %           skull-stripped version just to get some result.
            stime = cat_io_cmd(sprintf('SPM preprocessing 1 (estimate %d skull-stripped):',...
                2 + sampi),'caution','',job.extopts.verb-1,stime);
            if exist('Ybi','var') % use individual mask
              obj.image.dat = obj.image.dat .* (cat_vbdist(single(Ybi>0.5))<10);
            else % use template mask 
              VFa = VF; VFa.mat = Affine * VF.mat; %Fa.mat = res0(2).Affine * VF.mat;
              if isfield(VFa,'dat'), VFa = rmfield(VFa,'dat'); end
              [Vmsk,Yb] = cat_vol_imcalc([VFa,spm_vol(Pb)],Pbt,'i2',struct('interp',3,'verb',0,'mask',-1));  
              ds('d2sm','',1,Ym,Ym.*(Yb>0.5),round(size(Yb,3)*0.6))
              obj.image.dat = obj.image.dat .* (cat_vbdist(single(Yb>0.5))<10);
            end 
            res = cat_spm_preproc8(obj);
            if any(~isnan(res.ll))
              break
            else
              stime = cat_io_cmd(sprintf('SPM preprocessing 1 (estimate %d):',...
                2 + sampi),'caution','',job.extopts.verb-1,stime);
            end
          end
        end
        if ~exist('res','var')
          cat_io_printf('SPM preprocessing with default settings failed. Run backup settings. \n'); 
        end
        warning on; 
          
        if job.opts.redspmres
          res.image1 = image1; 
          clear reduce; 
        end
        
        % unknown BG detection problems in INDI_NHa > manual setting
        if ppe.affreg.skullstripped, res.mn(end) = 0; end 

      catch
        %%
        cat_io_addwarning([mfilename ':ignoreErrors'],'Run backup function (IN DEVELOPMENT).',1,[1 1]); 
        
        if isfield(obj.image,'dat')
          tmp = obj.image.dat; 
        else
          tmp = spm_read_vols(obj.image); 
          dt2 = obj.image.dt(1); 
          dts = cat_io_strrep(spm_type(dt2),{'float32','float64'},{'single','double'}); 
          obj.image.dat = eval(sprintf('%s(tmp);',dts)); 
          obj.image.pinfo = repmat([1;0],1,size(tmp,3));
        end
        if exist('Ybi','var')
          obj.image.dat = obj.image.dat .* (cat_vbdist(single(Ybi>0.5))<10);
        else
          VFa = VF; VFa.mat = Affine * VF.mat; %Fa.mat = res0(2).Affine * VF.mat;
          if isfield(VFa,'dat'), VFa = rmfield(VFa,'dat'); end
          [Vmsk,Yb] = cat_vol_imcalc([VFa,spm_vol(Pb)],Pbt,'i2',struct('interp',3,'verb',0,'mask',-1));  
          ds('d2sm','',1,Ym,Ym.*(Yb>0.5),round(size(Yb,3)*0.6))
          obj.image.dat = obj.image.dat .* (cat_vbdist(single(Yb>0.5))<10);
        end
        
        suc = 0;
        % try higher accuracy
        while obj.tol>10e-9 && suc == 0
          obj.tol = obj.tol / 10;
          try
            res = cat_spm_preproc8(obj);
            suc = 1;
          end
        end
        if suc == 0  
          % try lower accuracy
          while obj.tol<1 && suc == 0
            obj.tol = obj.tol * 10;
            try
              res = cat_spm_preproc8(obj);
              suc = 1;
            end
          end
        end
        
        if any( (vx_vol ~= vx_voli) ) || ~strcmp(job.extopts.species,'human')
          [pp,ff,ee] = spm_fileparts(job.channel(1).vols{subj});
          delete(fullfile(pp,[ff,ee]));
        end
        
        if suc==0
          %%
          mati = spm_imatrix(V.mat);
          
          error('cat_run_job:spm_preproc8',sprintf([
            'Error in spm_preproc8. Check image and orientation. \n'...
            '  Volume size (x,y,z):   %8.0f %8.0f %8.0f \n' ...
            '  Origin (x,y,z):        %8.1f %8.1f %8.1f \n' ...
            '  Rotation (deg):        %8.1f %8.1f %8.1f \n' ...
            '  Resolution:            %8.1f %8.1f %8.1f \n'],...
            V.dim,[mati(1:3),mati(4:6),mati(7:9),]));
        end
        
        %% set internal image
        if ~exist('dt2','var')
          %tmp = obj.image.dat;
          dt2 = obj.image.dt(1); 
          dts = cat_io_strrep(spm_type(dt2),{'float32','float64'},{'single','double'}); 
        end
        obj.image.dat   = eval(sprintf('%s(tmp);',dts)); 
        obj.image.pinfo = repmat([1;0],1,size(tmp,3));
        obj.image.dt(1) = dt2; 
        res.image.dat   = eval(sprintf('%s(tmp);',dts)); 
        res.image.pinfo = repmat([1;0],1,size(tmp,3));
        res.image.dt(1) = dt2; 
      end
      warning on 
      
      if job.extopts.expertgui>1
        %% print the tissue peaks 
        mnstr = sprintf('\n  SPM-US:  ll=%0.6f, Tissue-peaks: ',res.ll);
        for lkpi = 1:numel(res.lkp)
          if lkpi==1 || ( res.lkp(lkpi) ~= res.lkp(lkpi-1) )
            mnstr = sprintf('%s  (%d) ',mnstr,res.lkp(lkpi)); 
          end
          if lkpi>1 &&( res.lkp(lkpi) == res.lkp(lkpi-1) ), mnstr = sprintf('%s, ',mnstr); end
          if sum(res.lkp == res.lkp(lkpi))>1 && res.mg(lkpi)==max( res.mg( res.lkp == res.lkp(lkpi) )), mnstr = sprintf('%s*',mnstr); end
          mnstr = sprintf('%s%0.2f',mnstr,res.mn( lkpi )); 
        end
        cat_io_cprintf('blue',sprintf('%s\n',mnstr)); 
        cat_io_cmd(' ',' ');
      end
      fprintf('%5.0fs\n',etime(clock,stime));   

      %% check contrast (and convergence)
      %min(1,max(0,1 - sum( shiftdim(res.vr) ./ res.mn' .* res.mg ./ mean(res.mn(res.lkp(2))) ) ));  
        
      clsint = @(x) round( sum(res.mn(res.lkp==x) .* res.mg(res.lkp==x)') * 10^5)/10^5;
      Tgw = [cat_stat_nanmean(res.mn(res.lkp==1)) cat_stat_nanmean(res.mn(res.lkp==2))]; 
      Tth = [
        ... min(res.mn(res.lkp==6 & res.mg'>0.3)) ... % bg; ignore the background, because of MP2RGAGE, R1, and MT weighted images  
        max( min( clsint(3) ,  max(Tgw)+1.5*abs(diff(Tgw))) , min(Tgw)-1.5*abs(diff(Tgw)) ) ... % csf with limit for T2!
        clsint(1) ... gm
        clsint(2) ... wm 
        clsint(4) ... skull
        clsint(5) ... head tissue
        clsint(6) ... background
      ];
      
      res.Tth = Tth; 
      cat_err_res.res = res;   
      
      % inactive preprocessing of inverse images (PD/T2) ... just try
      if 0 % any(diff(Tth(1:3))<=0)
        error('cat_run_job:BadImageProperties', ...
        ['CAT12 is designed to work only on highres T1 images. \n' ...
         'T2/PD preprocessing can be forced on your own risk by setting \n' ...
         '''cat12.extopts.INV=1'' in the cat_default file. If this was a highres \n' ...
         'T1 image then the initial segmentation might be failed, probably \n' ...
         'because of alignment problems (please check image orientation).']);    
      end
      
      % RD202006: Throw warning/error?
      % Due to inaccuracies of the clsint function it is better to print 
      % this as intense warning.
      if any( Tth(2:3)<0 ) || job.test_warnings
        cat_io_addwarning([mfilename ':negVal'],sprintf( ...
         ['CAT12 was developed for images with positive values and \\\\n', ...
          'negative values can lead to preprocessing problems. The average \\\\n', ...
          'intensities of CSF/GM/WM are %0.4f/%0.4f/%0.4f. \\\\n', ...
          'If you observe problems, you can use the %s to scale your data.'], Tth(1:3), ...
          spm_file('Datatype-batch','link','spm_jobman(''interactive'','''',''spm.tools.cat.tools.spmtype'');')),2,[0 1],Tth);
      end
  end
  
  %% updated tpm information for skull-stripped data should be available for cat_main
  if isfield(obj.tpm,'bg1') && exist('ppe','var') && ( ppe.affreg.skullstripped || job.extopts.gcutstr<0 )
    fname = res.tpm(1).fname;
    res.tpm       = obj.tpm;
    res.tpm(1).fname = fname;
  end
  spm_progress_bar('Clear');
          
  % call main processing
  res.tpm     = obj.tpm.V;
  res.stime   = stime0;
  res.catlog  = catlog; 
  res.Affine0 = res.Affine; 
  res.image0  = spm_vol(job.channel(1).vols0{subj}); 
  if exist('ppe','var'), res.ppe = ppe; end
  
  if isfield(job.extopts,'affmod') && any(job.extopts.affmod)
    res.AffineUnmod = AffineUnmod; 
    res.AffineMod   = AffineMod;
  end
  
  if exist('Ylesion','var'), res.Ylesion = Ylesion; else res.Ylesion = false(size(res.image.dim)); end; clear Ylesion;
  if exist('redspmres','var'); res.redspmres = redspmres; res.image1 = image1; end
  job.subj = subj; 
  
  % call new pipeline in case of inverse images (PD/T2/FLAIR) 
  if exist('Tth','var') && ( any(diff(Tth(1:3))<=0) || job.test_warnings )
    % update SPM processing ?
    % use stronger bias correction ?
    
    % warning/error because we use the new pipeline 
    job.extopts.inv_weighting = 1; 
    cat_io_addwarning([mfilename ':nonT1contrast'],sprintf( ...
      ['A non-T1 contrast was detected and the NEW EXPERIMENTAL \\\\n' ...
       'PIPELINE is used! If it is was a high-resolution T1 image, \\\\n' ...
       'the initial segmentation might have failed, probably due \\\\n' ...
       'to alignment problems (please check the image orientation).']),1,[0 1],Tth);
    cat_main(res,obj.tpm,job);
  else
    cat_main1639(res,obj.tpm,job);
  end
  
  % delete denoised/interpolated image
  [pp,ff,ee] = spm_fileparts(job.channel(1).vols{subj});
  if exist(fullfile(pp,[ff,ee]),'file') 
    delete(fullfile(pp,[ff,ee]));
  end
  %%
  
  if usediary
    diary off;
  end
return

%=======================================================================
function [Ym,Yt,Ybg,WMth] = APPmini(obj,VF)
%% very simple affine preprocessing (APP)
%  ------------------------------------------------------------------------
%  Creates an intensity normalized image Ym by the average higher tissue
%  intensity WMth estimated in the mask Yt. Moreover, it estimates the
%  background region Ybg. 
%
%  [Ym,Yt,Ybg,WMth] = APPmini(obj,VF)
%  ------------------------------------------------------------------------
%  Robert Dahnke 2019/01

  Ysrc = single(obj.image.private.dat(:,:,:)); 

  % remove outlier and use SPM for intensity normalization to uint8 
  % empirical division by 200 to get WM values around 1.0
  Ysrc = cat_stat_histth(Ysrc,99.9);
  VF0  = cat_spm_smoothto8bit(VF,0.1); 
  Ym   = single(VF0.dat)/200; clear VG0 
  
  % find the larges object and estimate the averag intensity
  % keep in mind that this will may inlcude the head (and in MP2RAGE/MT/R1
  % images also the background), i.e. highest intensity is may head,
  % blood vessels or WM or CSF in T1/PD
  Yt   = cat_vol_morph(Ym>cat_stat_nanmean(Ym(Ym(:)>0.1)),'l',[100 1000])>0.5;
  WMth = cat_stat_kmeans( Ysrc(Yt(:)) , 1); 
  
  % rescale Ym and roughly estimate the background (not in MP2Rage/MT/R1)
  Ym   = Ysrc ./ WMth;
  Ybg  = cat_vol_morph(Ym<0.2,'l',[100 1000])>0;
  
return

function APP_RMSE = checkAPP(Ym,Ysrc) 
%% check Ym
%  ------------------------------------------------------------------------
%  Function to compare the normalized gradient maps of two input images 
%  that should be nearly identical.
%
%  APP_RMSE = checkAPP(Ym,Ysrc) 
%  ------------------------------------------------------------------------
%  Robert Dahnke 2019/01

  % remove strongest outlier
  Ym   = cat_stat_histth(Ym,99.9);
  Ysrc = cat_stat_histth(Ysrc,99.9);

  % avoid division by zeros
  Ym   = Ym   + min(Ym(:));
  Ysrc = Ysrc + min(Ysrc(:)); 
  
  % normalized gradient maps
  Ygm = cat_vol_grad(Ym)   ./ (Ym + eps);     
  Ygs = cat_vol_grad(Ysrc) ./ (Ysrc + eps);

  % use only the central region and values in the expected tissue range
  sYm  = round(size(Ym) / 5);
  Ymsk = false(size(Ym) ); Ymsk(sYm(1):end-sYm(1),sYm(2):end-sYm(2),sYm(3):end-sYm(3)) = true;  
  Ymsk = Ymsk & cat_vol_morph(Ygm<2 & Ygs<2 & Ym>0.5 & Ysrc>0.5,'e');
  
  % general error between both images within the mask
  APP_RMSE = cat_stat_nanmean( ( Ygm(Ymsk(:)) - Ygs(Ymsk(:)) ).^2 )^0.5;
  
return
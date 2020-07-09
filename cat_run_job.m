function cat_run_job(job,tpm,subj)
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
% Christian Gaser
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
    if job.extopts.subfolders
    
      folders = {'mri','report'};
      warning('off', 'MATLAB:MKDIR:DirectoryExists');
      for i=1:numel(folders)
        if ~exist(fullfile(pth,folders{i}),'dir')
          mkdir(fullfile(pth,folders{i}));
        end
      end
    
      if ~exist(fullfile(pth,'surf'),'dir') && job.output.surface
        mkdir(fullfile(pth,'surf'));
      end
    
      if ~exist(fullfile(pth,'label'),'dir') && job.output.ROI
        mkdir(fullfile(pth,'label'));
      end
      
      mrifolder    = 'mri';
      reportfolder = 'report';
    else
      mrifolder    = '';
      reportfolder = '';
    end
    
    % create subject-wise diary file with the command-line output
    [pp,ff,ee,ex] = spm_fileparts(job.data{subj});  %#ok<ASGLU>
    catlog = fullfile(pth,reportfolder,['catlog_' ff '.txt']);
    if exist(catlog,'file'), delete(catlog); end % write every time a new file, turn this of to have an additional log file
    
    % check if not another diary is already written that is not the default- or catlog-file. 
    if ~strcmpi(spm_check_version,'octave')
      olddiary = spm_str_manip( get(0,'DiaryFile') , 't');
      usediary = ~isempty(strfind( olddiary , 'diary' )) | ~isempty(strfind( olddiary , 'catlog_' )); 
      if usediary
        diary(catlog); 
        diary on; 
      else  
        cat_io_cprintf('warn',sprintf('External diary log is writen to "%s".\n',get(0,'DiaryFile'))); 
      end
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
       exist(fullfile(pp,[ff(3:end) '_seg8.mat']),'file')
       
        job.data{subj}          = fullfile(pp,[ff ee]); 
        job.channel.vols{subj}  = fullfile(pp,[ff ee]); 

        % prepare SPM preprocessing structure 
        images = job.channel(1).vols{subj};
        for n=2:numel(job.channel)
          images = char(images,job.channel(n).vols{subj});
        end

        [cv,cr]      = cat_version;
        obj.image    = spm_vol(images);
        obj.image.descrip = sprintf('%sR%s < %s',cv,cr,obj.image.descrip); % add CAT version 
        obj.fwhm     = job.opts.fwhm;
        obj.biasreg  = cat(1,job.opts.biasreg);
        obj.biasfwhm = cat(1,job.opts.biasfwhm);
        obj.tol      = job.opts.tol;
        obj.tpm      = tpm;
        obj.lkp      = [];
        spm_check_orientations(obj.image);
        
        if all(isfinite(cat(1,job.tissue.ngaus)))
            for k=1:numel(job.tissue)
                obj.lkp = [obj.lkp ones(1,job.tissue(k).ngaus)*k];
            end
        end
        
        obj.reg      = job.opts.warpreg;
        obj.samp     = job.opts.samp;              
        cfname  = fullfile(pp,[ff ee]);
        ofname  = fullfile(pp,[ff(3:end) ee]); 
        nfname  = fullfile(pp,mrifolder,['n' ff '.nii']); 
        copyfile(ofname,nfname); 

        Ysrc0    = single(spm_read_vols(obj.image)); 
        Ylesion  = single(isnan(Ysrc0) | isinf(Ysrc0) | Ysrc0==0); clear Ysrc0;
        
        res = load(fullfile(pp,[ff(3:end) '_seg8.mat']));
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
        %  -----------------------------------------------------------------
        for n=1:numel(job.channel) 
          V = spm_vol(job.channel(n).vols{subj});
          vx_vol = sqrt(sum(V.mat(1:3,1:3).^2));

          % maximum [ slice-thickness , volume^3 , anisotropy ]
          reslimits = [5 3 8]; 
         
          % too thin slices
          if any( vx_vol > reslimits(1) ) || job.test_warnings
            mid = 'cat_run_job:TooLowResolution'; 
            msg = sprintf(['Voxel resolution should be better than %d mm in any dimension for \\\\n' ...
              'reliable preprocessing! This image has a resolution of %0.2fx%0.2fx%0.2f mm%s. '], ... 
              reslimits(1),vx_vol,native2unicode(179, 'latin1'));
            if job.extopts.ignoreErrors < 2
              error(mid,msg); %#ok<SPERR>
            else
              cat_io_addwarning(mid,msg,2,[0 1],vx_vol);
            end
          end
          
          % too small voxel volume (smaller than 3x3x3 mm3)
          if prod(vx_vol) > reslimits(2)^3 || job.test_warnings
            mid = 'cat_run_job:TooHighVoxelVolume'; 
            msg = sprintf(['Voxel volume should be smaller than %d mm%s (around %dx%dx%d mm%s) for \\\\n' ...
                    'reliable preprocessing! This image has a voxel volume of %0.2f mm%s. '], ...
                    reslimits(2)^3,native2unicode(179, 'latin1'),reslimits(2),reslimits(2),reslimits(2),...
                    native2unicode(179, 'latin1'),prod(vx_vol),native2unicode(179, 'latin1'));
            if job.extopts.ignoreErrors < 2
              error(mid,msg); %#ok<SPERR>
            else
              cat_io_addwarning(mid,msg,2,[0 1],vx_vol);
            end
          end
          
          % anisotropy
          if max(vx_vol) / min(vx_vol) > reslimits(3) || job.test_warnings
            mid = 'cat_run_job:TooStrongIsotropy';
            msg = sprintf(['Voxel isotropy (max(vx_size)/min(vx_size)) should be smaller than %d for \\\\n' ...
                    'reliable preprocessing! This image has a resolution %0.2fx%0.2fx%0.2f mm%s \\\\nand a isotropy of %0.2f. '], ...
                    reslimits(3),vx_vol,native2unicode(179, 'latin1'),max(vx_vol)/min(vx_vol));
            if job.extopts.ignoreErrors < 2
              error(mid,msg);  %#ok<SPERR>
            else
              cat_io_addwarning(mid,msg,2,[0 1],vx_vol);
            end
          end
        end

        % save original file name 
        for n=1:numel(job.channel) 
          job.channel(n).vols0{subj} = job.channel(n).vols{subj};
        end
        
       
        % always create the n*.nii image because of the real masking of the
        % T1 data for spm_preproc8 that include rewriting the image!
        for n=1:numel(job.channel) 
          [pp,ff,ee] = spm_fileparts(job.channel(n).vols{subj}); 
          ofname  = fullfile(pp,[ff ee]); 
          nfname  = fullfile(pp,mrifolder,['n' ff '.nii']); 
          if strcmp(ee,'.nii')
            copyfile(ofname,nfname); 
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
          vx_vol    = round(vx_vol*10^2)/10^2; % avoid small differences 

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
           
            imat      = spm_imatrix(Vi.mat); 
            Vi.dim    = round(Vi.dim .* vx_vol./vx_voli);
            imat(7:9) = vx_voli .* sign(imat(7:9));
            Vi.mat    = spm_matrix(imat);

            Vn = spm_vol(job.channel(n).vols{subj}); 
            cat_vol_imcalc(Vn,Vi,'i1',struct('interp',2,'verb',0,'mask',-1));
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
   % 20181222       if ~strcmp(job.extopts.species,'human'), job.extopts.APP = 5; end
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
        if all(isfinite(cat(1,job.tissue.ngaus))),
          for k=1:numel(job.tissue),
            obj.lkp = [obj.lkp ones(1,job.tissue(k).ngaus)*k];
          end;
        end
        spm_check_orientations(obj.image);
        cat_err_res.obj = obj; 
        
        %% Initial affine registration.
        %  -----------------------------------------------------------------
        Affine  = eye(4);
        [pp,ff] = spm_fileparts(job.channel(1).vols{subj});
        Pbt     = fullfile(pp,mrifolder,['brainmask_' ff '.nii']);
        Pb      = char(job.extopts.brainmask);
        Pt1     = char(job.extopts.T1);
        
        
        % use affine transformation of given (average) data for longitudinal mode
        if ~isempty(job.useprior)
          priorname = job.useprior{1};
          [pp,ff,ee,ex] = spm_fileparts(priorname);  %#ok<ASGLU>
          catxml = fullfile(pp,reportfolder,['cat_' ff '.xml']);

          % check that file exists and get affine transformation
          if exist(catxml,'file') 
            cat_io_cmd(sprintf('Use affine transformation from %s',[ff ee]));
            xml             = cat_io_xml(catxml);
            Affine          = xml.parameter.spm.Affine;
            useprior        = 1; 
            job.useprior    = catxml; 
            job.opts.affreg = ''; 
          else
            cat_io_addwarning('cat_run_job:UseAffRegPrior',...
              sprintf('Affine prior file "%s" not found. \\\\nEstimate individual affine transformation. ',catxml),2,[1 2],catxml);
            useprior        = 0;
          end
          clear catxml; 
        else
          useprior          = 0;
        end
        if ~isempty(job.opts.affreg) && ~useprior 
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
              cat_io_addwarning('cat_run_job:skullStrippedInputWithSkullStripping',msg,[0 1],ppe.affreg.skullstrippedpara);
              
            elseif job.extopts.gcutstr<0 && ~ppe.affreg.skullstripped || job.test_warnings
              cat_io_addwarning('cat_run_job:noSkullStrippingButSkull',[...
                  'Skull-Stripping is deactivated but skull was detected. \\n' ...
                  'Go on without skull-stripping what possibly will fail.'],[0 1],ppe.affreg.skullstrippedpara);
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

          % APP step 1 rough bias correction 
          % --------------------------------------------------------------
          % Already for the rough initial affine registration a simple  
          % bias corrected and intensity scaled image is required, because
          % large head intensities can disturb the whole process.
          % --------------------------------------------------------------
          % ds('l2','',vx_vol,Ym, Yt + 2*Ybg,obj.image.private.dat(:,:,:)/WMth,Ym,60)
          if job.extopts.APP == 1070 || job.extopts.APP == 1144
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
              cat_io_addwarning('cat_run_job:APPerror','APP failed. Use simple scaling.',1,[0 0],apperr);
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
                cat_io_addwarning('cat_run_job:APPerror',...
                  sprintf('Detect problems in APP preprocessing (APPRMS: %0.4f). \\\\nDo not use APP results. ',APPRMS),1,[0 1],APPRMS);
              end 
            end
        
            if ~debug, clear Yt; end

            stime = cat_io_cmd('Affine registration','','',1,stime); 

            % use further data limitation and remove background for affreg 
            [Ym2,ths]     = cat_stat_histth(Ym);
            Ym2           = (Ym2 - ths(1)) ./ diff(ths) .* (1 - Ybg); 
            % write data to VF
            VF.dt         = [spm_type('UINT8') spm_platform('bigend')];
            VF.dat(:,:,:) = cat_vol_ctype(Ym2 * 200,'uint8'); 
            VF.pinfo      = repmat([1;0],1,size(Ym,3));
            clear WI ths; 

            % smoothing
            resa  = obj.samp*2; % definine smoothing by sample size
            VF1   = spm_smoothto8bit(VF,resa);
            VG1   = spm_smoothto8bit(VG,resa);

          else
            % standard approach with static resa value and no VG smoothing
            stime = cat_io_cmd('Coarse affine registration'); 
            resa  = 8;
            VF1   = spm_smoothto8bit(VF,resa);
            VG1   = VG; 
            [Ym,Yt,Ybg,WMth] = APPmini(obj,VF);
          end

          %% prepare affine parameter 
          aflags     = struct('sep',obj.samp,'regtype','subj','WG',[],'WF',[],'globnorm',1); 
          aflags.sep = max(aflags.sep,max(sqrt(sum(VG(1).mat(1:3,1:3).^2))));
          aflags.sep = max(aflags.sep,max(sqrt(sum(VF(1).mat(1:3,1:3).^2))));

          %{
          % use affine transformation of given (average) data for longitudinal mode
          if isfield(job,'useprior') && ~isempty(job.useprior)
            priorname = job.useprior{1};
            [pp,ff,ee,ex] = spm_fileparts(priorname);  %#ok<ASGLU>
            catxml = fullfile(pp,reportfolder,['cat_' ff '.xml']);
            
            % check that file exists and get affine transformation
            if exist(catxml,'file') || job.test_warnings
              cat_io_cmd(sprintf('Use affine transformation from %s',[ff ee]));
              xml = cat_io_xml(catxml);
              Affine   = xml.parameter.spm.Affine;
              affscale = 1;
              useprior = 1;
            else
              %%
              cat_io_addwarning('cat_run_job:UseAffRegPrior',...
                sprintf('Affine prior file "%s" not found. \\\\nEstimate individual affine transformation. ',catxml),2);
              useprior = 0;
            end
          else
            useprior = 0;
          end
          %}
          
          if strcmp('human',job.extopts.species) %&& ~useprior
            % affine registration
            try
              spm_plot_convergence('Init','Coarse affine registration','Mean squared difference','Iteration');
            catch
              spm_chi2_plot('Init','Coarse affine registration','Mean squared difference','Iteration');
            end

            warning off
            try 
              [Affine0, affscale]  = spm_affreg(VG1, VF1, aflags, eye(4)); Affine = Affine0; 
            catch
              affscale = 0; 
            end
            % RD202007: Unimportant information if maff8 works 
            if job.extopts.expertgui 
              if affscale>3 || affscale<0.5
                stime  = cat_io_cmd('Coarse affine registration failed. Try fine affine registration.','','',1,stime);
                Affine = eye(4); 
              end
            end
            warning on
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
          if strcmp('human',job.extopts.species) %&& ~useprior 
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
        %  We further discussed to use a separate mask images but finally desided
        %  to keep this as simple as possible using no additional options!
        obj.image0 = spm_vol(job.channel(1).vols0{subj});
        Ysrc0      = spm_read_vols(obj.image0); 
        Ylesion    = single(Ysrc0==0); clear Ysrc0; 
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
        if exist('Ybg','var'), Ylesion(Ybg)=0; end % denoising in background
        
        
        
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

            % mask the eroded background
            % RD202006: masking of distant background
            % This has strong effects for some images but I found no good 
            % explanation how to use the mask. However, it seems that it is
            % useful to mask unclear and/or bad background voxels but not 
            % all of them. So we use the eroded background segment mask and 
            % remove also regions with 0 and no gradient that are often the
            % result of defacing, skull-stripping and reslicing.
            % RD202006: SVE 32 dataset 
            % We use a noisy corona here to avoid that SPM try to fit a
            % head class into it.
            % RD202006: thickness phantom problems
            % Masking causes general problems in SPM US with Christian's 
            % thickness phantom (brain PVE voxels were aligned to class 5) 
            % that required further correction in cat_main_updateSPM. 
            if exist('Ybg','var')
              Ymsk          = cat_vol_morph( ~Ybg ,'dd',15,vx_vol) & ...        % remove voxels far from head
                              ~( Ybg & rand(size(Ybg))>0.5) & ...               % have a noisy corona
                              ~( cat_vol_grad( Ysrc , vx_vol)==0  &  Ysrc==0 ); % remove voxel that are 0 and have no gradient
              Ybge          = cat_vol_morph(Ybg,'de',10,vx_vol) | ...           % define background for cat_main_updateSPM
                              ( cat_vol_grad( Ysrc , vx_vol)==0  &  Ysrc==0 );  % RD202006: set value arbitrary to 10 mm 
              obj.msk       = VF; 
              obj.msk.pinfo = repmat([255;0],1,size(Ybg,3));
              obj.msk.dt    = [spm_type('uint8') spm_platform('bigend')];
              obj.msk.dat   = uint8( Ymsk );  % RD202006 background corona to have save background values
              obj.msk       = spm_smoothto8bit(obj.msk,0.1); 
            end            
            clear Ysrc; 
        else
          % defintion of basic variables in case of no APP 
            obj.image.dat(:,:,:)  = single(obj.image.private.dat(:,:,:));
            obj.image.dt          = [spm_type('FLOAT32') spm_platform('bigend')];
            obj.image.pinfo       = repmat([255;0],1,size(obj.image.dat,3));
            obj.msk               = VF; 
            obj.msk.pinfo         = repmat([255;0],1,size(Ybg,3));
            obj.msk.dt            = [spm_type('uint8') spm_platform('bigend')];
            obj.msk.dat           = zeros(size(obj.image.dat),'uint8');  
            obj.msk               = spm_smoothto8bit(obj.msk,0.1); 
        end

        

        
        %% Fine affine Registration with automatic selection in case of multiple TPMs. 
        %  This may not work for non human data (or very small brains).
        %  This part should be an external (coop?) function?
        stime = cat_io_cmd('SPM preprocessing 1 (estimate 1):','','',1,stime);
        if ~isempty(job.opts.affreg) && strcmp('human',job.extopts.species) 
          if ~useprior 
            % turn rand warning off
            wo = warning('QUERY','MATLAB:RandStream:ActivatingLegacyGenerators'); wo = strfind( wo.state , 'on');
            if wo, warning('OFF','MATLAB:RandStream:ActivatingLegacyGenerators'); end

            if numel(job.opts.tpm)>1
              %% merging of multiple TPMs
              obj2 = obj; obj2.image.dat(:,:,:) = max(0.0,Ym);
              [Affine,obj.tpm,res0] = cat_run_job_multiTPM(job,obj2,Affine,ppe.affreg.skullstripped,1); %#ok<ASGLU>
            elseif strcmp(job.extopts.species,'human')
              %% only one TPM (old approach); 
              spm_plot_convergence('Init','Fine affine registration','Mean squared difference','Iteration');
              
              % first we start with the given affine registration and affreg parameter (e.g. mni) and a very low resolution
              % RD202007: also here different maskings could be tested - however, it looks quite stable now 
              [Affine2,ppe.spm_maff8.ll(1)] = spm_maff8(obj.image(1),obj.samp,(obj.fwhm+1)*16,obj.tpm,Affine ,job.opts.affreg,80); 
              ppe.spm_maff8.ll_help = ['ll(1) with affreg result, ll(2) without spm_affreg init, ' ...
                'll(3) without spm_affreg and with opts.affreg=none; only test further cases if ll(i)<0.9'];   
              if ppe.spm_maff8.ll(1)<0.9
                % if there was no high overlap than we try if maff8 support better results without affreg initialization 
                [Affine2o,ppe.spm_maff8.ll(2)] = spm_maff8(obj.image(1),obj.samp,(obj.fwhm+1)*16,obj.tpm,eye(4),job.opts.affreg,80); 
                Affine2 = Affine2o; 
                if ppe.spm_maff8.ll(2)<0.9
                  % especially for very small heads the mni definiton is not good 
                  % we start here with the maff8 that is more robust to varying contrasts
                  [Affine2n,ppe.spm_maff8.ll(3)] = spm_maff8(obj.image(1),obj.samp,(obj.fwhm+1)*16,obj.tpm,eye(4),'none',80); 
                  if ppe.spm_maff8.ll(3) > ppe.spm_maff8.ll(2)
                    cat_io_addwarning('cat_run_job:spm_maff8','Use affreg=none due to better results.',1,[1 1],ppe.spm_maff8);
                    job.opts.affreg = 'none'; % in this case we have to update the affreg parameter
                    Affine2 = Affine2n;
                  else
                    if ppe.spm_maff8.ll(1) < ppe.spm_maff8.ll(2)
                      Affine2 = Affineo; 
                    end
                  end
                end
              end
              
              % after this initial step we do some refined registration with less smoothing 
              if any(any(isnan(Affine2(1:3,:))))
                Affine2 = spm_maff8(obj.image(1),obj.samp,(obj.fwhm+1)*4,obj.tpm,Affine ,job.opts.affreg,80);
                if any(any(isnan(Affine2(1:3,:)))) 
                  Affine2 = Affine; 
                end
              end
              Affine3 = spm_maff8(obj.image(1),obj.samp,obj.fwhm,obj.tpm,Affine2,job.opts.affreg,80);
              warning on  
              if ~any(any(isnan(Affine3(1:3,:)))), Affine = Affine3; else, Affine = Affine2; end
    
              % turn warning on 
              if wo, warning('ON','MATLAB:RandStream:ActivatingLegacyGenerators'); end
            end
          end
          if 0
            %% visual control for development and debugging
            VFa = VF; VFa.mat = Affine * VF.mat; %Fa.mat = res0(2).Affine * VF.mat;
            if isfield(VFa,'dat'), VFa = rmfield(VFa,'dat'); end
            [Vmsk,Yb] = cat_vol_imcalc([VFa,spm_vol(Pb)],Pbt,'i2',struct('interp',3,'verb',0,'mask',-1));  
            %[Vmsk,Yb] = cat_vol_imcalc([VFa;obj.tpm.V(1:3)],Pbt,'i2 + i3 + i4',struct('interp',3,'verb',0));  
            %[Vmsk,Yb] = cat_vol_imcalc([VFa;obj.tpm.V(5)],Pbt,'i2',struct('interp',3,'verb',0));  
            ds('d2sm','',1,Ym,Ym.*(Yb>0.5),round(size(Yb,3)*0.6))
          end
          
         
          if ~ppe.affreg.skullstripped 
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
            if 1
              for k=1:3
                obj.tpm.dat{k}     = spm_read_vols(obj.tpm.V(k));
                obj.tpm.V(k).dt(1) = 64;
                obj.tpm.V(k).dat   = double(obj.tpm.dat{k});
                obj.tpm.V(k).pinfo = repmat([1;0],1,size(Ybg,3));
              end
            end

            obj.tpm.V(4).dat = Ybg;
            obj.tpm.dat{4}   = Ybg; 
            obj.tpm.V(4).pinfo = repmat([1;0],1,size(Ybg,3));
            obj.tpm.V(4).dt(1) = 64;
            obj.tpm.dat(5:6) = []; 
            obj.tpm.V(5:6)   = []; 
            obj.tpm.bg1(4)   = obj.tpm.bg1(6);
            obj.tpm.bg2(4)   = obj.tpm.bg1(6);
            obj.tpm.bg1(5:6) = [];
            obj.tpm.bg2(5:6) = [];
            %obj.tpm.V = rmfield(obj.tpm.V,'private');
            
            % tryed 3 peaks per class, but BG detection error require manual 
            % correction (set 0) that is simple with only one class  
            job.opts.ngaus = [([job.tissue(1:3).ngaus])';1]; % 3*ones(4,1);1; 
            obj.lkp        = [];
            for k=1:numel(job.opts.ngaus)
              job.tissue(k).ngaus = job.opts.ngaus(k);
              obj.lkp = [obj.lkp ones(1,job.tissue(k).ngaus)*k];
            end
          end
        end
        
        % adpation parameter for affine registration? 0.98 and 1.02?
        %imat = spm_imatrix(Affine2); imat(7:9)=imat(7:9)*1.02; Affine2 = spm_matrix(imat); 
        
        obj.Affine = Affine;
        cat_err_res.obj = obj; 
      
        
        %% SPM preprocessing 1
        %  RD202006: Although affine registration works now quite good
        %  there are still severe problems in the SPM US in some cases 
        %  that of course also depends on our settings. Hence, we have 
        %  to try different settings, e.g. with/without masking.  In 
        %  some cases the US may run but produces a bad result that will 
        %  cause many problems in cat_main and its subfunctions, so we 
        %  should test the result, may run alternative settings and 
        %  choose the best trial. 
        %  Possible settings: 
        %    * with / without masking (seams to be quite powerful) 
        %    . optimized image (not so easy to handle)
        %      maybe use the low resolution option 
        %    . modify SPM accuracy (?)
        %    . apply skull-stripping (not optimal > last option)
        % 
        %     ds('l2','a',0.5,Ym,Ybg,Ym,Ym,140);
        %     ds('l2','a',0.5,Ysrc/WMth,Yb,Ysrc/WMth,Yb,140);

        old = 0; 
        if old
% #### keep this part temporary ... RD20200619.begin        
%{
          warning off 
          try 
            % This is the first level where we try to process the image with 
            % a mask that removes voxel far from the head and a noisy corona
            % to have some random background voxels. 
            % The mask helps in many cases to avoid bad voxels and also
            % increase processing speed but in some cases it cause also
            % severe problems. 


            % inital estimate
            stime = cat_io_cmd('SPM preprocessing 1 (estimate 2):','','',job.extopts.verb-1,stime);
            obj.tol = job.opts.tol;  
            if job.opts.redspmres==0 
              warning off; % turn off "Warning: Using 'state' to set RANDN's internal state causes RAND ..."
              res = cat_spm_preproc8(obj);
              warning on; 
            else
              % Use a low resolution version to speed up preprocessing.
              % This averages voxels and reduce noise and artifacts. 
              % RD202006: rarely used
              image1 = obj.image; 
              [obj.image,redspmres]  = cat_vol_resize(obj.image,'interpv',1);
              res = cat_spm_preproc8(obj);
              res.image1 = image1; 
              clear reduce; 
            end

            % unknown BG detection problems in INDI_NHa > manual setting
            if ppe.affreg.skullstripped, res.mn(end) = 0; end 

            % RD202006: Add test concept to catch problems and use another version. 
            % This should be a separate function that can be called from the 
            % other attempts to get a good initial segmentation. 
            % One criteria is to have a reasonable contrast between the brain
            % tissue peaks that we also used to create the error message later. 
          catch
            try 
              % This is the second level that uses no mask at all.
              stime = cat_io_cmd('SPM preprocessing 1 (estimate 3):','','',job.extopts.verb-1,stime);

              obj = rmfield(obj,'msk');  

              warning off; % turn off "Warning: Using 'state' to set RANDN's internal state causes RAND ..."
              res = cat_spm_preproc8(obj);
              warning on; 

              % RD202006: add test of the result 
            catch
              stime = cat_io_cmd('SPM preprocessing 1 (estimate 4):','','',job.extopts.verb-1,stime);

              % save the masked default image for the error report
              tmp = obj.image.dat; 
              tmp(obj.msk.dat<1) = NaN;           

              % Ignore larger parts of the image that a far from the brain mask. 
              % RD202006: masking of distant non-brain voxels 
              % Because the masking (obj.msk) was a bit unclear I previously
              % set voxel 10 mm outside the brain mask to zero. However, I now 
              % use the msk field and also increase the radius to 20 mm to have
              % more skull.
              if ~exist('Ybi','var') && exist(Pb,'file')
                % If there is no mask yet then load a affine registered brain mask 
                % version of it. 
                VFa = VF; VFa.mat = Affine * VF.mat; 
                if isfield(VFa,'dat'), VFa = rmfield(VFa,'dat'); end
                [Vmsk,Ybi] = cat_vol_imcalc([VFa,spm_vol(Pb)],Pbt,'i2',struct('interp',3,'verb',0,'mask',-1));  
                % ds('d2sm','',1,Ym,Ym.*(Ybi>0.5),round(size(Ybi,3)*0.6))
                obj.msk = cat_vbdist(single(Ybi>0.5),true(size(Ybi)),vx_vol)<20;
              end
              obj.msk.dat = cat_vbdist(single(Ybi>0.5),true(size(Ybi)),vx_vol)<20;


              % in some cases lower! accuracy worked better - I don't know why
              suc = 0;
              while obj.tol<1
                obj.tol = obj.tol * 10;
                try
                  res = cat_spm_preproc8(obj);
                  suc = 1;
                end
              end

              if any( (vx_vol ~= vx_voli) ) || ~strcmp(job.extopts.species,'human')
                [pp,ff,ee] = spm_fileparts(job.channel(1).vols{subj});
                delete(fullfile(pp,[ff,ee]));
              end


              if suc==0
                %% If it was not possible to run SPM then print an error message. 
                %  RD202006: Alternatively we could try the kmeans here but I  
                %  think it is better to give an error here and use manual
                %  selection by the job.extopts.init

                mati = spm_imatrix(V.mat);
                emsg = sprintf([
                  'Error in spm_preproc8 that possibly happens due to bad orientation or untypical contrast. \n' ...
                  '%s and/or apply %s with intensity normalization.  \n' ...
                  'If the image is a normal T1 image with good contrast and correct orientation that you \n' ...
                  'can share with us than please contact us (%s). \n' ...
                  '  Volume size (x,y,z):   %8.0f %8.0f %8.0f \n' ...
                  '  Origin (x,y,z):        %8.1f %8.1f %8.1f \n' ...
                  '  Rotation (deg):        %8.1f %8.1f %8.1f \n' ...
                  '  Resolution:            %8.1f %8.1f %8.1f \n'],...
                  spm_file('Check image orientation/contrast','link',['spm_image(''Display'', ''' ofname ''')']), ...
                  spm_file('head trimming','link','spm_jobman(''interactive'','''',''spm.tools.cat.tools.datatrimming'');'), ...
                  spm_file('vbmweb@gmail.com','link','web(''mailto:vbmweb@gmail.com'');'), ...
                  V.dim,[mati(1:3),mati(4:6),mati(7:9)]); 

            %   fprintf('%s%4.0fs - Display %s\n',nstr,etime(clock,stime),spm_file(Pgmv{gmvi},'link','cat_surf_display(''%s'')'));  

                if exist('Ybi','var') && exist('Ybg','var')
                  % estimate some thresholds with kmeans3D in specific ROIs
                  try %#ok<TRYNC> % this has to be save although a variable is miss in some cases 
                    Tcgw = kmeans3D(tmp(Ybi(:)>0.5),5); Tcgw([2,4]) = [];
                    Thd  = kmeans3D(tmp(~Ybg(:) & Ybi(:)<0.5),2);
                    Tbg  = kmeans3D(tmp(obj.msk.dat(:)>0 & Ybi(:)<0.5),1);
                    Tth  = [Tcgw Thd Tbg]; 
                    emsg = [emsg, sprintf([
                      '  Brain tissues:         %8.2f %8.2f %8.2f \n' ...
                      '  Head tissues / BG:     %8.2f %8.2f %8.2f \n'], Tth)];
                    cat_err_res.res.Tth = Tth;   
                  end 
                end

                if ~exist('res','var') || job.extopts.ignoreErrors < 2
                  error('cat_run_job:spm_preproc8',emsg);
                end
              end

              % udpate the last image version
              res.image.dat = tmp;
              obj.image.dat = tmp; 
              clear tmp;
            end
          end
          warning on 
          fprintf('%5.0fs\n',etime(clock,stime));
%}          
% #### keep this part temporary ... RD20200619.end        
        else
          %% Run SPM preprocessing with different settings. 
          % Start with the most reasonable model. 
          % I am not really sure how many levels are useful - but less than a hand
          % RD202006: Maybe we could run this with faster settings 
          % (lower samp & tol) to select the best mask definition? >> CAT12.8         
          stime  = cat_io_cmd('SPM preprocessing 1 (estimate 2):','','',job.extopts.verb-1,stime); 
          casei  = 0;     % iteration counter
          rerun  = 0;     % rerun in case of TPM problems 
          acccon = 0.33;  % acceptable contrast (optimal is 0.5, default maybe 0.35-0.45 in T1 )
          runcas = inf;     % stop for acceptable contrast (inf = test all, 2 = only casei<3) ... for test we start with inf
          if job.extopts.ignoreErrors > 2
            verbs  = 1;   % show results
          else
            verbs  = 1;   % show results
          end
          resi   = cell(1,4); modelname = cell(1,4); mincontrast = zeros(1,4); 
          
          % save the masked default image for the error report
          tmp = obj.image.dat; 
          tmp(obj.msk.dat<1) = NaN;    
          if verbs, fprintf('\n'); end
          while 1 
            casei = casei + 1; 

            switch casei
              case 1
                % use the coronal mask - nothing to do
                if ppe.affreg.skullstripped, continue, end
                obj2 = obj; 
                modelname{casei} = 'with coronal background mask'; 
              case 2
                % use now mask 
                obj2 = rmfield(obj,'msk'); 
                modelname{casei} = 'without background mask'; 
              case 3
                % use extrended brain mask
                obj2 = obj; 
                obj2.msk.dat = uint8(1 - Ybg); 
                modelname{casei} = 'with full background mask'; 
                %{
              case 4
                obj2 = obj; 
                if ~exist('Ybi','var') && exist(Pb,'file')
                  % If there is no mask yet then load a affine registered brain mask 
                  % version of it. 
                  VFa = VF; VFa.mat = Affine * VF.mat; 
                  if isfield(VFa,'dat'), VFa = rmfield(VFa,'dat'); end
                  [Vmsk,Ybi] = cat_vol_imcalc([VFa,spm_vol(Pb)],Pbt,'i2',struct('interp',3,'verb',0,'mask',-1));  
                  % ds('d2sm','',1,Ym,Ym.*(Ybi>0.5),round(size(Ybi,3)*0.6))
                  obj2.msk.dat = cat_vbdist(single(Ybi>0.5),true(size(Ybi)),vx_vol)<20;
                  modelname{casei} = 'with extended brain mask'; 
                end
                obj2.msk.dat = cat_vbdist(single(Ybi>0.5),true(size(Ybi)),vx_vol)<20;
                %}
              otherwise
                break
            end
            if verbs
              if casei == 1
                stime  = cat_io_cmd(sprintf('  Run %s',modelname{casei}),'g5',' ',job.extopts.verb-1);
              else
                stime  = cat_io_cmd(sprintf('  Run %s',modelname{casei}),'g5',' ',job.extopts.verb-1,stime); 
              end
            end
            
            % run SPM 
            warning off; % turn off "Warning: Using 'state' to set RANDN's internal state causes RAND ..."
            try
              resi{casei} = cat_spm_preproc8(obj2);
            catch
              % if something failed than define the tissue contrast parameter with the worst case 
              resi{casei}.mn  = zeros(size(obj.lkp)); 
              resi{casei}.mg  = zeros(size(obj.lkp))'; 
              resi{casei}.lkp = obj.lkp; 
              resi{casei}.ll  = 0; 
            end
            warning on; 
            
            % test result
            % - may use a weighted contrast, where a good WM-GM is more important? 
            clsint = @(resx,x) round( sum(resx.mn(resx.lkp==x) .* resx.mg(resx.lkp==x)') * 10^5)/10^5;
            clscon = @(resx,c1,c2) abs( diff([clsint(resx,c1),clsint(resx,c2)]) ) ./ ...  % specific tissue contrast between class c1 and c2
              (max([clsint(resx,1),clsint(resx,2),clsint(resx,3)]) - ...
               min([clsint(resx,1),clsint(resx,2),clsint(resx,3)])) .* ...                % signal intensity (maximum brain tissue contrast)
              (1 - any(isnan([clsint(resx,1),clsint(resx,2),clsint(resx,3)])));           % NaN as worst case
            mincontrast(casei) = min( [ clscon(resi{casei},1,2) , clscon(resi{casei},2,3) , clscon(resi{casei},1,3) ] ); 
            
            % test for TPM problems
            % this problem occured for the long TPM and no other parameter changed something 
            % RD202007 - version 0.1 - try to clean up later
            if all( [ clsint(resi{casei},1) clsint(resi{casei},2) clsint(resi{casei},3) ] == 0 ) && ~exist('rerun','var') && rerun<5
              %if strcmp( job.opts.tpm , dtpm ) % it may also help with the default SPM TPM 
              if ~exist('rerun','var') && rerun < 3
                casei = casei - 1; 
                rerun = rerun + 1; 
                cat_io_addwarning('cat_run_job:ownTPM',sprintf('Estimation problems with selected TPM. Try smoothing %d.',rerun),2,[1 1],resi{casei}); 
                for i=1:numel(obj.tpm.dat) 
                  spm_smooth( obj.tpm.dat{i} , obj.tpm.dat{i} , rerun * [1 1 1]); 
                end
                continue
              elseif strcmp( job.opts.tpm , dtpm )
                cat_io_addwarning('cat_run_job:ownTPMerror','The new selected TPM does not work. Use SPM TPM.',2,[1 1],resi{casei}); 
                casei = casei - 1; 
                rerun = 1; 
                dtpm  = fullfile(spm('dir'),'TPM','TPM.nii'); 
                obj.tpm = spm_load_priors8(dtpm);
                continue
              else
                rerun = inf; 
              end
            end
            
            ll(casei) = resi{casei}.ll; 
            if verbs
              txt = sprintf(' (min tissue contrast: %0.3f,%0.2f) ',mincontrast(casei),resi{casei}.ll ./ (obj.samp^3 * 1000)); 
              cat_io_cprintf('g5',sprintf('%s',txt)); %,repmat('\b',1,numel(txt)) 
            end
            if ( mincontrast(casei) > acccon )  && ( casei >= runcas )
              res   = resi{casei}; 
              maxci = casei; 
              if verbs, fprintf('\n'); end
              break
            end
          end
          if ~exist('res','var') || ~exist('maxci','var')
            % if no model was optimal than chouse the best
            [maxc,maxci] = max(mincontrast(1:min(numel(mincontrast)      ,find(cellfun('isempty',modelname))))); 
         %  [maxl,maxli] = max(ll(1:min(numel(ll ./ (obj.samp^3 * 1000) ),find(cellfun('isempty',modelname))))); 
     
            if ~isempty(maxci)
              res = resi{maxci}; 
              if verbs
                stime  = cat_io_cmd(sprintf('  Select %s',modelname{maxci}),'g5','',job.extopts.verb-1,stime); 
              end
            end
          else
            if verbs, cat_io_cmd(' ','g5','',job.extopts.verb-1); end
          end
          
          if isempty(maxci) || (mincontrast(maxci) < 0.05 && job.extopts.ignoreErrors < 2) || job.test_warnings
            %% If it was not possible to run SPM then print an error message. 
            %  RD202006: Alternatively we could try the kmeans here but I  
            %  think it is better to give an error here and use manual
            %  selection by the job.extopts.init

            mati = spm_imatrix(V.mat);
            emsg = sprintf([
              'Error in spm_preproc8 that possibly happens due to \\\\n ' ...
              'bad orientation or untypical contrast. %s \\\\n' ...
              'and/or apply %s with intensity normalization. If the image is a  \\\\n' ...
              'normal MRI image with good contrast and correct orientation that you \\\\n' ...
              'can share with us than please contact us (%s). \\\\n' ...
              '  Volume size (x,y,z):   %8.0f %8.0f %8.0f \\\\n' ...
              '  Origin (x,y,z):        %8.1f %8.1f %8.1f \\\\n' ...
              '  Rotation (deg):        %8.1f %8.1f %8.1f \\\\n' ...
              '  Resolution:            %8.1f %8.1f %8.1f \\\\n'],...
              spm_file('Check image orientation/contrast','link',['spm_image(''Display'', ''' ofname ''')']), ...
              spm_file('head trimming','link','spm_jobman(''interactive'','''',''spm.tools.cat.tools.datatrimming'');'), ...
              spm_file('vbmweb@gmail.com','link','web(''mailto:vbmweb@gmail.com'');'), ...
              V.dim,[mati(1:3),mati(4:6),mati(7:9)]); 
% ### 
% write AC voxel coordinate???
% ####
            if exist('Ybi','var') && exist('Ybg','var')
              % estimate some thresholds with kmeans3D in specific ROIs
              try %#ok<TRYNC> % this has to be save although a variable is miss in some cases 
                Tcgw = kmeans3D(tmp(Ybi(:)>0.5),5); Tcgw([2,4]) = [];
                Thd  = kmeans3D(tmp(~Ybg(:) & Ybi(:)<0.5),2);
                Tbg  = kmeans3D(tmp(obj.msk.dat(:)>0 & Ybi(:)<0.5),1);
                Tth  = [Tcgw Thd Tbg]; 
                emsg = [emsg, sprintf([
                  'Normalized tissue intensities: \\\\n' ...
                  '  Brain tissues:         %8.2f %8.2f %8.2f \\\\n' ...
                  '  Head tissues / BG:     %8.2f %8.2f %8.2f'], Tth)];
                cat_err_res.res.Tth = Tth;   
              end 
            end
            
            if ~exist('res','var')
              error('cat_run_job:spm_preproc8',emsg);
            else
              if exist('Tth','var')
                cat_io_addwarning('cat_run_job:spm_preproc8',emsg,2,[1 2],Tth); 
              else
                cat_io_addwarning('cat_run_job:spm_preproc8',emsg,2,[1 2]); 
              end  
            end
         
            % udpate the last image version
            res.image.dat = tmp;
            obj.image.dat = tmp; 
            clear tmp;
          end
          
          
          fprintf('%5.0fs\n',etime(clock,stime));
        end

        
        
        %% check contrast (and convergence)
        %  RD202006: SPM peak averaging 
        %  To get one single tissue value the following definition is correct 
        %  in principle but outliers can have strong effect on mean estimation: 
        %    clsint = @(x) round( sum(res.mn(res.lkp==x) .* res.mg(res.lkp==x)') * 10^5)/10^5;
        %  So we have to be careful by using these values.
        if ~isempty(res)  
          clsint = @(x) round( sum(res.mn(res.lkp==x) .* res.mg(res.lkp==x)') * 10^5)/10^5;
          Tgw = [cat_stat_nanmean(res.mn(res.lkp==1)) cat_stat_nanmean(res.mn(res.lkp==2))]; 
          Tth = [
            ... min(res.mn(res.lkp==6 & res.mg'>0.3)) ... % bg; ignore the background, because of MP2RGAGE, R1, and MT weighted images  
            max( min( clsint(3) ,  max(Tgw) + 1.5*abs(diff(Tgw))) , min(Tgw) - 1.5*abs(diff(Tgw)) ) ... % csf with limit for T2!
            clsint(1) ... gm
            clsint(2) ... wm 
            clsint(4) ... skull
            clsint(5) ... head tissue
            clsint(6) ... background
          ];

          res.Tth = Tth; 
          cat_err_res.res = res;   

          % RD202006: Throw warning/error?
          % Due to inaccuracies of the clsint function it is better to print 
          % this as intense warning.
          if any( Tth(2:3)<0 ) || job.test_warnings
            cat_io_addwarning('cat_run_job:negVal',sprintf( ...
             ['CAT12 was developed for images with positive values and \\\\n', ...
              'negative values can lead to preprocessing problems. The average \\\\n', ...
              'intensities of CSF/GM/WM are %0.4f/%0.4f/%0.4f. \\\\n', ...
              'If you observe problems, you can use the %s to scale your data.'], Tth(1:3), ...
              spm_file('Datatype-batch','link','spm_jobman(''interactive'','''',''spm.tools.cat.tools.spmtype'');')),2,[0 1],Tth);
          end
        end

    end
    
    % updated tpm information for skull-stripped data should be available for cat_main
    if isfield(obj.tpm,'bg1') && exist('ppe','var') && ( ppe.affreg.skullstripped || job.extopts.gcutstr<0 )
      fname = res.tpm(1).fname;
      res.tpm       = obj.tpm;
      res.tpm(1).fname = fname;
    end
    spm_progress_bar('Clear');
            
    %% call main processing
    res.tpm    = obj.tpm.V;
    res.stime  = stime0;
    res.catlog = catlog; 
    if exist('Ybge','var')
      % If the background was estimated we want to save it to improve the 
      % SPM segmentation in regions outside the TPM volume. 
      res.bge = Ybge; 
    end
    res.image0 = spm_vol(job.channel(1).vols0{subj}); 
    if exist('Ylesion','var'), res.Ylesion = Ylesion; else, res.Ylesion = false(size(res.image.dim)); end; clear Ylesion;
    if exist('redspmres','var'); res.redspmres = redspmres; res.image1 = image1; end
    job.subj = subj; 
    cat_main(res,obj.tpm,job);
    
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
  WMth = kmeans3D( Ysrc(Yt(:)) , 1); 
  
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
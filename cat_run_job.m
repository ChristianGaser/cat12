function cat_run_job(job,tpm,subj)
% run CAT 
% ______________________________________________________________________
%
% Initialization function of the CAT preprocessing. 
%  * creation of the subfolder structure (if active)
%  * check of image resolution (avoid scans with very low resolution)
%  * noise correction (ISARNLM)
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

    global cat_err_res; % for CAT error report

    stime = clock;

    % print current CAT release number and subject file
    [n,r] = cat_version;
    str  = sprintf('CAT12 r%s: %d/%d',r,subj,numel(job.channel(1).vols));
    str2 = spm_str_manip(job.channel(1).vols{subj}(1:end-2),['a' num2str(70 - length(str))]);
    cat_io_cprintf([0.2 0.2 0.8],'\n%s\n%s: %s%s\n%s\n',...
          repmat('-',1,72),str,...
          repmat(' ',1,70 - length(str) - length(str2)),str2,...
          repmat('-',1,72));
    clear r str str2

    
    % create subfolders if not exist
    pth = spm_fileparts(job.channel(1).vols{subj}); 
    if job.extopts.subfolders
      folders = {'mri','report'};
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
  
    
    %  -----------------------------------------------------------------
    %  separation of full CAT preprocessing and SPM segmentation
    %  preprocessing (running DARTEL and PBT with SPM segmentation)
    %  -----------------------------------------------------------------
    [pp,ff,ee,ex] = spm_fileparts(job.data{subj}); 
    if exist(fullfile(pp,['c1' ff(3:end) ee]),'file') && ...
       exist(fullfile(pp,['c2' ff(3:end) ee]),'file') && ...
       exist(fullfile(pp,['c3' ff(3:end) ee]),'file') && ...
       exist(fullfile(pp,[ff(3:end) '_seg8.mat']),'file');
       
        job.data{subj}          = fullfile(pp,[ff ee]); 
        job.channel.vols{subj}  = fullfile(pp,[ff ee]); 

        % prepare SPM preprocessing structure 
        images = job.channel(1).vols{subj};
        for n=2:numel(job.channel)
          images = char(images,job.channel(n).vols{subj});
        end

        obj.image    = spm_vol(images);
        spm_check_orientations(obj.image);

        obj.fwhm     = job.opts.fwhm;
        obj.biasreg  = cat(1,job.opts.biasreg);
        obj.biasfwhm = cat(1,job.opts.biasfwhm);
        obj.tpm      = tpm;
        obj.lkp      = [];
        if all(isfinite(cat(1,job.tissue.ngaus))),
            for k=1:numel(job.tissue),
                obj.lkp = [obj.lkp ones(1,job.tissue(k).ngaus)*k];
            end;
        end
        obj.reg      = job.opts.warpreg;
        obj.samp     = job.opts.samp;              

        cfname  = fullfile(pp,[ff ee]);
        ofname  = fullfile(pp,[ff(3:end) ee]); 
        nfname  = fullfile(pp,mrifolder,['n' ff '.nii']); 
        copyfile(ofname,nfname); 

        res = load(fullfile(pp,[ff(3:end) '_seg8.mat']));
        job.channel(1).vols{subj}  = [nfname ex];
        job.channel(1).vols0{subj} = [ofname ex];
        res.image  = spm_vol([nfname ex]);
        res.image0 = spm_vol([ofname ex]);
        res.imagec = spm_vol([cfname ex]);
        res.spmpp  = 1; 
    else

        %  -----------------------------------------------------------------
        %  check resolution properties
        %  -----------------------------------------------------------------
        %  There were some images that should not be processed. So we have  
        %  to check for high slice thickness and low resolution.
        %  -----------------------------------------------------------------
        for n=1:numel(job.channel) 
          V = spm_vol(job.channel(n).vols{subj});
          vx_vol = sqrt(sum(V.mat(1:3,1:3).^2));

          if any(vx_vol>5)  % too thin slices
            error('CAT:cat_main:TooLowResolution', sprintf(...
                 ['Voxel resolution has to be better than 5 mm in any dimension \n' ...
                  'for reliable CAT preprocessing! \n' ...
                  'This image has a resolution %0.2fx%0.2fx%0.2f mm%s. '], ... 
                    vx_vol,char(179))); %#ok<SPERR>
          end
          if prod(vx_vol)>27  % too small voxel volume (smaller than 3x3x3 mm3)
            error('CAT:cat_main:TooHighVoxelVolume', ...
                 ['Voxel volume has to be smaller than 10 mm%s (around 3x3x3 mm%s) to \n' ...
                  'allow a reliable CAT preprocessing! \n' ...
                  'This image has a voxel volume of %0.2f mm%s. '], ...
                  char(179),char(179),prod(vx_vol),char(179));
          end
          if max(vx_vol)/min(vx_vol)>8 % isotropy 
            error('CAT:cat_main:TooStrongIsotropy', sprintf(...
                 ['Voxel isotropy (max(vx_size)/min(vx_size)) has to be smaller than 8 to \n' ...
                  'allow a reliable CAT preprocessing! \n' ...
                  'This image has a resolution %0.2fx%0.2fx%0.2f mm%s and a isotropy of %0.2f. '], ...
                  vx_vol,char(179),max(vx_vol)/min(vx_vol))); %#ok<SPERR>
          end
        end


        % save original file name 
        for n=1:numel(job.channel) 
          job.channel(n).vols0{subj} = job.channel(n).vols{subj};
        end


        % allways create the n*.nii image because of the real masking of the
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


            % noise correction
            if job.extopts.NCstr>0
              if job.extopts.sanlm==1
                stime = cat_io_cmd(sprintf('SANLM denoising (NCstr=%0.2f)',job.extopts.NCstr));
                cat_vol_sanlm(struct('data',nfname,'verb',0,'prefix','')); 
              elseif job.extopts.sanlm==2
                stime = cat_io_cmd(sprintf('ISARNLM denoising (NCstr=%0.2f)',job.extopts.NCstr));
                if job.extopts.verb>1, fprintf('\n'); end
                cat_vol_isarnlm(struct('data',nfname,'verb',(job.extopts.verb>1)*2,'prefix','')); 
                if job.extopts.verb>1, cat_io_cmd(' ','',''); end
              end
              %V = spm_vol(job.channel(n).vols{subj});
              if job.extopts.sanlm>0
                fprintf('%4.0fs\n',etime(clock,stime));   
              end
            end
        end




        %% Interpolation
        %  -----------------------------------------------------------------
        %  The interpolation can help to reduce problems for morphological
        %  operations for low resolutions and strong isotropic images. 
        %  For interpolation from 2 to 1 mm the segmentation shows more 
        %  anatomical details.
        for n=1:numel(job.channel) 

          % prepare header of resampled volume
          Vi        = spm_vol(job.channel(n).vols{subj}); 
          vx_vol    = sqrt(sum(Vi.mat(1:3,1:3).^2));
          vx_vol    = round(vx_vol*10^2)/10^2; % avoid small differences 

          % we have to look for the name of the field due to the GUI job struct generation! 
          restype   = char(fieldnames(job.extopts.restypes));
          switch restype
            case 'native'
              vx_voli  = vx_vol;
            case 'fixed', 
              vx_voli  = min(vx_vol ,job.extopts.restypes.(restype)(1) ./ ...
                         ((vx_vol > (job.extopts.restypes.(restype)(1)+job.extopts.restypes.(restype)(2)))+eps));
              vx_voli  = max(vx_voli,job.extopts.restypes.(restype)(1) .* ...
                         ( vx_vol < (job.extopts.restypes.(restype)(1)-job.extopts.restypes.(restype)(2))));
            case 'best'
              vx_voli  = min(vx_vol ,job.extopts.restypes.(restype)(1) ./ ...
                         ((vx_vol > (job.extopts.restypes.(restype)(1)+job.extopts.restypes.(restype)(2)))+eps));
              %vx_voli  = min(vx_vold,vx_voli); % guarantee Dartel resolution
            otherwise 
              error('cat_run_job:restype','Unknown resolution type ''%s''. Choose between ''fixed'',''native'', and ''best''.',restype)
          end


          % interpolation 
          if any( (vx_vol ~= vx_voli) )  

            stime = cat_io_cmd(sprintf('Internal resampling (%4.2fx%4.2fx%4.2fmm > %4.2fx%4.2fx%4.2fmm)',vx_vol,vx_voli));

            Vi        = rmfield(Vi,'private'); 
            imat      = spm_imatrix(Vi.mat); 
            Vi.dim    = round(Vi.dim .* vx_vol./vx_voli);
            imat(7:9) = vx_voli .* sign(imat(7:9));
            Vi.mat    = spm_matrix(imat);

            Vn = spm_vol(job.channel(n).vols{subj}); 
            Vn = rmfield(Vn,'private'); 
            cat_vol_imcalc(Vn,Vi,'i1',struct('interp',6,'verb',0));
            vx_vol = vx_voli;

            fprintf('%4.0fs\n',etime(clock,stime));    
          else
            vx_vol = sqrt(sum(Vi.mat(1:3,1:3).^2));
          end
          clear Vi Vn;
        end


        %  prepare SPM preprocessing structure 
        images = job.channel(1).vols{subj};
        for n=2:numel(job.channel)
            images = char(images,job.channel(n).vols{subj});
        end

        obj.image    = spm_vol(images);
        spm_check_orientations(obj.image);

        obj.fwhm     = job.opts.fwhm;
        obj.biasreg  = cat(1,job.opts.biasreg);
        obj.biasfwhm = cat(1,job.opts.biasfwhm);
        obj.tpm      = tpm;
        obj.lkp      = [];
        if all(isfinite(cat(1,job.tissue.ngaus))),
            for k=1:numel(job.tissue),
                obj.lkp = [obj.lkp ones(1,job.tissue(k).ngaus)*k];
            end;
        end
        obj.reg      = job.opts.warpreg;
        obj.samp     = job.opts.samp;              



        %% Initial affine registration.
        %  -----------------------------------------------------------------
        %  APP option with subparameter
        %  Skull-stripping is helpful for correcting affine registration of neonates and other species. 
        %  Bias correction is important for the affine registration.
        %  However, the first registation can fail and further control is required   
        % 
        %  bias = 0-5 = none, light, light threshold, light apply, fine apply (light=only for registration) 
        %  msk  = 0-4 = none, head msk, head hard, brain msk, brain hard (msk=mask only, hard=remove nonmsk)
        %  aff  = 0-1 = no affreg, affreg

        if ~strcmp(job.extopts.species,'human'), job.extopts.APP='nonhuman'; end

        if ischar(job.extopts.APP) || (isnumeric(job.extopts.APP) && numel(job.extopts.APP)==1)
          switch job.extopts.APP
            case {0,'none'},     app.bias=0; app.msk=0; app.aff=1; % old default
            case {1,'light'},    app.bias=1; app.msk=1; app.aff=1; % affreg with BC; thresholding and head masking for SPM
            case {2,'medium'},   app.bias=2; app.msk=1; app.aff=1; % no-affreg; BC and head masking for SPM  
            case {3,'strong'},   app.bias=2; app.msk=1; app.aff=0; % no-affreg; BC and head masking for SPM  
            case {4,'heavy'},    app.bias=4; app.msk=3; app.aff=0; % no-affreg; BC and brain masking for SPM  
            case {5,'nonhuman'}, app.bias=4; app.msk=4; app.aff=0; % no-affreg; BC and brain masking for SPM  
            otherwise
              app.bias  = max(0,min(4,round(job.extopts.APP,-2)/100));
              app.msk   = max(0,min(4,round(mod(job.extopts.APP,100),-1)/10)); 
              app.msk   = max(0,min(1,mod(job.extopts.APP,10))); 
          end
        end
          
          

        Affine  = eye(4);
        [pp,ff] = spm_fileparts(job.channel(1).vols{subj});
        Pbt = fullfile(pp,mrifolder,['brainmask_' ff '.nii']);
        Pb  = char(job.extopts.brainmask);
        Pt1 = char(job.extopts.T1);
        if ~isempty(job.opts.affreg)      

            %% first affine registration (with APP)
            try 
                VG = spm_vol(Pt1);
            catch
                pause(rand(1))
                VG = spm_vol(Pt1);
            end
            VF = spm_vol(obj.image(1));

            % Rescale images so that globals are better conditioned
            VF.pinfo(1:2,:) = VF.pinfo(1:2,:)/spm_global(VF);
            VG.pinfo(1:2,:) = VG.pinfo(1:2,:)/spm_global(VG);


            % APP step 1 rough bias correction 
            % --------------------------------------------------------------
            % Already for the rought initial affine registration a simple  
            % bias corrected and intensity scaled image is required, because
            % high head intensities can disturb the whole process.
            % --------------------------------------------------------------
            % ds('l2','',vx_vol,Ym, Yt + 2*Ybg,obj.image.private.dat(:,:,:)/WMth,Ym,60)
            if  app.bias ||  app.msk %job.extopts.APP  
                stime = cat_io_cmd('APP: Rough bias correction'); 
                
                try
                  [Ym,Yt,Ybg,WMth,bias,Tth] = cat_run_job_APP_init(single(obj.image.private.dat(:,:,:)),vx_vol,job.extopts.verb);
                catch
                  cat_io_cprintf('err','\n  APP failed, try APP=0! \n'); 
                  job2             = job; 
                  job2.extopts.APP = 0;
                  job2.channel(1).vols{subj} = job.channel(1).vols0{subj};
                  cat_run_job(job2,tpm,subj); 
                  return
                end
                zeroBG = cat_stat_nanmean(Ym(Ybg))<1/3;
                
                % zero background is required for images with high intensity background
                Ymc = Ym; if ~zeroBG, Ymc(Ybg) = 0; end
                
                % correct AC if it is to far away from the image center 
                % ... this code does not work - 20160303 ... remove it after 201803 
                %{
                VFimat = spm_imatrix(VF.mat); nf = 8;
                VFacvx = VFimat(1:3) ./ VFimat(7:9);
                if app.aff==1 && (any(VFacvx<VG.dim/nf) || any(VFacvx>VG.dim * (nf-1)/nf))
                  YM  = false(size(Yt)); 
                  YM(round(VF.dim(1)/nf):round(VF.dim(1) * (nf-1)/nf),...
                     round(VF.dim(2)/nf):round(VF.dim(2) * (nf-1)/nf),...
                     round(VF.dim(3)/nf):round(VF.dim(3) * (nf-1)/nf)) = true; 
                  Ytt = Yt .* YM;
                  [indx,indy,indz] = ind2sub(size(Ytt),find(Ytt==1));
                  VFimat(1:3) = mean([indx,indy,indz]) .* VFimat(7:9);
                  VF.mat = spm_matrix(VFimat);
                  stime = cat_io_cmd('Coarse affine registration with COM','','',1,stime); 
                else
                  stime = cat_io_cmd('Coarse affine registration','','',1,stime); 
                end
                %}

                % update SPM parameter - only increasing of resolution paramter 
                % experimental 20161021
                if job.extopts.experimental
                  bias = max(0,bias - 0.1);                                       % there will allways be some bias
                  bias = double(bias);                                            % SPM need double!
                  obj.biasreg  = min(0.01,obj.biasreg * 10^round(2*min(1,bias))); % less regularisation in case of bias
                  obj.biasfwhm = max(30,obj.biasfwhm  * min(1,max(0.5,1-bias)));  % reduce bias fwhm in case of bias
                  obj.samp     = obj.samp             * min(1,max(0.5,1-bias));   % increase sample distance in case of bias
                  job.opts.biasreg  = obj.biasreg;
                  job.opts.biasfwhm = obj.biasfwhm;
                  job.opts.samp     = obj.samp;
                  job.opts.bias     = bias; 
                  cat_io_cmd(sprintf('  bias~%0.2f >> biasreg=%0.0e; biasfwhm=%0.2f; samp=%0.2f',...
                    bias,obj.biasreg,obj.biasfwhm,obj.samp),'','',1,stime); 
                  fprintf('\n');
                  stime = cat_io_cmd('Coarse affine registration:','','',1); 
                else 
                  stime = cat_io_cmd('Coarse affine registration:','','',1,stime); 
                end
                

                % write data to VF
                VF.dt         = [spm_type('UINT8') spm_platform('bigend')];
                VF.dat(:,:,:) = cat_vol_ctype(Ymc * max(1/2,Tth.Tmax) * 255,'uint8'); 
                VF.pinfo      = repmat([1;0],1,size(Ymc,3));
                clear WI; 

                % smoothing
                resa  = obj.samp*2; % definine smoothing by sample size
                VF1   = spm_smoothto8bit(VF,resa);
                VG1   = spm_smoothto8bit(VG,resa);


            else
                % standard approach with static resa value and no VG smoothing
                stime = cat_io_cmd('Coarse affine registration:'); 
                resa  = 8;
                VF1   = spm_smoothto8bit(VF,resa);
                VG1   = VG; 
            end


            % prepare affine parameter 
            aflags     = struct('sep',obj.samp,'regtype','subj','WG',[],'WF',[],'globnorm',1); 
            aflags.sep = max(aflags.sep,max(sqrt(sum(VG(1).mat(1:3,1:3).^2))));
            aflags.sep = max(aflags.sep,max(sqrt(sum(VF(1).mat(1:3,1:3).^2))));

            %% affine registration
            try
                spm_plot_convergence('Init','Coarse affine registration','Mean squared difference','Iteration');
            catch
                spm_chi2_plot('Init','Coarse affine registration','Mean squared difference','Iteration');
            end
            if app.aff % job.extopts.APP~=4
                warning off 
                try 
                  [Affine0, affscale]  = spm_affreg(VG1, VF1, aflags, eye(4)); Affine = Affine0; 
                catch
                  Affine = eye(4); affscale = 0; cat_io_cmd(' ','','',1);
                end
                if affscale>3 || affscale<0.5
                  stime  = cat_io_cmd('Coarse affine registration failed. Try fine affine registration.','','',1,stime);
                  Affine = eye(4); 
                end
                warning on
            else
              Affine = eye(4); affscale = 1;
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
            aflags.sep = obj.samp/2; 
            aflags.sep = max(aflags.sep,max(sqrt(sum(VG(1).mat(1:3,1:3).^2))));
            aflags.sep = max(aflags.sep,max(sqrt(sum(VF(1).mat(1:3,1:3).^2))));
            if app.bias>2 || app.msk>2 
                %% apply (first affine) registration on the default brain mask
                VFa = VF; 
                if app.aff, VFa.mat = Affine * VF.mat; else Affine = eye(4); affscale = 1; end
                if isfield(VFa,'dat'), VFa = rmfield(VFa,'dat'); end
                [Vmsk,Yb] = cat_vol_imcalc([VFa,spm_vol(Pb)],Pbt,'i2',struct('interp',3,'verb',0)); Yb = Yb>0.5 & ~Ybg; 
                stime = cat_io_cmd('APP: Fine bias correction and skull-stripping:','','',1,stime); 

                if sum(Yb(:))==0 || mean(Ym(Yb(:)>0))<0.5 || mean(Ym(Yb(:)>1.05))<0.5 
                  %%
                  [Ymr,Ytr,resT3] = cat_vol_resize({Ym,single(Yt)},'reduceV',vx_vol,min(1.2,cat_stat_nanmean(vx_vol)*2),64,'meanm'); 
                  Ybr  = cat_vol_morph(Ytr>0.5 & Ymr>0.7 & Ymr<1.2,'o',1); 
                  Ybsr = cat_vol_smooth3X(Ybr,8);  
                  Ybsr = Ybsr./max(Ybsr(:)); 
                  Ybr  = cat_vol_morph(Ytr>0.3 & Ymr>0.7 & Ymr<1.2 & Ybsr>0.1,'o',1); 
                  Ybr  = cat_vol_morph(Ybr,'lc',6); 
                  Yb   = cat_vol_resize(Ybr,'dereduceV',resT3); 
                  Affine = eye(4); affscale = 1;
                end
  
                % fine APP
                [Ym,Yp0,Yb] = cat_run_job_APP_final(single(obj.image.private.dat(:,:,:)),...
                    Ymc,Yb,Ybg,vx_vol,job.extopts.gcutstr,job.extopts.verb);
                stime = cat_io_cmd('Affine registration:','','',1,stime); 

                % zero background is required for images with high intensity background
                Ymc = Ym; if ~zeroBG, Ymc(Ybg) = 0; end

                %% smooth data
                VF.dat(:,:,:) =  cat_vol_ctype(Ymc * max(1/2,Tth.Tmax) * 255); 
                VF1 = spm_smoothto8bit(VF,aflags.sep);
                VG1 = spm_smoothto8bit(VG,aflags.sep);

                if 1 % brain masking for affine registration 
                  VB  = spm_vol(Pb);
                  Ybt = spm_read_vols(VB); 
                  VG1.dat(:,:,:) =  cat_vol_ctype(single(VG1.dat(:,:,:)) .* smooth3(Ybt));
                  VF1.dat(:,:,:) =  cat_vol_ctype(single(VF1.dat(:,:,:)) .* smooth3(Yb));
                end   

                % using brain volume for affine scaling???
                %   cat.extopts.brainscale   = 200; % non-human brain volume in cm3 (from literature) or scaling in mm (check your data)
            elseif app.bias || app.msk 
                % smooth data
                stime = cat_io_cmd('Affine registration','','',1,stime); 
                VF.dat(:,:,:) =  cat_vol_ctype(Ymc * max(1/2,Tth.Tmax) * 255); 
                VF1 = spm_smoothto8bit(VF,aflags.sep);
                VG1 = spm_smoothto8bit(VG,aflags.sep);
            else
                % standard approach 
                stime = cat_io_cmd('Affine registration','','',1,stime); 
                VF1 = spm_smoothto8bit(VF,aflags.sep);
                VG1 = spm_smoothto8bit(VG,0.5); 
            end


            %% fine affine registration 
            try
                spm_plot_convergence('Init','Affine registration','Mean squared difference','Iteration');
            catch
                spm_chi2_plot('Init','Affine registration','Mean squared difference','Iteration');
            end
            warning off
            try
              [Affine1,affscale1] = spm_affreg(VG1, VF1, aflags, Affine, affscale); 
            catch
              Affine1 = eye(4,4); affscale1 = 0;
            end
            warning on
            
            % if also the fine affine scaling failed, use original AC!
            if any(any(isnan(Affine1(1:3,:)))) || affscale1<0.5 || affscale1>3,
              if job.extopts.APP<4
                cat_io_cprintf('err','\n  Affine registration failed, try APP=4 with brainmask! \n'); 
                job2             = job; 
                job2.extopts.APP = 5; 
                job.opts.samp    = job.opts.samp*2;
                job2.channel(1).vols{subj} = job.channel(1).vols0{subj};
                cat_run_job(job2,tpm,subj); 
                return
              elseif job.extopts.APP<5
                cat_io_cprintf('err','\n  Affine registration failed, try APP=5 with applied brainmask! \n'); 
                job2             = job; 
                job2.extopts.APP = 5;
                job.opts.samp    = job.opts.samp*2;
                job2.channel(1).vols{subj} = job.channel(1).vols0{subj};
                cat_run_job(job2,tpm,subj); 
                return
              else
                cat_io_cprintf('err','\n  Affine registration failed, use no affine transformation! \n'); 
                Affine = Affine1; 
              end
            end
            clear VG1 VF1
        end


        %% APP for spm_maff8
        %  optimize intensity range
        %  we have to rewrite the image, because SPM reads it again 
        if job.extopts.APP>0
            % WM threshold
            Ysrc = single(obj.image.private.dat(:,:,:)); 
            Ysrc(isnan(Ysrc) | isinf(Ysrc)) = min(Ysrc(:));

            if exist('Yb','var')
                Yb = cat_vol_morph(cat_vol_morph(Yb,'d',1),'lc',1);
                th = cat_stat_nanmean(Ysrc(Yb(:) & Ysrc(:)>cat_stat_nanmean(Ysrc(Yb(:))))) / ...
                     cat_stat_nanmean(Ymc(Yb(:)  & Ymc(:)>cat_stat_nanmean(Ymc(Yb(:)))));
            else % only initial bias correction
                th = WMth;
            end
            if zeroBG
              bth = mean(single(Ysrc(Ybg(:)))) - std(single(Ysrc(Ybg(:))));
            else
              bth = 0; 
            end

            % add temporary skull-stripped images
            if app.msk>2 % use brain mask
                obj.msk       = VF; 
                obj.msk.pinfo = repmat([255;0],1,size(Yb,3));
                obj.msk.dt    = [spm_type('uint8') spm_platform('bigend')];
                obj.msk.dat(:,:,:) = cat_vol_ctype(Yb); 
                obj.msk       = spm_smoothto8bit(obj.msk,0.1); 
            elseif app.msk>0 % use head mask
                obj.msk       = VF; 
                obj.msk.pinfo = repmat([255;0],1,size(Ybg,3));
                obj.msk.dt    = [spm_type('uint8') spm_platform('bigend')];
                obj.msk.dat(:,:,:) = cat_vol_ctype(1-Ybg); 
                obj.msk       = spm_smoothto8bit(obj.msk,0.1); 
            else 
                if isfield(obj,'msk'), obj = rmfield(obj,'msk'); end
            end

            % add and write bias corrected (, skull-stripped) image
            if app.bias<=1 && zeroBG
             % app.bias=1 is just a simple bias correction for affreg and will cause errors in the BWP cerebellum, if used further! 
                Ymc = single(max(bth,min( max( th + max(3*diff([bth,th])) , th / max(1/2 - 3/8*Tth.inverse,Tth.Tmax) ) ,Ysrc ))); % just limit the image intensities
            else 
                Ymc = single(max(bth,min( max( th + max(3*diff([bth,th])) , th / max(1/2 - 3/8*Tth.inverse,Tth.Tmax) ) ,Ym * th))); % use the bias corrected image
            end

            % hard masking
            if ~zeroBG,    Ymc = Ymc .* (~Ybg); end % 20161229 - simple to do this than to change all things in cat_main
            if app.msk==2, Ymc = Ymc .* (~Ybg); end
            if app.msk==4, Ymc = Ymc .* cat_vol_morph(cat_vol_morph(Yb,'d',1),'lc',1); end

            % set variable and write image
            obj.image.dat(:,:,:)         = Ymc;  
            obj.image.private.dat(:,:,:) = Ymc; 

            obj.image.dt    = [spm_type('FLOAT32') spm_platform('bigend')];
            obj.image.pinfo = repmat([1;0],1,size(Ysrc,3));
            clear Ysrc; 
        end



        %  Fine Affine Registration with 3 mm sampling distance
        %  This does not work for non human (or very small brains)
        stime = cat_io_cmd('SPM preprocessing 1:','','',1,stime);
        if strcmp('human',job.extopts.species) 
            spm_plot_convergence('Init','Fine affine registration','Mean squared difference','Iteration');
            warning off 
            Affine2 = spm_maff8(obj.image(1),obj.samp,(obj.fwhm+1)*16,obj.tpm,Affine ,job.opts.affreg); 
            Affine3 = spm_maff8(obj.image(1),obj.samp,obj.fwhm,       obj.tpm,Affine2,job.opts.affreg);
            warning on  
            if ~any(any(isnan(Affine3(1:3,:)))), Affine = Affine3; end
        end
        obj.Affine = Affine;

        % set original non-bias corrected image
        if job.extopts.APP==1
          obj.image = spm_vol(images);
        end
        cat_err_res.obj = obj; 


        %% SPM preprocessing 1
        %  ds('l2','a',0.5,Ym,Ybg,Ym,Ym,140);
        %  ds('l2','a',0.5,Ysrc/WMth,Yb,Ysrc/WMth,Yb,140);
        warning off 
        try 
          try
            res = spm_preproc8(obj);
          catch
            obj = rmfield(obj,'msk'); % try without mask ... there was an datatype error ... 
            res = spm_preproc8(obj);
          end
        catch
            if (job.extopts.sanlm && job.extopts.NCstr) || any( (vx_vol ~= vx_voli) ) || ~strcmp(job.extopts.species,'human') 
                [pp,ff,ee] = spm_fileparts(job.channel(1).vols{subj});
                delete(fullfile(pp,[ff,ee]));
            end
            error('CAT:cat_run_job:spm_preproc8','Error in spm_preproc8. Check image and orientation. \n');
        end
        warning on 


        if job.extopts.experimental
            % save information for debuging and OS test
            [pth,nam] = spm_fileparts(job.channel(1).vols0{subj}); 
            tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s_%s.mat',nam,'runjob','postpreproc8')); 
            save(tmpmat,'obj','res','Affine','Affine0','Affine1','Affine3');     
        end 
        cat_err_res.res = res;   

        fprintf('%4.0fs\n',etime(clock,stime));   

        %error('Affine Registration test error')


        %% check contrast
        Tgw = [cat_stat_nanmean(res.mn(res.lkp==1)) cat_stat_nanmean(res.mn(res.lkp==2))]; 
        Tth = [
          ... min(res.mn(res.lkp==6 & res.mg'>0.3)) ... % bg; ignore the background, because of "R" weighted images  
          max( min( mean(res.mn(res.lkp==3 & res.mg'>0.2)) , ...
            max(Tgw)+abs(diff(Tgw))) , min(Tgw)-abs(diff(Tgw)) ) ... % csf with limit for T2!
          cat_stat_nanmean(res.mn(res.lkp==1 & res.mg'>0.3)) ... gm
          cat_stat_nanmean(res.mn(res.lkp==2 & res.mg'>0.3)) ... wm 
        ];
        if isfield(obj,'msk'), res.msk = obj.msk; end

        % inactive preprocessing of inverse images (PD/T2) 
        if job.extopts.INV==0 && any(diff(Tth)<=0)
          error('CAT:cat_main:BadImageProperties', ...
          ['CAT12 is designed to work only on highres T1 images.\n' ...
           'T2/PD preprocessing can be forced on your own risk by setting \n' ...
           '"cat12.extopts.INV=1" in the cat default file. If this was a highres \n' ...
           'T1 image than the initial segmentation seemed to be corrupded, maybe \n' ...
           'by alignment problems (check image orientation).']);    
        end

    end
    
    %% call main processing
    res.stime  = stime;
    res.image0 = spm_vol(job.channel(1).vols0{subj}); 
    cat_main(res,obj.tpm,job);
    
    % delete denoised/interpolated image
    [pp,ff,ee] = spm_fileparts(job.channel(1).vols{subj});
    if exist(fullfile(pp,[ff,ee]),'file'); 
      delete(fullfile(pp,[ff,ee]));
    end
%%
return
%=======================================================================


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

    % if there is a breakpoint in this file set debug=1 and do not clear temporary variables 
    dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

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
            if job.extopts.NCstr~=0 && job.extopts.sanlm>0
              if job.extopts.sanlm==1
                stime = cat_io_cmd(sprintf('SANLM denoising (NCstr=%0.2f)',job.extopts.NCstr));
                cat_vol_sanlm(struct('data',nfname,'verb',0,'prefix','')); 
              elseif job.extopts.sanlm==2
                stime = cat_io_cmd(sprintf('ISARNLM denoising (NCstr=%0.2f)',job.extopts.NCstr));
                if job.extopts.verb>1, fprintf('\n'); end
                cat_vol_isarnlm(struct('data',nfname,'verb',(job.extopts.verb>1)*2,'prefix','')); 
                if job.extopts.verb>1, cat_io_cmd(' ','',''); end
              end
              fprintf('%4.0fs\n',etime(clock,stime));   
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
              best_vx  = max( min( 1.0 , vx_vol) ,job.extopts.restypes.(restype)(1)); 
              vx_voli  = min(vx_vol ,best_vx ./ ((vx_vol > (best_vx + job.extopts.restypes.(restype)(2)))+eps));
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
        obj.biasreg  = job.opts.biasreg;
        obj.biasfwhm = job.opts.biasfwhm; 
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
        %  bias = 0-5: none, light, light threshold, light apply, fine apply (light=only for registration) 
        %  msk  = 0-4: none, head msk, head hard, brain msk, brain hard (msk=mask only, hard=remove nonmsk)
        %  aff  = 0-1: no affreg, affreg

        if ~strcmp(job.extopts.species,'human'), job.extopts.APP='nonhuman'; end

        if ischar(job.extopts.APP) || (isnumeric(job.extopts.APP) && numel(job.extopts.APP)==1)
          switch job.extopts.APP
            case {0,'none'},     app.bias=0; app.msk=0; app.aff=1; % old default
            case {1,'light'},    app.bias=1; app.msk=0; app.aff=1; % affreg with BC; thresholding and head masking for SPM
            case {2,'medium'},   app.bias=2; app.msk=1; app.aff=1; % no-affreg; BC and head masking for SPM  
            case {3,'strong'},   app.bias=2; app.msk=2; app.aff=1; % no-affreg; BC and head masking for SPM  
            case {4,'heavy'},    app.bias=4; app.msk=3; app.aff=1; % affreg with BC; BC and brain masking for SPM  
            case {5,'nonhuman'}, app.bias=4; app.msk=4; app.aff=0; % no-affreg; BC and brain masking for SPM  
            otherwise
              app.bias  = max(0,min(4,roundx(job.extopts.APP,-2)/100));
              app.msk   = max(0,min(4,roundx(mod(job.extopts.APP,100),-1)/10)); 
              app.aff   = max(0,min(1,mod(job.extopts.APP,10))); 
          end
          if app.aff==0 && app.bias==1, app.bias=0; end
        end
          

        Affine  = eye(4);
        [pp,ff] = spm_fileparts(job.channel(1).vols{subj});
        Pbt     = fullfile(pp,mrifolder,['brainmask_' ff '.nii']);
        Pb      = char(job.extopts.brainmask);
        Pt1     = char(job.extopts.T1);

        
        if ~isempty(job.opts.affreg) && (app.aff>0 || app.bias>0 || app.msk>0)  

          % read image
          VF = spm_vol(obj.image(1).fname); 
          YF = spm_read_vols(VF); 
          
          % Detect skull-stripping, because this strongly effects the registration!
          % If a brain mask was used than we expect 
          % - many zeros (50% for small background - 80-90% for large backgrounds)
          % - a brain like volume (below 2500 cm3)
          % - only on object (the masked regions)
          % - only on background (not in very case?)
          % - less variance of thissue intensity (only 3 brain classes)
          Oth   = cat_stat_nanmean(YF(YF(:)~=0 & YF(:)>cat_stat_nanmean(YF(:)))); 
          F0vol = cat_stat_nansum(YF(:)~=0) * prod(vx_vol) / 1000; 
          F0std = cat_stat_nanstd(YF(YF(:)>0.5*Oth & YF(:)>0)/Oth); 
          YFC = YF~=0; 
          if sum(YFC(:)>0)<numel(YFC)*0.9 && sum(YFC(:)>0)>numel(YFC)*0.1  % if there is a meanful background
            YFC = ~cat_vol_morph(YF~=0,'lc',1);                            % close noisy background
          end
          [YL,numo] = spm_bwlabel(double(YF~=0),26);  clear YL;            % number of objects
          [YL,numi] = spm_bwlabel(double(YFC==0),26); clear YL;            % number of background regions 
          ppe.affreg.skullstrippedpara = [sum(YF(:)==0)/numel(YF) numo numi F0vol F0std]; 
          ppe.affreg.skullstripped = ...
            ppe.affreg.skullstrippedpara(1)>0.5 && ...                     % many zeros
            ppe.affreg.skullstrippedpara(2)<5  && ...                      % only few object
            ppe.affreg.skullstrippedpara(3)<10 && ...                      % only few background regions 
            F0vol<2500 && F0std<0.5;                                       % many zeros  and  not to big
          ppe.affreg.skullstripped = ppe.affreg.skullstripped || ...
            sum([ppe.affreg.skullstrippedpara(1)>0.8 F0vol<1500 F0std<0.4])>1; % or 2 extrem values
          if ~debug, clear YFC F0vol F0std numo numi; end 
          
          
          % Estimate the COM and initial affine matrix that is used in case 
          % of failed initial affine registration.
          xsum  = shiftdim(sum(sum(YF/Oth>0.5,2),3),1); 
          ysum  = sum(sum(YF/Oth>0.5,1),3);             
          zsum  = shiftdim(sum(sum(YF/Oth>0.5,1),2),1); 
          COM   = [find(cumsum(xsum)/sum(xsum)>0.5,1,'first') ...
                   find(cumsum(ysum)/sum(ysum)>0.5,1,'first') ...
                   find(cumsum(zsum)/sum(zsum)>0.5,1,'first')];
          COMmm = (eye(4) * VF.mat) * [COM';1];
          AffineCOM = eye(4); AffineCOM(13:15) = -COMmm(1:3); 
          %ACvox  = round(inv(Affine * VF.mat) * [ 0; 0; 0; 1]); ACvox = ACvox(1:3)';
          COMvox = round(inv(AffineCOM * VF.mat) * [0;0;0;1]);  COMvox = COMvox(1:3)';
          if ~debug, clear YF xsum ysum zsum Oth; end 
         
          
          % load template and remove the sull if the image is skull-stripped
          try 
            VG = spm_vol(Pt1);
          catch
            pause(rand(1))
            VG = spm_vol(Pt1);
          end
          if ppe.affreg.skullstripped
            % print a warning for all users because processing of
            % skull-stripped data is not standard!
            cat_io_cprintf('warn',[...
              'Detected skull-stripping or strongly masked image.\n' ...
              'Use skull-stripped initial affine registration template!\n']);
            if job.extopts.verb>1
              cat_io_cprintf('warn',sprintf(...
                '  %0.2f%%%% zeros, %d object(s), %d background region(s), %4.0f cm%s, norm. Obj. SD %0.2f.\n',...
                ppe.affreg.skullstrippedpara(1:4),char(179),ppe.affreg.skullstrippedpara(5))); 
            end
            
            % skull-stripp the template
            VB = spm_vol(Pb);
            [VB2,YB] = cat_vol_imcalc([VG,VB],Pbt,'i1 .* i2',struct('interp',3,'verb',0)); 
            VB2.dat(:,:,:) = eval(sprintf('%s(YB/max(YB(:))*255);',spm_type(VB2.dt))); 
            VB2.pinfo      = repmat([1;0],1,size(YB,3));
            VG             = spm_smoothto8bit(VB2,0.5);
            VG.dat         = cat_vol_ctype(VG.dat);
            VG.dt          = [spm_type('UINT8') spm_platform('bigend')];
            clear VB2 YB; 
          end
          
          
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
          resa  = 8; % definine smoothing by sample size
          if  app.bias ||  app.msk %job.extopts.APP  
              stime = cat_io_cmd('APP: Rough bias correction'); 

              try
                [Ym,Yt,Ybg,WMth,bias,Tth,ppe.APPi] = cat_run_job_APP_init(single(obj.image.private.dat(:,:,:)+0),vx_vol,struct('verb',job.extopts.verb));
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

              % update SPM parameter - only increasing of resolution parameter 
              % experimental 20161021
              if job.extopts.experimental 
                bias = max(0,bias - 0.1);                                       % there will allways be some bias
                bias = double(bias);                                            % SPM need double!
                obj.biasreg  = min(0.01,obj.biasreg * 10^round(2*min(1,bias))); % less regularisation in case of bias
                obj.biasfwhm = max(30,obj.biasfwhm  * min(1,max(0.5,1-bias)));  % reduce bias fwhm in case of bias
                obj.samp     = obj.samp             * min(1,max(min(vx_vol),1-bias));   % increase sample distance in case of bias
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
              VF.dt    = [spm_type('UINT8') spm_platform('bigend')];
              VF.pinfo = repmat([1;0],1,size(Ymc,3));
              VF.dat   = cat_vol_ctype(Ymc * max(1/2,Tth.Tmax) * 255); 
              clear WI; 

              % smoothing
              VF1   = spm_smoothto8bit(VF,resa);
              VG1   = spm_smoothto8bit(VG,0.5);

          else
              % standard approach with static resa value and no VG smoothing
              stime = cat_io_cmd('Coarse affine registration'); 
              
              % smoothing
              VF1   = spm_smoothto8bit(VF,resa);
              VG1   = spm_smoothto8bit(VG,0.5); 
          end

          % prepare affine parameter 
          aflags     = struct('sep',obj.samp,'regtype','subj','WG',[],'WF',[],'globnorm',1); 
          aflags.sep = max(aflags.sep,max(sqrt(sum(VG(1).mat(1:3,1:3).^2))));
          aflags.sep = max(aflags.sep,max(sqrt(sum(VF(1).mat(1:3,1:3).^2))));

          % In some cases spm_smoothto8bit gives uint16 data that has the
          % correctly intensity range (0 to 255) but lead to a "dt" error
          % in cat_spm_affreg.
          VG1.dat = cat_vol_ctype(VG1.dat);
          VF1.dat = cat_vol_ctype(VF1.dat);
   
          % masking rather than changing the template did not work yet
          if 0%ppe.affreg.skullstripped
            aflags.WG = spm_vol(Pb); 
            YF        = spm_read_vols(VF); 
            VF.dt     = [spm_type('UINT8') spm_platform('bigend')];
            VF.pinfo  = repmat([1;0],1,size(YF,3));
            VF.dat    = YF~=0; 
            aflags.WF = spm_vol(Pb); 
          end
          
          % affine registration
          spm_plot_convergence('Init','Coarse affine registration','Mean squared difference','Iteration');
          if app.aff % job.extopts.APP~=4
              %% default affreg
              
              warning off 
              try 
                [Affine0, affscale]  = cat_spm_affreg(VG1, VF1, aflags,Affine, 1); Affine = Affine0; 
              catch
                Affine0 = eye(4); affscale = 1; cat_io_cmd(' ','','',1);
              end
              
              
              %  If this registration failed and the AC is somewhere outside 
              %  than use an initial affine registion that use the COM as AC. 
              %  Use the more successfull mapping
              mat0 = spm_imatrix(Affine0);
              ACvox  = round(inv(Affine * VF.mat) * [ 0; 0; 0; 1]); ACvox = ACvox(1:3)';
              if any(any(isnan(Affine0(1:3,:)))) || affscale<0.5 || affscale>3 || (all(all(roundx(Affine0,6) == eye(4))) &&  affscale == 1) || ...
                (max(mat0(7:9))/min(mat0(7:9)))>1.2 || any(ACvox<VF.dim/6) || any(ACvox>VF.dim - VF.dim/6)
                
                % do affreg
                try 
                  [AffineCOM2, affscaleCOM]  = cat_spm_affreg(VG1, VF1, aflags, AffineCOM, 1);
                catch
                  AffineCOM2 = eye(4); affscaleCOM = 1; 
                end
                matCOM  = spm_imatrix(AffineCOM2);
                COMvox2 = round(inv(AffineCOM2 * VF.mat) * [0;0;0;1]);  COMvox2 = COMvox2(1:3)';
          
                % check affreg and use the affine registration with less different scaling factors 
                if (( (max(matCOM(7:9))/min(matCOM(7:9))) < (max(mat0(7:9))/min(mat0(7:9)))*1.05 ) && ...
                   (  abs(affscaleCOM-1) < abs(affscale-1)*1.05 )) || ...
                   (any(COMvox2>VF.dim/6) && any(COMvox2<VF.dim - VF.dim/6))
                  Affine = AffineCOM2; affscale = affscaleCOM; 
                  
                  % give another (silent) warning
                  if job.extopts.verb>1
                    cat_io_cprintf('warn','\nInitial registration failed use center of mass as AC! \n');
                    cat_io_cmd(' ','','',1);  
                  end
                end
              end
              
              
              % final check of the new affine matrix 
              if any(any(isnan(Affine(1:3,:)))) || affscale<0.5 || affscale>3 || (all(all(roundx(Affine,6) == eye(4))) &&  affscale == 1)
                stime  = cat_io_cmd('Coarse affine registration failed. Try fine affine registration.','','',1,stime);
                cat_io_cmd(' ','','',1);  Affine = eye(4); affscale = 1; cat_io_cmd(' ','','',1);
              end
              warning on
          else
            % no affine registration 
            Affine = eye(4); affscale = 1; Affine0 = eye(4); %#ok<NASGU> % Affine 0 for debuging
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
              VFa = VF; if app.aff, VFa.mat = Affine * VF.mat; else Affine = eye(4); affscale = 1; end
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
              [Ym,Yp0,Yb] = cat_run_job_APP_final(single(obj.image.private.dat(:,:,:)+0),...
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
          % guaranty uint8
          VG1.dat = cat_vol_ctype(VG1.dat);
          VF1.dat = cat_vol_ctype(VF1.dat);
          
          

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
            Affine1 = eye(4,4); affscale1 = 1;
          end
          warning on
                   
          if any(any(isnan(Affine1(1:3,:)))) || affscale1<0.5 || affscale1>3,
            Affine1 = eye(4,4); affscale1 = 1;
          end
          Affine = Affine1; 
          
          if ~debug, clear VG1 VF1; end
        end
        ppe.affreg.Affine1 = Affine; ppe.affreg.mat1 = spm_imatrix(Affine); ppe.affreg.affscale1 = affscale1;

        if job.extopts.verb>2
          cat_io_cprintf(mean(mat1(7:9)),'  Initial affine scalling: %0.2f %0.2f %0.2f\n',ppe.affreg.mat1(7:9)); 
          cat_io_cmd(' ','','',1);
        end


        % APP for spm_maff8
        %  optimize intensity range
        %  we have to rewrite the image, because SPM reads it again 
        if job.extopts.APP>0 && (app.aff>0 || app.bias>0 || app.msk>0)  
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
        if ppe.affreg.skullstripped 
          obj.msk       = VF; 
          YF            = spm_read_vols(VF); 
          obj.msk.pinfo = repmat([255;0],1,size(YF,3));
          obj.msk.dt    = [spm_type('uint8') spm_platform('bigend')];
          obj.msk.dat   = cat_vol_ctype(YF~=0); 
          obj.msk       = spm_smoothto8bit(obj.msk,0.1); 
          
          
          % update number of SPM gaussian classes ... 
          % this did not solve the main problem, but may keeps things easier 
          job.opts.ngaus(4:6) = min([2 2 2],job.opts.ngaus(4:6)); 
          obj.lkp        = [];
          for k=1:numel(job.opts.ngaus)
            job.tissue(k).ngaus = job.opts.ngaus(k);
            obj.lkp = [obj.lkp ones(1,job.tissue(k).ngaus)*k];
          end
          
          clear YF; 
        end
       

        
        %%  Fine Affine Registration with 3 mm sampling distance
        %  This does not work for non human (or very small brains)
        stime = cat_io_cmd('SPM preprocessing 1:','','',1,stime); 
        if strcmp('human',job.extopts.species) 
          % sampling >6 mm seams to fail in some cases, but higher sampling (<3 mm) did not improve the result 
          % also the variation of smoothing did not work well in all cases and lower smoothing (<=4 mm) lead to problems in some cases 
          fwhm = [0 16 12 8]; samp = [0 6 6 6]; % entry 1 is the initial affine registration 
          ppe.affreg.Affine2   = cell(size(fwhm));  ppe.affreg.Affine2{1}   = Affine1; 
          ppe.affreg.affscale2 = cell(size(fwhm));  ppe.affreg.affscale2{1} = ppe.affreg.affscale1; 
          ppe.affreg.mat2      = cell(size(fwhm));  ppe.affreg.mat2{1}      = ppe.affreg.mat1; 
          ppe.affreg.useaff2   = zeros(size(fwhm)); ppe.affreg.useaff2(1)   = 1; 
          ppe.affreg.stopiter  = zeros(size(fwhm));
          ppe.affreg.affll     = inf(size(fwhm));
          ppe.affreg.sc1{1}    = (max(ppe.affreg.mat2{1}(7:9)) / min(ppe.affreg.mat2{1}(7:9)));
          ppe.affreg.sc2{1}    = abs(ppe.affreg.affscale2{1} - 1);
          %
          for fwhmi=2:numel(fwhm)
            %%
            spm_plot_convergence('Init','Fine affine registration','Mean squared difference','Iteration');
            regtime = clock;

            % spm_maff8 registration
            warning off
            [ppe.affreg.Affine2{fwhmi}, ppe.affreg.affll(fwhmi), ppe.affreg.affh{fwhmi}] = ...
              spm_maff8(obj.image(1) , samp(fwhmi) , fwhm(fwhmi) , obj.tpm,ppe.affreg.Affine2{fwhmi} , job.opts.affreg);
            warning on  

            ppe.affreg.mat2{fwhmi}       = spm_imatrix(ppe.affreg.Affine2{fwhmi}); 
            ppe.affreg.affscale2{fwhmi}  = mean(ppe.affreg.mat2{fwhmi}(7:9)); 

            % check registration and reset values to previous results
            ppe.affreg.sc1{fwhmi} = (max(ppe.affreg.mat2{fwhmi}(7:9)) / min(ppe.affreg.mat2{fwhmi}(7:9)));
            ppe.affreg.sc2{fwhmi} = abs(ppe.affreg.affscale2{fwhmi} - 1);
            ppe.affreg.useaff2(fwhmi) = ...
              (ppe.affreg.sc1{fwhmi} < ppe.affreg.sc1{fwhmi-1}*1.05 || ...            % check scaling properties
               ppe.affreg.sc2{fwhmi} < ppe.affreg.sc2{fwhmi-1}*1.05) && ...            % ckeck total scaling
              ppe.affreg.affll(fwhmi) < ppe.affreg.affll(fwhmi-1)*1.05 && ...
              ppe.affreg.affll(fwhmi) > 0.5 && ppe.affreg.affll(fwhmi) < 1.5 && ...  % check SPM convergence criteria
              ppe.affreg.affscale2{fwhmi}>0.5 && ppe.affreg.affscale2{fwhmi}<3;      % check principle scaling range
            ppe.affreg.stopiter(fwhmi) = fwhmi>1 && (...
              ppe.affreg.useaff2(fwhmi)==0 || ... 
              ppe.affreg.sc1{fwhmi}/ppe.affreg.sc1{fwhmi-1}>1.05 || ... % stop if values get worse
              ppe.affreg.sc1{fwhmi}/ppe.affreg.sc1{fwhmi-1}>1.05 || ...
              ppe.affreg.affll(fwhmi)/ppe.affreg.affll(fwhmi-1)>1.05 || ... % stop for low changes
              abs(ppe.affreg.sc1{fwhmi} - ppe.affreg.sc1{fwhmi-1})<0.01 || ...
              abs(ppe.affreg.sc1{fwhmi} - ppe.affreg.sc1{fwhmi-1})<0.01 || ...
              ((ppe.affreg.affll(fwhmi) - ppe.affreg.affll(fwhmi-1))<0.01 && fwhmi>2));
     
            %% some information in case of debugging
            if 0 || debug || job.extopts.verb > 2
              if fwhmi==2, fprintf('\n'); end; 
              fprintf('  sc: %5.3f >> %5.3f; sc2: %5.3f >> %5.3f; conv: %5.3f > %5.3f, time: %3.0fs - use: %d, stop: %d\n',...
                (max(ppe.affreg.mat2{fwhmi-1}(7:9)) / min(ppe.affreg.mat2{fwhmi-1}(7:9))), ...
                (max(ppe.affreg.mat2{fwhmi}(7:9))   / min(ppe.affreg.mat2{fwhmi}(7:9))), ...
                abs(ppe.affreg.affscale2{fwhmi-1} - 1), abs(ppe.affreg.affscale2{fwhmi} - 1), ...
                ppe.affreg.affll(fwhmi-1), ppe.affreg.affll(fwhmi) , etime(clock,regtime) , ...
                ppe.affreg.useaff2(fwhmi) ,ppe.affreg.stopiter(fwhmi) ); 
              if ppe.affreg.stopiter(fwhmi)
                cat_io_cmd(' ','','',1); 
              end
            end

            if ppe.affreg.useaff2(fwhmi)==0
              ppe.affreg.Affine2{fwhmi}   = ppe.affreg.Affine2{fwhmi-1};
              ppe.affreg.mat2{fwhmi}      = ppe.affreg.mat2{fwhmi-1}; 
              ppe.affreg.affll(fwhmi)     = ppe.affreg.affll(fwhmi-1); 
              ppe.affreg.affscale2{fwhmi} = ppe.affreg.affscale2{fwhmi-1}; 
            end
             Affine = ppe.affreg.Affine2{fwhmi}; 
            
            %% stop iteration if the results bring no advante 
            if ppe.affreg.stopiter(fwhmi), break; end 
 
          end
          
          %%
          if job.extopts.verb>2
            cat_io_cprintf(mean(mat(7:9)),'\n  Fine affine scalling: %0.2f %0.2f %0.2f\n',ppe.affreg.mat2{fwhmi}(7:9)); 
            cat_io_cmd(' ','','',1);
          end
          
        end
        obj.Affine = Affine;
        
        
        if 0
          %% just for debugging!
          Ysrc = single(obj.image.private.dat(:,:,:)); Ym = Ysrc/cat_stat_nanmean(Ysrc(Ysrc(:)>cat_stat_nanmean(Ysrc(:)))); 
          VFa  = VF; if isfield(VFa,'dat'), VFa = rmfield(VFa,'dat'); end
          if 1 
            if 0 & isfield(ppe.affreg,'Affine2') %&& ~isempty(ppe.affreg.Affine2{end})
              VFa.mat = ppe.affreg.Affine2{3} * VF.mat; 
            else
              VFa.mat = Affine * VF.mat;
            end
          else % old registration
            VFa.mat = Affine0 * VF.mat;
          end
          [Vmsk,Yb] = cat_vol_imcalc([VFa,spm_vol(Pb)],Pbt,'i2',struct('interp',3,'verb',0)); 
          AC = round(inv(Affine * VF.mat) * [ 0; 0; 0; 1]);
          if 1
            %%
            switch 2
              case 1, ds('l2','',vx_vol,Ym,Yb>0.5,Ym,Ym,AC(3))
              case 2, ds('l2','a',vx_vol,Ym,Yb>0.5,Ym,Ym,AC(2))
              case 3, ds('l2','m',vx_vol,Ym,Yb>0.5,Ym,Ym,AC(1))
            end
          end
        end
        
        
        % set original non-bias corrected image
        if job.extopts.APP==1
          obj.image = spm_vol(images);
        end
        cat_err_res.obj = obj; 

        
        if 1
          %% Adaptive biascorrection by modification of the biasfwhm
          
          % load image and brainmask
          VFa  = VF; if isfield(VFa,'dat'), VFa = rmfield(VFa,'dat'); end; VFa.mat = Affine * VF.mat;
          [Vmsk,Yb] = cat_vol_imcalc([VFa,spm_vol(Pb)],Pbt,'i2',struct('interp',3,'verb',0)); 
          Ysrc = single(obj.image.private.dat(:,:,:));
          
          % estimate some rought tissue tresholds although there was no bias correction 
          [Ysrcr,Ybr] = cat_vol_resize({Ysrc,Yb},'reduceV',vx_vol,2,32,'meanm'); 
          Tth0 = kmeans3D(Ysrcr(Ybr(:)>0.5),3);
          Ym   = (Ysrc-Tth0(1)) / diff(Tth0(1:2:3)); 
          Yg   = cat_vol_grad(Ym,vx_vol)./max(0.1,Ym);
          
          % use the larges part of continues edge free tissue (typically the WM, if old maybe the ventricle)
          % to estimate the normalized coefficent of variation (CV) 
          Ymsk = cat_vol_morph(Yg<max(0.05,min(0.2,cat_stat_nanmean(Yg(Yb(:)>0.5 & Ym(:)>0.5)/2 ))) & Yb>0.5 & Ym<3,'l',[2 0.5])>0.5; 
          Ymb = cat_vol_median3(Ym.*Ymsk,Ymsk,Ymsk);  
          Ymr  = cat_vol_resize(Ymb,'reduceV',vx_vol,cat_stat_nanmean(vx_vol)*4,32,'max'); 
          for i=1:1, Ymr = cat_vol_localstat(Ymr,Ymr>0,2,3); end % only one iteration!
          for i=1:2, Ymr = cat_vol_localstat(Ymr,Ymr>0,2,1); end
          %bias = double(std(Ymr(Ymr(:)>0)) ./ min(abs(diff(Tth0)/diff([min(Tth0),max(Tth0)])))); 
          Ymb = cat_vol_localstat(Ymr,Ymr>0,5,4) - cat_vol_localstat(Ymr,Ymr>0,1,4);
          bias = double(mean(Ymb(Ymr(:)>0)*5) ./ min(abs(diff(Tth0)/diff([min(Tth0),max(Tth0)])))); 
          if ~debug, clear Ymr Ymb; end
          
          % prepare parameters - simple version that only reduce the biasfwhm in case of strong bias
          %obj.biasreg       = max(0.0001,min(0.01,obj.biasreg * 10^(2*min(1,bias*10))));  % less regularisation in case of bias
          obj.biasfwhm      = max(min(obj.biasfwhm,40),min(90,obj.biasfwhm * min(1.0,max(0.5,1.5-bias))));  
          job.opts.biasfwhm = obj.biasfwhm;
          job.opts.bias     = bias; 
          if bias>0.5, fprintf('\n  strong bias~%0.2f >> biasfwhm=%0.2f mm\n',bias,obj.biasfwhm); cat_io_cmd(' ','','',1); end
          
          if ~debug, clear YFa Ysrc Ym Yg Ymsk Ymb Ym; end
        end
        
        

        %% SPM preprocessing 1
        %  ds('l2','a',0.5,Ym,Ybg,Ym,Ym,140);
        %  ds('l2','a',0.5,Ysrc/WMth,Yb,Ysrc/WMth,Yb,140);
        warning off 
        lowngaus = 0; 
        try 
          
          try
            % spm preprocessing
            res = spm_preproc8(obj);
          catch
            lowngaus = 1; 
            
            %if job.extopts.verb, cat_io_cprintf('warn','  Error: spm_preproc8 failed, try spm_preproc8 with brainmask.\n'); end
            %  obj = rmfield(obj,'msk'); % try without mask ... there was an datatype error ... 
            %  res = spm_preproc8(obj);
            %end
          end
          
          if lowngaus && ~ppe.affreg.skullstripped 
            job.opts.ngaus = max([3 3 2 3 4 2],job.opts.ngaus); 
            obj.lkp        = [];
            for k=1:numel(job.opts.ngaus)
              job.tissue(k).ngaus = job.opts.ngaus(k);
              obj.lkp = [obj.lkp ones(1,job.tissue(k).ngaus)*k];
            end
            
            % spm preprocessing with increased number of gaussian classes
            res = spm_preproc8(obj);
          end
            
        catch

          if job.extopts.APP<1
            cat_io_cprintf('err','\n  Affine registration failed, try APP=1! \n'); 
            job2             = job; 
            job2.extopts.APP = 1; 
            job.opts.samp    = job.opts.samp*2;
            job2.channel(1).vols{subj} = job.channel(1).vols0{subj};
            cat_run_job(job2,tpm,subj); 
            return
          elseif job.extopts.APP<4
            cat_io_cprintf('err','\n  Affine registration failed, try APP=4 with brainmask! \n'); 
            job2             = job; 
            job2.extopts.APP = 4; 
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
          end
          
          
          if (job.extopts.sanlm && job.extopts.NCstr) || any( (vx_vol ~= vx_voli) ) || ~strcmp(job.extopts.species,'human')
            [pp,ff,ee] = spm_fileparts(job.channel(1).vols{subj});
            delete(fullfile(pp,[ff,ee]));
            error('CAT:cat_run_job:spm_preproc8','Error in spm_preproc8. Check image and orientation. \n');
          end
        end
        
        % eigentlich schon , au?er SPM geht schief
        obj.Affine = res.Affine;
        
        warning on 
        if exist('ppe','var'), res.ppe = ppe; end

        if job.extopts.experimental
            % save information for debuging and OS test
            [pth,nam] = spm_fileparts(job.channel(1).vols0{subj}); 
            tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s_%s.mat',nam,'runjob','postpreproc8')); 
            save(tmpmat,'obj','res','Affine','Affine0','Affine1','Affine2');     
        end 
        cat_err_res.res = res;   

        fprintf('%4.0fs\n',etime(clock,stime));   

        %error('Affine Registration test error')


        %% check contrast
        clsint = @(x) round( sum(res.mn(res.lkp==x) .* res.mg(res.lkp==x)') * 10^5)/10^5;
        Tgw = [cat_stat_nanmean(res.mn(res.lkp==1)) cat_stat_nanmean(res.mn(res.lkp==2))]; 
        Tth = [
          ... min(res.mn(res.lkp==6 & res.mg'>0.3)) ... % bg; ignore the background, because of "R" weighted images  
          max( min( clsint(3) ,  max(Tgw)+abs(diff(Tgw))) , min(Tgw)-abs(diff(Tgw)) ) ... % csf with limit for T2!
          clsint(1) ... gm
          clsint(2) ... wm 
        ];
        if isfield(obj,'msk'), res.msk = obj.msk; end
        
        if ppe.affreg.skullstripped 
           % this is important for the second SPM processing!
           res.mn(res.lkp==4) = 0; 
           res.mn(res.lkp==5) = 0; 
           res.mn(res.lkp==6) = 0; 
        end
           
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
function r = roundx(r,rf)
  r(:) = round(r(:) * rf) / rf;
return
%=======================================================================

function cat_run_job(job,estwrite,tpm,subj)
%#ok<*WNOFF,*WNON>

    stime = clock;

    %% print current CAT release number and subject file
    A = ver; r = 0;
    for i=1:length(A)
        if strcmp(A(i).Name,'Computational Anatomy Toolbox')
            r = str2double(A(i).Version);
        end
    end
    clear A 

    str  = sprintf('CAT12 r%d',r);
    str2 = spm_str_manip(job.channel(1).vols{subj},['a' num2str(70 - length(str))]);
    cat_io_cprintf([0.2 0.2 0.8],'\n%s\n%s: %s%s\n%s\n',...
          repmat('-',1,72),str,...
          repmat(' ',1,70 - length(str) - length(str2)),str2,...
          repmat('-',1,72));
    clear r str str2


    %  -----------------------------------------------------------------
    %  check resolution properties
    %  -----------------------------------------------------------------
    %  There were some images that should not be processed. So we have  
    %  to check for high slice thickness, low resolution.
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
    
    % noise-correction
    if job.cat.sanlm && job.extopts.NCstr

        switch job.cat.sanlm
          case {1,2,3,4}, stime = cat_io_cmd('NLM-filter with multi-threading');
          case {5},       stime = cat_io_cmd('Temporary NLM-filter with multi-threading');
        end


        for n=1:numel(job.channel) 
            V = spm_vol(job.channel(n).vols{subj});
            Y = single(spm_read_vols(V));
            Y(isnan(Y)) = 0;
            switch job.cat.sanlm
              case {1,2,3,4,5},  
                % use isarnlm-filter only if voxel size <= 0.7mm
                if any(round(vx_vol*100)/100<=0.70)
                  if job.extopts.verb>1, fprintf('\n'); end
                  Y = cat_vol_isarnlm(Y,V,job.extopts.verb>1);   % use iterative multi-resolution multi-threaded version
                  if job.extopts.verb>1, cat_io_cmd(' '); end
                else
                  cat_sanlm_mex(Y,3,1,0);          % use multi-threaded version
                  fprintf(sprintf('%s',repmat('\b',1,numel('Using 8 processors '))));
                end
            end
            Vn = cat_io_writenii(V,Y,'n','noise corrected','float32',[0,1],[1 0 0]);
            job.channel(n).vols{subj} = Vn.fname;
            clear Y V Vn;
        end

        fprintf('%4.0fs\n',etime(clock,stime));     
    else
       if job.extopts.APP>0 || ~strcmp(job.cat.species,'human')
         % this is necessary because of the real masking of the T1 data 
         % for spm_preproc8 that include rewriting the image!
         for n=1:numel(job.channel) 
           [pp,ff,ee] = spm_fileparts(job.channel(n).vols{subj}); 
           ofname  = fullfile(pp,[ff ee]); 
           nfname  = fullfile(pp,['n' ff '.nii']); 
           copyfile(ofname,nfname); 
           job.channel(n).vols{subj} = nfname;
         end
       end
    end
   
      
    
    %% Interpolation
    % The interpolation can help to reduce problems for morphological
    % operations for low resolutions and strong isotropic images. 
    % Especially for Dartel registration a native resolution higher than the Dartel 
    % resolution helps to reduce normalization artifacts of the
    % deformations. Furthermore, even if artifacts can be reduced by the final smoothing
    % it is much better to avoid them.  
    Vt      = tpm.V(1); 
    vx_vold = min(cat_get_defaults('extopts.vox'),sqrt(sum(Vt.mat(1:3,1:3).^2))); clear Vt; % Dartel resolution 
    for n=1:numel(job.channel) 

      % prepare header of resampled volume
      Vi        = spm_vol(job.channel(n).vols{subj}); 
      vx_vol    = sqrt(sum(Vi.mat(1:3,1:3).^2));
      switch job.cat.restype 
        case 'native'
          vx_voli  = vx_vol;
        case 'fixed', 
          vx_voli  = min(vx_vol ,job.cat.resval(1) ./ ((vx_vol > (job.cat.resval(1)+job.cat.resval(2)))+eps));
          vx_voli  = max(vx_voli,job.cat.resval(1) .* (vx_vol < (job.cat.resval(1)-job.cat.resval(2))));
        case 'best'
          vx_voli  = min(vx_vol ,job.cat.resval(1) ./ ((vx_vol > (job.cat.resval(1)+job.cat.resval(2)))+eps));
          vx_voli  = min(vx_vold,vx_voli); % guarantee Dartel resolution
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
        if ~(job.cat.sanlm && job.extopts.NCstr)
          % if no noise correction we have to add the 'n' prefix here
          [pp,ff,ee] = spm_fileparts(Vn.fname);
          Vi.fname = fullfile(pp,['n' ff ee]);
          job.channel(n).vols{subj} = Vi.fname;
        end
        if job.cat.sanlm==0
          [pp,ff,ee,dd] = spm_fileparts(Vn.fname); 
          Vi.fname = fullfile(pp,['n' ff ee dd]);
          job.channel(n).vols{subj} = Vi.fname;
        end
        cat_vol_imcalc(Vn,Vi,'i1',struct('interp',6,'verb',0));
        vx_vol = vx_voli;

        fprintf('%4.0fs\n',etime(clock,stime));    
      end
      clear Vi Vn;
    end
    
    
    %%
    if estwrite % estimate and write segmentations            

        % 
        images = job.channel(1).vols{subj};
        for n=2:numel(job.channel)
            images = char(images,job.channel(n).vols{subj});
        end

        obj.image    = spm_vol(images);
        spm_check_orientations(obj.image);

        obj.fwhm     = job.cat.fwhm;
        obj.fudge    = 5;
        obj.biasreg  = cat(1,job.biasreg);
        obj.biasfwhm = cat(1,job.biasfwhm);
        obj.tpm      = tpm;
        obj.lkp      = [];
        if all(isfinite(cat(1,job.tissue.ngaus))),
            for k=1:numel(job.tissue),
                obj.lkp = [obj.lkp ones(1,job.tissue(k).ngaus)*k];
            end;
        end
        obj.reg      = job.cat.reg;
        obj.samp     = job.cat.samp;              
       


        %% Initial affine registration.
        
        % APP option with subparameter
        % Skull-stripping is helpful for correcting affine registration of neonates and other species. 
        % Bias correction is important for the affine registration.
        % However, the first registation can fail,
        if ~strcmp(job.cat.species,'human'), job.extopts.APP=2; end
        switch job.extopts.APP
          case 0 % no APP
            doskullstripping = 0;
            dobiascorrection = 0;
            doregistration   = 1;
          case 1 % APP light with simple bias correction
            doskullstripping = 0;
            dobiascorrection = 1;
            doregistration   = 1; % not available - allways active
          case 2 % APP with full bias correction, but without brain mask
            doskullstripping = 0;
            dobiascorrection = 1;
            doregistration   = 1;
          case 3 % full APP without initial affine registration (AC-PC has to be correct) for other species
            doskullstripping = 1;
            dobiascorrection = 1;
            doregistration   = 0;
          case 4 % full APP (full bias correction and brain masking)
            doskullstripping = 1;
            dobiascorrection = 1;
            doregistration   = 1;
        end
       
        Affine  = eye(4);
        [pp,ff] = spm_fileparts(job.channel(1).vols{subj});
        Pbt = fullfile(pp,['brainmask_' ff '.nii']);
        Pb  = char(cat_get_defaults('extopts.brainmask'));
        Pt1 = char(cat_get_defaults('extopts.T1'));
        if ~isempty(job.cat.affreg)
         

          
          
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
          %   ds('l2','',vx_vol,Ym, Yt + 2*Ybg,obj.image.private.dat(:,:,:)/WMth,Ym,60)
          if job.extopts.APP>0  
            stime = cat_io_cmd('APP1: rough bias correction'); 
            [Ym,Yt,Ybg,WMth] = APP_initial_bias_correction(single(obj.image.private.dat(:,:,:)),vx_vol,job.extopts.verb);;
            
            stime = cat_io_cmd('Coarse Affine registration','','',1,stime); 
            
            % write data to VF
            VF.dt         = [spm_type('UINT8') spm_platform('bigend')];
            VF.dat(:,:,:) = uint8(Ym * 200); 
            VF.pinfo      = repmat([1;0],1,size(Ym,3));
            clear WI; 
            
            % smoothing
            resa = obj.samp*3; % definine smoothing by sample size
            VF1  = spm_smoothto8bit(VF,resa);
            VG1  = spm_smoothto8bit(VG,resa);
          else
            stime = cat_io_cmd('Coarse Affine registration'); 
            
          % old approach with static resa value and no VG smoothing
            resa = 8;
            VF1  = spm_smoothto8bit(VF,resa);
            VG1  = VG; 
          end
          
          %% prepare affine parameter 
          aflags     = struct('sep',resa,'regtype',job.cat.affreg,'WG',[],'WF',[],'globnorm',0);
          aflags.sep = max(aflags.sep,max(sqrt(sum(VG(1).mat(1:3,1:3).^2))));
          aflags.sep = max(aflags.sep,max(sqrt(sum(VF(1).mat(1:3,1:3).^2))));

          % affine resistration
          try
            spm_plot_convergence('Init','Coarse Affine Registration','Mean squared difference','Iteration');
          catch
            spm_chi2_plot('Init','Coarse Affine Registration','Mean squared difference','Iteration');
          end
          if doregistration
            warning off 
            [Affine0, scale]  = spm_affreg(VG1, VF1, aflags, eye(4)); Affine = Affine0; 
            warning on
          end
          
          % improve sampling rate and use less smoothing and by brain masking
          try 
            aflags.WG  = spm_vol(Pb);
          catch
            pause(rand(1))
            aflags.WG  = spm_vol(Pb);
          end
          aflags.sep = aflags.sep/2;
          
          
          %% APP step 2 - brainmasking and second tissue separated bias correction  
          if job.extopts.APP>1
            %  ---------------------------------------------------------
            %  The second part of APP maps a brainmask to native space 
            %  and refines this mask by morphologic operations and
            %  region-growing to adapt for bad affine alignments. 
            %  The mask has to be a little bit wider and it is important
            %  that the brain is complete!
            %  ---------------------------------------------------------
            %  ds('l2','',vx_vol,Ym,Yb,Ym,Yp0,90)
            
            % apply (first affine) registration on the default brain mask
            VFa = VF; 
            if doregistration, VFa.mat = Affine0 * VF.mat; else Affine = eye(4); end
            if isfield(VFa,'dat'), VFa = rmfield(VFa,'dat'); end
            [Vmsk,Yb] = cat_vol_imcalc([VFa,spm_vol(Pb)],Pbt,'i2',struct('interp',3,'verb',0)); Yb = Yb>0.5; 
       
            stime = cat_io_cmd(sprintf('APP%d: fine bias correction',job.extopts.APP),'','',1,stime); 
            [Ym,Yp0,Yb] = APP_final_bias_correction(single(obj.image.private.dat(:,:,:)),Ym,Yb,Ybg,vx_vol,job.extopts.verb);;
            stime = cat_io_cmd('Affine registration','','',1,stime); 
                                  
            %% msk T1 & TPM
            VF.dat(:,:,:) =  cat_vol_ctype(Ym*200 .* Yb); 
            VF1 = spm_smoothto8bit(VF,aflags.sep/2);
            
            VG1 = spm_smoothto8bit(VG,0.1);
            VG1.dat = VG1.dat .* uint8(spm_read_vols(spm_vol(Pb))>0.5); 
            VG1 = spm_smoothto8bit(VG1,aflags.sep/2);
          elseif job.extopts.APP==1
            %% msk T1 & TPM
            stime = cat_io_cmd('Affine registration','','',1,stime); 
            VF.dat(:,:,:) =  cat_vol_ctype(Ym*200); 
            VF1 = spm_smoothto8bit(VF,aflags.sep/2);
            VG1 = spm_smoothto8bit(VG1,aflags.sep/2);
          else
          % old approach ... only smoothing of the VF with 2 mm
            stime = cat_io_cmd('Affine registration','','',1,stime); 
            VF1 = spm_smoothto8bit(VF,aflags.sep/2);
            VG1 = VG; 
          end

          
          % fine affine registration 
          try
            spm_plot_convergence('Init','Coarse Affine Registration 2','Mean squared difference','Iteration');
          catch
            spm_chi2_plot('Init','Coarse Affine Registration 2','Mean squared difference','Iteration');
          end
          warning off
          Affine1 = spm_affreg(VG1, VF1, aflags, Affine, scale);   
          warning on
          if ~any(isnan(Affine1(1:3,:))), Affine = Affine1; end
            
          clear VG1 VF1
        end
        if job.extopts.APP && dobiascorrection
          Ysrc = single(obj.image.private.dat(:,:,:)); 
          obj.image.dt    = [spm_type('FLOAT32') spm_platform('bigend')];
          obj.image.pinfo = repmat([1;0],1,size(Ysrc,3));
        
          % rewrite bias correctd, but not skull-stripped image
          %obj.image.private.dat(:,:,:) = single(max(-WMth*0.1,min(4*WMth,Ym * th))); 
          % add temporary skull-stripped images
          if doskullstripping
             th = cat_stat_nanmean(Ysrc(Yb(:) & Ysrc(:)>cat_stat_nanmean(Ysrc(Yb(:))))) / ...
                 cat_stat_nanmean(Ym(Yb(:)   & Ym(:)>cat_stat_nanmean(Ym(Yb(:)))));
             obj.image.dat(:,:,:) = single(max(-WMth*0.1,min(4*WMth,Ym * th))); % .* Yb))); 
          else
            if exist('Yb','var')
              th = cat_stat_nanmean(Ysrc(Yb(:) & Ysrc(:)>cat_stat_nanmean(Ysrc(Yb(:))))) / ...
                   cat_stat_nanmean(Ym(Yb(:)   & Ym(:)>cat_stat_nanmean(Ym(Yb(:)))));
            else % only initial bias correction
              th = WMth;
            end
            obj.image.dat(:,:,:) = single(max(-WMth*0.1,min(4*WMth,Ym * th))); 
            VFa = VF; 
            if doregistration, VFa.mat = Affine0 * VF.mat; else Affine = eye(4); end
            if isfield(VFa,'dat'), VFa = rmfield(VFa,'dat'); end
            [Vmsk,Yb] = cat_vol_imcalc([VFa,spm_vol(Pb)],Pbt,'i2',struct('interp',3,'verb',0)); Yb = Yb>0.5; 
          end
          obj.msk       = VF; 
          obj.msk.pinfo = repmat([255;0],1,size(Yb,3));
          obj.msk.dt    = [spm_type('uint8') spm_platform('bigend')];
          obj.msk.dat(:,:,:) = uint8(Yb); 
          obj.msk       = spm_smoothto8bit(obj.msk,0.1); 
          clear Ysrc; 
        end

        if job.extopts.APP
          stime = cat_io_cmd(sprintf('SPM preprocessing 1 (APP=%d):',job.extopts.APP),'','',1,stime);
        else
          stime = cat_io_cmd('SPM preprocessing 1:','','',1,stime);
        end
        
        
        
        %% Fine Affine Registration with 3 mm sampling distance
        % This does not work for non human (or very small brains)
        if strcmp('human',cat_get_defaults('extopts.species')) % job.extopts.APP==0
          spm_plot_convergence('Init','Fine Affine Registration','Mean squared difference','Iteration');
          warning off 
          Affine3 = spm_maff8(obj.image(1),obj.samp,obj.fudge,obj.tpm,Affine,job.cat.affreg);
          warning on  
          if ~any(isnan(Affine3(1:3,:))), Affine = Affine3; end
          
          if ~doskullstripping
            VFa = VF; 
            if doregistration, VFa.mat = Affine3 * VF.mat; else Affine = eye(4); end
            if isfield(VFa,'dat'), VFa = rmfield(VFa,'dat'); end
            [Vmsk,Yb] = cat_vol_imcalc([VFa,spm_vol(Pb)],Pbt,'i2',struct('interp',3,'verb',0)); Yb = Yb>0.5; 
            obj.msk       = VF; 
            obj.msk.pinfo = repmat([255;0],1,size(Yb,3));
            obj.msk.dt    = [spm_type('uint8') spm_platform('bigend')];
            obj.msk.dat(:,:,:) = uint8(Yb); 
            obj.msk = spm_smoothto8bit(obj.msk,0.1); 
          end
        end
        obj.Affine = Affine;
        
        
        %% SPM preprocessing 1
        % ds('l2','a',0.5,Ysrc/WMth,Yb,Ysrc/WMth,Yb,140);
        warning off 
        try 
          res = spm_preproc8(obj);
        catch
          if (job.cat.sanlm && job.extopts.NCstr) || any( (vx_vol ~= vx_voli) ) || ~strcmp(job.cat.species,'human') 
            [pp,ff,ee] = spm_fileparts(job.channel(1).vols{subj});
            delete(fullfile(pp,[ff,ee]));
          end
          error('CAT:cat_run_job:spm_preproc8','Error in spm_preproc8. Check image and orientation. \n');
        end
        warning on 
        
        if cat_get_defaults('extopts.debug')==2
          % save information for debuging and OS test
          [pth,nam] = spm_fileparts(job.channel(1).vols0{subj}); 
          tmpmat = fullfile(pth,sprintf('%s_%s_%s.mat',nam,'runjob','postpreproc8')); 
          save(tmpmat,'obj','res','Affine','Affine0','Affine1','Affine3');      
        end 
        
        fprintf('%4.0fs\n',etime(clock,stime));   
    end
    
    %% check contrast
    Tgw = [mean(res.mn(res.lkp==1)) mean(res.mn(res.lkp==2))]; 
    Tth = [
      max( min(mean(res.mn(res.lkp==3)) , max(Tgw)+abs(diff(Tgw))),min(Tgw)-abs(diff(Tgw)) ) ... % csf with limit for T2!
      mean(res.mn(res.lkp==1)) ... gm
      mean(res.mn(res.lkp==2)) ... wm 
    ];
    
    % inactive preprocessing of inverse images (PD/T2) 
    if cat_get_defaults('extopts.INV')==0 && any(diff(Tth)<=0)
      error('CAT:cat_main:BadImageProperties', ...
      ['CAT12 is designed to work only on highres T1 images.\n' ...
       'T2/PD preprocessing can be forced on your own risk by setting \n' ...
       '''cat.extopts.INV=1'' in the cat default file. If this was a highres \n' ...
       'T1 image than the initial segmentation seemed to be corrupded, maybe \n' ...
       'by alignment problems (check image orientation).']);    
    end
            
    %% Final iteration, so write out the required data.
    tc = [cat(1,job.tissue(:).native) cat(1,job.tissue(:).warped)];
    bf = job.bias;
    df = job.cat.warps;
    lb = job.label;
    jc = job.jacobian;
    res.stime = stime;
    res.image0 = spm_vol(job.channel(1).vols0{subj}); 
    cat_main(res, tc, bf, df, lb, jc, job.cat, obj.tpm, job);
    
    % delete denoised/interpolated image
    if (job.cat.sanlm && job.extopts.NCstr) || any( (vx_vol ~= vx_voli) ) || ~strcmp(job.cat.species,'human') 
      [pp,ff,ee] = spm_fileparts(job.channel(1).vols{subj});
      delete(fullfile(pp,[ff,ee]));
    end
%%
return
%=======================================================================

%=======================================================================
function r = roundx(r,rf)
  r(:) = round(r(:) * rf) / rf;
return
%=======================================================================

%=======================================================================
function  [Ym,Yp0,Yb] = APP_final_bias_correction(Ysrco,Ym,Yb,Ybg,vx_vol,verb)
%  ---------------------------------------------------------------------
%  rough scull-stripping:
%  The affine registration, especially spm_preproc8 requires a very good masking!
%  Because we it is also required for the Unified Segmenation
%  a wider mask with a complete brain is important
%    ds('l2','m',0.5,Ym*0.7+0.3,Yb,Ysrc/WMth,Ym,80)
%  ---------------------------------------------------------------------
  if verb, fprintf('\n'); end
  stime = cat_io_cmd('  Initialize','g5','',verb);
  msize = 222; %round(222 ./ max(size(Ysrco).*vx_vol) .* min(size(Ysrco).*vx_vol));  
  
  [Ysrc,Ym,resT3] = cat_vol_resize({Ysrco,Ym},'reduceV',vx_vol,min(1.5,min(vx_vol)*2),msize,'meanm'); 
  [Yb,Ybg]        = cat_vol_resize({single(Yb),single(Ybg)},'reduceV',vx_vol,min(1.5,min(vx_vol)*2),msize,'meanm'); 
  Ybg = Ybg>0.5;
  
  Yg   = cat_vol_grad(Ym,resT3.vx_volr) ./ max(eps,Ym); 
  Ydiv = cat_vol_div(Ym,resT3.vx_volr)./Ym;
  Ygs  = smooth3(Yg);
  
  % greater mask
  [dilmsk,resT2] = cat_vol_resize(Yb,'reduceV',resT3.vx_volr,mean(resT3.vx_volr)*2,32); 
  dilmsk  = cat_vbdist(dilmsk,true(size(dilmsk)),resT2.vx_volr);
  dilmsk  = dilmsk - cat_vbdist(single(dilmsk>0),true(size(dilmsk)),resT2.vx_volr);
  dilmsk  = cat_vol_resize(smooth3(dilmsk),'dereduceV',resT2); 
  headradius = -min(dilmsk(:)); 

  % thresholds
  rf   = 10^6; 
  Hth  = roundx(cat_stat_nanmean(Ym(Ym(:)>0.4 & Ym(:)<1.2  & Ygs(:)<0.2 & ~Yb(:) & Ydiv(:)<0.05 & Ydiv(:)>-0.5 & dilmsk(:)>0 & dilmsk(:)<10)),rf); % average intensity of major head tissues
  GMth = roundx(cat_stat_nanmean(Ym(Ym(:)>0.2 & Ym(:)<0.9  & Ygs(:)<0.2 & ~Yb(:) & Ydiv(:)<0.1 & Ydiv(:)>-0.1)),rf);  % first guess of the GM intensity
  CMth = roundx(cat_stat_nanmean(Ym(Ym(:)>0.1 & Ym(:)<GMth*0.7 & Ygs(:)<0.2 & ~Yb(:) & Ydiv(:)>-0.05)),rf);  % first guess of the CSF intensity
  %WMth = cat_stat_nanmean(Ym(Ym(:)>0.8 & Ym(:)<1.2 & Ygs(:)<0.2 & ~Yb(:) & Ydiv(:)>-0.05)); 
  BGth = roundx(cat_stat_nanmean(Ym(Ybg(:))),rf); 
  
  
  %% Skull-Stripping
  stime = cat_io_cmd('  Skull-Stripping','g5','',verb,stime);
  Yb = (dilmsk<20 & Ym<1.2) & (Ym>(GMth*0.7+0.3)) & Yg<0.5 & Ydiv<0.2; Yb(smooth3(Yb)<0.5)=0; 
  Yb = single(smooth3(cat_vol_morph(Yb,'l'))>0.2); 
  [dilmsk2,resT2] = cat_vol_resize(single(Yb),'reduceV',resT3.vx_volr,mean(resT3.vx_volr)*4,32); 
  dilmsk2  = cat_vbdist(dilmsk2,true(size(dilmsk)),resT2.vx_volr);
  dilmsk2  = cat_vol_resize(smooth3(dilmsk2),'dereduceV',resT2); 
  %% WM growing
  Yb(Yb<0.5 & (dilmsk2>20 | Ym>1.1 | Ym<mean([1,GMth]) | (Yg.*Ym)>0.5))=nan;
  [Yb1,YD] = cat_vol_downcut(Yb,Ym,0.05); 
  Yb(isnan(Yb))=0; Yb((YD.*dilmsk2/10)<400/mean(resT3.vx_volr))=1; Yb(isnan(Yb))=0;
  Yb = smooth3(Yb)>0.2; 
  % ventricle closing
  [Ybr,resT2] = cat_vol_resize(single(Yb),'reduceV',resT3.vx_volr,mean(resT3.vx_volr)*4,32); 
  Ybr = cat_vol_morph(Ybr>0.5,'lc',4);
  Ybr = cat_vol_resize(smooth3(Ybr),'dereduceV',resT2)>0.5; 
  Yb = single(Yb | (Ym<1.2 & Ybr));
  % GWM growing
  Yb(Yb<0.5 & (dilmsk2>20 | Ym>1.1 | Ym<GMth | (Yg.*Ym)>0.5))=nan;
  [Yb1,YD] = cat_vol_downcut(Yb,Ym,0.01); 
  Yb(isnan(Yb))=0; Yb((YD.*dilmsk2/10)<400/mean(resT3.vx_volr))=1; Yb(isnan(Yb))=0;
  Yb = smooth3(Yb)>0.5; 
  Yb = single(Yb | (Ym>0.3 & Ym<1.2 & cat_vol_morph(Yb,'lc',2)));
  % GM growing
  Yb(Yb<0.5 & (dilmsk2>20 | Ym>1.1 | Ym<CMth | (Yg.*Ym)>0.5))=nan;
  [Yb1,YD] = cat_vol_downcut(Yb,Ym,0.00);
  Yb(isnan(Yb))=0; Yb((YD.*dilmsk2/10)<200/mean(resT3.vx_volr))=1; Yb(isnan(Yb))=0; clear Yb1 YD; 
  Yb(smooth3(Yb)<0.5)=0;
  Yb = single(Yb | (Ym>0.1 & Ym<1.1 & cat_vol_morph(Yb,'lc',4)));
  % CSF growing (add some tissue around the brain)
  Yb(Yb<0.5 & (dilmsk2>30 | Ym<0.1 | Ym>1.1 | (Yg.*Ym)>0.5))=nan;
  [Yb1,YD] = cat_vol_downcut(Yb,Ym,-0.01); Yb(isnan(Yb))=0; 
  Yb(YD<100/mean(resT3.vx_volr))=1; Yb(isnan(Yb))=0; clear Yb1 YD; 
  Yb(smooth3(Yb)<0.7)=0;
  Yb  = single(cat_vol_morph(Yb,'lab'));
  Yb  = Yb | (Ym>0.2 & Ym<0.6 & cat_vol_morph(Yb,'lc',2));
  Yb  = cat_vol_smooth3X(Yb,2)>0.3;
  Ymo = Ym; 


  %% ---------------------------------------------------------
  %  improve bias correction:
  %  Also here it is important that the bias field is very smooth
  %  to avoid overcorrections. In contrast to the first
  %  correction we will try to separate between tissues...
  %  ---------------------------------------------------------
  stime = cat_io_cmd('  Tissue classification','g5','',verb,stime);
  Ym   = Ymo;
  Yg   = cat_vol_grad(Ym,vx_vol);%./max(eps,Ym)
  Ygs  = smooth3(Yg);
  Ydiv = cat_vol_div(Ym,vx_vol);

      % WM Skeleton:  Ydiv./Yg<-1    
% nicht WM: Ydiv.*Yg<-0.02

  %% tissue classes WM, GM, subcortical GM (Ybm), CSF, head tissue (Yhm) 
  Ywm  = ((Ym-max(0,(Ydiv+0.01))*10)>(GMth*0.1+0.9) | Ym>(GMth*0.1+0.9) | ...
          (Ydiv./Yg<0.5 & ((Ym>0.9 & Yg>0.1 & Ydiv<0) | (Ym>0.9)) & Ym<1.2) ) & ...
          Ym<1.3 & Yg<0.6 & Ygs<0.9 & Yb & Ydiv<0.1 & Ydiv>-0.5; 
  Ywm  = smooth3(Ywm)>0.5;      
  % subcotical GM 
  Ybm  = ((Ym-max(0,(Ydiv+0.01))*10)<0.98) & Ym<0.98 & dilmsk<-headradius*0.3 & ... 
         Ym>(GMth*0.6+0.4) & Yb & Ygs<0.2 & Yg<0.2 & ... Ydiv<0.1 & Ydiv>-0.02 & ... Yg<0.1 & 
         ~(Ydiv./Yg<0.5 & ((Ym>0.9 & Yg>0.1 & Ydiv<0) | (Ym>0.95)) & Ym<1.2);
       %& ~Ywm;  
  Ybm  = smooth3(Ybm)>0.5;
  Ybm  = cat_vol_morph(Ybm,'o',1);
  % cortical GM 
  Ygm  = Ym<(GMth*0.3+0.7) & Ym>(CMth*0.6+0.4*GMth) & Yg<0.4 & Yb & Ydiv<0.4 & Ydiv>-0.3 & ~Ywm & ~Ybm & ~Ywm; % & (Ym-Ydiv*2)<GMth;  
  Ygm(smooth3(Ygm)<0.3)=0;
  % CSF
  Ycm  = Ym<(CMth*0.5+0.5*GMth) & Yg<0.1 & Yb & ~Ygm & dilmsk<-headradius*0.3; 
  Ycm  = smooth3(Ycm)>0.5;
  % head tissue
  Yhm  = Ym>max(mean([CMth,GMth]),Hth*0.2) & Ym<1.2 & Yg<0.8 & cat_vol_smooth3X(Yb,2)<0.1 & Ydiv<0.6 & Ydiv>-0.6;
  Yhm  = smooth3(Yhm)>0.5;

  %% refine
  %Ygm  = Ygm | (cat_vol_morph(Ywm,'d',3) & ~Ybm & ~Ywm & (Ym-Ydiv*2)<GMth & ...
  %   ~Ycm & smooth3((Ym + Yg)<(CMth*0.8+0.2*GMth))<0.5) & Ym>(CMth*0.9+0.1*GMth) & Ym<(GMth*0.2+0.8);;
  
  %% masking of the original values and local filtering
  stime = cat_io_cmd('  Filtering','g5','',verb,stime);
  fi   = round(max(3,min(resT3.vx_volr)*3)/3); 
  Ywm  = Ysrc .* Ywm; Ywm  = cat_vol_localstat(Ywm,Ywm>0,1,3); % PVE
  Ycm  = Ysrc .* Ycm; Ycm  = cat_vol_localstat(Ycm,Ycm>0,1,2);
  for i=1:fi-1, Ywm = cat_vol_localstat(Ywm,Ywm>0,2,1); end
  for i=1:fi-1, Ycm = cat_vol_localstat(Ycm,Ycm>0,2,1); end
  Ybm  = Ysrc .* Ybm; for i=1:fi, Ybm = cat_vol_localstat(Ybm,Ybm>0,2,1); end
  Ygm  = Ysrc .* Ygm; for i=1:fi, Ygm = cat_vol_localstat(Ygm,Ygm>0,2,1); end
  Yhm  = Ysrc .* Yhm; for i=1:fi, Yhm = cat_vol_localstat(Yhm,Yhm>0,2,1); end

  % estimate intensity difference bettween the tissues
  Ywmr = cat_vol_resize(Ywm,'reduceV',resT3.vx_volr,8,16,'meanm'); 
  Ybmr = cat_vol_resize(Ybm,'reduceV',resT3.vx_volr,8,16,'meanm'); 
  Ygmr = cat_vol_resize(Ygm,'reduceV',resT3.vx_volr,8,16,'meanm'); 
  Ycmr = cat_vol_resize(Ycm,'reduceV',resT3.vx_volr,8,16,'meanm'); 
  bmth = mean(Ybmr(Ybmr(:)>0 & Ywmr(:)>0))/mean(Ywmr(Ybmr(:)>0 & Ywmr(:)>0));
  gmth = mean(Ygmr(Ygmr(:)>0 & Ywmr(:)>0))/mean(Ywmr(Ygmr(:)>0 & Ywmr(:)>0));
  cmth = mean(Ycmr(Ycmr(:)>0 & Ywmr(:)>0))/mean(Ywmr(Ycmr(:)>0 & Ywmr(:)>0));
  Ywmr = cat_vol_resize(Ywm,'reduceV',resT3.vx_volr,16,8,'meanm'); 
  Yhmr = cat_vol_resize(Yhm,'reduceV',resT3.vx_volr,16,8,'meanm'); 
  hmth = mean(Yhmr(Yhmr(:)>0 & Ywmr(:)>0))/mean(Ywmr(Ywmr(:)>0 & Ywmr(:)>0));
  hmth = min(max(hmth,mean([cmth,gmth])),mean([gmth,1]));
  % if something failed use global thresholds
  if isnan(bmth), bmth = GMth; end
  if isnan(gmth), gmth = GMth; end
  if isnan(cmth), cmth = CMth; end
  if isnan(hmth), hmth = GMth; end
  clear Ywmr Ybmr Ygmr Yhmr Ycmr; 

  Yhm = Yhm .* (dilmsk>20);  % to avoid near skull tissue
  
  %% estimate bias fields
  stime = cat_io_cmd('  Bias correction','g5','',verb,stime);
  Ywi = sum( cat(4,Ywm,Ygm/gmth,Ybm/bmth,Ycm/cmth,Yhm/hmth),4) ./ sum( cat(4,Ywm>0,Ybm>0,Ygm>0,Ycm>0,Yhm>0),4 );
  [Ywi,resT2]  = cat_vol_resize(Ywi,'reduceV',resT3.vx_volr,min(4,min(resT3.vx_volr)*2),32,'meanm'); 
  for i=1:4, Ywi=cat_vol_localstat(Ywi,Ywi>0,2,1); end
  Ywi   = cat_vol_approx(Ywi,'nn',resT2.vx_volr,2);
  Ywi   = cat_vol_smooth3X(Ywi,4.*mean(vx_vol)); 
  Ywi   = cat_vol_resize(Ywi,'dereduceV',resT2);

  %% background noise
  stime = cat_io_cmd('  Background correction','g5','',verb,stime);
  %Ybc  = cat_vol_morph(smooth3(Ym<mean([BGth,CMth]) & Ym<CMth & Ygs<0.05 & ~Yb & dilmsk2>8)>0.5,'lo',3); 
  [Ybc,resT2] = cat_vol_resize(Ysrc .* Ybg,'reduceV',resT2.vx_volr,max(8,max(16,cat_stat_nanmean(resT2.vx_volr)*4)),16,'min'); 
  Ybc  = cat_vol_localstat(Ybc,Ybc>0,2,2);
  for i=1:1, Ybc  = cat_vol_localstat(Ybc,Ybc>0,2,1); end
  %Ybc2 = cat_vol_approx(Ybc,'nn',resT2.vx_volr,4); % no aproximation to correct only in the backgound!  
  %Ybc2 = cat_vol_smooth3X(Ybc2,4);
  Ybc  = cat_vol_smooth3X(Ybc,2);
  Ybc  = cat_vol_resize(Ybc,'dereduceV',resT2); 
  %Ybc2 = cat_vol_resize(Ybc2,'dereduceV',resT2); 

  %% back to original size
  stime = cat_io_cmd('  Final scaling','g5','',verb,stime);
  [Ywi,Ybc] = cat_vol_resize({Ywi,Ybc},'dereduceV',resT3); 
  %Ybc2 = cat_vol_resize({Ybc2},'dereduceV',resT3); 
  [Yg,Ygs]  = cat_vol_resize({Yg,Ygs},'dereduceV',resT3); 
  Yb   = cat_vol_resize(Yb,'dereduceV',resT3)>0.5; 
  Yp0 = cat_vol_resize(((Ywm>0)*3 + (Ygm>0)*2 + (Ybm>0)*2.3 + (Ycm>0)*1 + (Yhm>0)*2.7 + (Ybg>0)*0.5)/3,'dereduceV',resT3);
  Ysrc = Ysrco; clear Ysrco;

  %%  Final intensity scaling
  Ym   = (Ysrc - Ybc) ./ (Ywi - Ybc); % correct for noise only in background
 % Ym   = (Ysrc - Ybc) ./ (Ywi - Ybc2 + Ybc); % correct for noise only in background
  Wth  = single(cat_stat_nanmedian(Ym(Ygs(:)<0.2 & Yb(:) & Ym(:)>0.9))); 
  [WIth,WMv] = hist(Ym(Ygs(:)<0.1 &  Yb(:) & Ym(:)>GMth & Ym(:)<Wth*1.2),0:0.01:2);
  WIth = find(cumsum(WIth)/sum(WIth)>0.90,1,'first'); WIth = roundx(WMv(WIth),rf);  
  [BIth,BMv] = hist(Ym(Ym(:)<mean([BGth,CMth]) & Yg(:)<0.2),-1:0.01:2);
  BIth = find(cumsum(BIth)/sum(BIth)>0.02,1,'first'); BIth = roundx(BMv(BIth),rf);  
  Ym   = (Ym - BIth) ./ (WIth - BIth); 
  
  cat_io_cmd(' ','','',verb,stime); 
 
return
%=======================================================================

%=======================================================================
function [Ym,Yt,Ybg,WMth] = APP_initial_bias_correction(Ysrco,vx_vol,verb)
%% ---------------------------------------------------------------------
%  rough bias correction:
%  All tissues (low gradient) should have a similar intensity.
%  A strong smoothing of this approximation is essential to 
%  avoid anatomical filtering between WM and GM that can first 
%  be seen in overfitting of the subcortical structures!
%  However, this filtering will overcorrect head tissue with
%  a typical intensity around GM.
%    ds('l2','',0.5,Yo/WMth,Yg<0.2,Yo/WMth,Ym,80)
%  ---------------------------------------------------------------------
  rf = 10^6; 
  bfsmoothness = 3; 
  if verb, fprintf('\n'); end
  
  stime = cat_io_cmd('  Initialize','g5','',verb);
  msize = 222; %round(222 ./ max(size(Ysrco).*vx_vol) .* min(size(Ysrco).*vx_vol));  

  [Ysrc,resT3] = cat_vol_resize(Ysrco,'reduceV',vx_vol,min(1.2,cat_stat_nanmean(vx_vol)*2),msize,'meanm'); 

  % correction for negative backgrounds (MT weighting)
  WMth = roundx(single(cat_stat_nanmedian(Ysrc(Ysrc(:)>cat_stat_nanmean( ...
          Ysrc(Ysrc(:)>cat_stat_nanmean(Ysrc(:))))))),rf); 
  [BGth,BGv] = hist(Ysrc(Ysrc(:)<WMth*0.5)/WMth,min(Ysrc(:)/WMth):0.05:max(Ysrc(:)/WMth));
  BGth = find(cumsum(BGth)/sum(BGth)>0.05,1,'first'); BGth = roundx(BGv(BGth),rf); 

  Ysrc = Ysrc - BGth; Ysrco = Ysrco - BGth;
  Yg   = cat_vol_grad(Ysrc,resT3.vx_volr) ./ max(eps,Ysrc); 
  Ydiv = cat_vol_div(Ysrc,resT3.vx_volr) ./ Ysrc;

  WMth = roundx(single(cat_stat_nanmedian(Ysrc(Yg(:)<0.2 & Ysrc(:)>cat_stat_nanmean( ...
           Ysrc(Yg(:)<0.2 & Ysrc(:)>cat_stat_nanmean(Ysrc(:))))))),rf); 
  [BGth,BGv] = hist(Ysrc(Ysrc(:)<WMth*0.2)/WMth,min(Ysrc(:)/WMth):0.05:max(Ysrc(:)/WMth));
  BGth = find(cumsum(BGth)/sum(BGth)>0.05,1,'first'); BGth = roundx(BGv(BGth),rf); 
  BGth = max(BGth*WMth*4,WMth*0.2); % * 2 to get the CSF
  Ym   = (Ysrc - BGth) ./ (WMth - BGth);

  
  % background
  stime = cat_io_cmd('  Estimate background','g5','',verb,stime);
  Ybg = ((Yg.*Ym)<cat_vol_smooth3X(Ym,2)*1.2) & Ym>0.1; 
  Ybg([1,end],:,:)=0; Ybg(:,[1,end],:)=0; Ybg(:,:,[1,end])=0; Ybg(smooth3(Ybg)<0.5)=0;
  [Ybg,resT2] = cat_vol_resize(single(Ybg),'reduceV',resT3.vx_volr,2,32,'meanm'); 
  Ybg([1,end],:,:)=0; Ybg(:,[1,end],:)=0; Ybg(:,:,[1,end])=0; Ybg = Ybg>0.5;
  Ybg  = cat_vol_morph(Ybg,'lc',8);
  Ybg  = cat_vol_smooth3X(Ybg,2); 
  Ybg  = cat_vol_resize(Ybg,'dereduceV',resT2)<0.5;    
  BGth = roundx(mean(Ysrc(Ybg(:))),rf);
  Ym   = (Ysrc - BGth) ./ (WMth - BGth);
  
  %% first WM inhomogeneity with low tissue boundary (may include CSF > strong filtering for IXI175)
  stime = cat_io_cmd('  Initial correction','g5','',verb,stime);
  Yms  = cat_vol_smooth3X( min(2 .* ~Ybg,Ym .* (Ydiv>-0.2) .* ~Ybg .* (Ym>0.1)),16*mean(vx_vol));     % this map is to avoid CSF in the mask!
  Yms  = (Yms ./ mean(Yms(~Ybg(:)))) * WMth;
  Yms  = cat_vol_smooth3X( min(Yms*1.5 .* ~Ybg,Ysrc .* ~Ybg),16*mean(vx_vol));
  Yms  = (Yms ./ mean(Yms(~Ybg(:)))) * WMth;
  Yt   = Ysrc>max(BGth,Yms*0.3) & Ysrc<Yms*2 & Ysrc<WMth*(1+Yms/WMth*2) & Yg<0.9 & Ydiv<0.2 & ...
         Ydiv>-0.6 & smooth3(Ysrc./Yms.*Yg.*Ydiv<-0.2)<0.3 & ~Ybg; Yt(smooth3(Yt)<0.5)=0;
  Ywi  = (Ysrc .* Yt) ./ max(eps,Yt);  
  [Ywi,resT2] = cat_vol_resize(Ywi,'reduceV',resT3.vx_volr,cat_stat_nanmean(resT3.vx_volr)*2,32,'max'); 
  for i=1:1, Ywi = cat_vol_localstat(Ywi,Ywi>0,2,3); end % only one iteration!
  for i=1:4, Ywi = cat_vol_localstat(Ywi,Ywi>0,2,1); end
  Ywi  = cat_vol_approx(Ywi,'nn',resT2.vx_volr,4);
  Ywi  = cat_vol_smooth3X(Ywi,bfsmoothness.*mean(vx_vol)); % highres data have may stronger inhomogeneities 
  Ywi  = cat_vol_resize(Ywi,'dereduceV',resT2);    
  Ybc  = Ysrc./Ywi;
  WMt2 = roundx(cat_stat_nanmedian(Ybc(Yg(:)<0.2 & Ybc(:)>0.9)),rf); 
  Ywi  = Ywi * WMt2;
  
  %% background update
  stime = cat_io_cmd('  Refine background','g5','',verb,stime);
  Ybg = ((Yg.*(Ysrc./Ywi))<cat_vol_smooth3X(Ysrc./Ywi,2)*1.2) & Ysrc./Ywi>0.2; 
  Ybg([1,end],:,:)=0; Ybg(:,[1,end],:)=0; Ybg(:,:,[1,end])=0; Ybg(smooth3(Ybg)<0.5)=0;
  [Ybg,resT2] = cat_vol_resize(single(Ybg),'reduceV',resT3.vx_volr,2,32,'meanm'); 
  Ybg([1,end],:,:)=0; Ybg(:,[1,end],:)=0; Ybg(:,:,[1,end])=0; Ybg = Ybg>0.5;
  Ybg  = cat_vol_morph(Ybg,'lc',8);
  Ybg  = cat_vol_smooth3X(Ybg,2); 
  Ybg  = cat_vol_resize(Ybg,'dereduceV',resT2)<0.5;    

  %% second WM inhomogeneity with improved Yt with higher lower threshold (avoid CSF and less filtering)
  stime = cat_io_cmd('  Final correction','g5','',verb,stime);
  Yt   = Ysrc>max(BGth,Yms*0.3)  & Ysrc./Ywi>0.2 & Ysrc./Ywi<1.2 & Ysrc./Ywi<Yms/WMth*2 & Yg<0.9 & Ydiv<0.2 & Ydiv>-0.6 & ...
         smooth3(Ysrc./Yms.*Yg.*Ydiv<-0.1)<0.1 & ~Ybg; Yt(smooth3(Yt)<0.5)=0;
  Yt   = Yt | (~Ybg & Ysrc>BGth/2 & Ysrc>Yms*0.5 & Ysrc<Yms*1.2 & Ydiv./Yg<0.5 & ((Ysrc./Ywi>0.3 & Yg>0.1 & Ydiv<0) | (~Ybg & Ysrc./Ywi>0.6)) & Ysrc./Ywi<1.2); 
  Yt(smooth3(Yt)<0.7)=0;
  Ywi  = (Ysrc .* Yt) ./ max(eps,Yt);  
  [Ywi,resT2] = cat_vol_resize(Ywi,'reduceV',resT3.vx_volr,cat_stat_nanmean(resT3.vx_volr)*2,32,'max'); 
  for i=1:1, Ywi = cat_vol_localstat(Ywi,Ywi>0,2,3); end % only one iteration!
  for i=1:4, Ywi = cat_vol_localstat(Ywi,Ywi>0,2,1); end
  Ywi  = cat_vol_approx(Ywi,'nn',resT2.vx_volr,4);
  Ywi  = cat_vol_smooth3X(Ywi,bfsmoothness.*mean(vx_vol)); %.*mean(vx_vol)); % highres data have may stronger inhomogeneities 
  Ywi  = cat_vol_resize(Ywi,'dereduceV',resT2);    
  Ybc  = Ysrc./Ywi;
  WMt2 = roundx(cat_stat_nanmedian(Ybc(Yg(:)<0.2 & Ybc(:)>0.9)),rf); 
  Ywi  = Ywi * WMt2;

  
  %% BG inhomogeneity (important for normalization of the background noise)
  %[Ybc,Ygr,resT2] = cat_vol_resize({Ysrc./Ywi,Yg},'reduceV',resT3.vx_volr,cat_stat_nanmean(resT3.vx_volr)*4,16,'meanm'); 
  %Ybc  = cat_vol_morph(Ybc<BGth/WMth*2 & Ygr<0.05,'lc',2);
  %Ybc  = cat_vol_resize(smooth3(Ybc),'dereduceV',resT2)>0.5; 
  stime = cat_io_cmd('  Background correction','g5','',verb,stime);
  [Ybc,resT2] = cat_vol_resize(single(Ysrc .* Ybg),'reduceV',resT3.vx_volr,max(8,min(16,cat_stat_nanmean(resT3.vx_volr)*4)),16,'min'); 
  Ybc  = cat_vol_localstat(Ybc,Ybc>0,2,2);
  Ybc  = cat_vol_localstat(Ybc,Ybc>0,2,1);
  %Ybc  = cat_vol_approx(Ybc,'nn',resT2.vx_volr,4); % no aproximation to correct only in the background! 
  Ybc  = cat_vol_smooth3X(Ybc,4);
  Ybc  = cat_vol_resize(Ybc,'dereduceV',resT2); 


  %% back to original size
  stime = cat_io_cmd('  Final scaling','g5','',verb,stime);
  [Ywi,Ybc] = cat_vol_resize({Ywi,Ybc},'dereduceV',resT3); 
  Yg        = cat_vol_resize(Yg,'dereduceV',resT3); 
  [Yt,Ybg]  = cat_vol_resize({single(Yt),single(Ybg)},'dereduceV',resT3); Yt = Yt>0.5; Ybg = Ybg>0.5;
  Ysrc      = Ysrco; clear Ysrco;

  %% intensity normalization (Ybc is the average background noise)
  % in data with strong inhomogeneities (7T) the signal can trop below the noise level 
  Ym   = (Ysrc - min(Ybc/2,Ywi/20)) ./ (Ywi - min(Ybc/2,Ywi/20)); 
  Wth  = single(cat_stat_nanmedian(Ym(Yg(:)<0.2 & Ym(:)>cat_stat_nanmean( Ym(Yg(:)<0.2 & Ym(:)>cat_stat_nanmean(Ym(:))))))); 
  [WIth,WMv] = hist(Ym(Yg(:)<0.2 & Ym(:)>Wth*0.5 & Ym(:)<Wth*1.5),0:0.01:2);
  WIth = find(cumsum(WIth)/sum(WIth)>0.8,1,'first'); WIth = roundx(WMv(WIth),rf); 
  Ym   = Ym ./ WIth; 

  cat_io_cmd(' ','','',verb,stime); 
return
%=======================================================================


%=======================================================================
function Yg = cat_vol_grad(Ym,vx_vol)
% ----------------------------------------------------------------------
% gradient map for edge description
% ----------------------------------------------------------------------
  [gx,gy,gz] = cat_vol_gradient3(Ym); 
  Yg = abs(gx./vx_vol(1))+abs(gy./vx_vol(2))+abs(gz./vx_vol(3)); 
  %Yg = Yg ./ (Ym+eps);
return
%=======================================================================

%=======================================================================
function Ydiv = cat_vol_div(Ym,vx_vol)
% ----------------------------------------------------------------------
% Divergence helps to identify all gyri that should not be in the GM, but 
% helps to improve the WM. Divergence estimation is very memory intensive 
% so it is better to limit the resolution.
% ----------------------------------------------------------------------
  [Ymr,resT2] = cat_vol_resize(Ym,'reduceV',vx_vol,1.5,32);
  [gx,gy,gz]  = cat_vol_gradient3(max(1/3,Ymr)); 
  Ydivr = smooth3(divergence(gy./vx_vol(1),gx./vx_vol(1),gz./vx_vol(3))); clear gx gy gz Ymr;
  Ydiv  = cat_vol_resize(Ydivr,'dereduceV',resT2); 
return
%=======================================================================

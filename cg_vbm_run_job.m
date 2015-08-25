function cg_vbm_run_job(job,estwrite,tpm,subj)
    
    stime = clock;

    %% print current VBM release number and subject file
    A = ver; r = 0;
    for i=1:length(A)
        if strcmp(A(i).Name,'Voxel Based Morphometry Toolbox')
            r = str2double(A(i).Version);
        end
    end
    clear A 

    str  = sprintf('VBM12 r%d',r);
    str2 = spm_str_manip(job.channel(1).vols{subj},['a' num2str(70 - length(str))]);
    vbm_io_cprintf([0.2 0.2 0.8],'\n%s\n%s: %s%s\n%s\n',...
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

      if any(vx_vol>5)  % to high slice thickness
        error('VBM:cg_vbm_write:ToLowResolution', sprintf(...
             ['Voxel resolution has to be better than 3.5 mm in any dimention \n' ...
              'for save VBM preprocessing and a reasonable anatomical analysis! \n' ...
              'This image has got a resolution %0.2fx%0.2fx%0.2f mm%s. '], ... 
                vx_vol,char(179))); %#ok<SPERR>
      end
      if prod(vx_vol)>27  % to low voxel volume (smaller than 3x3x3 mm3)
        error('VBM:cg_vbm_write:ToHighVoxelVolume', ...
             ['Voxel volume has to be smaller than 10 mm%s (around 2x2x2 mm%s) to \n' ...
              'allow a save VBM preprocessing and reasonable anatomical analysis! \n' ...
              'This image has got a voxel volume of %0.2f mm%s. '], ...
              char(179),char(179),prod(vx_vol),char(179));
      end
      if max(vx_vol)/min(vx_vol)>8 % isotropy 
        error('VBM:cg_vbm_write:ToStrongIsotropy', sprintf(...
             ['Voxel isotropy (max(vx_size)/min(vx_size)) has to be smaller 8 to \n' ...
              'allow a save VBM preprocessing and reasonable anatomical analysis! \n' ...
              'This image has got a resolution %0.2fx%0.2fx%0.2f mm%s and a isotropy of %0.2f. '], ...
              vx_vol,char(179),max(vx_vol)/min(vx_vol))); %#ok<SPERR>
      end
    end


    % save original file name 
    for n=1:numel(job.channel) 
      job.channel(n).vols0{subj} = job.channel(n).vols{subj};
    end
    % noise-correction
    if job.vbm.sanlm && job.extopts.NCstr
        % for windows always disable multi-threading
        if ispc
            if (job.vbm.sanlm == 2) || (job.vbm.sanlm == 4)
                job.vbm.sanlm = job.vbm.sanlm - 1;
            end
        end

        switch job.vbm.sanlm
          case {1,3}, stime = vbm_io_cmd('NLM-Filter'); 
          case {2,4}, stime = vbm_io_cmd('NLM-Filter with multi-threading');
          case {5},   stime = vbm_io_cmd('Temporary NLM-Filter with multi-threading');
        end


        for n=1:numel(job.channel) 
            V = spm_vol(job.channel(n).vols{subj});
            Y = single(spm_read_vols(V));
            Y(isnan(Y)) = 0;
            switch job.vbm.sanlm
              case {1,3,5}, sanlmMex_noopenmp(Y,3,1,0); % use single-threaded version
              case {2,4},   sanlmMex(Y,3,1,0);          % use multi-threaded version
            end
            Vn = vbm_io_writenii(V,Y,'n','noise corrected','float32',[0,1],[1 0 0]);
            job.channel(n).vols{subj} = Vn.fname;
            clear Y V Vn;
        end

        fprintf('%4.0fs\n',etime(clock,stime));     
    end
    
    
    %% Interpolation
    % The interpolation can help to reduce problems for morphological
    % operations for low resolutions and strong isotropic images. 
    % Especially for Dartel a native resolution higher than the Dartel 
    % resolution helps to reduce normalization artifacts of the
    % deformation. Also this artifacts were reduce by the final smoothing
    % it is much better to avoid them.  
    Vt      = tpm.V(1); 
    vx_vold = min(cg_vbm_get_defaults('extopts.vox'),sqrt(sum(Vt.mat(1:3,1:3).^2))); clear Vt; % Dartel resolution 
    for n=1:numel(job.channel) 

      % prepare header of resampled volume
      Vi        = spm_vol(job.channel(n).vols{subj}); 
      vx_vol    = sqrt(sum(Vi.mat(1:3,1:3).^2));
      switch job.vbm.restype 
        case 'native'
          vx_voli  = vx_vol;
        case 'fixed', 
          vx_voli  = min(vx_vol ,job.vbm.resval(1) ./ ((vx_vol > (job.vbm.resval(1)+job.vbm.resval(2)))+eps));
          vx_voli  = max(vx_voli,job.vbm.resval(1) .* (vx_vol < (job.vbm.resval(1)-job.vbm.resval(2))));
        case 'best'
          vx_voli  = min(vx_vol ,job.vbm.resval(1) ./ ((vx_vol > (job.vbm.resval(1)+job.vbm.resval(2)))+eps));
        otherwise 
          error('cg_vbm_run_job:restype','Unknown resolution type ''%s''. Choose between ''fixed'',''native'', and ''best''.',restype)
      end
      vx_voli   = max(0.5,min(vx_vold,vx_voli)); % guarantee Dartel resolution
      
      
      % interpolation 
      if any( (vx_vol ~= vx_voli) )  
       
        stime = vbm_io_cmd(sprintf('Internal resampling (%4.2fx%4.2fx%4.2fmm > %4.2fx%4.2fx%4.2fmm)',vx_vol,vx_voli));
       
        Vi        = rmfield(Vi,'private'); 
        imat      = spm_imatrix(Vi.mat); 
        Vi.dim    = round(Vi.dim .* vx_vol./vx_voli);
        imat(7:9) = vx_voli .* sign(imat(7:9));
        Vi.mat    = spm_matrix(imat);

        Vn = spm_vol(job.channel(n).vols{subj}); 
        Vn = rmfield(Vn,'private'); 
        if ~(job.vbm.sanlm && job.extopts.NCstr)
          % if no noise correction we have to add the 'n' prefix here
          [pp,ff,ee] = spm_fileparts(Vn.fname);
          Vi.fname = fullfile(pp,['n' ff ee]);
          job.channel(n).vols{subj} = Vi.fname;
        end
        if job.vbm.sanlm==0
          [pp,ff,ee,dd] = spm_fileparts(Vn.fname); 
          Vi.fname = fullfile(pp,['n' ff ee dd]);
          job.channel(n).vols{subj} = Vi.fname;
        end
        vbm_vol_imcalc(Vn,Vi,'i1',struct('interp',6,'verb',0));

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

        obj.fwhm     = job.vbm.fwhm;
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
        obj.reg      = job.vbm.reg;
        obj.samp     = job.vbm.samp;              
       


        %% Initial affine registration.
        stime = vbm_io_cmd('Affine Registration'); 
 
        Affine  = eye(4);
        [pp,ff] = spm_fileparts(job.channel(1).vols{subj});
        Pbt = fullfile(pp,['brainmask_' ff '.nii']);
        Pb  = char(cg_vbm_get_defaults('extopts.brainmask'));
        Pt1 = char(cg_vbm_get_defaults('extopts.T1'));
        if ~isempty(job.vbm.affreg)
          
          warning off %#ok<WNOFF>
          
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

          % smooth source with 8mm
          VF1 = spm_smoothto8bit(VF,8);
         
          %stime = vbm_io_cmd('Initial Coarse Affine Registration'); 
          aflags     = struct('sep',8,'regtype',job.vbm.affreg,'WG',[],'WF',[],'globnorm',0);
          aflags.sep = max(aflags.sep,max(sqrt(sum(VG(1).mat(1:3,1:3).^2))));
          aflags.sep = max(aflags.sep,max(sqrt(sum(VF(1).mat(1:3,1:3).^2))));

          try
            spm_plot_convergence('Init','Coarse Affine Registration','Mean squared difference','Iteration');
          catch
            spm_chi2_plot('Init','Coarse Affine Registration','Mean squared difference','Iteration');
          end
          [Affine0, scale]  = spm_affreg(VG, VF1, aflags, eye(4)); Affine = Affine0; 

          % 4 mm
          try 
            aflags.WG  = spm_vol(Pb);
          catch
            pause(rand(1))
            aflags.WG  = spm_vol(Pb);
          end
          
          aflags.sep = aflags.sep/2;
          try
            spm_plot_convergence('Init','Coarse Affine Registration 2','Mean squared difference','Iteration');
          catch
            spm_chi2_plot('Init','Coarse Affine Registration 2','Mean squared difference','Iteration');
          end
          Affine1 = spm_affreg(VG, VF1, aflags, Affine, scale);   
          if ~any(isnan(Affine1(1:3,:))), Affine = Affine1; end
            
          clear VG1 VF1

          %fprintf('%4.0fs\n',etime(clock,stime));
        end
        
        
        % Fine Affine Registration with 3 mm sampling distance
        % Especially for non-human TPMs a brain mask is important
        % to avoid 'SingularMatrix' errors!
        %stime = vbm_io_cmd('Fine Affine Registration');
        spm_plot_convergence('Init','Fine Affine Registration','Mean squared difference','Iteration');
        Affine3 = spm_maff8(obj.image(1),obj.samp,obj.fudge,obj.tpm,Affine,job.vbm.affreg);
        if ~any(isnan(Affine3(1:3,:))), Affine = Affine3; end
        fprintf('%4.0fs\n',etime(clock,stime));
        
        warning on  %#ok<WNON>
        obj.Affine = Affine;
        
        
        
        %% SPM preprocessing 1
        
        stime = vbm_io_cmd('SPM preprocessing 1');
        
        % species
        % --------------------------------------------------------------
        if ~exist(Pbt,'file') && ~strcmp(job.vbm.species,'human')
          % Similar to the fine Affine registration this VBM12 
          % function requries brain masking for non-human subjects.
          % Old SPM segment will runs much simpler, but does not 
          % work directly for VBM12 due to the change data structure.
          Vb = spm_vol(Pb);
          vbm_vol_imcalc([VF,Vb],Pbt,'i2',struct('interp',6,'verb',0)); 
          obj.msk = Pbt;
        end
        % --------------------------------------------------------------
        
        warning off %#ok<WNOFF>
        res = spm_preproc8(obj);
        warning on  %#ok<WNON>

        if cg_vbm_get_defaults('extopts.debug')==2
          % save information for debuging and OS test
          [pth,nam] = spm_fileparts(job.channel(1).vols0{subj}); 
          tmpmat = fullfile(pth,sprintf('%s_%s_%s.mat',nam,'runjob','postpreproc8')); 
          save(tmpmat,'obj','res','Affine','Affine0','Affine1','Affine3');      
        end 
        
        if exist(Pbt,'file'), delete(Pbt); end
       
        fprintf('%4.0fs\n',etime(clock,stime));   
            

        try %#ok<TRYNC>
            [pth,nam] = spm_fileparts(job.channel(1).vols{subj});
            if job.vbm.sanlm>0
              nam = nam(2:end);
            end
        end
    end
    
    %% check contrast
    Tgw = [mean(res.mn(res.lkp==1)) mean(res.mn(res.lkp==2))]; 
    Tth = [
      ... mean(res.mn(res.lkp==6)) ... background
      max( min(mean(res.mn(res.lkp==3)) , max(Tgw)+abs(diff(Tgw))),min(Tgw)-abs(diff(Tgw)) ) ... % csf with limit for T2!
      mean(res.mn(res.lkp==1)) ... gm
      mean(res.mn(res.lkp==2)) ... wm 
    ];
    
    % inactive preprocessing of inverse images (PD/T2) 
    if cg_vbm_get_defaults('extopts.INV')==0 && any(diff(Tth)<=0)
      error('VBM:cg_vbm_write:BadImageProperties', ...
      ['VBM12 is designed to work only on highres T1 images.\n' ...
       'T2/PD preprocessing can be forced on your own risk by setting \n' ...
       '''vbm.extopts.INV=1'' in the vbm default file. If this was a highres \n' ...
       'T1 image than the initial segmentation seemed to be corrupded, maybe \n' ...
       'by alignment problems (check image orientation).']);    
    end
            
    %% Final iteration, so write out the required data.
    tc = [cat(1,job.tissue(:).native) cat(1,job.tissue(:).warped)];
    bf = job.bias;
    df = job.vbm.write;
    lb = job.label;
    jc = job.jacobian;
    res.stime = stime;
    res.image0 = spm_vol(job.channel(1).vols0{subj}); 
    cg_vbm_write(res, tc, bf, df, lb, jc, job.vbm, obj.tpm, job);

    % delete denoised/interpolated image
    if (job.vbm.sanlm && job.extopts.NCstr) || any( (vx_vol ~= vx_voli) ) 
      [pp,ff,ee] = spm_fileparts(job.channel(1).vols{subj});
      delete(fullfile(pp,[ff,ee]));
    end

return

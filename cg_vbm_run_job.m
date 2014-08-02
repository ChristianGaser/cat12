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

      if any(vx_vol>3.5)  % to high slice thickness
        error('VBM:cg_vbm_write:ToLowResolution', sprintf(...
             ['Voxel resolution has to be better than 3.5 mm in any dimention \n' ...
              'for save VBM preprocessing and a reasonable anatomical analysis! \n' ...
              'This image has got a resolution %0.2fx%0.2fx%0.2f mm%s. '], ... 
                vx_vol,char(179))); %#ok<SPERR>
      end
      if prod(vx_vol)>10  % to low voxel volume (smaller than 2x2x2 mm3)
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
      if sqrt(mean([2 2 2].^2))>8 % resolution RMS value
        error('VBM:cg_vbm_write:ToLowRMSResolution', sprintf(...
             ['Voxel RMS value has to be smaller 2 to allow a save VBM preprocessing and  \n' ...
              'reasonable anatomical analysis! \n' ...
              'This image has got a resolution %0.2fx%0.2fx%0.2f mm%s and a RMS value of %0.2f. '], ...
              vx_vol,char(179),sqrt(mean([2 2 2].^2)))); %#ok<SPERR>
      end
    end


    % noise-correction
    if job.warp.sanlm
        % for windows disable multi-threading
        if strcmp(mexext,'mexw32') || strcmp(mexext,'mexw64')
            job.warp.sanlm = min(1,job.warp.sanlm);
        end

        switch job.warp.sanlm
          case {1,3}, stime = vbm_io_cmd('NLM-Filter'); 
          case {2,4}, stime = vbm_io_cmd('NLM-Filter with multi-threading');
        end


        for n=1:numel(job.channel) 
            V = spm_vol(job.channel(n).vols{subj});
            Y = single(spm_read_vols(V));
            Y(isnan(Y)) = 0;
            switch job.warp.sanlm
              case {1,3}, sanlmMex_noopenmp(Y,3,1); % use single-threaded version
              case {2,4}, sanlmMex(Y,3,1);          % use multi-threaded version
            end
            Vn = vbm_io_writenii(V,Y,'n','noise corrected','float32',[0,1],[1 0 0],0);
            job.channel(n).vols{subj} = Vn.fname;
            clear Y V Vn;
        end

        fprintf('%4.0fs\n',etime(clock,stime));     
    end

    if estwrite % estimate and write segmentations            

        % 
        images = job.channel(1).vols{subj};
        for n=2:numel(job.channel)
            images = char(images,job.channel(n).vols{subj});
        end

        obj.image    = spm_vol(images);
        spm_check_orientations(obj.image);

        obj.fwhm     = job.warp.fwhm;
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
        obj.reg      = job.warp.reg;
        obj.samp     = job.warp.samp;              
        M = eye(4);


        %{
        %%  -------------------------------------------------------------
        %  Correct orientation, if the the AC is 2 SD outside of the
        %  center of mass of the major object.
        %  -------------------------------------------------------------
        % estimate average object intensity 
        % (typically a value something between GM and WM)
        V   = spm_vol(obj.image(1));
        vol = spm_read_vols(V);
        avg = mean(vol(:));
        avg = mean(vol(vol>avg));

        % don't use background values
        [x,y,z] = ind2sub(size(vol),find(vol>avg));
        com = [mean(x) mean(y) mean(z)];
        cSD = [std(x)  std(y)  std(z)];

        vmat = spm_imatrix(obj.image(1).mat);
        if any((vmat(1:3).*vmat(7:9))<(com - 2*cSD)) || any((-vmat(1:3).*vmat(7:9))>(com + 2*cSD))
          VG  = spm_vol(fullfile(spm('Dir'),'toolbox','OldNorm','T1.nii'));
          mat = spm_matrix([-com.*vmat(7:9) 0 0 0 vmat(7:9) 0 0 0]);
          M   = eye(4);
          M   = VG.mat\M*mat;
        end
        %}


        %% Initial affine registration.
        Affine  = eye(4);
        if ~isempty(job.warp.affreg),
            try
              VG = spm_vol(fullfile(spm('Dir'),'toolbox','OldNorm','T1.nii'));
            catch
              pause(rand(1))
              VG = spm_vol(fullfile(spm('Dir'),'toolbox','OldNorm','T1.nii'));
            end
            VF = spm_vol(obj.image(1));

            % smooth source with 8mm
            VF1 = spm_smoothto8bit(VF,8);

            % Rescale images so that globals are better conditioned
            VF1.pinfo(1:2,:) = VF1.pinfo(1:2,:)/spm_global(VF1);
            VG.pinfo(1:2,:)  = VG.pinfo(1:2,:)/spm_global(VG);

            %fprintf('Initial Coarse Affine Registration..\n');
            stime = vbm_io_cmd('Initial Coarse Affine Registration'); 
            aflags    = struct('sep',8, 'regtype',job.warp.affreg,...
                        'WG',[],'WF',[],'globnorm',0);
            aflags.sep = max(aflags.sep,max(sqrt(sum(VG(1).mat(1:3,1:3).^2))));
            aflags.sep = max(aflags.sep,max(sqrt(sum(VF(1).mat(1:3,1:3).^2))));

            spm_plot_convergence('Init','Coarse Affine Registration','Mean squared difference','Iteration');
            warning off;
            [Affine, scale]  = spm_affreg(VG, VF1, aflags,M);
            warning on;

            aflags.WG  = spm_vol(fullfile(spm('Dir'),'toolbox','FieldMap','brainmask.nii'));
            aflags.sep = aflags.sep/2;
            spm_plot_convergence('Init','Fine Affine Registration','Mean squared difference','Iteration');
            warning off;
            Affine = spm_affreg(VG, VF1, aflags, Affine, scale);
            warning on;
            fprintf('%4.0fs\n',etime(clock,stime));


            % Fine Affine Registration with 3 mm sampling distance
            stime = vbm_io_cmd('Fine Affine Registration');
            warning off;
            Affine = spm_maff8(obj.image(1),3,obj.fudge,  tpm,Affine,job.warp.affreg);
            warning on;
            fprintf('%4.0fs\n',etime(clock,stime)); 
        end;
        obj.Affine = Affine;


        stime = vbm_io_cmd('SPM-Preprocessing 1');
        warning off; 
        res = spm_preproc8(obj);
        warning on;
        fprintf('%4.0fs\n',etime(clock,stime));   

        try
            [pth,nam] = spm_fileparts(job.channel(1).vols{subj});
            if job.warp.sanlm>0
              nam = nam(2:end);
            end
            save(fullfile(pth,['vbm12_' nam '.mat']),'-struct','res', spm_get_defaults('mat.format'));
        end

    else % only write segmentations

        [pth,nam] = spm_fileparts(job.channel(1).vols{subj});
        if job.warp.sanlm>0
          nam = nam(2:end);
        end
        seg12_name = fullfile(pth,['vbm12_' nam '.mat']);

        if exist(seg12_name,'file')
            res = load(seg12_name);

            % check for spm version
            if ~isfield(res,'wp')
                error([fullfile(pth,['vbm12_' nam '.mat']) ' was not processed using SPM12. Use Estimate&Write option.']);
            end

            % load original used tpm, which is save in seg12.mat file
            try
                tpm    = spm_load_priors8(res.tpm);
            catch
                % or use default TPM
                fprintf('Original TPM image %s was not found. use default TPM image instead.\n',res.tpm(1).fname);
                for i=1:6
                    job.tissue(i).tpm = fullfile(spm('dir'),'tpm',['TPM.nii,' num2str(i)]);
                end
                tpm    = char(cat(1,job.tissue(:).tpm));
                tpm    = spm_load_priors8(tpm);
            end

            % use path of mat-file in case that image was moved
            [image_pth,image_nam,image_ext] = spm_fileparts(job.channel(1).vols{subj});
            res.image(1).fname = fullfile(image_pth, [image_nam, image_ext]);
        else
            error(['Can''t load file ' seg12_name]);  
        end
    end

    % Final iteration, so write out the required data.
    tc = [cat(1,job.tissue(:).native) cat(1,job.tissue(:).warped)];
    bf = job.bias;
    df = job.warp.write;
    lb = job.label;
    jc = job.jacobian;
    res.stime = stime;
    cg_vbm_write(res, tc, bf, df, lb, jc, job.warp, tpm, job);

    % delete denoised image
    if job.warp.sanlm>0
      delete(job.channel(1).vols{subj});
    end

return

function [out,outs] = cat_vol_mp2rage(job)
%cat_vol_mp2rage. Optimization of MP2Rage weighted scans for CAT12. 
% The function uses the unified segmenation of SPM to further optimize a 
% set of input images (MP2RAGE but also other T1-weighteds scans).
%
%  out = cat_vol_mp2rage(job)
%
%  job
%   .files            .. list of MP2Rage images
%   .biascorrection   .. biascorrection (0-no,*1-yes*,2-yes-extended) 
%   .skullstripping   .. skull-stripping 
%                        (0-no, 1-SPM, 2-optimized, *3-backgroundremoval*)
%   .intscale         .. general parameter to use logscale and intnorm 
%                        (for default GUI, 0-no, 1-yes, -1-unset)   
%   .logscale         .. use log/exp scaling for more equally distributed
%                        tissues (0-none, 1-log, -1-exp, inf-*auto*);
%   .intnorm          .. contrast normalization using the tan of GM normed
%                        values with values between 1.0 - 2.0 for light to 
%                        strong adaptiong (0-none, inf-*auot*)
%   .restoreLCSFnoise .. restore CSF (noise) values below zero (0-no,1-yes)   
%   .spm_preprocessing.. do SPM preprocessing 
%                        (0-no, 1-yes (if required), *2-always*)
%   .spm_cleanupfiles .. remove SPM files (0-no, *1-yes*)
%   .report           .. create report and xml/mat file (0-no, *1-yes*)
%   .prefix           .. filename prefix (strong with PARA for parameter
%                        depending naming, e.g. ... default, 'PARA_)
%   .verb             .. be verbose (0-no,*1-yes*,2-details)
%
%  out.files          .. resulting images (for SPM batch dependencies)
%  outs               .. structure of major results
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$


  def.files             = {};       % list of MP2Rage images
  def.biascorrection    = 1;        % biascorrection (0-no,1-light,2-extended) 
  def.skullstripping    = 3;        % skull-stripping (0-no, 1-SPM, 2-optimized, *3-backgroundremoval*)
  def.intscale          = -1;       % general parameter to use logscale and intnorm (for default GUI, 0-no, 1-yes, -1-unset )
  def.logscale          = inf;      % use log/exp scaling for more equally distributed
                                    % tissues (0-none, 1-log, -1-exp, inf-*auto*);
  def.intnorm           = -.5;      % contrast normalization using the tan of GM normed
                                    % values with values between 1.0 - 2.0 for light to 
                                    % strong adaptiong (0-none, 1..2-manuel, -0..-2-*auto*)
  def.restoreLCSFnoise  = 1;        % restore values below zero (lower CSF noise)    
  def.prefix            = 'PARA_';  % filename prefix (strong with PARA for parameter
                                    % depending naming, e.g. ... ) 
  def.spm_preprocessing = 2;        % do SPM preprocessing (0-no, 1-yes (if required), 2-always)
  def.spm_cleanupfiles  = 1;        % remove SPM files (0-no, *1-yes*)
  def.report            = 1;        % create a report (0-no, *1-yes*)
  def.verb              = 1;        % be verbose (0-no,*1-yes*,2-details)
  %def.outputBF          = 1;        % write correction map?
  %def.bloodvesselscorrection = 0;   % not implemented >> should be done by CAT 
  def.quicktest         = 0;        % low quality setting for quick test of general running of settings
  
  job = cat_io_checkinopt(job,def);


  [job.resdir,job.prefix] = spm_fileparts(job.prefix);
  

  

  % this is a super parameter for default users to control logscale and  
  % intnorm in a simple manner (do / not do correction)
  if job.intscale == 0
    job.logscale = 0; job.intnorm = 0; 
  elseif job.intscale == 1 % use defaults
    job.logscale = def.logscale; job.intnorm = def.intnorm; 
  end
  if isinf(job.logscale), lg = 'A'; else, lg = sprintf('%d',job.logscale); end
  if isinf(job.intnorm),  in = 'A'; else, in = sprintf('%d',job.intnorm);  end
  
  if 0 %job.biascorrection == 0 
    % RD20241010: no sure if this is useful
    cat_io_addwarning('Warning:NoBiasCorrection','No bias correction > no intensity scaling!')
    job.logscale = 0; 
    job.intnorm  = 0;
  end

  if job.quicktest
  % Quicktest setting for faster SPM processing to test the principle running for various parameter
    job.spm_preprocessing = 2;
    job.spm_cleanupfiles  = 0;
    job.verb              = 2; 
  end

  % update prefix
  if cat_io_contains(job.prefix,'PARA')
    job.prefix = strrep(job.prefix,'PARA',sprintf('MP2R_bc%d_lg%s_in%s_rn%d_ss%d',...
       job.biascorrection, lg, in, job.restoreLCSFnoise, job.skullstripping));
  end
  parastr = sprintf('(bc%d;lg%s;in%s;rn%d;ss%d)',...
    job.biascorrection, lg, in, job.restoreLCSFnoise, job.skullstripping);

  % prepare output
  out.files = job.files; 


  % TODO: 
  % * Bugs: 
  %   - CJV values in T2/PD MP2R
  %   - PD/T2 MP2R with to strong background remval 
  %   - intensity setting need further tests 
  %
  % * Challanging data with extrem bias field: 
  %   - 7T-TRT inv2-T1w data >> strong GM and even stronger thickness 
  %     understimation in CAT
  % * Improve bone correction in case of MP2R use BG intensity as minimum  
  %   and fat intensity as maximum
  % * Add QC measures (as we have the SPM segmentation)
  % * Bias correction warning if brain tissues show high variation or  
  %   corrections of the segmentation (closing) where done? 
  % * Report (make changes transparent)
  %    > colors to support faster review
  %    > add histogram !!!
  %    > colormap for images?
  %    > intensity peaks org / cor
  %    > resolution 
  % - variable/memory cleanup
  %
  % x What to do in case of worse results? 
  %   - Just use the original? 
  %   - No, problems should be visible to be fixed
  % x blood vessel correction ? 
  %    - No, should be done by PP, but may aboid precorrection. 
  %    - Maybe indirectly as part of the skull-stripping (reinclude as brain or not) 

  %#ok<*CLOCK,*DETIM>
    
  %% main processing 
  spm_progress_bar('Init',numel(job.files),'MP2Rage Optimization','Volumes Complete');
  for si = 1:numel(job.files)
    clear Yseg Ym Yo; 

    % be verbose
    stime = clock; 
    if job.verb
      [pp,ff,ee] = spm_fileparts(job.files{si});
      fprintf('%4d) %80s: ',si,spm_str_manip( fullfile(pp,[ff,ee]), 'a80'));
    end

  
    %% estimate if MP2R to set up the SPM processing 
    %  MP2R is here defined as an image with high intensity background, ie. 
    %  there is no low intensity background 
    %  In low-int background data a threshold bit below the median will 
    %  separate the object (head/brain) from the background, whereas in 
    %  an MP2R this is not the case and the thresholded regions covers 
    %  the whole image.
    Vm = spm_vol( job.files{si} ); 
    Ym = spm_read_vols( Vm ); 
    [Ymr,res] = cat_vol_resize(Ym,'reduceV',1,4,16,'meanm'); 
    Ymo    = cat_vol_morph(smooth3(Ymr) > 0.9*median(Ymr(:)),'lc',4) & Ymr~=0; 
    res.isMP2R = sum(Ymo(:) | Ymr(:)==0) > 0.98*numel(Ymo); 
    clear Ymo Ymr Ym Vm; 
    

    %% check if processed data is available 
    SPMdata = [
      spm_file(job.files(si),'suffix','_seg8','ext','.mat')
      spm_file(job.files(si),'prefix','m' ,'ext','.nii')
      spm_file(job.files(si),'prefix','c1','ext','.nii')
      spm_file(job.files(si),'prefix','c2','ext','.nii')
      spm_file(job.files(si),'prefix','c3','ext','.nii')
      spm_file(job.files(si),'prefix','c4','ext','.nii')
      spm_file(job.files(si),'prefix','c5','ext','.nii')
      ];
    SPMdatae = zeros(numel(SPMdata),1);
    for fi = 1:numel(SPMdata)
      SPMdatae(fi) = exist( SPMdata{fi} , 'file'); 
    end
       
    if job.spm_preprocessing==2 || ... allways run SPM segmentation 
      (job.spm_preprocessing==1 && any(~SPMdatae)) ... only if the segments are missing

      %% SPM segmentation
      if job.verb>1, fprintf('\n'); end
      stime2 = cat_io_cmd('      SPM segmentation','g5','',job.verb-1);

      spmiter = max(1,2 - 1*res.isMP2R);
      % if no bias-correction is used also the additional iterations are not useful
      if job.biascorrection <= 1, spmiter = 1; end 
      % for quick tests we not need full quality ?
      if job.quicktest > 1, spmiter = 1; end
      for i = 1:spmiter
        % *** biasreg 1e-6 caussed problems in the head tissue *** 
        % *** even 2 mm are not always better ***
        if job.verb>1 && i>1, fprintf('.'); end
        if job.biascorrection && i~=1
          matlabbatch{1}.spm.spatial.preproc.channel.vols   = spm_file(job.files(si),'prefix','m');
        else
          matlabbatch{1}.spm.spatial.preproc.channel.vols   = job.files(si);
        end
        if i==1
          % the first iteration should be smooth
          if res.isMP2R % only one iteration
            matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
            matlabbatch{1}.spm.spatial.preproc.warp.samp        = 3;
          else
            matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 120;
            matlabbatch{1}.spm.spatial.preproc.warp.samp        = 6;
          end
        elseif i==spmiter
          % final iterations could do further corrections
          matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 30;
          matlabbatch{1}.spm.spatial.preproc.warp.samp        = 3; % 3 - (spmiter>2 && spmiter==i); % 2 mm only for 3 iters
        else % intermediate 
          matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
          matlabbatch{1}.spm.spatial.preproc.warp.samp        = 4; % 3 - (spmiter>2 && spmiter==i); % 2 mm only for 3 iters
        end
        matlabbatch{1}.spm.spatial.preproc.channel.biasreg    = 1e-4; % biasreg 1e-6 caussed problems on the bottom of the TPM 
        matlabbatch{1}.spm.spatial.preproc.channel.write      = [0 job.biascorrection>0]; % write bias corrected image
        ngaus = [1 1 2-res.isMP2R 3 4 2+res.isMP2R]; % use 3 to have more stable MP2R backgrounds
        for ti = 1:6
          matlabbatch{1}.spm.spatial.preproc.tissue(ti).tpm    = ...
            {fullfile(spm('dir'),'tpm',sprintf('TPM.nii,%d',ti))};
          matlabbatch{1}.spm.spatial.preproc.tissue(ti).ngaus  = ngaus(ti);
          matlabbatch{1}.spm.spatial.preproc.tissue(ti).native = ...
            [((ti<6 & spmiter==i) | ti==2 | ti==5) 0]; 
          matlabbatch{1}.spm.spatial.preproc.tissue(ti).warped = [0 0];
        end
        matlabbatch{1}.spm.spatial.preproc.warp.mrf           = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.cleanup       = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.reg           = [0 0.001 0.5 0.05 0.2];
        matlabbatch{1}.spm.spatial.preproc.warp.affreg        = 'mni';
        matlabbatch{1}.spm.spatial.preproc.warp.fwhm          = 2 * (i==1);
        matlabbatch{1}.spm.spatial.preproc.warp.write         = [0 0];
        matlabbatch{1}.spm.spatial.preproc.warp.vox           = NaN;
        matlabbatch{1}.spm.spatial.preproc.warp.bb            = [NaN NaN NaN; NaN NaN NaN];

        if job.quicktest
          matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
          matlabbatch{1}.spm.spatial.preproc.channel.biasreg  = 1e-4;
          matlabbatch{1}.spm.spatial.preproc.warp.samp        = 5; 
        end
        
        % run SPM 
        sii = 0; spmfailed = 0;
        while sii<3
          sii = sii + 1; 
          try
            evalc('spm_jobman(''run'',matlabbatch)'); 
          catch
            spmfailed = 1;
            continue
          end
  
          % use bias corrected in next iteration
          [pp,ff,ee] = spm_fileparts(job.files{si}); mfile = fullfile(pp,[ff,ee]); 
          if job.biascorrection  &&  i > 1
            movefile( spm_file(mfile,'prefix','mm'), spm_file(mfile,'prefix','m') );
          
            [pp,ff,ee] = spm_fileparts(job.files{si}); mfile = fullfile(pp,[ff,ee]); 
            for ti = 1:6
              if exist(  spm_file(mfile,'prefix',sprintf('c%dm',ti)) , 'file' )
                movefile( spm_file(mfile,'prefix',sprintf('c%dm',ti)), ...
                  spm_file(mfile,'prefix',sprintf('c%d',ti)) );
              end
            end
          end

          % basic test 
% test spm8 values          
          if exist( spm_file(mfile,'prefix','c5') ,'file')
            Vseg(5) = spm_vol(spm_file(mfile,'prefix','c5'));  
            Yseg{5} = single( spm_read_vols(Vseg(5)) ); 
            if std(Yseg{5}(Yseg{5}(:)>0)) > 0
              matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = ...
                max(60,round(matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm / 1.5,-1));
              matlabbatch{1}.spm.spatial.preproc.warp.samp        = ...
                max(3, round(matlabbatch{1}.spm.spatial.preproc.warp.samp * .8,1));  
              continue
            else
              break 
            end
          end
% optimize cls parameters
% .. in later iterations test for fatsup value
        end
        if spmfailed
          cat_io_cprintf('err','SPM failed - check origin');
          continue; 
        end 

        %% basic informations
        spm8       = load(spm_file(job.files{si},'suffix','_seg8','ext','.mat'));
        res.isT1   = (spm8.mn(spm8.lkp==2)*spm8.mg(spm8.lkp==2)) > ...
                     (spm8.mn(spm8.lkp==1)*spm8.mg(spm8.lkp==1)); 
        res.isMP2R = (spm8.mn(spm8.lkp==6)*spm8.mg(spm8.lkp==6)) > ...
          min( 0.5 * (spm8.mn(spm8.lkp==2)*spm8.mg(spm8.lkp==2)) , ...
               2.0 * (spm8.mn(spm8.lkp==3)*spm8.mg(spm8.lkp==3)) );
        res.fatsup = ((spm8.mn(spm8.lkp==5 & spm8.mg'>.2)*spm8.mg(spm8.lkp==5 & spm8.mg'>.2)) / ...
                      (spm8.mn(spm8.lkp<3)*spm8.mg(spm8.lkp<3)) ) < .3;
        res.fatsup = (mean(spm8.mn(spm8.lkp==5 & spm8.mg'>.25)) / mean(spm8.mn(spm8.lkp<3))) < .9;
        vx_vol     = sqrt(sum(spm8.image.mat(1:3,1:3).^2)); 
        


        %%
        if job.biascorrection>1 && ~res.isMP2R && i~=spmiter %&& i>1
        % Using only SPM for bias correction was not successful in super 
        % extrem cases as the SPM model was stable even the bias was still 
        % included. Hence, we use the most important segments the WM and 
        % head to do a very quick/simple/smooth additional correction. 
        % Skip the first iteration as it was worse with the SPM8 estimates. 
        
          % load SPM bias correct image and the WM and head segment
          Pm = spm_file(job.files{si},'prefix','m');
          Vm = spm_vol(Pm); 
          Ym = single( spm_read_vols(Vm) );
          
          Yseg = cell(1,6);
          Pseg = cell(1,6); Tth = zeros(1,6); 
          for ci = [1,2,5]
            Pseg{ci} = spm_file(job.files{si},'prefix',sprintf('c%d',ci),'ext','.nii'); 
            if exist(Pseg{ci},'file')
              Vseg(ci) = spm_vol(Pseg{ci}); %#ok<AGROW> 
              Yseg{ci} = single( spm_read_vols(Vseg(ci)) ) & Ym~=0; 
            end
            Tth(ci)  = spm8.mn(spm8.lkp==ci) * spm8.mg(spm8.lkp==ci);
          end
          if res.isT1 && std(Yseg{5}(Yseg{5}(:)>0))==0 % correct SPM bug
            Yseg{5}(Ym<.1*Tth(2)) = 0;
          end
          Yd = cat_vol_div(Ym./Tth(2),vx_vol,2);

          %% fatsup
          [Ymr,Yhd,Ywm,resV] = cat_vol_resize( {Ym,Yseg{5} & smooth3(Ym)>.2*Tth(2),Yseg{2}},'reduceV',vx_vol,1.2,4,'meanm'); %clear Ypx#
          Ywm  = cat_vbdist( single( smooth3(Ywm) ) , Ymr>0.15*Tth(2), resV.vx_volr ); 
          Yhd  = cat_vol_morph(Yhd & Ywm>5 & Ywm<30,'ldo',2,vx_vol);
          Yhd  = cat_vbdist( single( smooth3(1-Yhd) ) , Yhd>0, resV.vx_volr ); 
          Ygt  = cat_vol_pbtp( single(3-(smooth3(Yhd>0))) , Yhd , 0*Yhd );
          Yhdppr = Yhd./max(eps,Ygt); Yhdppr(Yhdppr==1 & smooth3(Yhdppr)<.5 | Yhdppr>1) = 0;
          isfatsup0 = [ mean( Ymr( Yhdppr(:)>0.5) ) mean( Ymr( Yhdppr(:)>0 & Yhdppr(:)<0.5) ) ]; % inner
          isfatsup = isfatsup0(1) / isfatsup0(2) < 1.2; 
          Yhdpp = cat_vol_resize(Yhdppr,'dereduceV',resV);
          if isfatsup, cat_io_cprintf([0 .5 0],'_'); else, cat_io_cprintf('warn','^'); end
          res.fatsup = isfatsup; 

          %% use the WM and head segment to estimate another smooth bias corrected imags
          mth  = max(min(Tth(1:2)), min( Tth(5) , mean(Tth(1:2)))); 
          Yg   = cat_vol_grad(Ym)./Ym; 
          Ypx  = (Yseg{2}>.125)*Tth(2) + ... 
                 cat_vol_morph( smooth3((Yseg{5})>.2 & (Ym)>0.05*Tth(2))>.7,'lo',3)*mth;
          Ypx( smooth3(Ypx>0)<.1 | Ym==0 | (Yg>.5 & Yseg{2}<.1)) = 0;  
          Ypx(Ypx==0 & Yhdpp>0.1 & Yhdpp<.5 & Ym>min(Tth(1:2))) = mth;
          Ympx = (Ypx~=0) .* cat_vol_localstat(Ym,Ypx>0,1,2+res.isT1); 
          if isfatsup==0
            [Ymr,Ypxr,resV] = cat_vol_resize( {max(Ym,Ympx),Ypx},'reduceV',vx_vol,1.5,4,'meanm'); %clear Ypx
            Ywi0 = cat_vol_approx(Ymr ./ Ypxr); 
            Ywi0 = cat_vol_smooth3X(Ywi0,16);
            Ywi0 = cat_vol_resize(Ywi0,'dereduceV',resV);
            Ywi0 = Ywi0 ./ cat_stat_nanmedian(Ywi0(Yseg{2}>.5));
            %Ypx( Ypx>0 & (Ym ./ Ywi0) > Tth(2)*0.9 & (Ym ./ Ywi0) < Tth(2)*0.9 & Yseg{5}>.1) = 0;
            Ypx( Ypx>0 & (Ym ./ Ywi0) > Tth(2)*0.9 &  Yseg{5}>.1) = max(spm8.mn(spm8.lkp==5));
            Ypx(Ypx==0 & Yhdpp>.5 & Ym>.5*Tth(2)) = max(spm8.mn(spm8.lkp==5));
          end
          %%
          [Ymr,Ypxr,resV] = cat_vol_resize( {max(Ym,Ympx),Ypx},'reduceV',vx_vol,1.5,4,'meanm'); %clear Ypx
          Ywr  = cat_vol_localstat(Ymr ./ Ypxr,Ypxr>0,2,2+res.isT1); 
          Ywir = cat_vol_approx(Ywr); clear Ywr Ymr; 
         % Ymxr = cat_vol_smooth3X(Ypxr~=0,8); clear Ypxr; 
         % Ywir = cat_vol_smooth3X(Ywir,2).*Ymxr + (1-Ymxr).*cat_vol_smooth3X(Ywir,2); 
          Ywir = cat_vol_smooth3X(Ywir,2 + 2*(i==1));
          Ywi  = cat_vol_resize(Ywir,'dereduceV',resV); clear Ywir;
          Ywi  = Ywi ./ cat_stat_nanmedian(Ywi(Yseg{2}>.5));
          
          % quantify bias field to interupt unnecessary loops?
          BF(i) = std(Ywi(:)) / min(abs(diff(Tth(1:2)/max(Tth(1:2)) ))); 

          %% write the result   
          spm_write_vol(Vm,Ym ./ Ywi);
          clear Ym Yseg Ywi

% document values for XML report          
        end
      end
      
    else
      stime2 = clock; 
      % if SPM data is incomplete but we are not allow to do it
      if any(~SPMdatae)
        cat_io_addwarning('Warning:BadSPMdata',['Cannot find all required preprocessed SPM files \\n' ...
            '(*seg8.mat, m*.nii, c#*.nii). Continue with next subject!']);
        continue
      end
    end




    %% load SPM mat (with tissue thresholds)
    if etime(clock,stime2) < 10 
      if job.verb>1, fprintf('\n'); end % no 
      stime2 = cat_io_cmd('      Load segmentation','g5','',job.verb-1);
    else
      stime2 = cat_io_cmd('      Load segmentation','g5','',job.verb-1,stime2);
    end
    spm8 = load(spm_file(job.files{si},'suffix','_seg8','ext','.mat'));
    
    % load original for skull-stripping
    if isfield(job,'ofiles')
      Po = job.ofiles{si};
      Vo = spm_vol(Po); 
      Yo = single( spm_read_vols(Vo));
    else 
      Po = job.files{si}; 
      Vo = spm_vol(Po); 
      Yo = single( spm_read_vols(Vo));
    end

    % load bias corrected 
% ############# use the original if no bias corrected was written? ########## 
    if job.biascorrection
      Pm = spm_file(job.files{si},'prefix','m');
    else
      Pm = job.files{si};
    end
    Vm = spm_vol(Pm); 
    Ym = single( spm_read_vols(Vm));
    vx_vol = sqrt(sum(Vm.mat(1:3,1:3).^2)); 
    
    % load segmentation (for brain masking)
    Yseg = zeros(Vm.dim,'single');
    Pseg = cell(1,6); Tth = zeros(1,6); Sth = Tth; 
    for ci = 1:5
      Pseg{ci}       = spm_file(job.files{si},'prefix',sprintf('c%d',ci)); 
      Vseg(ci)       = spm_vol(Pseg{ci}); 
      Yseg(:,:,:,ci) = single( spm_read_vols(Vseg(ci)) );
      Ymsk           = cat_vol_morph( Yseg(:,:,:,ci)>.9 ,'e') & Ym~=0; 
      Tth(ci)        = cat_stat_nanmedian(Ym(Ymsk(:))); 
      if isnan(Tth(ci))
        Ymsk           = cat_vol_morph( Yseg(:,:,:,ci)>.5 ,'e') & Ym~=0; 
        Tth(ci)        = cat_stat_nanmedian(Ym(Ymsk(:))); 
      end
      if isnan(Tth(ci))
        Ymsk           = Yseg(:,:,:,ci)>.5 & Ym~=0; 
        Tth(ci)        = cat_stat_nanmedian(Ym(Ymsk(:))); 
      end
      Sth(ci)        = cat_stat_nanstd(Ym(Ymsk(:))); 
    end
    clear Ymsk Vseg;
    Yseg(:,:,:,6) = 1 - cat_stat_nansum( Yseg, 4);  
    Tth(6) = cat_stat_nanmedian(Ym(Yseg(:,:,:,6) & Ym~=0));
    res.isMP2R = Tth(6) > min(Tth(2)/2,Tth(3)*2);                  % MP2Rage
    res.isT1   = Tth(3) < Tth(1) & Tth(1) < Tth(2);  % T1 defintion  
    if 0 % ~res.isT1
      cat_io_cprintf('err','No T1w contrast or bad segmenation!\n')
      %continue
    end
    if any( isnan(Tth) )
      cat_io_cprintf('err','Bad contrast or segmenation!\n')
    end

   

%% %##########################
% detect critical segmentations / skull-strippings
% >> stop case processing and continue with next subject
% >> don't repair! but maybe suggest to check and to what aparameter could
%    be changed
% Case 1 - Buchert 125849) failed GM/WM segmentation but BET is ok
% Case 2 - Buchert 150601) extrem low contrast (remove/ignore/bullshit in=bullshit out)  
% Case 3 - Buchert 155553) low contrast (high-res, low-freq-noise) but good segmentation
% Case 4 - Buchert 155802) failed bias correction ...
% Yp0  = Yseg(:,:,:,1)*2 + Yseg(:,:,:,2)*3 + Yseg(:,:,:,3);
    iCon = @(x) (exp(x) - 1) / (exp(1)-1);

% ######################## WMHs ?! #############  
% ### errors with BVs and MAs ... test for local high variance in tissues
% #############

    
    if res.isMP2R
      norm2  = @(Yx, Ycmm,Ywmm) (Yx - median(Yx(Ycmm(:)))) ./ ...
        ( median(Yx(Ywmm(:))) - median(Yx(Ycmm(:)))) * 2/3 + 1/3;
      norm2i = @(Yx,Yy,Ycmm,Ywmm) (Yx - 1/3) .* ... 
        ( median(Yy(Ywmm(:))) - median(Yy(Ycmm(:))) )*2/3 + median(Yy(Ycmm(:)));
    else
      norm2  = @(Yx,   Ycmm,Ywmm) (Yx - median(Yx(Ycmm(:))) + std(Yx(Ycmm(:))) ) ./ ...
        ( median(Yx(Ywmm(:))) - median(Yx(Ycmm(:))) + std(Yx(Ycmm(:))) );
      norm2i = @(Yx,Yy,Ycmm,Ywmm) (Yx - 1/3) .* ...
        ( median(Yy(Ywmm(:))) - median(Yy(Ycmm(:))) + std(Yx(Ycmm(:))) ) + median(Yy(Ycmm(:))) - std(Yx(Ycmm(:))) ;
    end
    if res.isT1
      exp1   = exp(1) - 1; 
      Ywmm   = cat_vol_morph( Yseg(:,:,:,2)>.9 ,'e'); 
      Ycmm   = smooth3( (Yseg(:,:,:,3)>.95 | ( (Yseg(:,:,:,6)>.95) .* (Tth(6) < Tth(3))) ) & Ym/spm8.mn(2)<.5)>.5; 
      Tthm   = norm2( Tth , 3, 2); 
      Ymm    = norm2( Ym, Ycmm, Ywmm); log1str = 'Y'; 
      for i = 1:3
        TX  = real(...
              [ norm2(  exp(  norm2(Tthm, 3, 2)*2/3+1/3 - 1) / exp1 , 3, 2);
                norm2(        norm2(Tthm, 3, 2)*2/3+1/3 , 3, 2); 
                norm2( log2(  norm2(Tthm, 3, 2)*2/3+1/3 + 1) , 3, 2); 
                norm2( log(  (norm2(Tthm, 3, 2)*2/3+1/3) * exp1 + 1), 3, 2);
                norm2( log10((norm2(Tthm, 3, 2)*2/3+1/3) * 9 + 1), 3, 2) ]);
        [mx,cci] = min( abs(TX(1:5) - 0.7)); 
        Tthm = TX(cci,:);
      
        if     cci == 1
          Ymm2 = exp  (  norm2(Ymm,  Ycmm, Ywmm)*2/3+1/3 - 1) / exp1;    lstr = 'exp';   %% without / exp1 ???
        elseif cci == 2
          Ymm2 =         norm2(Ymm,  Ycmm, Ywmm)*2/3+1/3;                lstr = 'non';
        elseif cci == 3
          Ymm2 = log2 (  norm2(Ymm,  Ycmm, Ywmm)*2/3+1/3 + 1);           lstr = 'log2';
        elseif cci == 4
          Ymm2 = log  ( (norm2(Ymm,  Ycmm, Ywmm)*2/3+1/3) * exp1 + 1);   lstr = 'log';
        elseif cci == 5
          Ymm2 = log10( (norm2(Ymm,  Ycmm, Ywmm)*2/3+1/3) * 9 + 1);      lstr = 'log10';
        end
        if cci == 2 || mx < 0.02; break; end
        log1str = sprintf('%s(%s)',lstr,log1str);
      
        if res.isMP2R
          Ymm  = norm2( real(Ymm2) , Ycmm , Ywmm);
        else
          Ymm  = norm2( real(Ymm2) , Yseg(:,:,:,6)>.9 & Ym<Tth(2)*0.3 , Ywmm);
        end
        Tthm   = norm2( Tth , 3, 2); 
      end
      Ym2  = real(Ymm)  * spm8.mn(2);
      Tthm = Tthm * spm8.mn(2);
      Ym2(Yo==0)=0; 
      clear Ymm 
    else
      %Ymm  = norm2( Ym, Yo==0 | Ym<Tth(2), Yseg(:,:,:,3)>.5); 
      Ym2  = Ym;
      Tthm = Tth; 
      Ym2(Yo==0)=0; 
      clear Ymm 
    end



    %% get (MP2R) background
    [Ybgr,resV] = cat_vol_resize( Yseg(:,:,:,6)>.1 & Ym<1.5*Tth(2),'reduceV',vx_vol,2,4,'meanm');
    Ybgr = cat_vol_morph(cat_vol_morph(Ybgr,'ldo',2),'l',[10 0.2]); 
    Ybgr = cat_vol_smooth3X(Ybgr,2); 
    Ybg  = cat_vol_resize(Ybgr,'dereduceV',resV);
    Ybg  = Ybg>.5  & smooth3(Yo)~=0; 
 


    % estimate if fat supression is used
    %  get tissue surrounding the upper skull to estimate if fat supression 
    [Ymr,resV] = cat_vol_resize(  Ym/Tth(2),'reduceV',vx_vol,2,4,'meanm');
    Ysegr = nan([resV.sizeTr,size(Yseg,4)],'single'); 
    for ti = 1:size(Yseg,4)
      Ysegr(:,:,:,ti) = cat_vol_resize( Yseg(:,:,:,ti),'reduceV',vx_vol,2,4,'meanm');
    end

    % estimate distance maps to identify the upper skull
    Ybgdistr = cat_vbdist( Ybgr , sum(Ysegr(:,:,:,1:3) ,4) <.5, resV.vx_volr); 
    Ybdistr  = cat_vbdist( sum(Ysegr(:,:,:,1:3) ,4) , Ybgr<.5 , resV.vx_volr); 
    Ysdistr  = cat_vbdist( sum(Ysegr(:,:,:,1:4) ,4) , Ybgr<.5 , resV.vx_volr); 
    Yuskinr  = Ysegr(:,:,:,5)>.5 & Ybdistr<20 & Ybgdistr<20 & Ybdistr>Ysdistr & ...
               cat_vol_morph(1 - Ysegr(:,:,:,5)>.2,'dc',10,resV.vx_volr); 
    Yhdppr    = min(Ysdistr,Ybgdistr) ./ (Ybgdistr/2 + Ysdistr/2) .* Yuskinr; 
  %  clear Ysegr Ybdistr
    
    % final evaluation 
    %{
    fatsup = [ median( Ymr(Yhdppr(:)<1/2 & Yhdppr(:)>0) ) > median( Ymr(Yhdppr(:)>1/2) ) , ... low-int central area
                   abs(median( Ymr(Yhdppr(:)<1/2 & Yhdppr(:)>0 & Ybgdistr(:)<Ysdistr(:) ) ) - ... symetric sides
                     median( Ymr(Yhdppr(:)<1/2 & Yhdppr(:)>0 & Ybgdistr(:)>=Ysdistr(:)) )) < .1 , ...
                   ( spm8.mn(spm8.lkp==5) * spm8.mg(spm8.lkp==5) / Tth(2) ) > 1 ]; 
   % clear Ymr Ybgdistr Ysdistr
    res.fatsup = sum(fatsup)>=1; 
    %}
    Yhdpp    = cat_vol_resize(Yhdppr,'dereduceV',resV);



    %% extended bias correction
    if job.biascorrection > 1
      stime2 = cat_io_cmd('      Bias correction','g5','',job.verb-1,stime2);

      %  Use also the GM area to get the closes WM value by a maximum operation (i.e. T1 only). 
      %  We are not using the CSF as we are not knowing much about it and and cannot trust the segmentation too much.  
      Yg = cat_vol_grad(Ym2)./Ym2; Ygs = smooth3(Yg); 
      if 0
        %% old version

        % define bias correction field for each tissue  
        Ywi = cat_vol_localstat(Ym2,cat_vol_morph(Yseg(:,:,:,2)>.5,'lc') ,1,2+res.isT1);
        Ywi = cat_vol_median3(Ywi,Ywi>0,Ywi>0);
        Ywi(Ybg | Yo==0) = Tth(2); % avoid overcorrections outside the brain
        % remove outlier voxels       res.fatsup = sum(fatsup)>=2;  
        Ywiv  = cat_vol_localstat(Ywi,Ywi~=0,1,4);
        wivth = cat_stat_nanmedian(Ywiv(Ywiv(:)~=0));
        Ywi(Ywiv>(wivth*2)) = 0; Ywi(smooth3(Ywi~=0)<.5) = 0;  
        clear Ywiv wivth; 

        % mean here to handle noise .. .* (Ym2>0 | Ym2<max(Tth)*1.2)
        [Ym2r,resV] = cat_vol_resize(Ym2 .* (Ygs<.2),'reduceV',vx_vol,1.5,4,'meanm'); 
        [Ygi1,Ygi2] = cat_vol_resize({ ...
          single(smooth3(cat_stat_nansum(Yseg(:,:,:,1:2),4))>.5 & Ygs<.5 ), ...& ~Ysc
          single(smooth3(cat_stat_nansum(Yseg(:,:,:,1:2),4))>.5 & Ygs<.5 )},'reduceV',vx_vol,1.5,4,'meanm');
        Ygi1  = cat_vol_localstat(Ym2r,Ym2r~=0 & Ygi1>.5,4,2 + res.isT1); % closer to WM
        Ygi2  = cat_vol_localstat(Ym2r,Ym2r~=0 & Ygi2>.5,2,2 + res.isT1); % more distant to WM
        Ygi1(Ygi2>0 & Ygi1==0) = Ygi2(Ygi2>0 & Ygi1==0);
        Ygi   = cat_vol_approx( Ygi1 ); 
        % clear Ym2r Ygi1 Ygi2
        Ygi   = cat_vol_resize(Ygi,'dereduceV',resV); % .* (Yseg(:,:,:,1)>.2);
        % adaption for high noise bias in low intensity regions
        Ygi   = Ygi*.5 + 0.5*(iCon(Ygi / spm8.mn(2))*spm8.mn(2)); 
        Ywi(Ygi>0 & Ywi==0) = Ygi(Ygi>0 & Ywi==0);
        Ywi   = cat_vol_approx( Ywi ); 
      else
        %% new version 

        % create values for brain and muscles 
        Ypx  = Yseg(:,:,:,1)*Tthm(1) + Yseg(:,:,:,2)*Tthm(2) + Yseg(:,:,:,3)*Tthm(3)*.9;
       
        if res.isMP2R
          % in MP2R only add background and ignore other things
          Ypx(Ybg & Yo~=0) = cat_stat_nanmean( Ym2(Ybg(:) & Yo(:)~=0));
          Yw   = Ym2./Ypx .* (Ypx>0 & (Ygs<.4 | Yseg(:,:,:,5)>0)); Yw(Yw>2)=0;% clear Ypx
        else
          % in MPR we need to use the skull to handle extreme cases
          % depending on the fat suppression ...
          Ypx  = Yseg(:,:,:,1)*Tthm(1) + Yseg(:,:,:,2)*Tthm(2) + Yseg(:,:,:,3)*Tthm(3)*.9;
          
          Ypx = Ypx + (Yseg(:,:,:,5)>.125) * mean(Tthm(1)) .* (cat_vol_morph(Ygs<0.5,'o') & ~Ybg & Yo~=0); 
          %Ypx  = (Yseg{2}>.125)*Tth(2) + ... 
          %       cat_vol_morph( smooth3((Yseg{5})>.2 & (Ym)>0.05*Tth(2))>.7,'lo',3)*mth;
          Ypx( smooth3(Ypx>0)<.1 | Ym==0 | (Yg>.5 & Yseg(:,:,:,2)<.1)) = 0;  
          Ypx(Ypx==0 & Yhdpp>0.1 & Yhdpp<.5 & Ym>min(Tth(1:2))) = mean(Tthm(1));
          %%
          if res.fatsup
            Ylb = (smooth3(Ypx>0)<0.1  &  smooth3(Yseg(:,:,:,5))>0.5 ); 
            Ypx = max(Ypx, (Ylb & Ym>Tth(3)) .* mean(Tthm(1:2)) );
            Ypx = max(Ypx, (Ylb & smooth3(Yhdpp)>Ym/3 & Ym>Tth(3)) .* min(1,1 - Yhdpp) * mean(Tthm(1:2))); % two skin layer
          else
            %%
            Ypx = Ypx + (Ypx==0 & Yhdpp>.1 & Ym>Tth(3) & Yg>0.1) .* (Yhdpp .* max(spm8.mn(spm8.lkp==5))); 
          
            Ympx = (Ypx~=0) .* cat_vol_localstat(Ym,Ypx>0,2,2+res.isT1); 
            %{
            Ywi0 = cat_vol_approx(Ympx ./ Ypx); 
            Ywi0 = cat_vol_smooth3X(Ywi0,32);
            Ywi0 = Ywi0 ./ cat_stat_nanmedian(Ywi0(Yseg(:,:,:,2)>.5));
            Ypx( cat_vol_morph( Ypx>0 & (Ym ./ Ywi0) > Tth(2)*.9 & (Ym ./ Ywi0) < Tth(2)*1.2 & Yseg(:,:,:,5)>.1 ,'d' )) = 0;
            Ypx( cat_vol_morph( Ypx>0 & (Ym ./ Ywi0) > Tth(2)*1.2 &  Yseg(:,:,:,5)>.1 ,'d' )) = ...
              max(spm8.mn(spm8.lkp==5));
            %}
            [Ymr,Ypxr,resV] = cat_vol_resize( {max(Ym,Ympx),Ypx},'reduceV',vx_vol,1.5,4,'meanm'); %clear Ypx
            Ywi0 = cat_vol_approx(Ymr ./ Ypxr); 
            Ywi0 = cat_vol_smooth3X(Ywi0,16);
            Ywi0 = cat_vol_resize(Ywi0,'dereduceV',resV);
            Ywi0 = Ywi0 ./ cat_stat_nanmedian(Ywi0(Yseg(:,:,:,2)>.5));
            Ypx( (Ym ./ Ywi0) > Tth(2)*0.9 & Yseg(:,:,:,5)>.1) = max(spm8.mn(spm8.lkp==5));
            Ypx( Yhdpp>.4 & Ym>.5*Tth(2)) = max(spm8.mn(spm8.lkp==5));

          end
          Yw   = Ym2./Ypx .* (Ypx>0 & (Ygs<.4 | Yseg(:,:,:,5)>0) & ~Ybg); Yw(Yw>2)=0;% clear Ypx
        end
        
        %Ypx  = Ypx .* cat_vol_morph(Ypx>0 & ~Ybg,'l',[10 .1]); 
        Yw   = cat_vol_localstat(Yw,Yw~=0,2,2 + res.isT1); 
        %
        [Ywi,Ywm,resV] = cat_vol_resize({Yw,single(Yw>0)},'reduceV',vx_vol,2,4,'meanm'); 
        Ywi(Ywm<.25) = 0; clear Ywm;  
        Ywi  = cat_vol_approx(Ywi,'rec');
        Ywi  = cat_vol_smooth3X(Ywi,2);
        Ywi  = cat_vol_resize(Ywi,'dereduceV',resV);
        % extended smoothing
        Ybs  = cat_vol_smooth3X(sum(Yseg(:,:,:,1:3),4),8);
        Ywis = cat_vol_smooth3X(Ywi,4);
        Ywi  = Ywi .* Ybs + (1-Ybs) .* Ywis; clear Ybs Ywis;
        % scaling 
        Ywi  = Ywi ./ median(Ywi(Yseg(:,:,:,2)>.9)); 
      end
      Ym3   = Ym2./Ywi;  
      Ywi   = Ywi .* median(Ym3(Yseg(:,:,:,2)>.5));
      clear Ygi Ym3


      %% mix tissues and approximate bias field
      if res.isMP2R
        Ygw = sum(Yseg(:,:,:,2),4)>.5; 
        Yw = Ywi / mean(Ywi(Ygw(:))) * mean(Ym2(Ygw(:))); 
      else 
        %% for the inpaint method the use of GM information was less optimal 
        % - first estimation for background 
        Yi = Ywi; %max(Ywi,Ybi .* ( cat_vol_morph(Yseg(:,:,:,5)>.5,'l') & (Ym./spm8.mn(3))>.5) * spm8.mn(2) ); 
        Yw = cat_vol_inpaint(Yi,10,40,2,1); Yw = Yw*0.9 + 0.1 * iCon(Yw./spm8.mn(2)).*spm8.mn(2);
        Yw = Yw ./ cat_stat_nanmedian(Yw(Yseg(:,:,:,2)>.95)) .* cat_stat_nanmedian(Ym2(Yseg(:,:,:,2)>.95));
        % - second estimation for 
        Yi = max(Ywi,Yw .* ( cat_stat_nansum(Yseg(:,:,:,6),4)>.5 )); 
        Yw = cat_vol_inpaint(Yi,10,6 / job.biascorrection,2,1);  
        Yw = Yw*0.9 + 0.1 * iCon(Yw./spm8.mn(2)).*spm8.mn(2);
        %Yw = Yw ./ cat_stat_nanmedian(Yw(Yseg(:,:,:,2)>.95)); % .* cat_stat_nanmedian(Ym2(Yseg(:,:,:,2)>.95));
        
      end
  
      % apply bias correction
      Ymm  = Ym2 ./ Yw;
    else
      % no bias correction
      Ymm  = Ym2; 
    end  
    Ymbc = Ymm;
  




    %% general exp/log scaling 
    if res.isMP2R && res.isT1
      norm2  = @(Yx, Ycmm,Ywmm) (Yx - median(Yx(Ycmm(:)))) ./ ...
        ( median(Yx(Ywmm(:))) - median(Yx(Ycmm(:)))) * 2/3 + 1/3;
      norm2i = @(Yx,Yy,Ycmm,Ywmm) (Yx - 1/3) .* ... 
        ( median(Yy(Ywmm(:))) - median(Yy(Ycmm(:))) ) * 2/3 + median(Yy(Ycmm(:)));
    else
      norm2  = @(Yx,   Ycmm,Ywmm) (Yx - median(Yx(Ycmm(:))) + std(Yx(Ycmm(:))) ) ./ ...
        ( median(Yx(Ywmm(:))) - median(Yx(Ycmm(:))) + std(Yx(Ycmm(:))) );
      norm2i = @(Yx,Yy,Ycmm,Ywmm) (Yx - 1/3) .* ...
        ( median(Yy(Ywmm(:))) - median(Yy(Ycmm(:))) + std(Yx(Ycmm(:))) ) + median(Yy(Ycmm(:))) - std(Yx(Ycmm(:))) ;
    end
    if res.isT1  
      if job.logscale>0 || job.intnorm>0
        stime2 = cat_io_cmd('      Scale intensities','g5','',job.verb-1,stime2);
      end
   %Ymm = Ym ./ Yw; % only for debugging
      exp1 = exp(1) - 1; 
      if job.logscale == 1 
        Ymm = log(  Ymm/spm8.mn(2) * exp1 + 1) * spm8.mn(2);
      elseif job.logscale == 2
        Ymm = log2( (Ymm/spm8.mn(2)) + 1) * spm8.mn(2);
      elseif job.logscale == 10
        Ymm = log10( (Ymm/spm8.mn(2)) * 9 + 1) * spm8.mn(2);
      elseif job.logscale == -1
        Ymm = (exp(  Ymm/spm8.mn(2) ) - 1)  * spm8.mn(2) / exp1;
      elseif isinf(job.logscale) % auto
        % we typically expect the GM value at about 60% of the WM peak
        opt    = 0.7; % BG
        
        % some images/classes that we need later
        Ygm = Yseg(:,:,:,1)>.1 & ~cat_vol_morph( Yseg(:,:,:,1)>.5 ,'o',4); 
        Ywm = cat_vol_morph( Yseg(:,:,:,2)>.9 ,'e'); 
  
        Ywmm   = Ywm; %Yseg(:,:,:,2)>.95;
        Ycmm   = smooth3( (Yseg(:,:,:,3)>.95 | ( (Yseg(:,:,:,6)>.95) .* (Tth(6) < Tth(3))) ) & Ym/spm8.mn(2)<.5)>.5; 
  
        %%
        %Ymm = norm2(Ymm,Ycmm,Ywmm); 
        res.logscaleres = ''; 
        res.log2str = 'Y';
        for opti = 1:3
    % ################ local adaption or mixture of values?
          for ii = 1:2
            if ii == 1, cli = 1:5; else, ccm = abs(cc - opt); [ccmin,cci] = min( ccm ); cli = cci; end
            for cci = cli
              if     cci == 1, Ymm2 = exp  (  norm2(Ymm, Ycmm, Ywmm)*2/3+1/3 - 1) / exp1;     lstr = 'exp';
              elseif cci == 2, Ymm2 =         norm2(Ymm, Ycmm, Ywmm)*2/3+1/3;                 lstr = 'non';
              elseif cci == 3, Ymm2 = log2 (  norm2(Ymm, Ycmm, Ywmm)*2/3+1/3 + 1);            lstr = 'log2';
              elseif cci == 4, Ymm2 = log  ( (norm2(Ymm, Ycmm, Ywmm)*2/3+1/3) * exp1 + 1);    lstr = 'log';
              elseif cci == 5, Ymm2 = log10( (norm2(Ymm, Ycmm, Ywmm)*2/3+1/3) * 9 + 1);       lstr = 'log10';
              end
              if ii == 1 
                Ymm2     = norm2i(Ymm2,Ymm2,Ycmm,Ywmm); 
                Ymmt     = norm2(Ymm2,Ycmm,Ywmm); 
                cc(cci) = median( Ymmt(Ygm(:)) ); clear Ymmt
              end
            end
          end
          %Ymm = norm2i(Ymm2,Ymm,Ycmm,Ywmm); 
          Ymm = norm2(real(Ymm2),Ycmm,Ywmm); 
          
     %  ########   ccm = max(0,ccm + [0.04 0 0.01 0.02 0.04] - 0.04); 
  
          ccistr = {'exponential','orgiginal','log2','log','log10'};
          scstr  = sprintf('%12s-scaling (exp|org|log2|log|log10 = %s\b)',ccistr{cci},sprintf('%0.2f|',ccm)); 
          if job.verb > 1 && opti==1, fprintf('\n'); end 
          if job.verb > 1, cat_io_cprintf('g5','        %s\n',scstr); end
          if cci == 2 || ccmin < 0.02; break; end
          res.logscaleres = [res.logscaleres '; ' scstr];
        end
        res.log2str = sprintf('%s(%s)',lstr,res.log2str);
        if job.verb > 1, cat_io_cmd(' ','g5','',job.verb-1); end
        %res.logscale    = {ccistr;ccm};
        res.logscaleres = scstr;
      else
        Ygm = Yseg(:,:,:,1)>.1 & ~cat_vol_morph( Yseg(:,:,:,1)>.5 ,'o',4); 
        Ywm = cat_vol_morph( Yseg(:,:,:,2)>.9 ,'e'); 
      end
      if 0% job.logscale>0 || job.intnorm>0
        stime2 = cat_io_cmd(sprintf('      Scale intensities (%s)',ccistr{cci}),'g5','',job.verb-1,stime2); %#ok<UNRCH>
      end
  
      
     
      %% scale intensities 
      %  This is the critical part. 
      % log .. optimize the tissue peaks 
      % fx  .. scaling function to reduce the 
      % Ygs .. add some left-sided noise of CSF values
      
      %fx  = @(x,y,z) min(20,max(0,tan(((x-.5) * z + 0.5) * pi - pi/2)/pi/y/z + .2 + 0.5/y));
      fx2 = @(x,con,sqr) min(inf,max(-inf, ( abs(tan( max(-pi/2, min(pi/2, x * pi * con)))).^sqr .* sign(x)) ));
    
      if job.intnorm ~= 0 && res.isMP2R
        if isinf(job.intnorm) || job.intnorm < 0
          if  isinf(job.intnorm), intnorm = 0.5; else, intnorm = abs(job.intnorm); end
          %Ym2 = fx( log( Ym/spm8.mn(2)* 1.71 + 1),1.4,.3); % this was worse
          %Ym2 = fx2(  (Ymm - median(Ymm(Ygm(:))) ) / spm8.mn(2) ,1.8,0.6);  %works ok 
          %Ym2 = fx2(  (Ymm - median(Ymm(Ygm(:))) ) / spm8.mn(2) ,.6,1.4); % stronger
          Ym2 = fx2(  (Ymm - median(Ymm(Ygm(:))) ) / spm8.mn(2) ,0.1,1 + intnorm/2 .* .3 * std(Ymm(Ygm(:))) / std(Ymm(Ywm(:))) ); % stronger
        else
          Ym2 = fx2(  (Ymm - median(Ymm(Ygm(:))) ) / spm8.mn(2) ,0.1,job.intnorm); % manual
        end
        Ywmm = Yseg(:,:,:,2)>.95; Ycmm = smooth3(Yseg(:,:,:,3)>.95 & Ym/spm8.mn(2)<.5)>.5; 
        Ym2 = (Ym2 - median(Ym2(Ycmm(:)))) ./ (median(Ym2(Ywmm(:))) - median(Ym2(Ycmm(:))));
        Ym2 = Ym2 * 2/3 + 1/3 .* (Ymm>0);
        clear Ywmm Ycmm;
      else
        %Ym2 = (Ymm / spm8.mn(2)) * 2/3 + 1/3 .* (Ym>0); 
        if res.isMP2R
          Ym2 = norm2(Ymm,Ycmm & (Yseg(:,:,:,3)>.5),Ywmm); 
        else 
          Ym2 = norm2(Ymm, Ybg & Ymm<Tth(2)*0.3 , Ywmm);
          % scale by BG and WM 
       %   Ybg   = cat_vol_morph( Ym/Tth(2) < .3 & Yseg(:,:,:,6)>.5 ,'e');
       %   Ymstd = cat_vol_localstat(single(Ym / Tth(2)),Ybg,1,4);
       %   noise = cat_stat_nanmean(Ymstd(Ybg(:))); 
  
       %   Ym2 = (norm2(Ymm,Ycmm  & (Yseg(:,:,:,3)>.5),Ywmm) - 1/3 + 4*noise ) * 1 / (2/3 + 8*noise);  
          clear Ymstd; 
  %        Ym2 = (norm2(Ymm,Ycmm  & (Yseg(:,:,:,6)>.5),Ywmm) - 1/3) * 3/2; 
        end
      end
    else
      Ym2 = norm2(Ym2,Yo==0 | Yo<Tth(2)/4,Yseg(:,:,:,3)>.5); 
    end



    %% restore left-side CSF noise (value below zero that were just cutted)
    if job.restoreLCSFnoise && res.isT1
      if res.isMP2R
        stime2 = cat_io_cmd('      Restore CSF noise','g5','',job.verb-1,stime2);
        
        % define noise
        pnoise = 0.03; % percentage of noise
        Ymsk = Ym~=0 & Yseg(:,:,:,3)>0.8; % Ym !
        Ygs  = Ymsk .* cat_vol_smooth3X( Yseg(:,:,:,3) .^2 ,2) .* ...
          pnoise .* randn(size(Ym2)) .* min(1,1 + (Ym2 - 1/3)*3).^2; 
        % add noise
        Ym2(Ymsk) = Ym2(Ymsk) + Ygs(Ymsk);
        clear Ygs Ymsk;
      %else % inactive because it is only an developer option
      %  stime2 = cat_io_cmd('      No MP2Rage - no CSF noise restoration.','g8','',job.verb-1,stime2);
      end
    end



    
    %% additional skull-stripping
    if job.skullstripping == 2 || ( job.skullstripping == 3 && res.isMP2R)
      stime2 = cat_io_cmd('      Optimized skull-stripping','g5','',job.verb-1,stime2);

      %% denoising based correction
      % In MP2Rage, background/bone/head tissues are quite noisy and strongly 
      % corrected by the SANLM. Hence, denoised regions are often background
      % and should be not part of the GM class

      % get noisy regions
      Yngw = cat_vol_morph( smooth3( cat_stat_nansum(Yseg(:,:,:,1:2),4) ) > .99 ,'do',3,vx_vol);
      Yngw = smooth3( Yngw ) < .5;
      
      % gradient to get the local noise before/after correction (Yog/Ymg)
      if exist('Yo','var')
        Ymm  = Ym; 
        Yog  = cat_vol_grad(Yo);
        Ymg  = cat_vol_grad(Ym);
      else
        Ymm  = single(Ym) + 0; cat_sanlm(Ymm,1,3); 
        Yog  = cat_vol_grad(Ym);
        Ymg  = cat_vol_grad(Ymm);
      end
      Ymm  = max(-1,min(10,Ymm)); 
      Ybgx = cat_vol_morph( abs(Ym)<10*eps, 'd',1); 
      Yog(Ybgx) = mean(Yog(:)); Ymg(Ybgx) = mean(Yog(:));

      % get correction area
      Ynog = max( Ybgx , Yngw .* (Ymm/spm8.mn(2)) .* ((Yog ./ Ymg) .* abs(Yog-Ymg) ./ spm8.mn(2)) ); %clear Yog Ymg Ymm; 
      Ynog = Yseg(:,:,:,1) .* smooth3(Ynog);
      % do correction
      Yseg(:,:,:,4) = Yseg(:,:,:,4) + Ynog; 
      Yseg(:,:,:,1) = Yseg(:,:,:,1) - Ynog; 
      clear Yngw %Ynog


      % cleanup
      Ygw = smooth3( cat_stat_nansum(Yseg(:,:,:,1:2),4) ) > .9;
      Ygw = single(cat_vol_morph(Ygw,'do',3,vx_vol)); 
      Ygw(cat_stat_nansum(Yseg(:,:,:,1:3),4)<.2) = nan;  
      [~,Yd] = cat_vol_downcut(Ygw,Ym/spm8.mn(3) .* (1-sum( Yseg(:,:,:,5:6), 4)),0.01);
      Ygw2 = Ym/spm8.mn(2)<1.2 & smooth3(Yd<.5)>.5; 


      % edge detection to improve skull-stripping
      %Yg  = cat_vol_grad(Ym/spm8.mn(2));
      %Yd  = cat_vol_div(Ym/spm8.mn(2));
      % optimize skull-stripping
      Yb  = smooth3( cat_stat_nansum( Yseg(:,:,:,1:3),4) ); 
      Yb  = max(Yb,cat_vol_morph(Yb > .5 ,'ldo',1,vx_vol));
      Yb  = max(Yb,cat_vol_morph(Yb > .5 ,'ldc',2,vx_vol));
      Yb  = max(Yb,Ygw2);
      Yb  = smooth3( Yb )>.8; 

    elseif job.skullstripping
      % simple SPM skull-stripping as CSF+GM+WM 
      stime2 = cat_io_cmd('      SPM skull-stripping','g5','',job.skullstripping && job.verb>1,stime2);
      Yb  = smooth3( sum( Yseg(:,:,:,1:3),4) )>.5; 

    else
      Yb  = true(size(Ym)); 

    end
  
    if job.skullstripping == 1 || job.skullstripping == 2
      % full classical skull-stripping
      Ym2 = Ym2 .* Yb; 
    elseif job.skullstripping == 3 && res.isMP2R
      % background-stripping and skull modification 
      % we keep here some minimum amount of noise and also some skull values
      Ym2 = max(eps,Ym2 .* smooth3( max( (.4*cat_vol_smooth3X(Yseg(:,:,:,4),2).^4), ...
         max(0.06, smooth3( 1 - Ybg - Yseg(:,:,:,4)) ))) - Yseg(:,:,:,4)*.05  +  ...
         (Ybg*0.02 + 0.04 .* Ybg .* randn(size(Ybg),'single'))); 
    end
    % restore defacing values
    Ym2(Yo==0) = 0; 

    
    %% evaluation
    % evaluate also this ... 
    % * compare the difference between Yp0 and Ym2 => RMSE
    % * compare the histogram of the tissues => RMSE
    % * run iterative optimization ...
    if res.isT1
      Yp0   = Yseg(:,:,:,1)*2 + Yseg(:,:,:,2)*3 + Yseg(:,:,:,3);
    else
      Yp0   = Yseg(:,:,:,1)*2 + Yseg(:,:,:,2)*1 + Yseg(:,:,:,3)*3;
    end
    % intensity-based measures
    % - the correction should optimize the coefficient of joint variation (CJV) 
    fcjvx = @(x,y) ( cat_stat_nanstd(x(round(y(:))==1)) +  cat_stat_nanstd(x(round(y(:))==2)) +  cat_stat_nanstd(x(round(y(:))==3)) ) ./ ...
                   (cat_stat_nanmean(x(round(y(:))==1)) + cat_stat_nanmean(x(round(y(:))==2)) + cat_stat_nanmean(x(round(y(:))==3)));
    fcjv  = @(x,y) ( cat_stat_nanstd(x(round(y(:))==2)) +  cat_stat_nanstd(x(round(y(:))==3)) ) ./ ...
                   (cat_stat_nanmean(x(round(y(:))==2)) + cat_stat_nanmean(x(round(y(:))==3)));
    res.CJVorg      = fcjv(Yo,Yp0);
    res.CJVcor      = fcjv(Ym2,Yp0);
    res.CJVCGWorg   = fcjvx(Yo,Yp0);
    res.CJVCGWcor   = fcjvx(Ym2,Yp0);
    % Peak width similarity
    res.sdCGWorg    = [ cat_stat_nanstd(  Yo(round(Yp0(:)*3)==3)) cat_stat_nanstd(  Yo(round(Yp0(:))==2))  cat_stat_nanstd(  Yo(round(Yp0(:))==3))] / spm8.mn(2);
    res.sdCGWcor    = [ cat_stat_nanstd( Ym2(round(Yp0(:)*3)==3)) cat_stat_nanstd( Ym2(round(Yp0(:))==2))  cat_stat_nanstd( Ym2(round(Yp0(:))==3))];
    res.sdCGorg     = std(res.sdCGWorg(2:3));
    res.sdCGcor     = std(res.sdCGWcor(2:3));
    % gradient test
    Ygo   = cat_vol_grad(Yo)   ./ max(eps,Yo);
    Ygm   = cat_vol_grad(Ymbc) ./ max(eps,Ymbc); 
    Ygw   = max(0,1 - Ygo  - Ygm); 
    res.gradbias_tot = cat_stat_nanmean( (Ygm(:) - Ygo(:)) .*  Ygw(:)); 
    res.gradbias_imp = cat_stat_nanmean( max(0,Ygm(:) - Ygo(:) ) .* Ygw(:) ); 
    res.gradbias_wor = cat_stat_nanmean( max(0,Ygo(:) - Ygm(:) ) .* Ygw(:) ); 
    clear Ygo ym Ygw; 
    
    % volumetric measures 
    % - the resuling map should fit better to the existing segmentation 
    res.vol_TIV     = cat_stat_nansum( Yp0(:)>.5)  .* prod(vx_vol) / 1000;
    res.vol_abs_CGW = [ cat_stat_nansum(round(Yp0(:))==1) cat_stat_nansum(round(Yp0(:))==2) cat_stat_nansum(round(Yp0(:))==3) ] .* prod(vx_vol) / 1000;
    res.vol_rel_CGW = res.vol_abs_CGW ./ res.vol_TIV; 
    res.vol_abs_CGWorg = [ 
      cat_stat_nansum(Yb(:) & Yo(:)<spm8.mn(1) ) ...
      cat_stat_nansum(Yb(:) & Yo(:)>spm8.mn(1) & Yo(:)<spm8.mn(2) ) ...
      cat_stat_nansum(Yb(:) & Yo(:)>spm8.mn(2) ) ] .* prod(vx_vol) / 1000;
    res.vol_rel_CGWorg = res.vol_abs_CGWorg ./ res.vol_TIV; 
    res.vol_abs_CGWcor = [ cat_stat_nansum(round(Yb(:).*Ym2(:)*3)==1) cat_stat_nansum(round(Yb(:).*Ym2(:)*3)==2) cat_stat_nansum(round(Yb(:).*Ym2(:)*3)==3) ] .* prod(vx_vol) / 1000;
    res.vol_rel_CGWcor = res.vol_abs_CGWcor ./ res.vol_TIV; 
    res.vol_fitorg     = mean( abs(res.vol_rel_CGWorg - res.vol_rel_CGW ) );
    res.vol_fitcor     = mean( abs(res.vol_rel_CGWcor - res.vol_rel_CGW ) );
    res.totalorg       = mean([ res.CJVorg res.sdCGorg res.vol_fitorg ]);
    res.totalcor       = mean([ res.CJVcor res.sdCGcor res.vol_fitcor ]);
    clear Yp0; 


    %% write output
    stime2 = cat_io_cmd('      Write output','g5','',job.verb-1,stime2);
    resdir = fullfile(spm_fileparts(job.files{si}),job.resdir); 
    if ~exist(resdir,'dir'), mkdir(resdir); end
    Pcm = spm_file(job.files{si},'path',resdir,'prefix',job.prefix,'ext','.nii');
    Vm2 = Vm; Vm2.fname = Pcm;
    Vm2.descrip = sprintf('%s >> CAT-MP2Rage %s', Vm.descrip,parastr); %####################
    if ~res.isMP2R && res.isT1
      spm_write_vol(Vm2,max(-.5,min(2,Ym2)));
    else
      spm_write_vol(Vm2,Ym2);
    end

% write bias-field / correction map?
    % 

    %% run SPM postbiascorrection 
    if 0 %~res.isMP2R  &&  job.biascorrection > 1  &&  ~job.quicktest
      stime2 = cat_io_cmd('      Post correction','g5','',job.verb-1,stime2);
      matlabbatch{1}.spm.spatial.preproc.channel.vols     = {Pcm};
      for ti = 1:6
        matlabbatch{1}.spm.spatial.preproc.tissue(ti).native = [0 0];
      end
      evalc('spm_jobman(''run'',matlabbatch)');
      movefile( spm_file(Pcm,'prefix','m'), spm_file(Pcm,'prefix','') );
    end

    % save XML/mat
    if job.report
      jobxml = job; jobxml.files = jobxml.files(si); 
      jobxml.result = {Pcm};
      jobxml.res    = res; 
      if exist('scstr','var')
        jobxml.scstr  = scstr; 
      end
      cat_io_xml(spm_file(Pcm,'prefix','catmp2r_','ext','.xml'),jobxml);
    end

    out.files{si} = Pcm;
    if nargout>1
      outs(si)    = res; %#ok<AGROW> 
    end



    % report ? >> subfunction 
    % - input parameter 
    % - estimated parameters
    % - original image + segmentation + new image + change image ???
    % - histogram tissues?
    if job.report 
      res.time = etime(clock,stime);
      report_figure(Vo,Vm2, Yo,Ym2,Yseg, job,spm8,res);
    end

    if res.totalorg > res.totalcor, ccol = [0 0.5 0]; else, ccol = [0.5 0 0]; end
    imgstr  = {'MP1R','MP2R'};
    imgseq  = {'PDw','T1w','T2w'};
    imgfat  = {'fs0','fs1','fsu'};
    if job.verb > 1
      fprintf('%5.0fs\n        Average change rating:  ',etime(clock,stime2));
      cat_io_cprintf( ccol , sprintf('%5.2g - %0.2f > %0.2f\n', res.gradbias_tot*100, res.totalorg, res.totalcor) );
      cat_io_cmd(' ','g5','',job.verb-1);
      fprintf('%5.0fs\n',etime(clock,stime));
    elseif job.verb == 1
      fprintf('%s-%s-%s: ',imgstr{res.isMP2R+1},imgseq{res.isT1+1},imgfat{res.fatsup+1}); 
      cat_io_cprintf( ccol , sprintf('%5.2g - %0.2f > %0.2f', res.gradbias_tot*100, res.totalorg, res.totalcor) ); 
      fprintf('%5.0fs\n',etime(clock,stime));
    end

    if job.spm_cleanupfiles
      [pp,ff,ee] = spm_fileparts(job.files{si}); mfile = fullfile(pp,[ff,ee]); 
      for ci = 1:5
        delete( spm_file(mfile, 'prefix', sprintf('c%d',ci)) ); 
      end
      if job.biascorrection
        delete( spm_file(mfile, 'prefix', 'm' ) ); 
      end
      delete( spm_file(mfile, 'suffix', '_seg8' , 'ext', '.mat' ) ); 
    end

    spm_progress_bar('Set',si); 
  end
end
function report_figure(V,V2, Ym,Ym2,Yc, opt,spm8,res)


  % in case of updates
  spm_orthviews('Reset')
  %clear -global st; global st %#ok<GVMIS> 
  
  % fontsettings
  fontsize  = 10; fontcolor = [0 0 0]; fontname  = 'monospace';
 
  % setup SPM figure
  fg = spm_figure('FindWin','Graphics'); 
  set(0,'CurrentFigure',fg)
  spm_figure('Clear',fg);
  if isempty(fg)
    fg = spm_figure('Create','Graphics','visible','on'); 
  end
  colormap gray; 
  

  % main text report box with header (filename)
  ax = axes('Position',[0.01 0.75 0.98 0.245],'Visible','off','Parent',fg);
  text(0,0.99, ['CAT MP2Rage Precessing: ' ...
    strrep( strrep( spm_str_manip(V.fname,'k80d'),'\','\\'), '_','\_') '       '],...
    'FontSize',fontsize+1,'FontWeight','Bold','Interpreter','tex','Parent',ax);


  % write parameters
  str{1}      = [];
  imgseq      = {'T2w','T1w','PDw'};
  imgstr      = {'MPR','MP2R'};
  imgfat      = {'','fatsup','unkown'};
  % fatshift 
  biasstr     = {'no','yes-basic','yes-extend'};
  ssstr       = {'no','SPM','optimized','background-removal'}; % 0-3
  % logstr      = {'no','exp','log2','log','log10','auto'};
  csfnoise    = {'no','yes'};
  % parameters
  str{1}(end+1).name  = 'Image type:'; 
  str{1}(end).value   = sprintf('%s %s ', ...
    imgseq{res.isT1+1}, imgstr{res.isMP2R+1}, imgfat{res.fatsup+1});
  if cat_get_defaults('extopts.expertgui')
    str{1}(end+1).name  = 'Bias-correction / skull-stripping:'; 
    str{1}(end).value   = sprintf('%s / %s ', ...
      biasstr{opt.biascorrection+1}, ssstr{opt.skullstripping+1});
  else
    str{1}(end+1).name  = 'Bias-correction / int-harmonization /skull-stripping:'; 
    str{1}(end).value   = sprintf('%s / %s / %s ', ...
      biasstr{opt.biascorrection+1}, csfnoise{opt.intscale+1}, ssstr{opt.skullstripping+1});
  end
  if cat_get_defaults('extopts.expertgui')
    str{1}(end+1).name  = 'log-norm / contrast-norm / CSF noise:'; 
    str{1}(end).value   = sprintf('%d / %0.2f / %s ', ...
      opt.logscale, opt.intnorm, csfnoise{opt.restoreLCSFnoise+1}); 
  end
  % SPM intensity thresholds ?
  % SPM volumes
  str{1}(end+1).name  = 'Absolute brain volumes (CSF + GM + WM = TIV, in mm): ';
  str{1}(end).value   = sprintf('%0.0f + %0.0f + %0.0f = %0.0f', res.vol_abs_CGW,res.vol_TIV); 
  str{1}(end+1).name  = 'Relative brain volumes (CSF + GM + WM, in %): ';
  str{1}(end).value   = sprintf('%0.3f + %0.3f + %0.3f', res.vol_rel_CGW); 
  % paras
  str{1}(end+1).name  = 'Intensity-volume fit / average-fit: ';
  if res.vol_fitorg > res.vol_fitcor, ccor = '\color[rgb]{0  0.5   0}'; else, ccor = '\color[rgb]{0.5  0   0}'; end
  str{1}(end).value   = sprintf('%s%0.3f > %0.3f /', ccor, res.vol_fitorg,res.vol_fitcor ); 
  if res.totalorg > res.totalcor, ccor = '\color[rgb]{0  0.5   0}'; else, ccor = '\color[rgb]{0.5  0   0}'; end
  str{1}(end).value   = [str{1}(end).value sprintf('%s%0.3f > %0.3f ', ccor, res.totalorg,res.totalcor )]; 
  %
  str{1}(end+1).name  = 'Grad-change / CJV(org) > CJV(cor) / CJVCGW(org) > CJVCGW(cor): ';
  if res.gradbias_tot < 0,          ccor0 = '\color[rgb]{0  0.5   0}'; else, ccor0 = '\color[rgb]{0.5  0   0}'; end
  if res.CJVorg    > res.CJVcor,    ccor1 = '\color[rgb]{0  0.5   0}'; else, ccor1 = '\color[rgb]{0.5  0   0}'; end
  if res.CJVCGWorg > res.CJVCGWcor, ccor2 = '\color[rgb]{0  0.5   0}'; else, ccor2 = '\color[rgb]{0.5  0   0}'; end
  str{1}(end).value   = sprintf('%s%8g / %s%0.3f > %0.3f / %s%0.3f > %0.3f', ...
    ccor0, res.gradbias_tot, ccor1, res.CJVorg,res.CJVcor, ccor2, res.CJVCGWorg,res.CJVCGWcor); 
  str{1}(end+1).name  = 'Procesing time: ';
  str{1}(end).value   = sprintf('%0.0fs',res.time); 
% +++ add simple/complex color rating (simple=just better, compex=marks)  
  

  htext = zeros(5,2,2);
  for i=1:size(str{1},2)   % main parameter
    htext(1,i,1) = text(0.01,0.98-(0.055*i), str{1}(i).name  ,'FontName',fontname,'FontSize',fontsize,'color',fontcolor,'Interpreter','none','Parent',ax);
    htext(1,i,2) = text(0.51,0.98-(0.055*i), str{1}(i).value ,'FontName',fontname,'FontSize',fontsize,'color',fontcolor,'Interpreter','tex','Parent',ax);
  end


  
  %% == images ==
  spm_orthviews('Reset')
  pos = {[0.008 0.375 0.486 0.35]; [0.506 0.375 0.486 0.35]; [0.008 0.01 0.486 0.35];};
  % T1 + SPM segmentation
  %colormap( [[0 0.02 0.07]; repmat([0.05 0.15 .35],round(59/(crange+2)),1); repmat([ 0 .3 .6],round(59/(crange+2)),1); jet(59 - 2*round(59/(crange+2)))]);
  if res.isT1
    V0 = V; V0.dat(:,:,:) = single(Ym/spm8.mn(2)); V0.dt(1) = 16;
  else
    V0 = V; V0.dat(:,:,:) = single(Ym/max(spm8.mn(1:3))); V0.dt(1) = 16;
  end  
  V0.pinfo  = repmat([1;0],1,size(Ym,3));
  V0.mat    = spm8.Affine * V0.mat; % Vo.mat;   \seg8t.tpm(1).mat
  hh0       = spm_orthviews('Image',spm_vol(fullfile(spm('dir'),'tpm','TPM.nii,1')),pos{1}); % avoid problems with high resolutions
  spm_orthviews('window',hh0,[0 10000]); % just make it black
  spm_orthviews('BB', [-85 -120 -90; 85 95 105]); % this has to be set in the low-resolution image
  hh0       = spm_orthviews('Image',V0,pos{1}); % add correct image after the other settings! 
  spm_orthviews('Caption',hh0,sprintf('%s','original'));
  spm_orthviews('window',hh0,[0 1.2]); 
% #### overlay skull-stripping  

  % print bone marrow
  if res.isT1
    V1 = V; V1.dat(:,:,:) = Ym2;
  else
    V1 = V; V1.dat(:,:,:) = single(Ym2/max(max(max(Ym2(sum(Yc(:,:,:,1:3),4)>.5))))); V1.dt(1) = 16;
  end
  V1.pinfo  = repmat([1;0],1,size(Ym2,3));
  V1.mat    = spm8.Affine * V1.mat; 
  V1.dt(1)  = 16;
  hh1       = spm_orthviews('Image',V1,pos{2}); 
  spm_orthviews('Interp',1);
  spm_orthviews('window',hh1,[0 1.1]);
  spm_orthviews('Caption',hh1,'Optimized');

  % print bone marrow
  V1 = V; V1.dat(:,:,:) = Yc(:,:,:,1)*2 + Yc(:,:,:,2)*3 + Yc(:,:,:,3) + Yc(:,:,:,4)*0.5 + Yc(:,:,:,5)*3.3;
  V1.pinfo  = repmat([1;0],1,size(Ym2,3));
  V1.mat    = spm8.Affine * V1.mat; 
  V1.dt(1)  = 16;
  hh1       = spm_orthviews('Image',V1,pos{3}); 
  spm_orthviews('window',hh1,[0 3.3]);
  spm_orthviews('Caption',hh1,'SPM Label Map');
  spm_orthviews('Reposition',[-25 0 0]);
  spm_orthviews('redraw'); 

% ### thresholded optimized ? 


  % save image
  try % does not work in headless mode without java
    figfs10 = [ findobj(fg,'FontSize',fontsize+1); findobj(fg,'FontSize',fontsize); findobj(fg,'FontSize',fontsize-1);  ...
      findobj(fg,'FontSize',fontsize*0.85); findobj(fg,'FontSize',fontsize*.6);  findobj(fg,'FontSize',fontsize*.4); ]; 
    for fsi = 1:numel(figfs10), try figfs10(fsi).FontSize = figfs10(fsi).FontSize * .75; end; end %#ok<TRYNC> 
    saveas(fg,spm_file(V2.fname,'prefix','report','ext','.png'));
    for fsi = 1:numel(figfs10), try figfs10(fsi).FontSize = figfs10(fsi).FontSize / .75; end; end %#ok<TRYNC> 
  catch
    cat_io_cprintf('err','Error while saving report figure using "saveas".\n');
  end

end

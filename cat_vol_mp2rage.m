function [out,outs] = cat_vol_mp2rage(job)
%cat_vol_mp2rage. Optimization of MP2Rage weighted scans for CAT12. 
% The function uses the unified segmenation of SPM to further optimize a 
% set of input images (MP2RAGE but also other T1-weighteds scans).
%
%  out = cat_vol_mp2rage(job)
%
%  job
%   .files            .. list of MP2Rage images
%   .headtrimming     .. trimming to brain or head (*0-none*,1-brain,2-head)
%   .biascorrection   .. biascorrection (0-no,1-*yes*) 
%   .skullstripping   .. skull-stripping (0-no, 1-SPM, 2-*optimized*)
%   .logscale         .. use log/exp scaling for more equally distributed
%                        tissues (0-none, 1-log, -1-exp, inf-*auto*);
%   .intnorm          .. contrast normalization using the tan of GM normed
%                        values with values between 1.0 - 2.0 for light to 
%                        strong adaptiong (0-none, inf-*auot*)
%   .restoreLCSFnoise .. restore CSF (noise) values below zero (0-no,1-yes)   
%   .report           .. create report file 
%   .prefix           .. filename prefix (strong with PARA for parameter
%                        depending naming, e.g. ... )
%   .verb             .. be verbose (0-no,1-yes,2-details)
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
  %def.bloodvesselscorrection = 0;   % not implemented >> should be done by CAT 
  def.headtrimming      = 0;        % trimming to brain or head (*0-none*,1-brain,2-head)
  def.biascorrection    = 1;        % biascorrection (0-no,1-light(SPM60mm),2-average(SPM60mm+X,3-strong(SPM30+X)) #######
  def.skullstripping    = 2;        % skull-stripping (0-no, 1-SPM, 2-*optimized*)
  def.logscale          = inf;      % use log/exp scaling for more equally distributed
                                    % tissues (0-none, 1-log, -1-exp, inf-*auto*);
  def.intnorm           = -.5;      % contrast normalization using the tan of GM normed
                                    % values with values between 1.0 - 2.0 for light to 
                                    % strong adaptiong (0-none, 1..2-manuel, -0..-2-*auto*)
  def.restoreLCSFnoise  = 1;        % restore values below zero (lower CSF noise)    
  def.prefix            = 'PARA_';  % filename prefix (strong with PARA for parameter
                                    % depending naming, e.g. ... ) 
  def.spm_preprocessing = 1;        % do SPM preprocessing (0-no, 1-yes (if required), 2-always)
  def.report            = 1;        % create a report
  def.verb              = 1;        % be verbose (0-no,1-yes,2-details)
  
  job = cat_io_checkinopt(job,def);


  [job.resdir,job.prefix] = spm_fileparts(job.prefix);
  

  % update prefix
  if isinf(job.logscale), lg = 'A'; else, lg = sprintf('%d',job.logscale); end
  if isinf(job.intnorm),  in = 'A'; else, in = sprintf('%d',job.intnorm);  end
  if cat_io_contains(job.prefix,'PARA')
    job.prefix = strrep(job.prefix,'PARA',sprintf('MP2R_hd%d_bc%d_lg%s_in%s_rn%d_ss%d',...
       job.headtrimming, job.biascorrection, lg, in, job.restoreLCSFnoise, job.skullstripping));
  end
  parastr = sprintf('(hd%d;bc%d;lg%s;in%s;rn%d;ss%d)',...
    job.headtrimming, job.biascorrection, lg, in, job.restoreLCSFnoise, job.skullstripping);

  % prepare output
  out.files = job.files; 


  % TODO: 
  % * optimization criteria
  %    > histogram equation > harmonization ?
  %    > CJV, VOL ... use it for model selection / test
  %    > automatic recursive application such as exp( exp( x ) )
  % * delete temporary files
  % * use side independent tan-scaling? 
  %    - but how?
  % * what to do in case of worse results? - just use the original?
  % * add T2/PD cases !
  % - integrate subdirectory definition in the prefix definition
  % - bias correction warning if brain tissues show high variation or  
  %   corrections of the segmentation (closing) where done? 
  % - add full tissue bias correction case
  % - variable/memory cleanup
  % - faster BET
  % - report (make changes transparent)
  %    > colors to support faster review
  %    > add histogram
  %    > colormap for images?
  %    > intensity peaks org / cor
  %    > resolution, 
  % - cmd line report (add colored result)
  % . additional information about the contrast settings (auto details) and time   
  % . blood vessel correction ? 
  %    - No, should be done by PP, but may aboid precorrection. 
  %    - Maybe indirectly as part of the skull-stripping (reinclude as brain or not) 


  %% main processing 
  spm_progress_bar('Init',numel(job.files),'MP2Rage Optimization','Volumes Complete');
  for si = 1:numel(job.files)
    
    % be verbose
    stime = clock; 
    if job.verb
      fprintf('%4d) %80s: ',si,spm_str_manip( job.files{si}, 'a80'));
    end


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
       
    if job.spm_preprocessing==2 || ... allways run SPM 
      (job.spm_preprocessing==1 && any(~SPMdatae)) ... only if data is missing

      %% SPM segmentation
      if job.verb>1, fprintf('\n'); end
      stime2 = cat_io_cmd('      SPM segmentation','g5','',job.verb-1);

      matlabbatch{1}.spm.spatial.preproc.channel.vols       = job.files(si);
      matlabbatch{1}.spm.spatial.preproc.channel.biasreg    = 0.001;
      matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm   = 60;
      matlabbatch{1}.spm.spatial.preproc.channel.write      = [0 1]; % write bias corrected image
      ngaus = [1 1 2 4 3 2]; 
      for ti = 1:6
        matlabbatch{1}.spm.spatial.preproc.tissue(ti).tpm    = ...
          {fullfile(spm('dir'),'tpm',sprintf('TPM.nii,%d',ti))};
        matlabbatch{1}.spm.spatial.preproc.tissue(ti).ngaus  = ngaus(ti);
        matlabbatch{1}.spm.spatial.preproc.tissue(ti).native = [ti<6 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(ti).warped = [0 0];
      end
      matlabbatch{1}.spm.spatial.preproc.warp.mrf           = 1;
      matlabbatch{1}.spm.spatial.preproc.warp.cleanup       = 1;
      matlabbatch{1}.spm.spatial.preproc.warp.reg           = [0 0.001 0.5 0.05 0.2];
      matlabbatch{1}.spm.spatial.preproc.warp.affreg        = 'mni';
      matlabbatch{1}.spm.spatial.preproc.warp.fwhm          = 0;
      matlabbatch{1}.spm.spatial.preproc.warp.samp          = 3;
      matlabbatch{1}.spm.spatial.preproc.warp.write         = [0 0];
      matlabbatch{1}.spm.spatial.preproc.warp.vox           = NaN;
      matlabbatch{1}.spm.spatial.preproc.warp.bb            = [NaN NaN NaN; NaN NaN NaN];


      % run SPM 
      evalc('spm_jobman(''run'',matlabbatch)');
    else
      stime2 = clock; 
      % if SPM data is incomplete but we are not allow to do it
      if any(~SPMdatae)
        cat_io_cprintf('err','Cannot find all required preprocessed SPM files (*seg8.mat, m*.nii, c#*.nii). Continue with next subject. \n');
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
    
    % load bias corrected 
% ############# use the original if no bias corrected was written? ##########    
    Pm = spm_file(job.files{si},'prefix','m');
    Vm = spm_vol(Pm); 
    Ym = single( spm_read_vols(Vm));
    vx_vol = sqrt(sum(Vm.mat(1:3,1:3).^2)); 
    
    % load segmentation (for brain masking)
    Yseg = zeros(Vm.dim,'single');
    Pseg = cell(1,6); Tth = zeros(1,6); Sth = Tth; 
    for ci = 1:5
      Pseg{ci}       = spm_file(job.files{si},'prefix',sprintf('c%d',ci)); 
      Vseg(ci)       = spm_vol(Pseg{ci}); %#ok<AGROW> 
      Yseg(:,:,:,ci) = single( spm_read_vols(Vseg(ci)) );
      Ymsk           = cat_vol_morph( Yseg(:,:,:,ci)>.9 ,'e') & Ym~=0; 
      Tth(ci)        = cat_stat_nanmedian(Ym(Ymsk(:))); 
      Sth(ci)        = cat_stat_nanstd(Ym(Ymsk(:))); 
    end
    clear Ymsk Vseg;
    Yseg(:,:,:,6) = 1 - sum( Yseg, 4);  
    Tth(6) = cat_stat_nanmedian(Ym(Yseg(:,:,:,6) & Ym~=0));
    isMP2R = Tth(6) > Tth(3)*4;                  % MP2Rage
    isT1   = Tth(3) < Tth(1) & Tth(1) < Tth(2);  % T1 defintion  



    if ~isT1
      cat_io_cprintf('err','No T1w contrast or bad segmenation!\n')
    end
    


    if job.headtrimming
      %% limit image to brain  
      Vm = spm_vol(Pm);   Ym = single( spm_read_vols(Vm)); Vm.mat
      [Ym,BB] = cat_vol_resize(Ym,'reduceBrain',Vm,10,sum(Yseg(:,:,:,1:3),4)>0.5); %Vm = BB.V; 
      Vm.mat
  %    Yseg0   = zeros([size(Ym),size(Yseg,4)],'single');    
  %    for ci = 1:size(Yseg,4)
  %      Yseg0(:,:,:,ci) = cat_vol_resize(Yseg(:,:,:,ci) ,'reduceBrain',Vm,10,sum(Yseg(:,:,:,1:3),4)>0.5);
  %    end
      % Update affine registration for trimming
      %{
      imat = spm_imatrix(spm8.Affine);
      imat(1:3) = imat(1:3) - BB.trans*2; 
      spm8.Affine = spm_matrix(imat);
      Yseg = Yseg0; clear Yseg0; 
      %}
      imat = spm_imatrix(Vm.mat);

      Vm.dim = size(Ym);
      
      stime2 = cat_io_cmd('      Write output','g5','',job.verb-1,stime2);
      Pcm = spm_file(job.files{si},'prefix',job.prefix,'ext','.nii');
      Vm2 = Vm; Vm2.fname = Pcm;
      Vm2.descrip = sprintf('%s >> CAT-MP2Rage %s', Vm.descrip,parastr); %####################
      spm_write_vol(Vm2,Ym);
  
      
    end
 

%##########################
% detect critical segmentations / skull-strippings
% >> stop case processing and continue with next subject
% >> don't repair! but maybe suggest to check and to what aparameter could
%    be changed
% Case 1 - Buchert 125849) failed GM/WM segmentation but BET is ok
% Case 2 - Buchert 150601) extrem low contrast (remove/ignore/bullshit in=bullshit out)  
% Case 3 - Buchert 155553) low contrast (high-res, low-freq-noise) but good segmentation
% Case 4 - Buchert 155802) failed bias correction ...
Yp0  = Yseg(:,:,:,1)*2 + Yseg(:,:,:,2)*3 + Yseg(:,:,:,3);
iCon = @(x) (exp(x) - 1) / (exp(1)-1);

% ######################## WMHs ?! #############  
% ### errors with BVs and MAs ... test for local high variance in tissues
% #############
if 1
  %%
   norm2  = @(Yx,   Ycmm,Ywmm) (Yx - median(Yx(Ycmm(:)))) ./ (median(Yx(Ywmm(:))) - median(Yx(Ycmm(:)))) * 2/3 + 1/3;
   exp1   = exp(1) - 1; 
   Ywmm   = cat_vol_morph( Yseg(:,:,:,2)>.9 ,'e'); 
   Ycmm   = smooth3( (Yseg(:,:,:,3)>.95 | ( (Yseg(:,:,:,6)>.95) .* (Tth(6) < Tth(3))) ) & Ym/spm8.mn(2)<.5)>.5; 
   Tthm   = norm2(Tth,3,2); 
   Ymm    = norm2( Ym, Ycmm, Ywmm); log1str = 'Y'; 
   for i = 1:10
     TX  = [ norm2(  exp(  norm2(Tthm, 3, 2)*2/3+1/3 - 1) , 3, 2);
             norm2(        norm2(Tthm, 3, 2)*2/3+1/3 , 3, 2); 
             norm2( log2(  norm2(Tthm, 3, 2)*2/3+1/3 + 1) , 3, 2); 
             norm2( log(  (norm2(Tthm, 3, 2)*2/3+1/3) * exp1 + 1), 3, 2);
             norm2( log10((norm2(Tthm, 3, 2)*2/3+1/3) * 9 + 1), 3, 2) ];
     [mx,cci] = min( abs(TX(1:5) - 0.7)); 
     Tthm = TX(cci,:);

     if     cci == 1, Ymm2 = exp  (  norm2(Ymm, Ycmm, Ywmm)*2/3+1/3 - 1) / exp1;    lstr = 'exp';
     elseif cci == 2, Ymm2 =         norm2(Ymm, Ycmm, Ywmm)*2/3+1/3;                lstr = 'non';
     elseif cci == 3, Ymm2 = log2 (  norm2(Ymm, Ycmm, Ywmm)*2/3+1/3 + 1);           lstr = 'log2';
     elseif cci == 4, Ymm2 = log  ( (norm2(Ymm, Ycmm, Ywmm)*2/3+1/3) * exp1 + 1);   lstr = 'log';
     elseif cci == 5, Ymm2 = log10( (norm2(Ymm, Ycmm, Ywmm)*2/3+1/3) * 9 + 1);      lstr = 'log10';
     end
     if cci == 2 || mx < 0.01; break; end
     log1str = sprintf('%s(%s)',lstr,log1str);
     Ymm = norm2( real(Ymm2) , Ycmm , Ywmm);
   end
   Ym2 = Ymm * spm8.mn(2);
%Ymmt    = norm2(Ymm2,Ycmm,Ywmm); 
 %             cc(cci) = median( Ymmt(Ygm(:)) 
end

    if job.biascorrection
    %% bias correction
    %  the core 

      stime2 = cat_io_cmd('      Bias correction','g5','',job.verb-1,stime2);

      % define bias correction field for each tissue  
      Ysc = cat_vol_morph( cat_vol_morph( Yseg(:,:,:,1)>.8 & Ym>spm8.mn(1),'do', 3 ), 'dd', 6); 
      Ywi = cat_vol_localstat(Ym2,cat_vol_morph(Yseg(:,:,:,2)>.5,'lc') ,1,2+isT1);
      % remove outlier voxels      
      Ywiv  = cat_vol_localstat(Ywi,Ywi~=0,1,4);
      wivth = cat_stat_nanmedian(Ywiv(Ywiv(:)~=0));
      Ywi(Ywiv>(wivth*2)) = 0; Ywi(smooth3(Ywi~=0)<.5) = 0;  
      clear Ywiv wivth; 
%%
      if job.biascorrection>0 
      %% Use also the GM area to get the closes WM value by a maximum operation (i.e. T1 only). 
      %  We are not using the CSF as we are not knowing much about it and and cannot trust the segmentation too much.  
        [Ym2r,resV] = cat_vol_resize(Ym2,'reduceV',vx_vol,2,4,'mean'); % mean here to handle noise 
        [Ygi1,Ygi2] = cat_vol_resize({ ...
          single(smooth3(sum(Yseg(:,:,:,1:2),4))>.95 & ~Ysc & Ym2>Tth(1)), ...
          single(smooth3(sum(Yseg(:,:,:,1:2),4))>.5  & ~Ysc)},'reduceV',vx_vol,2,4,'meanm'); 
        Ygi1  = cat_vol_localstat(Ym2r,Ygi1>.8,2,2 + isT1); % closer to WM
        Ygi2  = cat_vol_localstat(Ym2r,Ygi2>.8,4,2 + isT1); % more distant to WM
        Ygi   = cat_vol_approx( max(Ygi1,Ygi2) ); clear Ym2r Ygi1 Ygi2
        Ygi   = cat_vol_resize(Ygi,'dereduceV',resV) .* (Yseg(:,:,:,1)>.2);
        % adaption for high noise bias in low intensity regions
        Ygi   = Ygi*.5 + 0.5*(iCon(Ygi / spm8.mn(2))*spm8.mn(2)); 
        Ywi   = max(Ywi,Ygi); 
        clear Ygi
      end
      

      %% estimate bias in the background in MP2Rage or the head MPRage 
      if isMP2R
        %% MP3RAGE - use simple background 
        Ybi = cat_vol_morph( Yseg(:,:,:,6)>.5, 'e',3) & ...
          Ym~=0 & Ym<spm8.mn(2) & Ym>mean(spm8.mn(1:2:3)); Ybi2=Ybi;
        [Ybi,Ymsk,resV] = cat_vol_resize({Ym2 .* Ybi,Ybi},'reduceV',1,4,4,'meanm'); Ybi(Ymsk<.99)=0;  
        Ybi = cat_vol_approx(Ybi,'nh',1,2);
        Ybi = cat_vol_resize(Ybi,'dereduceV',resV) .* Ybi2;
      elseif 0
        %% create smooth head based bias correction map
        %Yp0 = Yseg(:,:,:,1) + Yseg(:,:,:,2)*3 + Yseg(:,:,:,3)*2;
        Yg  = cat_vol_grad(Ym2)./Ym2; 
        Ybi = (cat_vol_morph(Yseg(:,:,:,5)>.5,'l') | ...
               cat_vol_morph(Yseg(:,:,:,2)>.5,'lc')) & (Ym./spm8.mn(1))>.5 & Yg<0.5;% Ybi2=Ybi;
        [Ybi,Ymsk,Yp0r,resV] = cat_vol_resize({Ym2 .* Ybi,Ybi,Yp0},'reduceV',1,4,4,'meanm'); Ybi(Ymsk<.99)=0;
        %clear Yp0; 
        Ybi = cat_vol_morph(Ybi,'gc',9); % higher size to avoid overcorrection 
        Ybi = Ybi ./ median(Ybi(Ybi(:)>0));
        Ybi = cat_vol_inpaint(Ybi,10,10,1,1);
        Ybi = cat_vol_resize(Ybi,'dereduceV',resV); 
      else
        Ybi = zeros(size(Ym)); 
      end
      

      %% mix tissues and approximate bias field
      if isMP2R
        %% simpler faster correction 
        Yi = Ywi / max(eps,cat_stat_nanmedian(Ywi(Ywi(:)>0))) + ...
             Ybi / max(eps,cat_stat_nanmedian(Ybi(Ybi(:)>0))) + ...
             mean(cat(4,Ywi>0,Ywi/max(eps,cat_stat_nanmedian(Ywi(Ywi(:)>0)))),4); % use lower GM weighting
        Yw = cat_vol_approx(Yi,'nn',1,4); 
      else 
        %% for the inpaint method the use of GM information was less optimal 
        % - first estimation for background 
        Yi = max(Ywi,Ybi .* ( cat_vol_morph(Yseg(:,:,:,5)>.5,'l') & (Ym./spm8.mn(3))>.5) * spm8.mn(2) ); 
        Yw = cat_vol_inpaint(Yi,10,40,2,1); Yw = Yw*0.9 + 0.1 * iCon(Yw./spm8.mn(2)).*spm8.mn(2);
        Yw = Yw ./ cat_stat_nanmedian(Yw(Yseg(:,:,:,2)>.95)) .* cat_stat_nanmedian(Ym2(Yseg(:,:,:,2)>.95));
        % - second estimation for 
        Yi = max(Ywi,Yw .* ( sum(Yseg(:,:,:,6),4)>.5 )); 
        Yw = cat_vol_inpaint(Yi,10,6 / job.biascorrection,2,1);  Yw = Yw*0.9 + 0.1 * iCon(Yw./spm8.mn(2)).*spm8.mn(2);
        Yw = Yw ./ cat_stat_nanmedian(Yw(Yseg(:,:,:,2)>.95)); % .* cat_stat_nanmedian(Ym2(Yseg(:,:,:,2)>.95));
        
      end
      
      % apply bias correction
      Ymm = Ym2 ./ Yw;
    else
      % no bias correction
      Ymm = Ym2; 
    end  
  




    %% general exp/log scaling 
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
    elseif 0 %isinf(job.logscale) % auto
      % we typically expect the GM value at about 60% of the WM peak
      opt    = 0.7; % BG
      
      % some images/classes that we need later
      Ygm = Yseg(:,:,:,1)>.1 & ~cat_vol_morph( Yseg(:,:,:,1)>.5 ,'o',4); 
      Ywm = cat_vol_morph( Yseg(:,:,:,2)>.9 ,'e'); 

      Ywmm   = Ywm; %Yseg(:,:,:,2)>.95;
      Ycmm   = smooth3( (Yseg(:,:,:,3)>.95 | ( (Yseg(:,:,:,6)>.95) .* (Tth(6) < Tth(3))) ) & Ym/spm8.mn(2)<.5)>.5; 
   %   norm2  = @(Yx,   Ycmm,Ywmm) (Yx - min(Yx(Ycmm(:)))) ./ (median(Yx(Ywmm(:))) - min(Yx(Ycmm(:))));
  %    norm2i = @(Yx,Yy,Ycmm,Ywmm) Yx .* ( median(Yy(Ywmm(:))) - min(Yy(Ycmm(:))) ) + min(Yy(Ycmm(:)));
      norm2  = @(Yx,   Ycmm,Ywmm) (Yx - median(Yx(Ycmm(:)))) ./ (median(Yx(Ywmm(:))) - median(Yx(Ycmm(:)))) * 2/3 + 1/3;
      norm2i = @(Yx,Yy,Ycmm,Ywmm) (Yx - 1/3) .* ( median(Yy(Ywmm(:))) - median(Yy(Ycmm(:))) )*2/3 + median(Yy(Ycmm(:)));

      %%
      %Ymm = norm2(Ymm,Ycmm,Ywmm); 
      res.logscaleres = ''; 
      res.log2str = 'Y';
      for opti = 1:5
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
              Ymm2    = norm2i(Ymm2,Ymm2,Ycmm,Ywmm); 
              Ymmt    = norm2(Ymm2,Ycmm,Ywmm); 
              cc(cci) = median( Ymmt(Ygm(:)) ); clear Ymmt
            end
          end
        end
        %Ymm = norm2i(Ymm2,Ymm,Ycmm,Ywmm); 
        Ymm = norm2(real(Ymm2),Ycmm,Ywmm); 
        
        ccistr = {'exponential','orgiginal','log2','log','log10'};
        scstr  = sprintf('%12s-scaling (exp|org|log2|log|log10|=%s\b)',ccistr{cci},sprintf('%0.2f|',ccm)); 
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
      stime2 = cat_io_cmd(sprintf('      Scale intensities (%s)',ccistr{cci}),'g5','',job.verb-1,stime2);
    end

    
   
    %% scale intensities 
    %  This is the critical part. 
    % log .. optimize the tissue peaks 
    % fx  .. scaling function to reduce the 
    % Ygs .. add some left-sided noise of CSF values
    
    %fx  = @(x,y,z) min(20,max(0,tan(((x-.5) * z + 0.5) * pi - pi/2)/pi/y/z + .2 + 0.5/y));
    fx2 = @(x,con,sqr) min(inf,max(-inf, ( abs(tan( max(-pi/2, min(pi/2, x * pi * con)))).^sqr .* sign(x)) ));

    if job.intnorm ~= 0
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
      Ym2 = norm2(Ymm,Ycmm,Ywmm); 
    end
% evaluate this ... 
% * compare the difference between Yp0 and Ym2 => RMSE
% * compare the histogram of the tissues => RMSE
% * run iterative optimization ...



    %% restore left-side CSF noise (value below zero that were just cutted)
    if job.restoreLCSFnoise
      if isMP2R
        stime2 = cat_io_cmd('      Restore CSF noise','g5','',job.verb-1,stime2);
        
        % define noise
        pnoise = 0.03; % percentage of noise
        Ymsk = Ym~=0 & Yseg(:,:,:,3)>0.8; % Ym !
        Ygs  = Ymsk .* cat_vol_smooth3X( Yseg(:,:,:,3) .^2 ,2) .* ...
          pnoise .* randn(size(Ym2)) .* min(1,1 + (Ym2 - 1/3)*3).^2; 
        % add noise
        Ym2(Ymsk) = Ym2(Ymsk) + Ygs(Ymsk);
        clear Ygs Ymsk;
      else
        stime2 = cat_io_cmd('      No MP2Rage - no CSF noise restoration.','g8','',job.verb-1,stime2);
      end
    end


    if 0
    %% test figures
      fx = figure; clf; fx.Position(3:4) = [600 250];
      subplot(1,2,1); histogram(Ym(Yb(:))/spm8.mn(2),0:.02:1.5); 
      ylim(ylim.*[1 1.2]); title('original (normed for WM)'); hold on
      [mn,sd] = cat_stat_kmeans(Ym(Yb(:))/spm8.mn(2),3);
      plot(repmat(0.02 +(0:1/4:1) ,2,1),repmat( ylim' ,1,5),'b');
      plot(repmat(mn,2,1),repmat( ylim' ,1,numel(mn)),'r');
      xlabel(sprintf('%0.2f %0.2f %0.2f',sd ));

      subplot(1,2,2); histogram(min(1.5,Ym2(Yb(:))),0:.02:1.5);  
      ylim(ylim.*[1 1.2]); title('optimized'); hold on
      [mn,sd] = cat_stat_kmeans(Ym2(Yb(:)),3)
      plot(repmat( 2/6:1/6:1,2,1),repmat( ylim' ,1,5),'b'); 
      plot(repmat(mn,2,1),repmat( ylim' ,1,numel(mn)),'r');
    elseif false
      %% another test figure
      tisn = {'GM','WM','CSF'};
      for ci = 1:3
        Ymsk = Yseg(:,:,:,ci)>.5; 
        [segmn{ci},segsd{ci}] = cat_stat_kmeans(Ym(Ymsk)/spm8.mn(2),3);
        fx = figure(3949); fx.Position(3:4) = [600 400]; 
        subplot(2,3,ci);   histogram( Ym(Ymsk(:)) / spm8.mn(2),0:.02:1.5); title(tisn{ci})
        subplot(2,3,3+ci); histogram(Ym2(Ymsk(:)),0:.02:1.5); title(tisn{ci})

      end    
    end



    
    %% skull-stripping
    if job.skullstripping == 2 || ( job.skullstripping == 3 && isMP2R)
      stime2 = cat_io_cmd('      Optimized skull-stripping','g5','',job.verb-1,stime2);
    
      %% cleanup
      Ygw = smooth3( sum(Yseg(:,:,:,1:2),4) ) > .9;
      [~,Yd] = cat_vol_downcut(single(Ygw),Ym/spm8.mn(3) .* (1-sum( Yseg(:,:,:,5:6), 4)),0.01);
      Ygw2 = Ym/spm8.mn(2)<1.2 & smooth3(Yd<.5)>.5; 
      % edge detection to improve skull-stripping
      Yg  = cat_vol_grad(Ym/spm8.mn(2));
      Yd  = cat_vol_div(Ym/spm8.mn(2));
      % optimize skull-stripping
      Yb  = sum( Yseg(:,:,:,1:3),4); 
      Yb  = max(Yb,cat_vol_morph(Yb > .5 ,'ldc',5));
      Yb  = max(Yb,Ygw2);
      %
      Yb  = smooth3( Yb )>.5; 
      Yb  = Yb | cat_vol_smooth3X( cat_vol_morph(Yb,'dd',2.5) & ...
        Ym<(spm8.mn(2)*0.4+0.6*spm8.mn(1)) & Ym>mean(spm8.mn(3)*0.4+0.6*spm8.mn(1)) & ...
        Yd>-0.05 & Yg<.7 , 2)>.5;
    else
      % simple SPM skull-stripping as CSF+GM+WM 
      stime2 = cat_io_cmd('      SPM skull-stripping','g5','',job.skullstripping && job.verb>1,stime2);
      Yb  = smooth3( sum( Yseg(:,:,:,1:3),4) )>.5; 
    end
  


    %% prepare output
    if job.skullstripping
      % full classical skull-stripping
      Ym2 = Ym2 .* Yb; 
    else
      % background-stripping and skull modification 
      Ym2 = Ym2 .* smooth3( max(1/3 * smooth3(Yseg(:,:,:,4)>.5) , ...
        smooth3(sum( max( Yb ,Yseg(:,:,:,5)),4))>.5)) ; 
    end

    
    %% evaluation 
    Yp0   = Yseg(:,:,:,1)*2 + Yseg(:,:,:,2)*3 + Yseg(:,:,:,3);
    % intensity-based measures
    % - the correction should optimize the coefficient of joint variation (CJV) 
    fcjvx = @(x,y) ( cat_stat_nanstd(x(round(y(:))==1)) +  cat_stat_nanstd(x(round(y(:))==2)) +  cat_stat_nanstd(x(round(y(:))==3)) ) ./ ...
                   (cat_stat_nanmean(x(round(y(:))==1)) + cat_stat_nanmean(x(round(y(:))==2)) + cat_stat_nanmean(x(round(y(:))==3)));
    fcjv  = @(x,y) ( cat_stat_nanstd(x(round(y(:))==2)) +  cat_stat_nanstd(x(round(y(:))==3)) ) ./ ...
                   (cat_stat_nanmean(x(round(y(:))==2)) + cat_stat_nanmean(x(round(y(:))==3)));
    res.CJVorg      = fcjv(Ym,Yp0);
    res.CJVcor      = fcjv(Ym2,Yp0);
    res.CJVCGWorg   = fcjvx(Ym,Yp0);
    res.CJVCGWcor   = fcjvx(Ym2,Yp0);
    % Peak width similarity
    res.sdCGWorg    = [ cat_stat_nanstd(  Ym(round(Yp0(:)*3)==3)) cat_stat_nanstd(  Ym(round(Yp0(:))==2))  cat_stat_nanstd(  Ym(round(Yp0(:))==3))] / spm8.mn(2);
    res.sdCGWcor    = [ cat_stat_nanstd( Ym2(round(Yp0(:)*3)==3)) cat_stat_nanstd( Ym2(round(Yp0(:))==2))  cat_stat_nanstd( Ym2(round(Yp0(:))==3))];
    res.sdCGorg     = std(res.sdCGWorg(2:3));
    res.sdCGcor     = std(res.sdCGWcor(2:3));
    % volumetric measures 
    % - the resuling map should fit better to the existing segmentation 
    res.vol_TIV     = cat_stat_nansum( Yp0(:)>.5)  .* prod(vx_vol) / 1000;
    res.vol_abs_CGW = [ cat_stat_nansum(round(Yp0(:))==1) cat_stat_nansum(round(Yp0(:))==2) cat_stat_nansum(round(Yp0(:))==3) ] .* prod(vx_vol) / 1000;
    res.vol_rel_CGW = res.vol_abs_CGW ./ res.vol_TIV; 
    res.vol_abs_CGWorg = [ 
      cat_stat_nansum(Yb(:) & Ym(:)<spm8.mn(1) ) ...
      cat_stat_nansum(Yb(:) & Ym(:)>spm8.mn(1) & Ym(:)<spm8.mn(2) ) ...
      cat_stat_nansum(Yb(:) & Ym(:)>spm8.mn(2) ) ] .* prod(vx_vol) / 1000;
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
    spm_write_vol(Vm2,Ym2);

    % save XML/mat
    jobxml = job; jobxml.files = jobxml.files(si); 
    jobxml.result = {Pcm};
    jobxml.res    = res; 
    if exist('scstr','var')
      jobxml.scstr  = scstr;
    end
    cat_io_xml(spm_file(Pcm,'prefix','catmp2r_','ext','.xml'),jobxml);

    out.files{si} = Pcm;
    if nargout>1
      outs(si)    = res; %#ok<AGROW> 
    end

    % report ? >> subfunction 
    % - input parameter 
    % - estimated parameters
    % - original image + segmentation + new image + change image ???
    % - histogram tissues?
    report_figure(Vm,Vm2, Ym,Ym2,Yseg, job,spm8,res);

    if res.totalorg > res.totalcor, ccol = [0 0.5 0]; else, ccol = [0.5 0 0]; end
    if job.verb > 1
      fprintf('%5.0fs\n        Average change rating:  ',etime(clock,stime2));
      cat_io_cprintf( ccol , sprintf('%0.2f > %0.2f\n', res.totalorg, res.totalcor) );
      cat_io_cmd(' ','g5','',job.verb-1);
      fprintf('%5.0fs\n',etime(clock,stime));
    elseif job.verb == 1
      cat_io_cprintf( ccol , sprintf('%0.2f > %0.2f', res.totalorg, res.totalcor) ); 
      fprintf('%5.0fs\n',etime(clock,stime));
    end


    spm_progress_bar('Set',si); 
  end
end
function report_figure(V,V2, Ym,Ym2,Yc, opt,spm8,res)


  % in case of updates
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
  % - trimming - bias - skull-stripping
  % - log - contrast - noise
  % - CJV
  % - vols?
  str{1} = [];
  trimmingstr = {'no','brain','head'};
  biasstr     = {'no','yes','yes'};
  ssstr       = {'no','SPM','optimized'};
  % logstr      = {'no','exp','log2','log','log10','auto'};
  csfnoise    = {'no','yes'};
  % parameters
  str{1}(end+1).name  = 'trimming / bias-correction / skull-stripping:'; 
  str{1}(end).value   = sprintf('%s / %s / %s ', ...
    trimmingstr{opt.headtrimming+1}, biasstr{opt.biascorrection+1}, ssstr{opt.skullstripping+1}); 
  str{1}(end+1).name  = 'log-norm / contrast-norm / CSF noise:'; 
  str{1}(end).value   = sprintf('%d / %0.2f / %s ', ...
    opt.logscale, opt.intnorm, csfnoise{opt.restoreLCSFnoise+1}); 
  % SPM intensity thresholds ?
  % SPM volumes
  str{1}(end+1).name  = 'Absolute brain volumes (CSF + GM + WM = TIV): ';
  str{1}(end).value   = sprintf('%0.0f + %0.0f + %0.0f = %0.0f', res.vol_abs_CGW,res.vol_TIV); 
  str{1}(end+1).name  = 'Relative brain volumes (CSF + GM + WM): ';
  str{1}(end).value   = sprintf('%0.3f + %0.3f + %0.3f', res.vol_rel_CGW); 
  % paras
  str{1}(end+1).name  = 'Intensity-volume fit: ';
  if res.vol_fitorg > res.vol_fitcor, ccor = '\color[rgb]{0  0.5   0}'; else, ccor = '\color[rgb]{0.5  0   0}'; end
  str{1}(end).value   = sprintf('%s%0.3f > %0.3f ', ccor, res.vol_fitorg,res.vol_fitcor ); 
  str{1}(end+1).name  = 'CJV(org) > CJV(cor) / CJVCGW(org) > CJVCGW(cor): ';
  if res.CJVorg    > res.CJVcor,    ccor1 = '\color[rgb]{0  0.5   0}'; else, ccor1 = '\color[rgb]{0.5  0   0}'; end
  if res.CJVCGWorg > res.CJVCGWcor, ccor2 = '\color[rgb]{0  0.5   0}'; else, ccor2 = '\color[rgb]{0.5  0   0}'; end
  str{1}(end).value   = sprintf('%s%0.3f > %0.3f / %s%0.3f > %0.3f', ...
    ccor1, res.CJVorg,res.CJVcor, ccor2, res.CJVCGWorg,res.CJVCGWcor); 
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
  V0 = V; V0.dat(:,:,:) = single(Ym/spm8.mn(2)); V0.dt(1) = 16;
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
  V1 = V; V1.dat(:,:,:) = Ym2;
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
  %try % does not work in headless mode without java
    figfs10 = [ findobj(fg,'FontSize',fontsize+1); findobj(fg,'FontSize',fontsize); findobj(fg,'FontSize',fontsize-1);  ...
      findobj(fg,'FontSize',fontsize*0.85); findobj(fg,'FontSize',fontsize*.6);  findobj(fg,'FontSize',fontsize*.4); ]; 
    for fsi = 1:numel(figfs10), try figfs10(fsi).FontSize = figfs10(fsi).FontSize * .75; end; end %#ok<TRYNC> 
    saveas(fg,spm_file(V2.fname,'prefix','report','ext','.png'));
    for fsi = 1:numel(figfs10), try figfs10(fsi).FontSize = figfs10(fsi).FontSize / .75; end; end %#ok<TRYNC> 
 % end

end

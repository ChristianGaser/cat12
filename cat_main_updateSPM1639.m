function [Ysrc,Ycls,Yb,Yb0,job,res,T3th,stime2] = cat_main_updateSPM1639(Ysrc,P,Yy,tpm,job,res,stime,stime2)
% ______________________________________________________________________
%  Update SPM preprocessing. 
%  Subfunction of cat_main.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$


  %dbs = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

  global cat_err_res; 

  clsint = @(x) round( sum(res.mn(res.lkp==x) .* res.mg(res.lkp==x)') * 10^5)/10^5;
  
  [pth,nam] = spm_fileparts(res.image0(1).fname); %#ok<ASGLU> % original   
  
  % RD202211: Added SPM based detection of high intensity backgronds of MP2Rage scans
  res.isMP2RAGE = min( res.mn(res.lkp==3 & res.mg'>0.2) ) < min( res.mn(res.lkp==max(res.lkp) & res.mg'>0.2) ); 
    
  % voxel size parameter
  vx_vol  = sqrt(sum(res.image(1).mat(1:3,1:3).^2));    % voxel size of the processed image
  vx_vol0 = sqrt(sum(res.image0(1).mat(1:3,1:3).^2));
  vx_volp = prod(vx_vol)/1000;
  
  d = res.image(1).dim(1:3);

  %% RD20221102: Characterize usefulness of SPM classification.
  %  Analyse and store basic information of the SPM segmentation to improve
  %  handling of different protocols such as MP2rage. In general SPM gives 
  %  some clear sharp tissue output but in some cases (e.g., registration/
  %  segmenation errors) the tissues are quite smooth. See also help entry.
  %  I use 20 buckets to get the number of voxels below 25%, 50%, and 75%, 
  %  but also get 90% etc. 
  hbuckets = 20; 
  if isfield(res,'spmP0'), res = rmfield(res,'spmP0'); end
  for ti = 1:size(P,4)
    Pti = P(:,:,:,ti); 
    res.spmP0.hstbuckets    = 256/(hbuckets*2):256/hbuckets:256; 
    res.spmP0.hst(ti,:)     = hist(Pti(Pti(:)>0),256/(hbuckets*2):256/hbuckets:256) / sum(Pti(:)>0); %#ok<HIST> 
    res.spmP0.hstsumi(ti,:) = 1 - cumsum(res.spmP0.hst(ti,:)); 
    res.spmP0.mn(ti)        = cat_stat_nanmean(   single(Pti(Pti(:)>0)) / 255 ); 
    res.spmP0.sd(ti)        = cat_stat_nanstd(    single(Pti(Pti(:)>0)) / 255 ); 
    res.spmP0.md(ti)        = cat_stat_nanmedian( single(Pti(Pti(:)>0)) / 255 ); 
    %res.spmP0.Q1(ti)       = cat_stat_nanmedian( single(Pti( Pti(:)<res.spmP0.md(ti)*255 ) ) / 255); 
    %res.spmP0.Q2(ti)       = cat_stat_nanmedian( single(Pti( Pti(:)>res.spmP0.md(ti)*255 ) ) / 255); 
  end
  res.spmP0.help = { ...
    'This structure contains some parameters to characterize the initial SPM segmentation to detect possible problems. ' 
    sprintf('It contains a normalized histogram (hst) of %d buckets sampled at hstbuckets (center point with borders inbetween!). ',hbuckets) 
    'The inverse cumultative sum of the histogram (hstsumi) defines how many values are above the histogram bucket. ' 
    'A good segment should include a high proposion of high probable voxels. ' 
    'Moreover, we saved the mean(mn), standard deviation (sd), median(md) of the individual tissue desnity maps of SPM. ' 
    'In general, higher mn, md, Q1/2 and lower sd values (<.5) are expected otherwise segmentation is maybe inproper. ' 
     };
  clear Pti;
  % RD20221102: Add warnings 
  % These threshold are arbitrary choosen and updates are maybe required!
  if job.extopts.expertgui>1 && ...
    ( any( res.spmP0.md(1:3) < [.5 .9 .3],2) || any( res.spmP0.mn(1:3) < [.4 .7 .2],2) || ...
      any( any( res.spmP0.hstsumi(1:3,round(hbuckets*[.25,.50,.75])) < [0.6 0.5 0.4; 0.8 0.7 0.6; 0.5 0.2 0.05 ] )) )
      cat_io_addwarning('cat_main_updateSPM:lowSPMseg',...
        sprintf(['Low SPM segmentation quality with larger areas of mixed tissues. \\\\n' ...
                 '  Median: GM=%0.2f, WM=%0.2f, CSF=%0.2f \\\\n' ...
                 '  Mean:   GM=%0.2f, WM=%0.2f, CSF=%0.2f'], ...
                 res.spmP0.md(1:3),res.spmP0.mn(1:3)),2,[1 2]);
  end
  
  %% RD20221102: Intensity overlap in brain tissue classes 
  %  Besides the probabilty we may check also the intensities. If two
  %  classes show a high overlap of intensities the segmentation is mabe 
  %  inaccurate.
  minprob = 32; hbuckets = 30;
  trange = [min(res.mn(res.lkp<6)) max(res.mn(res.lkp<6))]; 
  trange = trange + [ -diff(trange) +diff(trange) ] / 4; 
  trange = trange(1) : diff(trange)/(hbuckets-2) : trange(2);
  res.spmP0.hstibuckets = trange;
  for ti=1:6
    htmp = hist(Ysrc(P(:,:,:,ti)>minprob),trange); %#ok<HIST> 
    res.spmP0.hsti(ti,:) = htmp; %.Values;
  end
  res.spmP0.hstidiffCGW = [ 
    mean( abs(diff(res.spmP0.hsti(1:2:3,:))) ./ max(eps,sum(res.spmP0.hsti(1:2:3,:)))), ...
    mean( abs(diff(res.spmP0.hsti(1:2,:)))   ./ max(eps,sum(res.spmP0.hsti(1:2,:)))), ...
    mean( abs(diff(res.spmP0.hsti(2:3,:)))   ./ max(eps,sum(res.spmP0.hsti(2:3,:))))];
  res.spmP0.hstidiff = mean( min([ 
    abs(diff(res.spmP0.hsti(1:2,:)))   ./ max(eps,sum(res.spmP0.hsti(1:2,:))); 
    abs(diff(res.spmP0.hsti(2:3,:)))   ./ max(eps,sum(res.spmP0.hsti(2:3,:))); 
    abs(diff(res.spmP0.hsti(1:2:3,:))) ./ max(eps,sum(res.spmP0.hsti(1:2:3,:))); 
    ]));
  res.spmP0.help = [res.spmP0.help; {
    sprintf('The hsti is the histogram of the T1 image for each class (prob>%0.2%%,%d buckets) ',minprob/255,hbuckets)};
    'that are combined into the average intensity difference between classes histidiff with 1 for the ideal case without overlap and 0 for high overlap. '];
  if job.extopts.expertgui>1 && res.spmP0.hstidiff<0.5
      cat_io_addwarning('cat_main_updateSPM:highIntOverlapBetweenClasses',...
        sprintf('High overlap of image intensies in different classes\\\\n  CG=%0.2f,GW=%0.2f,CW=%0.2f,CGW=%0.2f',res.spmP0.hstidiffCGW,res.spmP0.hstidiff),2,[1 2]);
  end
  if 0
    %%
    figure, hold on; histogram(Ysrc(P(:,:,:,1)>minprob),trange); histogram(Ysrc(P(:,:,:,2)>minprob),trange); histogram(Ysrc(P(:,:,:,3)>minprob),trange);

  end
  %%
  stime2 = cat_io_cmd('  Update Segmentation','g5','',job.extopts.verb-1,stime2); 
  

  % Create brain mask based on the the TPM classes
  % cleanup with brain mask - required for ngaus [1 1 2 4 3 2] and R1/MP2Rage like data 
  %YbA = zeros(d,'single');
  Vb = tpm.V(1); Vb.pinfo(3) = 0; Vb.dt=16; 
  Vb.dat = single(exp(tpm.dat{1}) + exp(tpm.dat{2}) + exp(tpm.dat{3})); 
  YbA = cat_vol_sample(res.tpm(1),Vb,Yy,1);
  %for z=1:d(3)
   % YbA(:,:,z) = spm_sample_vol(Vb,double(Yy(:,:,z,1)),double(Yy(:,:,z,2)),double(Yy(:,:,z,3)),1); 
  %end
  if (isfield(job,'useprior') && ~isempty(job.useprior) ), bth = 0.5; else, bth = 0.1; end
  if round(max(YbA(:))/Vb.pinfo(1)), YbA=YbA>bth*Vb.pinfo(1); else, YbA=YbA>bth; end
  % add some distance around brainmask (important for bias!)
  YbA = YbA | cat_vol_morph(YbA & sum(P(:,:,:,1:2),4)>4 ,'dd',2.4,vx_vol);
  
  if size(P,4)==3
    % create background class that is required later
    Yb = cat_vol_morph( smooth3( sum(P(:,:,:,1:2),4) ) > 192 , 'ldc' , 5 , vx_vol ); 
    Yb = cat_vol_morph( Yb , 'ldo' ); 
    for i=1:3, P(:,:,:,i) = P(:,:,:,i) .* cat_vol_ctype(Yb); end 
    P(:,:,:,4) = cat_vol_ctype((1 - Yb) * 255); 
  end
  
  % transfer tissue outside the brain mask to head  ... 
  % RD 201807: I am not sure if this is a good idea. Please test this with children! 
  for i=1:3
    P(:,:,:,4) = cat_vol_ctype(single(P(:,:,:,4)) + single(P(:,:,:,i)) .* single(~YbA)); 
    P(:,:,:,i) = cat_vol_ctype(single(P(:,:,:,i)) .* single(YbA)); 
  end
  
  
  % Cleanup for high resolution data
  % Alghough the old cleanup is very slow for high resolution data, the   
  % reduction of image resolution removes spatial segmentation information. 
  if job.opts.redspmres==0 % already done in case of redspmres
    if max(vx_vol)<1.5 && mean(vx_vol)<1.3
      % RD202102: resizing adds maybe too much blurring that can trouble other functions
      %for i=1:size(P,4), [Pc1(:,:,:,i),RR] = cat_vol_resize(P(:,:,:,i),'reduceV',vx_vol,job.extopts.uhrlim,32); end %#ok<AGROW>
      for i=1:size(P,4), [Pc1(:,:,:,i),BB] = cat_vol_resize(P(:,:,:,i),'reduceBrain',vx_vol,4,YbA); end %#ok<AGROW>
      Pc1 = cat_main_clean_gwc1639(Pc1,max(1,min(2,job.extopts.cleanupstr*2)));
      Ybb = ones(size(YbA),'uint8'); Ybb(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) = uint8(1); 
      for i=1:size(P,4), P(:,:,:,i) = Ybb.*P(:,:,:,i) + cat_vol_resize(Pc1(:,:,:,i),'dereduceBrain',BB); end 
      %for i=1:size(P,4), P(:,:,:,i)   = cat_vol_resize(Pc1(:,:,:,i),'dereduceV',RR); end 
      clear Pc1 Pc2 Ybb;
    end
  end
  if ~( (isfield(job,'useprior') && ~isempty(job.useprior) ) && ... 
        (isfield(res,'ppe') && ~res.ppe.affreg.highBG) )
    clear YbA;
  end
  
  % some reports
  for i=1:size(P,4), Pt = P(:,:,:,i); res.ppe.SPMvols0(i) = cat_stat_nansum(single(Pt(:)))/255 .* prod(vx_vol) / 1000; end; clear Pt; 
  
  
  % garantee probability 
  sP = (sum(single(P),4)+eps)/255;
  for k1=1:size(P,4), P(:,:,:,k1) = cat_vol_ctype(single(P(:,:,:,k1))./sP); end
  clear sP;

  
  % Use median for WM threshold estimation to avoid problems in case of WMHs!
  WMth = double(max( clsint(2) , cat_stat_nanmedian(Ysrc(P(:,:,:,2)>192)) )); 
  if clsint(3)>clsint(2) % invers
    CMth = clsint(3); 
  else
    CMth = min( [  clsint(1) - diff([clsint(1),WMth]) , clsint(3) ]);
  end
  T3th = double([ CMth , clsint(1) , WMth]);


  %% Some error handling
  %    ds('l2','',vx_vol,Ysrc./WMth,Yp0>0.3,Ysrc./WMth,Yp0,80)
  Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
  if isfield(res,'Ylesion') && sum(res.Ylesion(:)>0)
    res.Ylesion = cat_vol_ctype( single(res.Ylesion) .* (Yp0>0.2) ); 
    for k=1:size(P,4), Yl = P(:,:,:,k); Yl(res.Ylesion>0.5) = 0; P(:,:,:,k) = Yl; end  
    Yl = P(:,:,:,3); Yl(res.Ylesion>0.5) = 255; P(:,:,:,3) = Yl; clear Yl; 
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
  end
  if sum(Yp0(:)>0.3)<100 
    % this error often depends on a failed affine registration, where SPM
    % have to find the brain in the head or background
    BGth  = min(cat_stat_nanmean(Ysrc( P(:,:,:,end)>128 )),clsint(6));
    HDHth = clsint(5);
    HDLth = clsint(4);
    clsvol = nan(1,size(P,4)); for ci=1:size(P,4), Yct = P(:,:,:,ci)>128; clsvol(ci) = sum(Yct(:))*vx_volp; end; clear Yct; 
    if size(P,4)==6
        error('CAT:cat_main:SPMpreprocessing:emptySegmentation', ...
         sprintf(['Empty Segmentation: \n ' ...
          'Possibly the affine registration failed. Please check image orientation.\n' ...
          ' Tissue class:           %10s%10s%10s%10s%10s%10s\n' ...
          ' Rel. to image volume:   %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n' ...
          ' Rel. to brain volume:   %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n' ...
          ' Tissue intensity:       %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f'],...
          'BG','CSF','GM','WM','HDH','HDL', ...
          [ clsvol([6 3 1 2 4 5])/cat_stat_nansum(clsvol)*100, clsvol([6 3 1 2 4 5])/cat_stat_nansum(clsvol(1:3))*100, BGth,T3th,HDHth,HDLth]));  %#ok<SPERR>
    elseif size(P,4)==4 % skull-stripped
        error('CAT:cat_main:SPMpreprocessing:emptySegmentation', ...
         sprintf(['Empty Segmentation: \n ' ...
          'Possibly the affine registration failed. Please check image orientation.\n' ...
          ' Tissue class:           %10s%10s%10s%10s\n' ...
          ' Rel. to image volume:   %10.2f%10.2f%10.2f%10.2f\n' ...
          ' Rel. to brain volume:   %10.2f%10.2f%10.2f%10.2f\n' ...
          ' Tissue intensity:       %10.2f%10.2f%10.2f%10.2f'],...
          'BG','CSF','GM','WM', ...
          [ clsvol([4 3 1 2])/cat_stat_nansum(clsvol)*100, clsvol([4 3 1 2])/cat_stat_nansum(clsvol(1:3))*100, BGth,T3th]));  %#ok<SPERR>
    else
        error('CAT:cat_main:SPMpreprocessing:emptySegmentation', ['Empty Segmentation: ' ...
           'Possibly the affine registration failed. Please check image orientation.\n']); 
    end
  end

  
  
  %%
  Yp0(smooth3(cat_vol_morph(Yp0>0.3,'lo'))<0.5)=0; % not 1/6 because some ADNI scans have large "CSF" areas in the background 
  Yp0     = Yp0 .* cat_vol_morph(Yp0 & (Ysrc>WMth*0.05),'lc',2);
  Yp0toC  = @(Yp0,c) 1-min(1,abs(Yp0-c));

  
  % values are only used if errors occur
  cat_err_res.init.T3th = T3th; 
  cat_err_res.init.subjectmeasures.vol_abs_CGW = [prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),1)), ... CSF
                                                  prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),2)), ... GM 
                                                  prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),3)), ... WM
                                                  prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),4))];  % WMH
  cat_err_res.init.subjectmeasures.vol_TIV     =  sum(cat_err_res.init.subjectmeasures.vol_abs_CGW); 
  cat_err_res.init.subjectmeasures.vol_rel_CGW =  cat_err_res.init.subjectmeasures.vol_abs_CGW ./ ...
                                                  cat_err_res.init.subjectmeasures.vol_TIV;
  [cat_err_res.init.Yp0,cat_err_res.init.BB] = cat_vol_resize(Yp0,'reduceBrain',vx_vol,2,Yp0>0.5); 
  cat_err_res.init.Yp0 = cat_vol_ctype(cat_err_res.init.Yp0/3*255);
  clear Yp0; 

  % ### This can not be reached because the mask field is removed by SPM! ###
  if isfield(res,'msk') 
    Ybg = ~res.msk.dat; 
    P4  = cat_vol_ctype( single(P(:,:,:,6)) .* (Ysrc<T3th(2))  .* (Ybg<0.5) + single(P(:,:,:,4)) .* (Ybg<0.5) ); % remove air in head
    P5  = cat_vol_ctype( single(P(:,:,:,6)) .* (Ysrc>=T3th(2)) .* (Ybg<0.5) + single(P(:,:,:,5)) .* (Ybg<0.5) ); % remove air in head
    P6  = cat_vol_ctype( single(sum(P(:,:,:,4:5),4)) .* (Ybg>0.5) + single(P(:,:,:,6)) .* (Ybg>0.5) ); % add objects/artifacts to background
    P(:,:,:,4) = P4;
    P(:,:,:,5) = P5;
    P(:,:,:,6) = P6;
    clear P4 P5 P6 Ybg; 
  end
  

  
  
  %% Skull-Stripping
  %  ----------------------------------------------------------------------
  %  Update Skull-Stripping 1
  %  ----------------------------------------------------------------------
  stime2 = cat_io_cmd('  Update Skull-Stripping','g5','',job.extopts.verb-1,stime2); 
  if (isfield(job,'useprior') && ~isempty(job.useprior) && strcmp(job.opts.affreg,'prior') ) && ... 
     (isfield(res,'ppe') && ~res.ppe.affreg.highBG) 
    % RD202010: use longitudinal skull-stripping 
    [pp,ff,ee] = spm_fileparts(char(job.useprior));
    if isfield(job.output.BIDS,'BIDSyes') % I am not sure if separation is needed or if we simply try with/without mri-dir
      Pavgp0 = fullfile(pp,[strrep(ff,'avg_','p0avg_'),ee]);
      if ~exist(Pavgp0,'file')
        Pavgp0 = fullfile(pp,'mri',[strrep(ff,'avg_','p0avg_'),ee]);
      end      
    else
      Pavgp0 = fullfile(pp,'mri',[strrep(ff,'avg_','p0avg_'),ee]);
      if ~exist(Pavgp0,'file')
        Pavgp0 = fullfile(pp,[strrep(ff,'avg_','p0avg_'),ee]);
      end
    end

% RD20220213: 
%  For the development model with longitudinal TPM you may have to add the affine registration. 
%  However it seems that the adaption of the brainmask works quite well ... 
%  but maybe it is better to full deactive the skull-stripping in the 
%  plasticity/aging case

    % get gradient and divergence map (Yg and Ydiv)
    [Ytmp,Ytmp,Yg,Ydiv] = cat_main_updateSPM_gcut0(Ysrc,P,vx_vol,T3th); clear Ytmp;  %#ok<ASGLU>
    if exist(Pavgp0,'file')
      % the p0avg should fit optimal  
      if any(vx_vol0 ~= vx_vol) % if the data was internaly resampled we have to load it via imcalc
        [Vb,Yb] = cat_vol_imcalc(spm_vol(Pavgp0),spm_vol(res.image.fname),'i1',struct('interp',3,'verb',0,'mask',-1)); clear Vb;  %#ok<ASGLU>
      else
        Yb = spm_read_vols(spm_vol(Pavgp0));
      end
      Yb = Yb > 0.5; 
    else
      % otherwise it would be possible to use the individual TPM 
      % however, the TPM is more smoothed and is therefore only second choice  
      cat_io_addwarning('cat_main_updateSPM:miss_p0avg',sprintf('Cannot find p0avg use TPM for brainmask: \n  %s\n',Pavgp0),2,[1 2]);
      Yb = YbA > 0.5;
      clear YbA
    end
    Ybb = cat_vol_ctype(cat_vol_smooth3X(Yb,0.5)*256); 
    
    %% correct tissues
    %  RD20221224: Only the brainmask wasn't enough and we need to cleanup 
    %              the segmentation also here (only for long pipeline)
    % move brain tissue to head tissues or vice versa
    for ti = 1:3
      if ti == 1 % GM with soft boundary to reduce meninges
        Ynbm = cat_vol_ctype( single(P(:,:,:,ti)) .* (1 - max(0,2 * smooth3(Yb) - 1) ) ); 
        Ybm  = cat_vol_ctype( single(P(:,:,:,5))  .* (    max(0,2 * smooth3(Yb) - 1) ) ); 
      elseif ti == 2 % WM with very soft boundary because we exptect no WM close to the skull
        Ynbm = cat_vol_ctype( single(P(:,:,:,ti)) .* (1 - max(0,2 * single(Ybb)/255 - 1) ) ); 
        Ybm  = cat_vol_ctype( single(P(:,:,:,5))  .* (    max(0,2 * single(Ybb)/255 - 1) ) ); 
      else % CSF with hard boundary
        Ynbm = cat_vol_ctype( single(P(:,:,:,ti)) .* (1 - Yb) ); 
        Ybm  = cat_vol_ctype( single(P(:,:,:,5))  .* (    Yb) ); 
      end
      P(:,:,:,ti) = P(:,:,:,ti) - Ynbm + Ybm; 
      P(:,:,:,5)  = P(:,:,:,5)  + Ynbm - Ybm; 
      clear Ynbm Ybm;
    end
    % some extra GM cleanup for meninges
    Yngm = P(:,:,:,1) .* uint8( Ybb<255 & (P(:,:,:,1)>64) & (smooth3( single(P(:,:,:,1)>64) )<0.5) );
    P(:,:,:,1) = P(:,:,:,1) - Yngm; P(:,:,:,5) = P(:,:,:,5) + Yngm; %if ~debug, clear Yngm; end
    % some further hard GM cleanup ? 
    %{
    Yp0avg = spm_read_vols(spm_vol(Pavgp0));
    Yngm = P(:,:,:,1) .* uint8( cat_vol_morph( Yp0avg < 1.75 , 'de' , 3, vx_vol) & Yp0yvg>0 );
    P(:,:,:,1) = P(:,:,:,1) - Yngm; P(:,:,:,5) = P(:,:,:,5) + Yngm; %if ~debug, clear Yngm; end
    %}
    
  elseif size(P,4)==4 || size(P,4)==3 % skull-stripped
    [Yb,Ybb,Yg,Ydiv,P] = cat_main_updateSPM_skullstriped(Ysrc,P,res,vx_vol,T3th);
  elseif job.extopts.gcutstr==0 
    [Yb,Ybb,Yg,Ydiv] = cat_main_updateSPM_gcut0(Ysrc,P,vx_vol,T3th);
  elseif job.extopts.gcutstr==2
    [Yb,Ybb,Yg,Ydiv] = cat_main_APRG(Ysrc,P,res,T3th);
  else
    [Yb,Ybb,Yg,Ydiv] = cat_main_updateSPM_gcutold(Ysrc,P,res,vx_vol,T3th);
  end

  
  % RD202101: cleanup+ remove small unconnected components
  %{
    Pgm = P(:,:,:,1); 
    Ybb = cat_vol_ctype(cat_vol_morph(Yb,'de',4,vx_vol)); 
    Ym  = Pgm>192; Ygc = cat_vol_morph(Ym | Ybb,'l',[inf,27])>0; Pgm(Ym(:)) = Pgm(Ym(:)) .* cat_vol_ctype(Ygc(Ym(:))); 
    Ym  = Pgm>128; Ygc = cat_vol_morph(Ym | Ybb,'l',[inf,27])>0; Pgm(Ym(:)) = Pgm(Ym(:)) .* cat_vol_ctype(Ygc(Ym(:))); 
    Ym  = Pgm> 64; Ygc = cat_vol_morph(Ym | Ybb,'l',[inf,27])>0; Pgm(Ym(:)) = Pgm(Ym(:)) .* cat_vol_ctype(Ygc(Ym(:))); 
    Ym  = Pgm>  8; Ygc = cat_vol_morph(Ym | Ybb,'l',[inf,27])>0; Pgm(Ym(:)) = Pgm(Ym(:)) .* cat_vol_ctype(Ygc(Ym(:))); 
    Ygc = cat_vol_morph(Pgm> 64,'l',[inf,27])>0; Pgm = Pgm .* cat_vol_ctype(Ygc | Ybb); 
    P(:,:,:,2) = P(:,:,:,2) + ( (P(:,:,:,1) - Pgm) .* cat_vol_ctype(Ysrc>T3th(2))); % add to WM
    P(:,:,:,3) = P(:,:,:,3) + ( (P(:,:,:,1) - Pgm) .* cat_vol_ctype(Ysrc<T3th(2))); % add to CSF
    P(:,:,:,1) = Pgm; 
  %} 
  
  
  
  %% save brainmask using SPM12 segmentations for later use
  if ~exist('Ym0','var')
    Ym0 = single(P(:,:,:,3))/255 + single(P(:,:,:,1))/255 + single(P(:,:,:,2))/255;
  end
  Yb0 = (Ym0 > min(0.5,max(0.25, job.extopts.gcutstr))); clear Ym0
  Yb0 = cat_vol_morph(cat_vol_morph(Yb0,'lo'),'c');

  
  
% RD202010: In some images SPM selects the image BG and brain tisssue as class 4  
%{
  volcls4 = sum(sum(sum( single(P(:,:,:,4)>64) .* (Yb>0.5) ))) .* prod(vx_vol)/1000; 
  volcls5 = sum(sum(sum( single(P(:,:,:,5)>64) .* (Yb>0.5 & (Ysrc>=mean(T3th(1:2)) & Ysrc<T3th(3) + diff(T3th(2:3))) ) ))) .* prod(vx_vol)/1000; 
  volcls6 = sum(sum(sum( single(P(:,:,:,6)>64) .* (Yb>0.5 & (Ysrc>=mean(T3th(1:2)) & Ysrc<T3th(3) + diff(T3th(2:3))) ) ))) .* prod(vx_vol)/1000; 
  scls = setdiff( unique( [ 4*(volcls4>10) 5*(volcls5>10) 6*(volcls6>10) ] ) , 0); 
  if volcls4 > 10 || volcls5 > 10 || volcls6 > 10
    PN{1} = cat_vol_ctype( single(P(:,:,:,1)) + single(sum(P(:,:,:,scls),4)) .* (Yb>0.5) .* ...
           (Ysrc>=mean(T3th(1:2)) & Ysrc<mean(T3th(2:3))) ); % GM
    PN{2} = cat_vol_ctype( single(P(:,:,:,2)) + single(sum(P(:,:,:,scls),4)) .* (Yb>0.5) .* ...
           (Ysrc>=mean(T3th(2:3)) & Ysrc<T3th(3) + 0.25*diff(T3th(2:3)))  ); % WM ... be carefull with blood vessels and tumors
    PN{3} = cat_vol_ctype( single(P(:,:,:,3)) + single(sum(P(:,:,:,scls),4)) .* (Yb>0.5) .* ...
           (Ysrc< mean(T3th(1:2))) ); % CSF
    PN{4} = cat_vol_ctype( single( P(:,:,:,4) ) .* (Ysrc<T3th(3) + 0.25*diff(T3th(2:3))) .* (Yb<0.5) );
    PN{5} = cat_vol_ctype( single( P(:,:,:,5) ) .* (Ysrc<T3th(3) + 0.25*diff(T3th(2:3))) .* (Yb<0.5) );
    PN{6} = cat_vol_ctype( single( P(:,:,:,6) ) .* (Ysrc<T3th(3) + 0.25*diff(T3th(2:3))) .* (Yb<0.5) );
    %%
    for pi=1:6, P(:,:,:,pi) = PN{pi}; end
    clear PN pi; 
  end 
%}




  %% RD202110: Background correction in longitidunal mode
  %  We observed some problems in the SPM background segmentation for
  %  longidutidnal processing that detected the volume of the boundary box 
  %  whereas the real background was miss-aligned to class 5 that caused 
  %  further problems in the LAS function that were solved too. Although 
  %  it would be possible to adapt the SPM segmentation, eg. by adapting 
  %  the number of gaussians per class, we decided that it is simpler and 
  %  maybe saver to add further test in the longitudinal case, where the 
  %  TPM should be close to the segmentation outcome.  
  if (isfield(job,'useprior') && ~isempty(job.useprior) && strcmp(job.opts.affreg,'prior') ) && ... 
     (isfield(res,'ppe') && ~res.ppe.affreg.highBG)
    %% sum of all TPM classes without background
    Vall = tpm.V(end); Vall.pinfo(3) = 0; Vall.dt=16; 
    Vall.dat = zeros(size(tpm.dat{1})); for k1 = 1:numel(tpm.dat)-1, Vall.dat = Vall.dat + single(exp(tpm.dat{k1})); end 
    Yall = cat_vol_sample(res.tpm(1),Vall,Yy,1);

    % backgound class
    Ybg = 1 - Yall; clear Yall Vall; 

    %% estimate error and do correction 
    rmse = @(x,y) mean( (x(:) - y(:)).^2 ).^0.5; 
    % the TPM BG may be smaller due to the limited overlap and we need a
    % higher threshold to avoid unnecessary background corrections
    TPisSmaller = ( sum(sum(sum(single(P(:,:,:,end))/255))) - sum(Ybg(:))) < 0;  
    if rmse(Ybg,single(P(:,:,:,end))/255) > 0.3 + 0.2*TPisSmaller % just some threshold (RD20220103: adjusted by TPisSmaller)
      % setup new background
      Ynbg = Ybg>0.5 | ( P(:,:,:,end)>128 & Ysrc < mean( T3th(1:2) ) );
      Ynbg = cat_vol_morph(Ynbg,'dc',5,vx_vol); 
      Ynbg = uint8( 255 .* smooth3(Ynbg) ); 
      
      % correct classes
      for k1 = 1:size(P,4)-1, P(:,:,:,k1) = P(:,:,:,k1) - min(P(:,:,:,k1),Ynbg); end
      P(:,:,:,end) = max( Ynbg , P(:,:,:,end) ); 
      clear Ynbg; 
      
      % normalize all classes
      sP = (sum(single(P),4)+eps)/255;
      for k1=1:size(P,4), P(:,:,:,k1) = cat_vol_ctype(single(P(:,:,:,k1))./sP); end
      clear sP; 
    
      cat_io_addwarning('cat_main_updateSPM:ReplacedLongBackground','Detected and corrected inadequate background \\nsegmentation in longitudinal mode.',0,[1 2]);
    end
    clear Ybg; 
  end

  
  
  
  stime2 = cat_io_cmd('  Update probability maps','g5','',job.extopts.verb-1,stime2);
  if ~(any(sign(diff(T3th))==-1)) && ...
     ~( (isfield(job,'useprior') && ~isempty(job.useprior) && strcmp(job.opts.affreg,'prior') ) && ... % no single longitudinal timepoint
        (isfield(res,'ppe') && ~res.ppe.affreg.highBG) )
    %% Update probability maps
    % background vs. head - important for noisy backgrounds such as in MT weighting
    if size(P,4)==4 || size(P,4)==3 % skull-stripped
      Ybg = ~Yb;
    else
      if sum(sum(sum(P(:,:,:,6)>240 & Ysrc<cat_stat_nanmean(T3th(1:2)))))>10000
        Ybg = P(:,:,:,6); 
        [Ybgr,Ysrcr,resT2] = cat_vol_resize({Ybg,Ysrc},'reduceV',vx_vol,2,32); 
        Ybgrth = max(cat_stat_nanmean(Ysrcr(Ybgr(:)>128)) + 2*std(Ysrcr(Ybgr(:)>128)),T3th(1));
        Ybgr = cat_vol_morph(cat_vol_morph(cat_vol_morph(Ybgr>128,'d') & Ysrcr<Ybgrth,'lo',1),'lc',1);
        Ybg  = cat_vol_resize(cat_vol_smooth3X(Ybgr,1),'dereduceV',resT2); 
        clear Ysrcr Ybgr; 
      else
        Ybg = ~Yb;
      end
    end
    %%
    P4   = cat_vol_ctype( single(P(:,:,:,6)) .* (Ysrc<T3th(2))  .* (Ybg<0.5) + single(P(:,:,:,4)) .* (Ybg<0.5) ); % remove air in head
    P5   = cat_vol_ctype( single(P(:,:,:,6)) .* (Ysrc>=T3th(2)) .* (Ybg<0.5) + single(P(:,:,:,5)) .* (Ybg<0.5) ); % remove air in head
    P6   = cat_vol_ctype( single(sum(P(:,:,:,4:5),4)) .* (Ybg>0.5) + single(P(:,:,:,6)) .* (Ybg>0.5) ); % add objects/artifacts to background
    P(:,:,:,4) = P4;
    P(:,:,:,5) = P5;
    P(:,:,:,6) = P6;
    clear P4 P5 P6;
    
    sP = (sum(single(P),4)+eps)/255;
    for k1=1:size(P,4), P(:,:,:,k1) = cat_vol_ctype(single(P(:,:,:,k1))./sP); end
    
    %% correct probability maps to 100% 
    sumP = cat_vol_ctype(255 - sum(P(:,:,:,1:6),4));
    P(:,:,:,1) = P(:,:,:,1) + sumP .* uint8( Ybg<0.5  &  Yb & Ysrc>cat_stat_nanmean(T3th(1:2)) & Ysrc<cat_stat_nanmean(T3th(2:3)));
    P(:,:,:,2) = P(:,:,:,2) + sumP .* uint8( Ybg<0.5  &  Yb & Ysrc>=cat_stat_nanmean(T3th(2:3)));
    P(:,:,:,3) = P(:,:,:,3) + sumP .* uint8( Ybg<0.5  &  Yb & Ysrc<=cat_stat_nanmean(T3th(1:2)));
    P(:,:,:,4) = P(:,:,:,4) + sumP .* uint8( Ybg<0.5  & ~Yb & Ysrc<T3th(2));
    P(:,:,:,5) = P(:,:,:,5) + sumP .* uint8( Ybg<0.5  & ~Yb & Ysrc>=T3th(2));
    P(:,:,:,6) = P(:,:,:,6) + sumP .* uint8( Ybg>=0.5 & ~Yb );
    clear Ybg sumP;

    
    %% head to WM 
    % Undercorrection of strong inhomogeneities in high field scans 
    % (>1.5T) can cause missalignments of the template and therefore 
    % miss classifications of the tissues that finally avoid further 
    % corrections in by LAS. 
    % Typically the alginment failed in this cases because the high 
    % intensities next to the head that were counted as head and not
    % corrected by SPM.
    % e.g. HR075, Magdeburg7T, SRS_SRS_Jena_DaRo81_T1_20150320-191509_MPR-08mm-G2-bw330-nbc.nii, ...
    Ywm = single(P(:,:,:,2)>128 & Yg<0.3 & Ydiv<0.03); Ywm(Ybb<128 | (P(:,:,:,1)>128 & abs(Ysrc/T3th(3)-2/3)<1/3) | Ydiv>0.03) = nan;
    [Ywm1,YD] = cat_vol_downcut(Ywm,1-Ysrc/T3th(3),0.02); Yb(isnan(Yb))=0; Ywm(YD<300)=1; Ywm(isnan(Ywm))=0; clear Ywm1 YD; %#ok<ASGLU>
    Ywmc = uint8(smooth3(Ywm)>0.7);
    Ygmc = uint8(cat_vol_morph(Ywmc,'d',2) & ~Ywmc & Ydiv>0 & Yb & cat_vol_smooth3X(Yb,8)<0.9 & Ysrc>mean(T3th(1:2)));
    P(:,:,:,[1,3:6]) = P(:,:,:,[1,3:6]) .* repmat(1-Ywmc,[1,1,1,5]);
    P(:,:,:,2:6)     = P(:,:,:,2:6)     .* repmat(1-Ygmc,[1,1,1,5]);
    P(:,:,:,1)       = max(P(:,:,:,1),255*Ygmc);
    P(:,:,:,2)       = max(P(:,:,:,2),255*Ywmc);
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
    clear Ygmc Ywmc; 

    
    %% head to GM ... important for children
    [Ywmr,Ybr,resT2] = cat_vol_resize({Ywm,Yb},'reduceV',vx_vol,2,32); 
    Ygm = cat_vol_morph(Ywmr>0.5,'d',3) & (cat_vol_morph(~Ybr,'d',3) | cat_vol_morph(Ybr,'d',1)); clear Ybr Ywmr;  % close to the head
    Ygm = cat_vol_resize(single(Ygm),'dereduceV',resT2)>0.5;
    Ygm = Ygm & Yp0<2/3 & Yb & Yg<cat_stat_nanmean(Yg(P(:,:,:,1)>64)) & Ydiv<cat_stat_nanmean(Ydiv(P(:,:,:,1)>64)); % add GM with low SPM prob ... 
    Ygm = Ygm & (Ysrc>cat_stat_nansum(T3th(1:2).*[0.5 0.5])) & (Ysrc<cat_stat_nansum(T3th(2:3).*[0.2 0.8])); % but good intensity
    Ygm(smooth3(Ygm)<0.5)=0; 
    clear Ydiv;
    Ygm = uint8(Ygm); 
    P(:,:,:,5) = P(:,:,:,5) .* (1-Ygm);
    P(:,:,:,3) = P(:,:,:,3) .* (1-Ygm);
    P(:,:,:,2) = P(:,:,:,2) .* (1-Ygm);
    P(:,:,:,1) = cat_vol_ctype(single(P(:,:,:,1)) + 255*single(Ygm));
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
    clear Ywm Ygm;

    %% RD20221108: update CSF to reduce problems in MP2rage
    if res.isMP2RAGE
      Ycsf = uint8( (Ysrc < min( res.mn(res.lkp==3)) * 0.7 + 0.3 * min( res.mn(res.lkp==1)) )  &  Yb  & Yg<.3);
      for ti = setdiff(1:size(P,4),3), P(:,:,:,ti) = P(:,:,:,ti) .* (1-Ycsf); end
      P(:,:,:,3) = cat_vol_ctype(single(P(:,:,:,3)) + 255*single(Ycsf));
      Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
    end    

    %% remove brain tissues outside the brainmask ...
    %  tissues > skull (within the brainmask)
    Yhdc = uint8(smooth3( Ysrc/T3th(3).*(Ybb>cat_vol_ctype(0.2*255)) - Yp0 )>0.5); 
    sumP = sum(P(:,:,:,1:3),4); 
    P(:,:,:,4)   =  cat_vol_ctype( single(P(:,:,:,4)) + sumP .* ((Ybb<=cat_vol_ctype(0.05*255)) | Yhdc ) .* (Ysrc<T3th(2)));
    P(:,:,:,5)   =  cat_vol_ctype( single(P(:,:,:,5)) + sumP .* ((Ybb<=cat_vol_ctype(0.05*255)) | Yhdc ) .* (Ysrc>=T3th(2)));
    P(:,:,:,1:3) =  P(:,:,:,1:3) .* repmat(uint8(~(Ybb<=cat_vol_ctype(0.05*255)) | Yhdc ),[1,1,1,3]);
    clear sumP Yp0 Yhdc; 
  end
  clear Ybb;
  
  
  

  %% MRF
  % Used spm_mrf help and tested the probability TPM map for Q without good results.         
  nmrf_its = 0; % 10 iterations better to get full probability in thin GM areas 
  spm_progress_bar('init',nmrf_its,['MRF: Working on ' nam],'Iterations completed');
  if isfield(res,'mg'), Kb = max(res.lkp); else, Kb = size(res.intensity(1).lik,2); end
  G   = ones([Kb,1],'single');
  vx2 = single(sum(res.image(1).mat(1:3,1:3).^2));
  % P = zeros([d(1:3),Kb],'uint8');
  % P = spm_mrf(P,Q,G,vx2); % init: transfer data from Q to P 
  if 0
    %% use TPM as Q
    Q = zeros(size(P),'uint8');
    for di=1:6
      vol = cat_vol_ctype(spm_sample_vol(tpm.V(di),...
        double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0)*255,'uint8');
      Q(:,:,:,di) = reshape(vol,d);
    end
  end
  for iter=1:nmrf_its
      P = spm_mrf(P,single(P),G,vx2); % spm_mrf(P,Q,G,vx2);
      spm_progress_bar('set',iter);
  end

  %% update segmentation for error report
  Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
  [cat_err_res.init.Yp0,cat_err_res.init.BB] = cat_vol_resize(Yp0,'reduceBrain',vx_vol,2,Yp0>0.5); 
  cat_err_res.init.Yp0 = cat_vol_ctype(cat_err_res.init.Yp0/3*255);
 
  spm_progress_bar('clear');
  for k1=1:size(P,4)
      Ycls{k1} = P(:,:,:,k1); %#ok<AGROW>
  end
  clear Q P q q1 Coef b cr N lkp n wp M k1


  if job.extopts.verb>2
    % save information for debugging and OS test
    % input variables + bias corrected, bias field, class image
    % strong differences in bias fields can be the result of different 
    % registration > check 'res.image.mat' and 'res.Affine'
    [pth,nam] = spm_fileparts(res.image0(1).fname); 
    tpmci  = 1;
    tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',tpmci,'postbias'));
    save(tmpmat,'res','tpm','job','Ysrc','Ybf','Ycls');
  end
 
  clear Ybf
  stime2 = cat_io_cmd(' ','g5','',job.extopts.verb-1,stime2); 
  fprintf('%5.0fs\n',etime(clock,stime));

  
  % some reports 
  for i=1:numel(Ycls), res.ppe.SPMvols1(i) = cat_stat_nansum(single(Ycls{i}(:)))/255 .* prod(vx_vol) / 1000; end
  
  % display  some values for developers
  if job.extopts.expertgui > 1
    % ... I want to add the intensities later
    %cat_io_cprintf('blue',sprintf('    SPM  volumes (CGW=TIV; in mm%s):      %6.2f + %6.2f + %6.2f = %4.0f\n',...
    %  native2unicode(179, 'latin1'),res.ppe.SPMvols0([3 1 2]),sum(res.ppe.SPMvols0(1:3))));    
    if isfield(job.extopts,'spm_kamap') && job.extopts.spm_kamap 
      cat_io_cprintf('blue',sprintf('    SPM  volumes (CGW=TIV; in mm%s):     %7.2f +%7.2f +%7.2f = %4.0f\n',...
        native2unicode(179, 'latin1'),res.ppe.SPMvols0([3 1 2]),sum(res.ppe.SPMvols0(1:3))));    
      cat_io_cprintf('blue',sprintf('    AMAP volumes (CGW=TIV; in mm%s):     %7.2f +%7.2f +%7.2f = %4.0f\n',...
        native2unicode(179, 'latin1'),res.ppe.SPMvols1([3 1 2]),sum(res.ppe.SPMvols1(1:3))));    
    else
      cat_io_cprintf('blue',sprintf('    SPM volumes pre  (CGW=TIV; in mm%s): %7.2f +%7.2f +%7.2f = %4.0f\n',...
        native2unicode(179, 'latin1'),res.ppe.SPMvols0([3 1 2]),sum(res.ppe.SPMvols0(1:3)))); 
      cat_io_cprintf('blue',sprintf('    SPM volumes post (CGW=TIV; in mm%s): %7.2f +%7.2f +%7.2f = %4.0f\n',...
        native2unicode(179, 'latin1'),res.ppe.SPMvols1([3 1 2]),sum(res.ppe.SPMvols1(1:3)))); 
    end
  end
  


end
function [Yb,Ybb,Yg,Ydiv,P] = cat_main_updateSPM_skullstriped(Ysrc,P,res,vx_vol,T3th)
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
    Yb   = Yp0>=0.5/3 & Ysrc>0; 
    Ybb  = cat_vol_ctype(Yb)*255; 

    P(:,:,:,6) = P(:,:,:,4); 
    P(:,:,:,4) = zeros(size(Ysrc),'uint8');
    P(:,:,:,5) = zeros(size(Ysrc),'uint8'); 
    res.lkp = [res.lkp 5 6];
    res.mn  = [res.mn(1:end-1),0,0,0];
    res.mg  = [res.mg(1:end-1);1;1;1];
    res.vr(1,1,numel(res.lkp)-1:numel(res.lkp)) = 0;
     
    [Ysrcb,BB] = cat_vol_resize(Ysrc,'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yp0>1/3); clear Yp0;
    Yg   = cat_vol_grad(Ysrcb/T3th(3),vx_vol);
    Ydiv = cat_vol_div(Ysrcb/T3th(3),vx_vol);
    Yg   = cat_vol_resize(Yg ,'dereduceBrain',BB);
    Ydiv = cat_vol_resize(Ydiv ,'dereduceBrain',BB);
end
function [Yb,Ybb,Yg,Ydiv] = cat_main_updateSPM_gcut0(Ysrc,P,vx_vol,T3th)
    % brain mask
    Ym   = single(P(:,:,:,3))/255 + single(P(:,:,:,1))/255 + single(P(:,:,:,2))/255;
    Yb   = (Ym > 0.5);
    Yb   = cat_vol_morph(cat_vol_morph(Yb,'lo'),'c');
    Ybb  = cat_vol_ctype(cat_vol_smooth3X(Yb,2)*256); 

    [Ysrcb,BB] = cat_vol_resize({Ysrc},'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yb);
    Yg   = cat_vol_grad(Ysrcb/T3th(3),vx_vol);
    Ydiv = cat_vol_div(Ysrcb/T3th(3),vx_vol);
    Yg   = cat_vol_resize(Yg   ,'dereduceBrain',BB);
    Ydiv = cat_vol_resize(Ydiv ,'dereduceBrain',BB);
end
function [Yb,Ybb,Yg,Ydiv] = cat_main_updateSPM_gcutold(Ysrc,P,res,vx_vol,T3th)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T1 only > remove in future if gcut is removed too!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clsint = @(x) round( sum(res.mn(res.lkp==x) .* res.mg(res.lkp==x)') * 10^5)/10^5;
    voli   = @(v) (v ./ (pi * 4./3)).^(1/3);   % volume > radius
    Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));

    % old skull-stripping
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
    brad = voli(sum(Yp0(:)>0.5).*prod(vx_vol)/1000); 
    [Ysrcb,Yp0,BB] = cat_vol_resize({Ysrc,Yp0},'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yp0>1/3);
    %Ysrcb = max(0,min(Ysrcb,max(T3th)*2));
    BGth = min(cat_stat_nanmean(Ysrc( P(:,:,:,6)>128 )),clsint(6));
    Yg   = cat_vol_grad((Ysrcb-BGth)/diff([BGth,T3th(3)]),vx_vol);
    Ydiv = cat_vol_div((Ysrcb-BGth)/diff([BGth,T3th(3)]),vx_vol);
    Ybo  = cat_vol_morph(cat_vol_morph(Yp0>0.3,'lc',2),'d',brad/2/mean(vx_vol)); 
    BVth = diff(T3th(1:2:3))/abs(T3th(3))*1.5; 
    RGth = diff(T3th(2:3))/abs(T3th(3))*0.1; 
    Yb   = single(cat_vol_morph((Yp0>1.9/3) | (Ybo & Ysrcb>mean(T3th(2)) & ...
           Ysrcb<T3th(3)*1.5 & Yg<0.5),'lo',max(0,0.6/mean(vx_vol)))); 
    
    %% region-growing GM 1
    Yb(~Yb & (~Ybo | Ysrcb<cat_stat_nanmean(T3th(2)) | Ysrcb>cat_stat_nanmean(T3th(3)*1.2) | Yg>BVth))=nan;
    [Yb1,YD] = cat_vol_downcut(Yb,Ysrcb/T3th(3),RGth); clear Yb1; Yb(isnan(Yb))=0; Yb(YD<400/mean(vx_vol))=1; clear YD; %#ok<ASGLU> 
    Yb(smooth3(Yb)<0.5)=0; Yb = single(Yb | (Ysrcb>T3th(1) & Ysrcb<1.2*T3th(3) & cat_vol_morph(Yb,'lc',4)));
    
    %% region-growing GM 2
    Yb(~Yb & (~Ybo | Ysrcb<T3th(1) | Ysrcb>cat_stat_nanmean(T3th(3)*1.2) | Yg>BVth))=nan;
    [Yb1,YD] = cat_vol_downcut(Yb,Ysrcb/T3th(3),RGth/2); clear Yb1; Yb(isnan(Yb))=0; Yb(YD<400/mean(vx_vol))=1; clear YD; %#ok<ASGLU>
    Yb(smooth3(Yb)<0.5)=0; Yb = single(Yb | (Ysrcb>T3th(1) & Ysrcb<1.2*T3th(3) & cat_vol_morph(Yb,'lc',4)));
    
    %% region-growing GM 3
    Yb(~Yb & (~Ybo | Ysrcb<mean([BGth,T3th(1)]) | Ysrcb>cat_stat_nanmean(T3th(3)*1.2) | Yg>BVth))=nan; clear Ybo;
    [Yb1,YD] = cat_vol_downcut(Yb,Ysrcb/T3th(3),RGth/10); clear Yb1; Yb(isnan(Yb))=0; Yb(YD<400/mean(vx_vol))=1; clear YD; %#ok<ASGLU>
    Yb(smooth3(Yb)<0.5)=0; Yb(Yp0toC(Yp0*3,1)>0.9 & Yg<0.3 & Ysrcb>BGth & Ysrcb<T3th(2)) = 1; 
    
    %% ventricle closing
    [Ybr,Ymr,resT2] = cat_vol_resize({Yb>0,Ysrcb/T3th(3)},'reduceV',vx_vol,2,32); clear Ysrcb
    Ybr = Ybr | (Ymr<0.8 & cat_vol_morph(Ybr,'lc',6)); clear Ymr;  % large ventricle closing
    Ybr = cat_vol_morph(Ybr,'lc',2);                 % standard closing
    Yb  = Yb | cat_vol_resize(cat_vol_smooth3X(Ybr,2),'dereduceV',resT2)>0.7; clear Ybr
    Yb  = smooth3(Yb)>0.5; 
    Ybb = cat_vol_ctype(cat_vol_smooth3X(Yb,2)*255); 
    Yb   = cat_vol_resize(Yb   , 'dereduceBrain' , BB);
    Ybb  = cat_vol_resize(Ybb  , 'dereduceBrain' , BB);
    Yg   = cat_vol_resize(Yg   , 'dereduceBrain' , BB);
    Ydiv = cat_vol_resize(Ydiv , 'dereduceBrain' , BB);
    clear Ysrcb Ybo;
end
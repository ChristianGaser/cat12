function cat_tst_prepareManualSegmentation
% _________________________________________________________________________
% This script prepare data for semi-manual tissue segmenation. 
% As far as manual segmenation is a very time consuming process, we want to 
%   (a) focus only on a limit numer of slices and 
%   (b) we need use automatic preprocessing routines to give some intial 
%       segmentation that is not related to any standard preprocessing to 
%       avoid a bias. 
%
% Therefore, we used SPM to create a bias corrected images, create a tissue 
% segmenation (that is used for brain masking) and to estimate the average
% tissue intensities that were used for intensity normalization of the bias
% corrected image. Furthermore, a noise correction (ISARNLM) is applyed. 
% 
% There are manual control variables (mymask,myint,mypoints) to overwrite
% the standard/estimated values.
%
% _________________________________________________________________________
% Robert Dahnke
% $Id$


  % select images and result directory
  Pm     = cellstr(spm_select([1 inf],'image','Select bias corrected SPM images','','','^m.*'));
  resdir = '/Volumes/catDB/Tracing2/';
  if isempty(Pm), return; end

  job.isarnlm = 1; 
  
  
  %% Manual paraemter
  slicepoints = [ 27 -45   5]; % default points as nx3 matrix ...
  nslicepoints = size(slicepoints,1);
  
  % use smoothing to create softer segments in case of noisy data
  mysmooth = {
    'mWB01_T1' 1,
  };
  
  % skull-stripped data may failed
  mymask   = { % ff 
    'mINDI_HC_LPZ_sub00321_T1_SD000000-RS00' 
    'mINDI_HC_NHa_sub10033_T1_SD000000-RS00' 
    'mINDI_HC_NHb_sub01183_T1_SD000000-RS00' 
    'mINDI_HC_OGB_sub05191_T1_SD000000-RS00' 
    'mINDI_HC_PIT_sub01891_T1_SD000000-RS00' 
    };
  
  % manual intensity values
  myint    = { % filename [background CSF lowGM highGM WM]; ignore if [x,y,z] is empty
    'mBUSS_2002_1YO_t1'                       [0.03  0.40  0.70  0.80  1.00] *   214.9821;
    'mINDI_HC_AAa_sub04111_T1_SD000000-RS00'  [0.03  0.30  0.65  0.80  1.00] *   570.7427;
    'mINDI_HC_AAb_sub00306_T1_SD000000-RS00'  [0.03  0.30  0.60  0.80  1.00] *   757.8686;
    'mINDI_HC_LPZ_sub00321_T1_SD000000-RS00'  [0.03  0.30  0.65  0.83  1.00] *   203.9154;
    'mINDI_HC_DAL_sub04288_T1_SD000000-RS00'  [0.03  0.30  0.55  0.85  1.05] * 3.8388e+03;
    'mINDI_HC_MLb_sub00917_T1_SD000000-RS00'  [0.06  0.20  0.60  0.70  1.00] * 1.2759e+03;
    'mINDI_HC_NHa_sub10033_T1_SD000000-RS00'  [0.00  0.40  0.65  0.80  1.00] *   354.9029;
    'mINDI_HC_NHb_sub01183_T1_SD000000-RS00'  [0.00  0.20  0.50  0.75  1.00] *   302.5603;
    'mINDI_HC_NWK_sub13411_T1_SD000000-RS00'  [0.00  0.15  0.55  0.75  1.00] *   263.5865;
    'mINDI_HC_OGB_sub05191_T1_SD000000-RS00'  [0.00  0.30  0.80  1.00  1.15] * 2.7279e+03;
    'mINDI_HC_OUL_sub01077_T1_SD000000-RS00'  [0.00  0.15  0.45  0.75  1.00] *   895.8569;
    'mINDI_HC_PAL_sub04856_T1_SD000000-RS00'  [0.06  0.25  0.45  0.80  1.00] * 4.5968e+03;
    'mINDI_HC_PIT_sub01891_T1_SD000000-RS00'  [0.00  0.30  0.70  0.88  1.00] *   338.0128;
    }; 
  
  % manual AC correction 
  mypoints = { % filename [x y z]; ignore if [x,y,z] is empty
    'mHR075_MPRAGE'                           []; 
    'm4397-tfl'                               slicepoints + repmat([  0   0  20],nslicepoints,1); % [ 27 -45  25];
    'mINDI_HC_AAa_sub04111_T1_SD000000-RS00'  slicepoints + repmat([  0   0   5],nslicepoints,1);
    'mINDI_HC_BNG_sub00031_T1_SD000000-RS00'  slicepoints + repmat([  0   0   5],nslicepoints,1); % [ 27 -45  10];
    'mINDI_HC_CAM_sub00156_T1_SD000000-RS00'  slicepoints + repmat([  0   0   5],nslicepoints,1); % [ 27 -45   5];
    'mINDI_HC_DAL_sub04288_T1_SD000000-RS00'  slicepoints + repmat([ -2   5  15],nslicepoints,1); % [ 25 -40  20]; 
    'mINDI_HC_LDa_sub01553_T1_SD000000-RS00'  slicepoints + repmat([  0   0   5],nslicepoints,1); % [ 27 -45  10];  
    'mINDI_HC_LDb_sub01787_T1_SD000000-RS00'  slicepoints + repmat([ -3   0   4],nslicepoints,1); % [ 24 -45  10]; 
    'mINDI_HC_LPZ_sub00321_T1_SD000000-RS00'  slicepoints + repmat([ -4   3  10],nslicepoints,1); % [ 23 -42  15];
    'mINDI_HC_MLb_sub00917_T1_SD000000-RS00'  slicepoints + repmat([ -5   0  10],nslicepoints,1); % [ 22 -45  15];
    'mINDI_HC_NHa_sub10033_T1_SD000000-RS00'  slicepoints + repmat([ -2  -5   5],nslicepoints,1); % [ 25 -50  10];
    'mINDI_HC_NHb_sub01183_T1_SD000000-RS00'  slicepoints + repmat([-37  25 -55],nslicepoints,1); % [-10 -20 -50];
    'mINDI_HC_NWK_sub13411_T1_SD000000-RS00'  slicepoints + repmat([ -5   5  20],nslicepoints,1); % [ 22 -40  25];
    'mINDI_HC_NYa_sub01912_T1_SD000000-RS00'  slicepoints + repmat([  0   0   5],nslicepoints,1); % [ 27 -45  10]; 
    'mINDI_HC_OGB_sub05191_T1_SD000000-RS00'  slicepoints + repmat([ -5   0  60],nslicepoints,1); % [ 22 -45  65]; 
    'mINDI_HC_STL_sub02115_T1_SD000000-RS00'  slicepoints + repmat([  0   0  05],nslicepoints,1); % [ 22 -45  65]; 
    };
  
  
  
  %% Processing
  for pmi = 1:numel(Pm)
    
    %% -- load data -------------------------------------------------------
    [pp,ff] = spm_fileparts(Pm{pmi});
    [~,ppp] = spm_fileparts(pp);  
    rpp  = fullfile(resdir,ppp,ff(2:end)); 
    if ~exist(rpp,'dir'), mkdir(rpp); end
    Pmat = fullfile(pp,sprintf('%s_seg8.mat',ff(2:end)));
    Pp0  = fullfile(pp,sprintf('p0%s.nii',ff(2:end)));

    res = load(Pmat); 

    Vm  = spm_vol(Pm{pmi});
    Vp0 = spm_vol(Pp0);

    Ysrc = single(spm_read_vols(Vm)); 
    Yp0  = single(spm_read_vols(Vp0)); 

    vx_vol  = sqrt(sum(Vm.mat(1:3,1:3).^2));   

    fprintf('%s: ',ff);

    %% -- intensity normalization -----------------------------------------
    id = find(cellfun('isempty',strfind(myint(:,1),ff))==0);
    if isempty(id) || isempty(myint(id,2))
      T3ths = [ ...
             min(Ysrc(:)) ...
             cat_stat_nanmean(res.mn(res.lkp==6 & res.mg'>0.3)) ...
             min([  cat_stat_nanmean(res.mn(res.lkp==1 & res.mg'>0.3)) - ...
              diff([cat_stat_nanmean(res.mn(res.lkp==1 & res.mg'>0.3)),cat_stat_nanmean(res.mn(res.lkp==2 & res.mg'>0.3))]),...
             cat_stat_nanmean(res.mn(res.lkp==3 & res.mg'>0.3))]) ... CSF
             cat_stat_nanmean(res.mn(res.lkp==1 & res.mg'>0.3)) ... GMth
             cat_stat_nanmean(res.mn(res.lkp==2 & res.mg'>0.3)) ... WMth
             cat_stat_nanmean(res.mn(res.lkp==2 & res.mg'>0.3)) + ... WM+
             abs(diff(max(...
              [cat_stat_nanmean(res.mn(res.lkp==1 & res.mg'>0.3)),cat_stat_nanmean(res.mn(res.lkp==2 & res.mg'>0.3))],...
              [cat_stat_nanmean(res.mn(res.lkp==1 & res.mg'>0.3)),cat_stat_nanmean(res.mn(res.lkp==3 & res.mg'>0.3))]/2)*1.5)) ...
             max(Ysrc(:))];
      T3thx = [0,0.05,1,2,3,4,5];
    else
      T3ths = [min(Ysrc(:)),myint{id,2}(1:5),myint{id,2}(5) + abs(diff(myint{id,2}(3:2:5))),max(Ysrc(:))];
      T3thx = [0,0.05,1,1.66,2.33,3,4,5];
    end
    
    [T3ths,si] = sort(T3ths);
    T3thx     = T3thx(si);
    Ym = Ysrc+0; 
    for i=numel(T3ths):-1:2
      M = Ysrc>T3ths(i-1) & Ysrc<=T3ths(i);
      Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3ths(i-1))/diff(T3ths(i-1:i))*diff(T3thx(i-1:i));
    end
    M  = Ysrc>=T3ths(end); 
    Ym(M(:)) = numel(T3ths)/6 + (Ysrc(M(:)) - T3ths(i))/diff(T3ths(end-1:end))*diff(T3thx(i-1:i));    

    
    % create brainmask
    id = find(cellfun('isempty',strfind(mymask(:,1),ff))==0); %#ok<EFIND>
    if isempty(id) 
      Yb = smooth3(Yp0>0.5 | (Yp0>0.2 & Ym<1.5))>0.1; 
    else
      Yb = Ym>0; 
    end

    
    % copy of the original bias corrected
    copyfile(Vm.fname,fullfile(rpp,sprintf('b%s.nii',ff(2:end)))); 
    
    
    % ISARNLM noise correction and creation of output image
    fprintf('ISAR1 .. ',ff);
    if job.isarnlm, Ysrc = cat_vol_isarnlm(Ysrc,Vm,0,inf); end
    Vm2 = Vm; Vm2.fname = fullfile(rpp,sprintf('m%s.nii',ff(2:end)));
    spm_write_vol(Vm2,Ysrc); 

    
    % ISARNLM noise correction and creation of output image
    fprintf('ISAR2 .. ',ff);
    if job.isarnlm, Ym = cat_vol_isarnlm(Ym,Vm,0,inf); end
    Vm2 = Vm; Vm2.fname = fullfile(rpp,sprintf('n%s.nii',ff(2:end)));
    spm_write_vol(Vm2,Ym); 

    
    
    
    %% -- create slice mask ----------------------------------------------- 
    mati = spm_imatrix(Vm.mat); %*res.Affine); % .*sign(mati(7:9)) res.Affine
    %slicepointsubject = round((slicepoint - mati(1:3))./mati(7:9)); 
    id = find(cellfun('isempty',strfind(mypoints(:,1),ff))==0);
    if isempty(id) || isempty(mypoints{id,2})
      vmat = res.Affine * Vm.mat; vmat = inv(vmat); 
      slicepointsubject = round(vmat * [slicepoints .* sign(mati(7:9)) ,1]'); slicepointsubject = slicepointsubject(1:3)';
    else
      vmat = res.Affine * Vm.mat; vmat = inv(vmat); 
      slicepointsubject = round(vmat * [mypoints{id,2} .* sign(mati(7:9)) ,1]'); slicepointsubject = slicepointsubject(1:3)';
    end
    
    Yslicemask = false(size(Ym));
    for si=1:size(slicepointsubject,1) 
      Yslicemask(slicepointsubject(1),:,:) = true;
      Yslicemask(:,slicepointsubject(2),:) = true;
      Yslicemask(:,:,slicepointsubject(3)) = true;
    end

    % intensity-based segmentation 
    id = find(cellfun('isempty',strfind(mysmooth(:,1),ff))==0);
    if isempty(id) || isempty(mysmooth{id,2})
      Yp0pm = cat_vol_smooth3X(Ym,mysmooth{id,2}/mean(vx_vol));
      Yp0pm = round( max(1, Yp0pm ) ); 
    else
      Yp0pm = round( max(1, Ym ) ); 
    end
    Yp0pm(Yp0pm>3.2 | ~Yb) = 0; 
    Vp0m = Vm; Vp0m.fname = fullfile(rpp,sprintf('p0m%s.nii',ff(2:end)));
    spm_write_vol(Vp0m,Yp0pm); 
    Yp0pm(~Yslicemask ) = 0; 
    Vp0m = Vm; Vp0m.fname = fullfile(rpp,sprintf('p0s%s.nii',ff(2:end)));
    spm_write_vol(Vp0m,Yp0pm.*Yslicemask); 
    
    % display
    ds('l2','',vx_vol,Ym/3,Yb+6*Yslicemask,Ysrc/T3ths(end-2),round(Ym)/3,slicepointsubject(3)+1);
    T3ths/T3ths(end-2), T3ths(end-2)
    
    %%
    ds('l2','',vx_vol,Ym/3,Yb+6*Yslicemask,Ym/3,round( max(1,Ym) )/3,slicepointsubject(3)+1);
    fprintf('\n ');
    % slicemask 
    %system('/opt/local/lib/cmtk/bin/cmtk convertx 

  end
end







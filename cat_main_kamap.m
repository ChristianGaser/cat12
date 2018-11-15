function [P,res,stime2] = cat_main_kamap(Ysrc,Ycls,job,res,vx_vol,stime2)
% ______________________________________________________________________
%
% (initial) k-means AMAP segmentation
%
% This is an alternative pipeline in case of failed SPM brain tissue
% classification in datasets with abnormal anatomy, i.e. superlarge 
% ventricle. However, SPM is used for head tissue classification and
%  bias correction. It is (only) called from the cat_main function. 
%
%   [P,res,stime2] = cat_main_amap(Ysrc,Ycls,job,res,stime2)
%
%   Ysrc    .. cell structure of 3D probability maps of the segmentation 
%              [GM,WM,CSF,SK,HD,BG] 
%   P       .. 4D probability map of the segmentation
%   res     .. SPM segmentation structure
%   vx_vol  .. voxel_resolution
%   stime2  .. just a time stamp
% ______________________________________________________________________
%
% Robert Dahnke, Christian Gaser 
% Structural Brain Mapping Group
% University Jena 
% ______________________________________________________________________
% $Id$ 


% ______________________________________________________________________
% Possible developments:
% (1) Update for PD/T2 weightings?
%     Yes, in future!
% (2) Add head class model? 
%     Maybe later.
% (-) More complex skull-stripping? 
%     No, this is given by the input segmenation  
% ______________________________________________________________________


  % if there is a breakpoint in this file set debug=1 and do not clear temporary variables 
  dbs  = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

  
  stime2 = cat_io_cmd('  K-means AMAP bias-correction','g5','',job.extopts.verb-1,stime2);

  % simple skull-stripping
  Yb = Ycls{1} + Ycls{2} + Ycls{3};
  
  % In case of skull-stripping with have to combine SPM segmentation and 
  % the possible skull-stripping information 
  if (max(res.lkp)==4), Yb = Yb .* uint8(Ysrc>0); end % if skullstripped
  Yb = cat_vol_morph( smooth3(Yb)>64 , 'ldc', 8,vx_vol);

  
  %% intensity normalization
  Ym = Ysrc; Ym(~Yb) = nan; 
  [Ymic,th] = cat_stat_histth(Ym,0.98); 
  Ym(isnan(Ym) | Ym>(th(2) + diff(th)*2)) = 0;
  [T3th,T3sd,T3md] = kmeans3D(Ymic((Ym(:)>0)),5); clear Ymic;  %#ok<ASGLU>
  if T3md(end)<0.1, T3th(end)=[]; end
  Ym = (Ym - th(1)) ./ (T3th(end) - th(1)); 
  clear th; 

  
  %% bias correction for white matter (Yw) and ventricular CSF areas (Yv)
  Yw  = Ym>0.8 & Ym<1.5; Yw(smooth3(Yw)<0.6) = 0; Yw(smooth3(Yw)<0.6) = 0; 
  Yv  = Ym<0.5 & Yb; Yv  = cat_vol_morph(Yv,'e',4); 

  
  %% estimate value vth to mix CSF and WM information
  Yvw = cat_vol_smooth3X(Yv,6)>0.05 & cat_vol_morph(Yw,'e',2); 
  Ywv = cat_vol_smooth3X(Yw,6)>0.05 & cat_vol_morph(Yv,'e',2); 
  if sum(Yvw(:))>10 && sum(Yvw(:))>10
    vth = cat_stat_nanmedian(Ysrc(Yvw(:))) ./ cat_stat_nanmedian(Ysrc(Ywv(:)));
  elseif  sum(Yvw(:))>0 && sum(Yvw(:))>0
    vth = cat_stat_nanmedian(Ysrc(Yv(:))) ./ cat_stat_nanmedian(Ysrc(Yw(:)));
  else
    vth = 1; 
  end
  if debug==0, clear Ywv Yvw; end

  
  %% some precorrections and bias field approximation 
  Yi   = cat_vol_localstat(Ysrc .* Yw,Yw,1,3); 
  Yiv  = cat_vol_localstat(Ysrc .* Yv,Yv,1,1); 
  Yi   = Yi + Yiv * vth; if debug==0, clear Yiv; end
  Yi   = cat_vol_approx(Yi,'nn',vx_vol,4); 
  Yi   = cat_vol_smooth3X(Yi,4); 
  Ysrc = Ysrc ./ (Yi ./ median(Yi(Yw(:))));
  Ymi  = (Ysrc .* Yb) ./ Yi;

  
  %% remove of high intensity structures
  Ymi = cat_vol_median3(Ymi/median(Ymi(Yw(:))),Yb,Yb,0.2)*median(Ymi(Yw(:)));
  if debug==0, clear Ym Yi Yw Yv; end


  %%  similar to later AMAP call
  %   -------------------------------------------------------------------

  % correct for harder brain mask to avoid meninges in the segmentation
  Ymib = min(1.5,Ymi); Ymib(~Yb) = 0; 
  rf = 10^4; Ymib = round(Ymib*rf)/rf;

  %  prepare data for segmentation
  % more direct method ... a little bit more WM, less CSF
  Yp0 = cat_vol_ctype(max(1,min(3,round(Ymi * 3)))); Yp0(~Yb) = 0;

  % use index to speed up and save memory
  sz = size(Yb);
  [indx, indy, indz] = ind2sub(sz,find(Yb>0));
  indx = max((min(indx) - 1),1):min((max(indx) + 1),sz(1));
  indy = max((min(indy) - 1),1):min((max(indy) + 1),sz(2));
  indz = max((min(indz) - 1),1):min((max(indz) + 1),sz(3));

  % Yb source image because Amap needs a skull stripped image
  % set Yp0b and source inside outside Yb to 0
  Yp0b = Yp0(indx,indy,indz);  %#ok<NASGU>
  Ymib = Ymib(indx,indy,indz); 

  % adaptive mrf noise 
  if job.extopts.mrf>=1 || job.extopts.mrf<0; 
    % estimate noise
    [Yw,Yg] = cat_vol_resize({Ymi.*(Ycls{1}>240),Ymi.*(Ycls{2}>240)},'reduceV',vx_vol,3,32,'meanm');
    Yn = max(cat(4,cat_vol_localstat(Yw,Yw>0,2,4),cat_vol_localstat(Yg,Yg>0,2,4)),[],4);
    job.extopts.mrf = double(min(0.15,3*cat_stat_nanmean(Yn(Yn(:)>0)))) * 0.5; 
    clear Yn Yg
  end

  % display something
  stime2 = cat_io_cmd(sprintf('  Amap using k-means segmentations (MRF filter strength %0.2f)',...
    job.extopts.mrf),'g5','',job.extopts.verb-1,stime2);

  
  %% Amap parameters  - default sub=16 caused errors with highres data!
  % don't use bias_fwhm, because the Amap bias correction is not that efficient
  % and also changes intensity values
  Ymib = double(Ymib); n_iters = 50; sub = round(32/min(vx_vol)); %#ok<NASGU>
  n_classes = 3; pve = 5; bias_fwhm = 120; init_kmeans = 1;  %#ok<NASGU>
  if job.extopts.mrf~=0, iters_icm = 50; else iters_icm = 0; end %#ok<NASGU>

  
  % do segmentation and rep
  evalc(['prob = cat_amap(Ymib, Yp0b, n_classes, n_iters, sub, pve, init_kmeans, ' ...
    'job.extopts.mrf, vx_vol, iters_icm, bias_fwhm);']);
  clear Ymib Yp0b;
 
  
  % reorder probability maps according to spm order
  prob = prob(:,:,:,[2 3 1]);  %#ok<NODEF>

  
  % cleanup
  prob = cat_main_clean_gwc(prob,1);

  
  % update probability maps
  P = zeros([size(Ycls{1}) numel(Ycls)],'uint8');
  for i=1:numel(Ycls), P(:,:,:,i) = Ycls{i}; end
  for i=1:3
     P(:,:,:,i) = 0; P(indx,indy,indz,i) = prob(:,:,:,i);
  end
  clear prob;
  Ys = cat_vol_ctype(255 - sum(P(:,:,:,1:3),4));
  for i=4:size(P,4)
    P(:,:,:,i) = P(:,:,:,i) .* Ys; 
  end
  clear Ys;

  sP = (sum(single(P),4) + eps)/255;
  for k1=1:size(P,4), P(:,:,:,k1) = cat_vol_ctype(single(P(:,:,:,k1))./sP); end

  
  % update SPM segmentation information 
  for i=1:3
    [res.mn(res.lkp==i),tmp,res.mg(res.lkp==i)] = ...
      kmeans3D(Ysrc(P(:,:,:,i)>64),sum(res.lkp==i)); %#ok<ASGLU>
  end

end

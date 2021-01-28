function [Ym2,Ymi2,Tthm,Tthmi] = cat_main_update_intnorm(Ym,Ymi,Yb,Ycls,job,verb,indx,indy,indz)
%cat_main_update_intnorm. Fine intensity normalization in optimal images. 
% This is a temporary function to test the some concepts of improving the 
% intensity scaling in prepared images. This function expect optimal pre-
% processed data that is corrected for noise and inhomogeneities!
%
% - after cat_main_LAS:
%   [Ym2,Ymi2] = cat_main_update_intnorm(Ym,Ymi,Yb,Ycls,verb)
%
% - after AMAP with own peak estimation 
%   [Ym2,Ymi2] = cat_main_update_intnorm(Ym,Ymi,Yb,prob,verb,indx,indy,indz)
%
% - after AMAP with AMAP thresholds
%   [Ym2,Ymi2] = cat_main_update_intnorm(Ym,Ymi,AMAPths)
%
%   Ym[i][2] .. global/local intensity normalized input/output images
%   Yb       .. brain mask
%   Ycls     .. tissue classes as cell
%   prob     .. tissue classes as 4D matrix
%   ind[xyz] .. subregions defined by prob
%   verb     .. display some details (only for tests/development)
%   AMAPths  .. do not estimate peaks and just use the given 
%               (1x3 cell with 1x2 matrix with mean and std given by AMAP)
% _________________________________________________________________________
% Robert Dahnke, Structural Brain Mapping Group, University Jena, 202101
% $Id: cat_main_amap.m 1694 2020-08-19 13:42:00Z dahnke $ 

%
% ToDo: Test with PD/T2/FLAIR/MP2
%

  % be verbose for tests
  fprintf(' --- Update Ym and Ymi --- \n');

  def.extopts.ignoreErrors = 0; 
  job = cat_io_checkinopt(job,def); 
  if ~exist('verb','var'), verb = 0; end
  pve  = 0;  % estimate peaks also for the PVE boudary region (IN DEVELOPMENT)
  
  if numel(Yb) == 3
  % Final intensity scaling for local intensity normalized map Ymi that 
  % uses the AMAP threshholds.
    Tthamap = struct('T3th',{0:1/3:2},'T3thx',{3*[0 Yb{1}(1) Yb{2}(1) Yb{3}(1) 4/3:1/3:2] });
    Ymi2 = cat_main_gintnormi(Ymi,Tthamap);
    Ym2  = cat_main_gintnormi(Ym ,Tthamap);
    return
  end
  
  % thresholds variables and smaller brainmask (to avoid PVE to background)
  mith = zeros(1,3); mth = mith; 
  Ybe = cat_vol_morph(Yb,'e',2); 
  
  % The AMAP gives us a 4D-matrix, whereas SPM/cat_main uses a cell
  % structure with each tissue class.
  if ndims(Ycls)==4
    prob = Ycls; clear Ycls;  
    for i = 1:3
      Ycls{i} = zeros(size(Ym),'uint8'); 
      Ycls{i}(indx,indy,indz) = prob(:,:,:,i);
    end
    clear prob; 
  end
  
  % Different handling in case of other contrasts.
  % **** This needs more tests ****
  invers = cat_stat_nanmean(Ym(Ycls{2}(:)>0)) < cat_stat_nanmean(Ym(Ycls{3}(:)>0)); 
  if invers, thsel  = [0 1 2]; else, thsel = [0 2 1]; end 
  
  
  
  
  %% final intensity scaling for global intensity normalized map Ym
  %  The main aim is to estimate the best value that represents a tissue class. 
  %  This is in general the most typical value that occures most times and 
  %  that is the peak of the histogram. To avoid side effects by the PVE 
  %  it is usefull to use a mask for tissue with large/thick/compact volume
  %  such as the CSF and WM but not for the GM! In addition we expect a much
  %  smaller side peak, especially in the unmasked WM, where 5 peaks are 
  %  useful becuase there are also real local intensity differences in the 
  %  subcortical structures but also by the differently myelinated cortex. 
  %  For our highly corrected input the curve of the WM and CSF are not
  %  symetric and only include PVE torwards the CSF. 
  %
  %  Moreover, it is more save to use allways the peak with the larger GM
  %  distance, i.e., in T1 images the smalles/highest peak in the CSF/WM. 
  %  .. However, right now we are just using the highest one that is easier
  %     in case of different contrast 
  for i = 1:3 
    % masking
    if i==1
      Yclse = Ycls{i}>16; % just a little bit in the GM
    else
      if job.extopts.ignoreErrors < 3 
        Yclse = cat_vol_morph(Ycls{i}>240,'e'); % much more in the GM/CSF
      else
        Yclse = cat_vol_morph(Ycls{i}>16,'e'); % much more in the GM/CSF
      end
    end
    % assure that there are enough values (e.g. in children/olderly the
    % CSF/WM can be quite small 
    if sum(Yclse(:))<100, Yclse      = cat_vol_morph(Ycls{i}>128);     end
    if sum(Yclse(:))<50,  Yclse      = cat_vol_morph(Ycls{i}>0);       end
   
    % estimate the peaks
    if job.extopts.ignoreErrors < 3 
      [tth,sd,h] = cat_stat_kmeans(Ymi( Yclse & Ybe ),2 + 3*(i==1));  % 2 peaks for CSF/WM and 5 for GM
      if i==1, mith(i) = tth(h==max(h)); else, mith(i) = tth( thsel(i) ); end
      [tth,sd,h] = cat_stat_kmeans(Ym ( Yclse & Ybe ),2 + 3*(i==1)); 
      if i==1, mth(i) = tth(h==max(h)); else, mth(i) = tth( thsel(i) ); end
    else 
      % inoptimal preprocessing with some noise, bias, and more PVE
      % more classes but still inbalanced for CSF/WM (4 rather than the 5 GM classes)
      [tth,sd,h] = cat_stat_kmeans(Ymi( Yclse & Ybe ),4 + 1*(i==1));  % 2 peaks for CSF/WM and 5 for GM
      if i==1, mith(i) = tth(h==max(h)); else, mith(i) = tth( thsel(i) ); end
      [tth,sd,h] = cat_stat_kmeans(Ym ( Yclse & Ybe ),4 + 1*(i==1)); 
      if i==1, mth(i) = tth(h==max(h)); else, mth(i) = tth( thsel(i) ); end
    end
  end
  
  
  if pve
  %% PVE peaks (IN DEVELOPMENT)
  %  Peaks of PVE regions are difficult but the edge map of the image and 
  %  segmenation may allow the definition of this area ... 
  %  the matlab isosurface may also helps to extract a sublayer ...
    if ~exist('S2','var')
      S1 = isosurface( Ycls{1} + Ycls{2} , 128);
      S2 = isosurface( Ycls{2}           , 128);
    end
    mipveth(1) = cat_stat_kmeans( cat_surf_fun('isocolors',S1,Ymi ),1);
    mipveth(2) = cat_stat_kmeans( cat_surf_fun('isocolors',S2,Ymi ),1);
    mpveth(1)  = cat_stat_kmeans( cat_surf_fun('isocolors',S1,Ym  ),1);
    mpveth(2)  = cat_stat_kmeans( cat_surf_fun('isocolors',S2,Ym  ),1);
  end
  
  
  % final scaling
  if pve
    Tthmi = struct('CGW',mith,'PVEth',mipveth,'T3th',{[0 1/3:1/6:5/6 1:1/3:2]},'T3thx',{3*[0 sort( [mith mipveth]) 4/3:1/3:2] });
    Tthm  = struct('CGW',mth ,'PVEth',mpveth ,'T3th',{[0 1/3:1/6:5/6 1:1/3:2]},'T3thx',{3*[0 sort( [mth  mpveth])  4/3:1/3:2] });
  else
    Tthmi = struct('CGW',mith,'T3th',{0:1/3:2},'T3thx',{3*[0 sort( mith )  4/3:1/3:2] });
    Tthm  = struct('CGW',mth ,'T3th',{0:1/3:2},'T3thx',{3*[0 sort( mth  )  4/3:1/3:2] });
  end
  
  
  % the final scaling
  Ymi2  = cat_main_gintnormi(Ymi,Tthmi);
  Ym2   = cat_main_gintnormi(Ym ,Tthm);
  
  
  % display a result histogram figure
  if verb
    cg_hist2d(Ym2(Yb(:))); set(gcf,'Menubar','figure'); 
    hold on; for i=1:0.5:3, plot([i/3 i/3;],get(gca,'ylim'),'r'); end
    cg_hist2d(Ymi2(Yb(:))); set(gcf,'Menubar','figure'); 
    hold on; for i=1:0.5:3, plot([i/3 i/3;],get(gca,'ylim'),'r'); end
  end
return
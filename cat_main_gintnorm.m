function [Ym,T3th3,Tth,inv_weighting,noise] = cat_main_gintnorm(Ysrc,Ycls,Yb,vx_vol,res,extopts)
% This is a subfunction of cat_main.
% ______________________________________________________________________
% Global intensity normalization based on tissue thresholds estimated as 
% median intensity in the SPM tissue maps refined by edge (gradient) 
% information. Class propability should be higher than 50% (=128) to 
% avoid problems by the PVE or bias regions like basal ganglia or the CSF.
% Especialy correct CSF estimation can be problematic, because it is
% strongly influenced by the PVE and other tissues like blood vessels 
% and meninges. This structures with GM like intensity will cause a to 
% high global CSF value.
% For CSF, and WM we can use low gradient thesholds to avoid the PVE, but
% for GM this can lead to strong problems because to low thresholds will
% only give large GM areas like the basal ganlia, that have often a to high
% intensity. 
%
%   [Ym,Yb,T3th3,Tth,inv_weighting,noise] =
%     cat_main_gintnorm(Ysrc,Ycls,Yb,vx_vol,res,extopts)
%
%   Ym      .. intensity normalized image
%   Yb      .. brain mask
%   T3th3   .. [CSF,GM,WM] peak intensity
%   Tth     .. structure for inverse function cat_main_gintnormi
%   inv_weighting .. true in T2/PD images
%   noise   .. first guess of the noise level
%   cat_waring .. structure with warings
%
%   Ysrc    .. the original (noise/bias corrected) image
%   Ycls    .. SPM classification [GM,WM,CSF,HD1,HD2,BG]
%              (6 cells with uint8 classes images)
%   Yb      .. brain mask
%   vx_vol  .. voxel resolution of the images
%   res     .. SPM segmentation structure
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$
  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

  nwarnings = cat_io_addwarning; 
  if ~exist('extopts','var')
    extopts = cat_get_defaults('extopts');
  end
  if isstruct(Ycls)
    
    %% final peaks and intensity scaling
    %  -----------------------------------------------------------------
    %  RD202006: 
    %  This case is used for intensity normalized denoising.
    %  This case was created to run cat_main_gintnorm on low resolution and
    %  then update the high resolution case but I forgot that the full call
    %  of cat_main_gintnorm also includes a bias correction in some cases. 
    %  I.e. in the worst case of a high resolution scan the bias correction
    %  was not used and also the estimated peaks were biased by the missing
    %  correction. Hence, we skip the high resolution thing in case of bias
    %  correction (mostly inverse contrast).
    
    T3th  = Ycls.T3th;
    T3thx = Ycls.T3thx;

    % intensity scaling
    Ym    = Ysrc; 
    
    if all(T3th==T3thx), return; end

    isc   = 1;
    %T3th  = interp1(T3th,1:1/isc:numel(T3th)*isc,'spline');  %pchip');
    %T3thx = interp1(T3thx,1:1/isc:numel(T3th)*isc,'spline'); %pchip');

    for i=2:numel(T3th)
      M = Ysrc>T3th(i-1) & Ysrc<=T3th(i);
      Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3th(i-1))/diff(T3th(i-1:i))*diff(T3thx(i-1:i));
    end
    M  = Ysrc>=T3th(end); 
    Ym(M(:)) = numel(T3th)/isc/6 + (Ysrc(M(:)) - T3th(i))/diff(T3th(end-1:end))*diff(T3thx(i-1:i));    
    Ym = Ym / 3; 

    return
  end
    
  try
    if extopts.ignoreErrors > 2 
      error('cat_main_gintnorm:runbackup','Test backup function.');
    end
    
    %clsint  = @(x) cat_stat_nanmedian(Ysrc(Ycls{x} > 128)); 
    clsint  = @(x) round( sum(res.mn(res.lkp==x) .* res.mg(res.lkp==x)') * 10^5)/10^5; % better for ADNI ??? 
    clsints = @(x,y) [round( res.mn(res.lkp==x) * 10^5)/10^5; res.mg(res.lkp==x-((y==0)*8))']; 

    vxv    = 1/cat_stat_nanmean(vx_vol);
    res.mn = round(res.mn*10^5)/10^5; 

    % initial thresholds and intensity scaling
    T3th3 = [clsint(3) clsint(1) clsint(2)];
    BGth  = min(mean(Ysrc(Ycls{end}(:)>192)),clsint(end));
 
    inv_weighting = ~(T3th3(1)<T3th3(2) && T3th3(2)<T3th3(3));
   


    %% -------------------------------------------------------------------
    %  intensity checks and noise contrast ratio (contrast part 1)
    %  -------------------------------------------------------------------
    % relation between the GM/WM and CSF/GM and CSF/WM contrast has to be
    % greater than 3 times of the maximum contrast (max-min).
    clear Yn
    checkcontrast = @(T3th,minContrast) ...
      abs(diff(T3th([1,3]))) < (max(T3th(:))-min(T3th(:)))*minContrast || ...
      abs(diff(T3th(1:2)))   < (max(T3th(:))-min(T3th(:)))*minContrast || ...
      abs(diff(T3th(2:3)))   < (max(T3th(:))-min(T3th(:)))*minContrast;

    if checkcontrast(T3th3,1/9) && exist('cat_warnings','var') % contrast relation
      cat_io_addwarning([mfilename ':LowContrast'],...
        sprintf(['The contrast between different tissues is relative low! \n' ...
             '  (BG=%0.2f, CSF=%0.2f, GM=%0.2f, WM=%0.2f)\n'],BGth,T3th3),2,[1 0]);
% use backup pipeline?           
    end

    
% RD202006: Only process T1 data and call backup function otherwise. 
%           Try to use the neonate pipeline below (case 3) in the backup 
%           function to separate CSF from WM in the SPM segmentation in 
%           a future release (CAT12.8?) if you have more test and valiation 
%           data.
    if T3th3(1)>T3th3(3) || T3th3(2)>T3th3(3) || T3th3(1)>T3th3(2) % inverse (T2 / PD)
      % RD202006: WM < GM < CSF      
      % Use backup funtion in case of inverse contrast.
      error('cat_main_gintnorm:runbackup','Run PD/T2 processing.');

    elseif T3th3(1)<T3th3(2) && T3th3(2)<T3th3(3) % T1
      %%
      BGmin = min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))); 
      T3th3(1) = min( min(clsints(3,0)) , mean(Ysrc(Ycls{3}(:)>240))); 
      if isfield(res,'ppe') && isfield(res.ppe,'affreg') && isfield(res.ppe.affreg,'highBG') && res.ppe.affreg.highBG
        BGminl = BGmin - 8 * diff( [BGmin T3th3(1)] ); % compensate  BGminl*0.1+0.9*T3th3(1) the minimum is close to CSF here
        BGcon  = max([BGmin*1.1,T3th3(1) - cat_stat_nanmean(diff(T3th3))]);
      else
        %BGminl = max([BGmin*1.1,T3th3(1) - cat_stat_nanmean(diff(T3th3)),cat_stat_nanmedian(Ysrc(Ycls{end}(:)>128))]);
%       BGminl = BGmin - 8 * diff( [BGmin T3th3(1)] );
        BGminl = BGmin; 
        %BGcon  = max([BGmin*1.1,T3th3(1) - cat_stat_nanmean(diff(T3th3)),cat_stat_nanmedian(Ysrc(Ysrc(:)>BGminl & Ycls{end}(:)>128))]);
        BGcon = max([BGmin*1.1,T3th3(1) - cat_stat_nanmean(diff(T3th3)),median(Ysrc(Ycls{end}(:)>128))]);
      end
      %T3th3 = [max( min(res.mn(res.lkp==3 & res.mg'>0.3/sum(res.lkp==3)))*.05 + .95*max(res.mn(res.lkp==2 & res.mg'>0.3/sum(res.lkp==2))) , ...
      %              min(res.mn(res.lkp==3 & res.mg'>0.3/sum(res.lkp==3)))) ...
      %         max(res.mn(res.lkp==1 & res.mg'>0.1)) ...
      %         max(res.mn(res.lkp==2 & res.mg'>0.1))];
      T3th  = [BGmin ... minimum
               min(BGcon,BGminl*0.1+0.9*T3th3(1)) ... cat_stat_nanmean background (MT contrast with strong background noise)
               T3th3 ... csf gm wm 
               max(T3th3) + abs(diff(T3th3([1,numel(T3th3)])/2)) ... higher
               max(T3th3(end) + abs(diff(T3th3([1,numel(T3th3)])/2)) , ... maximum
                max(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:))))) ];
      T3thx = [0,0.05,1:5];
      Ysrcr = round(Ysrc*10^5)/10^5; 



    elseif T3th3(1)>T3th3(3) && T3th3(2)<T3th3(3) && T3th3(1)>T3th3(2)
      error('cat_main_gintnorm:runbackup','Run backup pipeline due to non T1 contrast or neonate dataset.');
%{   
% begin.neonatePipeline      
% #########################################################################
% RD202006: This part is the special neonate pipeline with problematic 
            thissue contrast and heavy processing.  It includes some useful
            ideas but need much further evaluation and tests. 
            It includes an extra bias correction and optimization of tissue
            classes when the WM intensity is quite similar to the CSF
            intensity.
      
    % This is a very special case of T2 weighting (of neonates) with
    % BG < GM < WM < CSF that need special correction for the sulcal 
    % CSF that mostly have WM like intensity due to the PVE. Hence, 
    % there are a lot of miss-labeled voxel that need further correction 
    % that follow in the second part of this file.
    %
    % Because, I observed a lot of problems by the SPM bias correction, I
    % increased the biasfwhm and the biasreg. Although, this did not work, 
    % the indroduce bias was smoother and could be corrected here very well. 

      cat_warnings = cat_io_addwarning(cat_warnings,...
        'cat_main:InverseContrast',...
        sprintf(['Inverse tissue contrast that requires strong modifications! \n' ...
             'In case of "BG<GM<WM<CSF", the CSF in sulci got WM-like intensities \n' ...
             'due to the PVE and require severe correction that may fail! \n' ...
             '  (BG=%0.2f, CSF=%0.2f, GM=%0.2f, WM=%0.2f)\n'],BGth,T3th3(1:3)),numel(cat_warnings)==0);

      Sth   = clsint(2);      
      Ym    = (Ysrc - BGth) / ( Sth - BGth); 
      Yp0   = single(Ycls{3})/255/3 + single(Ycls{1})/255*2/3 + single(Ycls{2})/255;
      Yg    = cat_vol_grad(Ym,vx_vol);
      Ydiv  = cat_vol_div(Ym,vx_vol);
      [Yp0r,resTh] = cat_vol_resize(Yp0,'reduceV',vx_vol,max(vx_vol)*3,16,'meanm'); 
      Yhd   = cat_vbdist(min(max(0,1-Yp0r*3),1)); 
      Yhd   = cat_vol_resize(Yhd,'dereduceV',resTh);   
      Yhd   = Yhd ./ max(Yhd(Yhd(:)>0)); 
      Ysrco = Ysrc+0;

      T3th2 = [clsint(3) ...
               min( [clsint(1) , ...
                     cat_stat_nanmean(Ysrco( Ysrco(:)>(BGth*0.8+0.2*Sth) & Ysrco(:)<cat_stat_nanmean(Sth) & ...
                      Yp0(:)>1.9/3 & Yp0(:)<2.2/3 & Ydiv(:)<0.05 & Yg(:)>0.05 & Yg(:)<0.15))  ]) ... 
               max( [clsint(2) cat_stat_nanmedian( Ysrco(Ycls{2}(:)>192) ) ]) ...
              ];


      %%      
      if ~exist('Yy','var') && isfield(res,'Twarp'), Yy = res.Twarp; end
      LAB = extopts.LAB;
      if ~exist('Yy','var')
        PA  = extopts.cat12atlas;
        VA  = spm_vol(PA{1});
        YA  = cat_vol_ctype(spm_sample_vol(VA,double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0));
        YA  = reshape(YA,size(Ym));
        YA(mod(YA,2)==0 & YA>0)=YA(mod(YA,2)==0 & YA>0)-1;  
      else
        YA = ones(size(Yg)); 
      end
      if ~debug, clear Yy; end


      %% bias correction
      Ym   = (Ysrco - BGth) / ( Sth - BGth); 
      Ycbg = Yg<0.1 & (Ym<1 | Ycls{1}>Ycls{2} & Ydiv>=0) & YA==LAB.CB; 
      Ycbw = Yg<0.1 & Ycls{2}>Ycls{1} & Ydiv<0 & YA==LAB.CB & ~Ycbg; 
      Ybg  = cat_vol_morph(cat_vol_morph(cat_vol_smooth3X(Ym<1.02 & Yp0>2/3 & Yp0<1 & Yg<cat_stat_nanmean(Yg(:))/2,1.2)>0.8 & ...
        Yhd>0.5,'o',1),'d',4) & Yp0>1.5/3 & Ym<1.02 & Yg<cat_stat_nanmean(Yg(:)) & cat_vol_morph(YA==LAB.BG | YA==LAB.TH,'d',2); % subcortical structures
      Ygw  = (Yp0>1.9/3 | (Yp0>0 & Ym>0.7 & Ym<1.4)) & cat_vol_morph(cat_vol_morph(smooth3(Ycls{1}>240 & Yg<0.1 & abs(Ydiv)<0.01)<0.3,'c',1),'e',3) & ...
        Ym<(clsint(3)/clsint(2)*0.8 + 0.2) & Ym>(clsint(1)/clsint(2)*0.6) & Yg<0.4 & ~Ybg & ~Ycbg; 
      Ygw = cat_vol_morph(Ygw,'l',[inf 0.05])>0;
      Ygw  = Ygw | (Yg<0.1 & Ycls{2}>Ycls{1} & Ycls{2}>Ycls{3} & (Yp0>1.5/3 | Ym>0.95 & Ym<1.5) & ~Ybg & ~Ycbg) | YA==LAB.BS | Ycbw;
      Ygw = cat_vol_morph(Ygw,'l',[inf 0.1])>0;
      Ygw  = Ym .* Ygw;
      Ygw  = cat_vol_median3(Ygw,Yp0>0);
      Ygm  = Ym .* (cat_vol_morph(smooth3(~Ygw & Yg<0.10 & abs(Ydiv)<0.05 & Ycls{1}>Ycls{2} & Ycls{1}>Ycls{3} & Yhd<0.4)>0.5,'o',2) | Ycbg )/ ...
              (T3th2(2) / T3th2(3)); 
      %Ycm  = Ym .* cat_vol_morph(smooth3(~Ygw & ~Ygm & Yg<0.1 & (Ycls{3}>Ycls{2} & Ycls{3}>Ycls{1})>0.5) | (Ym>1.2 & Yb)); 
      %Ycm  = cat_vol_localstat(Ycm,Ycm>0,1,3);
      %Ycm  = Ycm / mean(Ycm(Ycm(:)>0)); %  * (T3th2(1) / T3th2(3)); 
      Ygw  = Ygw + (Ygm .* (Ygw==0)) + ((Ym .* Ybg .* (Ygw==0)) / (cat_stat_nanmean(T3th2(2)) / T3th2(3)) ); % + (Ycm .* (Ygm==0)); 
      Ygw2 = Ym .* (Yg<max(0.1,min(0.2,cat_stat_nanmean(Yg(:)))) & ( ((Yp0>2.9/3 | Ym>0.9) & Ym<1.5 & Ygw>0) | Ycbw | ...
                    (abs(Ydiv)<0.1 & Ycls{2}/2>Ycls{1} & Ycls{2}/2>Ycls{3}) ) & ~Ybg & ~Ycbg);
      %% field approximation
      [Ywi,Ywi2,resT2] = cat_vol_resize({Ygw,Ygw2},'reduceV',vx_vol,max(vx_vol)*2,16,'max'); 
      for i=1:1, Ywi2 = cat_vol_localstat(Ywi2,Ywi2>0,1,3); end
      for i=1:1, Ywi = cat_vol_localstat(Ywi,Ywi>0,2,3); end % only one iteration!
      for i=1:1, Ywi = cat_vol_localstat(Ywi,Ywi>0,1,1); end
      Ywi  = cat_vol_approx(Ywi,'nn',resT2.vx_volr,2);
      Ywi  = cat_vol_smooth3X(Ywi,1); Ywi(Ywi2>0)=Ywi2(Ywi2>0);
      Ywi  = cat_vol_smooth3X(Ywi,1); % highres data have may stronger inhomogeneities 
      Ywi  = cat_vol_resize(Ywi,'dereduceV',resT2);    
      Ywi  = Ywi / cat_stat_nanmean(Ywi(Ygw2>0)); 
      Ysrc  = Ysrco ./ Ywi * ( cat_stat_nanmean(Ysrco(Ygw2>0)/T3th2(3)) / cat_stat_nanmean(Ywi(Ygw2>0)) );



      %% first initial scaling for gradients and divergence
      %T3th3 = [cat_stat_nanmean( res.mn(res.lkp==3 & res.mg'>0.1)) ...
      %         min( [res.mn(res.lkp==1 & res.mg'>0.1) , ...
      %               cat_stat_nanmean(Ysrc( Ysrc(:)>(BGth*0.8+0.2*Sth) & Ysrc(:)<cat_stat_nanmean(Sth) & ...
      %                Yp0(:)>1.9/3 & Yp0(:)<2.2/3 & Ydiv(:)<0.05 & Yg(:)>0.05 & Yg(:)<0.15))  ]) ... 
      %         max( [sum(res.mn(res.lkp==2) .* res.mg(res.lkp==2)') median( Ysrc(Ycls{2}(:)>192) ) ]) ...
      %        ];
      clear T3th T3thx; 
      T3th = [ min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) ...
               min( T3th3(2)*0.5+0.5*min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) , BGth ) ...
               ...
               min(Ysrc( Ysrc(:)>(BGth*0.8+0.2*Sth) & Ysrc(:)<cat_stat_nanmean(Sth) & ...
                      Yp0(:)>1.9/3 & Yp0(:)<2.2/3 & Ydiv(:)<0.05 & Yg(:)>0.05 & Yg(:)<0.15)) ... head
               T3th3 ...
               ...
               (T3th3(3)*0.8 + 0.2*max(res.mn(res.lkp==3 & res.mn>max(res.mn(res.lkp==2)))) ) ...
               (T3th3(3)*0.5 + 0.5*max(res.mn(res.lkp==3 & res.mn>max(res.mn(res.lkp==2)))) ) ...
               (T3th3(3)*0.2 + 0.8*max(res.mn(res.lkp==3 & res.mn>max(res.mn(res.lkp==2)))) ) ...
               ...
               max(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:))))]; ...

      T3thx = [0,0.05,1.75, 1,2,3,  3.1,2.0,1.1,  1.0]; 

      [T3th,si] = sort(T3th);
      T3thx     = T3thx(si);

      Ysrcr = round(Ysrc*10^5)/10^5; 

      inv_weighting = 1;

      if debug
          Ym = Ysrcr+0;
          for i=2:numel(T3th)
            M = Ysrc>T3th(i-1) & Ysrc<=T3th(i);
            Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3th(i-1))/diff(T3th(i-1:i))*diff(T3thx(i-1:i));
          end
          M  = Ysrc>=T3th(end); 
          Ym(M(:)) = numel(T3th)/6 + (Ysrc(M(:)) - T3th(i))/diff(T3th(end-1:end))*diff(T3thx(i-1:i));    
          Ym = Ym / 3; 
          ds('l2','',vx_vol,Ysrc/T3th3(3), round(Ym*3),Ysrc/T3th3(3),Ym,60)
      end
% #########################################################################
% end.neonatePipeline      
%}

    else
      error('cat_main_gintnorm:runbackup','Test backup function due to unknow contrast.');
    end

    
% RD202006: The following part is only for the standard T1 case.

    %% intensity scaling for gradient estimation
    Tth.T3th  = T3th;
    Tth.T3thx = T3thx;

    Ym = Ysrcr+0; 
    for i=2:numel(T3th)
      M = Ysrcr>T3th(i-1) & Ysrcr<=T3th(i);
      Ym(M(:)) = T3thx(i-1) + (Ysrcr(M(:)) - T3th(i-1))/diff(T3th(i-1:i))*diff(T3thx(i-1:i));
    end
    M  = Ysrcr>=T3th(end); 
    Ym(M(:)) = numel(T3th)/6 + (Ysrcr(M(:)) - T3th(i))/diff(T3th(end-1:end))*diff(T3thx(i-1:i));    
    Ym = Ym / 3; 



    %% new initial segment threshold
    if ~exist('Yg','var'), Yg  = cat_vol_grad(Ym,vx_vol)./max(eps,Ym); end
    T3th  = [cat_stat_nanmedian(Ysrcr(Ycls{3}(:)>192 & Yg(:)<0.20 & Ym(:)<0.45)) ...
             cat_stat_nanmedian(Ysrcr(Ycls{1}(:)>192 & Yg(:)<0.20)) ...
             cat_stat_nanmedian(Ysrcr(Ycls{2}(:)>192 & Yg(:)<0.10))];
    Ynw   = cat_vol_localstat(Ysrc,Ycls{3}>192,2,4);
    Ync   = cat_vol_localstat(Ysrc,Ycls{2}>192,2,4); 
    noise = round(min(cat_stat_nanmean(Ynw(Ynw(:)>0)),cat_stat_nanmean(Ync(Ync(:)>0))) / min(abs(diff(T3th(1:3)))) * 10^6)/10^6; 
    clear Ynw Ync;

    if debug==2
      [mrifolder, reportfolder] = cat_io_subfolders(res.image0(1).fname,struct('extopts',extopts));
      [pth,nam] = spm_fileparts(res.image0(1).fname);
      tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',1,'gintnorm00'));
      save(tmpmat,'Ysrc','Ycls','Yb','vx_vol','res','T3th','T3thx','Yg','Ym','noise');
    end






    %% -------------------------------------------------------------------
    %  check modality (contrast part 2)
    %  -------------------------------------------------------------------
    %  It is possible to invert T2 and PD images based on the SPM class 
    %  information, but actual there is no time to develope and proof this 
    %  function in detail, due to the most other functions ...
    %  -------------------------------------------------------------------
    if T3th(1)<T3th(2) && T3th(2)<T3th(3)
    %  -------------------------------------------------------------------
    %  standard T1 contrast
    %  -------------------------------------------------------------------
    %  For T1 data SPM mean tissue values were not always correct. 
    %  Especially, high and low contrast images or images with incomplete
    %  inhomogeneity correction can have bad peaks (ADHD200/..NYC..14). 
    %  So it is better to use the SPM segments and add some further 
    %  knowledge (gradient & divergence) to refine these segments and 
    %  estimate the median value of the segment that is typcialy more 
    %  stable than the mean value. 
    %  -------------------------------------------------------------------

      % check SPM segmentation
      if exist('cat_warnings','var')
        Ymx = single(Ycls{1})/255*2/3 + single(Ycls{2})/255+ single(Ycls{3})/255*1/3;  
        Ygw = Yb & ((single(Ycls{1})+single(Ycls{2}))>128);
        Ymp0diff = sqrt(cat_stat_nanmean(Ym(Ygw(:)) - Ymx(Ygw(:)))^2); 
        if Ymp0diff>0.10 && debug
          cat_io_addwarning([mfilename ':badSPMsegment'],sprintf(...
            ['SPM segmentation does not fit to the image (RMS(Ym,Yp0)=%0.2f).\n'...
             'This can be an alignment problem (check origin), ' ...
             'untypical subjects (neonates, non-human),\n'...
             'bad image contrast (C=%0.2f,G=%0.2f,W=%0.2f), \n'...
             'low image quality (NCR~%0.2f), or something else ...'],...
             Ymp0diff,T3th,noise),2,[1 0]); 
        end
        clear Ymx;
      end


      %% skull-stripping warning
% #########      
% RD202006: Skull-Stripping warning in backup case?  
%           Move to updateSPM ???
% #########      
      if numel(Ycls)>4
        skulltest = (cat_stat_nanmedian(Ysrc(Ycls{5}(:)>192 & Ysrc(:)>T3th(2))) < ... 
           cat_stat_nanmedian(Ysrc(Ycls{3}(:)>192 & Ysrc(:)>0))); 
        if exist('cat_warnings','var') &&  (isnan(skulltest) || skulltest)

          % Skull-Stripped images can of course lead to problems with to strong
          % brain masks, but the bigger problem here is that the CSF intensity 
          % threshold were maybe affected. 

          % If a skull-stripping was used, we will use this as initial mask 
          % that we close and dilate a little bit. 
          % Now, the original image can be corrected in the stripped area, 
          % because some images have missing points (slicewise). Becuase of 
          % the gaussian functions a hard boundary is better.
          if Ymp0diff<0.05 && numel(Ysrc>0)/numel(Ysrc)<0.8
            Yb    = smooth3(cat_vol_morph(cat_vol_morph(Ysrc>0,'lc',3),'d'))>0.5;
            CSFth = min([cat_stat_nanmedian(Ysrc(Ycls{3}(:)>240 & Ysrc(:)>0)), ... 
                         cat_stat_nanmedian(Ysrc(Ycls{3}(:)>192 & Ysrc(:)>0)), ... 
                         cat_stat_nanmedian(Ysrc(Ycls{3}(:)>128 & Ysrc(:)>0)), ...
                         cat_stat_nanmean(Ysrc(Ysrc>0))*0.5])*0.9; % 
            Ysrc  = cat_vol_laplace3R(max(CSFth,Ysrc),Yb & Ysrc==0,0.2) .* Yb;
            cat_io_addwarning([mfilename ':SkullStripped'],...
               'Skull-stripped input image detected! Try boundary cleanup.',1,[1 0]);  
          else
            cat_io_addwarning([mfilename ':SkullStripped'],...
               'Skull-stripped input image?',1,[1 0]); 
          end
        end
      end


      %% segment refinement and median peak estimation 
      %  -----------------------------------------------------------------
      Yg    = cat_vol_grad(Ym,vx_vol);
      Ydiv  = cat_vol_div(Ym,vx_vol);
      %noise = estimateNoiseLevel(Ym,Ycls{2}>192); 

      Yb2   = cat_vol_morph(Yb & Ym>0.5,'e',2*vxv); 
      gth   = max(0.06,min(0.3,noise*6));
      %Ybm   = cat_vol_morph(Ycls{6}>240 & Ysrc<min(T3th),'lc'); 
      BGmin = min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))); 
      BGcon = max([BGmin*1.1,T3th3(1) - cat_stat_nanmean(diff(T3th3)),cat_stat_nanmedian(Ysrc(Ycls{end}(:)>128))]);
      if isfield(res,'ppe') && isfield(res.ppe,'affreg') && isfield(res.ppe.affreg,'highBG') && ~res.ppe.affreg.highBG
        BGcon = max([BGmin,T3th3(1) - cat_stat_nanmean(abs(diff(T3th3)))]);
      end
      %BMth  = max(BGmin,mean([BGcon,T3th(1)])); % - diff(T3th(1:2)))); %max(0.01,cat_stat_nanmedian(Ysrc(Ybm(:))));
      BMth  = max(BGmin,min(BGcon,T3th(1) - diff(T3th(1:2)))); %max(0.01,cat_stat_nanmedian(Ysrc(Ybm(:))));
      Ywm   = (Ycls{2}>128  & Yg<gth) | ((Ym-Ydiv*2)>(1-0.05*cat_stat_nanmean(vx_vol)) & Yb2); % intensity | structure (neonate contast problem)
      if isfield(res,'ppe') && isfield(res.ppe,'affreg') && isfield(res.ppe.affreg,'highBG') && res.ppe.affreg.highBG
        Ycm   = smooth3((Ycls{3}>240 | Ym<0.4) & Yg<gth*3 & Yb & ~Ywm & Ycls{1}<8 & Ym<0.7)>0.5; % important to avoid PVE!
      else
        Ycm   = smooth3((Ycls{3}>240 | Ym<0.4) & Yg<gth*3 & Yb & ~Ywm & Ycls{1}<8 & Ysrc>BMth & Ym<0.7)>0.5; % important to avoid PVE!
      end
      
      % If SPM get totaly wrong maps due to bad image orientations our 
      % segment were incorrect too (or empty) and peak estimation fail.
      % I try to use the kmeans, but in WM it is affected by WMHs, in 
      % CSF by blood vessels and meninges and in GM noise and subcortical
      % structures were problematic. In ADHD/..NYC..14 the basal structes 
      % get the average peak and the cortex was detected as CSF. There 
      % were much more images with smaller problems ...
      Ysrcr  = round( Ysrc.*10^5 ) / 10^5;
      WMth   = cat_stat_nanmedian(Ysrcr(Ywm(:))); % cat_stat_kmeans(Ysrc(Ycls{2}(:)>192 & Yg(:)<gth),1); % GM/WM WM  
      CSFth  = cat_stat_nanmedian(Ysrcr(Ycm(:))); % cat_stat_kmeans(Ysrc(Ycls{3}(:)>64 & Yg(:)>gth & Yb(:)),2); % CSF CSF/GM
        %  0.05 <<<<< BMth + 4*cat_stat_nanstd(Ysrc(Ybm(:)))
      if isfield(res,'ppe') && isfield(res.ppe,'affreg') && isfield(res.ppe.affreg,'highBG') && res.ppe.affreg.highBG
        Ybg  = cat_vol_morph(Ycls{6}>128 & ~Yb,'o',2); 
      else
        Ybg  = cat_vol_morph(Yg<0.10 & Yb & Ysrc<WMth*(1-0.03*cat_stat_nanmean(vx_vol)) & Ysrc>CSFth*1.5 & Ycls{3}<64,'o',2);
      end
      Ygm    = ~Ybg & Yg<0.4 & Ysrc<min(clsint(2)*0.8+clsint(1)*0.2,WMth+0.5*diff([CSFth,WMth])) & Yg<gth*2 & Ycls{1}>32 & ~Ywm & Ycls{2}<64 & ...
        Ysrc>(CSFth+0.1*diff([CSFth,WMth])) & ~Ywm & ~Ycm & Yb & abs(Ydiv)<0.2; 
      %Ygm   = Ygm | (Ycls{1}>64 & Ybg & ~Ywm);
      GMth   = cat_stat_nanmedian(Ysrcr(Ygm(:))); %cat_stat_kmeans(Ysrc(Ygm(:)),3); % CSF/GM GM GM/WM
      T3th_cls  = round([CSFth(1) GMth(1) WMth(1)]*10^4)/10^4;
      %clear Ybg
     %
      if any(isnan(T3th_cls)) 
        fprintf('\n');
        error('cat_main:cat_main_gintnorm:nobrain',...
          'Bad SPM-Segmentation. %s!\n',...
          spm_file('Please check image orientation and quality','link',['spm_image(''Display'', ''' res.image0.fname ''')']) );
      end
      % median tissue peaks


      if debug==2
        tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',1,'gintnorm01'));
        save(tmpmat,'Ysrc','Ycls','Yb','vx_vol','res','T3th','Yg','Ydiv','Ym',...
         'Yb2','gth','Ybm','BMth','Ywm','Ygm','Ycm','Ybg','T3th_cls','T3th','noise');
      end


      %% final peaks and intensity scaling
      %  -----------------------------------------------------------------
      T3th3 = T3th_cls;
      if isfield(res,'ppe') && isfield(res.ppe,'affreg') && isfield(res.ppe.affreg,'highBG') && res.ppe.affreg.highBG
        BMth    = min(Ysrc(Yb(:)))*0.99 + 0.01*T3th3(1);
        BMCSFth = BMth*0.5 + 0.5*(T3th3(1)); % - mean(abs(diff(T3th(:)))));
        T3thx   = [2/5, 3/5, 4/5, 1.1, 2.2, 3:5];
      else 
        BMCSFth = min(BGth,mean([BMth,T3th3(1)]));
        T3thx   = [0,0.02,0.05,1:5];
      end
      T3th  = [min(Ysrcr(~isnan(Ysrcr(:)) & ~isinf(Ysrcr(:)))) BMth BMCSFth T3th3 ...
                T3th3(end) + diff(T3th3([1,numel(T3th3)])/2) ... WM+
                max(T3th3(end)+diff(T3th3([1,numel(T3th3)])/2) , ... max
                max(Ysrcr(~isnan(Ysrcr(:)) & ~isinf(Ysrcr(:))))) ];
      


      % intensity scaling
      Ym    = Ysrc; 
      isc   = 1;
      %T3th  = interp1(T3th,1:1/isc:numel(T3th)*isc,'spline');  %pchip');
      %T3thx = interp1(T3thx,1:1/isc:numel(T3th)*isc,'spline'); %pchip');

      for i=2:numel(T3th)
        M = Ysrc>T3th(i-1) & Ysrc<=T3th(i);
        Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3th(i-1))/diff(T3th(i-1:i))*diff(T3thx(i-1:i));
      end
      M  = Ysrc>=T3th(end); 
      Ym(M(:)) = numel(T3th)/isc/6 + (Ysrc(M(:)) - T3th(i))/diff(T3th(end-1:end))*diff(T3thx(i-1:i));    
      Ym = Ym / 3; 

      Tth.T3th  = T3th;
      Tth.T3thx = T3thx;
      
      
    elseif T3th3(1)>T3th3(3) && T3th3(2)>T3th3(3)
       error('runbackup')
%{       
% RD202006:  GM < WM < CSF - part 2 of the neonate contrast?
% Do not remove this special contrast case to early.
      
      %% reestimation of brain mask
      Yb  = Ym>0.8 & Ym<1.2 & (Ycls{5}<64); Yb  = single(cat_vol_morph(Yb,'lo',1));
      [Ybr,Ymr,Ycls5,resT2] = cat_vol_resize({single(Yb),Ym,single(Ycls{5})/255},'reduceV',vx_vol,2,32); 
      Ybr(~Ybr & (Ymr<2.5/3 | Ymr>3.2/3 | Ycls5>0.5))=nan; 
      [Ybr1,Ydr] = cat_vol_downcut(Ybr,Ymr,0.03); Ybr(Ydr<100)=1; Ybr(isnan(Ybr))=0;
      Ybr(~Ybr & (Ymr<1.9/3 | Ymr>3.2/3 | Ycls5>0.5))=nan; 
      [Ybr1,Ydr] = cat_vol_downcut(Ybr,Ymr,0.01); Ybr(Ydr<100)=1; Ybr(isnan(Ybr))=0;
      Ybr(~Ybr & (Ymr<1/3 | Ymr>2.5/3 | Ycls5>0.5))=nan; 
      [Ybr1,Ydr] = cat_vol_downcut(Ybr,Ymr,-0.01); Ybr(Ydr<100)=1; Ybr(isnan(Ybr))=0;
      Ybr = Ybr>0 | (Ymr<0.8 & cat_vol_morph(Ybr,'lc',6) & Ycls5<0.02); % large ventricle closing
      Ybr = cat_vol_morph(Ybr,'lc',2);                 % standard closing
      Yb  = cat_vol_resize(cat_vol_smooth3X(Ybr,2),'dereduceV',resT2)>0.4; 
      clear Ybr Ymr;
      %% filtering
      %YM  = Ysrc<Tth.T3th(5)/1.2; 
      %Ym(YM(:)) = Ysrc(YM(:)) / (Tth.T3th(5)/1.2);    
      YM  = (smooth3(Ysrc<Tth.T3th(5)/1.2) & smooth3(Ysrc>Tth.T3th(4))) | Ym>2; 
      Ym = cat_vol_median3(Ym,YM,Ym<1.5,0.1); 
      cat_sanlm(Ym,1,3)
%}
       
    elseif T3th3(1)>T3th3(3) && T3th3(2)<T3th3(3) && T3th3(1)>T3th3(2)
       error('runbackup')
    end
  catch e 
    if extopts.ignoreErrors < 1
      rethrow(e)
    else
      % [Ym,T3th3,Tth,inv_weighting,noise,cat_warnings] = cat_main_gintnorm(Ysrc,Ycls,Yb,vx_vol,res,Yy,extopts)
      if extopts.ignoreErrors > 2
        cat_io_addwarning([mfilename ':runGeneral'], ...
          'IgnoreErrors: Run generalized function.',3,[1 0]);
      elseif isfield(extopts,'inv_weighting') && extopts.inv_weighting==0 % if not already active
        if exist('BGth','var') && exist('T3th3','var') 
          cat_io_addwarning([mfilename ':InverseContrast'],...
            sprintf(['No T1 contrast detected. Run generalized function. \\\\n' ...
               '  (BG=%0.2f, CSF=%0.2f, GM=%0.2f, WM=%0.2f)'],BGth,T3th3(1:3)),1,[1 0]);
        else
          cat_io_addwarning([mfilename ':InverseContrast'],...
            sprintf('  No T1 contrast detected. Run generalized function.'),1,[1 0]);
        end
      end
      
      % if this is a untypical error then print it as warning
      if ~strcmp( e.identifier , [mfilename ':runbackup'] ),  warning(e.message); end
      
      Ysrc    = cat_stat_histth(Ysrc,extopts.histth); % remove outlier
      clsints = @(x,y) [round( res.mn(res.lkp==x) * 10^5)/10^5; res.mg(res.lkp==x-((y==0)*8))']; % SPM peak definition  
      clsint  = @(x) cat_stat_kmeans(Ysrc(cat_vol_morph(Ycls{x}>128,'e') | Ycls{x}>250));                                      % median within the tissue label 
      clsintv = [clsint(1) clsint(2) clsint(3) clsint(4) clsint(5) clsint(6)];
      for i=1:6, if isnan(clsintv(i)), clsintv(i) = sum(prod(clsints(i,1))); end; end % use median if possible
      
      bc = 1; 
      if bc
        % This is a quite simple bias corrections that should remove small 
        % local reminding inhomogeneities. It must consider WMHs to avoid
        % wrong corrections and get the correct WM peak. 

% #########       
% RD202006: Use bias correction in general?
%           + it sesms to work also with standard data
%           - in standard T1 data LAS already includes a BC but this may
%             also benefits by corrections here?
%           - bias correction will alter the intensity values and the 
%             spm estimation is not longer correct 
%           > this is super in this backup case (that do not support LAS 
%             inclusive the bias correction there)
% #########       
        Yg      = cat_vol_grad(Ysrc,vx_vol) ./ ... 
                min( max([clsintv(1) clsintv(2) clsintv(3)]) , Ysrc); % normalize extrem values
        gth     = max(0.10,min( 1/3 , cat_stat_nanmean( Yg( cat_vol_morph( smooth3( Ycls{2}>240 )>0.5 ,'e') )) ));  
        clsintg = @(x) cat_stat_nanmedian(Ysrc(Ycls{x}>128 & Yg<gth*2));   
        
              
        % remove strong outlier
        Ysrcs = cat_vol_median3(Ysrc,Yb,Yb,0.1); 
        
        % detect WMHs
        Ygmp = max(0,Yb - max(0,abs( Ysrc - clsintg(1))) / abs(diff([clsintg(2),clsintg(1)])));    % GM int prob map
        Ywmh = (Ygmp - single(Ycls{1})/255); Ywmh(smooth3(Ywmh)<0.5) = 0;                          % larger WMHs 
        Ywmh = (Ygmp - single(Ycls{1})/255) .* cat_vol_morph(Ywmh>0.2,'d',2,vx_vol);               % also smaller WMHs around the large ones
        Ywmh(smooth3(Ywmh)<0.3) = 0;                                                               % remove small dots
        clsintg2  = @(x) cat_stat_nanmedian(Ysrcs(Ycls{x}>128 & Yg<gth*2 & ~Ywmh));                % awoid WMH areas
        for i=1:6
          clsintv(i) = clsintg2(i);  
          if isnan(clsintv(i)), clsintv(i) = clsint(i);  end
          if isnan(clsintv(i)), clsintv(i) = sum(prod(clsints(i,1))); end
        end
          
        
        % WM region
        Ybe  = cat_vol_morph( (Ycls{2} + Ycls{1})>8 | cat_vol_morph(Yb,'e',5,vx_vol) ,'de',2,vx_vol);  % we have to be carefull with WM close to the skull
        Ywm  = Ycls{2}>8 & smooth3(Ycls{3})<8 & ~cat_vol_morph(Ycls{1}>32,'o',2,vx_vol) & ~Ywmh & Ybe; % WM with close GM/WM PVE voxel
        Ywm(smooth3(Ywm)<0.3) = 0;                                                                     % remove small dots
        Ygwm = Ycls{1}>8 & ~Ywmh & ~Ywm & Ycls{2}>Ycls{3} & ~cat_vol_morph(Ycls{1}>128,'o',2,vx_vol) & Ybe;  % WM with more distant GM/WM PVE voxel
        Ygwm = Ygwm | (cat_vol_morph(Ywm,'d') & Ybe);                                                  % add neigbors to get a more stable area (use large min-max distance)
        %
        Ywmi = cat_vol_localstat(Ysrcs,Ywm,2,2 + (clsintv(2)>clsintv(1))) + ...                        % filter size 1 (1vx-more aggressive,2vx-more-save)
               Ygwm .* cat_vol_localstat(Ysrcs,Ygwm,3,2 + (clsintv(2)>clsintv(1)));                    % filter size 2 (2vx-more aggressive,3vx-more-save)
        Ywmm = Ywm | Ygwm; Ywmi(Ywmm) = Ywmi(Ywmm) ./ (Ywm(Ywmm) + Ygwm(Ywmm)); Ywm = Ywmm;            % combine both areas
        clear Ygwm Ywmm Ybe;                                                                           % cleanup   
        
        % CSF region - for CSF we cannot expect to have large CSF in the neighborhood  
        %  CSF is not a good region in some images ...
        %    DB/ADNI_HC_100_0015_T1_sd000105-rs00_S8833_I33046
        con  = mean(abs(diff([clsintv(1),clsintv(2),clsintv(3),clsintv(1)]))); 
        Ycm  = Ycls{3}>192 & Ycls{2}<8 & smooth3(Ycls{1})<8 & ... Yg<min(0.5,gth*4) & ... Yg is instable
                Ysrc>( clsintv(3) - con/4 ) & Ysrc<( clsintv(3) + con/2 );                      % CSF with PVE
        Ycsd = cat_vol_localstat(Ysrc,Ycm,1,4);
        Ycm(Ycsd/clsintv(3)>0.2) = 0;  
        Ycm(smooth3(Ycm)<0.9) = 0; 
        % remove dots
        Ycmi = cat_vol_localstat(Ysrc,Ycm,3,1); %3 - (clsintg(1)>clsintg(3))); % minimum does not work in all cases
        if sum(Ycm) < 1000
          Ycm  = Ycls{3}>192; 
          Ycmi = Ycm .* clsintv(3); 
        end
        %% GM without PVE
        Ygm = (cat_vol_morph(Ycls{1}>128,'e',1) | Ycls{1}>240) & Yg<gth * 2; 
        Ygm = cat_vol_morph(Ygm,'l',[0.01 100 ]); Ygm(smooth3(Ygm)<0.5) = 0; 
        Ygmi = cat_vol_localstat(Ysrc,Ygm,2,1);
        Ygmi = cat_vol_localstat(Ygmi,Ygm,1,1);

% ########        
% RD202006: Use atlas data to improve/save subcortical areas > CAT12.?
%           Use divergence to find gyri/sulci > CAT12.8
%           Add skull-tissue bias correction  > CAT12.8
% ########        
        
        
        % get background area
% ########        
% RD202006: Add skull-stripped special case! > CAT12.8
% ########        
        if isfield(res,'bge') 
          if all(size(res.bge) == size(Ysrc))
            Ybge = res.bge; 
          else
            % low resolution update
            Ybge = cat_vol_resize(single(res.bge),'reduceV', vx_vol .* size(Ysrc)./size(res.bge), min(vx_vol*2, 1.4), 32, 'meanm')>0.5;
          end
        else
          Ybge = cat_vol_morph(Ycls{6}>128,'de',15,vx_vol); 
        end
        if isfield(res,'ppe') && isfield(res.ppe,'affreg') && isfield(res.ppe.affreg,'highBG') && res.ppe.affreg.highBG
          Ybge = cat_vol_localstat(Ysrc,Ybge,2,1) / cat_stat_nanmean(Ysrc(Ybge(:)>0));
        end
        %% final correction map with correction for double counts
        Yi   =  Ywmi/clsintv(2) + Ygmi/clsintv(1) + Ycmi/cat_stat_nanmean([clsintv(3),cat_stat_nanmean(Ycmi(Ycmi(:)>0))]) + single(Ybge); % and also the eroded background  
        if ~debug, clear Ygmi Ywmi Ycmi Ybge; end 
        Yim  =  Ycm &  Ygm &  Ywm; Yi(Yim) = Yi(Yim) ./ ( Ywm(Yim) + Ycm(Yim) + Ywm(Yim));
        Yim  =  Ycm &  Ygm & ~Ywm; Yi(Yim) = Yi(Yim) ./ ( Ycm(Yim) + Ygm(Yim));
        Yim  =  Ycm & ~Ygm &  Ywm; Yi(Yim) = Yi(Yim) ./ ( Ycm(Yim) + Ywm(Yim));
        Yim  = ~Ycm &  Ygm &  Ywm; Yi(Yim) = Yi(Yim) ./ ( Ygm(Yim) + Ywm(Yim));
        Yi(isnan(Yi) | isinf(Yi))=0;
        if ~debug, clear Yim Ycm Ygm Ywm; end 
        %{ 
        % RD202401: incorrect setup > try to simplify
        if isfield(res,'ppe') && isfield(res.ppe,'affreg') && isfield(res.ppe.affreg,'highBG') && res.ppe.affreg.highBG
          Yiw   = cat_vol_approx(Yi,'nn',4 * (res.ppe.affreg.highBG * 3 + 1)); 
        else
          Yiw   = cat_vol_approx(Yi,'nn',4); 
        end
        %}
        Yiw   = cat_vol_approx(Yi,'rec'); 
        %% correction 
        %  ds('d2sm','a',1,Ysrc/clsint(3),Ysrc ./ Yi / clsint(3),140)
        Ysrc = Ysrc ./ max(eps,Yiw); 
      end
      
      
      %% T3th3
      if ~exist('T3th3','var') || any(isnan(T3th3))
        T3th3 = [clsintv(3) clsintv(1) clsintv(2)];
      end
      
      % T3th, T3thx - keep it simple and avoid inversion
      if isfield(res,'ppe') && isfield(res.ppe,'affreg') && isfield(res.ppe.affreg,'highBG') && res.ppe.affreg.highBG
        BGth = round( clsint(6) / max([clsintv(3) clsintv(1) clsintv(2)]) * 6 ) / 2;
      else
        BGth = 0.1; 
      end
%% ##################
%  RD20200428:
%  Peak estimation needs further improvement (for the PD/T2 processing)
%  Noisy data and PVE lead no biased WM and CSF peaks that are to close to 
%  the GM values and cause GM volume/thickness underestimation.
%  ##################
      minYsrc = min(Ysrc(:)); 
      if 1%~exist('T3th','var') || ~exist('T3thx','var') || any(isnan(T3th3))
        if 0
          T3th  = [min(Ysrc(:)) sort( [ clsint(6) clsintv(3) clsintv(1) clsintv(2) clsint(5) ] ) max(Ysrc(:)) ];
        else
          order = [1 2 3 4 6];
          if extopts.gcutstr == -2
            clsn  = [1 1 1 1 1];
          else
            clsn  = [5 2 5 1 1]; % RD202108: only useful for T1?
          end
          T3th4 = zeros(1,5);
          for i = 1:5
            % try to use lower resolution to reduce noise
            Yclsr = Ycls{ order(i) }>255 | cat_vol_morph( Ycls{ order(i) }>192 , 'e' ); 
            %[Ysrcr,Yclsr] = cat_vol_resize({(Ysrc .* Yclsr) -
            %minYsrc,single(Yclsr)},'reduceV',vx_vol,1.1,16,'meanm'); %
            %RD202108 cause problems in PD ex-vivo
            try
              [Tmn,Tsd,Tn] = cat_stat_kmeans( Ysrcr(Yclsr(:)>0.9)  + minYsrc , clsn( i ) );
              T3th4(i) = Tmn(Tn==max(Tn));
            catch
              try
                [Tmn,Tsd,Tn] = cat_stat_kmeans( Ysrc(Ycls{order(i)}(:)>128) , clsn( i ) );
                T3th4(i) = Tmn(Tn==max(Tn));
              catch
                T3th4(i) = clsintv(order(i)); 
              end
            end
          end
          T3th4(4) = min(T3th4(4:5)); T3th4(5) = []; 
          T3th  = [min(Ysrc(:)) sort( T3th4 ) max(Ysrc(:)) ];
        end
        T3thx = [0,sort([BGth,1,2,3]),4];        
      end
      
      %% noise
      if ~exist('noise','var')
        Ynw   = cat_vol_localstat(Ysrc,Ycls{3}>128,2,4);
        Ync   = cat_vol_localstat(Ysrc,Ycls{2}>128,2,4); 
        noise = round(min(cat_stat_nanmean(Ynw(Ynw(:)>0)),cat_stat_nanmean(Ync(Ync(:)>0))) / min(abs(diff(T3th(1:3)))) * 10^6)/10^6; 
        clear Ynw nc; 
      end
      
% RD202006: Add contrast-rating to control WMHC and AMAP (CAT12.8?)

      if extopts.gcutstr == -2
        % ex-vivo CSF~BG
        T3thx(3) = 0.2; 
      end

      % Tth
      if ~exist('Tth','var')
        Tth.T3th  = T3th;
        Tth.T3thx = T3thx;
      end
  
      
      %% Ym 
      Ym = Ysrc - minYsrc;
      Tth.T3th = Tth.T3th - minYsrc; 
      for i=2:numel(T3th)
        M = Ysrc>T3th(i-1) & Ysrc<=T3th(i);
        Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3th(i-1))/diff(T3th(i-1:i))*diff(T3thx(i-1:i));
      end
      M  = Ysrc>=T3th(end); 
      Ym(M(:)) = numel(T3th)/6 + (Ysrc(M(:)) - T3th(i))/diff(T3th(end-1:end))*diff(T3thx(i-1:i));    
      Ym = Ym / 3; 
      clear M
      Tth.T3th = Tth.T3th + minYsrc; 
      
      inv_weighting = ~(T3th3(1)<T3th3(2) && T3th3(2)<T3th3(3)); % ~T1 
      
      
    end
  end

  
  if numel(nwarnings) < numel(cat_io_addwarning) 
    cat_io_cmd(' ','','',1); 
  end
end
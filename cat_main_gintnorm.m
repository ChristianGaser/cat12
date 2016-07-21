function [Ym,Yb,T3th3,Tth,inv_weighting,noise,cat_warnings] = cat_main_gintnorm(Ysrc,Ycls,Yb,vx_vol,res)
% This is an exclusive subfunction of cat_main.
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
%   [Ym,Yb,T3th3,Tth,inv_weighting,noise,cat_warnings] =
%     cat_main_gintnorm(Ysrc,Ycls,Yb,vx_vol,res)
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
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id$

  if isstruct(Ycls)
    
    %% final peaks and intesity scaling
    %  -----------------------------------------------------------------
    T3th  = Ycls.T3th;
    T3thx = Ycls.T3thx;


    % intensity scalling
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

    return
  end
    
  debug = cat_get_defaults('extopts.debug');
  inv_weighting = 0;
  if nargout==7
    cat_warnings = struct('identifier',{},'message',{});
  end
  vxv    = 1/mean(vx_vol);
  res.mn = round(res.mn*10^5)/10^5; 
  
  if cat_get_defaults('extopts.subfolders')
    reportfolder  = 'report';
  else
    reportfolder  = '';
  end
  
  %% initial thresholds and intensity scaling
  T3th3 = [mean(res.mn(res.lkp==3 & res.mg'>0.3)) ...
           mean(res.mn(res.lkp==1 & res.mg'>0.1)) ...
           mean(res.mn(res.lkp==2 & res.mg'>0.2))];
  T3th3 = round(T3th3*10^5)/10^5; 
  
  if T3th3(1)>T3th3(3) && T3th3(2)>T3th3(3) && T3th3(1)>T3th3(2) % invers (T2 / PD)
    cat_warnings = cat_io_addwarning(cat_warnings,...
      'CAT:cat_main:InverseContrast',...
      sprintf(['Inverse tissue contrast! ' ...
           '(Tissue Peaks: CSF=%0.2f, GM=%0.2f, WM=%0.2f)\n'],T3th3(1),T3th3(2),T3th3(3)),numel(cat_warnings)==0);
    
    % first initial scaling for gradients and divergence
    T3th3 = [max( res.mn(res.lkp==3 & res.mg'>0.3)) ...
             mean(res.mn(res.lkp==1 & res.mg'>0.1)) ...
             min( res.mn(res.lkp==2 & res.mg'>0.1))];
    T3th  = [min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))) ...
             max(res.mn(res.lkp==6 & res.mg'>0.1)) T3th3 ...
             mean([T3th3(3) max(res.mn(res.lkp==6 & res.mg'>0.1))]) ...
             min(res.mn(res.lkp==2 & res.mg'>0.1))*0.8 ...
             max(T3th3) + abs(diff(T3th3([1,numel(T3th3)])/2)) ...
             max(T3th3(end) + abs(diff(T3th3([1,numel(T3th3)])/2)) , ...
             max(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:))))) ];
    T3thx = [0,0.05, 1.4,2,3, 3.1, 2.9, 1.0, 0.7];
    
    [T3th,si] = sort(T3th);
    T3thx     = T3thx(si);
    
    Ym = Ysrc+0; 
    for i=2:numel(T3th)
      M = Ysrc>T3th(i-1) & Ysrc<=T3th(i);
      Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3th(i-1))/diff(T3th(i-1:i))*diff(T3thx(i-1:i));
    end
    M  = Ysrc>=T3th(end); 
    Ym(M(:)) = numel(T3th)/6 + (Ysrc(M(:)) - T3th(i))/diff(T3th(end-1:end))*diff(T3thx(i-1:i));    
    Ym = Ym / 3; 
    
    Yg    = cat_vol_grad(Ym,vx_vol);
    Ydiv  = cat_vol_div(Ym,vx_vol);
    
    %% tissues for bias correction
    Ycm   = (Ym + Yg - Ydiv + single(Ycls{2})/255)<2/3 | ...
            (Ycls{3} + Ycls{6} + Ycls{4} + Ycls{5})>128; 
    Ycd   = cat_vbdist(single(Ycm));
    Ybd   = cat_vbdist(cat_vol_morph(single((Ycls{6} + Ycls{4} + Ycls{5})>128),'lo',1));
    
    Ywm  = (single(Ycls{2})/255 - Yg - Ydiv - max(0,3-Ycd-Ybd/40)/2)>0.7 | ... 
           (Ym-Yg-Ydiv-max(0,3-Ycd-Ybd/40)/2)>0.8 & Ycls{1}+Ycls{2}>240;
    Ywm(smooth3(Ywm)<0.3)=0;
    Ygm  = (single(Ycls{1})/255 - abs(Ydiv)*8)>0.5 | ...
           (single(Ycls{2}+Ycls{1})>240 & max(0,2-Ycd - max(0,Ybd/2-10))>0 & abs(Ydiv)<0.1);

    %% bias correction
    [Yi,resT2] = cat_vol_resize(Ysrc.*Ywm,'reduceV',vx_vol,1,16,'min');
    Yig        = cat_vol_resize(Ysrc.*Ygm./median(Ysrc(Ygm(:)))*median(Ysrc(Ywm(:))),'reduceV',vx_vol,1,16,'meanm');
    Yi = max(Yi,Yig); Yi(Yig>0) = min(Yig(Yig>0),Yi(Yig>0));
    Yi = cat_vol_localstat(Yi,Yi>0,1,2);
    for xi=1:2, Yi = cat_vol_localstat(Yi,Yi>0,2,1); end
    Yi = cat_vol_approx(Yi,'nh',resT2.vx_volr,2); Yi = cat_vol_smooth3X(Yi,2); 
    Yi = cat_vol_resize(Yi,'dereduceV',resT2)./median(Yi(Ycls{2}>192));  
    Ysrcr = round(Ysrc ./ Yi * 10^5)/10^5; % * T3th3(3) * 1.05
    if debug==0, clear Yg Ydiv Yn Yi; end
    
    %% final thresholds
    T3th3 = [max( res.mn(res.lkp==3 & res.mg'>0.3)) ...
             mean(res.mn(res.lkp==1 & res.mg'>0.1)) ...
             min( res.mn(res.lkp==2 & res.mg'>0.1))];
    T3th  = [min(Ysrcr(~isnan(Ysrcr(:)) & ~isinf(Ysrcr(:)))) ...
             max(res.mn(res.lkp==6 & res.mg'>0.1)) ...
             T3th3 ...
             mean([T3th3(3) max(res.mn(res.lkp==6 & res.mg'>0.1))]) ...
             min(res.mn(res.lkp==2 & res.mg'>0.1))*0.8 ...
             max(T3th3) + abs(diff(T3th3([1,numel(T3th3)])/2)) ...
             max(T3th3(end) + abs(diff(T3th3([1,numel(T3th3)])/2)) , ...
             max(Ysrcr(~isnan(Ysrcr(:)) & ~isinf(Ysrcr(:))))) ];
    T3thx = [0,0.05,1.4,2,3,3.1, 2.9,1.0, 0.7];
    
    [T3th,si] = sort(T3th);
    T3thx     = T3thx(si);
    
    Tth.T3th  = T3th;
    Tth.T3thx = T3thx;
    
    inv_weighting = 1;
    
  elseif T3th3(1)<T3th3(2) && T3th3(2)<T3th3(3) % T1
    %%
    BGmin = min(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:)))); 
    BGcon = max([BGmin*1.1,T3th3(1) - mean(diff(T3th3)),median(Ysrc(Ycls{6}(:)>128))]);
    
    T3th3 = [max(max(res.mn(res.lkp==2 & res.mg'>0.1))*0.05,min(res.mn(res.lkp==3 & res.mg'>0.3))) ...
             max(res.mn(res.lkp==1 & res.mg'>0.1)) ...
             max(res.mn(res.lkp==2 & res.mg'>0.1))];
    T3th  = [BGmin ... minimum
             BGcon ... mean background (MT contrast with strong background noise)
             T3th3 ... csf gm wm 
             max(T3th3) + abs(diff(T3th3([1,numel(T3th3)])/2)) ... higher
             max(T3th3(end) + abs(diff(T3th3([1,numel(T3th3)])/2)) , ... maximum
              max(Ysrc(~isnan(Ysrc(:)) & ~isinf(Ysrc(:))))) ];
    T3thx = [0,0.05,1:5];
    Ysrcr = round(Ysrc*10^5)/10^5; 
    
    Tth.T3th  = T3th;
    Tth.T3thx = T3thx;
  else    
    error('CAT:cat_main:badTissueContrast',...
      sprintf('Bad tissue contrast (C=%0.2f, G=%0.2f, W=%0.2f)\n',...
        T3th3(1),T3th3(2),T3th3(3)),numel(cat_warnings)==0); %#ok<SPERR>
  end

  %% intensity scalling for gradient estimation
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
  T3th  = [median(Ysrcr(Ycls{3}(:)>192 & Yg(:)<0.20 & Ym(:)<0.45)) ...
           median(Ysrcr(Ycls{1}(:)>192 & Yg(:)<0.20)) ...
           median(Ysrcr(Ycls{2}(:)>192 & Yg(:)<0.10))];
  Yn    = cat_vol_localstat(Ysrc,Ycls{1}>192,2,4) + cat_vol_localstat(Ysrc,Ycls{2}>192,2,4); 
  noise = round(cat_stat_nanmean(Yn(Yn(:)>0)) / min(abs(diff(T3th(1:3)))) * 10^6)/10^6; 
  
 
  if debug==2
    [pth,nam] = spm_fileparts(res.image0(1).fname);
    tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',1,'gintnorm00'));
    save(tmpmat,'Ysrc','Ycls','Yb','vx_vol','res','T3th','T3thx','Yg','Ym','noise');
  end
  
  
  
  %% -------------------------------------------------------------------
  %  intensity checks and noise contrast ratio (contrast part 1)
  %  -------------------------------------------------------------------
  % relation between the GM/WM and CSF/GM and CSF/WM contrast has to be
  % greater that 3 times of the maximum contrast (max-min).
  clear Yn
  checkcontrast = @(T3th,minContrast) ...
    abs(diff(T3th([1,3]))) < (max(T3th(:))-min(T3th(:)))*minContrast || ...
    abs(diff(T3th(1:2)))   < (max(T3th(:))-min(T3th(:)))*minContrast || ...
    abs(diff(T3th(2:3)))   < (max(T3th(:))-min(T3th(:)))*minContrast;
  if checkcontrast(T3th3,1/9) && exist('cat_warnings','var') % contrast relation
    cat_warnings = cat_io_addwarning(cat_warnings,...
      'CAT:cat_main:LowContrast',...
      sprintf(['The contrast between the tissues is extremely low! ' ...
           '(C=%0.2f, G=%0.2f, W=%0.2f)\n'],T3th(1),T3th(2),T3th(3)),numel(cat_warnings)==0);
  end
  

  %  -------------------------------------------------------------------
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
 
    % spm tissue peaks
    T3th_spm = [min(res.mn(res.lkp==3 & res.mg'>0.3)) ...
                max(res.mn(res.lkp==1 & res.mg'>0.1)) ...
                max(res.mn(res.lkp==2 & res.mg'>0.1))];
    
    
    
    % check SPM segmentation
    if exist('cat_warnings','var')
      Ymx = single(Ycls{1})/255*2/3 + single(Ycls{2})/255+ single(Ycls{3})/255*1/3;  
      Ygw = Yb & ((Ycls{1}+Ycls{2})>128);
      Ymp0diff = sqrt(mean(Ym(Ygw(:)) - Ymx(Ygw(:)))^2); 
      if Ymp0diff>0.10 && debug
        cat_warnings = cat_io_addwarning(cat_warnings,...
          'CAT:cat_main:badSPMsegment',sprintf(...
          ['SPM segmentation does not fit to the image (RMS(Ym,Yp0)=%0.2f).\n'...
           'This can be an alignment problem (check origin), ' ...
           'untypical subjects (neonates, non-human),\n'...
           'bad image contrast (C=%0.2f,G=%0.2f,W=%0.2f), \n'...
           'low image quality (NCR~%0.2f), or something else ...'],Ymp0diff,T3th,noise),numel(cat_warnings)==0); 
      end
      clear Ymx;
    end
    
    
    %% skull-stripping warning
    skulltest = (median(Ysrc(Ycls{5}(:)>192 & Ysrc(:)>T3th(2))) < ... 
       median(Ysrc(Ycls{3}(:)>192 & Ysrc(:)>0))); 
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
        CSFth = min([nanmedian(Ysrc(Ycls{3}(:)>240 & Ysrc(:)>0)), ... 
                     nanmedian(Ysrc(Ycls{3}(:)>192 & Ysrc(:)>0)), ... 
                     nanmedian(Ysrc(Ycls{3}(:)>128 & Ysrc(:)>0)), ...
                     mean(Ysrc(Ysrc>0))*0.5])*0.9; % 
        Ysrc  = cat_vol_laplace3R(max(CSFth,Ysrc),Yb & Ysrc==0,0.2) .* Yb;
         cat_warnings = cat_io_addwarning(cat_warnings,...
           'CAT:cat_main:SkullStripped',...
           'Skull-stripped input image detected! Try boundary cleanup.',numel(cat_warnings)==0);  
      else
         cat_warnings = cat_io_addwarning(cat_warnings,...
           'CAT:cat_main:SkullStripped',...
           'Skull-stripped input image?',numel(cat_warnings)==0); 
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
    BGcon = max([BGmin*1.1,T3th3(1) - mean(diff(T3th3)),median(Ysrc(Ycls{6}(:)>128))]);
    BMth  = BGcon; %max(0.01,cat_stat_nanmedian(Ysrc(Ybm(:))));
    Ywm   = (Ycls{2}>128  & Yg<gth) | ((Ym-Ydiv*2)>(1-0.05*mean(vx_vol)) & Yb2); % intensity | structure (neonate contast problem)
    Ycm   = smooth3((Ycls{3}>240 | Ym<0.4) & Yg<gth*3 & Yb & ~Ywm & Ycls{1}<8 & Ysrc>BMth & Ym<0.5)>0.5; % important to avoid PVE!

    % If SPM get totaly wrong maps due to bad image orientations our 
    % segment were incorrect too (or empty) and peak estimation fail.
    % I try to use the kmeans, but in WM it is affected by WMHs, in 
    % CSF by blood vessels and meninges and in GM noise and subcortical
    % structures were problematic. In ADHD/..NYC..14 the basal structes 
    % get the average peak and the cortex was detected as CSF. There 
    % were much more images with smaller problems ...
    Ysrcr  = round( Ysrc.*10^5 ) / 10^5;
    WMth   = cat_stat_nanmedian(Ysrcr(Ywm(:))); % kmeans3D(Ysrc(Ycls{2}(:)>192 & Yg(:)<gth),1); % GM/WM WM  
    CSFth  = cat_stat_nanmedian(Ysrcr(Ycm(:))); % kmeans3D(Ysrc(Ycls{3}(:)>64 & Yg(:)>gth & Yb(:)),2); % CSF CSF/GM
      %  0.05 <<<<< BMth + 4*cat_stat_nanstd(Ysrc(Ybm(:)))
    Ybg    = cat_vol_morph(Yg<0.10 & Yb & Ysrc<WMth*(1-0.03*mean(vx_vol)) & Ysrc>CSFth*1.5 & Ycls{3}<64,'o',2);
    Ygm    = ~Ybg & Yg<0.4 & Ysrc<(WMth+0.9*diff([CSFth,WMth])) & Yg<gth*2 & ...
      Ysrc>(CSFth+0.1*diff([CSFth,WMth])) & ~Ywm & ~Ycm & Yb & abs(Ydiv)<0.1; 
    %Ygm   = Ygm | (Ycls{1}>64 & Ybg & ~Ywm);
    GMth   = cat_stat_nanmedian(Ysrcr(Ygm(:))); %kmeans3D(Ysrc(Ygm(:)),3); % CSF/GM GM GM/WM
    T3th_cls  = round([CSFth(1) GMth(1) WMth(1)]*10^4)/10^4;
    %clear Ybg
   %
    if any(isnan(T3th_cls)) 
      fprintf('\n');
      error('CAT:cat_main:cat_pre_gintnorm:nobrain',...
        'Bad SPM-Segmentation. Check image orientation!');
    end
    % median tissue peaks
    
   
    % print a warning for strong variation in the peaks
    % deactivated 2016-04-06 ... remove in final version
    %{
    T3th_diff = min(T3th_spm,T3th_cls) ./ max(T3th_spm,T3th_cls);
    if any(T3th_diff < [0.5 0.5 0.95]) 
      if exist('cat_warnings','var')
         cat_warnings = cat_io_addwarning(cat_warnings,...
          'CAT:cat_main:DiffTissuePeaks',...
          sprintf(['Peaks of SPM segmentation do not fit to average segmentation intensity!\\n '...
            'spm[%3.2f,%3.2f,%3.2f] ~= [%3.2f,%3.2f,%3.2f]median'],...
            T3th_spm/max(T3th_spm),T3th_cls/max(T3th_cls)),numel(cat_warnings)==0);
      end
      %T3th3 = T3th_cls;
    else
      if max(res.mn(res.lkp==5 & res.mg'>0.1)) < mean(res.mn(res.lkp==3 & res.mg'>0.3)), fprintf('\n'); end
      %T3th3 = T3th_spm;
    end
    %} 
   
    if debug==2
      tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',1,'gintnorm01'));
      save(tmpmat,'Ysrc','Ycls','Yb','vx_vol','res','T3th','Yg','Ydiv','Ym',...
       'Yb2','gth','Ybm','BMth','Ywm','Ygm','Ycm','Ybg','T3th_cls','T3th_spm','noise');
    end
    
    
    %% final peaks and intesity scaling
    %  -----------------------------------------------------------------
    T3th3 = T3th_cls;
    T3th  = [min(Ysrcr(~isnan(Ysrcr(:)) & ~isinf(Ysrcr(:)))) BMth T3th3 ...
              T3th3(end) + diff(T3th3([1,numel(T3th3)])/2) ... WM+
              max(T3th3(end)+diff(T3th3([1,numel(T3th3)])/2) , ... max
              max(Ysrcr(~isnan(Ysrcr(:)) & ~isinf(Ysrcr(:))))) ];
    T3thx = [0,0.02,1:5];


    % intensity scalling
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
  end
    
  
  %% if there was a warning we need a new line 
  if nargout==7 && numel(cat_warnings)>1, fprintf('\n'); cat_io_cmd(' ','','',1); end

end
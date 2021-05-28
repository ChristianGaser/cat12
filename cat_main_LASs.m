function [Yml,Ymg,Ycls] = cat_main_LASs(Ysrc,Ycls,Ym,Yb,Yy,Tth,res,vx_vol,extopts) 
% This is an exclusive subfunction of cat_main.
% ______________________________________________________________________
%
% Local Adaptive Segmentation (LAS):
%
% This version of the local adaptive intensity correction includes a  
% bias correction that based on a maximum filter for the WM and a mean
% filter of the GM to stabilize the correction in region with less WM.
%
% The extension based mostly on the assumption that the tissue next to 
% the CSF (and high divergence sulci) has to be WM (maximum, high 
% divergence) or GM. For each tissue a refined logical map is generated 
% and used to estimate the local intensity threshold.
%
% It is important to avoid high intensity blood vessels in the process, 
% because they will push down local WM and GM intensity. 
%
% There are further regionwise correction, e.g. , to avoid overfitting in 
% cerebellum, or adapt for age specific changes, e.g. enlarged ventricle.
%
% Based on this values a intensity transformation is used. Compared to 
% the global correciton this has to be done for each voxel. To save time
% only a rough linear transformation is used.
% ______________________________________________________________________
%
%   [Yml,Ymg,Yclsu] = ...
%     cat_main_LAS(Ysrc,Ycls,Ym,Yb0,Yy,T3th,res,vx_vol,extopts,Tth)
%
%   Yml     .. local intensity correct image
%   Ymg     .. global intensity correct image
%   Yclsu   .. corrected SPM tissue class map
%
%   Ysrc    .. (bias corrected) T1 image
%   Ycls     .. SPM tissue class map
%   Ym      .. intensity corrected T1 image (BG=0,CSF=1/3,GM=2/3,WM=1)
%   Yb      .. brain mask
%   Yy      .. deformation map for the cat atlas map
%   Tth     .. structure with tissue thresholds of CSF, GM, and WM in Ysrc
%   res     .. SPM segmentation structure
%   vx_vol  .. voxel dimensions
%   extopts .. cat options
% ______________________________________________________________________
% 
% internal maps:
%
%   Yg   .. gradient map   - edges between tissues
%   Ydiv .. divergence map - sulci, gyris pattern, and blood vessels
%   Yp0  .. label map      - tissue classes (BG=0,CSF=1,GM=2,WM=3) 
%
%   Ycm  .. CSF
%   Ygm  .. GM
%   Ywm  .. WM 
%
%   Yvt  .. WM next to the ventricle map 
%   Ygmt .. cortical thickness map
%   Ypp  .. cortical percentage position map
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  % set this variable to 1 for simpler debuging without reduceBrain
  % function (that normally save half of processing time)
  verb    = extopts.verb-1;
  NS      = @(Ys,s) Ys==s | Ys==s+1;          % function to ignore brain hemisphere coding
  LASstr  = max(eps,min(1,extopts.LASstr));   % LAS strenght (for GM/WM threshold)3
  LAB     = extopts.LAB;                      % atlas labels
  mvx     = mean(vx_vol);                     % mean voxel volume to correct for morphological operations  
  Tth.T3th(1) = min(0,Tth.T3th(1));           % correction of the background value
  
  % set debug = 1 and do not clear temporary variables if there is a breakpoint in this file 
  dbs   = dbstatus; debug = 0; 
  for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,'cat_main_LASs'); debug = 1; break; end; end
  
  
%% ---------------------------------------------------------------------
%  First, we have to optimize the segments using further information that 
%  SPM do not use, such as the gradient, divergence and distance maps. 
%  The gradient map (average of the first derivate of the T1 map) is an 
%  edge map and independent of the image intensity. It helps to avoid PVE 
%  regions and meninges. 
%  The divergence (second derivate of the T1 map) help to identfiy sulcal
%  and gyral pattern and therefore to find WM and CSF regions for further 
%  corrections and to avoid meninges and blood vessels. 
%  Furhtermore, special assumption can be used. 
%  The first one is the maximum property of the WM in T1 data that allows
%  using of a maxim filter for the GM/WM region. 
%  The second is the relative stable estimation of CSF/BG that allows to 
%  estimat a distance map. Because, most regions have a thin layer of 
%  GM around the WM we can avoid overestimation of the WM by the other 
%  maps (especially the divergence). 
%  ---------------------------------------------------------------------
  fprintf('\n'); 
  stime = cat_io_cmd('  Prepare maps','g5','',verb); 

  % map atlas to RAW space
  % -------------------------------------------------------------------
  % avoid read error in parallel processing
  for i=1:5
    try
      Vl1A = spm_vol(extopts.cat12atlas{1}); break
    catch 
      pause(1)
    end
  end
  Yl1 = cat_vol_ctype( cat_vol_sample(res.tpm(1),Vl1A,Yy,0) );
  if ~debug, clear Yy; end

  
  
  % lower resolution to save time and space
  rres = 0.9; 
  [Ym,resTb] = cat_vol_resize(Ym        ,'reduceV',vx_vol,rres,32,'meanm'); 
  if any(resTb.vx_vol ~= resTb.vx_volr), Ysrco = Ysrc+0; end % save orignal image for later correction
  Ysrc       = cat_vol_resize(Ysrc      ,'reduceV',vx_vol,rres,32,'meanm'); 
  Yb         = cat_vol_resize(single(Yb),'reduceV',vx_vol,rres,32,'meanm')>0.5; 
  Yl1        = cat_vol_resize(Yl1       ,'reduceV',vx_vol,rres,32,'nearest'); 
  for i=1:6, Ycls{i} = cat_vol_ctype(cat_vol_resize(single(Ycls{i}),'reduceV',vx_vol,rres,32,'meanm')); end
  if debug, Ymo = Ym; Ybo = Yb; end  %#ok<NASGU>
  vx_vol = resTb.vx_volr;
  
  % help maps to detect edges (Yg) and sulci/gyris (Ydiv)
  Yg    = cat_vol_grad(Ym,vx_vol);                                                  % mean gradient map  ...  ./max(0.1,Ym)
  Ydiv  = cat_vol_div(Ym,vx_vol);                                                   % divergence map 
  Yp0   = single(Ycls{1})/255*2 + single(Ycls{2})/255*3 + single(Ycls{3})/255;      % tissue label map
  Yp0   = max( cat_vol_morph(Yb,'lc',2),Yp0); 
  noise = std(Ym(Ycls{2}>192)); 
  
  noisef = max(0,min(0.3,noise./mvx)); 
  Ym  = cat_vol_smooth3X(Ym,noisef); 
  Yp0 = cat_vol_smooth3X(Yp0,noisef); 
  for i=1:6, Ycls{i} = cat_vol_ctype(cat_vol_smooth3X(single(Ycls{i}),noisef)); end
 
  Yb = smooth3(Yb | (cat_vol_morph(Yb,'d',2/mvx) & Ym<0.8 & Yg<0.3 & Ym>0 & Yp0>0.2))>0.5;  % increase brain mask, for missing GM 

  % correction for negative (and positive) values
  srcmin = min(Ysrc(:)); Ysrc = Ysrc - srcmin; Tth.T3th = Tth.T3th - srcmin;
  Tthc = Tth; Tthc.T3th = Tth.T3th - srcmin;
  
  
  %% Optimization of the tissue segments
  %  ----------------------------------------------------------------------
  %  The correction of the SPM tissue classification is required especially
  %  in cases of failed bias correction, that strongly based on the T1
  %  intensities that were corrected before. So this is main part of LAS.
  %  ----------------------------------------------------------------------
  if 1;; % also this if is just for the debuging mode and contain the second block the correction of the segmentation 
    %% brain segmentation can be restricted to the brain to save time 
    stime = cat_io_cmd('  Prepare partitions','g5','',verb,stime); 
    [Ymr,Yb,BB] = cat_vol_resize({Ym,Yb},'reduceBrain',vx_vol,round(10/mvx),Yb);
    [Ygr,Ydivr,Yp0]   = cat_vol_resize({Yg,Ydiv,Yp0},'reduceBrain',vx_vol,BB.BB);
    Yl1               = cat_vol_resize(Yl1          ,'reduceBrain',vx_vol,round(4/mvx),BB.BB);
    Yclsr = cell(size(Ycls)); for i=1:6, Yclsr{i} = cat_vol_resize(Ycls{i},'reduceBrain',vx_vol,BB.BB); end
    

    % adaption of the LASstr depending on average basal values ... not in T2/PD/FLAIR because of unknown variation 
    % LASmod  = min(2,max(0,cat_stat_kmeans((Ymr( NS(Yl1,LAB.BG) & Ygr<0.1 & Ydivr>-0.05  & Yclsr{1}>4)) - 2/3) * 8));
    % LASstr  = min(1,max(0.05,LASstr * LASmod)); clear LASmod                 % adaption by local BG variation
    LASfs   = 2 / max(0.05,LASstr);                                          % smoothing filter strength 
    

    %  GM thickness (Ygmt) and percentage position map (Ypp) estimation
    %  The Ypp and Ygmt maps are used to refine the GM especially to correct
    %  highly myelinated GM regions. Using allow to avoid overcorrections.
    fastppi = 1; % 1 mm: 15 vs. 40 seconds ... it thing fast is ok here 20161014 
    if fastppi
      [Ygmt,Ypp] = cat_vol_pbt( (Yp0 + (Ymr*3 .* (Yp0>0)))/2 , ...
        struct('resV',vx_vol,'verb',0,'dmethod','vbdist','method','pbt2x') );
    else
      [Ygmt,Ypp] = cat_vol_pbt( (Yp0 + (Ymr*3 .* (Yp0>0)))/2 , ...
        struct('resV',vx_vol,'verb',0) ); 
    end
    Ygmtroi = Ygmt>0 & Ygmt<6 & NS(Yl1,LAB.CT); 
    [GMTstat(1), GMTstat(2)] = cat_stat_kmeans(Ygmt(Ygmtroi(:))); 
    [D,I]   = cat_vbdist(single(Ygmt>0.1),Yp0>0,vx_vol); Ygmt = Ygmt(I); clear D I;  %#ok<ASGLU> % full GMT map for the whole brain
    if ~debug, clear Ygmtroi; end


    %% As far as SPM segmentation is not optimal we need some refinements.
    %  ----------------------------------------------------------------------
    stime = cat_io_cmd('  Prepare segments','g5','',verb,stime);

    % use SPM segment 
    Ygm = Yclsr{1}>128;
    Ywm = Yclsr{2}>128;
    Ycm = Yclsr{3}>128 | (Yp0==1 & Yclsr{3}>8);
    
    % remove small elements / noise
    Ywm = cat_vol_morph(Ywm,'l',[100 0.01])>0; 
    Ycm = cat_vol_morph(Ycm,'l',[100 0.01])>0; 

    % update subcortical structures 
    % (this can introduce some bias i
    Yss = (NS(Yl1,LAB.BG) | NS(Yl1,LAB.TH)) & Ygr<0.1 & Ydivr>-0.05; 
    Yss = smooth3(Yss)>0.5;
    Ygm = Ygm | Yss;
    Ywm = Ywm & ~Yss & ~(NS(Yl1,LAB.BG) | NS(Yl1,LAB.TH));

    % blood vessel correction
    Ybv = (Ymr - Ydivr)>1.2 & Yclsr{2}<192 & Yclsr{3}<192 & Yclsr{1}<192;
    Ybv(smooth3(Ybv)<0.3) = 0; 
    Ywm = Ywm & ~Ybv; 
    Ygm = Ygm & ~Ybv; 
    Ycm = Ycm & ~Ybv; 

    % tissue thresholds
    [gmm,gms,gmn] = cat_stat_kmeans(Ymr(Ygm),6); %#ok<ASGLU> % we need 3 GM types to include GM/WM and GM/CSF PVE voxels 
    [wmm,wms,wmn] = cat_stat_kmeans(Ymr(Ywm),2); %#ok<ASGLU> % here we need more to exclude other tissues
    [cmm,cms,cmn] = cat_stat_kmeans(Ymr(Ycm),2); %#ok<ASGLU> % here we need more to exclude other tissues
    TSth = [mean(gmm(2:4)) mean(wmm(wmn>0.5)) mean(cmm(cmn>0.5))]; 
    clear gms wms cms gmn; 
      
      
    %% correction for WMHs
    Ypm  = single(Yclsr{1})/255*TSth(1) + single(Yclsr{2})/255*TSth(2) + single(Yclsr{3})/255*TSth(3); 
    Yvt  = cat_vol_morph( NS(Yl1,LAB.CB) | NS(Yl1,LAB.BS) | NS(Yl1,LAB.PH) | ...
                          NS(Yl1,LAB.BG) | NS(Yl1,LAB.TH) | NS(Yl1,LAB.HC)  ,'d',5/mvx); 
    Ywmh = smooth3( (Ypm/3 - Ymr)>0.1 & ~Yss & ~Yvt)>0.7;
    Ywmh = cat_vol_morph(Ywmh,'l',[100 50]); 
    Ywm  = Ywm & ~Ywmh; 
    Ygm  = Ygm & ~Ywmh; 
    if ~debug, clear Yss Ypm; end

    % parahippocampale gyrus and hippocampus
    Yhcg = (NS(Yl1,LAB.PH) | NS(Yl1,LAB.HC)) & (Ymr - Ydivr)>2.5/3;
    Ywm  = Ywm | Yhcg; 
    Ygm  = Ygm & ~Yhcg;
    if ~debug, clear Yhcg; end

    % Ycsf
    Yngm = Ygm & Ypp<0.1 & Ymr<0.5 & Ygmt<=GMTstat(1); 
    Ygm  = Ygm & ~Yngm;
    if ~debug, clear Yngm Ygmt Ypp; end

    % remove sinus (rectus) and brainstem
    Ysr = cat_vbdist(single(NS(Yl1,LAB.CB)),NS(Yl1,LAB.CT) | NS(Yl1,LAB.BV),vx_vol)<30 & ...
          cat_vol_morph(NS(Yl1,LAB.BV),'d',10/mvx) & ...
          (NS(Yl1,LAB.CT)) & Ymr>0.5/3 & Ymr<1.8/3 & ~cat_vol_morph(NS(Yl1,LAB.CT) & Ymr>2/3,'d',1);
    Ysr(smooth3(Ysr)<0.5)=0; 
    Ygm = Ygm & ~Ysr & ~cat_vol_morph(NS(Yl1,LAB.BS),'d',2/mvx); 
    if ~debug, clear Ysr; end

    % non CSF
    Ycm = Ycm & abs(Ydivr)<0.2 & Ygr<0.2; 
    %Ycm(smooth3(Ycm)<0.5)=0;  
    if ~debug, clear Ymr Ydivr Ygr Yclsr; end 

    
    %% back to original resolution for full bias field estimation
    [Ycm,Ygm,Ywm] = cat_vol_resize({Ycm,Ygm,Ywm},'dereduceBrain',BB); 
    [Yp0,Yl1]     = cat_vol_resize({Yp0,Yl1},'dereduceBrain',BB); %#ok<ASGLU>
    [Yb]          = cat_vol_resize({Yb},'dereduceBrain',BB);
  
    if ~debug, clear  Yg; end  %Ybg
  end
  
  
  
  
  %% Estimation of local peaks and creation of normalized T1 maps
  %  --------------------------------------------------------------------- 
  %  Estimation of the local WM threshold with "corrected" GM voxels to
  %  avoid overfitting (see BWP cerebellum). 
  %  CSF is problematic in high contrast or skull-stripped image should 
  %  not be used here, or in GM peak estimation
  %  ---------------------------------------------------------------------
  if 1;; % the last if block for debuging mode 
    mres  = 1.1; 
    stime = cat_io_cmd('  Estimate local tissue thresholds','g5','',verb,stime);  if debug, fprintf('\n'); end
    stime2 = cat_io_cmd('    WM intensity','g5','',debug);  
    %%
    if     TSth(2)>TSth(1) && TSth(2)>TSth(3), wmtht = 3; % T1 - WM is maximum
    elseif TSth(2)<TSth(1) && TSth(2)<TSth(3), wmtht = 2; % T2 - WM is minimum
    else                                       wmtht = 1; 
    end
    
    %% tissue thresholds
    % you need to know which are the extrem classes to define there PVE
    [gmm,gms,gmn] = cat_stat_kmeans(Ysrc(cat_vol_morph(Ygm & Ycls{1}>128,'e')),3); % we need 3 GM types to include GM/WM and GM/CSF PVE voxels 
    [wmm,wms,wmn] = cat_stat_kmeans(Ysrc(cat_vol_morph(Ywm & Ycls{2}>128,'e')),2); % here we need more to exclude other tissues
    [cmm,cms,cmn] = cat_stat_kmeans(Ysrc(cat_vol_morph(Ycm & Ycls{3}>128,'e')),2); % here we need more to exclude other tissues (T2 with high CSF vals)
    TSth = [mean(gmm(gmn>0.4)) mean(wmm(wmn>0.5)) mean(cmm(cmn>0.5))];
    TSmx = (TSth == max(TSth)) - (TSth == min(TSth)); 
    TSth = TSth + TSmx .* [mean(gms(wmn>0.5)) mean(wms(wmn>0.5)) mean(cms(cmn>0.5))]; 

    %%
    if 1
      Ysrcm  = cat_vol_localstat(Ysrc .* Ywm,Ywm,1,wmtht); 
      Ysrcs  = cat_vol_localstat(Ysrcm,Ysrcm>0,2,4);
      Ysrcm(Ysrcs>mean(Ysrcs(Ysrcs(:)>0) + std(Ysrcs(Ysrcs(:)>0) ))) = 0;
      Ysrcm  = cat_vol_noPVE(Ysrcm,vx_vol,2); 
      Ysrca  = Ysrcm; 
      Ysrcmw = cat_vol_approx(Ysrcm,'nh',vx_vol,4); Ysrcmw = cat_vol_smooth3X(Ysrcmw,2); 

      Ysrcm  = cat_vol_localstat(Ysrc ./ TSth(1) .* TSth(2) .* (Ycls{1}>4),(Ygm & Ycls{1}>4),1,1);%Ygm & 
      Ysrcm  = cat_vol_noPVE(Ysrcm,vx_vol,2); 
      Ysrcs  = cat_vol_localstat(Ysrcm,Ysrcm>0,4,4);
      Ysrcm(Ysrcs>mean(Ysrcs(Ysrcs(:)>0)+std(Ysrcs(Ysrcs(:)>0) ))) = 0;
      Ysrca  = Ysrca + Ysrcm; 
      Ysrcmg = cat_vol_approx(Ysrcm,'nh',vx_vol,4); Ysrcmg = cat_vol_smooth3X(Ysrcmg,LASfs); 

      Ysrcm  = cat_vol_localstat(Ysrc ./ TSth(3) .* TSth(2) .* (Ycm & Ycls{3}>4),(Ycm & Ycls{3}>4),1,1); 
      Ysrcm  = cat_vol_noPVE(Ysrcm,vx_vol,2); 
      Ysrcs  = cat_vol_localstat(Ysrcm,Ysrcm>0,2,4);
      Ysrcm(Ysrcs>mean(Ysrcs(Ysrcs(:)>0)+std(Ysrcs(Ysrcs(:)>0) ))) = 0;
      Ysrca  = Ysrca + Ysrcm; 
      Ysrcmc = cat_vol_approx(Ysrcm,'nh',vx_vol,4); Ysrcmc = cat_vol_smooth3X(Ysrcmc,LASfs); 
   
      Ysrcma = cat_vol_approx(Ysrca,'nh',vx_vol,2); Ysrcma = cat_vol_smooth3X(Ysrcma,LASfs); 
      Ysrcs  = cat_vol_localstat(Ysrcma,Ysrcma>0,2,4);
      Ysrcma(Ysrcs>mean(Ysrcs(Ysrcs(:)>0)+std(Ysrcs(Ysrcs(:)>0) ))) = 0;
      
      Yi = cat_vol_approx(Ysrcma,'nh',vx_vol,2); Yi = cat_vol_smooth3X(Yi,2); 
      Ylab{2} = Yi; 
    else
     
       Ysrcm  = cat_vol_localstat(Ysrc .* Ywm,Ywm,1,wmtht) + ...
              cat_vol_localstat(Ysrc ./ TSth(1) .* TSth(2) .* (Ygm & Ycls{1}>128),(Ygm & Ycls{1}>128),1,1) + ...
              cat_vol_localstat(Ysrc ./ TSth(3) .* TSth(2) .* (Ycm & Ycls{3}>240),(Ycm & Ycls{3}>240),1,1); 
      Ysrcs = cat_vol_localstat(Ysrcm,Ysrcm>0,2,4);
      Ysrcm(Ysrcs>mean(Ysrcs(Ysrcs(:)>0)+std(Ysrcs(Ysrcs(:)>0) ))) = 0;
      Ysrcm2 = Ysrcm; 
      YM = cat_vol_morph(Ysrcm>0,'d') & Yp0>1 & Ysrcm==0; Ysrcm2(YM) = Ysrc(YM);
      Ysrcm2 = cat_vol_localstat(Ysrcm2,Ysrcm2>0,2,wmtht); Ysrcm(YM) = Ysrcm2(YM);
      Ysrcm  = cat_vol_noPVE(Ysrcm,vx_vol,2); 

      % major correction outside the brain 
      [Yi ,resT2] = cat_vol_resize(Ysrcm,'reduceV',vx_vol,mres,32,'meanm'); % maximum reduction for the WM
      if ~debug, clear Ysrcm Ydiv; end
      Yi(smooth3(Yi>0)<0.5)=0; Yi(smooth3(Yi>0)<0.5)=0;

      % Yi(Yi==0 & Ygi>0)=Ygi(Yi==0 & Ygi>0); if ~debug, clear Ygi; end
      Yi = cat_vol_approx(Yi,'nh',resT2.vx_volr,4); Yi = cat_vol_smooth3X(Yi,LASfs); 
      Ylab{2} = max(eps,cat_vol_resize(Yi,'dereduceV',resT2)); 
    end
    if ~debug, clear Yi; end


    %% GM
    stime2 = cat_io_cmd('    GM intensity','g5','',debug,stime2);  
    Yi = Ysrc ./ max(eps,Ylab{2}) .* (Ycls{1}>128) .* (Yp0>0);
    Yi = cat_vol_noPVE(Yi,vx_vol,1);
    Yi  = cat_vol_localstat(Yi,Yi>0,4,4);
    Yi(Yi>mean(Yi(Yi(:)>0)+std(Yi(Yi(:)>0) ))) = 0;
    if ~debug, clear Yhdm; end
    [Yir,resT2] = cat_vol_resize(Yi,'reduceV',vx_vol,mres,32,'meanm'); if ~debug, clear Yi; end
    Yir = cat_vol_noPVE(Yir,vx_vol,1);
    Yir = cat_vol_approx(Yir,'nh',resT2.vx_volr,2); 
    Yir = cat_vol_smooth3X(Yir,LASfs * 2); 
    Ylab{1} = cat_vol_resize(Yir,'dereduceV',resT2) .* Ylab{2};  
    if mean(Ylab{1}(:))<mean(Ylab{2}(:))
      Ylab{1} = min(Ylab{1},Ylab{2}*0.99); 
    else
      Ylab{1} = max(Ylab{1},Ylab{2}/0.99); 
    end
    if ~debug, clear Yir Ybd Yl1; end


    %% CSF & BG
    % no adaption yet
    Yi = Ysrc ./ max(eps,Ylab{2}) .* (smooth3(Ycm)>.6) .* (Yp0>0);
    Yi = cat_vol_noPVE(Yi,vx_vol,2,3,2);
    if ~debug, clear Yhdm; end
    [Yir,resT2] = cat_vol_resize(Yi,'reduceV',vx_vol,mres/2,32,'meanm'); if ~debug, clear Yi; end
    Yir = cat_vol_noPVE(Yir,vx_vol,2);
    Yir = cat_vol_approx(Yir,'nh',resT2.vx_volr,2); 
    Yir = cat_vol_smooth3X(Yir,LASfs * 4); 
    Ylab{3} = cat_vol_resize(Yir,'dereduceV',resT2) .* Ylab{2};  
    
    %Ylab{3} = cmm2(cmn2>0.5) .* Ylab{2}; %cat_vol_smooth3X(cat_vol_resize(Yca ,'dereduceV',resT2).*Ylab{2},LASfs*2);  
    if mean(Ylab{3}(:)) < mean(Ylab{1}(:)) 
      Ylab{3} = min(Ylab{3},Ylab{1}*0.99); 
    else
      Ylab{3} = max(Ylab{3},Ylab{1}/0.99); 
    end
    if mean(Ylab{3}(:)) < mean(Ylab{2}(:))
      Ylab{3} = min(Ylab{3},Ylab{2}*0.99); 
    else
      Ylab{3} = max(Ylab{3},Ylab{2}/0.99); 
    end
    %}
    Ysrcbc  = Ysrc ./ Ylab{2}; 
    Ylab{6} = min( cat_stat_kmeans( Ysrcbc( Ycls{6}(:)>128 ) ) , cat_stat_kmeans( Ysrcbc( Ycls{4}(:)>128 ) ) ) .* Ylab{2};
    clear Ynba Yca; 
    
    %% back to original resolution
    if any(resTb.vx_vol ~= resTb.vx_volr) 
      Ysrc = Ysrco; clear Ysrco; 
      for i=1:6, if ~isempty(Ylab{i}), Ylab{i} = cat_vol_resize(Ylab{i},'dereduceV',resTb); end; end
      Yb  = cat_vol_resize(Yb ,'dereduceV',resTb); %#ok<NASGU>
      Yp0 = cat_vol_resize(Yp0,'dereduceV',resTb); %#ok<NASGU>
      Ywm = cat_vol_resize(Ywm,'dereduceV',resTb); %#ok<NASGU>
    end 
    %%
    labmn = zeros(1,6); for i=[1 2 3 6], if ndims(Ylab{i})==3, labmn(i) = mean(Ylab{i}(Ycls{i}(:)>128)); else,  labmn(i) = mean(Ylab{i}); end; end
    [labmn,lo] = sort(labmn,'descend'); %clear labmn;  %#ok<ASGLU>
    
    % local intensity modification of the original image
    % --------------------------------------------------------------------
    stime2 = cat_io_cmd('    Intensity mapping','g5','',debug,stime2);  
    Yml = zeros(size(Ysrc));  mc = mean(abs(diff(Tthc.T3th)))/100; 
    Yml = Yml + ( (Ysrc>=Ylab{lo(1)}                    ) .* (3.0 + (Ysrc-Ylab{lo(1)}) ./ max(mc,abs(Ylab{lo(1)}-Ylab{lo(2)})) ));
    Yml = Yml + ( (Ysrc>=Ylab{lo(2)} & Ysrc<Ylab{lo(1)} ) .* (2.0 + (Ysrc-Ylab{lo(2)}) ./ max(mc,abs(Ylab{lo(1)}-Ylab{lo(2)})) ));
    Yml = Yml + ( (Ysrc>=Ylab{lo(3)} & Ysrc<Ylab{lo(2)} ) .* (1.0 + (Ysrc-Ylab{lo(3)}) ./ max(mc,abs(Ylab{lo(2)}-Ylab{lo(3)})) ));
    Yml = Yml + ( (Ysrc>=Ylab{lo(4)} & Ysrc<Ylab{lo(3)} ) .* (0.1 + (Ysrc-Ylab{lo(4)}) ./ max(mc,abs(Ylab{lo(3)}-Ylab{lo(4)}))*0.9));
    Yml = Yml + ( (Ysrc< Ylab{lo(4)} )                    .* (0.0 + (Ysrc)             ./ max(mc,abs(Ylab{lo(3)}-Ylab{lo(4)}))*0.1));
    Yml(isnan(Yml) | Yml<0)=0; Yml(Yml>10)=10;
    Yml = Yml / 3;

    %% update of the global intensity normalized map
    % ######## RD202105: this is not optimal #########
    Ymg = Ysrc ./ max(eps,Ylab{2}) * Tth.T3th(5);
    %Ymg = Ysrc ./ max(eps,Ylab{2}); 
    %Ymg = Ymg * max( Tthc.T3th(Tthc.T3thx>0 & Tthc.T3thx<4) ) / max( max(Tthc.T3thx(Tthc.T3thx>0 & Tthc.T3thx<4))/3); % RD202004: corrected srcmin-correction
    Ymg = cat_main_gintnorm(Ymg,Tth); 
    if ~debug, clear Ylab Ysrc; end
   
    
    if ~debug, clear Yngm Ynwm Yncm; end 

    cat_io_cmd(' ','','',debug,stime2);  
  end
  cat_io_cmd('','','',verb,stime);


end
function Ygi = cat_vol_noPVE(Ygi,vx_vol,mres,dist,iter)
  if ~exist('dist','var'), dist=2; end
  if ~exist('iter','var'), iter=1; end
  if ~exist('mres','var'), mres=1; end
  
  if any(vx_vol<0.5)
    YM = Ygi>0;
    [Ygi,resT2] = cat_vol_resize(Ygi,'reduceV',vx_vol,mres,32,'meanm');
  end
  for i=1:iter
    Ygi      = cat_vol_median3(Ygi,Ygi>0,Ygi>0);
    [Ygistd,Ygimean] = cat_vol_localstat(Ygi,Ygi>0,dist,4); 
    Ygx      = Ygi<(Ygimean-Ygistd/4) | Ygi>(Ygimean+Ygistd/2); 
    Ygi(Ygx) = Ygimean(Ygx);
    clear Ygx;
  end
  if any(vx_vol<0.5)
    Ygi = cat_vol_approx(Ygi,'nh',resT2.vx_volr,4); 
    Ygi = cat_vol_resize(Ygi,'dereduceV',resT2);
    Ygi(~YM)=0;
  end
end
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
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id$

  % set this variable to 1 for simpler debuging without reduceBrain
  % function (that normally save half of processing time)
  %debug   = extopts.debug; % debug = 1;
  verb    = extopts.verb-1;
  dsize   = size(Ysrc);
  NS      = @(Ys,s) Ys==s | Ys==s+1;          % function to ignore brain hemisphere coding
  LASstr  = max(eps,min(1,extopts.LASstr));   % LAS strenght (for GM/WM threshold)3
  LAB     = extopts.LAB;                      % atlas labels
  mvx     = mean(vx_vol);                     % mean voxel volume to correct for morphological operations  
  T3th    = Tth.T3th(3:6);                    % CSF, GM, and WM peak
  Tth.T3th(1) = min(0,Tth.T3th(1));           % correction of the background value
  
  % set debug = 1 and do not clear temporary variables if there is a breakpoint in this file 
  dbs   = dbstatus; debug = 0; 
  for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,'cat_main_LASs'); debug = 1; break; end; end
  if debug, Ymo = Ym; Ybo = Yb; end 
  
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

  
  % help maps to detect edges (Yg) and sulci/gyris (Ydiv)
  Yg    = cat_vol_grad(Ym,vx_vol);                                                  % mean gradient map  ...  ./max(0.1,Ym)
  Ydiv  = cat_vol_div(Ym,vx_vol);                                                   % divergence map 
  Yp0   = single(Ycls{1})/255*2 + single(Ycls{2})/255*3 + single(Ycls{3})/255;      % tissue label map
 
  
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
  Yl1  = cat_vol_ctype(round(spm_sample_vol(Vl1A,double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0)));
  Yl1  = reshape(Yl1,dsize);
  if ~debug, clear Yy; end

  
  
  %% initial bias correction 
  %  ----------------------------------------------------------------------
  %  required especial in case of strong bias (7 Tesla) for the redefinition 
  %  of tissue segments for the LAS correction
  %  ----------------------------------------------------------------------
  stime = cat_io_cmd('  Initial bias correction','g5','',verb,stime); 
  if debug, Ym = Ymo; Yb = Ybo; end 
  Yb   = smooth3(Yb | (cat_vol_morph(Yb,'d',2/mvx) & Ym<0.8 & Yg<0.3 & Ym>0 & Yp0>0.2))>0.5;  % increase brain mask, for missing GM 

  %  background and brain distance map
  Ybg = cat_vol_smooth3X(Ycls{6}>128 & Yg<0.1 & Ym<1/6 & Ysrc<T3th(1)/3,4/mvx)>0.9 & Ysrc>mean(Tth.T3th(1:2));  % background
  [Ybr,Ybgr,resT3] = cat_vol_resize({single(Yb),single(Ybg)},'reduceV',vx_vol,2,32,'meanm');                 % brain distance Ybd
  Ybd = cat_vbdist(single(Ybr>0.1),Ybgr<=0.05,resT3.vx_volr); 
  Ybd = max(eps,cat_vol_resize(Ybd,'dereduceV',resT3)); clear Ybgr Ybr; 
  Ybv = (Ym - 2*Ydiv)>1.1 & Ym>1.1 & Ycls{2}<64 & Yp0>0;                                                     % blood vessels
  
  % brain tissues WM, GM, CSF (Ywi, Ygi, Yci)
  stime = cat_io_cmd('    Brain Tissues','g5','',debug,stime);  
  Ywmh = smooth3(Yp0/3-Ym)>0.2 & Ym<2.5/3 & Ym>1.5/3 & NS(Yl1,LAB.CT)>0.5;
  Ywh = cat_vol_smooth3X(max(0,Ym-1) .* (Yp0>2.75 & Ym>1.05 & ~Ybv & ~Ywmh),4)*4;                            % biased WM
  Ym  = Ym - Ywh; clear Ywh; % hard but good working correction for uncorrected high intensity bias (WM with intensity over 1)
  Yss = (NS(Yl1,LAB.BG) | NS(Yl1,LAB.TH)) & Ym<0.95 & Ym>0.6;                                                % subcortical structures
  % WM 
  Ywi = Ysrc .* (Yp0>2.5 | ((Ym - 2*Ydiv)>0.8 & (Yp0 - 6*Ydiv)>2)) .* (~Yss & ~Ywmh & ~Ybv & Yb);            % WM    
  Ywi = cat_vol_noPVE(Ywi,1);       
  % GM 
  Ygi = (Yp0 - 8*abs(Ydiv./Ym))>1 & (Yp0 + 4*abs(Ydiv./Ym))<2.5 & Yb;                  % GM 
  Ygi = (Ygi | Yss) & ~Ywi & ~Ywmh & ~cat_vol_morph( NS(Yl1,LAB.VT) | NS(Yl1,LAB.BS) ,'d',2/mvx);            % no ventricle/brainstem
  Ygi = Ygi .* Ysrc * T3th(3)/T3th(2);  
  Ygi = cat_vol_noPVE(Ygi);                                                                                  % correct for PVE
  % CSF
  Yci = (Ycls{3}>16 & NS(Yl1,LAB.VT))  | (Ycls{3}>240 & Ym>0.9/3 );                                          % CSF
  Yci = Yci & (Ym - Ydiv)<1.1/3 & ~Ygi & ~Ywi & (Yg./Ym)<0.1 & Ydiv<0.1 & ~Ywmh & Ysrc>mean(Tth.T3th(1:2));        % only save CSF!
  Yci(smooth3(Yci)<0.5) = 0; 
  Yci = Yci .* Ysrc; Yci = cat_vol_noPVE(Yci);                                                               % correct for WM intensity
  Ygi = Ygi + Yci * T3th(3)/median(Ysrc(Yci(:)>0));                                                          % save memory
  if ~debug, clear Yss Ybv Yci; end

  % head tissues 
  stime = cat_io_cmd('    Head Tissues','g5','',debug,stime);
  Yhg = Yp0==0 & Yg<0.2 & Ym>max(0.1,min(0.2,mean([0.2,median(Ym(Ybg(:)))]))) & (Ym + 2*abs(Ydiv))<0.9 & abs(Ydiv./Ym)<0.1 & Ybd>5;            % head GM like (mussels)
  Yhg(smooth3(Yhg)<0.5) = 0;                                                                                 % remove small dots
  Yhi = (Yp0==0 & ~Yhg & (Ym - 2*Ydiv)>0.8 & smooth3(Yg)>0.2 & ~Yhg & ~Ybg & Ybd>5 & (Ym>1.5 | Ym>max(0,(Ybd-5)/10))) & smooth3(Ydiv<-0.3)==0;    % higher head tissue 
  Yhi(smooth3(Yhi)<0.5) = 0;                                                                                 % remove small dots
  Yhi(cat_vol_smooth3X(Yhi,1/mvx)<0.5)=0;
  Yhi = Yhi .* Ysrc; Yhi = cat_vol_noPVE(Yhi);                                                               % correct for PVE
  Yhi = Yhi * T3th(3)/median(Yhi(Yhi(:)>0));  
  % correct for WM intensity
  Yhg = Yhg | ((smooth3(Yhg | Ycls{6}>128)<0.2) & ~Ygi & Ycls{2}<8 & (Ym - 2*Ydiv)>0.5 & Ym<0.8 & Yp0<2);    % further low intensity head tissue 
  Yhg = cat_vol_morph(Yhg & ~Yhi & Ybd>5 & Yg<0.15 & abs(Ydiv)<0.15,'o') .* Ysrc; 
  Ygi = cat_vol_noPVE(Ygi,2);
  Yhg = Yhg * T3th(3)/T3th(2); %median(Yhg(Yhg(:)>0 & Ybd(:)<8));                                            % expect GM like intensity
  Ygi = Ygi + Yhg;                                                                                     % save memory
  if ~debug, clear Yhi Yhg; end  
  Ygi = cat_vol_noPVE(Ygi,2);

  %% background to avoid conflict by bad approximation 
  stime = cat_io_cmd('    Background','g5','',debug,stime);
  if median(Ysrc(Ybg(:)))/T3th(3)>0.02 && median(Ysrc(Ybg(:)))/T3th(3)<T3th(1)
    Ybgi = (Ybg & Ym>0.01) .* (Ysrc * T3th(3)/max(eps,median(Ysrc(Ybg(:))))); 
    [Ybgi,Yd,resT2] = cat_vol_resize({Ybgi,single(Ygi>0)},'reduceV',vx_vol,min(vx_vol*3,6),32,'meanm'); 
    Yd   = cat_vbdist(Yd,Yd<1,resT2.vx_volr); 
    Ybgi = cat_vol_noPVE(Ybgi);
    Ybgi = cat_vol_approx(Ybgi,'nh',resT2.vx_volr,2); Ybgi = cat_vol_smooth3X(Ybgi,2); 
    Ybgi = max(eps,cat_vol_resize(Ybgi,'dereduceV',resT2));
    Yd   = max(eps,cat_vol_resize(Yd,'dereduceV',resT2));
    Ygi = Ygi + Ybgi .* (Yd>30); clear Ybgi Yd;
  else
    Ygi = Ygi + Ybg * T3th(3); 
  end
  
  %% estimate correction on lower resolution
  stime = cat_io_cmd('    Bias Field','g5','',debug,stime);
  [Ywi,resT2] = cat_vol_resize(Ywi,'reduceV',vx_vol,min(vx_vol*2,2),32,'max');        % maximum reduction for the WM
  [Ygi,resT2] = cat_vol_resize(Ygi,'reduceV',vx_vol,min(vx_vol*2,2),32,'meanm');      % mean for the rest
  [Ybr,resT2] = cat_vol_resize(Yp0,'reduceV',vx_vol,min(vx_vol*2,2),32,'meanm');      % mean for the rest
  Ywi = cat_vol_localstat(Ywi,Ywi>0,1,3);                                             % maximum for stabilization of small WM structures
  Ywi(Ywi==0 & Ygi>0)=Ygi(Ywi==0 & Ygi>0); clear Ygir;                                % mixing of both maps
  
  %% Ywi = cat_vol_noPVE(Ywi);   
  Yis = cat_vol_approx(Ywi,'nh',resT2.vx_volr,4); Yis = cat_vol_smooth3X(Yis,2);      % first bias field map 
  Ywi(Ywi>Yis*1.1 | (Ywi>Yis*1.01 & ~Ybr) | (Ywi<Yis*0.5 & Ywi>0)) = 0; 
  Yd  = cat_vbdist(single(Ybr>1),Ybr<=1,resT2.vx_volr);                                 % distance map the brain tissue 
  Ywi(Yd>8) = Yis(Yd>8); clear Yis Yd Ybr;                                            % use the first bias field far away from the brain
  Ywi = cat_vol_approx(Ywi,'nh',resT2.vx_volr,2); Ywi = cat_vol_smooth3X(Ywi,2);      % ... to use less filtering here
  Ywi = max(eps,cat_vol_resize(Ywi,'dereduceV',resT2));
  Ybf = Ywi; 
  Ym  = (Ysrc./max(eps,Ywi))/mean((Ysrc(Yp0(:)>2.9 & ~Ywmh(:))./max(eps,Ywi(Yp0(:)>2.9 & ~Ywmh(:)))))*T3th(3); 
  Ym  = cat_main_gintnorm(Ym,Tth); 
  if ~debug, clear Ywi; end
  % hard but good working correction for uncorrected high intensity bias (WM with intensity over 1)
  Ybv = (Ym - 2*Ydiv)>1.05 & Ym>1.05 & Ycls{2}<192;              
  Ywh = cat_vol_smooth3X(max(0,Ym-1) .* (Yp0>2.75 & Ym>1.05 & ~Ybv),4)*4;
  Ym = Ym - Ywh; 
  clear Ybv Ywh;
  cat_io_cmd('','','',debug,stime);   
  
  
  
  %% brain segmentation can be restricted to the brain to save time 
  stime = cat_io_cmd('  Prepare partitions','g5','',verb,stime); 
  [Ymr,Yb,BB] = cat_vol_resize({Ym,Yb},'reduceBrain',vx_vol,round(10/mvx),Yb);
  [Ygr,Ydivr,Yp0]   = cat_vol_resize({Yg,Ydiv,Yp0},'reduceBrain',vx_vol,BB.BB);
  Yl1               = cat_vol_resize(Yl1          ,'reduceBrain',vx_vol,round(4/mvx),BB.BB);
  Yclsr = cell(size(Ycls)); for i=1:6, Yclsr{i} = cat_vol_resize(Ycls{i},'reduceBrain',vx_vol,BB.BB); end
  %Ywtpm = cat_vol_resize(Ywtpm,'reduceBrain',vx_vol,round(4/mvx),BB.BB);
  
  
  % adaption of the LASstr depending on average basal values 
  LASmod  = min(2,max(0,mean((Ymr( NS(Yl1,LAB.BG) & Ygr<0.1 & Ydivr>-0.05  & Yclsr{1}>4)) - 2/3) * 8));
  LASstr  = min(1,max(0.05,LASstr * LASmod)); clear LASmod                 % adaption by local BG variation
  LASfs   = 1 / max(0.05,LASstr);                                          % smoothing filter strength 
  LASi    = min(8,round(LASfs));                                           % smoothing interation (limited)
  
  
  %  GM thickness (Ygmt) and percentage possition map (Ypp) estimation
  %  The Ypp and Ygmt maps are used to refine the GM especially to correct
  %  highly myelinated GM regions. Using allow to avoid overcorrections.
  fastppi = 1; % 1 mm: 15 vs. 40 seconds ... it thing fast is ok here 20161014 
  if fastppi
    [Ygmt,Ypp,Ywmd,Ycsfd] = cat_vol_pbt( (Yp0 + (Ymr*3 .* (Yp0>0)))/2 , ...
      struct('resV',vx_vol,'verb',0,'dmethod','vbdist','method','pbt2x') );
  else
    [Ygmt,Ypp,Ywmd,Ycsfd] = cat_vol_pbt( (Yp0 + (Ymr*3 .* (Yp0>0)))/2 , ...
      struct('resV',vx_vol,'verb',0) );
  end
  Ygmtroi = Ygmt>0 & Ygmt<6 & NS(Yl1,LAB.CT); 
  GMTstat = [cat_stat_nanmedian(Ygmt(Ygmtroi(:))) cat_stat_nanstd(Ygmt(Ygmtroi(:)))]; 
  [D,I]   = cat_vbdist(single(Ygmt>0.1),Yp0>0,vx_vol); Ygmt = Ygmt(I); clear D I;  % full GMT map for the whole brain
  if ~debug, clear Ywmd Ygmtroi; end
  
  
  
  %% As far as SPM segmentation is not optimal we need some refinements.
  %  ----------------------------------------------------------------------
  stime = cat_io_cmd('  Prepare segments','g5','',verb,stime);
 
  % use SPM segment 
  Ycm = Yp0>0.5 & Yp0<=1.5;
  Ygm = Yp0>1.5 & Yp0<=2.5;
  Ywm = Yp0>2.5;
  
  % corrections for original intensies 
  Ywm = Ywm & Ymr>2/3 & Ymr<1.3; 
  Ygm = Ygm & Ymr>1/3 & Ymr<1.0; 
  Ycm = Ycm & Ymr>0.1 & Ymr<1.5/3;
   
  % add subcortical structures
  Yss = (NS(Yl1,LAB.BG) | NS(Yl1,LAB.TH)) & Ymr>2/3 & Ymr<2.75/3 & Ygr<0.1 & Ydivr>-0.05; 
  Yss = smooth3(Yss)>0.5;
  Ygm = Ygm | Yss;
  Ywm = Ywm & ~Yss;
 
  % blood vessel correction
  Ybv = (Ymr-Ydivr)>1.2 & Ymr>0.9 & Yclsr{2}<192;
  Ywm = Ywm & ~Ybv; 
  Ygm = Ygm & ~Ybv; 
  Ycm = Ycm & ~Ybv; 
  
  % correction by divergence to avoid gyri in the GM segment 
  Ymlwm = Ygm & (Ymr - 2*Ydivr + (Ygmt-3)/10 + 0.1*NS(Yl1,LAB.CB))>0.9 & ...
          ~Yss & Ymr>(2.25 + min(0.25,max(0,mvx)))/3  & ...
          cat_vol_morph((Ymr - 2*Ydivr + 0.2*NS(Yl1,LAB.CB))>1.75/3,'e',1/mvx);  
  Ymlwm = Ymlwm | (smooth3(Yp0>1.5 & Ymr>2.5/3 & Ymr<3.2/3 & (Ydivr<0 | Ygr<0.1) & ~Ybv)>0.3 & Ymr>2.5/3);
  Ymlwm = (Ymlwm & (Ymr - Ydivr)>0.9) | (cat_vol_morph(Ymlwm,'e') & (Ymr - Ydivr)>0.8); 
  Ygm   = Ygm & ~Ymlwm; 
  Ywm   = Ywm | Ymlwm; 
  if ~debug, clear Ymlwm; end
  
  % added undetected GM around WM 
  Ymlgm = ~Ywm & Yb & cat_vbdist(single(Ywm),Yb,vx_vol)<3 & ... 
    Ymr>(1.25 + min(0.25,max(0,mvx)))/3 & Ymr<2.5/3 &  ...
    abs(Ydivr)<0.1 & cat_vol_morph(NS(Yl1,LAB.CT),'d',2/mvx); 
  Ygm   = Ygm | Ymlgm; 
  Ycm   = Ycm & ~Ymlgm;
  if ~debug, clear Ymlgm; end
  
  % correct for uncorrected WM bias
  Ywh  = cat_vol_smooth3X(Ywm & Ymr>1.05 & ~Ybv,3)>0.05;
  Ygmc = Ywm & (Ywh & Ymr<0.98); 
  Ycmc = Ygm & (Ywh & Ymr<0.75); 
  if ~debug, clear Ywh; end
  Ygm  = (Ygm & ~Ycmc) | Ygmc; 
  Ywm  = Ywm & ~Ygmc; 
  Ycm  = Ycm | Ycmc;
  if ~debug, clear Ygmc Ycmc; end
  
  % correction by divergence to avoid sulci in the GM segment
  Ypvec = Ygm & (Ymr - 2*Ydivr - (Ygmt)/10 - Ypp/10)<0.5 & Ymr<1.25/3;  
  Ygm   = Ygm & ~Ypvec; 
  Ycm   = Ycm | Ypvec; 
  if ~debug, clear Ypvec; end
  
  % correction for WMHs
  Yvt  = cat_vol_morph( NS(Yl1,LAB.CB) | NS(Yl1,LAB.BS) | NS(Yl1,LAB.PH) | ...
                        NS(Yl1,LAB.BG) | NS(Yl1,LAB.TH) | NS(Yl1,LAB.HC)  ,'d',5/mvx); 
  Ywmh = smooth3( (Yp0/3-Ymr)>0.1 & Ymr<2.5/3 & Ymr>1.5/3 & ~Yss & ~Yvt)>0.5;
  Ywm  = Ywm & ~Ywmh; 
  Ygm  = Ygm & ~Ywmh; 
  if ~debug, clear Yss; end
  
  % myelinated GM  -  this is the major goal of LAS !!!
  Ymgm = ~(cat_vbdist(single(Ypp<0.1),Yb,vx_vol)>(min(2.5,max(1,GMTstat(1)+GMTstat(2))))) & ~Yvt & ...
            cat_vol_smooth3X(Ygmt,4)<=min(2,max(1,GMTstat(1)-GMTstat(2)/2)); 
  Ymgm = Ymgm | (Yclsr{1}>0 & Ycsfd<min(2,max(1,GMTstat(1)-GMTstat(2)/2)) & ~Ycm & ~Yvt);
  if ~debug, clear Ycsfd; end
  Ymgm = Ymgm & Yp0>2.0 & Ymr<0.95 & Ymr>0.8 & Ymr<1 & Ywm & ~Ybv & ~Ywmh & Ymr>2/3 & ~Yvt & ~Ycm & NS(Yl1,LAB.CT); 
  Yhd  = cat_vbdist(single(~Yb),Yb,vx_vol);
  Ymgm = Ymgm | (Yp0>2.0 & Ymr>0.5 & Yhd<min(2,max(1,GMTstat(1)-GMTstat(2)/2)) & ~Yvt & ~Ybv);
  if ~debug, clear Yvt Yhd Ywmh Ybv; end
  Ygm  = Ygm | Ymgm; 
  Ywm  = Ywm & ~Ymgm;
  Ycm  = Ycm & ~Ymgm;
  if ~debug, clear Ymgm; end
  
  % parahippocampale gyrus and hippocampus
  Yhcg  = (NS(Yl1,LAB.PH) | NS(Yl1,LAB.HC)) & (Ymr - Ydivr)>2.5/3;
  Ywm  = Ywm | Yhcg; 
  Ygm  = Ygm & ~Yhcg;
  if ~debug, clear Yhcg; end
  
  % Ycsf
  Yngm = Ygm & Ypp<0.2 & Ymr<0.5 & Ygmt<=GMTstat(1) & Ypp<0.1; 
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
  Ycm = Ycm & (Ymr-Ydivr)<.35 & Ymr>0.01 & Yclsr{3}>240; 
  Ycm(smooth3(Ycm)<0.5)=0;  
  Ygm = Ygm | (~Ycm & (Ymr-Ydivr)>.4 & Ygr<0.2 & Ymr<0.5 & Yb); 
  Ygm(smooth3(Ygm | Ywm)<0.5)=0;  
  if ~debug, clear Ymr Ydivr Ygr Yclsr; end 
  
  
  %% back to original resolution for full bias field estimation
  [Ycm,Ygm,Ywm] = cat_vol_resize({Ycm,Ygm,Ywm},'dereduceBrain',BB); 
  [Yp0,Yl1]     = cat_vol_resize({Yp0,Yl1},'dereduceBrain',BB);
  [Yb]          = cat_vol_resize({Yb},'dereduceBrain',BB);
  

  %% new head tissues
  Ynb  = smooth3(Ycls{6})>128 & Ysrc<mean(Tth.T3th(2:3)) & Ym<1/6; % ... zero background
  Yhdh = Yp0==0 & (Ym - Ydiv./Ym)>0.5 & (Ym - Yg)>0.8 & (Ym>1.2 | smooth3(Yg)>0.1) & ~Ybg & Ybd>2 & (Ym>1.2 | Ym>max(0.4,(1-Ybd/300)));
  Yhdm = Yp0==0 & Ym<1.0 & Ym>0.5; 
  Yhdm = Yhdm | (cat_vol_morph(Yp0>1.5,'d',3) & Yp0<1.5 & cat_vol_smooth3X(Yhdh>0.5,8)>0.1 & ~Yhdh & Ym>0.5 & Ym<0.9);
  Yhdm = Yhdm | ((smooth3(Yhdh | Ycls{6}>128)<0.2) & ~Ygm & ~Ywm & ~Ycm & (Ym - 2*Ydiv)>0.5 & Ym>0.6);   
  Yhdm = Yhdm & Yg<0.1 & abs(Ydiv)<0.2 & ~Ybg & ~Yhdh;
  Yhdm((smooth3(Yhdm) - smooth3(Yhdh)/1.5)<0.6) = 0;
  if ~debug, clear Ybg Ydiv Yg; end 

  
%% --------------------------------------------------------------------- 
%  Now, we can estimate the local peaks 
%  ---------------------------------------------------------------------
  % Estimation of the local WM threshold with "corrected" GM voxels to
  % avoid overfitting (see BWP cerebellum). 
  % CSF is problematic in high contrast or skull-stripped image should 
  % not be used here, or in GM peak estimation
  mres  = 1.1; 
  stime = cat_io_cmd('  Estimate local tissue thresholds','g5','',verb,stime); 
  Ysrcm = cat_vol_median3(Ysrc.*Ywm,Ywm,Ywm); 
  rf    = [10^5 10^4];
  T3th3 = max(1,min(10^6,rf(2) / (round(T3th(3)*rf(1))/rf(1))));
  Ysrcm = round(Ysrcm*T3th3)/T3th3;
  % major correction outside the brain 
  Ygi = cat_vol_noPVE(Ysrc .* Ygm  * T3th(3)/T3th(2)) + ...
        ... cat_vol_noPVE(Ysrc .* Ycm  * T3th(3)/T3th(1)) + ...
        cat_vol_noPVE(Ysrc .* Yhdm * T3th(3)/T3th(2)) + ...
        cat_vol_noPVE(Ysrc .* Yhdh * T3th(3)/median(Ysrc(Yhdh(:)>0 & Ybd(:)<15))); 
  Ysrcm(Ybd>20) = Ybf(Ybd>20); 
  % fine correction
  [Yi ,resT2] = cat_vol_resize(Ysrcm,'reduceV',vx_vol,mres,32,'max'); % maximum reduction for the WM
  if ~debug, clear Ysrcm; end
  Yi = cat_vol_localstat(Yi,Yi>0,1,3); % one maximum for stabilization of small WM structures
  Yi(Yi==0 & Ygi>0)=Ygi(Yi==0 & Ygi>0); if ~debug, clear Ygi; end
  Yi = cat_vol_approx(Yi,'nh',resT2.vx_volr,2); Yi = cat_vol_smooth3X(Yi,LASfs); 
  Ylab{2} = max(eps,cat_vol_resize(Yi,'dereduceV',resT2)); 
  if ~debug, clear Yi; end

  
  %% GM
  Yi = Ysrc ./ max(eps,Ylab{2}) .* (Ygm | Yhdm);
  Ybs = NS(Yl1,LAB.BS) & Ym<1.1 & Ym>0.9 & Yp0>2.5 & Ym>0.9;
  Yi(Ybs)  = Ysrc(Ybs) ./ max(eps,Ylab{2}(Ybs))   .* T3th(2)/T3th(3); clear Ybs;
  Yi = cat_vol_noPVE(Yi);
  Yi(Ybd>20 & Yi==0) = 2/3; 
  if ~debug, clear Yhdm; end
  [Yir,Ygmr,resT2] = cat_vol_resize({Yi,Ygm},'reduceV',vx_vol,mres,32,'meanm'); if ~debug, clear Yi; end
  Yir = cat_vol_noPVE(Yir);
  Yir = cat_vol_approx(Yir,'nh',resT2.vx_volr,2); 
  Yir = cat_vol_smooth3X(Yir,LASfs); 
  Ylab{1} = cat_vol_resize(Yir,'dereduceV',resT2).*Ylab{2};   
  if ~debug, clear Yir Ybd Yl1; end
  
  
  %% CSF & BG 
  Ynb = Ynb | smooth3( (Ycls{4}>128 | Ycls{6}>128) & Ym<median(Ym(Ycls{4} & Ym<0.4)))>0.5; if ~debug, clear Ym; end
  Ynb = cat_vol_morph(Ynb,'e',2/mvx);
  Ynb = Ynb & Ysrc~=0; 
  [Yc,resT2] = cat_vol_resize(round(Ysrc ./ max(eps,Ylab{2}) .* (smooth3(Ycm)>0.5) * rf(2))/rf(2),...
     'reduceV',vx_vol,8,16,'min');% only pure CSF !!!
  [Ynb,resT2] = cat_vol_resize(round(Ysrc ./ max(eps,Ylab{2}) .* Ynb * rf(2))/rf(2),...
     'reduceV',vx_vol,8,16,'meanm');
  Ynb(Yc>0)=0; Yc(Ynb>0)=0;
  for xi=1:2*LASi, Ynb = cat_vol_localstat(Ynb,Ynb>0,2,1); end
  for xi=1:2*LASi, Yc  = cat_vol_localstat(Yc,Yc>0,2,1); end
  Ynba = cat_vol_approx(Ynb ,'nh',resT2.vx_volr,2); 
  Yca  = cat_vol_approx(Yc ,'nh',resT2.vx_volr,2); % + min(max( meanYnb + stdYbc , meanYc - stdYbc ),...
  clear Yc Ynb; 
  
  %Yca  = max(Yca,Ynba*1.5); 
  Ynba  = min(Yca/1.5,Ynba); 
  Yca   = max(Yca,Ynba*1.5); 
  
  Ynba = cat_vol_smooth3X(Ynba,LASfs*2); 
  Yca  = cat_vol_smooth3X(Yca,LASfs*2); % * meanYc/mean(Yca(:)); 
  Ylab{3} = cat_vol_smooth3X(cat_vol_resize(Yca,'dereduceV',resT2).*Ylab{2},LASfs*2);  
  Ylab{6} = cat_vol_smooth3X(cat_vol_resize(Ynba,'dereduceV',resT2).*Ylab{2},LASfs*2);
  clear Ynba Yca; 
  
  
  %% local intensity modification of the original image
  % --------------------------------------------------------------------
  Yml = zeros(size(Ysrc));  
  Yml = Yml + ( (Ysrc>=Ylab{2}                ) .* (3.0 + (Ysrc-Ylab{2}) ./ max(eps,Ylab{2}-Ylab{1})    ) );
  Yml = Yml + ( (Ysrc>=Ylab{1} & Ysrc<Ylab{2} ) .* (2.0 + (Ysrc-Ylab{1}) ./ max(eps,Ylab{2}-Ylab{1})    ) );
  Yml = Yml + ( (Ysrc>=Ylab{3} & Ysrc<Ylab{1} ) .* (1.0 + (Ysrc-Ylab{3}) ./ max(eps,Ylab{1}-Ylab{3})    ) );
  Yml = Yml + ( (Ysrc>=Ylab{6} & Ysrc<Ylab{3} ) .* (1/3 + (Ysrc-Ylab{6}) ./ max(eps,Ylab{3}-Ylab{6})*2/3) );
  Yml = Yml + ( (Ysrc< Ylab{6}                ) .* (0.0 + (Ysrc        ) ./ max(eps,Ylab{6}*3          )) );
  Yml(isnan(Yml) | Yml<0)=0; Yml(Yml>10)=10;
  Yml = Yml/3;
  
  %% global
  Ymg = Ysrc ./ max(eps,Ylab{2}) * Tth.T3th(5);
  Ymg = cat_main_gintnorm(Ymg,Tth); 
  if ~debug, clear Ylab Ysrc; end
  
  % fill up CSF in the case of a skull stripped image 
  if max(res.mn(res.lkp==5 & res.mg'>0.1)) < mean(res.mn(res.lkp==3 & res.mg'>0.3))
    YM   = cat_vol_morph(Yb,'d'); 
    Ymls = smooth3(max(Yml,YM*0.5));
    Yml(YM & Yml<0.5)=Ymls(YM & Yml<0.5); 
    clear Ymls YM
  end
  
  
  %% prepare class correction 
  Yp0n = Yml .* Yb .* (Yml<1.1);
  Yp0n = min(Yp0n,2-(smooth3(Ywm)>0.3));
  Yp0n(Yp0n>1.1)=0;
  Yp0n(Yp0<1.5 & Yml<1.5/3 & Yp0>1/3) = 1/3; 
  Yp0n(smooth3(Yp0n>1.5/3)<0.5 & Yp0n>1/3) = 1/3;
  Yp0n(cat_vol_morph(Yp0n>1.25/3,'labopen',1/mvx)==0 & Yp0n>1/3)=1/3;
  Yp0n(smooth3((Yml<1/6 | Yml>0.45) & Yp0n<0.34 & Yp0<1.5)>0.5)=0;
  Yp0n(cat_vol_morph(Yp0n>0.5/3,'labopen',1/mvx)==0 & Yp0n<=0.34)=0;
  Yp0n(cat_vol_morph(Yp0n>0.5/3,'labclose',2/mvx) & Yp0n<=0.34)=1/3;
  if ~debug, clear Yp0; end
  
  %% update classes
  Yp0toC  = @(Yp0,c) 1-min(1,abs(Yp0-c));
  Ycls{1} = cat_vol_ctype(Yp0toC(Yp0n,2)*255,'uint8');
  Ycls{2} = cat_vol_ctype(Yp0toC(Yp0n,3)*255,'uint8');
  Ycls{3} = cat_vol_ctype(Yp0toC(Yp0n,1)*255,'uint8');
  Ycls{6} = cat_vol_ctype(cat_vol_morph(smooth3(1-max(0,Ymg*3) - Yp0n)>0.75,'lc')*255,'uint8'); 
  Ycls{5} = cat_vol_ctype((Ymg<4/6 & ~Yp0n & ~Ycls{6})*255,'uint8'); 
  Ycls{4} = cat_vol_ctype((Ymg>5/6 & ~Yp0n & ~Ycls{6})*255,'uint8'); 
  cat_io_cmd('','','',verb,stime);


end
function Ygi = cat_vol_noPVE(Ygi,dist,iter)
  if ~exist('dist','var'), dist=2; end
  if ~exist('iter','var'), iter=1; end

  for i=1:iter
    Ygi      = cat_vol_median3(Ygi,Ygi>0,Ygi>0);
    Ygimean  = cat_vol_localstat(Ygi,Ygi>0,dist,1); 
    Ygistd   = cat_vol_localstat(Ygi,Ygi>0,dist,4); 
    Ygx      = Ygi<(Ygimean-Ygistd/4) | Ygi>(Ygimean+Ygistd/2); 
    Ygi(Ygx) = Ygimean(Ygx);   
  end
end

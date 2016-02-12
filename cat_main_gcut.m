function [Yb,Yl1] = cat_main_gcut(Ysrc,Yb,Ycls,Yl1,YMF,vx_vol,opt)
% This is an exclusive subfunction of cat_main.
% ______________________________________________________________________
%
% gcut+: skull-stripping using graph-cut
% ----------------------------------------------------------------------
% This routine use morphological, region-growing and graph-cut methods. 
% It starts from the WM segment and grow for lower tissue intensities.
% Atlas knowledge is used to for separate handling of the cerebellum.
% Because its high frequency structures in combination with strong
% noise or other artifacts often lead to strong underestimations.
%
% There are 4 major parameters:
%   gcutstr - strengh of skull-stripping with str=0 - more tissue, str=1 less tissue
%   vx_res - resolution limit for skull-stripping (default 1.5)
%   gcutCSF 
% Especialy the second parameter controls many subparameters for  
% different tissue thresholds, region-growing, closing and smoothing
% parameters.
% This routine have a strong relation to the previous estimated main
% partition map l1, and the blood vessel correction. Therefore, it is
% maybe useful to move it...
%
%   [Yb,Yl1] = cat_main_gcut(Ysrc,Yb,Ycls,Yl1,YMF,vx_vol,opt)
% 
%   Yb   .. updated brain mask
%   Yl1  .. updated label map
% 
%   Ysrc .. anatomical image
%   Yb   .. initial brain mask
%   Ycls .. SPM tissue classification
%   Yl1  .. CAT atlas map
%   YMF  .. subcortical/ventricular regions (for filling in surf. recon.)
%   vx_vol .. image resolutino
%   opt  .. further options
% 
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id$

  LAB.CB =  3; % Cerebellum
  LAB.VT = 15; % Ventricle
  LAB.HI = 23; % WM hyperintensities
  
  if ~exist('opt','var'), opt=struct(); end
  def.debug   = cat_get_defaults('extopts.debug');
  def.gcutstr = cat_get_defaults('extopts.gcutstr');
  def.verb    = cat_get_defaults('extopts.verb')-1;
  opt = cat_io_checkinopt(opt,def); 
  opt.gcutstr = max(0,min(1,opt.gcutstr));
  %noise   = cat_stat_nanstd(Ym(cat_vol_morph(cat_vol_morph(Ym>0.95 & Ym<1.05,'lc',1),'e')));

  NS = @(Ys,s) Ys==s | Ys==s+1;
  
  voli = @(v) (v ./ (pi * 4./3)).^(1/3);           % volume > radius
  brad = double(voli(sum(Yb(:)>0).*prod(vx_vol))); % distance and volume based brain radius (brad)
  Yp0  = single(Ycls{3})/255/3 + single(Ycls{1})/255*2/3 + single(Ycls{2})/255;
  rvol = [sum(round(Yp0(:)*3)==1), sum(round(Yp0(:)*3)==2), sum(round(Yp0(:)*3)==3)]/sum(round(Yp0(:)*3)>0);
  
  %% set different paremeters to modifiy the stength of the skull-stripping 
  %gc.n = max(0.05,min(0.1,noise));
  % intensity parameter
  gc.h = 3.5  - 0.25*opt.gcutstr; % upper tissue intensity (WM vs. blood vessels)     - higher > more "tissue" (blood vessels)
  gc.l = 1.7  - 0.40*opt.gcutstr; % lower tissue intensity (WM vs. blood vessels)     - higher > more "tissue" (blood vessels)
  gc.o = 0.25 + 0.50*opt.gcutstr; % BG tissue intensity (for high contrast CSF=BG=0!) - lower value > more "tissue"
  % distance parameter
  gc.d = (brad*4 - brad*2*opt.gcutstr)/mean(vx_vol);          % distance  parameter for downcut - higher > more tissue
  gc.c = (0.005 - 0.02*opt.gcutstr)*mean(vx_vol);              % growing   parameter for downcut - higher > more tissue
  gc.f = (brad/10 * opt.gcutstr * rvol(1)/0.10)/mean(vx_vol); % closing   parameter             - higher > more tissue ... 8
  % smoothing parameter
  gc.s = 0.4  + 0.20*min(0.55,opt.gcutstr); % smoothing parameter                     - higher > less tissue
  
  if opt.verb, fprintf('\n'); end
  stime = cat_io_cmd('  WM initialisation','g5','',opt.verb); dispc=1;
  %% init: go to reduces resolution 
  [Ym,Yl1,YMF,BB] = cat_vol_resize({Ysrc,Yl1,YMF},'reduceBrain',vx_vol,round(4/mean(vx_vol)),Yb);
  [Ywm,Ygm,Ycsf,Ymg,Yb] = cat_vol_resize({single(Ycls{2})/255,single(Ycls{1})/255,...
    single(Ycls{3})/255,single(Ycls{5})/255,Yb},'reduceBrain',vx_vol,round(4/mean(vx_vol)),Yb);
  vxd  = max(1,1/mean(vx_vol)); 
  Ymg = Ymg>0.05 & Ym<0.45; 
  
  clear Ycls
  
  %% initial WM+ region
  YHDr = cat_vol_morph(Yl1>20 | Yl1<=0,'e',vxd*2);
  Yb  = Yb>0.25 & Ym>2.5/3 & Ym<gc.h/3 & Yl1<21 & Yb;  % init WM 
  Yb  = cat_vol_morph(Yb,'l'); 
  
  % if no largest object could be find it is very likeli that initial normalization failed
  if isempty(Yb)
          error('cat:cat_main:largestWM','No largest WM cluster could be found: Please try to set origin (AC) and run preprocessing again because it is very likeli that spatial normalization failed.');
  end
  
  Yb  = Yb | (Ym>2.5/3  & Ym<gc.h/3 & Yb) | NS(Yl1,LAB.HI) | NS(Yl1,LAB.VT);          % init further WM 
  Yb  = smooth3(Yb)>gc.s;
  Yb(smooth3(single(Yb))<0.5)=0;                          % remove small dots
  [Ybr,resT2] = cat_vol_resize(single(Yb),'reduceV',vx_vol,mean(vx_vol)*4,32); 
  Ybr = cat_vol_morph(Ybr>0,'lc',gc.f*mean(resT2.vx_vol)/mean(resT2.vx_volr));
  Ybr = cat_vol_resize(smooth3(Ybr),'dereduceV',resT2)>gc.s; 
  Yb = single(Yb | (Ym<1.1 & Ybr));
  %Yb  = single(cat_vol_morph(Yb,'labclose',gc.f));         % one WM object to remove vbs
  
  
  %% region growing GM/WM (here we have to get all WM gyris!)
  stime = cat_io_cmd('  GM region growing','g5','',opt.verb,stime); dispc=dispc+1;
  Yb(~Yb & (YHDr | Ym<1.9/3 | Ym>gc.h/3 | (Ywm + Ygm)<0.5))=nan; %clear Ywm Ygm; 
  [Yb1,YD] = cat_vol_downcut(Yb,Ym,0.01+gc.c); % this have to be not to small... 
  Yb(isnan(Yb) | YD>gc.d*vxd*2)=0; Yb(Yb1>0 & YD<gc.d*vxd*2)=1;
  Yb(smooth3(single(Yb))<gc.s)=0;
  Yb = single(Yb | (cat_vol_morph(Yb,'labclose',vxd) & Ym<1.1));

  
  %% region growing CSF/GM 
  stime = cat_io_cmd('  GM-CSF region growing','g5','',opt.verb,stime); dispc=dispc+1;
  Yb(~Yb & (YHDr | Ym<gc.l/3 | Ym>gc.h/3) | Ymg)=nan;
  [Yb1,YD] = cat_vol_downcut(Yb,Ym,0.00+gc.c);
  Yb(isnan(Yb) | YD>gc.d/2)=0; Yb(Yb1>0 & YD<gc.d)=1; 
  for i=1:2, Yb(smooth3(single(Yb))<gc.s)=0; end
  Yb  = single(cat_vol_morph(Yb,'o',max(1,min(3,4 - 0.2*gc.f* (rvol(1)/0.4) ))));
  Yb  = single(Yb | (cat_vol_morph(Yb ,'labclose',1) & Ym<1.1));
  
  %% region growing - add CSF
  Yb(~Yb & (YHDr | Ym<1/3 | Ym>gc.h/3) | Ymg)=nan;
  [Yb1,YD] = cat_vol_downcut(Yb,Ym,-0.02+gc.c);
  Yb(isnan(Yb) | YD>gc.d/2)=0; Yb(Yb1>0 & YD<gc.d)=1; 
  for i=1:2, Yb(smooth3(single(Yb))<gc.s)=0; end
  Yb  = single(Yb | (cat_vol_morph(Yb ,'labclose',1) & Ym<1.1));
  
  %% region growing - add CSF regions   
  stime = cat_io_cmd('  CSF region growing','g5','',opt.verb,stime); dispc=dispc+1;
  Ygr = cat_vol_grad(Ym,vx_vol);
  Yb(~Yb & smooth3(cat_vol_morph(smooth3(Ym<0.75/3 | (Ym>1.25/3 & ~Yb) | ...
    (Ygr>0.05 & ~Yb))>0.5,'lc',vxd*2) | Ymg )>0.5)=nan; 
  [Yb1,YD] = cat_vol_downcut(Yb,Ym,-0.02+gc.c); 
  Yb(isnan(Yb) | YD>gc.d/2)=0; Yb(Yb1>0 & YD<gc.d*2 & YD>0)=1;
  for i=1:2, Yb(cat_vol_smooth3X(Yb,2)<(gc.s - 0.25))=0; end
  Yb = Yb | YMF; 
  
  % smooth / low dilated boundary 
  Ybs = single(Yb); spm_smooth(Ybs,Ybs,4*gc.s./vx_vol); Yb   = Yb | (Ybs>(gc.s-0.25) & Ym<1.25/3);
  
  %% filling of ventricles and smooth mask
  stime = cat_io_cmd('  Ventricle closing','g5','',opt.verb,stime); dispc=dispc+1; %#ok<*NASGU>
  Yb  = Yb | (cat_vol_morph(Yb ,'labclose',vxd*gc.f) & ...
    Ym>=gc.o/3 & Ym<1.25/3 & ~Ymg & Ycsf>0.75);
  Yb  = single(cat_vol_morph(Yb,'o',max(1,min(3,4 - 0.2*gc.f* (rvol(1)/0.4) ))));
  Yb  = Yb | (cat_vol_morph(Yb ,'labclose',vxd) & Ym<1.1);
  Ybs = single(Yb); spm_smooth(Ybs,Ybs,0.6./vx_vol); Yb = max(Yb,Ybs)>gc.s;
 
  %%
  Yb   = cat_vol_resize(Yb  ,'dereduceBrain',BB)>0.5;
  Yl1  = cat_vol_resize(Yl1 ,'dereduceBrain',BB);
    
  %% update Yl1 with Yb
  Yl1(~Yb)  = 0;
  [tmp0,tmp1,Yl1] = cat_vbdist(single(Yl1),Yl1==0 & Yb); clear tmp0 tmp1;

  if opt.debug
    cat_io_cmd(' ','','',opt.verb,stime); 
  else
    cat_io_cmd(' ','','',opt.verb,stime);   
%    cat_io_cmd('cleanup',dispc,'',opt.verb);
  end

end
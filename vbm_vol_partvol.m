function [vol,Yl1,YB,YMF] = vbm_vol_partvol(Yl1A,Yp0,Ym,Yl0,opt)
% ______________________________________________________________________
% Use a segment map Yp0, the global intensity normalized T1 map Ym and 
% the atlas label map Yl1 to create a individual label map Yl1. 
% The atlas contain main regions like cerebrum, brainstem, midbrain,
% cerebellum, ventricle, and regions with blood vessels. 
%
% This function try to solve the following problems:
%  1) Finding of the cerebrum, the cerebellum, the head, blood vessels, 
%     brain skin and other mayor structures based on atlas (Yl1A) and 
%     tissue class information (Yp0). 
%     To do this it is important to use data from the T1-map (Ym) that
%     use the same intensity scalling as the segment map Yp0, but have 
%     more informations about partial volume regions.
%  2) Set Partions:
%     2.1) Find biggest WM part of each region.
%     2.2) Align the nearest region class for other voxel
%     2.3) Finding and Filling of the ventricle and the Basalganglia
%     2.4) Find blood vessels
% ______________________________________________________________________
%
%   Was ist neu im Vergleich zu anderen?
%   - Zuweisung durch Dartel mit hoher Genauigkeit
%   - Erweiterung von SPM/VBM durch MainROIs (Seiten, Lappen, ...)
%   - Verbesserung der SPM/VBM durch bessere Enfernung von unerwünschtem
%     Gewebe (ON, Blutgefäße ...)
%   - Blutgefäße könnnen als erweitere Masken für fMRI genutzt werden um
%     Seiteneffekte besser ausblenden zu können.
%  [- Beliebige Atlanten können genutzt werden.]
%
%  Todo:
%   - Entfernen von unsicheren Regionen
%   - stärkerer Abgleich mit Yp0
%   - stärkere Einbindung der LAS
% ______________________________________________________________________
%
% Structure:
%
%   [Yl1,pfT,mfT] = vbm_vol_partvol(Yl1A,Yp0,Ym,opt)
%
%   INPUT:  Yl1A = 3D-volume with brain regions (altas map)
%           Yp0 = 3D-volume with tissue propability map (CSF=1,GM=2;WM=3)
%           Ym = intensity normalized T1 image
%           opt
%            .res    = resolution for mapping
%            .vx_vol = Voxelsize
%            .LAB    = Label of Yl1 map (see def.LAB definition below)
%            
%
%   OUTPUT: Yl1 = individual label map 
%           pfT = Corrected and filled Yp0 map 
%           mfT = Corrected and filled Ym map
%
% ______________________________________________________________________
% Structural Brain Mapping Group, University Jena, Germany
% Robert Dahnke
% 2013/04
%
% $Id$


  if ~exist('opt','var'), opt=struct(); end
  def.res    = 2;
  def.vx_vol = [1 1 1];
  
  % definition of ROIs with: ID = [L R]
  def.LAB.CT = [ 1  2]; % cortex
  def.LAB.MB = [13 14]; % MidBrain
  def.LAB.BS = [13 14]; % BrainStem
  def.LAB.CB = [ 3  4]; % Cerebellum
  def.LAB.ON = [11 12]; % Optical Nerv
  def.LAB.BG = [ 5  6]; % BasalGanglia 
  def.LAB.TH = [ 9 10]; % Hypothalamus 
  def.LAB.HC = [19 20]; % Hippocampus 
  def.LAB.VT = [15 16]; % Ventricle
  def.LAB.NV = [17 18]; % no Ventricle
  def.LAB.BV = [ 7  8]; % Blood Vessels
  def.LAB.NB = [ 0  0]; % no brain 
  def.LAB.HD = [21 22]; % head
  
  opt = checkinopt(opt,def);
    
  vx_vol = opt.vx_vol;
    
  %% OPTIMIZATION:
  % ds('l2','',vx_vol,Ym,Yl1A,Ym,Yp0/3,130)
  Ym = Ym*3; mgTO = Ym;
  [Ym,Yp0,Yl1A,Yl0,BB] = vbm_vol_resize({Ym,Yp0,Yl1A,Yl0},'reduceBrain',vx_vol,10,Yp0>0);   % removing of background
  [Ym,Yp0,Yl0,resTr]   = vbm_vol_resize({Ym,Yp0,Yl0},'reduceV',vx_vol,opt.res,64); 
  Yl1A                 = vbm_vol_resize(Yl1A   ,'reduceV',vx_vol,opt.res,64,'nearest'); 

  Yl1Ans = round(Yl1A/2)*2-1;

  %% gradients and divergence
  [gx,gy,gz]=vbm_vol_gradient3(Ym); Yg = abs(gx)+abs(gy)+abs(gz); Yg=Yg./Ym;  clear gx gy gz;
  [gx,gy,gz]=vbm_vol_gradient3(max(2,Ym)); Ydiv = smooth3(divergence(gy,gx,gz)); clear gx gy gz;
  %CSFD = smooth3(vbdist(single(Ym<2 | Yp0==0))); 
  %[gx,gy,gz]=vbm_vol_gradient3(CSFD); div2=divergence(gy,gx,gz); clear gx gy gz;
  
  %% alignment of high intensity structures 
  Yl1=zeros(size(Ym),'single'); Yl1(vbm_vol_morph(Yp0<0.5 & Ym<0.4,'ldo',2,resTr.vx_volr)==1)=-1; Yl1(Yl1A>0)=0;       % backgound
  Yl1(Ym>2.25 & (Yl0==2 | Yl0==1) & ~vbm_vol_morph(Yp0>0,'lc')) = opt.LAB.HD(1);                                                           % head
  for s=1:2
    Yl1(vbm_vol_morph((Ym>2.8 & Yp0>2.8) & Ym<3.25 & Yl1A==opt.LAB.CT(s),'lab')) = opt.LAB.CT(1); % cortex
    Yl1(vbm_vol_morph((Ym>2.5 | Yp0>2.2) & Ym<3.25 & Yl1A==opt.LAB.BS(s),'lab')) = opt.LAB.BS(1); % brainstem
    Yl1(vbm_vol_morph((Ym>2.5 | Yp0>2.2) & Ym<3.25 & Yl1A==opt.LAB.MB(s),'lab')) = opt.LAB.MB(1); % midbrain
    Yl1(vbm_vol_morph((Ym>2.5 | Yp0>2.2) & Ym<3.25 & Yl1A==opt.LAB.CB(s),'lab')) = opt.LAB.CB(1); % cerebellum
    Yl1(vbm_vol_morph((Ym>2.5 | Yp0>2.2) & Ym<3.50 & Yl1A==opt.LAB.ON(s),'lab')) = opt.LAB.ON(1); % optical nerv
    Yl1(vbm_vol_morph(Yl1A==opt.LAB.HC(s) & Ym>1.75 & Ym<2.25,'lab'))  = opt.LAB.HC(1);           % hippocampus
  end
  
  % claustrum
  BGD = vbdist(single(Yl1Ans==opt.LAB.BG(1)),Yp0>0);
  Yl1(smooth3(BGD<6 & BGD>3 & Ydiv<0.1 & Yg>0.02 & Ym>2.1 & Ym<3.1 & Yl1==0 & Yl1Ans==opt.LAB.CT(1))>0.5)=opt.LAB.CT(1);
  
  % region-growing for special high intensity regions
  Yl1(((Ym<=2.85 | (Yl1Ans~=opt.LAB.CT(1) & Ydiv>0.05 & Ym<=2.9)) & Yl1==0) | Yl0==3)=-inf;
  [Yl1,D] = vbm_vol_simgrow(Yl1,Ym,0.05); Yl1(isinf(Yl1) | D>0.2)=0; 
  
  % alignment of medium and low intensity structures
  for s=1:2
    Yl1(vbm_vol_morph(Yl1A==opt.LAB.BG(s) & Yl1==0 & Ym>1.75 & Ym<2.75 & Yg<0.1 & Ydiv>0,'lab')) = opt.LAB.BG(1);  % basal ganglia
    Yl1(vbm_vol_morph(Yl1A==opt.LAB.TH(s) & Yl1==0 & Ym>1.75 & Ym<2.75 & Yg<0.1,'lo',0)) = opt.LAB.TH(1);          % hypothalamus
    Yl1(vbm_vol_morph(Yl1A==opt.LAB.HC(s) & Yl1==0 & Ym>1.50 & Ym<2.50 & Yg<0.1,'lo',0)) = opt.LAB.HC(1);          % hippocampus
  end
  
  % region growing
  Yl1(((Ym<=2 | Yl0==3 | Yg>0.15 | (Yl1Ans==opt.LAB.CT(1) & Ydiv<0.05) | ...
    Yl1Ans==opt.LAB.BG(1) & Ydiv<0.05 & Yg>0.1) & Yl1==0) | Yl0==3) = -inf;
  [Yl1,D] = vbm_vol_simgrow(Yl1,Ym,0.05); Yl1(isinf(Yl1) | D>0.15)=0;
  Yl1(((Ym<=2 | Yl0==3 | Yg>0.15 | (Yl1Ans==opt.LAB.CT(1) & Ydiv<0.02) | ...
    Yl1Ans==opt.LAB.BG(1) & Ydiv<0 & Yg>0.05) & Yl1==0) | Yl0==3) = -inf;
  [Yl1,D] = vbm_vol_simgrow(Yl1,Ym,0.02); Yl1(isinf(Yl1) | D>0.15)=0;     

  % refinement of brainstem and midbrain
  Yl1(Yl1==opt.LAB.BS(1) & ~vbm_vol_morph(Yl1==opt.LAB.BS(1),'lo',1))=0;
  Yl1(Yl1==opt.LAB.MB(1) & ~vbm_vol_morph(Yl1==opt.LAB.MB(1),'lo',1))=0;
  
  % refinement of hypothamlamus
  YM = Yl1==opt.LAB.TH(1) & ~vbm_vol_morph(Yl1==opt.LAB.TH(1),'o',2);
  Yl1(YM & Ym>2.75)=1;  Yl1(YM & Ym<=2.5)=0;
  
  % refinement of hippocampus
  Yl1((Yl1Ans~=opt.LAB.HC(1) | Ydiv<-0.02 | Ym>2.25) & Yl1==opt.LAB.HC(1)) = 1; 
  Yl1 = vbm_vol_median3c(Yl1,Yl1==opt.LAB.HC(1));
  Yl1(vbm_vol_morph(Yl1==opt.LAB.HC(1) & Yl1Ans==opt.LAB.HC(1),'lc',2))=opt.LAB.HC(1);
  for s=0:1, 
    Yl1(mod(Yl1A,2)==s & Yl1==opt.LAB.HC(1) & ...
      ~vbm_vol_morph(mod(Yl1A,2)==s & Yl1==opt.LAB.HC(1),'l')) = 0; 
  end
  
 
  
  %% refinement of the basal ganglia
  Yl1( ( (Ydiv<0.00 & Yl1Ans~=opt.LAB.BG(1)) | (Ydiv<-0.05 & Yl1Ans==opt.LAB.BG(1)) | Yg>0.1) & ...
    Yl1==opt.LAB.BG(1) & Ym>1.75 & Ym<2.5) = 0;  % basasl ganglia corr
  Yl1 = vbm_vol_median3c(Yl1,Ym>2 & Yl1>0,Ym>2 & Yl1>0);                                    % region smoothing
  for i=1:8, Yl1 = vbm_vol_median3c(Yl1,Yl1>1 & Ym<2.7,Yl1>0 & Ym<2.9); end                 % region smoothing
  Yl1((Ym<=2.25 | (Yl1Ans~=opt.LAB.CT(1) & Ydiv>-0.1 & Ym<=2.9)) & Yl1==0)=-inf;
  
  Yl2 = Yl1; Yl2(Yl2==opt.LAB.BG(1) | (Yl1Ans==opt.LAB.BG(1) & Yl1==0)| Ym<2.1 ) = -inf; 
  Yl2 = vbm_vol_simgrow(Yl2,Ym,0.1); Yl2 = vbm_vol_simgrow(Yl2,Ym,0.1);  Yl2(isinf(Yl2))=0;
  Yl1 = max(Yl1,Yl2); clear Yl2;
  Yl1((Ym<=2.0 | (Ym<1.5 & Yl1Ans==opt.LAB.BG(1))) & Yl1==0) = -inf;
  [Yl1,D] = vbm_vol_simgrow(Yl1,Ym,0.05); Yl1(isinf(Yl1) | D>0.1)=0;                        % region growing
  
  Yl1 = vbm_vol_median3c(Yl1,Ym>2 & Yg>0.05 & Yl1>4);                                       % region smoothing
  
  % regino growing (GM)
  Yl2=Yl1; Yl2( Yl2==0 & ( Ym<2.0 |  Ym>3.0 | Yl2>1))=-inf; 
  [Yl2,D] = vbm_vol_downcut(Yl2,Ym,0.04); Yl2(Yl2==-inf | D>50)=0;

  Yl2( Yl2==0 & ( Ym<1.7 |  Ym>3.0 | Yl2>1))=-inf; 
  [Yl2,D] = vbm_vol_downcut(Yl2,Ym,0.02); Yl2(Yl2==-inf | D>50)=0;

  Yl2( Yl2==0 & ( Ym<1.5 |  Ym>3.0 | Yl2>1))=-inf; 
  [Yl2,D] = vbm_vol_downcut(Yl2,Ym,-0.01); Yl2(Yl2==-inf | D>50)=0;

  Yl1((Yl2==opt.LAB.CT(1))>0.55 & Ym<2.9 & Yl1==0)=opt.LAB.CT(1);
  Yl1((Yl2==opt.LAB.CB(1))>0.55 & Ym<2.9 & Yl1==0)=opt.LAB.CB(1);
  %YM = Yl1==opt.LAB.BG(1) & ~vbm_vol_morph(Yl1==opt.LAB.BG(1),'o',1);
  %Yl1(YM & Ym>2.85 & Yl1Ans~=opt.LAB.BG(1))=1;  Yl1(YM & Ym<=2.5 & Yl1Ans~=opt.LAB.BG(1))=1;
  
  Yl1(vbm_vol_morph(Yl1==opt.LAB.BG(1),'lc') & Ym>1.75 & Ym<2.75 & Yg<0.2)=opt.LAB.BG(1);
  Yl1(vbm_vol_morph(Yl1==opt.LAB.TH(1),'lc') & Ym>1.75 & Ym<2.75 & Yg<0.2)=opt.LAB.TH(1);
  
  % head
  Yl1(Yl0==2 & Yl1==0 & Ym>1 & Yp0==0 & Yl1Ans==opt.LAB.HD(1)) = opt.LAB.HD(1);
  Yl1(Yl1==0 & Ym>1.2 & Yl0==1 & Yp0==0 &...
    (Yl1Ans==opt.LAB.HD(1) | Yl1Ans==opt.LAB.BV(1) | Yl1Ans==opt.LAB.NV(1))) = opt.LAB.HD(1);
  Yl1(Yl1==opt.LAB.HD(1) & smooth3(Yl1==opt.LAB.HD(1))<0.5)=0;
  
  %% blood vessel detection
  YWMF = single(smooth3(vbm_vol_morph(Ym>2.5 & Ym<3.5 & Yp0>1.5 & Yl1Ans~=opt.LAB.BV(1),'lo',0))>0.5);
  YWMF(YWMF==0 & (Ym<=2.5 | Ym>=3.5)) = -inf; [YWMF,D] = vbm_vol_simgrow(YWMF,Ym,0.1); YWMF(isinf(YWMF) | D>50)=0; 
  YWMF = vbm_vol_morph(YWMF,'labclose',1); YWMF = single(smooth3(YWMF)>0.5);
  YWMF(~YWMF & (Ym<2.2 | ~Yp0>0.5 | Ym>3.5))=nan; YWMF=vbm_vol_downcut(YWMF,Ym*0.9,+0.02); YWMF(isnan(YWMF))=0;
  YWMF = single(smooth3(YWMF)>0.5); 
  YWMF(~YWMF & (Ym<2.0 | ~Yp0>0.5 | Ym>3.5))=nan; YWMF=vbm_vol_downcut(YWMF,Ym*0.9,-0.02); YWMF(isnan(YWMF))=0;
  YWMF = YWMF | vbm_vol_morph( Yl1Ans==opt.LAB.HC(1) | Yl1Ans==opt.LAB.BS(1) | Yl1Ans==opt.LAB.ON(1) | ...
         Yl1Ans==opt.LAB.BG(1) | Yl1Ans==opt.LAB.TH(1) | Yl1Ans==opt.LAB.MB(1) | ...
         Yl1Ans==opt.LAB.CB(1) | Yl1Ans==opt.LAB.CB(1) | Yl1Ans==opt.LAB.VT(1),'dilate',1);
  YWMF = vbm_vol_morph(YWMF,'labclose',0); YWMF = single(smooth3(YWMF)>0.5); 
  
  % high intensity structures
  YBV = zeros(size(Ym),'uint8'); 
  YBV(YBV==0 & Yl1==0 & Ym>2.3 & Yp0>2.2 & ~YWMF & ((Yl1==0 & Yl1Ans==opt.LAB.CT(1)) | ...
    Yp0==0 | Yl1Ans==opt.LAB.NV(1) | Yl1Ans==opt.LAB.BV(1) )) = 3;
  YBV(vbm_vol_morph(YBV==3,'d',1) & Ym<2.8 & Ym>2.3) = 3; 
  YBV(vbm_vol_morph(YBV==3,'d',2) & Yl1==0 & Ym>2.1) = 3; 
  
  % low intensity structures
  YBV(Yl1==0 & smooth3(YBV==0 & Yl1==0 & Ym>1.1 & Yp0<1.2 & ~YWMF & ...
    (Ym<1.6 | Yp0==0 | Yl1Ans==opt.LAB.BV(1)) & ...
    (Yp0==0 | Yl1Ans==opt.LAB.BV(1) | Yl1Ans==opt.LAB.NV(1)))>0.55) = 2;
  YBV(vbm_vol_morph(smooth3(Ym>1.25 & Ym<1.75)>0.5,'o',2)) = 2;
  
  Yl1(YBV>0) =  opt.LAB.BV(1);  
  
  
  %% alignment of the ventricels
  for s=1:2
    Yl1(vbm_vol_morph(Yl1A==opt.LAB.VT(s) & Ym<1.8,'lab') & Yp0<1.8) = opt.LAB.VT(1);        % ventricle
  end
  VT=single(Yl1==opt.LAB.VT(1)); VT((Yl1Ans==opt.LAB.NV(1) & Yp0<2) | (Yl1Ans==opt.LAB.BV(1) & Yp0<2))=2;
  VT((Yp0>2.0 | Yp0<1 | Yl1~=0) & ~VT) = -inf; VT = vbm_vol_simgrow(VT,Ym,1.5); 
  Yl1(smooth3(VT==1)>0.2 & (Ym<1.75 | Yp0<1.75 | (Yl1==0 & Ym<1.75)))=opt.LAB.VT(1);
  
  % alignment of HD
  Yl1(Yl1==0 & Ym>2.5 & vbm_vol_morph(Yp0==0,'lc')) = opt.LAB.HD(1); 
  
    
  %% region-growing in GM only for non-blood vessels regions
  YM = Yl1==0 & Ym>0.9 & Yp0>0.5 & vbm_vol_morph(Yl1>0 & Yl1<20,'lc'); Yl1(YM)=Yl1(YM);
  Yl2=Yl1; Yl2(Yl2>1 | Ym<1.2 | Yl0==3)=-inf; 
  Yl2=smooth3(vbm_vol_downcut(Yl2,Ym,-.1,resTr.vx_volr))>0.5;
  Yl1(Yl1==0 & Yl2 & Ym>0.9)=opt.LAB.CT(1); Yl1(isinf(Yl1))=0;

  Yl1((Yl1==0 | Yl1==opt.LAB.BV(1)) & Yp0<0.5 & vbm_vol_morph(Yl1==opt.LAB.HD(1),'lc',2))=opt.LAB.HD(1);
  Yl1(Yl1==0 & ~vbm_vol_morph(Yl1>0,'lc',2))=-inf; Yl1=vbm_vol_downcut(Yl1,Ym,-.1,resTr.vx_volr);
  Yl1=vbm_vol_median3c(Yl1,Yl1==0,Yl1>0);
  
  %% prepare mask for filling of subcortical regions
  M3 = Yl1Ans==opt.LAB.BG(1)|Yl1Ans==opt.LAB.VT(1)|Yl1Ans==opt.LAB.TH(1);
  M2 = vbm_vol_morph(M3,'dilate',2,resTr.vx_volr);
  M  = 3*vbm_vol_smooth3X(single(vbm_vol_morph(vbm_vol_morph(vbm_vol_morph(M3 | ...
       (M2 & Ym>2.7 & Yl1~=7 & Yl1~=8),'labclose',1),'labopen',1),'erode',0)),0.5);  
  M2 = vbm_vol_morph(M3,'distdilate',4,resTr.vx_volr);
  YMF = M2.*max(Ym,M); clear M2 M;  
  mf3T = max(Ym,YMF); TI3FS = vbm_vol_smooth3X(mf3T,1.5); mf3T(mf3T>3 & YMF>3)=TI3FS(mf3T>3 & YMF>3);
  %Yl1(isinf(Yl1))=0; Yl1=vbm_vol_median3c(Yl1,Yl1>=0 & Ym<3,Yl1>=0); Yl1(Yl1==0 & ~Yp0 & Ym<0.75)=-inf;
  
  
  %% corrections
  % correct CSF in special brain regions
  SR = (Yl1==opt.LAB.BG(1) | Yl1==opt.LAB.TH(1) | Yl1==opt.LAB.HC(1));
  Yl1(SR & Yl1Ans==opt.LAB.CT(1) & Yp0<1.5) = opt.LAB.CT(1);    
  Yl1(SR & Yl1Ans==opt.LAB.VT(1) & Yp0<1.5) = opt.LAB.VT(1);
  Yl1(Yl1==opt.LAB.BV(1) & Yp0==0) = opt.LAB.HD(1);
  clear SR; 
    
  % corrections for non brain parts
  Yl1(smooth3(Yp0<1 & Ym>1.5 & (Yl1==0 | Yl1==opt.LAB.BV(1)))>0.5) = opt.LAB.HD(1); 
    
  % cleanup for labels that are not so excact in this version of Yl1A
  Yl1(Yl1==opt.LAB.HC(1))=opt.LAB.CT(1);
  Yl1=vbm_vol_median3c(Yl1,Yl1>=0 & Yl1<5,Yl1>=0 & Yl1<5);
  Yl1=vbm_vol_median3c(Yl1,Yl1>=0 & Yl1<5,Yl1>=0 & Yl1<5);
  
  
  %% brain mask
  YB  = Yl1>0 & Yl1<opt.LAB.HD(1) & Yl1~=opt.LAB.BV(1);
  YB  = vbm_vol_morph(smooth3(YB)>0.5,'o',3);
  YB  = vbm_vol_morph(smooth3(YB)>0.5,'lc');
  Yl1B = Yl1>0; 
  
  % fine correction based on SPM mask
  Yl1(~YB & Yl1~=opt.LAB.HD(1) & Yp0>0 ) = 0;
  Yl1( YB & Yl1==opt.LAB.HD(1) & Yp0>0 ) = 0; 
  Yl1(~YB & Yl1~=opt.LAB.HD(1) & Yp0==0 & Yl1>0) = opt.LAB.HD(1);   
  Yl1( YB & Yl1==opt.LAB.HD(1) & Yp0>0  & Yl1>0) = opt.LAB.HD(1);
  [D,I]=vbdist(Yl1,Yl1B); Yl1=Yl1(I);
  
  
  %% side aligment using laplace to correct for missalignments due to the normalization
  d = 5; M = vbm_vol_smooth3X(single(vbm_vol_morph(mod(Yl1A,2)==0,'distdilate',d,resTr.vx_volr)) & ...
      single(vbm_vol_morph(mod(Yl1A,2)==1,'distdilate',d,resTr.vx_volr)==1),20);
  S = 2*single(mod(Yl1A,2)==0 & Ym>2.5 & M<max(M(:))*0.9) +single(mod(Yl1A,2)==1 & ...
      Ym>2.5 &  M<max(M(:))*0.9); S(mf3T<=2.5)=-inf; S(S==0)=1.5;

  rS=vbm_vol_resize(S,'reduce'); rS=round(rS*2)/2; 
  rS=vbm_vol_laplace3R(rS,rS==1.5,0.001); 
  rS=vbm_vol_resize(single(rS),'dereduce',size(Ym)); 
  S(S==1.5)=round(rS(S==1.5)); 
  S(isinf(S) & mf3T>2)=0; S=vbm_vol_downcut(S,mf3T,3,resTr.vx_volr); 
  S(isinf(S) & mf3T>0)=0; S=vbm_vol_downcut(S,mf3T,3,resTr.vx_volr); 
  S(S<=0)=2-mod(Yl1A(S<=0),2);
  S=round(vbm_vol_smooth3X(S,2));
  Yl1(Yl1>0)=Yl1(Yl1>0)+(S(Yl1>0)==2); 
  Yl1(Yl1<0.5)=0;
  
  
  %% back to original resolution and full size
  YMF = vbm_vol_resize(YMF,'dereduceV',resTr);
  Yp0 = vbm_vol_resize(Yp0,'dereduceV',resTr);
  YB  = vbm_vol_resize(YB ,'dereduceV',resTr);
  Yl1 = vbm_vol_resize(Yl1,'dereduceV',resTr,'nearest');
  
  YMF = vbm_vol_resize(YMF,'dereduceBrain',BB);
  Yl1 = vbm_vol_resize(Yl1,'dereduceBrain',BB);
  Yp0 = vbm_vol_resize(Yp0,'dereduceBrain',BB);
  YB  = vbm_vol_resize(YB ,'dereduceBrain',BB);
  
  Ym  = mgTO; clear mgTO;
  
  
  %% setting head label for voxel outside the bounding box
  Yl1(smooth3(Yp0<1 & Ym>1.5 & (Yl1==0 | Yl1==opt.LAB.BV(1)))>0.5) = opt.LAB.HD(1);
  [HDr,resTr] = vbm_vol_resize(Yl1>0,'reduceV',vx_vol,3,64);
  HDr = smooth3(vbm_vol_morph(HDr>0.1,'lc',2))>0.5;
  HD  = smooth3(Yl1>0 | vbm_vol_resize(HDr,'dereduceV',resTr))>0.5;
  [D,I,Yl1] = vbdist(single(Yl1),HD); clear D I;
  Yl1(HD & Yl1<=0) = opt.LAB.HD(1); 
  Yl1(~HD) = 0;

  
  %% subvolumes
  fn = setdiff(fieldnames(opt.LAB),{'NV','HD','NB','ON','HC'}); vol=struct();
  for fni=1:numel(fn)
    eval(sprintf(['vol.vol_abs_%s = prod(vx_vol)/1000 .* ' ...
      '[sum(Yl1(:)==%d) sum(Yl1(:)==%d)];'],fn{fni},def.LAB.(fn{fni}))); 
  end 
  vol.vol_LRP = sum( mod(Yl1(:),2)==1 & Yl1(:)>0 & Yl1(:)<20) ./ ...
                sum( Yl1(:)>0 & Yl1(:)<20);

end

% old code ... delete later
  %{
  ( Ym>1.3 & (Yp0==0 | Yl1Ans==opt.LAB.NV(1) | Yl1Ans==opt.LAB.BV(1)))) ; 
  BV = smooth3(BV)>0.5;
  %%
  BV = (Yl1==0 & Ym>2.5 & ~WM & Yl1Ans~=opt.LAB.HD(1) & vbm_vol_morph(Yp0>0,'labclose')) | ...
       (Yl1==0 & Ym>2.3 & ~WM & Yl1Ans==opt.LAB.BV(1) & vbm_vol_morph(Yp0==0,'labclose')) | ...
       (Yl1==0 & Ym>2.3 & ~WM & Yl1Ans~=opt.LAB.HD(1))| ...
        vbm_vol_smooth3X(Ym>1.25 & Ym<1.75 & Yp0<0.1 & ...
        (Yl1Ans==opt.LAB.HD(1) | Yl1Ans==opt.LAB.BV(1) | Yl1Ans==opt.LAB.NV(1)) )>0.8 | ...
        vbm_vol_smooth3X((Yl1==0 & (Yl1Ans==opt.LAB.BV(1)) & ...
       ((Ym>1.25 & Ym<1.75 & Yp0<1.75) | (Ym>2.3)) & ~WM))>0.5; BV=single(smooth3(BV)>0.2);
  
  BV = BV | ...
      (~WM & vbm_vol_morph(Yp0<1.75,'e',2) & Ym>1.25 & (Yl1Ans==opt.LAB.NV(1) | ...
        Yl1Ans==opt.LAB.BV(1)|  Yl1Ans==opt.LAB.HD(1) | Yl0==1));
  
    %% 
  BV = BV | (vbm_vol_morph(BV,'dilate',1) & Ym>2.3 & ~WM); 
  BV = BV | (vbm_vol_morph(BV,'labclose',1) & (Ym>2.2 | (Ym>1.2 & Ym<1.6)) & ~WM);
  %BV = BV & vbm_vol_morph(Yp0>0,'lc');
  Yl1(BV) =  opt.LAB.BV(1);  
 % clear WM;
 %%
   Yl1(Yl1==0 & Ym<1.5)=-inf; [Yl1,D] = vbm_vol_downcut(Yl1,Ym,-0.01); Yl1(Yl1==-inf | D>100)=0;
  %}

  
  %{
  
  %[Yl1,D] = vbm_vol_simgrow(Yl1,Ym,0.1); Yl1(D>50)=0; %Yl1 = vbm_vol_median3c(Yl1,Ym>1.5);
  %Yl1(BV) = opt.LAB.BV(1); Yl1(BG) =  opt.LAB.BG(1);
  %Yl1((Yl1==0 & Ym>1.3 & ~WM & Yl1Ans==opt.LAB.BV(1) & vbm_vol_morph(Yp0==0,'labclose')))= opt.LAB.BV(1);
  %%
  Yl1=vbm_vol_downcut(Yl1,Ym,0.2,resTr.vx_volr);
%  Yl1(isinf(Yl1))=0; Yl1=vbm_vol_median3c(Yl1,Yl1>=0 & Ym<2.5 & ...
%    ~(Yl1==opt.LAB.BV(1) & Ym>2.25),Yl1>=0); Yl1(Yl1==0 & ~Yp0 & Ym<0.75)=-inf;
  Yl1(Yl1==opt.LAB.BV(1) & (Yl1Ans==opt.LAB.HD(1) | Yl1<1.5)) = opt.LAB.HD(1); Yl1(isinf(Yl1))=0;
  %Yl1((Ym<=1.00) & Yl1==0) = -inf; 
  %[Yl1,D] = vbm_vol_simgrow(Yl1,Ym,0.1); Yl1(D>50)=0; Yl1 = vbm_vol_median3c(Yl1,Ym>1.5);
   Yl1(YBV>0) =  opt.LAB.BV(1);  
  %}

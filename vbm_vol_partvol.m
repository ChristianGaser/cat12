function [l1T,MF] = vbm_vol_partvol(l1A,p0T,mgT,opt)
% ______________________________________________________________________
% Use a segment map p0T, the global intensity normalized T1 map mgT and 
% the atlas label map l1T to create a individual label map l1T. 
% The atlas contain main regions like cerebrum, brainstem, midbrain,
% cerebellum, ventricle, and regions with blood vessels. 
%
% This function try to solve the following problems:
%  1) Finding of the cerebrum, the cerebellum, the head, blood vessels, 
%     brain skin and other mayor structures based on atlas (l1A) and 
%     tissue class information (p0T). 
%     To do this it is important to use data from the T1-map (mgT) that
%     use the same intensity scalling as the segment map p0T, but have 
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
%   - stärkerer Abgleich mit p0T
%   - stärkere Einbindung der LAS
% ______________________________________________________________________
%
% Structure:
%
%   [l1T,pfT,mfT] = vbm_vol_partvol(l1A,p0T,mgT,opt)
%
%   INPUT:  l1A = 3D-volume with brain regions (altas map)
%           p0T = 3D-volume with tissue propability map (CSF=1,GM=2;WM=3)
%           mgT = intensity normalized T1 image
%           opt
%            .res    = resolution for mapping
%            .vx_vol = Voxelsize
%            .LAB    = Label of l1T map (see def.LAB definition below)
%            
%
%   OUTPUT: l1T = individual label map 
%           pfT = Corrected and filled p0T map 
%           mfT = Corrected and filled mgT map
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
  % ds('l2','',vx_vol,mgT,l1A,mgT,p0T/3,130)
  mgT = mgT*3; mgTO = mgT;
  [mgT,p0T,l1A,BB] = vbm_vol_resize({mgT,p0T,l1A},'reduceBrain',vx_vol,5,p0T>0);   % removing of background
  [mgT,p0T,resTr]  = vbm_vol_resize({mgT,p0T},'reduceV',vx_vol,opt.res,64); 
  l1A              = vbm_vol_resize(l1A   ,'reduceV',vx_vol,opt.res,64,'nearest'); 

  l1ANS = round(l1A/2)*2-1;

  
  %% alignment of high intensity structures 
  l1T=zeros(size(mgT),'single'); l1T(vbm_vol_morph(p0T<0.5 & mgT<0.4,'ldo',2,resTr.vx_volr)==1)=-1; l1T(l1A>0)=0;       % backgound
  l1T(vbm_vol_morph(p0T<0.5,'open')==1 & mgT>2.5) = opt.LAB.HD(1);                                                    % head
  for s=1:2
    l1T(vbm_vol_morph(mgT>2.5 & mgT<3.25 & l1A==opt.LAB.CT(s),'lab')) = opt.LAB.CT(1); % cortex
    l1T(vbm_vol_morph(mgT>2.5 & mgT<3.25 & l1A==opt.LAB.BS(s),'lab')) = opt.LAB.BS(1); % brainstem
    l1T(vbm_vol_morph(mgT>2.5 & mgT<3.25 & l1A==opt.LAB.MB(s),'lab')) = opt.LAB.MB(1); % midbrain
    l1T(vbm_vol_morph(mgT>2.5 & mgT<3.25 & l1A==opt.LAB.CB(s),'lab')) = opt.LAB.CB(1); % cerebellum
    l1T(vbm_vol_morph(mgT>2.5 & mgT<3.25 & l1A==opt.LAB.ON(s),'lab')) = opt.LAB.ON(1); % optical nerv
  end
  
  
  %% region-growing for special high intensity regions
  l1T((mgT<=2.75 & l1T==0)) = -inf; l1T = vbm_vol_simgrow(l1T,mgT,1);l1T(isinf(l1T))=0; 
  l1T(vbm_vol_morph(l1T==opt.LAB.CT(1),'labclose')) = opt.LAB.CT(1);
  l1T((mgT<=2.75 & l1T==0)) = -inf; l1T = vbm_vol_simgrow(l1T,mgT,1);l1T(isinf(l1T))=0; 
  
  
  %% blood vessel detection
  WM = single(smooth3(vbm_vol_morph(mgT>2.5 & mgT<3.5 & p0T>1.5 & l1ANS~=opt.LAB.BV(1),'lo',0))>0.5);
  WM(WM==0 & (mgT<=2.5 | mgT>=3.5)) = -inf; [WM,D] = vbm_vol_simgrow(WM,mgT,0.1); WM(isinf(WM) | D>50)=0; 
  WM = vbm_vol_morph(WM,'labclose',1); WM = single(smooth3(WM)>0.5);
  WM(~WM & (mgT<2.2 | ~p0T>0.5 | mgT>3.5))=nan; WM=vbm_vol_downcut(WM,mgT*0.9,+0.02); WM(isnan(WM))=0;
  WM = single(smooth3(WM)>0.5); 
  WM(~WM & (mgT<2.0 | ~p0T>0.5 | mgT>3.5))=nan; WM=vbm_vol_downcut(WM,mgT*0.9,-0.02); WM(isnan(WM))=0;
  WM = WM | vbm_vol_morph( l1ANS==opt.LAB.HC(1) | l1ANS==opt.LAB.BS(1) | l1ANS==opt.LAB.ON(1) | ...
         l1ANS==opt.LAB.BG(1) | l1ANS==opt.LAB.TH(1) | l1ANS==opt.LAB.MB(1) | ...
         l1ANS==opt.LAB.CB(1) | l1ANS==opt.LAB.CB(1) | l1ANS==opt.LAB.VT(1),'dilate',1);
  WM = vbm_vol_morph(WM,'labclose',0); WM = single(smooth3(WM)>0.5); 
  
  BV = (l1T==0 & mgT>2.2 & ~WM) | (l1T==0 & l1ANS==opt.LAB.BV(1) & ...
       ((mgT>1.2 & mgT<1.5) | (mgT>2.2)) & ~WM) ; BV=single(smooth3(BV)>0.3);
  BV = BV | (vbm_vol_morph(BV,'dilate',1) & mgT>2.2 & ~WM); 
  BV = BV | (vbm_vol_morph(BV,'labclose',1) & (mgT>2.2 | (mgT>1.2 & mgT<1.6)) & ~WM);
  BV = BV & p0T>0;
  l1T(BV) =  opt.LAB.BV(1);  
  clear WM;
  
  
  %% alignment of medium and low intensity structures
  l1T(l1ANS==opt.LAB.BG(1) & l1T==0 & mgT>1.75 & mgT<2.75 & mgT<2.9) = opt.LAB.BG(1);          % basal ganglia
  l1T(l1ANS==opt.LAB.TH(1) & mgT>1.75 & mgT<2.75 & mgT<2.9) = opt.LAB.TH(1);                   % hypothalamus
  l1T(l1ANS==opt.LAB.HC(1) & mgT>1.75 & mgT<2.75 & mgT<2.9) = opt.LAB.HC(1);                   % hippocampus
  % alignment of the ventricels
  for s=1:2
    l1T(vbm_vol_morph(l1A==opt.LAB.VT(s) & mgT<1.75,'lab') & p0T<1.75) = opt.LAB.VT(1);        % ventricle
  end
  VT=single(l1T==opt.LAB.VT(1));
  VT(l1ANS==opt.LAB.NV(1) & (p0T>=1 & p0T<1.75))=2;
  VT((p0T>1.5 | p0T<1 | l1T~=0) & ~VT) = -inf; VT = vbm_vol_simgrow(VT,mgT,1); 
  l1T(smooth3(VT==1)>0.5)=opt.LAB.VT(1);
  % alignment of HD
  l1T(p0T==0 & mgT>2) = opt.LAB.HD(1); 
  
  
  %% refinement of WM 
  l1T(l1ANS==opt.LAB.CT(1) & l1T==0 & vbm_vol_morph(mgT>2 & l1T==1,'labopen',1)) = opt.LAB.CT(1);  
  l1T(l1ANS==opt.LAB.BG(1) & ~vbm_vol_morph(l1ANS==opt.LAB.BG(1),'erode',1) & ...
    l1T==opt.LAB.BG(1) & (mgT>2.75 & mgT>2.5)) = opt.LAB.CT(1);                                % basal ganglia
  l1T=vbm_vol_median3c(l1T,l1T>=0,l1T>0);  

  
  %% regino-growing for all high intensity regions
  l1T((mgT<=2.25 & l1T==0))=-inf; l1T = vbm_vol_simgrow(l1T,mgT,0.5); l1T(isinf(l1T))=0; 
  l1T(l1T==0 & vbm_vol_morph(mgT>2.75 & p0T,'labopen',1)==0 & ...
    mgT>2.75 & l1T==0 & vbm_vol_morph(mgT<2,'labclose',1,resTr.vx_volr))=opt.LAB.BV(1);

  
  %% region-growing in GM only for non-blood vessels regions
  l1T(l1T==opt.LAB.BV(1) | isinf(l1T))=0; l1T((mgT<=1.9 | ~p0T) & l1T==0) = -inf; 
  l1T = vbm_vol_downcut(l1T,mgT,0.02,resTr.vx_volr); 
  l1T(l1T==opt.LAB.BV(1) | isinf(l1T))=0; l1T((mgT<=1.0 | ~p0T) & l1T==0) = -inf; 
  l1T = vbm_vol_downcut(l1T,mgT,-0.01,resTr.vx_volr); 
  %l1T((l1ANS==opt.LAB.BV(1) | D>500) & mgT>1.5 & ~vbm_vol_morph(mgT>2.5,'lab') & ...
  %  (l1ANS~=opt.LAB.CB(1) | l1ANS~=opt.LAB.HD(1)))=opt.LAB.BV(1);                         % adding blood vessels again
  l1T(l1T==opt.LAB.BG(1) & ~(l1ANS==opt.LAB.BG(1)|l1ANS==opt.LAB.VT(1)))=opt.LAB.CT(1);   % basal ganglia (avoid overgrowing)
  l1T(l1T==0 & mgT>3.5)=opt.LAB.BV(1); l1T=vbm_vol_downcut(l1T,mgT,0.2,resTr.vx_volr);
  l1T(isinf(l1T))=0; l1T=vbm_vol_median3c(l1T,l1T>=0 & mgT<2.5 & ...
    ~(l1T==opt.LAB.BV(1) & mgT>2.25),l1T>=0); l1T(l1T==0 & ~p0T & mgT<0.75)=-inf;

  
  %% prepare mask for filling of subcortical regions
  M3 = l1ANS==opt.LAB.BG(1)|l1ANS==opt.LAB.VT(1)|l1ANS==opt.LAB.TH(1);
  M2 = vbm_vol_morph(M3,'dilate',2,resTr.vx_volr);
  M  = 3*vbm_vol_smooth3X(single(vbm_vol_morph(vbm_vol_morph(vbm_vol_morph(M3 | ...
       (M2 & mgT>2.7 & l1T~=7 & l1T~=8),'labclose',1),'labopen',1),'erode',0)),0.5);  
  M2 = vbm_vol_morph(M3,'distdilate',4,resTr.vx_volr);
  MF = M2.*max(mgT,M); clear M2 M;  
  mf3T = max(mgT,MF); TI3FS = vbm_vol_smooth3X(mf3T,1.5); mf3T(mf3T>3 & MF>3)=TI3FS(mf3T>3 & MF>3);
  l1T(isinf(l1T))=0; l1T=vbm_vol_median3c(l1T,l1T>=0 & mgT<3,l1T>=0); l1T(l1T==0 & ~p0T & mgT<0.75)=-inf;
  
  
  %% corrections
  % correct CSF in special brain regions
  SR = (l1T==opt.LAB.BG(1) | l1T==opt.LAB.TH(1) | l1T==opt.LAB.HC(1));
  l1T(SR & l1ANS==opt.LAB.CT(1) & p0T<1.5) = opt.LAB.CT(1);    
  l1T(SR & l1ANS==opt.LAB.VT(1) & p0T<1.5) = opt.LAB.VT(1);
  clear SR; 
    
  % corrections for non brain parts
  l1T(BV) = opt.LAB.BV(1);  
  l1T(smooth3(p0T<1 & mgT>1.5 & (l1T==0 | l1T==opt.LAB.BV(1)))>0.5) = opt.LAB.HD(1); 
    
  % cleanup for labels that are not so excact in this version of l1A
  l1T(l1T==opt.LAB.HC(1))=opt.LAB.CT(1);
  
  
  %% side aligment using laplace to correct for missalignments due to the normalization
  d = 5; M = vbm_vol_smooth3X(single(vbm_vol_morph(mod(l1A,2)==0,'distdilate',d,resTr.vx_volr)) & ...
      single(vbm_vol_morph(mod(l1A,2)==1,'distdilate',d,resTr.vx_volr)==1),20);
  S = 2*single(mod(l1A,2)==0 & mgT>2.5 & M<max(M(:))*0.9) +single(mod(l1A,2)==1 & ...
      mgT>2.5 &  M<max(M(:))*0.9); S(mf3T<=2.5)=-inf; S(S==0)=1.5;

  rS=vbm_vol_resize(S,'reduce'); rS=round(rS*2)/2; 
  rS=vbm_vol_laplace3R(rS,rS==1.5,0.001); 
  rS=vbm_vol_resize(single(rS),'dereduce',size(mgT)); 
  S(S==1.5)=round(rS(S==1.5)); 
  S(isinf(S) & mf3T>2)=0; S=vbm_vol_downcut(S,mf3T,3,resTr.vx_volr); 
  S(isinf(S) & mf3T>0)=0; S=vbm_vol_downcut(S,mf3T,3,resTr.vx_volr); 
  S(S<=0)=2-mod(l1A(S<=0),2);
  S=round(vbm_vol_smooth3X(S,2));
  l1T(l1T>0)=l1T(l1T>0)+(S(l1T>0)==2); 
  l1T(l1T<0.5)=0;

  
  %% back to original resolution
  MF  = vbm_vol_resize(MF,'dereduceV',resTr);
  p0T = vbm_vol_resize(p0T,'dereduceV',resTr);
  l1T = vbm_vol_resize(l1T,'dereduceV',resTr,'nearest');
  [MF,l1T,p0T] = vbm_vol_resize({MF,l1T,p0T},'dereduceBrain',BB);
  mgT = mgTO;
  
  
  %% setting head label for voxel outside the bounding box
  l1T(smooth3(p0T<1 & mgT>1.5 & (l1T==0 | l1T==opt.LAB.BV(1)))>0.5) = opt.LAB.HD(1);
  [HDr,resTr] = vbm_vol_resize(l1T>0,'reduceV',vx_vol,3,64);
  HDr = smooth3(vbm_vol_morph(HDr>0.1,'lc',2))>0.5;
  HD  = smooth3(l1T>0 | vbm_vol_resize(HDr,'dereduceV',resTr))>0.5;
  [D,I,l1T] = vbdist(single(l1T),HD); clear D I;
  l1T(HD & l1T<=0) = opt.LAB.HD(1); 
  l1T(~HD) = 0;
end

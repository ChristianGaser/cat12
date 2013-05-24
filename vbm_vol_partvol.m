function [vol,Yl1,Yb,YMF] = vbm_vol_partvol(Yl1A,Yp0,Ym,Yl0,opt)
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
%     2.5) Brain extraction
%     2.6) Side alignment
% ______________________________________________________________________
%
% Structure:
%
%   [vol,Yl1,Yb,YMF] = vbm_vol_partvol(Yl1A,Yp0,Ym,Yl0,opt)
%
%   INPUT:  Yl1A = 3D-volume with brain regions (altas map)
%           Yp0  = 3D-volume with tissue propability map (CSF=1,GM=2;WM=3)
%           Ym   = intensity normalized T1 image (BG=0,CSF=1/3,GM=2/3,WM=1)
%           Yl0  = spm-classes 4-6 (intracranial=1,skull=2,background=3)
%           opt
%            .res    = resolution for mapping
%            .vx_vol = voxelsize
%            .LAB    = label of Yl1 map (see def.LAB definition below)
%            
%
%   OUTPUT: vol = structure with volumes
%           Yl1 = individual label map 
%           Yb  = brain mask
%           YMF = filling mask for ventricle and subcortical structures
%
% ______________________________________________________________________
% Structural Brain Mapping Group, University Jena, Germany
% Robert Dahnke
% 2013/05
%
% $Id$

% ______________________________________________________________________
%
% Development comments:
%
%   Was ist neu im Vergleich zu anderen?
%   - Zuweisung durch Dartel mit hoher Genauigkeit möglich 
%   - Erweiterung von SPM/VBM durch MainROIs (Seiten, Lappen, ...)
%   - Verbesserung der SPM/VBM durch bessere Enfernung von unerwünschtem
%     Gewebe (ON, Blutgefäße ...)
%   - Blutgefäße könnnen als erweitere Masken für fMRI genutzt werden um
%     Seiteneffekte besser ausblenden zu können.
%  [- Beliebige Atlanten können genutzt werden.]
%
%  Todo:
%   - Besserer Atlas
%   - BV vs. HD - glätten in dilated HD region
%   - Füllen von CSF lücken bei LAB~=BV und Ym<1.2 und LAB==NV?
%
% ______________________________________________________________________


  if ~exist('opt','var'), opt=struct(); end
  def.res    = 1;
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
  def.LAB.HI = [23 24]; % WM hyperintensities
  
  def.color.commend = [0 0 0.8];
  def.color.warning = [0.8 0 0];
  
  opt = checkinopt(opt,def);
  
  
%  opt.res = max(min(1.2,opt.res*2),opt.res);
  
  vx_vol = opt.vx_vol; %vx=mean(vx_vol);
    
  
  %% Initialization:
  % ds('l2','',vx_vol,Ym,Yl1A,Ym,Yp0/3,130)
  
  % Optimization:
  Ym = Ym*3; mgTO = Ym;
  [Ym,Yp0,Yl1A,Yl0,BB] = vbm_vol_resize({Ym,Yp0,Yl1A,Yl0},'reduceBrain',vx_vol,10,Yp0>0);   % removing of background
  [Ym,Yp0,Yl0,resTr]   = vbm_vol_resize({Ym,Yp0,Yl0},'reduceV',vx_vol,opt.res,64); 
  Yl1A                 = vbm_vol_resize(Yl1A   ,'reduceV',vx_vol,opt.res,64,'nearest'); 

  Yl1Ans = round(Yl1A/2)*2-1;

  % gradients and divergence
  [gx,gy,gz]=vbm_vol_gradient3(Ym);        Yg   = abs(gx)+abs(gy)+abs(gz); Yg=Yg./Ym;  clear gx gy gz;
  [gx,gy,gz]=vbm_vol_gradient3(max(2,Ym)); Ydiv = smooth3(divergence(gy,gx,gz));       clear gx gy gz;

  warning off 'MATLAB:vbm_vol_morph:NoObject'; 
  
  %% Find major structures  
  % alignment of high intensity structures 
  Yl1 = zeros(size(Ym),'single'); 
  YM  = Yp0<0.5 & Ym<0.4;
  Yl1(vbm_vol_morph(YM,'ldo',2,resTr.vx_volr)==1)=-1; Yl1(Yl1A>0)=0;              % backgound
 
  YM=Ym>2.25 & (Yl0==2 | Yl0==1) & ~vbm_vol_morph(Yp0>0,'lc');                    % head
  Yl1(YM) = opt.LAB.HD(1); ER.HD(1:2)=sum(YM(:))==0; 
  for s=1:2
    YM=vbm_vol_morph((Ym>2.5 & Yp0>2.5) & Ym<3.25 & Yl1A==opt.LAB.CT(s),'lc');    % cerebrum
    Yl1(YM)=opt.LAB.CT(1); ER.CT(s)=sum(YM(:))==0; 
    YM=vbm_vol_morph((Ym>2.5 | Yp0>2.2) & Ym<3.25 & Yl1A==opt.LAB.BS(s),'lc');    % brainstem
    Yl1(YM)=opt.LAB.BS(1); ER.BS(s)=sum(YM(:))==0;  
    YM=vbm_vol_morph((Ym>2.5 | Yp0>2.2) & Ym<3.25 & Yl1A==opt.LAB.MB(s),'lc');    % midbrain
    Yl1(YM)=opt.LAB.MB(1); ER.MB(s)=sum(YM(:))==0;  
    YM=vbm_vol_morph((Ym>2.5 | Yp0>2.2) & Ym<3.25 & Yl1A==opt.LAB.CB(s),'lc');    % cerebellum
    Yl1(YM)=opt.LAB.CB(1); ER.CB(s)=sum(YM(:))==0; 
    YM=vbm_vol_morph((Ym>2.1 | Yp0>2.1) & Ym<3.50 & Yl1A==opt.LAB.ON(s),'lab');   % optical nerv
    Yl1(YM)=opt.LAB.ON(1); %ER.ON(s)=sum(YM(:))==0; 
    YM=vbm_vol_morph(Yl1A==opt.LAB.HC(s) & Ym>1.75 & Ym<2.25,'lab');              % hippocampus
    Yl1(YM)=opt.LAB.HC(1); %ER.HC(s)=sum(YM(:))==0; 
  end
  warning off 'MATLAB:vbm_vol_morph:NoObject';
  
  % claustrum - to reduce BG overgrowing
  YBGD = vbdist(single(Yl1Ans==opt.LAB.BG(1)),Yp0>0,vx_vol);
  Yl1(smooth3(YBGD==0 & Ym>2.8 & Ym<3.2)>0.5)=opt.LAB.CT(1); 
  Yl1((YBGD<4 & YBGD>2 & Ydiv<0.1 & Yg>0.02 & Ym>2.3 & Ym<2.8 & Yl1==0 & Yl1Ans==opt.LAB.CT(1))>0.5)=opt.LAB.CT(1);
  clear YBGD;
  
  % region-growing for special high intensity regions
  Yl1(((Ym<=2.85 | (Yl1Ans~=opt.LAB.CT(1) & Ydiv>0.05 & Ym<=2.9)) & Yl1==0) | Yl0==3)=-inf;
  [Yl1,YD] = vbm_vol_simgrow(Yl1,Ym,0.05); Yl1(isinf(Yl1) | YD>0.1)=0; 
  
  
  % alignment of medium and low intensity structures
  for s=1:2
    YM=vbm_vol_morph(Yl1A==opt.LAB.BG(s) & Yl1==0 & Ym>1.75 & Ym<2.75 & Yg<0.1 & Ydiv>0,'lab'); % basal ganglia
    YM=vbm_vol_morph(YM,'lc',2); Yl1(YM)=opt.LAB.BG(1); ER.BG(s)=sum(YM(:))==0;  
    YM=vbm_vol_morph(Yl1A==opt.LAB.TH(s) & Yl1==0 & Ym>1.75 & Ym<2.75 & Yg<0.1,'lo',0);         % hypothalamus
    YM=vbm_vol_morph(YM,'lc',2); Yl1(YM)=opt.LAB.TH(1); ER.TH(s)=sum(YM(:))==0;          
%     YM=vbm_vol_morph(Yl1A==opt.LAB.HC(s) & Yl1==0 & Ym>1.50 & Ym<2.50 & Yg<0.1,'lo',0);         Yl1(YM)=opt.LAB.HC(1);          % hippocampus
%     if sum(YM(:))==0, vbm_io_cprintf([0.8 0 0],'\nWARNING: miss %s hippocampus. ',sides{s}); end 
  end
  

  
  
  %% Error managment, if we miss a major structure
  if any(struct2array(ER))
    sides = {'left' 'right'}; 
    ROI = {'HD','head';'CT','cerebrum';'CB','cerebellum';'MB','midbrain'; ...
           'HC','hippocampus';'BG','basal ganglia';'TH','thalamus'; ...
           'BS','brainstem'};
    for s=1:2
       fn = fieldnames(ER); fprintf(1,''); 
       for fni=1:numel(fn)
         if ER.(fn{fni})(s)
           fprintf(1,'\n');
           if strcmp('HD',fn{fni})
             if s==1
               str = sprintf('COMMEND: Did not find the head - already sklull-stripped?');
               col = opt.color.commend;
             end
           else
             str = sprintf('WARNING: Did not find %s %s.',sides{s},...
              ROI{find(cellfun('isempty',strfind(ROI(:,1),fn{fni}))==0,1),2});
             col = opt.color.warning;
           end
           str = sprintf('%s ',str,repmat(' ',1,67-length(str))); 
           vbm_io_cprintf(col,'%s',str); 
         end
       end
    end
  end
  
  
  %% region growing (WM)
  YM = Yl1==0 & smooth3(Yl1==0 & ((Ym<=2 | Yl0==3 | Yg>0.15 | (Yl1Ans==opt.LAB.CT(1) & Ydiv<0.05) | ...
    Yl1Ans==opt.LAB.BG(1) & Ydiv<0.05 & Yg>0.1)) | Yl0==3)>0.5; 
  for i=1:3, Yl1(YM) = -inf; [Yl1,YD] = vbm_vol_simgrow(Yl1,Ym,0.05); Yl1(isinf(Yl1) | YD>0.1)=0; end
  
  
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
    
  
  % refinement of the basal ganglia
  Yl1( ( (Ydiv<0.00 & Yl1Ans~=opt.LAB.BG(1)) | (Ydiv<-0.05 & Yl1Ans==opt.LAB.BG(1)) | Yg>0.1) & ...
    Yl1==opt.LAB.BG(1) & Ym>1.75 & Ym<2.5) = 0;                                             % basasl ganglia corr
  Yl1 = vbm_vol_median3c(Yl1,Ym>2 & Yl1>0,Ym>2 & Yl1>0);                                    % region smoothing
  for i=1:8, Yl1 = vbm_vol_median3c(Yl1,Yl1>1 & Ym<2.7,Yl1>0 & Ym<2.9); end                 % region smoothing
  
  
  % further region-growing (GM-WM) without ON!
  Yl2 = Yl1; Yl2(Yl2==opt.LAB.ON(1) | smooth3(Yl1==0 & (Yl2==opt.LAB.BG(1) | (Yl1Ans==opt.LAB.BG(1) & Yl1==0) | ...
    Ym<=2.25 | (Yl1Ans~=opt.LAB.CT(1) & Ydiv>-0.1 & Ym<=2.9) | Ym>3.2))>0.5) = -inf; 
  Yl2 = vbm_vol_simgrow(Yl2,Ym,0.1); [Yl2,YD] = vbm_vol_simgrow(Yl2,Ym,0.1);  Yl2(isinf(Yl2) | YD>0.1)=0;
  Yl1 = max(Yl1,Yl2); clear Yl2;
  
  
  % region growing (GM)
  %YM  = smooth3((Ym> 2.6 | Ym<=1.8 | (Ym<1.5 & Yl1Ans==opt.LAB.BG(1))) & Yl1==0)>0.5;
  %for i=1:3, Yl1(YM) = -inf; [Yl1,YD] = vbm_vol_simgrow(Yl1,Ym,1); Yl1(isinf(Yl1) | YD>0.05)=0; end 
  
  YM = Yp0>0 & Yl1~=opt.LAB.BG(1) & Yl1~=opt.LAB.ON(1) & Ym<2.5 & Ym>1.5; 
  Yl1 = vbm_vol_median3c(Yl1,YM,YM);                                     
  
  
  %% Hyperintensities
  YH = single(vbm_vol_morph(Yl1A==opt.LAB.CT(1) & Ym>1.5 & (Yp0 - Ym)>0.5,'o',1));
  YH(smooth3(Yl1~=0 | (Yp0 - Ym)<0.3 | Ym<1.5)>0.5) = -inf;
  [YH,YD] = vbm_vol_simgrow(YH,Ym,0.04);
  Yl1(YH==1) = opt.LAB.HI(1); clear YH;
  
  
  % regino growing (GM)
  Yl2=Yl1; Yl2(smooth3(Yl2==0 & ( Ym<2.0 |  Ym>3.0 | Yl2>1))>0.5)=-inf; 
  [Yl2,YD] = vbm_vol_downcut(Yl2,Ym,0.02); Yl2(Yl2==-inf | YD>50)=0;
  Yl2 = vbm_vol_median3c(Yl2,YM,YM); 
  
  Yl2(smooth3(Yl2==0 & ( Ym<1.5 |  Ym>3.0 | Yl2>1))>0.5)=-inf; 
  [Yl2,YD] = vbm_vol_downcut(Yl2,Ym,-0.00); Yl2(Yl2==-inf | YD>50)=0;
  Yl2 = vbm_vol_median3c(Yl2,YM,YM);

  Yl1((Yl2==opt.LAB.CT(1))>0.55 & Ym<2.9 & Yl1==0)=opt.LAB.CT(1);
  Yl1((Yl2==opt.LAB.CB(1))>0.55 & Ym<2.9 & Yl1==0)=opt.LAB.CB(1);
  
  %YM = Yl1==opt.LAB.BG(1) & ~vbm_vol_morph(Yl1==opt.LAB.BG(1),'o',1);
  %Yl1(YM & Ym>2.85 & Yl1Ans~=opt.LAB.BG(1))=1;  Yl1(YM & Ym<=2.5 & Yl1Ans~=opt.LAB.BG(1))=1;
  
  Yl1(vbm_vol_morph(Yl1==opt.LAB.BG(1),'lc') & Ym>1.75 & Ym<2.75 & Yg<0.2)=opt.LAB.BG(1);
  Yl1(vbm_vol_morph(Yl1==opt.LAB.TH(1),'lc') & Ym>1.75 & Ym<2.75 & Yg<0.2)=opt.LAB.TH(1);
  
  
  % head
  Yl1(Yl0==2 & Yl1==0 & Ym>1.2 & Yp0==0 & Yl1Ans==opt.LAB.HD(1)) = opt.LAB.HD(1);
  Yl1(Yl1==0 & Ym>1.2 & Yl0==1 & Yp0==0 &...
    (Yl1Ans==opt.LAB.HD(1) | Yl1Ans==opt.LAB.BV(1) | Yl1Ans==opt.LAB.NV(1))) = opt.LAB.HD(1);
  Yl1(Yl1==opt.LAB.HD(1) & smooth3(Yl1==opt.LAB.HD(1))<0.5)=0;
  
  % ON
  Yl1(Yl1==opt.LAB.ON(1) & smooth3(Yl1Ans==opt.LAB.ON(1))<0.5)=0;
  
  
  % CB
  Yl2=Yl1; Yl2(Yl1==0 & (Yl1Ans~=opt.LAB.CB(1) | Ym<1.9))=-inf;
  [Yl2,YD] = vbm_vol_downcut(Yl2,Ym,0.1); Yl2(Yl2==-inf | YD>100)=0; 
  Yl2(Yl1==0 & (Yl1Ans~=opt.LAB.CB(1) | Ym<1.25))=-inf;
  [Yl2,YD] = vbm_vol_downcut(Yl2,Ym,-0.2); Yl2(Yl2==-inf | YD>100)=0;
  Yl1(Yl1<=0 & smooth3(Yl2==opt.LAB.CB(1) & Yl1Ans~=opt.LAB.BV(1))>0.5)=opt.LAB.CB(1);
  

  %% blood vessel detection
  YWMF = single(smooth3(vbm_vol_morph(Ym>2.5 & Ym<3.5 & Yp0>1.5 & Yl1Ans~=opt.LAB.BV(1),'lo',0))>0.5);
  YWMF(YWMF==0 & (Ym<=2.5 | Ym>=3.5)) = -inf; [YWMF,YD] = vbm_vol_simgrow(YWMF,Ym,0.1); YWMF(isinf(YWMF) | YD>50)=0; 
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
  
  
  %% low intensity structures
  % WM Hyperintensities can be problematic (OASIS0031), so you have to
  % use the Yl1A!
  YM = Yl1==0 & smooth3(YBV==0 & Yl1==0 & (Ym>1.25 & Ym<1.75) & ...
    smooth3(Yl1A==opt.LAB.HD(1) | Yl1A==opt.LAB.BV(1) | Yl1A==opt.LAB.NV(1))>0.5 & ... 
    ((Yp0>1.25 & Yp0<1.75) | Yp0<0.5)  & ~YWMF & ...
    (Ym<1.6 | Yp0==0 | Yl1Ans==opt.LAB.BV(1)) & ...
    (Yp0==0 | Yl1Ans==opt.LAB.BV(1) | Yl1Ans==opt.LAB.NV(1)))>0.55; 
  YBV(YM) = 2;
  
  YM = vbm_vol_morph(smooth3(Ym>1.25 & Ym<1.75)>0.5 & ...
    smooth3(Yl1A==opt.LAB.HD(1) | Yl1A==opt.LAB.BV(1) | Yl1A==opt.LAB.NV(1))>0.5,'o',2);
  YBV(YM) = 2;
  
  Yl1(YBV>0) =  opt.LAB.BV(1);  
    
  
  
  %% alignment of the ventricels
  for s=1:2
    Yl1(vbm_vol_morph(Yl1A==opt.LAB.VT(s) & Ym<1.8,'lab') & Yp0<1.8) = opt.LAB.VT(1);        % ventricle
  end
  VT=single(Yl1==opt.LAB.VT(1)); VT((Yl1Ans==opt.LAB.NV(1) & Yp0<2) | (Yl1Ans==opt.LAB.BV(1) & Yp0<2))=2;
  VT((Yp0>2.0 | Yp0<0.5 | Yl1>0) & ~VT) = -inf; VT = vbm_vol_simgrow(VT,Ym,1.5); VT = vbm_vol_simgrow(VT,Ym,1.5); 
  Yl1(smooth3(VT==1)>0.2 & (Ym<1.75 | Yp0<1.75 | (Yl1==0 & Ym<1.75)))=opt.LAB.VT(1);
  
  % alignment of HD
  Yl1(Yl1==0 & Ym>2.5 & vbm_vol_morph(Yp0==0,'lc')) = opt.LAB.HD(1); 
  
  
    
  %% region-growing in GM only for non-blood vessels regions
  YM = Yl1==0 & Ym>0.9 & Yp0>0.5 & vbm_vol_morph(Yl1>0 & Yl1<20,'lc'); Yl1(YM)=Yl1(YM); clear YM;
  Yl2=Yl1; Yl2(Yl2>1 | Ym<2.1 | Yl0==3)=-inf; 
  Yl2=single(smooth3(vbm_vol_downcut(Yl2,Ym,0.01,resTr.vx_volr))>0.5); 
  Yl2(Yl2>1 | Ym<1.2 | Yl0==3)=-inf; 
  Yl2=smooth3(vbm_vol_downcut(Yl2,Ym,-.1,resTr.vx_volr))>0.5;
  Yl1(Yl1==0 & Yl2 & Ym>0.9)=opt.LAB.CT(1); Yl1(isinf(Yl1))=0;

  Yl1((Yl1==0 | Yl1==opt.LAB.BV(1)) & Yp0<0.5 & vbm_vol_morph(Yl1==opt.LAB.HD(1),'lc',2))=opt.LAB.HD(1);
  Yl1(Yl1==0 & (Ym<0.9 | ~vbm_vol_morph(Yl1>0,'lc',2)))=-inf; Yl1=vbm_vol_downcut(Yl1,Ym,-.1,resTr.vx_volr);
  Yl1=vbm_vol_median3c(Yl1,Yl1==0,Yl1>0);
  
  Yl1 = vbm_vol_median3c(Yl1,Yp0>0 & Yl1~=opt.LAB.BV(1));
  
  
  
  %% prepare mask for filling of subcortical regions
  YB0 = vbm_vol_morph(Yp0>0,'lc',2); 
  YMF = Yl1==opt.LAB.BG(1) | Yl1==opt.LAB.VT(1) | Yl1==opt.LAB.TH(1) | ...
        Yl1Ans==opt.LAB.BG(1) | Yl1Ans==opt.LAB.VT(1) | Yl1Ans==opt.LAB.TH(1) | Yl1Ans==opt.LAB.HI(1);
  Ymf = max(Ym,YMF*3); Ymf(YMF)=min(3,Ymf(YMF)); Ymfs = vbm_vol_smooth3X(Ymf,1); 
  YM  = vbm_vol_morph(YMF,'d',3) & Ymfs>2.3;
  Ymf(YM) = max(min(Ym(YM),3),Ymfs(YM)); 
  clear Ymfs YM; 
  
  
  
  %% corrections
  % correct CSF in special brain regions
  YM  = (Yl1==opt.LAB.BG(1) | Yl1==opt.LAB.TH(1) | Yl1==opt.LAB.HC(1));
  Yl1(YM & Yl1Ans==opt.LAB.CT(1) & Yp0<1.5) = opt.LAB.CT(1);    
  Yl1(YM & Yl1Ans==opt.LAB.VT(1) & Yp0<1.5) = opt.LAB.VT(1);
  Yl1(Yl1==opt.LAB.BV(1) & ~YB0) = opt.LAB.HD(1);
  clear YM; 
    
  
  % corrections for non brain parts
  Yl1(smooth3(~YB0 & Ym>1.5 & (Yl1==0 | Yl1==opt.LAB.BV(1)))>0.5) = opt.LAB.HD(1); 
   
  
  % cleanup for labels that are not so excact in this version of Yl1A
  Yl1(Yl1==opt.LAB.HC(1))=opt.LAB.CT(1);
  Yl1=vbm_vol_median3c(Yl1,Yl1>=0 & Yl1<5,Yl1>=0 & Yl1<5);
  Yl1=vbm_vol_median3c(Yl1,Yl1>=0 & Yl1<5,Yl1>=0 & Yl1<5);
  Yl1(isinf(Yl1))=0;
  
  
  %% brain mask
  % we have to close now some CSF areas (hull mask) and we have to label
  % high intensity tissues as BV
  Yb = Yl1>0 & Yl1<opt.LAB.HD(1) & Yl1~=opt.LAB.BV(1) & Ym>1.2;
  Yb = vbm_vol_morph(smooth3(Yb)>0.4,'o',1);
  Yb = vbm_vol_morph(smooth3(Yb)>0.4,'ldc',8);

  % new brain areas (previous head)
  YM1 = smooth3((Yl1<=0 | Yl1==opt.LAB.HD(1)) & Yb & (Ym>1.25))>0.5; 
  YM2 = (Yl1<=0 | Yl1==opt.LAB.HD(1)) & Yb & ~YM1; 
  [YD,YI] = vbdist(single(Yl1>0 & Yl1~=opt.LAB.HD(1)),YM2); 
  Yl1(YM2(:)) = Yl1(YI(YM2(:))); Yl1(YM1) = opt.LAB.BV(1);
  for i=1:3, Yl1 = vbm_vol_median3c(Yl1,YM1 | YM2);  end
  clear YM1 YM2;
    
  % new/old head areas
  Yl1(Yl1~=opt.LAB.HD(1) & ~YB0 & ~Yb & Yl1>0) = opt.LAB.HD(1);   
  
  % if there was no head at the beginning (skull-stripped input) than we
  % want no head in the labelmap (the brainmask still removes these areas)
  if ER.HD(1), Yl1(Yl1==opt.LAB.HD(1))=0; end
  

  
  %% side aligment using laplace to correct for missalignments due to the normalization
  d = 5; 
  YM = vbm_vol_smooth3X(single(vbm_vol_morph(mod(Yl1A,2)==0,'distdilate',d,resTr.vx_volr)) & ...
      single(vbm_vol_morph(mod(Yl1A,2)==1,'distdilate',d,resTr.vx_volr)==1),20);
  Ys = 2*single(mod(Yl1A,2)==0 & Yl1A>0 & Ym>2.5 & YM<max(YM(:))*0.9) + ...
         single(mod(Yl1A,2)==1 & Yl1A>0 & Ym>2.5 & YM<max(YM(:))*0.9);
  Ys(Ymf<=2.5)=-inf; Ys(Ys==0)=1.5; clear YM;

  Ysr=vbm_vol_resize(Ys,'reduce'); Ysr=round(Ysr*2)/2; 
  Ysr=vbm_vol_laplace3R(Ysr,Ysr==1.5,0.001); 
  Ysr=vbm_vol_resize(single(Ysr),'dereduce',size(Ym)); 
  
  Ys(Ys==1.5)=round(Ysr(Ys==1.5)); clear Ysr; 
  Ys(isinf(Ys) & Ymf>2)=0; Ys=vbm_vol_downcut(Ys,Ymf,3,resTr.vx_volr); 
  Ys(isinf(Ys) & Ymf>0)=0; Ys=vbm_vol_downcut(Ys,Ymf,3,resTr.vx_volr); 
  Ys(Ys<=0)=2-mod(Yl1A(Ys<=0),2);
  Ys=round(vbm_vol_smooth3X(Ys,2));
  Yl1(Yl1>0)=Yl1(Yl1>0)+(Ys(Yl1>0)==2); 
  Yl1(Yl1<0.5)=0;
  clear Ys d;

  
  
  %% back to original resolution and full size
  Yp0 = vbm_vol_resize(Yp0,'dereduceV',resTr);
  Yl1 = vbm_vol_resize(Yl1,'dereduceV',resTr,'nearest');
  Yb  = vbm_vol_resize(single(Yb) ,'dereduceV',resTr);
  YMF = vbm_vol_resize(single(YMF),'dereduceV',resTr);
  
  Yl1 = vbm_vol_resize(Yl1,'dereduceBrain',BB);
  Yp0 = vbm_vol_resize(Yp0,'dereduceBrain',BB);
  Yb  = vbm_vol_resize(Yb ,'dereduceBrain',BB)>0.5;
  YMF = vbm_vol_resize(YMF,'dereduceBrain',BB)>0.5;
  
  Ym  = mgTO; clear mgTO;
  
 
  
  %% setting head label for voxel outside the bounding box
  Yl1(smooth3(Yp0<1 & Ym>1.5 & (Yl1==0 | Yl1==opt.LAB.BV(1)))>0.5) = opt.LAB.HD(1);
  [HDr,resTr] = vbm_vol_resize(Yl1>0,'reduceV',vx_vol,3,64);
  HDr = smooth3(vbm_vol_morph(HDr>0.1,'lc',2))>0.5;
  HD  = smooth3(Yl1>0 | vbm_vol_resize(HDr,'dereduceV',resTr))>0.5;
  Yl1(HD & Yl1<=0) = opt.LAB.HD(1); 
  Yl1(~HD) = 0;

  
  % if there was no head at the beginning (skull-stripped input) than we
  % want no head in the labelmap (the brainmask still removes these areas)
  if ER.HD(1), Yl1(Yl1==opt.LAB.HD(1))=0; end
  
  
  
  %% subvolumes
  fn = setdiff(fieldnames(opt.LAB),{'NV','HD','NB','ON','HC'}); vol=struct();
  for fni=1:numel(fn)
    eval(sprintf(['vol.vol_abs_%s = prod(vx_vol)/1000 .* ' ...
      '[sum(Yl1(:)==%d) sum(Yl1(:)==%d)];'],fn{fni},def.LAB.(fn{fni}))); 
  end 
  vol.vol_LRP = sum( mod(Yl1(:),2)==1 & Yl1(:)>0 & Yl1(:)<20) ./ ...
                sum( Yl1(:)>0 & Yl1(:)<20);

end

function [vol,Yl1,Yb,YMF] = vbm_vol_partvol(Yl1A,Yp0,Ym,Yl0,Yb,opt)
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
  def.res    = 2.2;
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
  def.action = 'full';
  
  opt = checkinopt(opt,def);
  
  

  
  vx_vol = opt.vx_vol; %vx=mean(vx_vol);
    

  %% Initialization:
  % ds('l2','',vx_vol,Ym,Yl1A,Ym,Yp0/3,130)
  
  % Optimization:
  Ym = Ym*3; mgTO = Ym;
  [Yl1A,Yl0,BB] = vbm_vol_resize({Yl1A,Yl0} ,'reduceBrain',vx_vol,10,Yp0>0);   % removing of background
  [Ym,Yp0,Yb]   = vbm_vol_resize({Ym,Yp0,Yb},'reduceBrain',vx_vol,10,Yp0>0);  % removing of background
  [Ym,Yp0,Yb,resTr]    = vbm_vol_resize({Ym,Yp0,Yb},'reduceV',vx_vol,opt.res,64); Yb=Yb>0.5; 
  [Yl1A,Yl0]           = vbm_vol_resize({Yl1A,Yl0} ,'reduceV',vx_vol,opt.res,64,'nearest'); 
  
  % gradients and divergence 
  [gx,gy,gz]=vbm_vol_gradient3(Ym);        Yg   = abs(gx)+abs(gy)+abs(gz); Yg=Yg./Ym;  clear gx gy gz;

  

  
  %% Find major structures  
  % alignment of high intensity structures 
  warning off 'MATLAB:vbm_vol_morph:NoObject'; 

 %tic
  [Ymr,Yp0r,Ygr,resTrr] = vbm_vol_resize({Ym,Yp0,Yg},'reduceV',vx_vol,max(2.2,opt.res),64); 
  [Yl1Ar,Yl0r]          = vbm_vol_resize({Yl1A,Yl0} ,'reduceV',vx_vol,max(2.2,opt.res),64,'nearest'); 

  Yl1Ar = uint8(Yl1Ar); Yl0r = uint8(Yl0r); Yl1Ansr = uint8(round(single(Yl1Ar)/2)*2-1);

  [gx,gy,gz]=vbm_vol_gradient3(max(2,Ymr)); Ydivr = smooth3(divergence(gy,gx,gz)); clear gx gy gz;

  Yl1r = zeros(size(Ymr),'single'); 
  YMr  = vbm_vol_morph(vbm_vol_morph(smooth3(Yl0r==3 & Ymr<0.4)>0.5,'lo',1),'lc',1);
  Yl1r(YMr)=255;                                                                       % backgound

  YMr=Ymr>2.25 & (Yl0r==2 | Yl0r==1) & ~vbm_vol_morph(Yp0r>0.5,'lc');                  % head
  Yl1r(YMr) = opt.LAB.HD(1); ER.HD(1:2)=sum(YMr(:))==0; 
  for s=1:2
    YMr=vbm_vol_morph((Ymr>2.5 & Yp0r>2.5) & Ymr<3.25 & Yl1Ar==opt.LAB.CT(s),'lc');    % cerebrum
    Yl1r(YMr)=opt.LAB.CT(1); ER.CT(s)=sum(YMr(:))==0; 
    YMr=vbm_vol_morph((Ymr>2.5 | Yp0r>2.2) & Ymr<3.25 & Yl1Ar==opt.LAB.BS(s),'lc');    % brainstem
    Yl1r(YMr)=opt.LAB.BS(1); ER.BS(s)=sum(YMr(:))==0;  
    YMr=vbm_vol_morph((Ymr>2.5 | Yp0r>2.5) & Ymr<3.25 & Yl1Ar==opt.LAB.MB(s),'lc');    % midbrain
    Yl1r(YMr)=opt.LAB.MB(1); ER.MB(s)=sum(YMr(:))==0;  
    YMr=vbm_vol_morph((Ymr>2.5 | Yp0r>2.5) & Ymr<3.25 & Yl1Ar==opt.LAB.CB(s),'lc');    % cerebellum
    Yl1r(YMr)=opt.LAB.CB(1); ER.CB(s)=sum(YMr(:))==0; 
    YMr=vbm_vol_morph((Ymr>2.1 | Yp0r>2.1) & Ymr<3.50 & Yl1Ar==opt.LAB.ON(s),'lab');   % optical nerv
    Yl1r(YMr)=opt.LAB.ON(1); %ER.ON(s)=sum(YMr(:))==0; 
    YMr=vbm_vol_morph(Yl1Ar==opt.LAB.HC(s) & Ymr>1.75 & Ymr<2.25,'lab');               % hippocampus
    Yl1r(YMr)=opt.LAB.HC(1); %ER.HC(s)=sum(YMr(:))==0; 
  end
  warning off 'MATLAB:vbm_vol_morph:NoObject';

  % claustrum - to reduce BG overgrowing
  YBGD = vbdist(single(Yl1Ansr==opt.LAB.BG(1)),Yp0r>0,resTrr.vx_vol);
  Yl1r(smooth3(YBGD==0 & Ymr>2.8 & Ymr<3.2)>0.5)=opt.LAB.CT(1); 
  Yl1r((YBGD<4 & YBGD>2 & Ydivr<0.1 & Ygr>0.02 & Ymr>2.3 & Ymr<2.8 & Yl1r==0 & Yl1Ansr==opt.LAB.CT(1))>0.5)=opt.LAB.CT(1);
  clear YBGD;

  % region-growing for special high intensity regions
  Yl1r(((Ymr<=2.85 | (Yl1Ansr~=opt.LAB.CT(1) & Ydivr>0.05 & Ymr<=2.9)| Yl0r==3) & Yl1r==0) )=-inf;
  [Yl1r,YD] = vbm_vol_simgrow(Yl1r,Ymr,0.05); Yl1r(isinf(Yl1r) | YD>0.1)=0; 

  % alignment of medium intensity structures
  for s=1:2
    YMr=vbm_vol_morph(Yl1Ar==opt.LAB.BG(s) & Yl1r==0 & Ymr>1.75 & Ymr<2.75 & Ygr<0.2 & Ydivr>0,'lab'); % basal ganglia
    YMr=vbm_vol_morph(YMr,'lc',2) & Ymr<2.85; Yl1r(YMr)=opt.LAB.BG(1); ER.BG(s)=sum(YMr(:))==0;  
    YMr=vbm_vol_morph(Yl1Ar==opt.LAB.TH(s) & Yl1r==0 & Ymr>1.75 & Ymr<2.75 & Ygr<0.2,'lo',0);         % hypothalamus
    YMr=vbm_vol_morph(YMr,'lc',2) & Ymr<2.85; Yl1r(YMr)=opt.LAB.TH(1); ER.TH(s)=sum(YMr(:))==0;          
  end



  % Error managment, if we miss a major structure
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



  % refinement of brainstem and midbrain
  Yl1r(Yl1r==opt.LAB.BS(1) & ~vbm_vol_morph(Yl1r==opt.LAB.BS(1) & Ymr>2.25,'lo',1))=0;
  Yl1r(Yl1r==opt.LAB.MB(1) & ~vbm_vol_morph(Yl1r==opt.LAB.MB(1) & Ymr>2.25,'lo',1))=0;


  % refinement of hypothamlamus
  YMr = Yl1r==opt.LAB.TH(1) & ~vbm_vol_morph(Yl1r==opt.LAB.TH(1),'o',1);
  Yl1r(YMr & Ymr>2.75)=1;  Yl1r(YMr & Ymr<=2.5)=0;


  % refinement of hippocampus
  Yl1r((Yl1Ansr~=opt.LAB.HC(1) | Ydivr<-0.02 | Ymr>2.25) & Yl1r==opt.LAB.HC(1)) = 1; 
  Yl1r = vbm_vol_median3c(Yl1r,Yl1r==opt.LAB.HC(1));
  Yl1r(vbm_vol_morph(Yl1r==opt.LAB.HC(1) & Yl1Ansr==opt.LAB.HC(1),'lc',1))=opt.LAB.HC(1);
  for s=0:1, 
    Yl1r(mod(Yl1Ar,2)==s & Yl1r==opt.LAB.HC(1) & ...
      ~vbm_vol_morph(mod(Yl1Ar,2)==s & Yl1r==opt.LAB.HC(1),'l')) = 0; 
  end


  % Hyperintensities
  YH = single(vbm_vol_morph(Yl1Ar==opt.LAB.CT(1) & Ymr>1.5 & (Yp0r - Ymr)>0.5,'o',1));
  YH(smooth3(Yl1r~=0 | (Yp0r - Ymr)<0.3 | Ymr<1.5)>0.5) = -inf;
  YH = vbm_vol_simgrow(YH,Ymr,0.04);
  Yl1r(YH==1) = opt.LAB.HI(1); clear YH;


  % refinement of the basal ganglia
  Yl1r( ( (Ydivr<0.00 & Yl1Ansr~=opt.LAB.BG(1)) | (Ydivr<-0.05 & Yl1Ansr==opt.LAB.BG(1)) | Ygr>0.1) & ...
    Yl1r==opt.LAB.BG(1) & Ymr>1.75 & Ymr<2.5 & Yp0r>0.5) = 0;                                     % basasl ganglia corr


 
  Yl1  = vbm_vol_resize(Yl1r,'dereduceV',resTrr,'nearest');
  Ydiv = vbm_vol_resize(Ydivr,'dereduceV',resTrr);
  clear Yl0r Yp0r Ym0r Yg0r Yl1Ar Yl1nsr Yl1r;
  %toc, tic  


  %% original resolution
  Yl1(Yl1>0 & Yp0>0 & ((Ym>2.60 & Ym<2.9) | Ym>3.3))=0; 
  Yl1(Yl1>0 & Yp0>0 & Ym<2.60 & ~(Yl1==255 | Yl1==opt.LAB.CT(1) | Yl1==opt.LAB.TH(1) | Yl1==opt.LAB.BG(1) | Yl1==opt.LAB.HC(1)))=0; 

  Yl1(Yl1==0 & (Ym<2.1 | Ym>2.9))=-inf;
  [Yl1,YD] = vbm_vol_simgrow(Yl1,Ym,0.1); Yl1(isinf(Yl1) | YD>0.2)=0;
  Yl1  = vbm_vol_median3c(Yl1,Yp0>1,Yp0>1);

  Yl1Ans = round(Yl1A/2)*2-1;


  % further region-growing (GM-WM) without ON!
  Yl2 = Yl1; Yl2(Yl2==opt.LAB.ON(1) | smooth3(Yl1==0 & (Yl2==opt.LAB.BG(1) | (Yl1Ans==opt.LAB.BG(1) & Yl1==0) | ...
    Ym<=1.9 | (Yl1Ans~=opt.LAB.CT(1) & Ydiv>-0.1 & Ym<=2.9) | Ym>3.5))>0.5) = -inf; 
  Yl2 = vbm_vol_simgrow(Yl2,Ym,0.1); [Yl2,YD] = vbm_vol_simgrow(Yl2,Ym,0.1);  Yl2(isinf(Yl2) | YD>0.1)=0;
  Yl1 = max(Yl1,Yl2); clear Yl2;

  % regino growing (GM)
  YM = Yp0>0 & Yl1~=opt.LAB.BG(1) & Yl1~=opt.LAB.ON(1) & Ym<2.5 & Ym>1.5; 
  Yl2=Yl1; Yl2(smooth3(Yl2==0 & ( Ym<2.0 |  Ym>3.2 | Yl2>1 | Yl1Ans==opt.LAB.BV(1)))>0.5)=-inf; 
  [Yl2,YD] = vbm_vol_downcut(Yl2,Ym,0.02); Yl2(Yl2==-inf | YD>50)=0;
  Yl2 = vbm_vol_median3c(Yl2,YM,YM); 

  Yl2(smooth3(Yl2==0 & ( Ym<1.5 |  Ym>3.2 | Yl2>1 | Yl1Ans==opt.LAB.BV(1)))>0.5)=-inf; 
  [Yl2,YD] = vbm_vol_downcut(Yl2,Ym,-0.00); Yl2(Yl2==-inf | YD>50)=0;
  Yl2 = vbm_vol_median3c(Yl2,YM,YM);

  Yl1((Yl2==opt.LAB.CT(1))>0.55 & Ym<3.1 & Yl1==0)=opt.LAB.CT(1);
  Yl1((Yl2==opt.LAB.CB(1))>0.55 & Ym<3.1 & Yl1==0)=opt.LAB.CB(1);

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
  Yl2=Yl1; Yl2(Yl1==0 & (Yl1Ans~=opt.LAB.CB(1) | Ym<1.8 | Yp0==0))=-inf;
  [Yl2,YD] = vbm_vol_downcut(Yl2,Ym,0.1); Yl2(Yl2==-inf | YD>100)=0;
  Yl2(Yl2==opt.LAB.CT(1) & Yl1~=opt.LAB.CT(1))=0; 
  Yl2(Yl1==0 & (Ym<1.8 | Yp0==0 | Yl1Ans==opt.LAB.BV(1) ))=-inf;
  [Yl2,YD] = vbm_vol_downcut(Yl2,Ym,-0.05); Yl2(Yl2==-inf | YD>100)=0;
  Yl1(Yl1<=0 & smooth3(Yl2==opt.LAB.CB(1) & Yl1Ans~=opt.LAB.BV(1))>0.5 & Ym<2.8)=opt.LAB.CB(1);
  Yl1(Yl1<=0 & smooth3(Yl2==opt.LAB.CT(1) & Yl1Ans~=opt.LAB.BV(1))>0.5 & Ym<2.8)=opt.LAB.CT(1);
%toc, tic


  %% alignment of the ventricels
  for s=1:2
    Yl1(vbm_vol_morph(Yl1A==opt.LAB.VT(s) & Ym<1.8,'lab') & Yp0<1.8) = opt.LAB.VT(1);        % ventricle
  end

  [Yt,resT3] = vbm_vol_resize(Yb & Yl1Ans~=opt.LAB.CB(1),'reduceV',vx_vol,6,32); Yt = vbm_vol_morph(Yt,'e',2);
  Yt = vbm_vol_resize(vbm_vol_smooth3X(Yt),'dereduceV',resT3)>0.5; clear resT3;

  VT=single(Yl1==opt.LAB.VT(1)); 
  VT(smooth3((Yl1Ans==opt.LAB.NV(1) | Yl1Ans==opt.LAB.BV(1)) & Ym>1.25 & Ym<1.9 )>0.5 | ~Yt)=2; clear Yt;
  VT(((Yp0>2.0 & Ym>2.0) | Yl1>0) & ~VT) = -inf; VT = vbm_vol_simgrow(VT,Ym,1.5); 
  Yl1(smooth3(VT==1)>0.2 & (Ym<1.75 | Yp0<1.75 | (Yl1==0 & Ym<1.75)))=opt.LAB.VT(1); clear VT;

  % alignment of HD
  Yl1(Yl1==0 & Ym>2.5 & vbm_vol_morph(Yp0==0,'lc')) = opt.LAB.HD(1); 


  % blood vessel detection
  YWMF = vbm_vol_morph(Yl1>0,'e') | vbm_vol_morph( Yl1Ans==opt.LAB.HC(1) | Yl1Ans==opt.LAB.BS(1) | Yl1Ans==opt.LAB.ON(1) | ...
          Yl1Ans==opt.LAB.BG(1) | Yl1Ans==opt.LAB.TH(1) | Yl1Ans==opt.LAB.MB(1) | ...
          Yl1Ans==opt.LAB.CB(1) | Yl1Ans==opt.LAB.CB(1) | Yl1Ans==opt.LAB.VT(1),'dilate',ceil(1/mean(vx_vol)));

  % high intensity structures
  YBV = zeros(size(Ym),'uint8'); 
  YBV(Yp0==0 & vbm_vol_morph(Yp0>0,'lc') & Ym>2.2)=3;
  YBV(YBV==0 & Yl1==0 & Ym>2.3 & Yp0>2.3 & ~YWMF & ((Yl1==0 & Yl1Ans==opt.LAB.CT(1)) | ...
    Yp0==0 | Yl1Ans==opt.LAB.NV(1) | Yl1Ans==opt.LAB.BV(1) )) = 3;
  YBV(vbm_vol_morph(YBV==3,'d',1) & Ym<2.8 & Ym>2.3) = 3; 
  YBV(vbm_vol_morph(YBV==3,'d',2) & ((Yl1==0 & Ym>2.1) | Ym>3.2)) = 3; 


  % low intensity structures
  % WM Hyperintensities can be problematic (OASIS0031), so you have to
  % use the Yl1A!
  YM = Yl1==0 & smooth3(YBV==0 & Yl1==0 & (Ym>1.25 & Ym<1.75) & ...
    smooth3(Yl1A==opt.LAB.HD(1) | Yl1A==opt.LAB.BV(1) | Yl1A==opt.LAB.NV(1))>0.5 & ... 
    ((Yp0>1.25 & Yp0<1.75) | Yb)  & ~YWMF & ...
    (Ym<1.6 | Yp0==0 | Yl1Ans==opt.LAB.BV(1)) & ...
    (Yp0==0 | Yl1Ans==opt.LAB.BV(1) | Yl1Ans==opt.LAB.NV(1)))>0.55; 
  YBV(YM) = 2;

  YM = vbm_vol_morph(smooth3(Ym>1.25 & Ym<1.75)>0.5 & ...
    smooth3(Yl1A==opt.LAB.HD(1) | Yl1A==opt.LAB.BV(1) | Yl1A==opt.LAB.NV(1))>0.5,'o',2);
  YBV(YM) = 2;

  Yl1(YBV>0) =  opt.LAB.BV(1);  
%toc    



  %% region-growing in GM only for non-blood vessels regions
%tic
  Yl1((Yl1==0 | Yl1==opt.LAB.BV(1)) & Yb & vbm_vol_morph(Yl1==opt.LAB.HD(1),'lc',2))=opt.LAB.HD(1);
  Yl1(Yl1==0 & (Ym<0.75 & ~Yb))=-inf; 
  [Yl1,YD]=vbm_vol_downcut(Yl1,Ym,-.1,resTr.vx_volr); Yl1(isinf(Yl1) | YD>20)=0;
  Yl1=vbm_vol_median3c(Yl1,Yl1~=opt.LAB.BV(1) & Yp0>0,Yp0>0);


  %% prepare mask for filling of subcortical regions
  [Yt,Yt2,resT3] = vbm_vol_resize({Yb & Yl1Ans~=opt.LAB.CB(1),Yl1==opt.LAB.VT(1)},'reduceV',vx_vol,6,32); 
  Yt = vbm_vol_morph(Yt,'e',2); Yt2 = vbm_vol_morph(Yt2,'d',2);
  [Yt,Yt2] = vbm_vol_resize({vbm_vol_smooth3X(Yt),vbm_vol_smooth3X(Yt2)},'dereduceV',resT3); clear resT3;
  Yt=Yt>0.5; Yt2=Yt2>0.5;

  res=2;
  [Ymr,Yp0r,resT3] = vbm_vol_resize({Ym,Yp0},'reduceV',vx_vol,res,32);
  [YMF1r,YMF2r] = vbm_vol_resize({single(((~Yt | ~Yt2) & Ym<2 & Yp0<2) | ~Yb), ...
      single(Yl1==opt.LAB.BG(1) | Yl1==opt.LAB.VT(1) | Yl1==opt.LAB.TH(1) | ...
      Yl1Ans==opt.LAB.BG(1) | Yl1Ans==opt.LAB.VT(1) | ...
      Yl1Ans==opt.LAB.TH(1) | Yl1Ans==opt.LAB.HI(1))}, ...
      'reduceV',vx_vol,res,32); YMF1r=YMF1r>0.5; YMF2r=YMF2r>0.5; %clear Yt Yt2;
  YMFr=1.5 * ones(size(Ymr),'single'); YMFr(Ymr>2.5 & Yp0r>2.5)=-inf; 
  YMFr(YMF1r)=1; YMFr(YMF2r)=2; % clear YMF1r YMF2r;
  YMFr = vbm_vol_laplace3R(YMFr,YMFr==1.5,0.01); YMFr(isinf(YMFr))=1.5;  
  YMF  = vbm_vol_resize(YMFr,'dereduceV',resT3)>1.5; clear resT3;
  Yl1(smooth3(Yl1==opt.LAB.CT(1) & YMF & Ym<2.5 & Ym>1.5)>0.5)=opt.LAB.HI(1);
  %%   

  Ymf = max(Ym,YMF*3); Ymf(YMF)=min(3,Ymf(YMF)); Ymfs = vbm_vol_smooth3X(Ymf,1); 
  YM  = vbm_vol_morph(YMF,'d',3) & Ymfs>2.3;
  Ymf(YM) = max(min(Ym(YM),3),Ymfs(YM)); 
  clear Ymfs YM; 


  % corrections
  % correct CSF in special brain regions
  YM  = (Yl1==opt.LAB.BG(1) | Yl1==opt.LAB.TH(1) | Yl1==opt.LAB.HC(1));
  Yl1(YM & Yl1Ans==opt.LAB.CT(1) & Yp0<1.5) = opt.LAB.CT(1);    
  Yl1(YM & Yl1Ans==opt.LAB.VT(1) & Yp0<1.5) = opt.LAB.VT(1);
  Yl1(Yl1==opt.LAB.BV(1) & Yb) = opt.LAB.HD(1); %%%%%%
  clear YM; 


  % corrections for non brain parts
  Yl1(smooth3(Yp0==0 & Ym>1.5 & (Yl1==0 | Yl1==opt.LAB.BV(1)))>0.5) = opt.LAB.HD(1); 


  % cleanup for labels that are not so excact in this version of Yl1A
  Yl1(Yl1==opt.LAB.HC(1))=opt.LAB.CT(1);
  Yl1=vbm_vol_median3c(Yl1,Yp0>0 & Yl1<5 & Ym<2.2,Yp0>0 & Yl1<5 & Ym<2.2);
  Yl1(isinf(Yl1))=0;
%toc 

  %% brain mask
  % we have to close now some CSF areas (hull mask) and we have to label
  % high intensity tissues as BV
  Yb = Yl1>0 & Yl1<opt.LAB.HD(1) & Yl1~=opt.LAB.BV(1) & Ym>1.2;
  Yb = vbm_vol_morph(smooth3(Yb)>0.4,'o',1);
  [Yb,resT2] = vbm_vol_resize(Yb,'reduceV',vx_vol,2,32); Yb = vbm_vol_morph(Yb,'lc',4);
  Yb = vbm_vol_resize(vbm_vol_smooth3X(Yb),'dereduceV',resT2)>0.4; 

  % new brain areas (previous head)
  YM1 = smooth3((Yl1<=0 | Yl1==opt.LAB.HD(1)) & Yb & (Ym>1.25))>0.5; 
  YM2 = (Yl1<=0 | Yl1==opt.LAB.HD(1)) & Yb & ~YM1; 
  [YD,YI] = vbdist(single(Yl1>0 & Yl1~=opt.LAB.HD(1)),YM2); 
  Yl1(YM2(:)) = Yl1(YI(YM2(:))); Yl1(YM1) = opt.LAB.BV(1);
  for i=1:3, Yl1 = vbm_vol_median3c(Yl1,YM1 | YM2);  end
  clear YM1 YM2;

  % new/old head areas
  Yl1(Yl1~=opt.LAB.HD(1) & ~Yb & Yl1>0 & Yl1<255) = opt.LAB.HD(1);   

  % if there was no head at the beginning (skull-stripped input) than we
  % want no head in the labelmap (the brainmask still removes these areas)
  if ER.HD(1), Yl1(Yl1==opt.LAB.HD(1))=0; end



  %% side aligment using laplace to correct for missalignments due to the normalization
  [Ymr,Ymfr,resT3] = vbm_vol_resize({Ym,Ymf},'reduceV',vx_vol,max(1.2,opt.res),64); 
  [Yl1Ar]           = vbm_vol_resize({Yl1A,Yl1} ,'reduceV',vx_vol,max(1.2,opt.res),64,'nearest'); 

  d = 3; 
  YMr = vbm_vol_smooth3X(single(vbm_vol_morph(mod(Yl1Ar,2)==0,'distdilate',d,resTr.vx_volr)) & ...
      single(vbm_vol_morph(mod(Yl1Ar,2)==1,'distdilate',d,resTr.vx_volr)==1),10);
  Ysr = 2*single(mod(Yl1Ar,2)==0 & Yl1Ar>0 & Ymr>2.5 & YMr<max(YMr(:))*0.9) + ...
          single(mod(Yl1Ar,2)==1 & Yl1Ar>0 & Ymr>2.5 & YMr<max(YMr(:))*0.9);
  Ysr(Ymfr<=2.5)=-inf; Ysr(Ysr==0)=1.5; clear YMr;

  Ysrr=vbm_vol_resize(Ysr,'reduce'); Ysr=round(Ysr*2)/2; 
  Ysrr=vbm_vol_laplace3R(Ysrr,Ysrr==1.5,0.001); 
  Ysrr=vbm_vol_resize(single(Ysrr),'dereduce',size(Ym)); 
  Ysr(Ysr==1.5)=round(Ysrr(Ysr==1.5)); clear Ysrr; 

  Ysr(isinf(Ysr) & Ymfr>2)=0; Ysr=vbm_vol_downcut(Ysr,Ymfr,3,resTrr.vx_volr); 
  Ysr(isinf(Ysr) & Ymfr>0)=0; Ysr=vbm_vol_downcut(Ysr,Ymfr,3,resTrr.vx_volr); 
  Ysr(Ysr<=0 & Yl1Ar>0)=2-mod(Yl1Ar(Ysr<=0 & Yl1Ar>0),2);
  Ysr=round(vbm_vol_smooth3X(Ysr,2));
  Ys = vbm_vol_resize(Ysr,'dereduceV',resT3); 
  clear Ymr Ymfr Yl1Ar Ysr; 

  Yl1(Yl1>0)=Yl1(Yl1>0)+((Ys(Yl1>0)>1.5)); 
  Yl1(Yl1<0.5 | Yl1>=255)=0;
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
  [HDr,resTr] = vbm_vol_resize(Yl1>0,'reduceV',vx_vol,4,64);
  HDr = smooth3(vbm_vol_morph(HDr>0.1,'lc',3))>0.5;
  HD  = smooth3(Yl1>0 | vbm_vol_resize(HDr,'dereduceV',resTr))>0.5;
  Yl1(HD & Yl1<=0) = opt.LAB.HD(1); 
  Yl1(~HD | Yl1==255) = 0;
  Yl1 = uint8(Yl1);
  
  
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

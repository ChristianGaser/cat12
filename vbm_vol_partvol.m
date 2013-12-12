function [Ya1,Ycls,YBG,YMF] = vbm_vol_partvol(Ym,Ycls,Yb,Yy,vx_vol,opt)
% ______________________________________________________________________
% Use a segment map Ycls, the global intensity normalized T1 map Ym and 
% the atlas label map YA to create a individual label map Ya1. 
% The atlas contain main regions like cerebrum, brainstem, midbrain,
% cerebellum, ventricle, and regions with blood vessels. 
%
% This function try to solve the following problems:
%  1) Finding of the cerebrum, the cerebellum, the head, blood vessels, 
%     brain skin and other mayor structures based on atlas (YA) and 
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
%   [vol,Ya1,Yb,YMF] = vbm_vol_partvol(YA,Yp0,Ym,Yl0,opt)
%
%   INPUT:  YA = 3D-volume with brain regions (altas map)
%           Yp0  = 3D-volume with tissue propability map (CSF=1,GM=2;WM=3)
%           Ym   = intensity normalized T1 image (BG=0,CSF=1/3,GM=2/3,WM=1)
%           Yl0  = spm-classes 4-6 (intracranial=1,skull=2,background=3)
%           opt
%            .res    = resolution for mapping
%            .vx_vol = voxelsize
%            .LAB    = label of Ya1 map (see def.LAB definition below)
%            
%
%   OUTPUT: vol = structure with volumes
%           Ya1 = individual label map 
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



% ----------------------------------------------------------------------
% fast partitioning for B3C[, and LAS]
% ----------------------------------------------------------------------
% VBM atlas atlas map to find important structures for the LAS and the
% skull-stripping, which are the subcortical GM regions and the cerebellum.
% Maybe also WM hyperintensity have to be labeled here as a region without
% local correction - actual clear WMHs are handeled as GM.
% ----------------------------------------------------------------------
  if ~exist('opt','var'), opt=struct(); end
  def.res    = 1.6;
  def.vx_vol = [1 1 1];
  
  % definition of ROIs 
  def.LAB.CT =  1; % cortex
  def.LAB.MB = 13; % MidBrain
  def.LAB.BS = 13; % BrainStem
  def.LAB.CB =  3; % Cerebellum
  def.LAB.ON = 11; % Optical Nerv
  def.LAB.BG =  5; % BasalGanglia 
  def.LAB.TH =  9; % Hypothalamus 
  def.LAB.HC = 19; % Hippocampus 
  def.LAB.VT = 15; % Ventricle
  def.LAB.NV = 17; % no Ventricle
  def.LAB.BV =  7; % Blood Vessels
  def.LAB.NB =  0; % no brain 
  def.LAB.HD = 21; % head
  def.LAB.HI = 23; % WM hyperintensities
  
  opt = checkinopt(opt,def);


  %% map atlas to RAW space
  PA = cg_vbm_get_defaults('extopts.atlas'); PA = PA{1,1};
  VA = spm_vol(PA);
  YA = uint8(round(spm_sample_vol(VA,double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0)));
  YA = reshape(YA,size(Ym));
  clear Yy;
  
  Yp0  = (single(Ycls{1})*2/255 + single(Ycls{2})*3/255 + single(Ycls{3})/255) .* Yb; 

  
  % work on average resolution
  Ym0 = Ym; 
  [Ym,YA,Yp0,Yb,BB] = vbm_vol_resize({Ym,YA,Yp0,Yb},'reduceBrain',vx_vol,10,Yb);
  [Ym,Yp0,Yb,resTr] = vbm_vol_resize({Ym,Yp0,Yb},'reduceV',vx_vol,opt.res,64);
  [YA]              = vbm_vol_resize(YA ,'reduceV',vx_vol,opt.res,64,'nearest'); 
 
  vx_vol = resTr.vx_volr; 
  vxd    = 1/mean(vx_vol); 
  
  % prepare maps
  [~,~,YS] = vbdist(single(mod(YA,2)) + single(YA>0)); YS=(YS>1.5);     % side map
  YA(mod(YA,2)==0 & YA>0)=YA(mod(YA,2)==0 & YA>0)-1;                    % ROI map without side
  YA   = uint8(round(vbm_vol_median3c(single(YA),Yp0>0)));
  Yg   = vbm_vol_grad(Ym,vx_vol);
  Ydiv = vbm_vol_div(Ym,vx_vol);
  Ym   = Ym*3 .* (Yb);
  Yb   = Yb>0.5;
  
  %% Create individual mapping:
  noise = double(max(0.01,min(0.1,mean(Yg(Yp0>2.9))/3)));
  
  % Major structure mapping:
  % Major structure mapping with downcut to have a better alginment for 
  % the CB and CT. Simple setting of BG and TH as GM structures. 
  Ya1 = zeros(size(Ym),'single');
  Ya1(YA==opt.LAB.BG & Ym>1.9 & Ym<2.85 & Ydiv>-0.1)=opt.LAB.BG;                        % basal ganglia
  Ya1(YA==opt.LAB.TH & Ym>1.9 & Ym<2.85 & Ydiv>-0.1)=opt.LAB.TH;                        % thalamus
  Ybg = Ya1==0 & vbm_vol_morph(Ya1,'d',3) & ~vbm_vol_morph(YA==opt.LAB.VT,'d',4);       % VT correction area 
  Ya1(((Yp0>2.5 & Ym>2.5) & YA==opt.LAB.CT & Ya1==0) | Ybg)=opt.LAB.CT;                 % cerebrum
  Ya1((Yp0>2.0 & Ym>2.0) & YA==opt.LAB.CB)=opt.LAB.CB;                                  % cerebellum
  Ya1((Yp0>2.0 & Ym>2.0) & YA==opt.LAB.BS)=opt.LAB.BS;                                  % brainstem
  Ya1((Yp0>2.0 & Ym>2.0) & YA==opt.LAB.ON)=opt.LAB.ON;                                  % optical nerv
  Ya1((Yp0>2.0 & Ym>2.0) & YA==opt.LAB.MB)=opt.LAB.MB;                                  % midbrain
  Ya1(Ya1==0 & Yp0<1.5)=nan; Ya1 = vbm_vol_downcut(Ya1,Ym,2*noise); Ya1(isinf(Ya1))=0;  % region-growing
  Ya1=vbm_vol_median3c(Ya1,Yb);                                                         % smoothing
  Ya1((Yp0>1.75 & Ym>1.75 & Yp0<2.5 & Ym<2.5) & Ya1==opt.LAB.MB)=0;                       % midbrain correction
  
  
  %% Blood vessels
  Ywm = Yp0>2.8 & Ym>2.8 & Yp0<3.1 & Ym<4; % init WM 
  Ywm = single(vbm_vol_morph(Ywm,'lc',2));                 
  Ywm(smooth3(single(Ywm))<0.5)=0;           % remove small structures
  Ywm(~Ywm & (Yp0<1.2 | Ym<1.2 | Ym>4))=nan; 
  [Ywm1,YDr] = vbm_vol_downcut(Ywm,Ym,2*noise); Ywm(Ywm==-inf | YDr>20)=0; Ywm(Ywm1>0)=1; clear Ywm1
  Ywms=smooth3(single(Ywm)); Yms=smooth3(Ym);
  Ywm(Ywms<0.5)=0; Ywm(Ywms>0.5 & Yb & (Ym-Yms)<0.6)=1; Ywm(Ywms<0.5 & Yb & (Ym-Yms)>0.6)=0; clear Ywms
  % 
  Ybv=vbm_vol_morph( (Ym>3.3 & Yp0<2.2) | ((Ym-Yp0)>0.5 & Ym>2.2 & Yp0<2.2) | ...
    (Yms>2.5 & (Ym-Yms)>0.6) | (Ym>2.2 & Ywm==0) | ...
    (Ym>2.2 & Yp0<2.2 & Ya1==0 & YA==opt.LAB.CT),'c',1) & vbm_vol_morph(Ya1==opt.LAB.CT,'d',2) & ~Ywm;   
  Ybvs=smooth3(Ybv);
  Ybv(Ybvs>0.3 & Ym>2.5 & Yp0<2.5)=1; Ybv(Ybvs>0.3 & Ym>3.5 & Yp0<2.9)=1;
  Ybv(Ybvs<0.2 & Ym<4)=0;  clear Yvbs;
  Ya1(Ybv)=opt.LAB.BV; 
  clear Ybv Ywm
  
  
  %% Ventricle:
  % Ventricle estimation with a previous definition of non ventricle CSF
  % to have a second ROI in the region-growin. Using only the ventricle
  % ROI can lead to overgrowing. Using of non ventrilce ROI doesn't work
  % because dartle failed for large ventricle. 
  Ynv  = single(~Yb | YA==opt.LAB.CB | YA==opt.LAB.BS);
  Ynv  = vbm_vol_morph(Ynv,'d',8*vxd) | (YA==opt.LAB.NV);
  Ynv  = Ynv & Ym<2 & ~vbm_vol_morph(Yp0<1.8 & (YA==opt.LAB.VT) & Yg<0.3,'d',4*vxd);
  Ynv  = smooth3(Ynv)>0.5;
  %
  Yvt  = single(smooth3(Yp0<1.8 & (YA==opt.LAB.VT) & Yg<0.3 & Ydiv<0.01)>0.5); 
  Yvt(Yvt==0 & Ynv)=2; Yvt(Yvt==0 & Ym>2.5)=nan;
  Yvt = vbm_vol_downcut(Yvt,1-Ym,noise*2,vx_vol); Ynv = Yvt==2; 
  Ya1(Yvt==1 & Yp0<1.5)=opt.LAB.VT; 
 % clear Yvt
  
  %% WMH (White Matter Hyperintensities):
  % WMHs can be found as GM next to the ventricle (A) that do not belong 
  % to a subcortical structure (A) or there must be a big difference 
  % between the tissue SPM expect and the real intensity 'Yp0 - Ym' (C).
  % Furthermore no other Sulic (=near other CSF) should be labeld (D).
  Yvtd = vbm_vol_morph(Ya1==opt.LAB.VT,'d',2*vxd) & Ym>1.5 & Ym<2.5 & ...
        ~vbm_vol_morph(Ynv,'d',2*vxd) & Ya1~=opt.LAB.VT; % (A) - around the ventricle
  Yvtd = vbm_vol_downcut(single(Yvtd),1-Ym,-.1,vx_vol) & Ym>1.5 & Ym<2.5;
  Ybg  = vbm_vol_morph(Ya1>1 & Ya1<20 & Ya1~=opt.LAB.VT & Ya1~=opt.LAB.BV,'d',2); % (B) no deep GM
  Ywmh = single(vbm_vol_morph(( ~Ybg & ~Ynv & (Yp0 - Ym)>0.6) | ... % (C) 
    ( Yvtd & ~Ybg & ~Ynv & Ym>1.5 & Ym<2.5 & ( (Yp0 - Ym)>0.2) ),'c',2*vxd) & Ym>1.5 & Ym<2.5); 
%  clear Yvtd Ybg  
  % no a lit bit region growing for similiar values
  Ywmh = single(Ywmh); spm_smooth(Ywmh,Ywmh,vx_vol); Ywmh=round(Ywmh); Ywmh(Yp0<1.5 | Ym>2.75)=nan; 
  [Ywmh,YD]=vbm_vol_downcut(Ywmh,1-Ym,-.01,vx_vol); Ywmh(isinf(Ya1) | YD>100)=0; clear YD
  spm_smooth(Ywmh,Ywmh,vx_vol); Ywmh=round(Ywmh); 
  Ya1(Ywmh==1 & Ya1~=opt.LAB.VT)=opt.LAB.HI;
%  clear Ywmh;

  
  %% Closing of gaps between diffent structures:
  Yvtd2 = vbm_vol_morph(Ya1==opt.LAB.VT,'d',2*vxd) & Ya1~=opt.LAB.VT;
  % CT and VT
  Yt = vbm_vol_morph(Ya1==opt.LAB.VT,'d',2*vxd) & ...
       vbm_vol_morph(Ya1==opt.LAB.CT,'d',2*vxd) & Ya1==0 ;
  Ya1(Yt & Yp0<=1.5 & ~Ynv)=opt.LAB.VT; Ya1(Yt & Yp0>1.5)=opt.LAB.CT; % avoid other tissues in the ventricle
  % WMH and VT
  Yt = vbm_vol_morph(Ya1==opt.LAB.HI,'d',3*vxd) & Yvtd2 & ~Ynv & Ya1==0;
  Ya1(Yt &  Ym<=1.5)=opt.LAB.VT; Ya1(Yt & Ym>1.5 & Ym<2.5)=opt.LAB.HI; 
  % TH and VT
  Yt = vbm_vol_morph(Ya1==opt.LAB.TH,'d',3*vxd) & Yvtd2;
  Ya1(Yt & Ym<=1.5)=opt.LAB.VT; Ya1(Yt & Ym>1.5 & Ym<2.85)=opt.LAB.TH; 
  % BG and VT
  Yt = vbm_vol_morph(Ya1==opt.LAB.BG,'d',2*vxd) & Yvtd2;
  Ya1(Yt & Ym<=1.5)=opt.LAB.VT; Ya1(Yt & Ym>1.5 & Ym<2.85)=opt.LAB.BG;
  % BG and WMH
  Yt = vbm_vol_morph(Ya1==opt.LAB.BG,'d',4*vxd); 
  Yh = vbm_vol_morph(Ya1==opt.LAB.HI,'d',4*vxd);
  Ya1(Yt & Yh & YA~=opt.LAB.BG & Ym>1.5 & Ym<2.5 & ~Ynv)=opt.LAB.HI; 
  % no bloodvessels next to the ventricle, because for strong atrophy
  % brains the WM structures can be very thin and may still include strong
  % bias
  Ya1(Ya1==opt.LAB.BV & vbm_vol_morph(Ya1==opt.LAB.VT,'d',3*vxd))=0;
  clear Yt Yh Yvtd2 
  
 
  %% complete map
  [~,~,Ya1] = vbdist(Ya1,Yb);
  
  
  %% side aligment using laplace to correct for missalignments due to the normalization
  YBG  = Ya1==opt.LAB.BG | Ya1==opt.LAB.TH;
  YMF  = Ya1==opt.LAB.VT | Ya1==opt.LAB.HI | Ya1==opt.LAB.BG | Ya1==opt.LAB.TH; 
  YMF2 = vbm_vol_morph(YMF,'d',2*vxd) | Ya1==opt.LAB.CB | Ya1==opt.LAB.BS | Ya1==opt.LAB.MB;
  Ymf  = max(Ym,smooth3(single(YMF2*3))); 

  Yt = vbm_vol_smooth3X(YS==0,6)<0.9 & vbm_vol_smooth3X(YS==1,6)<0.9 & ~YMF2 & Yp0>0 & (Yp0<2.5 | Ya1==opt.LAB.BV);
  Ys = (single(YS)+1) .* ~Yt;
  Ys=vbm_vol_downcut(Ys,Ymf,0.1,vx_vol);
  clear YMF2 Yt YS;
  
  %% back to original size
  Ya1 = vbm_vol_resize(Ya1,'dereduceV',resTr,'nearest'); Ya1 = vbm_vol_median3c(Ya1,Ya1>0 & Ya1~=opt.LAB.BV);
  Ys  = vbm_vol_resize(Ys ,'dereduceV',resTr,'nearest'); Ys  = vbm_vol_median3c(Ys ,Ya1>0);
  YMF = vbm_vol_resize(YMF,'dereduceV',resTr);
  YBG = vbm_vol_resize(YBG,'dereduceV',resTr);
  
  Ya1 = vbm_vol_resize(Ya1,'dereduceBrain',BB); Ya1 = uint8(round(Ya1));
  Ys  = vbm_vol_resize(Ys ,'dereduceBrain',BB); [~,~,Ys] = vbdist(Ys); 
  YMF = vbm_vol_resize(YMF,'dereduceBrain',BB); 
  YBG = vbm_vol_resize(YBG,'dereduceBrain',BB); 
  Ym  = Ym0; clear Ym0;

  % final side alignment
  Ya1(Ya1>0)=Ya1(Ya1>0)+(Ys(Ya1>0)-1);
 
  
  % class correction
  YBGs = min(min(255-uint8(round(vbm_vol_smooth3X(Ya1==1 & Ycls{2}>250,0.8))).*Ycls{2},...
    uint8(round(255*vbm_vol_smooth3X(YBG,0.5) .* (Ym<2.9/3)))),255-Ycls{3});
  Ycls{2} = min(Ycls{2},255-YBGs);
  Ycls{1} = max(Ycls{1},YBGs);
  Ysum = zeros(size(Ym),'uint8'); for i=1:numel(Ycls), Ysum = Ysum + Ycls{i}; end;
  Ycls{2} = Ycls{2} + (255-Ysum).*uint8(Ym>2.75/3);
  Ycls{1} = Ycls{1} + (255-Ysum).*uint8(Ym<2.75/3 & YBGs==0);
  clear YBGs Ysum; 

end
%=======================================================================
function Yg = vbm_vol_grad(Ym,vx_vol)
% ----------------------------------------------------------------------
% gradient map for edge description
% ----------------------------------------------------------------------
  [gx,gy,gz] = vbm_vol_gradient3(Ym); 
  Yg = abs(gx./vx_vol(1))+abs(gy./vx_vol(2))+abs(gz./vx_vol(3)); 
  Yg = Yg ./ (Ym+eps);
end
%=======================================================================

%=======================================================================
function Ydiv = vbm_vol_div(Ym,vx_vol)
% ----------------------------------------------------------------------
% Diverence helps to identify all gyri that should not be in the GM, but 
% helps to improve the WM. Divergence estimation is very memory intensive 
% so it is better to limit the resolution.
% ----------------------------------------------------------------------
  [Ymr,resT2] = vbm_vol_resize(Ym,'reduceV',vx_vol,1.5,32);
  [gx,gy,gz]  = vbm_vol_gradient3(max(2/3,Ymr)); 
  Ydivr = smooth3(divergence(gy./vx_vol(1),gx./vx_vol(1),gz./vx_vol(3))); clear gx gy gz Ymr;
  Ydiv  = vbm_vol_resize(Ydivr,'dereduceV',resT2); 
end
%=======================================================================
function [Ya1,Ycls,YBG,YMF] = vbm_vol_partvol(Ym,Ycls,Yb,Yy,vx_vol)
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
%            .LAB    = label of Ya1 map (see LAB definition below)
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
  
  % definition of ROIs 
  LAB.CT =  1; % cortex
  LAB.MB = 13; % MidBrain
  LAB.BS = 13; % BrainStem
  LAB.CB =  3; % Cerebellum
  LAB.ON = 11; % Optical Nerv
  LAB.BG =  5; % BasalGanglia 
  LAB.TH =  9; % Hypothalamus 
  LAB.HC = 19; % Hippocampus 
  LAB.VT = 15; % Ventricle
  LAB.NV = 17; % no Ventricle
  LAB.BV =  7; % Blood Vessels
  LAB.NB =  0; % no brain 
  LAB.HD = 21; % head
  LAB.HI = 23; % WM hyperintensities
  
  WMHCstr = cg_vbm_get_defaults('extopts.WMHCstr'); 
  vx_res  = cg_vbm_get_defaults('extopts.vx_res')*1; %*2
  verb    = cg_vbm_get_defaults('extopts.verb')-1;
  debug   = cg_vbm_get_defaults('extopts.debug');

  %% map atlas to RAW space
  if verb, fprintf('\n'); end
  stime = vbm_io_cmd('  Atlas 2 Subjectspace','g5','',verb);
  PA = cg_vbm_get_defaults('extopts.atlas'); PA = PA{1,1};
  VA = spm_vol(PA);
  YA = vbm_vol_ctype(spm_sample_vol(VA,double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0));
  YA = reshape(YA,size(Ym));
 % clear Yy;
  
  Yp0  = (single(Ycls{1})*2/255 + single(Ycls{2})*3/255 + single(Ycls{3})/255) .* Yb; 

  
  % work on average resolution
  Ym0 = Ym; 
  [Ym,YA,Yp0,Yb,BB] = vbm_vol_resize({Ym,YA,Yp0,Yb},'reduceBrain',vx_vol,2,Yb);
  [Ym,Yp0,Yb,resTr] = vbm_vol_resize({Ym,Yp0,Yb},'reduceV',vx_vol,vx_res,64);
  [YA]              = vbm_vol_resize(YA ,'reduceV',vx_vol,vx_res,64,'nearest'); 
 
  vx_vol = resTr.vx_volr; 
  vxd    = 1/mean(vx_vol); 
  
  % prepare maps
  [tmp0,tmp1,YS] = vbdist(single(mod(YA,2)) + 2*single(YA>0)); YS=~mod(YS,2); clear tmp0 tmp1;  % side map
  YA(mod(YA,2)==0 & YA>0)=YA(mod(YA,2)==0 & YA>0)-1;                    % ROI map without side
  YA   = vbm_vol_ctype(vbm_vol_median3c(single(YA),Yp0>0));
  Yg   = vbm_vol_grad(Ym,vx_vol);
  Ydiv = vbm_vol_div(Ym,vx_vol);
  Ym   = Ym*3 .* (Yb);
  Yb   = Yb>0.5;

  
  
  %% Create individual mapping:
  stime = vbm_io_cmd('  Major Struktures','g5','',verb,stime);
  noise = double(max(0.02,min(0.1,mean(Yg(Yp0>2.9))/3)));
  
  % Major structure mapping:
  % Major structure mapping with downcut to have a better alginment for 
  % the CB and CT. Simple setting of BG and TH as GM structures. 
  Ya1 = zeros(size(Ym),'single');
  Ya1(YA==LAB.BG & Ym>1.9 & Ym<2.85 & Ydiv>-0.1)=LAB.BG;                % basal ganglia
  Ya1(YA==LAB.TH & Ym>1.9 & Ym<2.85 & Ydiv>-0.1)=LAB.TH;                % thalamus
  Ya1(YA==LAB.HC & Ym>1.9 & Ym<2.85 & Ydiv>-0.1)=LAB.HC;                % hippocampus
  Ybg = Ya1==0 & vbm_vol_morph(Ya1,'d',3) & ~vbm_vol_morph(YA==LAB.VT,'d',4); % VT correction area 
  Ya1(((Yp0>2.5 & Ym>2.5) & YA==LAB.CT & Ya1==0) | Ybg)=LAB.CT;         % cerebrum
  Ya1((Yp0>2.0 & Ym>2.0) & YA==LAB.CB)=LAB.CB;                          % cerebellum
  Ya1((Yp0>2.0 & Ym>2.0) & YA==LAB.BS)=LAB.BS;                          % brainstem
  Ya1((Yp0>2.0 & Ym>2.0) & YA==LAB.ON)=LAB.ON;                          % optical nerv
  Ya1((Yp0>2.0 & Ym>2.0) & YA==LAB.MB)=LAB.MB;                          % midbrain
  % region-growing
  Ya1(Ya1==0 & Yp0<1.5)=nan; 
  Ya1 = vbm_vol_downcut(Ya1,Ym,2*noise); Ya1(isinf(Ya1))=0; 
  Ya1 = vbm_vol_median3c(Ya1,Yb);                                       % smoothing
  Ya1((Yp0>1.75 & Ym>1.75 & Yp0<2.5 & Ym<2.5) & Ya1==LAB.MB)=0;         % midbrain correction
  
  
  
  %% Blood vessels
  % For this we may require the best resolution!
  % first a hard regions growing have to find the real WM-WM/GM region
  stime = vbm_io_cmd('  Blood Vessel detection','g5','',verb,stime);
  Ywm = Yp0>2.8 & Ym>2.8 & Yp0<3.1 & Ym<4;                              % init WM 
  Ywm = single(vbm_vol_morph(Ywm,'lc',2));                              % closing WM               
  Ywm(smooth3(single(Ywm))<0.5)=0;                                      % remove small dots
  Ywm(~Ywm & (Yp0<1.2 | Ym<1.2 | Ym>4))=nan;                            % set regions growing are
  [Ywm1,YDr] = vbm_vol_downcut(Ywm,Ym,2*noise);                         % region growing
  Ywm(Ywm==-inf | YDr>20)=0; Ywm(Ywm1>0)=1; clear Ywm1                  % set regions growing
  % smoothing
  Ywms = smooth3(single(Ywm)); Yms=smooth3(Ym);                         
  Ywm(Ywms<0.5)=0; Ywm(Ywms>0.5 & Yb & (Ym-Yms)<0.6)=1;                 
  Ywm(Ywms<0.5 & Yb & (Ym-Yms)>0.6)=0; clear Ywms                       
  % set blood vessels
  Ybv=vbm_vol_morph( (Ym>3.3 & Yp0<2.2) | ((Ym-Yp0)>0.5 & Ym>2.2 & Yp0<2.2) | ...
    (Yms>2.5 & (Ym-Yms)>0.6) | (Ym>2.2 & Ywm==0) | ...
    (Ym>2.2 & Yp0<2.2 & Ya1==0 & YA==LAB.CT),'c',1) & ...
    vbm_vol_morph(Ya1==LAB.CT,'d',2) & ~Ywm;  clear Ywm 
  % smoothing
  Ybvs = smooth3(Ybv);
  Ybv(Ybvs>0.3 & Ym>2.5 & Yp0<2.5)=1; Ybv(Ybvs>0.3 & Ym>3.5 & Yp0<2.9)=1;
  Ybv(Ybvs<0.2 & Ym<4)=0; clear Yvbs;
  Ya1(Ybv)=LAB.BV; clear Ybv 
  
  
  
  %% Ventricle:
  % Ventricle estimation with a previous definition of non ventricle CSF
  % to have a second ROI in the region-growin. Using only the ventricle
  % ROI can lead to overgrowing. Using of non ventrilce ROI doesn't work
  % because dartle failed for large ventricle. 
  stime = vbm_io_cmd('  Ventricle detection','g5','',verb,stime);
  Ynv = vbm_vol_morph(vbm_vol_morph(~Yb,'d',6*vxd) | (YA==LAB.CB | YA==LAB.BS),'d',2) | vbm_vol_morph(YA==LAB.NV,'e',1);
  Ynv = single(Ynv & Ym<2 & ~vbm_vol_morph(Yp0<2 & (YA==LAB.VT) & Yg<0.2,'d',4*vxd));
  Ynv = smooth3(round(Ynv))>0.5; 
  % between thamlamus
  Ynv = Ynv | (vbm_vol_morph(Ya1==LAB.TH,'c',10) & Yp0<2) | YA==LAB.CB | YA==LAB.BS;
  Ynv = smooth3(Ynv)>0.8;
  Yvt = single(smooth3(Yp0<1.5 & (YA==LAB.VT) & Yg<0.25 & ~Ynv)>0.7); 
  Yvt(Yvt==0 & Ynv)=2; Yvt(Yvt==0 & Ym>1.5)=nan; Yvt(Yvt==0)=1.5;
  Yvt2 = vbm_vol_laplace3R(Yvt,Yvt==1.5,0.005);
  Yvt = smooth3(round(Yvt2)==1 & Yp0<1.5)>0.5; 
  Ya1(Yvt)=LAB.VT; 

 
  %% WMH (White Matter Hyperintensities):
  % WMHs can be found as GM next to the ventricle (A) that do not belong 
  % to a subcortical structure (A) or there must be a big difference 
  % between the tissue SPM expect and the real intensity 'Yp0 - Ym' (C).
  % Furthermore no other Sulic (=near other CSF) should be labeld (D).
  % ####################################################################
  % There can also be deep GM Hyperintensities! 
  % ####################################################################
  stime = vbm_io_cmd('  WMH detection','g5','',verb,stime);
  Ywmh = single(smooth3(vbm_vol_morph(Yvt2<1.45,'d',1) & Ym<2.25 &...
    ~vbm_vol_morph(YA==LAB.HC,'d',4*vxd))>0.75-WMHCstr/2); % ventricle
  Ywmh(smooth3((Yp0 - Ym)>0.75-WMHCstr/2 & Ym<2.75 & Ym>1.5)>0.8)=1; % WMH
  Ywmh(vbm_vol_morph(Yp0>2.5,'c',0) & ~(Yp0>2.5))=1; % WMH wholes
  Ywmh(Yvt2>1.75 & Yvt2<3 | (Ywmh==0 & Ym<1.5))=2;
  Ywmh((Ywmh==0 & Ym>2.25) | Ya1==LAB.BG | Ya1==LAB.TH)=-inf;
  Ywmh = vbm_vol_downcut(Ywmh,(3-Ym)/3,max(-0.02,noise*WMHCstr - (0.5-WMHCstr)),vx_vol); 
  Ywmh(Ywmh==2 & smooth3(Ywmh==2)<0.25+WMHCstr/2)=0;
  Ywmh(Ywmh==1 & smooth3(Ywmh==1)<0.25+WMHCstr/2)=0;
  %
  Ywmh(isinf(Ywmh))=0;
  Ywmh((Ywmh==0 & Ym>2.75) | Ya1==LAB.BG | Ya1==LAB.TH)=nan; Ywmh(Ywmh==0)=1.5;
  Ywmh = vbm_vol_laplace3R(Ywmh,Ywmh==1.5,0.005);
  Ywmh(vbm_vol_morph(YS==1,'d',3*vxd) & vbm_vol_morph(YS==0,'d',3*vxd))=2; % not for the CC
  Ynwmh = ~smooth3(vbm_vol_morph(Ya1==LAB.VT | Ya1==LAB.TH | Ya1==LAB.BG,'c',4))>0.5;
  Ywmh = smooth3(Ywmh<1.25+WMHCstr/4 & Ya1~=LAB.VT & Ynwmh)>0.75-WMHCstr/2;
  Ywmh = smooth3(vbm_vol_morph(vbm_vol_morph(Ywmh,'o',1),'c',1) & Ym<2.75)>(0.75-WMHCstr/2);
  %%
  Ya1(Ywmh)=LAB.HI;
  %{
   Yvt2(Yvt2>1.45 & Yvt2<1.55)=inf; Yvt2=round(Yvt2);
  Yvt2(Yvt2>3 & Ym<2.5 & Ym>1.5 & YA~=LAB.CB & YA~=LAB.BS & Ya1~=LAB.BG & Ya1~=LAB.TH)=1.5;
  Yvt2(Yvt2==1.5 & smooth3((Yp0 - Ym)>0.6 & Ym<2.75 & Ym>1.75)>0.7)=1; 
  Yvt2(Yvt2>3)=1.5;
  Yvt2 = vbm_vol_laplace3R(Yvt2,Yvt2==1.5 & Ym<2.75 & YA~=LAB.BS & Ya1~=LAB.BG & Ya1~=LAB.TH,0.005); %Yvt2 = round(Yvt2);
  %}
  if ~debug, clear Ywmh Yvt Ynwmh; end

  
  
  %% Closing of gaps between diffent structures:
  stime = vbm_io_cmd('  Closing for deep structures','g5','',verb,stime);
  Yvtd2 = vbm_vol_morph(Ya1==LAB.VT,'d',2*vxd) & Ya1~=LAB.VT;
  % CT and VT
  Yt = vbm_vol_morph(Ya1==LAB.VT,'d',2*vxd) & ...
       vbm_vol_morph(Ya1==LAB.CT,'d',2*vxd) & Ya1==0 ;
  Ya1(Yt & Yp0<=1.5 & ~Ynv)=LAB.VT; Ya1(Yt & Yp0>1.5)=LAB.CT; 
  % WMH and VT
  Yt = vbm_vol_morph(Ya1==LAB.HI,'d',1*vxd) & Yvtd2 & ~Ynv & Ya1==0;
  Ya1(Yt &  Ym<=1.25)=LAB.VT; Ya1(Yt & Ym>1.25 & Ym<2.5)=LAB.HI; 
  % TH and VT
  Yt = vbm_vol_morph(Ya1==LAB.TH,'d',1*vxd) & Yvtd2;
  Ya1(Yt & Ym<=1.5)=LAB.VT; Ya1(Yt & Ym>1.5 & Ym<2.85)=LAB.TH; 
  % BG and VT
  Yt = vbm_vol_morph(Ya1==LAB.BG,'d',1*vxd) & Yvtd2;
  Ya1(Yt & Ym<=1.5)=LAB.VT; Ya1(Yt & Ym>1.5 & Ym<2.85)=LAB.BG;
  % no bloodvessels next to the ventricle, because for strong atrophy
  % brains the WM structures can be very thin and may still include 
  % strong bias
  Ya1(Ya1==LAB.BV & vbm_vol_morph(Ya1==LAB.VT,'d',3*vxd))=0;
  clear Yt Yh Yvtd2 Yw
 
  
  
  %% complete map
  [tmp0,tmp1,Ya1] = vbdist(Ya1,Yb); clear tmp0 tmp1;
  
  
  
  %% side aligment using laplace to correct for missalignments due to the normalization
  stime = vbm_io_cmd('  Side Alignment','g5','',verb,stime);
  YBG  = Ya1==LAB.BG | Ya1==LAB.TH;
  YMF  = Ya1==LAB.VT | Ya1==LAB.BG | Ya1==LAB.TH | Ya1==LAB.HI; 
  YMF2 = vbm_vol_morph(YMF,'d',2*vxd) | Ya1==LAB.CB | Ya1==LAB.BS | Ya1==LAB.MB;
  Ymf  = max(Ym,smooth3(single(YMF2*3))); 
  Yt = vbm_vol_smooth3X(YS==0,6)<0.9 & vbm_vol_smooth3X(YS==1,6)<0.9 & ~YMF2 & Yp0>0 & Ym<3.1 & (Yp0<2.5 | Ya1==LAB.BV);
  Ys = (2-single(YS)) .* single(smooth3(Yt)<0.4);
  Ys(Ys==0 & (Ym<1 | Ym>3.1))=nan; Ys = vbm_vol_downcut(Ys,Ymf,0.1,vx_vol); 
  [tmp0,tmp1,Ys] = vbdist(Ys,Ys==0);
  clear YMF2 Yt YS tmp0 tmp1;
  
  % YMF for FreeSurfer fsaverage
  Ysm  = vbm_vol_morph(Ys==2,'d',1.75*vxd) & vbm_vol_morph(Ys==1,'d',1.75*vxd);
  YMF  = vbm_vol_morph(Ya1==LAB.VT | Ya1==LAB.BG | Ya1==LAB.HI | (Ya1==LAB.TH & smooth3(Yp0)>2),'c',3) & ~Ysm; 
  YMF  = smooth3(YMF)>0.5;
  clear Ysm; 
  
  
  %% back to original size
  stime = vbm_io_cmd('  Final corrections','g5','',verb,stime);
  Ya1 = vbm_vol_resize(Ya1,'dereduceV',resTr,'nearest'); Ya1 = vbm_vol_median3c(Ya1,Ya1>1 & Ya1~=LAB.BV);
  Ys  = vbm_vol_resize(Ys ,'dereduceV',resTr,'nearest'); Ys  = 1 + single(smooth3(Ys)>1.5);
  YMF = vbm_vol_resize(YMF,'dereduceV',resTr);
  YBG = vbm_vol_resize(YBG,'dereduceV',resTr);
  
  Ya1 = vbm_vol_resize(Ya1,'dereduceBrain',BB); Ya1 = vbm_vol_ctype(Ya1);
  Ys  = vbm_vol_resize(Ys ,'dereduceBrain',BB); [tmp0,tmp1,Ys] = vbdist(Ys,Ya1>0); clear tmp0 tmp1;
  YMF = vbm_vol_resize(YMF,'dereduceBrain',BB); 
  YBG = vbm_vol_resize(YBG,'dereduceBrain',BB); 
  Ym  = Ym0; clear Ym0;

  % final side alignment
  Ya1(Ya1>0)=Ya1(Ya1>0)+(Ys(Ya1>0)-1);
 
  
  % class correction
  % YBG is smoothed a little bit and (B) reset all values that are related
  % with GM/WM intensity (Ym<2.9/3) (A)
  Yclssum = single(Ycls{1})+single(Ycls{2})+single(Ycls{3});
  YBGs    = min( max(0,min(255, 255 - vbm_vol_smooth3X(Ya1==1 & Ycls{2}>round(2.9/3),0.8) .* single(Ycls{2}) )), ... (A)
                 max(0,min(255, 255 * vbm_vol_smooth3X(YBG .* (Ym<=2.9/3 & Ym>2/3) ,0.5) )) ); % (B)
  Ycls{1} = vbm_vol_ctype(single(Ycls{1}) + YBGs .* (single(Ycls{2})./max(eps,Yclssum)));
  Ycls{2} = vbm_vol_ctype(single(Ycls{2}) - YBGs .* (single(Ycls{2})./max(eps,Yclssum)));
  clear YBGs Yclssum; 
 
  if verb, vbm_io_cmd(' ','','',verb,stime); end

end
%=======================================================================
function Yg = vbm_vol_grad(Ym,vx_vol)
% ----------------------------------------------------------------------
% gradient map for edge description
% ----------------------------------------------------------------------
  [gx,gy,gz] = vbm_vol_gradient3(Ym); 
  Yg = abs(gx./vx_vol(1))+abs(gy./vx_vol(2))+abs(gz./vx_vol(3)); 
  Yg = Yg ./ max(eps,Ym);
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
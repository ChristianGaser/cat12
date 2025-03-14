function [Ya1,Ycls,YMF,Ycortex] = cat_vol_partvol1639(Ym,Ycls,Yb,Yy,vx_vol,extopts,Vtpm,noise,job,Ylesionmsk,Ydt,Ydti)
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
%     use the same intensity scaling as the segment map Yp0, but have 
%     more information about partial volume regions.
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
%   [vol,Ya1,Yb,YMF] = cat_vol_partvol(YA,Yp0,Ym,Yl0,opt)
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
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

% ______________________________________________________________________
%
% Development comments:
%   ToDo:
%   - WMHs werden bei geringer Aufl?sung ?berschaetzt
%   - mehr Kommentare bei WMHC und SLC
%
%   Was ist neu im Vergleich zu anderen?
%   - Zuweisung durch Dartel mit hoher Genauigkeit moeglich 
%   - Erweiterung von SPM/VBM durch MainROIs (Seiten, Lappen, ...)
%   - Verbesserung der SPM/VBM durch bessere Enfernung von unerwuenschtem
%     Gewebe (ON, Blutgefaesse ...)
%   - Blutgefaesse koennnen als erweitere Masken fuer fMRI genutzt werden 
%     um Seiteneffekte besser ausblenden zu koennen.
%  [- Beliebige Atlanten koennen genutzt werden.]
%
%  Todo:
%   - Besserer Atlas
%   - BV vs. HD - glaetten in dilated HD region
%   - Fuellen von CSF Luecken bei LAB~=BV und Ym<1.2 und LAB==NV?
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
  
%   LAB.CT =  1; % cortex
%   LAB.MB = 13; % MidBrain
%   LAB.BS = 13; % BrainStem
%   LAB.CB =  3; % Cerebellum
%   LAB.ON = 11; % Optical Nerv
%   LAB.BG =  5; % BasalGanglia 
%   LAB.TH =  9; % Hypothalamus 
%   LAB.HC = 19; % Hippocampus 
%   LAB.VT = 15; % Ventricle
%   LAB.NV = 17; % no Ventricle
%   LAB.BV =  7; % Blood Vessels
%   LAB.NB =  0; % no brain 
%   LAB.HD = 21; % head
%   LAB.HI = 23; % WM hyperintensities
%   LAB.PH = 25; % Gyrus parahippocampalis

  % if there is a breakpoint in this file set debug=1 and do not clear temporary variables 
  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

  def.uhrlim = 0.7; 
  extopts = cat_io_checkinopt(extopts,def); 

  LAB     = extopts.LAB;
  BVCstr  = mod(extopts.BVCstr,1) + (extopts.BVCstr==1 || extopts.BVCstr==2); 
  verb    = extopts.verb-1;
  PA      = extopts.cat12atlas;
  vx_res  = max( extopts.uhrlim , max( [ max(vx_vol) min(vx_vol) ] )); 
  noise   = double(noise);

  Yp0toC  = @(Yp0,c) 1-min(1,abs(Yp0-c));
  
  %% map atlas to RAW space
  if verb, fprintf('\n'); end
  stime = cat_io_cmd('  Atlas -> subject space','g5','',verb); 
  % CAT atlas
  YA = cat_vol_ctype( cat_vol_sample(Vtpm(1),PA{1},Yy,0) );
  
  
  % template map
  Yp0A = single( cat_vol_sample(Vtpm(1),Vtpm(1),Yy,1) )*2 + ...
         single( cat_vol_sample(Vtpm(1),Vtpm(2),Yy,1) )*3 + ...
         single( cat_vol_sample(Vtpm(1),Vtpm(3),Yy,1) )*1;
  
  
  % WMH atlas
  watlas = 3; 
  switch watlas
    case 1, PwmhA = strrep(PA{1},'cat.nii','cat_wmh_soft.nii');
    case 2, PwmhA = strrep(PA{1},'cat.nii','cat_wmh.nii');
    case 3
      if isfield(job.extopts,'SLtpm')
        PwmhA = job.extopts.WMHtpm{1};
      else
        PwmhA = strrep(PA{1},'cat.nii','cat_wmh_miccai2017.nii');
      end
  end
  if ~isempty(PwmhA) && exist(PwmhA,'file') && ~strcmp(PwmhA,PA{1}) 
    YwmhA = cat_vol_ctype( cat_vol_sample(Vtpm(1),PwmhA,Yy,0) );
  else
    YwmhA = max(0,min(1,Yp0A-2)); 
  end
  switch watlas
    case 2, YwmhA = min(1,max(0,YwmhA - 0.1) * 0.8);
    case 3, YwmhA = min(1,max(0,cat_vol_smooth3X(YwmhA,1) - 0.01) * 10); 
  end

  
  % Stroke lesion atlas
  if isfield(job.extopts,'SLtpm')
    PslA = job.extopts.SLtpm{1}; 
  else
    PslA = strrep(PA{1},'cat.nii','cat_strokelesions_ATLAS303.nii');
  end
  if ~isempty(PslA) && exist(PslA,'file') && ~strcmp(PslA,PA{1}) 
    YslA = cat_vol_ctype( cat_vol_sample(Vtpm(1),PslA,Yy,0) );
    YslA = YslA./max(YslA(:)); 
  else
    YslA = max(0,min(1,Yp0A-2)); 
  end
  

  % Blood vessel probability map (added 202305):
  % Uses a blood vessel map created using MRA scans of the IXI and ICBM  
  % databases in combination with CSF and WM probability maps as far as 
  % larger blood vessels are typically located along the brainstem, corpus 
  % callosum and within the insula.
  if isfield(job.extopts,'BVtpm')
    Pbv = job.extopts.BVtpm{1}; 
  else
    Pbv = strrep(PA{1},'cat.nii','cat_bloodvessels.nii');
  end
  YwmA  = single(cat_vol_sample(Vtpm(1),Vtpm(2),Yy,1)); 
  YcsfA = single(cat_vol_sample(Vtpm(1),Vtpm(3),Yy,1)); 
  if ~isempty(Pbv) && exist(Pbv,'file') 
    Vbv  = spm_vol(Pbv);
    YbvA = single(cat_vol_sample(Vbv,Vbv,Yy,1));
    YbvA = 1 - YwmA + max(YcsfA * 0.1 , YbvA); 
  else
    YbvA = 1 - YwmA + YcsfA * 0.1;
  end
  clear YwmA YcsfA; 
  clear Yy; 


  % use addition FLAIR images
  if exist('job','var') && isfield(job,'data_wmh') && ~isempty(job.data_wmh) && isfield(job,'subj') && numel(job.data_wmh)>=job.subj
    %%
    [pp,ff,ee] = spm_fileparts(job.data_wmh{job.subj}); 
    Pflair = fullfile(pp,[ff ee]); 
    if ~isempty(Pflair) && exist(Pflair,'file')
      stime = cat_io_cmd('  FLAIR corregistration','g5','',verb,stime); 
    
      % coreg
      Vflair = spm_vol(job.data_wmh{job.subj}); 
      Vm     = spm_vol(job.data{job.subj}); 
      evalc('R = spm_coreg(Vm,Vflair,struct(''graphics'',0));'); 
      R      = spm_matrix(R);  %#ok<NODEF>

      % load
      Vflair.dat   = cat_vol_sanlm(struct('verb',0),Vflair,1,spm_read_vols(Vflair));
      Vflair.pinfo = repmat([1;0],1,size(Vflair.dat,3));
      Vflair.dt(1) = 16;
      Yflair = zeros(Vm.dim,'single');
      for i=1:Vm.dim(3)
        Yflair(:,:,i) = single( spm_slice_vol(Vflair, R \ Vflair.mat \ Vm.mat  * spm_matrix([0 0 i]) ,Vm.dim(1:2),[1,NaN])); 
      end    
    end
  end
  
  
  %% resize data
  if ~debug; clear Yy; end
  
  
  Yp0  = (single(Ycls{1})*2/255 + single(Ycls{2})*3/255 + single(Ycls{3})/255) .* Yb; 
  if isfield(job.extopts,'inv_weighting') && job.extopts.inv_weighting
    Ym = Yp0/3;
  end

  % work on average resolution
  Ym0 = Ym; 
  % remove background
  [Ym,BB]  = cat_vol_resize(Ym      ,'reduceBrain',vx_vol,2,Yb);
  YA       = cat_vol_resize(YA      ,'reduceBrain',vx_vol,2,Yb);
  Yp0      = cat_vol_resize(Yp0     ,'reduceBrain',vx_vol,2,Yb);
  if exist('Ydt','var')
    Ydt    = cat_vol_resize(Ydt     ,'reduceBrain',vx_vol,2,Yb);
  end
  if exist('Ydti','var')
    Ydti   = cat_vol_resize(Ydti    ,'reduceBrain',vx_vol,2,Yb);
  end
  Yp0A     = cat_vol_resize(Yp0A    ,'reduceBrain',vx_vol,2,Yb);
  YslA     = cat_vol_resize(YslA    ,'reduceBrain',vx_vol,2,Yb);
  YbvA     = cat_vol_resize(YbvA    ,'reduceBrain',vx_vol,2,Yb);
  YwmhA    = cat_vol_resize(YwmhA   ,'reduceBrain',vx_vol,2,Yb);
  if exist('Yflair','var')
    Yflair = cat_vol_resize(Yflair  ,'reduceBrain',vx_vol,2,Yb);
  end
  if exist('Ylesionmsk','var')
    Ylesionmsk = cat_vol_resize(Ylesionmsk  ,'reduceBrain',vx_vol,2,Yb);
  end
  Yb       = cat_vol_resize(Yb      ,'reduceBrain',vx_vol,2,Yb);
  % use lower resolution 
  [Ym,resTr] = cat_vol_resize(Ym    ,'reduceV',vx_vol,vx_res,64);
  YA         = cat_vol_resize(YA    ,'reduceV',vx_vol,vx_res,64,'nearest'); 
  Yp0        = cat_vol_resize(Yp0   ,'reduceV',vx_vol,vx_res,64);
  if exist('Ydt','var')
    Ydt      = cat_vol_resize(Ydt   ,'reduceV',vx_vol,vx_res,64);
  end
  if exist('Ydti','var')
    Ydti     = cat_vol_resize(Ydti  ,'reduceV',vx_vol,vx_res,64);
  end
  Yp0A       = cat_vol_resize(Yp0A  ,'reduceV',vx_vol,vx_res,64);
  YslA       = single(cat_vol_resize(YslA  ,'reduceV',vx_vol,vx_res,64));
  YbvA       = cat_vol_resize(YbvA  ,'reduceV',vx_vol,vx_res,64);
  YwmhA      = cat_vol_resize(YwmhA ,'reduceV',vx_vol,vx_res,64); 
  if exist('Yflair','var')
    Yflair   = cat_vol_resize(Yflair,'reduceV',vx_vol,vx_res,64); 
  end    
  if exist('Ylesionmsk','var')
    Ylesionmsk   = cat_vol_resize(Ylesionmsk,'reduceV',vx_vol,vx_res,64); 
  end    
  Yb         = cat_vol_resize(Yb    ,'reduceV',vx_vol,vx_res,64);
  vx_vol     = resTr.vx_volr; 
  
  
  % noise reduction
  spm_smooth(Ym,Ym,0.6./vx_vol);
  
  
  % prepare maps
  YA   = cat_vol_ctype(cat_vol_median3c(single(YA),Yp0>0));                % noise filter of atlas map
  [~,~,YS] = cat_vbdist(single(mod(YA,2)) + single(YA>0)); YS=~mod(YS,2);  % side map
  YA(mod(YA,2)==0 & YA>0)=YA(mod(YA,2)==0 & YA>0)-1;                       % ROI map without side
  YS   = cat_vol_smooth3X(YS,4) > 0.5;                                     % RD202501: side smoothing 
  Yg   = cat_vol_grad(Ym,vx_vol);                                          % gadient map (edge map)
  Ydiv = cat_vol_div(Ym,vx_vol);                                           % divergence map (edge map)
  Ym   = Ym*3 .* (Yb);
  Yb   = Yb>0.5;

  
  
  %% Create individual mapping:
  stime = cat_io_cmd('  Major structures','g5','',verb,stime); 
  
  % Mapping of major structure:
  % Major structure mapping with downcut to have a better alginment for 
  % the CB and CT. Simple setting of BG and TH as GM structures. 
  Ya1  = zeros(size(Ym),'single');
  Ybg  = zeros(size(Ym),'single');
  Ybgd = cat_vbdist(single(YA==LAB.BG),Yb,vx_vol); 
  Yosd = cat_vbdist(single(YA==LAB.TH | YA==LAB.VT | YA==LAB.HC | ...
          YA==LAB.BS | (YA==LAB.CT & Ym>2.9)),Yb,vx_vol); 
  Ybg(smooth3(Yosd>3 & Ybgd<5  & Ym>1.9 & Ym<2.85 & Yg<4*noise & ...
          ((Ybgd<1 & Ydiv>-0.01) | (Ydiv>-0.01+Ybgd/100)))>0.7) = 1;
  Ybg(smooth3((Ybg==0 & Yp0>2.8 & Ym>2.8 & YA==LAB.CT) | Ym>2.9 | ...
          YA==LAB.TH | YA==LAB.HC | Yosd<2 | (Ybg==0 & Yp0<1.25) | ...
          (Ybg==0 & Ybgd>8) | (Ybg==0 & Ydiv<-0.01+Ybgd/200))>0.3) = 2;
  Ybg(Ybg==0 & Ybgd>0 & Ybgd<10) = 1.5; 
  Ybg = cat_vol_laplace3R(Ybg,Ybg==1.5,0.005)<1.5 & Ym<2.9 & Ym>1.8 & Ydiv>-0.02;
  Ya1(Ybg & Ym>1.5 & Yp0>1.5) = LAB.BG;                                      % basal ganglia
  Ya1(YA==LAB.TH & Ym>1.9 & Ym<2.85 & Ydiv>-0.1) = LAB.TH;                   % thalamus
  % RD202501: hippocampus definition was not optimal here as we need a clear parahippocampus for surface reconstruction 
  Ya1( cat_vol_morph(YA==LAB.HC,'dd',3) & Ym>1.5 & Ym<2.5 & ~(cat_vol_morph(YA==LAB.PH,'dd',2) | ...
    YA==LAB.TH | cat_vol_morph(YA==LAB.NV,'dd',5,vx_vol))) = LAB.HC;         % hippocampus 
  Ya1(Yp0>1.9 & Ym>1.7 & YA==LAB.HC) = LAB.HC;                               % hippocampus 
  NBVC = BVCstr > 0 && BVCstr <= 1; % RD202306: new BV correction by default
  %%
  if NBVC
    % define some high intensity BVs and the neocortex
    % RD202501: extened the BV detection 
    Ybv = ( (smooth3(Yp0)-Yp0 )>.7 | (Ym-Ydiv)>3.4 | max(0,Ym - Yp0)>Ydiv+.8 & Ym>2 & Ydiv<0) & YA==LAB.CT & smooth3(Ya1>0)<.5; 
    Ybv = single( cat_vol_morph( Ybv , 'l', [inf 1])>0 ); % remove single voxels
    Ybv( Ybv==0 & Yp0<=2 & YA==LAB.CT) = nan; Ybv(YA~=LAB.CT ) = nan;  
    Ybv( cat_vol_morph( cat_vol_morph( Yp0>2.1 & Ym<3.2 & (YA==LAB.CT | Ybgd<7),'o'),'lc')  ) = 2; % add neocortex region
    Ybv( cat_vol_morph( YA==LAB.HC | YA==LAB.VT | YA==LAB.PH | YA==LAB.TH,'dd', 4) ) = 2;          % add "neocortex" region 
    Ybv = cat_vol_downcut(Ybv,Ym - Ydiv,noise); 
    Ybv(Yp0>2.1 & Ym>2.3 & Ybv<2 & YA==LAB.CT) = 1; 
    Ybv = single( cat_vol_morph( Ybv==1 , 'l', [inf 2])>0 ); % remove single voxels
    Ya1( cat_vol_morph( Ybv==1,'dc') ) = LAB.BV; clear Ybv
    Ya1( cat_vol_morph( ((Yp0>2.5 & Ym>2.5 & Ym<3 -(YbvA.^4*BVCstr - 1) & (YA==LAB.CT | YA==LAB.BG)) | ...
         (Yp0>1.5 & Ym<3.01 & Yp0A>2.5 & Ybgd>1 & Ybgd<8 & (Ybgd>4 | ...
         Ydiv<-0.02+Ybgd/200))) & Ya1==0 , 'l',[.1 10 ])) = LAB.CT;        % cerebrum
    Yhbv = max(0,8*-Ydiv).^4 .* (2*Yg).^4 .* Ym .* 2.*smooth3(Yp0<2.9) .* (YbvA - 1).^2 .* (YA~=LAB.CB); %   .* cat_vol_morph(Ya1==0,'e');
    Ya1(Ya1==0 & Yhbv>.2*(1-BVCstr)) = LAB.BV;
  else
    % older version
    Ya1(((Yp0>2.5 & Ym>2.5 & (YA==LAB.CT | YA==LAB.BG )) | ...
       (Yp0>1.5 & Ym<3.5 & Ybgd>1 & Ybgd<8 & (Ybgd>4 | ...
       Ydiv<-0.02+Ybgd/200))) & Ya1==0)=LAB.CT;                            % cerebrum
  end
  % use region-growing in case of the cerebellum to compensate for fine structure
  Ycb  = single(YA==LAB.CB) + 2*single(cat_vol_morph( YA==LAB.CT,'de',6)); 
  Ycb(Yp0<1 & Ym<2) = nan; Ycb(cat_vol_morph( YA~=LAB.CB ,'de',12)) = nan; 
  [Ycb,Yd] = cat_vol_downcut(Ycb,Ym,noise); 
  Ycb = single(YA==LAB.CB) | cat_vol_smooth3X(Ycb==1 & Yd<50,2)>.7; 
  % cortext growing without BVs
  Yct  = single( (Ya1==LAB.CT | (Ybgd>2 & Ybgd<7)) & (Ym - Ydiv)<3.1 ) + 2*single(cat_vol_morph( YA==LAB.CB,'de',6)); 
  Yct(Yp0<1 & Ym<1.7) = nan; Yct(cat_vol_morph( YA~=LAB.CT,'de',12)) = nan; 
  [Yct,Yd] = cat_vol_downcut(Yct,Ym - Ydiv,noise/4); 
  Yct = (Yct==1 & smooth3(Yd)<20); 
  
  
  % set other regions
  Ya1(Yp0>1.9 & Ym>1.7 & Ym<3.1 & Ycb)=LAB.CB;                             % cerebellum
  Ya1(Yp0>1.9 & Ym>1.7 & YA==LAB.BS)=LAB.BS;                               % brainstem
  Ya1(Yp0>1.9 & Ym>1.7 & YA==LAB.ON)=LAB.ON;                               % optical nerv
  Ya1(Yp0>1.9 & Ym>1.7 & Yp0<3 & YA==LAB.TH)=LAB.TH;                       % thalamus
  Ya1(Yp0>1.9 & Ym>1.7 & cat_vol_smooth3X(YA==LAB.MB,2)>.1)=LAB.MB;        % midbrain
  Ya1(Yp0>1.9 & Ym>1.7 & YA==LAB.PH & Ya1==0)=LAB.CT;                      % added PH as cortex
  Ya1(Yp0>1.9 & Ym>1.7 & Ym<3.0 & Yct & Ya1==0)=LAB.CT;                    % added CT
  Ya1(Yp0>1.9 & Ym>1.7 & cat_vol_morph(YA==LAB.VT | YA==LAB.HC | YA==LAB.PH ,'dd',2) & Ya1==0)=LAB.CT;                    % added CT
  % denoising for regions but not BV
  Ya1 = cat_vol_median3c(Ya1,Ya1>0 & ~YA==LAB.BV);    
      %Ya1(Ybv>1 & ~cat_vol_morph(Ya1>1,'l',[inf 2]))=0;
  %%
  if NBVC
    % create a cerebellar mask to avoid corrections there
    Ycb3 = cat_vol_morph(YA==LAB.CB,'d',1) | (cat_vol_morph(YA==LAB.CB,'d',3) & YbvA<0.5 & Ym<2.5); 
    Ysc5 = cat_vol_morph(YA==LAB.HC | YA==LAB.BG | YA==LAB.TH,'dd',3); 
    Yhc3 = smooth3(cat_vol_morph(YA==LAB.HC | YA==LAB.PH | YA==LAB.VT | Ysc5,'d'))>.5;
    Ya1(Ysc5 & Ym>2.6 & Ym<2.95 & Ya1==0) = LAB.CT;  
    
    % extend the neocortical area
    Ya1(Ya1<1 & Ym<1.8) = nan; 
    [Ya1,Yd] = cat_vol_downcut(Ya1,Ym,noise); Ya1(isinf(Ya1) | Yd>15 ) = 0; 

    % get blood vessels as intensity-based eikonal-distance map with
    % different grow rates, where the first one is more limited by the
    % local intensities and the second one allows to remove the general 
    % distance aspect
    Ya0 = single(Ya1>0 & Ya1~=LAB.BV); Ya0(Ym<1.7) = nan; 
    [~,Yd] = cat_vol_downcut(Ya0,Ym,noise  ); clear Ya0 
    Ya0 = single(Ya1>0 & Ya1~=LAB.BV); Ya0(Ym<1.7) = nan; 
    [~,Yd2] = cat_vol_downcut(Ya0,Ym,noise * 16); clear Ya0
    Ya1( ((Yd - Yd2)>10./YbvA) & Ym>2.4 & YA==LAB.CT & ~Ycb3 & ~Yhc3 & YbvA>.7) = LAB.BV; % add highly distant voxels
    Ya1(Ya1==0 & Ym>2.4 & Yd>70 & Yd2>50 & ~cat_vol_morph(YA>1 & YA~=LAB.HD,'d') & ~Ycb3 & ~Yhc3) = LAB.BV; % mask further BV
    Ya1(Ya1==0 & Ym>2.2 & Yd<20 & ~cat_vol_morph(YA>1 & YA~=LAB.HD,'d') & ~Ycb3) = LAB.CT; % add also save brain voxels
    Ya1(Ya1==LAB.BV & cat_vol_morph(Ya1==LAB.BV,'l',[inf 16])==0) = 0;
    clear Ycb3 Yhc3
    
    %% modified old lines
    Ya1((Ya1==0 & Yp0<1.5 & Ym<1.5 & Yp0>1.3 & Ym>1.3) & YA==LAB.BV) = LAB.BV; % low-int BV  (updated due to cerebellar errors RD20190929)
    Ya1((cat_vol_morph(Ya1==0 & YA~=LAB.CB,'e') & (Ym>2.5 | Ym<1.7) & YA==LAB.CT & ...
      Ym>2 - (YbvA-1).^4*BVCstr) & (YA==LAB.BV | YbvA>1))=LAB.BV;       % high-int BV (updated due to cerebellar errors RD20190929)
    
    %% some light region growing
    Ya1(Ya1<1 & Ym<2.2) = nan; 
    [Ya1,Yd] = cat_vol_downcut(Ya1,Ym,2*noise); Ya1(isinf(Ya1) | (Yd>15 & Ya1==LAB.BV) ) = 0; 
    Ya1(Ya1<1 & Ym<1.9) = nan; 
    [Ya1,Yd] = cat_vol_downcut(Ya1,Ym,2); Ya1(isinf(Ya1) | (Yd>15 & Ya1==LAB.BV) ) = 0; 
    Ya1(Ya1<1 & Ym<1.7) = nan; 
    [Ya1,Yd] = cat_vol_downcut(Ya1,Ym,2); Ya1(isinf(Ya1) | (Yd>15 & Ya1==LAB.BV ) ) = 0; 
  
    %% cleanup?
    Ymsk = single(smooth3(Ym)>2.2 | Ya1==LAB.BV) .*  ((Ya1==LAB.CT) + 2*(Ya1==LAB.BV));
    Ymsk = cat_vol_localstat(Ymsk,Ymsk>0,1,1,8); 
    %Ya1a = Ya1; Ya1a(Ymsk>0 & Ymsk<1.5) = LAB.CT; Ya1a(Ymsk>=1.5) = LAB.BV;
    Ya1(Ymsk>0 & Ymsk<1.5) = LAB.CT; Ya1(Ymsk>=1.5) = LAB.BV; clear Ymsk

  else
    Ya1((Ya1==0 & Yp0<1.5 & Ym<1.5 & Yp0>1.3 & Ym>1.3) & YA==LAB.BV)=LAB.BV; % low-int VB  (updated due to cerebellar erros RD20190929)
    Ya1((Ya1==0 & Yp0>3.0 & Ym>3.2) & YA==LAB.BV)=LAB.BV;                    % high-int VB (updated due to cerebellar erros RD20190929)
    Ya1(Ya1==0 & Ym<1.9)=nan; 
    Ya1 = cat_vol_downcut(Ya1,Ym,4*noise); Ya1(isinf(Ya1))=0; 
  end
  Ya1(YA==LAB.TH & Ym>1.75 & Ym<2.85 & Ydiv>-0.1 & Yg<.2 & ~cat_vol_morph(Ya1==LAB.MB | Ya1==LAB.VT | Ya1==LAB.HC | Ya1==LAB.PH,'dd',4,vx_vol))=LAB.TH;                   % thalamus
  Ya1 = cat_vol_median3c(Ya1,Ya1>0 & ~YA==LAB.BV);                                          % smoothing
  clear Ybg Ybgd; 
  Ya1((Yp0>1.75 & Ym>1.75 & Yp0<2.5 & Ym<2.5) & Ya1==LAB.MB) = 0;          % midbrain correction
  Ya1(Ya1==LAB.CT & ~cat_vol_morph(Ya1==LAB.CT,'do',1.4)) = 0; 
  if 0
    %% display for development
    [YD,YI] = cat_vbdist(single(Ya1)); Yax = Ya1(YI);
    figure, isosurface(smooth3(Ya1 & Ym>2.1),.5,Yax), axis equal off; colormap cool, clim([1,8])
    %%
    figure, isosurface(smooth3((Yd<50 | Ym<2.1) .* (Ym>2.1)),.5,Yax), axis equal off; colormap cool, clim([1,8])
    
  end

  % correction of structures that should be compact
  Ya1(Ya1==LAB.BS & ~cat_vol_morph(cat_vol_morph(Ya1==LAB.BS,'o'),'c')) = 0; 
  Ya1(Ya1==LAB.MB & ~cat_vol_morph(cat_vol_morph(Ya1==LAB.MB,'o'),'c')) = 0; 
  Ya1(Ya1==LAB.PH) = 0; 
  Ya1 = cat_vol_median3c(Ya1,Ya1>0 & ~YA==LAB.BV);      


  %% Mapping of ventricles:
  %  Ventricle estimation with a previous definition of non ventricle CSF
  %  to have a second ROI in the region-growin. Using only the ventricle
  %  ROI can lead to overgrowing. Using of non ventrilce ROI doesn't work
  %  because dartle failed for large ventricle. 
  %  It is important to use labopen for each side!
  stime = cat_io_cmd('  Ventricle detection','g5','',verb,stime); 
  %% RD202501: added parahypocampus here to avoid overgrowing ventricles and filling issues and defects 
  Yph = Ym>2 & cat_vol_smooth3X(YA==LAB.PH | YA==LAB.HC,2)>.1; 
  Ynv = cat_vol_morph(~Yb,'d',4,vx_vol) | cat_vol_morph(YA==LAB.CB,'dd',10,vx_vol) | ...% RD202501: extend CB 
        cat_vol_morph(YA==LAB.BS,'dd',2,vx_vol) | ...
        (cat_vol_morph(smooth3(YA==LAB.NV | YA==LAB.TH)>.8,'dd',2,vx_vol) & Ym<1.8) | ...
        cat_vol_morph(YA==LAB.TH & Ym<2.5,'e',1,vx_vol); 
  Ynv = Ynv | Yph | (~cat_vol_morph((Yp0>0 & Yp0<1.5 & (YA==LAB.VT)),'dd',20,vx_vol) & Yp0<1.5);
  Ynv = single(Ynv & Ym<2 & ~cat_vol_morph(Yp0<2 & (YA==LAB.VT) & Yg<0.2,'d',4,vx_vol) & ~Yph);
  Ynv = single(Ynv | (cat_vol_morph(smooth3(YA==LAB.NV | YA==LAB.TH)>.8,'dd',1,vx_vol) & Ym<1.8)); 
  Ynv = Ynv & (~cat_vol_morph(Ya1==LAB.HC | Ya1==LAB.PH,'d',8) | Ya1==LAB.NV);
  Ynv = smooth3(round(Ynv))>0.5;
  % between thamlamus
  Ynv = Ynv | (cat_vol_morph(Ya1==LAB.TH,'c',10) & Yp0<2) | YA==LAB.CB | YA==LAB.BS;
  Ynv = smooth3(Ynv)>0.8;
  Yvt = single(smooth3(Yp0<1.5 & (YA==LAB.VT) & Yg<0.25 & ~Ynv)>0.7); 
  Yvt(cat_vol_morph( (YA==LAB.HC | YA==LAB.PH) & ~(Ya1==LAB.HC | Ya1==LAB.PH) & Yg<noise*2 & Ydiv>-.1 & Ym<1.2 & Yp0<1.2,'do',3)) = 1; 
  Yvt(Yvt==0 & ~(smooth3(Ya1==LAB.HC)>.1 & Ym>1.1) & cat_vol_morph( smooth3(Ya1==LAB.HC)>.7 | Yvt,'dc',7,vx_vol) & Ym<1.25 & Yg<0.3)=1;
  Yvt(cat_vol_morph(Yvt,'l',[inf,10])==0 & Yvt>0) = 0; 
  Yvt(Yvt==0 & Ynv)=2; Yvt(Yvt==0 & Ym>1.8)=nan; Yvt(Yvt==0)=1.5;
  Yvt(Yvt==0 & Ydiv<-0.1) = nan; Yvt(~Yb) = nan;

  %% subcortical stroke lesions
  if exist('Ylesionmsk','var'), Yvt(Ylesionmsk>0.5) = nan; end
  Yvt( cat_vol_morph(YslA>0.6 & Ym<2 & Ydiv./(Ym+eps)>0 & Yp0A>2.5   , 'do', 1 ) ) = 2; % WM
  Yvt( cat_vol_morph(YslA>0.2 & Ym<2 & Ydiv./(Ym+eps)>0 & YA==LAB.BG , 'do', 1 ) ) = 2; % BG lesions
  if exist('Ydt','var') && exist('Ydti','var')
    % by deformation
    Ystsl = (( cat_vol_smooth3X( 1-single(cat_vol_morph(YA==LAB.VT,'dd',10,vx_vol)) ,4 ) .* ...
               cat_vol_smooth3X( single((Ydt - Ydti)>0.7),4) )>0.1 & Ym>0.5 & Ym<2.5 ); 
%%            
    Ystsl = Ystsl | ...
            ((cat_vol_smooth3X(1-single(cat_vol_morph(YA==LAB.VT,'dd',5,vx_vol)),4) .* ...
              cat_vol_smooth3X(single((Ydt - Ydti)>0.7),4))>0.05 & Ym>1.5 & Ym<2.8 & ...
              ~cat_vol_morph(YA==LAB.BG | YA==LAB.TH,'d',2));
    Ystsl = Ystsl | ...
            (cat_vol_smooth3X(single((1/(Ydt+eps) - 1/(Ydti+eps))>0.7),4)>0.4 & Ym>1.5); 
    Ystsl =cat_vol_morph(Ystsl,'dc',4,vx_vol); 
    Yvt(  cat_vol_morph(Ystsl,'e',2)) = 2; 
  else 
    Ystsl = false(size(Ynv)); 
  end
       
  %% bottleneck
  Yvt2 = cat_vol_laplace3R(Yvt,Yvt==1.5,0.01); % first growing for large regions
  Yvt(cat_vol_morph(Yvt2<1.4,'o',2) & ~isnan(Yvt) & Yp0<1.5) = 1; 
  Yvt(cat_vol_morph(Yvt2>1.6,'o',2) & ~isnan(Yvt) & Yp0<1.5) = 2; 
  Ygx = Ym/3 ./ cat_vol_localstat(Ym/3,Ym>1.125,1,3,1); Yvt(Yvt==1.5 & Ygx<.5)=nan; 
  Ygx = Ym/3 ./ cat_vol_localstat(Ym/3,Ym>1.25,1,3,1); Yvt(Yvt==1.5 & Ygx<.8)=nan; 
  Yvt2 = cat_vol_laplace3R(Yvt,Yvt==1.5,0.001);
  % remove small objects
  warning('off','MATLAB:cat_vol_morph:NoObject');
  Yvt = cat_vol_morph(Yvt2<1.5, 'l', [10 0.1]);
  Yvt = cat_vol_morph(Yvt     , 'l', [10 50]);
  warning('on','MATLAB:cat_vol_morph:NoObject');
  Yvt = smooth3((Yvt | (YA==LAB.VT & Ym<1.7)) & Yp0<1.5 & Ym<1.5)>0.5; 
  %%
  Ya1(Yvt) = LAB.VT; Yvtp = cat_vol_morph( Yvt ,'dd', 2) ; 
  Ya1( Ya1 == 0 & YA==LAB.HC & Yvtp & Yp0>1.9 & Ym<3.1 & Ym>1.9 ) = LAB.HC; 
  Ya1( Ya1 == 0 & YA==LAB.PH & Yvtp & Yp0>1.9 & Ym<3.1 & Ym>1.9 ) = LAB.PH; 
  Ya1( Ya1 == 0 & YA==LAB.CT & Yvtp & Yp0>1.9 & Ym<3.1 & Ym>1.9 ) = LAB.CT; 
  if ~debug, clear Yvts1; end

  
  
  %% Mapping of blood vessels
  %  For this we may require the best resolution!
  %  first a hard regions growing have to find the real WM-WM/GM region
  if BVCstr
    stime = cat_io_cmd('  Blood vessel detection','g5','',verb,stime); 
    Ywm = Yp0>2.25 & Ym>2.25 & Yp0<3.1 & Ym<4;                               % init WM 
    Ywm = Ywm | (cat_vol_morph(Ywm,'dd') & Ym<3.5); 
    %%
    Ywm = single(cat_vol_morph(Ywm,'lc',2,vx_vol));                        % closing WM               
    Ywm(smooth3(single(Ywm))<0.5)=0;                                       % remove small dots
    Ywm(~Ywm & (Yp0<0.5 | Ym<1.2 | Ym>4))=nan;                             % set regions growing are
    [Ywm1,YDr] = cat_vol_downcut(Ywm,Ym,2*noise*(1-BVCstr/2));             % region growing
    Ywm(Ywm==-inf | YDr>20)=0; Ywm(Ywm1>0)=1; clear Ywm1                   % set regions growing
    % smoothing
    Ywms = smooth3(single(Ywm)); Yms=smooth3(Ym);                         
    Ywm(Ywms<0.5)=0; Ywm(Ywms>0.5 & Yb & (Ym-Yms)<0.5)=1;                 
    Ywm(Ywms<0.5 & Yb & (Ym-Yms)>0.5)=0; clear Ywms                       
    %% set blood vessels
    Ybv=cat_vol_morph( (Ym>3.75-(0.5*BVCstr) & Yp0<2+(0.5*BVCstr)) | ... % high intensity, but not classified as WM (SPM)
      (Yms>2.5 & (Ym-Yms)>0.6) | ...                                     % regions that strongly change by smoothing
      (Ym>2.5-(0.5*BVCstr) & Ywm==0) | ...                               % high intensity, but not classified as WM (SPM)
      0 ...(Ym>2.5-(0.5*BVCstr) & Yp0<2+(0.5*BVCstr) & Ya1==0 & YA==LAB.CT),'c',1,vx_vol) & ... RD 201901 ADNI 128S0216 error
      ,'c',1,vx_vol) & ...
      cat_vol_morph(Ya1==LAB.CT,'d',2,vx_vol) & ~cat_vol_morph(Ya1==LAB.HC,'d',2,vx_vol) & ...
      cat_vol_morph((Ya1==0 | Ya1==LAB.CB | Ya1==LAB.CT | Ya1==LAB.BV | Ym>1.5) & Ya1~=LAB.VT & Yp0<2.5,'e',1,vx_vol) & ... avoid subcortical regions
      ~Ywm; if ~debug, clear Ywm; end
    Ybb = cat_vol_morph(Yp0>0.5,'lc',1,vx_vol); 
    
    %% RD 201901 ADNI 128S0216 error
 %   Ycenter  = cat_vol_smooth3X(YS==0,2)<0.95 & cat_vol_smooth3X(YS==1,2)<0.95 & Yp0>0;
    Yb2 = Ybv; %smooth3( Ydiv<0.1 & Ym>1.5 & (Ym-Yp0)>0.5 & (Ycenter | cat_vol_morph(~Yb,'dd',3)) & Ya1==0 )>0.5;
    
    %%
    Ybv = ((Ybv | Yb2) & Ybb) | smooth3(Yp0<0.5 & Ybb)>0.4; clear Ybb; 
    %% smoothing
    Ybvs = smooth3(Ybv);
    Ybv(Ybvs>0.3 & Ym>2.5 & Yp0<2.5)=1; Ybv(Ybvs>0.3 & Ym>3.5 & Yp0<2.9)=1;
    Ybv(Ybvs<0.2 & Ym<4-2*BVCstr)=0; clear Yvbs;
    Ya1(Ybv)=LAB.BV; clear Ybv 
  end

  

  %% WMH (White Matter Hyperintensities):
  % WMHs can be found as GM next to the ventricle (A) that do not belong 
  % to a subcortical structure (A) or there must be a big difference 
  % between the tissue SPM expect and the real intensity 'Yp0 - Ym' (C).
  % Furthermore no other Sulic (=near other CSF) should be labeld (D).
  % ####################################################################
  % There can also be deep GM Hyperintensities! 
  % ####################################################################
  % ds('l2','',vx_vol,Ym,Ywmh,Ym/3,Ym/3,90)
  % ####################################################################
  % ToDo: Separate detection of ventricular lesion and subventriculars 
  % ####################################################################
  if extopts.WMHC>=0 && extopts.WMHCstr>=0 && ~extopts.inv_weighting
    % T1 bias correction
    Yi     = Ym .* (Yp0>2.5 & Ym>2.5 & cat_vol_morph(Ya1~=LAB.BG | Ya1~=LAB.VT | Ya1~=LAB.TH,'d',2)); 
    Yi     = cat_vol_median3(Yi,Yi>0,Yi>0); Yi = cat_vol_localstat(Yi,Yi>0,1,3);
    for i=1:2, Yi = cat_vol_localstat(Yi,Yi>0,1,1); end
    Yi     = cat_vol_approx(Yi,'nh',vx_vol,3);
    Ymi    = Ym ./ Yi * 3; 
    
    
    % estimate relative CSF volume 
    Yp0e    = Yp0.*cat_vol_morph(Yb,'e',2); 
    vols    = mean([sum(round(Yp0e(:))==1) sum(round(Yp0e(:))==1 & Yvt(:))] / ...
              sum(round(Yp0e(:))>0.5)); clear Yp0e; 
    noisew  = cat_vol_localstat(smooth3(Ymi)*2,cat_vol_morph(Yp0>2.5,'e') & Ya1==LAB.CT & ~(YwmhA>0.5 & Ymi<2.7),3,4);
    noisew  = cat_stat_nanmean(noisew(noisew(:)>0)); 
    noisec  = cat_vol_localstat(smooth3(Ymi)*2,cat_vol_morph((cat_vol_morph(Yp0>0.5,'e') & ....
              cat_vol_morph(Yp0<2.5,'e')) | cat_vol_morph(Ya1==LAB.VT,'e'),'o'),3,4);
    noisec  = cat_stat_nanmean(noisec(noisec(:)>0));
    if sum(noisec(:)>0)>100, noisew = min(noisew,cat_stat_nanmean(noisec(noisec(:)>0))); end
    
    % control variables 
    % only if there is a lot of CSF and not too much noise
    extopts.WMHCstr = min( 1, max( 0, extopts.WMHCstr ./ max(1,mean(vx_vol)) )); % adaptiv for resolution
    csfvol  = max(eps,min(1.0, (vols  - 0.05) * 10 ));                     % relative CSF volume weighting
    WMHCstr = max(eps,min(1.0, extopts.WMHCstr .* csfvol ));               % normalized WMHCstr 
    wmhvols = 40 - 30 * (1 - extopts.WMHCstr);                             % absolute WMH volume threshold
    mth     = [ min( 1.2 , 1.1 + noisew * (2 - extopts.WMHCstr) ) , ...  % lower tissue threshold
                max( 2.5 , 2.9 - noisew * (2 - extopts.WMHCstr) ) ];     % upper tissue threshold
    %ath     = 2.85 - 0.1 * WMHCstr;                                         % tissue probability threshold
    
    stime   = cat_io_cmd(sprintf('  WMH detection (WMHCstr=%0.02f > WMHCstr''=%0.02f)',...
              extopts.WMHCstr,WMHCstr),'g5','',verb,stime); 

    
    %% creation of helping masks with *YwmhL* and *Ycortex* as most important mask!
    % Ybgth:        subcortical GM regions with WMHs with intensities 
    %               between CSF and GM 
    % Ycenter:      area between the hemispheres that we ignore to avoid 
    %               clossing of the small sulci (eg. close to the CC) 
    % Ycortex1:     initial map to classify cortical GM         
    % Ycortex2:     laplace filtered map to classify cortical GM 
    % Ycortex:      final cortex map that also include CSF and parts of the
    %               normal subcortical structures
    % YwmhL:        extrem large WMs (important to use age-based threshold!)
    % Yinsula:      insula and claustrum 
    Ybgth    = cat_vol_morph( YA==LAB.BG | Ya1==LAB.BS | Ya1==LAB.TH | YA==LAB.PH | YA==LAB.CB ,...
                'dd',1.5 + max(0,2-WMHCstr*4),vx_vol);
    Ybgth    = cat_vol_morph( Ybgth | Ya1==LAB.VT ,'dc',10,vx_vol); 
    Ycenter  = cat_vol_smooth3X(YS==0,8)<0.95 & cat_vol_smooth3X(YS==1,8)<0.95 & Yp0>0;
   
    Ycortex1 = single(1.5 + 0.5*(Yvt2>1.5 & Yvt2<3 & Yp0<=2 & Ym<2.1) - 0.5*Yvt); Ycortex1(Ym>2.1) = nan; 
    %%
    YwmhL    = cat_vol_morph( smooth3((Yp0A + YwmhA)>1.9 & Ym>1.7 & Ym<2.2 & ~Ybgth & ~Ycenter & ...
               cat_vol_morph( YwmhA>0 ,'dd',8) & YA==LAB.CT & ... % use Ywmh atlas
               ~cat_vol_morph(Yvt2>1.6 & Yvt2<3 & Ym<1.8 ,'dd',3,vx_vol))>0.5,'do',4-3*WMHCstr,vx_vol); % age adaptation
             %%
    Ycortex1(cat_vol_morph(YwmhL,'dd',3) & Ycortex1==2) = 0; Ycortex1(YwmhL) = 1; 
    Ycortex2 = cat_vol_laplace3R(Ycortex1,Ycortex1==1.5,0.005);
    Ycortex1( smooth3(Ycortex1==1.5 & Ycortex2>1.5 & Yp0<=2 & Ym<2.1)>0.6) = 2; 
    Ycortex1(isnan(Ycortex1) & Ym>2.8) = 1; 
    Ycortex2 = cat_vol_laplace3R(Ycortex1,Ycortex1==1.5,0.002);
    Ycortex  = cat_vol_morph( (Ycortex2>1.5 & Ycortex2<3) | Ya1==LAB.HC,'dd',2) & Yp0<2.8;
    if ~debug, clear Ycortex1 Ycortex2; end
    Yinsula  = cat_vol_morph( Ycortex,'dd',2,vx_vol) & Yp0<2.8 & ...
               cat_vol_morph( Ya1==LAB.BG | Ya1==LAB.BG,'dd',12,vx_vol); 
    Ycortex  = Ycortex | Yinsula;
    if ~debug, clear Yinsula; end
    
    
    %% create initial WMH map
    %  (A1) classical WMHs + (A2) subcortical (GM) WMHs + (A3) large WMHs
    %  (B1) cortical GM + (B2) deep ventriclal CSF (to avoid the mapping 
    %       of CSF-WM-PVE close to the ventricle)
    %  (C1 & C2) WM regions as boundary for the bottleneck region growing 
    Ywmh = single( smooth3( ...
      cat_vol_morph(~cat_vol_morph(Ycortex,'dd',2 - WMHCstr) & ~Ycenter & (Yp0A+YwmhA)>2.1 & ...
      Ym<2.8 & Ym>1.2 & ~Ybgth & ~cat_vol_morph(Yvt,'d'),'l',[inf,15 - 10*WMHCstr]) > 0 ) > 0.4 ); % (A1)
    Ywmh( smooth3(Ym<1.9 & Ybgth & Yp0A>2)>0.5 & ~cat_vol_morph(Yvt,'o',3) & ...
      ~Ycenter & ~Yvt & ~Ya1==LAB.TH & ...
      ~cat_vol_morph( Ya1==LAB.HC | YA==LAB.PH  ,'dc',10,vx_vol)) = 1; % (A2)
    Ywmh( cat_vol_morph(Ywmh | Yvt,'dc',4,vx_vol) & Ym>1 & Ym<2.8 & Yp0<2.8 & ~Ycenter & ~Yvt ) = 1; 
    Ywmh( YwmhL ) = 1; if ~debug, clear YwmhL; end % (A3)
    Ywmh( smooth3( Ycortex)>0.8 & ~Yvt ) = 2; % (B1)
    Ywmh( cat_vol_morph( Yvt, 'de', 2) ) = 2; % (B2) 
    Ywmh( Ywmh==0 & Ymi>2.8 ) = nan; % (C1)
    Ywmh( Ywmh==0 & Ym>1.9 & cat_vol_morph( Ya1==LAB.BG | Ya1==LAB.TH | ... 
      Ya1==LAB.HC | YA==LAB.PH  ,'dc',10,vx_vol) ) = nan; % (C2)
    
    
    %% remove small dots
    Ywmh(Ywmh==2 & smooth3(Ywmh==2)<0.1 + WMHCstr/2 + noisew) = 0;
    Ywmh(Ywmh==1 & smooth3(Ywmh==1)<0.1 + WMHCstr/2 + noisew) = 0;
    
    
    % Ywmhp:  tiny WMHs that were found by closing the WM 
    if noisew<0.07
      % tiny WMHs as regions with low T1 intensity that can be described 
      % by closing of WM-like regions
      Ywmhp = (Ymi<mth(2) & cat_vol_morph(Ymi>2.2 & Yp0>2.2,'de')) | ...
              (cat_vol_morph(Ymi>mth(2) | cat_vol_morph(Ya1==LAB.VT,'dd',1,vx_vol),'lc',1) & Ymi<mth(2)) | ...
              (cat_vol_morph(Ymi>2.3 | cat_vol_morph(Ya1==LAB.VT,'dd',1,vx_vol),'lc',1) & ...
                Ymi<mth(2) & ~cat_vol_morph(Yp0>2.9,'dd',1.4) & ~Ycortex);
      Ywmhp = Ywmhp & (YwmhA>0 | Yp0A>2.6) & ~Ybgth & Ymi>mth(1) & Ymi<mth(2) & ~Ycenter & ~Ycortex & ...
              ~cat_vol_morph(Ymi<1.7 & Yp0<1.7,'dd',1.8,vx_vol) & ~cat_vol_morph(Ymi<1.5 & Yp0<2.3,'dd',2.5,vx_vol); 
      Ywmhp(cat_vol_morph(Yp0>2.9,'e') & Ym<2.9 & ~Ycortex & ~Ybgth & Ya1~=LAB.BS & Ya1~=LAB.CB) = 1; 
      % no WMHs
      Ywmhp(cat_vol_morph(Ya1==LAB.VT,'dd',1,vx_vol)) = 0;  % not close to the ventricle (PVE range)
      Ywmhp(smooth3(Ywmhp)<0.2 + noisew) = 0;               % avoid WMHs caused by noise
      Ywmhp = cat_vol_morph(Ywmhp,'l',[inf 4 * max(1,min(4,noisew*20))])>0; % remove small WMHs
      Ywmh(smooth3(Ywmhp)>0.5 & Ymi<mth(2)) = 1; 
      if ~debug, clear Ywmhlp; end
    else
      Ywmhp = zeros(size(Ywmh),'single');
    end
    Ywmh( Ywmh==1 & cat_vol_morph(Ywmh==2,'dd',2) ) = 0; 
    
    
    if exist('Yflair','var')
    %% If a flair image is available that we try to use it. 
    %  Although no noise correction is is necessary a bias correction is 
    %  required. Due to our T1 segmentation this is quite easy. Next, we 
    %  can create a FLAIR lesion mask. 
    
      % FLAIR tissue thresholds
      T3thf = [cat_stat_nanmean( Yflair( Yp0toC(Yp0,1)>0.5 )), ...
               cat_stat_nanmean( Yflair( Yp0toC(Yp0,2)>0.5 )), ...
               cat_stat_nanmean( Yflair( Yp0toC(Yp0,3)>0.5 ))]; 
      Tstd   = cat_stat_nanstd( Yflair(Yflair(:)>T3thf(1) & Yflair(:)<max(T3thf(2))) ); 
      
      % FLAIR bias corrections
      Yi     = Yflair > cat_stat_nanmean(T3thf(3)) - Tstd*2 & ...
               Yflair < cat_stat_nanmean(T3thf(3)) + Tstd*2 & ...
               (Yflair/T3thf(3)./(Ymi/3))<2 & Yp0>2.5; 
      Yi     = Yi | (Yp0>2.8 & Ymi>2.8);
      Yg     = Yflair > cat_stat_nanmean(T3thf(2)) - Tstd*3 & Yp0>1.5 & Yp0<2.5 & ~Yi;
      Yi     = cat_vol_localstat(Yflair .* Yi,Yi>0,1,2) + ...
               cat_vol_localstat(Yflair .* Yg,Yg>0,2,3)*T3thf(3)/T3thf(2);  
      Yi = cat_vol_median3(Yi,Yi>0,Yi>0);
      for i=1:2, Yi = cat_vol_localstat(Yi,Yi>0,1,1); end
      Yi     = cat_vol_approx(Yi,'nn',vx_vol,8);
      
      %%
      Yflairn = Yflair./(Yi+eps); 
      T3thf = [cat_stat_nanmean( Yflairn( Yp0toC(Yp0,1)>0.8 )), ...
               cat_stat_nanmean( Yflairn( Yp0toC(Yp0,2)>0.9 & Yflairn>1.1 )), ...
               cat_stat_nanmean( Yflairn( Yp0toC(Yp0,3)>0.9 & Yflairn<1.1 & Ym>2.8 ))]; 
      [T3thfs,T3this] = sort([0, T3thf, max(T3thf) + 1*abs(diff(T3thf(2:3)))]); T3this = [0 1/3 2/3 2/3 2];
      Yflairn2 = zeros(size(Yflairn),'single'); 
      for i=numel(T3thfs):-1:2
        M = Yflairn>T3thfs(i-1) & Yflairn<=T3thfs(i);
        Yflairn2(M(:)) = T3this(i-1) + (Yflairn(M(:)) - T3thfs(i-1))/diff(T3thfs(i-1:i))*diff(T3this(i-1:i));
      end
      M  = Yflairn>=T3thfs(end); 
      Yflairn2(M(:)) = max(T3this) + (Yflairn(M(:)) - T3thfs(end))/diff(T3thfs(end-1:end))*diff(T3this(i-1:i));  
      Yflairn2 = cat_vol_median3(Yflairn2,Yp0>0,Yp0>0,0.1); 
      
      %% create FLAIR mask
      Yflairl = Yflairn2>0.8 & ~Ycenter & ... % flair intensiy 
              cat_vol_morph(cat_vol_morph(Yp0>1 | Yvt,'dc',1.5),'do',1.5); % brain tissue or ventricular area!
              % & ... Yflair./Yi > (T3thf(2)/T3thf(3) + Tstd/T3thf(3)*1.5 - (Tstd/T3thf(3)*extopts.WMHCstr)) & ...  % flair intensiy 
                %Yp0>1.1 & Ymi<2.9 & Yp0A>1.1  & ~Ycenter & ... & YwmhA>eps              % T1 intensities & altas limits 
               % (Yflair./Yi + Ymi/3)>max(2.1,3*(1-YwmhA)); % & ...                    % another flair intensity limit
                %cat_vol_morph(Ya1~=LAB.TH & Ya1~=LAB.BG & Ya1~=LAB.HC,'d',1);     % avoid some regions 
%%                
      %Yflairl = Yflairl | (Yflair./Yi + Ymi/3)>max(2.3,3*(1-YwmhA)) & Ymi<2.9;      % add some regions
      Yflairl = cat_vol_morph(Yflairl,'l',[inf (4 - extopts.WMHCstr)^3])>0;
      Yflairr = smooth3(Yflairl) .* min( 1, 2 * max(0, Yflair./(Yi+eps) - ...
        ( (T3thf(2)/T3thf(3) + Tstd/T3thf(3)*1.5 - (Tstd/T3thf(3)*extopts.WMHCstr)) ) )); 
      if ~debug, clear Yi Yflair; end
      
      %% add FLAIR mask
      Ywmh(Yflairl & (Ywmh<2 | isnan(Ywmh))) = 1; 
    else
      Yflairl = false(size(Ymi)); 
    end
    if ~debug, clear Ycenter; end
    
    % lesions as a regions :
    % - that did not fit to the expected tissue? - not stong enough
    % - with low self/mirror similarity? > Shooting+
    % cat_vol_morph( Ycortex & Yp0<Yp0A-1 & Ym>1.1 & Yp0A>1.8 & ~cat_vol_morph(Yp0<1.2,'do',4) ,'do',1.8); % lesions without FLAIR
    % cat_vol_morph( (Yflairn - (Ym-1)) .* (Ycortex & Yp0<Yp0A-0.2 & Ym>1.1),'do',0.9); % lesion with FLAIR 
    
    %% bottleneck region growing [Manjon:1995]
    Ywmh(Ywmh==0) = 1.5; Ywmh(Ym>2.2 & Ywmh==1.5) = nan;  % harder mask with low T1 threshold
    Ywmh2 = cat_vol_laplace3R(Ywmh, Ywmh==1.5, 0.001);    % bottleneck
    Ywmh(Ywmh==1.5 & Ywmh2>1.7 & Ywmh<3)       = 2;       % add cortex
    Ywmh(Ywmh==1.5 & (Ywmh2<1.2 | Ywmh2==1.5)) = 1;       % add WMHs
    Ywmh(isnan(Ywmh) & Ym<2.5 & Ym>1.5)        = 1.5;     % soft mask with higher T1 threshold
    Ywmh2 = cat_vol_laplace3R(Ywmh, Ywmh==1.5, 0.001);    % bottleneck
    
    % final mask and remove small WMHs
    Ywmh = ~Yvt & Ymi>mth(1) & Ymi<mth(2) & (Ywmh2 < 1.3 | Ywmh2 == 1.5) & ~(Ybgth & Ym>1.8);
    if ~debug, clear Ywmh2; end
    Ywmh = cat_vol_morph(Ywmh , 'l', [inf wmhvols/2])>0;
    
    % final mask and add tiny WMHs and FLAIR WMHs
    Ywmh = ( ( Ywmh | Ywmhp | Yflairl ) & Ymi>mth(1) & Ym<mth(2) & ~cat_vol_morph(Ybgth,'dd',1) );
    Ywmh( Ya1==LAB.CB | Ya1==LAB.BS ) = 0; 
    
    % final intensity (PVE-like segmentation)
    Ywmhr = min( 1, max( 0, Ywmh .* max( Yp0toC(Ym,2) , max( Yp0toC(Ym*8,1.125*8), (Ym>1.125 & Ym<2) ) ) )); 
    if exist('Yflairr','var'), Ywmhr = max( Yflairr .* Ywmh, Ywmhr ); clear Yflairr; end
    if ~debug, clear Ywmhp Ybgth Yflairl Yt; end 
    
    
    %% apply to atlas
    Ya1(Ywmh) = LAB.HI;
  else
    Ywmhr   = false(size(Ym)); 
    Ycortex = false(size(Ym)); 
  end
  if ~debug, clear Ywmh Ynwmh Yvt2; end


  
  %% stroke lesion detection

  % 1. Detection of manually masked regions (zeros within brainmask)
  if exist('Ylesionmsk','var')
    stime = cat_io_cmd('  Manual stroke lesion detection','g5','',verb,stime); 
    Ylesion = Ylesionmsk>0.5; % add manual lesions
  end
  
  % 2. Automatic stroke lesion detection
  %    * use prior maps to identify reginos that should be tissue but
  %      are in stroke areas and have CSF intensity 
  %    * use distance properies to differentiate between normal brain 
  %      atrophy and stroke lesions as local CSF areas that differ stongly
  %      from the average values
  if extopts.SLC
    stime = cat_io_cmd('  Stroke lesion detection','g5','',verb,stime); 
    
    %% large CSF regions without lesion prior
    Ysd       = cat_vbdist(single(cat_vol_morph(~Yb | Yp0>2 | cat_vol_morph(Ya1>0,'ldc',1),'do',3,vx_vol))); 
    mdYsd     = max(3.0,median(Ysd(Ysd(:)>0))); 
    sdYsd     = max(1.5,std(Ysd(Ysd(:)>0))); 
    Yclesion  = smooth3(Ysd>(mdYsd + 3*sdYsd) & Yp0<1.5 & Yp0A>1.25 & Ydiv>-0.1 & (YslA>0 | Yp0A>1.25))>0.5;
    % large CSF regions with lesion prior
    Ysd2      = Ysd .* YslA .* Yp0A/3; 
    mdYsd2    = max(3.0,median(Ysd2(Ysd2(:)>0))); 
    sdYsd2    = max(1.5,std(Ysd2(Ysd2(:)>0))); 
    Yclesion( smooth3(Ysd>(mdYsd + 3*sdYsd) & Yp0<1.5 & Yp0A>1.25 & Ydiv>-0.1 & (YslA>0 | Yp0A>1.25))>0.5 ) = 1;
    if ~debug, clear Ysd; end 
    Yclesion  = cat_vol_morph(cat_vol_morph(Yclesion | ~(Yb & Yp0A>0),'dc',4,vx_vol) & Ym<2 & Yp0>0.5,'do',2); 
    Yclesion  = cat_vol_morph(Yclesion,'l',[inf 400],vx_vol)>0; 
    
    % further lesions
    Ywlesion  = smooth3(abs(Yp0-Ym)/3 .* (1-max(0,abs(3-Yp0A))) .* YslA * 10 )>0.3 & Ym<2.2 & Yp0>1 & Ya1~=LAB.HI; % WM-GM leson
    Ywlesion( smooth3(Yp0A/3 .* YslA .* (3-Ym) .* (Ya1~=LAB.HI))>0.5 ) = 1; % WM lesions
    %
    Ywlesion( cat_vol_morph(YslA>0.6 & Yp0<2.0 & Ydiv./(Ym+eps)>0 & YslA>0.2 & Yp0A>2.5 , 'do', 1 ) ) = 1; % WM
    Ywlesion( cat_vol_morph(YslA>0.2 & Yp0<1.6 & Ydiv./(Ym+eps)>0 & YslA>0.2 & (YA==LAB.BG | YA==LAB.TH) , 'do', 1 ) ) = 1; % BG lesions
    % closing and opening
    Ywlesion  = cat_vol_morph(cat_vol_morph(Ywlesion | ~(Yb & Yp0A>0),'dc',4,vx_vol) & Ym<2 & Yp0>0.5,'do',1); 
    Ywlesion  = cat_vol_morph(Ywlesion,'l',[inf 200],vx_vol)>0; 
    
    %%
    Ystsl     = Ystsl & Ym<2.8 & Ya1~=LAB.VT;
    Yilesion  = single(Yclesion | Ywlesion | Ystsl); Yilesion(Yilesion==0  & ( Yp0<0.5 | Yp0>2.1 | Ydiv<-0.1 | Ya1==LAB.VT ) ) = -inf;
    [Yilesion,Ydd]  = cat_vol_downcut(Yilesion,3-Ym,-0.0001); Yilesion(Ydd>200) = 0;  
    Yilesion  = cat_vol_morph(Yilesion,'dc',1,vx_vol); 
    Yilesion  = cat_vol_morph(Yilesion,'do',4,vx_vol)>0 | Ywlesion| Yclesion; 
    Yilesion  = cat_vol_morph(Yilesion,'l',[inf 200],vx_vol)>0; 
    Ynlesion  = smooth3(~Yilesion & Ysd2<(mdYsd2 + 2*sdYsd2) & Ym<2.0 & Yp0A<1.9 & Ydiv>0)>0.5;
    if ~debug, clear Ysd2; end 
    
    % region-growing  
    Ysl = single(Yilesion) + 2*single(Ynlesion | Yp0<0.5 | Yvt | ~Yb); if ~debug, clear Yilesion Ynlesion; end
    Ysl(Yp0>2.8)=nan; Ysl(Ysl==0) = 1.5; % harder mask with low T1 threshold
    Ysl2 = cat_vol_laplace3R(Ysl, Ysl==1.5, 0.05);  % bottleneck
    Ysl(Ysl2<1.45 & ~isnan(Ysl))=1; Ysl(Ysl2>1.55 & ~isnan(Ysl))=2; 
    Ysl2 = cat_vol_laplace3R(Ysl, Ysl==1.5, 0.01); if ~debug, clear Ysl; end    % bottleneck
    Ylesion  = cat_vol_morph(Ysl2<1.45 & Ysl2>0 & Ym<2.8 & Ya1~=LAB.HI,'do',1,vx_vol); if ~debug, clear Ysl2; end
    %%
    
    %%
    %Ylesion = single( Yp0A./(Ym+2)>1.3 & Yp0>0.5 & Yp0<2 & Ya1~=LAB.VT & Ya1~=LAB.HI & Ym<1.5 ); % , 'do', 3);
    %Ylesion = single(cat_vol_morph(Ylesion,'l',[inf 50])>0); 
    %{
    Ylesion( cat_vol_morph(YslA>0.6 & Ym<2 & Ydiv./Ym>0 & Yp0A>2.5   , 'do', 1 ) ) = 2; % WM
    Ylesion( cat_vol_morph(YslA>0.2 & Ym<2 & Ydiv./Ym>0 & YA==LAB.BG , 'do', 1 ) ) = 2; % BG lesions
    Ylesion(smooth3(Yp0A>2.9 & Ym<2)>0.7) = 1; 
    
    Ylesion(Ym>2 | (Ya1==LAB.VT & Yp0A<2.9) | Ya1==LAB.HI | Yp0==0) = nan;
    [Ylesion,Ydx] = cat_vol_simgrow(Ylesion,max(1,Ym),1); Ylesion(Ydx>0.1)=0; 
    % different volume boundaries depending on the position of the lesion
    Ylesion = cat_vol_morph(cat_vol_morph(Ylesion,'do',2) & Yp0A<2.5,'l',[inf 200])>0 | ... % maybe just atrophy
              cat_vol_morph(cat_vol_morph(Ylesion,'do',1) & Yp0A>2.5,'l',[inf 100])>0 | ...
              cat_vol_morph(Ylesion & Yp0A>2.9,'l',[inf 10])>0;
    %}
    % add manual leisons
    if exist('Ylesionmsk','var'), Ylesion(Ylesionmsk>0.5) = 1; end             
  end
  Ya1(Ylesion>0) = LAB.LE;
  
  
  %% Closing of gaps between diffent structures:
  stime = cat_io_cmd('  Closing of deep structures','g5','',verb,stime); 
  Yvtd2 = cat_vol_morph(Ya1==LAB.VT,'dd',2,vx_vol) & Ya1~=LAB.VT;
  % CT and VT
  Ycenter = cat_vol_morph(Ya1==LAB.VT,'dd',2,vx_vol) & ...
       cat_vol_morph(Ya1==LAB.CT,'d',2,vx_vol) & Ya1==0 ;
  Ya1(Ycenter & Yp0<=1.5 & ~Ynv)=LAB.VT; Ya1(Ycenter & Yp0>1.5)=LAB.CT; 
  % WMH and VT
  Ycenter = cat_vol_morph(Ya1==LAB.HI,'dd',2,vx_vol) & Yvtd2 & ~Ynv & Ya1==0;
  Ya1(Ycenter &  Ym<=1.25)=LAB.VT; Ya1(Ycenter & Ym>1.25 & Ym<2.5)=LAB.HI; 
  % TH and VT
  if 0 % RD202501
    Ycenter = cat_vol_morph(Ya1==LAB.TH,'dd',2,vx_vol) & Yvtd2;
    Ya1(Ycenter & Ym<=1.5)=LAB.VT; Ya1(Ycenter & Ym>1.5 & Ym<2.85)=LAB.TH; 
  end
  % BG and VT
  Ycenter = cat_vol_morph(Ya1==LAB.BG,'dd',2,vx_vol) & Yvtd2;
  Ya1(Ycenter & Ym<=1.5)=LAB.VT; Ya1(Ycenter & Ym>1.5 & Ym<2.85)=LAB.BG;
  % no bloodvessels next to the ventricle, because for strong atrophy
  % brains the WM structures can be very thin and may still include 
  % strong bias
  Ya1(Ya1==LAB.BV & cat_vol_morph(Ya1==LAB.VT,'dd',3,vx_vol))=0;
  if ~debug, clear Yt Yh Yvtd2 Yw; end

  
  
  %% complete map
  Ya1(Ya1==0 & Yp0<1.75) = nan; 
  Ya1 = cat_vol_downcut(Ya1,Ym,noise); Ya1(isinf(Ya1)) = 0; 
  Ybv = Ya1==LAB.BV; Ya1(Ya1==LAB.BV) = 0; 
  [~,~,Ya1x] = cat_vbdist(Ya1,Yb); 
  Ya1(Ya1x==LAB.CB | Ya1x==LAB.BS | Ya1x==LAB.MB | Ya1x==LAB.CT) = Ya1x(Ya1x==LAB.CB | Ya1x==LAB.BS | Ya1x==LAB.MB | Ya1x==LAB.CT);
  Ya1(Ya1x>0 & Ya1==0) = 1; 
  Ya1(Ybv) = LAB.BV; clear Ybv; 
  
  % consider gyrus parahippocampalis | cat_vol_morph(Ya1==LAB.HC,'e')
  % RD202501: add defintion of parahippocampal gyrus (with some extensive 
  %           processing to really have a whole free thing but we will try to first keep it simple)
  Yph = cat_vol_morph((YA==LAB.PH ),'dd',1.9) & cat_vol_morph(Ya1==LAB.VT | Ya1==LAB.HC,'dd',4,vx_vol) & Ym>2.125 & Ydiv<0.05; 
  Ya1(Yph>0) = LAB.PH; clear Yph; % parahippocampus
  
  %% side aligment using laplace to correct for missalignments due to the normalization
  stime = cat_io_cmd('  Side alignment','g5','',verb,stime); 
  YBG  = Ya1==LAB.BG | Ya1==LAB.TH;
  YMF  = Ya1==LAB.VT | Ya1==LAB.BG | Ya1==LAB.HI | Ya1==LAB.PH | (Ya1==LAB.TH & cat_vol_smooth3X(Ym>1.9)); % add the thalamus (RD20190913)
  YMF(smooth3(YMF)<.5) = 0; 
  YMF  = cat_vol_morph(cat_vol_morph(YMF,'dc',2),'do'); 
  YMF2 = cat_vol_morph(YMF,'dd',2,vx_vol) | Ya1==LAB.CB | Ya1==LAB.BS | Ya1==LAB.MB;
  Ymf  = max(Ym,smooth3(single(YMF2*3))); 
  Ycenter = cat_vol_smooth3X(YS==0,6)<0.9 & cat_vol_smooth3X(YS==1,6)<0.9 & ~YMF2 & Yp0>0 & Ym<3.1 & (Yp0<2.5 | Ya1==LAB.BV);
  Ys = (2-single(YS)) .* single(smooth3(Ycenter)<0.4);
  Ys(Ys==0 & (Ym<1 | Ym>3.1))=nan; Ys = cat_vol_downcut(Ys,Ymf,0.1,vx_vol); 
  [~,~,Ys] = cat_vbdist(Ys,Ys==0); 
  if ~debug, clear YMF2 Yt YS; end
  
  % YMF for FreeSurfer fsaverage
  Ysm  = cat_vol_morph(Ys==2,'d',1.5,vx_vol) & cat_vol_morph(Ys==1,'d',1.5,vx_vol); 
  Ynn  = 0 * cat_vol_morph(Ysm & Ya1==LAB.MB,'dd',10);
  YMF  = Ya1==LAB.VT | (cat_vol_morph(Ya1==LAB.PH | Ya1==LAB.BG | Ya1==LAB.HI | (Ya1==LAB.TH & cat_vol_smooth3X(Ym)>1.9),'dc',2,vx_vol) & ~Ysm); % changed thalamus (RD20190913)
  YMF  = Ya1~=LAB.CB &  ~Ynn & Ym<=2.75  & cat_vol_morph(YMF | Ym>2.3,'c',1) & cat_vol_morph(YMF,'dd',2,vx_vol);
  YMF  = smooth3(YMF)>0.5;
  Ycenter = cat_vol_morph(Ya1==LAB.TH | Ya1==LAB.VT,'dc',4,vx_vol) & ~(Ya1==LAB.TH | Ya1==LAB.VT | cat_vol_morph(Ya1==LAB.HC | Ya1==LAB.PH,'dd',2,vx_vol));
  YMF(Ycenter)=1; 
  clear Ysm; 
  
  
  %% back to original size
  stime = cat_io_cmd('  Final corrections','g5','',verb,stime); 
  Ya1   = cat_vol_resize(Ya1,'dereduceV',resTr,'nearest'); Ya1 = cat_vol_median3c(Ya1,Ya1>0 & Ya1~=LAB.BV);
  Ys    = cat_vol_resize(Ys ,'dereduceV',resTr,'nearest'); Ys  = 1 + single(smooth3(Ys)>1.5);
  YMF   = cat_vol_resize(single(YMF),'dereduceV',resTr)>0.5;
  YBG   = cat_vol_resize(single(YBG),'dereduceV',resTr)>0.5;
  Ywmhr = cat_vol_resize(Ywmhr,'dereduceV',resTr);
  Ycortex = cat_vol_resize(single(Ycortex),'dereduceV',resTr)>0.5;
  
  Ya1   = cat_vol_resize(Ya1,'dereduceBrain',BB); Ya1 = cat_vol_ctype(Ya1);
  Ys    = cat_vol_resize(Ys ,'dereduceBrain',BB); [~,~,Ys] = cat_vbdist(Ys,Ya1>0);
  YMF   = cat_vol_resize(YMF,'dereduceBrain',BB); 
  YBG   = cat_vol_resize(YBG,'dereduceBrain',BB); 
  Ywmhr = cat_vol_resize(Ywmhr,'dereduceBrain',BB); 
  Ycortex = cat_vol_resize(Ycortex,'dereduceBrain',BB); 
  Ym    = Ym0; clear Ym0;

  % final side alignment
  Ya1(Ya1>0)=Ya1(Ya1>0)+(Ys(Ya1>0)-1);
 
  
  %% correction of tissue classes 
  
  % add WMH class
  Ywmhrd  = cat_vol_morph(Ywmhr,'dd');
  Yclssum = (single(Ycls{1}) + single(Ycls{3})) .* (Ywmhrd);
  Ycls{7} = cat_vol_ctype(Yclssum); 
  Ycls{1} = cat_vol_ctype(single(Ycls{1}) .* (~Ywmhrd)); 
  Ycls{3} = cat_vol_ctype(single(Ycls{3}) .* (~Ywmhrd)); 
  clear Ywmhrd
  
  % set possible blood vessels to class 4
  NS = @(Ys,s) Ys==s | Ys==s+1;
  if ~NBVC
    % this should be done later with consideration of the BV intensity and WM distance! 
    Ybv     = NS(Ya1,LAB.BV); 
    Yclssum = (single(Ycls{1}) + single(Ycls{2})) .* (Ybv);
    Ycls{5} = cat_vol_ctype(single(Ycls{5}) + Yclssum); 
    Ycls{1} = cat_vol_ctype(single(Ycls{1}) .* (~Ybv)); 
    Ycls{2} = cat_vol_ctype(single(Ycls{2}) .* (~Ybv)); 
    clear Ybv; 
  end
  
  
  % YBG is smoothed a little bit and (B) reset all values that are related with GM/WM intensity (Ym<2.9/3) (A)
  Yclssum = single(Ycls{1}) + single(Ycls{2}) + single(Ycls{3});
  YBGs    = min( max(0,min(255, 255 - cat_vol_smooth3X(Ya1==1 & Ycls{2}>round(2.9/3),0.8) .* single(Ycls{2}) )), ... (A)
                 max(0,min(255, 255 * cat_vol_smooth3X(YBG .* (Ym<=2.9/3 & Ym>2/3) ,0.5) )) ); % (B)
  Ycls{1} = cat_vol_ctype(single(Ycls{1}) + YBGs .* (single(Ycls{2})./max(eps,Yclssum)));
  Ycls{2} = cat_vol_ctype(single(Ycls{2}) - YBGs .* (single(Ycls{2})./max(eps,Yclssum)));
  clear YBGs Yclssum; 
 
  % assure that the sum of all tissues is 255 
  Yclss = zeros(size(Ym),'single'); 
  for ci=1:numel(Ycls), Yclss = Yclss + single(Ycls{ci}); end
  for ci=1:numel(Ycls), Ycls{ci} = cat_vol_ctype(single(Ycls{ci}) ./ max(eps,Yclss) * 255); end
  
  
  cat_io_cmd(' ','','',verb,stime); 
  
end
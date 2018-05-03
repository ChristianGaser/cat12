function [Ya1,Ycls,YBG,YMF] = cat_vol_partvol(Ym,Ycls,Yb,Yy,vx_vol,extopts,Vtpm,noise,job)
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
%   - Zuweisung durch Dartel mit hoher Genauigkeit m??glich 
%   - Erweiterung von SPM/VBM durch MainROIs (Seiten, Lappen, ...)
%   - Verbesserung der SPM/VBM durch bessere Enfernung von unerw??nschtem
%     Gewebe (ON, Blutgef????e ...)
%   - Blutgef????e k??nnnen als erweitere Masken f??r fMRI genutzt werden um
%     Seiteneffekte besser ausblenden zu k??nnen.
%  [- Beliebige Atlanten k??nnen genutzt werden.]
%
%  Todo:
%   - Besserer Atlas
%   - BV vs. HD - gl??tten in dilated HD region
%   - F??llen von CSF l??cken bei LAB~=BV und Ym<1.2 und LAB==NV?
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
  BVCstr  = extopts.BVCstr; 
  verb    = extopts.verb-1;
  PA      = extopts.cat12atlas;
  vx_res  = max( extopts.uhrlim , max( [ max(vx_vol) min(vx_vol) ] )); 
  noise   = double(noise);

  %% map atlas to RAW space
  if verb, fprintf('\n'); end
  stime = cat_io_cmd('  Atlas -> subject space','g5','',verb); dispc=1;
  % CAT atlas
  VA = spm_vol(PA{1});
  YA = cat_vol_ctype(spm_sample_vol(VA,double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0));
  YA = reshape(YA,size(Ym));
  
   % template map
  Yp0A = single(spm_sample_vol(Vtpm(1),double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),1))*2 + ...
         single(spm_sample_vol(Vtpm(2),double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),1))*3 + ...
         single(spm_sample_vol(Vtpm(3),double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),1))*1;
  Yp0A = reshape(Yp0A,size(Ym)); 
  
  %% WMH atlas
  watlas = 3; 
  switch watlas
    case 1, PwmhA = strrep(PA{1},'cat.nii','cat_wmh_soft.nii');
    case 2, PwmhA = strrep(PA{1},'cat.nii','cat_wmh.nii');
    case 3, PwmhA = strrep(PA{1},'cat.nii','cat_wmh_miccai2017.nii');
  end
  if exist(PwmhA,'file') && ~strcmp(PwmhA,PA{1}) 
    VwmhA = spm_vol(PwmhA);
    YwmhA = spm_sample_vol(VwmhA,double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0);
    YwmhA = reshape(YwmhA,size(Ym));
  else
    YwmhA = max(0,min(1,Yp0A-2)); 
  end
  switch watlas
    case 2, YwmhA = min(1,max(0,YwmhA - 0.1) * 0.8);
    case 3, YwmhA = min(1,max(0,cat_vol_smooth3X(YwmhA,1) - 0.01) * 10); 
  end
  clear Yy2; 

  if exist('job','var') && isfield(job,'data_wmh') && ~isempty(job.data_wmh)
    stime = cat_io_cmd('  FLAIR corregistration','g5','',verb,stime); dispc=dispc+1;
    
    % coreg
    Vflair = spm_vol(job.data_wmh{job.subj}); 
    Vm     = spm_vol(job.data{job.subj}); 
    evalc('R = spm_coreg(Vm,Vflair,struct(''graphics'',0));'); 
    R      = spm_matrix(R);  %#ok<NODEF>
    
    %% load
    for i=1:Vm.dim(3)
      Yflair(:,:,i) = single( spm_slice_vol(Vflair, R \ Vflair.mat \ Vm.mat  * spm_matrix([0 0 i]) ,Vm.dim(1:2),[1,NaN])); 
    end    

    
  end
  
  %%
  if ~debug; clear Yy; end
  
  Yp0  = (single(Ycls{1})*2/255 + single(Ycls{2})*3/255 + single(Ycls{3})/255) .* Yb; 
  
  % work on average resolution
  Ym0 = Ym; 
  % remove background
  [Ym,BB]  = cat_vol_resize(Ym      ,'reduceBrain',vx_vol,2,Yb);
  YA       = cat_vol_resize(YA      ,'reduceBrain',vx_vol,2,Yb);
  Yp0      = cat_vol_resize(Yp0     ,'reduceBrain',vx_vol,2,Yb);
  Yp0A     = cat_vol_resize(Yp0A    ,'reduceBrain',vx_vol,2,Yb);
  YwmhA    = cat_vol_resize(YwmhA   ,'reduceBrain',vx_vol,2,Yb);
  if exist('Yflair','var')
    Yflair = cat_vol_resize(Yflair  ,'reduceBrain',vx_vol,2,Yb);
  end
  Yb       = cat_vol_resize(Yb      ,'reduceBrain',vx_vol,2,Yb);
  % use lower resolution 
  [Ym,resTr] = cat_vol_resize(Ym    ,'reduceV',vx_vol,vx_res,64);
  YA         = cat_vol_resize(YA    ,'reduceV',vx_vol,vx_res,64,'nearest'); 
  Yp0        = cat_vol_resize(Yp0   ,'reduceV',vx_vol,vx_res,64);
  Yp0A       = cat_vol_resize(Yp0A  ,'reduceV',vx_vol,vx_res,64);
  YwmhA      = cat_vol_resize(YwmhA ,'reduceV',vx_vol,vx_res,64); 
  if exist('Yflair','var')
    Yflair   = cat_vol_resize(Yflair,'reduceV',vx_vol,vx_res,64); 
  end    
  Yb         = cat_vol_resize(Yb    ,'reduceV',vx_vol,vx_res,64);
  vx_vol     = resTr.vx_volr; 
  
  
  % noise reduction
  spm_smooth(Ym,Ym,0.6./vx_vol);
  
  
  % prepare maps
  [tmp0,tmp1,YS] = cat_vbdist(single(mod(YA,2)) + single(YA>0)); YS=~mod(YS,2); clear tmp0 tmp1;  % side map
  YA(mod(YA,2)==0 & YA>0)=YA(mod(YA,2)==0 & YA>0)-1;                       % ROI map without side
  YA   = cat_vol_ctype(cat_vol_median3c(single(YA),Yp0>0));                % noise filter
  Yg   = cat_vol_grad(Ym,vx_vol);                                          % gadient map (edge map)
  Ydiv = cat_vol_div(Ym,vx_vol);                                           % divergence map (edge map)
  Ym   = Ym*3 .* (Yb);
  Yb   = Yb>0.5;

  
  
  %% Create individual mapping:
  stime = cat_io_cmd('  Major structures','g5','',verb,stime); dispc=dispc+1;
  
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
  Ya1(Ybg)=LAB.BG;                                                         % basal ganglia
  Ya1(YA==LAB.TH & Ym>1.9 & Ym<2.85 & Ydiv>-0.1)=LAB.TH;                   % thalamus
  Ya1(YA==LAB.HC & Ym>1.9 & Ym<2.85 & Ydiv>-0.1)=LAB.HC;                   % hippocampus
  Ya1(((Yp0>2.5 & Ym>2.5 & (YA==LAB.CT | YA==LAB.BG)) | ...
       (Yp0>1.5 & Ym<3.5 & Ybgd>1 & Ybgd<8 & (Ybgd>4 | ...
       Ydiv<-0.02+Ybgd/200))) & Ya1==0)=LAB.CT;                            % cerebrum
  Ya1((Yp0>2.0 & Ym>2.0) & YA==LAB.CB)=LAB.CB;                             % cerebellum
  Ya1((Yp0>2.0 & Ym>2.0) & YA==LAB.BS)=LAB.BS;                             % brainstem
  Ya1((Yp0>2.0 & Ym>2.0) & YA==LAB.ON)=LAB.ON;                             % optical nerv
  Ya1((Yp0>2.0 & Ym>2.0) & YA==LAB.MB)=LAB.MB;                             % midbrain
  clear Ybg Ybgd; 
  %% region-growing
  Ya1(Ya1==0 & Yp0<1.9)=nan; 
  Ya1 = cat_vol_downcut(Ya1,Ym,4*noise); Ya1(isinf(Ya1))=0; 
  Ya1 = cat_vol_median3c(Ya1,Yb);                                          % smoothing
  Ya1((Yp0>1.75 & Ym>1.75 & Yp0<2.5 & Ym<2.5) & Ya1==LAB.MB)=0;            % midbrain correction
  
  
  
  %% Mapping of ventricles:
  %  Ventricle estimation with a previous definition of non ventricle CSF
  %  to have a second ROI in the region-growin. Using only the ventricle
  %  ROI can lead to overgrowing. Using of non ventrilce ROI doesn't work
  %  because dartle failed for large ventricle. 
  %  It is important to use labopen for each side!
  stime = cat_io_cmd('  Ventricle detection','g5','',verb,stime); dispc=dispc+1;
  Ynv = cat_vol_morph(cat_vol_morph(~Yb,'d',4,vx_vol) | ...
        (YA==LAB.CB | YA==LAB.BS),'d',2,vx_vol) | cat_vol_morph(YA==LAB.NV,'e',1,vx_vol);
  Ynv = Ynv | ~cat_vol_morph(Yp0>0 & Yp0<1.5 & (YA==LAB.VT),'d',8,vx_vol);
  Ynv = single(Ynv & Ym<2 & ~cat_vol_morph(Yp0<2 & (YA==LAB.VT) & Yg<0.2,'d',4,vx_vol));
  Ynv = smooth3(round(Ynv))>0.5; 
  % between thamlamus
  Ynv = Ynv | (cat_vol_morph(Ya1==LAB.TH,'c',10) & Yp0<2) | YA==LAB.CB | YA==LAB.BS;
  Ynv = smooth3(Ynv)>0.8;
  Yvt = single(smooth3(Yp0<1.5 & (YA==LAB.VT) & Yg<0.25 & ~Ynv)>0.7); 
  Yvt(Yvt==0 & Ynv)=2; Yvt(Yvt==0 & Ym>1.8)=nan; Yvt(Yvt==0)=1.5;
  Yvt2 = cat_vol_laplace3R(Yvt,Yvt==1.5,0.005);
  % remove small objects
  warning('off','MATLAB:cat_vol_morph:NoObject');
  Yvt2 = cat_vol_morph(Yvt2<1.5, 'l', [10 0.1]);
  Yvt2 = cat_vol_morph(Yvt2    , 'l', [10 50]);
  warning('on','MATLAB:cat_vol_morph:NoObject');
  Yvt = smooth3((Yvt2 | (YA==LAB.VT & Ym<1.7)) & Yp0<1.5 & Ym<1.5)>0.5; 
  Ya1(Yvt) = LAB.VT; 
  if ~debug, clear Yvts1; end

  
  
  %% Mapping of blood vessels
  %  For this we may require the best resolution!
  %  first a hard regions growing have to find the real WM-WM/GM region
  if BVCstr
    stime = cat_io_cmd('  Blood vessel detection','g5','',verb,stime); dispc=dispc+1;
    Ywm = Yp0>2.5 & Ym>2.5 & Yp0<3.1 & Ym<4;                               % init WM 
    Ywm = Ywm | (cat_vol_morph(Ywm,'d') & Ym<3.5); 
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
      (Ym>2.5-(0.5*BVCstr) & Yp0<2+(0.5*BVCstr) & Ya1==0 & YA==LAB.CT),'c',1,vx_vol) & ...
      cat_vol_morph(Ya1==LAB.CT,'d',2,vx_vol) & ~cat_vol_morph(Ya1==LAB.HC,'d',2,vx_vol) & ...
      cat_vol_morph((Ya1==0 | Ya1==LAB.CT | Ya1==LAB.BV | Ym>1.5) & Ya1~=LAB.VT & Yp0<2.5,'e',1,vx_vol) & ... avoid subcortical regions
      ~Ywm;  clear Ywm 
    Ybb = cat_vol_morph(Yp0>0.5,'lc',1,vx_vol); 
    Ybv = (Ybv & Ybb) | smooth3(Yp0<0.5 & Ybb)>0.4; clear Ybb; 
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
  % ToDo: Seperate detection of ventricular lession and subventriculars 
  % ####################################################################
  if extopts.WMHC>0 && extopts.WMHCstr>0
    %%
    % estimate relative CSF volume 
    Yp0e    = Yp0.*cat_vol_morph(Yb,'e',2); 
    vols    = mean([sum(round(Yp0e(:))==1) sum(round(Yp0e(:))==1 & Yvt(:))] / ...
              sum(round(Yp0e(:))>0.5)); clear Yp0e; 

    % control variables 
    % only if there is a lot of CSF and not too much noise
    csfvol  = max(eps,min(1.5, (vols  - 0.05) * 10 ));                     % relative CSF volume weighting
    WMHCstr = max(eps,min(1.0, extopts.WMHCstr .* csfvol ));               % normalized WMHCstr 
    wmhvols = max([0.01 50],min([0.1,100], ...
                  [0.1 200] - [0.08 100] .* (1 - WMHCstr)));               % WMH volume thresholds
    mth     = max([1.1 + min(0.3 , noise) , 2.6 ] , ...
              min([1.4                    , max(2.7 , 2.8 - noise)], ...
                  [1.4 2.8] + [-0.2 0.2] .* WMHCstr));                     % tissue thresholds
    ath     = max(2.0,min(2.8,2.85 - 0.1 * WMHCstr));                      % tissue probability threshold
    vtd     = max(4,min(12,4 + 4 * WMHCstr));                              % ventricle distance threshold

    stime   = cat_io_cmd(sprintf('  WMH detection (WMHCstr=%0.02f > WMHCstr''=%0.02f)',...
              extopts.WMHCstr,WMHCstr),'g5','',verb,stime); dispc=dispc+1;

    YBG2 = cat_vol_morph(Ya1==LAB.BG,'d',1); 

    Ywmhp = min(2,single(1.5 - 0.5*(Yp0<1.5) +  0.5*(Yp0>2.5) + (Ya1==LAB.BG | Ya1==LAB.TH | Ya1==LAB.VT)));
    Ywmhp = cat_vol_laplace3R(Ywmhp, Ywmhp==1.5, 0.01); 
    
    % set pre WMHs
    Ywmh = smooth3( cat_vol_morph(Yvt,'dd',vtd,vx_vol) & ~YBG2 & Yp0A>ath & ...
                    Ya1==LAB.CT & YwmhA>(0.8 - WMHCstr*0.8) & Ym>mth(1) & Ym<mth(2) )>0.5 ;
    Ywmh(cat_vol_morph(Yvt,'dd',vtd*1.5,vx_vol) & Ywmhp>1.8 & Yp0A>ath & YwmhA>(0.8 - WMHCstr*0.8) & Ym>mth(1) & Ym<mth(2) & ~Yvt & ~YBG2 &  YA==LAB.CT) = 1;
    Ywmh(Ywmhp>1.9 & Ym<min(2.9,mth(2)+0.1) & Ym>2 & cat_vol_morph(Yp0>2.5,'e')) = 1;
    % rwmhvol = sum(Ywmh(:)==1) ./ sum(Yp0(:)>2.5 & Ym(:)>2.5); % ... not yet ... 
    Ywmh = cat_vol_morph(Ywmh,'l',[inf wmhvols(1)/2])>0;
    Ywmh = cat_vol_morph(Ywmh,'l',[inf wmhvols(2)/2])>0;
    % set no-go area
    Ywmh = single(Ywmh); 
    Ywmh((Ywmh==0 & Ym>2.5) | YBG2 | Ya1==LAB.BG | Ya1==LAB.TH | Ya1==LAB.TH | Ya1==LAB.PH) = nan;
    % set non-WMH
    Ywmh(cat_vol_morph(Ya1==LAB.HC,'dd',5) & Ywmh==0) = 2; 
    Ywmh( cat_vol_morph( ( (Yvt2>1.99 & Yvt2<2.05) | (Ywmh==0 & Ym<1.5) ) & ...
         ~cat_vol_morph(Yvt,'dd',10) , 'l',[20 100])>0 ) = 2;
    Ywmh(cat_vol_morph(Ya1==LAB.VT,'e',2)) = 2; 
    Ywmh(Ywmhp<1.5 & ~cat_vol_morph(Yvt,'dd',vtd*1.5,vx_vol) & Yp0<1.6 & Ym<1.6) = 2; 
    
    % remove small dots
    Ywmh(Ywmh==2 & smooth3(Ywmh==2)<0.1 + WMHCstr/2) = 0;
    Ywmh(Ywmh==1 & smooth3(Ywmh==1)<0.1 + WMHCstr/2) = 0;
    
    %%
    if exist('Yflair','var')
      Yp0toC  = @(Yp0,c) 1-min(1,abs(Yp0-c));
      
      % FLAIR tissue thresholds
      T3thf = [cat_stat_nanmean( Yflair( Yp0toC(Yp0,1)>0.5 )), ...
               cat_stat_nanmean( Yflair( Yp0toC(Yp0,2)>0.5 )), ...
               cat_stat_nanmean( Yflair( Yp0toC(Yp0,3)>0.5 ))]; 
      Tstd   = cat_stat_nanstd( Yflair(Yflair(:)>T3thf(1) & Yflair(:)<max(T3thf(2:3))) ); 
      
      % FLAIR bias corrections
      Yi     = Yflair > cat_stat_nanmean(T3thf(2:3)) - Tstd*2 & ...
               Yflair < cat_stat_nanmean(T3thf(2:3)) + Tstd*2 & ...
               (Yflair/T3thf(3)./(Ym/3))<2 & Yp0>1.5; 
      Yi     = Yflair .* Yi;        
      for i=1:2, Yi = cat_vol_localstat(Yi,Yi>0,1,1); end
      Yi     = cat_vol_approx(Yi,'nh',vx_vol,4);
      Yt = cat_vol_smooth3X(YS==0,6)<0.9 & cat_vol_smooth3X(YS==1,6)<0.9 & Yp0>0;

      % create FLAIR mask
      Yflairl = Yflair./Yi > 1.3-(0.1*WMHCstr) & Ym>0.5 & Ym<2.8 & Yp0A>2.3 & YwmhA>eps & ~Yt & ...
                cat_vol_morph(Ya1~=LAB.TH & Ya1~=LAB.BG & Ya1~=LAB.HC,'d',1);
      if ~debug, clear Yt Yi; end
      Yflairl = cat_vol_morph(Yflairl,'l',[inf (4 - WMHCstr)^3])>0;
      
      % add FLAIR mask
      Ywmh(Yflairl & Ywmh<2) = 1; 
    end
    
    % region-growing / bottleneck
    %Ywmh = cat_vol_downcut(Ywmh, 1+(3-Ym)/3, WMHCstr/2, vx_vol); 
    Ywmh(Ywmh==0) = 1.5; 
    Ywmh = cat_vol_laplace3R(Ywmh, Ywmh==1.5, 0.001); if debug, Ywmh2 = Ywmh; end
    %
    Ywmh = ...YwmhA>=(1 - WMHCstr) & 
      Yp0A>2 & ~Yvt & Ym>mth(1) & Ym<mth(2) & (Ywmh <= 1.5);
    % remove small WMHs
    Ywmh = cat_vol_morph(Ywmh==1, 'l', [inf wmhvols(1)])>0;
    Ywmh = cat_vol_morph(Ywmh   , 'l', [inf wmhvols(2)])>0;

    %% apply to atlas
    Ya1(Ywmh) = LAB.HI;
  end
  if ~debug, clear Ywmh Yvt Ynwmh Yvt2; end

  
  
  %% Closing of gaps between diffent structures:
  stime = cat_io_cmd('  Closing of deep structures','g5','',verb,stime); dispc=dispc+1;
  Yvtd2 = cat_vol_morph(Ya1==LAB.VT,'d',2,vx_vol) & Ya1~=LAB.VT;
  % CT and VT
  Yt = cat_vol_morph(Ya1==LAB.VT,'d',2,vx_vol) & ...
       cat_vol_morph(Ya1==LAB.CT,'d',2,vx_vol) & Ya1==0 ;
  Ya1(Yt & Yp0<=1.5 & ~Ynv)=LAB.VT; Ya1(Yt & Yp0>1.5)=LAB.CT; 
  % WMH and VT
  Yt = cat_vol_morph(Ya1==LAB.HI,'d',1,vx_vol) & Yvtd2 & ~Ynv & Ya1==0;
  Ya1(Yt &  Ym<=1.25)=LAB.VT; Ya1(Yt & Ym>1.25 & Ym<2.5)=LAB.HI; 
  % TH and VT
  Yt = cat_vol_morph(Ya1==LAB.TH,'d',1,vx_vol) & Yvtd2;
  Ya1(Yt & Ym<=1.5)=LAB.VT; Ya1(Yt & Ym>1.5 & Ym<2.85)=LAB.TH; 
  % BG and VT
  Yt = cat_vol_morph(Ya1==LAB.BG,'d',1,vx_vol) & Yvtd2;
  Ya1(Yt & Ym<=1.5)=LAB.VT; Ya1(Yt & Ym>1.5 & Ym<2.85)=LAB.BG;
  % no bloodvessels next to the ventricle, because for strong atrophy
  % brains the WM structures can be very thin and may still include 
  % strong bias
  Ya1(Ya1==LAB.BV & cat_vol_morph(Ya1==LAB.VT,'d',3,vx_vol))=0;
  if ~debug, clear Yt Yh Yvtd2 Yw; end
 
  
  
  %% complete map
  [tmp0,tmp1,Ya1] = cat_vbdist(Ya1,Yb); clear tmp0 tmp1;
  
  % consider gyrus parahippocampalis
  Ya1(YA==LAB.PH) = LAB.PH;
  
  %% side aligment using laplace to correct for missalignments due to the normalization
  stime = cat_io_cmd('  Side alignment','g5','',verb,stime); dispc=dispc+1;
  YBG  = Ya1==LAB.BG | Ya1==LAB.TH;
  YMF  = Ya1==LAB.VT | Ya1==LAB.BG | Ya1==LAB.TH | Ya1==LAB.HI; 
  YMF2 = cat_vol_morph(YMF,'d',2,vx_vol) | Ya1==LAB.CB | Ya1==LAB.BS | Ya1==LAB.MB;
  Ymf  = max(Ym,smooth3(single(YMF2*3))); 
  Yt = cat_vol_smooth3X(YS==0,6)<0.9 & cat_vol_smooth3X(YS==1,6)<0.9 & ~YMF2 & Yp0>0 & Ym<3.1 & (Yp0<2.5 | Ya1==LAB.BV);
  Ys = (2-single(YS)) .* single(smooth3(Yt)<0.4);
  Ys(Ys==0 & (Ym<1 | Ym>3.1))=nan; Ys = cat_vol_downcut(Ys,Ymf,0.1,vx_vol); 
  [tmp0,tmp1,Ys] = cat_vbdist(Ys,Ys==0); clear tmp0 tmp1;
  if ~debug, clear YMF2 Yt YS; end
  
  %% YMF for FreeSurfer fsaverage
  Ysm  = cat_vol_morph(Ys==2,'d',1.75,vx_vol) & cat_vol_morph(Ys==1,'d',1.75,vx_vol);
  YMF  = cat_vol_morph(Ya1==LAB.VT | Ya1==LAB.BG | Ya1==LAB.HI | (Ya1==LAB.TH & smooth3(Yp0)>2),'c',3,vx_vol) & ~Ysm; 
  %YMF  = YMF | (cat_vol_morph(YA==LAB.CT & YBG,'c',6) & ~Ysm); 
  YMF  = Ym<=2.5  & cat_vol_morph(YMF | Ym>2.3,'c',1) & cat_vol_morph(YMF,'d',2,vx_vol);
  YMF  = smooth3(YMF)>0.5;
  clear Ysm; 
  
  
  %% back to original size
  stime = cat_io_cmd('  Final corrections','g5','',verb,stime); dispc=dispc+1;
  Ya1 = cat_vol_resize(Ya1,'dereduceV',resTr,'nearest'); Ya1 = cat_vol_median3c(Ya1,Ya1>1 & Ya1~=LAB.BV);
  Ys  = cat_vol_resize(Ys ,'dereduceV',resTr,'nearest'); Ys  = 1 + single(smooth3(Ys)>1.5);
  YMF = cat_vol_resize(YMF,'dereduceV',resTr);
  YBG = cat_vol_resize(YBG,'dereduceV',resTr);
  
  Ya1 = cat_vol_resize(Ya1,'dereduceBrain',BB); Ya1 = cat_vol_ctype(Ya1);
  Ys  = cat_vol_resize(Ys ,'dereduceBrain',BB); [tmp0,tmp1,Ys] = cat_vbdist(Ys,Ya1>0); clear tmp0 tmp1;
  YMF = cat_vol_resize(YMF,'dereduceBrain',BB); 
  YBG = cat_vol_resize(YBG,'dereduceBrain',BB); 
  Ym  = Ym0; clear Ym0;

  % final side alignment
  Ya1(Ya1>0)=Ya1(Ya1>0)+(Ys(Ya1>0)-1);
 
  
  % class correction
  % YBG is smoothed a little bit and (B) reset all values that are related
  % with GM/WM intensity (Ym<2.9/3) (A)
  Yclssum = single(Ycls{1})+single(Ycls{2})+single(Ycls{3});
  YBGs    = min( max(0,min(255, 255 - cat_vol_smooth3X(Ya1==1 & Ycls{2}>round(2.9/3),0.8) .* single(Ycls{2}) )), ... (A)
                 max(0,min(255, 255 * cat_vol_smooth3X(YBG .* (Ym<=2.9/3 & Ym>2/3) ,0.5) )) ); % (B)
  Ycls{1} = cat_vol_ctype(single(Ycls{1}) + YBGs .* (single(Ycls{2})./max(eps,Yclssum)));
  Ycls{2} = cat_vol_ctype(single(Ycls{2}) - YBGs .* (single(Ycls{2})./max(eps,Yclssum)));
  clear YBGs Yclssum; 
 
  if debug
    cat_io_cmd(' ','','',verb,stime); 
  else
    cat_io_cmd(' ','','',verb,stime); 
    %cat_io_cmd('cleanup',dispc,'',verb); 
  end
  
end
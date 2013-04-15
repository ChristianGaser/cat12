function [Vp4T,VpfT,VmfT,VepT,time]=vbm_vol_partvol(VmxT,Vp0T,opt)
% __________________________________________________________________________________________________
% Use a segment map p0T, the T1 data (with an intensity scale given by the p0T map) and the atlas 
% label map p4T to create a individual label map p4T. Furthermore, a map with filled ventricle and
% subcortical GM structures is generated. The atlas contain main regions like cerebrum, brainstem,
% midbrain, cerebellum and ventricle. Furthermore, it try to detect blood vessels by a combination
% of (atlas, thickness,) intensity, and distance information. 
%
% This function try to solve the following problems:
%  1) Finding of the cerebrum, the cerebellum, the head, blood vessels, brain skin and other mayor
%     structures based on atlas (p4T) and tissue class information (p0T). 
%     To do this it is important to use data from the T1-map (mxT) that use the same intensity 
%     scalling as the segment map p0T. 
%  2) Set Partions:
%     2.1) Find biggest WM part of each region... works
%     2.2) Align the nearest region class for other voxel
%     2.3) Finding and Filling of the ventricle and the Basalganglia
%     3.1) Removing of blood vessels and brain skin and the optical nerv
%          (set to 1.45 near the GM, else 1.00)
%     3.2) Removing of head structures (set to 0.00)
%
% _________________________________________________________________________
% ToDo: - use of interpolated data: actual segment maps has to have the original resolution, 
%         because of the downcut functions allow to get to higher values (based on the tissue class) 
%         => update of downcut functions
% _________________________________________________________________________
%
%   Was ist neu im Vergleich zu anderen?
%   - Zuweisung durch Dartel mit hoher Genauigkeit
%   - Erweiterung von SPM/VBM durch MainROIs (Seiten, Lappen, ...)
%   - Verbesserung der SPM/VBM durch bessere Enfernung von unerwünschtem
%     Gewebe (ON, Blutgefäße ...)
%   - Blutgefäße könnnen als erweitere Masken für fMRI genutzt werden um
%     Seiteneffekte besser ausblenden zu können.
%   - Beliebige Atlanten können genutzt werden.
%
% _________________________________________________________________________
%
% Structure:
%
%   [VpfT,Vp4T,VmfT,time]=vbm_vol_partvol(VmxT,Vp0T,opt)
%
%   INPUT:  p4T   = 3D-volume with brain regions (altas map)
%           p0T    = 3D-volume with tissue propability map (CSF=1,GM=2;WM=3)
%           opt
%            .resV = Voxelsize
%            .LAB  = Label of p4T map (see def.LAB definition below)
%            
%
%   OUTPUT: PA   = individual label map 
%           p0PF  = Corrected p0T map based on opt.PF
%                     without optical nerv   (opt.PF = {...,'-ON')
%                     without bloodvessels   (opt.PF = {...,'-BV')
%                     filled ventricle       (opt.PF = {...,'+VmxT')
%                     filled Basalganglias   (opt.PF = {...,'+BG')
%                     filled Thalamus        (opt.PF = {...,'+TH')
%
% _________________________________________________________________________
% Center of Neuroimaging, University Jena, Germany
% Robert Dahnke, Christian Gaser
% 2010/09

  if nargin<3, error('MATLAB:vbm_vol_partvol','ERROR: Need the intensity correctd T1, the p0T segmenation and the template atlas map!\n'); end

  if ~exist('opt','var'), opt=struct(); end
  % Templates
  def.Fp0A            = fullfile(spm('Dir'),'toolbox','vbm12+','templates_1.50mm','p0A.nii');
  def.Fp4A            = fullfile(spm('Dir'),'toolbox','vbm12+','templates_1.50mm','p4A.nii');
  def.res             = 2; 
  % Normalization
  def.norm.smoref     = 0;     % no smoothing for template image 
  def.norm.smosrc     = 8;     % similar smoothing like the template image
  def.norm.regtype    = 'mni'; 
  def.norm.cutoff     = 15;    % need high resolution
  def.norm.nits       = 3;     % need only few interations
  def.norm.reg        = 0.1;
  def.write.interp    = 0;     % nearest neighbor 
  def.write.preserve  = 0;     % no modulation
  def.write.prefix    = '';    % no prefix
  def.write.vox       = [1 1 1]; % set resolution later 
  def.write.bb        = [[-inf -inf -inf];[inf inf inf]]; % no bb
  % definition of ROIs with: ID = [L R]
  def.LAB.CT = [ 1  2]; % cortex
  def.LAB.MB = [13 14]; % MidBrain
  def.LAB.BS = [13 14]; % BrainStem
  def.LAB.CB = [ 3  4]; % Cerebellum
  def.LAB.ON = [11 12]; % Optical Nerv
  def.LAB.BG = [ 5  6]; % BasalGanglia 
  def.LAB.TH = [ 9 10]; % Hypothalamus 
  def.LAB.HC = [19 20]; % Hippocampus 
  def.LAB.VmxT = [15 16]; % Ventricle
  def.LAB.NV = [17 18]; % no Ventricle
  def.LAB.BV = [ 7  8]; % BlodsodVessels
  def.LAB.NB = [ 0  0]; % no brain 
  def.LAB.HD = [21 22]; % head
  % output files
  def.FpfTpre       = 'pf';
  def.FmfTpre       = 'mf';
  def.Fp4Tpre       = 'l1';
  def.FepTpre       = 'ep';
  opt               = checkinopt(opt,def);
  opt.FpfT          = vbm_io_handle_pre(VmxT.fname,opt.FpfTpre);
  opt.FmfT          = vbm_io_handle_pre(VmxT.fname,opt.FmfTpre);
  opt.Fp4T          = vbm_io_handle_pre(VmxT.fname,opt.Fp4Tpre);
  opt.FepT          = vbm_io_handle_pre(VmxT.fname,opt.FepTpre);
  opt.mat           = vbm_io_handle_pre(VmxT.fname,''); opt.mat = [opt.mat(1:end-3) 'mat'];
  
  V.mxT = VmxT; clear VmxT; 
  V.p0T = Vp0T; clear Vp0T;
  
  stime = clock;
    
% check if p4T is not only a copy of the p4T file (maybe because the previos partioning was interrupted).
  if any([~exist(opt.FpfT,'file'),~exist(opt.FmfT,'file'),~exist(opt.Fp4T,'file'),~exist(opt.FepT,'file')]) || opt.recalc==1
    
    % map template label map to individual space
    spm_normalise(V.p0T,spm_vol(opt.Fp0A),opt.mat,'','',opt.norm); 
    fprintf(sprintf('%s',repmat('\b',1,85 + 33 + 43*opt.norm.nits))); 
    opt.write.vox = sqrt(sum(V.p0T.mat(1:3,1:3).^2));
    copyfile(opt.Fp4A,opt.Fp4T);
    spm_write_sn(opt.Fp4T,opt.mat,opt.write);
    V.p4T = V.mxT; V.p4T.fname = opt.Fp4T; 
    spm_imcalc(spm_vol(char({V.mxT.fname,V.p4T.fname})),V.p4T,'i2');
    
    % first load the images if files are given
    p0T    = single(spm_read_vols(V.p0T));
    mxT    = single(spm_read_vols(V.mxT)); 
    p4A    = single(round(spm_read_vols(V.p4T))); [D,D,p4A]=vbdist(p4A); clear D; p4A=uint8(p4A);  
    vx_vol = sqrt(sum(V.p0T.mat(1:3,1:3).^2)); 
    
    % OPTIMIZATION:
    % ds('l2','',vx_vol,mxT,p4A,mxT,p0T/3,130)
    [mxT,p0T,p4A,BB] = vbm_vol_resize({mxT,p0T,p4A},'reduceBrain',vx_vol,5,p0T>0);   % removing of background
    [mxT,p0T,resTr]  = vbm_vol_resize({mxT,p0T},'reduceV',vx_vol,opt.res,64); 
    p4A              = vbm_vol_resize(p4A   ,'reduceV',vx_vol,opt.res,64,'nearest'); 
    
    p4ANS = round(p4A/2)*2-1;
    mi3T  = mxT*3; %vbm_vol_iscale(mxT,'gCGW',resTr.vx_volr,p0T)*3;
    
    % alignment of high intensity structures 
    p4T=zeros(size(mxT),'single'); p4T(vbm_vol_morph(~p0T & mi3T<0.3,'ldo',2,resTr.vx_volr)==1)=-1; p4T(p4A)=0;       % backgound
    p4T(vbm_vol_morph(~p0T,'open')==1 & mi3T>2.5) = opt.LAB.HD(1);                                                    % head
    for s=1:2
      p4T(vbm_vol_morph(mi3T>2.5 & p4A==opt.LAB.CT(s),'lab')) = opt.LAB.CT(1); % cortex
      p4T(vbm_vol_morph(mi3T>2.5 & p4A==opt.LAB.BS(s),'lab')) = opt.LAB.BS(1); % brainstem
      p4T(vbm_vol_morph(mi3T>2.5 & p4A==opt.LAB.MB(s),'lab')) = opt.LAB.MB(1); % midbrain
      p4T(vbm_vol_morph(mi3T>2.5 & p4A==opt.LAB.CB(s),'lab')) = opt.LAB.CB(1); % cerebellum
      p4T(vbm_vol_morph(mi3T>2.5 & p4A==opt.LAB.ON(s),'lab')) = opt.LAB.ON(1); % optical nerv
    end
    p4T(p4ANS==opt.LAB.BV(1) & mi3T>3 & p0T>0 & vbm_vol_morph(mi3T<2,'labclose',1,resTr.vx_volr) & ...
      (~vbm_vol_morph(mi3T>2.5,'lab') | mi3T>4.0)) = opt.LAB.BV(1);                                 % blood vessels

    % region-growing for special high intensity regions
    p4T((mi3T<=2.9 & p4T==0)) = -inf; p4T = vbm_vol_simgrow(p4T,mi3T,1);p4T(isinf(p4T))=0; 
  %  M = mi3T>2.5 & (p4ANS==opt.LAB.CT(1) | p4ANS==opt.LAB.CB(1) | p4ANS==opt.LAB.BV(1) | p4ANS==opt.LAB.HD(1) | p4ANS==opt.LAB.ON(1));
  %  p4T(M)=p4ANS(M); p4T(isinf(p4T))=0; 

    % alignment of medium and low intensity structures
    p4T(p4ANS==opt.LAB.BG(1) & p4T==0 & mi3T>1.75 & mi3T<2.75 & mi3T<2.9) = opt.LAB.BG(1);          % basal ganglia
    p4T(p4ANS==opt.LAB.TH(1) & mi3T>1.75 & mi3T<2.75 & mi3T<2.9) = opt.LAB.TH(1);                   % hypothalamus
    p4T(p4ANS==opt.LAB.HC(1) & mi3T>1.75 & mi3T<2.75 & mi3T<2.9) = opt.LAB.HC(1);                   % hippocampus
    for s=1:2
      p4T(vbm_vol_morph(p4A==opt.LAB.VmxT(s) & mi3T<1.75,'lab') & mi3T<1.25) = opt.LAB.VmxT(1);     % ventricle
    end
      
    % refinement of WM 
    p4T(p4ANS==opt.LAB.CT(1) & p4T==0 & vbm_vol_morph(mi3T>2 & p4T==1,'labopen',1)) = opt.LAB.CT(1);  
    p4T(p4ANS==opt.LAB.BG(1) & ~vbm_vol_morph(p4ANS==opt.LAB.BG(1),'erode',1) & p4T==opt.LAB.BG(1) & (mi3T>2.75 & mi3T>2.5)) = opt.LAB.CT(1);                               % basal ganglia
    p4T=vbm_vol_median3(p4T,p4T>=0,p4T>0);  

    % regino-growing for all high intensity regions
    p4T((mi3T<=2.25 & p4T==0))=-inf; p4T = vbm_vol_simgrow(p4T,mi3T,0.5); p4T(isinf(p4T))=0; 
    p4T(p4T==0 & vbm_vol_morph(mi3T>2.75 & p0T,'labopen',1)==0 & mi3T>2.75 & p4T==0 & vbm_vol_morph(mi3T<2,'labclose',1,resTr.vx_volr))=opt.LAB.BV(1);

    % region-growing in GM only for non-blood vessels regions
    p4T(p4T==opt.LAB.BV(1) | isinf(p4T))=0; p4T((mi3T<=1 | ~p0T) & p4T==0) = -inf;  [p4T,D]=vbm_vol_downcut(p4T,mi3T,0.1,resTr.vx_volr); 
    p4T((p4ANS==opt.LAB.BV(1) | D>500) & mi3T>1.5 & ~vbm_vol_morph(mi3T>2.5,'lab') & ...
      (p4ANS~=opt.LAB.CB(1) | p4ANS~=opt.LAB.HD(1)))=opt.LAB.BV(1);                     % adding blood vessels again
    p4T(p4T==opt.LAB.BG(1) & ~(p4ANS==opt.LAB.BG(1)|p4ANS==opt.LAB.VmxT(1)))=opt.LAB.CT(1); % basal ganglia (avoid overgrowing)
    p4T(p4T==0 & mi3T>3.5)=opt.LAB.BV(1); p4T=vbm_vol_downcut(p4T,mi3T,0.2,resTr.vx_volr);
    p4T(isinf(p4T))=0; p4T=vbm_vol_median3(p4T,p4T>=0 & mi3T<2.5 & ~(p4T==opt.LAB.BV(1) & mi3T>2.25),p4T>=0); p4T(p4T==0 & ~p0T & mi3T<0.75)=-inf;

    % filling of subcortical regions
    M3 = p4ANS==opt.LAB.BG(1)|p4ANS==opt.LAB.VmxT(1)|p4ANS==opt.LAB.TH(1);
    M2 = vbm_vol_morph(M3,'d',2,resTr.vx_volr);
    M  = 3*vbm_vol_smooth3X(single(vbm_vol_morph(vbm_vol_morph(vbm_vol_morph(M3 | (M2 & mi3T>2.7 & p4T~=7 & p4T~=8),'lc',1),'lo',1),'e',0)),0.5);  
    M2 = vbm_vol_morph(M3,'dd',4,resTr.vx_volr);
    MF = M2.*max(mi3T,M); clear M2 M;  
    mf3T = max(mi3T,MF); TI3FS = vbm_vol_smooth3X(mf3T,1.5); mf3T(mf3T>3 & MF>3)=TI3FS(mf3T>3 & MF>3);
    
    % removement of bv & meninges ... NEED MORE WORK
   % BV = vbm_vol_morph((p4ANS==opt.LAB.HD(1) | p4ANS==opt.LAB.BV(1) | p4T<=0) & ~p0T,'labclose')==1; D=vbdist(single(~BV))*2;
   % pfT(BV)=max(1.9-D(BV),1); mf3T(BV)=max(min(1.9-D(BV),mf3T(BV)),1);

    % side aligment using laplace to correct for missalignments due to the normalization
    d = 5; M = vbm_vol_smooth3X(single(vbm_vol_morph(mod(p4A,2)==0,'dd',d,resTr.vx_volr)) & single(vbm_vol_morph(mod(p4A,2)==1,'dd',d,resTr.vx_volr)==1),20);
    S = 2*single(mod(p4A,2)==0 & mi3T>2.5 & M<max(M(:))*0.9) +single(mod(p4A,2)==1 & mi3T>2.5 &  M<max(M(:))*0.9); S(mf3T<=2.5)=-inf; S(S==0)=1.5;

    rS=vbm_vol_resize(S,'reduce'); rS=round(rS*2)/2; rS=vbm_vol_laplace3R(rS,rS==1.5,0.001); rS=vbm_vol_resize(single(rS),'dereduce',size(mi3T)); S(S==1.5)=round(rS(S==1.5)); 
    S(isinf(S) & mf3T>2)=0; S=vbm_vol_downcut(S,mf3T,3,resTr.vx_volr); S(isinf(S) & mf3T>0)=0; S=vbm_vol_downcut(S,mf3T,3,resTr.vx_volr); S(S<=0)=2-mod(p4A(S<=0),2);
    S=round(vbm_vol_smooth3X(S,2));
    p4T(p4T>0)=p4T(p4T>0)+(S(p4T>0)==2); 

    MF       = vbm_vol_resize(MF,'dereduceV',resTr);
    p4T      = vbm_vol_resize(p4T,'dereduceV',resTr,'nearest');
    [MF,p4T] = vbm_vol_resize({MF,p4T},'dereduceBrain',BB);
    
    p0T = single(spm_read_vols(V.p0T)); pfT = max(p0T,min(3,MF)); 
    mfT = single(spm_read_vols(V.mxT)); mfT = vbm_vol_iscale(mfT,'gCGW',resTr.vx_volr,p0T);
    mfT = max(mfT,min(1,MF/3)); smfT = vbm_vol_smooth3X(mfT,1.5); mfT(mfT>1 & MF>3)=smfT(mfT>1 & MF>3);
    
    VpfT = vbm_io_write_nii(pfT,V.mxT,'pf','filled','float32',[0,1]);
    VmfT = vbm_io_write_nii(mfT,V.mxT,'mf','filled','float32',[0,1]);
    Vp4T = vbm_io_write_nii(p4T,V.mxT,'p4','atlas' ,'uint8'  ,[0,1]);
    VepT = [];
    
    time = dp('|P',opt,stime);   
  else 
    VpfT = spm_vol(opt.FpfT);
    Vp4T = spm_vol(opt.Fp4T);
    VmfT = spm_vol(opt.FmfT);
    VepT = [];
    
    time = 0;
  end
end

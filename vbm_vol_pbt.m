function [GMT,PP]=vbm_vol_pbt(SSEG,opt)
% ______________________________________________________________________
%
% Use a graph-based thickness estimation (gbdist) to estimate the distance
% from the WM/GM boundary and a projection scheme to transfer the values at
% GM/CSF boundary over the whole GM.
%
%   [GMT,PP]=vbm_vol_pbt(SSEG,opt)
%  
%   GMT:      GM thickness map 
%   PP:       percentage possition map
%   SSEG:     tissue segment image
%   opt.resV  voxel resolution (only isotropic)
%
%   Needed functions:
%   - qedist.cpp (fast estimation of a submask for speedup - only isotropic)
%   - gbdist.m   (only isotropic)
%   - WMDskeleton.m
%
%   See also BWT, VBT, GBT, GBTX, GBCL, LBT, LBTS, PBTV.
%   
% ______________________________________________________________________
%
% More details and test:
%   ... Paper
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of psychiatry and psychotherapy 
%   University Jena
%
%   Version: 1.00 © 2010/08
% ______________________________________________________________________
% $Id$ 

% default variables and check/set function  
  if ~exist('opt','var'), opt=struct(); end
  def.resV   = 1;
  def.fsize  = @(x) 1+2*max(round(x * 1.5),1); 
  def.method = 'pbt2x';
  %cond = {'resV(1)==resV(2)';'resV(1)==resV(3)'}; % only isotropic (see qedist and gbdist)
  opt  = checkinopt(opt,def);%,cond);
  mvxs = mean(opt.resV);
  SSEG = round(SSEG*10)/10;
  
  % additional re-estimation of the boundarys... in development:
  % funkt nicht ... der einschnitt erzeugt einen fehler der vergleichbar
  % hoch ist - zwar etwas netter (vor allem dünner) - aber nicht sinnvoll
  % begründbar
  if  0 && strcmp(opt.method,'pbt2x')
    % Refinement of sulcal areas: Although pbt do not need this to find
    % the sulcus it helps to find the CSF boundary a little bit more exact.
    % It it not usefull for the WM boundary because it may remove the
    % sulcus!
    
%     M = max(0,min(1,(SSEG-2.9)*10)); M(SSEG<=2)=-inf; WMD  = eidist_r04(M,ones(size(M),'single'));
%     M = max(0,min(1,(2.1-SSEG)*10)); M(SSEG>=3)=-inf; CSFD = eidist_r04(M,ones(size(M),'single'));
%     [~,PP] = vbm_vol_pbtp( (1.1*(2.5-SSEG)) + 2 , CSFD , WMD);  % 1.25 , 1.1
%     PP = 1-vbm_vol_median3(1-PP,0.5,0,0.5);
%     M = SSEG>2 & SSEG<3; SSEG(M) = max(SSEG(M),2.75-0.75*PP(M)); %max(SSEG(M),min(2.75,3-PP(M)));
%     clear WMD CSFD PP;

    M = max(0,min(1,(SSEG-1.9)*10)); M(SSEG<=1)=-inf; WMD  = vbm_vol_eidist(M,ones(size(M),'single')); 
    M = max(0,min(1,(1.1-SSEG)*10)); M(SSEG>=2)=-inf; CSFD = vbm_vol_eidist(M,ones(size(M),'single')); 

    try %#ok<TRYNC>
      [GMT,PP] = vbm_vol_pbtp( (1.1*(SSEG-1.5) ) + 2 , WMD , CSFD );  % 1.25 , 1.1
      PP = vbm_vol_median3(PP,0.5,0,0.5); PP=1-PP; PP=PP.*(GMT*0.375); 
      M = SSEG>1 & SSEG<2; SSEG(M) = min(SSEG(M),max(1,2-PP(M))); %min(SSEG(M),max(1.25,1+PP)); %min(SSEG(M),1.45+PP(M).*(0.55*(SSEG(M)-1)));
    end
    clear PP SSEGO WMD CSFD;
  end
  
  

  % estimate WM distance WMD and the non-corrected CSF distance CSFD (not correct in sulcal areas)  
  %L = laplace3(single(SSEG==2)*0.5 + single(SSEG>2),0,1,0.01);
  M = max(0,min(1,(SSEG-2))); M(SSEG<1)=-inf; WMD  = vbm_vol_eidist(M,max(0,( max(0,min(1,SSEG-1))))); % L + min.../2
  M = max(0,min(1,(SSEG-1))); M(SSEG<1)=-inf; CSFD = vbm_vol_eidist(M,max(0,  max(0,min(1,SSEG-1))));   
  M = SSEG>1 & SSEG<1.5; WMD(M) = WMD(M) - CSFD(M); clear CSFD;
  WMDM = vbm_vol_median3(WMD,WMD>mvxs & M,WMD>mvxs & M,opt.resV/8); WMD(M) = WMDM(M); clear WMDM;

  if strcmp(opt.method,'pbt2x')  
    M = max(0,min(1,(2-SSEG))); M(SSEG>=3)=-inf; CSFD = vbm_vol_eidist(M,max(0,(min(1,3-SSEG)))); %(1-L + min(1,3-SSEG))/2))
    M = max(0,min(1,(3-SSEG))); M(SSEG>=3)=-inf; WMDC = vbm_vol_eidist(M,max(0, min(1,3-SSEG)));   
    M = SSEG>2.5 & SSEG<3.0; CSFD(M) = CSFD(M) - WMDC(M); clear WMDC;
    CSFDM = vbm_vol_median3(CSFD,CSFD>mvxs & M,CSFD>mvxs & M,opt.resV/8); CSFD(M) = CSFDM(M); clear CSFDM;
    clear L M; 

    WMD(isnan(WMD) | isinf(WMD) | isinf(-WMD))=0; CSFD(isnan(CSFD) | isinf(CSFD) | isinf(-CSFD))=0; % das nicht nötig sein...
    
    GMT1 = vbm_vol_pbtp(  SSEG,WMD,CSFD); %########################## TIM
    try GMT2 = vbm_vol_pbtp(4-SSEG,CSFD,WMD); catch, GMT2=inf(size(SSEG),'single'); end %#ok<CTCH>
    GMT2(GMT2<=(2/mean(opt.resV)) | SSEG<2)=inf; % try this correction only in thick regions
    
    [GMT,I] = min(cat(4,GMT1,GMT2+0.25*mean(opt.resV)),[],4); 
    PP=zeros(size(SSEG),'single');
    M=SSEG>=1.5 & SSEG<=2.0; PP(M) = ((GMT1(M) - WMD(M))) ./ (GMT1(M)); 
    M=SSEG> 2.0 & SSEG<=2.5; PP(M) = ((GMT1(M) - WMD(M)).*(I(M)==1) + CSFD(M).*(I(M)==2)) ./ (GMT(M) + eps); 
    PP(SSEG>2.5)=1;
  else
    M = max(0,min(1,(2-SSEG))); M(SSEG>=3)=-inf; CSFD = vbm_vol_eidist(M,max(0,(min(1,3-SSEG))/2)); % (1-L + min(1,3-SSEG))/2))

    clear L M; %CSFD(DG>(sqrt(3)/opt.resV) | SSEG<1.5) = 0; 
    WMD(isnan(WMD) | isinf(WMD) | isinf(-WMD))=0; CSFD(isnan(CSFD) | isinf(CSFD) | isinf(-CSFD))=0; % das nicht nötig sein...

   [GMT,PP] = vbm_vol_pbtp(SSEG,WMD,CSFD);%########################## TIM
  end
  clear WMD CSFD;

  if opt.resV<1
    PP  = vbm_vol_median3(PP,PP>0.125 & PP<0.175,PP>0 & PP<1,0.25);
    PP  = vbm_vol_median3(PP,PP>0.125 & PP<0.175,PP>0 & PP<1,0.00); % 0.10
    GMT = vbm_vol_median3(GMT,GMT>eps,GMT>eps,opt.resV/4);
    GMT = vbm_vol_median3(GMT,GMT>eps,GMT>eps,0);
  end
  
  GMT = GMT*mean(opt.resV); % darf ich erst hier machen wegen der projektion die aktuell nur für 1mm läuft!
  GMT(SSEG<1.5 | SSEG>2.5)=nan; GMT=vbm_vol_nanmean3(GMT); GMT(isnan(GMT))=eps; % erweitern
  %GMTS = smooth3(GMT,'gaussian',3,0.9); GMT(SSEG>=2 & SSEG<2.5)=GMTS(SSEG>=2 & SSEG<2.5); clear GMTs % smoothen allerdings nur für sicheren bereich

  GMT(GMT>10)=0; PP(isnan(PP))=0; 
  YM  = vbm_vol_morph(GMT>0,'d');
  GMT = vbm_vol_median3(GMT,YM);
  PP  = vbm_vol_median3(PP ,YM);
  clear Ymt;

  % removing non brain objects
  Ymm = PP>=0.5 & ~vbm_vol_morph(PP>=0.5,'labopen',1); 
  PP(Ymm) = 0.5-eps; GMT(Ymm) = 0;
 
end
function [Ygmt,Ypp] = cat_vol_pbt(Ymf,opt)
% ______________________________________________________________________
%
% Cortical thickness and surface position estimation. 
% 
%   [Ygmt,Ypp]=cat_vol_pbt(Ymf,opt)
%  
%   Ygmt:      GM thickness map 
%   Ypp:       percentage position map
%   Ymf:       tissue segment image or better the noise, bias, and 
%              intensity corrected 
%   opt.resV   voxel resolution (only isotropic)
%   opt.method choose of method {'pbt2x','pbt2'} with default=pbt2x as 
%              the method that is described in the paper.
% ______________________________________________________________________
%
%   Dahnke, R; Yotter R; Gaser C.
%   Cortical thickness and central surface estimation.
%   NeuroImage 65 (2013) 226-248.
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
%
%   Version: 1.10 � 2014/02
% ______________________________________________________________________
% $Id$ 


% default variables and check/set function  
  if ~exist('opt','var'), opt=struct(); end
  def.resV   = 1;
  def.method = 'pbt2x';
  def.debug  = 0; %cat_get_defaults('extopts.debug');
  def.verb   = cat_get_defaults('extopts.verb')-1;
  opt  = cat_io_checkinopt(opt,def);
  opt.resV = mean(opt.resV);

   
  %% Estimate WM distance Ywmd and the outer CSF distance Ycsfdc to correct
  %  values in CSF area to limit the Ywmd to the maximum value that is 
  %  possible within the cortex.
  if opt.verb, fprintf('\n'); end
  stime = cat_io_cmd('    WM distance: ','g5','',opt.verb);
  YMM = ~(cat_vol_morph(Ymf>0.5,'d',1) | cat_vol_morph(Ymf>0.5,'d',1));
  YM = max(0,min(1,(Ymf-2))); YM(YMM)=nan; Ywmd   = cat_vol_eidist(YM,max(0,min(1,Ymf/2)),[1 1 1],1,1,0,opt.debug); 
  YM = max(0,min(1,(Ymf-1))); YM(YMM)=nan; Ycsfdc = cat_vol_eidist(YM,max(0,min(1,Ymf/2)),[1 1 1],1,1,0,opt.debug); 
  YM = cat_vol_morph(Ymf>2.5,'d'); Ywmd(YM) =  Ywmd(YM)/2;
  YM = cat_vol_morph(Ymf>1.5,'d'); Ycsfdc(YM) =  Ycsfdc(YM)/2;
  YM = Ymf>0 & Ymf<1.5; Ywmd(YM) = Ywmd(YM) - Ycsfdc(YM); Ywmd(isinf(Ywmd)) = 0; clear Ycsfdc;
  YwmdM = cat_vol_median3(Ywmd,Ywmd>opt.resV & YM,Ywmd>1,0.25); Ywmd(YM) = YwmdM(YM); clear YwmdM YM;

  
  %% Estimate CSF distance Ycsfd and the CSF distance Ycsfdc to correct 
  %  values in WM area to limit the Ycsfd to the maximum value that
  %  is possible within the cortex.
  stime = cat_io_cmd('    CSF distance: ','g5','',opt.verb,stime);
  Ywmm = cat_vol_morph(Ymf>2.5,'e'); % this was dilate???
  YM = max(0,min(1,(2-Ymf))); YM(Ywmm)=nan; Ycsfd = cat_vol_eidist(YM,max(0,min(1,(4-Ymf)/2)),[1 1 1],1,1,0,opt.debug); 
  YM = max(0,min(1,(3-Ymf))); YM(Ywmm)=nan; Ywmdc = cat_vol_eidist(YM,max(0,min(1,(4-Ymf)/2)),[1 1 1],1,1,0,opt.debug);   
  YM = cat_vol_morph(Ymf<1.5,'d'); Ycsfd(YM) =  Ycsfd(YM)/2;
  YM = cat_vol_morph(Ymf<2.5,'d'); Ywmdc(YM) =  Ywmdc(YM)/2;
  YM = Ymf>2.5 & ~Ywmm; Ycsfd(YM) = Ycsfd(YM) - Ywmdc(YM); Ycsfd(isinf(-Ycsfd)) = 0; clear Ywmdc;
  YcsfdM = cat_vol_median3(Ycsfd,Ycsfd>opt.resV & YM,Ycsfd>1,0.25); Ycsfd(YM) = YcsfdM(YM); clear YcsfdM YM;
  
  
  %% PBT is the default thickness estimation, but PBT2x is the optimized
  %  version that use both sulci and gyris refinements, because not only 
  %  thin sulci can blurred. PBT2x is furthermore the methode that is
  %  described in the paper, but PBT is faster.
  if strcmp(opt.method,'pbt2x')  
    % Estimation of the cortical thickness with sulcus (Ygmt1) and gyri 
    % correction (Ygmt2) to create the final thickness as the minimum map
    % of both.
    stime = cat_io_cmd('    PBT2x thickness: ','g5','',opt.verb,stime);
    %Ywmd  = cat_vol_median3(Ywmd,Ywmd>0,Ywmd>0);
    %Ycsfd = cat_vol_median3(Ycsfd,Ycsfd>0,Ycsfd>0);
    Ygmt1 = cat_vol_pbtp(Ymf,Ywmd,Ycsfd);   Ygmt1 = cat_vol_median3(Ygmt1,Ygmt1>0,Ygmt1>0,opt.resV/4); %,opt.resV/4
    Ygmt2 = cat_vol_pbtp(4-Ymf,Ycsfd,Ywmd); Ygmt2 = cat_vol_median3(Ygmt2,Ygmt2>0,Ygmt2>0,opt.resV/4); %,opt.resV/4
    Ygmt1(Ygmt1<=0.5 & Ygmt2>0) = Ygmt2(Ygmt1<=0.5 & Ygmt2>0);
    Ygmt2(Ygmt2<=0.5 & Ygmt1>0) = Ygmt1(Ygmt2<=0.5 & Ygmt1>0);
    [Ygmt,Yi] = min(cat(4,Ygmt1,Ygmt2+0.25*mean(opt.resV)),[],4); 
    Ygmt  = cat_vol_median3(Ygmt,Ygmt>eps,Ygmt>eps);    
    
    %% Estimation of a mixed percentual possion map Ypp.
    Ypp=zeros(size(Ymf),'single');
    YM=Ymf>1.5 & Ymf<=2;  Ypp(YM) =   (Ygmt1(YM) - Ywmd(YM)) ./ (Ygmt1(YM) + eps); 
    YM=Ymf>2.0 & Ymf<2.5; Ypp(YM) = ( (Ygmt1(YM) - Ywmd(YM)).*(Yi(YM)==1) + Ycsfd(YM).*(Yi(YM)==2) ) ./ (Ygmt(YM) + eps); 
    Ypp(Ymf>2.5 | (Ymf>2 & Ygmt<=1))=1;
    Ypp = cat_vol_median3(Ypp,Ymf>0 & Ygmt>2/opt.resV,Ymf>0 & Ygmt>2/opt.resV);
  else
    % Estimation of thickness map Ygmt and percentual position map Ypp.
    stime = cat_io_cmd('    PBT2 thickness: ','g5','',opt.verb,stime);
    [Ygmt,Ypp] = cat_vol_pbtp(Ymf,Ywmd,Ycsfd);
  end
  clear Ywmd Ycsfd;

  
  %% Final corrections for position map with removing of non brain objects.
  %   ds('l2','',1,Ymf/3,Ypp>0.5,Ypp,Ygmt/opt.resV,250)
  stime = cat_io_cmd('    Final Corrections: ','g5','',opt.verb,stime);
  Ypp(isnan(Ypp)) = 0; 
  Ypp = cat_vol_median3(Ypp,Ymf>0 & Ymf<3 & Ygmt>2/opt.resV,Ymf>0 & Ymf<3 & Ygmt>2/opt.resV,0.3);
  YM  = Ypp>=0.5 & ~cat_vol_morph(Ypp>=0.5,'labopen',1);  
  Ypp(YM) = 0.5-eps; clear YM
  
  
  % Final corrections for thickness map with thickness limit of 10 mm. 
  % Resolution correction of the thickness map after all other operations, 
  % because PBT actually works only with the voxel-distance (isotropic 1 mm)
  [tmp0,Yi] = cat_vbdist(single(Ygmt>eps),cat_vol_morph(Ygmt>eps | (Ymf>1.5 & Ymf<2.5),'d',2)); Ygmt=Ygmt(Yi); clear Yi tmp0; %#ok<ASGLU>
  Ygmt = cat_vol_median3(Ygmt,Ymf>0 & Ymf<3,Ygmt>eps,0.25);
  Ygmt = cat_vol_median3(Ygmt,Ymf>0 & Ymf<3,Ygmt>eps);
  Ygmt = Ygmt*opt.resV; 
  Ygmt(Ygmt>10) = 10; 
  
  
  if opt.debug, cat_io_cmd(' ','','',opt.debug,stime); end
end
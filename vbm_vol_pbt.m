function [Ygmt,Ypp] = vbm_vol_pbt(Ymf,opt)
% ______________________________________________________________________
%
% Cortical thickness and surface possition estimation. 
% 
%   [Ygmt,Ypp]=vbm_vol_pbt(Ymf,opt)
%  
%   Ygmt:      GM thickness map 
%   Ypp:       percentage possition map
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
%   Department of neurology
%   University Jena
%
%   Version: 1.10 © 2014/02
% ______________________________________________________________________
% $Id$ 


% default variables and check/set function  
  if ~exist('opt','var'), opt=struct(); end
  def.resV   = 1;
  def.method = 'pbt2x';
  def.debug  = 0; %cg_vbm_get_defaults('extopts.debug');
  def.verb   = cg_vbm_get_defaults('extopts.verb')-1;
  opt  = checkinopt(opt,def);
  opt.resV = mean(opt.resV);

  
  % rounding of values for simpler projection
  %Ymf  = round(Ymf*20)/20; 
  
   
  % Estimate WM distance Ywmd and the outer CSF distance Ycsfdc to correct
  % values in CSF area to limit the Ywmd to the maximum value that is 
  % possible within the cortex.
  if opt.verb, fprintf('\n'); end
  stime = vbm_io_cmd('    WM distance: ','g5','',opt.verb);
  YM = max(0,min(1,(Ymf-2))); YM(Ymf<=0)=nan; Ywmd   = vbm_vol_eidist(YM,max(0,min(1,Ymf/2)),[1 1 1],1,1,0,opt.debug); 
  YM = max(0,min(1,(Ymf-1))); YM(Ymf<=0)=nan; Ycsfdc = vbm_vol_eidist(YM,max(0,min(1,Ymf/2)),[1 1 1],1,1,0,opt.debug);   
  YM = Ymf>0 & Ymf<1.5; Ywmd(YM) = Ywmd(YM) - Ycsfdc(YM); Ywmd(isinf(Ywmd)) = 0; clear Ycsfdc;
  YwmdM = vbm_vol_median3(Ywmd,Ywmd>opt.resV & YM,Ywmd>1,0.25); Ywmd(YM) = YwmdM(YM); clear YwmdM YM;

  
  % Estimate CSF distance Ycsfd and the CSF distance Ycsfdc to correct 
  % values in WM area to limit the Ycsfd to the maximum value that
  % is possible within the cortex.
  stime = vbm_io_cmd('    CSF distance: ','g5','',opt.verb,stime);
  Ywmm = vbm_vol_morph(Ymf>2.5,'d');
  YM = max(0,min(1,(2-Ymf))); YM(Ywmm)=nan; Ycsfd = vbm_vol_eidist(YM,max(0,min(1,(4-Ymf)/2)),[1 1 1],1,1,0,opt.debug); 
  YM = max(0,min(1,(3-Ymf))); YM(Ywmm)=nan; Ywmdc = vbm_vol_eidist(YM,max(0,min(1,(4-Ymf)/2)),[1 1 1],1,1,0,opt.debug);   
  YM = Ymf>2.5 & ~Ywmm; Ycsfd(YM) = Ycsfd(YM) - Ywmdc(YM); Ycsfd(isinf(-Ycsfd)) = 0; clear Ywmdc;
  YcsfdM = vbm_vol_median3(Ycsfd,Ycsfd>opt.resV & YM,Ycsfd>1,0.25); Ycsfd(YM) = YcsfdM(YM); clear YcsfdM YM;
  
  
  % PBT is the default thickness estimation, but PBT2x is the optimized
  % version that use both sulci and gyris refinements, because not only 
  % thin sulci can blurred. PBT2x is furthermore the methode that is
  % described in the paper, but PBT is faster.
  if strcmp(opt.method,'pbt2x')  
    % Estimation of the cortical thickness with sulcus (Ygmt1) and gyri 
    % correction (Ygmt2) to create the final thickness as the minimum map
    % of both.
    stime = vbm_io_cmd('    PBT2x thickness: ','g5','',opt.verb,stime);
    Ygmt1 = vbm_vol_pbtp(Ymf,Ywmd,Ycsfd);   Ygmt1 = vbm_vol_median3(Ygmt1,Ygmt1>0,Ygmt1>0,opt.resV/4);
    Ygmt2 = vbm_vol_pbtp(4-Ymf,Ycsfd,Ywmd); Ygmt2 = vbm_vol_median3(Ygmt2,Ygmt2>0,Ygmt2>0,opt.resV/4);
    [Ygmt,Yi] = min(cat(4,Ygmt1,Ygmt2+0.25*mean(opt.resV)),[],4); 
        
    
    % Estimation of a mixed percentual possion map Ypp.
    Ypp=zeros(size(Ymf),'single');
    YM=Ymf>=1.5 & Ymf<=2.0; Ypp(YM) =   (Ygmt1(YM) - Ywmd(YM)) ./ (Ygmt1(YM) + eps); 
    YM=Ymf> 2.0 & Ymf<2.5;  Ypp(YM) = ( (Ygmt1(YM) - Ywmd(YM)).*(Yi(YM)==1) + Ycsfd(YM).*(Yi(YM)==2) ) ./ (Ygmt(YM) + eps); 
    Ypp(Ymf>=2.5 | (Ymf>2 & Ygmt<=1))=1;
  else
    % Estimation of thickness map Ygmt and percentual possion map Ypp.
    stime = vbm_io_cmd('    PBT2 thickness: ','g5','',opt.verb,stime);
    [Ygmt,Ypp] = vbm_vol_pbtp(Ymf,Ywmd,Ycsfd);
  end
  clear Ywmd Ycsfd;

  
  % Final corrections for position map with removing of non brain objects.
  stime = vbm_io_cmd('    Final Corrections: ','g5','',opt.verb,stime);
  Ypp(isnan(Ypp)) = 0; 
  Ypp = vbm_vol_median3(Ypp,Ymf>0 & Ymf<3,Ymf>0 & Ymf<3,2);
  YM  = Ypp>=0.5 & ~vbm_vol_morph(Ypp>=0.5,'labopen',1);  
  Ypp(YM) = 0.5-eps; clear YM
  
  
  % Final corrections for thickness map with thickness limit of 10 mm. 
  % Resolution correction of the thickness map after all other operations, 
  % because PBT actually works only with the voxel-distance (isotropic 1 mm)
  [tmp0,Yi] = vbdist(single(Ygmt>eps),vbm_vol_morph(Ygmt>eps | (Ymf>1.5 & Ymf<2.5),'d',2)); Ygmt=Ygmt(Yi); clear Yi tmp0;
  Ygmt = vbm_vol_median3(Ygmt,Ymf>0 & Ymf<3,Ygmt>eps,0.25);
  Ygmt = vbm_vol_median3(Ygmt,Ymf>0 & Ymf<3,Ygmt>eps);
  Ygmt = Ygmt*opt.resV; 
  Ygmt(Ygmt>10) = 10; 
  
  
  if opt.debug, vbm_io_cmd(' ','','',opt.debug,stime); end
end
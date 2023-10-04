function [Ygmt,Ypp] = cat_vol_pbtsimple(Yp0,vx_vol,method)
%cat_vol_pbtsimple. Simple cortical thickness/position estimation.  
% Voxel-based distance estimation and projection-based thickness (PBT) 
% and distance-based surface position estimation. Uses a label map as 
% input. Required a isotropic input map. 
%
%   [Ygmt,Ypp] = cat_vol_pbtsimple(Yp0[,vx_vol,method])
%
%   Ygmt    .. GM thickness map 
%   Ypp     .. percentage position map 
%   Yp0     .. tissue label map (1-CSF, 2-GM, 3-WM)
%   vx_vol  .. voxel-size (in mm, default = 1)
%   method  .. use voxel- or eikonal distance (0-voxel,1-eikonal)
%              voxel is more robust and faster but less accurate
%
%   See also cat_vol_pbt, cat_vol_createCS3.
% ______________________________________________________________________
%
%   Dahnke, R; Yotter R; Gaser C.
%   Cortical thickness and central surface estimation.
%   NeuroImage 65 (2013) 226-248.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$
  
  if ~exist('method','var'), method = 1; end
  



  % (X) Interpolation:
  % --------------------------------------------------------------------
  % Extra interpolation especially for the simple grid-based interpolation 
  % (that does not support partial volume effects) supports some  
  % improvements if the deinterpolation includes enough smoothing.
  % However, it is not clear if the increase in processing time is worth it.
  % For the simple vbdist (method=0) it seems to improves lower values that
  % where otherwise light underestimated. For the eidist I am not sure yet.
  %
  % This parameter is only for manual tests (especially vbdist)
  %   interp .. 1.0 - no interpolation, 0.5 - half resolution,
  %             0.8 - a bit interpolation (doulbed voxel size) 
  if method == 0, interp = 1.0; else; interp = 1; end % ########
  if interp ~= 1 
    V = struct('mat', eye(4), 'dim', size(Yp0));
    [Yp0,resI] = cat_vol_resize(max(1,min(3,Yp0)),'interp',V,interp,'cubic'); 
    vx_vol = vx_vol * interp;
  end




  % (1) Distance estimation:
  % --------------------------------------------------------------------
  if method == 0
  % Simple voxel-based function without considering partial volume effect 
  % (PVE) or asymetries but also without projection issues!
  % Due to it simpicity this is more robust.
    
    % CSF and WM distance maps (based on a binary boundary without PVE)
    Ycd = cat_vbdist(single(Yp0 < 1.5), Yp0 < 3); 
    Ywd = cat_vbdist(single(Yp0 > 2.5), Yp0 > 1); 

  else
  % More accurate Eikonal-based Euclidean distance estimation. 
  % But I am not really happy with the result, the projection seems to
  % include a lot of problems. 
  % Advantages should be vissible in the motorcortex or the insula.

    % CSF and WM speed maps to handle asymetries 
    Ycf  = max(eps, min(1, ((Yp0 - 1) / 1.1) .^4 )); 
    Ywf  = max(eps, min(1, ((4 - Yp0) / 2.0) .^2 )); 

    % CSF and WM boundary maps with PVE 
    Ycb  = max(0 , min(1,  2 - Yp0 )); Ycb(Yp0>2.5) = nan; 
    Ywb  = max(0 , min(1,  Yp0 - 2 )); Ywb(Yp0<1.5) = nan; 

    % CSF and WM distance maps with PVE
    Ywd  = cat_vol_eidist(Ywb, Ywf); 
    Ycd  = cat_vol_eidist(Ycb, Ycf);

    clear Ycf Ywf Ywb Ycb;
  end




  % (2) Thickness estimation:
  % --------------------------------------------------------------------
  if 1 % method == 0
  % Using the PBT apporach only to reconstruct the sulci.  

    % projection-based thickness mapping
    Ygmt = cat_vol_pbtp( round(Yp0) , Ywd, Ycd);

    % minimum to reduce issues with meninges
    Ygmt = min(Ygmt, Ycd + Ywd); 

  else
  % Using PBT to reconstruct the sulci and the gyri. The combination with 
  % the minimum is in priciple correct but the update of the CSF distance
  % is not optimal. As the later FS-based thickness estimation is used now
  % the benifits of this approach are maybe less usefull. 
  % Seems to be better for the eidist results 

    % reconstruct sulci as well as gyri 
    Ygmt1 = cat_vol_pbtp(round(Yp0)  , Ywd, Ycd);  
    Ygmt2 = cat_vol_pbtp(round(4-Yp0), Ycd, Ywd);

    % avoid meninges !
    Ygmt1 = min(Ygmt1, Ycd + Ywd);
    Ygmt2 = min(Ygmt2, Ycd + Ywd); 

    % average GM thickness maps
    Ygmt  = min(cat(4, Ygmt1, Ygmt2),[],4);

    clear Ygmt1 Ygmt2
  end




  % (3) Smoothing of thickness:
  % --------------------------------------------------------------------
  % To use a simple gaussian smoothing in an extend thickness map. 
  % A smoothing of 2 (1 mm) was some good compromise between details and 
  % less topology-correction artifacts. 
  % - cat_vol_approx: this function should fit a bit better but the other 
  %                   one is simpler and more direct 
  %Ygmt = cat_vol_approx(Ygmt,1); 
  % - need mask smoothing due to interpolation artifacts and to get 
  %   enlarge the filter area (saver but slower is to use the hole 
  %   volume by removing the last msk parameter)  
  Ygmt = simple_approx(Ygmt, 1, smooth3( Yp0>1 & Yp0<3 )>0 );  
   



  % (4) Estimate percentage position map:
  % --------------------------------------------------------------------
  % We first create a corrected CSF distance map with reconstructed sulci.
  % If gyri were reconstructed too than also the WMD would have to be
  % corrected to avoid underestimation of the position map with surfaces 
  % running to close to the WM.
  YM      = Yp0>=1.5 & Yp0<2.5 & Ygmt>eps;
  Ycdc    = Ycd; Ycdc(YM) = min(Ycd(YM), Ygmt(YM) - Ywd(YM)); 
  Ypp     = zeros(size(Yp0),'single'); Ypp(Yp0>=2.5) = 1;
  Ypp(YM) = Ycdc(YM) ./ (Ygmt(YM) + eps); 
  Ypp(Ypp>1) = 0;
  clear Ycdc YM;




  % Voxel-size resolution correction:
  % --------------------------------------------------------------------
  % Because Ycd and Ywd measure a grid-based distance (defined as the 
  % center of a voxel), we have to correct 1 voxel (2x.5 voxel, and finally 
  % correct for the _isotropic_ size of our voxel-grid. This is not
  % required for the Eikonal-based distance estimation! 
  Ygmt = (Ygmt - (method == 0) ) * mean(vx_vol); 




  if interp ~= 1 
    % back to original resolution
    Ygmt = cat_vol_resize(Ygmt,'deinterp',resI);                         
    Ypp  = cat_vol_resize(Ypp ,'deinterp',resI);                         
  end




  % final smoothing (for the folling downsampling) to reduce artifacts
  % --------------------------------------------------------------------
  spm_smooth(Ygmt , Ygmt , repmat(.5,1,3)); 
  spm_smooth(Ypp  , Ypp  , repmat(.5,1,3)); 

end
function Yo = simple_approx(Y,s,Ymsk)
%simple_approx. Simple approximation by the closest Euclidean value.
%
%  Y = simple_approx(Y[,s,Ymsk])
%  Y .. in/output image
%  s .. smoothing filter size
%

  if ~exist('s', 'var'), s = 1; end
  if ~exist('Ymsk', 'var'), Ymsk = true(size(Y)); end

  % estimate closest object point
  [~,I] = cat_vbdist(single(Y~=0),Ymsk > 0); 
  
  % align (masked) non-object voxels with closest object value
  Yo = Y(I);

  % smooth the result - correct for average to avoid smoothing boundary issues
  mnYo = median(Yo(Yo(:)~=0)); Yo = Yo - mnYo; Yo(~Ymsk) = 0; 
  spm_smooth(Yo , Yo , repmat(s,1,3)); 
  Yo = Yo + mnYo; Yo(~Ymsk) = 0; 
end

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
  
  
  if method == 0
  % Simple voxel-based function without considering partial volume effect 
  % (PVE) or asymetries but also without projection issues!

    % extra interpolation only for simple grid-based interpolation supports
    % light imoprovements although the deinterpolation is not working well
    % yet (would need the voxel binning)
    interp = 0.5;
    if interp > 0
      V = struct('mat', eye(4), 'dim', size(Yp0));
      [Yp0,resI] = cat_vol_resize(max(1,min(3,Yp0)),'interp',V,interp,'cubic'); 
      vx_vol = vx_vol * interp;
    end

    % CSF and WM distance maps (based on a binary boundary without PVE)
    Ycd = cat_vbdist(single(Yp0 < 1.5), Yp0 < 3); 
    Ywd = cat_vbdist(single(Yp0 > 2.5), Yp0 > 1); 

    % projection-based thickness mapping
    [Ygmt,Ypp] = cat_vol_pbtp(round(Yp0), Ywd, Ycd); 

    % Because Ycd and Ywd measure a grid-based distance (defined as the 
    % center of a voxel), we have to correct 1 voxel, and finally correct
    % for the _isotropic_ size of our voxel-grid.
    Ygmt = (Ygmt - 1) * mean(vx_vol); 

    if interp > 0
      % smoothing to integrate neighborhood to use "simple" interpolation
      spm_smooth(Ypp , Ypp , repmat(1.5,1,3)); %Ypp  = smooth3(Ypp);
      spm_smooth(Ygmt, Ygmt, repmat(1.5,1,3)); %Ygmt = smooth3(Ygmt);

      % back to original resolution
      Ygmt = cat_vol_resize(Ygmt,'deinterp',resI);                         
      Ypp  = cat_vol_resize(Ypp ,'deinterp',resI);                         
    end
  else
  % More accurate Eikonal-based Euclidean distance estimation. 
  % But I am not really happy with the result, the projection seems to
  % include a lot of problems. 
    
    % CSF and WM speed maps to handle asymetries 
    Ycf  = max(eps, min(1, ((Yp0 - 1) / 1.1) .^4 )); 
    Ywf  = max(eps, min(1, ((4 - Yp0) / 2.0) .^2 )); 

    % CSF and WM boundary maps with PVE
    Ycb  = max(0 , min(1,  2 - Yp0 ));
    Ywb  = max(0 , min(1,  Yp0 - 2 ));

    % CSF and WM distance maps with PVE
    Ywd  = cat_vol_eidist(Ywb, Ywf); 
    Ycd  = cat_vol_eidist(Ycb, Ycf);

    % as the values are further rising beyond the boundary we have to cut 
    % them of (in pbtx the values where corrected by another distance) 
    Ywd(Yp0<1.5) = 0; 
    Ycd(Yp0>2.5) = 0; 

    % projection-based thickness mapping 
    % - it is important to uses integer labels, e.g., round(Yp0) !
    [Ygmt,Ypp] = cat_vol_pbtp( round(Yp0) , Ywd, Ycd);

    % correction for voxel size
    Ygmt = Ygmt * mean(vx_vol); 
  end
  
 
end
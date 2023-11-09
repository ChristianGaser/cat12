function [Ygmt,Ypp] = cat_vol_pbtsimple(Yp0,vx_vol,method,interp,enhance)
%cat_vol_pbtsimple. Simple cortical thickness/position estimation.  
% Voxel-based distance estimation and projection-based thickness (PBT) 
% and distance-based surface position estimation. Uses a label map as 
% input. Required a isotropic input map. 
%
%   [Ygmt, Ypp] = cat_vol_pbtsimple(Yp0[, vx_vol, method, interp, enhance])
%
%   Ygmt    .. GM thickness map 
%   Ypp     .. percentage position map 
%   Yp0     .. tissue label map (1-CSF, 2-GM, 3-WM)
%   vx_vol  .. voxel-size (in mm, default = 1)
%   method  .. use voxel- or eikonal distance (0-voxel,1-eikonal)
%              voxel is more robust and faster but less accurate
%   interp  .. aditional interpolation resolution (1-no, 0.5-half res)
%   enhance .. Use a fast thickness estimation (method0 witout further
%              interpolation) to optimzate the PVE between GM and WM by
%              aligning PVE voxels to WM if the GMT is lower and to WM if 
%              the GMT is higher than the average thickness. 
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
  
%#ok<*UNRCH> 

  if ~exist('method', 'var'), method  = 0;   end 
  if ~exist('interp', 'var'), interp  = nan; end % default defintion below
  if ~exist('enhance','var'), enhance = 1;   end


  %% closing for subjects with less/thin WM or with severe SVD/WMHs
    %  - correction by morphops
    % - NISALS
  if 1  
    Yp0o = Yp0; 
    %% laplacian-based closing
    % get smoother lower GM/WM boundary
    Yp0s  = smooth3(Yp0);
    Yp0   = Yp0o;
    Ymsk  = ~cat_vol_morph( Yp0s<2.25 ,'ldo',1.5); 
    Yp0   = max(2.255*Ymsk,Yp0); 
    YmskD = cat_vbdist( single(~Ymsk) , Yp0s<2.75 , repmat(vx_vol(1),1,3)); 
    Yp0   = max(Yp0,min(3, YmskD / vx_vol)); 
    Yp0   = min(Yp0,2.24 + (~cat_vol_morph(Yp0<2.25,'lc'))); 
    % part 2
    Ylp  = single(0.5 - ...
      0.5 * cat_vol_morph( Yp0s<2.25, 'ldo', 0.75, vx_vol) + ...
      0.5 * cat_vol_morph( Yp0s>2.75, 'ldc', 1.50, vx_vol)); 
    Ylp = cat_vol_laplace3R(Ylp,Ylp==0.5,0.01); 
    Yp0(Yp0>2.75) = max( Yp0(Yp0>2.75), 2.5 + Ylp(Yp0>2.75)/2 ); 
    Yp0 = max( Yp0, 2.76 * cat_vol_morph( Yp0>2.25 & ...
      cat_vol_morph(Ylp>0.5 , 'ldc', 3) ,'l') ); 
    
  end



  % (X) GM/WM boundary enhancement: 
  % Use a fast thickness estimation (method0 witout further interpolation) 
  % to optimzate the PVE between GM and WM by aligning PVE voxels to WM  
  % if the GMT is lower and to WM if the GMT is higher than the average 
  % thickness. The effects are larger if intensity normalized rather  
  % than label images are used. 
  % --------------------------------------------------------------------
  if enhance 
    % use fast thickness approximation without enhancement/correction
    Ygmt = cat_vol_pbtsimple(Yp0,vx_vol,method,1,enhance - 1);
   
    % limit and normalize thickness
    % To keep things simple we start with the global median thickness but 
    % it would be also possible to use a very smoothe map here to support
    % more local adjustments and just filter local outliers. 
    Ygmt  = max(1, min(5, Ygmt));
    Ygmt  = cat_vol_approx(Ygmt,1); 
    if 1;; % ##### MANUAL SETTING #####
      % global correction (I would prefere this one as it is simple and stronger)
      gmtmn = max(1.5, min(3.0, median(Ygmt(round(Yp0(:))==2))));
      Ygmtn = ( Ygmt - gmtmn ) / gmtmn;
    else
      % regional corrrection - not optimal as we smooth within the volume
      % but otherwise it would be to slow (+ 10 seconds in Collins) 
      %   gmtmn = max(1.5, min(3.0, median(Ygmt(round(Yp0(:))==2))));
      %   Ygmts = Ygmt - gmtmn; spm_smooth(Ygmts,Ygmts,32 * 4); Ygmts = Ygmts + gmtmn; 
      Ygmts = cat_vol_smooth3X(Ygmt,8); % values smaller <8 have issues (suclal closing) 
      Ygmtn = ( Ygmt - Ygmts ) ./ Ygmts;
    end
    spm_smooth(Ygmtn,Ygmtn,4); % avoid to local results (between 1 and 4)

    % create a corrected T1/segmenation map that will be partialy applied later
    % .4 is already close to equal thickness
    % .2 is quite natural and mostly avoid the underestimations
    Yp0e = ( (Yp0 - 2.5) + (0.3 * Ygmtn) ... 2.5 is our expected boundary threshold
      .* cat_vol_morph( Ygmtn<0 | cat_vol_morph(Yp0>2.5,'lc',1)<0.5 ,'d',2)  ...
      .* smooth3( cat_vol_morph(cat_vol_morph(Yp0>2.25,'do',1),'dd',2) ) ...
      ) + 2.5;
    spm_smooth(Yp0e,Yp0e,1); % keep this map smooth
 
    Ymsk = Yp0>2.1 & Yp0<2.9;    % general GM/WM area 
    Ymsk(smooth3(Ymsk)<0.8) = 0; % remove small real PVE areas that should not be changed!
    Yp0(Ymsk) = Yp0e(Ymsk);      % applay corrected map in GM/WM PVE area

    % post corrections of the GM/WM area 
    if sum(Yp0(:) > 2.5) / sum(Yp0(:) > 0.5) > .3 % critical in case of thin WM! 
      Yp0 = min(2.49 + cat_vol_morph(Yp0>2.5,'ldo',1.5),Yp0); % open WM 
    end
    Yp0 = max(2.51 .* cat_vol_morph(Yp0>2.5,'ldc',1.5),Yp0); % close WM
  end


  % pre-correction of the CSF/GM partial volume area 
  % - Tested in ADHD200NYC with quite similar intensity and position values
  %   but reduction of WM self intersections (9.37% > 8.86%).
  if 1;; % ##### MANUAL SETTING #####
    Ytissue = cat_vol_morph( Yp0>1.75 , 'ldc', .75, vx_vol); 
    Ytissue = cat_vol_morph( Ytissue  , 'dd' , .75, vx_vol);
    Yp0(~Ytissue) = min(Yp0(~Ytissue) , 1.49);
  end



  % (X) Interpolation:
  % --------------------------------------------------------------------
  % Extra interpolation especially for the simple grid-based interpolation 
  % (that does not support partial volume effects) allow improvements if 
  % the deinterpolation includes enough smoothing.
  % However, it is not clear if the increase in processing time is worth it.
  % For the simple vbdist (method=0) it seems to support better thickness,  
  % intensity, and especially position values. Without interpolation small
  % thickness values are otherwise a bit underestimated althought this is 
  % partially helpful for sulcus reconstruction. 
  % For the Eiknonal distance the improvement are less obvious. 
  %
  % This parameter is only for manual tests (especially vbdist)
  %   interp .. 1.0 - no interpolation, 0.5 - half resolution,
  %             0.8 - a bit interpolation (doulbed voxel size) 
  if isnan(interp)
    if method == 0, interp = 1.0; else; interp = 1; end % ########
  end
  if interp ~= 1 
    V = struct('mat', eye(4), 'dim', size(Yp0));
    [Yp0,resI] = cat_vol_resize(max(1,min(3,Yp0)),'interp',V,interp,'cubic'); 
    vx_vol = vx_vol * interp;
  end




  % (1) Distance estimation:
  % --------------------------------------------------------------------
  if method == 0;; % ##### MANUAL SETTING #####
  % Simple voxel-based function without considering partial volume effect 
  % (PVE) or asymetries but also without projection issues!
  % Due to it simpicity this is more robust.
    
    % CSF and WM distance maps (based on a binary boundary without PVE)
    Ycd = cat_vbdist(single(Yp0 < 1.51), Yp0 < 3); 
    Ywd = cat_vbdist(single(Yp0 > 2.51), Yp0 > 1); 

  else
  % More accurate Eikonal-based Euclidean distance estimation. 
  % But I am not really happy with the result, the projection seems to
  % include a lot of problems. 
  % Advantages should be vissible in the motorcortex or the insula.

    % CSF and WM speed maps to handle asymetries 
    Ycf  = max(eps, min(1, ((Yp0 - 1) / 1.1) .^4 )); 
    Ywf  = max(eps, min(1, ((4 - Yp0) / 2.0) .^2 )); 

    % CSF and WM boundary maps with PVE 
    Ycb  = max(0 , min(1,  2.01 - Yp0 )); Ycb(Yp0>2.51) = nan; 
    Ywb  = max(0 , min(1,  Yp0 - 2.01 )); Ywb(Yp0<1.51) = nan; 

    % CSF and WM distance maps with PVE
    Ywd  = cat_vol_eidist(Ywb, Ywf); 
    Ycd  = cat_vol_eidist(Ycb, Ycf);

    if 1
      Ywd = Ywd*.5 + .5*cat_vbdist(single(Yp0 > 2.51), Yp0 > 1); 
      Ycd = Ycd*.5 + .5*cat_vbdist(single(Yp0 < 1.51), Yp0 < 3); 
    end

    clear Ycf Ywf Ywb Ycb;
  end




  % (2) Thickness estimation:
  % --------------------------------------------------------------------
  % Typically, only the sulci are reconstructed but also thins gyral
  % structures suffer from blurring by low resolution or artificats. 
  % The idea is to reconstruct both sulci as well as gyris and just use
  % the minimum thickness. The missing update of the CSF distance is not 
  % optimal. 
  % However, there is some small evidence (ie. better intensity/position 
  % RMSE) and surface look (i.e. less errors in topology correction) that
  % this still support some small benefits.  
  % - Basic tests in ADHD200NYC and Collins. 
  if 0;; % ##### MANUAL SETTING ##### 
  % Using the PBT apporach only to reconstruct the sulci.  
  % As this is the more simple apprach we should keep and test it in 
  % case of pipeline changes.

    % projection-based thickness mapping
    Ygmt = cat_vol_pbtp( round(Yp0) , Ywd, Ycd);

    % minimum to reduce issues with meninges
    Ygmt = min(Ygmt, Ycd + Ywd); 

  else
  % Using PBT to reconstruct the sulci and the gyri. 

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
  %                   Ygmt = cat_vol_approx(Ygmt,1);
  %
  % - need mask smoothing due to interpolation artifacts and to get 
  %   enlarge the filter area (saver but slower is to use the hole 
  %   volume by removing the last msk parameter)  
  Ygmt = simple_approx(Ygmt, 2, smooth3( Yp0>1.05 & Yp0<2.95 )>0 );  % adapt for interpolation?                     
   



  % (4) Estimate percentage position map:
  % --------------------------------------------------------------------
  % We first create a corrected CSF distance map with reconstructed sulci.
  % If gyri were reconstructed too than also the WMD would have to be
  % corrected to avoid underestimation of the position map with surfaces 
  % running to close to the WM.
  YM      = Yp0>=1.51 & Yp0<2.51 & Ygmt>eps;
  Ycdc    = Ycd; Ycdc(YM) = min(Ycd(YM), Ygmt(YM) - Ywd(YM)); 
  Ypp     = zeros(size(Yp0),'single'); Ypp(Yp0>=2.51) = 1;
  Ypp(YM) = Ycdc(YM) ./ (Ygmt(YM) + eps); 
  Ypp(Ypp>1) = 0;
  clear Ycdc YM;


  % pp cleanup (this may comes also later)
  % - Tested on ADHD200NYC with light improvements of intensity and 
  %   position but also less self-intersections (12%>8%)
  if 1;; % ##### MANUAL SETTING #####
    Ypp = max( Ypp , 0.75 * cat_vol_morph( Ypp > .75  , 'ldc' , 1.5)); % avoid holes
    Ypp = min( Ypp , 0.49 * (1 + cat_vol_morph( Ypp > .49  , 'ldo', 2.5 , vx_vol ) )); % avoid islands
  end



  % Voxel-size resolution correction:
  % --------------------------------------------------------------------
  % Because Ycd and Ywd measure a grid-based distance (defined as the 
  % center of a voxel), we have to correct 1 voxel (2x.5 voxel, and finally 
  % correct for the _isotropic_ size of our voxel-grid. This is not
  % required for the Eikonal-based distance estimation! 
  Ygmt = (Ygmt - (method == 0) ) * mean(vx_vol); 




  if interp ~= 1 
    % back to original resolution
    Ygmt = cat_vol_resize(Ygmt, 'deinterp', resI);                         
    Ypp  = cat_vol_resize(Ypp , 'deinterp', resI);                          
  end

 
  % Final smoothing to reduce (topology) artifacts what support lower 
  % intensity and position RMSE (ADHD200NYC, method=0, interp=1):
  %     1.2:  0.1244 / 0.1583
  %     1.0:  0.1237 / 0.1549
  %     0.9:  0.1229 / 0.1517 ***
  %     0.8:  0.1239 / 0.1543
  %     0.6:  0.1259 / 0.1611
  %     0.4:  0.1262 / 0.1614
  % This also expect that the original volume resolution is worse as
  % the PBT resolution. 
  % --------------------------------------------------------------------
  if 1;; % ##### MANUAL SETTING #####
    spm_smooth(Ygmt , Ygmt , 0.9 ); 
    spm_smooth(Ypp  , Ypp  , 0.9 ); 
  end
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

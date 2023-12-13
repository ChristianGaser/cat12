function [Ygmt,Ypp] = cat_vol_pbtsimple(Yp0,vx_vol,opt)
%cat_vol_pbtsimple. Simple cortical thickness/position estimation.  
% Voxel-based distance estimation and projection-based thickness (PBT) 
% and distance-based surface position estimation. Uses a label map as 
% input. Required a isotropic input map. 
%
%   [Ygmt, Ypp] = cat_vol_pbtsimple( Yp0 , vx_vol, opt )
%
%   Ygmt    .. GM thickness map 
%   Ypp     .. percentage position map 
%   Yp0     .. tissue label map (1-CSF, 2-GM, 3-WM)
%   vx_vol  .. voxel-size (in mm, default = 1)
%   opt     .. parameter structure
%    .levels    .. 
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

  if ~exist('opt','var'), opt = struct(); end

  def.levels          = 4; % 0 .. 4 
  def.refine          = 1; % refined WM distance based on CSF distance (myelination correction)
  def.gyrusrecon      = 1; % additional PBT gyri reconstruction 
  def.correctoffeset  = 2;
  def.extendedrange   = 1; 
  opt = cat_io_checkinopt(opt,def);


  % (1) Distance estimation:
  % --------------------------------------------------------------------
  if opt.levels == 0   
  % Classic CSF and WM distance maps based on a binary boundary without 
  % partial volume effect (only for tests).
    Ycd = cat_vbdist(single(Yp0 < 1.5), Yp0 < 3); 
    Ywd = cat_vbdist(single(Yp0 > 2.5), Yp0 > 1);

  elseif opt.levels > 0  &&  opt.refine == 0
  % Simple CSF and WM distance maps with partial volume effect by multiple  
  % distance estimations but without further corrections.
    Ycd = cat_vol_eudist( max(0,min(1,2.0 - Yp0)), Yp0 <= 2.50, opt.levels, opt.correctoffeset); 
    Ywd = cat_vol_eudist( max(0,min(1,Yp0 - 2.0)), Yp0 >= 1.50, opt.levels, opt.correctoffeset);
    opt.extendedrange = 0;

  elseif opt.levels > 0  &&  opt.refine > 0
  % Enhanced CSF and WM distance maps with partial volume effect by   
  % multiple distance estimations but without further corrections.
    [Ycd, Ywd] = cat_vol_cwdist(Yp0,opt);

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
  if  ~opt.gyrusrecon
  % Using the PBT apporach only to reconstruct the sulci.  
  % As this is the more simple apprach we should keep and test it in 
  % case of pipeline changes.

    % projection-based thickness mapping
    Ygmt = cat_vol_pbtp( round(Yp0) , Ywd, Ycd);
    
    if opt.levels >= 0
      Ycd = max(0,Ycd - 0.5); Ywd = max(0,Ywd - 0.5); % now correct also this values
    end

    % minimum to reduce issues with meninges
    Ygmt = min(Ygmt, Ycd + Ywd); 

  else
  % Using PBT to reconstruct the sulci and the gyri. 

    % reconstruct sulci as well as gyri 
    Ygmt1 = cat_vol_pbtp(round(Yp0)  , Ywd, Ycd);  
    Ygmt2 = cat_vol_pbtp(round(4-Yp0), Ycd, Ywd);
    Ygmt2 = max(Ygmt2 , 1.5 .* (Ygmt2>0)); %only in thick regions

    if opt.levels >= 0
      Ycd = max(0,Ycd - 0.5); Ywd = max(0,Ywd - 0.5); % now correct also this values
    end

    % avoid meninges !
    Ygmt1 = min(Ygmt1, Ycd + Ywd);
    Ygmt2 = min(Ygmt2, Ycd + Ywd); 

    % average GM thickness maps
    Ygmt  = min(cat(4, Ygmt1, Ygmt2),[],4);

    clear Ygmt1 Ygmt2
  end




  % (3) Approximation of non GM voxels for resampling:
  % --------------------------------------------------------------------
  Ygmt = cat_vol_approx(Ygmt,'rec-test');                   
   



  % (4) Estimate percentage position map:
  % --------------------------------------------------------------------
  % We first create a corrected CSF distance map with reconstructed sulci.
  % If gyri were reconstructed too than also the WMD would have to be
  % corrected to avoid underestimation of the position map with surfaces 
  % running to close to the WM.
  YM      = Yp0 > 1.5-0.45*opt.extendedrange & Yp0 < 2.5+0.45*opt.extendedrange & Ygmt>eps;
  Ycdc    = Ycd; Ycdc(YM) = min(Ycd(YM), Ygmt(YM) - Ywd(YM)); 
  Ypp     = zeros(size(Yp0),'single'); Ypp(Yp0 >= 2.5+0.45*opt.extendedrange) = 1;
  Ypp(YM) = Ycdc(YM) ./ (Ygmt(YM) + eps); 
  Ypp(Ypp>1) = 0;
  clear Ycdc YM;


  % (5) Voxel-size resolution correction:
  % --------------------------------------------------------------------
  Ygmt = Ygmt * mean(vx_vol); 

end
function [Ycd, Ywd] = cat_vol_cwdist(Yp0,opt)

    % CSF and WM distance 
    Ycd = zeros(size(Yp0),'single'); 
    Ywd = zeros(size(Yp0),'single'); 

    % multi-level distance estimation
    hss = opt.levels; % number of opt.levels (as pairs)
    for si = 1:hss
      offset = max(0,min(.5,si * .5/hss)) / 2; 

      Ycdl  = cat_vbdist(single(Yp0 < ( 1.5 - offset)), Yp0 < 2.5 + 0.45*opt.extendedrange );  
      Ycdh  = cat_vbdist(single(Yp0 < ( 1.5 + offset)), Yp0 < 2.5 + 0.45*opt.extendedrange );  
      
      if opt.extendedrange
        Ycdlc  = max(0,cat_vbdist(single(Yp0 < ( 2.5 - offset)), Yp0 < 2.95 ) - (opt.levels<8));  
        Ycdhc  = max(0,cat_vbdist(single(Yp0 < ( 2.5 + offset)), Yp0 < 2.95 ) - (opt.levels<8));  
        Ycdl   = Ycdl - Ycdlc; 
        Ycdh   = Ycdh - Ycdhc; 
      end

      if opt.correctoffeset
        if opt.correctoffeset==2
          offsetc = cat_stat_nanmedian(Ycdl(Ycdl>0) - Ycdh(Ycdl>0))/2; 
        else
          offsetc = offset; 
        end
        Ycdl(Ycdl>0) = max(eps, Ycdl(Ycdl>0) - (.5 - offsetc )); 
        Ycdh(Ycdh>0) = max(eps, Ycdh(Ycdh>0) + (.5 - offsetc )); 
      end
      Ycd = Ycd + .5/hss  .* Ycdl  +  .5/hss .* Ycdh;
      Ycw = max(0.1,.5/hss - opt.refine * Ycd/10 );
      clear Ycdl Ycdh; 

      % WM distances
      Ywdl  = cat_vbdist(single(Yp0 > ( 2.5 - offset)), Yp0 > 1.5 - 0.45*opt.extendedrange );  
      Ywdh  = cat_vbdist(single(Yp0 > ( 2.5 + offset)), Yp0 > 1.5 - 0.45*opt.extendedrange); 

      if opt.extendedrange
        Ywdlc  = max(0,cat_vbdist(single(Yp0 > ( 1.5 - offset)), Yp0 > 1.05 ) - (opt.levels<8));  
        Ywdhc  = max(0,cat_vbdist(single(Yp0 > ( 1.5 + offset)), Yp0 > 1.05 ) - (opt.levels<8));  
        Ywdl   = Ywdl - Ywdlc; 
        Ywdh   = Ywdh - Ywdhc; 
      end

      if opt.correctoffeset
        if opt.correctoffeset==2
          offsetc = median(Ywdl(Ywdl>0) - Ywdh(Ywdl>0))/2; 
        else
          offsetc = offset; 
        end
        Ywdl(Ywdl>0) = max(eps, Ywdl(Ywdl>0) + (.5 - offsetc )); 
        Ywdh(Ywdh>0) = max(eps, Ywdh(Ywdh>0) - (.5 - offsetc )); 
      end
      Ywd = Ywd + Ycw .* Ywdl  +  (.5/hss*2 - Ycw) .* Ywdh;
    end

    Ycd(Ycd>100) = 0; 
    Ywd(Ywd>100) = 0; 
  
end
function Yd = cat_vol_eudist(Yb,Ymsk,levels,correctoffeset)
%cat_vol_eudist. Euclidean distance estimation to mixed boundaries.
% ...
%
%  Yd = cat_vol_eudist(Yb,Ymsk,levels,correctoffeset)


% add voxel size?


  if ~exist('Ymsk','var'), Ymsk = ~Yb; else, Ymsk = Ymsk>.5; end
  if ~exist('hss','var')
    levels = 4; 
  elseif levels < 0 
    error('cat_vol_eudist:BadLevelValue','Levels must be larger equal 0.') 
  end
  if ~exist('correctoffeset','var'), correctoffeset = 2; end

  % do not estimate the distance for NAN
  Ymsk( isnan(Yb) | (isinf(Yb) & Yb<0) ) = 0; 

  if levels == 0
    % simple single distance estimation 
    Yd  = cat_vbdist(single(Yb > .5), Ymsk );  
  else

    Yd = zeros(size(Yb),'single'); 
    for si = 1:levels
      % estimate the offset of the boundary 
      offset = max(0,min(.5,si * .5/levels)) / 2; 
  
      % estimate the distance to paired sublevels
      Ydl  = cat_vbdist(single(Yb > ( 0.5 - offset)), Ymsk );  
      Ydh  = cat_vbdist(single(Yb > ( 0.5 + offset)), Ymsk );  
  
      % correct for possible outliers
      Ydl(Ydl>max(size(Yb))) = 0; 
      Ydh(Ydh>max(size(Yb))) = 0; 

      % it is possible to correct for the theoretical offset by the
      % voxel-wise partial volume effect if two pure tissues are mixed
      % (typical tissue boundary vs. the myelinated regions)
      if correctoffeset
        if correctoffeset==2
          offsetc = cat_stat_nanmedian(Ydl(Ydl > 0) - Ydh(Ydl > 0))/2; 
        else
          offsetc = offset; 
        end
        Ydl(Ydl>0) = max(eps, Ydl(Ydl>0) - (.5 - offsetc )); 
        Ydh(Ydh>0) = max(eps, Ydh(Ydh>0) + (.5 - offsetc )); 
      end
  
      % add the distance from this level
      Yd = Yd + .5/levels  .* Ydl  +  .5/levels .* Ydh;
  
    end

  end

  % correct for possible outliers
  Yd(Yd > max(size(Yb))) = 0; 
    
end

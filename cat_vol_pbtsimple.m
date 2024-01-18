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

  def.WMTC            = 0; % Partial voxel-based topology correction of small 
                           % WM-defects based on the volume output of cat_vol_genus0 
  def.levels          = 4; % Number of dual distance estimate (requires refinement) with 0 
                           % for the classic approach with 1.5 and 2.5 as CSF and WM boudnary.
                           % Good values are between 2 and 4, whease higher values (>6)
                           % run into problems for interpolation artifacts of close to the 
                           % full tissue class.
  def.refine          = 1; % Refined WM distance based on CSF distance (myelination correction)
  def.gyrusrecon      = 1; % additional PBT gyri reconstruction (TRUE/FALSE)
                           % this reduce thickness overestimations but maybe underestimates  
                           % the outer surface position in gyral regions (running to much inside)  
  def.correctoffeset  = 2; % not really required if no refinement is used
  def.extendedrange   = 1; % may causes closed gyri
  def.sharpening      = 0; % sharpening the Ypp map to support more better resampling to lower resolution for the 0.5 layer - worse
  opt = cat_io_checkinopt(opt,def);


% == preparations ==
% should maybe be better part of create CS

  if 0
    % INITIAL SHARPENING?
    % RD202401: remove this part when the rest is finished
    % not so powerfull as you use already the multple distance estimates
    % should also not work in a label map!

    % sharpen WM
    Yp0x = min(3,max(2,Yp0));
    Yp0x = min(3,max(2,Yp0x + .5*(Yp0x - cat_vol_smooth3X(Yp0x,1)) + .5*(Yp0x - cat_vol_smooth3X(Yp0x,2)) ));
    Yp0(Yp0(:)>2) = Yp0x(Yp0(:)>2);

    % sharpen CSF
    Yp0x = min(2,max(1,Yp0));
    Yp0x = min(2,max(1,Yp0x + .5*(Yp0x - cat_vol_smooth3X(Yp0x,1)) + .5*(Yp0x - cat_vol_smooth3X(Yp0x,2)) ));
    Yp0(Yp0(:)<2 & Yp0(:)>1) = Yp0x(Yp0(:)<2 & Yp0(:)>1);
  end


  if opt.WMTC
  % RD20231224: for WM topology smoothing
  % the idea is to apply the voxel-based correction and _close_ small WM holes 
  % (<10 voxel) for the 2.75 and 2.25 boundaries for values above 2.

    % indicate isolated holes and replace by median of the neighbors
    Yp0(Yp0<0.35 & ~cat_vol_morph(Yp0<1,'l')) = 1;  % close major wholes in the WM 
    Ymsk = Yp0==-1 & cat_vol_morph(Yp0>0.9,'d',1);  % filter small wholes close to the WM
    Yp0 = cat_vol_median3(single(Yp0)+1,Ymsk,~Ymsk)-1; 

    % indicate isolated objects and replace by median of the neighbors
    Yp0(Yp0>0.65 & cat_vol_morph(Yp0==-1,'l')) = -1;
    Ymsk = Yp0>0.95 & cat_vol_morph(Yp0<-0.9,'d',1); 
    Yp0 = cat_vol_median3(single(Yp0)+1,Ymsk,~Ymsk)-1;
  
  % RD20210401: but there are some new background dots
    Ymsk = Yp0==-1 & cat_vol_morph(Yp0>-1,'lc');  % close major wholes in the WM 
    Yp0 = cat_vol_median3(single(Yp0)+1,Ymsk,~Ymsk)-1; 
  end 



  if opt.WMTC
  % topology correction for WM based on the volume output of cat_vol_genus0 

    % Closing of the WM object to avoid smaller holes.
    % Opening is here not useful as the WM !
    tclevels = [2.75, 2.25];
    for ii = 1:1
      for tci = 1:numel(tclevels)
        evalc(sprintf('Yppc = cat_vol_genus0(Yp0, %0.2f,0);',tclevels(tci)));
        %if tclevels(tci) > 2.5
          Yholes = Yppc==1 & Yp0<tclevels(tci); 
          Yholes = Yholes & ~cat_vol_morph(Yholes,'l',[1e4,8]);  % allow only small corrections
          Yp0(Yholes) = tclevels(tci) + 0.01; % close
          clear Yholes;
        %else
          Ybridges = Yppc==0 & Yp0>=tclevels(tci); 
          Ybridges = Ybridges & ~cat_vol_morph(Ybridges,'l',[1e4,8]); % allow only small corrections  
          Yp0(Ybridges) = tclevels(tci) - 0.01; % open
          clear Ybridges;
        %end
      end
    end

    
    % object correction in the background (opening of GM-WM objects)
    tclevels = [1.75, 1.25, 1.01];
    for tci = 1:numel(tclevels)
      Yp0(Yp0>tclevels(tci) & ~cat_vol_morph(Yp0>tclevels(tci), ...
        'ldo',2.5 - tclevels(tci))) = tclevels(tci) - 0.01;  
    end
  end


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
  Ygmt = cat_vol_approx(Ygmt,'nh');                   
   



  % (4) Estimate percentage position map:
  % --------------------------------------------------------------------
  % We first create a corrected CSF distance map with reconstructed sulci.
  % If gyri were reconstructed too than also the WMD would have to be
  % corrected to avoid underestimation of the position map with surfaces 
  % running to close to the WM.
  overrange = .0;  % range extention of the Ypp was resulting in worse
                   % T1 and position values and more topology defects
  YM      = Yp0 > 1.5-0.45*opt.extendedrange & Yp0 < 2.5+0.45*opt.extendedrange & Ygmt>eps;
  Ycdc    = Ycd; Ycdc(YM) = min(Ycd(YM), Ygmt(YM) - Ywd(YM)); 
  Ypp     = zeros(size(Yp0),'single'); Ypp(Yp0 >= 2.5+0.45*opt.extendedrange) = 1;
  Ypp(YM) = Ycdc(YM) ./ (Ygmt(YM) + eps); 
  Yppl    = -2 + 2*smooth3(Ypp>0); Yppl(Ypp>0) = 0;
  Ypph    = 2*smooth3(Ypp>=1);     Ypph(Ypp<1) = 0;
  Ypp     = min(1 + overrange,max(0 - overrange, Ypp + Yppl + Ypph ));
  clear Ycdc YM;


  % (5) Voxel-size resolution correction:
  % --------------------------------------------------------------------
  Ygmt = Ygmt * mean(vx_vol); 


  if 0 
    % FINAL TOPOLOGY CORRECTION
    % final voxel-based topology correction of some levels 
    % however, the topology correction in the surface create should be more powerfull
    tclevels = 1:-0.25:0;
    for tci = 1:numel(tclevels)
      evalc(sprintf('Yppc = cat_vol_genus0(Ypp, %0.2f,0);',tclevels(tci)));
      if tclevels(tci) > 0.5
        Ypp(Yppc==1 & Ypp<tclevels(tci))  = tclevels(tci) + 0.01; % close
      else
        Ypp(Yppc==0 & Ypp>=tclevels(tci)) = tclevels(tci) - 0.01; % open
      end
    end
  end


  % sharpening - important to compensate the downsampling
  if opt.sharpening
    Ypp = Ypp + 1*(Ypp - cat_vol_smooth3X(Ypp,1/mean(vx_vol))) + 2*(Ypp - cat_vol_smooth3X(Ypp,2/mean(vx_vol))) + 4*(Ypp - cat_vol_smooth3X(Ypp,4/mean(vx_vol)));
   
    % final smoothing to prepare surface reconstruction that also correct for WM topology issues 
    spm_smooth(Ypp,Ypp,.5/mean(vx_vol));
    Ypp = min(1 + overrange,max(0 - overrange,Ypp)); 
  
  end
end
% ======================================================================
function [Ycd, Ywd] = cat_vol_cwdist(Yp0,opt)

    % CSF and WM distance 
    Ycd = zeros(size(Yp0),'single'); 
    Ywd = zeros(size(Yp0),'single'); 

    % multi-level distance estimation
    hss = opt.levels; % number of opt.levels (as pairs)
    for si = 1:hss
      offset = max(0,min(.5, si * (.5)/hss)) / 2; 

      Ycdl  = cat_vbdist(single(Yp0 < ( 1.5 - offset)), Yp0 < 2.5 + 0.45*opt.extendedrange ); Ycdl(Ycdl > 1000) = 0; 
      Ycdh  = cat_vbdist(single(Yp0 < ( 1.5 + offset)), Yp0 < 2.5 + 0.45*opt.extendedrange ); Ycdh(Ycdh > 1000) = 0;

      if opt.extendedrange
        % additiona correction for levels>8
        Ycdlc  = max(0,cat_vbdist(single(Yp0 < ( 2.5 - offset)), Ycdl>0 ) - .5 ); Ycdlc(Ycdlc > 1000) = 0;  
        Ycdhc  = max(0,cat_vbdist(single(Yp0 < ( 2.5 + offset)), Ycdh>0 ) - .5 ); Ycdhc(Ycdhc > 1000) = 0;  
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

      % mixing
      Ycd = Ycd + .5/hss  .* Ycdl  +  .5/hss .* Ycdh; 
      Ycw = max(0.1,.5/hss - opt.refine * Ycd/5 );
      %%
      clear Ycdl Ycdh; 


      % WM distances
      Ywdl  = cat_vbdist(single(Yp0 > ( 2.5 - offset)), Yp0 > 1.5 - 0.45*opt.extendedrange ); Ywdl(Ywdl > 1000) = 0; 
      Ywdh  = cat_vbdist(single(Yp0 > ( 2.5 + offset)), Yp0 > 1.5 - 0.45*opt.extendedrange ); Ywdh(Ywdh > 1000) = 0; 

      if opt.extendedrange
        Ywdlc  = max(0,cat_vbdist(single(Yp0 > ( 1.5 - offset)), Yp0 > 1.05 ) - .5); Ywdlc(Ywdlc > 1000) = 0; 
        Ywdhc  = max(0,cat_vbdist(single(Yp0 > ( 1.5 + offset)), Yp0 > 1.05 ) - .5); Ywdhc(Ywdhc > 1000) = 0;  
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

      % mixing
      Ywd = Ywd + Ycw .* Ywdl  +  (.5/hss*2 - Ycw) .* Ywdh;
    end

  
end
% ======================================================================
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

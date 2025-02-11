function [Ygmt,Ypp] = cat_vol_pbtsimple(Yp0,vx_vol,opt)
%cat_vol_pbtsimple. Simple cortical thickness/position estimation.  
% Estimation of voxel-based distances and projection-based thickness (PBT) 
% and surface position based on a label map. Required isotropic input. 
% 
% After further optimisation the function becomes more complex again. 
% However, the supersimple parameter can be used to switch off extra steps.
% Although, the values look quite similar difference in the surface are 
% good detectable for small structures, in particular occipital sulci.
%
%   [Ygmt, Ypp] = cat_vol_pbtsimple( Yp0 , vx_vol, opt )
%
%   Ygmt    .. GM thickness map 
%   Ypp     .. percentage position map 
%   Yp0     .. tissue label map (1-CSF, 2-GM, 3-WM)
%   vx_vol  .. voxel-size (in mm, default=1)

%   opt     .. parameter structure 
%
%    .supersimple (0-no, 1-yes; default=1)
%      WARNING: Activation will run only the basic routines 
%               without enhanced refinements!
%               Although surface intensity and position values are not bad
%               the number of self-intersections indicate further problems.
%               In particular, small sulci are often affected by blurring. 
%      
%    .levels (integer; default=8)
%      Number of dual distance estimations to reduce sampling effects.
%      With logarithmic improvement and good results between 2 and 16. 
%
%    .extendedrange (0-no, 1-yes; default=1)
%      Uses also the PVE range to estimate the distance there. 
%      Important to avoid thickness underestimations.
%
%    .correctoffeset (0-none, 1-fixed, 2-adaptive; default=2)
%      Correction for the offset of boundaries beyond the paired concept.
%      Only minor effects. 
%
%    .distcleanup (0-no, 1-yes; default=1); 
%      Remove distance outliers, e.g., blood vessels, by the histogram but
%      more important by a first very smooth thickness approximation that 
%      define local outliers.  Although, this gives only tiny global 
%      improvements, the local differences are essential, e.g. occipital
%      close to the cerebellum and major veins. 
%
%    .gyrusrecon (0-none, 1-simple, 2-complex; default=3)
%      Additional calls of PBT to reconstruct also the gyri to reduce 
%      thickness overestimations. 
%       0 - none
%           Less accurate thickness with gyral overestimation but more 
%           accurate outer surface compared to method 1.
%       1 - simple approach 
%           Independent reconstruction of sulci and gyri that simply takes
%           the minimum thickness. This causes light underestimations of 
%           the outer surface position in gyral regions (i.e. running too 
%           much inside) but better thickness values. 
%       2 - complex approach (default)
%           Stepwise reconstruction of sulci and gyri (i.e. dependent). 
%           This is sensitive for (blood vessel) artefacts that is handled
%           by distcleanup step.
%
%    .range (real value <=0.5, good between 0.2 and 0.4; default=0.3)
%      Limitation of the offset of multiple thickness levels to avoid 
%      running into partial volume effects with thickness overestimation.
%
%    .keepdetails (0-off, 1-sulci, 2-sulci+gyri; default=1)
%      Enhance thin sulci (and gyri) to avoid blurring. 
%      Small global differences but important in occipial regions. 
%      Gyri enhancement is not optimal yet. 
%
%    .sharpening (0-no, 1-yes; default=1)
%      Further optimization of the maps that can help to save small sulci.
%      Sharpening the Ypp map to support more better resampling to lower 
%      resolution for the 0.5 layer adopted for slight changes that improves
%      especially the position RMSE value by ~0.01 (ie, unclear if this is 
%      really better). 
%
%   See also cat_vol_pbt, cat_vol_createCS[123].
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

%#ok<*UNRCH,*HIST>

  if ~exist('opt','var'), opt = struct(); end
  
  if 0 
    % test old function that resulted in better intensity and position values
    [Ygmt,Ypp] = cat_vol_pbtsimple_goodold(Yp0,vx_vol,opt);
    return
  end

  % Default parameter settings with short impression for some Collins  
  % processing to roughly quantify the global effect. 
  %
  % Values: 
  %  thickness (mn±sd), surface intensity / position RMSE, SI=self-intersections
  %
  % Final tests with the following subjects are important : 
  %  Collins, HR075, 4397, OASIS031 (WMHs), BUSS02 (Child), 
  %  BWPT (Phantom), NISALS_UTR (PVEs), ADHD200NYC (LowQuality/Motion)
  %  
  % Overall, it is important to eep an eye on 
  %  * thickness values >> BWPT
  %  * unterestimation of very thick young cortices >> BUSS02 
  %  * local breaks/defects with thickness underestimation in older subjects with WMHs >> OASIS31, NISALS_UTR
  %  * or in case of motion artefact >> ADHD200NYC 
  %
  def.supersimple     = 0;   % Refined WM distance based on CSF distance (myelination correction)
                             % This internal options is to turn of all feature functions.
                             %   Collins-1:   2.5713 ± 0.6525 mm, 0.0723 / 0.0934,  SI=16.15%
                             %   Collins-0:   2.5633 ± 0.6876 mm, 0.0728 / 0.0805,  SI= 7.67%
                             %                                   -0.0005 / 0.0129,     +8.48%
                             %
  def.WMTC            = 0;   % Voxel-based topology correction of small WM-defects 
                             % based on the volume output of cat_vol_genus0 
                             % ... not really working ... just keep it now for sentimental reasons before it will be remove in future  
                             %
  def.levels          = 8;   % Number of dual distance estimates.
                             % Larger values are more accurate but need more time (log-change).
                             % Good values are between 2 and 16 with limited range parameter.
                             %  
  def.correctoffeset  = 2;   % correct for layer intensity (0-no, 1-yes-simple, 2-yes-complex; tiny improvement) 
                             %   Collins-0:  2.5647 ± 0.6883 mm, 0.0726 / 0.0814,  SI=7.78%
                             %   Collins-1:  2.5633 ± 0.6876 mm, 0.0728 / 0.0805,  SI=7.67%
                             %                                  -0.0002 / 0.0009,    +0.11%  
                             %
  def.extendedrange   = 1;   % may causes closed gyri (0-no, 1-yes)
  def.distcleanup     = 1;   % remove distance outliers, e.g., blood vessels 
                             % (0-no, 1-yes; tiny global improvement but locally pretty important)     
                             %   Collins-1:  2.6693 ± 0.6199 mm, 0.0734 / 0.0746,  SI=4.975%
                             %   Collins-0:  2.6493 ± 0.6147 mm, 0.0728 / 0.0730,  SI=4.530%
                             %                                  +0.0006 /+0.0016,    +0.445%
                             %
  def.gyrusrecon      = 2;   % Additional PBT gyri reconstruction 
                             % this reduce thickness overestimations but maybe underestimates  
                             % the outer surface position in gyral regions (running too much inside)  
                             %   0 - none                     0.0666 / 0.0856, SI=10.19%
                             %   1 - simple old approach      0.0677 / 0.0838, SI=17.09%
                             %   2 - more complex approach    0.0729 / 0.0745, SI= 4.30% *** 
                             %
  def.range           = 0.3; % Default value for range extension (should be between 0.2-0.4)
                             %
  def.keepdetails     = 1;   % enhance thin (occipital) sulci (and gyri) to avoid blurring 
                             % (0-no, 1-yes; 2-yes(+gyri) worse values but pretty important!
                             % (0-off, 1-sulci, 2-sulci+gyri)
                             %   Collins-0:  2.6413 ± 0.6082 mm, 0.0708 / 0.0722,  SI=4.585%    
                             %   Collins-1:  2.6509 ± 0.6151 mm, 0.0728 / 0.0730,  SI=5.365%
                             %                                  -0.0020 /-0.0012,    -1.220% 
                             %
  def.sharpening      = 1;   % sharpening the Ypp map to avoid blurring while 
                             % resampling to lower resolution
                             % (0-no, 1-yes; reduce collisions)
                             %   Collins-0:  2.5654 ± 0.6841 mm, 0.0730 / 0.0845,  SI=8.22% 
                             %   Collins-1:  2.5633 ± 0.6876 mm, 0.0728 / 0.0805,  SI=7.67%
                             %                                   0.0002 / 0.0040,  SI=0.55%
                             %

  opt = cat_io_checkinopt(opt,def);

  % just a fast option to switch of the extra functions
  if opt.supersimple % avoid extras
    opt.WMTC          = 0; % WM topology correction
    opt.levels        = 0; % number of dual distance measurements
    opt.correctoffset = 0; % correct of dual displacements (in principle not required)
    opt.extendedrange = 0; % use further PVE area
    opt.distcleanup   = 0; % important correction of local outliers like blood vessels
    opt.gyrusrecon    = 0; % reconstruct gyri (in very young/old subjects with thin WM)
    opt.keepdetails   = 0; % important for small (occipital) sulci
    opt.sharpening    = 0; % reduce blurring of the Ypp map 
  end

  % display paras for developers
  if 0 % cat_get_defaults('extopts.expertgui') > 1
    fprintf('\nPBTsimple parameters:\n'); 
    disp(opt);
  end


  
  %  WM topology correction (relevant in thin WM subjects such as OASIS31):
  if opt.WMTC 
    Yp0 = WMTC(Yp0);
  end



  % (1) Distance estimation:
  % --------------------------------------------------------------------
  if opt.levels == 0 
  % Classic CSF and WM distance maps based on a binary boundary without 
  % partial volume effect (only as simplest approach and for tests).
    Ycd = cat_vbdist( single(Yp0 < 1.5), Yp0 < 2.5 + opt.range*opt.extendedrange ); 
    Ywd = cat_vbdist( single(Yp0 > 2.5), Yp0 > 1.5 - opt.range*opt.extendedrange );

  else 
    if opt.supersimple
    % Simple CSF and WM distance maps with partial volume effect by multiple  
    % distance estimations but without further corrections.
      Ycd = cat_vol_eudist( max(0,min(1,2.0 - Yp0)), ...
        Yp0 <= 2.5 + opt.range*opt.extendedrange, ...
        opt.levels, opt.correctoffeset); 
      Ywd = cat_vol_eudist( max(0,min(1,Yp0 - 2.0)), ...
        Yp0 >= 1.5 - opt.range*opt.extendedrange, ...
        opt.levels, opt.correctoffeset);

    else
      % Enhanced CSF and WM distance maps with partial volume effect by   
      % multiple distance estimations but without further corrections.
      [Ycd, Ywd] = cat_vol_cwdist(Yp0, vx_vol, opt);
      
    end
  end
 


  % (2) Thickness estimation:
  % --------------------------------------------------------------------
  % Typically, only the sulci are reconstructed but also thins gyral
  % structures suffer from blurring by low resolution or artefacts. 
  % The idea is to reconstruct both sulci as well as gyri and just use
  % the minimum thickness. The missing update of the CSF distance is not 
  % optimal. 
  % However, there is some small evidence (ie. better intensity/position 
  % RMSE) and surface look (i.e. less errors in topology correction) that
  % this still support some small benefits.  
  % - Basic tests in ADHD200NYC and Collins. 
  if  opt.gyrusrecon == 0
  % Using the PBT approach only to reconstruct the sulci.
    distcorval = 0.5; % in theory 0.5 

    % remove highly distant outliers
    [Yp0, Ywd, Ycd] = cleanupVessels(Yp0, Ywd, Ycd, opt.distcleanup);
  
    % projection-based thickness mapping
    Ygmt = cat_vol_pbtp( round(Yp0) , Ywd, Ycd);
  
    % now correct also these values (it is a bit better this way)
    Ycd = max(0,Ycd - distcorval); Ywd = max(0,Ywd - distcorval); 
    
    % minimum to reduce issues with meninges
    Ygmt = min( Ycd + Ywd , max(0,Ygmt - distcorval - distcorval * (Ygmt < (Ycd+Ywd)) ) ); 
    Ygmt = min( Ygmt , cat_vol_median3(Ygmt,Ygmt>0,Ygmt>0) );

  elseif opt.gyrusrecon == 1 % EXTRA 
  % Using PBT to reconstruct the sulci and the gyri to get the minimal 
  % thickness. Classic quite simple approach. 
    distcorval = 0.5; % in theory 0.5 

    % remove highly distant outliers
    [Yp0, Ywd, Ycd] = cleanupVessels(Yp0, Ywd, Ycd, opt.distcleanup);
  
    % reconstruct sulci as well as gyri 
    Ygmt1 = cat_vol_pbtp(round(Yp0)  , Ywd, Ycd);  
    Ygmt2 = cat_vol_pbtp(round(4-Yp0), Ycd, Ywd);
    Ygmt2 = max(Ygmt2 , 1.75 / mean(vx_vol) .* (Ygmt2>0)); %only in thick regions

    % now correct also this values
    Ycd = max(0,Ycd - distcorval); Ywd = max(0,Ywd - distcorval); 
  
    % avoid meninges !
    Ygmt1 = min( max(0,Ygmt1 - distcorval - distcorval * (Ygmt1 < (Ycd+Ywd)) ) , Ycd + Ywd); 
    Ygmt2 = min( max(0,Ygmt2 - distcorval - distcorval * (Ygmt2 < (Ycd+Ywd)) ) , Ycd + Ywd); 

    % average GM thickness maps
    Ygmt  = min(cat(4, Ygmt1, Ygmt2),[],4); clear Ygmt1 Ygmt2
    Ygmt  = cat_vol_median3(Ygmt,Ygmt>0,Ygmt>0);

  elseif opt.gyrusrecon == 2 % EXTRA
    % Using PBT to reconstruct the sulci and then the gyri.
    distcorval = 0.5; % in theory 0.5 

    % correct general voxel offset (vbdist quantifies the distance to the 
    % object grid points rather than the boundary between object and background. 
    Ycd = max(0,Ycd - distcorval); Ywd = max(0,Ywd - distcorval); 

    % remove highly distant outliers
    [Yp0c, Ywdc, Ycdc] = cleanupVessels(Yp0, Ywd, Ycd, opt.distcleanup); 

    % Optimized correction factor to avoid defects that is caused by the PBT 
    % mapping that maps the maximum distance so without any distnace between
    % sulcal banks. The correction should be between 0.0 and 0.25 being as 
    % small as possible (for a correction of 0.5 the voxel would be already
    % outside of the GM). Arbitrary correction by the half of "-0.01 - 0.15" 
    pbtsulccor = @(Ygmtx, Ycdx, Ywdx) max(0,Ygmtx -0.01 - 0.15 .* (Ygmtx < (Ycdx + Ywdx))); 

    % 1. sulcus-reconstruction: 
    %    first, we just use the classic estimation to get a corrected CSF distance
    Ygmt1 = cat_vol_pbtp( round(Yp0c) , Ywdc, Ycdc);
    Ygmt1 = pbtsulccor(Ygmt1, Ycdc, Ywdc);
    Ygmt1 = cat_vol_approx(Ygmt1, 'rec');            
    Ycdc  = min(Ygmt1,Ycdc); Ywdc = min(Ygmt1,Ywdc); 
    YM    = Ywdc>0 & Ycdc>0; Ycdc(YM) = min(Ycdc(YM), Ygmt1(YM) - Ywdc(YM)); % update
    
    % 2. gyrus reconstruction: 
    %    now, we can process the inverse/dual case go get
    Ygmt2 = cat_vol_pbtp(round(4 - Yp0), Ycdc, Ywdc);
    Ygmt2 = pbtsulccor(Ygmt2, Ycdc, Ywdc);
    Ygmt2 = cat_vol_approx(Ygmt2, 'rec');
    YM    = YM & Ycdc > median(Ygmt2(:));
    Ywdc  = Ywd; Ywdc(YM) = min(Ywdc(YM), Ygmt2(YM) - Ycdc(YM)); clear Ygmt2 YM; % Ywdc=Ywd  is required to avoid meninges !
    
    % 3. final estimation: 
    %    having now the correct distance values the finale thickness estimations are applied
    Ygmt  = cat_vol_pbtp( round(Yp0) , Ywdc, Ycdc);
    Ygmt  = pbtsulccor(Ygmt, Ycdc, Ywdc); 

    % remove thick outliers to avoid meninges  
    Ygmt( Ygmt > Ygmt1 .* (median(Ygmt1(Ygmt1(:)>0)) / median(Ygmt(Ygmt(:)>0))) ) = 0; clear Ygmt1;
   
    % final update of distance functions
    Ycd = Ycdc; Ywd = Ywdc;
  end



  % (3) Approximation of non GM voxels for resampling:
  % --------------------------------------------------------------------
  Ygmt = cat_vol_approx(Ygmt,'nh');   



  % (*) EXTRA: optimize map - enhancement of thin areas
  if opt.keepdetails 
    [Ywd, Ycd, Ygmt] = keepdetails(Yp0, Ywd, Ycd, Ygmt, vx_vol, ...
      opt.range * opt.extendedrange, opt.keepdetails);
  end      



  % (4) Estimate percentage position map:
  % --------------------------------------------------------------------
  % We first create a corrected CSF distance map with reconstructed sulci.
  % If gyri were reconstructed too than also the WMD would have to be
  % corrected to avoid underestimation of the position map with surfaces 
  % running to close to the WM.
  YM       = Yp0 > 1.5 - opt.range * opt.extendedrange & ...
             Yp0 < 2.5 + opt.range * opt.extendedrange & ...
             Ygmt > eps;
  Ycdc     = Ycd; 
  Ycdc(YM) = min(Ycd(YM), Ygmt(YM) - Ywd(YM)); 
  Ypp      = zeros(size(Yp0),'single'); 
  Ypp(Yp0 >= 2.5 + opt.range * opt.extendedrange) = 1;
  Ypp(YM)  = max(0,Ycdc(YM) ./ (Ygmt(YM) + eps)); 
  Ypp(Ypp<.9 & Yp0>2.5) = 1; 
  clear Ycdc YM;



  % (5) Voxel-size resolution correction:
  % --------------------------------------------------------------------
  Ygmt = Ygmt * mean(vx_vol); 
  

  % (*) EXTRA: magic sharpending function - enhancement of fine structures
  if opt.sharpening 
    [Ypp, Ygmt] = sharpening(Ypp, Ygmt, vx_vol);
  end
end
% ======================================================================
function Yp0 = WMTC(Yp0)
% WMTC. WM topology correction.

  numvx = 8; % size of correction 
  
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
  Ymsk = Yp0==-1 & cat_vol_morph(Yp0>-1,'ldc',1.5);  % close major wholes in the WM 
  Yp0 = cat_vol_median3(single(Yp0)+1,Ymsk,~Ymsk)-1; 


  % topology correction for WM based on the volume output of cat_vol_genus0 
  % Closing of the WM object to avoid smaller holes.
  % Opening is here not useful as the WM !
  tclevels = [2.75, 2.25]; 
  for ii = 1:1 % iterative ? 
    for tci = 1:numel(tclevels)
      evalc(sprintf('Yppc = cat_vol_genus0(Yp0, %0.2f,0);',tclevels(tci)));
      if 1 %tclevels(tci) > 2.5 % only closing
        Yholes = Yppc==1 & Yp0<tclevels(tci); 
        Yholes = Yholes & ~cat_vol_morph(Yholes,'l',[1e4,numvx]);  % allow only small corrections
        Yp0(Yholes) = tclevels(tci) + 0.01; % close
        clear Yholes;
      else
        Ybridges = Yppc==0 & Yp0>=tclevels(tci); 
        Ybridges = Ybridges & ~cat_vol_morph(Ybridges,'l',[1e4,numvx]); % allow only small corrections  
        Yp0(Ybridges) = tclevels(tci) - 0.01; % open
        clear Ybridges;
      end
    end
  end

  
  % object correction in the background (opening of GM-WM objects)
  tclevels = [1.75, 1.25, 1.01];
  for tci = 1:numel(tclevels)
    Yp0(Yp0>tclevels(tci) & ~cat_vol_morph(Yp0>tclevels(tci), ...
      'ldo',2.5 - tclevels(tci))) = tclevels(tci) - 0.01;  
  end
end
% ======================================================================
function [Ywd, Ycd, Ygmt] = keepdetails(Yp0, Ywd, Ycd, Ygmt, vx_vol, extendedrange,level)   
% Although distances and thickness are quite good, PBT slightly trend to  
% over-estimate the thickness in sulcal regions without CSF as the middle 
% voxel is counted for both sides (simplified). Moreover, initial surface 
% are partially created just on the original internal resolution (~1 mm), 
% what result in blurring of small sucli. To keep these structures, we 
% optimize regions where blurring/closing is happening by making them a 
% bit thinner and correcting the CSF distance (keepdetails>0). The dual
% operation (of avoiding of blurring small gyri) can also be applied but
% is expected to have less effects as these structures are already quite
% save by the WM (keepdetails>1). 
%
% 1         -1.0*ppth ... 0           2.5235 ± 0.5959 mm      0.0686 / 0.0605     0.74% (9.16 mm²)    6.0 / 1.0 / 0.42% **** 
% 3         -0.8*ppth ... 0           2.5403 ± 0.5872 mm      0.0682 / 0.0608     0.81% (10.01 mm²)   6.0 / 1.0 / 0.42% **** 
% 2         -0.5*ppth ... 0           2.5671 ± 0.5727 mm      0.0677 / 0.0621     0.85% (10.49 mm²)   10.0 / 2.0 / 0.46%  
% 2         -0.5*ppth ... .1          2.5515 ± 0.5863 mm      0.0690 / 0.0620     0.80% (9.87 mm²)    6.0 / 1.0 / 0.42% 
 
  Ygmto = Ygmt; 

  % estimate current percentage map (same as bellow) to identify
  % and correct problematic areas
  YM      = Yp0 > 1.5 - extendedrange & ...
            Yp0 < 2.5 + extendedrange & ...
            Ygmt > eps;
  Ycdc    = Ycd; Ycdc(YM) = min(Ycd(YM), Ygmt(YM) - Ywd(YM)); 
  Ypp     = zeros(size(Yp0),'single'); Ypp(Yp0 >= 2.5 + extendedrange) = 1;
  Ypp(YM) = Ycdc(YM) ./ (Ygmt(YM) + eps); 

  % reduce thickness and CSF distance if this results in smoothing
  % do this only in thin regions
  for ppth = 1:-.1:0.1
    Yblurredsulci = Ypp<ppth & cat_vol_morph(Ypp>=ppth,'c',1) & Ygmt < median(Ygmt(:))/2;  
    Yblurredsulci = cat_vol_smooth3X(Yblurredsulci, 2);
    Ygmt = max( ...
      max(0,Ygmto - 1/mean(vx_vol)), ...
      Ygmt .* max(0.5,1 - 0.8 .* ppth*Yblurredsulci) );  
    Ycd(Yblurredsulci>0) = max(0,min(Ycd(Yblurredsulci>0), Ygmt(Yblurredsulci>0) - Ywd(Yblurredsulci>0)));
  end   

  % reduce thickness and WM distance if this results in smoothing
  % ... slow, worse values, 
  % - helpful to remove blue outliers?
  if level == 2 
    for ppth = 0.9:-0.1:0.1 
      Yblurredgyri = Ypp>=ppth & cat_vol_morph(Ypp<ppth,'c',1) & Ygmt > median(Ygmt(:))/2; 
      Yblurredgyri = cat_vol_smooth3X(Yblurredgyri, 2);
      Ygmt = max( ...
        min(0,Ygmto + 1/mean(vx_vol)), ...
        Ygmt .* min(1.5,1 + .1 .* ppth*Yblurredgyri));   
      Ywd(Yblurredgyri>0)  = max(0,min(Ywd(Yblurredgyri>0), Ygmt(Yblurredgyri>0) + Ycd(Yblurredgyri>0)));
    end 
  end
end
% ======================================================================
function [Ypp, Ygmt] = sharpening(Ypp, Ygmt, vx_vol)
% sharpening. Reduce  
%
%  [Ypp, Ygmt] = sharpening(Ypp, Ygmt, vx_vol)
%

  Ypp0 = Ypp; 
 
  smoothfactor = .02; 
  smoothoffset = 1; 
  modthickness = 2; 

  % create a sharpend version of the Ypps that enphalize regions that were smoothed 
  Ypps = Ypp + smoothfactor*1*(Ypp - cat_vol_smooth3X(Ypp,1/mean(vx_vol)) + smoothoffset * 0.000) + ...
               smoothfactor*2*(Ypp - cat_vol_smooth3X(Ypp,2/mean(vx_vol)) + smoothoffset * 0.002) + ... 
               smoothfactor*4*(Ypp - cat_vol_smooth3X(Ypp,4/mean(vx_vol)) + smoothoffset * 0.004);

  % final smoothing to prepare surface reconstruction that also correct for WM topology issues 
  spm_smooth(Ypps,Ypps,.5/mean(vx_vol)); 
  Ypps = min(Ypp>0, max(0,Ypps)); % remove background artifacts  

  % New version that only change mid-position values relevant for surface
  % creation. Although this avoid the old binary-like better thickness 
  % cannot be expected
  Ypp  = min(Ypp,max(0.49,Ypps));
  Ypp  = max(Ypp,min(0.51,Ypps)); 

  % adopt thicnkess
  Ymt  = smooth3(Ypp - Ypp0) * modthickness; 
  Ygmt = Ygmt + Ymt;  
end
% ======================================================================
function [Yp0, Ywd, Ycd] = cleanupVessels(Yp0, Ywd, Ycd, doit)
% Filtering the Ywd to avoid mapping of artifacts by blood vessels or  
% other uncleared tissues and regions. 

  if doit
    % == (1) cleanup distance values ==
    % estimate histogram
    [hh,hr] = hist(Ywd(Ywd>eps) ./ Yp0(Ywd>eps),0:.1:50); 
    maxdist = max( hr( cumsum(hh)/sum(hh) < .99 ) );
    Ymsk = (Ywd ./ Yp0) > maxdist;
    % correction
    Yp0(Ymsk) = 1; 
    Ycd(Ymsk) = 0; 
    Ywd(Ymsk) = 0; 
  
    
    % == (2) cleanup by thickness values == 
    % estimate thickness
    Ygmt  = cat_vol_pbtp( round(Yp0) .* (Yp0 > 1.5) , Ywd, Ycd);
    % remove outliers
    gmtmn = cat_stat_nanmean(Ygmt(:)); 
    gmtsd = cat_stat_nanstd(Ygmt(:)); 
    Ymsk  = Ywd < (Ygmt + gmtsd)  &  Ywd < (gmtmn + gmtsd);
    Ygmt  = Ygmt .* Ymsk; 
    gmtsd = cat_stat_nanstd(Ygmt(:)); 
    % approximate and strong smoothing
    Ygmt  = cat_vol_approx(Ygmt,'nh');
    Ygmt  = cat_vol_smooth3X(Ygmt,16);
    
    % final correction
    Ymsk  = Ywd > (Ygmt + gmtsd)  &  Ywd > (gmtmn + gmtsd);
    Ywd(Ymsk) = min(Ywd(Ymsk), max( Ygmt(Ymsk) + gmtsd , Ygmt(Ymsk)*1.25 )); 
    Yp0(Ymsk) = 1.5; 
    Ycd(Ymsk) = 0.5; 
  end
end
% ======================================================================
function Yd = cat_vol_eudist(Yb, Ymsk, levels, correctoffeset)
%cat_vol_eudist. Euclidean distance estimation to mixed boundaries.
%
%  Yd = cat_vol_eudist(Yb, Ymsk, levels, correctoffeset)
% 
%  Yd              .. distance map
%  Yb              .. boundary map (with PVE)
%  Ymsk            .. masked regions for distance estimation 
%
%  levels          .. number of dual distance measurements
%                     optimimum between 2 and 4  
%  correctoffeset  .. use generalized correction for the additional distance
%                     estimations, eg., for a more WM like value of 2.75 all
%                     distance values are assumed to be over
%                     (0 - none, 1 - default difference, 2 - estimated difference)
%
%  see also cat_vol_pbtsimple:cat_vol_cwdist
%

%
% Todo: 
%  * Add voxel size?
%  * Use as esternal function?     
%

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
% ======================================================================
function [Ycd, Ywd] = cat_vol_cwdist(Yp0,vx_vol,opt)
%cat_vol_cwdist. Estimation of CSF and WM distance in a label map Yp0.
% 
% [Ycd, Ywd] = cat_vol_cwdist(Yp0,opt)
%
% Ycd, Ywd         .. CSF and WM distance maps
% opt              .. parameter structure
%  .levels         .. number of dual distance measurements
%  .extendedrange  .. estimate values also beyond the boundary to improve
%                     thickness mapping
%  .correctoffeset .. use generalized correction for the additional distance
%                     estimations, eg., for a more WM like value of 2.75 all
%                     distance values are assumed to be over
%                     (0 - none, 1 - default difference, 2 - estimated difference)
%  .range          .. limitation to avoid bias by interpolation overshoot 
% 

    % CSF and WM distance 
    Ycd = zeros(size(Yp0),'single'); 
    Ywd = zeros(size(Yp0),'single'); 

    opt.extendedrange = opt.extendedrange * opt.range; 

    % additional correction map for values behind tissue boundary, e.g., 
    % for the WMD we estimate the distance from the GM/CSF boundary to 
    % limit WMD values to the maximal thickness value
    % same idea as below
    if opt.extendedrange
      Ycdlc = Ycd; Ycdhc = Ycd; 
      Ywdlc = Ywd; Ywdhc = Ywd; 
      hss = opt.levels; % number of opt.levels (as pairs)
      for si = 1:hss
        offset = max(0,min(opt.range, opt.range * si/(hss+1))); 

        Ycdlc  = Ycdlc + 1/hss * max(0,cat_vbdist(single(Yp0 < ( 2.5 - offset)), Yp0 < 2.5 + opt.extendedrange ) - .5); 
        Ycdhc  = Ycdhc + 1/hss * max(0,cat_vbdist(single(Yp0 < ( 2.5 + offset)), Yp0 < 2.5 + opt.extendedrange ) - .5);  

        Ywdlc  = Ywdlc + 1/hss * max(0,cat_vbdist(single(Yp0 > ( 1.5 - offset)), Yp0 > 1.5 - opt.extendedrange ) - .5); 
        Ywdhc  = Ywdhc + 1/hss * max(0,cat_vbdist(single(Yp0 > ( 1.5 + offset)), Yp0 > 1.5 - opt.extendedrange ) - .5); 
      end
      Ycdlc(Ycdlc > 1000) = 0; Ycdhc(Ycdhc > 1000) = 0; 
      Ywdlc(Ywdlc > 1000) = 0; Ywdhc(Ywdhc > 1000) = 0;  
    end


    % multi-level distance estimation
    hss = opt.levels; % number of opt.levels (as pairs)
    for si = 1:hss
      offset = max(0,min(opt.range, opt.range * si/(hss+1))); 

      Ycdl  = cat_vbdist(single(Yp0 < ( 1.5 - offset)), Yp0 < 2.5 + opt.extendedrange ); Ycdl(Ycdl > 1000) = 0; 
      Ycdh  = cat_vbdist(single(Yp0 < ( 1.5 + offset)), Yp0 < 2.5 + opt.extendedrange ); Ycdh(Ycdh > 1000) = 0;

      if opt.extendedrange
        Ycdl   = Ycdl - Ycdlc; 
        Ycdh   = Ycdh - Ycdhc; 
      end

      if opt.extendedrange
        if opt.correctoffeset==2
          offsetc = offset/mean(vx_vol) + (cat_stat_nanmedian(Ycdl(Ycdl>0 & Ycdh>0) - Ycdh(Ycdl>0 & Ycdh>0)))/2; 
        else
          offsetc = offset/mean(vx_vol); 
        end
        Ycdl(Ycdl>0) = max(eps, Ycdl(Ycdl>0) - (.5 + offsetc) ); 
        Ycdh(Ycdh>0) = max(eps, Ycdh(Ycdh>0) + (.5 + offsetc) ); 
      end

      % mixing
      Ycd = Ycd + .5/hss  .* Ycdl  +  .5/hss .* Ycdh; 
      %Ycw = max(0.1,.5/hss - opt.refine * Ycd/5 );
      % idea was to could the boundaries different depending on the CSF distance
      clear Ycdl Ycdh; 


      % WM distances
      Ywdl  = cat_vbdist(single(Yp0 > ( 2.5 - offset)), Yp0 > 1.5 - opt.extendedrange ); Ywdl(Ywdl > 1000) = 0; 
      Ywdh  = cat_vbdist(single(Yp0 > ( 2.5 + offset)), Yp0 > 1.5 - opt.extendedrange ); Ywdh(Ywdh > 1000) = 0; 

      if opt.extendedrange
        Ywdl   = Ywdl - Ywdlc; 
        Ywdh   = Ywdh - Ywdhc; 
      end
      if opt.correctoffeset
        if opt.correctoffeset==2
          offsetc = offset + (cat_stat_nanmedian(Ywdl(Ywdl>0 & Ywdh>0) - Ywdh(Ywdl>0 & Ywdh>0)))/2; 
        else
          offsetc = offset; 
        end
        Ywdl(Ywdl>0) = max(eps, Ywdl(Ywdl>0) + (.5 + offsetc) ); 
        Ywdh(Ywdh>0) = max(eps, Ywdh(Ywdh>0) - (.5 + offsetc) ); 
      end

      % mixing
      Ywd = Ywd + .5/hss .* Ywdl  +  .5/hss .* Ywdh;
      %Ywd = Ywd + Ycw .* Ywdl  +  (1/hss - Ycw) .* Ywdh;
    end

  
end
% ======================================================================
% 20240217
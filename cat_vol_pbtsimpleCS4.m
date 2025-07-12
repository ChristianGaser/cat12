function [Ygmt,Ypp,Yp0] = cat_vol_pbtsimpleCS4(Yp0,vx_vol,opt)
%cat_vol_pbtsimple. Simple cortical thickness/position estimation.  
% Estimation of voxel-based distances and projection-based thickness (PBT) 
% and surface position based on a label map. Required isotropic input. 
% 
% After further optimisation the function becomes more complex again. 
% However, the supersimple parameter can be used to switch off extra steps.
% Although, the values look quite similar difference in the surface are 
% good detectable for small structures, in particular occipital sulci.
%
%   [Ygmt, Ypp, Yp0] = cat_vol_pbtsimple( Yp0 , vx_vol, opt )
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
%    .NBV (0-no, 1-yes; default=1)
%      Additional, new blood vessel correction based on a WM region growing.
%
%    .myelinCorrection (0-no, 0.25-light, 0.33-default, 0.5-strong, 1.0-maximum)
%      Correction of myelinated GM areas that often result in severe GM 
%      unterestimations. The code is experimental code from the (pre) LAS
%      correction that estimated the GM and WM thickness to extend thin GM 
%      areas if there is a lot of GM-WM PVE behind.
%      This correction is not fully correct but the introduced error is 
%      expected to be smaller than before. 
%
%   See also cat_vol_pbt, cat_vol_createCS[23].
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
                             %
  def.levels          = 4;   % Number of dual distance estimates. 
                             % Larger values are more accurate but need more time (log-change).
                             % Good values are between 1 and 4.
                             %  
  def.correctoffeset  = 0;   % correct for layer intensity 
                             % (0-no, 1-yes-simple, 2-yes-complex; tiny improvement) 
                             %
  def.extendedrange   = 1;   % Estimate the distance from boundary A to boundary B 
                             % also for voxels beyond B with a correction of the 
                             % extra distance, to stabilize the values and avoid
                             % underestimations (0-no, 1-yes)
                             %
  def.range           = 0.3; % Default value for range extension (should be between 0.2-0.4)
                             %
  def.keepdetails     = 1;   % enhance thin (occipital) sulci (and gyri) to avoid blurring 
                             % (0-no, 1-yes (sulci); 2-yes(sulci+gyri) 
                             % worse values but pretty important!
                             %
  def.sharpening      = 1;   % sharpening the Ypp map to avoid blurring while 
                             % resampling to lower resolution
                             % (0-no, 1-yes; reduce collisions)
                             %
  def.NBVC             = 1;  % new blood vessel correction (RD202503)
                             %
  def.eidist           = 0;  % distance metric (0-vbdist,1-eidist)
                             %  vbdist is faster but do not support to consider other thissue boundaries 
                             %  eidist is much slower but support other tissue boudaries
                             %
  def.myelinCorrection = .3; % correction for large cortical myelination artefacts that 
                             % cause strong local underestimation that is similar to the 
                             % extended LAS correction (cat_main_correctmyelination) (RD202503)
                             %
  def.gyrusrecon       = 1;  % use PBT also to reconstruct gyri  
  def.PVErefinement    = 1;  % 
  def.verb             = 0;  % be verbose

  opt = cat_io_checkinopt(opt,def);

  % just a fast option to switch of the extra functions
  if opt.supersimple % avoid extras
    opt.levels        = 0; % number of dual distance measurements
    opt.correctoffset = 0; % correct of dual displacements (in principle not required)
    opt.extendedrange = 1; % use further PVE area
    opt.keepdetails   = 0; % important for small (occipital) sulci
    opt.sharpening    = 0; % reduce blurring of the Ypp map 
    opt.NBVC          = 0; % new additional blood vessel correction (RD202503)
    opt.myelinCorrection = 0; % myelination correction
  end
  

  if opt.extendedrange
    % extend hard cuts, by a low value just to have a broader estimate
    % - could be improved by a basic WMD and GMT estimate to assure a
    %   local minimum thickness related to the GMT distance
    Yp0 = max(Yp0,(Yp0==1 & cat_vol_morph(Yp0==2,'d')) * 1.2);
  end

  % close holes (important for SPM with unsufficient WM correction)
if 0  
  Yp0 = max(Yp0, 3.0 * smooth3(cat_vol_morph(Yp0>2.75,'ldc',1.5)) ); 
  Yp0 = max(Yp0, 2.5 * smooth3(cat_vol_morph(Yp0>2.25,'ldc',1.5)) ); 
end

  %% RD202503: new blood vessel correction 
  if opt.NBVC, Yp0 = NBVC(Yp0,vx_vol); end
  if opt.myelinCorrection, Yp0 = myelincorrection(Yp0,vx_vol,opt); end
 
  c = clock; %#ok<*CLOCK>

  % estimation of distance measures ... 
  % - the fancy estimation does not show useful advantage
  % - the eikonal distance (more correct) is much slower (one level eidist ~60s vs. 4 level vbdist ~15s) 
  %   but support assymetrical mapping  and not really worth it
  %[Ycd0,Ywd0] = cat_vol_PVEdist(Yp0, opt.PVErefinement ); % this function was designed to optimize the GM-WM PVE but is not fully working yet
  [Ycd0, Ywd0] = cat_vol_cwdist(Yp0,opt,vx_vol); 
  

  % raw thickness maps
  Ygmtw0 = cat_vol_pbtp( round(Yp0)   , Ywd0, Ycd0); Ygmtw0(Ygmtw0>1000) = 0; 
  Ygmtc0 = cat_vol_pbtp( 4-round(Yp0) , Ycd0, Ywd0); Ygmtc0(Ygmtc0>1000) = 0; 


  % correct reconstrution overestimation 
  % - important to keep small sulci open, eg. rh.BWPT central and CC sulcus 
  % - correct thinner areas stronger to improve reconstruction
  pbtsulccor = @(Ygmtx, Ycdx, Ywdx) max(0,Ygmtx - 0.125 .* (Ygmtx < (Ycdx + Ywdx))); 
  Ygmtw0 = pbtsulccor(Ygmtw0, Ycd0, Ywd0); 
  Ygmtc0 = pbtsulccor(Ygmtc0, Ycd0, Ywd0); 


  % minimum tickness map and cleanup (removal of extrem outliers and approximation) 
  % - not useful for Ygmtw0/Ygmtc0!
  Ygmt0  = min(Ygmtw0,Ygmtc0); 
  Ygmt0  = cleanupPBT(Ygmt0, 1, 0); % filter limits has only minor effects
  

  % update distance information
  Ycd0  = min(Ygmt0,Ycd0); Ywd0 = min(Ygmt0,Ywd0); 
 
  % this functions emphasize fine structurs, ie., to avoid blurring of small sulci
  if opt.keepdetails 
    [Ywd0, Ycd0, Ygmt0] = keepdetails(Yp0, Ywd0, Ycd0, Ygmt0, vx_vol, ...
      opt.range * opt.extendedrange, opt.keepdetails);
  end  

  if opt.verb, fprintf('\n    Thickness estimation:           %0.3fs\n', etime(clock,c)); c = clock; end %#ok<*DETIM>


  %% CSF/WM blurring/reconstruction maps
  % - define the blurred sulcal/CSF areas as the area of WM closing and 
  %   'overestimated' thickness (i.e., where the PBT thickness is sign. 
  %   smaller as the simple sum of the WM and CSF distance)
  % - blurred gyri/WM is defined vite versa
  % - undefined values GM values are defined by neighbours
  % - neutral regions are defined by CSF
  % - next both areas are extend by intensity emphasized values 
  if opt.gyrusrecon
    % define sulci/gyri as areas closed WM/CSF regions
    Ygsr   = cat_vol_morph( cat_vol_morph( Yp0<2.5 & cat_vol_morph(Yp0>2.5,'dc',10,vx_vol) & (Ygmtw0 < Ygmtc0) & (Ygmt0*1.05 <= Ycd0+Ywd0) , 'do',1.5), 'dd',1); % large sculci
    Ygsr   = max(Ygsr,cat_vol_morph( cat_vol_morph( Yp0<2.75 & cat_vol_morph(Yp0>2.75,'dc',3,vx_vol) & (Ygmtw0 < Ygmtc0) & (Ygmt0*1.05 <= Ycd0+Ywd0) , 'do',1), 'dd',2)); % small sulci
    Ygsr   = Ygsr*.5 + max(Yp0>2.5,cat_vol_morph(Yp0>=2.25 & ~Ygsr & cat_vol_morph(Yp0<2.25,'dc',10,vx_vol) & (Ygmtc0 < Ygmt0*1.1) & (Ygmt0*1.05 <= Ycd0+Ywd0),'do',1.5));
    % basic extension by distance function
    [~,I]  = cat_vbdist( single(Ygsr>.25),Yp0>1.1);  Ygsr = single(Ygsr(I)); clear I;
    
  
    % general (global) relation between sulci and gyri
    %  to avoid gyrus reconstruction in case of many sulci that comes with high risk of bridge defects
    gsr    = max(.7,min(1.3,nnz(Ygsr(:)==.5) ./ nnz(Ygsr(:)==1  & Yp0(:)<2.5))); 
    % in case of very thick cortices as in children (about 3 mm, e.g., BUSS## dataset) 
    % it is better to avoid the reconstruction whereas in case of low thickness 
    % (about 2 mm) the need is even higher!   
    % Besides thickness also high variance of thickness (> 0.75 mm) the risk 
    % of sulcul blurring increases stongly
    mdGMT  = median(Ygmt0(round(Yp0(:))==2)) * mean(vx_vol); 
    iqrGMT = iqr(Ygmt0(round(Yp0(:))==2)) * mean(vx_vol); 
    gsr    = max(0.3, min(1.7, gsr  .*  max(0.5,min(2,1 + (mdGMT - 2.5))) .* max(.5,min(2,(iqrGMT-.5))))); 
   
    % local refinement
    Ygsr(Yp0<1.5)                                         = max(0.55,min(0.80, 0.75  / gsr));      % extend blurred sulci
    Ygsr(Yp0<1.5 & ~cat_vol_morph(Yp0<1.5,'do',2,vx_vol)) = max(0.50,min(0.75, 0.625 / gsr.^.5));  % pro sulci
    Ygsr(cat_vol_morph(Yp0<1.5,'do',2,vx_vol))            = max(0.75,min(0.95, 0.875 / gsr.^.5));  % pro gyri
    %Ygsr  = 1 - max(1-Ygsr, max(0,min(1,2-Yp0)).^.75);            % limited emphasization, as thin gyri often have lower GM probability ... 
    Ygsr   =     max(Ygsr,   max(0,min(1,Yp0-2)).^.75);            % emphasize WM 
    Ygsr   = Yppsmooth(Ygsr,Ygmt0,vx_vol,[0,-1]);                  % outlier correction 
    Ygsr   = max(0,min(1,cat_vol_smooth3X(Ygsr,2) * 2 - 1));       % final smoothing  .^ (1.5 - Ygsr); ... too small cause edges and bridges (eg. Collins) 
  
    
    %% percentage blurred sulcal/gyral volume (regions that need reconstrution
    srecon = sum(Ygsr(:)<.25 & round(Yp0(:))==2 & (Ygmtw0(:) < Ygmt0(:)*1.1) & (Ygmt0(:)*1.05 <= Ycd0(:)+Ywd0(:)) ) / sum( round(Yp0(:))==2 ); 
    grecon = sum(Ygsr(:)>.75 & round(Yp0(:))==2 & (Ygmtc0(:) < Ygmt0(:)*1.1) & (Ygmt0(:)*1.05 <= Ycd0(:)+Ywd0(:)) ) / sum( round(Yp0(:))==2 ); 
  
    % position map 
    Yppg = min(1, max(  Ygmt0 .* (Yp0>2.5) , max( Ycd0 .* Ygsr.^.1, (Yp0>1.5) .* (Ygmt0-Ywd0) .* Ygsr.^2)) ./ Ygmt0); 
    Ypps = 0.5 .* min(1, max( 0 , max(Yp0>2.5, (Yp0>1.5) .* max(eps,Ygmtw0-Ywd0) ./ max(eps,Ygmtw0) .* Ygsr ))) + ... 
           0.5 .* min(1, max( 0 , max(Yp0>2.5, (Yp0>1.5) .* min(Ycd0,Ygmt0-Ywd0) ./ max(eps,Ygmt0))));
    Yppg = Yppg .* Ypps.^max(0.05, min(.5, .1 * (gsr.^4)));
    Ypp  = Yppg.*Ygsr + (1-Ygsr).*Ypps; 
  else
    gsr  = 0; srecon = 1; grecon = 0; 
    Ypp  = 0.5 .* min(1, max( 0 , max(Yp0>2.5, (Yp0>1.5) .* max(eps,Ygmtw0-Ywd0) ./ max(eps,Ygmtw0)))) + ... 
           0.5 .* min(1, max( 0 , max(Yp0>2.5, (Yp0>1.5) .* min(Ycd0,Ygmt0-Ywd0) ./ max(eps,Ygmt0))));
  end 

  % cleanup of position map
  % - a topology correction is not really useful as we might use different thresholds later (cost ~15 seconds) !
  Ypp(cat_vol_morph(Yp0>2.75,'l')) = 1; Ypp(Yp0<1.5) = 0; 
  Ypp = oneObject(Ypp,vx_vol); 
  Ypp = Yppsmooth(Ypp,Ygmt0,vx_vol,[0.5,-1]);  % correction of strong outliers (eg. from mixing)
  Ypp = max( -.1 , min( 1.1, Ypp + (Ypp>.95) .* max(0,Yp0 - 2.5) + (Ypp<.05) .* min(0,Yp0 - 1.5)  )); % PVE for interpolation and deformation 

  % final scaling
  Ygmt = Ygmt0 * mean(vx_vol); 

  % (*) EXTRA: magic sharpending function - enhancement of fine structures
  if opt.sharpening, [Ypp, Ygmt] = sharpening(Ypp, Ygmt, vx_vol, 0); end

  % evaluation
  if opt.verb
    fprintf('    PP preparation (gsr=%0.3f):     %0.3fs\n', gsr, etime(clock,c)); c = clock; 
    fprintf('    Sulcus / gyrus reconstruction:  %5.2f%% / %4.2f%%\n', srecon*100, grecon*100);  
    fprintf('    Median thickness + IQR:        %5.2f ± %4.2f mm\n', ...
      median( Ygmt( Ypp(:)>.4 & Ypp(:)<.6 )) , iqr( Ygmt( Ypp(:)>.3 & Ypp(:)<.7 ))); 
    x = 0:0.01:10;  
    h = smooth( hist( Ygmt( Ypp(:)>.45 & Ypp(:)<.55 ) , x),2); h = h/sum(h);
    try
      hi = find(h==max(h),1); hil = find(h==max(h(100:hi-30)),1); hih = find(h==max(h(hi+30:end)),1); 
      fprintf('    Peak (x:y):                    %5.2f:%4.4f | %5.2f:%4.4f | %5.2f:%4.4f\n', x(hil), h(hil), x(hi), h(hi), x(hih), h(hih));
    end
  end

  Yp0o = Yp0; 
  %% Update Yp0
  % Yp0=Yp0o; Yp0 = max(Yp0,(Ypp>.5) .* 2 + (max(0,min(1,0-Ywd0))));  Yp0 = min(Yp0,2*(Ypp>.5) + 1 + max(0,min(1,max(0,((Ygmt0 - Ywd0 + .75))).^4)));
  Yp0=Yp0o; Yp0 = max(Yp0,(Ypp>.5) .* 2 + (max(0,cat_vol_smooth3X(Ypp,.5)*2-1)).^4); Yp0 = min(Yp0,2*(Ypp>.5) + 2 - max(0,1 - cat_vol_smooth3X(Ypp,.5)*.5).^64); 
  %ds('d2sm','',1,abs(Yp0fs-Ymfs),abs(Yp0e - Ymfs),150)
  %ds('d2sm','',1,Yp0o/3,Yp0/3,140)


end
% ======================================================================
function Ypp = Yppsmooth(Ypp,Ygmt,vx_vol,th)
  if ~exist('th','var'), th = 0.1; end
  if isscalar(th), th = repmat(th,1,2); end

  % filter in large areas
  Yflt    = .5; %abs(smooth3(Ypp)*2 - 1) .* max(0,min(1,Ygmt - median(Ygmt(Ypp(:)>0 & Ypp(:)<1))/2 )); 
  if th(1) >= 0
    Ypp   = cat_vol_median3(Ypp, cat_vol_morph(Ypp>0 & Ypp<1,'d') , true(size(Ypp)), th(1));
    Ypp   = cat_vol_smooth3X(Ypp,0.5); 
    Ypp   = max(0,min(1, Ypp + .5 * (1-th(1))*(Ypp - smooth3(Ypp)) )); 
  end
  %Ypp    = Ypp.*(1-Yflt) + Yflt.*Yppm; 
  if th(2) >= 0
    [Yppr,resYp0] = cat_vol_resize(Ypp,'reduceV',vx_vol,1,32,'median');
    Yppr  = cat_vol_median3(Yppr, cat_vol_morph(Yppr>0 & Yppr<1,'d') , true(size(Yppr)), th(2)); 
    Yppr  = cat_vol_resize(Yppr,'dereduceV',resYp0); 
    Yppr  = max(0,min(1,  Yppr + .5 * (1-th(2))*(Yppr - smooth3(Yppr)) )); 
    Yppr  = max(0,min(1,tan(Yppr - .5) + .5)); 
    Ypp   = Ypp.*(1-Yflt) + Yflt.*Yppr;
    Ypp   = max(0,min(1, Ypp + .5 * (1-th(1))*(Ypp - smooth3(Ypp)) )); 
  end
end
% ======================================================================
function Ygmtc = cleanupPBT(Ygmt,lim,lim2)
% need resolution to handle filter size better !
  if ~exist('lim','var'), lim = .05; end
  if ~exist('lim2','var'), lim2 = 1; end
  Ygmta = cat_vol_approx(Ygmt,'rec') ; 
  if lim2 > 0
    Ygmtc = cat_vol_approx(Ygmt .* ((Ygmt - Ygmta .* (Ygmt>0))<=lim  &  (Ygmt - Ygmta .* (Ygmt>0))<=lim*4  &  Ygmt>0), 'rec'); 
    Ygmt  = Ygmt .* (Ygmt>0 & abs( Ygmt - Ygmtc ) < 1); 
    Ygmta = cat_vol_smooth3X( Ygmta .* (Ygmt==0) + Ygmt , .5);
  end
  if lim > 0
    Ygmtc  = cat_vol_approx(Ygmt .* ((Ygmt - Ygmta .* (Ygmt>0))<=lim  &  (Ygmt - Ygmta .* (Ygmt>0))<=lim*4  &  Ygmt>0), 'rec'); 
    Ygmtc  = cat_vol_smooth3X( min(Ygmta,Ygmtc),.5); % better
  else
    Ygmtc  = Ygmta; 
  end
end
% ======================================================================
function Yp0 = oneObject(Yp0, vx_vol, n)
  if ~exist('n','var'), n = 10; end 

  for i = 1:n
    if max(Yp0(:))>1.5 % Yp0 with values from 1 to 3
      th  = 1 + ((n - i) / n) * 2;
    else % Ypp with values from 0 to 1
      th  = 0 + ((n - i) / n) * 1;
    end
    
    if th > 0.5 + 1.5*(max(Yp0(:))>1.5)
      Yl  = single(cat_vol_morph( Yp0 > th ,'l')); 
    else % additional opening
      Yl  = single(cat_vol_morph( Yp0 > th ,'ldo',(3 - th)/2,vx_vol)); 
    end
    Yp0 = min(Yp0,th) + Yl .* max(0,Yp0-th);  
  end
end
% ======================================================================
function Yp0 = NBVC(Yp0,vx_vol)
%% RD202503: new blood vessel correction 
% In principle this is not really new and should be done by other functions before. 
% However, it is so essential for PBT and not to difficult to include it here (too). 
% First the save WM area is estimated and then extended by a region growing method
% that focus on lower intensities in the Yp0 label mapf for GM and CSF. 
  F      = max(0.0,Yp0-1); F(Yp0<=1.1) = inf; 
  Ywm    = cat_vol_morph(Yp0>2.5,'ldo',2,vx_vol) | cat_vol_morph(Yp0>2.75,'ldo',1); 
  [~,Yd] = cat_vol_downcut(single(Ywm), F,-0.001); Yd(isnan(Yd))=inf; clear Ywm; 
  Ymsk   = Yd > 1000000 & Yp0>2; 
 
  % Ymsk as regions that we want to correct 
  %[~,I]  = cat_vbdist(single(~Ymsk),Yp0>1.5); Yp0(Ymsk) = Yp0(I(Ymsk));
  Yp0(Ymsk) = 2; 
  Ymsk   = cat_vol_morph(Ymsk,'dd',1,vx_vol) & Yp0>1.5;
  Yp0s   = cat_vol_median3(Yp0,Ymsk);
  Yp0(Ymsk) = Yp0s(Ymsk); 
  
  % CSF enhancement (by WM distnace) for better estimation of asymmetry and suppression of artefact ?
  %Yd(Ymsk) = 0; 
  %Yp0    = max( min(1.5,Yp0) , Yp0 - max(0,Yd/10000 .* (Ymsk)) );
end
% ======================================================================
function [Ypp, Ygmt] = sharpening(Ypp, Ygmt, vx_vol, extended)
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
  %if ~extended
  %  Ypps = min(Ypp>0, max(0,Ypps)); % remove background artifacts  
  %end

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
function [Ycd,Ywd] = cat_vol_PVEdist(Yp0,PVErefinement)
% function to estimate the raw distance maps

  if PVErefinement == 0
  % simple distance estimation with two levels to be robust also in case 
  % of the 5-cls AMAP segmentation with a subclass at 1.5 and 2.5 (although
  % this is blurred in general by the interpolation)
    Ycd = cat_vbdist( single(Yp0 <  1.5), Yp0 < 3 )/2 + ...
          cat_vbdist( single(Yp0 <= 1.5), Yp0 < 3 )/2 - 0.5; 
    Ywd = cat_vbdist( single(Yp0 >  2.5), Yp0 > 1 )/2 + ...
          cat_vbdist( single(Yp0 >= 2.5), Yp0 > 1 )/2 - 0.5;
  elseif PVErefinement == 1
    % Ycci/Ywwi is correction for voxel behind the boundary to avoid 
    % overestimation by PBT mapping (using the maximum values)
    Ycci = max(0,cat_vbdist( single(Yp0 <  2.5), Yp0 < 3 )/2 + ...
                 cat_vbdist( single(Yp0 <= 2.5), Yp0 < 3 )/2 - .5); 
    Ywwi = max(0,cat_vbdist( single(Yp0 >  1.5), Yp0 > 1 )/2 + ...
                 cat_vbdist( single(Yp0 >= 1.5), Yp0 > 1 )/2 - .5); 
    Ycd = (cat_vbdist( single(Yp0 <  1.5), Yp0 < 3 ) - Ycci)/2 + ...
          (cat_vbdist( single(Yp0 <= 1.5), Yp0 < 3 ) - Ycci)/2 - .5; 
    Ywd = (cat_vbdist( single(Yp0 >  2.5), Yp0 > 1 ) - Ywwi)/2 + ...
          (cat_vbdist( single(Yp0 >= 2.5), Yp0 > 1 ) - Ywwi)/2 - .5;
  elseif PVErefinement == 3
    Ycci = max(0,cat_vbdist( single(Yp0 < 2.325), Yp0 < 3 )/2 + ...
                 cat_vbdist( single(Yp0 < 2.675), Yp0 < 3 )/2 - 1); 
    Ywwi = max(0,cat_vbdist( single(Yp0 > 1.325), Yp0 > 1 )/2 + ...
                 cat_vbdist( single(Yp0 > 1.675), Yp0 > 1 )/2 - 1); 
    Ycd = (cat_vbdist( single(Yp0 < 1.325), Yp0 < 3 ) - Ycci)/2 + ...
          (cat_vbdist( single(Yp0 < 1.675), Yp0 < 3 ) - Ycci)/2 - 1; 
    Ywd = (cat_vbdist( single(Yp0 > 2.325), Yp0 > 1 ) - Ywwi)/2 + ...
          (cat_vbdist( single(Yp0 > 2.675), Yp0 > 1 ) - Ywwi)/2 - 1;
  
  elseif PVErefinement == 2
    Ycd = cat_vbdist( single(Yp0 < 1.5), Yp0 < 2.5 ); 
    Ywd = cat_vbdist( single(Yp0 > 2.5), Yp0 > 1.5 );

  else
  % complex distance estimation with multiple levels (but not to close to 
  % the standard classes that are more prone to interpolation errors)

    %level = [0.25 0.75];
    level = [0.125 0.25 0.375   0.5 + [0.125 0.25 0.375]]; % in SPM no real improvent in Colins

    % distance maps for all levels 
    Ywda = zeros([size(Yp0),numel(level)],'single'); 
    Ycda = zeros([size(Yp0),numel(level)],'single'); 
    for li = 1:numel(level)
      % raw distance measure
      Ycda(:,:,:,li) = cat_vbdist( single(Yp0 < 1 + level(li)), Yp0 < 2.5); 
      Ywda(:,:,:,li) = cat_vbdist( single(Yp0 > 2 + level(li)), Yp0 > 1.5);
      
      % simple offset correction 
      % **** this could be maybe improved later as the level intensity 
      %      offset is only a rough approximation of the real distance 
      %      offest 
      % **** use/correct for negative values?
      Ycda(:,:,:,li) = Ycda(:,:,:,li) + (0.5 - level(li)) .* (Yp0>1 & Yp0<3); 
      Ywda(:,:,:,li) = Ywda(:,:,:,li) + (0.5 - level(li)) .* (Yp0>1 & Yp0<3); 
    end
    % simple average thickness
    Ycd0 = single(cat_stat_nanmean(Ycda,4));
    Ywd0 = single(cat_stat_nanmean(Ywda,4));
    
    % raw thickness maps that we use to select the best single distance maps
    Ygmtw0 = cat_vol_pbtp( round(Yp0)   , Ywd0, Ycd0); Ygmtw0(Ygmtw0>1000) = 0; 
    Ygmtc0 = cat_vol_pbtp( 4-round(Yp0) , Ycd0, Ywd0); Ygmtc0(Ygmtc0>1000) = 0; 
% ***** add cleanupPBT of Ygmtw0 and Ygmtc0 ?    
    Ygmt0  = min(Ygmtw0,Ygmtc0); 
    Ygmt0  = cleanupPBT(Ygmt0); 
    Ycd0c  = (Ygmt0 - Ywd0) .* (Yp0>1 & Yp0<3); 
    Ywd0c  = (Ygmt0 - Ycd0) .* (Yp0>1 & Yp0<3); 
    
    %% we assume a continuous thickness pattern without extrem values 
    mdgmt   = median(Ygmt0(round(Yp0(:))==2)); 
    liqrgmt = prctile(Ygmt0(round(Yp0(:))==2),20) - mdgmt; % thinner values are expected (sulci) and ok  
    hiqrgmt = prctile(Ygmt0(round(Yp0(:))==2),60) - mdgmt; % larger values are more often outliers
% **** maybe the recon maps could be handy here?

    %% estimation of the deviation of the diance from the local expected value 
    Ywde = Ywd0 * 0; 
    for li = 1:numel(level)
      % estimate thickness - avg. exp. thickness 
      %  - pos. value means too thin, i.e., good in gyri bad in sulci
      %  - neg. value means too thick, i.e., 
      Ywdel = (Ycd0c + Ywda(:,:,:,li)) - Ygmt0;
      Ywdel(Ywdel > 0  &  Ywdel < liqrgmt) = max(0,Ywdel(Ywdel > 0  &  Ywdel < liqrgmt) + liqrgmt); 
      Ywdel(Ywdel < 0  &  Ywdel > hiqrgmt) = min(0,Ywdel(Ywdel < 0  &  Ywdel > hiqrgmt) + hiqrgmt); 
      Ywde(:,:,:,li) = Ywdel; 
    end

    % estimate and correct for minimum error 
    Ywdem = min(Ywde,[],4);
    for li = 1:numel(level)
      Ywde(:,:,:,li)  = cat_vol_smooth3X( Ywde(:,:,:,li) - Ywdem , 2 ); 
    end
    clear Ywdem

    % include only good distance estimates
    Ywd = Ywd0*0; Ywdn = Ywd; 
    for li = 1:numel(level)
      Ywdn = Ywdn + Ywde(:,:,:,li); 
      Ywd  = Ywd  + Ywde(:,:,:,li) .* Ywda(:,:,:,li); 
    end
    Ywd = Ywd ./ max(eps,Ywdn); 
    clear Ywde; 
    

    %% estimation of the deviation of the diance from the local expected value 
    Ycde = Ycd0 * 0; 
    for li = 1:numel(level)
      Ycde(:,:,:,li) = abs( Ygmt0 - (Ywd0c + Ycda(:,:,:,li)) );
      Ycdel = (Ywd0c + Ycda(:,:,:,li)) - Ygmt0;
      Ycdel(Ycdel > 0  &  Ycdel < liqrgmt) = max(0,Ycdel(Ycdel > 0  &  Ycdel < liqrgmt) + liqrgmt); 
      Ycdel(Ycdel < 0  &  Ycdel > hiqrgmt) = min(0,Ycdel(Ycdel < 0  &  Ycdel > hiqrgmt) + hiqrgmt); 
      Ycde(:,:,:,li) = Ycdel; 
    end

    % estimate the minimum error include only good distance estimates
    Ycdem = min(Ycde,[],4); 
    for li = 1:numel(level)
      Ycde(:,:,:,li)  = cat_vol_smooth3X( Ycde(:,:,:,li) - Ycdem , 2 ); 
    end
    clear Ycdem; 

    Ycd = Ycd0*0; Ycdn = Ycd;  
    for li = 1:numel(level)
      Ycdn = Ycdn + Ycde(:,:,:,li); 
      Ycd  = Ycd  + Ycde(:,:,:,li) .* Ycda(:,:,:,li); 
    end
    clear Ycde; 
    Ycd = Ycd ./ max(eps,Ycdn); 
    
  end

  % final correction for the voxel distance error (simplified half grid distance)
  % .. thicky ... should be .5 but this results in worse results
 % Ycd = max(0,Ycd - 1);
 % Ywd = max(0,Ywd - 1);
end
% ======================================================================
function [Ycd, Ywd] = cat_vol_cwdist(Yp0,opt,vx_vol)
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
    YMM = cat_vol_morph(Yp0 > 2.5 + opt.extendedrange,'e',1) | isnan(Yp0); 
    YMC = cat_vol_morph(Yp0 < 1.5 - opt.extendedrange,'e',1) | isnan(Yp0); 
  
    % estimation of the extend range to correct values beyond the other tissue boundary
    if opt.extendedrange > 0
      Ycdlc = Ycd; Ycdhc = Ycd; 
      Ywdlc = Ywd; Ywdhc = Ywd; 
      hss = opt.levels; % number of opt.levels (as pairs)
      for si = 1:hss
        offset = max(0,min(opt.range, opt.range * si/(hss+1))); 

        % CSF distance correction beyond WM boudary
        if opt.eidist == 3
          % Eikonal-based Euclidean distance
          F          = double( max(0,min(1, 3 - Yp0)) ); 
          [~,Ycdlct] = cat_vol_eidist3(Yp0 < 2.5 - offset, F); Ycdlc = Ycdlc + 1/hss * max(0,Ycdlct - .5);
          [~,Ycdhct] = cat_vol_eidist3(Yp0 < 2.5 + offset, F); Ycdhc = Ycdhc + 1/hss * max(0,Ycdhct - .5);
        elseif opt.eidist == 2
          % Eikonal-based Euclidean distance
          F      = max(eps,min(1, 3 - Yp0 )); 
          YM     = max(0,min(1,(3 - Yp0 - offset))); YM(YMM) = nan; Ycdlc = Ycdlc + 1/hss * max(0, cat_vol_eidist(YM,F,[1 1 1],1,1,0,0) - .5); 
          YM     = max(0,min(1,(3 - Yp0 + offset))); YM(YMM) = nan; Ycdhc = Ycdhc + 1/hss * max(0, cat_vol_eidist(YM,F,[1 1 1],1,1,0,0) - .5); 
        else
          % simple Euclidean distance 
          Ycdlc  = Ycdlc + 1/hss * max(0,cat_vbdist(single(Yp0 < ( 2.5 - offset)), ~YMM ) - .5);
          Ycdhc  = Ycdhc + 1/hss * max(0,cat_vbdist(single(Yp0 < ( 2.5 + offset)), ~YMM ) - .5); 
         
        end
        

        % WM distance correction beyond CSF boudary
        if opt.eidist == 3
          %%
          F          = double( max(0,min(1, Yp0 - 1)) ); 
          [~,Ywdlct]  = cat_vol_eidist3(Yp0 > 1.5 - offset,F); Ycdlc = Ycdlc + 1/hss * max(0,Ywdlct - .5);
          [~,Ywdhct]  = cat_vol_eidist3(Yp0 > 1.5 + offset,F); Ywdhc = Ywdhc + 1/hss * max(0,Ywdhct - .5);
        elseif opt.eidist == 2
          %%
          F      = max(eps,min(1, Yp0 - 1 )); 
          YM     = max(0,min(1,(Yp0 - 1 - offset))); YM(YMC) = nan; Ywdlc = Ywdlc + 1/hss * max(0,cat_vol_eidist(YM,F,[1 1 1],1,1,0,0) - .5); 
          YM     = max(0,min(1,(Yp0 - 1 + offset))); YM(YMC) = nan; Ywdhc = Ywdhc + 1/hss * max(0,cat_vol_eidist(YM,F,[1 1 1],1,1,0,0) - .5); 
        else
          Ywdlc  = Ywdlc + 1/hss * max(0,cat_vbdist(single(Yp0 > ( 1.5 - offset)), ~YMC) - .5);  
          Ywdhc  = Ywdhc + 1/hss * max(0,cat_vbdist(single(Yp0 > ( 1.5 + offset)), ~YMC) - .5); 
        end
      end
      
      Ywdlc(Ywdlc > 1000) = 0; Ywdhc(Ywdhc > 1000) = 0;  
      Ycdlc(Ycdlc > 1000) = 0; Ycdhc(Ycdhc > 1000) = 0; 

    end


    % multi-level distance estimation
    hss = opt.levels; % number of opt.levels (as pairs)
    for si = 1:hss
      offset = max(0,min(opt.range, opt.range * si/(hss+1))); 

      % CSF dist
      if opt.eidist == 3
        F        = max(0,min(1,((3 - Yp0)))); 
        [~,Ycdl] = cat_vol_eidist3(Yp0<1.5 - offset,double(F)); Ycdl = max(0, Ycdl - .5); 
        [~,Ycdh] = cat_vol_eidist3(Yp0<1.5 + offset,double(F)); Ycdh = max(0, Ycdh - .5); 
      elseif opt.eidist > 0 
        F     = max(eps,min(1,((3 - Yp0)))); 
        YM    = max(0,min(1,(2 - Yp0 - offset))); YM(YMM) = nan; Ycdl = max(0, cat_vol_eidist(YM,F,[1 1 1],1,1,0,0) - .5); 
        YM    = max(0,min(1,(2 - Yp0 + offset))); YM(YMM) = nan; Ycdh = max(0, cat_vol_eidist(YM,F,[1 1 1],1,1,0,0) - .5); 
      else
        Ycdl  = max(0,cat_vbdist(single(Yp0 < ( 1.5 - offset)), ~YMM ) -.5); Ycdl(Ycdl > 1000) = 0; 
        Ycdh  = max(0,cat_vbdist(single(Yp0 < ( 1.5 + offset)), ~YMM ) -.5); Ycdh(Ycdh > 1000) = 0;
      end

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
      % idea was to could the boundaries different depending on the CSF distance
      clear Ycdl Ycdh; 


      % WM distances
      if opt.eidist == 3
        F      = max(0,min(1,((Yp0-1 ))));
        [~,Ywdl]   = cat_vol_eidist3(Yp0>2.5 - offset,double(F));
        [~,Ywdh]   = cat_vol_eidist3(Yp0>2.5 + offset,double(F));
      elseif opt.eidist > 0
        F     = max(0,min(1,((Yp0-1 ))));
        YM    = max(0,min(1,(Yp0 - 2 - offset))); YM(YMM) = nan; Ywdl = max(0, cat_vol_eidist(YM,F,[1 1 1],1,1,0,0) -.5); 
        YM    = max(0,min(1,(Yp0 - 2 + offset))); YM(YMM) = nan; Ywdh = max(0, cat_vol_eidist(YM,F,[1 1 1],1,1,0,0) -.5); 
      else
        Ywdl  = max(0,cat_vbdist(single(Yp0 > ( 2.5 - offset)), Yp0 > 1.5 - opt.extendedrange ) -.5); Ywdl(Ywdl > 1000) = 0; 
        Ywdh  = max(0,cat_vbdist(single(Yp0 > ( 2.5 + offset)), Yp0 > 1.5 - opt.extendedrange ) -.5); Ywdh(Ywdh > 1000) = 0; 
      end

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
function Yp0 = myelincorrection(Yp0,vx_vol,opt)
  if opt.myelinCorrection > 0
    % quick estimation of the cortical thickness
    opt.verb   = 0; 
    opt.levels = 1; 
    opt.eidist = 0; 
    [Ycd, Ywd] = cat_vol_cwdist(Yp0, opt, vx_vol);
  
    % projection-based thickness mapping
    Ygmt0 = cat_vol_pbtp( round(Yp0) , Ywd, Ycd);
    Ygmt0 = cat_vol_approx(Ygmt0); 
  
    % reestimation of the CSF distance 
    Ypp   = min(1,min(Ygmt0,Ycd) ./ max(eps,Ygmt0)); Ypp(Yp0>2.5 & Ypp==0) = 1; 
    Ycdc2 = cat_vbdist( single( max(Yp0<=1, 1 - Ycd - Ypp) ), true(size(Ycd)) );
    Ycdc2(Ycdc2 > 6 / mean(vx_vol)) = 0; 
   
    % estimate the full tissue thickness (we needed the GM thickness and WM to reconstruct the sulcus)
    Ybmt  = cat_vol_pbtp( min(3,4 - min(2,Yp0)), Ycdc2, Ycdc2*inf); 
    Ybmt  = cat_vol_approx(Ybmt); 

    % estimate correction area
    medgmt = median(Ygmt0(:)); 
    try iqrgmt = iqr(Ygmt0(:)); catch, iqrgmt = std(Ygmt0(:)); end
    YenoughWM  = Ycdc2 < Ybmt - 1.5;
    YthinnerGM = max(0,medgmt - 1.5*iqrgmt - Ygmt0); 
    Yclose2CSF = Ycdc2>0 & Ycdc2<(medgmt - 1.5*iqrgmt); 
    Ygmwmpve   = cat_vol_morph(Yp0>2 & Yp0<2.9,'do',1); % | smooth3(Yp0>2 & Yp0<2.9)>.7);
    Ycor = YenoughWM & YthinnerGM & Yclose2CSF & Ygmwmpve;
    clear YenoughWM YthinnerGM Yclose2CSF Ygmwmpve;

    Yp0  = max(min(Yp0,2),max(Yp0>=2.95,Yp0 - smooth3( Ycor ) * opt.myelinCorrection));
    clear Ycdc2 Ybmt Ycor; 
  end
end
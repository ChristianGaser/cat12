function [Ygmt,Ypp,Yp0] = cat_vol_pbtsimpleCS4(Yp0,vx_vol,opt)
%cat_vol_pbtsimple. Simple cortical thickness/position estimation.  
% Estimation of voxel-based distances and projection-based thickness (PBT) 
% and surface position based on a label map. Required isotropic input. 
% 
% After further optimisation the function becomes more complex again. 
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
%    .levels (integer; default=8)
%      Number of dual distance estimations to reduce sampling effects.
%      With logarithmic improvement and good results between 2 and 16. 
%
%    .extendedrange (0-no, 1-yes; default=1)
%      Uses also the PVE range to estimate the distance there. 
%      Important to avoid thickness underestimations.
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
  def.levels          = 4;    % Number of dual distance estimates. 
                              % Larger values are more accurate but need more time (log-change).
                              % 1 means estimation for .25 and .75 boundaries
                              % Good values are between 1 and 8.
                              % RD202601: less was partially better in phantoms
                              %  
  def.extendedrange   = 0;    % Estimate the distance from boundary A to boundary B 
                              % also for voxels beyond B with a correction of the 
                              % extra distance, to stabilize the values and avoid
                              % underestimations (0-no, 1-yes)
                              % RD2026: difficult
                              %
  def.range           = 0.4;  % Default value for range extension (should be between 0.2-0.4)
                              %
  def.keepdetails     = 2;    % enhance thin (occipital) sulci (and gyri) to avoid blurring 
                              % (0-no, 1-yes (sulci); 2-yes(sulci+gyri) 
                              % worse values but pretty important!
                              %
  def.NBVC             = 1;   % new blood vessel correction (RD202503)
                              %
  def.myelinCorrection = .3;  % correction for large cortical myelination artefacts that 
                              % cause strong local underestimation that is similar to the 
                              % extended LAS correction (cat_main_correctmyelination) (RD202503)
                              %
  def.gyrusrecon       = 1;   % use PBT also to reconstruct gyri 
                              % key aspect for better surfaces but also risk for blurred sulci
                              %
  def.verb             = 0;   % be verbose
                              %
  def.smoothPVE        = 2;   % smooth and sharpen partial volume effects
                              % 0-none, 1-median, 2-median+sharpen
                              %
  def.usemedian        = 1;   % GMT filtering setting
                              % 
  def.corpve           = 3;   % correct PVE effect
                              % 1-pve*, 2-dist, 3-both* ... + not working

  opt = cat_io_checkinopt(opt,def);


  if 1 %opt.extendedrange > 0 
    % extend hard cuts, by a low value just to have a broader estimate
    % - could be improved by a basic WMD and GMT estimate to assure a
    %   local minimum thickness related to the GMT distance
    % - RD202601: important for small distances but less for large
    % - Ok for phantom.
    Yp0 = max(Yp0,(Yp0==1 & cat_vol_morph(Yp0==2,'d')) * 1.2);
  end


  % close holes 
  % Helpful for SPM with unsufficient WM correction but minor/no-effects for CAT.
  % Ok for phantom.
  if 1
    Yp0 = max(Yp0, 2.75 * smooth3(cat_vol_morph(Yp0>2.75,'ldc',1.5)) ); 
    Yp0 = max(Yp0, 2.50 * smooth3(cat_vol_morph(Yp0>2.25,'ldc',1.5)) ); 
  end


  % Estimate the size/thickness of the partial volume effect area to
  % Scale the partial volume area to use the intensity information to
  % approximate the boundary correction.
  % The fast version (second parameter = 1) could cause more issues. 
  opt.pvet = estimatePVEsize( Yp0 , 0); 


  % corrections for blood vessels, partial volume effects and myelination
  % RD202503: new blood vessel correction 
  if opt.NBVC, Yp0 = NBVC(Yp0, vx_vol); end 
  
  % RD202601: use median filter to reduce (interpolation/stair) artifacts!
  Yp0 = cleanPVE(Yp0, opt.smoothPVE); % 0-no, 1-median, 2+sharpen
  
  % RD2025: correct myelinated areas (using quick thickness estimation)
  if opt.myelinCorrection, Yp0 = myelincorrection(Yp0, vx_vol, opt); end
  c = clock; %#ok<*CLOCK>


  

  %% == Estimation of distance measures ==
  % It is important to use the PVE for multiple (=levels) distance estimations 
  % but also to extened (range) and correct (corpve) the values at the oposite  
  % boundary by the PVE or distance. 
  tic 
  if 0
    % quick manuell settings
    opt.levels = 4; 
    opt.corpve = 3;    % 1-pve*, 2-dist, 3-both* ... + not working
    opt.range  = 0.45; %0.125 / 2; % addition range (half distance between (PVE) classes, i.e., .25, what is optimal; 
    opt.range  = min(.4,opt.range);
  end
  % tested here multple estimations
  switch 1
    case 1 % current default
      opt.extendedrange = 0;  
      [Ycd0, Ywd0] = cat_vol_cwdist( Yp0 ,opt);;;
    case 2 
    % with basic extention and filtering .. improvement after mapping
    % but their are odd values outside 
      opt.extendedrange = 0.5;    
      [Ycd0, Ywd0] = cat_vol_cwdist_old(Yp0,opt);
    case 3     
      opt.extendedrange = 1;  
      [Ycd0, Ywd0] = cat_vol_cwdist_old(Yp0,opt);
    case 4  %%
      opt.extendedrange = 2;  
      [Ycd0, Ywd0] = cat_vol_cwdist_old(Yp0,opt);
  end
  toc, disttoc = toc; 

if 0
  % RD202602: DELETE ME LATER
  % just a stupid test if the round thickness levels phantom might prefere 
  % the less accurate measurements
  Ycd0 = round(Ycd0*2)/2;
  Ywd0 = round(Ywd0*2)/2;
end


  % == Estimation of thickness maps == 
  % Adjustment for the unnderestimation in case of the half the PBT maps with tighter boundary.
  % For PBT we have to consider that we slighly underestimate the ribbon with a correct distance value.
  % plusx = @(Y) Y + (.5 * (opt.extendedrange==.5)).*(Y>0); 
  tic; plusx = @(Y) Y + .5 .* (Yp0>1.5-opt.range & Yp0<2.5+opt.range);
  Ygmtw0 = cat_vol_pbtp( single(1 + (Yp0>1.5-opt.range) + (Yp0>2.5+opt.range)) , plusx(Ywd0), plusx(Ycd0));
  Ygmtc0 = cat_vol_pbtp( single(3 - (Yp0>1.5-opt.range) - (Yp0>2.5+opt.range)) , plusx(Ycd0), plusx(Ywd0));
  
  % basic thickness corrections
  Ygmtw0(Ygmtw0>1000 | Ygmtw0<0) = 0; Ygmtw0 = min(Ygmtw0,max(0,Ycd0+Ywd0));
  Ygmtc0(Ygmtc0>1000 | Ygmtc0<0) = 0; Ygmtc0 = min(Ygmtc0,max(0,Ycd0+Ywd0));

  % basic outlier correction 
  if 1
    Ygmtw0 = outierfilter(Ygmtw0,Ywd0,Ycd0);
    Ygmtc0 = outierfilter(Ygmtc0,Ycd0,Ywd0);
  else
    % in principle thickness corrections need also to adapt the distance maps
    [Ygmtw0,Ywd0,Ycd0] = outierfilter(Ygmtw0,Ywd0,Ycd0);
    [Ygmtc0,Ycd0,Ywd0] = outierfilter(Ygmtc0,Ycd0,Ywd0);
  end
  toc, pbttoc = toc; 

  Ymsk = Ycd0<=0 | Ywd0<=0 | Ygmtw0<=0 | Ygmtc0<=0; 
  Ycd0(Ymsk) = 0; Ywd0(Ymsk) = 0; Ygmtw0(Ymsk) = 0; Ygmtc0(Ymsk) = 0; 
  

  % Internal tests
  if 0
    %% median filter
    %  red > 2 create artefacts but this can be used as a feature to see how good the raw data is
    tic; Ygmtw0a = min(Ygmtw0,Ygmtc0); %Ygmtw0a = mixgmt(Ygmtw0,Ygmtc0,Ywd0+Ycd0,.5); 
    Ygmtw0as = Ygmtw0a; red = [4]; fs = [1]; ol=0; 
    for redi = red, Ygmtw0as = medfilter(Ygmtw0as, Ygmtc0, Ywd0, Ycd0, Yp0, redi, fs, ol, 0); end;;;
    
    if 1
      %% create histogram
      %  * we need here a very fine sampling (many bins) to see artifacts 
      binsfac = 10; 
      Ymsk = Yp0>1.5 & Yp0<2.5 & Ycd0>0 & Ywd0>0 & Ygmtw0<10 & Ygmtw0>0; 
      Ypp  = min(1,max(0, min(Ygmtw0-Ywd0,Ywd0) .* (Ygmtw0>1)) ./ max(eps,Ygmtw0/2));
      Ypp  = min(1,max(0, (Ypp ./ cat_vol_approx(Ypp .* Ypp>.5)).^2 * 1.5 ));
      %
      data = {}; dataname = {}; datacolor = []; 
      %data{end+1} = Ycd0(Ypp(:)>.01)/2 + Ywd0(Ypp(:)>.01)/2; dataname{end+1} = 'cd0+wd0 (GM)'; datacolor(end+1,:) = [1 0.6 1];
      data{end+1} = Ycd0(Ypp(:)>.56)/2 + Ywd0(Ypp(:)>.56)/2; dataname{end+1} = 'cd0+wd0 (CL)'; datacolor(end+1,:) = [1 0.7 1];
      data{end+1} = Ygmtw0a(Ypp(:)>.54)/2; dataname{end+1} = 'gmtw0 (CL)';  datacolor(end+1,:) = [0 0.4 0.7];
      %data{end+1} = Ygmtw0(Ypp(:)<.74 & Ypp(:)>.01)/2; dataname{end+1} = 'gmtw0 (OL)';  datacolor(end+1,:) = [0.1 0.6 0];
      %data{end+1} = Ygmtc0(Ypp(:)>.73)/2; dataname{end+1} = 'gmtc0 (CL)';  datacolor(end+1,:) = [0.7 0.3 0];
      data{end+1} = Ygmtw0as(Ypp(:)>.53)/2; dataname{end+1} = 'gmtw0a (CL)'; datacolor(end+1,:) = [.8 0 0];
      %data{end+1} = Ygmt0(Ypp(:)>.73)/2; dataname{end+1} = 'gmt0 (CL)';  datacolor(end+1,:) = [0.1 0.6 0];
      %data{end+1} = Ygmt0m(Ypp(:)>.72)/2; dataname{end+1} = 'gmt0m (CL)';  datacolor(end+1,:) = [1 0 0];
      cat_plot_histogram(data,struct('color',datacolor,'bins',100*binsfac));
      xlim([-.1 round(prctile(Ygmtw0(Ygmtw0(:)>0)/2,90)+1)]); ylim([0 .5/binsfac]); grid on; 
      ah = gcf; ah.Position(3:4) = [300,200]*1.5;
      title( [char(datetime) sprintf(' - ex/pve/r=%0.1f/%0.0f/%0.2f - D/T/S=%0.2fs/%0.2fs/%0.2fs',...
        opt.extendedrange,opt.corpve,opt.range,  disttoc,pbttoc,toc)]); 
      legend(dataname,'Location','Northwest'); 
    end
    if 0
      %% create surface
      Ygmts = Ygmtw0as; 
      %Ygmts = Ygmt0; 
      Ygmts = smooth3( Ygmts + cat_vol_approx(Ygmts) .* (Ygmts==0) );
      Ypp   = max(Yp0>2.5, max(0,(Ygmts/2 .* (Ywd0>0)) - Ywd0/2) ./ max(eps,Ygmts)); % basic map without mixing
      CS    = isosurface(smooth3(1-Ypp),.5,min(5,Ygmts/2)); cat_surf_render2(CS); 
      title( [char(datetime) sprintf(' - ex/pve/r=%0.1f/%0.0f/%0.2f - D/T/S=%0.2fs/%0.2fs/%0.2fs',...
        opt.extendedrange,opt.corpve,opt.range,  disttoc,pbttoc,toc)]);
      %if median(Ygmts(Ygmts(:))) < 3.5, cat_surf_render2('clim',[1 4]); end 
    end
  end

  %% median filter
  % red and fs can pre define as matrixes to be used as loop
  % red:   1 is very slow, 2 works quite well (fast and save), 
  %        >3 - allow strong filtering but could causse problems (multiple peaks)
  % fs:   internal filter size with fixed loops
  % ol:   outlier correction (of high values)
  %       2 was ok, but stronger correction biased the results
  % good combinations for the phantom where  red = [2 2]; fs = [2 1]; ol=2;
  if opt.usemedian
    red = [2]; fs = [3]; ol=0; % ol>=2! 
  
    Ygmtw0m = Ygmtw0; for redi = red, Ygmtw0m = medfilter(Ygmtw0m, Ygmtc0, Ywd0, Ycd0, Yp0, redi, fs, ol, 1); end
    Ygmtc0m = Ygmtc0; for redi = red, Ygmtc0m = medfilter(Ygmtc0m, Ygmtw0, Ywd0, Ycd0, Yp0, redi, fs, ol, 1); end
    
    if 0
    % also here we would have to adapt the distance maps too
      Ygmtmd = Ygmtw0 - Ygmtw0m.*(Ygmtw0>0); 
      Ywd0   = Ywd0 + max(0,Ygmtmd)/4; Ywd0 = min(Ywd0,Ygmtmd); 
      Ycd0   = Ycd0 + max(0,Ygmtmd)/4; Ycd0 = min(Ycd0,Ygmtmd); 
  
      Ygmtmd = Ygmtc0 - Ygmtc0m.*(Ygmtc0>0); 
      Ywd0   = Ywd0 + max(0,Ygmtmd)/4; Ywd0 = min(Ywd0,Ygmtmd); 
      Ycd0   = Ycd0 + max(0,Ygmtmd)/4; Ycd0 = min(Ycd0,Ygmtmd); 
    end

    Ygmtw0 = Ygmtw0m; 
    Ygmtc0 = Ygmtc0m; 
  end

  % only within hull
  [Yp0r,resYp0] = cat_vol_resize(Yp0,'reduceV',vx_vol,1,32,'meanm'); % 1mm
  [Ywr,resYp0w] = cat_vol_resize( cat_vol_morph( Yp0r>1.5 ,'de',2,resYp0.vx_volr) ,'reduceV',resYp0.vx_volr,2,32,'meanm'); clear Yp0r;
  Ywr = cat_vol_morph(Ywr>.5,'ldc',8);
  Ywr = cat_vol_smooth3X(Ywr,2); 
  Ywr = cat_vol_resize(Ywr,'dereduceV',resYp0w); 
  Yhull = cat_vol_resize(Ywr,'dereduceV',resYp0); clear Ywr; 
  Yhull = max(0,min(1,atan(Yhull*pi/2))); 


  % avoid blurring of small sulci ... this is critical but works quite well
  if 1 % RD202601 ok for phantom? - not really but could be worse, however it helps in creal data
    
    % This is critical for the phantom ! 
    % But also for small/blurring sulci (see msk defintion) in real data Defined on HR075. 
    if 1
      fac          = .5; 
      Ymsk         = Yhull > .5  &  Ygmtw0 < Ygmtc0  &  Ygmtw0 < 4  &  cat_vol_smooth3X(Yp0,2)>=2; 
      Ygmtc0(Ymsk) = max(0,Ygmtw0(Ymsk) - .25*fac); 
      Ygmtw0(Ymsk) = max(0,Ygmtw0(Ymsk) - .50*fac); 
      Ywd0(Ymsk)   = Ywd0(Ymsk) + .5*fac; 
      Ycd0(Ymsk)   = min(Ycd0(Ymsk), max(0,Ygmtw0(Ymsk) - Ywd0(Ymsk)));  
      Ygmtw0       = min(Ygmtw0,Ycd0+Ywd0); Ygmtc0 = min(Ygmtc0,Ycd0+Ywd0); 
    end
   
    % save sucli that would blur 
    if 1
      Ymsk       = Yhull > .5  &  Ygmtw0 < Ygmtc0  &  Yp0 < cat_vol_smooth3X(Yp0,2)  &  Yp0>=2  &  cat_vol_smooth3X(Yp0,2)>2.1  &  Ycd0>4;
      Ycd0(Ymsk) = min( Ycd0(Ymsk), max(0,Ygmtw0(Ymsk) - Ywd0(Ymsk))); 
      Ygmtw0     = min(Ygmtw0,Ycd0+Ywd0); Ygmtc0 = min(Ygmtc0,Ycd0+Ywd0); 
    end

    if 1
      % save gyri that would blur ... could cause sulcal blurring ! 
      Ymsk       = Yhull > .5  &  Ygmtc0 < Ygmtw0 - .5  &  Yp0 > cat_vol_smooth3X(Yp0,2)  & smooth3(Yp0)>2  &  Ywd0>4  &  Ygmtc0<4  & Ycd0<4;
      Ywd0(Ymsk) = min(Ywd0(Ymsk), max(0,Ygmtc0(Ymsk) - Ycd0(Ymsk)));
      Ygmtw0     = min(Ygmtw0,Ycd0+Ywd0); Ygmtc0 = min(Ygmtc0,Ycd0+Ywd0); 
    end

  end



  % minimum tickness map and cleanup (removal of extrem outliers and approximation) 
  [Ygmt0,Ywd0,Ycd0] = outierfilter(min(Ygmtw0,Ygmtc0), Ywd0, Ycd0);
  if opt.usemedian 
    red = [2]; fs = [2]; ol=2;
    Ygmt0m = Ygmt0; for redi = red, Ygmt0m = medfilter(Ygmt0m, Ygmt0m, Ywd0, Ycd0, Yp0, redi, fs, ol, 1); end
    Ygmt0m = smooth3( Ygmt0m  +  cat_vol_approx(Ygmt0m) .* (Ygmt0m==0) ); % rm outliers
    if 0
      % adjust distance measures for correction - minimal effect?
      Ygmtmd = Ygmt0 - Ygmt0m.*(Ygmt0>0); 
      Ywd0   = Ywd0 + max(0,Ygmtmd)/4;
      Ycd0   = Ycd0 + max(0,Ygmtmd)/4;
      clear Ygmtmd
    end
    Ygmt0 = Ygmt0m; 
  else
    Ygmt0 = smooth3(max(Ygmt0, cleanupPBT(Ygmt0 .* (Ygmt0>0), 2, 0) .* (Ygmt0<0.5) ));
  end
  % update distance information
  Ycd0 = min(Ygmt0,Ycd0); Ywd0 = min(Ygmt0,Ywd0); % limit
  clear Ycmsk;


  % final adjustment based no the phantom ?
  if 1
    cf = 0.1; shift = +.1;
    Ycd0   = max(0,Ycd0 + cf/2.*(Ycd0>0) + shift.*(Ycd0>0));
    Ywd0   = max(0,Ywd0 + cf/2.*(Ywd0>0) - shift.*(Ycd0>0));
    Ygmt0  = max(0,Ygmt0  - cf); 
    Ygmtw0 = max(0,Ygmtw0 - cf); 
    Ygmtc0 = max(0,Ygmtc0 - cf); 
  end

  

  % this functions emphasize fine structurs, ie., to avoid blurring of small sulci
  if 1 %opt.keepdetails 
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
    if opt.gyrusrecon == 1 
    % complex version ... 
    
      % define sulci/gyri as areas closed WM/CSF regions (1=sulcus, 0=gyrus)
      Ygsr   = cat_vol_morph( cat_vol_morph( Yp0<2.5 & cat_vol_morph(Yp0>2.5,'dc',10,vx_vol) & ...
               (Ygmtw0 < Ygmtc0) & (Ygmt0*1.05 <= Ycd0+Ywd0) , 'do',1.5), 'dd',1); % large atrophic sculci
      Ygsr   = max(Ygsr, cat_vol_morph( cat_vol_morph( Yp0<2.75 & cat_vol_morph(Yp0>2.75,'dc',3,vx_vol) & ...
               (Ygmtw0 < Ygmtc0) & (Ygmt0*1.05 <= Ycd0+Ywd0) , 'do',1), 'dd',2)); % small sulci
      Ygsr   = Ygsr*.5 + max(Yp0>2.5, cat_vol_morph(Yp0>=2.25 & ~Ygsr & cat_vol_morph(Yp0<2.25,'dc',10,vx_vol) & ...
               (Ygmtc0 < Ygmt0*1.1) & (Ygmt0*1.05 <= Ycd0+Ywd0),'do',1.5));
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

      Ymsk3  = max(Yp0>2.5,min(1,max(0,min(Ycd0,Ygmt0-Ywd0)) ./ Ygmt0)).^(0.25 ./ (Ygmt0/median(Ygmt0(:)))); 
      Ygsr   = min(Ygsr, Ymsk3); 
      Ygsr   = min(Ygsr, Yhull); 
      
    else 
    % simpler version based on the minimum thickness selector
      
      % local gyrus-sulcus definition 
      [~,Ygmt0I] = min(cat(4,Ygmtw0,Ygmtc0),[],4);
      Ymsk1  = max(0,min(1,(Ycd0>1/vx_vol) .* Ycd0./max(eps,Ygmtw0) .* (5 ./ (Ygmt0*mean(vx_vol))) + (Yp0>2.5))); % use the uncorrected maps to outline the sulcus .. not working
      Ymsk2  = max(0,cat_vol_smooth3X(Yp0-1,2)); % up-weight WM and down-weight CSF regions
      Ymsk3  = max(Yp0>2.5,min(1,max(0,min(Ycd0,Ygmt0-Ywd0)) ./ Ygmt0)).^.1; 

      Ygsr   = max(0,min(1,cat_vol_smooth3X(Ygmt0I .* Ymsk1 .* Ymsk2 .* Ymsk3, 2 ) - 1)) .^ 1.25; % .^ x with x>1 to prefere sulci
    
      %% general (global) relation between sulci and gyri
      mdGMT  = median(Ygmt0(round(Yp0(:))==2)) * mean(vx_vol); 
      iqrGMT = iqr(Ygmt0(round(Yp0(:))==2)) * mean(vx_vol); 
      gsr    = max(0.5, min(1.7, max(0.5,min(2,1 + (mdGMT - 2.5))) .* max(.5,min(2,(iqrGMT-.5))))); 
    end
    
    %% percentage blurred sulcal/gyral volume (regions that need reconstrution as evaluation parameter)
    srecon = sum(Ygsr(:)<.25 & round(Yp0(:))==2 & (Ygmtw0(:) < Ygmt0(:)*1.1) & (Ygmt0(:)*1.05 <= Ycd0(:)+Ywd0(:)) ) / sum( round(Yp0(:))==2 ); 
    grecon = sum(Ygsr(:)>.75 & round(Yp0(:))==2 & (Ygmtc0(:) < Ygmt0(:)*1.1) & (Ygmt0(:)*1.05 <= Ycd0(:)+Ywd0(:)) ) / sum( round(Yp0(:))==2 ); 
  
    % position maps
    % Yppg - gyrus map with further weighting to avoid bridges
    % Ypps - suclus map with two defintions based on the minimum thickness Ygmt0 and the WM driven thicknes Ygmtw0. 
    %        The Ygmtw0 is better in gyrus reconstruction but also more noisy and prone to bridges
    Yppg = min(1, max(  Ygmt0 .* (Yp0>2.5) , max( Ycd0 .* Ygsr.^.1, (Yp0>1.5) .* (Ygmt0-Ywd0) .* Ygsr.^2)) ./ Ygmt0); 
    Ypps = 0.5 .* min(1, max( 0 , max(Yp0>2.5, (Yp0>1.5) .* max(eps,Ygmtw0-Ywd0) ./ max(eps,Ygmtw0) .* Ygsr ))) + ... 
           0.5 .* min(1, max( 0 , max(Yp0>2.5, (Yp0>1.5) .* min(Ycd0,Ygmt0-Ywd0) ./ max(eps,Ygmt0))));
    % add the global weighting to avoid bridges 
    Yppg = Yppg .* Ypps.^max(0.05, min(.5, .1 * (gsr.^4)));
    % final combination 
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

 
  % evaluation
  if opt.verb
    fprintf('    PP preparation (gsr=%0.3f):     %0.3fs\n', gsr, etime(clock,c)); c = clock; 
    fprintf('    Sulcus / gyrus reconstruction:  %5.2f%% / %4.2f%%\n', srecon*100, grecon*100);  
    fprintf('    Median thickness + IQR:        %5.2f ± %4.2f mm\n', ...
      median( Ygmt( Ypp(:)>.3 & Ypp(:)<.7 )) , iqr( Ygmt( Ypp(:)>.3 & Ypp(:)<.7 ))); 
    x = 0:0.01:10;
    % smooth gives an error here @Robert
    %h = smooth( hist( Ygmt( Ypp(:)>.45 & Ypp(:)<.55 ) , x),2); h = h/sum(h);
    try
      hi = find(h==max(h),1); hil = find(h==max(h(100:hi-30)),1); hih = find(h==max(h(hi+30:end)),1); 
      fprintf('    Peak (x:y):                    %5.2f:%4.4f | %5.2f:%4.4f | %5.2f:%4.4f\n', x(hil), h(hil), x(hi), h(hi), x(hih), h(hih));
    end
  end

  Yp0o = Yp0; 
  
  %% Update Yp0
  % Yp0=Yp0o; Yp0 = max(Yp0,(Ypp>.5) .* 2 + (max(0,min(1,0-Ywd0))));  Yp0 = min(Yp0,2*(Ypp>.5) + 1 + max(0,min(1,max(0,((Ygmt0 - Ywd0 + .75))).^4)));
  Yp0=Yp0o; Yp0 = max(Yp0,(Ypp>.5) .* 2 + (max(0,cat_vol_smooth3X(Ypp,.5)*2-1)).^4); Yp0 = min(Yp0,2*(Ypp>.5) + 2 - max(0,1 - cat_vol_smooth3X(Ypp,.5)*.5).^64); 

  if 0 
    %ds('d2sm','',1,abs(Yp0fs-Ymfs),abs(Yp0e - Ymfs),150)
    %ds('d2sm','',1,Yp0o/3,Yp0/3,140)
    %%
      CS    = isosurface(smooth3(1-Ypp),.5,min(5,Ygmt)); cat_surf_render2(CS); 
      cat_surf_render2('clim',[0 6]); 
  end

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
%cleanupPBT. Upper limit for thickness outliers. 
%  Remove higher outliers (eg. blood vessels & meninges). No filtering of 
%  low values as these also represent small sulci. 
%  
%  Ygmtc = cleanupPBT(Ygmt,lim,lim2)
% 
%  lim1 .. main threshhold
%  lim2 .. additional prefiltering
%

% possible extention:  need resolution to handle filter size better !
  
  if ~exist('lim','var'),  lim  = .05; end
  if ~exist('lim2','var'), lim2 = 1; end

  % overall this causes some systematic underestimation that we try compensate
  thoffsetcor = 0.0;
  
  % basic approximation & filtering
  Ygmta = cat_vol_approx(Ygmt,'rec') + thoffsetcor; 

  % prelimitation to avoid more noisy data with stronger outliers
  if lim2 > 0
    Ygmtc = cat_vol_approx(Ygmt .* ( (Ygmta .* (Ygmt>0) - Ygmt)<=lim*2 & (Ygmt - Ygmta .* (Ygmt>0))<=lim & Ygmt>0), 'rec'); 
    Ygmt  = Ygmt .* (Ygmt>0 & abs( Ygmt - Ygmtc ) < 1); 
    Ygmta = cat_vol_smooth3X( Ygmta .* (Ygmt==0) + Ygmt , .5);
  end

  % main limitation 
  if lim > 0
    % #### the masking looks more complicated then necessary >> test simplification 
    %Ygmtc  = cat_vol_approx(Ygmt .* ( (Ygmt .* (Ygmt>0) - Ygmt)<=2*lim*Ygmta & (Ygmt - Ygmta .* (Ygmt>0))<=lim*Ygmta & Ygmt>0), 'rec'); 
    Ygmtc  = cat_vol_approx(Ygmt .* ( (Ygmta .* (Ygmt>0) - Ygmt)<=lim*2 & (Ygmt - Ygmta .* (Ygmt>0))<=lim & Ygmt>0), 'rec'); % abs
    %Ygmtc  = cat_vol_approx(Ygmt .* ( Ygmt < .5*(Ygmta .* (Ygmt>0)) & Ygmt > 1.25*(Ygmta .* (Ygmt>0)) & Ygmt>0), 'rec'); % rel
    % #### maybe the filtering is here still a bit too strong >> .25 ? 
    Ygmtc  = cat_vol_smooth3X( min(Ygmta,Ygmtc + thoffsetcor),.25); % better
  else
    Ygmtc  = Ygmta; 
  end
end
% ======================================================================
function Yp0 = oneObject(Yp0, vx_vol, n)
% oneObject. Remove addition objects for multiple threshold levels.

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
function [Ycd, Ywd] = cat_vol_cwdist(Yp0,opt)
%cat_vol_cwdist. Estimation of CSF and WM distance in a label map Yp0.
% 
% [Ycd, Ywd] = cat_vol_cwdist(Yp0,opt)
%
% Ycd, Ywd         .. CSF and WM distance maps
% opt              .. parameter structure
%  .levels         .. number of dual distance measurements
%  .extendedrange  .. estimate values also beyond the boundary to improve
%                     thickness mapping
%  .range          .. limitation to avoid bias by interpolation overshoot 
%  .corpve         .. 1-pve, 2-dist, 3-both
%

    corpve = opt.corpve;   % 1-pve, 2-dist, 3-both
    range  = min(.4,opt.range);    % defined by the PVE boundary of the AMAP, i.e range of .25 - .75  
    vxcor  = 0.5;          % correct for voxel to boundary (between voxels) rather than voxel to voxel distance
    pvet   = max(1,min(4,opt.pvet / 2)); % divide by 2 

    
    % The idea is that the we use here the full range of of the 1.5 and 2.5 
    % AMAP class to define the full thickness. However, we measure still
    % from the .5er boundary and we have to handle the 
    YMM = Yp0 > 1.5 - range | isnan(Yp0); 
    YMC = Yp0 < 2.5 + range | isnan(Yp0); 
    if range>0
      YMM = YMM | cat_vol_morph(Yp0 > 1.5,'d',1); 
      YMC = YMC | cat_vol_morph(Yp0 < 2.5,'d',1); 
    end

    % estimation of the extend range to correct values beyond the other tissue boundary
    % additional correction map for values behind tissue boundary, e.g., 
    % for the WMD we estimate the distance from the GM/CSF boundary to 
    % limit WMD values to the maximal thickness value
    % same idea as below
    % CSF/WM limit, ie extended WM/CSF boundary 
    pvemincor = -.25;; % -.5 to 0 ... with more robust result for -.35 
    if corpve == 0
      Ycdlc  = zeros(size(Yp0),'single'); 
      Ywdlc  = zeros(size(Yp0),'single');  
    elseif corpve == 1
      % pure PVE correction  
      Ycdlc = (max(pvemincor,min(.5,Yp0 - 2.5)) * pvet) .* (YMM & YMC);
      Ywdlc = (max(pvemincor,min(.5,1.5 - Yp0)) * pvet) .* (YMM & YMC); 
    elseif corpve >= 2
      % (mixed PVE and) distance map correction
      vxcorc =  .5; % use also minimum for correction  
      mincor = -.5; 
      Ycdlc  = zeros(size(Yp0),'single'); 
      Ywdlc  = zeros(size(Yp0),'single');  
      hss = ceil( opt.levels ); % number of opt.levels (as pairs)
      for si = 1:hss
        offset = max(0,min(.4, range * si/(hss+1))); 
        
        % CSF distance correction beyond WM boudary
       % Ycdlc  = Ycdlc + .5/hss * max(mincor,cat_vol_localstat(cat_vbdist(single(Yp0 < ( 2.5 - offset)), YMM & YMC ) - vxcorc .* (YMM & YMC), YMM & YMC, 1,3)) + ...
       %                  .5/hss * max(mincor,cat_vbdist(single(Yp0 < ( 2.5 + offset)), YMM & YMC ) - vxcorc .* (YMM & YMC)); 
        Ycdlc  = Ycdlc + .5/hss * max(mincor,cat_vbdist(single(Yp0 < ( 2.5 - offset)), YMM & YMC ) - vxcorc .* (YMM & YMC)) + ...
                         .5/hss * max(mincor,cat_vbdist(single(Yp0 < ( 2.5 + offset)), YMM & YMC ) - vxcorc .* (YMM & YMC)); 
      
        % WM distance correction beyond CSF boudary
      %  Ywdlc  = Ywdlc + .5/hss * max(mincor,cat_vbdist(single(Yp0 > ( 1.5 - offset)), YMC & YMM) - vxcorc .* (YMM & YMC)) + ...
      %                   .5/hss * max(mincor,cat_vol_localstat(cat_vbdist(single(Yp0 > ( 1.5 + offset)), YMC & YMM) - vxcorc .* (YMM & YMC), YMM & YMC, 1,3)); 
        Ywdlc  = Ywdlc + .5/hss * max(mincor,cat_vbdist(single(Yp0 > ( 1.5 - offset)), YMC & YMM) - vxcorc .* (YMM & YMC)) + ...
                         .5/hss * max(mincor,cat_vbdist(single(Yp0 > ( 1.5 + offset)), YMC & YMM) - vxcorc .* (YMM & YMC)); 

      end
      Ywdlc(Ywdlc > 1000) = 0;
      Ycdlc(Ycdlc > 1000) = 0; 

      if corpve == 3
      % mix PVE
        Ycdlc = Ycdlc*.5 + 0.5*( (max(pvemincor,min(.5,Yp0 - 2.5)) * pvet) .* (YMM & YMC) );
        Ywdlc = Ywdlc*.5 + 0.5*( (max(pvemincor,min(.5,1.5 - Yp0)) * pvet) .* (YMM & YMC) ); 
      elseif corpve == 4
      % use stronger correction
        Ycdlc = max( Ycdlc , (max(pvemincor,min(.5,Yp0 - 2.5)) * pvet) .* (YMM & YMC) );
        Ywdlc = max( Ywdlc , (max(pvemincor,min(.5,1.5 - Yp0)) * pvet) .* (YMM & YMC) ); 
      end  
    end



    % multi-level distance estimation
    Ycd = zeros(size(Yp0),'single'); 
    Ywd = zeros(size(Yp0),'single'); 
    hss = opt.levels; % number of opt.levels (as pairs)
    mincor2 = 0; 
    for si = 1:hss
      offset = max(0,min(.4, range * si/(hss+1))); 

      % CSF dist
      Ycdl  = max(mincor2,cat_vbdist(single(Yp0 < ( 1.5 - offset)), YMC ) - vxcor .* (YMM & YMC)); Ycdl(Ycdl > 1000) = 0; 
      Ycdh  = max(mincor2,cat_vbdist(single(Yp0 < ( 1.5 + offset)), YMC ) - vxcor .* (YMM & YMC)); Ycdh(Ycdh > 1000) = 0;
      % mixing
      Ycd = Ycd + .5/hss .* Ycdl  +  .5/hss .* Ycdh  -  Ycdlc/(hss*(1 + (corpve==4))); 

      % WM distances
      Ywdl  = max(mincor2,cat_vbdist(single(Yp0 > ( 2.5 - offset)), YMM ) - vxcor .* (YMM & YMC)); Ywdl(Ywdl > 1000) = 0; 
      Ywdh  = max(mincor2,cat_vbdist(single(Yp0 > ( 2.5 + offset)), YMM ) - vxcor .* (YMM & YMC)); Ywdh(Ywdh > 1000) = 0; 
      % mixing
      Ywd = Ywd + .5/hss .* Ywdl  +  .5/hss .* Ywdh  -  Ywdlc/(hss*(1 + (corpve==4)));
    end
    if corpve > 0% == 1 
      Ymsk = Yp0>2.1 & Yp0<2.5; Ywd(Ymsk) = min( max(pvemincor,2.5-Yp0(Ymsk))*pvet , Ywd(Ymsk) ); 
      Ymsk = Yp0>1.5 & Yp0<1.9; Ycd(Ymsk) = min( max(pvemincor,Yp0(Ymsk)-1.5)*pvet , Ycd(Ymsk) ); 
    elseif  corpve == 3
      Ywd = Ywd - (YMM & YMC & Yp0>1.05 & Yp0<1.95) .* max(-.5,1.5-Yp0)*pvet/2; 
      Ycd = Ycd - (YMM & YMC & Yp0>2.05 & Yp0<2.95) .* max(-.5,Yp0-2.5)*pvet/2; 
    elseif  corpve == 4
      Ywd = Ywd - (YMM & YMC & Yp0>1.05 & Yp0<1.95) .* max(-.5,1.5-Yp0)*pvet/2; 
      Ycd = Ycd - (YMM & YMC & Yp0>2.05 & Yp0<2.95) .* max(-.5,Yp0-2.5)*pvet/2; 
    end


    % == Filtering of values on the distant end / opposite boundary. ==
    % This shoud be relative save. As we have a thick subclass and use  
    % also correction before the boundary (with sligt overestimation) 
    % the minimum filters (within the ribbon) are helpful.
    if opt.range > 0;;
      Ymsk = YMM & Yp0<=1.75 & Ywd>0; 
      %Ywd  = Ywd .* (~Ymsk) + cat_vol_localstat(Ywd,Ymsk,2,8,2); 
      for i=1:4, Ywd = cat_vol_median3(Ywd,Ymsk,Ymsk,0); end
      
      Ymsk = YMC & Yp0>=2.25 & Ycd>0; 
      %Ycd  = Ycd .* (~Ymsk) + cat_vol_localstat(Ycd,Ymsk,2,8,2); 
      for i=1:4, Ycd = cat_vol_median3(Ycd,Ymsk,Ymsk,0); end
    end


    % == Final filtering of the full map. == 
    % Filtering here is dangerous and can close structures!
    % So we should include values outside the ribbon and only allow small
    % changes for the median filter.
    if 1 %  worse without - cat_vol_median better then cat_vol_localstat
      Ywd = cat_vol_median3(Ywd,Ywd>.10,true(size(Ywd)) & Ywd>0,0); 
      Ycd = cat_vol_median3(Ycd,Ycd>.10,true(size(Ywd)) & Ywd>0,0); 
    elseif false
      % much worse than the cat_vol_median3
      Ywd = cat_vol_localstat(Ywd,Ywd>0,1.9,8,2); 
      Ycd = cat_vol_localstat(Ycd,Ycd>0,1.9,8,2); 
    end

    
    % A bit general smoothing might help too
    mix = 0; 
    if mix > 0
      if 1
        Ywd = Ywd.*(1-mix | Ywd<=1) + mix.*cat_vol_localstat(Ywd,Ywd>1,1,1); 
        Ycd = Ycd.*(1-mix | Ycd<=1) + mix.*cat_vol_localstat(Ycd,Ycd>1,1,1); 
      else 
        Ywd = Ywd.*(1-mix) + mix.*smooth3(Ywd); 
        Ycd = Ycd.*(1-mix) + mix.*smooth3(Ycd); 
      end
    end
end
% ======================================================================
function [Ycd, Ywd] = cat_vol_cwdist_old(Yp0,opt)
%cat_vol_cwdist. Estimation of CSF and WM distance in a label map Yp0.
% 
% [Ycd, Ywd] = cat_vol_cwdist(Yp0,opt)
%
% Ycd, Ywd         .. CSF and WM distance maps
% opt              .. parameter structure
%  .levels         .. number of dual distance measurements
%  .extendedrange  .. estimate values also beyond the boundary to improve
%                     thickness mapping
%  .range          .. limitation to avoid bias by interpolation overshoot 
% 


    range = 0.5; % defined by the PVE boundary of the AMAP, i.e range of .25 - .75  
    vxcor = 0.5; % this should be .5 but .25 gives the right value - Why? Second AMAP boundary?
    extendedrange = opt.extendedrange;  

    %% additional correction map for values behind tissue boundary, e.g., 
    % for the WMD we estimate the distance from the GM/CSF boundary to 
    % limit WMD values to the maximal thickness value
    % same idea as below
    % CSF/WM limit, ie extended WM/CSF boundary 
    if extendedrange==0
      % The idea is that the we use here the full range of of the 1.5 and 2.5 
      % AMAP class to define the full thickness. However, we measure still
      % from the .5er boundary and we have to handle the 
      YMM = Yp0 > 1.25 | isnan(Yp0); 
      YMC = Yp0 < 2.75 | isnan(Yp0); 
  %    YMM = (cat_vol_morph(Yp0 > 1.25,'d') & Yp0>2.5) | Yp0 > 1.5 | isnan(Yp0); 
  %    YMC = (cat_vol_morph(Yp0 < 2.75,'d') & Yp0>2.5) | Yp0 < 2.5 | isnan(Yp0); 
    else
      YMM = Yp0 >= 1.5 - extendedrange | isnan(Yp0); 
      YMC = Yp0 <= 2.5 + extendedrange | isnan(Yp0); 
    end
    if extendedrange > 0
      YMM = YMM & cat_vol_morph(Yp0 > 1.5,'d',extendedrange); % limit
      YMC = YMC & cat_vol_morph(Yp0 < 2.5,'d',extendedrange); % limit
    end
    
    % estimation of the extend range to correct values beyond the other tissue boundary
    if extendedrange == 1
      Ywdlc = 0; Ywdhc = Ywdlc;
      Ycdlc = 0; Ycdhc = Ycdlc;
    elseif extendedrange == 2 
      vxcorc = 1.5; % use also minimum for correction  
      Ycdlc  = zeros(size(Yp0),'single'); Ycdhc = zeros(size(Yp0),'single'); 
      Ywdlc  = zeros(size(Yp0),'single'); Ywdhc = zeros(size(Yp0),'single'); 
      hss = opt.levels; % number of opt.levels (as pairs)
      for si = 1:hss
        %offset = max(0,min(extendedrange, extendedrange * si/(hss+1))); 
        offset = max(0,min(.4, range * si/(hss+1))); 
        
        % CSF distance correction beyond WM boudary
        Ycdlc  = Ycdlc + 1/hss * max(0,cat_vbdist(single(Yp0 < ( 2.5 - offset)), YMM & YMC ) - vxcorc);
        Ycdhc  = Ycdhc + 1/hss * max(0,cat_vbdist(single(Yp0 < ( 2.5 + offset)), YMM & YMC ) - vxcorc); 
      
        % WM distance correction beyond CSF boudary
        Ywdlc  = Ywdlc + 1/hss * max(0,cat_vbdist(single(Yp0 > ( 1.5 - offset)), YMC & YMM) - vxcorc);  
        Ywdhc  = Ywdhc + 1/hss * max(0,cat_vbdist(single(Yp0 > ( 1.5 + offset)), YMC & YMM) - vxcorc); 
      end
      corlim = -0;
      Ywdlc  = max(Ywdlc, corlim);
      Ywdlc  = max(Ywdlc, corlim);
      Ycdlc  = max(Ycdlc, corlim);
      Ycdhc  = max(Ycdhc, corlim);
      
      Ywdlc(Yp0   > 1.5)  = 0; Ywdhc(Yp0   > 1.5 ) = 0;  
      Ycdlc(Yp0   < 2.5)  = 0; Ycdhc(Yp0   < 2.5 ) = 0;  
      Ywdlc(Ywdlc > 1000) = 0; Ywdhc(Ywdhc > 1000) = 0;  
      Ycdlc(Ycdlc > 1000) = 0; Ycdhc(Ycdhc > 1000) = 0; 
    elseif extendedrange > 0
      % correct for half a voxel 
      corlim = 1;
      Ywdlc  = max(-corlim,min(corlim,(1.5 - Yp0) .* (Yp0>1.0 & Yp0<1.5) .* YMM * 2)); Ywdhc = Ywdlc;
      Ycdlc  = max(-corlim,min(corlim,(Yp0 - 2.5) .* (Yp0>2.5 & Yp0<3.0) .* YMC * 2)); Ycdhc = Ycdlc;
    else
      vxcorc = 0.5; % use also minimum for correction  
      Ycdlc  = zeros(size(Yp0),'single'); Ycdhc = zeros(size(Yp0),'single'); 
      Ywdlc  = zeros(size(Yp0),'single'); Ywdhc = zeros(size(Yp0),'single'); 
      hss = ceil( opt.levels ); % number of opt.levels (as pairs)
      for si = 1:hss
        %offset = max(0,min(extendedrange, extendedrange * si/(hss+1))); 
        offset = max(0,min(.4, range * si/(hss+1))); 
        
        % CSF distance correction beyond WM boudary
        Ycdlc  = Ycdlc + .5/hss * max(0,cat_vbdist(single(Yp0 < ( 2.25 - offset/2)), YMM & YMC ) - vxcorc) + ...
                         .5/hss * max(0,cat_vbdist(single(Yp0 < ( 2.25 + offset/2)), YMM & YMC ) - vxcorc); 
      
        % WM distance correction beyond CSF boudary
        Ywdlc  = Ywdlc + .5/hss * max(0,cat_vbdist(single(Yp0 > ( 1.75 - offset/2)), YMC & YMM) - vxcorc) + ...
                         .5/hss * max(0,cat_vbdist(single(Yp0 > ( 1.75 + offset/2)), YMC & YMM) - vxcorc); 
      end
      corlim = 0;
      Ywdlc  = max(Ywdlc, corlim);
      Ycdlc  = max(Ycdlc, corlim);
      
     % Ywdlc(Yp0   > 1.25) = 0; 
     % Ycdlc(Yp0   < 2.25) = 0; 
      Ywdlc(Ywdlc > 1000) = 0;
      Ycdlc(Ycdlc > 1000) = 0; 
    end


    % multi-level distance estimation
    Ycd = zeros(size(Yp0),'single'); 
    Ywd = zeros(size(Yp0),'single'); 
    hss = opt.levels; % number of opt.levels (as pairs)
    for si = 1:hss
      offset = max(0,min(.4, range * si/(hss+1))); 

      % CSF dist
      Ycdl  = max(-1,cat_vbdist(single(Yp0 < ( 1.5 - offset)), YMC ) - vxcor); Ycdl(Ycdl > 1000) = 0; 
      Ycdh  = max(-1,cat_vbdist(single(Yp0 < ( 1.5 + offset)), YMC ) - vxcor); Ycdh(Ycdh > 1000) = 0;

      if extendedrange > 0
        Ycdl  = Ycdl - Ycdlc;
        Ycdh  = Ycdh - Ycdhc;
      else
        Ycdl  = Ycdl - Ycdlc;
      end
    
      % mixing
      Ycd = Ycd + .5/hss  .* Ycdl  +  .5/hss .* Ycdh; 
      clear Ycdl Ycdh; 

      % WM distances
      Ywdl  = max(0,cat_vbdist(single(Yp0 > ( 2.5 - offset)), YMM ) - vxcor); Ywdl(Ywdl > 1000) = 0; 
      Ywdh  = max(0,cat_vbdist(single(Yp0 > ( 2.5 + offset)), YMM ) - vxcor); Ywdh(Ywdh > 1000) = 0; 

      if extendedrange > 0
        Ywdl   = Ywdl - Ywdlc; 
        Ywdh   = Ywdh - Ywdhc; 
      else
        Ywdl   = Ywdl - Ywdlc; 
      end

      % mixing
      Ywd = Ywd + .5/hss .* Ywdl  +  .5/hss .* Ywdh;
    end

   % Ywd=max(Ywd,0); 
   % Ycd=max(Ycd,0);

    % filtering in final band with equal values
    if extendedrange > 0 %&& extendedrange~=1
      Ymsk = YMM & Yp0<1.5 & Ywd>0; %+extendedrange; 
      if 0% extendedrange%>=1
        Ywd = Ywd .* (~Ymsk) + cat_vol_localstat(Ywd,Ymsk,1,2,extendedrange); 
      end
      Ywd = Ywd .* (~Ymsk) + cat_vol_localstat(Ywd,Ymsk,2,8,8); 
      Ywd = Ywd .* (~Ymsk) + cat_vol_localstat(Ywd,Ymsk,1,1,8); 
      
      Ymsk = YMC & Yp0>2.5 & Ycd>0; %-extendedrange;
      if 0%extendedrange%>=1
        Ycd = Ycd .* (~Ymsk) + cat_vol_localstat(Ycd,Ymsk,1,2,extendedrange); 
      end
      Ycd = Ycd .* (~Ymsk) + cat_vol_localstat(Ycd,Ymsk,2,8,8); 
      Ycd = Ycd .* (~Ymsk) + cat_vol_localstat(Ycd,Ymsk,1,1,8); 
    else
      mincor = 3;
      Ymsk = YMM & Yp0<1.75 & Ywd>0; 
      %Ymsk = YMM & (Yp0<1.5 | (Yp0>=1.5 & cat_vol_morph(Yp0<1.75,'d'))) & Ywd>0; %+extendedrange; 
      %Ywd = Ywd .* (~Ymsk) + Ymsk.*Ywd.*.25 + 0.75.*cat_vol_localstat(Ywd,Ymsk,1,2,1);
      if mincor
        Ymskmin = Ymsk & Yp0<(1.50+1/16); 
        Ywd = Ywd .* (~Ymskmin) + (~Ymskmin) .* cat_vol_localstat(Ywd,YMC & Yp0<=2,2,2,2); 
        clear Ymskmin
      elseif mincor == 2 
        Ywd = Ywd - Ymsk .* min(1,max(-.25,(1.5-Yp0) + .5));
      elseif mincor == 3

      else
        Ywd = Ywd - (Ymsk & Yp0<1.5-1/8) .* 0.35;
        Ywd = Ywd - (Ymsk & Yp0<1.5-1/8 & Yp0<1.5-1/16) .* 0.35/2;
      end
      Ywd = Ywd .* (~Ymsk) + cat_vol_localstat(Ywd,Ymsk,2,8,8); 
      Ywd = Ywd .* (~Ymsk) + cat_vol_localstat(Ywd,Ymsk,1,1,8); 
      
      %Ymsk = YMC & (Yp0>2.5 | (Yp0<=2.5 & cat_vol_morph(Yp0>2.25,'d'))) & Ycd>0; %-extendedrange;
      Ymsk = YMC & Yp0>2.25 & Ycd>0; 
      if mincor
        Ymskmin = Ymsk & Yp0>2.625; 
        Ycd = Ycd .* (~Ymskmin) + (~Ymskmin) .* cat_vol_localstat(Ycd,YMC & Yp0>=2,2,2,2); 
        clear Ymskmin;
      elseif mincor == 2
        Ycd = Ycd - Ymsk .* min(1,max(-.25,(Yp0-2.5) + .5));
      elseif mincor == 3

      else
        Ycd = Ycd - (Ymsk & Yp0>(2.5+1/8)) .* 0.35;
        Ycd = Ycd - (Ymsk & Yp0<(2.5+1/8) & (Yp0>2.5+1/4)) .* 0.35/2;
      end
      Ycd = Ycd .* (~Ymsk) + cat_vol_localstat(Ycd,Ymsk,2,8,8); 
      Ycd = Ycd .* (~Ymsk) + cat_vol_localstat(Ycd,Ymsk,1,1,8); 
    end

    if 0
      Ygmt = (Ycd+Ywd)/2  .* (Yp0>1.5 & Yp0<2.5); 
   
      figure, histogram(Ygmt(Ygmt(:)>0),'binlimit',[0 5], 'binwidth',0.01); xlim([0 5]); grid on
      title( [char(datetime) sprintf(' - %0.2fs',toc)]); 
    end
    
end
% ======================================================================
function [Ywd, Ycd, Ygmt] = keepdetails(Yp0, Ywd, Ycd, Ygmt, vx_vol, extendedrange,level)   
% Although distances and thickness are quite good, PBT slightly tend to  
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
    opt.levels = 2; 

    [Ycd, Ywd] = cat_vol_cwdist(Yp0, opt);
  
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
% ======================================================================
function Ygmt = medfilter( Ygmtw , Ygmtc, Ywd, Ycd, Yp0, red, fs, outliers, mask)
  if ~exist('outliers','var'), outliers = 2; end
  if ~exist('fs','var'), fs = 3; end
  if fs == 0, Ygmt = Ygmtw; return; end

  Ygmtwmsk = Ygmtw>0; 

  % only fully defined voxels
  Ymsk = Ywd>0 & Ycd>0 & Ygmtw>0 & Ygmtc>0; 
  Ygmtw(Ymsk) = max( Ygmtw(Ymsk) , eps );
  Ygmtc(Ymsk) = max( Ygmtc(Ymsk) , eps );

  % no negative values
  Ymsk = Ygmtw<0 | Ygmtc<0 | Ycd<0 | Ywd<0; 
  Ygmtw(Ymsk) = 0; Ycd(Ymsk) = 0; 
  Ygmtc(Ymsk) = 0; Ywd(Ymsk) = 0; 
  
  % remove outliers
  if outliers > 0 
    Ygmtwa = cat_vol_approx(Ygmtw); 
    Ymsk   = (Ygmtw > Ygmtwa * outliers | Ygmtw < Ygmtwa / outliers) & Ygmtw>0; % only positive outliers
    Ygmtw(Ymsk) = Ygmtwa(Ymsk); 
    clear Ygmtwa
  end

  %% cortical band 
  Ypp = min(1,max(Ygmtw>0 & Ygmtw<2,min(min(min(Ygmtw,Ygmtc)-Ywd,Ycd),Ywd) .* (min(Ygmtw,Ygmtc)>0.5)) ./ max(eps,min(Ygmtw,Ygmtc)/2));
  
  %% reduce resolution  
  [Ygmtr,resR] = cat_vol_resize(Ygmtw .* (Ypp>.5),'reduceV',1,red,16,'median'); 
  Ygmtwr       = cat_vol_resize(Ygmtw,'reduceV',1,red,16,'meanm'); 
  Ygmtcr       = cat_vol_resize(Ygmtc,'reduceV',1,red,16,'meanm'); 
  Ycsf         = cat_vol_resize( min(1,2 - Yp0),'reduceV',1,red,16,'meanm'); 
  Ywm          = cat_vol_resize( max(0,Yp0 - 2),'reduceV',1,red,16,'meanm'); 
  Yppr         = cat_vol_resize(Ypp  ,'reduceV',1,red,16,'meanm'); 
  Yppr         = max(0,min(1, (2 * Yppr ./ cat_vol_approx( cat_vol_localstat(Yppr,Yppr>0,1,3))).^2 / 3)) ;
  

  %Ygmtr(Yppr < .75) = 0; 
  mdgmt = median(Ygmtr(Ygmtr(:)>0)); 
  Ygmtr(Yppr < max(.7,min(.9,.8 .* Ygmtr/mdgmt))) = 0; 
  %Ygmtr(smooth3(Ygmtr>0)<.2 ) = 0; 
  Ygmtwr(Ygmtr==0) = 0;
  Ygmtcr(Ygmtr==0) = 0;

  %% main median filter
  Ygmtrm = Ygmtr; sth = [-.5  .1];
  for i = fs %( unique(max(1,min(10,round( ( 10:-1:1 ) / (mean(resR.vx_red / 3)) ))) ))
    if 1
      Ymsk   = Ygmtrm > (1.5./mean(resR.vx_red)); 
      Ygmtrm = Ygmtrm.*(~Ymsk & Ygmtrm>0) + cat_vol_localstat(Ygmtrm,Ymsk,i,8,1); 
      Ygmtrm = cat_vol_median3(  Ygmtrm, Ymsk, Ymsk ); 
    elseif true
      Ygmtrmd = cat_vol_localstat(Ygmtrm,Ygmtrm>0,i,8,1); 
      Ygmtdiv = Ygmtrm - Ygmtrmd; %Ygmtdiv(Ygmtdiv<.1) = 0; 
      Ygmtrm  = (Ygmtwr<=Ygmtcr - sth(1) ) .* min(Ygmtwr,Ygmtrm - Ygmtdiv)  +  ... & Ywm>Ycsf*2
                (Ygmtwr> Ygmtcr - sth(1) & Ygmtwr<=Ygmtcr + sth(2)) .* (Ygmtrm - Ygmtdiv) + ... & Ywm>Ycsf*2
                (Ygmtwr> Ygmtcr + sth(2) ) .* max(Ygmtwr,Ygmtrm - Ygmtdiv); %| Ywm<Ycsf*2
    else
      Ygmtrmd = cat_vol_localstat(Ygmtrm,Ygmtrm>0,i,8,1); 
      %Ygmtrmd = (Ygmtrmd - median(Ygmtrmd(Ygmtrmd(:)>0))) ./ iqr(Ygmtrmd(Ygmtrmd(:)>0)) .* iqr(Ygmtrm(Ygmtrm(:)>0)) + median(Ygmtrm(Ygmtrm(:)>0));
      Ygmtdiv = Ygmtrm - Ygmtrmd; %Ygmtdiv(Ygmtdiv<.1) = 0; 
      Ygmtrm  = ( Ygmtrmd <= median(Ygmtrmd(Ygmtrmd(:)>0))-sth(1) ) .* min(Ygmtwr,Ygmtrm - Ygmtdiv)  +  ... & Ywm>Ycsf*2
                ( (Ygmtrmd >  median(Ygmtrmd(Ygmtrmd(:)>0))-sth(1)) & (Ygmtrmd < median(Ygmtrmd(Ygmtrmd(:)>0))+sth(2)) ) .* (Ygmtrm - Ygmtdiv) + ... & Ywm>Ycsf*2
                ( Ygmtrmd >= median(Ygmtrmd(Ygmtrmd(:)>0))+sth(2) ) .* max(Ygmtwr,Ygmtrm - Ygmtdiv); %| Ywm<Ycsf*2
    end
  end 
  Ygmtr  = Ygmtrm; 
  
  % back to full resolution with approximation of undefined areas
  [~,I] = cat_vbdist(single(Ygmtr>0)); Ygmtr = Ygmtr(I);
  Ygmtr = cat_vol_median3(Ygmtr,Ygmtr>0,Ygmtr>0); 
  Ygmt  = cat_vol_resize(Ygmtr,'dereduceV',resR); 
  Ygmt  = smooth3(Ygmt); 

  % masking
  if mask, Ygmt = Ygmt .* (Ygmtwmsk); end

  % masking cut regions
  Ymsk = cat_vol_morph(cat_vol_morph( cat_vol_morph(Yp0>2.75,'d',2) & cat_vol_morph(Yp0<1.5,'d',2) ,'lo'),'dd',2);
  Ygmt = min(Ygmt, cat_vol_smooth3X(Ygmt .* (1-Ymsk),2)); 
end
% ======================================================================
function [Ygmt,Ywd,Ycd] = outierfilter(Ygmt, Ywd,Ycd, sdth,th)
%outierfilter. Clear local outliers. 
%
%  [Ygmt,Ywd,Ycd] = outierfilter(Ygmt, Ywd,Ycd, sdth,th)
%

  if ~exist('sdth','var'), sdth = 1; end
  if ~exist('th','var'), th = 1.9; end
  
  % define outliers as voxels as deviating voxels
  Ystd   = cat_vol_localstat(Ygmt, Ygmt > 0, 1, 4); 
  Ymsk   = (Ystd < sdth  &  Ygmt > th)  &  Ygmt > 0;

  % replace these voxels by neighbors
  [~,YI] = cat_vbdist( single(Ymsk), Ygmt>0); 
  Ygmt   = Ygmt(YI); clear YI; 
  Ygmt   = min(Ygmt, Ywd + Ycd);

  % filter neighbor voxels 
  Ymsk   = cat_vol_morph(Ystd >= sdth | ( Ygmt<=th & Ygmt>0),'d')  &  Ygmt > 0;
  Ygmt   = cat_vol_median3(  Ygmt, Ymsk > 3, Ygmt > 0 , 0); 

  % at least limit the distance values ...
  Ywd = min(Ywd,Ygmt);
  Ycd = min(Ycd,Ygmt);
  
  Ygmt = min( ...
    min(Ygmt, cat_vol_median3(  Ygmt, Ygmt > 0 & Ygmt < 3 , Ygmt > 0 & Ygmt < 3  , 0 )), ...
    cat_vol_median3(  Ygmt, Ygmt > 3, Ygmt > 0 , 0 )); 
  
  if 0
    Ymsk = (Ywd+Ycd) - Ygmt; 
    Ywd(Ymsk>0) = Ywd(Ymsk>0) - Ymsk(Ymsk>0)/2;  
    Ycd(Ymsk>0) = Ycd(Ymsk>0) - Ycsk(Ymsk>0)/2;  
  end

end
% ======================================================================
function Yp0s = cleanPVE(Yp0,level)
%cleanPVE. Use median (and local) filter to reduce interpolation artifacts.
%
%  Yp0s = cleanPVE(Yp0[,level])
%
%  Yp0   .. label map 
%  Yp0s  .. filtered label map
%  level .. 1-median, 2-add PVE sharpening, 3-add light final smoothing
% 
Yp0o=Yp0; 
%%
Yp0=Yp0o; 
  if ~exist('level','var'), level = 2; end
 % if level <= 0, Yp0s = Yp0; return; end
  
  % just a median filter in the main GM region
  Ymsk = cat_vol_morph( Yp0>1.25 & Yp0<2.75 , 'd',2); 
  Yp0s = cat_vol_median3(Yp0, Ymsk, Ymsk, max(0,1-level));  
%  if level <= 1, return; end
  
  % PVE boundary to enhance/sharpen (value +.5)
  for i = [1.0 2.0] 
    Ymski = Yp0>i+.1 & Yp0<i+1-.1; 
    Yp0ss = Yp0s + (Yp0s*2 - cat_vol_smooth3X(Yp0s,2) - cat_vol_smooth3X(Yp0s,1)); 
    Yp0l  = cat_vol_laplace3R(Yp0ss,Ymski,.01); 
    Yp0s  = Yp0s  +  Ymski .* max(-.5,min(.5, (Yp0 - Yp0l)));
  end 
  Yp0s = cat_vol_median3(Yp0s, Ymsk, Ymsk, max(0,2-level)); 
 %%
  if level <= 2, return; end
  
  % 0.5 was optimal to avoid artifacts in phantom but it should be as small as possible
  Yp0s = cat_vol_smooth3X(Yp0s, 0.25);  

end
% ======================================================================
function pvet = estimatePVEsize( Yp0 , fast ) 
%estimatePVEsize. Estimate the voxel-size of partial volume label map.
%
%  pvet = estimatePVEsize( Yp0 [, fast]) 
% 
%  Yp0  .. label map
%  pvet .. partial volume effect size in voxel
%  fast .. use fast approximation (default) or full estimation (=0)

  if ~exist('fast','var'), fast = 1; end

  if fast
  % fast approximation
    Ypvet = cat_vbdist(single( Yp0<1.1 | (Yp0>1.9 & Yp0<2.1) | Yp0>2.9 )) * 2.5; 
    Ypvet(Ypvet>100 | Ypvet<0) = 0;
  else
  % accurate estimation 
    Ycgd  = cat_vbdist(single(Yp0<1.1), Yp0<1.9); Ycgd(Ycgd>100 | Ycgd<0) = 0; 
    Ygcd  = cat_vbdist(single(Yp0>1.9), Yp0>1.1); Ygcd(Ygcd>100 | Ygcd<0) = 0; 
    Ygwd  = cat_vbdist(single(Yp0<2.1), Yp0<2.9); Ygwd(Ygwd>100 | Ygwd<0) = 0; 
    Ywgd  = cat_vbdist(single(Yp0>2.9), Yp0>2.1); Ywgd(Ywgd>100 | Ywgd<0) = 0; 

    Ypvet = max(Ycgd + Ygcd, Ygwd + Ywgd); 
  end

  % final evaluation
  pvet  = max(1,median(Ypvet(Ypvet(:)>0) - 1)); 
end
% ======================================================================
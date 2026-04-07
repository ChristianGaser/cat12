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
%    .range (real value <=0.5, good between 0.2 and 0.45; default=0.45)
%      Limitation of the offset of multiple thickness levels to avoid 
%      running into partial volume effects with thickness overestimation.
%
%    .NBVC (0-no, 1-yes; default=1)
%      Additional, new blood vessel correction based on a WM region growing.
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
  def.levels           = 2;   % Number of dual distance estimates. 
                              % Larger values are more accurate but need more time (log-change).
                              % 1 means estimation for .25 and .75 boundaries
                              % Good values are between 2 and 8.
                              %  
  def.range            = .4;  % Default value for range extension for *first* boundary (should be between 0.2-0.45)
                              %
  def.rangeE           = .4;  % Default value for range extension for *second* boundary(should be between 0.2-0.45)
                              %
  def.NBVC             = 1;   % new blood vessel correction 
                              %
  def.gyrusrecon       = 1;   % use PBT also to reconstruct gyri 
                              % key aspect for better surfaces but also risk for blurred sulci
                              %
  def.verb             = 0;   % be verbose
                              %
  def.usemedian        = 1;   % GMT filtering setting (0-approx,1-median)
                              %
  def.medfs            = [1 2]; % [processingResolution filterSize] default [1.5 2]; 
                              %
  def.eidist           = 1;   % eikonal distance (not so good yet)

  opt = cat_io_checkinopt(opt,def);
  vxs = mean(vx_vol); 

  opt.vxs    = vxs;
  opt.levels = max(1,min(8,opt.levels  * mean(vx_vol*2)^2)); 
    
  
  c = clock; %#ok<*CLOCK>

% TODO: 
% - remove tic-toc run-time evaluation blocks
  if opt.verb, fprintf('\nPrep:    '); tic; end  
  

  % estimate partial volume effect size in voxel (with fast option)
  % - used to controle filter size (eg. closing) 
  opt.pvet0 = estimatePVEsize( Yp0 , 0); 
  opt.pvet  = max(0,min(4 / vxs,opt.pvet0)); 

  % correct CSF interpolation artifacts
  Yp0 = CS4_desnoise(Yp0,opt);

  if mean(vx_vol) < .75 
    % extend hard cuts, by a low value just to have a broader estimate
    % - could be improved by a basic WMD and GMT estimate to assure a
    %   local minimum thickness related to the GMT distance
    % - RD202601: important for small distances but less for large
    % - Ok for phantom.
    Yp0 = max(Yp0,(Yp0==1 & cat_vol_morph(Yp0==2,'d')) * 1.2);

    Ymsk = cat_vol_smooth3X(Yp0>1.5 & Yp0<2.5,2); 
    Yp0  = max(Yp0, Ymsk .* 2.90 .* smooth3(cat_vol_morph(Yp0>2.90,'ldc',opt.pvet/1)) ); 
    Yp0  = max(Yp0, Ymsk .* 2.75 .* smooth3(cat_vol_morph(Yp0>2.75,'ldc',opt.pvet/2)) ); 
  end
  

  % corrections for blood vessels, partial volume effects and skull-stripping
  if opt.NBVC, Yp0 = NBVC(Yp0, vx_vol); end 
 

  % estimate hull to limit eg. gyrus reconstruction 
  Yhull = estimateHull(Yp0,vx_vol);


  %% == Estimation of distance measures ==
  % It is important to use the PVE for multiple (=levels) distance estimations 
  % but also to extened (range) and correct (corpve) the values at the oposite  
  % boundary by the PVE or distance. 

  % tested here multple estimations
  if opt.verb, toc; fprintf('Dist:    '); tic; end  
  [Ycd0, Ywd0] = cat_vol_cwdist( Yp0 ,opt);
  Ycd0(Ycd0>1000) = 0; Ywd0(Ywd0>1000) = 0; % lazy debugging

  %% == Estimation of thickness maps == 
  if opt.verb, toc; fprintf('Thick:   '); tic; end      
  Ygmtw0 = cat_vol_pbtp( single(1 + (Yp0>=1.5-opt.rangeE) + (Yp0>2.5+opt.rangeE)) , Ywd0, Ycd0);
  Ygmtc0 = cat_vol_pbtp( single(3 - (Yp0>=1.5-opt.rangeE) - (Yp0>2.5+opt.rangeE)) , Ycd0, Ywd0); 
  Ygmtw0(Ygmtw0>1000) = 0; Ygmtc0(Ygmtc0>1000) = 0; % lazy debugging

  %% extend areas that are too thin for regular quantification 
  %  This is typically the cutting area of the corpus callosum.
  %  That is so far not in our thickness maps but needs to be present
  %  in the following filter steps!
  minthick = 0.05/vxs; 
  Ycut   = Ycd0>0 & Ywd0>0 & Ycd0<vxs & Ywd0<vxs & Ygmtw0<vxs & Ygmtc0<vxs;
  % Ycut   = Yp0>1.25 & Yp0<2.75 & Ygmtw0<vxs & Ygmtc0<vxs;
  %Ycut   = cat_vol_morph(Ycut,'l',[10 .1]) > 0; 
  Ycut   = Ycut | ( cat_vol_morph(Ycut,'dd',1,vx_vol) & Ygmtw0>0); 
  Ycuts  = cat_vol_localstat( min(Ygmtw0,Ygmtc0) , Ycut, max(1,round(1/vxs)) , 2, 1); 
  Ycuts  = cat_vol_localstat( Ycuts , Ycut, max(1,round(1/vxs)) , 1, 1); 
  Ygmtw0(Ycut) = max(minthick, Ycuts(Ycut)); 
  %Ygmtc0(Ycut) = max(minthick, Ycuts(Ycut)); clear Ycuts; 
  Ygmtw0 = cat_vol_median3( Ygmtw0, cat_vol_morph(Ycut,'d',1) & Ygmtw0>0, Ygmtw0>0 ); 
  %Ygmtc0 = cat_vol_median3( Ygmtc0, cat_vol_morph(Ycut,'d',1) & Ygmtc0>0, Ygmtc0>0 ); 
 

  %% Correct thickness and distance in narrow areas of about a voxel
  % When the thickess is below about one voxel, we have to assume that it 
  % overestimated by a factor of 2. 
    % we assume bended cortical band adjacent to both sides. 
      %  4 corrections of the same thing ... not good   
  Ypve = (Ygmtw0>0 & Ywd0<1.5 & Ygmtw0<1.9 & Ygmtw0<Ygmtc0) | Ycut;
  Ypve(smooth3(Ypve)<.5) = 0; 
  Ygmtw0c = min(Ygmtw0,max(eps,Ygmtw0*3/4).^2); 
  Ygmtw0  = Ygmtw0.*(1-Ypve) + Ypve.*Ygmtw0c; clear Ygmtw0c;
  Ywd0    = Ywd0.*(1-Ypve) + (Ypve).*min(2,max(eps,Ywd0 * 2)); % 2 would be correct but 3 opens a bit better

  if opt.verb, toc; fprintf('Filter:  '); tic; end  


  %% minimum tickness map that represents both the sulcus as well as the gyrus reconstruction 
  Ygmt0 = min(Ygmtw0,Ygmtc0); 

  % Median smoothing or appoximation 
  % an update of the distance maps is not required (and gave worse results)
  if opt.usemedian 
    Ygmt0m = medfilter(Ygmt0, Yp0, opt.medfs(2), opt.medfs(2)/opt.medfs(1), 0, vxs); 
    Ygmt0  = smooth3( Ygmt0m  +  cat_vol_approx(Ygmt0m,'rec',2) .* (Ygmt0m<.25) ); 
  else
    Ygmt0 = cat_vol_approx(Ygmt0,'rec',2);  
  end
  if opt.verb, toc; end  
  if opt.verb, fprintf('  Thickness estimation:           %0.3fs\n', etime(clock,c)); end %#ok<*DETIM>
  


  %% CSF/WM blurring/reconstruction maps
  % - define the blurred sulcal/CSF areas as the area of WM closing and 
  %   'overestimated' thickness (i.e., where the PBT thickness is sign. 
  %   smaller as the simple sum of the WM and CSF distance)
  % - blurred gyri/WM is defined vite versa
  % - undefined values GM values are defined by neighbours
  % - neutral regions are defined by CSF
  % - next both areas are extend by intensity emphasized values 
  gmt = [cat_stat_nanmedian(Ygmtw0(Ygmtw0(:)>0)) cat_stat_nanstd(Ygmtw0(Ygmtw0(:)>0))]; 
  if opt.gyrusrecon 
    % run on 1 mm should be sufficient
    res = .5; 
    [Yp0r,resR]              = cat_vol_resize(Yp0,'reduceV',vx_vol,res,16,'meanm'); 
    [Ygmtw0r,Ygmtc0r,Ygmt0r] = cat_vol_resize({Ygmtw0,Ygmtc0,Ygmt0},'reduceV',vx_vol,res,16,'meanm'); 
    
    % low-res voxel-size, average thickness and deviation 
    vxsr = mean( resR.vx_volr ); 
    Ywm  = Yp0r > 2.125; 
    Ycsf = Yp0r < 1.875;

    % define blurred sulci/gyri by closing of WM and CSF
    Yscr = cat_vol_morph( Ywm  , 'dc', (2*gmt(1) + 1*gmt(2)) * vxs / vxsr ) + Ywm; 
    Ygcr = cat_vol_morph( Ycsf , 'dc', (2*gmt(1) + 1*gmt(2)) * vxs / vxsr ) + Ycsf; 
    
    % extend these regions to avoid WM/CSF bridges 
    [~,I] = cat_vbdist(single(cat_vol_morph(Yscr>0,'e')),~Ycsf);  Yscr = single(Yscr(I)); Ygmt0sr = Ygmt0r(I); clear I; 
    [~,I] = cat_vbdist(single(cat_vol_morph(Ygcr>0,'e')),~Ywm);   Ygcr = single(Ygcr(I)); Ygmt0gr = Ygmt0r(I); clear I;
    Ygmt1r = min(Ygmt0r,min(Ygmt0sr,Ygmt0gr)); 
    Ymskr  = (Ygmt1r<Ygmt0r | Ygmt1r<Ygmt0sr | Ygmt1r<Ygmt0gr) & Yp0r>1; clear Ygmt0sr Ygmt0gr;
    Ygmt1r = cat_vol_median3(Ygmt1r,Ymskr,Ymskr); 
    Ygmt0r = min(Ygmt0r, Ygmt1r) + 10*(1-Ymskr); clear Ygmt1r; 
    Ygmt0m = min(Ygmt0, cat_vol_resize(Ygmt0r,'dereduceV',resR)); % update main gmt map in mask area
    Ygmt0  = min(Ygmt0, smooth3( Ygmt0m ));

    % define gyrus (=1) / sulcus (=0) map
    Ygsr = (Yscr==2 | ( Ygcr==1 & Ygmtw0r > 0  &  Ygmtw0r > Ygmtc0r ))  &  ...  % gyrus
          ~(Ycsf | Yp0r==3 | Yscr==1 | cat_vol_smooth3X(Ygmtw0r < gmt(1) - 0.5 * gmt(2),1)>.125 );  % sulcus 
    Ygsr = cat_vol_smooth3X( cat_vol_resize(Ygsr,'dereduceV',resR) ,1) .* (Ygmt0 / gmt(1)) .* max(0,Yp0-1.5).^.1; 
    Ygsr = min(Ygsr, cat_vol_smooth3X(Yhull,1)); 
    Ygsr = Ygsr .* max(0.25,1 - max(0,vxs-.5)/2);
   
    % percentage blurred sulcal/gyral volume (regions that need reconstrution as evaluation parameter)
    % position maps
    % Yppg - gyrus map with further weighting to avoid bridges
    % Ypps - suclus map with two defintions based on the minimum thickness Ygmt0 and the WM driven thicknes Ygmtw0. 
    %        The Ygmtw0 is better in gyrus reconstruction but also more noisy and prone to bridges
    Yppg = min(1, max( 0 , max(Yp0>2.5, (Yp0>1.5) .* 1 - ( min(Ywd0,Ygmt0-Ycd0) ./ max(eps,Ygmt0)))));
    Ypps = min(1, max( 0 , max(Yp0>2.5, (Yp0>1.5) .*       min(Ycd0,Ygmt0-Ywd0) ./ max(eps,Ygmt0))));
 
    % final combination 
    Ypp  = Yppg.*Ygsr + (1-Ygsr).*Ypps; 
  else
    % only sulcus reconstruction 
    Ypp  = min(1, max( 0 , max(Yp0>2.5, (Yp0>1.5) .* min(Ycd0,Ygmt0-Ywd0) ./ max(eps,Ygmt0))));
  end 

  
  %% cleanup of position map
  Ypp = max((Yp0-2), max(Ypp, smooth3(Ypp) .* (Ygmt0>gmt(1)*.5) )); % stabize low res fine structures
  Ypp = max(0,min(Yp0-1,Ypp));
  Ypp(cat_vol_morph(Yp0>2.75,'l')) = 1; Ypp(Yp0<1.5) = 0; % close holes
  Ypp = min(Ypp,cat_vol_median3(Ypp,smooth3(Ypp)>0 & smooth3(Ypp)<.75, Ypp>-1, .5)); % prefere open
  Ypp = max(Ypp,cat_vol_median3(Ypp,smooth3(Ypp)>0 & Ygmt0>gmt(1), Ypp>-1, .5));     % prefere close
  % use PVE if better?
  Ypp = max(Ypp, max(0,min(.5,Yp0-2)*2).^2 );  % ^.5 stronger
  Ypp = min(Ypp, min(1,max(0,Yp0-1.5)*2).^2 );  % ^2  stronger

  % final scaling
  Ygmt = Ygmt0 * vxs; 
  
  % final evaluation
  if opt.verb
    fprintf('    Median thickness + IQR:        %5.2f ± %4.2f mm\n', ...
      median( Ygmt( Ypp(:)>.3 & Ypp(:)<.7 )) , iqr( Ygmt( Ypp(:)>.3 & Ypp(:)<.7 ))); 
  end
  
  if 0 
    %% only internal tests
    CS    = isosurface(smooth3(1-Ypp),.5,min(5,Ygmt)); cat_surf_render2(CS); 
    cat_surf_render2('clim',[0 6]); 
  end

end
% ======================================================================
function Yhull = estimateHull(Yp0,vx_vol)
% estimateHull. Create a smooth hull that inclose the brain tissue. 
%  For fast estimation we go to lower resolution (2 mm) and run their   
%  a closing operation. 
  [Yp0r,resp] = cat_vol_resize(Yp0,'reduceV',vx_vol,1,32,'meanm'); 
  Yp0r  = cat_vol_morph( Yp0r>1.25 ,'do',.5,resp.vx_volr); 
  [Ywr,resYp0w] = cat_vol_resize(Yp0r,'reduceV',resp.vx_volr,2,32,'meanm'); clear Yp0r;
  Ywr   = cat_vol_morph(Ywr>.5,'ldc',8);
  Ywr   = cat_vol_smooth3X(Ywr,1.2); 
  Ywr   = cat_vol_resize(Ywr,'dereduceV',resYp0w); 
  Yhull = cat_vol_resize(Ywr,'dereduceV',resp); clear Ywr; 
  Yhull = max(0,min(1,atan(Yhull*pi/1.5))); 
end
% ======================================================================
function Yp0 = CS4_desnoise(Yp0,opt)
% filter interolation artifacts (important for eikonal speedmap)
  inth  = cat_stat_nanmedian(Yp0(Yp0(:)>.5 & Yp0(:)<1)); 
  Yp0   = min(3,max(1,Yp0)); 
  if inth>.01
    Yp0 = cat_vol_median3(Yp0, Yp0>1.9 & Yp0<2.1,true(size(Yp0)),inth/2);
    Yp0 = cat_vol_median3(Yp0, Yp0>2.9 & Yp0<3.5,true(size(Yp0)),inth/2);
  end

  %% main filter depending on WM image _and_ anatomical noise to handle WMHs

  Yp0om = Yp0; 
  FECth = min(1,max(0,opt.pvet - 2)*2); 
  
  Ymm = cat_vol_morph( Yp0>2 & Yp0<3 , 'd'); 
  Yp0 = Yp0om; 
  opt.wmnoise = 0.02;
  for i = 1:opt.wmnoise*100/3 
    Yp0f = Yp0;% .* (Yp0<=2) + cat_vol_localstat(Yp0,Ymm,2, 8 );
    Yp0f = cat_vol_median3(Yp0f, Ymm, Yp0>0, .05 );
    
    % refinement
    Yp0d = abs(Yp0f-Yp0); 
    Yp0a = cat_vol_approx(Yp0d); 
    Yp0d = max( 0,Yp0 - Yp0a ) ./ Yp0a*2;
    Yp0d = min(1,Yp0d .* smooth3(Yp0>2 & Yp0<3)); 
    Yp0  = Yp0 .* (1-FECth.*Yp0d) + FECth .* Yp0d .* Yp0f;  
    Yw   = cat_vol_smooth3X(Yp0./Yp0om,4); 
    Yp0  = min(3,max(1,Yp0 ./ Yw));
    Ymm  = Ymm & Yp0d>.01; 
  end
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
function Yp0 = NBVC(Yp0,vx_vol)
%% RD202503: new blood vessel correction 
% In principle this is not really new and should be done by other functions before. 
% However, it is so essential for PBT and not to difficult to include it here (too). 
% First the save WM area is estimated and then extended by a region growing method
% that focus on lower intensities in the Yp0 label mapf for GM and CSF. 

  F      = min(1,max(0.0,Yp0-2)).^.01; F(Yp0<=1.9) = inf; 
  Ywm    = cat_vol_morph(cat_vol_morph(Yp0>2.25,'de',1,vx_vol),'ldo',1,vx_vol) | cat_vol_morph(Yp0>2.75,'l'); 
  % add closest voxels
  [~,Yd] = cat_vol_downcut(single(Ywm), F,0.01); Yd(isnan(Yd))=inf; 
  Ywm    = Yp0>2.1 & Yd<median(Yd(Yp0(:)>=2.25 & Yp0(:)<2.75 & Yd(:)>0));

  %% prepare mask
  [~,Yd] = cat_vol_downcut(single(Ywm), F,-0.001); Yd(isnan(Yd))=inf; %clear Ywm; 
  Ymsk   = Yd > 20000 * mean(vx_vol) & Yp0>2.25 & ~cat_vol_morph(Ywm,'dd',1); 
 
  % Ymsk as regions that we want to correct 
  %[~,I]  = cat_vbdist(single(~Ymsk),Yp0>1.5); Yp0(Ymsk) = Yp0(I(Ymsk));
  Yp0(Ymsk) = 2; 
  Ymsk      = cat_vol_morph(Ymsk,'dd',1,vx_vol) & Yp0>1.5;
  Yp0s      = cat_vol_median3(Yp0,Ymsk);
  Yp0(Ymsk) = Yp0s(Ymsk); 

  %% skull-stripping 2
  F      = min(1,max(0.0,Yp0-1).^2); F(Yp0<=1.01) = inf; 
  [~,Yd] = cat_vol_downcut(single(Ywm), F,0.1); Yd(isnan(Yd))=inf; 
  [~,Yd] = cat_vol_downcut(single(Yd<3*cat_stat_nanmedian(Yd(Yp0>1.9 & Yp0<2.25)) & Yp0>1.95), F,-0.01); Yd(isnan(Yd))=inf; 
  Ymsk   = smooth3(Yd > cat_stat_nanmedian(Yd(Yp0>1.25 & Yp0<1.75)))>.5;
  % need the core for the insula
  [Ycore,resR] = cat_vol_resize(Yp0>2.5,'reduceV',vx_vol,2,32,'median');
  Ycore = smooth3(cat_vol_morph(Ycore,'dc',4));
  Ycore = cat_vol_resize(Ycore,'dereduceV',resR); 
  Ymsk  = cat_vol_morph( Ymsk & Ycore<.5,'e');
  %%
  Yp0(Ymsk) = 1; 
  Ymsk      = cat_vol_morph(Ymsk,'dd',1,vx_vol) & Yp0>1 & Yp0<2;
  Yp0s      = cat_vol_median3(Yp0,Ymsk);
  Yp0(Ymsk) = Yp0s(Ymsk); 

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
%  .range          .. limitation to avoid bias by interpolation overshoot 
%
    
  def.levels = 2; 
  def.rangeE = 0.50; % extimation limit extention - smaller better
  def.range  = 0.50; % addition range - larger better but limited by interpolation artifacts
  def.vxs    = 1;
 
  % Only for quick tests!
  % There are two good combinations, with rangeE=.4: 
  %  1) voxelcorr/voxelcorrE = 0/1
  %  2) voxelcorr/voxelcorrE = 2/1
  % Other combinations are worse.
  voxelcorr  = 2; % 0-none, 1-simple vx-correction-worse, 2-PVE-correction (default=2)
  voxelcorrE = 1; % 0-none, 1-simple vx-correction, 2-PVE-correction-worse (default=1)
  
  opt = cat_io_checkinopt(opt,def); 

  range  = min(.4,opt.range);    % defined by the PVE boundary of the AMAP, i.e range of .25 - .75  
  rangeE = min(.4,opt.rangeE);   % defined by the PVE boundary of the AMAP, i.e range of .25 - .75  
  
  % The idea is that the we use here the full range of of the 1.5 and 2.5 
  % AMAP class to define the full thickness. However, we measure still
  % from the .5er boundary and we have to handle the 
  YMM = Yp0 >= 1.5 - rangeE; % range for WMD
  YMC = Yp0 <= 2.5 + rangeE; % range for CSFD
  
  % multi-level distance estimation
  Ycd = zeros(size(Yp0),'single'); 
  Ywd = zeros(size(Yp0),'single'); 
  hss = opt.levels; % number of opt.levels (as pairs)
  if opt.eidist 
    % speed/time map for eikonal distance
    Fc = smooth3(max(.001,min(1,cat_vol_smooth3X(min(1,max(eps,2 - Yp0)),2).^4))); Fc = max(Fc,cat_vol_morph(Yp0<2,'dd',1.9));
    Fw = smooth3(max(.001,min(1,cat_vol_smooth3X(min(1,max(eps,Yp0 - 1)),2).^4))); Fw = max(Fw,cat_vol_morph(Yp0>2,'dd',1.9));
  end
  for si = 1:hss
    offset = max(0,min(.5, range * si/(hss+1))); 

    % CSF dist
    if opt.eidist && exist('cat_vbdist3','file')
      [Ycdl,YIl] = cat_vbdist3(single(Yp0 < ( 1.5 - offset)), Fc .* YMC); 
      [Ycdh,YIh] = cat_vbdist3(single(Yp0 < ( 1.5 + offset)), Fc .* YMC); 
    else
      [Ycdl,YIl] = cat_vbdist(single(Yp0 < ( 1.5 - offset)), YMC ); 
      [Ycdh,YIh] = cat_vbdist(single(Yp0 < ( 1.5 + offset)), YMC ); 
    end
    % PVE correction 
    if voxelcorr==2
      vxc  = .5; 
      Ycdl = (Ycdl - min(vxc,max(0,- (Yp0(YIl) - (1.5-offset))*opt.pvet ))) .* (YMM & YMC); 
      Ycdh = (Ycdh - min(vxc,max(0,- (Yp0(YIh) - (1.5+offset))*opt.pvet ))) .* (YMM & YMC); 
    elseif voxelcorr==1
      Ycdl = max(0,Ycdl - 0.5) .* (YMM & YMC);
      Ycdh = max(0,Ycdh - 0.5) .* (YMM & YMC);
    end
    % classic correction (in interaction with second boundary correction!)
    Ycd  = Ycd + .5/hss .* Ycdl  +  .5/hss .* Ycdh; 

    % WM distances
    if opt.eidist && exist('cat_vbdist3','file')
      [Ywdl,YIl] = cat_vbdist3(single(Yp0 > ( 2.5 - offset)), Fw .* YMM);
      [Ywdh,YIh] = cat_vbdist3(single(Yp0 > ( 2.5 + offset)), Fw .* YMM); 
    else
      [Ywdl,YIl] = cat_vbdist(single(Yp0 > ( 2.5 - offset)), YMM );
      [Ywdh,YIh] = cat_vbdist(single(Yp0 > ( 2.5 + offset)), YMM ); 
    end
    % PVE correction 
    if voxelcorr==2
      Ywdl = (Ywdl - min(vxc,max(0, (Yp0(YIl) - (2.5-offset))*opt.pvet ))) .* (YMM & YMC); 
      Ywdh = (Ywdh - min(vxc,max(0, (Yp0(YIh) - (2.5+offset))*opt.pvet ))) .* (YMM & YMC); 
    elseif voxelcorr==1
      Ywdl = max(0,Ywdl - 0.5) .* (YMM & YMC);
      Ywdh = max(0,Ywdh - 0.5) .* (YMM & YMC);
    end
    Ywd  = Ywd + .5/hss .* Ywdl  +  .5/hss .* Ywdh;
  end
  Ycd = max(0,Ycd);
  Ywd = max(0,Ywd);


  %% endpoint PVE correction 
  for si = 1:hss
    offset = max(0,min(.5, range * si/(hss+1))); 

    % CSF dist
    if opt.eidist && exist('cat_vbdist3','file') 
      [Ycdl,YIl] = cat_vbdist3(single(Yp0 < ( 2.5 - offset)), YMC); 
      [Ycdh,YIh] = cat_vbdist3(single(Yp0 < ( 2.5 + offset)), YMC); 
    else 
      [Ycdl,YIl] = cat_vbdist(single(Yp0 < ( 2.5 - offset)), YMC ); 
      [Ycdh,YIh] = cat_vbdist(single(Yp0 < ( 2.5 + offset)), YMC ); 
    end
    if voxelcorrE==2
      Ycdl = (Ycdl - min(vxc,max(0,- (Yp0(YIl) - (2.5-offset))*opt.pvet ))) .* (YMM & YMC); 
      Ycdh = (Ycdh - min(vxc,max(0,- (Yp0(YIh) - (2.5+offset))*opt.pvet ))) .* (YMM & YMC); 
    elseif voxelcorrE==1
      Ycdl = max(0,Ycdl - 0.5) .* (YMM & YMC);
      Ycdh = max(0,Ycdh - 0.5) .* (YMM & YMC);
    end
    Ycd  = max(0,Ycd - (.5/hss .* Ycdl  +  .5/hss .* Ycdh)); 

    % WM distances
    if opt.eidist && exist('cat_vbdist3','file')
      [Ywdl,YIl] = cat_vbdist3(single(Yp0 > ( 1.5 - offset)), YMM);
      [Ywdh,YIh] = cat_vbdist3(single(Yp0 > ( 1.5 + offset)), YMM); 
    else
      [Ywdl,YIl] = cat_vbdist(single(Yp0 > ( 1.5 - offset)), YMM );
      [Ywdh,YIh] = cat_vbdist(single(Yp0 > ( 1.5 + offset)), YMM ); 
    end
    if voxelcorrE==2
      Ywdl = (Ywdl - min(vxc,max(0, (Yp0(YIl) - (1.5-offset))*opt.pvet ))) .* (YMM & YMC); 
      Ywdh = (Ywdh - min(vxc,max(0, (Yp0(YIh) - (1.5+offset))*opt.pvet ))) .* (YMM & YMC); 
    elseif voxelcorrE==1
      Ywdl = max(0,Ywdl - 0.5) .* (YMM & YMC);
      Ywdh = max(0,Ywdh - 0.5) .* (YMM & YMC);
    end
    Ywd  = max(0,Ywd - (.5/hss .* Ywdl  +  .5/hss .* Ywdh));
  end
end
% ======================================================================
function Ygmt = medfilter( Ygmtw , Yp0, red, fs, mask, vxs)
%medfilter. Simple median filter.

  if ~exist('vxs','var'), vxs = 1; end
  if ~exist('fs','var'), fs = 1; end
  if fs == 0, Ygmt = Ygmtw; return; end

  % save tissue mask for final masking later
  Ygmtwmsk = Ygmtw>0; 

  % fill-up
  [~,I] = cat_vbdist(single(Ygmtwmsk)); Ygmtw = Ygmtw(I); 

  % reduce resolution  
  [Ygmtr,resR] = cat_vol_resize(Ygmtw,'reduceV',vxs,red,16,'median'); 
  Yp0r         = cat_vol_resize(Yp0,  'reduceV',vxs,red,16,'median'); 
  Ymsk         = Yp0r>1 & Yp0r<3; 
     
  % main median filter
  Ygmtmd = Ygmtr; 
  for i = 1:fs
    Ygmtmd = cat_vol_median3(Ygmtmd,Ymsk); 
  end 
  Ygmtr = Ygmtmd; 

  % back to full resolution with approximation of undefined areas and slight smoothing
  Ygmt  = cat_vol_resize(Ygmtr,'dereduceV',resR); 

  % final tissue masking
  if mask, Ygmt = Ygmt .* (Ygmtwmsk); end
 
end
% ======================================================================

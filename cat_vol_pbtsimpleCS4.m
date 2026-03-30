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
  def.levels          = 2;    % Number of dual distance estimates. 
                              % Larger values are more accurate but need more time (log-change).
                              % 1 means estimation for .25 and .75 boundaries
                              % Good values are between 2 and 8.
                              %  
  def.range            = .4;  % Default value for range extension for first boundary (should be between 0.2-0.45)
                              %
  def.rangeE           = .4;  % Default value for range extension for second boundary(should be between 0.2-0.45)
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
  def.usemedian        = 2;   % GMT filtering setting (1-median)
  def.medfs            = [1 2]; % [processingResolution filterSize] default [1.5 2]; 
                              %
  def.denoise          = 1;   % correct for noise, artifacts and anatomical issues as WMHs
                              % 
  def.eidist           = 0;   % eikonal distance (not so good yet)

  opt = cat_io_checkinopt(opt,def);
  vxs = mean(vx_vol); 
  opt.vxs = vxs;

  opt.wmnoise = 0.02; 

  c = clock; %#ok<*CLOCK>


% TODO: 
% - remove tic-toc run-time evaluation blocks
  if opt.verb, fprintf('\nPrep:    '); tic; end  
  

  % estimate partial volume effect size in voxel (with fast option)
  opt.pvet0 = estimatePVEsize( Yp0 , 0); 
  opt.pvet  = max(0,min(4 / vxs,opt.pvet0)); 


  % denoising 
  if opt.pvet > 2  &&  mean(vx_vol) < .75  &&  opt.denoise
    [Yp0,opt.pvet] = CS4_desnoise(Yp0,opt);
  end

  if mean(vx_vol) < .75 
    % extend hard cuts, by a low value just to have a broader estimate
    % - could be improved by a basic WMD and GMT estimate to assure a
    %   local minimum thickness related to the GMT distance
    % - RD202601: important for small distances but less for large
    % - Ok for phantom.
    Yp0 = max(Yp0,(Yp0==1 & cat_vol_morph(Yp0==2,'d')) * 1.2);
  
    % close holes 
    % Helpful for SPM with unsufficient WM correction but minor/no-effects for CAT.
    % Ok for phantom.
    Yp0 = max(Yp0, 2.90 * smooth3(cat_vol_morph(Yp0>2.90,'ldc',opt.pvet/1)) ); 
    Yp0 = max(Yp0, 2.75 * smooth3(cat_vol_morph(Yp0>2.75,'ldc',opt.pvet/2)) ); 
    Yp0 = max(Yp0, 2.25 * smooth3(cat_vol_morph(Yp0>2.25,'ldc',0.50,vx_vol)) ); 
  end
  
  % corrections for blood vessels, partial volume effects and myelination
  % RD202503: new blood vessel correction 
  if opt.NBVC, Yp0 = NBVC(Yp0, vx_vol, opt.pvet); end 
  

  % RD2025: correct myelinated areas (using quick thickness estimation)
  if opt.myelinCorrection , Yp0 = myelincorrection(Yp0, vx_vol, opt); end  


  % only within hull
  [Yp0r,resYp0] = cat_vol_resize(Yp0,'reduceV',vx_vol,1,32,'meanm'); % 1mm
  [Ywr,resYp0w] = cat_vol_resize( cat_vol_morph( Yp0r>1.25 ,'do',.5,resYp0.vx_volr) ,'reduceV',resYp0.vx_volr,2,32,'meanm'); clear Yp0r;
  Ywr = cat_vol_morph(Ywr>.5,'ldc',8);
  Ywr = cat_vol_smooth3X(Ywr,1.2); 
  Ywr = cat_vol_resize(Ywr,'dereduceV',resYp0w); 
  Yhull = cat_vol_resize(Ywr,'dereduceV',resYp0); clear Ywr; 
  Yhull = max(0,min(1,atan(Yhull*pi/1.5))); 
  



  %% == Estimation of distance measures ==
  % It is important to use the PVE for multiple (=levels) distance estimations 
  % but also to extened (range) and correct (corpve) the values at the oposite  
  % boundary by the PVE or distance. 
  if 1
    % quick manuell settings
    % at least 1 but use more in low resolution to cope with various settings
    opt.levels = max(1,min(8,2  * mean(vx_vol*2)^2)); 
    opt.rangeE = 0.4; % extimation limit extention - smaller better
    opt.range  = 0.4; % addition range - larger better
    opt.range  = min(.4,opt.range);
  end

  
  %% tested here multple estimations
  if opt.verb, toc; fprintf('Dist:    '); tic; end  
  opt.eidist = 0; % not full ready 
  [Ycd0, Ywd0] = cat_vol_cwdist( Yp0 ,opt);;;


  %% == Estimation of thickness maps == 
  if opt.verb, toc; fprintf('Thick:   '); tic; end      
  if 0 % old version
    Ygmtw0 = cat_vol_pbtp( single(1 + (Yp0>=1.5-opt.rangeE) + (Yp0>2.5+opt.rangeE)) , Ywd0, Ycd0);
    Ygmtc0 = cat_vol_pbtp( single(3 - (Yp0>=1.5-opt.rangeE) - (Yp0>2.5+opt.rangeE)) , Ycd0, Ywd0);
  else
    Ygmtw0 = cat_vol_pbtp2( single(1 + (Yp0>=1.5-opt.rangeE) + (Yp0>2.5+opt.rangeE)) , Ywd0, Ycd0);
    Ygmtc0 = cat_vol_pbtp2( single(3 - (Yp0>=1.5-opt.rangeE) - (Yp0>2.5+opt.rangeE)) , Ycd0, Ywd0);         
  end
  Ygmtw0 = min(Ygmtw0, Ycd0+Ywd0);
  Ygmtc0 = min(Ygmtc0, Ycd0+Ywd0);


  %% extend areas that are too thin for quantification 
  minthick = 0.05/vx_vol; 
  Ycut   = cat_vol_morph(Yp0>2.5,'d',1)  &  cat_vol_morph(Yp0<1.5,'d',1)  &  Ygmtw0 < min(vxs , 1/opt.pvet0)/10 & Ycd0<.125 & Ywd0<.125;
  Ycut   = cat_vol_morph(Ycut,'l',[100 20]) > 0; 
  Ycut   = cat_vol_morph(Ycut,'dd',3); 
  Ycut   = cat_vol_morph(Ycut,'dd',5,vx_vol) & Ygmtw0>0; 
  Ycuts  = cat_vol_localstat( min(Ygmtw0,Ygmtc0) , Ycut, max(1,round(1/vxs)) , 2, 1); 
  Ycuts  = cat_vol_localstat( Ycuts , Ycut, max(1,round(1/vxs)) , 1, 1); 
  Ygmtw0(Ycut) = max(minthick, Ycuts(Ycut));
  Ygmtw0 = cat_vol_median3( Ygmtw0, Ycut ); 
  
  Ypve   = cat_vol_morph(Yp0>2.5,'d',1)  &  cat_vol_morph(Yp0<1.5,'d',1)  &  Ygmtw0 < min(vxs , 1/opt.pvet0)/10 & ~Ycut;
  Ygmtw0(Ypve) = max(minthick, Ygmtw0(Ypve) - mean(Ygmtw0(Ypve(:))) * 2); 
  clear Ycuts; 

  %% low dist correction
  pvefc = 1;
  if 0 %pvefc
    % thickness distribution 
    Ygmt0  = min(Ygmtw0,Ygmtc0); Ygmt0 = Ygmt0(Ygmt0>eps); 
    gmt    = [ cat_stat_nanmedian(Ygmt0) cat_stat_nanstd(Ygmt0) ]; clear Ygmt0;
    % correction term
    pvefc  = min( max(vxs,opt.pvet0) * pvefc , max(vxs, gmt(1) - 3*gmt(2)) );

    %% keep thin sulci 
    Ymsk   = Ygmtw0 > eps & Ygmtw0 < pvefc; 
    Ymsk   = cat_vol_morph(Ymsk,'l',[100 5])>0; 
    Ymsk   = cat_vol_morph(Ymsk,'dd',pvefc,vx_vol); 
    Ygmtw0(Ymsk) = max(minthick,min(pvefc, (Ygmtw0(Ymsk)-(pvefc/2))*2) ); 
    Ygmtw0 = cat_vol_median3(Ygmtw0,Ymsk,Ymsk);

    %% not sure about that
    if 0
      Ymsk   = Ygmtc0 > eps & Ygmtc0 < opt.pvet0*pvefc; 
      Ymsk   = cat_vol_morph(Ymsk,'l',[100 10])>0; 
      Ymsk   = cat_vol_morph(Ymsk,'dd',pvefc,vx_vol); 
      Ygmtc0(Ymsk) = max(minthick,min(pvefc, (Ygmtc0(Ymsk)-(pvefc/2))*2) ); 
      Ygmtc0 = cat_vol_median3(Ygmtc0,Ymsk,Ymsk); 
    end
  end
  
  
  %% final masking (after downsampling)
  th = .01;
  Ymsk = Ycd0<=th | Ywd0<=th | Ygmtw0<=th | Ygmtc0<=th | Yp0>2.5+opt.rangeE | Yp0<1.5-opt.rangeE | Ycd0>1000 | Ygmtw0>1000; 
  Ycd0(Ymsk) = 0; Ywd0(Ymsk) = 0; Ygmtw0(Ymsk) = 0; Ygmtc0(Ymsk) = 0; 

  % correct thickness in narrow areas (about a voxel)
  Ygmtmsk = Ygmtw0>0; 
  Ywm  = Yp0>2.75 & ~Ygmtmsk; 
  Ymsk = Ygmtw0>0 & Ygmtw0<1.25 & Ywd0<1.25 & cat_vol_morph(Ywm,'dc',1) & ~Ywm; clear Ywm;
  Ymsk = cat_vol_morph(Ymsk,'l',[10 .1])>0; 
  Ygmtw0(Ymsk) = max( min(Ygmtw0(Ymsk),minthick) ,Ygmtw0(Ymsk) .* (max(eps,Ygmtw0(Ymsk)./1.25).^2)); clear Ymsk; 
  Ygmtw0 = max( Ygmtmsk*minthick , cat_vol_median3(Ygmtw0,Ygmtmsk,Ygmtmsk)); clear Ygmtmsk; 

  % correction for thin sulci
  Ywd0   = Ywd0   .* ( 1 + (2./(Ygmtw0*vxs)) .* Yhull .* ( (Ygmtw0-Ywd0)./Ygmtw0 ).^(2) .* smooth3(Ygmtw0<Ygmtc0*.95 &  cat_vol_smooth3X(Yp0,4*vxs)>2 & (Ygmtw0*vxs<1.25 | Ygmtw0<1)) ) ; % sulci recon
  Ygmtw0 = Ygmtw0 ./ ( 1 + (2./(Ygmtw0*vxs)) .* Yhull .* ( (Ygmtw0-Ywd0)./Ygmtw0 ).^(2) .* smooth3(Ygmtw0<Ygmtc0*.95 &  cat_vol_smooth3X(Yp0,4*vxs)>2 & (Ygmtw0*vxs<1.25 | Ygmtw0<1)) ) ; % thick update
  Ygmtw0 = max( Ycut * minthick, Ygmtw0);
  Ygmtc0 = max( Ycut * minthick, Ygmtc0);

  % Internal tests
  if 0
    %% median filter
    %  red > 2 create artefacts but this can be used as a feature to see how good the raw data is
    tic; Ygmtw0a = min(Ygmtw0,Ygmtc0); %Ygmtw0a = mixgmt(Ygmtw0,Ygmtc0,Ywd0+Ycd0,.5); 
    Ygmtw0as = Ygmtw0a; red = 1; fs = 2; 
    if 1
      for redi = red, Ygmtw0as = medfilter2(Ygmtw0as, Yp0, redi, fs, 0, vxs); end;;;
    else
      Ygmtw0as = cat_vol_approx(Ygmtw0as .* (Ygmtw0as>0.5),'rec',2); 
    end

    if 1
      %% create histogram
      %  * we need here a very fine sampling (many bins) to see artifacts 
      binsfac = 1000; % .* vxs.^2; 
      Ymsk = Yp0>1.5 & Yp0<2.5 & Ycd0>0 & Ywd0>0 & Ygmtw0<5/vxs & Ygmtw0>0; 
      Ypp  = min(1,max(0, min(Ygmtw0-Ywd0,Ywd0) .* (Ygmtw0>1)) ./ max(eps,Ygmtw0/2));
      Ypp  = min(1,max(0, (Ypp ./ cat_vol_approx(Ypp .* Ypp>.5)).^2 * 1.5 ));
      %
      data = {}; dataname = {}; datacolor = []; 
      %data{end+1} = Ywd0(Ypp(:)<.25 & Ypp(:)>0)*vxs; dataname{end+1} = 'cd0+wd0 (GM)'; datacolor(end+1,:) = [1 0.6 .6];
      %data{end+1} = Ycd0(Ypp(:)>.56)*vxs + Ywd0(Ypp(:)>.56)/vxs; dataname{end+1} = 'cd0+wd0 (CL)'; datacolor(end+1,:) = [1 0.7 1];
      data{end+1} = Ygmtw0a(Ypp(:)>.54)*vxs; dataname{end+1} = 'gmtw0 (CL)';  datacolor(end+1,:) = [0 0.4 0.7];
      %data{end+1} = Ygmtw0(Ypp(:)<.74 & Ypp(:)>.01)*vxs; dataname{end+1} = 'gmtw0 (OL)';  datacolor(end+1,:) = [0.1 0.6 0];
      %data{end+1} = Ygmtc0(Ypp(:)>.73)/vxs; dataname{end+1} = 'gmtc0 (CL)';  datacolor(end+1,:) = [0.7 0.3 0];
      data{end+1} = Ygmtw0as(Ypp(:)>.53)*vxs; dataname{end+1} = 'gmtw0a (CL)'; datacolor(end+1,:) = [.8 0 0];
      %data{end+1} = Ygmt0(Ypp(:)>.73)*vxs; dataname{end+1} = 'gmt0 (CL)';  datacolor(end+1,:) = [0.1 0.6 0];
      %data{end+1} = Ygmt0m(Ypp(:)>.72)*vxs; dataname{end+1} = 'gmt0m (CL)';  datacolor(end+1,:) = [1 0 0];
      cat_plot_histogram(data,struct('color',datacolor,'bins',binsfac,'rawline',0,'xrange',[0 10]));
      xlim([0 round(prctile(Ygmtw0(Ygmtw0(:)>0)*vxs,90)+1)]); ylim([0 50/binsfac]); grid on; 
      ah = gcf; ah.Position(3:4) = [300,200]*1.5;
      title( [char(datetime) sprintf(' - R=%0.1f/r=%0.1f', vxs, opt.range)]); 
      legend(dataname,'Location','Northwest'); 
    end
    if 0
      %% create surface
      Ygmts = Ygmtw0as; 
      Ygmts = smooth3( Ygmts + cat_vol_approx(Ygmts) .* (Ygmts==0) );
      %Ypp   = min(1,max(Yp0>2.5, max(Ycd0,(Ygmts .* (Ywd0>0)) - Ywd0) ./ max(eps,Ygmts))); % basic map without mixing
      Ypp   = min(1,max(Yp0>2.5, max(0,(Ygmts .* (Ywd0>0)) - Ywd0) ./ max(eps,Ygmts))); % basic map without mixing
      CS    = isosurface(smooth3(1-Ypp),.5,min(5,Ygmts*vxs)); cat_surf_render2(CS); 
      title( [char(datetime) sprintf(' - R=%0.1f/r=%0.1f',vxs,opt.range)]);
      if median(Ygmts(Ygmts(:)>.5)*vxs) < 3.5, cat_surf_render2('Colorbar'); cat_surf_render2('clim',[1 4]); end 
    end
  end

  % This is critical for the phantom ! 
  if 1
  % avoid blurring of tiny sulci
    fac          = .5; 
    Ymsk         = Yhull > .5  &  Ygmtw0 < Ygmtc0*0.9  &  (Ygmtw0*vxs < 1.25 | Ygmtw0 < 0.75)  &  (cat_vol_smooth3X(Yp0,2) - Yp0)>2  &  Yhull>.5; 
    Ygmtc0(Ymsk) = max(eps,Ygmtw0(Ymsk) - .25/vxs*fac); 
    Ygmtw0(Ymsk) = max(eps,Ygmtw0(Ymsk) - .50/vxs*fac); 
    Ywd0(Ymsk)   = Ywd0(Ymsk) + .25/vxs*fac; 
    Ycd0(Ymsk)   = min(Ycd0(Ymsk), max(0,Ygmtw0(Ymsk) - Ywd0(Ymsk)));  
    Ygmtw0       = min(Ygmtw0,Ycd0+Ywd0); 
    Ygmtc0       = min(Ygmtc0,Ycd0+Ywd0); 
  end

  if opt.verb, toc; fprintf('MainMed: '); tic; end  

  % minimum tickness map and cleanup (removal of extrem outliers and approximation) 
  Ygmt0 = min(Ygmtw0,Ygmtc0); 
  Ygmt0 = min(Ygmt0,Ycd0+Ywd0);
  Ygmt0 = Ygmt0 .* ( (Ycd0>eps & Ywd0>eps) | Ycut); 
 
  
  if opt.usemedian 
    Ygmt0m = medfilter2(Ygmt0, Yp0, opt.medfs(2), opt.medfs(2)/opt.medfs(1), 0, vxs); 
    Ygmt0  = smooth3( Ygmt0m  +  cat_vol_approx(Ygmt0m,'rec',2) .* (Ygmt0m<.25) ); % rm outliers
  else
    Ygmt0 = cat_vol_approx(Ygmt0,'rec',2);  
  end
  Ycd0  = min(Ygmt0,Ycd0); 
  Ywd0  = min(Ygmt0,Ywd0); % limit
  
  %% low dist correction
  pvefc = 1; 
  if 0 %pvefc
    pvefc = min( max(vxs,opt.pvet0) * pvefc , max(vxs, gmt(1) - 2*gmt(2)) );
    Ymsk  = Ygmtw0 > eps & Ygmtw0 < pvefc; 
    Ymsk  = cat_vol_morph(Ymsk,'l',[100 5])>0; 
    Ymsk  = cat_vol_morph(Ymsk,'dd',pvefc,vx_vol); 
    Ygmt1 = Ygmt0; Ygmt1(Ymsk>0) = max(minthick,min(pvefc,(Ygmt1(Ymsk>0) - ((pvefc)/2))*2 )); 
    Ygmt1 = cat_vol_median3(Ygmt1,Ymsk,Ymsk); 
    Ygmt0 = min(Ygmt0,smooth3(Ygmt1)); clear Ygmt1;
  end
  if opt.verb, toc; end  
  
  

  %% CSF/WM blurring/reconstruction maps
  % - define the blurred sulcal/CSF areas as the area of WM closing and 
  %   'overestimated' thickness (i.e., where the PBT thickness is sign. 
  %   smaller as the simple sum of the WM and CSF distance)
  % - blurred gyri/WM is defined vite versa
  % - undefined values GM values are defined by neighbours
  % - neutral regions are defined by CSF
  % - next both areas are extend by intensity emphasized values 
  if opt.verb, fprintf('  Thickness estimation:           %0.3fs\n', etime(clock,c)); c = clock; end %#ok<*DETIM>
  gmt = [median(Ygmtw0(Ygmtw0(:)>0)) std(Ygmtw0(Ygmtw0(:)>0))]; 
  if opt.gyrusrecon 
    % run on 1 mm should be sufficient
    [Yp0r,resR]       = cat_vol_resize(Yp0,'reduceV',vx_vol,1,16,'meanm'); 
    [Ygmtw0r,Ygmtc0r] = cat_vol_resize({Ygmtw0,Ygmtc0},'reduceV',vx_vol,1,16,'meanm'); 
    
    % low-res voxel-size, average thickness and deviation 
    vxsr = mean( resR.vx_volr ); 
    Ywm  = Yp0r>2.125; 
    Ycsf = Yp0r<1.875;

    % define blurred sulci/gyri by closing of WM and CSF
    Yscr = cat_vol_morph( Ywm  , 'dc', (2*gmt(1) + 1*gmt(2)) * vxs / vxsr ) + Ywm; 
    Ygcr = cat_vol_morph( Ycsf , 'dc', (2*gmt(1) + 1*gmt(2)) * vxs / vxsr ) + Ycsf; 
    
    % extend these regions to avoid WM/CSF bridges (often
    [~,I] = cat_vbdist(single(cat_vol_morph(Yscr,'e')),~Ycsf);  Yscr = single(Yscr(I)); clear I; 
    [~,I] = cat_vbdist(single(cat_vol_morph(Ygcr,'e')),~Ywm);   Ygcr = single(Ygcr(I)); clear I;
    
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
    Ypp  = min(1, max( 0 , max(Yp0>2.5, (Yp0>1.5) .* min(Ycd0,Ygmt0-Ywd0) ./ max(eps,Ygmt0))));
  end 

  
  %% cleanup of position map
  Ypp = max((Yp0-2), max(Ypp, smooth3(Ypp) .* (Ygmt0>gmt(1)*.5) )); % stabize low res fine structures
  Ypp = max(0,min(Yp0-1,Ypp));
  Ypp(cat_vol_morph(Yp0>2.75,'l')) = 1; Ypp(Yp0<1.5) = 0; % close holes
  Ypp = oneObject(Ypp,vx_vol); 
  Ypp = min(Ypp,cat_vol_median3(Ypp,smooth3(Ypp)>0 & smooth3(Ypp)<.75, Ypp>-1, .5)); % prefere open
  Ypp = max(Ypp,cat_vol_median3(Ypp,smooth3(Ypp)>0 & Ygmt0>gmt(1), Ypp>-1, .5));     % prefere close
  % use PVE if better?
  Ypp = max(Ypp, max(0,min(.5,Yp0-2)*2).^2 );  % ^.5 stronger
  Ypp = min(Ypp, min(1,max(0,Yp0-1.5)*2).^2 );  % ^2  stronger

  % final scaling
  Ygmt = Ygmt0 * vxs; 
  
  % evaluation
  if opt.verb
    fprintf('    Median thickness + IQR:        %5.2f ± %4.2f mm\n', ...
      median( Ygmt( Ypp(:)>.3 & Ypp(:)<.7 )) , iqr( Ygmt( Ypp(:)>.3 & Ypp(:)<.7 ))); 
  end
  
  if 0 
    %%
    CS    = isosurface(smooth3(1-Ypp),.5,min(5,Ygmt)); cat_surf_render2(CS); 
    cat_surf_render2('clim',[0 6]); 
  end

end
% ======================================================================
function [Yp0,pvet] = CS4_desnoise(Yp0,opt)
  
  % filter interolation artifacts
  inth  = cat_stat_nanmedian(Yp0(Yp0(:)>.5 & Yp0(:)<1)); 
  Yp0   = min(3,max(1,Yp0)); 
  if inth>.01
    Yp0 = cat_vol_median3(Yp0, Yp0>1.9 & Yp0<2.1,true(size(Yp0)),inth/2);
  end
  
  Yp0om = Yp0; 
  FECth = min(1,max(0,opt.pvet - 2)*2); 
  
  %% main filter depending on WM image _and_ anatomical noise to handle WMHs
  Ymm = cat_vol_morph( Yp0>2 & Yp0<3 , 'd'); 
  Yp0 = Yp0om; 
  for i = 1:opt.wmnoise*100/3 
    Yp0f = Yp0;% .* (Yp0<=2) + cat_vol_localstat(Yp0,Ymm,2, 8 );
    Yp0f = cat_vol_median3(Yp0f, Ymm, Yp0>0, .05 );
    if 1 % refinement
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
  
  % udpate
  pvet = estimatePVEsize( Yp0 , 0); 
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
function Yp0 = NBVC(Yp0,vx_vol,pvet)
%% RD202503: new blood vessel correction 
% In principle this is not really new and should be done by other functions before. 
% However, it is so essential for PBT and not to difficult to include it here (too). 
% First the save WM area is estimated and then extended by a region growing method
% that focus on lower intensities in the Yp0 label mapf for GM and CSF. 
  F      = min(1,max(0.0,Yp0-2)).^.01; F(Yp0<=1.9) = inf; 
  Ywm    = cat_vol_morph(cat_vol_morph(Yp0>2.25,'de',1,vx_vol),'ldo',1,vx_vol) | cat_vol_morph(Yp0>2.75,'l'); 
  % add closest voxels
  [~,Yd] = cat_vol_downcut(single(Ywm), F,0.01); Yd(isnan(Yd))=inf; 
  Ywm    = Yp0>2.1 & Yd<median(Yd(Yp0(:)>=2.25 & Yp0(:)<2.7 & Yd(:)>0));
  
  % update Yp0 by sharpening the WM structures that blurr most
  if 1 % pvet > 2
    Ywmh   = double(cat_vol_morph(Yp0>2.1 & Yd<median(Yd(Yp0(:)>=2.5 & Yp0(:)<2.8 & Yd(:)>0)),'dc',1.5,vx_vol));
    Ywml   = double(Yp0>2.1 & Yd<median(Yd(Yp0(:)>=2.1 & Yp0(:)<2.3 & Yd(:)>0)));
    spm_smooth(Ywmh, Ywmh, double(repmat(.6,1,3)) );
    spm_smooth(Ywml, Ywml, double(repmat(.6,1,3)) ); 
    Yp0    = max(Yp0, Ywmh * 3.0); 
    Yp0    = max(Yp0, Ywml * 2.5); 
  end
  
  %% prepare mask
  [~,Yd] = cat_vol_downcut(single(Ywm), F,-0.001); Yd(isnan(Yd))=inf; %clear Ywm; 
  Ymsk   = Yd > 20000 * mean(vx_vol) & Yp0>2.25 & ~cat_vol_morph(Ywm,'dd',1); 
 
  %% Ymsk as regions that we want to correct 
  %[~,I]  = cat_vbdist(single(~Ymsk),Yp0>1.5); Yp0(Ymsk) = Yp0(I(Ymsk));
  Yp0(Ymsk) = 2; 
  Ymsk      = cat_vol_morph(Ymsk,'dd',1,vx_vol) & Yp0>1.5;
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
  def.eidist = 0; 

  opt = cat_io_checkinopt(opt,def); 

  vxc  = .5 - opt.range; % the more PVE the less we need to correct for
 
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
    Fc = smooth3(max(.001,min(1,cat_vol_smooth3X(min(1,max(eps,2 - Yp0)),2).^4))); 
    Fw = smooth3(max(.001,min(1,cat_vol_smooth3X(min(1,max(eps,Yp0 - 1)),2).^4))); 
  end
  for si = 1:hss
    offset = max(0,min(.5, range * si/(hss+1))); 

    % CSF dist
    if opt.eidist && exist('cat_vbdist3','file')
%      [Ycdl,YIl] = cat_vbdist2(single(Yp0 < ( 1.5 - offset)), Fc , YMC); 
%      [Ycdh,YIh] = cat_vbdist2(single(Yp0 < ( 1.5 + offset)), Fc , YMC); 
      [Ycdl,YIl] = cat_vbdist3(single(Yp0 < ( 1.5 - offset)), Fc .* YMC); 
      [Ycdh,YIh] = cat_vbdist3(single(Yp0 < ( 1.5 + offset)), Fc .* YMC); 
    else
      [Ycdl,YIl] = cat_vbdist(single(Yp0 < ( 1.5 - offset)), YMC ); 
      [Ycdh,YIh] = cat_vbdist(single(Yp0 < ( 1.5 + offset)), YMC ); 
    end
    % PVE correction 
    Ycdl = (Ycdl - min(vxc,max(0,- (Yp0(YIl) - (1.5-offset))*opt.pvet ))) .* (YMM & YMC); 
    Ycdh = (Ycdh - min(vxc,max(0,- (Yp0(YIh) - (1.5+offset))*opt.pvet ))) .* (YMM & YMC); 
    % classic correction (in interaction with second boundary correction!)
   % Ycdl = max(0,Ycdl - 0.5) .* (YMM & YMC);
   % Ycdh = max(0,Ycdh - 0.5) .* (YMM & YMC);
   % Ycdl = min(Ycdl,Ycdh + opt.pvet);
    Ycd  = Ycd + .5/hss .* Ycdl  +  .5/hss .* Ycdh; 

    % WM distances
    if opt.eidist && exist('cat_vbdist3','file')
    %  [Ywdl,YIl] = cat_vbdist2(single(Yp0 > ( 2.5 - offset)), Fw , YMM);
    %  [Ywdh,YIh] = cat_vbdist2(single(Yp0 > ( 2.5 + offset)), Fw , YMM); 
      [Ywdl,YIl] = cat_vbdist3(single(Yp0 > ( 2.5 - offset)), Fw .* YMM);
      [Ywdh,YIh] = cat_vbdist3(single(Yp0 > ( 2.5 + offset)), Fw .* YMM); 
    else
      [Ywdl,YIl] = cat_vbdist(single(Yp0 > ( 2.5 - offset)), YMM );
      [Ywdh,YIh] = cat_vbdist(single(Yp0 > ( 2.5 + offset)), YMM ); 
    end
    % PVE correction 
    Ywdl = (Ywdl - min(vxc,max(0, (Yp0(YIl) - (2.5-offset))*opt.pvet ))) .* (YMM & YMC); 
    Ywdh = (Ywdh - min(vxc,max(0, (Yp0(YIh) - (2.5+offset))*opt.pvet ))) .* (YMM & YMC); 
  %  Ywdl = max(0,Ywdl - 0.5) .* (YMM & YMC);
  %  Ywdh = max(0,Ywdh - 0.5) .* (YMM & YMC);
  %  Ywdh = min(Ywdh,Ywdl + opt.pvet);
    Ywd  = Ywd + .5/hss .* Ywdl  +  .5/hss .* Ywdh;
  end
  Ycd = max(0,Ycd);
  Ywd = max(0,Ywd);


  %% endpoint PVE correction 
  if 0
  % it would be possible to also correct for not achived points but this might increase the risk for defects
  % ... worse than thickness
    Ycdc = min(.5,max(-.5, (Yp0>2.5-0*opt.rangeE & Yp0<2.5+opt.rangeE) .* (Yp0-2.5) ));
    Ywdc = min(.5,max(-.5, (Yp0>1.5-opt.rangeE & Yp0<1.5+0*opt.rangeE) .* (Yp0-1.5) ));
    Ycdc(Ycdc~=0) = Ycdc(Ycdc~=0) ./ median(abs(Ycdc(Ycdc(:)~=0))) * opt.rangeE;
    Ywdc(Ywdc~=0) = Ywdc(Ywdc~=0) ./ median(abs(Ywdc(Ywdc(:)~=0))) * opt.rangeE;
    Ycd = max(0,Ycd - max(0,min( 1,Ycdc)));
    Ywd = max(0,Ywd + min(0,max(-1,Ywdc)));
  end
  if 1
    for si = 1:hss
      offset = max(0,min(.5, range * si/(hss+1))); 

      % CSF dist
      if opt.eidist 
        [Ycdl,YIl] = cat_vbdist2(single(Yp0 < ( 2.5 - offset)), Fc , YMC); 
        [Ycdh,YIh] = cat_vbdist2(single(Yp0 < ( 2.5 + offset)), Fc , YMC); 
      else
        [Ycdl,YIl] = cat_vbdist(single(Yp0 < ( 2.5 - offset)), YMC ); 
        [Ycdh,YIh] = cat_vbdist(single(Yp0 < ( 2.5 + offset)), YMC ); 
      end
      %Ycdl = (Ycdl - min(vxc,max(0,- (Yp0(YIl) - (2.5-offset))*opt.pvet ))) .* (YMM & YMC); 
      %Ycdh = (Ycdh - min(vxc,max(0,- (Yp0(YIh) - (2.5+offset))*opt.pvet ))) .* (YMM & YMC); 
      Ycdl = max(0,Ycdl - 0.5) .* (YMM & YMC);
      Ycdh = max(0,Ycdh - 0.5) .* (YMM & YMC);
      Ycd  = max(0,Ycd - (.5/hss .* Ycdl  +  .5/hss .* Ycdh)); 
  
      % WM distances
      if opt.eidist
        [Ywdl,YIl] = cat_vbdist2(single(Yp0 > ( 1.5 - offset)), Fw , YMM);
        [Ywdh,YIh] = cat_vbdist2(single(Yp0 > ( 1.5 + offset)), Fw , YMM); 
      else
        [Ywdl,YIl] = cat_vbdist(single(Yp0 > ( 1.5 - offset)), YMM );
        [Ywdh,YIh] = cat_vbdist(single(Yp0 > ( 1.5 + offset)), YMM ); 
      end
      %Ywdl = (Ywdl - min(vxc,max(0, (Yp0(YIl) - (1.5-offset))*opt.pvet ))) .* (YMM & YMC); 
      %Ywdh = (Ywdh - min(vxc,max(0, (Yp0(YIh) - (1.5+offset))*opt.pvet ))) .* (YMM & YMC); 
      Ywdl = max(0,Ywdl - 0.5) .* (YMM & YMC);
      Ywdh = max(0,Ywdh - 0.5) .* (YMM & YMC);
      Ywd  = max(0,Ywd - (.5/hss .* Ywdl  +  .5/hss .* Ywdh));
    end
  end
end
% ======================================================================
function Yp0 = myelincorrection(Yp0,vx_vol,opt)
  if opt.myelinCorrection > 0 %&& mean(vx_vol(:))<1
    % quick estimation of the cortical thickness
    opt.verb      = 0; 
    opt.levels    = 2; 
    opt.usemedian = 0; 
    opt.denoise   = 0; 
    opt.NBVC      = 0;
    %opt.eidist    = 0; 

    [Ycd, Ywd] = cat_vol_cwdist(Yp0, opt);
  
    % projection-based thickness mapping
    Ygmt0 = cat_vol_pbtp2( round(Yp0) , Ywd, Ycd);
    Ygmt0 = cat_vol_approx(Ygmt0,'rec',2); 
  
    % reestimation of the CSF distance 
    Ypp   = min(1,min(Ygmt0,Ycd) ./ max(eps,Ygmt0)); Ypp(Yp0>2.5 & Ypp==0) = 1; 
    Ycdc2 = cat_vbdist( single( max(Yp0<=1, 1 - Ycd - Ypp) ), true(size(Ycd)) );
    Ycdc2(Ycdc2 > 6 / mean(vx_vol)) = 0; 
   
    % estimate the full tissue thickness (we needed the GM thickness and WM to reconstruct the sulcus)
    Ybmt  = cat_vol_pbtp2( min(3,4 - min(2,Yp0)), Ycdc2, Ycdc2*inf); 
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

    % smooth resolution dependency
   %orgth = max(0,min(1,mean(vx_vol(:))*2 - 1));
   % Yp0   = Yp0.*orgth + (1-orgth).*Yp0c;
  end
end
% ======================================================================
function Ygmt = medfilter2( Ygmtw , Yp0, red, fs, mask, vxs)
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

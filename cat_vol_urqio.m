function cat_vol_urqio(job)
% Ultrahigh Resolution Quantitative Image Optimization
% ______________________________________________________________________
%
% Skript to reduce inhomogeneity, noise, and blood vessels in ultra-high
% resolution quantitative MR data. The R1, PD (A), and R2s images are 
% required. The function use a simple tissue classification to apply a
% bias correction, a gobal intensity normalization, a blood vessel 
% correction, and finally a noise reduction. 
% 
% WARNING: This function is in an early development stage (alpha)
%
% cat_vol_urqio(job)
%
%  job.data .. structure with the filenames of the input images
%   .r1   .. R1 images
%   .pd   .. PD images (A)
%   .r2s  .. R2s images
%  job.ouput .. structure to control writing of output images
%   .pd  .. write corrected PD images
%   .r1  .. write corrected R1 images
%   .r2s .. write corrected R2s images
%   .t1  .. write synthetic T1 images (invertation of the PD)
%   .bv  .. write detected blood vessels
%  job.opts .. structure of option parameter
%   .prefix .. filename prefix (default = 1)
%   .verb   .. display processing notes
%   .bc     .. apply bias correction [0|1]
%   .in     .. apply global intensity normalisation [0|1]
%   .bvc    .. apply blood vessel correction [0|1]
%   .nc     .. apply noise correction
%
% ______________________________________________________________________
% Robert Dahnke
% Structural Brain Mapping Group
% University Jena
%
% $Id: cat_io_writenii.m 1159 2017-08-04 14:08:14Z dahnke $


  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

  % default options
  if ~exist('job','var'), job = struct(); end
  % input data
  def.data.r1   = {};
  def.data.pd   = {};
  def.data.r2s  = {};
  % output data
  def.output.pd  = 1; 
  def.output.t1  = 1; 
  def.output.r1  = 1;
  def.output.r2s = 1; 
  def.output.bv  = 1; 
  % parameter
  def.opts.prefix = 'catsyn_';
  def.opts.verb = 2;
  def.opts.bc   = 1; 
  def.opts.in   = 1; 
  def.opts.bvc  = 1; 
  def.opts.nc   = 1; 
  % checkin
  job = cat_io_checkinopt(job,def); 
  
  % set empty data fields
  if isempty(job.data.r1) || isempty(job.data.r1{1})
    job.data.r1 = cellstr(spm_select(1,'image','Select R1 files'));
  end
  if isempty(job.data.pd) || isempty(job.data.pd{1})
    job.data.pd = cellstr(spm_select(1,'image','Select PD files'));
  end
  if  isempty(job.data.r2s) || isempty(job.data.r2s{1})
    job.data.r2s = cellstr(spm_select(1,'image','Select R2s files'));
  end
  job.data.r1   = cellstr(job.data.r1); 
  job.data.r2s  = cellstr(job.data.r2s); 
  job.data.pd   = cellstr(job.data.pd); 
  
  si = 1; %#ok<NASGU>
  
  
  
  %% main loop
  for si=1:numel(job.data.r1)
    stime2 = cat_io_cmd(sprintf('Preprocessing subject %d/%d:',si,numel(job.data.r1)),'','',job.opts.verb); fprintf('\n'); 
    stime  = cat_io_cmd('  Load images:','g5','',job.opts.verb-1);

    % get header
    Vd  = spm_vol(job.data.r1{si});
    Vr  = spm_vol(job.data.r2s{si}); 
    Vp  = spm_vol(job.data.pd{si});

    % voxel size
    vx_vol = sqrt(sum(Vd.mat(1:3,1:3).^2));

    % load images
    Yd = single(spm_read_vols(Vd));
    Yr = single(spm_read_vols(Vr));
    Yp = single(spm_read_vols(Vp));
   
    % some kind of brain masking
    Ypg  = cat_vol_grad(Yp)./Yp;
    Ydg  = cat_vol_grad(Yd)./Yd;
    pgth = max( mean(Ypg(Ypg(:)>0 & Ypg(:)<0.5)) , abs(mean(Ypg(Ypg(:)>0 & Ypg(:)<0.5))) );
    Ytis = Ypg>0 & Ypg<pgth*2 & Ydg<pgth*2 & Ydg>0; 
    Ytis = cat_vol_morph(cat_vol_morph(Ytis,'lo',1),'lc',10/mean(vx_vol)) & Ypg>0 & Ypg<pgth*4 & Ydg<pgth*4 & Ydg>0;
    
    % some simple thresholds 
    dth = max( mean(Yd(Ytis(:))) , abs(mean(Yd(Ytis(:)))));
    rth = max( mean(Yr(Ytis(:))) , abs(mean(Yr(Ytis(:)))));
    pth = max( mean(Yp(Ytis(:))) , abs(mean(Yp(Ytis(:)))));
    if ~debug, clear Ytis; end
    
    %  gradient maps and thresholds
    %  only use PD contrast that have best SNR here
    Ypg  = cat_vol_grad(Yp);
    Ypd  = cat_vol_div(Yp);
    pgth = max( mean(Ypg(Ypg(:)>0)) , abs(mean(Ypg(Ypg(:)<0))));
    
    
        
    
    %% tissue classification
    %  --------------------------------------------------------------------
    %  Yw = white matter map
    %  Yg = gray matter map
    %  Yc = cerebrospinal fluid map
    %  --------------------------------------------------------------------
    stime  = cat_io_cmd('  Tissue classification:','g5','',job.opts.verb-1,stime);
    
    % WM
    Yw = Yp<pth*1 & Yp>pth*0.6  & Yd>dth & Yd<dth*3.5  &  Yr>rth/2 & Yr<rth*2.5  & ...
         Ypg<pgth*4 & Ypd>-pgth/2 & Ypd<pgth/2; 
    Yw = Yw & Yr<mean(Yr(Yw(:)))*1.5; % remove
    Yw(smooth3(Yw)<0.4)=0; % small open
    Yw = cat_vol_morph(Yw,'l');
    
    % GM
    Yg = Yp>pth & Yp<pth*2  &  Yd<dth*1.6 & Yd>dth*0.5  &  Yr>rth/2 & Yr<rth*2.5  & ...
         Ypg<pgth*2 & Ypd>-pgth/2 & Ypd<pgth/2  & ~Yw; 
    Yg(smooth3(Yg)<0.4)=0; % small open 
    
    % CSF ... have to use it ...
    Yc = Yp>pth*1.5 & Yp<pth*3  &  Yd<dth*0.5 & Yd>dth*0.1  &  Yr<rth/2 & Yr<rth*2.5  & ...
         Ypg<pgth*8 & Ypd>-pgth*2  & ~Yw & ~Yg; 
    Yc(smooth3(Yc)<0.4)=0; % small open 
    
    % super WM (higher intensity in R1 and R2s) and the question if and how 
    % we need to correction it for normal preprocessing, because actual it 
    % lead to GM like brainstem that get lost ...
    %if debug, Yx = Yp; end  
    if ~debug, clear Ypd; end
    
    if 0
      %% display segmentation V1
      ds('l2','',vx_vol,Yp./mean(Yp(Yw(:))),(Yw*3 + Yg*2 + Yc)/3,Yp./mean(Yp(Yw(:))),(Yw*3 + Yg*2 + Yc)/3*2,40); colormap gray
      %% display segmentation V2
      ds('l2','',vx_vol,Yd/dth,Yw,Yr./rth,Yd/dth,200); colormap gray
    end
    
        
    
    
    %% bias correction
    %  --------------------------------------------------------------------
    %  Ybiw  = initial signal intensity map (of the white matter)
    %  Ybig  = intitla signal intensity map of the gray matter
    %  Ybiws = approximated bias field map
    %          normalization in the next step
    %  --------------------------------------------------------------------
    if job.opts.bc
      useGMandCSF = [1 0]; % use also the GM (and CSF) for bias correction (CSF is problematic) 

      % approx .. 2 is fast and smooth, 1 sh
      biasstr = { % name var WMfilter CSFfitler approxres approxsmooth
        'PD' ,'Yp', 2, 3, 1.2, 2/job.opts.bc;
        'R1' ,'Yd', 3, 2, 1.2, 2/job.opts.bc;
        'R2s','Yr', 3, 2, 1.2, 1/job.opts.bc;
      };
      for bi=1:size(biasstr,1)
        stime  = cat_io_cmd(sprintf('  Bias correction (%s):',biasstr{bi,1}),'g5','',job.opts.verb-1,stime);
        eval(['Ybi = ' biasstr{bi,2} '; bith = ' biasstr{bi,2}(2) 'th;']); 
        Ybiw = Ybi .* Yw; Ybiw = cat_vol_localstat(Ybiw,Ybiw>0,1,1); 
        Ybiw = Ybi .* (Yw & Ybiw~=Ybi & abs(Ybiw-Ybi)<bith*0.2); 
        Ybiw = cat_vol_localstat(Ybiw,Ybiw>0,1,biasstr{bi,3}); 
        for ii=1:round(8/mean(vx_vol)), Ybiw = cat_vol_localstat(Ybiw,Ybiw>0,1,1); end
        if useGMandCSF(1)
          % use additonal GM
          Ybig = Ybi .* Yg; Ybig = cat_vol_localstat(Ybig,Ybig>0,1,1); 
          Ybig = Ybi .* (Yg & Ybig~=Ybi & abs(Ybig-Ybi)<bith*0.1); 
          for ii=1:round(8/mean(vx_vol)), Ybig = cat_vol_localstat(Ybig,Ybig>0,1,1); end
          Ybiw = Ybiw + (Ybig / median(Ybig(Ybig(:)>0)) * median(Ybiw(Ybiw(:)>0)));
          clear Ybig; 
        end
        if useGMandCSF(2)
          % use additonal CSF
          Ybig = Ybi .* Yc; Ybig = cat_vol_localstat(Ybig,Ybig>0,1,2); 
          Ybig = Ybi .* (Yc & Ybig~=Ybi & abs(Ybig-Ybi)<bith*0.1); 
          Ybiw = cat_vol_localstat(Ybiw,Ybiw>0,1,biasstr{bi,4}); 
          for ii=1:round(8/mean(vx_vol)), Ybig = cat_vol_localstat(Ybig,Ybig>0,1,2); end
          Ybiw = Ybiw + (Ybig / median(Ybig(Ybig(:)>0)) * median(Ybiw(Ybiw(:)>0)));
          clear Ybig; 
        end
        Ybiws = cat_vol_approx(Ybiw,'nn',vx_vol,biasstr{bi,5},struct('lfO',biasstr{bi,6}));  %#ok<NASGU>
        eval([biasstr{bi,2} 'ws = Ybiws;']);
        if debug, eval([biasstr{bi,2} 'w = Ybiw;']); end
        clear Ybiw Ybiws;    
      end
      if 0
        % this is just for manual development & debugging
        %% display PD
        ds('d2','',vx_vol,Yp./mean(Yp(Yw(:))),Ypws./mean(Yp(Yw(:))),Yp./mean(Yp(Yw(:))),Yp./Ypws,100); colormap gray
        %% display R1
        ds('d2','',vx_vol,Yd./mean(Yd(Yw(:))),Ydws./mean(Yd(Yw(:))),Yd./mean(Yd(Yw(:))),Yd./Ydws,200); colormap gray
        %% display R2s
        ds('d2','',vx_vol,Yr./mean(Yr(Yw(:))),Yrws./mean(Yr(Yw(:))),Yr./mean(Yr(Yw(:))),Yr./Yrws,200); colormap gray
      end
    end
    
    
    
    
    %% intensity normalization
    %  --------------------------------------------------------------------
    %  Ydm  = bias corrected R1 map
    %  Ypm  = bias corrected PD map
    %  Yrm  = bias corrected R2s map
    %
    %  Ydmc = intensity normalized bias corrected R1 map
    %  Ypmc = intensity normalized bias corrected PD map
    %  Yrmc = intensity normalized bias corrected R2s map
    %  --------------------------------------------------------------------
    if job.opts.in
      stime  = cat_io_cmd('  Intensity Normalization:','g5','',job.opts.verb-1,stime);

      % apply bias correction, create inverted version of the pd image
      if job.opts.bc
        Ydm = Yd ./ Ydws; 
        Ypm = Yp ./ Ypws; 
        Yrm = Yr ./ Yrws; 
      else
        Ydm = Yd ./ mean(Yd(Yw(:))); 
        Ypm = Yp ./ mean(Yp(Yw(:))); 
        Yrm = Yr ./ mean(Yr(Yw(:))); 
      end
      if ~debug, clear Yd Yp Yr; end

      %% create synthetic T1 contrast based on the PD image
      Ytm = 1.5 - Ypm/2; 
      YM  = smooth3(~Yw & Ypm<0.8 & cat_vol_morph(cat_vol_morph(Ypm<0.5 | Yrm>1.5,'l'),'d',1.5))>0.6; 
      YM  = YM | smooth3( ~Yw & (abs(Ydm-Yrm)>1.1 | Yrm>2.1) )>0.7;  
      Ytm(YM)=0; 
      if debug, clear YM; end 

      %%  the region maps are only used to estimate a the mean intensity of a tissue 
      Ywm = Yw & Yrm<mean(Yrm(Yw(:)))*1.5 & Ypg<pgth*0.5;
      Ygm = Yg & Yrm<mean(Yrm(Yw(:)))*1.5 & Ypg<pgth; 
      Ycm = Yc & Ypg<pgth;
      Yhm = ~Yw & ~Yg & ~Yc & Ydm>1.2 & Ypm<0.5 & Yrm>1.2;
      Ylm = ~Yw & ~Yg & ~Yc & Ydm>0 & Ypm>0 & Yrm>0 & Ydm<mean(Ydm(Ycm(:)))*1.5 & Ypm>mean(Ypm(Ycm(:)))/1.5 & Yrm<mean(Yrm(Ycm(:)))*1.5;
      if sum(Ylm(:))<10
        Ylm = ~Yw & ~Yg & ~Yc & Ydm>min(Ydm(:)) & Ypm>min(Ypm(:)) & Yrm>min(Yrm(:)) & Ydm<mean(Ydm(Ycm(:)))*2 & Ypm>mean(Ypm(Ycm(:)))/2 & Yrm<mean(Yrm(Ycm(:)))*2;
      end  
      if ~debug, clear Yg Yc; end

      % real thresholds [bg lm cm gm wm hm max] 
      T3th(1:4,1) = {'Ydm';'Ypm';'Yrm';'Ytm'};
      T3th{1,2} = [ 0 min(Ydm(Ylm(:))) mean(Ydm(Ycm(:))) mean(Ydm(Ygm(:))) mean(Ydm(Ywm(:))) mean(Ydm(Yhm(:))) max(Ydm(:))];
      T3th{2,2} = [ 0 mean(Ypm(Ycm(:))) mean(Ypm(Ygm(:))) mean(Ypm(Ywm(:))) mean(Ypm(Yhm(:))) max(Ypm(:))];
      T3th{3,2} = [ 0 min(Yrm(Ylm(:))) mean(Yrm(Ycm(:))) mean(Yrm(Ygm(:))) mean(Yrm(Ywm(:))) mean(Yrm(Yhm(:))) max(Yrm(:))];
      T3th{4,2} = [ 0 min(Ytm(Ylm(:))) mean(Ytm(Ycm(:))) mean(Ytm(Ygm(:))) mean(Ytm(Ywm(:))) mean(Ytm(Yhm(:))) max(Ytm(:))];
      if ~debug, clear Ylm Yhm Ycm Ygm Ywm; end
      
      %% final thresholds
      T3th{1,3} = [ 0 1/12 1/3 2/3 1   T3th{1,2}(5)+diff(T3th{1,2}(5:6)) T3th{1,2}(6)+diff(T3th{1,2}(6:7))]; 
      T3th{2,3} = [ 0 1 2/3 1/3 T3th{2,2}(5)/(1/T3th{2,2}(2)) T3th{2,2}(6)/(1/T3th{2,2}(2))]; 
      T3th{3,3} = [ 0 1/12 1/3 2/3 1   T3th{3,2}(5)+diff(T3th{3,2}(5:6)) T3th{3,2}(6)+diff(T3th{3,2}(6:7))]; 
      T3th{4,3} = [ 0 1/12 1/3 2/3 1   T3th{4,2}(5)+diff(T3th{4,2}(5:6)) T3th{4,2}(6)+diff(T3th{4,2}(6:7))]; 

    
      
      % global intensity normalization 
      for bi=1:size(T3th,1)
        %%
        [T3ths,tsi] = sort(T3th{bi,2});
        T3thx       = T3th{bi,3}(tsi);
        eval(['Ym = ' T3th{bi,1} ' + 0; if ~debug, clear ' T3th{bi,1} '; end; Ysrc = Ym;']); 
        for i=numel(T3ths):-1:2
          M = Ysrc>T3ths(i-1) & Ysrc<=T3ths(i);
          Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3ths(i-1))/diff(T3ths(i-1:i))*diff(T3thx(i-1:i));
        end
        M  = Ysrc>=T3ths(end); 
        Ym(M(:)) = numel(T3ths)/6 + (Ysrc(M(:)) - T3ths(i))/diff(T3ths(end-1:end))*diff(T3thx(i-1:i));    
        eval([T3th{bi,1} 'c = Ym;']); 
        clear Ym Ysrc M;
      end
    else
      % no intensity normalization 
      Ydmc = Yd ./ mean(Yd(Yw(:))); 
      Ypmc = Yp ./ mean(Yp(Yw(:))); 
      Yrmc = Yr ./ mean(Yr(Yw(:))); 
      Ytmc = 1.5 - Ypmc/2; Ytmc(Ypmc<0.8 & cat_vol_morph(cat_vol_morph(Ypmc<0.5 | Yrmc>1.5,'l'),'d',1.5))=0;  
    end
    
    
    
    
    %% find and correct blood vessels
    %  --------------------------------------------------------------------
    %  Yrd, Ydd, Ytd = divergence maps
    %  --------------------------------------------------------------------
    if job.opts.bvc || job.output.bv 
      stime  = cat_io_cmd('  Blood Vessel Detection:','g5','',job.opts.verb-1,stime);

      % divergence update
      Yrd = cat_vol_div(Yrmc);
      Ydd = cat_vol_div(Ydmc);
      Ytd = cat_vol_div(Ytmc);


      %% first we remove high intensity blood vessels with large area 
      %  in the save WM region only a mean filter should be used to create an
      %  corrected map, whereas in the other regions a larger and minimum
      %  based filter is applied
      Yb    = cat_vol_morph(Ydmc>0 & Ydmc<1.1 & Yrmc>0 & Yrmc<1.1 & Ytmc>0 & Ytmc<1.1,'lc',2);  % brain mask for CSF minimum value
      Ywc   = cat_vol_morph(Yw,'lc',1) | ~cat_vol_morph(~cat_vol_morph(Yw,'o',1),'lc',2);       % wm mask with close microstructures
      Ywc   = Ywc | smooth3(Yrmc>0.8 & Yrmc<1.1)>0.5| cat_vol_morph(smooth3(Ydmc>0.7 & Ydmc<1.2)>0.6,'o',1);
      Ywc   = cat_vol_morph(Ywc,'lc',1);
     
      %% correction maps    
      Ymin  = max(Yb.*3/12,min(Ydmc,min(Ytmc,Yrmc)));   % minimum intensity in all images (CSF)
      Ymax  = max(Yb.*3/12,max(Ydmc,max(Ytmc,Yrmc)));   % maximum intensity in all images (BV)
      Ymn1  = cat_vol_localstat(Ymin,Ymin>0,1,1);       % low    short  correction
      Ymin1 = cat_vol_localstat(Ymin,Ymin>0,1,2);       % storng short  correction
      Ymax1 = cat_vol_localstat(Ymax,Ymax>0,1,3);       % very close to vessel
      Ymax2 = cat_vol_localstat(Ymax1,Ymax>0,1,3);      % close to vessel
      Ymin2 = cat_vol_localstat(Ymin1,Ymin>0,1,2);      % strong medium correction
      Ymin4 = cat_vol_localstat(Ymin2,Ymin>0,2,2);      % strong long   correction


      %% specific masks 
      %  - save regions that need no corrections
      %  - save blood vessels that requires maxim correction (long distance with min)
      %  - save blood vessels that requires carfull correction (short distance with mean) 
      %  - small blood vessels that requires less agressive correction 
      YminA = Ytmc  .* (                                     Ymax<1.1 &  Ywc) + ... low  values in the WM 
              Ymn1  .* (                         Ymax>=1.1 &             Ywc) + ... high values in the WM 
              Ymin  .* (                         Ymax< 0.7 &            ~Ywc) + ... low  good values
              Ymn1  .* (                         Ymax>=0.7 & Ymax<1.0 & ~Ywc) + ... affected voxel 
              Ymn1  .* (                         Ymax>=1.0 & Ymax<1.1 & ~Ywc) + ...
              Ymin  .* (Ymin< 0.9 & Ymax2< 1.3 & Ymax>=1.1 & Ymax<1.3 & ~Ywc) + ...
              Ymin1 .* (Ymin< 0.9 & Ymax2>=1.3 & Ymax>=1.1 & Ymax<1.3 & ~Ywc) + ...
              Ymin  .* (Ymin< 0.9 & Ymax2< 1.3 & Ymax>=1.3            & ~Ywc) + ...
              Ymin2 .* (Ymin< 0.9 & Ymax2>=1.3 & Ymax>=1.3            & ~Ywc) + ...
              Ymn1  .* (Ymin>=0.9              & Ymax>=1.1 & Ymax<1.3 & ~Ywc) + ...
              Ymin4 .* (Ymin>=0.9              & Ymax>=1.3            & ~Ywc);
      

      %% find blood vessels
      Ybv1 = (Yrmc>2 | Ytmc>1.5 | Ydmc>1.5) & ~Ywc; 
      Ybv  = max(0 , max(Yrmc.*Ydmc.*Ytmc,Ymax) - min(Yrmc.*Ydmc.*Ytmc,Ymin) - 2/3 + 2*max(Yrd.*Ydd.*Ytd,max(Yrd,max(Ydd,Ytd)))  - Ywc );
      
      if debug 
        %% display BV
        ds('d2','',vx_vol,Ydmc*1.5,Ytmc*1.5,Ybv,Ybv1,220); colormap gray
      end
    end
    if job.opts.bvc
      stime  = cat_io_cmd('  Blood Vessel Correction:','g5','',job.opts.verb-1,stime);

      %% first correction 
      Ybvd = cat_vol_div(Ybv);
      Ylow = smooth3(((Ymax1+1)./(Ymin1+1))>(Ymin+1) | Ybv>2 | Ybv1)>0.2;
      %%
      %Ypc  = min(1.1.*(Ypm>0.1),max( (1.25-Ymin2) .* (Ypm>0),Ypmc - (Ylow | Ybv1) .* ~Yw .* (min(0,Ybvd) - max(0,Ybv-0.5) - max(0,(Ydmc./Ypmc) - 1.5) )));
      Ydc  = max(Ymin2,Ydmc + (Ylow | Ybv1) .* ~Yw .* (min(0,Ybvd) - max(0,Ybv-0.5) - max(0,(Ydmc./Ypmc) - 1.5) ));
      Ytc  = max(Ymin2,Ytmc + (Ylow | Ybv1) .* ~Yw .* (min(0,Ybvd) - max(0,Ybv-0.5) - max(0,(Ydmc./Ypmc) - 1.5) ));
      Yrc  = max(Ymin2,Yrmc + (Ylow | Ybv1) .* ~Yw .* (min(0,Ybvd) - max(0,Ybv-0.5) - max(0,(Ydmc./Ypmc) - 1.5) ));

      %% median fitlering
      %Ydcm = cat_vol_median3(Ydc,Ydc>0,true(size(Ydc)),0.05); 
      Ytcm = cat_vol_median3(Ytc,Ytc>0,true(size(Ydc)),0.05);
      Yrcm = cat_vol_median3(Yrc,Yrc>0,true(size(Ydc)),0.05);

      %% second correction
      % r1
      YminAm = cat_vol_median3(YminA,YminA>0,true(size(Ydc)),0.1); 
      YM   = YminAm>0.2 & ~Ywc & ((Ydmc - YminAm)>0.2 | Ydmc>1.5 | Ybv>1); 
      Ydc2 = Ydmc .* (YminA>0.05); Ydc2(YM) = YminAm(YM); Ydc2 = min(Ydc2,Ydmc); 
      %% pd
      YM   = Ypmc>0.01 & Yb & ~Ywc & smooth3((Ypmc - (1.25-YminAm))>0.1 | Ypmc<0.3 | Ybv>1)>0.1; 
      YH   = 1.1 - 0.05*rand(size(YM)); 
      Ypc2 = Ypmc; Ypc2(YM) = min(YH(YM),1.25-YminAm(YM).*Yb(YM)); 
      % t1 = inverse pd
      YM   = Ytmc>0.2 & ~Ywc & ((Ytmc - Ytcm)>0.1 | Ytmc>1.1 | Ybv>0.5); 
      Ytc2 = Ytmc; Ytc2(YM) = Ytcm(YM); Ytc2 = min(Ytc2,Ytmc);
      YM   = (Ydc2 - Ytc2)<-0.4 & Ydc2<0.5; 
      Ytc2(YM) = min(Ytc2(YM) ,Ydc2(YM) );
      % r2s
      YM   = Yrmc>0.2 & ~Ywc & ((Yrmc - Yrcm)>0.2 | Yrmc>1.1 | Ybv>0.5); 
      Yrc2 = Yrmc; Yrc2(YM) = Yrcm(YM); Yrc2 = min(Yrc2,Yrmc);
      YM   = smooth3(Ydc2 - Yrc2)<-0.2 & (Ydc2<0.8 | Yrc2>1.3) & (~Ywc | (Yrd<-0.1 & Yrc2>1.3 & Ypg>0.05))>0.5; 
      Yrc2(YM) = min(Yrc2(YM) ,Ydc2(YM) );
      Yrc2 = cat_vol_median3(Yrc2,YM | Yrc2>1.3 & Ypg>0.2); % remove small reminding WM vessels and artifacts
      if ~debug, clear Ypg; end
      
      %% this is just for manual development & debugging
      if 0
        
        %% display BV
        ds('d2','',vx_vol,Ybv1,Ybv,Ydcm,abs(Ydc2-Yd),200); colormap gray
        %% display PD
        ds('d2','',vx_vol,Ypm/2.5*1.5,Ypws./mean(Yp(Yw(:))),Ypmc*1.5,Ypc2*1.5,220); colormap gray
        %% display PDi
        ds('d2','',vx_vol,2 - Ypc2*1.5,Yp/mean(Yp(Yw(:))),Ydc2*1.5,Ytc2*1.5,220); colormap gray
        %% display R1
        ds('d2','',vx_vol,Yd./mean(Yd(Yw(:)))*1.5,Ydws./mean(Yd(Yw(:))),Yd./mean(Yd(Yw(:)))*1.5,Ydc2*1.5,220); colormap gray
        %% display R2s
        ds('d2','',vx_vol,Yr./mean(Yr(Yw(:))),Yrws./mean(Yr(Yw(:))),Yr./mean(Yr(Yw(:))),Yrc2,220); colormap gray
        
        
        %% create test WM surface with correct blood vessel system - Warning very slow! 
        S   = isosurface(cat_vol_smooth3X(Ytc2,0.5),0.8); 
        Sb  = isosurface(cat_vol_smooth3X(Ybv .* Yb,0.5),0.25); 
        
        % reduce surface
        SR  = reducepatch(S,500000);
        SbR = reducepatch(Sb,500000);
        
        % add some color
        SR.facevertexcdata  = zeros(size(SR.vertices,1),1); 
        SbR.facevertexcdata = ones(size(SbR.vertices,1),1); 
        
        % combine surfaces
        SRX.vertices        = [SR.vertices;SbR.vertices];
        SRX.faces           = [SR.faces;SbR.faces+size(SR.vertices,1)];
        SRX.cdata = [SR.facevertexcdata;SbR.facevertexcdata];

        % dispaly surfaces
        cat_surf_render2(SRX);
        
      end
    else
      Ypc2 = Ypmc;
      Ytc2 = Ytmc;
      Ydc2 = Ydmc;
      Yrc2 = Yrmc;
    end
    if ~debug, clear Ypmc Ytmc Ydmc Yrmc; end
    
    
    
    
    %% write data
    %  --------------------------------------------------------------------
    stime  = cat_io_cmd('  Write Output:','g5','',job.opts.verb-1,stime);
    
    % PD 
    if job.output.pd  
      [pp,ff,ee,dd] = spm_fileparts(Vp.fname);
      Vpc = Vp;  Vpc.fname = fullfile(pp,[job.opts.prefix 'pd_' ff ee dd]); Vpc.dt(1) = 16; 
      spm_write_vol(Vpc,Ypc2); 
    end
    
    % inverted PD
    if job.output.t1 
      [pp,ff,ee,dd] = spm_fileparts(Vp.fname);
      Vtc = Vp;  Vtc.fname = fullfile(pp,[job.opts.prefix 't1_' ff ee dd]); Vtc.dt(1) = 16; 
      spm_write_vol(Vtc,Ytc2); 
    end
    
    % R1
    if job.output.r1% r1
      [pp,ff,ee,dd] = spm_fileparts(Vd.fname);
      Vdc = Vd;  Vdc.fname = fullfile(pp,[job.opts.prefix 'r1_' ff ee dd]); Vdc.dt(1) = 16; 
      spm_write_vol(Vdc,Ydc2); 
    end
    
    % R2s
    if job.output.r2s
      [pp,ff,ee,dd] = spm_fileparts(Vr.fname);
      Vrc = Vr;  Vrc.fname = fullfile(pp,[job.opts.prefix 'r2s_' ff ee dd]); Vrc.dt(1) = 16; 
      spm_write_vol(Vrc,Yrc2); 
    end
    
    % BV
    if job.output.bv 
      [pp,ff,ee,dd] = spm_fileparts(Vp.fname);
      Vbvc = Vp;  Vbvc.fname = fullfile(pp,[job.opts.prefix 'bv_' ff ee dd]); Vbvc.dt(1) = 16; 
      spm_write_vol(Vbvc,Ybv .* Yb); 
    end    
    if ~debug, clear Ypc2 Ytc2 Ydc2 Yrc2 Ybv Yb; end
    
    
    
    %% noise reduction
    if job.opts.nc
      stime  = cat_io_cmd('  Noise Correction:','g5','',job.opts.verb-1,stime); fprintf('\n')
      if job.output.pd,  cat_vol_sanlm(struct('data',cellstr(Vpc.fname),'prefix','')); end
      if job.output.t1,  cat_vol_sanlm(struct('data',cellstr(Vtc.fname),'prefix','')); end
      if job.output.r1,  cat_vol_sanlm(struct('data',cellstr(Vdc.fname),'prefix','')); end
      if job.output.r2s, cat_vol_sanlm(struct('data',cellstr(Vrc.fname),'prefix','')); end
      if job.output.bv,  cat_vol_sanlm(struct('data',cellstr(Vbvc.fname),'prefix','')); end
    end
    
    cat_io_cmd(' ','g5','',job.opts.verb-1); cat_io_cmd(' ','g5','',job.opts.verb-1,stime); % time of the last step
    cat_io_cmd(' ','','',job.opts.verb-1,stime2); fprintf('\n'); % time of the full subject
  end  
end
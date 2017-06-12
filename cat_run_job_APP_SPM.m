function  Ymc = cat_run_job_APP_SPM(Po,vout,vx_vol,verb,version,fwhm)
%  _____________________________________________________________________
%  The final bias correction is a subfunction of cat_run_job.
% 
%  The affine registration, especially spm_preproc8 requires a very good 
%  masking! Because this is also required for the Unified Segmenation
%  a wider mask with a complete brain is important.
%  _____________________________________________________________________
%  Robert Dahnke
%  $Id$


  % if there is a breakpoint in this file set debug=1 and do not clear temporary variables 
  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end
  
  resMth  = 1.4; % general resolution limitation 
  resTth  = 1.4; % tissue intensity estimation resolution limitation 
  
  if ~exist('verb','var'),     verb = 1; end
  if ~exist('version','var'),  version = 1; end
  if ~exist('fwhm','var'),     fwhm = 0.5; end
  fwhmx   = max(4,2 + 4*fwhm); % only lower filtering 
    
  %version = 2;
  if version == 2
    stime = cat_io_cmd('    Initialize','g5','',verb);
    Vo   = spm_vol(Po);
    Yo   = single(spm_read_vols(Vo));
    Ym   = vout.Ym; 
    Ycls = vout.Ycls; 
    res  = vout.res; 
    if ~debug; clear vout; end 
    
    % general resolution limitation 
    [Ym,resT1] = cat_vol_resize(Ym,'reduceV',vx_vol,resMth,32,'meanm');
    Yo         = cat_vol_resize(Yo,'reduceV',vx_vol,resMth,32,'meanm');
    for i=1:6, Ycls{i} = cat_vol_resize(Ycls{i},'reduceV',vx_vol,resMth,32); end
    vx_vol = resT1.vx_volr;
    
    Yp0  = single(Ycls{1})/255*2 + single(Ycls{2})/255*3 + single(Ycls{3})/255;
    Yb   = cat_vol_morph(cat_vol_morph(Yp0>1,'lo'),'lc',2); 
    
    
    %% global intensity normalization 
    %  ds('l2','',0.5,Ym,Yb,Ym./mean(Ym(Yp0(:)==3)),Yp0/3,80)
    %  ds('l2','',0.5,Ymi,Yp0,Ym/T3th(3),Ymi,80)
    if any( min(vx_vol*2,resTth)./vx_vol >= 2 )
      Ymr = cat_vol_resize(Ym,'reduceV',vx_vol,min(vx_vol*2,resTth),32,'meanm');
      Ybr   = cat_vol_resize(single(Yb),'reduceV',vx_vol,min(vx_vol*2,resTth),32,'meanm')>0.5;
      Yclsr = cell(size(Ycls)); for i=1:6, Yclsr{i} = cat_vol_resize(Ycls{i},'reduceV',vx_vol,min(vx_vol*2,resTth),32); end
      [Ymr,Ybr,T3th,Txth,inv_weighting,noise] = cat_main_gintnorm(Ymr,Yclsr,Ybr,vx_vol,res);
      clear Ymr Ybr Yclsr; 
      Ymi = cat_main_gintnorm(Ym,Txth); 
    else
      [Ymi,Yb,T3th,Tth,inv_weighting,noise] = cat_main_gintnorm(Ym,Ycls,Yb,vx_vol,res);
    end
    
    
    %% prepare data
    stime = cat_io_cmd('    Prepare measures','g5','',verb,stime);
    Yg    = cat_vol_grad(Ymi,vx_vol); 
    Ydiv  = cat_vol_div(Ymi,vx_vol/2);

    %% tissues
    stime = cat_io_cmd('    Optimize tissues','g5','',verb,stime);
    % WM
    Ywm   = cat_vol_morph(cat_vol_morph(Yp0>2,'l'),'lc');                  % for minimum thickness
    Ycsfd = cat_vbdist(single(Yp0<1.75 & ~Ywm),true(size(Yp0)),vx_vol);
    Ywm   = Yb & (Ywm | cat_vol_morph(Ymi>.9 & Ymi<1.2 & Yp0>2 & Ycsfd>1.5 & ...
            (Yp0<1.8 | Yp0>2.2) ,'lc')) & Yg<0.5 & abs(Ydiv)<0.5;
    Ywm   = Ywm & ~(Ymi<0.95 & Ydiv>0.02);                                 % avoid subcortical structures
    %% GM
    Ygm  = ~Ywm & Ycls{1}>8 & Ycls{2}<192 & Yg<0.6 & abs(Ydiv)<0.5 & Ymi<1.1 & Ymi>0.5;
    Ygmi = cat_vol_localstat(Ygm.*Ymi,Ygm,2,1);
    Ygm(abs(((Ygm.*Ymi)./Ygmi)-Ygm)>0.15 | abs(((Ygm.*Ymi)./Ygmi)-Ygm)==0)=0;
    Ygms = smooth3(Ygm); Ygm(Ygms<0.2)=0; Ygm(Ygms>0.6)=1; 
    %% CSF
    Ycm = cat_vol_morph(cat_vol_morph(Yp0>0.2,'lo'),'lc',1) & ~Ywm & ~Ygm & Ycls{3}>128 & Yg<0.3 & Ydiv>0 & Ymi>0 & Ymi<0.4; 
    Ycm = Ycm | (Ycls{3}>192 & Yg<0.3);
    Ycms = smooth3(Ycm); Ycm(Ycms<0.4)=0; Ycm(Ycms>0.6)=1; 
    Ycm = cat_vol_morph(Ycm,'l',[20 0.05])>0;
    
    %% background
    stime = cat_io_cmd('    Estimate background','g5','',verb,stime);
    Ybg   = Yg<(mean(Yg(Ymi(:)~=0))) | isnan(Yg);  
    % avoid edges on the image border
    bb = 1; Ybb1 = true(size(Ybg)); Ybb1(bb+1:end-bb,bb+1:end-bb,bb+1:end-bb) = 0; 
    bb = 4; Ybb2 = true(size(Ybg)); Ybb2(bb+1:end-bb,bb+1:end-bb,bb+1:end-bb) = 0; 
    Ybgc  = Ybb2 & cat_vol_morph(Ybg | (Ybb1 & Yg<0.5),'c',2);
    Ybg = Ybg | Ybgc | smooth3(Yg./Ydiv > 1000)>0.5; % >100
    Ybg = Ybg & Ycls{6}>128;
    if ~debug, clear Ybb1 Ybb2 bb; end
    % filling
    [Ybg,resT2] = cat_vol_resize(single(~Ybg),'reduceV',vx_vol,2.4,32,'meanm'); 
    Ybg  = Ybg>0.5;
    Ybg  = cat_vol_morph(Ybg,'lc',4);
    Ybg  = cat_vol_smooth3X(Ybg,1); 
    Ybg  = cat_vol_resize(Ybg,'dereduceV',resT2)<0.5;    
    zeroBG = cat_stat_nanmean(Ymi(Ybg(:)>0))<0.2;

    %% HD
    Yhm = Yg>0.01 & Yg<(mean(Yg(Yp0(:)>2.8))*2) | isnan(Yg);
    Yp0d = cat_vbdist(single(Yp0>0.1),~Ybg,vx_vol); 
    Ypd = cat_vbdist(single(Ybg),Yp0==0,vx_vol); 
    Yhm = Yhm & Ypd>10 & cat_vbdist(single(Yp0>0.1),~Ybg,vx_vol)>6 & Yp0<0.1 & Ymi<1.5 & abs(Ydiv)<0.1 & Ycls{5}>8; 
    for i=1:2, Yhds = smooth3(Yhm); Yhm(Yhds<0.4)=0; Yhm(Yhds>0.6)=1; end
    Yhm = cat_vol_morph(Yhm,'l',[100 0.05])>0;
    
    %%
    Yhh = (Ycls{5}>128 | Ydiv<-0.3) & ~Yb & Ydiv<0 & (Ymi>0.8 | Ydiv<-0.3) & ...
      (Yg>0.3 | Ymi>1) & ~Yhm & Ycls{2}<8 & Ycls{3}<8 & Ydiv<0.6; 
    Yhh = cat_vol_morph(Yhh,'l',[20 0.05])>0;
    Ywmi = Ywm.*Yo; Yhhi = Yhh.*Yo; Yhhi = cat_vol_localstat(Yhhi,Yhh,2,3);
    Yhhi = cat_vol_approx(Ywmi + Yhhi * ...
      (cat_stat_nanmean(Ywmi(Ywmi(:)>0))/cat_stat_nanmedian(Ym(Yhhi(:)>0 & Yp0d(:)<10))),'nearest',vx_vol,4); 
    Yhhi = cat_vol_smooth3X(Yhhi,4); 
    
    
    %%
    stime = cat_io_cmd('    Smooth values','g5','',verb,stime);
    if inv_weighting 
      Ywmi = cat_vol_localstat(Ywm.*Yo,Ywm,1,2);
    else
      Ywmi = cat_vol_localstat(Ywm.*Yo,Ywm,1,3);
    end
    for i=1:1, Ywmi = cat_vol_localstat(Ywmi,Ywm,1,1); end
    Ygmi = Ygm.*Yo; for i=1:1, Ygmi = cat_vol_localstat(Ygmi,Ygm,1,1); end
    Ycmi = Ycm.*Yo; for i=1:1, Ycmi = cat_vol_localstat(Ycmi,Ycm,1,1); end
    Yhmi = Yhm.*Yo; for i=1:1, Yhmi = cat_vol_localstat(Yhmi,Yhm,1,1); end
    Ycmi = cat_vol_approx(Ycmi,'nearest',vx_vol,4); Ycmi = cat_vol_smooth3X(Ycmi,4); 
    Yhmi = cat_vol_approx(Yhmi,'nearest',vx_vol,4); Yhmi = cat_vol_smooth3X(Yhmi,4); 
   
    %%
    stime = cat_io_cmd('    Estimate bias field','g5','',verb,stime);
    Yhmith = (cat_stat_nanmean(Ywmi(Ywm(:)))/ max( cat_stat_nanmean(Ygmi(Ygmi(:)>0)) , min( cat_stat_nanmean(Ywmi(Ywmi(:)>0)) ,...
             cat_stat_nanmedian(Yhmi(Yhmi(:)>0 & Ycls{5}(:)>225 & Yp0d(:)<15)))));
    Ytmi = Ywmi + ...
           Ygmi * (cat_stat_nanmean(Ywmi(Ywm(:)))/cat_stat_nanmean(Ygmi(Ygmi(:)>0))) + ...
           Ycm  .* Ycmi * (cat_stat_nanmean(Ywmi(Ywm(:)))/cat_stat_nanmedian(Ycmi(Ywm(:)))) + ...
           Yhm  .* Yhmi * Yhmith + ...
           (Yhh  & ~Yhm & ~Ywm & ~Ygm & ~Ycm) .* Yhhi * (cat_stat_nanmedian(Ywmi(Ywm(:)))/cat_stat_nanmean(Yhhi(Ywm(:)))) + ...
           (Ybgc & ~Yhm & ~Ywm & ~Ygm & ~Ycm & ~Yhh) .* (Ycmi * (cat_stat_nanmean(Ywmi(Ywm(:)))/cat_stat_nanmedian(Ycmi(Ywm(:)))) + ...
                    Yhmi * Yhmith + ...
                    Yhhi * (cat_stat_nanmedian(Ywmi(Ywm(:)))/cat_stat_nanmean(Yhhi(Ywm(:)))))/3;
    for i=1:1, Ytmi = cat_vol_localstat(Ytmi,Ytmi>0,1,1); end
    Ywi = cat_vol_approx(Ytmi,'nearest',vx_vol,2);
    Ywi = cat_vol_smooth3X(Ywi,4); 

    if debug
      %%
      Ymc = Yo ./ Ywi; Ymc = cat_main_gintnorm(Ymc * (cat_stat_nanmean(Ym(Ywm(:)))/cat_stat_nanmean(Ymc(Ywm(:)))) ,Txth); 
      %%
     % ds('l2','',0.5,Ymi,Yhm,Ymi,Ymc,80)
    end
    %%
    Ywi  = cat_vol_resize(Ywi,'dereduceV',resT1); 
    Ywm  = cat_vol_resize(single(Ywm),'dereduceV',resT1)>0.5; 
    Yo   = single(spm_read_vols(Vo));
    

    Ymc = Yo ./ Ywi; Ymc = Ymc * (cat_stat_nanmean(Ym(Ywm(:)))/cat_stat_nanmean(Ymc(Ywm(:))));
    cat_io_cmd(' ','','',verb,stime); 
    
    
    
    
    
  elseif version == 1
    % ds('l2','m',0.5,Ym*0.7+0.3,Yb,Ysrc/WMth,Ym,80)
    Vo = spm_vol(Po);
    Yo = single(spm_read_vols(Vo));
    Ym = vout.Ym; 
    Ycls = vout.Ycls; 
    res  = vout.res; 
    if ~debug; clear vout; end 

    %% prepare data
    stime = cat_io_cmd('    Prepare measures','g5','',verb);
    
    Yp0  = single(Ycls{1})/255*2 + single(Ycls{2})/255*3 + single(Ycls{3})/255;
    Yb   = cat_vol_morph(cat_vol_morph(Yp0>1,'lo'),'lc',1); 

    Ymr = cat_vol_resize(Ym,'reduceV',vx_vol,min(vx_vol*2,resTth),32,'meanm');
    Ybr   = cat_vol_resize(single(Yb),'reduceV',vx_vol,min(vx_vol*2,resTth),32,'meanm')>0.5;
    Yclsr = cell(size(Ycls)); for i=1:6, Yclsr{i} = cat_vol_resize(Ycls{i},'reduceV',vx_vol,min(vx_vol*2,resTth),32); end
    [Ymr,Ybr,T3th,Txth,inv_weighting,noise,cat_warnings] = cat_main_gintnorm(Ymr,Yclsr,Ybr,vx_vol,res,res.Twarp);
    clear Ymr Ybr Ysrcr Yclsr; 
    Ymi = cat_main_gintnorm(Ym,Txth); 

    %Yb   = cat_vol_morph(cat_vol_morph(Yp0>1,'lo'),'lc',1); 
    Yb2  = cat_vol_morph(cat_vol_morph(Yp0>2,'l'),'lc'); 
    Yg   = cat_vol_grad(Ymi); 
    Ydiv = cat_vol_div(Ymi);
    Ycsfd = cat_vbdist(single(Yp0<1.75 & ~Yb2))/mean(vx_vol);

    %% tissues
    stime = cat_io_cmd('    Optimize tissues','g5','',verb,stime);
    Ywm = Yb & (Yb2 | cat_vol_morph(Ymi>.9 & Ymi<1.2 & Yp0>2 & (Yp0<1.8 | Yp0>2.2) ,'lc')) & Yg<0.5 & abs(Ydiv)<0.5 & Ycsfd>2;
    Ywm = cat_vol_morph(Ywm,'l',[3 0.3])>0; Ywm(smooth3(Ywm)<0.5)=0;
    Ywm = Ywm & ~(Ymi<0.95 & Ydiv>0.02); 
    Ygm = Yb & ~Ywm & Yp0>1.5 & Yp0<2.5 & Yg<0.4 & Ydiv<0.1 & Ydiv>-0.1 & Ymi<0.9; 
    Ygmi = cat_vol_localstat(Ygm.*Ymi,Ygm,1,1);
    Ygm(abs(((Ygm.*Ymi)./Ygmi)-Ygm)>0.1 | abs(((Ygm.*Ymi)./Ygmi)-Ygm)==0)=0;
    Ycm = Yb & ~Ywm & ~Ygm & Yp0>0.5 & Yp0<1.25 & Yg<0.3; 
    %Yhm = Yg<0.2 & ~Yb & abs(Ydiv)<0.1 & Ym>Tth*0.1  & Ym<Tth*1.5;

    %
    stime = cat_io_cmd('    Smooth values','g5','',verb,stime);
     if inv_weighting 
      Ywmi = cat_vol_localstat(Ywm.*Yo,Ywm,1,2);
    else
      Ywmi = cat_vol_localstat(Ywm.*Yo,Ywm,1,3);
     end
    Ygmi = cat_vol_localstat(Ygm.*Yo,Ygm,1,1); 
    Ycmi = cat_vol_localstat(Ycm.*Yo,Ycm,1,1);
    %Yhmi = cat_vol_localstat(Yhm.*Yo,Yhm,1,1);

    %
    stime = cat_io_cmd('    Estimate bias field','g5','',verb,stime);
    Ytmi = Ywmi + Ygmi*(mean(Ywmi(Ywmi(:)>0))/mean(Ygmi(Ygmi(:)>0))) + Ycmi*(mean(Ywmi(Ywmi(:)>0))/cat_stat_nanmedian(Ycmi(Ycmi(:)>0 & Yg(:)<0.1)));
    for i=1:2, Ytmi = cat_vol_localstat(Ytmi,Ytmi>0,2,1); end
    Ywi = cat_vol_approx(Ytmi,'nearest',vx_vol,2);
    Ywi = cat_vol_smooth3X(Ywi,fwhmx); 

    %%
    %Ymc = Yo ./ Ywi; Ymc = cat_main_gintnorm(Ymc * (mean(Ym(Ywm(:)))/mean(Ymc(Ywm(:)))) ,Txth); 
    Ymc = Yo ./ Ywi; Ymc = Ymc * (mean(Ym(Ywm(:)))/mean(Ymc(Ywm(:)))); 

    cat_io_cmd(' ','','',verb,stime); 
  end
end


function  Ym = cat_run_job_APP_old()
  % skull-stripping
  Yb  = cat_vol_morph(cat_vol_morph(Yp0>1,'lo'),'lc',1); 

  %% background
  stime = cat_io_cmd('  Estimate background','g5','',verb);
  Ybg   = Yg<(mean(Yg(Ym(:)~=0))) | isnan(Yg);  
  % avoid edges on the image border
  bb = 1; Ybb1 = true(size(Ybg)); Ybb1(bb+1:end-bb,bb+1:end-bb,bb+1:end-bb) = 0; 
  bb = 4; Ybb2 = true(size(Ybg)); Ybb2(bb+1:end-bb,bb+1:end-bb,bb+1:end-bb) = 0; 
  Ybgc  = Ybb2 & cat_vol_morph(Ybg | (Ybb1 & Yg<0.5),'c',2);
  Ybg = Ybg | Ybgc | smooth3(Yg./Ydiv > 100)>0.5;
  if ~debug, clear Ybb1 Ybb2 bb; end
  % filling
  [Ybg,resT2] = cat_vol_resize(single(~Ybg),'reduceV',resT1.vx_volr,2,32,'meanm'); 
  Ybg  = Ybg>0.5;
  Ybg  = cat_vol_morph(Ybg,'lc',4);
  Ybg  = cat_vol_smooth3X(Ybg,2); 
  Ybg  = cat_vol_resize(Ybg,'dereduceV',resT2)<0.5;    
  zeroBG = cat_stat_nanmean(Ym(Ybg(:)>0))<0.2;
  
  
  if verb, fprintf('\n'); end
  stime = cat_io_cmd('  Initialize','g5','',verb,stime);
  msize = 222; %round(222 ./ max(size(Ysrco).*vx_vol) .* min(size(Ysrco).*vx_vol));  
  
  [Ysrc,Ym,resT3] = cat_vol_resize({Ysrco,Ym},'reduceV',vx_vol,min(1.5,min(vx_vol)*2),msize,'meanm'); 
  [Yb,Ybg]        = cat_vol_resize({single(Yb),single(Ybg)},'reduceV',vx_vol,min(1.5,min(vx_vol)*2),msize,'meanm'); 
  if debug, Ybo = Yb; end %#ok<NASGU>
  Ybg = Ybg>0.5;
  
  Yg   = cat_vol_grad(Ym,resT3.vx_volr) ./ max(eps,Ym); 
  Ydiv = cat_vol_div(Ym,resT3.vx_volr) ./ Ym;
  Ygs  = smooth3(Yg);
  
  % greater mask - distance based brain radius (brad)
  [dilmsk,resT2] = cat_vol_resize(Yb,'reduceV',resT3.vx_volr,mean(resT3.vx_volr)*2,32); 
  dilmsk  = cat_vbdist(dilmsk,true(size(dilmsk)))*mean(resT2.vx_volr); %resT2.vx_volr);
  dilmsk  = dilmsk - cat_vbdist(single(dilmsk>0),true(size(dilmsk)))*mean(resT2.vx_volr); %,resT2.vx_volr);
  dilmsk  = cat_vol_resize(smooth3(dilmsk),'dereduceV',resT2); 
  brad    = -min(dilmsk(:));
  dilmsk  = dilmsk / brad; 

  voli  = @(v) (v ./ (pi * 4./3)).^(1/3);                        % volume > radius
  brad  = double(mean([brad,voli(sum(Yb(:)>0).*prod(vx_vol))])); % distance and volume based brain radius (brad)
  
  
  % thresholds
  rf   = 6; 
  Hth  = round2(cat_stat_nanmean(Ym(Ym(:)>0.4 & Ym(:)<1.2  & Ygs(:)<0.2 & ~Yb(:) & Ydiv(:)<0.05 & Ydiv(:)>-0.5 & dilmsk(:)>0 & dilmsk(:)<10)),rf); % average intensity of major head tissues
  if isnan(Hth), Hth = 0.8; end
  GMth = round2(cat_stat_nanmean(Ym(Ym(:)>0.2  & Ym(:)<0.9      & Ygs(:)<0.2 & Yb(:) & Ydiv(:)<0.1 & Ydiv(:)>-0.1)),rf);  % first guess of the GM intensity
  CMth = round2(cat_stat_nanmean(Ym(Ym(:)>0.05 & Ym(:)<GMth*0.5 & Ygs(:)<0.2 & Yb(:) & Ydiv(:)>-0.10)),rf);  % first guess of the CSF intensity
  %WMth = cat_stat_nanmean(Ym(Ym(:)>0.8 & Ym(:)<1.2 & Ygs(:)<0.2 & ~Yb(:) & Ydiv(:)>-0.05)); 
  BGth = round2(cat_stat_nanmean(Ym(Ybg(:))),rf); 
  if isnan(CMth), CMth=mean([BGth,GMth]); end
  
  
  
  Ymo = Ym; 

  %% ---------------------------------------------------------
  %  improve bias correction:
  %  Also here it is important that the bias field is very smooth
  %  to avoid overcorrections. In contrast to the first
  %  correction we will try to separate between tissues...
  %  ---------------------------------------------------------
  stime = cat_io_cmd('  Tissue classification','g5','',verb,stime);
  Ym   = Ymo;
  Yg   = cat_vol_grad(Ym,vx_vol);%./max(eps,Ym)
  Ygs  = smooth3(Yg);
  Ydiv = cat_vol_div(Ym,vx_vol);

% WM Skeleton:  Ydiv./Yg<-1    
% nicht WM: Ydiv.*Yg<-0.02

  %% tissue classes WM, GM, subcortical GM (Ybm), CSF, head tissue (Yhm) 
  Ywm  = ((Ym-max(0,(Ydiv+0.01))*10)>(GMth*0.1+0.9) | Ym>(GMth*0.1+0.9) | ...
          (Ydiv./Yg<0.5 & ((Ym>0.9 & Yg>0.1 & Ydiv<0) | (Ym>0.9)) & Ym<1.2) ) & ...
          Ym<1.3 & Yg<0.6 & Ygs<0.9 & Yb & Ydiv<0.1 & Ydiv>-0.5; 
  CSFD = cat_vbdist(single(Ym<(CMth*0.5+0.5*GMth)),Yb,vx_vol);
  Ywm  = Ywm & CSFD>2;
  Ywm  = smooth3(Ywm)>0.5;      
  % subcotical GM 
  Ybm  = ((Ym-max(0,(Ydiv+0.01))*10)<0.98) & Ym<0.98 & dilmsk<-brad*0.3 & ... 
         Ym>(GMth*0.6+0.4) & Yb & Ygs<0.2 & Yg<0.2 & ... Ydiv<0.1 & Ydiv>-0.02 & ... Yg<0.1 & 
         ~(Ydiv./Yg<0.5 & ((Ym>0.9 & Yg>0.1 & Ydiv<0) | (Ym>0.95)) & Ym<1.2);
       %& ~Ywm;  
  Ybm  = smooth3(Ybm)>0.5;
  Ybm  = cat_vol_morph(Ybm,'o',1);
  % cortical GM 
  Ygm  = Ym<(GMth*0.3+0.7) & Ym>(CMth*0.6+0.4*GMth) & Yg<0.4 & Yb & Ydiv<0.4 & Ydiv>-0.3 & ~Ywm & ~Ybm; % & (Ym-Ydiv*2)<GMth;  
  Ygm(smooth3(Ygm)<0.3 | ~cat_vol_morph(Ywm,'d',3/mean(vx_vol)))=0;
  Ygm(CSFD<3 & Ym>(CMth*0.5+0.5*GMth) & ~Ywm & Ym<(CMth*0.5+0.5*GMth))=0; 
  % CSF
  Ycm  = Ym<(CMth*0.5+0.5*GMth) & Yg<0.1 & Yb & ~Ygm & dilmsk<-brad*0.3; 
  Ycm  = smooth3(Ycm)>0.5;
  % head tissue
  Yhm  = Ym>max(mean([CMth,GMth]),Hth*0.2) & Ym<1.2 & Yg<0.8 & cat_vol_smooth3X(Yb,2)<0.1 & Ydiv<0.6 & Ydiv>-0.6;
  Yhm  = smooth3(Yhm)>0.5;

  %% refine
  %Ygm  = Ygm | (cat_vol_morph(Ywm,'d',3) & ~Ybm & ~Ywm & (Ym-Ydiv*2)<GMth & ...
  %   ~Ycm & smooth3((Ym + Yg)<(CMth*0.8+0.2*GMth))<0.5) & Ym>(CMth*0.9+0.1*GMth) & Ym<(GMth*0.2+0.8);;
  
  %% masking of the original values and local filtering
  stime = cat_io_cmd('  Filtering','g5','',verb,stime);
  fi   = round2(max(3,min(resT3.vx_volr)*3)/3); 
  Ywm  = Ysrc .* Ywm; Ywm  = cat_vol_localstat(Ywm,Ywm>0,1,3); % PVE
  Ycm  = Ysrc .* Ycm; Ycm  = cat_vol_localstat(Ycm,Ycm>0,1,2);
  for i=1:fi-1, Ywm = cat_vol_localstat(Ywm,Ywm>0,2,1); end
  for i=1:fi-1, Ycm = cat_vol_localstat(Ycm,Ycm>0,2,1); end
  Ybm  = Ysrc .* Ybm; for i=1:fi, Ybm = cat_vol_localstat(Ybm,Ybm>0,2,1); end
  Ygm  = Ysrc .* Ygm; for i=1:fi, Ygm = cat_vol_localstat(Ygm,Ygm>0,2,1); end
  Yhm  = Ysrc .* Yhm; for i=1:fi, Yhm = cat_vol_localstat(Yhm,Yhm>0,2,1); end

  % estimate intensity difference bettween the tissues
  Ywmr = cat_vol_resize(Ywm,'reduceV',resT3.vx_volr,8,16,'meanm'); 
  Ybmr = cat_vol_resize(Ybm,'reduceV',resT3.vx_volr,8,16,'meanm'); 
  Ygmr = cat_vol_resize(Ygm,'reduceV',resT3.vx_volr,8,16,'meanm'); 
  Ycmr = cat_vol_resize(Ycm,'reduceV',resT3.vx_volr,8,16,'meanm'); 
  bmth = mean(Ybmr(Ybmr(:)>0 & Ywmr(:)>0))/mean(Ywmr(Ybmr(:)>0 & Ywmr(:)>0));
  gmth = mean(Ygmr(Ygmr(:)>0 & Ywmr(:)>0))/mean(Ywmr(Ygmr(:)>0 & Ywmr(:)>0));
  cmth = mean(Ycmr(Ycmr(:)>0 & Ywmr(:)>0))/mean(Ywmr(Ycmr(:)>0 & Ywmr(:)>0));
  Ywmr = cat_vol_resize(Ywm,'reduceV',resT3.vx_volr,16,8,'meanm'); 
  Yhmr = cat_vol_resize(Yhm,'reduceV',resT3.vx_volr,16,8,'meanm'); 
  hmth = mean(Yhmr(Yhmr(:)>0 & Ywmr(:)>0))/mean(Ywmr(Ywmr(:)>0 & Ywmr(:)>0));
  hmth = min(max(hmth,mean([cmth,gmth])),mean([gmth,1]));
  % if something failed use global thresholds
  if isnan(bmth), bmth = GMth; end
  if isnan(gmth), gmth = GMth; end
  if isnan(cmth), cmth = CMth; end
  if isnan(hmth), hmth = GMth; end
  clear Ywmr Ybmr Ygmr Yhmr Ycmr; 

  Yhm = Yhm .* (dilmsk>20);  % to avoid near skull tissue
  
  %% estimate bias fields
  stime = cat_io_cmd('  Bias correction','g5','',verb,stime);
  Ywi = sum( cat(4,Ywm,Ygm/gmth,Ybm/bmth,Ycm/cmth,Yhm/hmth),4) ./ sum( cat(4,Ywm>0,Ybm>0,Ygm>0,Ycm>0,Yhm>0),4 );
  if ~zeroBG
    Ybg2 = Ybg(:) & Yg(:)<(cat_stat_nanmean(Yg(Ybg(:))) + 2*(cat_stat_nanstd(Yg(Ybg(:))))); 
    Ywi(Ybg2) = Ysrc(Ybg2); clear Ybg2;
  end
  %%
  [Ywi,resT2]  = cat_vol_resize(Ywi,'reduceV',resT3.vx_volr,min(4,min(resT3.vx_volr)*2),32,'meanm'); 
  for i=1:4, Ywi=cat_vol_localstat(Ywi,Ywi>0,2,1); end
  Ywi   = cat_vol_approx(Ywi,'nn',resT2.vx_volr,2);
  Ywi   = cat_vol_smooth3X(Ywi,4.*mean(vx_vol)); 
  Ywi   = cat_vol_resize(Ywi,'dereduceV',resT2);

  %% background noise
  if zeroBG
    stime = cat_io_cmd('  Background correction','g5','',verb,stime);
    %Ybc  = cat_vol_morph(smooth3(Ym<mean([BGth,CMth]) & Ym<CMth & Ygs<0.05 & ~Yb & dilmsk2>8)>0.5,'lo',3); 
    [Ybc,resT2] = cat_vol_resize(Ysrc .* Ybg,'reduceV',resT2.vx_volr,max(8,max(16,cat_stat_nanmean(resT2.vx_volr)*4)),16,'min'); 
    Ybc  = cat_vol_localstat(Ybc,Ybc>0,2,2);
    for i=1:1, Ybc  = cat_vol_localstat(Ybc,Ybc>0,2,1); end
    %Ybc2 = cat_vol_approx(Ybc,'nn',resT2.vx_volr,4); % no aproximation to correct only in the backgound!  
    %Ybc2 = cat_vol_smooth3X(Ybc2,4);
    Ybc  = cat_vol_smooth3X(Ybc,2);
    Ybc  = cat_vol_resize(Ybc,'dereduceV',resT2); 
    %Ybc2 = cat_vol_resize(Ybc2,'dereduceV',resT2); 
  else
    % correction for negative backgrounds (MT weighting)
    [x,y]=hist(Ysrc(:),200); cx = cumsum(x)/sum(x);
    Ybc = y(find(cx<0.0001,1,'last')); 
  end

  %% back to original size
  stime = cat_io_cmd('  Final scaling','g5','',verb,stime);
  Ywi   = cat_vol_resize(Ywi,'dereduceV',resT3); 
  if zeroBG, Ybc = cat_vol_resize(Ybc,'dereduceV',resT3); end
  %Ybc2 = cat_vol_resize({Ybc2},'dereduceV',resT3); 
  [Yg,Ygs]  = cat_vol_resize({Yg,Ygs},'dereduceV',resT3); 
  Yb   = cat_vol_resize(Yb,'dereduceV',resT3)>0.5; 
  Yp0 = cat_vol_resize(((Ywm>0)*3 + (Ygm>0)*2 + (Ybm>0)*2.3 + (Ycm>0)*1 + (Yhm>0)*2.7 + (Ybg>0)*0.5)/3,'dereduceV',resT3);
  Ysrc = Ysrco; clear Ysrco;

  %%  Final intensity scaling
  Ym   = (Ysrc - Ybc) ./ (Ywi - Ybc); % correct for noise only in background
 % Ym   = (Ysrc - Ybc) ./ (Ywi - Ybc2 + Ybc); % correct for noise only in background
  Wth  = single(cat_stat_nanmedian(Ym(Ygs(:)<0.2 & Yb(:) & Ym(:)>0.95))); 
  [WIth,WMv] = hist(Ym(Ygs(:)<0.1 &  Yb(:) & Ym(:)>mean([GMth,Wth]) & Ym(:)<Wth*1.1),0:0.01:2);
  WIth = find(cumsum(WIth)/sum(WIth)>0.8,1,'first'); WIth = round2(WMv(WIth),rf);  
  %[BIth,BMv] = hist(Ym(Ym(:)<mean([BGth,CMth]) & Yg(:)<0.2),-1:0.01:2);
  %BIth = find(cumsum(BIth)/sum(BIth)>0.02,1,'first'); BIth = round2(BMv(BIth),rf);  
  Ym   = Ym ./ WIth; 
  
  cat_io_cmd(' ','','',verb,stime); 
end
function X = round2(X,N)
  if ~exist('N','var'), N = 0; end
  X = round(X*10^N)/10^N; 
end
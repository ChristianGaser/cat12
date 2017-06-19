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
  
  resMth  = 0.5; % general resolution limitation 
  resTth  = 1.0; % tissue intensity estimation resolution limitation 
  
  if ~exist('verb','var'),     verb = 1; end
  if ~exist('fwhm','var'),     fwhm = 0.5; end
    
  stime = cat_io_cmd('    Initialize','g5','',verb);
  Vo   = spm_vol(Po);
  Yo   = single(spm_read_vols(Vo));
  Ym   = vout.Ym; 
  Ymo  = Ym; 
  Ycls = vout.Ycls; 
  res  = vout.res; 
  if ~debug; clear vout; end 

  % general resolution limitation 
  [Ym,resTM] = cat_vol_resize(Ym,'reduceV',vx_vol,resMth,32,'meanm');
  Yo         = cat_vol_resize(Yo,'reduceV',vx_vol,resMth,32,'meanm');
  for i=1:6, Ycls{i} = cat_vol_resize(Ycls{i},'reduceV',vx_vol,resMth,32); end
  vx_vol = resTM.vx_volr;

 
  %% prepare data
  Yp0   = single(Ycls{1})/255*2 + single(Ycls{2})/255*3 + single(Ycls{3})/255;
  Yb    = cat_vol_morph(cat_vol_morph(Yp0>0.6,'lo'),'lc',1); 

  % global intensity normalization 
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

  
  % gradient & divergence maps
  stime = cat_io_cmd('    Prepare measures','g5','',verb,stime);
  Yb2   = cat_vol_morph(cat_vol_morph(Ycls{2}>8 & Ycls{5}<192,'l'),'lc'); 
  Yg    = cat_vol_grad(Ymi,vx_vol*2); gth = mean(Yg(Yb(:)));
  Ydiv  = cat_vol_div(Ymi,vx_vol/2);
  Ycsfd = cat_vbdist(single(Yp0<1.5 | Ycls{3}>128),true(size(Yp0)),vx_vol); % csf distance

  %% brain and head distances
  [Yp0r,Ymir,Ybr,resT1] = cat_vol_resize({Yp0,Ymi.*Yb,Yb},'reduceV',vx_vol,2,32,'meanm');
  Ybd   = cat_vbdist(min(Yp0r,1) ,true(size(Yp0r)),resT1.vx_volr); 
  Yhd   = cat_vbdist(single(~Ybr),true(size(Yp0r)),resT1.vx_volr); Yhd  = Yhd/max(Yhd(Yhd(:)<10^6));
  Yd    = cat_vbdist(single(Ymir<0.55)); Yd  = Yd/max(Yd(:));
  Ydi   = cat_vbdist(single(Yd>0.8),Ymir>0.5); Ydi = Ydi/max(Ydi(Ydi(:)<10^4));
  Ybd   = cat_vol_resize(Ybd,'dereduceV',resT1); % brain distance
  Yhd   = cat_vol_resize(Yhd,'dereduceV',resT1); % head distance
  Ydi   = cat_vol_resize(Ydi,'dereduceV',resT1); % head distance
  %Ybs   = cat_vol_smooth3X(Yb,16/mean(vx_vol)); Ybs = Ybs/max(Ybs(:));

  %% == tissues ==
  
  % subcortical structures
  Yss  = Ymi>0.6 & (Ymi + max(0,Ydi*30-3) < 0.98 ) & Ymi<0.98 & Yhd./(Yhd + Ydi)>0.7 & Ydi<0.2; 
  Ybgi = cat_vol_localstat(Yss.*Ymi,Yss,2,1);
  Yss(abs(((Yss.*Ymi)./Ybgi)-Yss)>0.1 | abs(((Yss.*Ymi)./Ybgi)-Yss)==0)=0;
  if ~debug, clear Ybgi; end 
  Yss(smooth3(Yss)<0.6)=0; 
  Yss  = cat_vol_morph(Yss,'l',[6 0.2])>0;
  [Yssc,resT1] = cat_vol_resize(Yss,'reduceV',vx_vol,4,32,'meanm');
  Yssc  = cat_vol_morph(Yssc,'c',8); 
  Yssc  = cat_vol_resize(Yssc,'dereduceV',resT1);
  Yss   = Yss | (Yssc & Ymi<2.8/3 & Ymi>0.5 & cat_vol_morph(Yss,'d',4));
  Yss(smooth3(Yss)<0.6)=0; 
  if ~debug, clear Yssc; end 
  
  % cortex close to head
  Yct = Ycsfd>0 & Ycsfd<3 & Yhd<0.2 & cat_vol_smooth3X(Ycls{6}>128,16)>0.01 & Ymi>0.5 & Yg<2*gth & Ycls{2}<128;
  Yct(smooth3(Yct)<0.5) = 0;

  % WM
  Ywm = Yb & (Yb2 | cat_vol_morph(Ymi>.95 & Ymi<1.2 & Yp0>2 & (Yp0<1.8 | Yp0>2.2) ,'lc')) &  ...
        Ycls{1}<240 & Ycls{3}<32 & Yg<0.5 & abs(Ydiv)<0.5 & Ycsfd>2 & ~Yss & ~Yct;
  Ywm = Ywm | (Yb & Ymi-Ydiv+((Ycsfd-3)/10)>0.9 & Ydiv<-0.05 & Ycls{1}<240 & Ycls{3}<32 & ~Yss & ~Yct); 
  Ywm = cat_vol_morph(Ywm,'l',[3 0.1])>0; Ywm(smooth3(Ywm)<0.5)=0;
  
  % GM
  Ygm  = Yb & ~Ywm & Ycls{1}>64 & Ycls{2}<128 & Ycls{3}<128 & Yg./Ymi<gth & abs(Ydiv)<gth*3;  
  Ygm  = Ygm | Yct | ...
         (Yss & abs(Ydiv)<gth & Yg./Ymi<gth*1 & Ymi<0.95 & Ycls{1}>4 & Ycls{2}<252 ); %  | ...
         %(Ycls{1}>64 & Ycls{2}<240 & Ycls{3}<128 & Ycsfd<3 & Yhd<3 & ~Ywm & abs(Ydiv)<gth*3); 
  Ygmi = cat_vol_localstat(Ygm.*Ymi,Ygm,2,1);
  Ygm(abs(((Ygm.*Ymi)./Ygmi)-Ygm)>gth*mean(vx_vol)*3 | abs(((Ygm.*Ymi)./Ygmi)-Ygm)==0)=0;
  Ygm(smooth3(Ygm)<0.4)=0; 
  
  % CM (this did not work 
  if 0
    Ycm = Yb & ~Ywm & ~Ygm & Yg<gth*2 & Ymi>0.1 & (Ymi-Ydiv)<0.5 & Ycls{3}>128; 
    Ycm(smooth3(Ycm)<0.6)=0;
  end
  
  % HM 
  % there is now way to use high intensity information from the scull, but
  % it is possbile to use the large areas of mussels 
  Yh2 = Ybd>2 & smooth3( Yg>gth*2 | Ymi>1.2 | Ydiv<-0.1)>0.5; 
  Yhm = Yg<gth & ~Yb & abs(Ydiv)<0.1 & Ycls{6}<128 & Ym<T3th(3)*1.2 & Ybd>3 & ...
    Ydiv>-0.1 & Ydiv<0.1 & ~Yh2 & cat_vol_smooth3X(Ycls{6}>128,4)<0.5;
  Yhmi = cat_vol_localstat(Yhm.*Ymi,Yhm,1,2);
  Yhm(abs(((Yhm.*Ymi)./Yhmi)-Yhm)>0.5 | abs(((Yhm.*Ymi)./Yhmi)-Yhm)==0)=0;
  Yhm(smooth3(Yhm)<0.4)=0; 
  
  %%
  stime = cat_io_cmd('    Smooth values','g5','',verb,stime);
  Ywmi = cat_vol_median3(Ywm.*Yo,Ywm,Ywm,0.2);
  if inv_weighting % PVE filting 
    Ywmi = cat_vol_localstat(Ywmi,Ywm,1,2);
  else
    Ywmi = cat_vol_localstat(Ywmi,Ywm,1,3);
  end
  Ygmi = cat_vol_median3(Ygm.*Yo,Ygm,Ygm); Ygmi = cat_vol_localstat(Ygmi,Ygm,1,1); 
  Yhmi = cat_vol_median3(Yhm.*Yo,Yhm,Yhm); Yhmi = cat_vol_localstat(Yhmi,Yhm,1,1);
  if exist('Ycm','var')
    Ycmi = cat_vol_median3(Ycm.*Yo,Ycm,Ycm); Ycmi = cat_vol_localstat(Ycmi,Ycm,1,1);
  end

  hdth = res.mn(res.lkp==5); hdth(hdth<T3th(2) | hdth>T3th(3)) = []; % threshold for head tissues 
  if isempty(hdth), hdth=T3th(2); end
  Ytmi = Ywmi + ...
         Ygmi*(T3th(3)/T3th(2)) + ...
         Yhmi*(T3th(3)/hdth); % + ... cat_stat_nanmedian(Yhmi(Yhmi(:)>0 & Yg(:)<0.1 & Ybd(:)<10)));
  if exist('Ycm','var')
    Ytmi = Ytmi + Ycmi*(cat_stat_nanmean(Ym(Ywm(:)))/cat_stat_nanmean(Ym(Ycm(:) & Yg(:)<gth/2 & Ycls{3}(:)>192)));
  end
  
  % final smoothing
  for i=1:2, Ytmi = cat_vol_localstat(Ytmi,Ytmi>0,2,1); end
  
  % bias field estimation 
  stime = cat_io_cmd('    Estimate bias field','g5','',verb,stime);
  Ywi = cat_vol_approx(Ytmi,'',vx_vol,2,struct('lfO',8*fwhm));

  if debug
    Ymc = Yo ./ Ywi; Ymc = cat_main_gintnorm(Ymc * (mean(Ym(Ywm(:)))/mean(Ymc(Ywm(:)))) ,Txth); 
  end

  %%
  Ywi = cat_vol_resize(Ywi,'dereduceV',resTM); 
  Ywm = cat_vol_resize(Ywm,'dereduceV',resTM)>0.9; 
  Yo  = single(spm_read_vols(Vo));
  Ymc = Yo ./ Ywi; Ymc = Ymc * (mean(Ymo(Ywm(:)))/mean(Ymc(Ywm(:))));
  cat_io_cmd(' ','','',verb,stime); 

end
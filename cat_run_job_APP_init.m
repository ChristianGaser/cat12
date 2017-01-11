function [Ym,Yt,Ybg,WMth,bias,Tth] = cat_run_job_APP_init(Ysrco,vx_vol,verb,icall)
%  _____________________________________________________________________
%  The rough bias correction is a subfunction of cat_run_rob.
% 
%  All tissues (low gradient areas) should have a similar intensity.
%  A strong smoothing of this approximation is essential to 
%  avoid anatomical filtering between WM and GM that can first 
%  be seen in overfitting of the subcortical structures.
%  However, this filtering will overcorrect head tissue with
%  a typical intensity around GM.
%  _____________________________________________________________________
%  Robert Dahnke
%  $Id$

  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,'cat_run_job_APP_init'); debug = 1; break; end; end
  if ~exist('icall','var'), icall = 1; end

%    ds('l2','',0.5,Yo/WMth,Yg<0.2,Yo/WMth,Ym,80)
  Ysrcmin = min(Ysrco(~isinf(Ysrco(:))));
  Ysrco   = Ysrco - Ysrcmin;

  rf = 10^9; 
  bfsmoothness = 3; 
  if verb, fprintf('\n'); end
  
  stime = cat_io_cmd('  Initialize','g5','',verb);
  msize = 222; %round(222 ./ max(size(Ysrco).*vx_vol) .* min(size(Ysrco).*vx_vol));  

  %%
  Ysrc = Ysrco + 0; spm_smooth(Ysrc,Ysrc,0.8./(vx_vol./mean(vx_vol)));
  [Ysrc,resT3] = cat_vol_resize(Ysrc,'reduceV',vx_vol,min(1.2,cat_stat_nanmean(vx_vol)*2),msize,'meanm'); 

  % initial thresholds
  % correction for negative backgrounds (MT weighting)
  Yg    = cat_vol_grad(Ysrc,resT3.vx_volr) ./ Ysrc; 
  Ymsk  = Ysrc>0 & Yg>0 & Yg<mean(Yg(Yg(:)>0)); 
  [x,y] = hist(Ysrc( Ymsk(:)  ),200); cx = cumsum(x)/sum(x);
  WMth0 = y(min([numel(cx),find(cx>0.99,1,'first')])); 
  BGth0 = y(max([1        ,find(cx<0.01,1,'last' )]));
  Tth0  = WMth0*0.2 + 0.8*BGth0; 
  [x,y] = hist(Ysrc( Ymsk(:) & Ysrc(:)<Tth0 ),200); cx = cumsum(x)/sum(x); BGth1 = y(max([1        ,find(cx<0.01,1,'last')])); 
  [x,y] = hist(Ysrc( Ymsk(:) & Ysrc(:)>Tth0 ),200); cx = cumsum(x)/sum(x); WMth1 = y(min([numel(cx),find(cx>0.90,1,'first')])); 
  Ym    = (Ysrc - BGth1) ./ (WMth1 - BGth1); 

  %% improved WM threshold
  Yg    = cat_vol_grad(Ym,resT3.vx_volr) ./ max(0.3,Ym); 
  Ydiv  = cat_vol_div(Ym,resT3.vx_volr/2) ./ (Ym+eps); % lower resolution is 8 times faster 
  Ymsk  = Yg>0 & Yg<mean(Yg(Yg(:)>0))/4 & Yg./abs(Ydiv)<0.5 & ...
           Ym>cat_stat_nanmean(Ym(Yg(:)<0.2 & Ym(:)>cat_stat_nanmean(Ym(:)))) & Ym<3; 
  Ymsk  = cat_vol_morph( Ymsk ,'lo',1); 
  WMth2 = roundx(single(cat_stat_nanmean( Ysrc( Ymsk(:) ) )),rf); if ~debug, clear WMth1 Ymsk, end
  Ym    = (Ysrc - BGth1) ./ (WMth2 - BGth1);
  
  
  %% background
  stime = cat_io_cmd('  Estimate background','g5','',verb,stime);
  Ybg   = Yg<(mean(Yg(Ym(:)~=0))) | isnan(Yg);  
  % avoid edges on the image border
  bb = 1; Ybb1 = true(size(Ybg)); Ybb1(bb+1:end-bb,bb+1:end-bb,bb+1:end-bb) = 0; 
  bb = 4; Ybb2 = true(size(Ybg)); Ybb2(bb+1:end-bb,bb+1:end-bb,bb+1:end-bb) = 0; 
  Ybgc  = Ybb2 & cat_vol_morph(Ybg | (Ybb1 & Yg<0.5),'c',2);
  Ybg = Ybg | Ybgc | smooth3(Yg./Ydiv > 100)>0.5;
  if ~debug, clear Ybb1 Ybb2 bb; end
  % filling
  [Ybg,resT2] = cat_vol_resize(single(~Ybg),'reduceV',resT3.vx_volr,2,32,'meanm'); 
  Ybg = Ybg>0.5;
  Ybg  = cat_vol_morph(Ybg,'lc',4);
  Ybg  = cat_vol_smooth3X(Ybg,2); 
  Ybg  = cat_vol_resize(Ybg,'dereduceV',resT2)<0.5;    
  zeroBG = cat_stat_nanmean(Ym(Ybg(:)>0))<0.2;
  if zeroBG, BGth1 = mean(Ysrc(Ybg(:))); end
  
  WMth3 = WMth2 * roundx(single(cat_stat_nanmedian(Ym(Yg(:)>0 & Yg(:)<mean(Yg(:))/4 & Yg(:)./abs(Ydiv(:))<0.5 & ~Ybg(:) & ...
           Ym(:)>cat_stat_nanmean(Ym(Yg(:)<0.2 & ~Ybg(:) & Ym(:)>cat_stat_nanmean(Ym(:))))))),rf); if ~debug, clear WMth2, end
  Ym    = (Ysrc - BGth1) ./ (WMth3 - BGth1);
  
   
  %% first WM inhomogeneity with low tissue boundary (may include CSF > strong filtering for IXI175)
  stime = cat_io_cmd('  Initial correction','g5','',verb,stime);
  Yms  = cat_vol_smooth3X( min(2 .* ~Ybg,Ym .* (Ydiv>-0.2) .* ~Ybg .* (Ym>0.1)),16*mean(vx_vol));     % this map is to avoid CSF in the mask!
  Yms  = (Yms ./ mean(Yms(~Ybg(:))));
  Yms  = cat_vol_smooth3X( min(Yms*1.5 .* ~Ybg,Ym .* ~Ybg),16*mean(vx_vol));
  Yms  = (Yms ./ mean(Yms(~Ybg(:))));
  Yt   = Ym>max(0,Yms*0.3) & Ym<Yms*2 & Ym<(1+Yms*2) & Yg<0.5 & Ydiv<0.2 & ~Ybg & ...
         Ydiv>-0.6 & smooth3(Ym./(Yms+eps).*Yg.*Ydiv<-0.2)<0.3 & ~Ybg; Yt(smooth3(Yt)<0.5)=0;
  Ywi  = (Ym .* Yt) ./ max(eps,Yt);  
  if ~zeroBG
    Ybg2 = Ybg(:) & Yg(:)<(cat_stat_nanmean(Yg(Ybg(:))) + 2*(cat_stat_nanstd(Yg(Ybg(:))))); 
    Ywi(Ybg2) = Ym(Ybg2) / cat_stat_nanmean(Ym(Ybg2(:))); clear Ybg2;
  else
    %Ybg2 = (Yms<0.1 & Ybg) | smooth3(Yg./Ydiv > 1000)>0.5; 
  end
  [Ywi,resT2] = cat_vol_resize(Ywi,'reduceV',resT3.vx_volr,cat_stat_nanmean(resT3.vx_volr)*2,32,'max'); 
  for i=1:1, Ywi = cat_vol_localstat(Ywi,Ywi>0,2,3); end % only one iteration!
  for i=1:4, Ywi = cat_vol_localstat(Ywi,Ywi>0,2,1); end
  Ywi  = cat_vol_approx(Ywi,'nn',resT2.vx_volr,4);
  Ywi  = cat_vol_smooth3X(Ywi,bfsmoothness.*mean(vx_vol)); % highres data have may stronger inhomogeneities 
  Ywi  = cat_vol_resize(Ywi,'dereduceV',resT2);    
  Ym   = Ym./Ywi;
  WMt2 = roundx(cat_stat_nanmedian(Ym(Yg(:)<0.2 & Ym(:)>0.9)),rf); 
  Ywi  = Ywi * WMt2;
  
  %% background update
  if zeroBG
    stime = cat_io_cmd('  Refine background','g5','',verb,stime);
    Ybg  = ((Yg.*Ym)<cat_vol_smooth3X(Ym,2)*1.2) & Ym>0.2;
    Ybg  = Ybg & ~Ybgc;
    Ybg  = Ybg & ~isnan(Yg) & ~isnan(Ym); 
    [Ybg,resT2] = cat_vol_resize(single(Ybg),'reduceV',resT3.vx_volr,2,32,'meanm'); 
    Ybg  = Ybg>0.5;
    Ybg  = cat_vol_morph(Ybg,'lc',8);
    Ybg  = cat_vol_smooth3X(Ybg,2); 
    Ybg  = cat_vol_resize(Ybg,'dereduceV',resT2)<0.5 & (Ym<0.2 | isnan(Ym));
    Ybg  = cat_vol_morph(Ybg,'lo');
  end
  if ~debug, clear Ybgc Ybb1 Ybb3 bb; end
  
  %% second WM inhomogeneity with improved Yt with higher lower threshold (avoid CSF and less filtering)
  stime = cat_io_cmd('  Final correction','g5','',verb,stime);
  Yt   = Ym>max(0,Yms*0.3) & Ym>0.2 & Ym<1.2 & Ym<Yms*2 & Yg<0.2 & Ydiv<0.2 & Ydiv>-0.6 & ...
         smooth3(Ym./(Yms+eps).*Yg.*Ydiv<-0.1)<0.1 & ~Ybg; Yt(smooth3(Yt)<0.5)=0;
  Yt   = Yt | (~Ybg & Ym>0.1 & Ydiv./(Yg+eps)<0.5 & (Ym>0.3 & Yg>0.1 & Ydiv<0) | (~Ybg & Ym>0.6) & Ym<1.2 & Yg<0.1); 
  Yt   = Yt & Ym>Yms*0.3 & Ym<Yms*1.2 & ~(-Ydiv.*Ym./Yms>0.15);
  Yt(smooth3(Yt)<0.5)=0;
  Ywi2 = ( Ym .* Yt) ./ max(eps,Yt);
  % it would be nice to use futher regions, but as far as we did not know
  % their average intensity in relation to the bias field it not so easy
  if ~zeroBG
    Ybg2 = Ybg(:) & Yg(:)<(cat_stat_nanmean(Yg(Ybg(:))) + 2*(cat_stat_nanstd(Yg(Ybg(:))))); 
    Ywi2(Ybg2)   = Ym(Ybg2) / cat_stat_nanmean(Ym(Ybg2(:))); %clear Ybg2;
  else
    %Ybg2 = (Yms<0.1 & Ybg) | smooth3(Yg./Ydiv > 1000)>0.5; 
    %Yhht = -Ydiv.*Ym./Yms>0.2; 
    %Ywi2(Yhht)   = Ym(Yhht) / cat_stat_nanmean(Ym(Yhht(:))); %clear Yhht;
  end
  [Ywi2,resT2] = cat_vol_resize(Ywi2,'reduceV',resT3.vx_volr,cat_stat_nanmean(resT3.vx_volr)*2,32,'max'); 
  for i=1:1, Ywi2 = cat_vol_localstat(Ywi2,Ywi2>0,2,3); end % only one iteration!
  for i=1:4, Ywi2 = cat_vol_localstat(Ywi2,Ywi2>0,2,1); end
  Ywi2 = cat_vol_approx(Ywi2,'nn',resT2.vx_volr,4);
  Ywi2 = cat_vol_smooth3X(Ywi2,bfsmoothness.*mean(vx_vol)); %.*mean(vx_vol)); % highres data have may stronger inhomogeneities 
  Ywi2 = cat_vol_resize(Ywi2,'dereduceV',resT2);    
  Ywi  = Ywi2 .* Ywi; % both bias fields
  bias = std(Ywi(:))/mean(Ywi(:)); 
  
  
  %% BG inhomogeneity (important for normalization of the background noise)
  %[Ybc,Ygr,resT2] = cat_vol_resize({Ysrc./Ywi,Yg},'reduceV',resT3.vx_volr,cat_stat_nanmean(resT3.vx_volr)*4,16,'meanm'); 
  %Ybc  = cat_vol_morph(Ybc<BGth/WMth*2 & Ygr<0.05,'lc',2);
  %Ybc  = cat_vol_resize(smooth3(Ybc),'dereduceV',resT2)>0.5; 
  %{
  if zeroBG
    stime = cat_io_cmd('  Background correction','g5','',verb,stime);
    [Ybc,resT2] = cat_vol_resize(single(Ysrc .* Ybg),'reduceV',resT3.vx_volr,max(8,min(16,cat_stat_nanmean(resT3.vx_volr)*4)),8,'meanm'); 
    Ybc  = cat_vol_localstat(Ybc,Ybc>0,2,1);
    Ybc  = cat_vol_approx(Ybc,'nn',resT2.vx_volr,4); % no aproximation to correct only in the background! 
    Ybc  = cat_vol_smooth3X(Ybc,4);
    Ybc  = cat_vol_resize(Ybc,'dereduceV',resT2) .* (~Ywi); 
  else
    Ybc  = ones(size(Ywi),'single')*BGth1;
  end
  %}
  
  %% back to original size
  %clear Ysrc; 
  stime = cat_io_cmd('  Final scaling','g5','',verb,stime);
  Ywi       = cat_vol_resize(Ywi,'dereduceV',resT3);   
  Yg        = cat_vol_resize(Yg,'dereduceV',resT3); 
  Ydiv      = cat_vol_resize(Ydiv,'dereduceV',resT3); 
  Yms       = cat_vol_resize(Yms,'dereduceV',resT3); 
  [Yt,Ybg]  = cat_vol_resize({single(Yt),single(Ybg)},'dereduceV',resT3); Yt = Yt>0.5; Ybg = Ybg>0.5;

  
  %% intensity normalization (Ybc is the average background noise)
  % in data with strong inhomogeneities (7T) the signal can trop below the noise level 
  Ym   = ( (Ysrco - BGth1) ./ (WMth3 - BGth1)) ./ Ywi;
  Ymw  = cat_vol_morph(Yg<mean(Yg(:)) & ~Ybg & Yg./abs(Ydiv)<0.5 & Yms>max(Yms(:))*0.5,'lo',1); 
  Wth  = single(cat_stat_nanmedian(Ym(Ymw(:)))); 
  Ymw  = cat_vol_morph(Yg<mean(Yg(:)) & Ym>Wth*0.5 & Ym<Wth*1.5 & ~Ybg & Yg./abs(Ydiv)<0.5 & Yms>max(Yms(:))*0.5,'lo',1); 
  [WIth3,WMv] = hist(Ym(Ymw(:)),1000);
  WItm = find(cumsum(WIth3)/sum(WIth3)>0.8,1,'first'); WItm = roundx(WMv(WItm),rf); 
  Ym   = Ym ./ WItm; 
  Yo   = ((Ysrco-BGth1)/(WMth3-BGth1))./WItm; 
  [WIth3,WMv] = hist(Yo(Ymw(:)),1000);
  WMto = find(cumsum(WIth3)/sum(WIth3)>0.7,1,'first'); WMto = roundx(WMv(WMto),rf); 
  Yo   = Yo / WMto;
 
  
  %% update WMth of the original image!
  Ysrco = Ysrco + Ysrcmin;
  [WIth3,WMv] = hist(Ysrco(Yg(:)<0.2 & Ym(:)>Wth*0.5 & Ym(:)<Wth*1.5 & ~Ybg(:)),1000);
  WMth = find(cumsum(WIth3)/sum(WIth3)>0.7,1,'first'); WMth = roundx(WMv(WMth),rf); 
  
  
  %% check if correction was successful
  Ymx = cat_vol_morph(Ymw,'lo'); 
  Ymm = Ym .* Ymx; Ymm = cat_vol_resize(Ymm,'reduceV',vx_vol,min(1.2,cat_stat_nanmean(vx_vol)*4),msize,'meanm'); mm = Ymm(Ymm(:)>0);
  Yom = Yo .* Ymx; Yom = cat_vol_resize(Yom,'reduceV',vx_vol,min(1.2,cat_stat_nanmean(vx_vol)*4),msize,'meanm'); mo = Yom(Yom(:)>0);
  Ymb = Ym .* Ymx; Ymb = cat_vol_resize(Ymb,'reduceV',vx_vol,min(1.2,cat_stat_nanmean(vx_vol)*4),msize,'meanm'); bm = Ymb(Ymb(:)>0);
  Yob = Yo .* Ymx; Yob = cat_vol_resize(Yob,'reduceV',vx_vol,min(1.2,cat_stat_nanmean(vx_vol)*4),msize,'meanm'); bo = Yob(Yob(:)>0);
  Tth.biascorr = mean([ (cat_stat_nanstd(mo)/cat_stat_nanmean(mo)) / (cat_stat_nanstd(mm)/cat_stat_nanmean(mm)) ...
                        (cat_stat_nanstd(bo)/cat_stat_nanmean(bo)) / (cat_stat_nanstd(bm)/cat_stat_nanmean(bm)) ]); 
  if ~debug, clear Ymw Ymm Yom Ymb Yob; end 

  Tth.inverse  = Wth<0.6; 
  Tth.Tmax     = WMth/WMth0;
  cat_io_cmd('','','',verb,stime); 
  
  %% value greater 1 describe a successful correction
  if Tth.biascorr<1.01
    
    % try inverse processing 
    if icall || Tth.inverse %&& Tth.biascorr>1.05 
      if verb>1
        cat_io_cprintf('warn',sprintf('    T1-Bias correction failed (CJVR=%0.2f), try inverse correction: ',Tth.biascorr)); 
      end
      
      % second call with inverted image (due to maximum correction theorem)  
      Ysrcmax = max(Ysrco(~isinf(Ysrco(:)))); 
      [Ym2,Yt2,Ybg2,WMth3,bias2,Tthi] = cat_run_job_APP_init(Ysrcmax - Ysrco,vx_vol,verb,0); % no further recursion
    
      Ym = (Ym2 - min(Ym2(:))) ./ (1 - min(Ym2(:))); Ym(Ybg2) = 0; % inverted map (T2/PD > T1
      Ym = Ym ./ mean(Ym(Ymx(:))); % well inverse map need another scaling otherwise we cut the intensity in affreg
      %Ym = (max(Ym2(Ybg2(:))) - Ym2) ./ ([max(Ym2(Ybg2(:))) - 1); % original map (same weighting as input)
      Yt = Yt2; Ybg = Ybg2; WMth = Ysrcmax - WMth3; bias = bias2; Tth = Tthi; 
      
      Tth.inverse = 1; 
    else
      if verb>1
        cat_io_cprintf('warn',sprintf('    Bias correction failed use only scaling (CJVR=%0.2f). \n',Tth.biascorr)); 
      end
    
      Ym = Yo;
    end
  else    
    if verb>1
      cat_io_cprintf([0 0.5 0],sprintf('    Bias correction successful (CJVR=%0.2f). \n',Tth.biascorr)); 
    end
  end
  if verb>1 && icall, cat_io_cmd(' ','','',verb); end
      
end
%=======================================================================
function r = roundx(r,rf)
  r(:) = round(r(:) * rf) / rf;
end
%=======================================================================

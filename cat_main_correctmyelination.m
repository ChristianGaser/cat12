function [Ym2,Ysrc2,Ycls,Ycor,cf] = cat_main_correctmyelination(Ym,Ysrc,Ycls,Yb,vx_vol,V,T3th,LASmyostr,Yy,cat12atlas,tpm,fname)
%cat_main_correctmyelination. Correction of GM myelination/artefacts.
% _________________________________________________________________________
% This function perform high frequency correction in the GM ribbon. We
% expect a cortical ribbon of about 2-3 mm characterized by increased 
% variance. A fast thickness estimation is used to dectect underestimated 
% areas and have en estimate for the CSF distance. 
%
% [Ym2,Ysrc2,Ycls,Ycor,glcor,cf] = cat_main_correctmyelination(...
%   Ym,Ysrc,Ycls,Yb,vx_vol,vx_volo,T3th,LASmyostr,Yy,cat12atlas)
%
%  Ym[2]      .. intensity normalized image (0 BG, 1/3 CSF, 2/3 GM, 1 WM)
%  Ysrc[2]    .. bias corrected image with original intensity scaling
%  Ycls       .. tissue classes (uint8)
%  Ycor       .. correction map
%  Yb         .. brain mask
%  vx_vol     .. voxel volume (internal interpolated)
%  vx_volo    .. voxel volume (original)
%  T3th       .. tissue treshholds in Ysrc
%  LASmyostr  .. correction value
%  Yy         .. deformation map for atlas
%  cat12atlas .. atlas map for cortex setting
% 
% Called only from cat_main[#]. 
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

%
% TODO: 
%  * test with PD/T2 weightings
%  * tests with different distance estimation resolution 
%  * e
%  . thickness setting should depend on brain size to support animal processing
%    (currently only partially relevant as larger marmals also have 1-1.5 mm) 

  dbs = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

  if LASmyostr == 2
    GMV = single(sum(Ycls{1}(:)))/255 * prod(vx_vol) / 1000; 
    WMV = single(sum(Ycls{1}(:)))/255 * prod(vx_vol) / 1000; 
    LASmyostr = min(1,max(0,WMV / GMV - 1)); 
  end
  
  Yp0  = single(Ycls{3})/255/3 + single(Ycls{1})/255*2/3 + single(Ycls{2})/255;

  % To save time & memory operations are performed in a timmed subspace 
  % that is given by the brainmask.
  if ~debug
    [Ybb,BB] = cat_vol_resize( Yb      , 'reduceBrain' , vx_vol , 3 , Yb );
    Ymb      = cat_vol_resize( Ym      , 'reduceBrain' , vx_vol , BB.BB );
    Yp0b     = cat_vol_resize( Yp0     , 'reduceBrain' , vx_vol , BB.BB );
    Ycls1b   = cat_vol_resize( Ycls{1} , 'reduceBrain' , vx_vol , BB.BB );
    Ycls2b   = cat_vol_resize( Ycls{2} , 'reduceBrain' , vx_vol , BB.BB );
  else
    Ybb = Yb; Ymb = Ym; Yp0b = Yp0; 
    Ycls1b = Ycls{1};  Ycls2b = Ycls{2}; 
  end
  % dilated brainmask to correct possible artefacts
  Ybb   = cat_vol_morph(Ybb,'d',2); 
 
  
  
  %% only work in the cortical areas
  %  ----------------------------------------------------------------------
  %  structures without GM ribbon
  % get some cortical ROIs to perform the corrections
  TIV   = sum(Yb(:) * prod(vx_vol)) / 1000; 
  drad  = 0.3 * TIV^(1/3) / mean(vx_vol) ; 
  if exist('Yy','var') 
  % if we use an atlas ... 
    if exist('cat12atlas','var')
      NS  = @(Ys,s) Ys==s | Ys==s+1; % function to ignore brain hemisphere coding
      LAB = cat_get_defaults('extopts.LAB'); 

      % map atlas to RAW space
      for i=1:5
        try
          VA = spm_vol(cat12atlas{1});
          break
        catch 
          % avoid read error in parallel processing
          pause(rand(1))
        end
      end
      YA = cat_vol_ctype( cat_vol_sample(tpm(1),VA,Yy,0) ); 
      %YA   = cat_vol_ctype(round(spm_sample_vol(VA,...
      %  double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0)));
    else
      YA = Yy; 
    end
    YA   = reshape(YA,size(Yp0));
    YA   = cat_vol_ctype(cat_vol_median3c(single(YA))); 

    if ~debug
      YA = cat_vol_resize( YA , 'reduceBrain' , vx_vol , BB.BB );
    end
    
    % extend atlas for all voxels within the brainmask
    [D,I] = cat_vbdist(single(YA>0),Ybb); YA = YA(I); clear D I;  %#ok<ASGLU>
    Yct   = NS(YA,LAB.CT) & ~cat_vol_morph( NS(YA,LAB.VT) | NS(YA,LAB.LE) | NS(YA,LAB.HC) | NS(YA,LAB.PH) | NS(YA,LAB.BG) | NS(YA,LAB.TH), 'dd' ,drad); 
  else
  % without atlas
    YA   = ones(size(Ybb),'single');
    Yvt  = cat_vol_morph(cat_vol_morph(round(Yp0b*3)==1,'o',2/mean(vx_vol)),'d',2/mean(vx_vol)) & ...
           cat_vol_morph(cat_vol_morph(round(Yp0b*3)==3,'o',2/mean(vx_vol)),'d',2/mean(vx_vol));
    Yct  = Ybb & ~Yvt; 
  end
  
  

  % estimate variance maps to identify the GM and its variance
  %Ysdg = cat_vol_localstat(max(1/3,Ymb),Ybb,round(3/mean(vx_vol)),4);         % values at the border between CSF and WM have high values
  %Ysdw = cat_vol_localstat(max(1/3,Ymb),Yp0b>2.8/3,round(2/mean(vx_vol)),4);  % estimate average WM variance
  %Ysdw = cat_vol_approx(Ysdw);                                       

  
  % gobal evaluation 
  % - thx: 
  % - rvth: relative volume threshold: 0 (no correction) to 1 (heavy correction)
  % - gwp:  gray-white matter proportion (>1 atypical)
  for i=1:3, vol(i) = cat_stat_nansum( single( Ycls{i}(:) )/255 ) * prod(vx_vol) / 1000; end
  LASmyo.volr = vol ./ sum(vol);
  LASmyo.thx  = min(2.2,3 - 4*LASmyo.volr(3));
  LASmyo.rvth = [  max(0,1/3 - LASmyo.volr(1))  max(0,LASmyo.volr(2) - 1/3)*3  max(0,1/3 - LASmyo.volr(3))*3 ];
  LASmyo.gwp  = LASmyo.volr(2) ./ LASmyo.volr(1); % 
  % cortical evaluation (maybe not allways defined)
  for i=1:3, ctvol(i) = cat_stat_nansum( single( Ycls{i}(Yct(:)) )/255 ) * prod(vx_vol) / 1000; end
  LASmyo.ctvolr = ctvol ./ sum(ctvol);
  LASmyo.ctthx  = min(2.2,3 - 4*LASmyo.ctvolr(3));
  LASmyo.rctvth =  max(0,[ 1/3 - LASmyo.ctvolr(1) , LASmyo.ctvolr(2) - 1/3 , 1/3 - LASmyo.ctvolr(3) ]) * 3;
  LASmyo.ctgwp  = LASmyo.ctvolr(2) ./ LASmyo.ctvolr(1);  
  
  LASmyo.clsrmse = cat_stat_nanmean( (Yp0b(:) - Ymb(:)).^2 .* Ybb(:) .* (Yp0b(:)<1 & Yp0b(:)>1/3)).^.5; % only GM
  LASmyo.noise   = min( cat_stat_nanstd( Ymb( cat_vol_morph( Yp0b>2.5/3 & Ymb  > 2.5/3,'e')) ), ...
                        cat_stat_nanstd( Ymb( cat_vol_morph( Yp0b<1.5/3 & Yp0b > 0.5/3,'e')) ) );

  %LASmyo.ctf = LASmyo.rctvth(2) .* LASmyo.rctvth(3) .* LASmyo.ctgwp * 10;
  LASmyo.ctf = min(2,sum(LASmyo.rctvth) .* LASmyo.ctgwp * 2);


  
  %% estimate GM and WM thickness 
  %  ----------------------------------------------------------------------
  %  This is just a fast approximation with many limitations!
  %  However, it should allow to detect problematic areas and roughly limit
  %  out corrections.
  %  ----------------------------------------------------------------------

  % resampling for speed and equal voxel size
  res       = 1.2;  %     % simple reduce ... hmm, is this required / useful any longer ? - yes for high res denoising ... need testing - first simple 
  resi      = 1;    %0.5; % distance/thickness estimation resolution >> higher resolution would need further adaptations 
%{ 
  [YA,redR] = cat_vol_resize(YA  ,'reduceV',vx_vol,res,32,'nearest');  
  Yp0b      = cat_vol_resize(Yp0b,'reduceV',vx_vol,res,32,'meanm'); 
  Ymb       = cat_vol_resize(Ymb ,'reduceV',vx_vol,res,32,'meanm');  
  Ybb       = cat_vol_resize(single(Ybb) ,'reduceV',vx_vol,res,32,'meanm');   
  Yct       = cat_vol_resize(single(Yct) ,'reduceV',vx_vol,res,32,'meanm');  
  Ycls1b    = cat_vol_resize(single(Ycls1b),'reduceV',vx_vol,res,32,'meanm');  
  Ycls2b    = cat_vol_resize(single(Ycls2b),'reduceV',vx_vol,res,32,'meanm');  
%}
  [YA,resI] = cat_vol_resize(YA  ,'interp',V,resi,'nearest');   
  Yp0b      = cat_vol_resize(Yp0b,'interp',V,resi,'linear');   
  Ymb       = cat_vol_resize(Ymb ,'interp',V,resi,'cubic');   
  Ybb       = cat_vol_resize(Ybb ,'interp',V,resi,'linear') > .5;   
  Yct       = cat_vol_resize(Yct ,'interp',V,resi,'linear') > .5;   
  Ycls1b    = cat_vol_resize(Ycls1b,'interp',V,resi,'linear') > .5;   
  Ycls2b    = cat_vol_resize(Ycls2b,'interp',V,resi,'linear') > .5;   

 
  %% estimate "blurred" CSF (with sulcual overestimation) and WM distance
  
  % there are some regions we will ignore
  Ydeepgm = NS(YA,LAB.VT) | NS(YA,LAB.BG) | NS(YA,LAB.TH) | NS(YA,LAB.HC); % | NS(YA,LAB.PH)
  Ydeepgm = cat_vol_smooth3X( cat_vol_morph( Ydeepgm ,'dd',2,resI.resN) , 4);
  
  % first distance maps 
  Ycd   = cat_vol_eudist(2-max(Yp0b,Ydeepgm)*3, Yct & Yp0b>1.5/3 & Yp0b<2.5/3, 2, 1); % overestimated in blurred sulci
  Ywd   = cat_vol_eudist(Ybb .* Yp0b .* 3 - 2 , Yct & Yp0b>1.5/3 & Yp0b<2.5/3, 2, 1);
  
  % estimate (filtered) thickness
  Ygmt  = cat_vol_pbtp(Yp0b * 3,Ywd,Ycd);  
  Ygmta = cat_vol_approx(Ygmt); Ygmt( abs(Ygmt - Ygmta ) > .5 ) = 0; clear Ygmta;
  Ygmt  = cat_vol_approx(Ygmt); 
  
  % estimate "real" CSF distance by using the thickness and WM distance map as Ygmt = Ycd + Ywd
  Ycdc  = Ygmt - Ywd; Ycdc(Ycd<Ycdc) = Ycd(Ycd<Ycdc); clear Ycd Ywd
  Ypp   = min(1,Ycdc ./ max(eps,Ygmt)); Ypp(Ycls2b>0 & Yct & Yp0b>2.5/3) = 1; 
  
  
  %% re-estimate the CSF distance (this time also going deep in the WM)
  %  the limitation is useful as we only need local gyral information!
  %  and want to avoid bias by deep Wm
  YdeepWM = cat_vol_morph(Ycls2b>.95,'de',5,resi); 
  Ycdc2 = cat_vol_eudist(smooth3(Ycdc<.5 & Ypp<.75 & ~Ydeepgm), Yct & ~YdeepWM, 2, 1); clear Ycdc
  Ycdc2(Ycdc2 > 6 * resi) = 0; 
 
  %% estimate the full tissue thickness (we needed the GM thickness and WM to reconstruct the sulcus)
  Ybmt  = cat_vol_pbtp( min(3,4-min(2,Yp0b*3)), Ycdc2, Ycdc2*inf); 
  Ybmt  = cat_vol_approx(Ybmt); 
  
  % estimate WM skeleleton distance map (to keep allways a central voxel)
  Ywid  = (Ybmt - Ycdc2); Ywid(Ywid>10 | Ycdc2==0) = 0; 

  % size scaling
  Ygmt = Ygmt*resi; Ybmt=Ybmt*resi; Ycdc2=Ycdc2*resi; Ywid=Ywid*resi; 
  
  % eval
  LASmyo.mnYgmt = cat_stat_nanmean(Ygmt(Ycls1b(:)>.5 & Yct(:))); 
  LASmyo.mnYbmt = cat_stat_nanmean(Ybmt(Ycls1b(:)>.5 & Yct(:))) - LASmyo.mnYgmt; 
  
    %% define regions with clear GM, where we don't want to correct
  Yisok = cat_vol_morph( cat_vol_morph(Ymb>1.5/3 & Ymb<2.5/3 ,'do',1,resI.resN),'dd',2,resI.resN) & Ymb>1.85/3 & Ymb<2.15/3 & ...
          cat_vol_morph( Ymb<2.5/3 ,'dd',2,resI.resN) & Ymb<2.5 & Ywid>2;
  Yisok = Yisok | (Ygmt>2.5 & Ycdc2>2.5 & Ymb>1.75/3 & Ymb<2.25/3) | (Ycdc2>3.5); 
  Yisok = smooth3(Yisok); 
  
  %% avoid high intenisty (WM) edges
  % define regions, where we want to correct
  Ybd = cat_vbdist( single( ~( (Ycls1b + Ycls2b)>0.5 &  cat_vol_morph(Ybb>.5,'e','3',resi) ) ) );
  Ybv = cat_vol_morph( (Ycls1b + Ycls2b)>0.5, 'd') & Ymb>2/3 & Ycdc2<2 & ~cat_vol_morph( (Ycls1b + Ycls2b)>0.5, 'e') & Yct & Ybd<5*resi; 
  Yisct = Yct & (Ywid>1 | Ybv) & Ydeepgm<.05 & (Ymb>1.5/3 | Yp0b>1.5/3) & ...
    (( Ybmt>1 & (Ycdc2 <= Ybmt-1) & (Ycdc2 < max(.5,min(2.5,LASmyo.ctf * LASmyostr*2)))) | ( (Ybmt>1 | Ybv) & Ycdc2 < 1.5));
  Yisct = min(Yisct,cat_vol_smooth3X(Yisct,.6)); % not extend but soften!
  %Yisct = max( 0, min( max(1.5, median(Ygmt(Yp0b(:)>1.5/3 & Yp0b(:)<2.5/3))) - Ycdc2/2 * LASmyo.ctf,Yisct));
  Ynpp  = max(0,1 - Ycdc2./max(eps,min(median(Ygmt(:)),Ybmt))); Ynpp(Ycdc2==0)=0; Ynpp(Ybv) = 1;  
%  Yisct = Yisct .* smooth3(max(0,min(1,Ywid/2))); % only thick structures .5 + 
  Yisct = Yisct .* (1 - Yisok/2); % only critical structures
  Yisct = Yisct .* cat_vol_smooth3X(Yct,4) .* (1 - Ydeepgm); % only neocortex
  %Yisct = Yisct .* max(smooth3(Ywd+.5),Ybv) .* Ynpp .* (2.5/Ygmt) .* min(1,LASmyo.ctf * LASmyostr*2) .* smooth3(min(1,Ybmt-1 + Ybv*10));  % important to avoid overcorrections
  Yisct = Yisct .* max(0,min(2,(Ybmt - Ygmt))) .*  Ynpp .* smooth3(median(Ygmt(:))./Ygmt) .* min(1,LASmyo.ctf * LASmyostr) .* smooth3(min(1,Ybmt-1 + Ybv*10));  % important to avoid overcorrections
 %Yisct = Yisct .* min(2 , smooth3( max(Ywd/2,Ybv)  .* max(0.5,LASmyo.ctf) )); % .* Ynpp  .* smooth3(min(1,Ybmt-1 + Ybv*10));  % important to avoid overcorrections
  Ym2   = max( min(Ymb   ,  min(2.25    ,Ymb  )*0.25 + 0.75*(2.25     - 0.25 * LASmyostr) / 3       ) , Ymb - Yisct / 3 );

  %% resample data to original space
  Yct   = cat_vol_resize(single(Yct),'deinterp' ,resI,'linear') > .5;   
  Yisct = cat_vol_resize(Yisct ,'deinterp' ,resI,'linear');   
  %Yct   = cat_vol_resize(Yct   ,'dereduceV',redR,'linear') > .5;   
  %Yisct = cat_vol_resize(Yisct ,'dereduceV',redR,'linear');  
  clear Yisok Ywid Ybmt Ygmt Ycdc2 Ydeepgm 

  %% get back to the orinal space
  if ~debug
    Yct   = cat_vol_resize( Yct   , 'dereduceBrain' , BB );
    Yisct = cat_vol_resize( Yisct , 'dereduceBrain' , BB );
  else
    Yclso = Ycls; 
  end  

  %% update segmentation 
  Ygmo    = single(Ycls{1}); 
  Ygm     = cat_vol_ctype( single( Ycls{2}) .* Yisct); 
  Ycls{2} = Ycls{2} - Ygm; 
  Ycls{1} = Ycls{1} + Ygm; 
  Ygm     = cat_vol_ctype( single( Ycls{3}) .* Yisct); 
  Ycls{3} = Ycls{3} - Ygm; 
  Ycls{1} = Ycls{1} + Ygm; 
  

  %% light bias correction
  if LASmyostr > .125 
    %%
    Tth = [cat_stat_nanmedian( Ym(Ycls{1}(:)>128) ) cat_stat_nanmedian( Ym(Ycls{2}(:)>128) ) cat_stat_nanmedian( Ym(Ycls{3}(:)>128) )]; 
    Yw = ( ((Ycls{1}+Ycls{2}+Ycls{3})>128) .* Ym) ./ (single(Ycls{1})/255*Tth(1) + single(Ycls{2})/255*Tth(2) + single(Ycls{3})/255*Tth(3));  
    Yw = Yw + (single(Ycls{6})>.5);
    Yw = cat_vol_approx(Yw); 
    Yw = cat_vol_smooth3X(Yw,4); 
    Yw = Yw ./ cat_stat_nanmedian(Yw((Ycls{1}+Ycls{2}+Ycls{3})>128));
    Ym = Ym ./ Yw; 
    Ysrc = Ysrc ./ Yw; 

    % additional intensity normalization 
  end

  %% final corrections
  if 0 % not here ... here only BC  LASmyostr >= .5
    Ym2   = max( min(Ym   ,  min(2.25    ,Ym  )*0.25 + 0.75*(2.25     - 0.25 * LASmyostr) / 3       ) , Ym - Yisct / 3 );
    Ysrc2 = max( min(Ysrc ,  min(T3th(2),Ysrc)*0.25 + 0.75*(T3th(3) - 0.5 * LASmyostr) / T3th(3) ) , ...
        Ysrc - Yisct * diff(T3th(2:3)) * 3 * (diff(T3th(1:2)) / diff(T3th(2:3)) * 0.5) ); % under/over-correction for high/low contrast
  else
    Ym2   = Ym; 
    Ysrc2 = Ysrc; 
  end  
  Ycor  = Yisct; 
  

  % eval
  for i=1:3, volc(i) = cat_stat_nansum( single( Ycls{i}(:) )/255 ) * prod(vx_vol) / 1000; end
  LASmyo.volcr = volc ./ sum(volc);
  for i=1:3, volc(i) = cat_stat_nansum( single( Ycls{i}(Yct(:)) )/255 ) * prod(vx_vol) / 1000; end
  LASmyo.volcr = volc ./ sum(volc);
  LASmyo.change = mean( Ym( Yct(:) & Ym2(:)>1.5/3 & Ym2(:)<2.5/3 )  -  Ym2( Yct(:) & Ym2(:)>1.5/3 & Ym2(:)<2.5/3 ) ) / prod(vx_vol) * 1000;
  try 
    cat_io_write2csv(fname,'LASmyo',LASmyo); 
  end
  
  
  % estiamte modification factor
  cf = abs( sum(Ygm(:)) - sum(Ygmo(:)) ) / sum(Ygmo(:)); 
  
  % display changes
  if cat_get_defaults('extopts.verb') > 2 
    if cf>0.04 || glcor>0.8
      cat_io_cprintf('warn'   ,sprintf('\n    Myelination correction of %0.2f%%%% of the GM volume (strength=%0.2f)! ', cf*100, glcor/3));  
    elseif cf>0.02 || glcor>0.4
      cat_io_cprintf('note'   ,sprintf('\n    Myelination correction of %0.2f%%%% of the GM volume (strength=%0.2f). ', cf*100, glcor/3));  
    else
      cat_io_cprintf([0 0 0.5],sprintf('\n    Myelination correction of %0.2f%%%% of the GM volume (strength=%0.2f). ', cf*100, glcor/3));  
    end
  
    % display voluminas
    for i=1:numel(Ycls), ppe.LASvols(i) = cat_stat_nansum(single(Ycls{i}(:)))/255 .* prod(vx_vol) / 1000; end
    cat_io_cprintf('blue',sprintf('\n    LAS  volumes (CGW=TIV; in mm%s):     %7.2f +%7.2f +%7.2f = %4.0f\n',...
          native2unicode(179, 'latin1'),ppe.LASvols([3 1 2]),sum(ppe.LASvols(1:3))));   

    cat_io_cmd(' ',' ','',job.extopts.verb);
  end
end
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
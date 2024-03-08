function [Ym2,Ysrc2,Ycls,Ycor,glcor,cf] = cat_main_correctmyelination(Ym,Ysrc,Ycls,Yb,vx_vol,vx_volo,T3th,LASstr,Yy,cat12atlas,tpm)
%cat_main_correctmyelination. Correction of GM myelination/artefacts.
% _________________________________________________________________________
% This function perform high frequency correction in the GM ribbon. We
% expect a cortical ribbon of about 2-3 mm characterized by increased 
% variance. A fast thickness estimation is used to dectect underestimated 
% areas and have en estimate for the CSF distance. 
%
% [Ym2,Ysrc2,Ycls,Ycor,glcor,cf] = cat_main_correctmyelination(...
%   Ym,Ysrc,Ycls,Yb,vx_vol,vx_volo,T3th,LASstr,Yy,cat12atlas)
%
%  Ym[2]      .. intensity normalized image (0 BG, 1/3 CSF, 2/3 GM, 1 WM)
%  Ysrc[2]    .. bias corrected image with original intensity scaling
%  Ycls       .. tissue classes (uint8)
%  Ycor       .. correction map
%  Yb         .. brain mask
%  vx_vol     .. voxel volume (internal interpolated)
%  vx_volo    .. voxel volume (original)
%  T3th       .. tissue treshholds in Ysrc
%  LASstr     .. correction value
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
% $Id: 2558 2024-02-28 $

%
% TODO: 
%  * thickness setting should depend on brain size to support animal processing
%  * test with PD/T2 weightings

  dbs = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

  LASmyostr = LASstr; 
  LASartstr = 1; %LASstr; 
  
  if LASmyostr == 2
    GMV = single(sum(Ycls{1}(:)))/255 * prod(vx_vol) / 1000; 
    WMV = single(sum(Ycls{1}(:)))/255 * prod(vx_vol) / 1000; 
    LASmyostr = min(1,max(0,WMV / GMV - 1)); 
  end
  
  Yp0  = single(Ycls{3})/255/3 + single(Ycls{1})/255*2/3 + single(Ycls{2})/255;

  % To save time & memory operations are performed in a timmed subspace 
  % that is given by the brainmask.
  if ~debug
    [Ybb,BB] = cat_vol_resize( Yb  , 'reduceBrain' , vx_vol , 3 , Yb );
    Ymb      = cat_vol_resize( Ym  , 'reduceBrain' , vx_vol , BB.BB );
    Yp0b     = cat_vol_resize( Yp0 , 'reduceBrain' , vx_vol , BB.BB );
    Ycls1b   = cat_vol_resize( Ycls{1} , 'reduceBrain' , vx_vol , BB.BB );
    Ycls2b   = cat_vol_resize( Ycls{2} , 'reduceBrain' , vx_vol , BB.BB );
  else
    Ybb = Yb; Ymb = Ym; Yp0b = Yp0; 
    Ycls1b = Ycls{1}; 
    Ycls2b = Ycls{2}; 
  end
  % dilated brainmask to correct possible artefacts
  Ybb   = cat_vol_morph(Ybb,'d'); 
 
  
  % estimate variance maps to identify the GM and its variance
  Ysdg = cat_vol_localstat(max(1/3,Ymb),Ybb,round(3/mean(vx_vol)),4);         % values at the border between CSF and WM have high values
  Ysdw = cat_vol_localstat(max(1/3,Ymb),Yp0b>2.8/3,round(2/mean(vx_vol)),4);  % estimate average WM variance
  Ysdw = cat_vol_approx(Ysdw);                                       
  
  
  
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
    
    if ~debug
      YA = cat_vol_resize( YA , 'reduceBrain' , vx_vol , BB.BB );
    end
    
    % extend atlas for all voxels within the brainmask
    [D,I] = cat_vbdist(single(YA>0),Ybb); YA = YA(I); clear D I;  %#ok<ASGLU>
    Yct   = NS(YA,LAB.CT) & ~cat_vol_morph( NS(YA,LAB.VT) | NS(YA,LAB.LE) | NS(YA,LAB.HC) | NS(YA,LAB.PH) | NS(YA,LAB.BG) | NS(YA,LAB.TH), 'dd' ,drad); 
  else
  % without atlas
    Yvt  = cat_vol_morph(cat_vol_morph(round(Yp0b*3)==1,'o',2/mean(vx_vol)),'d',2/mean(vx_vol)) & ...
           cat_vol_morph(cat_vol_morph(round(Yp0b*3)==3,'o',2/mean(vx_vol)),'d',2/mean(vx_vol));
    Yct  = Ybb & ~Yvt; 
  end
  
  
  
  %% estimate GM and WM thickness 
  %  ----------------------------------------------------------------------
  %  This is just a fast approximation with many limitations!
  %  However, it should allow to detect problematic areas and roughly limit
  %  out corrections.
  %  ----------------------------------------------------------------------
  
  %  we can only correct voxel to GM if there is enough space
  Ycd    = max(0,cat_vbdist(max(0,min(1,2 - max(Ymb,NS(YA,LAB.VT) | NS(YA,LAB.BG) | NS(YA,LAB.TH)) .* 3.* Ybb)),Yct) - 0.5); Ycd(Ycd>10^10) = 0; % uncorrected CSF distance estimate 
  Ywd    = max(0,cat_vbdist(max(0,min(1,Ybb .* Ymb .* 3 - 2)),Yct) - 0.5); Ywd(Ywd>10^10) = 0;      % WM distance 
  Ygmt   = cat_vol_pbtp(Yp0b*3,Ywd,Ycd) .* mean(vx_vol);                                            % GM thickness
  
  %  We need dynamic adaptation to handle different brain sizes in human
  %  but also non human data. 
  %  Denoising of the thickness estimate based on our expectations that 
  %  thickness in normally distributed but also linked to the TIV
  %  (Hofman, 1989) that is used here as lower limit to correct in bad
  %  quality cases. 
  if 1
    kn = 1; % number of neighbor peaks 
    [mnYgmt,sdYgmt] = cat_stat_kmeans(Ygmt(Ygmt(:)>0),kn*2+1); 
    mnYgmt = max( mnYgmt( kn+1 ) ,  ( TIV / 100 ).^0.2 ); 
    sdYgmt = max( sdYgmt( kn+1 ) * (kn*2+1) ,  ( TIV / 100 ).^0.2 / 4 ); 
  else
    % old aproach
    mnYgmt = cat_stat_nanmean(Ygmt(Ygmt(:)>0)); 
    sdYgmt = cat_stat_nanstd(Ygmt(Ygmt(:)>0)); 
  end
  %  denoising of the thickness estimate based on our expectations that thickness in normally distributed
  bd     = [mnYgmt mnYgmt] + [-sdYgmt sdYgmt]; 
  for  i = 1:1, Ygmt = cat_vol_median3(Ygmt,Ygmt>0,Ygmt>bd(1) & Ygmt<bd(2)); end    % remove outlier
  for  i = 1:2, Ygmt = cat_vol_localstat(Ygmt,Ygmt>0,1,1); end                      % smoothing within the GM  
  Ygmt   = cat_vol_approx(Ygmt,'nh',2); 
  if ~debug, clear bd; end
  
  % WM thickness
  [Ywdi,I] = cat_vbdist(max(0,min(1,2.6 - smooth3(Ymb) .* 3 .* Ybb)));                              % inner WM boundary distance
  Ywmt   = cat_vol_pbtp(max(2,min(3,5 - Ymb*3 .* Ybb)),Ywdi,ones(size(Ywdi),'single'))  .* mean(vx_vol); %if ~debug, clear Ywdi; end
  for  i = 1:2, Ywmt = cat_vol_localstat(Ywmt,Ywmt>0,1,1); end
  Ywmt   = cat_vol_approx(Ywmt,'nh',2);
  Ycdmm  = min( Ycd , (Ygmt/mean(vx_vol) - Ywd)) .* mean(vx_vol);                                                % final CSF/GM boundary distance map in mm
  Ycdmm(Ywdi>0) = 0; Ycdmm = Ycdmm(I) + Ywdi .* mean(vx_vol); 
  spm_smooth(Ycdmm,Ycdmm,0.3./vx_vol); 
   if ~debug, clear Ycd Ywd; end
  
  % correction factor based on the volumina and thickness values and on the size of the amount of higher GM PVE values
  clsvol     = zeros(1,3); for i=1:3, clsvol(i) = sum(single(Ycls{i}(:))/255); end %  * prod(vx_vol)/1000
  [mn,sd,nm] = cat_stat_kmeans(Ymb(Yp0b(:)>1.8/3 & Yct(:)),4);                       %#ok<ASGLU> % 4 classes [GM,HGM,LWM,WM] were we focus on the myelinated GM (LWM)  
  glcor      = min(1,max(0,( ( nm(3)/sd(3) ) / ( nm(4)/sd(4) ) - 0.1 ) * 2)) / mean(vx_volo);    % lower rate for low res -  prod was not enough correction at 2 mm
  glcor      = min(1,glcor * mean( vx_vol ./ vx_volo ) * clsvol(2)/clsvol(1)) * 3;               % lower correction in case of interpolation 
  if ~debug, clear clsvol mn sd nm; end
  
  
  
  %% combine measures
  %  ----------------------------------------------------------------------
% ######### this part need improvments and better explaination 
% ######### mark points were that allow easy modification/adaptations

  % Ywmp .. WM region and its PVE that should not or only partially be alterated 
  Ydiv  = cat_vol_div(Ywdi); 
  Ywmp  = single( Ycdmm>1.5 & ... tissue deeper than about 1.5 mm 
    ( ( Ymb - Ysdg ) > 2.5/3  | ... save WM tissue
    ( ( Ydiv < -0.05 .* Ymb ) & Ymb > 2.75 ) | ... WM divergence skeletton ... 
    ( max(0,((Ywdi + 1) ./ (Ywmt + 0.5))*4 - 3) > 0.25 & Ymb>2.75/3 ) ) ); 
  Ywmp(smooth3(Ywmp)<0.4) = 0;
  Ywmp  = cat_vol_morph(Ywmp,'l',[inf round( 5 / mean(vx_vol)) ]); % remove tiny dots
  Ywmps = min(1,Ywmp + max( max(0,-Ydiv) , 0.8*smooth3(Ywmp)) ); spm_smooth(Ywmps,Ywmps,0.5./vx_vol); 
  Ywmp  = min(1,max(Ywmp*0.95,Ywmps*1.5)); % use smoothing as PVE
%  Ywmp  = Ywmp .* max(single(Ycls2b)/255,cat_vol_smooth3X(cat_vol_morph(Ybb,'de',mnYgmt + sdYgmt + 1,vx_vol),1)); 
  if ~debug, clear Ywmps Ydiv; end
  
  
  %% Yct  .. cortical region
  if 0
    mnYgmt = 2.5;
    sdYgmt = 0.5;
  else 
    mnYgmt = mnYgmt; % + 0.5;   
  end
  
  Yct2  = Yct .* max(0,min(1,mnYgmt - Ycdmm)); spm_smooth(Yct2,Yct2,2./vx_vol);
  if ~debug, clear Yct; end
  
  %% Ycm  .. first correction map (range 0-1') + 
  Ycm   = (smooth3(Ymb<0.98 & Ymb>2.5/3 & Ycdmm<mnYgmt) + max(0,Ymb - 2/3)) .* ... GM-WM Voxel
          Yct2 .* (1 - Ywmp) .* 3 .* ... that belong to the cortex (Yct2) and not to the WM (Ywmp)
          min(2,max(0,mnYgmt - Ycdmm)) .* ... higher correction with increased WM distance 
          ( max(0,Ywmt - Ygmt + 1)>1 | Ycdmm<(mnYgmt+sdYgmt) ) .* ... correction only in regions with enough WM
          max(0,(mnYgmt+sdYgmt) - Ygmt) .* max(0,min(1, (Ywmt-2)/2)); 
  Ycm   = Ycm * 0.9 * LASmyostr.^0.2;  
  Ycm(smooth3(Ycm>0)<0.4) = 0; % remove PVE
  Ycm(Ycm>0 & cat_vol_morph(Ycm>0,'l',[inf round( 5 / mean(vx_vol)) ])==0) = 0; % remove PVE
  Ycm(smooth3(Ycm>0)<0.4) = 0; % remove PVE
  
  % Ycma .. artefact map (range 0-1')
  Yma   = max(0,min(2,3 * max(0,Ymb - 2/3) .* max(0,1.5 - min(0,Ycdmm-1)) .* Yct2 .* ...
          max( (Ymb-max(Yp0b,Ywmp))*3 , max( 1 - Ywmp , 1-single(Ycls1b)/255) .* max(0,min(1, (Ywmt-2)/2))) )); 
  Yma   = Yma * 2 * LASartstr.^0.2;  
  
  % global limitation of the correction 
  Ycma  = min(1,max(Yma,Ycm));
  Ycma  = cat_vol_approx(Ycma);
  Ycm(Ycm > Ycma*2) = Ycma(Ycm > Ycma*2);                                  % limit corrections ##############
  Ycm   = max(Yma,Ycm .* Ycma/max(Ycma(:)));                               % normalizes setting
  Ycm   = Ycm .^ 1.5;                                                      % exponential scaling to avoid small corrections ##################
 % Ycm   = Ycm + min(0, max(0,Ymb - 2/3) - Ycm/3);                          % correct overcorrections
  if ~debug, clear Ycma Ywmt Ygmt Ybb Yct2 Ymb Yp0b Ycdmm; end
  
  
  %% get back to the orinal space
  if ~debug
    Yma   = cat_vol_resize( Yma  , 'dereduceBrain' , BB );
    Ycm   = cat_vol_resize( Ycm  , 'dereduceBrain' , BB );
    Ysdg  = cat_vol_resize( Ysdg , 'dereduceBrain' , BB );
    Ysdw  = cat_vol_resize( Ysdw , 'dereduceBrain' , BB );
  end  
  
  
  
  
  %% correct maps
  Ycor  = max(0,min(1, max( Yma , (Ysdg - (2.5-LASmyostr) * mean(Ysdw(Yp0>2.5/3)))) .* Ycm .* glcor)); 
  %Ycor  = Ycor + min(0, max(0,Ym - 2/3) - Ycor/3) * 3; % correct overcorrections
  %Ycor  = max(0,Ym - max( min(Ym , (2.9 - 0.9 * LASmyostr.^0.5) / 3) , Ym - Ycor / 3 )) * 3;
  Ycor  = Ycor * LASmyostr.^0.5;            % correct only 90% to keep some real noise; use this function for fine tuning of the LASstr
  Ycor  = max(0,Ycor); 
  if ~debug, clear Ysdg Ysdw Ycm Yma; end
  
  % update original maps
  Ym2   = max( min(Ym   ,  min(2.5    ,Ym  )*0.25 + 0.75*(2.5     - 0.5 * LASmyostr) / 3       ) , Ym - Ycor / 3 );
  Ysrc2 = max( min(Ysrc ,  min(T3th(2),Ysrc)*0.25 + 0.75*(T3th(3) - 0.5 * LASmyostr) / T3th(3) ) , ...
      Ysrc - Ycor * diff(T3th(2:3)) * 3 * (diff(T3th(1:2)) / diff(T3th(2:3)) * 0.5) ); % under/over-correction for high/low contrast

  
  
  
  %% update segmentation 
  Ygmo    = single(Ycls{1}); 
  Ygm     = cat_vol_ctype( single( Ycls{1}) + single( Ycls{2} + Ycls{3} ) .* Ycor); 
  Ycls{2} = Ycls{2} - Ycls{2} .* cat_vol_ctype(Ycor); 
  Ycls{3} = Ycls{3} - Ycls{3} .* cat_vol_ctype(Ycor); 
  Ycls{1} = Ygm; 

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
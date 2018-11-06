function [Ysrc,Ycls,Yb,Yb0,job,res,T3th,stime2] = cat_main_updateSPM(Ysrc,P,Yy,tpm,job,res,stime,stime2)
% ______________________________________________________________________
%  Update SPM preprocessing. 
%  Subfunction of cat_main.
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id: cat_main_gcut.m 1315 2018-05-03 09:34:57Z dahnke $


  dbs = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

  global cat_err_res; 

  clsint = @(x) round( sum(res.mn(res.lkp==x) .* res.mg(res.lkp==x)') * 10^5)/10^5;
  
  VT  = res.image(1);  % denoised/interpolated n*.nii
  VT0 = res.image0(1); % original 
  [pth,nam] = spm_fileparts(VT0.fname);  %#ok<ASGLU>

  % voxel size parameter
  vx_vol  = sqrt(sum(VT.mat(1:3,1:3).^2));    % voxel size of the processed image
  vx_volp = prod(vx_vol)/1000;
  voli    = @(v) (v ./ (pi * 4./3)).^(1/3);   % volume > radius

  d = VT.dim(1:3);

  stime2 = cat_io_cmd('  Update Segmentation','g5','',job.extopts.verb-1,stime2); 
  % cleanup with brain mask - required for ngaus [1 1 2 4 3 2] and R1/MP2Rage like data 
  
  % tpm brain mask
  YbA = zeros(d,'single');
  Vb = tpm.V(1); Vb.pinfo(3) = 0; Vb.dt=16; Vb.dat = single(exp(tpm.dat{1}) + exp(tpm.dat{2}) + exp(tpm.dat{3})); 
  
  for z=1:d(3)
    YbA(:,:,z) = spm_sample_vol(Vb,double(Yy(:,:,z,1)),double(Yy(:,:,z,2)),double(Yy(:,:,z,3)),1); 
  end
  if round(max(YbA(:))/Vb.pinfo(1)), YbA=YbA>0.1*Vb.pinfo(1); else YbA=YbA>0.1; end
  

  % add some distance around brainmask (important for bias!)
  YbA = YbA | cat_vol_morph(YbA & sum(P(:,:,:,1:2),4)>4 ,'dd',2.4,vx_vol);
  % transfer tissue outside the brain mask to head  ... 
  % RD 201807: I am not sure if this is a good idea. Please test this with children! 
  for i=1:3
    P(:,:,:,4) = P(:,:,:,4) + P(:,:,:,i) .* uint8(~YbA); 
    P(:,:,:,i) = P(:,:,:,i) .* uint8(YbA); 
  end
  clear YbA;
  
  % cleanup for high resolution data
  % Alghough the old cleanup is very slow for high resolution data, the   
  % reduction of image resolution removes spatial segmentation information. 
  if job.opts.redspmres==0 % already done in case of redspmres
    if max(vx_vol)<1.5 && mean(vx_vol)<1.3
      for i=1:size(P,4), [Pc1(:,:,:,i),RR] = cat_vol_resize(P(:,:,:,i)  ,'reduceV',vx_vol,job.extopts.uhrlim,32); end %#ok<AGROW>
      for i=1:size(P,4), [Pc2(:,:,:,i),BB] = cat_vol_resize(Pc1(:,:,:,i),'reduceBrain',vx_vol,2,sum(Pc1,4)); end %#ok<AGROW>
      Pc2 = cat_main_clean_gwc(Pc2,1);
      for i=1:size(P,4), Pc1(:,:,:,i) = cat_vol_resize(Pc2(:,:,:,i),'dereduceBrain',BB); end
      for i=1:size(P,4), P(:,:,:,i)   = cat_vol_resize(Pc1(:,:,:,i),'dereduceV',RR); end 
      clear Pc1 Pc2;
    end
  end
  
  % garantee probability 
  sP = (sum(single(P),4)+eps)/255;
  for k1=1:size(P,4), P(:,:,:,k1) = cat_vol_ctype(single(P(:,:,:,k1))./sP); end
  clear sP;

  % median in case of WMHs!
  WMth = double(max( clsint(2) , cat_stat_nanmedian(Ysrc(P(:,:,:,2)>192)) )); 
  if clsint(3)>clsint(2) % invers
    CMth = clsint(3); 
  else
    CMth = min( [  clsint(1) - diff([clsint(1),WMth]) , clsint(3) ]);
  end
  T3th = [ CMth , clsint(1) , WMth];


  %    ds('l2','',vx_vol,Ysrc./WMth,Yp0>0.3,Ysrc./WMth,Yp0,80)
  Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
  if isfield(res,'Ylesion') && sum(res.Ylesion(:)>0)
    res.Ylesion = cat_vol_ctype( single(res.Ylesion) .* (Yp0>0.2) ); 
    for k=1:size(P,4), Yl = P(:,:,:,k); Yl(res.Ylesion>0.5) = 0; P(:,:,:,k) = Yl; end  
    Yl = P(:,:,:,3); Yl(res.Ylesion>0.5) = 255; P(:,:,:,3) = Yl; clear Yl; 
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
  end
  if sum(Yp0(:)>0.3)<100 
    % this error often depends on a failed affine registration, where SPM
    % have to find the brain in the head or background
    BGth  = min(cat_stat_nanmean(Ysrc( P(:,:,:,end)>128 )),clsint(6));
    HDHth = clsint(5);
    HDLth = clsint(4);
    clsvol = nan(1,size(P,4)); for ci=1:size(P,4), Yct = P(:,:,:,ci)>128; clsvol(ci) = sum(Yct(:))*vx_volp; end; clear Yct; 
    if size(P,4)==6
        error('CAT:cat_main:SPMpreprocessing:emptySegmentation', ...
         sprintf(['Empty Segmentation: \n ' ...
          'Possibly the affine registration failed. Please check image orientation.\n' ...
          ' Tissue class:           %10s%10s%10s%10s%10s%10s\n' ...
          ' Rel. to image volume:   %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n' ...
          ' Rel. to brain volume:   %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n' ...
          ' Tissue intensity:       %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f'],...
          'BG','CSF','GM','WM','HDH','HDL', ...
          [ clsvol([6 3 1 2 4 5])/cat_stat_nansum(clsvol)*100, clsvol([6 3 1 2 4 5])/cat_stat_nansum(clsvol(1:3))*100, BGth,T3th,HDHth,HDLth]));  %#ok<SPERR>
    elseif size(P,4)==4 % skull-stripped
        error('CAT:cat_main:SPMpreprocessing:emptySegmentation', ...
         sprintf(['Empty Segmentation: \n ' ...
          'Possibly the affine registration failed. Please check image orientation.\n' ...
          ' Tissue class:           %10s%10s%10s%10s\n' ...
          ' Rel. to image volume:   %10.2f%10.2f%10.2f%10.2f\n' ...
          ' Rel. to brain volume:   %10.2f%10.2f%10.2f%10.2f\n' ...
          ' Tissue intensity:       %10.2f%10.2f%10.2f%10.2f'],...
          'BG','CSF','GM','WM', ...
          [ clsvol([4 3 1 2])/cat_stat_nansum(clsvol)*100, clsvol([4 3 1 2])/cat_stat_nansum(clsvol(1:3))*100, BGth,T3th]));  %#ok<SPERR>
    else
        error('CAT:cat_main:SPMpreprocessing:emptySegmentation', ['Empty Segmentation: ' ...
           'Possibly the affine registration failed. Please check image orientation.\n']); 
    end
  end

  
  
  %%
  Yp0(smooth3(cat_vol_morph(Yp0>0.3,'lo'))<0.5)=0; % not 1/6 because some ADNI scans have large "CSF" areas in the background 
  Yp0     = Yp0 .* cat_vol_morph(Yp0 & (Ysrc>WMth*0.05),'lc',2);
  Yp0toC  = @(Yp0,c) 1-min(1,abs(Yp0-c));

  % values are only used if errors occur
  cat_err_res.init.T3th = T3th; 
  cat_err_res.init.subjectmeasures.vol_abs_CGW = [prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),1)), ... CSF
                                                  prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),2)), ... GM 
                                                  prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),3)), ... WM
                                                  prod(vx_vol)/1000 .* sum(Yp0toC(Yp0(:),4))];  % WMH
  cat_err_res.init.subjectmeasures.vol_TIV     =  sum(cat_err_res.init.subjectmeasures.vol_abs_CGW); 
  cat_err_res.init.subjectmeasures.vol_rel_CGW =  cat_err_res.init.subjectmeasures.vol_abs_CGW ./ ...
                                                  cat_err_res.init.subjectmeasures.vol_TIV;
  [cat_err_res.init.Yp0,cat_err_res.init.BB] = cat_vol_resize(Yp0,'reduceBrain',vx_vol,2,Yp0>0.5); 
  cat_err_res.init.Yp0 = cat_vol_ctype(cat_err_res.init.Yp0/3*255);

% ### This can not be reached because the mask field is removed by SPM! ###
  if isfield(res,'msk') 
    Ybg = ~res.msk.dat; 
    P4  = cat_vol_ctype( single(P(:,:,:,6)) .* (Ysrc<T3th(2))  .* (Ybg==0) + single(P(:,:,:,4)) .* (Ybg<1) ); % remove air in head
    P5  = cat_vol_ctype( single(P(:,:,:,6)) .* (Ysrc>=T3th(2)) .* (Ybg==0) + single(P(:,:,:,5)) .* (Ybg<1) ); % remove air in head
    P6  = cat_vol_ctype( single(sum(P(:,:,:,4:5),4)) .* (Ybg==1) + single(P(:,:,:,6)) .* (Ybg>0) ); % add objects/artifacts to background
    P(:,:,:,4) = P4;
    P(:,:,:,5) = P5;
    P(:,:,:,6) = P6;
    clear P4 P5 P6 Ybg; 
  end
  
  %%
  stime2 = cat_io_cmd('  Update Skull-Stripping','g5','',job.extopts.verb-1,stime2); 
  if size(P,4)==4 % skull-stripped
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
    Yb   = Yp0>=0.5/3 & Ysrc>0; 
    Ybb  = cat_vol_ctype(Yb)*255; 

    P(:,:,:,6) = P(:,:,:,4); 
    P(:,:,:,4) = zeros(size(Yp0),'uint8');
    P(:,:,:,5) = zeros(size(Yp0),'uint8'); 
    res.lkp = [res.lkp 5 6];
    res.mn  = [res.mn(1:end-1),0,0,0];
    res.mg  = [res.mg(1:end-1);1;1;1];
    res.vr(1,1,numel(res.lkp)-1:numel(res.lkp)) = 0;
     
    [Ysrcb,Yp0,BB] = cat_vol_resize({Ysrc,Yp0},'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yp0>1/3);
    Yg   = cat_vol_grad(Ysrcb/T3th(3),vx_vol);
    Ydiv = cat_vol_div(Ysrcb/T3th(3),vx_vol);
    Yg   = cat_vol_resize(Yg ,'dereduceBrain',BB);
    Ydiv = cat_vol_resize(Ydiv ,'dereduceBrain',BB);
  elseif (job.extopts.experimental || (job.extopts.INV && any(sign(diff(T3th))==-1))) && ...
      job.extopts.gcutstr>0 %&& job.extopts.gcutstr<=1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% not used for 2 years >> remove this path in future!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % (sum( abs( (Ysrc(:)==0) - (Yp0(:)<0.5) ) ) / sum(Ysrc(:)==0)) < 0.1  || ...
  % use gcut2

    %   brad = voli(sum(Yp0(:)>0).*prod(vx_vol)/1000); 
    %   noise = nanstd(Ysrc(cat_vol_morph(Yp0>0.8,'o') & Yp0>0.99)/diff(T3th(2:3))); 

    Vl1 = spm_vol(job.extopts.cat12atlas{1});
    Yl1 = cat_vol_ctype(spm_sample_vol(Vl1,double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0));
    Yl1 = reshape(Yl1,size(Ysrc)); [D,I] = cat_vbdist(single(Yl1>0)); Yl1 = Yl1(I); clear D I;  %#ok<ASGLU>

    %%
    Ybg = P(:,:,:,6)>128; 
    T3ths = [min(min(min(single(Ysrc(P(:,:,:,6)>128))))),...
             min( cat_stat_nanmean(Ysrc(Ybg(:))) + 2*cat_stat_nanstd(Ysrc(Ybg(:))) , ...
              mean([cat_stat_nanmean(Ysrc(Ybg(:))),min(T3th)])), ...
             T3th, T3th(3) + cat_stat_nanmean(diff(T3th))]; clear Ybg;
    T3thx = [0,0.05,1,2,3,4];
    if T3th(1)>T3th(3), T3thx = [0,0.05,1,2,3,2]; T3ths(end) = T3ths(2); end; 
    [T3ths,si] = sort(T3ths);
    T3thx      = T3thx(si);
    Ym = Ysrc+0; 
    for i=numel(T3ths):-1:2
      M = Ysrc>T3ths(i-1) & Ysrc<=T3ths(i);
      Ym(M(:)) = T3thx(i-1) + (Ysrc(M(:)) - T3ths(i-1))/diff(T3ths(i-1:i))*diff(T3thx(i-1:i));
    end
    M  = Ysrc>=T3ths(end); 
    Ym(M(:)) = numel(T3ths)/6 + (Ysrc(M(:)) - T3ths(i))/diff(T3ths(end-1:end))*diff(T3thx(i-1:i));    
    Ym = Ym / 3; 
    clear M; 
    
    for k1=1:size(P,4), Ycls{k1} = P(:,:,:,k1); end 
    Yb = (Ycls{1}+Ycls{2}+Ycls{3})>128; 
    %Yb = cat_main_gcut(Ym,Yp0>0.1,Ycls,Yl1,false(size(Ym)),vx_vol,...
    %  struct('gcutstr',0.1,'verb',0,'LAB',job.extopts.LAB,'LASstr',0,'red',1)); 
    clear Ycls; 
    Ybb  = cat_vol_ctype(cat_vol_smooth3X(Yb,2)*255); 
    
    [Ysrcb,Yp0,BB] = cat_vol_resize({Ysrc,Yp0},'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yp0>1/3);
    Ysrcb = max(0,min(Ysrcb,max(T3th)*2));
    Yg   = cat_vol_grad(Ysrcb/T3th(3),vx_vol);
    Ydiv = cat_vol_div(Ysrcb/T3th(3),vx_vol);

    Yb   = smooth3(Yb)>0.5; 
    Yg   = cat_vol_resize(Yg   ,'dereduceBrain',BB);
    Ydiv = cat_vol_resize(Ydiv ,'dereduceBrain',BB);
  elseif job.extopts.gcutstr==0 
    % brain mask
    Ym  = single(P(:,:,:,3))/255 + single(P(:,:,:,1))/255 + single(P(:,:,:,2))/255;
    Yb   = (Ym > 0.5);
    Yb = cat_vol_morph(cat_vol_morph(Yb,'lo'),'c');
    Ybb  = cat_vol_ctype(cat_vol_smooth3X(Yb,2)*256); 

    [Ysrcb,BB] = cat_vol_resize({Ysrc},'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yb);
    Yg   = cat_vol_grad(Ysrcb/T3th(3),vx_vol);
    Ydiv = cat_vol_div(Ysrcb/T3th(3),vx_vol);
    Yg   = cat_vol_resize(Yg   ,'dereduceBrain',BB);
    Ydiv = cat_vol_resize(Ydiv ,'dereduceBrain',BB);

  elseif job.extopts.gcutstr==2
    % adaptive probability region-growing
    if debug, Po=P; end
    
    %% tissue treshholds depending on MR modality
    if T3th(1) < T3th(3) % T1
      cth = min(res.mn(res.lkp==3));
      cth = min(cth,T3th(1));
      cth = min(cth,cat_stat_nanmedian(Ysrc( cat_vol_morph( smooth3(P(:,:,:,3))>200 & ...
              Ysrc<sum(T3th(1:2).*[0.8 0.2]) ,'de',2,vx_vol) ) ));
    else
      cth = max(res.mn(res.lkp==3));
      cth = max(cth,T3th(1));
      cth = max(cth,cat_stat_nanmedian(Ysrc( cat_vol_morph( smooth3(P(:,:,:,3))>200 & ...
              Ysrc>sum(T3th(1:2).*[0.8 0.2]),'de',2,vx_vol) ) ));
    end  
    if max(res.lkp)==4
      bth = 0;
    else
      bth = min( res.mn(res.lkp==4) );
    end
    if T3th(1) < T3th(3) % T1: CSF<GM<WM
      tth(1,:)  = [ mean(T3th(1:2))                     mean(T3th(2:3))                ]; % GM
      tth(2,:)  = [ mean(T3th(2:3))                     T3th(3) + 0.25*diff(T3th(2:3)) ]; % WM
      tth(3,:)  = [ T3th(1) - 0.25*diff([bth,T3th(1)])  sum(T3th(1:2) .* [0.75 0.25])  ]; % CSF
    elseif T3th(1) > T3th(3) % T2/PD: WM<GM<CSF ... not tested
      tth(1,:)  = [ mean(T3th(2:3))                 mean(T3th(1:2))                ];
      tth(2,:)  = [ T3th(3) - 0.25*diff(T3th(2:3))  mean(T3th(2:3))                ];
      tth(3,:)  = [ sum(T3th(1:2) .* [0.75 0.25])   sum(T3th(1:2) .* [0.75 0.25]) ];
    else % other contrast
      tth(1,:)  = [ T3th(1) - std(T3th(1:2))        T3th(1) + std(T3th(1:2)) ];
      tth(2,:)  = [ T3th(2) - std(T3th(1:2))        T3th(1) + std(T3th(2:3)) ];
      tth(3,:)  = [ T3th(3) - std(T3th(2:3))        T3th(1) + std(T3th(2:3)) ];
    end
    
    
    %% CSF mask
    %  Yc .. Combination of the CSF probability map and intensity map to 
    %        avoid meninges (especially in older subjectes) at the outer
    %        boundary. Due to failed registration we directly use the 
    %        brain mask Yb in the center of the brain, ie. we allow brain 
    %        tissue in the CSF far from the skull. 
    if T3th(1) < T3th(3)
      Yc   = single(P(:,:,:,3))/255 .* ...
        max(0,min(1,1 - ( max(0, Ysrc - cth) / abs( mean(res.mn(res.lkp(:)==1).*res.mg(res.lkp(:)==1)' ) - cth ) + ...
                          max(0,-Ysrc + cth) / abs( mean(res.mn(res.lkp(:)==4).*res.mg(res.lkp(:)==4)' ) - cth ) ) ));
    else
      Yc   = single(P(:,:,:,3))/255 .* ...
        max(0,min(1,1 - ( max(0,Ysrc - cth) / abs( min(res.mn(res.lkp==3)) - mean(res.mn(res.lkp==1)) )) ));
    end
    Ycg  = (single(P(:,:,:,3))/255 + single(P(:,:,:,1))/255 - Yc) .* ...
      max(0,min(1,1 - (0.5*abs(Ysrc - T3th(2)) / abs(diff(T3th(1:2))) ) ));
    
    
    % improved brain mask by region growing
    %  Yb .. improving the brain mask is necessary in case of missing
    %        structures (e.g. in children) or failed registration where 
    %        the TPM does not fit well and the CSF include brain tissue
    %        simply by chance.
    BGth = min(cat_stat_nanmean(Ysrc( P(:,:,:,6)>128 )),clsint(6));
    Yg   = cat_vol_grad((Ysrc - BGth)/diff([BGth,T3th(3)]),vx_vol);
    BVth = abs(diff(T3th(1:2:3)) / abs(T3th(3)) * 3);   % avoid blood vessels (high gradients) 
    RGth = abs(diff(T3th(2:3))   / abs(T3th(3)) * 0.1); % region growing threshold
    
    
    
    
    %% initial brain mask by region-growing
% Todo: T2/PD, skull-stripped
    %  as CSF/GM+GM+WM without blood vessels (Yg<0.5) 
    Yb   = min(1,single(P(:,:,:,1))/255 + single(P(:,:,:,2))/255 + Ycg/2); 
    Yb   = cat_vol_median3(Yb,Yb>0,true(size(Yb)));
    if T3th(1) < T3th(3) % T1 
      Yb   = min(1,Yb + cat_vol_morph(smooth3(sum(P(:,:,:,1:3),4))>253,'ldc',1));
    else % T2/PD ... more tolerant
      Yb   = min(1,Yb + cat_vol_morph(smooth3(sum(P(:,:,:,1:3),4))>200,'ldc',1));
    end
    Yb   = cat_vol_morph(Yb>0.8,'ldo',1.9,vx_vol);
    
    % mask for region growing and WM-GM region growing
    Yb2  = single(cat_vol_morph(Yb,'de',1.9,vx_vol)); 
    if T3th(1) < T3th(3) % T1 
      Yh   = (Yb2<0.5) & (Ysrc<sum(T3th(2:3).*[0.75 0.25]) | Ysrc>(T3th(3)*1.2) | Yg>BVth);
    else
      Yh   = (Yb2<0.5) & (Ysrc<mean([T3th(3),BGth]) | Ysrc>sum(T3th(2:3).*[0.75 0.25]) | Yg>BVth);
    end
    Yh   = cat_vol_morph(Yh,'ldc',1,vx_vol); 
    Yh   = cat_vol_morph(Yh,'de',1,vx_vol); Yb2(Yh) = nan; if ~debug, clear Yh; end
    if T3th(1) < T3th(3) % T1 
      [Yb2,YD] = cat_vol_downcut(Yb2,Ysrc/T3th(3),RGth/2); clear Yb2; %#ok<ASGLU>
    else
      [Yb2,YD] = cat_vol_downcut(Yb2,1 - Ysrc/T3th(3),RGth/2); clear Yb2; %#ok<ASGLU>
    end
    Yb(YD<400/mean(vx_vol)) = 1; clear YD; 
    Yb(smooth3(Yb)<0.5) = 0; 
    Yb   = cat_vol_morph(Yb,'ldo',1.9,vx_vol);
    
    % GM-CSF region
    Yb2  = single(cat_vol_morph(Yb,'de',1.9,vx_vol)); 
    if T3th(1) < T3th(3)
      Yh   = (Yb2<0.5) & (Ysrc<sum(T3th(1:2).*[0.9 0.1]) | sum(P(:,:,:,4:6),4)>250 | ...
              Ysrc>cat_stat_nanmean(T3th(3)*1.2) | Yg>BVth); 
    else
      Yh   = (Yb2<0.5) & (Ysrc>sum(T3th(1:2).*[0.9 0.1]) | sum(P(:,:,:,4:6),4)>250 | ...
              Ysrc<(T3th(3) - sum(T3th(2:3).*[0.5 0.5])) | Yg>BVth); 
    end
    Yh   = cat_vol_morph(Yh,'ldc',1) | cat_vol_morph(~Yb,'de',10,vx_vol); 
    Yh   = cat_vol_morph(Yh,'de',1,vx_vol);  Yb2(Yh) = nan; if ~debug, clear Yh; end
    if T3th(1) < T3th(3) % T1 
      [Yb2,YD] = cat_vol_downcut(Yb2,Ysrc/T3th(3),-RGth); clear Yb2; %#ok<ASGLU>
    else
      [Yb2,YD] = cat_vol_downcut(Yb2,1 - Ysrc/T3th(3),-RGth); clear Yb2; %#ok<ASGLU>
    end
    Yb(YD<400/mean(vx_vol)) = 1; clear YD; 
    Yb(smooth3(Yb)<0.5) = 0; Yb(smooth3(Yb)>0.5) = 1; 
    Yb   = cat_vol_morph(Yb,'ldo',1.9,vx_vol);
    Yb   = cat_vol_morph(Yb,'lc');
   
    
    
    %% CSF mask 2 with Yb
    %  Yc .. Combination of the CSF probability map and intensity map to 
    %        avoid meninges (especially in older subjectes) at the outer
    %        boundary. Due to failed registration we directly use the 
    %        brain mask Yb in the center of the brain, ie. we allow brain 
    %        tissue in the CSF far from the skull. 
    cth2 = min(cth,cat_stat_nanmedian( Ysrc( smooth3( ...
      cat_vol_morph(P(:,:,:,3)>200 & Yb & Yg<0.2,'de',1.5))>0.9 ) ));
    if T3th(1) < T3th(3)
      Yc   = single(P(:,:,:,3) )/255 .* ...
        max(Yb & Yg<BVth & Ysrc>tth(3,1) & Ysrc<tth(3,2), ...
        min(1,1 - ( max(0, Ysrc - cth2) / (2 * abs( mean(res.mn(res.lkp(:)==1).*res.mg(res.lkp(:)==1)' - cth2) )) + ...
                    max(0,-Ysrc + cth2) / (2 * abs( mean(res.mn(res.lkp(:)==4).*res.mg(res.lkp(:)==4)' - cth2) )) ) ));
    else
      Yc   = single(P(:,:,:,3))/255 .* ...
        max(Yb & Yg<BVth & Ysrc>tth(3,1) & Ysrc<tth(3,2), ...
        min(1,1 - ( min(0,Ysrc - cth2) / (2 * abs( mean(res.mn(res.lkp==1)) - cth2) ) ) ));
    end
    % single(P(:,:,:,3))/255 .* ...
    %  max(Yb & Ysrc>cth2,min(1,1 - ( abs(Ysrc - cth2) / ...
    %  (2 * abs( mean(res.mn(res.lkp==1)) - cth2) ) ) ));

     
    
    
    %% update GM map (include tissue that was previously labeled as head) 
    for i=1:3
      % smaller mask in case of WM
      if i==2, Ybt = cat_vol_morph(Yb,'de',1.5,vx_vol); else Ybt = Yb; end
      % smooth mask in case of GM
      if i==1
        Ytmp = single(P(:,:,:,4)) .* smooth3(Ybt & Ysrc>tth(i,1) & Ysrc<tth(i,2) & Yg<BVth); 
      else %if i==2
        Ytmp = single(P(:,:,:,4)) .*        (Ybt & Ysrc>tth(i,1) & Ysrc<tth(i,2) & Yg<BVth); 
      %{
      elseif i==3
        % add CSF around mask only if the intensity fit very well
        Ytmp = single(P(:,:,:,4)) .* Yc .* ...
          max(Yb & Ysrc>tth(i,1) & Ysrc<tth(i,2) & Yg<BVth, ...
          min(1,1 - ( abs(Ysrc - cth2) / (4 * abs( mean(res.mn(res.lkp==1)) - cth2) ) ) ));
      %}
     end
      % tissue transfer
      P(:,:,:,i) = cat_vol_ctype(single(P(:,:,:,i)) + Ytmp);
      P(:,:,:,4) = cat_vol_ctype(single(P(:,:,:,4)) - Ytmp);
    end
    % transfer tissue in case of bad TPM matching
    for i=[1,2,4]
      % smaller mask in case of WM
      if i==2, Ybt = cat_vol_morph(Yb,'de',1.5,vx_vol); else Ybt = Yb; end
      % smooth mask in case of GM
      if i==1
        Ytmp = single(P(:,:,:,3)) .* smooth3(Ybt & Ysrc>tth(i,1) & Ysrc<tth(i,2) & Yg<BVth);
      elseif i==2
        Ytmp = single(P(:,:,:,3)) .*        (Ybt & Ysrc>tth(i,1) & Ysrc<tth(i,2) & Yg<BVth); 
      elseif i==4
        Ytmp = single(P(:,:,:,3)) .* (1 - min(1,Ybt + Yc)); 
      end
      % tissue transfer
      P(:,:,:,i) = cat_vol_ctype(single(P(:,:,:,i)) + Ytmp);
      P(:,:,:,3) = cat_vol_ctype(single(P(:,:,:,3)) - Ytmp);
    end

    
    %% qa control parameter?
    % SPMsegfit = mean(P(Yb(:))-Po(Yb(:))); 
    
    
    %% create brain level set map
    %  Ym .. combination of brain tissue and CSF that is further corrected
    %        for noise (median) and smoothness (Laplace) an finally 
    %        threshholded 
    Ym  = min(1,Yc + single(P(:,:,:,1))/255 + single(P(:,:,:,2))/255 + Yb);
    Ym  = cat_vol_median3(Ym,Ym>0 & Ym<1);  % remove noise 
    % region-growing 
    Ym2 = Ym; Ym2(Ym2==0)=nan;
    [Ym2,YD] = cat_vol_downcut(single(Ym2>0.99),Ym2,0.01); clear Ym2; %#ok<ASGLU>
    Ym(YD>400/mean(vx_vol))=0; clear YD; 
    Ym(cat_vol_morph(Ym>0.95,'ldc',1)) = 1; 
    Ym(cat_vol_morph(Yb,'e') & Ym<0.9 & Yc<0.25) = 0;
    Ym  = Ym .* cat_vol_morph(Ym>0.5,'ldo',2);  % remove extensions (BV, eye)
    Ym = cat_vol_laplace3R(Ym,Ym>0.1 & Ym<0.9,0.2); % smooth mask
    Ym = cat_vol_laplace3R(Ym,Ym<0.25,0.2); 
    Ym(cat_vol_morph(Yb,'e') & Ym<0.9 & Yc<0.25) = 0;
    
    %% cutting parameter
    %  This is maybe a nice parameter to control the CSF masking.
    %  And even better we can use a surface to find the optimal value. :)
    cutstr    = 1.0; % 0.85; 
    cutstrs   = linspace(0.2,0.8,4); % 0.05,0.35,0.65,0.95]; 
    cutstrval = nan(1,4); 
    if debug, cutstrsa = zeros(0,8); end
    if cutstr == 1 % auto
      for l=1:3
        for i=1:numel(cutstrs)
          if isnan( cutstrval(i) )
            S = isosurface(Ym,cutstrs(i),Ysrc/T3th(3)); 
            cutstrval(i) = cutstrs(i)/20 + ... % litte offset to get more CSF
              abs(mean(S.facevertexcdata) - 2*T3th(1)/T3th(3)) + std(S.facevertexcdata);
          end
        end
        [tmp,cutstrid] = sort(cutstrval); clear tmp; %#ok<ASGLU>
        if debug, cutstrsa  = [cutstrsa; cutstrs, cutstrval]; end %#ok<AGROW>
        cutstrs   = linspace(cutstrs(max(1,cutstrid(1)-1)),cutstrs(min(4,cutstrid(1)+1)),4);
        cutstrval = [cutstrval(max(1,cutstrid(1)-1)),nan,nan,cutstrval(min(4,cutstrid(1)+1))];
        
      end
      cutstr = cutstrs(cutstrid(1));
    end
    
    
    %% normalize this map depending on the cutstr parameter 
    Yb  = cat_vol_morph(cat_vol_morph(Ym > cutstr,'lo'),'c');
    Yb  = cat_vol_morph(Yb,'e') | (Ym>0.9) | (Yb & Yc>0.5);
    Yb(smooth3(Yb)<0.5)=0;
    Ybb = cat_vol_ctype( max(0,min(1,(Ym - cutstr)/(1-cutstr))) * 256); 
    Ym0 = Ybb; 
    
    %% estimate gradient (edge) and divergence maps
    [Ysrcb,BB] = cat_vol_resize({Ysrc},'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yb);
    Yg   = cat_vol_grad(Ysrcb/T3th(3),vx_vol);
    Ydiv = cat_vol_div(Ysrcb/T3th(3),vx_vol);
    Yg   = cat_vol_resize(Yg   ,'dereduceBrain',BB);
    Ydiv = cat_vol_resize(Ydiv ,'dereduceBrain',BB);
  else
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % T1 only > remove in future if gcut is removed too!
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    % old skull-stripping
    brad = voli(sum(Yp0(:)>0.5).*prod(vx_vol)/1000); 
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
    [Ysrcb,Yp0,BB] = cat_vol_resize({Ysrc,Yp0},'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yp0>1/3);
    %Ysrcb = max(0,min(Ysrcb,max(T3th)*2));
    BGth = min(cat_stat_nanmean(Ysrc( P(:,:,:,6)>128 )),clsint(6));
    Yg   = cat_vol_grad((Ysrcb-BGth)/diff([BGth,T3th(3)]),vx_vol);
    Ydiv = cat_vol_div((Ysrcb-BGth)/diff([BGth,T3th(3)]),vx_vol);
    Ybo  = cat_vol_morph(cat_vol_morph(Yp0>0.3,'lc',2),'d',brad/2/mean(vx_vol)); 
    BVth = diff(T3th(1:2:3))/abs(T3th(3))*1.5; 
    RGth = diff(T3th(2:3))/abs(T3th(3))*0.1; 
    Yb   = single(cat_vol_morph((Yp0>1.9/3) | (Ybo & Ysrcb>mean(T3th(2)) & Ysrcb<T3th(3)*1.5 & Yg<0.5),'lo',max(0,0.6/mean(vx_vol)))); 
    %% region-growing GM 1
    Yb(~Yb & (~Ybo | Ysrcb<cat_stat_nanmean(T3th(2)) | Ysrcb>cat_stat_nanmean(T3th(3)*1.2) | Yg>BVth))=nan;
    [Yb1,YD] = cat_vol_downcut(Yb,Ysrcb/T3th(3),RGth); clear Yb1; Yb(isnan(Yb))=0; Yb(YD<400/mean(vx_vol))=1; clear YD; %#ok<ASGLU> 
    Yb(smooth3(Yb)<0.5)=0; Yb = single(Yb | (Ysrcb>T3th(1) & Ysrcb<1.2*T3th(3) & cat_vol_morph(Yb,'lc',4)));
    %% region-growing GM 2
    Yb(~Yb & (~Ybo | Ysrcb<T3th(1) | Ysrcb>cat_stat_nanmean(T3th(3)*1.2) | Yg>BVth))=nan;
    [Yb1,YD] = cat_vol_downcut(Yb,Ysrcb/T3th(3),RGth/2); clear Yb1; Yb(isnan(Yb))=0; Yb(YD<400/mean(vx_vol))=1; clear YD; %#ok<ASGLU>
    Yb(smooth3(Yb)<0.5)=0; Yb = single(Yb | (Ysrcb>T3th(1) & Ysrcb<1.2*T3th(3) & cat_vol_morph(Yb,'lc',4)));
    %% region-growing GM 3
    Yb(~Yb & (~Ybo | Ysrcb<mean([BGth,T3th(1)]) | Ysrcb>cat_stat_nanmean(T3th(3)*1.2) | Yg>BVth))=nan; clear Ybo;
    [Yb1,YD] = cat_vol_downcut(Yb,Ysrcb/T3th(3),RGth/10); clear Yb1; Yb(isnan(Yb))=0; Yb(YD<400/mean(vx_vol))=1; clear YD; %#ok<ASGLU>
    Yb(smooth3(Yb)<0.5)=0; Yb(Yp0toC(Yp0*3,1)>0.9 & Yg<0.3 & Ysrcb>BGth & Ysrcb<T3th(2)) = 1; 
    %% ventricle closing
    [Ybr,Ymr,resT2] = cat_vol_resize({Yb>0,Ysrcb/T3th(3)},'reduceV',vx_vol,2,32); clear Ysrcb
    Ybr = Ybr | (Ymr<0.8 & cat_vol_morph(Ybr,'lc',6)); clear Ymr;  % large ventricle closing
    Ybr = cat_vol_morph(Ybr,'lc',2);                 % standard closing
    Yb  = Yb | cat_vol_resize(cat_vol_smooth3X(Ybr,2),'dereduceV',resT2)>0.7; clear Ybr
    Yb  = smooth3(Yb)>0.5; 
    Ybb = cat_vol_ctype(cat_vol_smooth3X(Yb,2)*255); 
    Yb   = cat_vol_resize(Yb   , 'dereduceBrain' , BB);
    Ybb  = cat_vol_resize(Ybb  , 'dereduceBrain' , BB);
    Yg   = cat_vol_resize(Yg   , 'dereduceBrain' , BB);
    Ydiv = cat_vol_resize(Ydiv , 'dereduceBrain' , BB);
    clear Ybo;
  end
  clear Ysrcb Yp0; 

  
  
  %% save brainmask using SPM12 segmentations for later use
  if ~exist('Ym0','var'),
    Ym0 = single(P(:,:,:,3))/255 + single(P(:,:,:,1))/255 + single(P(:,:,:,2))/255;
  end
  Yb0 = (Ym0 > min(0.5,max(0.25, job.extopts.gcutstr))); clear Ym0
  Yb0 = cat_vol_morph(cat_vol_morph(Yb0,'lo'),'c');


%%
  stime2 = cat_io_cmd('  Update probability maps','g5','',job.extopts.verb-1,stime2);
  if ~(job.extopts.INV && any(sign(diff(T3th))==-1))
    %% Update probability maps
    % background vs. head - important for noisy backgrounds such as in MT weighting
    if size(P,4)==4 % skull-stripped
      Ybg = ~Yb;
    else
      if sum(sum(sum(P(:,:,:,6)>240 & Ysrc<cat_stat_nanmean(T3th(1:2)))))>10000
        Ybg = P(:,:,:,6); 
        [Ybgr,Ysrcr,resT2] = cat_vol_resize({Ybg,Ysrc},'reduceV',vx_vol,2,32); 
        Ybgrth = max(cat_stat_nanmean(Ysrcr(Ybgr(:)>128)) + 2*std(Ysrcr(Ybgr(:)>128)),T3th(1));
        Ybgr = cat_vol_morph(cat_vol_morph(cat_vol_morph(Ybgr>128,'d') & Ysrcr<Ybgrth,'lo',1),'lc',1);
        Ybg  = cat_vol_resize(cat_vol_smooth3X(Ybgr,1),'dereduceV',resT2); 
        clear Ysrcr Ybgr; 
      else
        Ybg = ~Yb;
      end
    end
    P4   = cat_vol_ctype( single(P(:,:,:,6)) .* (Ysrc<T3th(2))  .* (Ybg==0) + single(P(:,:,:,4)) .* (Ybg<1) ); % remove air in head
    P5   = cat_vol_ctype( single(P(:,:,:,6)) .* (Ysrc>=T3th(2)) .* (Ybg==0) + single(P(:,:,:,5)) .* (Ybg<1) ); % remove air in head
    P6   = cat_vol_ctype( single(sum(P(:,:,:,4:5),4)) .* (Ybg==1) + single(P(:,:,:,6)) .* (Ybg>0) ); % add objects/artifacts to background
    P(:,:,:,4) = P4;
    P(:,:,:,5) = P5;
    P(:,:,:,6) = P6;
    clear P4 P5 P6;

    %% correct probability maps to 100% 
    sumP = cat_vol_ctype(255 - sum(P(:,:,:,1:6),4));
    P(:,:,:,1) = P(:,:,:,1) + sumP .* uint8( Ybg<0.5  &  Yb & Ysrc>cat_stat_nanmean(T3th(1:2)) & Ysrc<cat_stat_nanmean(T3th(2:3)));
    P(:,:,:,2) = P(:,:,:,2) + sumP .* uint8( Ybg<0.5  &  Yb & Ysrc>=cat_stat_nanmean(T3th(2:3)));
    P(:,:,:,3) = P(:,:,:,3) + sumP .* uint8( Ybg<0.5  &  Yb & Ysrc<=cat_stat_nanmean(T3th(1:2)));
    P(:,:,:,4) = P(:,:,:,4) + sumP .* uint8( Ybg<0.5  & ~Yb & Ysrc<T3th(2));
    P(:,:,:,5) = P(:,:,:,5) + sumP .* uint8( Ybg<0.5  & ~Yb & Ysrc>=T3th(2));
    P(:,:,:,6) = P(:,:,:,6) + sumP .* uint8( Ybg>=0.5 & ~Yb );
    clear Ybg sumP;

    %% head to WM 
    % Undercorrection of strong inhomogeneities in high field scans 
    % (>1.5T) can cause missalignments of the template and therefore 
    % miss classifications of the tissues that finally avoid further 
    % corrections in by LAS. 
    % Typically the alginment failed in this cases because the high 
    % intensities next to the head that were counted as head and not
    % corrected by SPM.
    % e.g. HR075, Magdeburg7T, SRS_SRS_Jena_DaRo81_T1_20150320-191509_MPR-08mm-G2-bw330-nbc.nii, ...
    Ywm = single(P(:,:,:,2)>128 & Yg<0.3 & Ydiv<0.03); Ywm(Ybb<128 | (P(:,:,:,1)>128 & abs(Ysrc/T3th(3)-2/3)<1/3) | Ydiv>0.03) = nan;
    [Ywm1,YD] = cat_vol_downcut(Ywm,1-Ysrc/T3th(3),0.02); Yb(isnan(Yb))=0; Ywm(YD<300)=1; Ywm(isnan(Ywm))=0; clear Ywm1 YD; %#ok<ASGLU>
    Ywmc = uint8(smooth3(Ywm)>0.7);
    %%
    Ygmc = uint8(cat_vol_morph(Ywmc,'d',2) & ~Ywmc & Ydiv>0 & Yb & cat_vol_smooth3X(Yb,8)<0.9);
    P(:,:,:,[1,3:6]) = P(:,:,:,[1,3:6]) .* repmat(1-Ywmc,[1,1,1,5]);
    P(:,:,:,2:6)     = P(:,:,:,2:6)     .* repmat(1-Ygmc,[1,1,1,5]);
    P(:,:,:,1)       = max(P(:,:,:,1),255*Ygmc);
    P(:,:,:,2)       = max(P(:,:,:,2),255*Ywmc);
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
    clear Ygmc Ywmc; 

    %% head to GM ... important for children
    [Ywmr,Ybr,resT2] = cat_vol_resize({Ywm,Yb},'reduceV',vx_vol,2,32); 
    Ygm = cat_vol_morph(Ywmr>0.5,'d',3) & (cat_vol_morph(~Ybr,'d',3) | cat_vol_morph(Ybr,'d',1)); clear Ybr Ywmr;  % close to the head
    Ygm = cat_vol_resize(single(Ygm),'dereduceV',resT2)>0.5;
    Ygm = Ygm & Yp0<2/3 & Yb & Yg<cat_stat_nanmean(Yg(P(:,:,:,1)>64)) & Ydiv<cat_stat_nanmean(Ydiv(P(:,:,:,1)>64)); % add GM with low SPM prob ... 
    Ygm = Ygm & (Ysrc>cat_stat_nansum(T3th(1:2).*[0.5 0.5])) & (Ysrc<cat_stat_nansum(T3th(2:3).*[0.2 0.8])); % but good intensity
    Ygm(smooth3(Ygm)<0.5)=0; 
    clear Yg Ydiv;
    %%
    Ygm = uint8(Ygm); 
    P(:,:,:,5) = P(:,:,:,5) .* (1-Ygm);
    P(:,:,:,3) = P(:,:,:,3) .* (1-Ygm);
    P(:,:,:,2) = P(:,:,:,2) .* (1-Ygm);
    P(:,:,:,1) = cat_vol_ctype(single(P(:,:,:,1)) + 255*single(Ygm));
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
    clear Ywm Ygm;

    %% remove brain tissues outside the brainmask ...
    % tissues > skull (within the brainmask)
    Yhdc = uint8(smooth3( Ysrc/T3th(3).*(Ybb>cat_vol_ctype(0.2*255)) - Yp0 )>0.5); 
    sumP = sum(P(:,:,:,1:3),4); 
    P(:,:,:,4)   =  cat_vol_ctype( single(P(:,:,:,4)) + sumP .* ((Ybb<=cat_vol_ctype(0.05*255)) | Yhdc ) .* (Ysrc<T3th(2)));
    P(:,:,:,5)   =  cat_vol_ctype( single(P(:,:,:,5)) + sumP .* ((Ybb<=cat_vol_ctype(0.05*255)) | Yhdc ) .* (Ysrc>=T3th(2)));
    P(:,:,:,1:3) =  P(:,:,:,1:3) .* repmat(uint8(~(Ybb<=cat_vol_ctype(0.05*255)) | Yhdc ),[1,1,1,3]);
    clear sumP Yp0 Yhdc; 
  end
  clear Ybb;
  

  %% MRF
  % Used spm_mrf help and tested the probability TPM map for Q without good results.         
  nmrf_its = 0; % 10 interations better to get full probability in thin GM areas 
  spm_progress_bar('init',nmrf_its,['MRF: Working on ' nam],'Iterations completed');
  if isfield(res,'mg'), Kb = max(res.lkp); else Kb = size(res.intensity(1).lik,2); end
  G   = ones([Kb,1],'single');
  vx2 = single(sum(VT.mat(1:3,1:3).^2));
  % P = zeros([d(1:3),Kb],'uint8');
  % P = spm_mrf(P,Q,G,vx2); % init: transfer data from Q to P 
  if 0
    %% use TPM as Q
    Q = zeros(size(P),'uint8');
    for di=1:6
      vol = cat_vol_ctype(spm_sample_vol(tpm.V(di),...
        double(Yy(:,:,:,1)),double(Yy(:,:,:,2)),double(Yy(:,:,:,3)),0)*255,'uint8');
      Q(:,:,:,di) = reshape(vol,d);
    end
  end
  for iter=1:nmrf_its,
      P = spm_mrf(P,single(P),G,vx2); % spm_mrf(P,Q,G,vx2);
      spm_progress_bar('set',iter);
  end

  spm_progress_bar('clear');
  for k1=1:size(P,4)
      Ycls{k1} = P(:,:,:,k1);
  end
  clear Q P q q1 Coef b cr N lkp n wp M k1


  if job.extopts.verb>2
    % save information for debugging and OS test
    % input variables + bias corrected, bias field, class image
    % strong differences in bias fields can be the result of different 
    % registration > check 'res.image.mat' and 'res.Affine'
    [pth,nam] = spm_fileparts(res.image0(1).fname); 
    tpmci  = 1;
    tmpmat = fullfile(pth,reportfolder,sprintf('%s_%s%02d%s.mat',nam,'write',tpmci,'postbias'));
    save(tmpmat,'res','tpm','job','Ysrc','Ybf','Ycls');
  end

  clear Ybf
  stime2 = cat_io_cmd(' ','g5','',job.extopts.verb-1,stime2); 
  fprintf('%5.0fs\n',etime(clock,stime));



end
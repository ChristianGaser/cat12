function [Affine,Yb,Ymi,Ym0] = cat_run_job_APRGs(Ysrc,Ybg,VF,Pb,Pbt,Affine,vx_vol,obj,job)
% ______________________________________________________________________
%  Skull-stripping subfunction APRG (adaptive probability region-growing)
%  of cat_main_updateSPM.
%
%  [Yb,Ymi,Ym0] = cat_run_job_APRGs(Ysrc,Ybg,VF,Pb,Pbt,Affine,vx_vol)
%
%   Ysrc  .. original input images
%   Ybg   ..
%   res   .. SPM preprocessing structure
%   Yb    .. binary brain mask
%   Ym0   .. probability brain mask
%   Yg    .. absolute gradient map (eg. for tissue edges)
%   Ydiv  .. divergence maps (eg. for blood vessels)
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id: cat_main_APRG.m 1392 2018-11-15 20:24:07Z gaser $

  dbs = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end
 
  % load default brain mask
  VFa = VF; VFa.mat = Affine * VF.mat; %Fa.mat = res0(2).Affine * VF.mat;
  if isfield(VFa,'dat'), VFa = rmfield(VFa,'dat'); end
  [Vmsk,Yb0] = cat_vol_imcalc([VFa,spm_vol(Pb)],Pbt,'i2',...
    struct('interp',2,'verb',0,'mask',-1)); clear Vmsk %#ok<ASGLU>
  
  % remove areas far away from the brain and use low resolution
  Ymi  = Ysrc; 
  [Ysrc,BB] = cat_vol_resize(Ysrc,'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yb0,10);
  Ybg       = cat_vol_resize(Ybg ,'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yb0,10);
  Yb0       = cat_vol_resize(Yb0 ,'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yb0,10);
  [Ysrc,rV] = cat_vol_resize(Ysrc,'reduceV',vx_vol,2,32);
  Ybg       = cat_vol_resize(Ybg ,'reduceV',vx_vol,2,32);
  Yb0       = cat_vol_resize(Yb0 ,'reduceV',vx_vol,2,32);
  
  vx_vol    = rV.vx_volr;
  
  %% estimate tissue thresholds for intensity normalization 
  bgth = kmeans3D(Ysrc(Ybg)); 
  T5th = kmeans3D(Ysrc(Yb0>0.5),5); T3th = T5th(1:2:end);
  [Ysrc,th] = cat_stat_histth(Ysrc);
  Ym = (Ysrc - th(1)) ./ (T3th(3) - th(1)); 
  Yg = cat_vol_grad(Ym,vx_vol) ./ Ym;   
  
  %% create background if is not available (no APP)
  if isempty(Ybg)
    %% background
    Ydiv  = cat_vol_div(Ym,vx_vol,2*ones(1,3)) ./ (Ym+eps); Ydiv(Ydiv==0)=eps;
    
    %stime = cat_io_cmd('  Estimate background','g5','',opt.verb,stime);
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
  end
    

  %% improved brain mask by region growing
  %  Yb .. improving the brain mask is necessary in case of missing
  %        structures (e.g. in children) or failed registration where 
  %        the TPM does not fit well and the CSF include brain tissue
  %        simply by chance.
  BGth = double(cat_stat_nanmean(Ysrc( Ybg )));
  BVth = double(abs(diff(T3th(1:2:3)) / abs(T3th(3)) * 3));   % avoid blood vessels (high gradients) 
  RGth = double(abs(diff(T3th(2:3))   / abs(T3th(3)) * 0.1)); % region growing threshold


  %% initial brain mask by region-growing
  % mask for region growing and WM-GM region growing
  Yb2  = single(cat_vol_morph(Yb0>0.5 | (Yb0.*Ym)>0.25 & (Yb0.*Ym)<1.3,'de',2.9,vx_vol)); Yb=false(size(Yb0));
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

  %% GM-CSF region
  Yb2  = single(cat_vol_morph(Yb,'de',2.9,vx_vol)); 
  if T3th(1) < T3th(3)
    Yh   = (Yb2<0.5) & (Ysrc<sum(T3th(1:2).*[0.9 0.1]) | ...
            Ysrc>cat_stat_nanmean(T3th(3)*1.2) | Yg>BVth); 
  else
    Yh   = (Yb2<0.5) & (Ysrc>sum(T3th(1:2).*[0.9 0.1]) | ...
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

  
  
  
  
  %% CSF region
  if T3th(1) < T3th(3)
    %% create brain level set map
    %  Ym .. combination of brain tissue and CSF that is further corrected
    %        for noise (median) and smoothness (Laplace) an finally 
    %        threshholded 
    Ybd  = cat_vbdist(single(smooth3(Yb)>0.5),~Ybg,vx_vol);
    mnhd = kmeans3D( Ybd( cat_vol_morph(Ybg,'d') & Ybd<10 ) ); 
    Ybgd = cat_vbdist(single(Ybg | Ybd>mnhd*2),~Yb,vx_vol);
    Ymx  = Ybgd./(Ybd+Ybgd);
    %
    Yb2  = single(cat_vol_morph(Yb,'de',1.9,vx_vol)); 
    Yh   = (Yb2<0.5) & (Ymx<0.5 | (Ysrc>sum(T3th(1:2).*[0.5 0.5])) & cat_vol_morph(~Ybg,'de',2)); 
    Yh   = cat_vol_morph(Yh,'ldc',1) | cat_vol_morph(~Yb,'de',10,vx_vol); 
    Yh   = cat_vol_morph(Yh,'de',1,vx_vol);  Yb2(smooth3(Yh)>0.9) = nan; if ~debug, clear Yh; end
    %%
    [Yb2,YD] = cat_vol_downcut(Yb2,Ysrc/T3th(3),-RGth); clear Yb2; %#ok<ASGLU>
    Yb(YD<400/mean(vx_vol)) = 1; clear YD; 
    Yb(smooth3(Yb)<0.5) = 0; Yb(smooth3(Yb)>0.5) = 1; 
    Yb   = cat_vol_morph(Yb,'ldo',1.9,vx_vol);
    Yb   = cat_vol_morph(Yb,'lc');
  end

  
  %% create brain level set map
  %  Ym .. combination of brain tissue and CSF that is further corrected
  %        for noise (median) and smoothness (Laplace) an finally 
  %        threshholded 
  Ybd  = cat_vbdist(single(smooth3(Yb)>0.5),~Ybg,vx_vol);
  mnhd = kmeans3D( Ybd( cat_vol_morph(Ybg,'d') & Ybd<10 ) ); 
  Ybgd = cat_vbdist(single(Ybg | Ybd>mnhd),~Yb,vx_vol);
  Ymx  = Ybgd./(Ybd+Ybgd);
    
  
  %% cutting parameter
  %  This is maybe a nice parameter to control the CSF masking.
  %  And even better we can use a surface to find the optimal value. :)
  cutstr    = 1.0; % 0.85; 
  cutstrs   = linspace(0.3,0.95,4); % 0.05,0.35,0.65,0.95]; 
  cutstrval = nan(1,4); 
  if debug, cutstrsa = zeros(0,8); end
  Ysrc2 = (Ysrc>T3th(1)) .* (abs(Ysrc - T3th(1))/(T3th(2) - T3th(1))) + ...
          (Ysrc<T3th(1)) .* (abs(Ysrc - T3th(1))/(T3th(1) - BGth)) ;
        Ysrc2 = smooth3(Ysrc2);
  if cutstr == 1 % auto
    for l=1:3
      for i=1:numel(cutstrs)
        if isnan( cutstrval(i) )
          S = isosurface(Ymx,cutstrs(i),Ysrc2); 
          cutstrval(i) = cutstrs(i)/10 + cat_stat_nanmean(S.facevertexcdata.^2).^0.5; %... % litte offset to get more CSF
            %mean(S.facevertexcdata) + std(S.facevertexcdata);
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
  Yb  = cat_vol_morph(cat_vol_morph(Ymx > cutstr,'lo'),'c');
  Yb  = cat_vol_morph(Yb,'e') | (Ymx>0.9);
  Yb(smooth3(Yb)<0.5)=0;
  Ybb = cat_vol_ctype( max(0,min(1,(Ymx - cutstr)/(1-cutstr))) * 256); 
  Ym0 = Ybb; 

  
  %% estimate gradient (edge) and divergence maps
  Yb  = cat_vol_resize(Yb ,'dereduceV',rV);
  Ym0 = cat_vol_resize(Ym0,'dereduceV',rV);
  Yb  = cat_vol_resize(Yb ,'dereduceBrain',BB);
  Ym0 = cat_vol_resize(Ym0,'dereduceBrain',BB);
  
  %% intensity normalization
  Tth.T3thx = [bgth T5th(1) T5th(3) T5th(5) T5th(5)+diff(T5th(3:2:5))];
  Tth.T3th   = 0:1/3:4/3;
  Ymi = cat_main_gintnormi(Ymi/3,Tth);
  
  
  %% create obj
  obj2 = obj; 
  obj2.image.pinfo   = repmat([255;0],1,size(Yb,3));
  obj2.image.dt      = [spm_type('UINT8') spm_platform('bigend')];
  obj2.image.dat     = cat_vol_ctype(Ymi*255 .* Yb);  
  if isfield(obj2.image,'private'), obj2.image = rmfield(obj2.image,'private'); end
  obj2.samp = 1.5; 
  obj2.fwhm = 1;
  if isfield(obj2,'tpm'), obj2 = rmfield(obj2,'tpm'); end
  obj2.tpm = spm_load_priors8(obj.tpm.V(1:3));
  % do registration
  Affine   = spm_maff8(obj2.image,obj.samp,obj2.fwhm,obj2.tpm,Affine,job.opts.affreg,80);
  spm_progress_bar('Clear');
end
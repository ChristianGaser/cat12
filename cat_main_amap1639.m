function [prob,indx,indy,indz,th] = cat_main_amap1639(Ymi,Yb,Yb0,Ycls,job,res)
% ______________________________________________________________________
%
% AMAP segmentation:
% Most corrections were done before and the AMAP routine is used with 
% a low level of iterations and no further bias correction, because
% some images get tile artifacts. 
%
% [prob,indx,indy,indz] = cat_main_amap1639(Ymi,Yb,Yb0,Ycls,job,res)
%
% prob .. new AMAP segmenation (4D)
% ind* .. index elements to asign a subvolume
% Ymi  .. local intensity normalized source image
% Yb   .. brain mask
% Yb0  .. origina brain mask 
% Ycls .. SPM segmentation 
% job  .. SPM/CAT parameter structure
% res  .. SPM segmentation structure
% th   .. AMAP treshholds
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$ 


  % this function adds noise to the data to stabilize processing and we
  % have to define a specific random pattern to get the same results each time
  if exist('rng','file') == 2, rng('default'); rng(0); else, rand('state',0); randn('state',0); end
  if ~isfield(job.extopts,'AMAPframing'), job.extopts.AMAPframing = 0; end

  % if there is a breakpoint in this file set debug=1 and do not clear temporary variables 
  dbs  = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

  % correct for harder brain mask to avoid meninges in the segmentation
  Ymib   = Ymi; Ymib(~Yb) = 0; 
  rf     = 10^4; Ymib = round(Ymib*rf)/rf;
  d      = size(Ymi); 
  vx_vol = sqrt(sum(res.image(1).mat(1:3,1:3).^2));
  
  % use framing
  %{
    sz = size(Yb); 
    [indx, indy, indz] = ind2sub(sz,find(Yb>0));
    tx = [min(indx) szx - max(indx)];
    ty = [min(indx) szy - max(indx)];
    tz = [min(indx) szz - max(indx)];
    tb = min( [ tx ty tz ]; 
  %}
  framing.tissue = 4; 
  framing.pve    = 1; 
  
  %  prepare data for segmentation
  if 1
    %% classic approach, consider the WMH!
    Kb2 = 4;
    cls2 = zeros([d(1:2) Kb2]);
    Yp0  = zeros(d,'uint8');
    for i=1:d(3)
        for k1 = 1:Kb2, cls2(:,:,k1) = Ycls{k1}(:,:,i); end
        % find maximum for reordered segmentations
        [maxi,maxind] = max(cls2(:,:,[3,1,2,4:Kb2]),[],3);
        k1ind = [1 2 3 1 0 0 1 0]; 
        for k1 = 1:Kb2
          Yp0(:,:,i) = Yp0(:,:,i) + cat_vol_ctype((maxind == k1) .* (maxi~=0) * k1ind(k1) .* Yb(:,:,i)); 
        end
    end
    if ~debug, clear maxi maxind Kb k1 cls2; end
  else
    % more direct method ... a little bit more WM, less CSF
    Yp0 = uint8(max(Yb,min(3,round(Ymi*3)))); Yp0(~Yb) = 0;
  end  


  % use index to speed up and save memory
  sz = size(Yb);
  [indx, indy, indz] = ind2sub(sz,find(Yb>0));
  if job.extopts.AMAPframing
    bx   = (framing.tissue + framing.pve) * 3 * job.extopts.AMAPframing + 2;
    indx = [min(indx) max(indx)] + [-bx bx]; indy = [min(indy) max(indy)] + [-bx bx]; indz = [min(indz) max(indz)] + [-bx bx];
  end
  indx = max((min(indx) - 1),1):min((max(indx) + 1),sz(1));
  indy = max((min(indy) - 1),1):min((max(indy) + 1),sz(2));
  indz = max((min(indz) - 1),1):min((max(indz) + 1),sz(3));

  % Yb source image because Amap needs a skull stripped image
  % set Yp0b and source inside outside Yb to 0
  Yp0b = Yp0(indx,indy,indz);
  Ymib = Ymib(indx,indy,indz); 
  
  % remove non-brain tissue with a smooth mask and set values inside the
  % brain at least to CSF to avoid wholes for images with CSF==BG.
  if job.extopts.LASstr>0 
    Ywmstd = cat_vol_localstat(single(Ymib),Yp0b==3,1,4); 
    CSFnoise(1) = cat_stat_nanmean(Ywmstd(Ywmstd(:)>0))/mean(vx_vol); 
    Ywmstd = cat_vol_localstat(cat_vol_resize(single(Ymib),'reduceV',vx_vol,vx_vol*2,16,'meanm'),...
      cat_vol_resize(Yp0==3,'reduceV',vx_vol,vx_vol*2,16,'meanm')>0.5,1,4); 
    CSFnoise(2) = cat_stat_nanmean(Ywmstd(Ywmstd(:)>0))/mean(vx_vol); 
    Ycsf = double(0.33 * Yb(indx,indy,indz)); spm_smooth(Ycsf,Ycsf,0.6*vx_vol);
    Ycsf = Ycsf + cat_vol_smooth3X(randn(size(Ycsf)),0.5) * max(0.005,min(0.2,CSFnoise(1)/4)); % high-frequency noise
    Ycsf = Ycsf + cat_vol_smooth3X(randn(size(Ycsf)),1.0) * max(0.005,min(0.2,CSFnoise(2)*1)); % high-frequency noise
    Ymib = max(Ycsf*0.8 .* cat_vol_smooth3X(Ycsf>0,2),Ymib); 
    clear Ycsf; 
    % Yb is needed for surface reconstruction
  %     if ~job.output.surface, clear Yb; end
  end
  
  % adaptive mrf noise 
  if job.extopts.mrf>=1 || job.extopts.mrf<0 
    % estimate noise
    [Yw,Yg] = cat_vol_resize({Ymi.*(Ycls{1}>240),Ymi.*(Ycls{2}>240)},'reduceV',vx_vol,3,32,'meanm');
    Yn = max(cat(4,cat_vol_localstat(Yw,Yw>0,2,4),cat_vol_localstat(Yg,Yg>0,2,4)),[],4);
    job.extopts.mrf = double(min(0.15,3*cat_stat_nanmean(Yn(Yn(:)>0)))) * 0.5; 
    clear Yn Yg
  end

  % display something
  stime = cat_io_cmd(sprintf('Amap using initial SPM12 segmentations (MRF filter strength %0.2f)',job.extopts.mrf));       

  
  % intensity values
  Ymib = abs(double(Ymib)); 
  
  if job.extopts.AMAPframing
    Ybb  = Yb0(indx,indy,indz)>0; 
    Ybb  = cat_vol_morph(Ybb,'d',3); 
    pn   = max(double(Yp0b(:))) * 0.015; 
    pnx  = max(double(Yp0b(:))) * 0.015; 
    Yn   = randn(size(Yp0b)); 
    ex   = framing.tissue; 
    ep   = framing.pve; 
    BBe  = ~Ybb; BBe (2      :end-1         , 2      :end-1         , 3   :end-2        ) = false; 
    BBww = ~Ybb; BBww(ex*1   :end+1-ex*1    , ex*1   :end+1-ex*1    , ex*1:end+1-ex*1   ) = false; 
    BBgw = ~Ybb; BBgw(ex*1+0 :end+1-ex*1-0  , ex*1+ep:end+1-ex*1-ep , ex*1:end+1-ex*1-ep) = false; BBgw(BBww) = false; 
    BBgg = ~Ybb; BBgg(ex*2   :end+1-ex*2    , ex*2   :end+1-ex*2    , ex*2:end+1-ex*2   ) = false; BBgg(BBww | BBgw) = false; 
    BBgc = ~Ybb; BBgc(ex*2+0 :end+1-ex*2-0  , ex*2+ep:end+1-ex*2-ep , ex*2:end+1-ex*2-ep) = false; BBgc(BBww | BBgw | BBgg) = false; 
    BBcc = ~Ybb; BBcc(ex*3   :end+1-ex*3    , ex*3   :end+1-ex*3    , ex*3:end+1-ex*3   ) = false; BBcc(BBww | BBgw | BBgg | BBgc) = false; 
    BBcb = ~Ybb; BBcb(ex*3+0 :end+1-ex*3-0  , ex*3+ep:end+1-ex*3-ep , ex*3:end+1-ex*3-ep) = false; BBcb(BBww | BBgw | BBgg | BBgc | BBcc) = false; 
    % extra values + extra noise
    % all peaks have an offset of 0.05 that produces better results
    Ymib(BBcb) = 0.55/3 + pnx * Yn(BBcb); Yp0b(BBcb) = 0; 
    Ymib(BBcc) = 1.05/3 + pnx * Yn(BBcc); Yp0b(BBcc) = 1; 
    Ymib(BBgc) = 1.55/3 + pnx * Yn(BBgc); Yp0b(BBgc) = 1; 
    Ymib(BBgg) = 2.05/3 + pnx * Yn(BBgg); Yp0b(BBgg) = 2; 
    Ymib(BBgw) = 2.75/3 + pnx * Yn(BBgw); Yp0b(BBgw) = 2; 
    Ymib(BBww) = 3.05/3 + pnx * Yn(BBww); Yp0b(BBww) = 3; 
    Ymib(BBe ) = 3.55/3 + pnx * Yn(BBe ); Yp0b(BBe ) = 0; 
    clear BBww BBgw BBgg BBgc BBcc BBcb BBe; 
    
    % add noise 
    addnoise = 0; 
    if addnoise == 2
      % we add noise only in save regions
      WMe = cat_vol_morph( Yp0b==3 , 'e' , 1 )>0; 
      CMe = cat_vol_morph( Yp0b==1 , 'e' , 1 )>0; 
      Ymib( WMe ) = Ymib( WMe ) + pn * Yn( WMe ); 
      Ymib( CMe ) = Ymib( CMe ) + pn * Yn( CMe ); 
    elseif addnoise == 1
      % add noise
      Ymib = Ymib + pn * Yn;
    end
    Ymib = abs(Ymib); 
  end
  
  % Amap parameters  
  % - sub       .. size of sub-elementes is linked to the anatomy and needs adaptation for voxel size 
  %                in additition, test showed that 32 is quite optimal, whereas higher values >64 are worse 
  % - n_iters   .. for highly optimized data is about 10 iterations
  % - bias_fwhm .. the bias correction should be inactive 
  n_iters = 10; sub = round(64/mean(vx_vol));   
  n_classes = 3;  pve = 5; bias_fwhm = 0; init_kmeans = 0;           
  if job.extopts.mrf~=0, iters_icm = 50; else, iters_icm = 0; end    

  % remove noisy background for kmeans
  if init_kmeans, Ymib(Ymib<0.1) = 0; end %#ok<UNRCH>
  
  % do segmentation  
  [prob,amap_means,amap_stds] = cat_amap(Ymib, Yp0b, n_classes, n_iters, sub, pve, init_kmeans,  ...
    job.extopts.mrf, vx_vol, iters_icm, bias_fwhm, 0); 
  fprintf('%5.0fs\n',etime(clock,stime));
  if sum(prob(:)) == 0, error('cat_main:amap','AMAP output empty. '); end
  
  % analyse segmentation ... the input Ym is normalized an the tissue peaks should be around [1/3 2/3 3/3]
  th = {[amap_means(1) amap_stds(1)],[amap_means(2) amap_stds(2)],[amap_means(3) amap_stds(3)]}; 
  
  if job.extopts.AMAPframing
    for i=1:3, prob(:,:,:,i) = prob(:,:,:,i) .* uint8(Ybb); end
  end
  clear Ybb;
  
  if job.extopts.verb>1 
    if strcmpi(spm_check_version,'octave'), pm = '+/-'; else, pm = char(177); end
    fprintf('    AMAP peaks: [CSF,GM,WM] = [%0.2f%s%0.2f,%0.2f%s%0.2f,%0.2f%s%0.2f]\n',...
      th{1}(1),pm,th{1}(2),th{2}(1),pm,th{2}(2),th{3}(1),pm,th{3}(2));
  end
  if th{1}(1)<0 || th{1}(1)>0.6 || th{2}(1)<0.5 || th{2}(1)>0.9 || th{3}(1)<0.95-th{3}(2) || th{3}(1)>1.1 
    error('cat_main:amap',['AMAP estimated untypical tissue peaks that point to an \n' ...
                           'error in the preprocessing before the AMAP segmentation. ']);
  end
  % reorder probability maps according to spm order
  clear Yp0b Ymib; 
  prob = prob(:,:,:,[2 3 1]);  
  clear vol Ymib

  % finally use brainmask before cleanup that was derived from SPM12 segmentations and additionally include
  % areas where GM from Amap > GM from SPM12. This will result in a brainmask where GM areas
  % hopefully are all included and not cut 
  if job.extopts.gcutstr>0 && ~job.inv_weighting
    Yb0(indx,indy,indz) = Yb0(indx,indy,indz) | ((prob(:,:,:,1) > 0) & Yb(indx,indy,indz)); % & ~Ycls{1}(indx,indy,indz));
    for i=1:3
      prob(:,:,:,i) = prob(:,:,:,i) .* uint8(Yb0(indx,indy,indz));
    end
  end
  
  global cat_err_res
  
  % update segmentation for error report
  Yp0  = single(prob(:,:,:,3))/255/3 + single(prob(:,:,:,1))/255*2/3 + single(prob(:,:,:,2))/255;
  [cat_err_res.init.Yp0,cat_err_res.init.BB] = cat_vol_resize(Yp0,'reduceBrain',vx_vol,2,Yp0>0.5); 
  cat_err_res.init.Yp0 = cat_vol_ctype(cat_err_res.init.Yp0/3*255);

end
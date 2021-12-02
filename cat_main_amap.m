function [prob,indx,indy,indz,th,Yrep] = cat_main_amap(Ymi,Yb,Yb0,Ycls,job,res)
% ______________________________________________________________________
%
% AMAP segmentation:
% Most corrections were done before and the AMAP routine is used with 
% a low level of iterations and no further bias correction, because
% some images get tile artifacts. 
%
% [prob,indx,indy,indz,th] = cat_main_amap(Ymi,Yb,Yb0,Ycls,job,res)
%
% prob .. new AMAP segmenation (4D)
% ind* .. index elements to asign a subvolume
% Ymi  .. local intensity normalized source image
% Yb   .. brain mask
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

  global cat_err_res

  % this function adds noise to the data to stabilize processing and we
  % have to define a specific random pattern to get the same results each time
  if exist('rng','file') == 2, rng('default'); rng(0); else, rand('state',0); randn('state',0); end


  % if there is a breakpoint in this file set debug=1 and do not clear temporary variables 
  dbs  = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

  try 
    if job.extopts.ignoreErrors > 3
      error('cat_main_amap:useSPM','Use SPM segmentation.'); 
    end
    
    
    % correct for harder brain mask to avoid meninges in the segmentation
    % RD20200619: unclear severe MATLAB crashes calling cat_amap in the ignoreErrors pipeline
    Ymib   = real(Ymi); Ymib(~Yb) = 0; 
    Ymib(isnan(Ymib) | isinf(Ymib)) = 0;  
    Ymib   = max(0,min(Ymib,2)); 
    rf     = 10^4; Ymib = round(Ymib*rf)/rf;
    d      = size(Ymi); 
    vx_vol = sqrt(sum(res.image(1).mat(1:3,1:3).^2));

    framing.tissue = 4; 
    framing.pve    = 1; 
  
    %  prepare data for segmentation
    if 1 
      %% classic approach, consider the WMH!
      k1ind = [1 2 3 1 0 0 1 0];  
      Kb2   = min(numel(Ycls),numel(k1ind));
      cls2  = zeros([d(1:2) Kb2]);
      Yp0   = zeros(d,'uint8');
      for i=1:d(3)
          for k1 = 1:Kb2, cls2(:,:,k1) = Ycls{k1}(:,:,i); end
          % find maximum for reordered segmentations
          [maxi,maxind] = max(cls2(:,:,[3,1,2,4:Kb2]),[],3);
          for k1 = 1:Kb2
            Yp0(:,:,i) = Yp0(:,:,i) + cat_vol_ctype((maxind == k1) .* (maxi~=0) * k1ind(k1) .* Yb(:,:,i)); 
          end
      end
      if ~debug, clear maxi maxind Kb k1 cls2; end
  
      
      %% correct missing parts by using the intensity normalized map or the old Yp0 label map
% ########## NEW ###########
% RD202008: this may not working correctly 
defineMissingParts = 0;
if defineMissingParts
      if job.extopts.ignoreErrors < 3
        Yp0o = min(3,max(1,uint8(max(Yb,min(3,round(Ymi*3)))))); 
      else
        Yp0o = min(3,max(1,single(Ycls{3})/255*3 + single(Ycls{1})/255*2 + single(Ycls{2})/255));
      end
      Yp0(Yb & Yp0==0) = cat_vol_ctype( round( Yp0o(Yb & Yp0==0) ) ); 
end
% ########## NEW ###########

      if ~debug, clear maxi maxind Kb k1 cls2 Yp0o; end
      
    else
      % more direct method ... a little bit more WM, less CSF
      Yp0 = uint8(max(Yb,min(3,round(Ymi*3)))); Yp0(~Yb) = 0;
    end 
% ########## NEW ###########    
% RD202008: class 0 is required and ignored by the AMAP
if defineMissingParts
    Yp0o = min(3,max(0,uint8(max(Yb,min(3,round(Ymi*3)))))); % class 0 is required!
    Yp0(Yb & (Yp0<=0 | Yp0>3)) = min(3,max(1,cat_vol_ctype( round( Yp0o(Yb & (Yp0<=0 | Yp0>3)) ) ))); clear Yp0o
    if sum( Yp0(Yb(:)>0.5)==0 ) ~= 0 
      error('cat_main_amap:badYp0def', ...
        'Undefined tissues within the brain area may cause severe MATLAB errors while calling cat_amap.c.');
    end
end
% ########## NEW ###########


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
    if job.extopts.LASstr>0 && job.extopts.ignoreErrors < 3 && ...
      ~isfield(job.extopts,'inv_weighting') && ~job.extopts.inv_weighting
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
    end


    % adaptive mrf noise 
    if (job.extopts.mrf>=1 || job.extopts.mrf<0) %&& job.extopts.ignoreErrors < 3
      % estimate noise
      [Yw,Yg] = cat_vol_resize({Ymi.*(Ycls{1}>240),Ymi.*(Ycls{2}>240)},'reduceV',vx_vol,3,32,'meanm');
      Yn = max(cat(4,cat_vol_localstat(Yw,Yw>0,2,4),cat_vol_localstat(Yg,Yg>0,2,4)),[],4);
      job.extopts.mrf = double(min(0.15,3*cat_stat_nanmean(Yn(Yn(:)>0)))) * 0.5; 
      clear Yn Yg
    end

    % display something
    stime = cat_io_cmd(sprintf('Amap using initial SPM12 segmentations (MRF filter strength %0.2f)',job.extopts.mrf));       

    %% intensity values
    Ymib = abs(double(Ymib)); 
  
    if job.extopts.AMAPframing
      Ymib = double(Ymi(indx,indy,indz)) .* double(Yb(indx,indy,indz));  
      Yp0b = uint8(Yp0(indx,indy,indz)) .* uint8(Yb(indx,indy,indz));
      %%
      if ~job.extopts.inv_weighting
        tvals = [0 1 2 3];
      else 
       % [x,xy] = max( res.mg(:) .* (res.lkp(:) == 3) ); 
       % [x,xi] = sort([ mean(res.mn(res.lkp(:) == 1)) mean(res.mn(res.lkp(:) == 2)) mean(res.mn( xy )) ]);
        tvals = [0 cat_stat_kmeans(Ymib(Yp0b==1)) cat_stat_kmeans(Ymib(Yp0b==2)) cat_stat_kmeans(Ymib(Yp0b==3)) ];
      end
     %% 
      Ybb  = Yb0(indx,indy,indz)>0; 
      Ybb  = cat_vol_morph(Ybb,'d',3); 
      pn   = max(double(Yp0b(:))) * 0.015; 
      pnx  = max(double(Yp0b(:))) * 0.015; 
      Yn   = randn(size(Yp0b)); 
      ex   = framing.tissue; 
      ep   = framing.pve; 
      BBe  = ~Ybb; BBe (2      :end-1         , 2      :end-1         , 3      :end-2        ) = false; 
      BBww = ~Ybb; BBww(ex*1   :end+1-ex*1    , ex*1   :end+1-ex*1    , ex*1   :end+1-ex*1   ) = false; 
      BBgw = ~Ybb; BBgw(ex*1+0 :end+1-ex*1-0  , ex*1+ep:end+1-ex*1-ep , ex*1+ep:end+1-ex*1-ep) = false; BBgw(BBww) = false; 
      BBgg = ~Ybb; BBgg(ex*2   :end+1-ex*2    , ex*2   :end+1-ex*2    , ex*2   :end+1-ex*2   ) = false; BBgg(BBww | BBgw) = false; 
      BBgc = ~Ybb; BBgc(ex*2+0 :end+1-ex*2-0  , ex*2+ep:end+1-ex*2-ep , ex*2+ep:end+1-ex*2-ep) = false; BBgc(BBww | BBgw | BBgg) = false; 
      BBcc = ~Ybb; BBcc(ex*3   :end+1-ex*3    , ex*3   :end+1-ex*3    , ex*3   :end+1-ex*3   ) = false; BBcc(BBww | BBgw | BBgg | BBgc) = false; 
      BBcb = ~Ybb; BBcb(ex*3+0 :end+1-ex*3-0  , ex*3+ep:end+1-ex*3-ep , ex*3+ep:end+1-ex*3-ep) = false; BBcb(BBww | BBgw | BBgg | BBgc | BBcc) = false; 
      % extra values + extra noise
      % all peaks have an offset of 0.05 that produces better results
      if ~job.extopts.inv_weighting
        Ymib(BBcb) = 0.55/3 + pnx * Yn(BBcb); Yp0b(BBcb) = 0;  
        Ymib(BBcc) = 1.05/3 + pnx * Yn(BBcc); Yp0b(BBcc) = 1; 
        Ymib(BBgc) = 1.55/3 + pnx * Yn(BBgc); Yp0b(BBgc) = 2; 
        Ymib(BBgg) = 2.05/3 + pnx * Yn(BBgg); Yp0b(BBgg) = 2;
        Ymib(BBgw) = 2.75/3 + pnx * Yn(BBgw); Yp0b(BBgw) = 3;
        Ymib(BBww) = 3.05/3 + pnx * Yn(BBww); Yp0b(BBww) = 3;
        Ymib(BBe ) = 3.55/3 + pnx * Yn(BBe ); Yp0b(BBe ) = 0;
      else
        Ymib(BBcb) = tvals(1)/2        + pnx * Yn(BBcb); Yp0b(BBcb) = 0;  
        Ymib(BBcc) = tvals(2) + 0.05/3 + pnx * Yn(BBcc); Yp0b(BBcc) = 1; 
        Ymib(BBgc) = mean(tvals(2:3))  + pnx * Yn(BBgc); Yp0b(BBgc) = 2; 
        Ymib(BBgg) = tvals(3) + 0.05/3 + pnx * Yn(BBgg); Yp0b(BBgg) = 2;
        Ymib(BBgw) = mean(tvals(3:4))  + pnx * Yn(BBgw); Yp0b(BBgw) = 3;
        Ymib(BBww) = tvals(4) + 0.05/3 + pnx * Yn(BBww); Yp0b(BBww) = 3;
        %Ymib(BBe ) = tvals(1) + 0.55/3 + pnx * Yn(BBe ); Yp0b(BBe ) = 0;
      end
      clear BBww BBgw BBgg BBgc BBcc BBcb BBe; 

      %% add noise 
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
    
    % Amap parameters  - default sub=16 caused errors with highres data!
    % don't use bias_fwhm, because the Amap bias correction is not that efficient and also changes
    % intensity values 
    % RD202006 the bias_fwhm paraemter (and/or other) cause also MATLAB crashes in the ignoreError pipeline 
    % RD202104 sub has to be large because the AMAP bias corrections seams to be in-optimal and caused bad threshold 
    % RD202104 Ymib should only include positive values
    n_iters = 10; sub = round(64/mean(vx_vol));   %#ok<NASGU>
    n_classes = 3;  pve = 5; bias_fwhm = 0; init_kmeans = 0;           %#ok<NASGU>
    if job.extopts.mrf~=0, iters_icm = 50; else iters_icm = 0; end    %#ok<NASGU>
    if 0 %job.extopts.ignoreErrors > 2 || job.extopts.inv_weighting
      % init_kmeans = 0; % k-means was not stable working (e.g. HR075T2) 
      %                  % and it is better to use also here the previous intensity scaling  
      %
      % clsint = @(x) round( sum(res.mn(res.lkp==x) .* res.mg(res.lkp==x)') * 10^5)/10^5; % SPM peak definition  
      % pve = 6 - (clsint(1)<clsint(2)); % seems to be more stable in HR075T2 ...
      %                                  % but I will try the pve5 again after Christians comments  
      %
      % we need more iterations and also some bias correction here because LAS is missing
      n_iters = 200; %#ok<NASGU>
      % bias_fwhm = 60; % using bias_fwhm=60 caused MATLAB errors so we only use more iterations 
    end


    % remove noisy background for kmeans
    if init_kmeans && job.extopts.ignoreErrors < 3, Ymib(Ymib<0.1) = 0; end %#ok<NASGU>

    %% do segmentation  
    amapres = evalc(['[prob,mean] = cat_amap(Ymib, Yp0b, n_classes, n_iters, sub, pve, init_kmeans, ' ...
      'job.extopts.mrf, vx_vol, iters_icm, bias_fwhm);']);
    fprintf('%5.0fs\n',etime(clock,stime));

    %% analyse segmentation ... the input Ym is normalized an the tissue peaks should be around [1/3 2/3 3/3]
    if isempty(amapres)
      % in Ocatave evalc is not working properly and we have to skip the evaluation  
      th = {{nan nan},{nan nan},{nan nan}}; 
    else
      amapres = textscan(amapres,'%s'); amapres = amapres{1}; 
      th{1}   = cell2mat(textscan(amapres{11},'%f*%f')); 
      th{2}   = cell2mat(textscan(amapres{12},'%f*%f')); 
      th{3}   = cell2mat(textscan(amapres{13},'%f*%f')); 
    end

    if job.extopts.AMAPframing
      for i=1:3, prob(:,:,:,i) = prob(:,:,:,i) .* uint8(Ybb); end
    end
    clear Ybb;
    

    if job.extopts.verb>1 && ~isempty(amapres)
      if (th{1}(1) < th{2}(1)) && (th{2}(1) < th{3}(1)) % T1 
        fprintf('    AMAP peaks: [CSF,GM,WM] = [%0.2f%s%0.2f,%0.2f%s%0.2f,%0.2f%s%0.2f]\n',...
          th{1}(1),char(177),th{1}(2),th{2}(1),char(177),th{2}(2),th{3}(1),char(177),th{3}(2));
      else
        fprintf('    AMAP peaks: [%0.2f%s%0.2f,%0.2f%s%0.2f,%0.2f%s%0.2f]\n',...
          th{1}(1),char(177),th{1}(2),th{2}(1),char(177),th{2}(2),th{3}(1),char(177),th{3}(2));
      end
    end
    % if one of the peaks is NaN than create an error
    if ~isempty(amapres) && any( isnan( cell2mat(th) ) )
      error('cat_main_amap:nan',['AMAP estimated NaN tissue peaks that point to an error in the \\\\n' ...
                                 'preprocessing before the AMAP segmentation or inadequate input. ']);
    end
    % fine evaluation of tissue peaks
    if ~isempty(amapres) && (th{1}(1) < th{2}(1)) && (th{2}(1) < th{3}(1)) % T1 
      % in the T1 contrast we have a clear expectation for each tissue class 
      if th{1}(1)<0 || th{1}(1)>0.6 || th{2}(1)<0.5 || th{2}(1)>0.9 || th{3}(1)<0.95-th{3}(2) || th{3}(1)>1.1
        error('cat_main_amap:peaks',['AMAP estimated untypical tissue peaks that point to an \\\\n' ...
                                     'error in the preprocessing before the AMAP segmentation. ']);
      end
    elseif ~isempty(amapres)
      % if there is a very low contrast between two peaks then create an error  
      % because the intensity normalization unsed before was probably incorrect 
      con = [ abs(diff([th{1}(1) th{2}(1)])) , abs(diff([th{1}(1) th{3}(1)])) , abs(diff([th{2}(1) th{3}(1)]))];
      if any( con < 0.15 )
        %% We only create an servere warning here but go one with processing 
        cat_io_addwarning([mfilename ':lowTissueContrast'],sprintf( ...
          ['AMAP estimated quite low tissue contrast that point to problems \\\\n' ...
           'in the preprocessing before the AMAP segmentation or inadequate input \\\\n' ...
           '  [con(c1,c2),con(c1,c3),con(c2,c3)] = [%0.2f,%0.2f,%0.2f] '],con),1 + (min(con) < 0.1),[0 1]);
      end
    end
    
    if isfield(job.extopts,'inv_weighting') && job.extopts.inv_weighting % job.extopts.ignoreErrors > 1 &&
      % RD202006: catching of problems in low quality data - in development 
      probs = prob; 
      ap = [3 1 2]; if numel(Ycls)==7, ap(4)=7; end
      for i=1:numel(ap), probs(:,:,:,i) = Ycls{ap(i)}(indx,indy,indz); end
      % WMHC
      probs(:,:,:,2) = sum( probs(:,:,:,2:2:end) ,4); 
      
      Yp0s = (single(probs(:,:,:,1)) + single(probs(:,:,:,2))*2 + single( probs(:,:,:,3))*3)/255/3; % SPM  segmentation 
      Yp0  = (single(prob(:,:,:,1))  + single(prob(:,:,:,2))*2  + single(prob(:,:,:,3))*3)/255/3;     % AMAP segmentation
      if ( sum(abs(Yp0s(:) - Yp0(:))>0.4) / sum(Yp0(:)>0.5) ) > 0.02
        Yrep = min(1, abs(Yp0s - Yp0) * 3); % * 1 = low correction (less SPM) uint8( abs(Yp0s - Yp0) > 0.3  &  Yp0s>0.1); 
        prob2 = prob; 
        for i=1:3, prob2(:,:,:,i) = cat_vol_ctype( single(prob(:,:,:,i)) .* (1-Yrep) + Yrep .* single(probs(:,:,:,i)) ); end
      
        %rep = mean(Yrep(Yp0(:)))*100; 
        %cat_io_addwarning('cat_main_amap:mixSPMAMAP',sprintf( 'Mix SPM and AMAP segmentation. Use SPM in case of strong differences (%0.2f%%).',rep),1,[0 1]) 
      
        %Yp0c = (single(prob2(:,:,:,1)) + single(prob2(:,:,:,2))*2 + single(prob2(:,:,:,3))*3)/255/3;
      else
        Yrep = false(size(Ymib)); 
      end
      %% separation just for tests/debugging
      if ( sum(abs(Yp0s(:) - Yp0(:))>0.4) / sum(Yp0(:)>0.5) ) > 0.02
        prob = prob2; 
      end
      clear Yp0 Yp0s Yp0c probs prob2
    else
      Yrep = false(size(Ymib)); 
    end

    % reorder probability maps according to spm order
    clear Yp0b Ymib; 
    prob = prob(:,:,:,[2 3 1]);  
    clear vol Ymib

    % finally use brainmask before cleanup that was derived from SPM12 segmentations and additionally include
    % areas where GM from Amap > GM from SPM12. This will result in a brainmask where GM areas
    % hopefully are all included and not cut 
    if job.extopts.gcutstr>0 && ~isfield(job.extopts,'inv_weighting') && ~job.extopts.inv_weighting
      Yb0(indx,indy,indz) = Yb0(indx,indy,indz) | ((prob(:,:,:,1) > 0) & Yb(indx,indy,indz)); % & ~Ycls{1}(indx,indy,indz));
      for i=1:3
        prob(:,:,:,i) = prob(:,:,:,i).*uint8(Yb0(indx,indy,indz));
      end
    end

    
  catch e
    if job.extopts.ignoreErrors < 2
      % update segmentation for error report
      if exist('prob','var') 
        Yp0  = single(prob(:,:,:,3))/255/3 + single(prob(:,:,:,1))/255*2/3 + single(prob(:,:,:,2))/255;
      else
        Yp0  = single(Ycls{3})/255/3 + single(Ycls{1})/255*2/3 + single(Ycls{2})/255;
      end  
      vx_vol = sqrt(sum(res.image(1).mat(1:3,1:3).^2));
      [cat_err_res.init.Yp0,cat_err_res.init.BB] = cat_vol_resize(Yp0,'reduceBrain',vx_vol,2,Yp0>0.5); 
      cat_err_res.init.Yp0 = cat_vol_ctype(cat_err_res.init.Yp0/3*255);
    
      rethrow(e)
    end
    
    % use SPM
    if job.extopts.ignoreErrors < 3
      cat_io_addwarning(e.identifier,e.message,2,[0 1]);  
    end
    
    prob = zeros([size(Ymi),3],'uint8');
    for i = 1:3, prob(:,:,:,i) = Ycls{i}; end
    sz = size(Yb);
    [indx, indy, indz] = ind2sub(sz,find(Yb>0));
    indx = max((min(indx) - 1),1):min((max(indx) + 1),sz(1));
    indy = max((min(indy) - 1),1):min((max(indy) + 1),sz(2));
    indz = max((min(indz) - 1),1):min((max(indz) + 1),sz(3));
    prob = prob(indx,indy,indz,:);
    
    Yrep = false(size(Ymib)); 
  end
end
function [Ysrc,Ycls,Yb,Yb0,Yy,job,res,trans,T3th,stime2] = cat_main_updateSPM(Ysrc,P,Yy,tpm,job,res,stime,stime2)
% ______________________________________________________________________
%  Update SPM preprocessing. 
%  Subfunction of cat_main.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  global cat_err_res; 
  
  res.AffineSPM = res.Affine;
  
  clsint = @(x) round( sum(res.mn(res.lkp==x) .* res.mg(res.lkp==x)') * 10^5)/10^5;

  [pth,nam] = spm_fileparts(res.image0(1).fname); %#ok<ASGLU> % original   

  % voxel size parameter
  vx_vol  = sqrt(sum(res.image(1).mat(1:3,1:3).^2));    % voxel size of the processed image
  vx_volp = prod(vx_vol)/1000;

  %d = res.image(1).dim(1:3);

  % some reports
  for i=1:size(P,4), res.ppe.SPMvols0(i) = cat_stat_nansum(single(P(:,:,:,i)),0)/255 .* prod(vx_vol) / 1000; end

   
  try
    if job.extopts.ignoreErrors > 2 % && ~( ( clsint(3) < clsint(1) ) &&  ( clsint(1) < clsint(2) ) ) % ~T1
      error('cat_main_updateSPM:runbackup','Test backup function.');
    end
    
    stime2 = cat_io_cmd('  Update segmentation','g5','',job.extopts.verb-1,stime2);
  
    % Create brain mask based on the the TPM classes
    % cleanup with brain mask - required for ngaus [1 1 2 4 3 2] and R1/MP2Rage like data 
   % YbA = zeros(d,'single');
    Vb = tpm.V(1); Vb.pinfo(3) = 0; Vb.dt=16; 
    Vb.dat = single(exp(tpm.dat{1}) + exp(tpm.dat{2}) + exp(tpm.dat{3})); 
    YbA = cat_vol_sample(res.tpm(1),Vb,Yy,1);
   %for z=1:d(3)
   %   YbA(:,:,z) = spm_sample_vol(Vb,double(Yy(:,:,z,1)),double(Yy(:,:,z,2)),double(Yy(:,:,z,3)),1); 
   % end
    if (isfield(job,'useprior') && ~isempty(job.useprior) ), bth = 0.5; else, bth = 0.1; end
    if round(max(YbA(:))/Vb.pinfo(1)), YbA=YbA>bth*Vb.pinfo(1); else, YbA=YbA>bth; end
    % add some distance around brainmask (important for bias!)
    YbA = YbA | cat_vol_morph(YbA & sum(P(:,:,:,1:2),4)>4 ,'dd',2.4,vx_vol);


    if size(P,4)==3
      %% Skull-stripping of post mortem data (RD202108). 
      %  Here we only have 3 classes. GM, WM and a combined CSF/background 
      %  class and we use the GM/WM block to create a rough brain mask that
      %  is used to artificially seperate between "CSF" and background to  
      %  obtain some useful CSF volume values and avoid preprocessing 
      %  problems by other structures. 
      Yb = smooth3( sum(P(:,:,:,1:2),4) ) > 64;        % initial mask 
      Yb = cat_vol_morph( Yb , 'ldo' ,   2 , vx_vol );  % (light) opening to remove menignes and (large) unconnected components
      Yb = cat_vol_morph( Yb , 'dc'  ,  10 , vx_vol );  % major closing to include internal CSF
      Yb = cat_vol_morph( Yb , 'ldo' ,   5 , vx_vol );  % final opening to remove further menignes
      Yb = cat_vol_morph( Yb , 'dd'  ,   2 , vx_vol );  % add one mm to avoid to hard skull-strippings that could trouble thickness estimation  
      Yb = cat_vol_ctype(Yb);                           % convert to uint8 like P  
      for i=1:3, P(:,:,:,i) = P(:,:,:,i) .* Yb; end     % apply masking for major tissues ...
      P(:,:,:,4) = cat_vol_ctype((1 - Yb) * 255);       % ... and create a new background class
      Yb = Yb>0.5;                                      % convert to logical 
      postmortem = 1; 
    else 
      postmortem = 0; 
    end

    %% correction of CSF-GM-WM PVE voxels that were miss aligned to skull
    %  RD202006: This is a special cases observed for the thickness phantom. 
    %  With a coronal head-masking, SPM had problems with brain PVE values
    %  classes that were to close to some head peaks and therefore miss 
    %  classified in P(:,:,:,5).
    % GM-WM
    if size(P,4)>4
      pth = 128; 
      lth = min( clsint(1) , clsint(2) ); hth = max( clsint(1) , clsint(2) ); 
      Ycls5b = YbA & cat_vol_morph(P(:,:,:,1)>pth | P(:,:,:,2)>pth,'c') & Ysrc>lth & Ysrc<hth; 
      P5w = single(P(:,:,:,5)) .* Ycls5b .* max(0,min(1,abs(Ysrc - clsint(1)) ./ abs(diff([clsint(1),clsint(2)])))); 
      P5g = single(P(:,:,:,5)) .* Ycls5b  -  P5w; 
      P(:,:,:,1) = P(:,:,:,1) + uint8( P5w ); 
      P(:,:,:,2) = P(:,:,:,2) + uint8( P5g ); 
      P(:,:,:,5) = P(:,:,:,5) - uint8( P5g + P5w ); 
      clear P5g P5w
      % GM-CSF
      lth = min( clsint(1) , clsint(3) ); hth = max( clsint(1) , clsint(3) ); 
      Ycls5b = YbA & cat_vol_morph(P(:,:,:,1)>pth | P(:,:,:,3)>pth,'c') & Ysrc>lth & Ysrc<hth; 
      P5c = single(P(:,:,:,5)) .* Ycls5b .* max(0,min(1,abs(Ysrc - clsint(1)) ./ abs(diff([clsint(1),clsint(3)])))); 
      P5g = single(P(:,:,:,5)) .* Ycls5b  -  P5c; 
      P(:,:,:,1) = P(:,:,:,1) + uint8( P5g ); 
      P(:,:,:,3) = P(:,:,:,3) + uint8( P5c ); 
      P(:,:,:,5) = P(:,:,:,5) - uint8( P5g + P5c ); 
      clear P5g P5c
      % WM-CSF
      lth = min( clsint(2) , clsint(3) ); hth = max( clsint(2) , clsint(3) ); 
      Ycls5b = YbA & cat_vol_morph(P(:,:,:,2)>pth | P(:,:,:,3)>pth,'c') & Ysrc>lth & Ysrc<hth; 
      P5w = single(P(:,:,:,5)) .* Ycls5b .* max(0,min(1,abs(Ysrc - clsint(3)) ./ abs(diff([clsint(1),clsint(3)])))); 
      P5c = single(P(:,:,:,5)) .* Ycls5b  -  P5w; 
      P(:,:,:,2) = P(:,:,:,2) + uint8( P5w ); 
      P(:,:,:,3) = P(:,:,:,3) + uint8( P5c ); 
      P(:,:,:,5) = P(:,:,:,5) - uint8( P5w + P5c );
      clear P5w P5c
      
      if isfield(job.extopts,'inv_weighting') && job.extopts.inv_weighting
        % RD202105:  Cleanup for CSF artefacts between GM and WM (PD/T2 BWP data)
        %            there should be no CSF noise close to the WM 
        Ywd  = cat_vol_morph( P(:,:,:,2)>128 , 'd');
        Ycc  = P(:,:,:,3)>0 & smooth3(single(P(:,:,:,3)/255) < 0.5) & Ywd; 
        P(:,:,:,1) = P(:,:,:,1) + P(:,:,:,3) .* uint8( Ycc & abs(Ysrc-clsint(1))<=abs(Ysrc-clsint(2)) & ...
          abs(Ysrc-clsint(1))<=abs(Ysrc-clsint(3)) & abs(Ysrc-clsint(2))<=abs(Ysrc-clsint(3)) ); 
        P(:,:,:,2) = P(:,:,:,2) + P(:,:,:,3) .* uint8( Ycc & abs(Ysrc-clsint(1))> abs(Ysrc-clsint(2)) & ...
          abs(Ysrc-clsint(1))<=abs(Ysrc-clsint(3)) & abs(Ysrc-clsint(2))<=abs(Ysrc-clsint(3)) );  
        P(:,:,:,3) = P(:,:,:,3) - P(:,:,:,3) .* uint8( Ycc ); 
        clear Ywd Ycc; 
        % RD20215:  Cleanup for CSF artefacts between GM and WM (PD/T2 BWP data)
        %           there should be no head tissue inside the brain that can be explained by brain values
        tth = [clsint(1) clsint(2) clsint(3)]; 
        Ycc = cat_vol_morph(YbA,'e',3) & sum( P(:,:,:,4:6) , 4 ) & ...
          Ysrc>(min(tth) - mean( abs(diff(tth))) ) & Ysrc<(max(tth) + mean( abs(diff(tth)) ));
        [Ytmp,Yct] = min( cat( 4 , abs(Ysrc - tth(1)) ,  abs(Ysrc - tth(2)) ,  abs(Ysrc - tth(3)) ) , [], 4); clear Ytmp
        P(:,:,:,1) = P(:,:,:,1) + P(:,:,:,5) .* uint8( Ycc & Yct==1 );
        P(:,:,:,2) = P(:,:,:,2) + P(:,:,:,5) .* uint8( Ycc & Yct==2 );  
        P(:,:,:,3) = P(:,:,:,3) + P(:,:,:,5) .* uint8( Ycc & Yct==3 );  
        P(:,:,:,5) = P(:,:,:,5) - P(:,:,:,5) .* uint8( Ycc ); 
        clear Ycc tth Yct; 
      end
    end

    % clear WM segment
    Ywm = P(:,:,:,2) > 128; 
    Ybe = YbA & ~cat_vol_morph(YbA,'de',8/mean(vx_vol)); 
    Ywm = cat_vol_morph(Ywm,'l',[100 0.1 ]);
    P(:,:,:,min(size(P,4),5)) = P(:,:,:,min(size(P,4),5)) + P(:,:,:,2) .* uint8( Ybe & ~Ywm ) / 2; 
    P(:,:,:,2) = P(:,:,:,2) - P(:,:,:,2) .* uint8( Ybe & ~Ywm ) / 2;  
    clear Ywm be; 

    %% transfer tissue outside the brain mask to head  ... 
    % RD 201807: I am not sure if this is a good idea. Please test this with children! 
    for i=1:3
      P(:,:,:,4) = cat_vol_ctype(single(P(:,:,:,4)) + single(P(:,:,:,i)) .* single(~YbA)); 
      P(:,:,:,i) = cat_vol_ctype(single(P(:,:,:,i)) .* single(YbA)); 
    end
    

    % RD202006: Correct background (from cat_run_job)
    % RD202007: The noisy zeros background in resliced data (e.g. long avg)
    %           can possibly cause problems? - no 
    if isfield(res,'bge')
      P(:,:,:,end) = max( cat_vol_ctype( res.bge * 255 ) , P(:,:,:,end) ); 
      for i=1:size(P,4)-1, P(:,:,:,i) = P(:,:,:,i) .* cat_vol_ctype(1 - res.bge); end
    end

    %%
    % Cleanup for high resolution data
    % Although the old cleanup is very slow for high resolution data, the   
    % reduction of image resolution removes spatial segmentation information. 
    % RD202008: This operation has to be done in high-resolution and it is 
    %           maybe better to avoid the older gcw cleaning step the is
    %           very very slow. 
    if job.opts.redspmres==0 % already done in case of redspmres
      if max(vx_vol)<1.5 && mean(vx_vol)<1.3
        for i=1:size(P,4), [Pc1(:,:,:,i),BB] = cat_vol_resize(P(:,:,:,i),'reduceBrain',vx_vol,4,YbA); end %#ok<AGROW>
        Pc1 = cat_main_clean_gwc(Pc1,max(1,min(2,job.extopts.cleanupstr*2)) / mean(vx_vol));
        Ybb = ones(size(YbA),'uint8'); Ybb(BB.BB(1):BB.BB(2),BB.BB(3):BB.BB(4),BB.BB(5):BB.BB(6)) = uint8(1); 
        for i=1:size(P,4), P(:,:,:,i) = Ybb.*P(:,:,:,i) + cat_vol_resize(Pc1(:,:,:,i),'dereduceBrain',BB); end 
        clear Pc1 Ybb;
      end
    end
    clear YbA;

    

    %% guarantee probability 
    sP = (sum(single(P),4)+eps)/255;
    for k1=1:size(P,4), P(:,:,:,k1) = cat_vol_ctype(single(P(:,:,:,k1))./sP); end
    clear sP;


    % Use median for WM threshold estimation to avoid problems in case of WMHs!
    WMth = double(max( clsint(2) , cat_stat_nanmedian(Ysrc(P(:,:,:,2)>192)) )); 
    if clsint(3)>clsint(2) % invers
      CMth = clsint(3); 
    else
      CMth = min( [  clsint(1) - diff([clsint(1),WMth]) , clsint(3) ]);
    end
    T3th = double( [ CMth , clsint(1) , WMth]);


    %% Some error handling
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
            'Possibly the affine registration failed. %s.\n' ...
            ' Tissue class:           %10s%10s%10s%10s%10s%10s\n' ...
            ' Rel. to image volume:   %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n' ...
            ' Rel. to brain volume:   %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n' ...
            ' Tissue intensity:       %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f'],...
            spm_file('Please check image orientation and quality','link',['spm_image(''Display'', ''' res.image0.fname ''')']), ...
            'BG','CSF','GM','WM','HDH','HDL', ...
            [ clsvol([6 3 1 2 4 5])/cat_stat_nansum(clsvol)*100, clsvol([6 3 1 2 4 5])/cat_stat_nansum(clsvol(1:3))*100, BGth,T3th,HDHth,HDLth]));  %#ok<SPERR>
      elseif size(P,4)==4 % skull-stripped
          error('CAT:cat_main:SPMpreprocessing:emptySegmentation', ...
           sprintf(['Empty Segmentation: \n ' ...
            'Possibly the affine registration failed. %s.\n' ...
            ' Tissue class:           %10s%10s%10s%10s\n' ...
            ' Rel. to image volume:   %10.2f%10.2f%10.2f%10.2f\n' ...
            ' Rel. to brain volume:   %10.2f%10.2f%10.2f%10.2f\n' ...
            ' Tissue intensity:       %10.2f%10.2f%10.2f%10.2f'],...
            spm_file('Please check image orientation and quality','link',['spm_image(''Display'', ''' res.image0.fname ''')']), ...
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
    clear Yp0; 

    % ### This can not be reached because the mask field is removed by SPM! ###
    if isfield(res,'msk') 
      Ybg = ~res.msk.dat; 
      P4  = cat_vol_ctype( single(P(:,:,:,6)) .* (Ysrc<T3th(2))  .* (Ybg<0.5) + single(P(:,:,:,4)) .* (Ybg<0.5) ); % remove air in head
      P5  = cat_vol_ctype( single(P(:,:,:,6)) .* (Ysrc>=T3th(2)) .* (Ybg<0.5) + single(P(:,:,:,5)) .* (Ybg<0.5) ); % remove air in head
      P6  = cat_vol_ctype( single(sum(P(:,:,:,4:5),4)) .* (Ybg>0.5) + single(P(:,:,:,6)) .* (Ybg>0.5) ); % add objects/artifacts to background
      P(:,:,:,4) = P4;
      P(:,:,:,5) = P5;
      P(:,:,:,6) = P6;
      clear P4 P5 P6 Ybg; 
    end




    %% Skull-Stripping
    %  ----------------------------------------------------------------------
    %  Update Skull-Stripping 1
    %  ----------------------------------------------------------------------
    stime2 = cat_io_cmd('  Update skull-stripping','g5','',job.extopts.verb-1,stime2); 
    if (isfield(job,'useprior') && ~isempty(job.useprior) && strcmp(job.opts.affreg,'prior') ) && ... 
       (isfield(res,'ppe') && ~res.ppe.affreg.highBG)
      % RD202010: use longitudinal skull-stripping 
      [pp,ff,ee] = spm_fileparts(char(job.useprior));
      Pavgp0 = fullfile(pp,'mri',[strrep(ff,'avg_','p0avg_'),ee]);

      % get gradient and divergence map (Yg and Ydiv)
      [Ytmp,Ytmp,Yg,Ydiv] = cat_main_updateSPM_gcut0(Ysrc,P,vx_vol,T3th); clear Ytmp;  %#ok<ASGLU>
      if exist(Pavgp0,'file')
        % the p0avg should be optimal 
        if any(vx_vol0 ~= vx_vol) % if the data was internaly resampled we have to load it via imcalc
          [Vb,Yb] = cat_vol_imcalc(spm_vol(Pavgp0),spm_vol(res.image.fname),'i1',struct('interp',3,'verb',0,'mask',-1)); clear Vb;  %#ok<ASGLU>
        else
          Yb = spm_read_vols(spm_vol(Pavgp0));
        end
        Yb = Yb > 0.5; 
      else
        % otherwise it would be possible to use the individual TPM 
        % however, the TPM is more smoothed and is therefore only second choice  
        cat_io_cprintf('warn','Cannot find p0avg use TPM for brainmask: \n  %s\n',Pavgp0);
        Yb = YbA > 0.5;
        clear YbA
      end
      Ybb = cat_vol_ctype(cat_vol_smooth3X(Yb,0.5)*256); 
      
      %% correct tissues
      %  RD20221224: Only the brainmask wasn't enough and we need to cleanup 
      %              the segmentation also here (only for long pipeline)
      % move brain tissue to head tissues or vice versa
      for ti = 1:3
        if ti == 1 % GM with soft bounary to reduce meninges
          Ynbm = cat_vol_ctype( single(P(:,:,:,ti)) .* (1 - max(0,2 * smooth3(Yb) - 1) ) ); 
          Ybm  = cat_vol_ctype( single(P(:,:,:,5))  .* (    max(0,2 * smooth3(Yb) - 1) ) ); 
        elseif ti == 2 % WM with very soft boundary because we exptect no WM close to the skull
          Ynbm = cat_vol_ctype( single(P(:,:,:,ti)) .* (1 - max(0,2 * single(Ybb)/255 - 1) ) ); 
          Ybm  = cat_vol_ctype( single(P(:,:,:,5))  .* (    max(0,2 * single(Ybb)/255 - 1) ) ); 
        else % CSF with hard boundary
          Ynbm = cat_vol_ctype( single(P(:,:,:,ti)) .* (1 - Yb) ); 
          Ybm  = cat_vol_ctype( single(P(:,:,:,5))  .* (    Yb) ); 
        end
        P(:,:,:,ti) = P(:,:,:,ti) - Ynbm + Ybm; 
        P(:,:,:,5)  = P(:,:,:,5)  + Ynbm - Ybm; 
        clear Ynbm Ybm;
      end
      % some extra GM cleanup for meninges
      Yngm = P(:,:,:,1) .* uint8( Ybb<255 & (P(:,:,:,1)>64) & (smooth3( single(P(:,:,:,1)>64) )<0.5) );
      P(:,:,:,1) = P(:,:,:,1) - Yngm; P(:,:,:,5) = P(:,:,:,5) + Yngm; %if ~debug, clear Yngm; end
      % some further hard GM cleanup ? 
      %{
      Yp0avg = spm_read_vols(spm_vol(Pavgp0));
      Yngm = P(:,:,:,1) .* uint8( cat_vol_morph( Yp0avg < 1.75 , 'de' , 3, vx_vol) & Yp0yvg>0 );
      P(:,:,:,1) = P(:,:,:,1) - Yngm; P(:,:,:,5) = P(:,:,:,5) + Yngm; %if ~debug, clear Yngm; end
      %}
    elseif postmortem
      % already done 
    elseif size(P,4)==4 || size(P,4)==3 % skull-stripped
      [Yb,Ybb,Yg,Ydiv,P] = cat_main_updateSPM_skullstriped(Ysrc,P,res,vx_vol,T3th);
    elseif isfield(job.extopts,'inv_weighting') && job.extopts.inv_weighting
      %% estimate gradient (edge) and divergence maps
      if ~exist('Ym0','var')
        Ym0 = single(P(:,:,:,3))/255 + single(P(:,:,:,1))/255 + single(P(:,:,:,2))/255;
      end
      Ym0 = cat_vol_smooth3X(Ym0,4/mean(vx_vol)); 
      Yb  = Ym0 > min(0.5,max(0.25, job.extopts.gcutstr)); clear Ym0
      Yb  = cat_vol_morph(cat_vol_morph(Yb,'lo'),'c');
      
      %% create brain level set map (from cat_main_APRG)
      Yc   = single(P(:,:,:,3))/255;
      
      %  Ym .. combination of brain tissue and CSF that is further corrected
      %        for noise (median) and smoothness (Laplace) and finally 
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

      cutstr = min(0.5,max(0.25, job.extopts.gcutstr)); 
      Ybb  = cat_vol_ctype( max(0,min(1,(Ym - cutstr)/(1-cutstr))) * 256); 
      
      %%
      [Ysrcb,BB] = cat_vol_resize({Ysrc},'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yb);
      Yg   = cat_vol_grad(Ysrcb/T3th(3),vx_vol);
      Ydiv = cat_vol_div(Ysrcb/T3th(3),vx_vol);
      Yg   = cat_vol_resize(Yg   ,'dereduceBrain',BB);
      Ydiv = cat_vol_resize(Ydiv ,'dereduceBrain',BB);
    elseif job.extopts.gcutstr==0 
      [Yb,Ybb,Yg,Ydiv] = cat_main_updateSPM_gcut0(Ysrc,P,vx_vol,T3th);
    elseif job.extopts.gcutstr==2
      [Yb,Ybb,Yg,Ydiv] = cat_main_APRG(Ysrc,P,res,T3th,0);
    elseif job.extopts.gcutstr>2 && job.extopts.gcutstr<3
      [Yb,Ybb,Yg,Ydiv] = cat_main_APRG(Ysrc,P,res,T3th,job.extopts.gcutstr);
    else
      [Yb,Ybb,Yg,Ydiv] = cat_main_updateSPM_gcutold(Ysrc,P,res,vx_vol,T3th);
    end





    %% save brain mask using SPM12 segmentations for later use
    if ~exist('Ym0','var')
      Ym0 = single(P(:,:,:,3))/255 + single(P(:,:,:,1))/255 + single(P(:,:,:,2))/255;
    end
    Ym0 = cat_vol_smooth3X(Ym0,4/mean(vx_vol)); 
    Yb0 = (Ym0 > min(0.5,max(0.25, job.extopts.gcutstr))); clear Ym0
    Yb0 = cat_vol_morph(cat_vol_morph(Yb0,'ldo',1),'c',3);




    %% RD202110: Background correction in longitidunal mode
    %  We observed some problems in the SPM background segmentation for
    %  longidutidnal processing that detected the volume of the boundary box 
    %  whereas the real background was miss-aligned to class 5 that caused 
    %  further problems in the LAS function that were solved too. Although 
    %  it would be possible to adapt the SPM segmentation, eg. by adapting 
    %  the number of gaussians per class, we decided that it is simpler and 
    %  maybe saver to add further test in the longitudinal case, where the 
    %  TPM should be close to the segmentation outcome.  
    if (isfield(job,'useprior') && ~isempty(job.useprior) ) && ... 
       (isfield(job,'ppe') && ~job.ppe.affreg.highBG)
      % sum of all TPM classes without background
      Vall = tpm.V(end); Vall.pinfo(3) = 0; Vall.dt=16; 
      Vall.dat = zeros(size(tpm.dat{1})); for k1 = 1:numel(tpm.dat)-1, Vall.dat = Vall.dat + single(exp(tpm.dat{k1})); end 
      Yall = cat_vol_sample(res.tpm(1),Vall,Yy,1);

      % backgound class
      Ybg = 1 - Yall; clear Yall Vall; 

      % estimate error and do correction 
      rmse = @(x,y) mean( (x(:) - y(:)).^2 ).^0.5; 
      % the TPM BG may be smaller due to the limited overlap and we need a
      % higher threshold to avoid unnecessary background corrections
      TPisSmaller = ( sum(sum(sum(single(P(:,:,:,end))/255))) - sum(Ybg(:))) < 0;  
      if rmse(Ybg,single(P(:,:,:,end))/255) > 0.3 + 0.2*TPisSmaller % just some threshold (RD20220103: adjusted by TPisSmaller)
        % setup new background
        Ynbg = Ybg>0.5 | P(:,:,:,end)>128;
        Ynbg = cat_vol_morph(Ynbg,'dc',5,vx_vol); 
        Ynbg = uint8( 255 .* smooth3(Ynbg) ); 

        % correct classes
        for k1 = 1:size(P,4)-1, P(:,:,:,k1) = P(:,:,:,k1) - min(P(:,:,:,k1),Ynbg); end
        P(:,:,:,end) = max( Ynbg , P(:,:,:,end) ); 
        clear Ynbg; 

        % normalize all classes
        sP = (sum(single(P),4)+eps)/255;
        for k1=1:size(P,4), P(:,:,:,k1) = cat_vol_ctype(single(P(:,:,:,k1))./sP); end
        clear sP; 

        cat_io_addwarning('cat_main_updateSPM:ReplacedLongBackground','Detected and corrected inadequate background \\nsegmentation in longitudinal mode.',0,[1 2]);
      end
      clear Ybg; 
    end

  
  %%
    stime2 = cat_io_cmd('  Update probability maps','g5','',job.extopts.verb-1,stime2);
    if ~(any(sign(diff(T3th))==-1)) && ...
       ~( (isfield(job,'useprior') && ~isempty(job.useprior) ) && ... % no single longitudinal timepoint
        (isfield(res,'ppe') && ~res.ppe.affreg.highBG) )
      %% Update probability maps
      % background vs. head - important for noisy backgrounds such as in MT weighting
      if size(P,4)==4 || size(P,4)==3 % skull-stripped
        Ybg = ~Yb;
      else
        if sum(sum(sum(P(:,:,:,6)>240 & Ysrc<cat_stat_nanmean(T3th(1:2)))))>10000
          Ybg = P(:,:,:,6); 
          [Ybgr,Ysrcr,resT2] = cat_vol_resize({Ybg,Ysrc},'reduceV',vx_vol,2,32); 
          Ybgrth = max(cat_stat_nanmean(Ysrcr(Ybgr(:)>128)) + 2*std(Ysrcr(Ybgr(:)>128)),T3th(1));
          Ybgr = cat_vol_morph(cat_vol_morph(cat_vol_morph(Ybgr>128,'d') & Ysrcr<Ybgrth,'lo',1),'lc',1);
          Ybg  = cat_vol_resize(cat_vol_smooth3X(Ybgr,1),'dereduceV',resT2); 
          clear Ysrcr Ybgr; 
        elseif sum(sum(sum(P(:,:,:,6)>8 & Ysrc<cat_stat_nanmean(T3th(1:2)))))>10000
          % RD202010 bad SPM background
          Ybg = cat_vol_smooth3X( single(P(:,:,:,6)) * 240 ,2); 
          [Ybgr,Ysrcr,resT2] = cat_vol_resize({Ybg,Ysrc},'reduceV',vx_vol,2,32); 
          Ybgrth = max(cat_stat_nanmean(Ysrcr(Ybgr(:)>128)) + 2*std(Ysrcr(Ybgr(:)>128)),T3th(1));
          Ybgr = cat_vol_morph(cat_vol_morph(cat_vol_morph(Ybgr>128,'d') & Ysrcr<Ybgrth,'lo',1),'lc',1);
          Ybg  = cat_vol_resize(cat_vol_smooth3X(Ybgr,1),'dereduceV',resT2); 
          clear Ysrcr Ybgr; 
        else
          Ybg = P(:,:,:,6); %~Yb;
        end
      end
      if ~res.ppe.affreg.highBG
        P4   = cat_vol_ctype( single(P(:,:,:,6)) .* (Ysrc<T3th(2))  .* (Ybg<0.5) + single(P(:,:,:,4)) .* (Ybg<0.5) ); % remove air in head
        P5   = cat_vol_ctype( single(P(:,:,:,6)) .* (Ysrc>=T3th(2)) .* (Ybg<0.5) + single(P(:,:,:,5)) .* (Ybg<0.5) ); % remove air in head
        P6   = cat_vol_ctype( single(sum(P(:,:,:,4:5),4)) .* (Ybg>0.5) + single(P(:,:,:,6)) .* (Ybg>0.5) ); % add objects/artifacts to background
        P(:,:,:,4) = P4;
        P(:,:,:,5) = P5;
        P(:,:,:,6) = P6;
        clear P4 P5 P6;
      end

      % correct probability maps to 100% 
      sumP = cat_vol_ctype(255 - sum(P(:,:,:,1:6),4));
      P(:,:,:,1) = P(:,:,:,1) + sumP .* uint8( Ybg<0.5  &  Yb & Ysrc>cat_stat_nanmean(T3th(1:2)) & Ysrc<cat_stat_nanmean(T3th(2:3)));
      P(:,:,:,2) = P(:,:,:,2) + sumP .* uint8( Ybg<0.5  &  Yb & Ysrc>=cat_stat_nanmean(T3th(2:3)));
      P(:,:,:,3) = P(:,:,:,3) + sumP .* uint8( Ybg<0.5  &  Yb & Ysrc<=cat_stat_nanmean(T3th(1:2)));
      P(:,:,:,4) = P(:,:,:,4) + sumP .* uint8( Ybg<0.5  & ~Yb & Ysrc<T3th(2));
      P(:,:,:,5) = P(:,:,:,5) + sumP .* uint8( Ybg<0.5  & ~Yb & Ysrc>=T3th(2));
      P(:,:,:,6) = P(:,:,:,6) + sumP .* uint8( Ybg>=0.5 & ~Yb );
      clear Ybg sumP;
      
    
      %% head to WM 
      % Under-correction of strong inhomogeneities in high field scans 
      % (>1.5T) can cause miss-alignments of the template and therefore 
      % miss classifications of the tissues that finally avoid further 
      % corrections in by LAS. 
      % Typically the alignment failed in this cases because the high 
      % intensities next to the head that were counted as head and not
      % corrected by SPM.
      % e.g. HR075, Magdeburg7T, SRS_SRS_Jena_DaRo81_T1_20150320-191509_MPR-08mm-G2-bw330-nbc.nii, ...
      Ywm = single(P(:,:,:,2)>128 & Yg<0.3 & Ydiv<0.03); Ywm(Ybb<128 | (P(:,:,:,1)>128 & abs(Ysrc/T3th(3)-2/3)<1/3) | Ydiv>0.03) = nan;
      [Ywm1,YD] = cat_vol_downcut(Ywm,1-Ysrc/T3th(3),0.02); Yb(isnan(Yb))=0; Ywm(YD<300)=1; Ywm(isnan(Ywm))=0; clear Ywm1 YD; %#ok<ASGLU>
      Ywmc = uint8(smooth3(Ywm)>0.7);
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
      clear Ydiv;
      Ygm = uint8(Ygm); 
      P(:,:,:,5) = P(:,:,:,5) .* (1-Ygm);
      P(:,:,:,3) = P(:,:,:,3) .* (1-Ygm);
      P(:,:,:,2) = P(:,:,:,2) .* (1-Ygm);
      P(:,:,:,1) = cat_vol_ctype(single(P(:,:,:,1)) + 255*single(Ygm));
      Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
      clear Ywm Ygm;

      %% head to CSF
      Ycsf = (Ysrc > T3th(1) -  min( abs( [ diff(T3th(1:2)) diff(T3th(1:2:3)) ] ))/2 ) & ...
             (Ysrc > T3th(1) +  min( abs( [ diff(T3th(1:2)) diff(T3th(1:2:3)) ] ))/2 ) & ...
             Yb & sum(P(:,:,:,4:2:end),4)>0 & sum(P(:,:,:,1:3),4)<250 & ...
             Yg<cat_stat_nanmean(Yg(P(:,:,:,3)>64)*1.5); % & Ydiv<cat_stat_nanmean(Ydiv(P(:,:,:,3)>64));
      Ycsf(smooth3(Ycsf)<0.5)=0;
      for pi=4:2:size(P,4), P(:,:,:,pi) = P(:,:,:,pi) .* cat_vol_ctype(1-Ycsf); end
      P(:,:,:,3) = cat_vol_ctype(single(P(:,:,:,3)) + 255*single(Ycsf));
      clear Yg;
      
      %% remove brain tissues outside the brain mask ...
      %  tissues > skull (within the brain mask)
      Yhdc = uint8(smooth3( Ysrc/T3th(3).*(Ybb>cat_vol_ctype(0.2*255)) - Yp0 )>0.5); 
      sumP = sum(P(:,:,:,1:3),4); 
      P(:,:,:,4)   =  cat_vol_ctype( single(P(:,:,:,4)) + sumP .* ((Ybb<=cat_vol_ctype(0.05*255)) | Yhdc ) .* (Ysrc<T3th(2)));
      P(:,:,:,5)   =  cat_vol_ctype( single(P(:,:,:,5)) + sumP .* ((Ybb<=cat_vol_ctype(0.05*255)) | Yhdc ) .* (Ysrc>=T3th(2)));
      P(:,:,:,1:3) =  P(:,:,:,1:3) .* repmat(uint8(~(Ybb<=cat_vol_ctype(0.05*255)) | Yhdc ),[1,1,1,3]);
      clear sumP Yp0 Yhdc; 
    end
    clear Ybb;

    sP = (sum(single(P),4)+eps)/255;
    for k1=1:size(P,4), P(:,:,:,k1) = cat_vol_ctype(single(P(:,:,:,k1))./sP); end
    



    %% MRF
    % Used spm_mrf help and tested the probability TPM map for Q without good results.         
    nmrf_its = 0; % 10 iterations better to get full probability in thin GM areas 
    spm_progress_bar('init',nmrf_its,['MRF: Working on ' nam],'Iterations completed');
    if isfield(res,'mg'), Kb = max(res.lkp); else, Kb = size(res.intensity(1).lik,2); end
    G   = ones([Kb,1],'single');
    vx2 = single(sum(res.image(1).mat(1:3,1:3).^2));
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
    for iter=1:nmrf_its
        P = spm_mrf(P,single(P),G,vx2); % spm_mrf(P,Q,G,vx2);
        spm_progress_bar('set',iter);
    end

    % update segmentation for error report
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
    [cat_err_res.init.Yp0,cat_err_res.init.BB] = cat_vol_resize(Yp0,'reduceBrain',vx_vol,2,Yp0>0.5); 
    cat_err_res.init.Yp0 = cat_vol_ctype(cat_err_res.init.Yp0/3*255);

    spm_progress_bar('clear');
    Ycls = cell(1,size(P,4)); 
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
    
  catch e
  % just try to translate the input to the output
  %   [Ysrc,Ycls,Yb,Yb0,job,res,T3th,stime2] = cat_main_updateSPM(Ysrc,P,Yy,tpm,job,res,stime,stime2)

    if job.extopts.ignoreErrors < 2
      rethrow(e)
    else 
      % print that that backup function is used 
      % print the error message as warning in the developer mode
      cat_io_addwarning('cat_main_updateSPM:runbackup','IgnoreErrors: Run backup function',1)
      if ~strcmp(e.identifier, 'cat_main_updateSPM:runbackup')
        fprintf('\n'); 
        warning(e.message); 
        stime2 = cat_io_cmd(' ','g5','',job.extopts.verb-1); 
      end
    end
    
    % voxel size
    vx_vol = sqrt(sum(res.image(1).mat(1:3,1:3).^2)); 
    
    
    % correct background
    if isfield(res,'bge')
      P(:,:,:,end) = max( cat_vol_ctype( res.bge * 255 ) , P(:,:,:,end) ); 
      for i=1:size(P,4)-1, P(:,:,:,i) = P(:,:,:,i) .* cat_vol_ctype(1 - res.bge); end
    end


    % initial definition for threshold function 
    Ycls = cell(1,size(P,4)); 
    Yb0  = cat_vol_smooth3X( sum(P(:,:,:,1:3),4), 4/mean(vx_vol)) > 128; 
    Yb0  = cat_vol_morph( cat_vol_morph( Yb0 , 'do' ,3,vx_vol) , 'lc' , 2,vx_vol);  
    for k1 = 1:size(P,4)
      if k1<4
        Ycls{k1} = P(:,:,:,k1) .* uint8(Yb0);
      else
        Ycls{k1} = P(:,:,:,k1) .* uint8(~Yb0);
      end        
    end
    clear Yb0; 
    
    
    % T3th - brain tissue thresholds
    % Use median for WM threshold estimation to avoid problems in case of WMHs!
    if ~exist('T3th','var') || any( isnan(T3th) )
      clsint  = @(x) cat_stat_nanmedian(Ysrc(Ycls{x}>128));
      clsints = @(x) round( sum(res.mn(res.lkp==x) .* res.mg(res.lkp==x)') * 10^5)/10^5;
      T3th = [clsint(3) clsint(1) clsint(2)];
      for i=1:3, if isnan(T3th(i)), T3th(i) = clsints(i); end; end
    end      

    
    % Yb - skull-stripping
    if ~exist('Yb','var')
      if (isfield(job,'useprior') && ~isempty(job.useprior) ) && ... 
         (isfield(res,'ppe') && ~res.ppe.affreg.highBG)
        % use longitudinal TPM
        [Ytmp,Ytmp,Yg,Ydiv] = cat_main_updateSPM_gcut0(Ysrc,P,vx_vol,T3th); clear Ytmp;  %#ok<ASGLU>
        Yb  = YbA > 0.5;
        Ybb = Yb; 
        clear YbA
      elseif size(P,4)==4 || size(P,4)==3
        cat_io_cprintf('warn','\n  Skull-stripped input - refine original mask                                 ')
        [Yb,Ybb,Yg,Ydiv,P] = cat_main_updateSPM_skullstriped(Ysrc,P,res,vx_vol,T3th); %#ok<ASGLU>
        clear Ybb Yg Ydiv 
      elseif isfield(job.extopts,'inv_weighting') && job.extopts.inv_weighting
        Ym0 = single(P(:,:,:,3))/255 + single(P(:,:,:,1))/255 + single(P(:,:,:,2))/255;
        Ym0 = cat_vol_smooth3X(Ym0,4/mean(vx_vol)); 
        Yb  = Ym0 > min(0.5,max(0.25, job.extopts.gcutstr)); clear Ym0
        Yb  = cat_vol_morph(cat_vol_morph(Yb,'lo'),'c');
      else
        try
          Yb = cat_main_APRG(Ysrc,P,res,double(T3th)); 
        catch
          try
            cat_io_addwarning('cat_main_updateSPM:runbackup:APRGerr','IgnoreErrors: cat_main_updateSPM APRG failed - run GCUT function.\n',1); 
            Yb = cat_main_updateSPM_gcut0(Ysrc,P,vx_vol,T3th);
          catch   
            cat_io_addwarning('cat_main_updateSPM:runbackup:GCUTerr','IgnoreErrors: cat_main_updateSPM GCUT failed - use SPM brain mask.\n',1); 
            Yb = cat_vol_morph( sum(P(:,:,:,1:3) > 128,4) , 'lc' , 2);  
          end
        end
      end
    end

    % Yb0 - save brain mask using SPM12 segmentations for later use
    if ~exist('Yb0','var')
      Ym0 = single(P(:,:,:,3))/255 + single(P(:,:,:,1))/255 + single(P(:,:,:,2))/255;
      Yb0 = (Ym0 > min(0.5,max(0.25, job.extopts.gcutstr)));
      Yb0 = cat_vol_morph(cat_vol_morph(Yb0,'lo'),'c');
      clear Ym0
    end

    % final definition with corrected probability maps 
    Ycls = cell(1,size(P,4)); 
    Ysb  = zeros(size(Ysrc),'single'); 
    Ysnb = zeros(size(Ysrc),'single'); 
    for k1 = 1:size(P,4)
      if k1<4
        Ycls{k1} = P(:,:,:,k1) .* uint8(Yb);
        Ysb = Ysb + single(Ycls{k1}); 
      else
        Ycls{k1} = P(:,:,:,k1) .* uint8(~Yb);
        Ysnb = Ysnb + single(Ycls{k1}); 
      end
    end
    for k1 = 1:size(P,4)
      if k1<4
        Ycls{k1}(Yb>=0.5) = cat_vol_ctype(single( (Ycls{k1}(Yb>=0.5))) ./ Ysb(Yb>=0.5)  * 255 ); 
      else
        Ycls{k1}(Yb< 0.5) = cat_vol_ctype(single( (Ycls{k1}(Yb< 0.5))) ./ Ysnb(Yb< 0.5) * 255 ); 
      end
    end
    clear Ysb Ysnb;
    
  end

  
  %% prepared for improved partitioning - RD20170320, RD20180416
  %  Update the initial SPM normalization by a fast version of Shooting 
  %  to improve the skull-stripping, the partitioning and LAS.
  %  We need stong deformations in the ventricle for the partitioning 
  %  but low deformations for the skull-stripping. Moreover, it has to 
  %  be really fast > low resolution (3 mm) and less iterations. 
  %  The mapping has to be done for the TPM resolution, but we have to 
  %  use the Shooting template for mapping rather then the TPM because
  %  of the cat12 atlas map.
  %  
  %  #### move fast shooting to the cat_main_updateSPM function ####
  % 
  stime2 = cat_io_cmd(sprintf('  Update registration'),'g5','',job.extopts.verb,stime2); 

  res2 = res; 
  job2 = job;
  job2.extopts.bb             = 1; 
  job2.extopts.verb           = 0;      % do not display process (people would may get confused) 
  job2.extopts.vox            = abs(res.tpm(1).mat(1));  % TPM resolution to replace old Yy  
  if job.extopts.regstr>0
    job2.extopts.regstr       = 15;     % low resolution 
    job2.extopts.reg.nits     = 16;     % less iterations
    job2.extopts.reg.affreg   = 0;      % new affine registration
    job2.extopts.shootingtpms(3:end) = [];             % remove high templates, we only need low frequency corrections
    res2 = res; 
    res2.do_dartel            = 2;      % use shooting
  else
    fprintf('\n');
    job2.extopts.verb         = 0; 
    job2.extopts.vox          = abs(res.tpm(1).mat(1));  % TPM resolution to replace old Yy 
    job2.extopts.reg.iterlim  = 1;      % only 1-2 inner iterations
    job2.extopts.reg.affreg   = 0;      % new affine registration
    res2.do_dartel            = 1;      % use dartel
  end
  if isfield(res,'Ylesion') && sum(res.Ylesion(:)>0)
    [trans,res.ppe.reginitp,res.Affine] = cat_main_registration(job2,res2,Ycls(1:2),Yy,res.Ylesion); 
  else
    [trans,res.ppe.reginitp,res.Affine] = cat_main_registration(job2,res2,Ycls(1:2),Yy); 
  end
  Yy2  = trans.warped.y;
 
  % Shooting did not include areas outside of the boundary box
  %
  % ### add to cat_main_registration?
  %
  try
    Ybd = true(size(Ysrc)); Ybd(3:end-2,3:end-2,3:end-2) = 0; Ybd(~isnan(Yy2(:,:,:,1))) = 0; Yy2(isnan(Yy2))=0; 
    for k1=1:3
      Yy2(:,:,:,k1) = Yy(:,:,:,k1) .* Ybd + Yy2(:,:,:,k1) .* (1-Ybd);
      Yy2(:,:,:,k1) = cat_vol_approx(Yy2(:,:,:,k1),'nn',vx_vol,3); 
    end
    Yy = Yy2; 
    clear Yy2; 
  end
  stime2 = cat_io_cmd(' ','g5','',job.extopts.verb-1,stime2); 
  fprintf('%5.0fs\n',etime(clock,stime));
 
  
  % some reports 
  for i=1:numel(Ycls), res.ppe.SPMvols1(i) = cat_stat_nansum(single(Ycls{i}(:)))/255 .* prod(vx_vol) / 1000; end
  
  % display  some values for developers
  if job.extopts.expertgui > 1
    if isfield(job.extopts,'spm_kamap') && job.extopts.spm_kamap 
       cat_io_cprintf('blue',sprintf('    SPM  volumes (CGW = TIV in mm%s):%7.2f +%7.2f +%7.2f = %4.0f\n',...
        native2unicode(179, 'latin1'),res.ppe.SPMvols0([3 1 2]),sum(res.ppe.SPMvols0(1:3))));    
       cat_io_cprintf('blue',sprintf('    AMAP volumes (CGW = TIV in mm%s):%7.2f +%7.2f +%7.2f = %4.0f\n',...
        native2unicode(179, 'latin1'),res.ppe.SPMvols1([3 1 2]),sum(res.ppe.SPMvols1(1:3))));    
    else
      cat_io_cprintf('blue',sprintf('    SPM volumes pre  (CGW = TIV in mm%s): %7.2f +%7.2f +%7.2f = %4.0f\n',...
        native2unicode(179, 'latin1'),res.ppe.SPMvols0([3 1 2]),sum(res.ppe.SPMvols0(1:3)))); 
      cat_io_cprintf('blue',sprintf('    SPM volumes post (CGW = TIV in mm%s): %7.2f +%7.2f +%7.2f = %4.0f\n',...
        native2unicode(179, 'latin1'),res.ppe.SPMvols1([3 1 2]),sum(res.ppe.SPMvols1(1:3)))); 
    end
  end
  
  
end
function [Yb,Ybb,Yg,Ydiv,P] = cat_main_updateSPM_skullstriped(Ysrc,P,res,vx_vol,T3th)
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
    Yb   = Yp0>=0.5/3 & Ysrc>0; 
    Ybb  = cat_vol_ctype(Yb)*255; 

    P(:,:,:,6) = P(:,:,:,4); 
    P(:,:,:,4) = zeros(size(Ysrc),'uint8');
    P(:,:,:,5) = zeros(size(Ysrc),'uint8'); 
    res.lkp = [res.lkp 5 6];
    res.mn  = [res.mn(1:end-1),0,0,0];
    res.mg  = [res.mg(1:end-1);1;1;1];
    res.vr(1,1,numel(res.lkp)-1:numel(res.lkp)) = 0;
     
    [Ysrcb,BB] = cat_vol_resize(Ysrc,'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yp0>1/3); clear Yp0;
    Yg   = cat_vol_grad(Ysrcb/T3th(3),vx_vol);
    Ydiv = cat_vol_div(Ysrcb/T3th(3),vx_vol);
    Yg   = cat_vol_resize(Yg ,'dereduceBrain',BB);
    Ydiv = cat_vol_resize(Ydiv ,'dereduceBrain',BB);
end
function [Yb,Ybb,Yg,Ydiv] = cat_main_updateSPM_gcut0(Ysrc,P,vx_vol,T3th)
    % brain mask
    Ym   = single(P(:,:,:,3))/255 + single(P(:,:,:,1))/255 + single(P(:,:,:,2))/255;
    Yb   = (Ym > 0.5);
    Yb   = cat_vol_morph(cat_vol_morph(Yb,'lo'),'c');
    Ybb  = cat_vol_ctype(cat_vol_smooth3X(Yb,2)*256); 

    [Ysrcb,BB] = cat_vol_resize({Ysrc},'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yb);
    Yg   = cat_vol_grad(Ysrcb/T3th(3),vx_vol);
    Ydiv = cat_vol_div(Ysrcb/T3th(3),vx_vol);
    Yg   = cat_vol_resize(Yg   ,'dereduceBrain',BB);
    Ydiv = cat_vol_resize(Ydiv ,'dereduceBrain',BB);
end
function [Yb,Ybb,Yg,Ydiv] = cat_main_updateSPM_gcutold(Ysrc,P,res,vx_vol,T3th)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T1 only > remove in future if gcut is removed too!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clsint = @(x) round( sum(res.mn(res.lkp==x) .* res.mg(res.lkp==x)') * 10^5)/10^5;
    voli   = @(v) (v ./ (pi * 4./3)).^(1/3);   % volume > radius
    Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));

    % old skull-stripping
    Yp0  = single(P(:,:,:,3))/255/3 + single(P(:,:,:,1))/255*2/3 + single(P(:,:,:,2))/255;
    brad = voli(sum(Yp0(:)>0.5).*prod(vx_vol)/1000); 
    [Ysrcb,Yp0,BB] = cat_vol_resize({Ysrc,Yp0},'reduceBrain',vx_vol,round(6/mean(vx_vol)),Yp0>1/3);
    %Ysrcb = max(0,min(Ysrcb,max(T3th)*2));
    BGth = min(cat_stat_nanmean(Ysrc( P(:,:,:,6)>128 )),clsint(6));
    Yg   = cat_vol_grad((Ysrcb-BGth)/diff([BGth,T3th(3)]),vx_vol);
    Ydiv = cat_vol_div((Ysrcb-BGth)/diff([BGth,T3th(3)]),vx_vol);
    Ybo  = cat_vol_morph(cat_vol_morph(Yp0>0.3,'lc',2),'d',brad/2/mean(vx_vol)); 
    BVth = diff(T3th(1:2:3))/abs(T3th(3))*1.5; 
    RGth = double(diff(T3th(2:3))/abs(T3th(3))*0.1); 
    Yb   = single(cat_vol_morph((Yp0>1.9/3) | (Ybo & Ysrcb>mean(T3th(2)) & ...
           Ysrcb<T3th(3)*1.5 & Yg<0.5),'lo',max(0,0.6/mean(vx_vol)))); 
    
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
    clear Ysrcb Ybo;
end
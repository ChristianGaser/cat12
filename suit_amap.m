function P = suit_amap(job,varargin)
% SUIT extenstion 
% =========================================================================
% This function replace the SPM segmentation created by the SUIT TPM by 
% an AMAP segmentation used also in CAT12. It further includes a maximum 
% based bias correction, an intensity-scaling, and a (SA)NLM denosing step.
% =========================================================================
  
% #########################################################################
% TODO; 
%  - rewrite segmentation (interative separation of GM)
%  - keep it simple or at least transparent (parameter)

  if ischar(job)
    switch job
      case 'run'
        Yp0s = varargin{1};
        Ymi  = varargin{2};
        vx_vol = varargin{3};

        % Intensity normalization of Ym
        P.mm{1} = 'tmp'; 
        Ymi = suit_amap_quickscale(Ymi,Yp0s,vx_vol,P,1);
    
        [job,suit] = suit_amap_update_job(struct('gray',{{'tmp'}}));

        [Yg,gth] = suit_amap_getYg(Ymi,Yp0s);

        Ycbd = cat_vol_smooth3X(cat_vbdist(single(Yp0s>0)  , true(size(Yp0s)), vx_vol) - ... 
                                cat_vbdist(single(Yp0s==0) , true(size(Yp0s)), vx_vol),2);
 
        [~,~,~,Yc2h] = suit_amap_cleanup_cerebellum(Ymi,Yp0s,Ycbd,Yg,Yg<-inf,vx_vol,1);
        
    
        %  Local adaptive segmentation
        %  ====================================================================   
        %  * the local intensity normalization allows correction of artifacts 
        %    and further denoising 
        [~,Yml] = suit_amap_LAS2(Ymi,Yp0s,Yc2h,Yg,gth,job);
 
        Ycb = Yp0s>1; 

        Ywd = Ycb .* (cat_vbdist(single(Yml>0.5),Ycb) + cat_vbdist(single(Yml>0.55),Ycb) + cat_vbdist(single(Yml>0.6),Ycb)) / 3;
        Ycs = max(0,min(1,-cat_vol_div((Ywd).^2))); Ycs = (Ycs + cat_vol_smooth3X(Ycs,4).^4)/2; 
        % WM skeleton ...
        Ycd = Ycb .* (cat_vbdist(single(Yml<0.5),Ycb) + cat_vbdist(single(Yml<0.55),Ycb) + cat_vbdist(single(Yml<0.6),Ycb)) / 3;
        Yws = max(0,min(1,-cat_vol_div((Ycd)))).^.25; Yws = Yws + smooth3(Yws); Yws = Yws/2; Yws = Yws.^.5;
        
        %% mixed measure
        P.Yppi  = Ycb .* max(0,min(1,(Yml*4-1)/3  - Ycs/3)); % + Yws/3));
        P.Yml   = Yml; 
        


        [Ycba,Ybscb,cbac,cbas,suit] = suit_amap_parcelate(Yml,Yp0s,Yp0s,Ycb,Yg,gth,vx_vol,suit,1);
    
         Ycore  = cat_vol_morph( Ycba==cbac{cat_io_contains(cbac(:,2),'core'),1} ,'de',1.5);
        if sum(Ycore(:))<100, Ycore = Ycba==cbac{cat_io_contains(cbac(:,2),'core'),1}; end
        Yhull  = cat_vol_morph( ~Ycb ,'de',1.5); 
        speedlim = @(Y,x) min(1 - x,max(x, Y)) ./ Ycb; % div by 0 to get nan
        [~,~,Ycbhd] = cat_vol_eidist(single(Yhull), speedlim(1 - max(0.001,(Yml*4)-1)/3,.01).^16); 
        [~,~,Ycbcd] = cat_vol_eidist(single(Ycore), speedlim(    max(0.001,(Yml*4)-1)/3,.01).^16);  
      
        %% create a percentage map
        Yt  = min(smooth3(cat_vol_morph(Ycb,'de',2)),Ycbhd ./ max(eps,Ycbhd + Ycbcd)); % clear Ycbhd Ycbcd
        Yt  = Yt .* max(0,min(1,1 - (Yt - Yml).^2)) .* Ycb; % correct areas (close to the core) that have grown to fast
        Yta = Yt .^ (0.25 \ median(Yt(round(Yp0s(:))==2 & Yt(:)>0 & Yt(:)<1)));

         P.Yta   = Yta; 
        return
  
    end
  end


  % update job struction with default settings
  [job,suit] = suit_amap_update_job(job);

  % prepare filenames                           
  P = suit_amap_get_filenames(job);

  % load extracted volumes
  for fi = 1:numel( P.m )
    cat_io_cprintf([.5 .5 .5],'  Run %s.',spm_file(P.m{fi},'basename')); stime = datetime('now'); 


    % Load SPM data: Yp0s-segment, Ym-T1w, Vmsk-brainsask
    [Yp0s,Ym,Vmsk,Vgm,vx_vol] = suit_amap_load_segments(P,fi);
    

    % Intensity normalization of Ym
    Ym  = suit_amap_quickscale(Ym,Yp0s,vx_vol,P,fi);
    % correction of severe GM underestimates
    if 1
      Yp0x = (cat_vol_median3(Ym,Yp0s>1)*3);
      %Ymsk = smooth3(Yp0x<2.9 & Yp0x>1.5 & Yp0s>1.5 & Yp0s<2.95)>.4;  
      Ymsk = smooth3(Yp0x>1.5 & Yp0s>1.5 & ~cat_vol_morph(Yp0x>2.75 & Yp0s>2.75,'do',2,vx_vol))>.4;  
      Yp0s(Ymsk) = min(Yp0x(Ymsk),Yp0s(Ymsk)); 
    end

    % prepare basic maps: 
    % (1) General distance map from the cerebellar boundary based on the  
    %     SPM CSF/GM boundary with pos. outward and neg. inward values.
    Ycbd = cat_vol_smooth3X(cat_vbdist(single(Yp0s>0)  , true(size(Yp0s)), vx_vol) - ... 
                            cat_vbdist(single(Yp0s==0) , true(size(Yp0s)), vx_vol),2);
    % (2) get gradient map
    [Yg,gth] = suit_amap_getYg(Ym,Yp0s);
    % (3) MP2rage noise
    Ybv = smooth3(cat_vol_localstat(Ym,Ym>0,1,4).*(Ym + (Ym-smooth3(Ym)))*3.*(2-Yp0s))>.125; 


    % Cerebral Cleanup and update of input segmentation 
    % (3) percentage position map between core and hull 
    %     (need for bias correction, LAS, and CB extraction)
    [Ycb,~,Ycbd,Yc2h] = suit_amap_cleanup_cerebellum(Ym,Yp0s,Ycbd,Yg,Ybv,vx_vol,1);
    Yp0s = max(Yp0s .* Ycb, min(2,Ym .* Ycb)); % apply correction



    % Whatfore - skull-stripping???
    %{
      Ycore  = cat_vol_morph( Yp0s>2.5 ,'ldo',1,vx_vol);
      Yhull  = cat_vol_morph( ~cat_vol_morph( cat_vol_morph( Yp0s>0 ,'ldo',2,vx_vol) ,'ldc',4,vx_vol),'de',2,vx_vol); 
      speedlim = @(Y,x) min(1 - x,max(x, Y)) ./ (~Yhull); % div by 0 to get nan
      [~,~,Ycbhd] = cat_vol_eidist(single(Yhull), speedlim(1 - max(0.01,Ym),.01).^16); 
      [~,~,Ycbcd] = cat_vol_eidist(single(Ycore), speedlim(    max(0.01,Ym),.01).^16);  
      Ycb0 = Ycbcd./Ycbhd < 1; 
    
      %% create a percentage map
      Yt  = min(cat_vol_smooth3X(~Yhull,4),Ycbhd ./ max(eps,Ycbhd + Ycbcd)); % clear Ycbhd Ycbcd
      Yt  = Yt .* max(0,min(1,1 - (Yt - Ym).^2)); % correct areas (close to the core) that have grown to fast
      Yta = Yt .^ (0.25 \ median(Yt(round(Yp0s(:))==2 & Yt(:)>0 & Yt(:)<1)));
    %}



    % refine segmentation by AMAP
    %  ====================================================================
    % * AMAP might help here ?
    %   - when you have to stretch the GM then the bias correction might be no problem
    %   - would be better after the fine bias correction but this might profit by AMAP too
    % * maybe to use the AMAP here would be good to stablize the segments and avoid side defects 
    %   this might allow to avoid the distbias modulation 
    %   - as it is quick you might use it in the inner loop ? 
%##### ... not sure if this is good - I don't get the central line any longer ... 
    if 0 %job.amap 
      Yp0a = suit_amap_AMAP(Ym,Yp0s,vx_vol,job.sanlm,0);
    else
      Yp0a = Yp0s; % SPM underestimate the CSF in aging!
    end



    %%  Fine bias correction - assuming that GM is in average GM
    %  ====================================================================
    Ymi  = suit_amap_biascorr(Ym,Yp0a,Yg,Yc2h,Ycbd,gth,job,vx_vol); 
    
    %% basic evaluation 
    [Tth,suit] = suit_amap_eval(Ym,Ymi,Yp0a,Yp0s,vx_vol,suit,fi);



    %  Local adaptive segmentation
    %  ====================================================================   
    %  * the local intensity normalization allows correction of artifacts 
    %    and further denoising 
    [Ymg,Yml] = suit_amap_LAS2(Ymi,Yp0a,Yc2h,Yg,gth,job);
    if suit(fi).stripped > .1
      Ymg = max(Ymg,Yp0s/3); 
      Yml = max(Yml,(Yp0s+1)/4);
    end



    %%  Cerebral Cleanup
    %  ====================================================================
    %  - a good skull-stripping is highly important for the registration!
    %  - we therefore focus on GM and WM and the mask ignores external CSF
    %    that is not between structures
    Ycb = suit_amap_cleanup_cerebellum(Yml,Yp0s,Ycbd,Yg,Ybv,vx_vol,1);
    suit_amap_saveimg(Ycb, P.msk1{fi}, Vmsk, Vgm, job, job.mask, 'SUIT cerebellar mask', 'uint8');
    suit(fi).vol_cbmsk0 = sum(min(1,Yp0s(:)>0)) * prod(vx_vol) / 1000;
    suit(fi).vol_cbmsk1 = sum(min(1,Ycb (:)>0)) * prod(vx_vol) / 1000;
    


    %% final cerebellar models
    %  ====================================================================
    suit = suit_amap_write(Ymg,Yml,Yp0a,Ycb,Vmsk,Vgm,vx_vol,P,job,suit,fi);

    

    %% cerebellar structures (core, hull, brainstem, peduncles)
    %  needed for cerebellar coordinate system and cerebellar surface reconstruction
    %  ====================================================================
    [Ycba,Ybscb,cbac,cbas,suit] = suit_amap_parcelate(Yml,Yp0a,Yp0s,Ycb,Yg,gth,vx_vol,suit,fi);
    suit_amap_saveimg(Ycba, P.cba{fi}, Vmsk, Vgm, job, job.th0, cbas, 'uint8');
    


    %% cerebellar coordinates (x-transversal, y-longitudinal, z-central axis) 
    %  Unclear how much we need when ...
%########  * the potential maps depend on the image resolution 
    %  ====================================================================
    if job.axis==1 || job.axis==4
      %% x-transversal
      Ycoord{1} = suit_amap_get_transversal(Yml,Ybscb,Ycba,cbac);
      suit_amap_saveimg(Ycb .* (Ycoord{1}/2 + .5*(Ycoord{1}~=0)), P.ax1t{fi}, Vmsk, Vgm, job, job.axis==1 | job.axis==4,...
        'SUIT transversal cerebellar coordinates (x-axis)' ,'uint8');
    end
    if job.axis==2 || job.axis==4
      %% y-longitudinal
      Ycoord{2} = suit_amap_get_longitudinal(Yml,Ybscb,Ycba,cbac);
      suit_amap_saveimg(Ycb .* (Ycoord{2}/2 + .5*(Ycoord{1}~=0)), P.ax2l{fi}, Vmsk, Vgm,job, job.axis==2 | job.axis==4,...
        'SUIT longitudinal cerebellar coordinates (y-axis)','uint8');
    end
    % z-central axis
    [Yth0,Ycoord{3},suit] = suit_amap_cb0thickness( Ycba, cbac, vx_vol, suit, fi);
    suit_amap_saveimg(Yth0, P.th0{fi}, Vmsk, Vgm, job, job.th0, 'cerebellar thickness (core-hull-distance)', 'uint8');
    suit_amap_saveimg(Ycb .* Ycoord{3}, P.ax3c{fi}, Vmsk, Vgm, job, job.axis==3 | job.axis==4,...
      'SUIT central cerebellar coordinates (z-axis)','uint8');
   
    

    %% cortical percentage position estimate 
    %  ====================================================================
    %  - as an optimized map for registration and surface reconstruction
    %  - use of skull-stripping to clearly outline the CSF but also the 
    %    ongoing branches by a WM skeleton
    % CSF skeleton (this will enhance sulcal areas also in young subjects)
    Ywd = Ycb .* (cat_vbdist(single(Yml>0.5),Ycb) + cat_vbdist(single(Yml>0.55),Ycb) + cat_vbdist(single(Yml>0.6),Ycb)) / 3;
    Ycs = max(0,min(1,-cat_vol_div((Ywd).^2))); Ycs = (Ycs + cat_vol_smooth3X(Ycs,4).^4)/2; 
    % WM skeleton ...
    Ycd = Ycb .* (cat_vbdist(single(Yml<0.5),Ycb) + cat_vbdist(single(Yml<0.55),Ycb) + cat_vbdist(single(Yml<0.6),Ycb)) / 3;
    Yws = max(0,min(1,-cat_vol_div((Ycd)))).^.25; Yws = Yws + smooth3(Yws); Yws = Yws/2; Yws = Yws.^.5;
    %% mixed measure
    Yx  = Ycb .* max(0,min(1,(Yml*4-1)/3  - Ycs/3)); % + Yws/3));
    suit_amap_saveimg(Yx, P.cba{fi}, Vmsk, Vgm, job, 1,...
      'central deformation map','uint8');  
    %% CBA2
    if 0
      Yx  = Ycb .* max(0,min(1.25,(Yml*4-1)/3  - Ycs/3 + Yws/3));
      Yx  = min(1,max(Yx*4/5,min(Ycb .* Yml.^2, smooth3(Ycba==cbac{cat_io_contains(cbac(:,2),'brainstem'),1}) )));
      
      suit_amap_saveimg(Yx, strrep( P.cba{fi} ,'cba','cba2'), Vmsk, Vgm, job, 1,...
        'central deformation map 2','uint8');  
    end
    % use general core-hull-bias to weight WM structures
    %suit_amap_saveimg(Yx .* smooth3(max(max(eps,Ycb .* (Yp0s-2)),0.5.*Ycb + 0.5.*Ycoord{3}).^.5), ...
    %  strrep( P.cba{fi} ,'cba','cba2'), Vmsk, Vgm, job, 1,...
    %  'central deformation map','uint8');  

  % use of a real WM percentage map?
  % - it would be possible to do this for the CSF and WM range as some 
  %   more abstract model but the normalized and skeletonized map seems 
  %   to be more powerful 
  %  [Ywt,Ywp] = cat_vol_pbtp(max(2/3,1.5-Yml)*3, Ycd, 0*Ycd );
  %  Ywt = cat_vol_approx(Ywt); 


   
   %% cortical weighing map
   %  ====================================================================
   %  - use of the Eikonal distance estimation to define a more abstract 
   %    position map between the cerebellar core and hull
   if 1
     try
       %% estimate of weighed distance map (traveling time) ... processing time about 1 minute 
        Ycore  = cat_vol_morph( Ycba==cbac{cat_io_contains(cbac(:,2),'core'),1} ,'de',1.5);
        if sum(Ycore(:))<100, Ycore = Ycba==cbac{cat_io_contains(cbac(:,2),'core'),1}; end
        Yhull  = cat_vol_morph( ~Ycb ,'de',1.5); 
        speedlim = @(Y,x) min(1 - x,max(x, Y)) ./ Ycb; % div by 0 to get nan
        [~,~,Ycbhd] = cat_vol_eidist(single(Yhull), speedlim(1 - max(0.001,(Yml*4)-1)/3,.01).^16); 
        [~,~,Ycbcd] = cat_vol_eidist(single(Ycore), speedlim(    max(0.001,(Yml*4)-1)/3,.01).^16);  
      
        %% create a percentage map
        Yt  = min(smooth3(cat_vol_morph(Ycb,'de',2)),Ycbhd ./ max(eps,Ycbhd + Ycbcd)); % clear Ycbhd Ycbcd
        Yt  = Yt .* max(0,min(1,1 - (Yt - Yml).^2)) .* Ycb; % correct areas (close to the core) that have grown to fast
        Yta = Yt .^ (0.25 \ median(Yt(round(Yp0s(:))==2 & Yt(:)>0 & Yt(:)<1)));
      catch
        Yta = Yml; 
      end
      %%
      suit_amap_saveimg(Yx .* smooth3(Yta).^.5, ... (Yx .* smooth3(Yta) .* (Yta)).^1.2, ...
        strrep( P.cba{fi} ,'cba','cba3'), Vmsk, Vgm, job, 1, 'central deformation map','uint8');  
   end

    if 0
    %% CBA4
      levels = 0.9:-0.05:0.6; Ymis = smooth3(Ymi); 
      Ywd = Ymi*0; for i=levels,  Ywd = Ywd + cat_vbdist(single(Ymis>i & Ycb)) / numel(levels); end
      Yd  = cat_vol_div(Ywd,1,1); 
      %      Yxx = max(0,Ymi*4-1)/2.*Ycb - max(0,-(Yd)/2).^.8 - (max(0,Ywd)/40).^.8; 
      Yxx  = max(0,min(1-Ywd/200,0.5 - Ywd/40 + smooth3(Yd)/2 + min(1,max(0,smooth3(Ymi.*Ycb)*4-2).^4)/2)).^2;
      Yxxx = max(0,Yxx.*Yx.*Ymi);
      % use general core-hull-bias to weight WM structures
      suit_amap_saveimg( max(0,min(1,Yxxx)), ...
        strrep( P.cba{fi} ,'cba','cba4'), Vmsk, Vgm, job, 1, 'central deformation map','uint8');  
    end
    

    %% cerebellar surface measures
    %  ====================================================================
% #### - complexity measures? 
% #### - surface-based - area, folding
    %  

    if 0
      %% 
      %  ====================================================================
      %  (1) use cerebellar orientation map (x-,y-long,z-axial)
      %  (2) (abs) gradient maps >> gx,gy,gz
      st = cat_io_cmd('Create local coordinate systems','g5','',1);
      Ygcb  = cell(3,3);                      % cerebellar coordinate system
      Ygcbw = 1 - max(0,min(1,abs(Ym*2-1)));  % filter masking (1=filter)
      gcs   = [10 2 1];
      for ci = 1:3; [Ygcb{ci,1},Ygcb{ci,2},Ygcb{ci,3}] = cat_vol_gradient3( Ycoord{ci} ); end
      %  (3) orthogonalize
      Ygcb = orthogonalize(Ygcb);
      st    = cat_io_cmd('Next?','g5','',1,st);
  
      %%
      [Ymlr,Ycbr,redR]  = cat_vol_resize({Yml,Ycb} ,'reduceV',vx_vol,1,32,'meanm'); 
      Ygcbr = cell(size(Ygcb)); for i=1:numel(Ygcb), Ygcbr{i}  = cat_vol_resize(Ygcb{i},'reduceV',vx_vol,.9,32,'meanm'); end
      %  (4) filter cb field 
      st    = cat_io_cmd('Filter local coordinate systems','g5','',1,st);
      Yfmsk = min(1,1-abs((2*cat_vol_smooth3X(Ymlr,2))-1)) .* Ycbr; 
      Ygcbr2 = Ygcbr; for i=1, Ygcbr2 = smoothcoord(Ygcbr2,smooth3(single(Yfmsk)),2); Ygcbr2 = orthogonalize(Ygcbr2); end
      st    = cat_io_cmd('Next?','g5','',1,st);
      
      %%
      st   = cat_io_cmd('Filter data by local coordinate systems','g5','',1,st);
      Ymlrs = Ymlr; for i=1:10, Ymlrs = oriented_smooth(Ygcbr2, Ymlrs, single(Yfmsk), [1 1 1]); end  %max(Ycb*.2,smooth3(Yfmsk))
      st   = cat_io_cmd('Next?','g5','',1,st);
    end




    %% simple quality measures:  
    %  ====================================================================
    %  - NER = noise to edge ratio (larger worse)
    %  - NCR = noise to contrast ratio (larger worse)
    %  - add FEC = Fast Euler Characteristic ? 
    %  - add RMSE of diff Ymi and Yml
    Ymlsd = cat_vol_localstat( Ymi , Yml>3.5/4 & Ycb>.5, 2, 4); 
    Ymlgd = cat_vol_localstat( Ymi , Yml>1.5/4 & Yml<3.5/4 & Ycb>.5, 2, 4); 
    suit(fi).qc_NER = mean(Ymlsd(Ymlsd>0)) / mean(Ymlgd(Ymlgd>0)); 
    suit(fi).qc_NCR = mean(Ymlsd(Ymlsd>0)) / diff(Tth(2:3)); 
    clear Ymlsd Ymlgd; 
    
    % write XML/json/mat
    cat_io_xml(P.suit{fi},suit(fi)); 
    
    cat_io_cprintf([.5 .5 .5],sprintf(' %0.0fs.\n',seconds(datetime('now') - stime))); 
  end
end
% =========================================================================
function suit_amap_saveimg(Y,fname,Vmsk,Vgm,job,doit,desc,dt)
%suit_amap_saveimg. Save nifti volumes. 
  if doit
    Vci          = Vmsk; 
    Vci.dt(1)    = spm_type(strrep(dt,'x3',''));
    switch dt  
      case 'uint8x3',  Vci.pinfo(1) = 1/50;
      case 'int8',     Vci.pinfo(1) = 1/127; 
      case 'uint8',    Vci.pinfo(1) = 1/255;
      case 'float',    Vci.pinfo(1) = 1; 
      otherwise 
        error('unknow type dt=%s',dt); 
    end

    if job.cropped
      Vci.fname  = spm_file(fname,'prefix','c_');  
    else
      Vci.fname  = fname;
    end
    Vci.descrip  = desc; 

    % first write the cropped version 
    if isa(Y,'uint8') || max(Y(:))==255 || max(Y(:))==127
      spm_write_vol(Vci, double(Y) / double(max(Y(:))));
    else
      spm_write_vol(Vci, Y);
    end

    % now save the full ouput
    if job.full 
      Vgm.pinfo(1)  = Vci.pinfo(1); 
      Vgm.dt(1)     = spm_type(strrep(dt,'x3',''));
      [Vcir,Ycia]   = cat_vol_imcalc(Vci, Vgm, 'i1'); 
      Vcir.fname    = fname;  
    %  Vcir.pinfo(1) = Vci.pinfo(1); 
    %  Vcir.dt(1)    = spm_type(strrep(dt,'x3',''));
      spm_write_vol(Vcir, Ycia); 
    end
  end
end
% =========================================================================
% Preparation 
% =========================================================================
function [Yg,gth] = suit_amap_getYg(Ym,Yp0s)
% get normalized gradient map with tissue threshold

  Ygo = cat_vol_grad(Ym) ./ max(eps,abs(Ym)); 

  % correct for cerbellar GM area 
  Ygb = cat_vol_approx( real( max( min(1,abs(Yp0s-2)), max(0,1 - 2 * Ygo ./ median(Ygo(round(Yp0s(:))==2)))) ) , 'rec'); 
  Ygb = Ygb ./ median(Ygb(:));
  
  % correction for artifacts within GM
  Yg  = Ygo ./ Ygb; clear Ygb 

  % define tissue threshold 
  gth = min(.6, mean(Yg(Yp0s(:)>2.5 & Yg(:)<.5)) * 2);
end
% =========================================================================
function P = suit_amap_get_filenames(job)
  % prepare input filenmaes if not specified otherwise
  P.m    = spm_file(job.files, 'prefix', 'c_', 'ext', '.nii');             % original (not-bias-corrected!)
  P.mm   = spm_file(job.files, 'prefix', 'c_m', 'ext', '.nii');            % SPM bias-corrected
  P.gm   = spm_file(job.files, 'suffix', '_seg1', 'ext', '.nii');          % SPM GM
  P.wm   = spm_file(job.files, 'suffix', '_seg2', 'ext', '.nii');          % SPM WM
  P.msk  = spm_file(job.files, 'suffix', '_pcereb',  'prefix','c_', 'ext', '.nii'); % mask
  P.spm  = spm_file(job.files, 'prefix', '_seg8', 'ext', '.mat'); %        % SPM unified segmentation data
 
  % updated aps (prefix set later if required)
  P.msk1 = spm_file(job.files, 'suffix', '_pcereb1', 'ext', '.nii');       % updated mask
  P.bc   = spm_file(job.files, 'suffix', job.suffixbc, 'ext', '.nii');     % updatad bias corrected
  
  % 3-class AMAP segmenation maps (2 used shooting outputs as CSF=BG)
  % The 3 class amap model is not fitting to the cerebellum (as class 2 is 
  % not really existing) but presents the basis for comparison. Although an 
  % intensity based 3 class model would be possible it is not really suited
  % and we will go directly to the desired 4 class model.
  P.amap0  = spm_file(job.files, 'suffix', [job.suffixamap '0'], 'ext', '.nii'); % AMAP label map 
  P.amap1  = spm_file(job.files, 'suffix', [job.suffixamap '1'], 'ext', '.nii'); % AMAP GM
  P.amap2  = spm_file(job.files, 'suffix', [job.suffixamap '2'], 'ext', '.nii'); % AMAP WM
  P.amap3  = spm_file(job.files, 'suffix', [job.suffixamap '3'], 'ext', '.nii'); % AMAP CSF
  
  % 4-class cerebellar model (3 used shooting outputs)
  % The class number should play now role for shooting so we use here a
  % shell like structure CSF-GMlow-GMhigh-WM. 
  P.cb0  = spm_file(job.files, 'suffix', [job.suffixcb '0'], 'ext', '.nii'); % label map that respresents the range CSF to WM
  P.cb1  = spm_file(job.files, 'suffix', [job.suffixcb '1'], 'ext', '.nii'); % intensity based CSF
  P.cb2  = spm_file(job.files, 'suffix', [job.suffixcb '2'], 'ext', '.nii'); % intensity based lower GM (CSF/GM)
  P.cb3  = spm_file(job.files, 'suffix', [job.suffixcb '3'], 'ext', '.nii'); % intensity based higher GM (GM/WM)
  P.cb4  = spm_file(job.files, 'suffix', [job.suffixcb '4'], 'ext', '.nii'); % intensity based WM

  % 4-class cerebellar model (3 used shooting outputs)
  % The class number should play now role for shooting so we use here a
  % shell like structure CSF-GMlow-GMhigh-WM. 
  P.cc0  = spm_file(job.files, 'suffix', [job.suffixcc '0'], 'ext', '.nii'); % label map that respresents the range CSF to WM
  P.cc1  = spm_file(job.files, 'suffix', [job.suffixcc '1'], 'ext', '.nii'); % intensity based CSF
  P.cc2  = spm_file(job.files, 'suffix', [job.suffixcc '2'], 'ext', '.nii'); % intensity based lower GM (CSF/GM)
  P.cc3  = spm_file(job.files, 'suffix', [job.suffixcc '3'], 'ext', '.nii'); % intensity based higher GM (GM/WM)
 
  % cerebellar location maps
  P.ax1t = spm_file(job.files, 'suffix', '_ax1t', 'ext', '.nii'); % intensity based WM
  P.ax2l = spm_file(job.files, 'suffix', '_ax2l', 'ext', '.nii'); % intensity based WM
  P.ax3c = spm_file(job.files, 'suffix', '_ax3c', 'ext', '.nii'); % intensity based WM
  
  % atlas for cerebellar location and thickness maps
  P.cba = spm_file(job.files, 'suffix', '_cba', 'ext', '.nii'); 
  P.cbb = spm_file(job.files, 'suffix', '_cbb', 'ext', '.nii'); 

  % thickness maps
  P.pp  = spm_file(job.files, 'suffix', '_pp0', 'ext', '.nii'); 
  P.th0 = spm_file(job.files, 'suffix', '_th0', 'ext', '.nii'); 
  P.th1 = spm_file(job.files, 'suffix', '_th1', 'ext', '.nii'); 

  P.suit = spm_file(job.files, 'suffix', '_suit', 'ext', '.xml'); 
end
% =========================================================================
function [job,suit] = suit_amap_update_job(job)
  job.files           = job.gray; 

  % internal parameters
  def.suffixbc        = '_bc';   % bias corrected
  def.suffixamap      = '_amap'; % AMAP segmentation with CSF,GM,WM classes and CSF/GM and GM/WM PVE classes 
  def.suffixcb        = '_cggw'; % cereblleum intensity-scaling with 4 classes - CSF, GM/CSF, GM/WM and WM
  def.suffixcc        = '_cgw';  % cereblleum intensity-scaling with 3 classes - CSF, GM and WM#

  def.sanlm           = 1;   % run denoising
  def.cbdistbias      = .05; % correction to avoid bias by changing WM proportion between core and hull
  def.cbseg           = 1;   % run/output cerebellum segment (0-none, 1-write label, 2-write classes, 3-both)      
  def.mask            = 0;   % write updated cerebellar mask (0-no, 1-yes)
  def.amap            = 1;   % run/output AMAP (0-none, 1-write label, 2-write classes, 3-both)    
  def.biascorr        = 0;   % output bias corrected    

  % output maps 
  def.axis            = 0;   % write axis (0-none, 1-transversal, 2-longitudinal, 3-central, 4-all)
  def.th0             = 0;   % write cerebellar thickness map (0-none, 1-core-hull, 2-brachnes, 3-cortical, 4-all) 
  def.th1             = 0;   % write cerebellar thickness map (0-none, 1-core-hull, 2-brachnes, 3-cortical, 4-all) 
  def.atlas           = 1;   % write atlas map (0-no, 1-yes)
 
  % QC map
  def.cropped         = 1;   % save cropped images too (eg. only for evaluation)
  def.full            = 1;   % save full images (this is the input for other functions!)
  def.GMcls           = 0;   % 0 - no-GM-class (only PVE, basic-global-scaling), 1 - one-cls-model (classic), 2 - two-class-model (new)
                             % GM (real) vs. layersurf output (midline-correction)

  % update
  job = cat_io_checkinopt(job,def); 

  % output sturctures
  suit = struct();                    % XML/json output structure
  
end
% =========================================================================
% Main processing
% =========================================================================
function [Ymg,Yml] = suit_amap_LAS(Ymi,Yp0a,Yc2h,Yg,gth,T4th,vx_vol,job)
  Ymio = Ymi; 
  Yp0r = round(Yp0a); 
  Ytisdistprob = (1+job.cbdistbias - job.cbdistbias*Yc2h); 
  %%
  Ymi=Ymio; 
  iter = 3; 
  for si = 1:iter
    %% local intensity estimation 
    % the global value is more save here ! 

    % CSF
    %Ylab{1}  = median( Ymi( Yp0a(:)>0.5 & Yp0a(:)<1.25 & Ymi(:)<mean(Tth(1:2)) & Ymi(:)>Tth(1)/2 )); 
    Ylab{1}  = T4th(1);
    
    % WM
    % * maybe it is not necessary to handle the WM here as we did this already 
    if 1
      Ylab{4}  = Ymi .* (Yp0a>2.5 & Yp0a<3.5 & Yg>0 & Yg<gth*2 & Ymi>(T4th(1)*.25+.75*T4th(4)) & Yc2h>0.25); %.05; 
      Ylab{4}  = cat_vol_localstat( Ylab{4}, Ylab{4}>0, 1,3) * .95; % ### maybe avoid overestimation that we cannot represent in the segmentation
      Ylaba    = cat_vol_median3( Ylab{4}, Ylab{4}>0, Ylab{4}>0 ); % important to avoid bias by outliers!
      Ylab{4}(Yp0a==0) = .95; 
      Ylaba    = cat_vol_approx( Ylaba ); 
      Ylab{4}  = cat_vol_approx( Ylab{4} .* (Ymi ./ Ylaba > .9 & Ymi ./ Ylaba < 1.1) ); 
    else
      Ylab{4}  = T4th(4); 
    end
    
    % GM 
    Ylabg    = Ymi .* (Yp0a>1.5 & Yp0a<2.5 & Yg>0 & Yg<gth*4 & Ymi>mean(T4th(1:2)) & Ymi<mean(T4th(3:4))); 
    Ylabg    = Ylabg .* (Ylabg>prctile(Ylabg(Ylabg(:)>0),10) & Ylabg<prctile(Ylabg(Ylabg(:)>0),90)); % 
    %Ylabg    = cat_vol_median3( Ylabg, Ylabg>0, Ylabg>0 ); % important to avoid bias by outliers!
    Ylaba    = cat_vol_approx( Ylabg ); 
    Ylabg    = Ylabg + (Ylabg==0).*cat_vol_approx( Ylabg .* (Ymi ./ Ylaba > .9 & Ymi ./ Ylaba < 1.1) ); 
    Ylabg    = Ylabg .* Ytisdistprob; % reduced weighting to boundaries 
    Ylabg    = cat_vol_smooth3X(Ylabg,2); 
    
    % lower GM
    Ylab{2}  = ( Ymi .* (Yp0a>1.5 & Yp0a<2.5 & Ymi<Ylabg*.99 & Ymi>mean(T4th(1:2)) & Yg<gth*4)); 
    %Ylaba    = cat_vol_median3( Ylab{2}, Ylab{2}>0, Ylab{2}>0 ); % important to avoid bias by outliers!
    Ylaba    = cat_vol_approx( Ylaba ); 
    Ylab{2}  = Ylab{2} + (Ylab{2}==0).*cat_vol_approx( Ylab{2} .* (Ymi ./ Ylaba > .95 & Ymi ./ Ylaba < 1.05) ); 
    Ylab{2}  = Ylab{2} .* Ytisdistprob; 
    Ylab{2}  = cat_vol_smooth3X(Ylab{2},3); 
    
    % higher GM
    Ylab{3}  = ( Ymi .* (Yp0a>1.5 & Yp0a<2.5 & Ymi>Ylabg/.99 & Ymi<mean(T4th(3:4)) & Yg<gth*4)); 
    %Ylaba    = cat_vol_median3( Ylab{3}, Ylab{3}>0, Ylab{3}>0 ); % important to avoid bias by outliers!
    Ylaba    = cat_vol_approx( Ylaba ); 
    Ylab{3}  = Ylab{3} + (Ylab{3}==0).*cat_vol_approx( Ylab{3} .* (Ymi ./ Ylaba > .95 & Ymi ./ Ylaba < 1.05) ); 
    Ylab{3}  = Ylab{3} .* Ytisdistprob;
    Ylab{3}  = cat_vol_smooth3X(Ylab{3},3); 
    %   clear Ylabg Ylaba; 

    % local intensity normalization 
    Ymg = zeros(size(Ymi)); Ylab23 = (Ylab{2}+Ylab{3})/2; 
    Ymg = Ymg + ( (Ymi>=Ylab{4}               ) .* (3 + (Ymi - Ylab{4})  ./ max(eps,Ylab{4} - Ylab23)) ); % scale highest tissue (WM in T1w)
    Ymg = Ymg + ( (Ymi>=Ylab23  & Ymi<Ylab{4} ) .* (2 + (Ymi - Ylab23)   ./ max(eps,Ylab{4} - Ylab23)) ); % scale second highest tissue (GM in T1w)
    Ymg = Ymg + ( (Ymi>=Ylab{1} & Ymi<Ylab23  ) .* (1 + (Ymi - Ylab{1})  ./ max(eps,Ylab23  - Ylab{1})) ); % scale third highest tissue (CSF in T1w)
    Ymg = Ymg + (  Ymi< Ylab{1}               ) .* (    (Ymi - min(Ymi)) ./ max(eps,Ylab{1} - min(Ymi)) ); % scale background
    Ymg = Ymg / 3; 
   
%%
    if 0% si < 3
      %% magic corretion 
      Ywd  = zeros(size(Ymg)); for i=.1:.2:.9, Ywd = Ywd + 1/5 * cat_vbdist( single(Ymg*3 > 1.75 + i/2 ) , Yp0a>0 & Ymg>1.5/3, vx_vol); end
      Ycd  = zeros(size(Ymg)); for i=.1:.2:.9, Ycd = Ycd + 1/5 * cat_vbdist( single(Ymg*3 < 1.75 + i/2 ) , Yp0a>0 & Ymg<2.5/3, vx_vol); end
      Yth1 = cat_vol_approx( (Ywd + Ycd) .* (Yp0r==2 & Ywd<1e10 & Ycd<1e10) ); 
      Ytbc = cat_vol_smooth3X( (Yp0r==2 & Ycd<1e10 & Ywd<1e10) .* max(0,Ycd - 1.5*Yth1),4); %max(0,Yth1 - median(Yth1(Yp0r(:)==2)))) ); 
      Ymg  = Ymg .* (1+Ytbc); 
      Ytbc = cat_vol_smooth3X( (Yp0r==2 & Ycd<1e10 & Ywd<1e10) .* max(0,Ywd - 1.5*Yth1),4); %max(0,Yth1 - median(Yth1(Yp0r(:)==2)))) ); 
      Ymg  = Ymg .* (1+Ytbc); 
    end


    %% denoising
    if job.sanlm && si < iter
       cat_sanlm(Ymg,1,3); 

       Ymi  = Ymg; 
       T4th = [ prctile(Ymi(Yp0r(:)==1  & Ymi(:)>0),10), ...
                cat_stat_kmeans(Ymi(Yp0r(:)==2 & Ymg(:)<2.5/3 & Ymg(:)>1.5/3),2), ... % modeling two GM subregions by layer and PVE
                prctile(Ymi(Yp0r(:)==3  & Ymi(:)>0),80) ];
    else
      % local intensity normalization 
      % this function will also enhance the contrast between the tissues!
      Yml = zeros(size(Ymi)); 
      Yml = Yml + ( (Ymi>=Ylab{4}               ) .* (4 + (Ymi - Ylab{4})  ./ max(eps,Ylab{4} - Ylab{3})) ); % scale highest tissue (WM in T1w)
      Yml = Yml + ( (Ymi>=Ylab{3} & Ymi<Ylab{4} ) .* (3 + (Ymi - Ylab{3})  ./ max(eps,Ylab{4} - Ylab{3})) ); % scale second highest tissue (GM in T1w)
      Yml = Yml + ( (Ymi>=Ylab{2} & Ymi<Ylab{3} ) .* (2 + (Ymi - Ylab{2})  ./ max(eps,Ylab{3} - Ylab{2})) ); % scale second highest tissue (GM in T1w)
      Yml = Yml + ( (Ymi>=Ylab{1} & Ymi<Ylab{2} ) .* (1 + (Ymi - Ylab{1})  ./ max(eps,Ylab{2} - Ylab{1})) ); % scale third highest tissue (CSF in T1w)
      Yml = Yml + (  Ymi< Ylab{1}               ) .* (    (Ymi - min(Ymi)) ./ max(eps,Ylab{1} - min(Ymi)) ); % scale background
      Yml = Yml / 4; 

      if job.sanlm 
        cat_sanlm(Yml,1,3);
      end
    end
  end
end
% =========================================================================
function [Ymg,Yml] = suit_amap_LAS2(Ymi,Yp0a,Yc2h,Yg,gth,job)
  Ymio = Ymi; 
  Yp0r = round(Yp0a); 
  %Ytisdistprob = 1; % (1+job.cbdistbias - 2*job.cbdistbias*Yc2h); % not working!
  %%
  verb = 0; 
  Ymi  = Ymio; 
  T4th = [ prctile(Ymi(Yp0r(:)==1  & Ymi(:)>0),10), ...
          cat_stat_kmeans(Ymi(Yp0r(:)==2 & Ymi(:)<2.5/3 & Ymi(:)>0.5/3),4), ... % modeling two GM subregions by layer and PVE
          prctile(Ymi(Yp0r(:)==3  & Ymi(:)>0),80) ];
  T4th([2,5])=[]; % remove PVE helper
  T4vol = @(Yp0,Ymi,T4th) ...
      [ sum( Yp0(:)>0 & Ymi(:)<mean(T4th(1:2)) ), ...
        sum( Yp0(:)>0 & Ymi(:)>mean(T4th(1:2)) & Ymi(:)<mean(T4th(2:3)) ), ...
        sum( Yp0(:)>0 & Ymi(:)>mean(T4th(2:3)) & Ymi(:)<mean(T4th(3:4)) ), ...
        sum( Yp0(:)>0 & Ymi(:)>mean(T4th(3:4)) & Ymi(:)<(T4th(4) + mean(T4th(3:4))) )] ./ nnz(Yp0(:)>0); 
  if verb, T4vol(Yp0a,Ymi,T4th), end

  Ymi=Ymio; 
  iter = 3; 
  for si = 1:iter
    %% local intensity estimation 
    % the global value is more save here ! 

    % CSF
    %Ylab{1}  = median( Ymi( Yp0a(:)>0.5 & Yp0a(:)<1.25 & Ymi(:)<mean(Tth(1:2)) & Ymi(:)>Tth(1)/2 )); 
    Ylab{1}  = T4th(1);
    
    % WM
    % * maybe it is not necessary to handle the WM here as we did this already 
    if 1
      Ylab{4}  = Ymi .* (Yp0a>2.5 & Yp0a<3.5 & Yg>0 & Yg<gth*4); % & Ymi>(T4th(1)*.25+.75*T4th(4)) & Yc2h<0.25);  
      Ylab{4}  = cat_vol_localstat( Ylab{4}, Ylab{4}>0, 1,3); % ### maybe avoid overestimation that we cannot represent in the segmentation
      Ylaba    = cat_vol_median3( Ylab{4}, Ylab{4}>0, Ylab{4}>0 ); % important to avoid bias by outliers!
      Ylab{4}(Yp0a==0) = 1; 
      Ylaba    = cat_vol_approx( Ylaba ); 
      Ylab{4}  = cat_vol_approx( Ylab{4} .* (Ymi ./ Ylaba > .9 & Ymi ./ Ylaba < 1.2) ) * .95; 
    else
      Ylab{4}  = T4th(4); 
    end
    
    % GM 
    cx = .5; cx = 2*[cx 1-cx]; wx = 2 - cx; 
    Ylabg    = Ymi .* (Yp0a>1.5 & Yp0a<2.5 & Yg>0 & Yg<gth*8 & Ymi>mean(T4th(1:2).*cx) & Ymi<mean(T4th(3:4).*wx)); 
    Ylabg    = Ylabg .* (Ylabg>prctile(Ylabg(Ylabg(:)>0),5) & Ylabg<prctile(Ylabg(Ylabg(:)>0),95)); % 
    Ylaba    = cat_vol_approx( Ylabg ); 
    Ylabg    = Ylabg + (Ylabg==0).*cat_vol_approx( Ylabg .* (Ymi ./ Ylaba > .58 & Ymi ./ Ylaba < 1.52) ); 
    %Ylabg    = Ylabg .* Ytisdistprob; % reduced weighting to boundaries 
    Ylabg    = cat_vol_smooth3X(Ylabg,2); 
  
    % lower GM
    Ylab{2}  = ( Ymi .* (Yp0a>1.5 & Yp0a<2.5 & Ymi<Ylabg*.99 & Ymi>mean(T4th(1:2)) & Yg<gth*4)); 
    Ylaba    = cat_vol_approx( Ylaba ); 
    Ylab{2}  = Ylab{2} + (Ylab{2}==0).*cat_vol_approx( Ylab{2} .* (Ymi ./ Ylaba > .8 & Ymi ./ Ylaba < 1.2) ); 
    %Ylab{2}  = Ylab{2} .* Ytisdistprob; 
    Ylab{2}  = cat_vol_smooth3X(Ylab{2},2); 
    
    % higher GM
    Ylab{3}  = ( Ymi .* (Yp0a>1.5 & Yp0a<2.5 & Ymi>Ylabg/.99 & Ymi<mean(T4th(3:4)) & Yg<gth*4)); 
    Ylaba    = cat_vol_approx( Ylaba ); 
    Ylab{3}  = Ylab{3} + (Ylab{3}==0).*cat_vol_approx( Ylab{3} .* (Ymi ./ Ylaba > .8 & Ymi ./ Ylaba < 1.2) ); 
    %Ylab{3}  = Ylab{3} .* Ytisdistprob;
    Ylab{3}  = cat_vol_smooth3X(Ylab{3},2); 
    %   clear Ylabg Ylaba; 

    % local intensity normalization 
    maxi = 4;
    for i=1:maxi
      Ymg = zeros(size(Ymi)); Ylab23 = (Ylab{2}+Ylab{3})/2; 
      Ymg = Ymg + ( (Ymi>=Ylab{4}               ) .* (3 + (Ymi - Ylab{4})  ./ max(eps,Ylab{4} - Ylab23)) ); % scale highest tissue (WM in T1w)
      Ymg = Ymg + ( (Ymi>=Ylab23  & Ymi<Ylab{4} ) .* (2 + (Ymi - Ylab23)   ./ max(eps,Ylab{4} - Ylab23)) ); % scale second highest tissue (GM in T1w)
      Ymg = Ymg + ( (Ymi>=Ylab{1} & Ymi<Ylab23  ) .* (1 + (Ymi - Ylab{1})  ./ max(eps,Ylab23  - Ylab{1})) ); % scale third highest tissue (CSF in T1w)
      Ymg = Ymg + (  Ymi< Ylab{1}               ) .* (    (Ymi - min(Ymi)) ./ max(eps,Ylab{1} - min(Ymi)) ); % scale background
      Ymg = Ymg / 3; 

      tv = T4vol(Yp0a,Ymg,T4th); 
      if i < maxi
        Ylabd   = (Ylab{3} - Ylab{2}); 
        if tv(3) / sum(tv(2:3)) > .5
          Ylab{2} = Ylab{2} + Ylabd * .05;
        else
          Ylab{3} = Ylab{3} - Ylabd * .05;
        end
      end
    end
   
%%
    if 0 % si < 3
      %% magic corretion ... to slow
      Ywd  = zeros(size(Ymg)); for i=.1:.2:.9, Ywd = Ywd + 1/5 * cat_vbdist( single(Ymg*3 > 1.75 + i/2 ) , Yp0a>0 & Ymg>1.5/3, vx_vol); end
      Ycd  = zeros(size(Ymg)); for i=.1:.2:.9, Ycd = Ycd + 1/5 * cat_vbdist( single(Ymg*3 < 1.75 + i/2 ) , Yp0a>0 & Ymg<2.5/3, vx_vol); end
      Yth1 = cat_vol_approx( (Ywd + Ycd) .* (Yp0r==2 & Ywd<1e10 & Ycd<1e10) ); 
      Ytbc = cat_vol_smooth3X( (Yp0r==2 & Ycd<1e10 & Ywd<1e10) .* max(0,Ycd - 1.5*Yth1),4); %max(0,Yth1 - median(Yth1(Yp0r(:)==2)))) ); 
      Ymg  = Ymg .* (1+Ytbc); 
      Ytbc = cat_vol_smooth3X( (Yp0r==2 & Ycd<1e10 & Ywd<1e10) .* max(0,Ywd - 1.5*Yth1),4); %max(0,Yth1 - median(Yth1(Yp0r(:)==2)))) ); 
      Ymg  = Ymg .* (1+Ytbc); 
    end


    %% denoising
    if job.sanlm && si < iter
       Ymgs = Ymg .* (Yp0a>0); cat_sanlm(Ymgs,1,3); 
       Ymg(Yp0a>0) = Ymg(Yp0a>0).*.5 + .5.*Ymgs(Yp0a>0); clear Ymgs;
       Ymi  = Ymg; 
       Ymi  = Ymi.^(1./(Ylab{4}./(Ylab23/mean(Ylab23(:)))).^.5); %Ymi.^(1./(Ylab{4}*.5 + 0.5.*Ylab23/mean(Ylab23(:))).^2);
     
       T4th = [ prctile(Ymi(Yp0r(:)==1  & Ymi(:)>0),10), ...
                cat_stat_kmeans(Ymi(Yp0r(:)==2 & Ymg(:)<2.5/3 & Ymg(:)>1.5/3),4), ... % modeling two GM subregions by layer and PVE
                prctile(Ymi(Yp0r(:)==3  & Ymi(:)>0),80) ];
       T4th([2,5])=[]; % remove PVE helper

       if verb, T4vol(Yp0a,Ymi,T4th), end
    else
      % local intensity normalization 
      % this function will also enhance the contrast between the tissues!
      for i = 1:maxi
        Yml = zeros(size(Ymi)); 
        Yml = Yml + ( (Ymi>=Ylab{4}               ) .* (4 + (Ymi - Ylab{4})  ./ max(eps,Ylab{4} - Ylab{3})) ); % scale highest tissue (WM in T1w)
        Yml = Yml + ( (Ymi>=Ylab{3} & Ymi<Ylab{4} ) .* (3 + (Ymi - Ylab{3})  ./ max(eps,Ylab{4} - Ylab{3})) ); % scale second highest tissue (GM in T1w)
        Yml = Yml + ( (Ymi>=Ylab{2} & Ymi<Ylab{3} ) .* (2 + (Ymi - Ylab{2})  ./ max(eps,Ylab{3} - Ylab{2})) ); % scale second highest tissue (GM in T1w)
        Yml = Yml + ( (Ymi>=Ylab{1} & Ymi<Ylab{2} ) .* (1 + (Ymi - Ylab{1})  ./ max(eps,Ylab{2} - Ylab{1})) ); % scale third highest tissue (CSF in T1w)
        Yml = Yml + (  Ymi< Ylab{1}               ) .* (    (Ymi - min(Ymi)) ./ max(eps,Ylab{1} - min(Ymi)) ); % scale background
        Yml = Yml / 4; 
  
        tv = T4vol(Yp0a,Ymg,T4th); 
        if i < maxi
          Ylabd   = (Ylab{3} - Ylab{2}); 
          if tv(3) / sum(tv(2:3)) > .5
            Ylab{2} = Ylab{2} + Ylabd * .05;
          else
            Ylab{3} = Ylab{3} - Ylabd * .05;
          end
        end
      end

      if job.sanlm 
        Ymls = Yml .* (Yp0a>0); cat_sanlm(Ymls,1,3); 
        Yml(Yp0a>0) = Yml(Yp0a>0).*.25 + .75.*Ymls(Yp0a>0); clear Ymls;
      end
    end
  end
end
% =========================================================================
function Ymi = suit_amap_biascorr(Ym,Yp0a,Yg,Yc2h,Ycbd,gth,job,vx_vol)
%% quick bias correction 
  Ymo = Ym; 

  Ym = Ym ./ cat_vol_approx(smooth3(Ym) ./ max(eps,Yp0a/3) .* (Ycbd>0 & Yp0a>1.5)) ;  
  Ym = Ym ./ prctile(Ym(Yp0a(:)>2.5),60) * prctile(Ymo(Yp0a(:)>2.5),60);

  Ymo = Ym; 

  [Ymi,redR]  = cat_vol_resize(Ym  ,'reduceV',vx_vol,.9,32,'meanm');   % CSF thr. (minimum to avoid PVE)
  Yp0a       = cat_vol_resize(Yp0a,'reduceV',vx_vol,.9,32,'meanm');   % CSF thr. (minimum to avoid PVE)
  Yg         = cat_vol_resize(Yg  ,'reduceV',vx_vol,.9,32,'meanm');   % CSF thr. (minimum to avoid PVE)
  Yc2h       = cat_vol_resize(Yc2h,'reduceV',vx_vol,.9,32,'meanm');   % CSF thr. (minimum to avoid PVE)
  Ycbd       = cat_vol_resize(Ycbd,'reduceV',vx_vol,.9,32,'meanm');   % CSF thr. (minimum to avoid PVE)

  %%
  % correct a slight bias in more WM in and more CSF oudwards
  Ytisdistprob = (1+job.cbdistbias - job.cbdistbias*Yc2h); 
  for i = 1:3
    Ybs = cat_vol_morph( cat_vol_morph( Yp0a>2.75 ,'lo' ,2,vx_vol ),'lc',2);
    % GM
    Ylabgw1  = Ymi .* (Yp0a>1.5 & Yp0a<2.95 & Yg<gth*8 & Ycbd<0 & ~Ybs) .* Ytisdistprob; 
    Ylabgw1a = cat_vol_approx( Ylabgw1 );  
    Ylabgw1  = cat_vol_approx( Ylabgw1 .* (Ymi ./ Ylabgw1a > .8 & Ymi ./ Ylabgw1a < 1.3) ); clear Ylabgw1a;
    % WM
    Ylabgw2  = Ymi .* ( (Yp0a>2.15 & Yp0a<3.5 & Yg<gth*4 & Ycbd<-3) | Ybs )  .* Ytisdistprob; 
    Ylabgw2  = max( Ylabgw2, cat_vol_localstat(Ymi,Yp0a>1.5 & Yp0a<2.5 & Yg<1,3/mean(vx_vol),1) .* Ytisdistprob * 1.1); %1.1= noise 
    Ylabgw2a = cat_vol_approx( Ylabgw2 ); 
    Ylabgw2  = cat_vol_approx( Ylabgw2 .* (Ymi ./ Ylabgw2a > .8 & Ymi ./ Ylabgw2a < 1.3) ); clear Ylabgw2a;
    %% WM cerebrum
    Ycbr     = cat_vol_morph( (Ymi<1.2 & Yg<gth*2 & Ymi>mean(Ymi(Yp0a>2.5))*.95 & Ycbd>2), 'ldo'); 
    Ylabgw3  = cat_vol_localstat( Ymi , Yp0a==0  & Ycbr , 2/mean(vx_vol), 3); 
    Ylabgw3a = cat_vol_approx( Ylabgw3 );
    Ylabgw3  = cat_vol_approx( Ylabgw3 .* (Ymi ./ Ylabgw3a > .8 & Ymi ./ Ylabgw3a < 1.3) ); clear Ylabgw3a;
    %% mix
    Ylabgw1  = Ylabgw1 ./ cat_stat_nanmedian(Ylabgw1(Yp0a(:)>1.5)) .* cat_stat_nanmedian(Ylabgw2(Yp0a(:)>1.5));
    Ylabgw   = Ylabgw1 .* (1-Yc2h) + Ylabgw2 .* Yc2h; 
    Ylabgw   = Ylabgw .* cat_vol_smooth3X(Yp0a/3,2) + Ylabgw3 .* (1-cat_vol_smooth3X(Yp0a/3,2));
    Ylabw    = Ylabgw ./ cat_stat_nanmedian(Ylabgw(Yp0a(:)>2.5));

    Ymi = Ymi ./ cat_vol_smooth3X(Ylabw,4); 
    Ymi = Ymi ./ prctile( Ymi(Yp0a(:)>2.5) , 60);
  end
  
  %%
  Ymi  = cat_vol_resize(Ymi,'dereduceV',redR); 
  Yw   = cat_vol_approx( Ymo ./ Ymi , 'rec'); 
  Ymi  = Ymo ./ Yw; 

end
% =========================================================================
function [Ycb,Yml,Ycbd,Yc2h] = suit_amap_cleanup_cerebellum(Yml,Yp0s,Ycbd,Yg,Ybv,vx_vol,init)
%suit_amap_cleanup_cerebellum. Extract cerebellar structure.
% - We use here a quite hard approach, i.e. focusing on GM and WM but 
%   ignoring CSF, to avoid interaction with other tissues that can 
%   trouble the registration routine

  Ymlo = Yml; % Yml=Ym; init=1; Ycbd0=Ycbd; ii=1;

  [Yml,redR] = cat_vol_resize(Yml ,'reduceV',vx_vol,.8,32,'meanm');   % CSF thr. (minimum to avoid PVE)
  Yp0s       = cat_vol_resize(Yp0s,'reduceV',vx_vol,.8,32,'meanm');   % CSF thr. (minimum to avoid PVE)
  Yg         = cat_vol_resize(Yg  ,'reduceV',vx_vol,.8,32,'meanm');   % CSF thr. (minimum to avoid PVE)
  Ycbd       = cat_vol_resize(Ycbd,'reduceV',vx_vol,.8,32,'meanm');   % CSF thr. (minimum to avoid PVE)
  
  if init
    csfth = prctile( Yml(Yml(:)<.6 & Ycbd(:)<5) , 5); 
    gmth  = cat_stat_nanmedian( Yml(Yg(:)>1/9 & Yml(:)<.9 & Yml(:)>csfth*1.5 & Ycbd(:)<-2  & Ycbd(:)>-20 & Yml(:)<.95) ); 
  else
    csfth = 1/4; 
    gmth  = 3/4; 
  end
  % remove extremly high intensities that point to skull/blood vessels
  % especially blood vessels have to be removed rather then aligned to the 
  % skull as they would otherwise take too much in the region growing
  Yml(Yml > 1.2 & Ycbd>-5) = nan; 
%  Yml(smooth3(Yg ./ Yml) > 2 & Yml>.8 & Ycbd>-5) = nan; 

  %% define cerebellar and non-cerebellare areas for region-growing
  for ii = 1:3
    % second loop for small refinements
    Ybr = Ycbd .* max( Yml , 1.5 - Yml) > 3  |  smooth3( Ycbd>-1 & Yml<1.75/3 & Yml>1.25/3)>.5; 
    Ybr = cat_vol_morph(Ybr & Yml>(csfth*.25+.75*gmth) | cat_vol_morph(-Ycbd ./ Yml < -4,'ldo',5,redR.vx_volr) ,'dc',1,redR.vx_volr);
    Ybr = Ybr | (Yml > csfth & Yml<(csfth*.25+.75*gmth) & Ycbd>-3 & Yg<.2); 
    Ycb = -Ycbd .* Yml > 0 & Yml<-Ycbd/2; 
    %
    Ybr = cat_vol_morph( smooth3(Ybr)>.7 + 0.1*(3-ii) , 'de', 3/ii, redR.vx_volr); 
    if ii < 3
      Ycb = cat_vol_morph( smooth3(Ycb)>.5 + 0.1*(3-ii) , 'do', 2-ii, redR.vx_volr); 
      Ycb = cat_vol_morph( Ycb , 'ldc', 5, redR.vx_volr);  
    end
    %Ycb = smooth3(Ycb)>.95;  
    % run region growing with different stop criteria
    Ya  = single(Ybr + 2*Ycb); %Ya(Yml<(csfth*.25+0.75*gmth)) = nan; 
    ths = [-.01 0.01 0.05 .1 1] - 0.05*(ii>1); 
    for thi = 1:numel(ths)
      [Ya,Yd] = cat_vol_downcut(Ya,Yml + .01*Ycbd,ths(thi));
      Ya(Ya==2 & cat_vol_smooth3X(Ya)<1.125 & Yd>100) = 0; 
      Ya(Ya==2 & Ybv) = 1.5; 
      Ya(smooth3(Ya==1)<.35 & Ya==1) = 0; % 1.5
      Ya(smooth3(Ya==2)<.35 & Ya==2) = 0; % 1.5 
      %Ya(cat_vol_morph(Ya>1.5,'ldc',1)) = 2;
      Ya(isnan(Yml))=nan; 
    end
    Ya(isnan(Ya))=1.25; 
    %%
    Ynbv  = cat_vol_morph( smooth3( ~isnan(Yml) )>.75 ,'e',2) ;
    Ycb   = cat_vol_morph( cat_vol_morph( smooth3(Ya)>1.5 & Yml>(csfth*.5+.5*gmth) & ~isnan(Yml)  ...
      ... & (Yml<.95 | smooth3(Ycb)>.75) ...
      ,'ldo',.5,redR.vx_volr),'ldc',2)  & Ynbv;
    Ycb(Yp0s>2.5) = 1; 
    Ycb( isnan(Yml) & cat_vol_smooth3X(Ycb,2)<.75 & Yp0s<2.5 ) = 0;  
    Ycb(Ybv) = 0; 
    Ycb   = min(1,cat_vol_smooth3X( Ycb ,2) .* 1.2 .* (Ynbv | Yml<.9 | smooth3(Ycb)>.7 | Yp0s>2.5)); 
%%
    Ycbd = cat_vol_smooth3X(cat_vbdist(single(Ycb>0.5) , true(size(Ycb)), redR.vx_volr) - ... 
                            cat_vbdist(single(Ycb<0.5) , true(size(Ycb)), redR.vx_volr),2);    
  
  
    %% cleanup
    if 0
      Ycb = ( (Ycb .* Yml)>(csfth*.5+.5*gmth)  & Yml<(1*.75+0.25*gmth) ) | (Ycbd<-3) | ...
            cat_vol_morph(cat_vol_morph( Ycb>0 & Yml>mean([1,gmth]) & Yml<1.1,'ldo',2,redR.vx_volr),'d'); % brainstem
      % closing
      Ycb = smooth3( cat_vol_morph( cat_vol_morph( Ycb ,'do',0.5,redR.vx_volr) ,'ldc' , 1.5,redR.vx_volr)) > .25 & Ya>1.5;
      Ycb = Ycb .* (Yml<.9 | smooth3(Ycb)>.9); 
    end
  end
  %%
  if 0
    % high intensity blood vessels
    % - only relevant close to brainstam
    % - the sistance was not really working ... maybe add smoothing 
    % - to what do you correct ? 
    % - maybe just better to leave it his way?
    [~,YI] = cat_vbdist(single( Yml<1.05 ), Ycb); 
    Yml = Yml(YI); clear YI;  
  end

  if 0
  % CAT/SPM cleanup alternative
    cls = 3;
    
    Yprobt = zeros([size(Yml),5],'uint8'); 
    for ci=1:cls, Yprobt(:,:,:,ci) = uint8( max(0,min(255, Yp0toC(Yml .* cls .* Ycb,ci) * 255))); end
    Yprobt(:,:,:,cls) = Yprobt(:,:,:,cls) + uint8( max(0,min(255, Yp0toC(Yml .* cls .* Ycb,cls+1) * 255)));
    Yprobt(:,:,:,1) = 255.*Ycb - sum(Yprobt,4);
    Yprobt = cat_main_clean_gwc1639( Yprobt , 1, 1); 
    Yml2   = zeros( size(Yml) , 'single'); 
    for ci=1:cls, Yml2 = Yml2 + single(Yprobt(:,:,:,ci)) / 255 * ci/cls; end
    Yml2 = cat_vol_downcut(single(cat_vol_smooth3X(Yml2,4)>1/4), Yml + .01*max(0,Ycbd2),-0.005) & Ycbd2<3 & Yml>.125;
    Yml2 = cat_vol_morph( cat_vol_smooth3X(Yml2,2) > .5, 'ldc', 2).^4; 
    Ycb  = min(1,cat_vol_smooth3X(Yml2,1)).^4; 
  end

  %% estimate core-hull percentage position
  Ycc  = smooth3(cat_vol_morph( smooth3(Yml)>mean([gmth,1]) & Ycb,'ldo',2,redR.vx_volr));
  Yhh  = smooth3(cat_vol_morph( smooth3(Yml)>mean(gmth*.5 + 0.5*csfth) & Ycb ,'ldc',4,redR.vx_volr)); 
  Ycd  = zeros(size(Ycc),'single'); Yhd = Ycd; ni = 5;
  for i = 1:ni
    Ycd = Ycd + (1/ni) * cat_vbdist(single( Ycc>(i-.5)/ni), Yhh>0, redR.vx_volr);
    Yhd = Yhd + (1/ni) * cat_vbdist(single( Yhh<(i-.5)/ni), Ycc<1, redR.vx_volr);
  end
  Yc2h = max(Yhd ./ max(eps,Ycd + Yhd), Ycc ); 
 % clear Ycc Yhh Ycd Yhd;
    
  %%
  Ycb  = cat_vol_resize(single(Ycb), 'dereduceV', redR) > .5; 
  Ycbd = cat_vol_resize(Ycbd, 'dereduceV', redR); 
  Yc2h = cat_vol_resize(Yc2h, 'dereduceV', redR); 

  % general WM limit
  if init
    Yml = Ymlo; 
  else
    Yml = min(Ymlo,0.9 + 0.20*Yc2h) .* Ycb; 
  end
end
% =========================================================================
function Ymi = suit_amap_quickscale(Ym,Yp0s,vx_vol,P,fi)
%suit_amap_quickscale. Quick bias correction (if original image is used) & intensity normalization.

  Tth = [prctile(Ym(round(Yp0s(:))==1),10), prctile(Ym(round(Yp0s(:))==2),50), prctile(Ym(round(Yp0s(:))==3),90)];

  % estimate bias 
  Ywme = cat_vol_morph(Yp0s>2.5,'e');     % use save WM for estimation 
  Ywme(smooth3(Ywme)<.5) = 0; 
  Ygm  = Yp0s>1.5 & Yp0s<2.5 & Ym>min(Tth) & Ym<max(Tth);    % main GM area with CSF and WM PVE
  Ygm  = cat_vol_morph(Ygm,'dc',1); 
  Ymw  = cat_vol_localstat(Ym,Ywme,3,1);  % denoised WM
  Ymg  = cat_vol_localstat(Ym,Ygm,1,1);   % denoised GM
  bias = std(Ymw(Ywme(:))) + std(Ymg(Ygm(:))); % bias estimate (~CJV)
  Yex  = cat_vol_smooth3X(Yp0s<.5,16)>.9; % brainmask
  Yg   = cat_vol_grad(Ym) ./ Ym; 

  Yp0s(~Ywme) = min(2,Yp0s(~Ywme));
  
  Ymo = Ym; 
  %% initial bias correction 
  res = 1.5; % run on lower resolution for speed but also lower noise 
  for bci = 1:2
    if bci>1, Ym = Ymi; else, Ym = Ymo; end
    if ~exist(P.mm{fi},'file') && bias>.02
    %  important for later fine correction especially in case of strong bias and SPM masking issues 
      [Ymr,redR]    = cat_vol_resize(cat_vol_median3(Ym.*(Yp0s>1.25),Yp0s>1.25),'reduceV',vx_vol,res,32,'meanm');   % CSF thr. (minimum to avoid PVE)
      Yp0sr = cat_vol_resize(Yp0s,'reduceV',vx_vol,res,32,'meanm');   % CSF thr. (minimum to avoid PVE)
      Ywmer = cat_vol_resize(Ywme,'reduceV',vx_vol,res,32,'meanm') > .5;  
      Ygmr  = cat_vol_resize(Ygm,'reduceV',vx_vol,res,32,'meanm') > .5;  
      Yexr  = cat_vol_resize(Yex,'reduceV',vx_vol,res,32,'meanm') > .5;  
      Ygr   = cat_vol_resize(Yg,'reduceV',vx_vol,res,32,'meanm') > .5;  
      % get average WM and GM intensies
      Ywr = cat_vol_approx(max(Yexr * prctile(Ymr(Ywmer),95), ...
        cat_vol_localstat(Ymr,Ymr>median(Ymr(:)) & Yp0sr>2.5,2,2 + ( Tth(2) < Tth(3) )))); 
      Ygr = cat_vol_approx(max(Yexr * prctile(Ymr(Ygmr),50), ...
        cat_vol_localstat(Ymr,Ymr>median(Ymr(:)) & Yp0sr>1.5 & Yp0sr<2.5 & ...
        Ymr>min(Tth)*1.5 & Ymr<max(Tth)*.9 & Ygr>median(Ygr(round(Yp0sr(:))==2)),2,1))); 
      %% adapt and mix both tissue intensities  
      Ywr = (Ywr*.5 + .5*(Ygr / prctile(Ygr(Ywmer),50)) * prctile(Ywr(Ywmer),50) ); 
      % go back to original resolution
      Yw  = cat_vol_resize(Ywr,'dereduceV',redR); clear Ywr Ygr Yexr Ywmer Ywmer; 
      Ymi = real( (Ym ./ Yw) .^ max(eps,Yw./mean(Yw(Ywme(:)))) ); clear Yw;
      Ymw = cat_vol_localstat(Ymi,Ywme,1,1); 
      % final test if corretion was good otherwise take previous version
      biasi = std(Ymw(Ywme(:))) + std(Ymg(Ygm(:)));
      if bias < biasi, Ymi = Ym; break; end
    else
      Ymi = Ym;
    end
  
    % intensity normalization
    Ymi  = Ymi / prctile(Ymi(Ywme(:)),80);

    % further normalization ???
    Ymi = suit_quadintnorm(Ymi,Yp0s);
  end

  % normalization 
  Ymi = suit_quadintnorm(Ymi,Yp0s);

  %%

  %Ymi = max(-10,min(Ymi,10));
end
% =========================================================================
function [Ym2, coeff] = suit_quadintnorm(Ym,Yp0s)
%suit_quadintnorm. quadratic intensity correction.

  % define normalization intensities based on the segmentation 
  if 0
    x = [ ...
      prctile(Ym(round(Yp0s(:))==1),10);
      prctile(Ym(round(Yp0s(:))==2),50);
      prctile(Ym(round(Yp0s(:))==3),90);
      prctile(Ym(round(Yp0s(:))==3),99);
    ]; 
    y = [1; 2; 3; 4] / 3;
  else
    x = [ ...
      prctile(Ym(Ym(:)~=0),.5);
      prctile(Ym(round(Yp0s(:))==1),10);
      prctile(Ym(round(Yp0s(:))==2),50);
      prctile(Ym(round(Yp0s(:))==3),90);
      prctile(Ym(round(Yp0s(:))==3),99);
    ]; 
    y = [0; 1; 2; 3; 4] / 3;
  end

  % Build the Vandermonde matrix
  A = [x.^2, x, ones(numel(y),1)];
  coeff = A \ y;
  
  % apply to image
  Ym2 = coeff(1)*Ym.^2 + coeff(2)*Ym + coeff(3);
end
% =========================================================================
function [Yp0s,Ym,Vmsk,Vgm,vx_vol] = suit_amap_load_segments(P,fi)
%suit_amap_load_segments. Load SPM data.

  % load header
  Vgm  = spm_vol(P.gm{fi});
  Vwm  = spm_vol(P.wm{fi}); 
  Vmsk = spm_vol(P.msk{fi}); 
  Vmsk.fname = ''; % cropping properties
  if exist(P.mm{fi},'file')
    Vm = spm_vol(P.mm{fi});
  else
    Vm = spm_vol(P.m{fi});
  end

  % load data into subspace (via imcalc) and prepare CSF
  [~,Ygm] = cat_vol_imcalc(Vgm, Vmsk, 'i1'); Ygm = single(Ygm);
  [~,Ywm] = cat_vol_imcalc(Vwm, Vmsk, 'i1'); Ywm = single(Ywm); 
  Ycsf    = cat_vol_morph( (Ygm + Ywm) > 0,'ldc',4) - max(0,Ygm + Ywm);
  
  % final label map
  Yp0s    = Ycsf + 2*Ygm + 3*Ywm; 
  Ycb     = cat_vol_morph(Yp0s>.5,'ldo'); Yp0s(~Ycb) = 0; % quick cb stripping
  Yp0s    = cat_vol_median3(Yp0s,Yp0s>0); % denoising
  Ym      = single(spm_read_vols(Vm)); 

  % estimate voxel size
  vx_vol = sqrt(sum(Vm.mat(1:3,1:3).^2));
end
% =========================================================================
function [Tth,suit] = suit_amap_eval(Ym,Ymi,Yp0a,Yp0s,vx_vol,suit,fi)
%suit_amap_eval. 

  Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0 - c));

  % threshold for the tissues - use extrem values for CSF and WM 
  Yp0r = uint8(round(Yp0a));
  Tth  = [ prctile(Ymi(Yp0r(:)==1  & Ymi(:)>0),10), ...
           cat_stat_nanmedian(Ymi(Yp0r(:)==2  & Ymi(:)>0)), ...
           prctile(Ymi(Yp0r(:)==3  & Ymi(:)>0),80) ];
  % new 4 class model
  T4th = [ prctile(Ymi(Yp0r(:)==1  & Ymi(:)>0),10), ...
           cat_stat_kmeans(Ymi(Yp0r(:)==2),2), ... % modeling two GM subregions by layer and PVE
           prctile(Ymi(Yp0r(:)==3  & Ymi(:)>0),80) ];
  clear Ylabw Ylabgw Ylabgw1 Ylabgw2 Ylabgw3; 

  % test for skull-stripping problems (e.g., AHEAD dataset)
  suit(fi).stripped = nnz(Ym(Yp0s(:)>0)==0) ./ nnz(Yp0s(:)>0);
  if suit(fi).stripped > .1
    cat_io_cprintf('warn',...
      sprintf('Warning: Large empty areas (%0.2f%%%%) indicate skull-stripping problems! ', 100*suit(fi).stripped));
  end

  % basic SPM volume values
  suit(fi).vol_spm_BSCBV   = sum(min(1,Yp0s(:))) * prod(vx_vol) / 1000;
  suit(fi).vol_spm_abs_CGW = [sum(Yp0toC(Yp0s(:),1)) sum(Yp0toC(Yp0s(:),2)) sum(Yp0toC(Yp0s(:),3))] * prod(vx_vol) / 1000;
  suit(fi).vol_spm_rel_CGW = suit(fi).vol_spm_abs_CGW ./ suit(fi).vol_spm_BSCBV;

  % basic intensity values
  suit(fi).int_abs_CGW  = Tth; 
  suit(fi).int_rel_CGW  = Tth ./ max(Tth); 
  suit(fi).int_abs_CGGW = T4th; 
  suit(fi).int_rel_CGGW = T4th ./ max(Tth); 
end
% =========================================================================
function [Ycba,Ybscb,cba,cbas,suit] = suit_amap_parcelate(Yml,Yp0,Yp0s,Ycb,Yg,gth,vx_vol,suit,fi)
%suit_amap_parcelate. Define major cerregions.

  % basic regions
  cba = {
    1 'SUITmask'
    2 'brainstem'
    3 'upper left peduncle'
    4 'upper right peduncle'
    5 'lower left peduncle'
    6 'lower right peduncle'
    7 'cerebellar hull'
    8 'cerebellar core'
  };

  % helping function
  Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0 - c));
     
  % basic tissues
  Ywms = Yp0toC(smooth3(Yp0s),3); 
  Ygms = Yp0toC(smooth3(Yp0s),2);

  %% basic regions: brainstem, cerebellum, brain and head to describe up and down
  Ycb0 = cat_vol_morph( cat_vol_morph( Ygms>.5 , 'lo',1),'ldc',8,vx_vol);                % cerebellum
  Ybs  = cat_vol_morph( Ywms>.5 & ~Ycb0, 'ldo',1,vx_vol);                                % brainstem
  Ycbs = single(cat_vol_morph( Ycb0 .* (Yml>.625) + (Ygms>.5 & Yml<.5),'lo')) + ...
       2*single(cat_vol_morph( Ybs .* (Yml>.75),'lo')); Ycbs(~Ycb) = nan;  
  Ycbs = cat_vol_downcut(Ycbs,Yml,0);
  Ycb0 = cat_vol_morph(Ycbs>0 & Ycbs<1.5 & Ycb,'lc');
  Ybs  = cat_vol_morph(Ycbs>1.5,'lc');
  %
  Ybr  = cat_vol_morph( Yml>.6 & Yml<1.2 & Yg<gth*2 & ~Ywms & ~Ygms,'lo');               % brain
  Ybrd = cat_vbdist(single(Ybr),true(size(Yml)),vx_vol);                                 % brain distance
  Ybsd = cat_vbdist(single(Ybs),Yp0>0,vx_vol);                                           % brainstem distance
  Ycbd = cat_vbdist(single(Ycb0),Yp0>0,vx_vol);                                          % brainstem distance
  Yhd  = cat_vol_morph( Yml>.3 & Yml<1.3 & ~Ywms & ~Ygms & Ybrd>20 ,'ldo',3);            % head
  Ytd  = cat_vol_smooth3X(cat_vol_approx(cat_vol_approx(Yhd + 2*Ybr) .* Ybs),32);        % top-down (cerebrum-head)
%  clear Yhd Ybr Ybrd Ygms; 
  
  % core & hull
  Ycore  = cat_vol_morph( smooth3(Yml) > 3.5/4 & Ycb0, 'ldo', 2);
  Yhull  = cat_vol_morph( smooth3(Yml) > 1.5/4 & Ycb0, 'ldc', 2); 
  Ycba   = zeros(size(Ycb)); 
  Ycba(Ycb & Yhull) = cba{cat_io_contains(cba(:,2),'hull'),1};
  Ycba(Ycb & Ycore) = cba{cat_io_contains(cba(:,2),'core'),1};
  Ycba(Ycb & Ybs)   = cba{cat_io_contains(cba(:,2),'brainstem'),1};
 

  %% peduncles: 
  % brainstem-cerebrum map
  Ybscb = single( 1.5 + ...
    0.5*(Ycbd>5 | cat_vol_morph(cat_vol_morph(Ybs  & Ywms>.5, 'de',2.9,vx_vol))) - ...
    0.5*(Ybsd>5 | cat_vol_morph(cat_vol_morph(Ycb0 & Ywms>.5, 'de',2.9,vx_vol)))); Ybscb(Yp0s<=2 | Ycbd>6 | Ybsd>6) = nan;
  Ybscb = cat_vol_laplace3R(Ybscb, Ybscb==1.5,.001); 
%######## use V.mat orientation 
  Ycbc  = cat_vol_morph(Ybscb>1.2 & Ybscb<1.8 & Ybscb~=2 & Yp0s>1.5 & Ywms>.5, 'l', [4, .2]);
  for i = 0:2:4 % backups
    if max(Ycbc(:))<2, Ycbc  = cat_vol_morph(Ybscb>1.5+i/10 & Ybscb<1.6+i/10 & Ybscb~=2 & Yp0s>2 & Ywms>.5, 'l', [4, .2]); end
    if max(Ycbc(:))<2, Ycbc  = cat_vol_morph(Ybscb>1.5-i/10 & Ybscb<1.4-i/10 & Ybscb~=2 & Yp0s>2 & Ywms>.5, 'l', [4, .2]); end
  end
  %%
  if nnz(Ycbc(:)==1)==0 || nnz(Ycbc(:)==2)==0
    cat_io_cprintf('err','Error finding peduncles');
  end

  Ycba(Ycb & Ycba==0 & Ybscb>1.5) = cba{cat_io_contains(cba(:,2),'SUITmask'),1};
  
  clear Ywms; 

  % seed definition for a upper and lower region of both peduncles
  Ycbc1 = Ytd; Ycbc1(Ycbc~=1) = nan; 
  Ycbc2 = Ytd; Ycbc2(Ycbc~=2) = nan; 
  Ycbc1 = Ycbc1 - cat_stat_nanmedian(Ycbc1(:)); 
  Ycbc2 = Ycbc2 - cat_stat_nanmedian(Ycbc2(:)); 
  Ycbc1 = Ycbc1 / cat_stat_nanstd(Ycbc1(:)); 
  Ycbc2 = Ycbc2 / cat_stat_nanstd(Ycbc2(:)); 
  clear Ycbc

  %% update atlas 
  Ycba(Ycb & Ycbc1>=0) = cba{cat_io_contains(cba(:,2),'upper left peduncle'),1};
  Ycba(Ycb & Ycbc1<0)  = cba{cat_io_contains(cba(:,2),'upper right peduncle'),1};
  Ycba(Ycb & Ycbc2>=0) = cba{cat_io_contains(cba(:,2),'lower left peduncle'),1};
  Ycba(Ycb & Ycbc2<0)  = cba{cat_io_contains(cba(:,2),'lower right peduncle'),1};
  Ycba(Ycb & Ycba==0)  = cba{cat_io_contains(cba(:,2),'SUITmask'),1};
  clear Ycbc1 cbc2 th; 
 
  cbas  = 'Cerebelar ROI atlas (';
  for ci = 1:size(cba,1)
    cbas = sprintf('%d-%s,', cba{ci,1}, cba{ci,2}); 
    suit(fi).vols.(strrep(cba{ci,2},' ','_')) = sum(Ycba(:) == cba{ci,1}) .* prod(vx_vol); 
  end
  cbas = [cbas(1:end-1) ')']; 
end
% =========================================================================
function [Yth0,Ypp0,suit] = suit_amap_cb0thickness(Ycba,cbac,vx_vol,suit,fi)
%%  (4) cerebellar measures for level 0 (core vs. hull)
%  ====================================================================

%##### eikonal !    
  Ycore  = Ycba==cbac{cat_io_contains(cbac(:,2),'core'),1};
  Yhull  = Ycba>=cbac{cat_io_contains(cbac(:,2),'hull'),1};
  Ydcore = cat_vbdist( single( Ycore) ,  Yhull , vx_vol); 
  Ydhull = cat_vbdist( single(~Yhull) , ~Ycore , vx_vol);
  Yth0   = cat_vol_approx( (Ydcore + Ydhull) .* (~Ycore & Yhull & Ydcore<1e10 & Ydhull<1e10 ) ); 
  Ypp0   = min(1,((Yth0-Ydcore).*(Ycba>0)*.5 + .5*Ydhull) ./ max(eps,Yth0)); % mid surface visual
  suit(fi).thickness_level0 = median(Yth0(~Yhull(:) & ~Ycore(:) )); 
end
% =========================================================================
function Ycoord = suit_amap_get_transversal(Yml,Ybscb,Ycba,cbac)
%suit_amap_get_transversal. x-axis between left to right peduncles

  Ycbc   = 1*(Ycba==cbac{cat_io_contains(cbac(:,2),'upper left peduncle'),1}) + ...
           1*(Ycba==cbac{cat_io_contains(cbac(:,2),'lower left peduncle'),1}) + ...
           2*(Ycba==cbac{cat_io_contains(cbac(:,2),'upper right peduncle'),1}) + ...
           2*(Ycba==cbac{cat_io_contains(cbac(:,2),'lower right peduncle'),1});

  Ycbt = single( 1.5 + .5*(Ycbc==2) - .5*(Ycbc==1)); 
  Ycbt( cat_vol_morph(Yml<1.5/4 | Yml>3.5/4 | (Ybscb>1.5 & Ybscb<4)) & Ycbc<1 ) = nan;
  Ycbt( cat_vol_morph( ~isnan(Ycbt) ,'l')==0 & ~isnan(Ycbt) ) = nan; 
  Ycbt = cat_vol_laplace3R(Ycbt,Ycbt==1.5,.00001);
  Ycbt = (Ycbt-1.5)*2;
  Ycbt = sign(Ycbt) .* (abs(Ycbt).^(1/4));
  [~,YI] = cat_vbdist(single(Ycbt<3 & Ycbt~=1.5)); 
  Ycoord = Ycbt(YI); clear YI;

  %figure, isosurface(Yml .* Ycb0,.6,Ycoord{1}); axis equal off; colormap jet; lighting none
end
% =========================================================================
function Ycoord = suit_amap_get_longitudinal(Yml,Ybscb,Ycba,cbac)
%suit_amap_get_longitudinal: y-axis between inferior (neg.) and superior (pos.)

  % seed definition for a upper and lower region of both peduncles
  Ycbd   = 2*(Ycba==cbac{cat_io_contains(cbac(:,2),'upper left peduncle'),1}) + ...
           1*(Ycba==cbac{cat_io_contains(cbac(:,2),'lower left peduncle'),1}) + ...
           2*(Ycba==cbac{cat_io_contains(cbac(:,2),'upper right peduncle'),1}) + ...
           1*(Ycba==cbac{cat_io_contains(cbac(:,2),'lower right peduncle'),1});

  % filtering 
  Ycbl = single( 1.5 + .5*(Ycbd==2) - .5*(Ycbd==1)); 
  Ycbl( cat_vol_morph(Yml<1.5/4 | Yml>3/4 | (Ybscb>2 & Ybscb<4)) & Ycbd<1  ) = nan; %clear Ycbd;
  Ycbl( (Ycba<=cbac{cat_io_contains(cbac(:,2),'brainstem'),1}) ) = nan;
  Ycbl = cat_vol_laplace3R(Ycbl*1000,Ycbl==1.5,.00001)/1000;
  Ycbl = (Ycbl-1.5)*2;
  Ycbl = sign(Ycbl) .* (abs(Ycbl).^(1/3));
  [~,YI] = cat_vbdist(single(Ycbl<3 & Ycbl~=1.5)); 
  Ycoord = Ycbl(YI); clear YI;    
 
  %figure, isosurface(Yml .* Ycb0,.6,Ycbtd ); axis equal off; colormap jet; lighting none
end 
% =========================================================================
function [Yth1,Ypp1,suit] = suit_amap_cb1thickness(Yml,Ycb,Yp0,vx_vol,suit,fi) 
  %% cerebellar measure for level 1+ ~cortex
  sx = -.1:.02:.1; % distance estimation offset
  if 0
    % direct measurement from central layer - minimum concept
    Ywd  = zeros(size(Yml)); for i=sx, Ywd = Ywd + 1/numel(sx) * cat_vbdist( single( Yml*4 > 2+i ) , Ycb>0 & Yml>1.5/4, vx_vol); end
    Ycd  = zeros(size(Yml)); for i=sx, Ycd = Ycd + 1/numel(sx) * cat_vbdist( single( Yml*4 < 2+i ) , Ycb>0 & Yml<3.5/4, vx_vol); end
    Yth1 = cat_vol_approx( (Ywd + Ycd) .* (Ycb>0 & Ywd<1e10 & Ycd<1e10) ); 
    Ypp1 = max(0,min(Ycb,((Yth1-Ywd).*(Ycb>0)*.5 + .5*Ycd) ./ max(eps,Yth1))); 
  else
    % pbt measurement form central layer
    Ywd  = zeros(size(Yml)); for i=sx, Ywd = Ywd + 1/numel(sx) * cat_vbdist( single(Yml*4 > 2 + i ) , Ycb>0 & Yml>1.5/4, vx_vol); end
    Ycd  = zeros(size(Yml)); for i=sx, Ycd = Ycd + 1/numel(sx) * cat_vbdist( single(Yml*4 < 2 + i ) , Ycb>0 & Yml<3.5/4, vx_vol); end
    Ythc = cat_vol_pbtp(min(3,max(1,Ycb.*Yml*4 - 1)),Ywd,Ycd); 
    Ythw = cat_vol_pbtp(min(3,max(1,4 - Ycb.*Yml*4)),Ycd,Ywd); 
    Yth1 = cat_vol_approx( min(cat_vol_approx(Ythc),cat_vol_approx(Ythw)) .* (Ycb>0 & Ywd<1e10 & Ycd<1e10) ); 
    Ypp1 = max(0,min(Ycb,((Yth1-Ywd).*(Ycb>0)*.5 + .5*Ycd) ./ max(eps,Yth1))); 
  end
  suit(fi).level1depth = mean(Yth1(Yp0(:)>1.5 & Yp0(:)<2.5));
end
% =========================================================================
function [Yp0a,Yprob] = suit_amap_AMAP(Ym,Yp0s,vx_vol,denoise,dobc)
%suit_amap_amap. Run AMAP segmentation.

  % denoise input map and prepare segmentmap
  Ymib   = Ym + 0;
  if denoise 
    cat_sanlm(Ymib,1,3); 
  end
  Ymib   = double(min(1.1,Ymib) .* (Yp0s>0.25)); 
  Yp0sb  = uint8(round(Yp0s)); % use Ym alternatively ?

  % run amap (warning wrong parameters will cause fatal errors that crash matlab)
  % [prob,amap_means,amap_stds] = cat_amap(Ymib, Yp0sb, n_classes, n_iters, sub, pve, init_kmeans,  ...
  %     job.extopts.mrf, vx_vol, iters_icm, bias_fwhm); 
  Yprob  = cat_amap(Ymib,Yp0sb, 3,10 + 90 * dobc ,64,  5,0,   0, vx_vol, 0, 100 * dobc); 
  if size(Yprob,4)>3, Yprob(:,:,:,4:end) = []; end 

  for ci = 1:3, Yprob(:,:,:,ci) = Yprob(:,:,:,ci) .* uint8(Yp0s>0.25); end
  Yp0a   = single(Yprob(:,:,:,1))/255 + single(Yprob(:,:,:,2))/255*2 + single(Yprob(:,:,:,3))/255*3;
end
% =========================================================================
function suit = suit_amap_write(Ymg,Yml,Yp0a,Ycb,Vmsk,Vgm,vx_vol,P,job,suit,fi)

% - as subfunction
% - replace saveCB/AMAP by simple defaults
    Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0 - c));
    sum4   = @(x,d) sum(sum(sum( single(x(:,:,:,d))/255 ))); 

    % save skull-stripped bias-corrected
    suit_amap_saveimg(Ycb .* Ymg, P.bc{fi}, Vmsk, Vgm, job, job.biascorr,...
      'SUIT noise & bias-corrected cerebellum (0-BG, 1-WM)','uint8');

    %% 4-class setup 
    Yprobt  = zeros([size(Yml),4],'uint8'); 
    for ci=1:4, Yprobt(:,:,:,ci) = cat_vol_ctype( Yp0toC(Yml .* 4 .* Ycb,ci) * 255 ); end
    Yprobt(:,:,:,4) = cat_vol_ctype( single(Yprobt(:,:,:,4)) + Yp0toC(Yml .* 4 .* Ycb,5) * 255); % add too high values
    suit(fi).vol_cggw_abs_CGGW = [sum4(Yprobt,1) sum4(Yprobt,2) sum4(Yprobt,3) sum4(Yprobt,4)] * prod(vx_vol) / 1000;
    suit(fi).vol_cggw_rel_CGGW = suit(fi).vol_cggw_abs_CGGW ./ suit(fi).vol_cbmsk1;
    if 1
      suit_amap_saveimg(Yprobt(:,:,:,1), P.cb1{fi}, Vmsk, Vgm, job, job.cbseg, ...
        'SUIT segment 1/4 - CSF', 'uint8');
      suit_amap_saveimg(Yprobt(:,:,:,2), P.cb2{fi}, Vmsk, Vgm, job, job.cbseg, ...
        'SUIT segment 2/4 - GM - molecular layer', 'uint8');
      suit_amap_saveimg(Yprobt(:,:,:,3), P.cb3{fi}, Vmsk, Vgm, job, job.cbseg, ...
        'SUIT segment 3/4 - GM - granular layer', 'uint8');
      suit_amap_saveimg(Yprobt(:,:,:,4), P.cb4{fi}, Vmsk, Vgm, job, job.cbseg, ...
        'SUIT segment 4/4 - WM', 'uint8');
    end
    % label map
    suit_amap_saveimg(Ycb .* max(0,Yml*4-1)/3, P.cb0{fi}, Vmsk, Vgm, job, job.cbseg, ...
      'SUIT cerbellar layers with partial volumes (0-CSF,0.33-GMl,66-GMh,1-WM)', 'uint8');
    clear Yprobt;
  

    %% 3-class setup 
    % - although 2-GM classes are used in our segmentation concept
    %   a simpler 2 output for analysis migh be more better
    Yprobc  = zeros([size(Ymg),3],'uint8'); 
    for ci=1:3, Yprobc(:,:,:,ci) = cat_vol_ctype( Yp0toC(Ymg .* 3 .* Ycb,ci) * 255 ); end
    Yprobc(:,:,:,3) = cat_vol_ctype( single(Yprobc(:,:,:,3)) + Yp0toC(Ymg .* 3 .* Ycb,5) * 255); % add too high values
    suit(fi).vol_cggw_abs_CGW  = [sum4(Yprobc,1) sum4(Yprobc,2) sum4(Yprobc,3)] * prod(vx_vol) / 1000;
    suit(fi).vol_cggw_rel_CGW  = suit(fi).vol_cggw_abs_CGW ./ suit(fi).vol_cbmsk1;
    if 1
      % new general function
      suit_amap_saveimg(Yprobc(:,:,:,1), P.cc1{fi}, Vmsk, Vgm, job, job.cbseg, ...
        'SUIT segment 1/3 - CSF', 'uint8');
      suit_amap_saveimg(Yprobc(:,:,:,2), P.cc2{fi}, Vmsk, Vgm, job, job.cbseg, ...
        'SUIT segment 2/3 - GM', 'uint8');
      suit_amap_saveimg(Yprobc(:,:,:,3), P.cc3{fi}, Vmsk, Vgm, job, job.cbseg, ...
        'SUIT segment 3/3 - WM', 'uint8');
      % label map
      suit_amap_saveimg(Ycb .* Ymg, P.cc0{fi}, Vmsk, Vgm, job, job.cbseg, ...
        'SUIT cerbellar layers with partial volumes (0-CSF,0.5-GM,1-WM)', 'uint8');
      clear Yprobc;
    end

    %% final AMAP 3-class segmentation
    if job.amap 
      [Yp0,Yprob] = suit_amap_AMAP(Ymg,Yp0a .* Ycb,vx_vol,0,0);
      suit(fi).vol_amap_abs_CGW  = [sum4(Yprob,1) sum4(Yprob,2) sum4(Yprob,3)] * prod(vx_vol) / 1000;
      suit(fi).vol_amap_rel_CGW  = suit(fi).vol_amap_abs_CGW ./ suit(fi).vol_cbmsk1;
      % new general function
      suit_amap_saveimg(Yprob(:,:,:,2), P.amap1{fi}, Vmsk, Vgm, job, job.amap, ...
        'SUIT AMAP segment 1/3 - GM', 'uint8');
      suit_amap_saveimg(Yprob(:,:,:,3), P.amap2{fi}, Vmsk, Vgm, job, job.amap, ...
        'SUIT AMAP segment 2/3 - WM', 'uint8');
      suit_amap_saveimg(Yprob(:,:,:,1), P.amap3{fi}, Vmsk, Vgm, job, job.amap > 1, ...
        'SUIT AMAP segment 3/3 - CSF', 'uint8');
      % label map
      suit_amap_saveimg(Yp0, P.amap0{fi}, Vmsk, Vgm, job, job.amap, ...
        'SUIT cerbellar layers with partial volumes (1-CSF,2-GM,3-WM)', 'uint8x3');
      clear Yp0 Yprob; 
    end
end
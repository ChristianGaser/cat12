function cat_tst_qa_simerrBWP( datadir, qaversions, fasttest, rerun ) 
%% BWP resolutions
%  ------------------------------------------------------------------------
%  This function modifies an existing segmentation to test how segmentation
%  errors effects the quality measuresments by modifying the Yp0 label map.
%  (goal is the evaluation of the general mean effect)  
%
%   (1) Simulation of skull-stripping errors
%        - missing parts
%        - additional parts
%        - interpreation of intensity scaled Ym as Yp0 segment
%       >> this works quite well and the QC is stable even for low Kappas 
%
%   (2) Simulation of segmentation bias 
%        - WM overestimation by WM dilatation (> GM erosion) 
%        - WM underestimation by WM erosion (GM dilation)
%        - CSF overestimation by GM/WM erosion (CSF dilation) 
%        - CSF underestimation by GM/WM dilation (CSF erosion)
%        - use of min/max-filters rather than binar morphometric
%          operations 
%       >> this was a bit challenging but it is also ok
%
%   (3) Distortion of p0 by noise and smoothing or random deformations
%       >> noise worked well 
%
%   (4) Role of partial volume effects (PVE) - yes/no
%       >> not important 
%
%   Allows quasi-random cases and levels:
%    - size in percent of the TIV, e.g., 0:5:20 or 0,5,10,20,40
%    - dist in mm, e.g., 0:5:20 or 0,5,10,20,40;
%    - Kappa value as final outcome, e.g., range 0.5 - 1.0 
%
%  ------------------------------------------------------------------------

%#ok<*UNRCH>
%#ok<*SAGROW>  

cat_io_cprintf([0 0.5 0],'\n\n== Run cat_tst_qa_simerrBWP ==\n') 

% ### datadir ### 
if ~exist( 'datadir' , 'var' )
  Pddir  =  '/Volumes/SG5TB/MRData/202503_QA/BWP';
else
  Pddir  = fullfile(datadir,'BWP'); 
end
% ### QC version ### 
if ~exist( 'qaversions' , 'var')
  qaversions = {
    ...'cat_vol_qa201901';  % classic version (quite stable since 2016)
    'cat_vol_qa201901x'; % refined, debugged version of 201901 
    ...'cat_vol_qa202110';  % second classic version (successor of 201901)
    ...'cat_vol_qa202110x'; % refined, debugged version of 202110   
    ...'cat_vol_qa202205';  % last regular version before update (successor of 202110, stopped)
    ...'cat_vol_qa202310';  % redesigned version based on 201901 and 202110  * default *
    ...'cat_vol_qa202412';  % experimental version with internal segmentation >> qcseg
     };
end
if ~exist( 'fasttest', 'var'), fasttest = 0; end
if ~exist( 'rerun', 'var'), rerun = 0; end

outdir = {fullfile(fileparts(Pddir),'BWPE')};

% Brain Web Phantom (BWP): 
% * Using also cases with strong inhomogeneity causees severe changes of 
%   NCR and IQR (as expected), whereas contrast and ECR are quite stable.
P   = {fullfile(Pddir,'BWPC_HC_T1_pn1_rf020pA_vx100x100x100.nii')};
Pt  = {
     {
       fullfile(Pddir,'BWPC_HC_T1_pn1_rf020pA_vx100x100x100.nii');
       fullfile(Pddir,'BWPC_HC_T1_pn3_rf020pA_vx100x100x100.nii');
       fullfile(Pddir,'BWPC_HC_T1_pn5_rf020pA_vx100x100x100.nii');
       fullfile(Pddir,'BWPC_HC_T1_pn7_rf020pA_vx100x100x100.nii');
       fullfile(Pddir,'BWPC_HC_T1_pn9_rf020pA_vx100x100x100.nii');
       fullfile(Pddir,'BWPC_HC_T1_pn1_rf020pB_vx100x100x100.nii');
       fullfile(Pddir,'BWPC_HC_T1_pn3_rf020pB_vx100x100x100.nii');
       fullfile(Pddir,'BWPC_HC_T1_pn5_rf020pB_vx100x100x100.nii');
       fullfile(Pddir,'BWPC_HC_T1_pn7_rf020pB_vx100x100x100.nii');
       fullfile(Pddir,'BWPC_HC_T1_pn9_rf020pB_vx100x100x100.nii');
       fullfile(Pddir,'BWPC_HC_T1_pn1_rf020pC_vx100x100x100.nii');
       fullfile(Pddir,'BWPC_HC_T1_pn3_rf020pC_vx100x100x100.nii');
       fullfile(Pddir,'BWPC_HC_T1_pn5_rf020pC_vx100x100x100.nii');
       fullfile(Pddir,'BWPC_HC_T1_pn7_rf020pC_vx100x100x100.nii');
       fullfile(Pddir,'BWPC_HC_T1_pn9_rf020pC_vx100x100x100.nii');
       ... other bias level introduce a lot of variapCe in pCR as they localy change the noise pattern
       ... same for negative fields (probably even worse)
       ... but maybe your measure is now biased a bit > ########### may run an additional test later ##############
       ... maybe as some add data loop with other BWP settings ...
     }
  };
Pt{1}( cellfun(@(x) exist(x,'file'),Pt{1})==0 ) = [];

% 60% bias as average bias case (or 40% as old worst case) but with stronger differences for noise
replaceRF = 20; 
if replaceRF
  for pti = 1:numel(Pt)
    Pt{pti} = strrep(Pt{pti},'_rf020',sprintf('_rf%03d',replaceRF)); 
  end
end

% some other normalized Yp0 maps
%Px  = {
%  { fullfile( spm('dir'),'toolbox','cat12','templates_MNI152NLin2009cAsym','trimmed_Template_T1_masked.nii') }; 
%  }; 
Pp0 = spm_file(P,'prefix',['mri' filesep 'p0']); 
%res = repmat( ( 1:0.25:2 )', 1,3);
%ss  = repmat( ( 0.25:0.25:2 )', 1,3);

Pcolor1   = cat_io_colormaps('cold',4);
Pcolor2   = flip( cat_io_colormaps('hot',4) ,1);
Pcolor    = [Pcolor1(2:end-1,:); mean([Pcolor1(2:end-1,:);Pcolor2(2:end-1,:)]); Pcolor2(2:end-1,:)];
Pcolor    = repmat( Pcolor , max(1,numel(Pt{1}) / 5) , 1);
QR        = {'NCR','ICR','res_ECR','contrastr','IQR','SIQR'};
QRname    = {'NCR','ICR','ECR','CON','IQR','SIQR'};
markrange = [-.5 .5] * 2;
fastname  = {'full','fast'};
rerunqc    = rerun; 
rerunkappa = rerun; 

if fasttest && numel(Pt{1})>=15
  % 1%, 3%, and 9% noise of field A + 9% of field B and C 
  % (e.g. all worst cases that drives the outliers)
  Pt{1} = Pt{1}([1:2:5,10,15]); 
else
  Pt{1} = Pt{1}(1:end);
end
printdir  = fullfile(fileparts(Pddir),'+results',['BWPE_' fastname{fasttest+1} '_' datestr(clock,'YYYYmm')]);

f29 = figure(3235);  
f29.Visible = 'off';
f29.Interruptible = 'off';
set(f29,'Position',[0 0 1440 200],'Name','Skull-Stripping','color',[1 1 1]); 
clear valr2;

qais   = 1:numel(qaversions);
for qai = qais
  %%
  qafile = qaversions{qai}; 
  
  % (1) Simulation of skull-stripping errors
  for pi = 1:numel(P) 
    cat_io_cprintf('blue',sprintf(' Prepare subject %0.0f/%0.0f with %d cases\n',pi,numel(P),numel(Pt{pi}))); 
  
    [~,printname] = spm_fileparts(P{pi}); 
    
    V   = spm_vol(P{pi}); 
    Vp0 = spm_vol(Pp0{pi}); 
    
    Y   = spm_read_vols(V); 
    Yp0 = spm_read_vols(Vp0); 
    
    % simple global intensity normalization 
    bth = cat_stat_kmeans(Y(Yp0(:)==0)); 
    wth = cat_stat_kmeans(Y(round(Yp0(:))==3));
    Ym  = (Y - bth) ./ abs(diff([bth wth])) ; 
    Yp0e = max(Yp0,(Yp0==0) .* min(3,Ym*3));  
    
    % estimate brainmask distance as weighting 
    Yp0d  = cat_vbdist(single(Yp0>0.5)); 
    Yp0di = cat_vbdist(single(Yp0<0.5)); 
    
    fx = 1; dx = 40; % some paraemters to control the size of the modified area that I am not furhter using
    Yp0dw = max( 0 , (Yp0d>0) - Yp0d/dx) + max( 0 , (Yp0di>0) - Yp0di/dx);
    
    
    
  
    %% (1) Brain masking test
    %  ----------------------------------------------------------------------
    %  There are obvious differences for strong bias in the NCR (and ICR) 
    %  rating that are related to the BWP, i.e. real changes. 
    % 
    
    % create some quasi-random changes
    rng('default') 
    rng(33*pi); Ymsk = cat_vol_smooth3X( randn(size(Y)) , 8 * fx) * 16; 
    rng(48*pi); Ymsk = Ymsk + cat_vol_smooth3X( randn(size(Y)) , 4 * fx) * 8;
    rng(99*pi); Ymsk = Ymsk + cat_vol_smooth3X( randn(size(Y)) , 2 * fx) * .5;
    Ymsk = (Ymsk - mean(Ymsk(:))) ./ std(Ymsk(:));
    
    % bstr 16 is kappa ~0.12 and maybe to bad enough to be ingnored
    if fasttest
      bstr = [0 4 8 12 14 ]; 
    else
      bstr = [0 1 2 4 6 8 10 12 14 ]; 
    end
    ptix  = 0; 
    clear ktab Vmskb; 
    cat_io_cprintf('blue',sprintf(' Brain masking simulation of %d subcase of %d testcases\n',numel(bstr),numel(Pt{pi}))); 
    for bstri = 1:numel(bstr) 
      % create brainmask 
      %  I tried first just to use some distance-weighted cloudy field to create 
      %  some distortions but this was generally not enough. 
      %   bstr = 0.5; % strenght of distortion from 0 (none) to 8 (severe)
    
      Pp0BSE = spm_file(Pt{pi},'path',outdir,'prefix',sprintf('p0BSE%d',bstr(bstri))); 
      if ~all(cellfun(@exist,Pp0BSE)==2)
        Ymskb = cat_vol_smooth3X( (min(1,Yp0) + Ymsk * bstr(bstri).^2 .* Yp0dw.^4 )>=0.5 , ...
          min(8,4*bstr(bstri))) >= 0.5 - 0.01*min(2,bstr(bstri)); 
  
        %  Hence, I added some simplyfied Skull-Stripping like model of opening
        %  and closing to avoid strange holes or large outstanding parts.
        [Ymskbr,resr] = cat_vol_resize(Ymskb,'reduceV',1,4); 
        Ymskbr = cat_vol_morph(Ymskbr,'l');
        Ymskbr = cat_vol_morph(Ymskbr,'do',bstr(bstri)/2);
        Ymskbr = cat_vol_morph(Ymskbr,'ldc',bstr(bstri)/2);
        Ymskbr = cat_vol_morph(Ymskbr,'ldo',bstr(bstri)/2);
        Ymskb  = cat_vol_resize(Ymskbr,'dereduceV',resr); 
      end
      
      % write data
      for pti = 1:numel(Pt{pi})
        ptix = ptix + 1; 
        Pp0BSEpti = spm_file(Pt{pi}{pti},'path',outdir,'prefix',sprintf('p0BSE%d',bstr(bstri))); 
      
        if ~exist(outdir{1},'dir'), mkdir(outdir{1}); end
        if ~exist(Pp0BSEpti,'file')
          Vmskb(bstri,pti) = Vp0; Vmskb(bstri,pti).fname = Pp0BSEpti; 
          spm_write_vol(Vmskb(bstri,pti),Yp0e .* Ymskb);
        else
          Vmskb(bstri,pti) = spm_vol(Pp0BSEpti);
        end
        
        % T1 copies
        if ~exist(spm_file(Pt{pi}{pti},'path',outdir,'prefix',sprintf('BSE%d',bstr(bstri))),'file')
          copyfile( Pt{pi}{pti} , spm_file(Pt{pi}{pti},'path',outdir,'prefix',sprintf('BSE%d',bstr(bstri))) );
        end
        if ~exist(spm_file(Pt{pi}{pti},'path',outdir,'prefix',sprintf('mBSE%d',bstr(bstri))),'file')
          copyfile( Pt{pi}{pti} , spm_file(Pt{pi}{pti},'path',outdir,'prefix',sprintf('mBSE%d',bstr(bstri))) ); 
        end
      end
    end 
    
    
    %% estimate QC and Kappa
    qav = cat_vol_qa('p0',{Vmskb(:).fname},struct('prefix',[qafile '_'],'version',qafile,'rerun',rerunqc));
    if ~exist('valr','var') 
      [~,valr] = eva_vol_calcKappa({Vmskb(:).fname},Pp0,struct('recalc',rerunkappa)); 
    end


    %% Evaluation 
    set(f29,'Position',[0 0 1440 200],'Name','Skull-Stripping','color',[1 1 1]); 
    clf(f29); hold on; clear tab; 
    qav2r = reshape(qav,numel(bstr),numel(Pt{pi})); 
    RMSE  = @(x) mean(x.^2).^0.5; 
    for qri = 1:numel(QR) 
      sp29 = subplot('Position',[ (qri-1) * 1 / (numel(QR) + 2) + 0.035 , 0.17 , 0.73 / (numel(QR) + 2) , 0.73],'replace');  
      for pti = 1:numel(Pt{pi}) % BWP test case
        for bstri = 1:numel(bstr) % skull-stripping
          try
            tab.(QR{qri})(bstri,pti) = qav2r(bstri,pti).qualityratings.(QR{qri});
          catch
            tab.(QR{qri})(bstri,pti) = nan; 
          end
          ktab(bstri,pti) = valr{bstri+1,5};  
        end
      end
      
      tab.(QR{qri}) = tab.(QR{qri}) - repmat(tab.(QR{qri})(1,:),numel(bstr),1);
      hold on; pl=plot([0 10],[0 0]); pl.Color = [0.3 0.3 0.3]; 
      for pti = 1:numel(Pt{pi}), px = plot( ktab(:,1), tab.(QR{qri})(:,pti) ); set(px,'Color',Pcolor(pti,:) ); end
      pt = plot( ktab(:,1), mean(tab.(QR{qri}),2)); set(pt,'LineWidth',2,'Color',[0 0 0],'Marker','x');
      xlim([0.25,1]), ylim(markrange); grid on; box on; 
      title(strrep(QRname{qri},'_','\_'));
      xlabel('Kappa'); ylabel('mark error');
      set(sp29,'ytick',markrange(1):0.5:markrange(2),'yticklabel',num2str((markrange(1):0.5:markrange(2))','%0.2f '), ...
               'xtick',0:0.25:1,'xticklabel',num2str((0:0.25:1)','%0.2f ') );
    end

    % absolute error
    subplot('Position',[ (numel(QR)) / (numel(QR) + 2) + 0.035 ,  0.17 , 0.72 / (numel(QR) + 2) , 0.73]);  
    cat_plot_boxplot([ tab.NCR(:) , tab.ICR(:) , tab.res_ECR(:) , tab.contrastr(:) , ...
      tab.IQR(:) , tab.SIQR(:) ],struct('names',{QRname},'ygrid',0,'ylim',markrange,'style',4,'datasymbol','o','usescatter',1))
    ax=gca; ax.XTickLabelRotation = 0; ax.YGrid = 'on'; ax.FontSize = 9.5;
    set(ax,'ytick',markrange(1):0.5:markrange(2),'yticklabel',num2str((markrange(1):0.5:markrange(2))','%0.2f '));
    title('Absolute error'); xlabel('measure'); ylabel('mark error');

    % RMSE
    subplot('Position',[ (numel(QR) + 1) / (numel(QR) + 2) + 0.032 ,  0.17 , 0.72 / (numel(QR) + 2) , 0.73]);  
    rms     = @(a)   max(0,cat_stat_nanmean(a.^2).^(1/2));
    rmseval = [ rms(tab.NCR(:)) , rms(tab.ICR(:)) , rms(tab.res_ECR(:)) , rms(tab.contrastr(:)) , ...
      rms(tab.IQR(:)) , rms(tab.SIQR(:)) ]; 
    bar(rmseval(1:end)); 
    ylim([0,1]); xlim([.4 6.6]); xticklabels(QRname); 
    h = gca; h.YGrid = 'on'; h.XTickLabelRotation = 0; 
    for fi = 1:numel(rmseval)
      text(fi-.38, rmseval(fi) + .03, sprintf('%0.2f',rmseval(fi)),'FontSize',8,'Color',[.0 .2 .4]);
    end  
    ax=gca; ax.XTickLabelRotation = 0; ax.YGrid = 'on'; ax.FontSize = 9.5;
    set(ax,'ytick',markrange(1):0.25:markrange(2),'yticklabel',num2str((markrange(1):0.25:markrange(2))','%0.2f '));
    title(sprintf('RMSE (%0.3f)',mean(rmseval))); xlabel('measure'); ylabel('mark error');
  
    % caption
    capf = 1.3; 
    f29.Position(4) = f29.Position(4)*capf; ax = findobj(f29.Children,'type','Axes');
    for axi=1:numel(ax), ax(axi).Position(4) = ax(axi).Position(4)/1.3; ax(axi).Position(2) = (capf-1) * capf; end
    an = annotation('textbox',[0 0 1 ((capf-1)/capf + 0.02) ],'FitBoxToText','off','LineStyle','none','String',{['\bfFig. 2 (' ...
        strrep(qafile,'_','\_') '):\rm '...
        'Shown are the results of the quality ratings (as school marks) for simuated brain extraction problems ' ... 
        'quantified by the Kappa rating. Even in case of strong errors (Kappa<0.5) the measures stay relative stable ' ...
        'and it has to be considered that the worst Kappa rating of CAT12 was about 0.8. ' ...
        'Blue lines indicate low noise data (BWP noise 1 and 3%), whereas red lines represent noisy data (BWP noise 7 and 9%). ' ... 
        'For the\it noise contrast rating\rm (\bfNCR\rm) a half mark is equal to one BWP noise levels and indicate light (motion) artifacts. ' ...
        'For the\it inhomogeneity contrast rating\rm (\bfICR\rm) a 1/4 grade presents a change of 20% of the BWP inhomogeneity. ' ...
        'For\it edge contrast rating\rm (\bfECR\rm) a change of 1 grade correspond to a change of 0.5 mm that is somewhat similar like ' ...
        'Gaussian smoothing of 1.4 mm). ' ...
        'The\it structural image quality rating\rm (\bfSIQR\rm) is the weighted combination of the single measures. ' ...
        '']},'Fontsize',12);
  
    subdir = fullfile(printdir,'Fig2-skullstripping');
    if ~exist(subdir,'dir'), mkdir(subdir); end
    print(f29, '-djpeg', '-r300', fullfile(subdir,sprintf('%s_skullstrippingphantom_%s',printname,qafile))); 
    %delete(an);
    
    
    
     
    %% (2) Tissue Segment Modification  
    %  ----------------------------------------------------------------------
    %  Tissue changes affect the NCR and ECR but in oposit direction what  
    %  partially depends on the contrast rating. ICR is quite stable.
    %  Especially the dilatation of the WM was critical - probably because it
    %  strongly reduce the GM ribbon. 
    %  
    
    if fasttest
      levels = 2; 
    else
      levels = 1:2; %3:4; %1:2;
    end
    segmentation = cell(1,numel(levels)); seg = segmentation;
    for li = 1:numel(levels)
      seg{li} = sprintf('L%d',levels(li));
      segmentation{li} = sprintf('segmentation L%d',levels(li)); 
    end 
    clear Vmsks; 
    qav2      = cell(1,numel(levels)); 
    tabsegt   = cell(1,numel(levels));
    ktabsegt  = cell(1,numel(levels));
    tcaselong = {'original','dilated WM','eroded WM','dilated CSF','eroded CSF', ...
                 'dilated WM and CSF (~eroded GM)','eroded WM and CSF (dilated GM'};
    %%
    for li = 1:numel(levels)
      tcase = {'org','dw','ew','dc','ec','dd','ee'};
      for ti = 2:numel(tcase), tcase{ti} = sprintf('%s%d',tcase{ti},levels(li)); end 
  
      for ti = 1:numel(tcase)
        fprintf(' %s\n', tcase{ti}(1:3))
        deval  = str2double(tcase{ti}(3)) / 2;
        deval2 = mod(str2double(tcase{ti}(3)) + 1,2) / 2 + 0.5;
        Yp0d   = Yp0; 
        Pp0SEG = spm_file(Pt{pi},'path',outdir,'prefix',sprintf('p0SEG%s',tcase{ti}));
        if ~all(cellfun(@exist,Pp0SEG)==2)
          for it = 1:ceil( deval )
            Yp0o = Yp0d;
            switch tcase{ti}(1:2)
              case 'or'
                Yp0d = Yp0;
              case 'dw' % increased NCR, reduced ECR, but ok
                Yp0d = max( Yp0d  .* (Yp0d<=2)     , cat_vol_localstat(single(Yp0d) , Yp0d>=2,1,3));
              case 'ew' % strongly increased NCR, reduced ECR, NOT OK
                Yp0d = max( Yp0d  .* (Yp0d<=2)     , cat_vol_localstat(single(Yp0d) , Yp0d>=2,1,2));
              case 'dc' % increased NCR, reduced ECR, but ok (similar to dw)
                Yp0d = max( Yp0d  .* (Yp0d>=2.01)  , cat_vol_localstat(single(Yp0d) , Yp0d<=2.01,1,2));
              case 'ec' % strongly increased NCR, reduced ECR, NOT OK (better than ew) 
                Yp0d = max( Yp0d  .* (Yp0d>=2.01)  , cat_vol_localstat(single(Yp0d) , Yp0d<=2.01,1,3));
              case 'ee' % this is ok
                Yp0d = max( Yp0d  .* (Yp0d<=2)     , cat_vol_localstat(single(Yp0d) , Yp0d>=2,1,2));
                Yp0d = max( Yp0d .* (Yp0d>=2.01)   , cat_vol_localstat(single(Yp0d) , Yp0d<=2.01,1,3));
              case 'dd' % dilate WM and CSF result in very small GM ribbon and strong overestimation of NCR and unterstimtion of ECR
                Yp0d = max( Yp0d  .* (Yp0d<=2)     , cat_vol_localstat(single(Yp0d) , Yp0d>=2,1,3));
                Yp0d = max( Yp0d .* (Yp0d>=2.01)   , cat_vol_localstat(single(Yp0d) , Yp0d<=2.01,1,2));
              otherwise
                error('unknown case')
            end
            Yp0dc = Yp0d - Yp0o; 
            Yp0d  = Yp0o + min(deval,0.75) * min(min(deval2,1),max(-min(deval2,1),Yp0dc)); %min(min(deval2,1)/2,max(-min(deval2,1)/2,Yp0dc));
          end
          
          
          % write data
          for pti = 1:numel(Pt{pi})
            if ~exist(outdir{1},'dir'), mkdir(outdir{1}); end
            Pp0SEGpti = spm_file(Pt{pi}{pti},'path',outdir,'prefix',sprintf('p0SEG%s',tcase{ti}));
            if ~exist(Pp0SEGpti,'file')
              Vmsks(ti,pti) = Vp0; Vmsks(ti,pti).fname = Pp0SEGpti; 
              spm_write_vol(Vmsks(ti,pti),Yp0d);
            else
              Vmsks(ti,pti) = spm_vol(Pp0SEGpti);
            end
            
            % T1 copies
            if ~exist(spm_file(Pt{pi}{pti},'path',outdir,'prefix',sprintf('SEG%s',tcase{ti})),'file')
              copyfile( Pt{pi}{pti} , spm_file(Pt{pi}{pti},'path',outdir,'prefix',sprintf('SEG%s',tcase{ti})) ); 
            end
            if ~exist(spm_file(Pt{pi}{pti},'path',outdir,'prefix',sprintf('mSEG%s',tcase{ti})),'file')
              copyfile( Pt{pi}{pti} , spm_file(Pt{pi}{pti},'path',outdir,'prefix',sprintf('mSEG%s',tcase{ti})) ); 
            end
          end
        else
          for pti = 1:numel(Pt{pi})
            Pp0SEGpti = spm_file(Pt{pi}{pti},'path',outdir,'prefix',sprintf('p0SEG%s',tcase{ti}));
            Vmsks(ti,pti) = Vp0; Vmsks(ti,pti).fname = Pp0SEGpti;
          end
        end
      end
    
  
      %% estimate QC
  %    eval(sprintf('qav2{li} = %s(''p0'',{Vmsks(:).fname},qaopt);',qafile)); 
      qav2{li} = cat_vol_qa('p0',{Vmsks(:).fname},struct('prefix',[qafile '_'],'version',qafile,'rerun',rerunqc));
      if ~exist('valr2','var') || numel(valr2)<li || (size(valr2{li},1)-3) ~= numel({Vmsks(:).fname})
        [~,valr2{li}] = eva_vol_calcKappa({Vmsks(:).fname},Pp0,struct('recalc',rerunkappa)); 
      end    
  
      %% Evaluation  
      %  NCR and ECR are less stable, what partially depend non the contrast 
      %  but also other factors. ICR is quite stable (does this support some 
      %  further use?)
      clf(f29); hold on
      set(f29,'Position',[0 0 1440 200],'color',[1 1 1],'Name',sprintf('segmentation L%d',levels(li)));
      qav2r = reshape(qav2{li},numel(tcase),numel(Pt{pi})); 
      RMSE = @(x) mean(x.^2).^0.5; 
      for qri = 1:numel(QR) 
        subplot('Position',[ (qri-1) * 1 / (numel(QR) + 2) + 0.035 , 0.17 , 0.73 / (numel(QR) + 2) , 0.73],'replace');  
        for pti = 1:numel(Pt{pi}) % BWP test case
          for bstri = 1:numel(tcase) % segmentation case
            try 
              tabsegt{li}.(QR{qri})(bstri,pti) = qav2r(bstri,pti).qualityratings.(QR{qri});
            catch
              tabsegt{li}.(QR{qri})(bstri,pti) = nan;
            end
            ktabsegt{li}(bstri,pti)            = valr2{li}{bstri+1,5};  
          end
        end
        tabsegt{li}.(QR{qri}) = (tabsegt{li}.(QR{qri}) - repmat(tabsegt{li}.(QR{qri})(1,:),numel(tcase),1) );
        tcase2 = {'org','dw','ew','dc','ec','dd','ee'};
        tabtmp = (tabsegt{li}.(QR{qri})(2:size(tabsegt{li}.(QR{qri}),1),:))'; tabtmp = reshape(tabtmp,numel(tabtmp) / 6,6); 
        cat_plot_boxplot( tabtmp ,struct('ygrid',1,'names',{tcase2(2:end)},...
          'ylim',markrange,'style',4,'datasymbol','o','usescatter',1)); %num2str(ktabsegt(:,1),'%0.2f')
        hold on, plot([0 10],[0 0],'color',[0.3 0.3 0.3]);
  
        title(strrep(QRname{qri},'_','\_'));
        xlabel('Phantom');  ylabel('mark error')
      end
    
      % absolute error
      subplot('Position',[ numel(QR) / (numel(QR) + 2) + 0.035 ,  0.17 , 0.73 / (numel(QR) + 2) , 0.73],'replace'); 
      cat_plot_boxplot([ tabsegt{li}.NCR(:) , tabsegt{li}.ICR(:) , tabsegt{li}.res_ECR(:) , ...
        tabsegt{li}.contrastr(:) , tabsegt{li}.IQR(:) , tabsegt{li}.SIQR(:) ], ...
        struct('names',{QRname},'ygrid',1,'ylim',markrange,'style',4,'datasymbol','o','usescatter',1))
      ax=gca; ax.XTickLabelRotation = 0; ax.YGrid = 'on'; ax.FontSize = 9.5;
      title('Absolute error'); xlabel('measure'); ylabel('mark error')
  
      % RMSE
      subplot('Position',[ (numel(QR) + 1) / (numel(QR) + 2) + 0.032 ,  0.17 , 0.72 / (numel(QR) + 2) , 0.73]);  
      rmseval = [ rms(tabsegt{li}.NCR(:)) , rms(tabsegt{li}.ICR(:)) , rms(tabsegt{li}.res_ECR(:)) , ...
        rms(tabsegt{li}.contrastr(:)) , rms(tabsegt{li}.IQR(:)) , rms(tabsegt{li}.SIQR(:)) ]; 
      bar(rmseval(1:end)); 
      ylim([0,1]); xlim([.4 6.6]); xticklabels(QRname); 
      h = gca; h.YGrid = 'on'; h.XTickLabelRotation = 0; 
      for fi = 1:numel(rmseval)
        text(fi-.38, rmseval(fi) + .03, sprintf('%0.2f',rmseval(fi)),'FontSize',8,'Color',[.0 .2 .4]);
      end  
      ax=gca; ax.XTickLabelRotation = 0; ax.YGrid = 'on'; ax.FontSize = 9.5;
      set(ax,'ytick',markrange(1):0.25:markrange(2),'yticklabel',num2str((markrange(1):0.25:markrange(2))','%0.2f '));
      title(sprintf('RMSE (%0.3f)',mean(rmseval))); xlabel('measure'); ylabel('mark error');


      % caption
      capf = 1.5; 
      f29.Position(4) = f29.Position(4)*capf; ax = findobj(f29.Children,'type','Axes');
      for axi=1:numel(ax), ax(axi).Position(4) = ax(axi).Position(4)/capf; ax(axi).Position(2) = ax(axi).Position(2) + 0.8*(capf-1)/capf; end
      annotation('textbox',[0 0 1 (capf-1)/capf - 0.02 ],'FitBoxToText','off','LineStyle','none','String',{['\bfFig. 3 (' ...
        strrep(qafile,'_','\_') '):\rm '...
          'Shown are the deviations of different quality ratings in case of 6 different simulated segmentation problems: ' ...
          '\bf (dw)\rm dilatated WM (WM overestimation (GM unerstimation), with a too low T1 WM peak that include also the GM-WM partial volume ' ...
          'caused e.g. by inconsidering of WM hyperintesities); ' ...
          '\bf (ew)\rm erorded WM (WM underestimation that is less typical than the WM overestimation); ' ...
          '\bf (dc)\rm dilatated CSF (CSF overestimation/GM underestimation, e.g., ' ...
          'in image with very strong GM-WM constrast and therefore low BG/CSF-GM-contrast ' ...
          'where the CSF peak also include the CSF-GM partial volumes); ' ...
          '\bf (ec)\rm eroded CSF (CSF understimation/GM overestimation, inverse of the dc-case but for similar data); ' ...
          '\bf (dd)\rm dilated CSF and WM (overstimation of CSF and WM and underestimation of GM); and ' ...
          '\bf (ee)\rm eroded CSF and WM (understimation of CSF and WM and overestimation of GM). ' ...
          'For the\it noise contrast rating\rm (\bfNCR\rm) a half mark is equal to one BWP noise levels and indicate light (motion) artifacts, ' ...
          'where for the\it inhomogeneity contrast rating\rm (\bfICR\rm) a 1/4 grade presents a change of 20% of the BWP inhomogeneity. ' ...
          'For\it edge contrast rating\rm (\bfECR\rm) a change of 1 grade correspond to a change of 0.5 mm that is somewhat similar like ' ...
          'Gaussian smoothing of 1.4 mm). ' ...
          'The\it structural image quality rating\rm (\bfSIQR\rm) is the weighted combination of the single measures. ' ... 
          '']},'Fontsize',12);
  
      subdir = fullfile(printdir,'Fig3-classissues');
      if ~exist(subdir,'dir'), mkdir(subdir); end
      print(f29, '-djpeg', '-r300', fullfile(subdir,sprintf('%s_segmentationphantomL%d_%s',printname,li,qafile))); 
    end
     
   
    
    
    
    %% Final plot of all phantoms
    %  ----------------------------------------------------------------------
    %  It shows that the measures are quite stable for skull-stripping errors 
    %  (blue) but that noisy segmentations (yellow) but that tissue segmentation
    %  errors (red) are more challenging even if Kappa is quite good. 
    %  It also shows how strong NCR and ECR depend on the contrast estimation. 
    %  The scatter and the final boxplot also show that the measurement
    %  errors are typically within the range of one mark.
    %
    %  However, I need to show that it is below 0.5 for normal image quality
    %  as far this presents our outlier criteria (=5 rps) 
    %  ----------------------------------------------------------------------
    marker = 'o^v><sp+x'; clear msknoise mskkappa kappa data
    clf(f29); hold on
    set(f29,'Position',[0 800 1440 220],'color',[1 1 1],'name','conclusion'); 
    QRsubset = [1 2 3 5 6]; QRP = QR(QRsubset); QRPname = QRname(QRsubset); % without contrast and old IQR  
    dcolor = lines; dcolor = dcolor([1 3 2],:); 
    for qri = 1:numel(QRP) 
      sp32 = subplot('Position',[ (qri-1) * 1 / (numel(QRP) + 3) + 0.035 , 0.17 , 0.7 / (numel(QRP) + 3) , 0.73],'replace');  
      % prepare data
      for pti = 1:numel(Pt{pi}) % BWP test case
        msknoise(1,pti) = ~cellfun('isempty',strfind(Pt{pi}(pti),{'pn1'})) | ...
                          ~cellfun('isempty',strfind(Pt{pi}(pti),{'pn3'})); %#ok<STRCL1> 
      end
      data{qri}     = cell(1,numel(levels)+1);
      kappa{qri}    = cell(1,numel(levels)+1);
      data{qri}{1}  = tab.(QRP{qri})(:); 
      kappa{qri}{1} = repmat( ktab(:,1), numel(tab.(QRP{qri})) ./ size(ktab,1) ,1); 
      for li = numel(levels):-1:1 
        data{qri}{li+1}  = tabsegt{li}.(QRP{qri})(:);
        kappa{qri}{li+1} = repmat( ktabsegt{li}(:,1), numel(tabsegt{li}.(QRP{qri})) ./ size(ktabsegt{li},1) ,1 );
      end

      
      % scatterplot of the different phantoms
      hold(sp32,'on');
      for di = 1:numel(data{qri})
        hs = scatter(sp32, kappa{qri}{di}, data{qri}{di} , marker(di),'filled');
        set(hs,'MarkerEdgeAlpha',0.5 - 0.2,'MarkerFaceAlpha',0.5 - 0.3, ...
          'MarkerFaceColor', dcolor(di,:),'MarkerEdgeColor', dcolor(di,:));
      end
      pl=plot([0 1],[0 0]); pl.Color = [0.3 0.3 0.3];

      % ########### add lines (maybe not for the final print)
      if 1
        hold on; 
        tab.(QRP{qri}) = tab.(QRP{qri}) - repmat(tab.(QRP{qri})(1,:),numel(bstr),1);
        pl=plot([0 10],[0 0]); pl.Color = [0.3 0.3 0.3]; 
        for pti = 1:numel(Pt{pi}), px = plot( ktab(:,1), tab.(QRP{qri})(:,pti) ); set(px,'Color',min([.8 .8 .8],Pcolor(pti,:)+.5) ); end
        pt = plot( ktab(:,1), mean(tab.(QRP{qri}),2)); set(pt,'LineWidth',1.5,'Color',ones(1,3)*.3,'Marker','x');
      end
      
      % add a legend to the first plot with NCR that only have overestimation
      if qri == 1
        legend(sp32,[{'brain extraction (BE)'},segmentation],'Location','Southwest','box','off');
      end
      

      xlim([0.25,1]), ylim(sp32,markrange); grid(sp32); box on; 
      title(sp32,strrep(QRPname{qri},'_','\_')); xlabel(sp32,'Kappa'); ylabel(sp32,'mark error')
      set(sp32,'ytick',markrange(1):0.5:markrange(2),'yticklabel',num2str((markrange(1):0.5:markrange(2))','%0.2f '), ...
               'xtick',0:0.25:1,'xticklabel',num2str((0:0.25:1)','%0.2f ') );
    end
    
    %% create the final boxplot of all phantoms
    subplot('Position',[ numel(QRP) / (numel(QRP) + 3) + 0.035 ,  0.17 , 0.7 / (numel(QRP) + 3) , 0.73],'replace'); 
    pdata = cell(1,numel(data{1})); 
    for qri = 1:numel(QRP), for pti = 1:numel(data{qri}), pdata{pti} = [pdata{pti}; data{qri}{pti}]; end; end 
    cat_plot_boxplot( pdata , struct('names',{[{'BE'},seg]},'ygrid',0,'ylim',markrange,...
      'groupcolor',dcolor,'style',4,'datasymbol','o','usescatter',1));
    ax=gca; ax.XTickLabelRotation = 0; ax.YGrid = 'on';
    set(ax,'ytick',markrange(1):0.5:markrange(2),'yticklabel',num2str((markrange(1):0.5:markrange(2))','%0.2f '));
    title('Phantoms (all measures)'); xlabel('phantom'); ylabel('mark error')
    hold on, pl=plot([0 10],[0 0]); pl.Color = [0.3 0.3 0.3];
    
    % create the final boxplot of all measures
    subplot('Position',[ (numel(QRP) + 1) / (numel(QRP) + 3) + 0.035 ,  0.17 , 0.7 / (numel(QRP) + 3) , 0.73],'replace'); 
    fdata = cell(1,numel(QRP)); 
    for qri = 1:numel(QRP), fdata{qri} = [fdata{qri}; cell2mat(data{qri}(:)) ]; end 
    cat_plot_boxplot( fdata , struct('names',{QRPname},'ygrid',0,'ylim',markrange,'box',1,'style',4,'datasymbol','o','usescatter',1));
    ax=gca; ax.XTickLabelRotation = 0; ax.YGrid = 'on';
    set(ax,'ytick',markrange(1):0.5:markrange(2),'yticklabel',num2str((markrange(1):0.5:markrange(2))','%0.2f '));
    title('Measures (all phantoms)'); xlabel('measure'); ylabel('mark error')
    hold on, pl=plot([0 10],[0 0]); pl.Color = [0.3 0.3 0.3];
    
    % RMSE plot
    sp30 = subplot('Position',[ (numel(QRP) + 2) / (numel(QRP) + 3) + 0.035 ,  0.17 , 0.7 / (numel(QRP) + 3) , 0.73]);  
    rmseval = [ rms(fdata{1}(:)) , rms(fdata{2}(:)) , rms(fdata{3}(:)) , rms(fdata{4}(:)) , rms(fdata{5}(:)) ]; 
    bh = bar(rmseval(1:end)); 
    ylim([0,1]); xlim([.4 5.6]); xticklabels(QRPname); 
    h = gca; h.YGrid = 'on'; h.XTickLabelRotation = 0; 
    for fi = 1:numel(rmseval)
      text(fi-.38, rmseval(fi) + .03, sprintf('%0.2f',rmseval(fi)),'FontSize',8,'Color',[.0 .2 .4]);
    end  
    ax=gca; ax.XTickLabelRotation = 0; ax.YGrid = 'on'; ax.FontSize = 9.5;
    set(ax,'ytick',markrange(1):0.25:markrange(2),'yticklabel',num2str((markrange(1):0.25:markrange(2))','%0.2f '));
    title(sprintf('RMSE (%0.3f)',mean(rmseval))); xlabel('measure'); ylabel('mark error');


    % caption
    capf = 1.3; 
    f29.Position(4) = f29.Position(4)*capf; ax = findobj(f29.Children,'type','Axes');
    for axi=1:numel(ax), ax(axi).Position(4) = ax(axi).Position(4)/capf; ax(axi).Position(2) = (capf-1) * capf; end
    annotation('textbox',[0 0 1 ((capf-1)/capf + 0.04) ],'FitBoxToText','off','LineStyle','none','String',{['\bfFig. 1 (' ...
        strrep(qafile,'_','\_') '):\rm '...
        'Shown are the deviations of different quality ratings in case of simulated brain extraction (\bfBE\rm; blue circles) and ' ...
        'segmentation errors (with two levels with chnanges of a half (\bfL1\rm - yellow >) and full voxel (\bfL2\rm - red <)) ' ... 
        'in relation to the Kappa rating of the modified segmentation used for quality estimation. ' ...
        'Overall the ratings are quite robust against brain-extraction errors but more susceptible for segmentation errors ' ...
        'of erosion/dilation of the WM and CSF tissue that affect the measures as well as the contrast estimation. ' ...
        'All ratings are defined as school marks (1 grade = 10 rps of the percentage system). ' ... 
        'For the\it noise contrast rating\rm (\bfNCR\rm) a half mark is equal to one BWP noise levels and indicate light (motion) artifacts,' ...
        'whereas for the\it inhomogeneity contrast rating\rm (\bfICR\rm) a 1/4 grade presents a change of 20% of the BWP inhomogeneity. ' ...
        'For\it edge contrast rating\rm (\bfECR\rm) a change of 1 grade correspond to a change of 0.5 mm that is somewhat similar like ' ...
        'Gaussian smoothing of 1.4 mm). ' ...
        'The\it structural image quality rating\rm (\bfSIQR\rm) is the weighted combination of the single measures. ' ...
        '']},'Fontsize',12); 
  
    % save image
    subdir = fullfile(printdir,'Fig1-overview');
    if ~exist(subdir,'dir'), mkdir(subdir); end
    print(f29, '-djpeg', '-r300', fullfile(subdir,sprintf('%s_conclusion_%s',printname,qafile))); 
    
  end
end
fprintf('BWPE done.\n')


function cat_tst_qa_resizeBWP( datadir0, qaversions, rerun ) 
%% BWP resolutions
%  ------------------------------------------------------------------------
%  This function resizes and renames the BWP files to obtain some
%  additional resolution versions for further tests (especially of the 
%  image resolution).
%
%  ------------------------------------------------------------------------
% 
%  Main Theses: 
%   (1) Upsampling should not effect a good (resolution) rating. 
%       Using a pure voxel-based measure (RES) is of course biased, 
%       whereas the edge-based measure (edgeRES) should be not effected. 
%   (2) Downsampling should effect the (resolution) rating 
%       Both resolution measures (RES, edgeRES) should be affected in the 
%       same way.
%   (3) Smoothing an image should be reduce the spatial resolution and 
%       therefore also the overlap to the original images used as ground-
%       truth, e.g., to estimate Kappa or the RMSE.
%       A pure voxel-based measure (RES) is here unaffected but the new 
%       edge-based measure result in worse ratings similar to a reduction 
%       of image resoltion. 
% 
%  Side Aspects: 
%   (1) Using real rather than the ground-truth segmentation should not 
%       largely affect the results.
%   (2) Using simulated randomly distored segmentation should not bias
%       the results as far as these changes are not biased itself. I.e., 
%       (i) the skull-stripping could be inoptimal 
%       (ii) the segmentation could be 
%
%  Data: 
%   (1) BWP       - my standard, but default resolution of 1.0 mm with 
%                   partial volume effects (sufficent results)
%   (2) Tohoku    - only private, average 0.5 mm 
%                   (possible supplement)
%   (3) 28&me     - publicly availbe, should allow 0.5 mm average 
%                   (good supplement)
%   (4) Falk7T    - publicly availbe, should allow <0.5 mm average 
%                   (not suitable now, would need deforamtions)
%   (5) Buchert   - similar to Falk7T (need deformations)
%
%  ------------------------------------------------------------------------


% TODO:
%  * real data test
%  * anisotrophy test > add somehow to figure?
%
%  * add further case?  - Yes, to stabilize the scaling and values 
%                       - bias fields?   - easy, should be save
%                       - noise levels?  - maybe challanging, but good for variance
%                       - bias levels    - similar to noise but less effective 
%                       >> try with noise ... hmmm noise is a problem because it will change :-/ 
%                       >> use bias ABC and 40% to get more values and more stable results 
%   >> test just a noise case to check for problems and maybe add this later 
% +++++++++++ loop data
% +++++++++++ add 28&me
%
  
  cat_io_cprintf([0 0.5 0],'\n\n== Run cat_tst_qa_resizeBWP ==\n') 
  
  % ### data ###
  if ~exist( 'datadir0' , 'var' )
    datadir  = '/Volumes/SG5TB/MRData/202503_QA';
  else
    datadir  = fullfile(datadir0,'BWP'); 
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
  if ~exist( 'rerun', 'var'), rerun = 0; end
  
  resdir   = fullfile(fileparts(datadir),'BWPrestest');
  outdir   = fullfile(fileparts(datadir), '+results',['BWPrestest_' datestr(clock,'YYYYmm')]); 
  if ~exist(fullfile(resdir,'mri'),'dir'), mkdir(fullfile(resdir,'mri')); end  
  
  % testdata from BWP
  Po = {
    fullfile(datadir,'BWPC_HC_T1_pn1_rf020pA_vx100x100x100.nii')
    ...fullfile(datadir,'BWPC_HC_T1_pn1_rf020pB_vx100x100x100.nii')
    ...fullfile(datadir,'BWPC_HC_T1_pn1_rf020pB_vx100x100x100.nii')
    };
  P    = spm_file(Po,'path',resdir); 
  Pp0o = spm_file(Po,'prefix',['mri' filesep 'p0']); 
  Pp0  = spm_file(P ,'prefix',['mri' filesep 'p0']); 
  for fi = 1:numel(Po), copyfile(Po{fi}  ,P{fi}); end
  for fi = 1:numel(Po), copyfile(Pp0o{fi},Pp0{fi}); end
  
  testcase = 1; 
  switch testcase
    case 1 % BWP isotroph
      res = repmat( ( 1:0.25:3 )', 1,3);
      ss  = repmat( ( 0:0.25:3 )', 1,3);
    case 2 % BWP anisotroph
      res = 1:0.25:3; res = [ones(numel(res),2) res'];
      ss  = 0:0.25:3; ss  = [zeros(numel(ss) ,2) ss'];
  end
  
  % qc matlabbatch basic
  qcmatlabbatch{1}.spm.tools.cat.tools.iqe.opts.outdir  = outdir;
  qcmatlabbatch{1}.spm.tools.cat.tools.iqe.opts.verb    = 1;
  for qai = 1%6:numel(qaversions)
    qafile = qaversions{qai};
  
    % main variables
    recalc  = rerun; % force reprocessing of QC
    pi      = 1; % loop variable 
    FN      = {'NCR','ICR','SIQR','IQR','res_RMS','res_ECR'};
    Prs     = cell(numel(P),size(res,1)); 
    Prp     = cell(numel(P),size(ss,1));  
    Pss     = cell(numel(P),size(ss,1));  
    Psp     = cell(numel(P),size(ss,1));  
    Pp0real = cell(numel(P),1); 
    Pp0rs   = cell(numel(P),1); 
    Pp0ss   = cell(numel(P),1); 
    Preal   = cell(numel(P),1); 
    Pproc   = cell(numel(P),1);
    Prsp0   = cell(numel(P),size(res,1));
    Prsp0qc = cell(numel(P),size(res,1));
    Pssp0   = cell(numel(P),size(ss,1));
    Pssp0qc = cell(numel(P),size(ss,1));
    
    
    % print overview of parameters
    fprintf('\nCAT QC resolution measure test: \n')
    fprintf('========================================================================\n')
    fprintf('  Testcase:     %d\n' , testcase)
    fprintf('  File:     %s'  , sprintf('    %s\n',char(P{pi}) ) );
    fprintf('  Outdir:       %s\n' , outdir);
    fprintf('  QC-file:      %s\n' , qafile);
    fprintf('  Resolutions:\n%s'   , sprintf('    %0.2fx%0.2fx%0.2f\n',res') );
    fprintf('  Smoothings: \n%s'   , sprintf('    %0.2fx%0.2fx%0.2f\n',ss')  );
    fprintf('========================================================================\n')
    
    
    %% main loop
    for pi = 1:numel(P)
      [pp,ff,ee] = spm_fileparts(P{pi});
      fprintf('\n  %s:', ff); 
      
      
      % Resmapling resolution changes
      % -----------------------------------------------------------------------
      fprintf('\n    Resolution Reduction:');  
      rxlabel = cell(1,size(res,1));
      for ri = 1:size(res,1)
        rxlabel{ri} = sprintf('r%0.2fx%0.2fx%0.2fmm',res(ri,:)); 
        fprintf('\n      %s',rxlabel{ri});   
       
        % we just use the batch functionality of the cat_vol_resize function
        if recalc>0 || ~exist(Prp{pi,ri},'file')
          %% T1
          clear job; 
          job.data          = P(pi);
          job.restype.res   = res(ri,:); 
          job.interp        = -2005; % linear?,  cubic?-no!,  smooth-linear?,  smooth-cubic?
          job.prefix        = [rxlabel{ri} '_'];
          job.outdir        = {outdir}; 
          job.verb          = 0;
          job.lazy          = 1-recalc;
          rr                = cat_vol_resize(job); % (1) reduction
          Prs{pi,ri}        = rr.res{1};
          job               = rmfield(job,'restype');
          job.data          = Prs(pi,ri);
          job.restype.Pref  = P(pi);
          job.prefix        = ''; 
          cat_vol_resize(job); % (2) reinterpolate
    
          %% P0 
          % this simulates the case of optimal low resolution
          % segmentation and following interpolation.
          job               = rmfield(job,'restype');
          job.data          = Pp0(pi);
          job.restype.res   = res(ri,:); 
          job.interp        = -2002;
          job.prefix        = [rxlabel{ri} '_'];
          rr                = cat_vol_resize(job); % (1) reduction
          
          Prp{pi,ri}        = rr.res{1}; 
          job               = rmfield(job,'restype');
          job.data          = Prp(pi,ri);
          job.restype.Pref  = P(pi);
          job.prefix        = ''; 
          cat_vol_resize(job);
        end
    
        % define both output (the Prs for processing and the Prp and Pp0rs as its simulated and real result for evaluation)
        Prs{pi,ri}        = char(spm_file(P{pi},'path',fullfile(outdir),      'prefix',sprintf('r%0.2fx%0.2fx%0.2fmm_',res(ri,:)))); 
        Prp{pi,ri}        = char(spm_file(P{pi},'path',fullfile(outdir),      'prefix',sprintf('p0r%0.2fx%0.2fx%0.2fmm_',res(ri,:)))); 
        Pp0rs{pi,ri}      = char(spm_file(P{pi},'path',fullfile(outdir,'mri'),'prefix',sprintf('p0r%0.2fx%0.2fx%0.2fmm_',res(ri,:)))); 
      end
      fprintf('\n');
        
    
    
      % smooth
      % -----------------------------------------------------------------------
      fprintf('\n    Smoothing:');
      sxlabel = cell(1,size(ss,1));
      for si = 1:size(ss,1)
        sxlabel{si} = sprintf('s%0.2fx%0.2fx%0.2fmm',ss(si,:)); 
        % define both output (the Pss for processing and the Pp0ss as its result for evaluation)
        Pss(pi,si)  = spm_file(P(pi)  ,'path',outdir,'prefix',[sxlabel{si} '_']);
        Psp(pi,si)  = spm_file(Pp0(pi),'path',outdir,'prefix',[sxlabel{si} '_']);
        Pp0ss{pi,si} = char(spm_file(P{pi},'path',fullfile(outdir,'mri'),'prefix',sprintf('p0s%0.2fx%0.2fx%0.2fmm_',ss(si,:)))); 
    
        fprintf('\n      %s',sxlabel{si});   
        if ~exist(Pss{pi,si},'file')
          spm_smooth(P{pi},Pss{pi,si},ss(si,:)); 
        end
        % How does Kappa changes if we just smooth the original segmentation. 
        % However, this is not realistic because there should be no interpolated 
        % segmentation and we assume that intermediate values are from the PVE.
        if ~exist(Psp{pi,si},'file')
          spm_smooth(Pp0{pi},Psp{pi,si},ss(si,:)); 
        end
       
      end
      fprintf('\n');
    
    
      % run preprocessing for real segmentations 
      fprintf('\n    CAT Preprocessing:');  
      Preal{pi}   = [ P{pi} , Prs(pi,:) , Pss(pi,:) ];
      Pproc{pi}   = [ P{pi} , Prs(pi,:) , Pss(pi,:) ];
      Pp0real{pi} = spm_file(Pproc{pi},'prefix',['mri' filesep 'p0']); 
      for fi = numel(Pp0real{pi}):-1:1
        if exist(Pp0real{pi}{fi},'file'), Pproc{pi}(fi) = []; end 
      end
      if ~isempty(Pproc{pi})
        CATpreprocessing4qc;
        matlabbatch{1}.spm.tools.cat.estwrite.data = Pproc{pi}';
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.lazy = 1; 
        spm_jobman('run',matlabbatch); 
        fprintf(' .. done.\n');
      else
        fprintf(' .. is not required.\n');
      end
      
    
      % Kappa estimation
      % - estimate kappa resolution
      [~,valr] = eva_vol_calcKappa([Pp0;Pp0rs(pi,2:end)'],Pp0,struct('recalc',1 + recalc,'realign',0)); 
      % - estimate kappa smooth
      [~,vals] = eva_vol_calcKappa([Pp0;Pp0ss(pi,2:end)'],Pp0,struct('recalc',1 + recalc,'realign',0)); 
      
    
      %% Image quality estimation
      stepsize = 1; % this is only for fast test and should be 1 by default 
      
      % run with ground truth segmentation 
      prefixgt                                              = ['qc_gt_' qafile];
      qcmatlabbatch{1}.spm.tools.cat.tools.iqe.images       = Prs(pi,1:stepsize:end)';
      qcmatlabbatch{1}.spm.tools.cat.tools.iqe.model.catp0  = Pp0; 
      qcmatlabbatch{1}.spm.tools.cat.tools.iqe.opts.prefix  = prefixgt;
      qcmatlabbatch{1}.spm.tools.cat.tools.iqe.opts.verb    = 2;
      qcmatlabbatch{1}.spm.tools.cat.tools.iqe.opts.rerun   = 1;
      qcmatlabbatch{1}.spm.tools.cat.tools.iqe.opts.version = qafile;
      Prsqc = spm_file( Prs(pi,:) ,'path',fullfile(outdir,'report'),'prefix',prefixgt,'ext','.xml');
      spm_jobman('run',qcmatlabbatch); 
      xmlr = cat_io_xml(Prsqc); clear Xr;
      for fni=1:numel(FN)
        for xi=1:numel(xmlr)
          try 
            Xr.(FN{fni})(xi) = xmlr(xi).qualityratings.(FN{fni}); 
          catch
            Xr.(FN{fni})(xi) = nan;
          end
        end
      end
      
      %% run with ground truth segmentation 
      qcmatlabbatch{1}.spm.tools.cat.tools.iqe.images       = Psp(pi,1:stepsize:end)';
      qcmatlabbatch{1}.spm.tools.cat.tools.iqe.model.catp0  = Pp0; 
      Pspqc = spm_file( Psp(pi,:) ,'path',fullfile(outdir,'report'),'prefix',prefixgt,'ext','.xml');
      spm_jobman('run',qcmatlabbatch); 
      xmls = cat_io_xml(Pspqc); clear Xs
      for fni=1:numel(FN)
        for xi=1:numel(xmls)
          try 
            Xs.(FN{fni})(xi) = xmls(xi).qualityratings.(FN{fni});
          catch 
            Xs.(FN{fni})(xi) = nan;
          end
        end
      end
      
    
    
      %% run with real segmentation (resolution)
      prefixp0      = ['qc_p0_' qafile];
      Prs(pi,1)     = Pss(pi,1); % ############# the image is somehow modified 
      Prsp0(pi,:)   = spm_file( Prs(pi,:) ,'path',fullfile(outdir,'mri')   ,'prefix','p0');
      Prsp0qc(pi,:) = spm_file( Prs(pi,:) ,'path',fullfile(outdir,'report'),'prefix',prefixp0,'ext','.xml');
      qcmatlabbatch{1}.spm.tools.cat.tools.iqe.images       = Prs(pi,1:stepsize:end)';      
      qcmatlabbatch{1}.spm.tools.cat.tools.iqe.model.catp0  = Prsp0(pi,1:stepsize:end)';  
      qcmatlabbatch{1}.spm.tools.cat.tools.iqe.opts.prefix  = prefixp0;
      spm_jobman('run',qcmatlabbatch); 
      xmlr2 = cat_io_xml(Prsp0qc); clear Xr2
      for fni=1:numel(FN)
        for xi=1:numel(xmlr2)
          try
            Xr2.(FN{fni})(xi) = xmlr2(xi).qualityratings.(FN{fni}); 
          catch
            Xr2.(FN{fni})(xi) = nan; 
          end
        end
      end
      
      %% run with real segmentation (smoothing)
      Pssp0(pi,:)   = spm_file( Pss(pi,:) ,'path',fullfile(outdir,'mri')   ,'prefix','p0');
      Pssp0qc(pi,:) = spm_file( Pss(pi,:) ,'path',fullfile(outdir,'report'),'prefix',prefixp0,'ext','.xml');
      qcmatlabbatch{1}.spm.tools.cat.tools.iqe.images       = Pss(pi,1:stepsize:end)';
      qcmatlabbatch{1}.spm.tools.cat.tools.iqe.model.catp0  = Pssp0(pi,1:stepsize:end)';
      spm_jobman('run',qcmatlabbatch); 
      xmls2 = cat_io_xml(Pssp0qc(1:end)'); clear Xs2
      for fni=1:numel(FN)
        for xi=1:numel(xmls2)
          try
            Xs2.(FN{fni})(xi) = xmls2(xi).qualityratings.(FN{fni});
          catch
            Xs2.(FN{fni})(xi) = nan;
          end
        end
      end
  
  
  
      %% run interpolation test
      % interpolate image from 1 to 0.9:0.1:0.5 mm
      % run kappa
      % run QC
      % eval
     
    
    
    % +++++ estimate correlations ... just some code fragment for fast adaptions    
      if 0
        %%
        %ECR  = [.2700 .2568 .2406 .2292 0.2182 .2081 .2017 .1980 0.1965]; 
        ECR  = [.2901 .2756 .2480 .2303 0.2295 .2161 .2080 .2019 0.1947]; 
        ECRs = [.2699 .2483 .2288 .2096 .1922 .1810 .1731];
        K    = cell2mat([vals(2,5);valr(3:end-2,5)]); 
        Ks   = cell2mat(vals(2:end-2,5));
        EX   = @(x) log10( (max(0,x - 0.15)  )  ); 
        figure(3); plot(ECR    ,K,'x-'); hold on; plot(ECRs    ,Ks,'o-'); hold off; 
        figure(4); plot(EX(ECR),K,'x-'); hold on; plot(EX(ECRs),Ks,'o-'); hold off;
      end
      
    
      
      
      %% create some figure with subfigure#
      measures = {'res_RMS','res_ECR','NCR','ICR','IQR','SIQR'}; % print measures (see loop for limited tests) 
      mr       = [1 7]; % limitation of marks
      for mi = 1:numel(measures) % real measures and side effects ... 1:numel(measures) 
        figid   = 201; 
        measure = measures{mi};
      
        %% corr
        [rhor,pvalr] = corr( Xr.(measure)' , cell2mat(valr(2:end-2,5)')' );
        [rhos,pvals] = corr( Xs.(measure)' , cell2mat(vals(2:end-2,5)')' );
      
    
        % ---------------------------------------------------------------------
        % Figure 1 with just some bars for kappa (1st row), measure with ground
        % truth (2nd row), and measure with real segmentation (3rd row).
        % ---------------------------------------------------------------------
        figure(figid); fp = get(gcf,'Position'); fp(4) = 700; set(gcf,'Position',fp,'name',sprintf('%s-%s',qafile,measure)); 
    
        % plot kappa values
        figure(figid); subplot(3,2,1);  bar(cell2mat(valr(2:end-2,5))); grid on; ylim([0.85 1]); xlim([0 numel(rxlabel)+1])
        title('Resolution (processing quality)'); ylabel('Kappa')
        set(gca,'XTick',1:size(res,1),'xTickLabel',rxlabel,'TickLabelInterpreter','none','XTickLabelRotation',90);
        figure(figid); subplot(3,2,2);  bar(cell2mat(vals(2:end-2,5))); ylim([0.9 1]); xlim([0 numel(sxlabel)+1]); grid on; 
        title('Smoothing (processing quality)'); ylabel('Kappa')
        set(gca,'XTick',1:size(ss,1),'xTickLabel',sxlabel,'TickLabelInterpreter','none','XTickLabelRotation',90);
        
        % plot measures for ground turth segmentation 
        figure(figid); subplot(3,2,3);  bar(Xr.(measure)); ylim(mr); grid on; 
        title(sprintf('QC measure (%s; GT)',strrep(measure,'_',' '))); 
        ylabel('grad') %ylabel(sprintf('QC grad (%s;r=0.2f)',measure,kmr));
        set(gca,'XTick',1:size(res,1),'xTickLabel',rxlabel,'TickLabelInterpreter','none','XTickLabelRotation',90);
        figure(figid); subplot(3,2,4);  bar(Xs.(measure)); ylim(mr); grid on; 
        title(sprintf('QC measure (%s; GT)',strrep(measure,'_',' '))); 
        ylabel('grad') %ylabel(sprintf('QC grad (%s;r=0.2f)',measure,kmr));
        set(gca,'XTick',1:size(ss,1),'xTickLabel',sxlabel,'TickLabelInterpreter','none','XTickLabelRotation',90);
       
        % plot measures for real segmentation
        figure(figid); subplot(3,2,5);  bar(Xr2.(measure)); ylim(mr); grid on
        title(sprintf('QC measure (%s; DS); rho=%0.3f, p=%0.0e)',strrep(measure,'_',' '),rhor,pvalr));
        ylabel('grad') %ylabel(sprintf('QC grad (%s;r=0.2f)',measure,kmr));
        set(gca,'XTick',1:size(res,1),'xTickLabel',rxlabel,'TickLabelInterpreter','none','XTickLabelRotation',90);
        figure(figid); subplot(3,2,6);  bar(Xs2.(measure)); ylim(mr); grid on
        title(sprintf('QC measure (%s; DS; rho=%0.3f, p=%0.0e)',strrep(measure,'_',' '),rhos,pvals));
        ylabel('grad') %ylabel(sprintf('QC grad (%s;r=0.2f)',measure,kmr));
        set(gca,'XTick',1:size(ss,1),'xTickLabel',sxlabel,'TickLabelInterpreter','none','XTickLabelRotation',90);
        
        % save figures
        % update this here to have up to date value in case of batch mode 
        printoutdir   = fullfile(outdir,'results',datestr(clock,'YYYYmm')); 
        if ~exist(printoutdir,'dir'), mkdir(printoutdir); end
        [~,ff]        = spm_fileparts(P{pi}); 
        printname     = sprintf('testcase%d_%s',testcase,ff);
    
        % final print of the figure
        print(figid, '-djpeg', '-r300', fullfile(printoutdir,sprintf('%s_bars-%s_%s',printname,measure,qafile))); 
        
    
    
        %% ---------------------------------------------------------------------
        % Second figure with kappa vs. measure. 
        % ---------------------------------------------------------------------
        plotGT    = 0; % plot only real data (2 lines) or also GT cases (4 lines)
        plotIdeal = 0; % plot some ideal line 
        mark2rps  = @(mark) min(100,max(0,105 - mark*10)) + isnan(mark).*mark;
        figid     = 202; 
        colors    = lines(4); 
        
        fh = figure(figid); fh.Name = measure; fh.Position(3:4) = [260 260]; clf; hold on; 
        % line plots
        if plotGT
          plot( Xr.(measure)  , cell2mat(valr(2:end-2,5)') ,'Marker','<', 'Color', colors(3,:),'Markerfacecolor',[1 1 1]); 
          plot( Xs.(measure)  , cell2mat(vals(2:end-2,5)') ,'Marker','>', 'Color', colors(4,:),'Markerfacecolor',[1 1 1]); 
        end
        plot( Xr2.(measure) , cell2mat(valr(2:end-2,5)') ,'Marker','o', 'Color', colors(1,:),'Markerfacecolor',[1 1 1]);  
        plot( Xs2.(measure) , cell2mat(vals(2:end-2,5)') ,'Marker','s', 'Color', colors(2,:),'Markerfacecolor',[1 1 1]); 
        if plotIdeal
          plot( 2* mean(res.^2,2).^0.5 , cell2mat(valr(2:end-2,5)') ,'Marker','^', 'Color', [0.6 0.7 0.8],'Markerfacecolor',[1 1 1],'LineStyle','--'); 
        end
    
        % Mark the full resolution values with filled symbols (e.g. 1,2,3 mm) - this is hard coded
        if numel(Xr2.(measure))==9 
          for ri = [1,5,9]
            if plotGT
              plot( Xr.(measure)(ri) , cell2mat(valr(ri+1,5)') ,'Marker','<', 'Color', colors(3,:), 'MarkerFaceColor', colors(3,:) ); 
            end
            if plotIdeal
              plot( 2* mean(res(ri,:).^2,2).^0.5 , cell2mat(valr(ri+1,5)') ,'Marker','^', 'Color', [0.6 0.7 0.8], 'MarkerFaceColor', [0.6 0.7 0.8]); 
            end
            plot( Xr2.(measure)(ri)  , cell2mat(valr(ri+1,5)') ,'Marker','o', 'Color', colors(1,:), 'MarkerFaceColor', colors(1,:) ); 
          end
        end
        if numel(Xs2.(measure))==13
          for si = [1,5,9,13]
            if plotGT
              plot( Xs.(measure)(si) , cell2mat(vals(si+1,5)') ,'Marker','>', 'Color', colors(4,:), 'MarkerFaceColor', colors(4,:) );
            end
            plot( Xs2.(measure)(si)  , cell2mat(vals(si+1,5)') ,'Marker','s', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:) );
          end
          %plot( Xr2.(measure)(1) , cell2mat(valr(2,5)') ,'Marker','s', 'Color', colors(2,:), 'MarkerFaceColor', colors(2,:) ); 
        end
    
        % add legend
        hold off; box on; grid on;  ylim([0.795 1.005]); xlim([0.9 8.1]); 
        ylabel('Kappa'); xlabel(strrep([measure ' (grade)'],'_','\_')); 
        ax = gca; ax.XTick = [1:8]; ax.YTick = [.8:0.02:1.0]; ax.Position(3) = .83;
        if plotGT
          leg = {'ResamplingGT','SmoothingGT','Resampling','Smoothing'};
        else
          leg = {sprintf('Resicampling (\\it{}rho\\rm{}=%0.3f, \\it{}p\\rm{}=%0.0e)',rhor,pvalr), ...
                 sprintf('Smoothing (\\it{}rho\\rm{}=%0.3f, \\it{}p\\rm{}=%0.0e)',rhos,pvals)};
        end
        if plotIdeal
          leg = [leg,{'Expected for lower res.'}];
        end
        legend(leg,'Location','South')
        title('Quantification of anatomical resolution');
    
        % add further in figure text to sign out the marked elements - hard coded 
        if plotGT
          if numel(Xr.(measure))==9 
            % resolution data labels
            text(Xr.(measure)(1)   , valr{2,5}     , sprintf('\\leftarrow %s','1 mm') , 'Color', colors(3,:));
            text(Xr.(measure)(5)   , valr{6,5}     , sprintf('%s\\rightarrow','2 mm') , 'Color', colors(3,:),'HorizontalAlignment','right');
            text(Xr.(measure)(end) , valr{end-2,5} , sprintf('%s\\rightarrow','3 mm') , 'Color', colors(3,:),'HorizontalAlignment','right');
          end
          if numel(Xs.(measure))==13
            % smoothing data labels 
            text(Xs.(measure)(1)   , vals{2,5}     , sprintf('\\leftarrow %s','0 mm') , 'Color', colors(4,:));
            text(Xs.(measure)(5)   , vals{6,5}     , sprintf('\\leftarrow %s','1 mm') , 'Color', colors(4,:));
            text(Xs.(measure)(9)   , vals{10,5}    , sprintf('\\leftarrow %s','2 mm') , 'Color', colors(4,:));
            text(Xs.(measure)(13)  , vals{end-2,5} , sprintf('\\leftarrow %s','3 mm') , 'Color', colors(4,:));
          end
        end
        if numel(Xr2.(measure))==9 && plotIdeal
          % resolution data labels
          text(2* mean(res(1,:).^2,2).^0.5 , valr{2,5}     , sprintf('\\leftarrow %s','1 mm') , 'Color', [0.6 0.6 .8]);
          text(2* mean(res(round(end/2),:).^2,2).^0.5 , valr{6,5}     , sprintf('\\leftarrow %s','2 mm') , 'Color', [0.6 0.6 .8]);
          text(2* mean(res(end,:).^2,2).^0.5 , valr{end-2,5} , sprintf('\\leftarrow %s','3 mm') , 'Color', [0.6 0.6 .8]);
        end
        if numel(Xr2.(measure))==9 
          % resolution data labels
          text(Xr2.(measure)(1)   , valr{2,5}     , sprintf('%s\\rightarrow','1 mm') , 'Color', colors(1,:),'HorizontalAlignment','right');
          text(Xr2.(measure)(5)   , valr{6,5}     , sprintf('%s\\rightarrow','2 mm') , 'Color', colors(1,:),'HorizontalAlignment','right');
          text(Xr2.(measure)(end) , valr{end-2,5} , sprintf('%s\\rightarrow','3 mm') , 'Color', colors(1,:),'HorizontalAlignment','right');
        end
        if numel(Xs2.(measure))==13
          % smoothing data labels 
          text(Xs2.(measure)(1)   , vals{2,5}     , sprintf('\\leftarrow %s','0 mm') , 'Color', colors(2,:));
          text(Xs2.(measure)(5)   , vals{6,5}     , sprintf('\\leftarrow %s','1 mm') , 'Color', colors(2,:));
          text(Xs2.(measure)(9)   , vals{10,5}    , sprintf('\\leftarrow %s','2 mm') , 'Color', colors(2,:));
          text(Xs2.(measure)(13)  , vals{end-2,5} , sprintf('\\leftarrow %s','3 mm') , 'Color', colors(2,:));
        end
      
        % final print of the figure
        print(fh, '-djpeg', '-r300', fullfile(printoutdir,sprintf('%s_Kappa-%s_%s',printname,measure,qafile))); 
      end
    end
  end
end
% == this is the CAT job for preprocessing in developer mode! == 
function matlabbatch = catjob(nproc)
% There are some minor changes to avoid reprocessing and output of unused 
% volumes and especially surfaces.  

  if ~exist('nproc','var'), nproc = 0; end

  matlabbatch{1}.spm.tools.cat.estwrite.data = '<UNDEFINED>';
  matlabbatch{1}.spm.tools.cat.estwrite.data_wmh = {''};
  matlabbatch{1}.spm.tools.cat.estwrite.nproc = nproc;
  matlabbatch{1}.spm.tools.cat.estwrite.useprior = '';
  matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm = {'/Users/dahnke/Documents/MATLAB/spm12/tpm/TPM.nii'};
  matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg = 'mni';
  matlabbatch{1}.spm.tools.cat.estwrite.opts.ngaus = [1 1 2 3 4 2];
  matlabbatch{1}.spm.tools.cat.estwrite.opts.warpreg = [0 0.001 0.5 0.05 0.2];
  matlabbatch{1}.spm.tools.cat.estwrite.opts.bias.biasstr = 0.5;
  matlabbatch{1}.spm.tools.cat.estwrite.opts.acc.accstr = 0.5;
  matlabbatch{1}.spm.tools.cat.estwrite.opts.redspmres = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.restypes.optimal = [1 0.3];
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.setCOM = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.APP = 1070;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.affmod = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.NCstr = -Inf;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.spm_kamap = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.LASstr = 0.5;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.LASmyostr = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.gcutstr = 2;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.cleanupstr = 0.5;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.BVCstr = 0.5;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.WMHC = 2;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.SLC = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.mrf = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.T1 = ...
    {'/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_MNI152NLin2009cAsym/T1.nii'};
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.brainmask = ...
    {'/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_MNI152NLin2009cAsym/brainmask.nii'};
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.cat12atlas = ...
    {'/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_MNI152NLin2009cAsym/cat.nii'};
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.darteltpm = ...
    {'/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_MNI152NLin2009cAsym/Template_1_Dartel.nii'};
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shootingtpm = ...
    {'/Users/dahnke/Documents/MATLAB/spm12/toolbox/cat12/templates_MNI152NLin2009cAsym/Template_0_GS.nii'};
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.regstr = 0.5;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.bb = 12;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.vox = 1.5;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.pbtres = 0.5;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.pbtmethod = 'pbt2x';
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.SRP = 22;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.reduce_mesh = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.vdist = 2;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.scale_cortex = 0.7;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.add_parahipp = 0.1;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.close_parahipp = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.experimental = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.new_release = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.lazy = 1; % ########## changed
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.ignoreErrors = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.verb = 2;
  matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.print = 2;
  matlabbatch{1}.spm.tools.cat.estwrite.output.BIDS.BIDSno = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.surf_measures = 3;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 0; % ########## changed
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = 0; % ########## changed
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.thalamus = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.suit = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ibsr = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.aal3 = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.mori = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.anatomy3 = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.julichbrain = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_100Parcels_17Networks_order = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_200Parcels_17Networks_order = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_400Parcels_17Networks_order = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_600Parcels_17Networks_order = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ownatlas = {''};
  matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.GM.warped = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod    = 0; % ########## changed
  matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.WM.warped = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod    = 0; % ########## changed
  matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.native = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.warped = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.mod = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.dartel = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ct.native = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ct.warped = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.ct.dartel = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.pp.native = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.pp.warped = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.pp.dartel = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.native = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.warped = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.mod = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.dartel = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.SL.native = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.SL.warped = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.SL.mod = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.SL.dartel = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.native = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.warped = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.mod = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.dartel = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.native = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.warped = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.dartel = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.label.native = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.output.label.warped = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.label.dartel = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.labelnative = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.output.bias.native = 1;
  matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.bias.dartel = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.las.native = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.las.warped = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.las.dartel = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 0;
  matlabbatch{1}.spm.tools.cat.estwrite.output.warps = [0 0]; % ########## changed
  matlabbatch{1}.spm.tools.cat.estwrite.output.rmat = 0;
end


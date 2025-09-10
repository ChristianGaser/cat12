function cat_tst_qa_iqrRMS( datadir, qaversions, segment, fasttest)
% This script is about how to set up the average rating function. 
%
% - Correlation between Volume/Kappa and SIQR on BWP and MRART dataset 
%
% - Kappa is defined in relation to the zero-artifact case. 
%
% - My original orientation was: 
%   Image with severe artifacts should be rated with worse than 4 what 
%   mostly effected the defintion of the NCR estimation. 
%
% - Averaging with higher powers should reduce the impact of single good ratings.
%   This is important as the average score should be useful to detect outliers
%   that might cause issues in preprocessing.
%
% - In early versions of the QC we (had to) focused on NCR and res_RMS, 
%   while the ICR caused more problems then it was useful 
%   (e.g. highfield has strong bias but also better SNR)
%   As the additional later measures (FEC and res_ECR) where also sensitive 
%   for motion artifacts lower power values become more useful. 
%
% - We also have to consider that we only have two evaluation datasets and
%   the BWP and MRART that both have some limiations:
%     BWP:    artificial with 1/3 interpolated data 
%             (what troubles the res_RMS measure), and is overrepresented) 
%     MRART:  only one protocol (one resolution)
%
% - It is also possible to weight the tissue classes as WM and GM are much 
%   more relevant for processing and evaluation.
%
% - Different segmentations will also influence the result. 
%   Especially on the BWP SPM performend quite pour in some cases that 
%   are not fully representive for real data.  
%
% - The results should therefore not be overinterpretated.
%

cat_io_cprintf([0 0.5 0],'\n\n== Run cat_tst_qa_iqrRMS ==\n') 

%#ok<*AGROW>

if ~exist( 'fasttest', 'var'), fasttest = 0; end
fast            = {'full','fast'}; 

opt.type        = '-depsc';
opt.res         = '-r300'; 
opt.dpi         = 90; 
opt.closefig    = 0;
opt.maindir     = datadir; 
opt.resdir      = fullfile(opt.maindir,'+results','SIQRweighting'); 
opt.resdirBWP   = fullfile(opt.maindir,'+results', ...
  sprintf('%s_%s_%s', 'BWPmain', fast{fasttest+1}, '202508' )); %f datestr(clock,'YYYYmm')) );
opt.resdirART   = fullfile(opt.maindir,'+results', ...
  sprintf('%s_%s_%s', 'MR-ART', fast{fasttest+1}, '202508' )); %f datestr(clock,'YYYYmm')) );

% get BWP and MRART files
P.BWPfiles = [
    cat_vol_findfiles( fullfile( datadir, 'BWP'),   'BWP*.nii'  ,struct('depth',0));
    cat_vol_findfiles( fullfile( datadir, 'BWPr' ), 'rBWP*.nii' ,struct('depth',0));
    cat_vol_findfiles( fullfile( datadir, 'BWPr' ), 'irBWP*.nii',struct('depth',0))];
MRARTdir = 'ds004173-download'; 
P.ARTfiles = cat_vol_findfiles( fullfile( datadir, MRARTdir) , 'sub*.nii',struct('depth',3));


%% find XML files 
qai = 1; mi = 1; qcsi =1; 
segmentname  = cat_io_strrep(segment,{'CAT','SPM'},{'','spm_'}); 
BWPkappaname = cat_io_strrep(segment,{'CAT','SPM'},{'CAT12','SPM12'});  

% Defintion of quality ratings use for averaging
% It is possible to test different combinations of the quality measures
% However, for the article we focus only on the final one that use all QMs.
if 1 % ... not working need another cell layer
  QCs = {
    {
      {'NCR','res_RMS'};                                      % 2 para - old without ICR
      {'NCR','ICR','res_RMS'};                                % 3 para - old with ICR
      {'NCR','ICR','res_RMS','res_ECR','FEC'};  
      {'NCR','res_RMS','res_ECR','FEC'};  
    };
  };
elseif true
   QCs = {
    {
      {'NCR','res_RMS','res_ECR','FEC'};  
    }
  };
else
  QCs = {
    {
      {'ICR','res_RMS','res_ECR'}; {'NCR','res_RMS','res_ECR'}; {'NCR','ICR','res_ECR'}; {'NCR','ICR','res_RMS'};
      {'ICR','res_RMS','FEC'}; {'NCR','res_RMS','FEC'}; {'NCR','ICR','FEC'}; {'NCR','ICR','res_RMS'};
      {'ICR','FEC','res_ECR'}; {'NCR','FEC','res_ECR'}; {'NCR','ICR','res_ECR'}; {'NCR','ICR','FEC'}; 
      {'FEC','res_RMS','res_ECR'}; {'NCR','res_RMS','res_ECR'}; {'NCR','FEC','res_ECR'}; {'NCR','FEC','res_RMS'};
      {'ICR','res_RMS','res_ECR'}; {'FEC','res_RMS','res_ECR'}; {'FEC','ICR','res_ECR'}; {'FEC','ICR','res_RMS'};    
  
      {'NCR','ICR','res_RMS','res_ECR'}; {'NCR','ICR','res_RMS','FEC'}; {'NCR','ICR','FEC','res_ECR'}; 
      {'NCR','FEC','res_RMS','res_ECR'}; {'FEC','ICR','res_RMS','res_ECR'};  
      {'NCR','ICR','res_RMS','res_ECR','FEC'}; 

      {'NCR','res_RMS','res_ECR','FEC'};  
    }
  };
end

bwptestset = 1; 
tissue    = {'CSF','GM','WM','avg'}; % use the average also the CSF is quite bad
%CM        = cat_io_colormaps('set1',numel(QCs{qcsi})); CM(end,:) = [0 0 0];
MK        = '<>^vso+x';
%corrtype  = 'Spearman'; % use pearson to focus on linear relationship 
corrtype  = 'Pearson'; 


%%
clear kappa X QC SIQR rkappa name rtissue; 
for mi = 1:numel(segment)
  %% loading QC ratings
  for qai = 1:numel(qaversions)
    P.BWPfilesX{mi,qai} = spm_file(P.BWPfiles,'prefix',['report' filesep qaversions{qai} '_' segmentname{mi}],'ext','xml');
    %P.ARTfilesX{mi,qai} = cat_vol_findfiles( fullfile( datadir, 'ds004173-download','derivatives', 'CAT12.9') , ...
    %  sprintf('%s_sub*.xml',qaversions{qai}),struct('depth',3));
    for fi = 1:numel(P.ARTfiles) % 'derivatives' filesep 'CAT12.9' filesep 
      P.ARTfilesX{mi,qai}{fi} = spm_file(P.ARTfiles{fi},'prefix',[qaversions{qai} '_' segmentname{mi}],'ext','xml', ...
        'path', strrep(spm_file(P.ARTfiles{fi},'path'),MRARTdir,fullfile(MRARTdir,'derivatives','CAT12.9') ));
    end
  
    % load XML
    X{1,mi,qai} = cat_io_xml(P.BWPfilesX{mi,qai}); 
    X{2,mi,qai} = cat_io_xml(P.ARTfilesX{mi,qai}); 
  end
  
  %% BWP kappa (search for the results from the prevous script)
  %P.BWPkappamat{mi} = fullfile(opt.resdirBWP,sprintf('bwp_kappa_NIR%d_%s.mat',numel(P.BWPfilesX{mi}),BWPkappaname{mi}));
  tmp = cat_vol_findfiles(fileparts(fullfile(opt.resdirBWP,'BWPmain_full_202508')), sprintf('bwp_kappa_NIR%d_%s.mat',numel(P.BWPfilesX{mi}),BWPkappaname{mi}));
  if ~isempty(tmp), P.BWPkappamat{mi} = tmp{1}; else, P.BWPkappamat{mi} = ''; end
  if ~exist(P.BWPkappamat{mi},'file')
    error('Cannot find BWP kappa estimates. Run  cat_tst_qa_bwpmaintest  in  cat_tst_qa_main.')
  else
    S = load(P.BWPkappamat{mi}); 
  end
  kappa{1} = cell2mat(S.val(2:end-2,2:5));

  % only BWP test dataset (see cat_tst_qa_bwpmaintest)
  if bwptestset
    test = ~( mod(1:numel(P.BWPfiles),8)' < 4); 
    P.BWPfiles(test) = []; 
    for  qai = 1:numel(qaversions)
      P.BWPfilesX{mi,qai}(test) = []; 
      X{1,mi,qai}(test) = []; 
    end
    kappa{1}(test,:) = []; 
  end
 

  %% MRART kappa
  di = 2; % this is dataset two
  for si = 1:numel(X{di,1,1})
    subID{di}{si} = X{di,1,1}(si).filedata.file(5:10);
    BL{di}(si)    = contains(X{di,1}(si).filedata.file,'standard'); 
    P.ARTfilesp0{mi}{si} = X{di,1,1}(si).filedata.Fp0; %, ...
     % {'ds004173-download','mri'},{fullfile('ds004173-download','derivatives','CAT12.9'),''}); 
  end
  tmp = cat_vol_findfiles(fileparts(opt.resdirART), sprintf('mrart_kappa_NIR%d_%s.mat',numel(P.ARTfilesX{mi}),BWPkappaname{mi})); 
  if ~isempty(tmp), P.ARTkappamat{mi} = tmp{1}; else, P.ARTkappamat{mi} = ''; end
  if ~exist(P.ARTkappamat{mi},'file')
    P.ARTkappamat{mi} = fullfile(opt.resdirART,sprintf('mrart_kappa_NIR%d_%s.mat',numel(P.ARTfilesX{mi}),BWPkappaname{mi})); 
    for si = 1:numel( P.ARTfilesp0{mi}) % only first dim!
      BLsiID(si) = find( contains( subID{di} , subID{di}{si}) & BL{di} == 1);
    end

    [~,val] = eva_vol_calcKappa(P.ARTfilesp0{mi} ,P.ARTfilesp0{mi}(BLsiID) ,struct('recalc',0,'realign',2,'realignres',1));
    if ~exist(opt.resdirART,'dir'), mkdir(opt.resdirART); end
    save(P.ARTkappamat{mi},'val'); 
    kappa{2} = real(cell2mat(val(2:end-2,2:5)));
  else
    S = load(P.ARTkappamat{mi}); 
    kappa{2} = real(cell2mat(S.val(2:end-2,2:5)));
  end


  %% estimate correlations
  clear rkappa rtissue 
  for qai = 1:numel(qaversions)
    for qcsi = 1:numel(QCs) 
      %%
      clear SIQR name rBV; 
      for di = 1:size(X,1) % datasets
        for qi = 1:numel(QCs{qcsi}) % quality measure set
          for si = 1:numel(X{di,mi,qai}) % subject
            for qci = 1:numel(QCs{qcsi}{qi}) % quality measure
              try
                QC{di,qi}(si,qci) = X{di,mi,qai}(si).qualityratings.(QCs{qcsi}{qi}{qci}); 
              catch
                QC{di,qi}(si,qci) = nan; 
              end
              if qci == 1 
                name{qi} = sprintf('%d %s',qi,QCs{qcsi}{qi}{qci});
              else
                name{qi} = sprintf('%s+%s',name{qi},QCs{qcsi}{qi}{qci}); 
              end
            end
          end
          
          % estimate general SIQR measure and the correlation to Kappa
          % * genSIQR .. subjectwise SIQR values ( {dataset,QMset}( subject , averagemode ) 
          % * gkappa  .. correlation kappa-generalSIQR ( {segmentation,QCversion,dataset}(QMset,rms-potenial,tissue)      
          % * rkappa  .. correlation kappa-SIQR  ( {segmentation,QCversion,dataset}(QMset,rms-potenial,tissue)  
          genSIQRmod = {'median','mean','rms^2','rms^4','rms^8','rms^{16}','max'}; 
          genSIQR{di,qi}(:,1) = cat_stat_nanmedian( QC{di,qi} ,2); 
          genSIQR{di,qi}(:,2) = cat_stat_nanmean( QC{di,qi},2); 
          genSIQR{di,qi}(:,3) = cat_stat_nanmean( QC{di,qi}.^2 ,2).^(1/2); 
          genSIQR{di,qi}(:,4) = cat_stat_nanmean( QC{di,qi}.^4 ,2).^(1/4); 
          genSIQR{di,qi}(:,5) = cat_stat_nanmean( QC{di,qi}.^8 ,2).^(1/8); 
          genSIQR{di,qi}(:,6) = cat_stat_nanmean( QC{di,qi}.^16,2).^(1/16); 
          genSIQR{di,qi}(:,7) = max(  QC{di,qi},[],2); 
          % estimate RMS-based SIQR measure and the correlation to Kappa
          for rmsi = 1:16
            SIQR{di,qi}(:,rmsi) = cat_stat_nanmean( QC{di,qi}.^rmsi ,2).^(1/rmsi); 
            if bwptestset && di==1
              % to have similar values as before we use here directly the average as in cat_tst_qa_bwpmaintest
              rkappa{mi,qai,di}(qi,rmsi,1:size(kappa{1},2))   = repmat( ...
                corr(cat_stat_nanmean(kappa{di}(:,1:3),2), SIQR{di,qi}(:,rmsi), 'Type',corrtype) , 1, size(kappa{di},2));
              if rmsi <= size(genSIQR{di,qi},2)
                gkappa{mi,qai,di}(qi,rmsi,1:size(kappa{1},2)) = repmat( ...
                  corr(cat_stat_nanmean(kappa{di}(:,1:3),2), genSIQR{di,qi}(:,rmsi), 'Type',corrtype), 1, size(kappa{di},2));
              end
            else
              for ti = 1:size(kappa{1},2)
                rkappa{mi,qai,di}(qi,rmsi,ti)   = corr(kappa{di}(:,ti), SIQR{di,qi}(:,rmsi), 'Type',corrtype);
                if rmsi <= size(genSIQR{di,qi},2)
                  gkappa{mi,qai,di}(qi,rmsi,ti) = corr(kappa{di}(:,ti), genSIQR{di,qi}(:,rmsi),'Type',corrtype);
                end
              end
            end
          end
          
          if isfield( X{di,mi,qai}(si) , 'subjectmeasures' )
            for si = 1:numel(X{di}) % subject
              if di==1 
                % normalize by first scan of the BWP without distortions
                rBV{di}(si,1:3) = X{di,mi,qai}(si).subjectmeasures.vol_rel_CGW - ...
                  X{di,mi,qai}(1).subjectmeasures.vol_rel_CGW;
              else
                % normalize by standard scan (without motion artifacts)
                BLsi = contains( subID{di} , subID{di}{si}) & BL{di};
    
                rBV{di}(si,1:3) = X{di,mi,qai}(si).subjectmeasures.vol_rel_CGW - ... 
                                  X{di,mi,qai}(BLsi).subjectmeasures.vol_rel_CGW ; 
              end
            end
          else
            fprintf('ERR')
          end

          % average
          for ti = 1:size(rBV{di},2)
            for rmsi = 1:16
              rtissue{mi,qai,di}(qi,rmsi,ti)   = abs( corr(rBV{di}(:,ti), SIQR{di,qi}(:,rmsi),'Type',corrtype) );
              if rmsi <= size(genSIQR{di,qi},2)
                gtissue{mi,qai,di}(qi,rmsi,ti) = abs( corr(rBV{di}(:,ti), genSIQR{di,qi}(:,rmsi),'Type',corrtype) );
              end
            end
          end
        end
        ti = size(rBV{di},2)+1;
        rtissue{mi,qai,di}(:,:,ti) = cat_stat_nanmean(abs(rtissue{mi,qai,di}),3);
        gtissue{mi,qai,di}(:,:,ti) = cat_stat_nanmean(abs(gtissue{mi,qai,di}),3);
      end
    end
  end
end


%%
if numel(QCs{1})>1
  maindir   = fullfile(datadir,'ds004173-download'); 
  exprating = fullfile(maindir,'derivatives','scores.tsv');
  copyfile(exprating,strrep(exprating,'.tsv','.csv')); 
  tsv = cat_io_csv(exprating,'','',struct('delimiter','\t'));
  for ti = 2:size(tsv,1)  
    tid = find( contains( P.ARTfiles , tsv{ti,1}) == 1);
    group(tid,1) = tsv{ti,2};
  end
  
  di = 2;  %
  
  fx = figure(23); fx.Visible = 'on';  
  fx.Position(3:4) = [1700 500]; 
  gcolor = [0 0.8 0; 0.8 0.7 0; 0.8 0 0]; 
  for plotcase = 1:2
    if plotcase
      subpl = size(genSIQR{di,1},2); 
    else
      subpl = size(kappa{1},2);
    end
    tiledlayout(2, subpl, 'TileSpacing', 'compact', 'Padding', 'compact'); 
    for qi = [1,numel(QCs{1})]
      for ddi = 1:subpl
        nexttile; 
        tmp = genSIQR{di,qi}(:,ddi);
        if plotcase == 1
          cat_plot_boxplot( { tmp(group==1), tmp(group==2), tmp(group==3) }, ...
            struct('usescatter',1,'groupcolor',gcolor,'names', ...
            {{'no','light','strong'}}));
          ylim([1.5 6.5]); hold on
          plot([0.5 3.5],[4 4],'--','Color',[.8 0 0])
          plot([0.5 3.5],[3 3],'--','Color',[.5 .5 0])
          plot([0.5 3.5],[2 2],'--','Color',[0 .5 0])
        else
          for ti = 1:3
            sch = scatter( ...
              tmp(group==ti & kappa{di}(:,4)<1) , ...
              kappa{di}(group==ti & kappa{di}(:,4)<1,qi),'filled' ); hold on
            sch.MarkerEdgeAlpha = .5; 
            sch.MarkerFaceAlpha = .5; 
            sch.MarkerEdgeColor = gcolor(ti,:);
            sch.MarkerFaceColor = gcolor(ti,:);
          end   
          [curve1{ti}, goodness] = fit(  tmp(kappa{di}(:,qi)<1) , ...
            kappa{di}( kappa{di}(:,qi)<1 ,qi)  ,'poly2');%,'robust','LAR');
          ph = plot(curve1{ti});
          set(ph,'Color',gcolor(ti,:)); 
          legend(sprintf('fit(r^2=%0.6f)',goodness.rsquare))
        end
        title(genSIQRmod{ddi});
        if ddi==1, ylabel( name{end}(3:end),'Interpreter','none'); end
        grid on 
      end
    end
  
    % save result
    if plotcase == 1
      fx.Name = sprintf('maincorrelation%d_boxplotMRART_SIQR_%s_%s',qcsi,segment{mi},qaversions{qai});
    else
      fx.Name = sprintf('maincorrelation%d_kappaMRART_SIQR_%s_%s',qcsi,segment{mi},qaversions{qai});
    end
    if ~exist(opt.resdir,'dir'), mkdir(opt.resdir); end
    print(fx, '-djpeg', '-r300', fullfile(opt.resdir,fx.Name)); 
  end
end      



%% create a average variable over all segments and QC versions
qaisel = 1:numel(qaversions); 
clear tgtissue tgkappa
for di = 1:size(X,1)
  tgtissue{di} = gtissue{mi,qai,di}*0; 
  tgkappa{di}  = gkappa{mi,qai,di}*0;  
  for mi = 1:numel(segment)
    for qai = qaisel
      %qaversions{qai}
      tgtissue{di}  = tgtissue{di} + gtissue{mi,qai,di}(end,:,:) / (numel(segment) * numel(qaisel));
      tgkappa{di}   = tgkappa{di}  + gkappa{mi,qai,di}(end,:,:)  / (numel(segment) * numel(qaisel));
    end
  end
end

% main table 
% BWP and MRART  volume  values for each class
C1 = num2cell( abs([ shiftdim(tgtissue{1}(1,:,:),1) , shiftdim(tgtissue{2}(1,:,:),1)] )); 
cell2table( C1 , 'RowNames',genSIQRmod,'VariableNames',[strcat('BWP-',tissue),strcat('MRART-',tissue)] ); 
% BWP and MRART  kappa  values for each class
C2 = num2cell( abs([ shiftdim(tgkappa{1}(1,:,:),1) , shiftdim(tgkappa{2}(1,:,:),1)] )); 
cell2table( C2 , 'RowNames',genSIQRmod,'VariableNames',[strcat('BWP-',tissue),strcat('MRART-',tissue)]);
% BWP and MRART volume and kappa values (average tissues values) 
C3 = num2cell( abs([ ...
  shiftdim(tgtissue{1}(1,:,4),1) , shiftdim(tgtissue{2}(1,:,4),1) ...
  shiftdim(tgkappa{1}(1,:,4),1)  , shiftdim(tgkappa{2}(1,:,4),1) ...
  mean( abs( [shiftdim(tgtissue{1}(1,:,4),1) , shiftdim(tgtissue{2}(1,:,4),1), ...
    shiftdim(tgkappa{1}(1,:,4),1) , shiftdim(tgkappa{2}(1,:,4),1) ] ),2) ] ));
C4 = num2cell( abs([ ...
  gtissue{mi,qai,1}(end,:,4)', gtissue{mi,qai,2}(end,:,4)', ...
  gkappa{mi,qai,1}(end,:,4)', gkappa{mi,qai,2}(end,:,4)', ...
  mean( abs( [gtissue{mi,qai,1}(end,:,4)', gtissue{mi,qai,2}(end,:,4)', ...
  gkappa{mi,qai,1}(end,:,4)', gkappa{mi,qai,2}(end,:,4)' ] ) , 2)]));
CH = [strcat('Vol-BWP-',tissue(4)), strcat('Vol-MRART-',tissue(4)), ...
      strcat('Kappa-BWP-',tissue(4)), strcat('Kappa-MRART-',tissue(4)) ,'avg']; 
T = cell2table( C3 ,'RowNames',genSIQRmod,'VariableNames', CH);

% save result
if ~exist(opt.resdir,'dir'), mkdir(opt.resdir); end
fname = fullfile(opt.resdir,...
  sprintf('table_correlation_avgQRs_%0.0fsegmentations_%0.0fQAversion', ...
    numel(segment), numel(qaversions)));
cat_io_csv(fname, [{corrtype} CH; genSIQRmod' C3]);
fname = fullfile(opt.resdir,...
  sprintf('table_correlation_SIQR_%0.0fsegmentations_%0.0fQAversion', ...
    numel(segment), numel(qaversions)));
T = cell2table( C4 ,'RowNames',genSIQRmod,'VariableNames', CH)
cat_io_csv(fname, [{corrtype} CH; genSIQRmod' C4]);
%%
for mi = 1:numel(segment)
  for qai = 1:numel(qaversions)
    for qcsi = 1:numel(QCs) 
      CM = cat_io_colormaps('set1',numel(QCs{qcsi})); CM(end,:) = [0 0 0];
      
      
      %% table 
      C1 = num2cell( abs([ shiftdim(gtissue{mi,1}(1,:,:),1) , shiftdim(gtissue{mi,2}(1,:,:),1)] )); 
      cell2table( C1 , 'RowNames',genSIQRmod,'VariableNames',[strcat('BWP-',tissue),strcat('MRART-',tissue)] ); 

      C2 = num2cell( abs([ shiftdim(gkappa{mi,1}(1,:,:),1) , shiftdim(gkappa{mi,2}(1,:,:),1)] )); 
      cell2table( C2 , 'RowNames',genSIQRmod,'VariableNames',[strcat('BWP-',tissue),strcat('MRART-',tissue)])

      C3 = num2cell( abs([ ...
        shiftdim(gtissue{mi,1}(1,:,4),1) , shiftdim(gtissue{mi,2}(1,:,4),1) ...
        shiftdim(gkappa{mi,1}(1,:,4),1)  , shiftdim(gkappa{mi,2}(1,:,4),1) ...
        mean( abs( [shiftdim(gtissue{mi,1}(1,:,4),1) , shiftdim(gtissue{mi,2}(1,:,4),1), ...
          shiftdim(gkappa{mi,1}(1,:,4),1) , shiftdim(gkappa{mi,2}(1,:,4),1) ] ),2) ] )); 
      cell2table( C3 , 'RowNames',genSIQRmod,'VariableNames', ...
        [strcat('Vol-BWP-',tissue(4)), strcat('Vol-MRART-',tissue(4)), strcat('Kappa-BWP-',tissue(4)), strcat('Kappa-MRART-',tissue(4)) ,'avg'])
      
       
       %% correlation between tissue volume reduction and SIQR measures
       
       if 0
        for di = 1:4
          fx = figure(di+1); clf(di+1); fx.Position(3:4) = [1200 400];
  
          switch di 
            case 1, fx.Name = sprintf('correlation%d-Volume-SIQR_BWP_%s_%s'  ,qcsi,segment{mi},qaversions{qai});
            case 2, fx.Name = sprintf('correlation%d-Volume-SIQR_MRART_%s_%s',qcsi,segment{mi},qaversions{qai});
            case 3, fx.Name = sprintf('correlation%d-Kappa-SIQR_BWP_%s_%s'   ,qcsi,segment{mi},qaversions{qai});
            case 4, fx.Name = sprintf('correlation%d-Kappa-SIQR_MRART_%s_%s' ,qcsi,segment{mi},qaversions{qai});
          end
          
          % create subfigures for the tissues as well as the average
          tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact'); 
          for ti = 1:4 % 3 tissue + avg
            nexttile;
      
            % select figure
            if di>=3
              p = plot(abs(rkappa{mi,di-2}(:,:,ti)')); 
              title(sprintf('corr(%sKappa,SIQR)',tissue{ti}),'Interpreter','none')
              ylabel([tissue{ti} ' Kappa']); 
            else
              p = plot(abs(rtissue{mi,qai,di}(:,:,ti)')); 
              title(sprintf('corr(%sV,SIQR)', tissue{ti}),'Interpreter','none')
              ylabel([corrtype ' correlation']); 
            end
      
            % format lines
            for pi=1:numel(p), p(pi).Color = CM(pi,:); p(pi).LineWidth = 1.5; p(pi).Marker = MK( 1+mod(pi,numel(MK))); end  
      
            % label
            subtitle(sprintf('%s - %s',qaversions{qai}, segment{mi}),'Interpreter','none');
            xlabel('power used for SIQR')
            if ti == 4, legend(name,'Interpreter','none','Location','SouthEast'); end
            ax=gca; ax.XTick = [1 2 4 8 16]; grid on
            switch qcsi
              case 1
                switch di 
                  %case {1}, ylim([.4 .8]);
                  %case {2}, ylim([.4 .8]);
                  %case {3}, ylim([.5 1]);
                  case {4}, ylim([.82 .9]);
                end
              otherwise
                switch di 
                  case {1},   ylim([.4 .9]);
                  %case {2},   ylim([.4 .95]);
                  case {3},   ylim([.6  1]);
                  case {4},   ylim([.75 .92]);
                end
            end
          end
      
          % save result
          if ~exist(opt.resdir,'dir'), mkdir(opt.resdir); end
          print(fx, '-djpeg', '-r300', fullfile(opt.resdir,fx.Name)); 
        end
      end

      
      %% total figure only with averages
      tisnam = {'CSF','GM','WM','avg'}; 
      for ti = 1:4
        fy = figure(6); clf(6); fy.Position(3:4) = [1200 350];
        fy.Name = sprintf('maincorrelation%d-SIQR_%s_%s_%s',...
          qcsi,segment{mi},tisnam{ti},qaversions{qai});
        tiledlayout(1, 5, 'TileSpacing', 'compact', 'Padding', 'compact'); 
        dataset = {'BWP','MRART','BWP','MRART'};
        for di = 1:5
          nexttile; hold on
         
          % select figure
          if di >= 5
            p = plot( (abs(rkappa{mi,1}(:,:,ti)') + abs(rkappa{mi,2}(:,:,ti)') + ...
                       abs(rtissue{mi,qai,1}(:,:,ti)') + abs(rtissue{mi,qai,2}(:,:,ti)') ) / 4); 
            title(sprintf('corr(%sKappa+dV%s,SIQR) BWP+MRART',tissue{ti},tissue{ti}),'Interpreter','none')
            ylabel([tissue{ti} ' Kappa']); 
          elseif di >= 3
            p = plot(abs(rkappa{mi,di-2}(:,:,ti)')); 
            title(sprintf('corr(%sKappa,SIQR) %s',tissue{ti}, dataset{di}),'Interpreter','none')
            ylabel([tissue{ti} ' Kappa']); 
          else
            p = plot(abs(rtissue{mi,qai,di}(:,:,ti)')); 
            title(sprintf('corr(dV%s,SIQR) %s', tissue{ti}, dataset{di}),'Interpreter','none')
            ylabel([corrtype ' correlation']); 
          end
  
          % format lines
          for pi=1:numel(p), p(pi).Color = CM(1+mod(pi-1,size(CM,1)),:); p(pi).LineWidth = 1.5; p(pi).Marker = MK( 1+mod(pi,numel(MK) )); end  
        
          % label
          subtitle(sprintf('%s - %s',qaversions{qai}, segment{mi}),'Interpreter','none');
          xlabel('power used for SIQR')
          if di==5, legend(name,'Interpreter','none','Location','SouthWest'); end
          ax=gca; ax.XTick = [1 2 4 8 16]; grid on; box on; 
          if 1
            ylim(ax,[0 1])
          else
            switch qcsi 
              case 1
                switch di 
                  case {1 2}, ylim(ax,[.40 1]);
                  case {3 4}, ylim(ax,[.40 1]);
                end
              otherwise
                switch di 
                  case {1 2}, ylim(ax,[.45 .80]);
                  case {3 4}, ylim(ax,[.75 1]);
                end
            end
          end
        end
        % save result
        if ~exist(opt.resdir,'dir'), mkdir(opt.resdir); end
        print(fy, '-djpeg', '-r300', fullfile(opt.resdir,fy.Name)); 
      end
    end
  end
end

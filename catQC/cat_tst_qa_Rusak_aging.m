function cat_tst_qa_Rusak_aging( datadir0, qaversions, segment, fasttest, rerun )  
%% Rusak_aging 
%  ------------------------------------------------------------------------
%
%  Requirements: 
%   0. Download and install SPM and CAT
%   1. Download Rusak T1 data from: 
%        https://doi.org/10.25919/4ycc-fc11
%
%   2. Specify in this script: 
%      1) the data directory "datadir" 
%      2) the QC version you would like to tests (the file has to exist in the cat directory) 
%      3) the segmentation you would like to use
%
%  ------------------------------------------------------------------------

%#ok<*AGROW>

cat_io_cprintf([0 0.5 0],'\n\n== Run cat_tst_qa_Rusak_aging ==\n') 

% directories
runPP      = 1; 

% ### datadir ###
if ~exist( 'datadir0' , 'var' )
  datadir    = '/Volumes/SG5TB/MRData/202503_QA/Rusak2021';
  datadir   = '/Volumes/WDE18TB/MRData/Dahnke2025_QC/Rusak2021'; %20211122-SyntheticDataset 
else
  datadir    = fullfile(datadir0,'Rusak2021');
end
% ### segmention ###
if ~exist( 'segment' , 'var')
  segment = {'SPM'}; % {'SPM','CAT','qcseg'}; % qcseg requires cat_vol_qa2024012
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
fast = {'full','fast'};

resdir     = fullfile(fileparts(datadir), '+results',['Rusak2021_' fast{fasttest+1} '_' datestr(clock,'YYYYmm')]); 
if ~exist(resdir,'dir'), mkdir(resdir); end


if runPP
  for si = 1:numel(segment)
    clear matlabbatch; 
    Rusakfiles = cat_vol_findfiles( datadir , 'sub-ADNI*.nii', struct('maxdepth',3));
    %Rusakfiles( contains( Rusakfiles , 'err' ) ) = [];
    if isempty(Rusakfiles), fprintf('Cannot find Rusak''s files. Check the directory but also the depth parameter parameter in the line ahead this message.'); end
    switch segment{si}
      case 'CAT'
        CATpreprocessing4qc;
        matlabbatch{1}.spm.tools.cat.estwrite.data = Rusakfiles;
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.lazy = 1; % derivatives
        spm_jobman('run',matlabbatch);
      case 'SPM'
        SPMpreprocessing4qc;
        RusakfilesSPM = Rusakfiles; 
        RusakfilesSPM( cellfun(@(x) exist(x,'file'),spm_file(RusakfilesSPM,'prefix','c1'))>0 ) = [];
        if ~isempty( RusakfilesSPM )
          matlabbatch{1}.spm.spatial.preproc.channel.vols = RusakfilesSPM;
          spm_jobman('run',matlabbatch);
        end
      case 'synthseg'
        error('synthseg is not prepared in the public script ... use MATLAB help')
      case 'qcseg'
        fprintf('No preprocessing required.\n\n');
    end
  end
end

%%
qias = 1:numel(qaversions);
for segi = 1:numel(segment)
  switch segment{segi}
    case 'CAT'
      Pp0{segi} = cat_vol_findfiles( datadir ,'p0*',struct('depth',0));
      %Pp0{segi} = cat_vol_findfiles(fullfile(datadir,'derivatives','CAT12.9_2583'),'p0*.nii'); 
    case 'SPM'
      Pp0{segi} = cat_vol_findfiles( datadir , 'c1*',struct('depth',0));
    case 'synthseg'
      Pp0{segi} = cat_vol_findfiles( datadir , 'synthseg_p0*',struct('depth',0));
    case 'qcseg'
      Pp0{segi} = cat_vol_findfiles( datadir , 'sub*.nii',struct('depth',2));
  end
  if fasttest, Pp0{si} = Pp0{si}(1:20*3); Pp0{si} = Pp0{si}(1:2:end); end
  fprintf('Process Rusak2021: \n')
  for qai = qias
    switch segment{segi}
      case 'CAT'
        cat_vol_qa('p0',Pp0{segi},struct('prefix','VERSION_','version',qaversions{ qai },'rerun',rerun));
      case 'SPM'
        cat_vol_qa('p0',Pp0{segi},struct('prefix','VERSION_spm_','version',qaversions{ qai },'rerun',rerun));
      case 'synthseg'
        cat_vol_qa('p0',Pp0{segi},struct('prefix','VERSION_synthseg_','version',qaversions{ qai },'rerun',rerun));
      case 'qcseg'
        cat_vol_qa('p0',Pp0{segi},struct('prefix','VERSION_qcseg_','version',qaversions{ qai },'rerun',rerun));
    end
  end
  fprintf('Process Rusak2021 done. \n')
  
  
  
  
  %% evaluate test
  sub = cat_vol_findfiles(datadir,'sub*',struct('dirs',true,'depth',2));
  if fasttest, sub = sub(1:3); end 
  for qai = qias
    % get and read XMLs, set subjects and time points
    Pxml = {}; SID = []; TID = []; TP = []; 
    for si = 1:numel(sub)
      switch segment{segi}
        case 'CAT'
          Pxmls = sort(cat_vol_findfiles( fullfile( sub{si} ,'report') , sprintf('%s_sub*.xml',strrep(qaversions{qai},'_','')) )); 
        case 'SPM'
          Pxmls = sort(cat_vol_findfiles( sub{si} , sprintf('%s_spm_sub*.xml',strrep(qaversions{qai},'_','')) )); 
        case 'qcseg'
          Pxmls = sort(cat_vol_findfiles( fullfile( sub{si} ,'report') , sprintf('%s_qcseg_sub*.xml',strrep(qaversions{qai},'_','')) )); 
      end
      %fprintf('%d - %d - %s\n',si,numel(Pxmls),sub{si}); 
      if fasttest, Pxmls = Pxmls(1:2:end); end
      Pxml  = [Pxml; Pxmls]; 
      TPs   = char(Pxmls); TPs = str2num(TPs(:,end-9:end-6)); %#ok<ST2NM>
      TP    = [TP; TPs];
      SID(end+1:end+numel(Pxmls)) = si; 
      TID(end+1:end+numel(Pxmls)) = 1:numel(Pxmls); 
    end
    XML = cat_io_xml(Pxml); 
    
    %%
    clear Q N; 
    ptn = {'full','short'}; 
    for pt = 2 %:2
      % set QMs 
      if pt == 1
        mylabel = 'grad error'; 
        QM = {'NCR','FEC','res_RMS','res_ECR','ICR','contrastr','IQR','SIQR'}; 
        cl = max(0,flip(hsv(8)) - 0.2); cl(end,:) = [ 0 0 0]; 
        mylim = [-1 1]; loc = 'Southwest'; 
      else
        % Paper version with less values
        mylabel = 'error (rps)'; 
        QM = {'NCR','ICR', 'res_RMS','res_ECR','FEC','SIQR'};  % 
        cl = [.8 0 .2;   0.8 0.6 0;  0.2 0.5 0;   .0 .5 .9; 0.1 0 0.9;  0 0 0 ]; %
        mylim = [-1 1]*3; loc = 'Southwest'; 
      end
      QMs = cat_io_strrep(QM,{'res_RMS','res_ECR','contrastr'},{'RMS','ECR','CON'});
      cm = '<>^vosdph*+_.<>^vosdph*+_.';
      for qmi = 1:numel(QM)
        for xi = 1:numel(Pxml) 
          if isfield(XML(xi).qualityratings,QM{qmi})
            Q.(QM{qmi})(SID(xi),TID(xi)) = XML(xi).qualityratings.(QM{qmi}); 
            N.(QM{qmi})(SID(xi),TID(xi)) = Q.(QM{qmi})(SID(xi),TID(xi)) - Q.(QM{qmi})(SID(xi),1); 
          else
            Q.(QM{qmi})(SID(xi),TID(xi)) = XML(xi).subjectmeasures.vol_rel_CGW(2) * 10;
            N.(QM{qmi})(SID(xi),TID(xi)) = Q.(QM{qmi})(SID(xi),TID(xi))  - Q.(QM{qmi})(SID(xi),1) ;
          end
        end
        try
          [N.r.(QM{qmi}),    N.p.(QM{qmi})]              = corr(TP,N.(QM{qmi})(:));
        catch
          N.r.(QM{qmi}) = nan; 
          N.p.(QM{qmi}) = nan; 
        end
        try
          [N.fit.(QM{qmi}),  N.fit.([QM{qmi} '_stat'])]  = robustfit(TP,N.(QM{qmi})(:));
          [N.fitm.(QM{qmi}), N.fitm.([QM{qmi} '_stat'])] = robustfit(TPs,mean(N.(QM{qmi})));
        catch
          N.fit.(QM{qmi})  = nan(1,2); 
          N.fitm.(QM{qmi}) = nan(1,2); 
          N.fit.([QM{qmi} '_stat'])  = nan(1,2); 
          N.fitm.([QM{qmi} '_stat']) = nan(1,2);  
        end
      end
      
      % figure
      fh = figure(4); clf; 
      fh.Position(3:4) = [1000 400];
      fname = sprintf('Rusak2021_%s_%s',qaversions{qai},datestr(clock,'YYYYmmdd'));
      annotation('textbox',[0 0.95 1 0.05 ],'FitBoxToText','off','LineStyle','none','Fontsize',12,'String',...
        sprintf('Rusak2021: %s - %s', strrep(qaversions{qai},'_','\_'), datestr(clock,'YYYYmmdd')) );
      pos = [1 1 2 2 3 3 4 4 5 5 6 6; 
             2 1 2 1 2 1 2 1 2 1 2 1] - 1; 
      mpos = @(pos,qmi) [0.03+(0.85/max(1+pos(1,:)))*pos(1,qmi)   0.09+(.95/max(1+pos(2,:)))*pos(2,qmi)  .68/max(pos(1,:))  0.33]; 
      
    
      % print each QR
      for qmi = 1:numel(QM)
        subplot('Position',mpos(pos,qmi)); hold on;  
        for si=1:max(SID)
          try
            pl = scatter(TP(SID==si), N.(QM{qmi})(si,:)); pl.SizeData = 10; pl.Marker = cm(qmi); 
            pl.MarkerFaceColor = cl(qmi,:); pl.MarkerEdgeColor =  pl.MarkerFaceColor; pl.MarkerFaceAlpha = .3; pl.MarkerEdgeAlpha = .4; 
            pl = plot(TP(SID==si), N.(QM{qmi})(si,:)); pl.Color = [cl(qmi,:) .05]; pl.LineWidth = .5; 
          end
        end
        try
          pl = plot(TP(SID==1), mean(N.(QM{qmi}),1)); pl.Color = max(0,cl(qmi,:)-.2); pl.LineWidth = 1; 
        end
        a=gca; a.XTick = [0 0.1 0.5 1];
        if pos(1,qmi)==0, ylabel(mylabel); else, a.YTickLabel = {};  end
        xlabel('atrophy (in mm)');  
        QMsx{qmi} = sprintf('%s (%+0.2f,r=%0.2f,p=%0.0e)',QMs{qmi}, ...
          N.fitm.(QM{qmi})(2), N.r.(QM{qmi}), N.p.(QM{qmi})  );
        ylim(mylim); title(QMsx{qmi});  box on; grid on;  
        if strcmp(QM{qmi},'GMV'), ylim([-0.15 0.15]); end
      end
    
    
      % overview
      subplot('Position',(mpos(pos,9) + [.05 0 0 0]) .* [1 1 1 1]); hold on, 
      for qmi = 1:numel(QM)
        try
          pl = plot(TP(SID==1), mean(N.(QM{qmi}),1)); 
          pl.Color = [cl(qmi,:) 0.4 + 0.6*(qmi==numel(QM))]; pl.LineWidth = .5 + .5*(qmi > (numel(QM)-2)); 
          pl.Marker = cm(qmi); pl.MarkerSize = 4; pl.MarkerFaceColor = pl.Color;  
        end
      end
      for qmi = 1:numel(QM)
        try
          pl = scatter(TP(SID==1), mean(N.(QM{qmi}),1)); pl.SizeData = 10; pl.Marker = cm(qmi); 
          pl.MarkerFaceColor = cl(qmi,:); pl.MarkerEdgeColor =  pl.MarkerFaceColor; 
          pl.MarkerFaceAlpha = .2 + 0.8*(qmi > (size(QM,1)-2)); pl.MarkerEdgeAlpha = pl.MarkerFaceAlpha;
        end
      end
      a=gca; a.XTick = [0 0.1 0.5 1];
      legend(QMsx,'location',loc,'FontSize',6,'box','off'); 
      ylabel(mylabel); xlabel('atrophy (in mm)'); title(sprintf('overview QRs')); 
      ylim(mylim); box on; grid on; a=gca; %a.YTickLabel = {};
       
    
      % boxplot
      subplot('Position',(mpos(pos,10) + [.05 0 0 0]) .* [1 1 1 1]); hold on, clear boxval; 
      for qmi = 1:numel(QM), boxval{qmi} = single(N.(QM{qmi})(:)); end
      cat_plot_boxplot( boxval , struct('names',{QMs},'ylim',mylim,'style',0,'usescatter',1,'ygrid',1,'groupcolor',cl));
      ylabel(mylabel); title(sprintf('boxplot QRs')); a=gca; a.XTickLabelRotation = 90; 
    
      % change
      subplot('Position',(mpos(pos,11) + [.1 0 0 0]) .* [1 1 1 1]); hold on, clear boxval; 
      for qmi = 1:numel(QM), mcorr(qmi) = N.fitm.(QM{qmi})(2); end
      mcorr(qmi+1) = mean(mcorr); 
      bh = bar(mcorr(1:end-1)); bh.CData = cl; bh.FaceColor = 'flat';
      ylim(mylim); xlim([.4 qmi+0.6]); xticks(1:numel(QMs)); xticklabels([QMs {'avg'}]);
      a=gca; a.XTickLabelRotation = 90; a.YGrid = 'on'; 
      ylabel(mylabel); box on;
      title(sprintf('Change (avg=%0.3f)',mcorr(end)),'FontWeight','bold','Fontsize',10 - 2*(pt==1)); 
      for fi = 1:numel(mcorr)-1
        dt(fi) = text(fi-.45, mcorr(fi) + (1-2*(mcorr(fi)<0)) * .065*mylim(2), sprintf('%0.3f',mcorr(fi)),'FontSize',8,'Color',cl(fi,:)); %#ok<*SAGROW>
      end  
      
      % RMSE
      rmse = @(a)   max(0,cat_stat_nanmean(a.^2).^(1/2)); clear rmseval
      subplot('Position',(mpos(pos,12) + [.1 0 0 0]) .* [1 1 1 1]); hold on, clear boxval; 
      for qmi = 1:numel(QM), rmseval(qmi) = rmse(N.(QM{qmi})(:)); end
      rmseval(qmi+1) = mean(rmseval); 
      bh = bar(rmseval(1:end-1)); bh.CData = cl; bh.FaceColor = 'flat';
      ylim([0,mylim(2)]); xlim([.4 qmi+0.6]); xticks(1:numel(QMs)); xticklabels([QMs {'avg'}]);
      a=gca; a.XTickLabelRotation = 90; a.YGrid = 'on'; 
      ylabel(mylabel); box on;
      title(sprintf('RMSE (avg=%0.3f)',rmseval(end)),'FontWeight','bold','Fontsize',10 - 2*(pt==1)); 
      for fi = 1:numel(rmseval)-1
        dt(fi) = text(fi-.45, rmseval(fi) + .1/3*mylim(2), sprintf('%0.3f',rmseval(fi)),'FontSize',8,'Color',cl(fi,:)); %#ok<*SAGROW>
      end  
     
      % print
      print(fh, '-dpng', '-r600', fullfile(resdir,[ptn{pt} '_' segment{segi} '_' fname])); 
    end
  end
end
fprintf('Evaluate Rusak2021 done.\n'); 


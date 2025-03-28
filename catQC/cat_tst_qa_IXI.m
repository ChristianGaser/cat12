function cat_tst_qa_IXI( datadir0 , qaversions , segment, fasttest, rerun ) 
%% Evaluation of CATQC in IXI 
%  ------------------------------------------------------------------------
%
%  Requirements: 
%   0. Download and install SPM and CAT
%   1. Download IXI T1 data from: 
%      https://brain-development.org/ixi-dataset/
%      http://biomedic.doc.ic.ac.uk/brain-development/downloads/IXI/IXI-T1.tar
%
%   2. Specify in this script: 
%      1) the data directory "datadir" 
%      2) the QC version you would like to tests (the file has to exist in the cat directory) 
%      3) the segmentation you would like to use
%
%  ------------------------------------------------------------------------

  cat_io_cprintf([0 0.5 0],'\n\n== Run cat_tst_qa_IXI ==\n') 
  
  % ### datadir ###
  if ~exist( 'datadir0' , 'var' )
    datadir   = '/Volumes/SG5TB/MRData/202503_QA/IXI'; 
  else
    datadir   = fullfile(datadir0,'IXI-T1'); 
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

  % ### segmention ###
  if ~exist( 'segment' , 'var')
    segment = {'CAT'}; % {'SPM','CAT','qcseg'}; % qcseg requires cat_vol_qa2024012
  end

  if ~exist( 'fasttest', 'var'), fasttest = 0; end
  if ~exist( 'rerun', 'var'), rerun = 0; end
  fast = {'full','fast'}; 

  runPP      = 1; % run CAT/SPM preprocessing 
  useratings = 1; 
  printall   = 1; 

  resdir     = fullfile(fileparts(datadir), '+results',['IXI_' fast{fasttest+1} '_' datestr(clock,'YYYYmm')]); 
  if ~exist(resdir,'dir'), mkdir(resdir); end
  
  qias = 1:numel(qaversions);
  
  if runPP
    for si = 1:numel(segment)
      clear matlabbatch; 
      IXIfiles = cat_vol_findfiles( datadir , 'IXI*.nii',struct('depth',0));
      switch segment{si}
        case 'CAT'
          CATpreprocessing4qc;
          matlabbatch{1}.spm.tools.cat.estwrite.data = IXIfiles;
          matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.lazy = 1; 
          spm_jobman('run',matlabbatch);
        case 'SPM'
          SPMpreprocessing4qc;
          IXIfilesSPM = IXIfiles; 
          IXIfilesSPM( cellfun(@(x) exist(x,'file'),spm_file(IXIfilesSPM,'prefix','c1'))>0 ) = [];
          if ~isempty( IXIfilesSPM )
            matlabbatch{1}.spm.spatial.preproc.channel.vols = IXIfilesSPM;
            spm_jobman('run',matlabbatch);
          end
        case 'synthseg'
          error('synthseg is not prepared in the public script ... use MATLAB help')
        case 'qcseg'
          fprintf('No preprocessing required.\n\n');
      end
    end
  end
  
  
  
  
  %% (re)process QC values
  fprintf('Prepare IXI: \n')
  for si = 1:numel(segment)
    for qai = qias
      switch segment{si}
        case 'CAT'
          Pp0{si}{qai} = cat_vol_findfiles( fullfile( datadir , 'mri') , 'p0*',struct('depth',0));
        case 'SPM'
          Pp0{si}{qai} = cat_vol_findfiles( datadir , 'c1*',struct('depth',0));
        case 'synthseg'
          Pp0{si}{qai} = cat_vol_findfiles( datadir , 'synthseg_p0*',struct('depth',0));
        case 'qcseg'
          Pp0{si}{qai} = cat_vol_findfiles( datadir , 'IXI*.nii',struct('depth',0));
      end
      fasttestss = 10; 
      if fasttest, Pp0{si}{qai} = Pp0{si}{qai}(1:fasttestss:end); end
    end
  end
  fprintf('Prepare IXI done. \n')
  
  fprintf('Process IXI: \n')
  for si = 1:numel(segment)
    for qai = qias
      switch segment{si}
        case 'CAT'
          cat_vol_qa('p0',Pp0{si}{qai},struct('prefix',[qaversions{qai} '_'],'version',qaversions{ qai },'rerun',rerun));
        case 'SPM'
          cat_vol_qa('p0',Pp0{si}{qai},struct('prefix',[qaversions{qai} '_spm_'],'version',qaversions{ qai },'rerun',rerun));
        case 'synthseg'
          cat_vol_qa('p0',Pp0{si}{qai},struct('prefix',[qaversions{qai} '_synthseg_'],'version',qaversions{ qai },'rerun',rerun));
        case 'qcseg'
          cat_vol_qa('p0',Pp0{si}{qai},struct('prefix',[qaversions{qai} '_qcseg_'],'version',qaversions{ qai },'rerun',rerun));
      end  
    end
  end
  fprintf('Process IXI done. \n')
  
  fprintf('Print IXI: \n')
  Pcsv = fullfile( fileparts(datadir) , 'IXI_2019main.csv' ); 
  csv  = cat_io_csv(Pcsv,'','',struct('delimiter',';','csv','.')); 
   
  
  
  % create figure
  fh = figure(39);
  %if ~fasttest
    fh.Visible = 'off';
  %end
  fh.Interruptible = 'off';
  fh.Position(3:4) = [600,200]; 
  for si = 1:numel(segment)
    for qai = qias
      fprintf('  Print %s with %s: \n',qaversions{qai},segment{si})
    
      % find and load QC XML files 
      if exist('Pp0','var') % old
        clear Pxml;
        for pi = 1:numel(Pp0{si}{qai})
          [pp,ff] = fileparts(Pp0{si}{qai}{pi});
          switch segment{si}
            case 'CAT'
              Pxml{pi,1} = fullfile(fileparts(pp),'report',[qaversions{qai} '_' ff(3:end) '.xml']);
            case 'SPM'
              Pxml{pi,1} = fullfile(pp,'report',[qaversions{qai} '_spm_' ff(3:end) '.xml']);
            case 'synthseg'
              Pxml{pi,1} = fullfile(pp,'report',strrep([qaversions{qai} ff '.xml'],'synthseg_p0','_synthseg_'));
            case 'qcseg'
              Pxml{pi,1} = fullfile(pp,'report',[qaversions{qai} '_qcseg_' ff '.xml']);
          end
        end
      else
        Pxml = cat_vol_findfiles( fullfile( datadir , 'report') , sprintf('%s_IXI*.xml',qaversions{qai}) ); 
        if fasttest && numel(Pxml)>numel(Pp0{si}{qai}), Pxml = Pxml(1:fasttestss:end); end
      end
      xml  = cat_io_xml(Pxml);
      
      % read out CAT XML QC value and get age values
      clear siten; 
      if printall
        measures = {'NCR','ICR','ECR','SIQR','FEC','CON','IQR','TIV','rGMV','rWMV','rCSFV'}; 
        scaling  = {[.5 6.5],[.5 6.5],[.5 6.5],[.5 6.5],[.5 6.5],[.5 6.5],[.5 6.5],[900 2100],[.2 .6],[.2 .6],[.2 .6]}; 
      else
        measures = {'SIQR','rGMV'};
        scaling  = {[.5 6.5],[.2 .6]};
      end
    
      clear NCR IQR ICR SIQR ECR ECRmm FEC CON TIV rGMV ID csvID age sex; ni = -1; 
      for xi = 1:numel(xml)
        fname{xi} = xml(xi).filedata.file; %#ok<*SAGROW>
        CON(xi)   = xml(xi).qualityratings.contrastr;
        NCR(xi)   = xml(xi).qualityratings.NCR;
        ICR(xi)   = xml(xi).qualityratings.ICR;
        IQR(xi)   = xml(xi).qualityratings.IQR;
        try
          FEC(xi) = xml(xi).qualityratings.FEC;
        catch
          FEC(xi) = nan; 
        end
        SIQR(xi)  = xml(xi).qualityratings.SIQR;
        if isfield(xml(xi).qualityratings,'res_ECR')
          ECR(xi) = xml(xi).qualityratings.res_ECR; 
        end
        if isfield(xml(xi).qualityratings,'res_ECRmm')
          ECRmm(xi) = xml(xi).qualityratings.res_ECRmm; 
        end
        TIV(xi)   = xml(xi).subjectmeasures.vol_TIV;
        rGMV(xi)  = xml(xi).subjectmeasures.vol_rel_CGW(2);
       
        % get ID of the XML entry
        seps      = find( fname{xi} == '-' );
        ID(xi)    = str2double( fname{xi}(4:6));
        siten{xi} = fname{xi}(seps(1)+1:seps(2)-1);

        % get age of the subject with the ID 
        try
          csvID(xi) = find( cell2mat(csv(2:end,1)) == ID(xi) , 1 , 'first' ) + 1; 
          age(xi)   = csv{ csvID(xi) , end }; 
          sex(xi)   = csv{ csvID(xi) , 2 }==1; 
        catch
          csvID(xi) = ni; ni = ni - 1; 
          age(xi)   = 0;
          sex(xi)   = 0; 
          site(xi)  = 0; 
        end
 
      end

      
      [siteu,~,site] = unique(siten);
      
      if useratings
        mark2rps = @(mark) min(100,max(0,105 - mark*10));
        ratings  = {'NCR','IQR','ICR','SIQR','ECR','FEC','CON'};
        for ri = 1:numel(ratings)
          if exist(ratings{ri},'var')
            eval(sprintf('%s = mark2rps(%s);',ratings{ri},ratings{ri})); 
          end
        end
        for ri = find(contains(measures,ratings))
          scaling{ri} = [40 100]; 
        end
        lgpos = 'South';
      else
        lgpos = 'North'; 
      end
  
      avar  = {'age','TIV'}; % analysis variables (x-axis)
      avart = {'in years; '; 'in ml; '};
      avars = {[20 90];[900 1900]};
      acols = [[0.0 0.4 0];   [0.2 0.6 0.2]; [0.4 .8 0.4];]; % aging 
      gcols = {[0.5 0 0];     [0 0 0.5 ];    [.5 .5 0]}; % group
      scols = [[0.0 0.4 0.8]; [1 0.5 0.0];   [0.8 .0 0.1];];
  
      %%
      tables.agecorr   = {'measure','r1','r2','r3','ra','mean(ra)', 'p1','p2','p3','pa','mean(pa)', 'a1','a2','a3','aa','mean(aa)' }; 
      tables.tivcorr   = {'measure','r1','r2','r3','ra','mean(ra)', 'p1','p2','p3','pa','mean(pa)', 'a1','a2','a3','aa','mean(aa)' }; 
      tables.ageanova  = {'measure','F','Prob>F'}; 
      tables.siteanova = {'measure','F','Prob>F'}; 
      tables.tivanova  = {'measure','F','Prob>F'}; 
      tables.sexanova  = {'measure','F','Prob>F'}; 
      tables.sexmwwt   = {'measure','p','Z','U'}; 
      tables.sidestd   = {'measure','std(site1)','std(site2)','std(site3)'};
      for mi = 1:numel(measures)
        if exist(measures{mi},'var') && eval(sprintf('any(~isnan(%s))',measures{mi})) 
          for ai = 1:numel(avar)
    
            if ~strcmp( avar{ai} , measures{mi} )
              %% first subfigure with scatter plot of all values colored by density 
              clf
              subplot('position',[0.08 0.16 .40 .75]);
              marker = {'o'  's'  'd'  'v'  '^'  '>'  '<'  'p'  'h'  }; 
              hold on;
              eval(sprintf('scatter( %s(site==1)'' , %s(site==1)'' ,20, ''o'')', avar{ai},measures{mi}));
              eval(sprintf('scatter( %s(site==2)'' , %s(site==2)'' ,20, ''^'')', avar{ai},measures{mi}));
              eval(sprintf('scatter( %s(site==3)'' , %s(site==3)'' ,20, ''v'')', avar{ai},measures{mi}));
              ax  = gca; 
              ax.Position = [0.08 0.16 .40 .75];
              ax.FontSize = 10; 
              hsd = findobj( get(ax,'Children') ,'Type','Scatter'); 
              for hsdi = 1:numel(hsd)
                hsd(hsdi).MarkerEdgeColor = scols(numel(hsd) + 1 - hsdi,:);
                hsd(hsdi).MarkerFaceColor = hsd(hsdi).MarkerEdgeColor; 
                hsd(hsdi).MarkerFaceAlpha = 0.2; 
                hsd(hsdi).MarkerEdgeAlpha = 0.4; 
              end
              box on; 
              title(sprintf('%s density in aging (IXI, N=%d)',measures{mi},numel(NCR))); 
              myylim = scaling{mi}; % round( ax.YLim .* [0.9 1.1] * 4 ) / 4 ; 
               xlim(avars{ai}); 
              % fitting biased by outiers ... this is pearson ...
              eval(sprintf('[R1,P1] = corrcoef( %s(site==1) , %s(site==1));',avar{ai},measures{mi}));  
              eval(sprintf('[R2,P2] = corrcoef( %s(site==2) , %s(site==2));',avar{ai},measures{mi}));  
              eval(sprintf('[R3,P3] = corrcoef( %s(site==3) , %s(site==3));',avar{ai},measures{mi}));  
              eval(sprintf('[R4,P4] = corrcoef( %s , %s);',avar{ai},measures{mi}));  
              P = cat_stat_nanmean([P1(2),P2(2),P3(2)]); R = cat_stat_nanmean([R1(2),R2(2),R3(2)]); %,P4(2), ,R4(2)
              if mi==1 && ai==1  % licence call
                try
                  eval(sprintf('[curve1, goodness, output] = fit( %s(site==1)'', double(%s(site==1))'',''poly1'');', ...
                    avar{ai}, measures{mi})); pf(1) = plot(curve1); pf(1).Color = gcols{1}; 
                end
                pause(2); 
              end
              % correlation 
              eval(sprintf('[curve1, goodness, output] = fit( %s(site==1)'', double(%s(site==1))'',''poly1'');', ...
                avar{ai}, measures{mi})); pf(1) = plot(curve1); pf(1).Color = scols(1,:); 
              eval(sprintf('[curve2, goodness, output] = fit( %s(site==2)'', double(%s(site==2))'',''poly1'');', ...
                avar{ai}, measures{mi})); pf(2) = plot(curve2); pf(2).Color = scols(2,:);
              eval(sprintf('[curve3, goodness, output] = fit( %s(site==3)'', double(%s(site==3))'',''poly1'');', ...
                avar{ai}, measures{mi})); pf(3) = plot(curve3); pf(3).Color = scols(3,:);
              eval(sprintf('[curve4, goodness, output] = fit( %s'', double(%s)'',''poly1'');', ...
                avar{ai}, measures{mi})); pf(4) = plot(curve4); pf(4).Color = [ 0.2 0.2 0.2 ]; pf(4).LineStyle = '--'; 
              pan = cat_stat_nanmean([curve1.p1, curve2.p1, curve3.p1]); 
              legend off; 
              
              % anova
              tab = cell2table([num2cell(age'),num2cell(site),num2cell(TIV')],'Variablenames',{'age','site','TIV'},'rownames',cellstr(num2str(csvID','IXI%03d')));
              eval(sprintf('[an0] = anova(tab,%s'');',measures{mi})); 
              cat_io_csv(fullfile(resdir,'anovas',sprintf('IXI_anova_%s-age-site_%s_%s.csv', measures{mi}, qaversions{qai}, segment{si})),...
                [{''} an0.stats.Properties.VariableNames; an0.stats.Properties.RowNames table2cell(an0.stats)]); 
  
              % man-wikney            
              tables.sidestd = [ tables.sidestd ;
                [ measures(mi) , num2cell( eval(sprintf('[std(%s(site==1)),std(%s(site==2)),std(%s(site==3))];',measures{mi},measures{mi},measures{mi})))] ]; 
      
              if ai == 1 
                tables.agecorr = [ tables.agecorr ; 
                  [ measures(mi) num2cell( [R1(2),R2(2),R3(2),R4(2),R,  P1(2),P2(2),P3(2),P4(2),P,  curve1.p1,curve2.p1,curve3.p1,curve4.p1,pan ]) ] ]; 
              else
                tables.tivcorr = [ tables.tivcorr ; 
                  [ measures(mi) num2cell( [R1(2),R2(2),R3(2),R4(2),R,  P1(2),P2(2),P3(2),P4(2),P,  curve1.p1,curve2.p1,curve3.p1,curve3.p1,pan ]) ] ]; 
              end
              ylim( myylim ); xlim(avars{ai}); 
              set(gca,'Ygrid',1)
              if useratings && any(contains(measures{mi},ratings))
                set(gca,'Ytick',scaling{mi}(1)+5 : 10 : scaling{mi}(2))
              end
  
              if strcmp(measures{mi},'rGMV'), set(gca,'YTickLabel',num2str(get(gca,'YTick')','%0.2f') ); end
              
              xlabel(sprintf('%s (%sb=%0.3f, \\it{r}\\rm{}=%0.3f, \\it{p}\\rm{}=%0.0e)',avar{ai},avart{ai},pan,R,P)); 
              ylabel(measures{mi}); grid on; hold on; 
              legend({
                sprintf('Guys (b=%0.3f, \\it{r}\\rm{}=%0.3f, \\it{p}\\rm{}=%0.3f)',curve1.p1, R1(2), P1(2) ), ...
                sprintf('HH (b=%0.3f, \\it{r}\\rm{}=%0.3f, \\it{p}\\rm{}=%0.3f)'  ,curve2.p1, R2(2), P2(2)) , ...
                sprintf('IOP (b=%0.3f, \\it{rp}\\rm{}=%0.3f, \\it{p}\\rm{}=%0.3f)' ,curve3.p1, R3(2), P3(2)) , ...
                },'Location',lgpos)
    
              % second figure with poxplot of age groups
              subplot('position',[0.575 0.16 .16 .75]);
              if ai == 1
                eval(sprintf('[an0,anova1age{mi}] = anova1(%s,1 +(age>40) + (age>60),''off'');',measures{mi})); 
                tables.ageanova = [ tables.siteanova; [ measures(mi) anova1age{mi}(2,4) anova1age{mi}(2,5) ] ]; 
                cat_plot_boxplot( eval(sprintf('{ %s(age<40); %s(age>40 & age<60); %s(age>60) }',...
                  measures{mi},measures{mi},measures{mi})) , ...
                  struct('names',{{'20-40','40-60','60-90'}},'style',4,'ylim',myylim,'ygrid',0, ...
                  'datasymbol','o','usescatter',1,'groupcolor',acols) );
                title('age groups');  
                xlabel('range (years)'); ylabel(measures{mi});
              else
                eval(sprintf('[an0,anova1tiv{mi}] = anova1(%s,1 +(TIV>1300) + (TIV>1500),''off'');',measures{mi})); 
                tables.tivanova = [ tables.siteanova; [ measures(mi) anova1tiv{mi}(2,4) anova1tiv{mi}(2,5) ] ]; 
                cat_plot_boxplot( eval(sprintf('{ %s(TIV<1300); %s(TIV>1300 & TIV<1500); %s(TIV>1500) }',...
                  measures{mi},measures{mi},measures{mi})) , ...
                  struct('names',{{'<1.3','~1.4','>1.5'}},'style',4,'ylim',myylim,'ygrid',0, ...
                  'datasymbol','o','usescatter',1,'groupcolor',acols) ); 
                title('TIV groups'); 
                xlabel('TIV (liter)'); ylabel(measures{mi});
              end
              ax  = gca;  ax.FontSize = 10; hold on; 
              ylim( myylim ); hold off; set(gca,'Ygrid',1,'YTickLabelRotation',0,'XTickLabelRotation',0);
              if strcmp(measures{mi},'rGMV'), set(gca,'YTickLabel',num2str(get(gca,'YTick')','%0.2f') ); end
              if useratings && any(contains(measures{mi},ratings))
                set(gca,'Ytick',scaling{mi}(1)+5 : 10 : scaling{mi}(2))
              end
        
              %% second figure with poxplot of centers
              subplot('position',[0.83 0.16 .16 .75]);
              if ai == 1 
                eval(sprintf('[an0,anova1site{mi}] = anova1(%s,site,''off'');',measures{mi})); 
                tables.siteanova = [ tables.siteanova; [ measures(mi) anova1site{mi}(2,4)' anova1site{mi}(2,5)' ] ]; 
                cat_plot_boxplot( eval(sprintf('{ %s(site==1) ; %s(site==2); %s(site==3) }',...
                  measures{mi},measures{mi},measures{mi})) , ...
                  struct('names',{{'Guys','HH','IOP'}},'style',4,'ylim',myylim ,'ygrid',0, ... siteu
                  'datasymbol','o','usescatter',1,'groupcolor',scols) ); 
                title('scan cites');  
                xlabel('center'); ylabel(measures{mi});
              else
                if ~exist('mwwtest','file')
                  warning('Miss mwwtest function (https://github.com/dnafinder/mwwtest/blob/master/mwwtest.m) use anova1.')
                  eval(sprintf('[an0,anova1sex{mi}] = anova1(%s,sex,''off'');',measures{mi})); 
                  tables.sexanova = [ tables.sexanova; [ measures(mi) anova1sex{mi}(2,5:6) ] ]; 
                else 
                  % https://github.com/dnafinder/mwwtest/blob/master/mwwtest.m
                  txt = evalc( sprintf('mwwt = mwwtest( %s(sex==0) , %s(sex==1) );',measures{mi},measures{mi}));
                  if isfield(mwwt,'Z')
                    tables.sexmwwt = [ tables.sexmwwt; [ measures(mi) mwwt.p(2) mean(mwwt.U) mwwt.Z] ]; 
                  else
                    tables.sexmwwt = [ tables.sexmwwt; [ measures(mi) mwwt.p(2) mean(mwwt.U) nan] ]; 
                  end
                end
                cat_plot_boxplot( eval(sprintf('{ %s(sex==0) ; %s(sex==1); }',...
                  measures{mi},measures{mi})) , ...
                  struct('names',{{'female','male'}},'style',4,'ylim',myylim ,'ygrid',0, ... siteu
                  'datasymbol','o','usescatter',1) ); 
                title('sex');  
                xlabel('sex'); ylabel(measures{mi});
              end
              ax  = gca;  ax.FontSize = 10; hold on; 
              if strcmp(measures{mi},'IQR') || strcmp(measures{mi},'SIQR')
                ph = plot([0 5],[1.9 1.9]); ph.Color = [1 0 0 ]; ph.LineStyle = '--'; ph.LineWidth = 1.5;     
              end
              ylim( myylim ); 
             
              if useratings && any(contains(measures{mi},ratings))
                set(gca,'Ytick',scaling{mi}(1)+5 : 10 : scaling{mi}(2))
              end
              set(gca,'Ygrid',1,'XTickLabelRotation',0);
              if strcmp(measures{mi},'rGMV'), set(gca,'YTickLabel',num2str(get(gca,'YTick')','%0.2f') ); end
                
    
              % print figure
              pdir = fullfile(resdir,sprintf('%s',measures{mi})); 
              if ~exist(pdir,'dir'), mkdir(pdir); end
              print(fh,fullfile(pdir,sprintf('IXI_%s_%s_%s_%s',avar{ai}, measures{mi}, qaversions{qai}, segment{si})),'-r1200','-djpeg');
            end
          end
        end
      end
    end
    
    % print final table
    tf = fieldnames(tables); 
    for tfi = 1:numel(tf)
      cat_io_csv(fullfile(resdir,sprintf('IXI_%s_%s_%s.csv', tf{tfi}, qaversions{qai}, segment{si})),tables.(tf{tfi})); 
    end
  end
  fprintf('Print IXI done. \n')

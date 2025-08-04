function cat_tst_aging(Pmethod, Presdir, Praw, subsets)
  %% - prepare and add further test sites 
  % Aspects and expectations: 
  % - all results depend on the underlying segmentation
  % - evaluation of atrophy in (i) aging ... and (ii) disease (e.g. AD,MCI,HC) ... controled?
  %   higher correlation and clear decline over different protocols with 
  %   less outiers
  % - evaluation of field strength 
  %   less effects by scanner or protocol 
  % - use intensity values? 
 

  % TODO: 
  % - handling of many curvs >> additional subfigure plot (1-3 lines)
  % - save tables 
  % - regression estimation on healthy only
  % - run multiple fit functions and select the best fitting one 
  % - check steps: 
  %   - name missing files (for reprocessing) or just to count the failures
  %   - add SPM parameter tests >> to improve segmentation?
  % - use different age ranges: 
  %   - full 
  %   - development (<18)
  %   - adult (>=18)
  %   - elderly and Alzheimers (>50)
  %   
  % - integrate further variables (problem is - what do we expect? )
  %   - intensity 
  %   - curvature
  %   - global GI? 
  %   
  % - extend super young and super old
  % - add non-healty in ageing plot 
  % - add further cohorts plots
  % - need real tests > plot of stat histograms
  % - create subfunctions
  %

  % SIDE ISSUE
  %#ok<*AGROW>
  addpath(fullfile(fileparts(which('cat12')),'internal')); 
  



  %  Find preocessed data that is available for all methods. 
  %  ----------------------------------------------------------------------
  datas = 'pbt'; Pdata = cell(1,size(Pmethod,1),1);
  for pi = 1:size(Pmethod,1)
    if pi == 1
      if ~exist(Pmethod{pi,2},'dir') && ~exist(fullfile(Pmethod{pi,2},subsets{sdi}),'dir')
        error('Need to start with a SPM/CAT directory with subdirectories.\n'); 
      end
      for sdi = 1:numel(subsets)
        Pdata{pi} = [ Pdata{pi}; cat_vol_findfiles(fullfile(Pmethod{pi,2}, subsets{sdi}), ...
          sprintf('lh.%s.*',datas),struct('depth',1)) ]; % pbt
      end
    elseif exist( fullfile(Pmethod{pi,2},'01_Bestoff'),'dir')
      Pdata{pi} = cat_io_strrep( Pdata{1}, Pmethod{1,2}, Pmethod{pi,2}) ; 
    else
      Pdata{pi} = spm_file( Pdata{1}, 'path', fullfile(Pmethod{pi,2},'surf') ); 
    end
    % find missing files
    Pexist(:,pi) = cellfun(@(x) exist(x,'file'),Pdata{pi})  &  ...
      cellfun(@(x) exist(x,'file'),cat_io_strrep(Pdata{pi},[filesep 'lh.'],[filesep 'rh.']));
  end
  % remove missing files in all testsets
  for pi = 1:size(Pmethod,1)
    Pdata{pi}( sum(Pexist==0,2)>0 ) = [];
  end
  % print report
  fprintf('\nFound %d/%d cases for evaluation of %d methods: \n', numel(Pdata{1}), size(Pdata,1), size(Pmethod,1));
  for pi = 1:size(Pmethod,1)
    fprintf('%16s: %4d files\n', Pmethod{pi,1}, nnz(Pexist(:,pi)>0) ); 
  end
  fprintf('----------------------------\n'); 
  fprintf('%16s: %4d files\n\n', 'overlap', length(Pdata{pi}) ); 
  



  %% define further files
  %  ----------------------------------------------------------------------
  for pi = 1:size(Pmethod,1)
    
    Pcentral{pi} = cat_io_strrep( Pdata{pi}, datas, 'central'); 
    Pwhite{pi}   = cat_io_strrep( Pdata{pi}, datas, 'white'); 
    Ppial{pi}    = cat_io_strrep( Pdata{pi}, datas, 'pial'); 
    Pgmv{pi}     = spm_file( cat_io_strrep( Pdata{pi}, ...
      {[filesep 'surf'  filesep], [filesep sprintf('lh.%s.',datas)]}, ...
      {[filesep 'mri'   filesep],  filesep}), 'ext', '.nii','prefix','mwp1');
    Pp0{pi}     = strrep( Pgmv{pi} , [filesep 'mwp1'], [filesep 'p0']); 
    if all( cellfun( @(x) exist(x,'file') , Pgmv{pi} )==0 )
      Pgmv{pi} = strcat(Pgmv{pi},'.gz'); 
    end
    if all( cellfun( @(x) exist(x,'file')==0, Pp0{pi} ) )
      Pp0{pi} = strcat(Pp0{pi},'.gz'); 
    end
    Pxml{pi}     = spm_file( cat_io_strrep( Pdata{pi}, ...
      {[filesep 'surf'   filesep], [filesep sprintf('lh.%s.',datas)]}, ...
      {[filesep 'report' filesep],  filesep}), 'ext', '.xml','prefix','cat_');
    if any(contains(Pmethod{pi,1},{'S25','SPM'})), Pxml{pi} = cat_io_strrep(Pxml{pi},[filesep 'cat_'],[filesep 'cat_c1']); end
    if any(contains(Pmethod{pi,1},{'S25','SPM'})), Pgmv{pi} = cat_io_strrep(Pgmv{pi},[filesep 'mwp1'],[filesep 'mwp1c1']); end
    Pff{pi}      = strrep(spm_str_manip(Pcentral{pi}','tr'),'lh.central',''); 

    % get mini xml if no xml is available 
    cat_io_createXML( struct('data', {Pp0{pi}( ...
      cellfun(@(x) exist(x,'file')==0, Pxml{pi})  &  cellfun(@(x) exist(x,'file')==2, Pp0{pi}) )} ) ); 
    Pexistx(:,pi) = cellfun(@(x) exist(x,'file'),Pxml{pi});
  end
  for pi = 1:size(Pmethod,1)
    Pdata{pi}( sum(Pexistx==0,2)>0 ) = [];
    Pwhite{pi}( sum(Pexistx==0,2)>0 ) = [];
    Ppial{pi}( sum(Pexistx==0,2)>0 ) = [];
    Pgmv{pi}( sum(Pexistx==0,2)>0 ) = [];
    Pp0{pi}( sum(Pexistx==0,2)>0 ) = [];
    Pff{pi}( sum(Pexistx==0,2)>0 ) = [];
    Pxml{pi}( sum(Pexistx==0,2)>0 ) = [];
  end
  Pdemo = spm_file( cat_io_strrep( Pdata{1}, {sprintf('lh.%s.',datas); Pmethod{1,2}}, {'',Praw}), 'ext','.xml');




  %% extract basic demographic parameters
  %  ----------------------------------------------------------------------
  age = []; sex = []; siten = {}; project = {}; cohortn = []; scanner = {}; 
  for si = 1:numel(Pdemo)
    demo  = cat_io_xml( Pdemo{si} ); 
    FN    = fieldnames(demo); 
    FNp   = FN{ contains( FN ,'patternCog_') };
    try age(si,1) = demo.(FNp).Age; catch, age(si,1) = -1; end
    try sex(si,1) = demo.(FNp).Sex; catch, sex(si,1) = str2double(cat_io_strrep(demo.(FNp).Sex,{'''';'F';'M'},{'','0','1'})); end
    try cohortn{si,1} = strrep( demo.(FNp).Cohort ,'''',''); catch, cohort{si,1} = 'NA'; end
    try field(si,1)   = demo.(FNp).Fieldstrength; catch, field(si,1) = -1; end
    try scanner{si,1} = strrep( demo.(FNp).Scanner ,'''',''); catch, scanner{si,1} = ''; end
    try
      project{si,1} = cat_io_strrep(FNp,{'patternCog_';'_demographics_T1w'},{'',''}); 
      if isnumeric( demo.(FNp).Site )
        siten{si,1}   = [ project{si,1} '_' sprintf('%d',demo.(FNp).Site)];
      else
        siten{si,1}   = [ project{si,1} '_' strrep( demo.(FNp).Site , '''', '')];
      end
    catch
      project{si,1} = cat_io_strrep(FNp,{'patternCog_';'_demographics_T1w'},{'',''});
      siten{si,1}   = [ project{si,1} "_NA"];
    end
  end
  site = str2double(cat_io_strrep( siten , ...
    unique(siten), cellfun(@(x) {num2str(x)},num2cell(1:numel(unique(siten))))) );
  cohort = str2double(cat_io_strrep( cohortn , ...
    unique(cohortn), cellfun(@(x) {num2str(x)},num2cell(1:numel(unique(cohortn))))) ); %#ok<NASGU>
  HC = contains(cohortn,'HC'); 

  % ##### CSV import for special cases #####
  % MR-ART?
  if 0
    csv    = cat_io_csv(Pcsv{pi}{1}); 
    csv    = sortrows(csv,1);
    csvf   = contains(csv(2:end,1), files); 
    IDage  = contains(csv(1,:),'Age'); 
    IDsex  = contains(csv(1,:),'Sex'); 
    IDsite = contains(csv(1,:),'Site'); 
    age    = [age  cell2mat(csv(2:end,IDage)) ]; 
    sex    = [sex  cell2mat(csv(2:end,IDsex))==1 ]; 
    site   = [site (sdi*1000 + str2double(cat_io_strrep( csv(2:end,IDsite) , ...
      unique(csv(2:end,IDsite)), cellfun(@(x) {num2str(x)},num2cell(1:numel(unique(csv(2:end,IDsite))))) ) ) ) ]; 
  end
    

  % update result directory
  Presdir = fullfile(Presdir,sprintf('aging_%0.0f-%0.0f_N%d_M%d',min(age),max(age),numel(age),size(Pmethod,1))); 
  if ~exist(Presdir,'dir'), mkdir(Presdir); end


  % #### 
  %  What about folding depending thickness? - We would expect less variation (although sulci fund are probably sig. thinner) 
  % 
  tiv = zeros(numel(Pdata{pi}),size(Pmethod,1)); rGMV = tiv; tsa = tiv; 
  for pi = 1:size(Pmethod,1)
    fprintf('\n  Method %2d/%2d: %16s', pi, size(Pmethod,1), Pmethod{pi} ); 
    
    % estimate average thickness
    for si = 1:numel(Pdata{pi})
      txt = sprintf(' Load subjects %4d/%4d', si, numel(Pdata{pi}));
      if si>1, fprintf(repmat('\b',1,numel(txt))); end
      fprintf(txt');

      cdata{pi,si}   = cat_io_FreeSurfer('read_surf_data',Pdata{pi}{si}); 
      mnthick(si,pi) = cat_stat_nanmean(cdata{pi,si}); 
    end
  
    % CAT XML data
    xml{pi} = cat_io_xml( Pxml{pi} ); 
    for si = 1:numel(Pdata{pi})
      try
        tiv(si,pi)  = xml{pi}(si).subjectmeasures.vol_TIV;
      catch
        tiv(si,pi)  = nan;
      end
      try 
        if isfield( xml{pi}(si).subjectmeasures,'vol_rel_CGW')
          rGMV(si,pi)  = xml{pi}(si).subjectmeasures.vol_rel_CGW(2);
        elseif isfield( xml{pi}(si).subjectmeasures,'vol_CGW_rel')
          rGMV(si,pi)  = xml{pi}(si).subjectmeasures.vol_CGW_rel(2);
        end
      catch
        rGMV(si,pi) = nan;
      end
      try 
        tsa(si,pi)   = xml{pi}(si).subjectmeasures.surf_TSA;
      catch
        try
          tsa(si,pi) = xml{pi}(si).subjectmeasures.TSA / 100;
        catch
          tsa(si,pi) = nan;
        end
      end
    end
  end
  
  fprintf('\n'); 
  

%% Pdata{pi}
  onlyAllreadySmoothed = 0; %#####################
  for pi=1:size(Pmethod,1)
    
%    matlabbatch = cat_tst_agereg( Pgmv{pi}, Presdir, Pmethod{pi,1}, 'GMV', age, sex, cohort==1, tiv(:,pi), 1, onlyAllreadySmoothed);
%    if ~isempty(matlabbatch), spm_jobman('run',matlabbatch); end

% thickness not working
    for sm = {'pbt','thickness'}
      %%
      Pdatapi = strrep( Pdata{pi} ,sprintf('lh.%s.',datas), sprintf('lh.%s.',sm{1})); 
      matlabbatch = cat_tst_agereg( Pdatapi, Presdir, Pmethod{pi,1}, sm{1}, age, sex, HC, tiv(:,pi), 1, onlyAllreadySmoothed);
      if ~isempty(matlabbatch), spm_jobman('run',matlabbatch); end
    end
  end
  
  
  
  %% global statistics
  %  ----------------------------------------------------------------------
  fprintf('Print figures: \n'); 
  fh = figure(50); set(fh,'Visible',1,'Interruptible',0); 
  %%
  for mi = 1:4 
    % fitting models: 
    % - not working: cubicinterp
    % - not fitting: log10
    % - too simple:  logistic4
    % - too fitting: smoothingspline
    % >> test the poly# fit and use the best fitting one
    dfit   = {'poly1','poly2','poly3','poly4','poly5'}; 
    roundf = 1;  
    xrange = [ floor( prctile(age(:),.1)/5)*5 , ceil( prctile(age(:),99.9)/5)*5 ];
    switch mi
      case 1
        data   = mnthick; 
        ylab   = 'cortical thickness (mm)';
        tlab   = 'CT'; 
        range  = [1.5 4];
      case 2 
        data   = rGMV; 
        ylab   = 'relative gray matter volume'; 
        tlab   = 'rGMV';
        range  = [0.2 .7];
      case 3 
        data   = tsa; 
        ylab   = 'TSA (total surface area)'; 
        tlab   = 'TSA';
        %dfit   = {'poly2'};
        range  = [1000 3000]; 
        roundf = -3;
      case 4
        data   = (mnthick .* tsa) ./ tiv; 
        ylab   = 'relative surface volume';
        tlab   = 'rCTV ';
        range  = [2 5];
    end
    if isempty(range)
      range  = round( [ prctile(data(:),.1) , prctile(data(:),99.9) ], roundf-1);
    end
    fprintf('  Measure %d - %s:\n',mi,tlab); 
  
    
    for spi = 1:2
      clf(fh); 
    
      fh.Position(3:4) = [600*spi 600];
      if spi == 2
        tiledlayout( round(size(Pmethod,1)/3) , ...
          ceil(size(Pmethod,1) /  round(size(Pmethod,1)/3)), ...
          'TileSpacing', 'compact', 'Padding', 'compact'); 
      end
    
      % print just for legend
      if spi == 1 
        fhagex = axes(fh); hold(fhagex,'on'); %#ok<LAXES>
        for pi = 1:size(Pmethod,1)
          sh = scatter(fhagex, age(1),data(1,pi));
          sh.MarkerFaceColor = Pmethod{pi,3};
          sh.MarkerEdgeColor = Pmethod{pi,3};
          sh.MarkerFaceAlpha = .2; 
          sh.MarkerEdgeAlpha = .3; 
          sh.Marker          = Pmethod{pi,4};
        end
      end

      % main print
      for pi = 1:size(Pmethod,1)
        if spi == 2
          fhagex = nexttile; hold(fhagex,'on');
        end

        sh = scatter(fhagex, age(HC),data(HC,pi)); 
        sh.MarkerFaceColor = Pmethod{pi,3};
        sh.MarkerEdgeColor = Pmethod{pi,3};
        sh.MarkerFaceAlpha = .1; 
        sh.MarkerEdgeAlpha = .2; 
        sh.Marker          = Pmethod{pi,4};
      
        if ~all( isnan(data(:,pi)) ) 
          ndnan = ~isnan(data(:,pi)) & HC; 
          %%
          [afit{pi}, agof{pi}, afitr{pi}, agofr{pi}, bfit{pi}] = bestfit(age(ndnan(:)) , data(ndnan(:),pi),dfit);
          arsq(pi)   = agof{pi}.rsquare; 
          arsqr(pi)  = agofr{pi}.rsquare; 
          acor(pi)   = corr( age(ndnan) , data(ndnan,pi) ,'Type','Spearman');  
        
          if spi == 2  ||  size(Pmethod,1)<5
            xint = linspace(0, 100,100);
            CIF = predint(afit{pi},xint,0.95,'Functional');
            CIO = predint(afit{pi},xint,0.95,'obs');
            p = fill(fhagex, [xint'; flip(xint')],[CIO(:,1); flip(CIO(:,2))], ...
              Pmethod{pi,3}); p.FaceAlpha = .1; p.EdgeAlpha = .2; p.EdgeColor = Pmethod{pi,3};
            p = fill(fhagex, [xint'; flip(xint')],[CIF(:,1); flip(CIF(:,2))], ...
              Pmethod{pi,3}); p.FaceAlpha = .3; p.EdgeAlpha = .4; p.EdgeColor = Pmethod{pi,3}; 
          end

          fph           = plot( fhagex, afit{pi} );
          fph.Color     = Pmethod{pi,3};
          fph.LineWidth = 1;
          fph.MarkerFaceColor = fph.Color;
        else
          acor(pi)   = nan; 
          arsq(pi)   = nan;
        end
      
        if spi==2
          grid(fhagex, 'on'); box(fhagex, 'on'); 
          ylim(range); xlim(xrange); 
          xlabel(fhagex, 'age (years)'); ylabel(fhagex, ylab);
          title(fhagex, sprintf('%s (Nsub=%d)',Pmethod{pi,1} ,sum(ndnan)));
          lgd = legend( fhagex, strcat( ...
            ' (', bfit{pi} ,', r=', num2str(acor(pi),'%0.3f'), ', r^2=', num2str(arsq(pi),'%0.3f'), ')' ));  %, ', rr^2=', num2str(arsqr(pi),'%0.3f')
        end
        
      end
      % format
      if spi==1
        grid(fhagex, 'on'); box(fhagex, 'on'); 
        xlabel(fhagex, 'age (years)'); ylabel(fhagex, ylab);
        title(fhagex, sprintf('%s (Nsub=%d, Nsites=%d)',tlab,numel(age),numel(unique(site))));
        subtitle(sprintf('[ r=corr(Spearman), r^2-fit, rr^2-robustfit ]')); 
        lgd = legend( fhagex, strcat( Pmethod(:,1), ...
          ' (', bfit', ' ,r=', num2str(acor','%0.2f'), ', r^2=', num2str(arsq','%0.2f'), ')' )); %, ', rr^2=', num2str(arsqr','%0.2f')
        fhagex.FontSize = 12; lgd.FontSize = 12;
      end
      ylim(range); xlim(xrange); 
      if spi==2, spis = 'single'; else, spis = 'mixed'; end 
      print(fh,fullfile(Presdir,sprintf('cat_tst_aging_%s%s',spis,deblank(tlab))),'-r300','-dpng');
      fprintf('    %s\n',fullfile(Presdir,sprintf('cat_tst_aging_%s%s.png',spis,deblank(tlab)))); 
    end



    %matlabbatch = cat_tst_aging(age,data,); 
    %spm('run',matlabbatch); 




    %% fieldstrength diffs
    % regress out age from data
    for pi = 1:size(Pmethod,1)
      % Get regression coefficients
      b = regress(data(HC,pi), [ones(size(age(HC), 1), 1), age(HC)]);
      % Predicted Y from X
      d_hat(:,pi) = [ones(size(age, 1), 1), age] * b;
      % Residuals = part of Y not explained by X
      data_residual(:,pi) = data(:,pi) - d_hat(:,pi);
    end



    %% test of field strength difference
    %  - this is not really working as it does not consider the diffent slopes in aging
    %    (eg. a method that did not show aging might have better results here)
    clear name dstrength methodfielddiff
    dsi = 0; fields = [1.5 3.0]; subname = '';
    for pi = 1:size(Pmethod,1)
      for si = 1:numel(fields)
        dsi = dsi + 1; 
        dstrength{dsi} = data_residual( age>5 & age<95 & field==fields(si) ,pi); 
        name{dsi} = sprintf('%sT%0.1f',Pmethod{pi},fields(si));
        if si == numel(fields)
          [~,ttsp(pi),tts] = ttest2( dstrength{dsi-1} , dstrength{dsi} ); 
          methodfielddiff(pi) = diff(tts);
        end
      end
      subname = sprintf('%s   %s(tst p=%0.3f)',subname,Pmethod{pi,1},ttsp(pi));
    end
    mylim = round( prctile(abs(cell2mat(dstrength')),99),roundf) * 1.2;
    clf(fh); fh.Position(3:4) = [500 200 + 50*size(Pmethod,1)];
    fhfieldx = axes(fh); %#ok<LAXES>
    cat_plot_boxplot( dstrength , struct('names',{name},'ygrid',1, 'ylim',[-mylim mylim],...
      'subsets',mod( 1+round((1:size(Pmethod,1)*2)/2),2) , ...
      'groupcolor',min(1,max(0,cell2mat(reshape( [Pmethod(:,3)'; Pmethod(:,3)'] , ...
      size(Pmethod,1)*2, 1)) + repmat([-.2;+.2],size(Pmethod,1),1))))); 
    title(fhfieldx, sprintf('Method vs. fieldstrength (Nsub=%d, Nsites=%d, higher p-values are better)', ...
      numel(age), numel(unique(site))) ); subtitle(fhfieldx,subname);
    xlabel(fhfieldx,'Method + fieldstrength'); ylabel(fhfieldx,ylab,'Rotation',90); 
    print(fh,fullfile(Presdir,sprintf('cat_tst_fieldstrength_%s',deblank(tlab))),'-r300','-dpng');
    fprintf('    %s\n',fullfile(Presdir,sprintf('cat_tst_fieldstrength_%s.png',deblank(tlab)))); 



    %% group-diffs 
    %  - also not optimal
    clear name dstrength methodfielddiff tts ttsp
    projects = {'ADNI','AIBL','OASIS3'}; 
    dsi = 0; subname = ''; cohorts = {'HC','MCI','AD'}; %unique(cohortn);
    for pi = 1:size(Pmethod,1)
      % regress out age from data
      % Get regression coefficients
      for si = 1:numel(cohorts)
        dsi = dsi + 1; 
        dstrength{dsi} = data_residual( age>50 & contains(cohortn, cohorts{si}) & contains(project,projects),pi); 
        name{dsi} = sprintf('%sT%s',Pmethod{pi},cohorts{si});
        if si > 1
          [~,ttsp(pi,si-1),tts{pi,si-1}] = ttest2( dstrength{dsi-1} , dstrength{dsi} ); 
        end
        if si == numel(cohorts)
          [~,ttsp(pi,si),tts{pi,si}] = ttest2( dstrength{dsi+(1-numel(cohorts))} , dstrength{dsi} );
          methodfielddiff(pi) = mean( diff([tts{pi,:}]) );
        end
      end
      subname = sprintf('%s   %s(p=%0.3f)', subname, Pmethod{pi,1}, mean(ttsp(pi,end)));
    end
    if any( cellfun(@(x) ~isempty(x), dstrength )) 
      clf(fh); fh.Position(3:4) = [600 400];
      mylim = round( prctile(abs(cell2mat(dstrength')),99),roundf) * 1.2; 
      fhcohortx = axes(fh); %#ok<LAXES>
      cat_plot_boxplot( dstrength , struct('names',{name},'ygrid',1,'maxwhisker',1.5,'ylim',[-mylim mylim],...
        'handle',fh,'subsets',mod( floor( (0:size(Pmethod,1)*3-1)/3) ,2), ...
        'groupcolor',repmat( cat_io_colormaps('trafficlight',3), size(Pmethod,1), 1) ));
      title(fhcohortx, sprintf('Alzheimer''s in %s (tst HC vs. AD, Nsub=%d, Nsites=%d, lower p-values are better)', ...
        [projects{:}], numel(age), numel(unique(site))) ); subtitle(fhcohortx, subname);
      xlabel(fhcohortx, 'Method+Cohort'); ylabel(fhcohortx, ylab,'Rotation',90); hold on;
      % ... save ...
      print(fh, fullfile(Presdir,sprintf('cat_tst_Alzheimers_%s',deblank(tlab))),'-r300','-dpng');
      fprintf('    %s\n',fullfile(Presdir,sprintf('cat_tst_Alzheimers_%s.png',deblank(tlab)))); 
   end



    %% sex differences (in puberty)
    clf(fh); 
    fh.Position(3:4) = [1200 800]; sgcolor = [.5 0 0; 0 0 .5]; sgmarker = {'o','^'};
    tiledlayout( round(size(Pmethod,1)*2) , ...
      ceil(size(Pmethod,1) /  round(size(Pmethod,1)*2)), ...
      'TileSpacing', 'compact', 'Padding', 'compact'); 
    for pi = 1:size(Pmethod,1)
      nexttile; 
      hold on; 
      for si = 1:2
        sg = HC & sex==(si-1) & ~isnan(data(:,pi)); 
        if sum(sg)>0
          [afit{pi}, agof{pi}, afitr{pi}, agofr{pi}, bfit{pi}] = bestfit(age(sg(:)) , data(sg(:),pi),dfit);
          arsq(pi)   = agof{pi}.rsquare; 
          arsqr(pi)  = agofr{pi}.rsquare; 
          afitp1(pi) = afit{pi}.p1; 
          afitp2(pi) = afit{pi}.p2;
          acor(pi)   = corr( age(sg) , data(sg,pi) ,'Type','Spearman');  
        
          xint = linspace(0, 100,100);
          CIF = predint(afit{pi},xint,0.80,'Functional');
          CIO = predint(afit{pi},xint,0.80,'obs');
          p = fill([xint'; flip(xint')],[CIO(:,1); flip(CIO(:,2))], sgcolor(si,:)); p.FaceAlpha = .1; p.EdgeAlpha = .2; p.EdgeColor = sgcolor(si,:);
          p = fill([xint'; flip(xint')],[CIF(:,1); flip(CIF(:,2))], sgcolor(si,:)); p.FaceAlpha = .2; p.EdgeAlpha = .3; p.EdgeColor = sgcolor(si,:); 
     
          sh = scatter(age(sg),data(sg,pi));
          sh.Marker = sgmarker{si}; 
          sh.MarkerFaceColor = sgcolor(si,:); sh.MarkerFaceAlpha = .1;
          sh.MarkerEdgeColor = sgcolor(si,:); sh.MarkerEdgeAlpha = .1; 
          
          fph           = plot( afit{pi} );
          fph.Color     = sgcolor(si,:);
          fph.LineWidth = 1;
  
          legend off; 
  
          title(sprintf('sex difference %s',Pmethod{pi,1})); %subtitle(subname);
          xlim(xrange); ylim(range); 
          xlabel('age'); ylabel(ylab); box on; grid on;
        end
      end
      ylim(round([prctile(data(:),2) prctile(data(:),98)] .* [.8 1.1],1))
    end
    print(fh,fullfile(Presdir,sprintf('cat_tst_sex_%s',deblank(tlab))),'-r300','-dpng');
    fprintf('    %s\n',fullfile(Presdir,sprintf('cat_tst_sex_%s.png',deblank(tlab)))); 
    
  end
  fprintf('done.\n')
end
function [afit,agof,afitr,agofr,bfit] = bestfit(x,y,dfit) %#ok<INUSD>
%#ok<*NODEF>
  for fi = 1:numel(dfit)
    evalc('[afitr{fi}, agofr{fi}]  = fit( x , y, dfit{fi},''Robust'',''Bisquare'');');  %,struct('robust','LAR'));Bisquare
    evalc('[afit{fi},  agof{fi}]   = fit( x , y, dfit{fi});');
    if strcmp(lastwarn,'Equation is badly conditioned. Remove repeated data points or try centering and scaling.')
      lastwarn(''); 
      rs(fi)  = 0; 
    else
      rs(fi)  =  agofr{fi}.rsquare; %mean( [ agof{fi}.rsquare, agofr{fi}.rsquare ]); 
    end
  end   
  [~,bfid] = max(rs); 

  afit  = afit{bfid}; 
  agof  = agof{bfid}; 
  afitr = afitr{bfid}; 
  agofr = agofr{bfid}; 
  bfit  = dfit{bfid}; 
end
function matlabbatch = SPMpreprocessing4bc(BWPfilesSPM)
%-----------------------------------------------------------------------
% Job saved on 07-Mar-2025 14:46:45 by cfg_util (rev $Rev: 8183 $)
% spm SPM - SPM25 (25.01.02)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
  matlabbatch{1}.spm.spatial.preproc.channel.vols = BWPfilesSPM;
  matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.0001;
  matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
  matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
  matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {fullfile(spm('dir'),'tpm','TPM.nii,1')};
  matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
  matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [0 0];
  matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
  matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {fullfile(spm('dir'),'tpm','TPM.nii,2')};
  matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
  matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [0 0];
  matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
  matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {fullfile(spm('dir'),'tpm','TPM.nii,3')};
  matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
  matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [0 0];
  matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
  matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {fullfile(spm('dir'),'tpm','TPM.nii,4')};
  matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
  matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
  matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
  matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {fullfile(spm('dir'),'tpm','TPM.nii,5')};
  matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
  matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
  matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
  matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {fullfile(spm('dir'),'tpm','TPM.nii,6')};
  matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
  matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
  matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
  matlabbatch{1}.spm.spatial.preproc.warp.mrf = 0;
  matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 0;
  matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0 0.1 0.01 0.04];
  matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
  matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
  matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
  matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
  matlabbatch{1}.spm.spatial.preproc.warp.vox = NaN;
  matlabbatch{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
                                                NaN NaN NaN];
end
function CS1 = loadSurf(P)
  if ~exist(P,'file'), error('Surface file %s could not be found due to previous processing errors.',P); end 
  
  try
    CS = gifti(P);
  catch
    error('Surface file %s could not be read due to previous processing errors.',P);
  end
  
  CS1.vertices = CS.vertices; CS1.faces = CS.faces; 
  if isfield(CS,'cdata'), CS1.cdata = CS.cdata; end
end
function borderintensity( Praw , Pdata)
%% outsource into function 
  % native data gunzipping
  [~,~,ee] = fileparts(Pdata{pi}{si});
  Prawsi   = cat_vol_findfiles(Praw,sprintf('%s.nii',ee(2:end)));
  if isempty(Prawsi)
    Prawgz = cat_vol_findfiles(Praw,sprintf('%s.nii.gz',ee(2:end)));
    gunzip(Prawgz{1})
    delraw = 1; 
    Prawsi   = cat_vol_findfiles(Praw,sprintf('%s.nii',ee(2:end)));
  else
    delraw = 0; 
  end

  %% bias-corrected > but not intensity normalized :/ 
  Pbcsi = cat_vol_findfiles(Praw,sprintf('m%s.nii',ee(2:end)));
  if numel(Pbcsi) < numel(Prawsi) 
    BWPfilesSPM = Prawsi; 
    BWPfilesSPM( cellfun(@(x) exist(x,'file'),spm_file(BWPfilesSPM,'prefix','m'))>0 ) = [];
    if ~isempty( BWPfilesSPM )
      matlabbatchbc = SPMpreprocessing4bc( BWPfilesSPM ) ;
      if verbose, spm_jobman('run',matlabbatchbc); else, evalc( 'spm_jobman(''run'',matlabbatchbc);'); end
    end
  end

  %% 
  Pbcsi = cat_vol_findfiles(Praw,sprintf('m%s.nii',ee(2:end)));
  Pgsi  = spm_file(Pbcsi,'prefix','g'); 
  if ~exist(Pgsi{1},'file')
    %%
    Vbcsi = spm_vol(char(Pbcsi)); 
    Ybcsi = spm_read_vols(Vbcsi);
    Ygsi  = cat_vol_grad(Ybcsi,1) ./ max( prctile( Ybcsi(:),90) , cat_vol_smooth3X(Ybcsi,4) ) * 3; 
    Ygsi(Ygsi>1) = max(0,1 - 1/3*(Ygsi(Ygsi>1)));
    Ygsi  = log10(Ygsi * 9 + 1); 
    cat_sanlm(Ygsi,1,3); 
    Vgsi  = Vbcsi; Vgsi.fname = Pgsi{1};
    spm_write_vol(Vgsi,Ygsi); 
  end

  
  S = gifti(Pcentral{pi}{si}); 
  %%
  surfname = {'white','pial','layer4','cortex'}; 
  for bi = 1:3
    clear matlabbatch; 
    matlabbatch{1}.spm.tools.cat.stools.vol2surf.data_vol        = Pbcsi; % Pgsi
    matlabbatch{1}.spm.tools.cat.stools.vol2surf.datafieldname   = 'int';
    if bi == 1
      matlabbatch{1}.spm.tools.cat.stools.vol2surf.data_mesh_lh  = Pwhite{pi}(si);
    elseif bi == 2
      matlabbatch{1}.spm.tools.cat.stools.vol2surf.data_mesh_lh  = Ppial{pi}(si);
    else
      matlabbatch{1}.spm.tools.cat.stools.vol2surf.data_mesh_lh  = Pcentral{pi}(si);
      matlabbatch{1}.spm.tools.cat.stools.vol2surf.datafieldname = sprintf('int_%s',surfname{bi});
    end
    matlabbatch{1}.spm.tools.cat.stools.vol2surf.sample          = {'mean'};
    matlabbatch{1}.spm.tools.cat.stools.vol2surf.interp          = {'linear'};
    matlabbatch{1}.spm.tools.cat.stools.vol2surf.mapping.rel_equivol_mapping.class      = 'GM';
    if bi == 4
      matlabbatch{1}.spm.tools.cat.stools.vol2surf.sample          = {'weighted_avg'}; %'multi'
      matlabbatch{1}.spm.tools.cat.stools.vol2surf.mapping.rel_equivol_mapping.startpoint = -0.5;
      matlabbatch{1}.spm.tools.cat.stools.vol2surf.mapping.rel_equivol_mapping.steps      = 7;
      matlabbatch{1}.spm.tools.cat.stools.vol2surf.mapping.rel_equivol_mapping.endpoint   = 0.5;
    else
      matlabbatch{1}.spm.tools.cat.stools.vol2surf.mapping.rel_equivol_mapping.startpoint = 0;
      matlabbatch{1}.spm.tools.cat.stools.vol2surf.mapping.rel_equivol_mapping.steps      = 1;
      matlabbatch{1}.spm.tools.cat.stools.vol2surf.mapping.rel_equivol_mapping.endpoint   = 0;
    end

% request for gradient or STD?        
    if verbose, spm_jobman('run',matlabbatch); else, evalc( 'spm_jobman(''run'',matlabbatch);'); end

    %%
    if bi == 3
      Pint{pi}{si}{bi,1} = strrep( spm_str_manip( Pcentral{pi}{si}, 'r'), '.central.', ...
        sprintf('.%s_%s.', char(matlabbatch{1}.spm.tools.cat.stools.vol2surf.datafieldname), ...
        char(spm_str_manip(matlabbatch{1}.spm.tools.cat.stools.vol2surf.data_vol,'tr'))));
    else
      Pint{pi}{si}{bi,1} = strrep( spm_str_manip( Pcentral{pi}{si}, 'r'), '.central.', ...
        sprintf('.%s_%s_%s.', char(matlabbatch{1}.spm.tools.cat.stools.vol2surf.datafieldname), ...
        surfname{bi}, char(spm_str_manip(matlabbatch{1}.spm.tools.cat.stools.vol2surf.data_vol,'tr'))));
    end
    %cdata = cat_io_FreeSurfer('read_surf_data',  Pint{pi}{si}{bi} );
    Pintp{si}{bi}{pi,1} = Pint{pi}{si}{bi,1}; 
  end
  
%%
  if 0 % need intnorm and orient!
    Ybcsi = spm_read_vols(spm_vol(char(Pbcsi))); 
    rres = cat_surf_fun('evalCS', ...
      loadSurf(Pcentral{pi}{si}), cat_io_FreeSurfer('read_surf_data',Pdata{pi}{si}), cat_io_FreeSurfer('read_surf_data',Pdata{pi}{si}), ...
      Ybcsi,Ybcsi,Pcentral{pi}{si},[],2,cat_get_defaults('extopts.expertgui')>1);
  end
%%
  if delraw, delete(Prawsi{1}); end


    if 0
    %%
    datalim = 3; 
    [hst1_int,hst2_int] = cat_plot_histogram( char([Pintp{1}{1:datalim}]) , struct( ...
      'color', min(1,repmat( cat_io_colormaps('nejm',size(Pintp{1}{1},1)),datalim,1) + ...
      repmat( round(1:size(Pintp{1}{1},1)*datalim)' / size(Pintp{1}{1},1)*datalim / 20 -.05,1,3) ) ));
    %% thickness
    [hst1_data,hst2_data] = cat_plot_histogram( char([Pdata{1}(1); Pdata{2}(1); Pdata{3}(1)]) , struct( ...
      'color', repmat( cat_io_colormaps('nejm',size(Pintp{1}{1},1)),1,1)));
  end

end
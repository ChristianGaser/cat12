function cat_tst_qa_MRART_expertgroups( datadir, qaversions, segment, fasttest, rerun ) 
%% cat_tst_qa_MRART_expertgroups. QC evaluation script
%  ------------------------------------------------------------------------
%  Script to evaluate the CAT QC measures on the MR-ART dataset with motion
%  affected groups (1 - no MA, 2 - light MA, and 3 - strong MA).
%
%  Requirements: 
%   0. Download and install SPM and CAT
%   1. Download MR-ART data from: 
%        
%
%   2. Specify in this script: 
%      1) the data directory "datadir" 
%      2) the QC version you would like to tests (the file has to exist in the cat directory) 
%      3) the segmentation you would like to use
%
%  ------------------------------------------------------------------------

cat_io_cprintf([0 0.5 0],'\n\n== Run cat_tst_qa_MRART_expertgroups ==\n') 

% ### datadir ###
if ~exist( 'datadir' , 'var' )
  maindir  = '/Volumes/WDE18TB/MRData/Dahnke2025_QC/ds004173-download'; 
else
  maindir  = fullfile(datadir,'ds004173-download'); 
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


[cv,rn]   = cat_version;
catppdir  = fullfile(maindir,'derivatives',sprintf('%s',cv)); %sprintf('%s_R%s',cv,rn)); 
[~,CATver,ext] = fileparts(catppdir); CATver = [CATver ext]; clear ext
mriqcdir  = fullfile(maindir,'derivatives','mriqc-0.16.1');
resultdir = fullfile(maindir,'derivatives',['results_' CATver]);
printdir  = fullfile(fileparts(maindir),'+results',sprintf('MR-ART_%s_%s',fast{fasttest+1},datestr(clock,'YYYYmm')));
exprating = fullfile(maindir,'derivatives','scores.tsv');
partis    = fullfile(maindir,'participants.tsv');
FS        = [10 10];
if fasttest, pres = '-r300'; else, pres = '-r1200'; end
% translate marks into percentage ratingop
mark2rps  = @(mark) min(100,max(0,105 - mark*10));
verb      = 0; 

si = 1; 
qais = 1:numel(qaversions); 


% preprocessing
for si = 1:numel(segment)
  clear matlabbatch; 
  ARTfiles = cat_vol_findfiles( maindir , 'sub*.nii',struct('depth',3));
  switch segment{si}
    case 'CAT'
      CATpreprocessing4qc;
      ARTfilesCAT = ARTfiles;
      for fi = 1:numel(ARTfiles)
        p0files{fi} = spm_file(ARTfiles{fi},'path',fullfile(maindir,'derivatives','CAT12.9',spm_str_manip(ARTfiles{fi},'hht'),'anat'),'prefix','p0');
      end
      ARTfilesCAT( cellfun(@(x) exist(x,'file'),p0files)>0 ) = [];
      if ~isempty( ARTfilesCAT )
        matlabbatch{1}.spm.tools.cat.estwrite.data = ARTfilesCAT;
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.lazy = 1; 
        matlabbatch{1}.spm.tools.cat.estwrite.output.BIDS = ...
          struct('BIDSyes',struct('BIDSfolder',fullfile('..','derivatives','CAT12.9')));
        spm_jobman('run',matlabbatch);
      end
    case 'SPM'
      SPMpreprocessing4qc;
      ARTfilesSPM = ARTfiles;
      ARTfilesSPM( cellfun(@(x) exist(x,'file'),spm_file(ARTfiles,'prefix','c1'))>0 ) = [];
      if ~isempty( ARTfilesSPM )
        matlabbatch{1}.spm.spatial.preproc.channel.vols = ARTfilesSPM;
        spm_jobman('run',matlabbatch);
      end
    case 'synthseg'
      error('synthseg is not prepared in the public script ... use MATLAB help')
    case 'qcseg'
      fprintf('No preprocessing required.\n\n');
  end
end



%%
for qai = qais
  Pp0{qai} = cat_vol_findfiles( catppdir , 'p0*.nii');  
  Pc1{qai} = cat_vol_findfiles( maindir  , 'c1*.nii');  
  Pss{qai} = cat_vol_findfiles( maindir  , 'synthseg_p0*.nii');  
  Pqs{qai} = cat_vol_findfiles( maindir  , 'sub*.nii');  
  if fasttest, Pp0{qai} = Pp0{qai}(1:4:end); end
  if fasttest, Pc1{qai} = Pc1{qai}(1:4:end); end
  if fasttest, Pss{qai} = Pss{qai}(1:4:end); end
  if fasttest, Pqs{qai} = Pqs{qai}(1:4:end); end
  %% (re)processing of QC values
  if 1 %rerun
    for si = 1:numel(segment)
      switch segment{si}
        case 'CAT'
          cat_vol_qa('p0',Pp0{qai},struct('prefix',[qaversions{qai} '_'],'version',qaversions{ qai },'rerun',rerun*2));
        case 'SPM' 
          cat_vol_qa('p0',Pc1{qai},struct('model',struct('spmc1',1),'prefix',[qaversions{qai} '_spm_'],'version',qaversions{ qai },'rerun',rerun*2));
        case 'synthseg'
          cat_vol_qa('p0',Pss{qai},struct('prefix',[qaversions{qai} '_synthseg_'],'version',qaversions{ qai },'rerun',rerun*2));
        case 'qcseg'
          cat_vol_qa('p0',Pqs{qai},struct('prefix',[qaversions{qai} '_qcseg_'],'version',qaversions{ qai },'rerun',rerun*2));
      end
    end
  end
end


CSV = cat_io_csv(partis,struct('delimiter',' ')); 

%% create figure
%segment = 'SPM';%SPM'; %'synthseg';
for si = 1:numel(segment)
  for qai = qais
    fprintf('MRART %s:\n',qaversions{qai})
  
    resultdirqai = [resultdir '_'  fast{fasttest+1} '_' datestr(clock,'YYYYmm') ]; 
  
    if ~exist(resultdirqai,'dir'), mkdir(resultdirqai); end
    if ~exist(printdir    ,'dir'), mkdir(printdir); end
    
    % find CAT XML files - catppdir
    clear Pxml Pxmlspm; 
    if exist('Pp0','var')
      for pi=1:numel(Pp0{qai})
        [pp,ff] = fileparts(Pp0{qai}{pi});
        if contains( segment{si} ,'synthseg')
          [pp,ff] = fileparts(Pss{qai}{pi});
          Pxml{pi,1} = fullfile(pp,'report',[qaversions{qai} '_' strrep(ff,'synthseg_p0','synthseg_') '.xml']);
        elseif contains( segment{si} ,'SPM')
          [pp,ff] = fileparts(Pc1{qai}{pi});
          Pxml{pi,1} = fullfile(pp,[qaversions{qai} '_spm_' ff(3:end) '.xml']); % spm_c1
        elseif contains( segment{si} ,'qcseg')
          [pp,ff] = fileparts(Pqs{qai}{pi});
          Pxml{pi,1} = fullfile(pp,'report',[qaversions{qai} '_qcseg_' ff '.xml']); % spm_c1
        else
          Pxml{pi,1} = fullfile(maindir,spm_str_manip(pp,'ht'),spm_str_manip(pp,'t'),'report',[qaversions{qai} '_' ff(3:end) '.xml']);
          if ~exist(Pxml{pi},'file')
            Pxml{pi,1} = fullfile(maindir,spm_str_manip(pp,'ht'),[qaversions{qai} '_' ff(3:end) '.xml']);
          end
          if ~exist(Pxml{pi},'file')
            Pxml{pi,1} = fullfile(pp,'report',[qaversions{qai} '_' ff(3:end) '.xml']);
          end
        end
        [pp,ff] = fileparts(Pc1{qai}{pi});
        Pxmlspm{pi,1} = fullfile(pp,[qaversions{qai} '_spm_' ff(3:end) '.xml']);
      end
    else
      Pxml = cat_vol_findfiles(fullfile(maindir),sprintf('%s_sub*.xml',qaversions{qai})); % ,'derivatives' 
      if fasttest && numel(Pp0)>numel(Pxml), Pxml = Pxml(1:4:end); end
    end
    xml    = cat_io_xml(Pxml); 
    if contains( segment{si} ,'CAT')
      xmlspm = cat_io_xml(Pxmlspm); 
    end
  
    % get json and other files
    clear Q;
    %%
    for xi = 1:numel(Pxml)
  
      % get json and other files
      if any( contains(segment{si},'CAT' ) )
        [pp,ff,ee]     = fileparts(strrep(Pxml{xi},[qaversions{qai} '_'],''));
        SID{xi}        = spm_str_manip(Pxml{xi},'hhht');
      elseif any( contains(segment{si},'SPM' ) )
        [pp,ff,ee]     = fileparts(strrep(strrep(Pxml{xi},[qaversions{qai} '_'],''),'spm_',''));
        SID{xi}        = spm_str_manip(Pxml{xi},'hht');
      elseif any( contains(segment{si},'qcseg' ) )
        [pp,ff,ee]     = fileparts(strrep(strrep(Pxml{xi},[qaversions{qai} '_'],''),'qcseg_',''));
        SID{xi}        = spm_str_manip(Pxml{xi},'hhht');
      else
        [pp,ff,ee]     = fileparts(strrep(strrep(Pxml{xi},[qaversions{qai} '_'],''),'synthseg_',''));
        SID{xi}        = spm_str_manip(Pxml{xi},'hhht');
      end
      jsonfile{xi,1} = fullfile(mriqcdir,SID{xi},'anat',[ff '.json']); 
      
      % get age/sex from csv
      csvSIDi = find(contains(CSV(:,1),SID{xi}));
      xml(xi).age = CSV{ csvSIDi , 3}; 
      xml(xi).sex = CSV{ csvSIDi , 2}=='M'; 
      Q.age(1,xi) = xml(xi).age; 
      Q.sex(1,xi) = xml(xi).sex; 
      
      %% get SPM 
      [pp,ff] = fileparts(Pp0{qai}{xi}); pp = strrep(pp,catppdir,maindir); 
      spm_vol_mat = fullfile(pp,['spmvol_' ff(3:end) '.mat']); 
      if exist(spm_vol_mat,'file')
        load(spm_vol_mat,'spmvol'); 
      else
        for ci = 1:3
          P = fullfile(pp,sprintf('c%d%s.nii',ci,ff(3:end))); 
          V = spm_vol(P); 
          Y = spm_read_vols(V); 
          vx_vol = sqrt(sum(V.mat(1:3,1:3).^2));
          spmvol(ci) = sum(Y(:)) * prod(vx_vol)/1000; 
        end
        save(spm_vol_mat,'spmvol'); 
      end
      spmvols(xi,1:3) = spmvol / sum(spmvol); 
  
      if contains(segment{si},'CAT')
        xml(xi).subjectmeasures.SPMrGMV  = spmvols(xi,1); 
        xml(xi).subjectmeasures.SPMrCSFV = spmvols(xi,3); 
        xml(xi).subjectmeasures.SPMrWMV  = spmvols(xi,2); 
      else
        xml(xi).subjectmeasures.vol_rel_CGW(1:3)  = spmvols(xi,[3,1,2]); 
      end
  
      %% get synthseg 
      if contains(segment{si},'synthseg')
        [pp,ff] = fileparts(Pp0{qai}{xi}); pp = strrep(pp,catppdir,maindir); 
        synthseg_vol_mat = fullfile(pp,['synthsegvol_' ff(3:end) '.mat']); 
        Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));
        if exist(synthseg_vol_mat,'file')
          load(synthseg_vol_mat,'synthsegvol'); 
        else
          P = fullfile(pp,sprintf('synthseg_%s.nii',ff)); 
          V = spm_vol(P); 
          Y = spm_read_vols(V); 
          vx_vol = sqrt(sum(V.mat(1:3,1:3).^2));
          for ci = 1:3
            synthsegvol(ci) = sum(Yp0toC(Y(:),ci)) * prod(vx_vol)/1000; 
          end
          save(synthseg_vol_mat,'synthsegvol'); 
        end
        synthsegvols(xi,1:3) = synthsegvol / sum(synthsegvol);
  
        xml(xi).subjectmeasures.vol_rel_CGW(1:3)  = synthsegvols(xi,:); 
      end     
  
      
  
    end
    json = cat_io_json(jsonfile); 
    
    
    %% inverse measure
    evalmriqc = 0;
  
    mriqcQFNinv = {'cjv','efc','qi_1','inu','snr_CSF' ...
      'wm_p95_','wm_stdv_','wm_median_','wm_mean_','wm_mad_','bg_p95','bg_stdv' ...
      'bg_mad','bg_mean','bg_median','bg_p95','bg_stdv', ...
      'csf_mean','csf_median','csf_n', ...
      'gm_mad','gm_mean','gm_median','gm_p95','gm_stdv', ...
      'wm_mad','wm_mean','wm_median','wm_p95','wm_stdv', ...
      }; 
    if 0
      mriqcQFN    = {
        'snr_total'; 'snr_wm'; % excelent
        'cnr'; 'cjv'; 'fber'; % good  
        ...'summary_bg_p95'; 'summary_bg_stdv'; 
        ...'summery_gm_p95'; 'summery_gm_k'; 'summery_gm_mad'; 'summery_gm_mean';
        ...'summery_wm_p95'; 'summery_wm_k'; 'summery_wm_mad'; 'summery_wm_mean';
      };
    else
      mriqcQFN    = {
        'cnr'; 'cjv'; 'fber'; 'efc'; 'qi_1'; 'qi_2'; 'wm2max'; 
        'fwhm_avg'; 'fwhm_x'; 'fwhm_y'; 'fwhm_z';
        'icvs_csf'; 'icvs_gm'; 'icvs_wm';
        'inu_med'; 'inu_range';
        'rpve_csf'; 'rpve_gm'; 'rpve_wm'; 
        'snr_csf'; 'snr_gm'; 'snr_total'; 'snr_wm'; 
        'snrd_csf'; 'snrd_gm'; 'snrd_total'; 'snrd_wm'; 
        'summary_bg_k';  'summary_bg_mad';  'summary_bg_mean';  'summary_bg_median';  'summary_bg_n';  'summary_bg_p05';  'summary_bg_p95';  'summary_bg_stdv';
        'summary_csf_k'; 'summary_csf_mad'; 'summary_csf_mean'; 'summary_csf_median'; 'summary_csf_n'; 'summary_csf_p05'; 'summary_csf_p95'; 'summary_csf_stdv';
        'summary_gm_k';  'summary_gm_mad';  'summary_gm_mean';  'summary_gm_median';  'summary_gm_n';  'summary_gm_p05';  'summary_gm_p95';  'summary_gm_stdv';
        'summary_wm_k';  'summary_wm_mad';  'summary_wm_mean';  'summary_wm_median';  'summary_wm_n';  'summary_wm_p05';  'summary_wm_p95';  'summary_wm_stdv';
        'tpm_overlap_csf'; 'tpm_overlap_gm'; 'tpm_overlap_wm';
        }'; 
    end
    for fni = 1:numel(mriqcQFN)
      for xi = 1:numel(Pxml)
        Q.(mriqcQFN{fni})(xi,1) = json(xi).(mriqcQFN{fni}); 
      end
      if any(contains( mriqcQFN(fni) , mriqcQFNinv ))
        Q.(mriqcQFN{fni}) = -Q.(mriqcQFN{fni}); 
      end
    end
    
    % read expert rating (rename tsv as csv)
    copyfile(exprating,strrep(exprating,'.tsv','.csv')); 
    tsv = cat_io_csv(exprating,'','',struct('delimiter','\t'));
  
    
      
  
    % extract IQR values
    % QFN: { fieldname1 fieldname2 fieldindex outlier range name print }
    QFN = {
      'qualityratings'   'NCR'          1 1 [1 6]         'NCR'       1
      'qualityratings'   'ICR'          1 1 [1 6]         'ICR'       1
      ...'qualityratings'   'res_RES'      1 1 [1 6]         'res_RES'   0-fasttest
      'qualityratings'   'res_ECR'      1 1 [1 6]         'res_ECR'   1
      'qualityratings'   'FEC'          1 1 [1 6]         'FEC'       1
      ...'qualityratings'   'IQR'          1 1 [1 6]         'IQR'       1
      'qualityratings'   'SIQR'         1 1 [1 6]         'SIQR'      1
     ... 'qualityratings'   'contrastr'    1 1 [1 6]         'CON'       1
      'subjectmeasures'  'vol_rel_CGW'  2 0 [0 1]         'rGMV'      1 %-fasttest
      'subjectmeasures'  'vol_rel_CGW'  1 0 [0 1]         'rCSFV'      1-fasttest
      'subjectmeasures'  'vol_rel_CGW'  3 0 [0 1]         'rWMV'      1-fasttest
      ...'subjectmeasures'  'SPMrGMV'      1 0 [0 1]         'SPMrGMV'      1-fasttest
      ...'subjectmeasures'  'SPMrWMV'      1 0 [0 1]         'SPMrCSFV'      1-fasttest
      ...'subjectmeasures'  'SPMrCSFV'     1 0 [0 1]         'SPMrWMV'      1-fasttest
      ...'subjectmeasures'  'vol_abs_CGW'  1 0 [400 800]     'CSFV'      0
      ...'subjectmeasures'  'vol_abs_CGW'  2 0 [400 800]     'GMV'       0
      ...'subjectmeasures'  'vol_abs_CGW'  3 0 [400 800]     'WMV'       0
      ...'subjectmeasures'  'vol_TIV'      1 0 [800 2000]    'TIV'       1-fasttest
      };
    QS = Q; 
    for xi = 1:numel(xml)
        Q.name{xi,1} = xml(xi).filedata.file;
        QS.name{xi,1} = xml(xi).filedata.file;
        for fni1 = 1:size(QFN,1)
          if isfield(xml(xi),QFN{fni1,1})  &&  isfield(xml(xi).(QFN{fni1,1}),QFN{fni1,2})
            Q.(QFN{fni1,6})(xi,1)   = xml(xi).(QFN{fni1,1}).(QFN{fni1,2})(QFN{fni1,3});
          else
            Q.(QFN{fni1,6})(xi,1)   = nan;
          end
          if exist('xmlspm','var')
            if xi<=numel(xmlspm)
              if isfield(xmlspm(xi),QFN{fni1,1})  &&  isfield(xmlspm(xi).(QFN{fni1,1}),QFN{fni1,2})
                QS.(QFN{fni1,6})(xi,1)   = xmlspm(xi).(QFN{fni1,1}).(QFN{fni1,2})(QFN{fni1,3});
              else
                QS.(QFN{fni1,6})(xi,1)   = nan;
              end
            end
          end
        end
    end
    
  if 0
    Q = QS; 
  %  switch segment{si}
   %   case 'SPM', qaversions{ qai } = [qaversions{ qai } '_spm'];
    %qaversions{ qai } = [qaversions{ qai } '_synthseg'];
  end
  Q.site  = ones(size(Q.NCR));
    
  
    
    % get the motion groups by instruction 
    for xi = 1:numel(xml)
        Q.sub{xi,1}  = Q.name{xi}(1:10); 
        
        if ~isempty( strfind(Q.name{xi},'standard') )
          Q.group0(xi,1) = 1; 
        elseif ~isempty( strfind(Q.name{xi},'headmotion1_T1w') )
          Q.group0(xi,1) = 2; 
        elseif ~isempty( strfind(Q.name{xi},'headmotion2_T1w') )
          Q.group0(xi,1) = 3; 
        end
    end
   
     
    if evalmriqc == 2, QFN = {}; qaversions{qai} = 'mriQC'; end
    if evalmriqc
      for fni = 1:numel(mriqcQFN)
        QFN = [ QFN ; {'' mriqcQFN{fni} 1 1 ...
          [min(Q.(mriqcQFN{fni})) max(Q.(mriqcQFN{fni}))] mriqcQFN{fni} 1 } ]; 
      end
    end
  
  
    %%
    if 0
      clear QSD; 
      fh = figure(22);
      fh.Visible = 'off';
      fh.Interruptible = 'off';
      fh.Position(3:4) = [200 200];
      QM = {'NCR','ICR','res_RES','res_ECR','FEC','SIQR'};   
      cl = [.8 0 .2;   0.8 0.6 0;  0.2 0.5 0;  .0 .5 .9; 0. 0 0.9;  0 0 0 ]; 
      for fni1 = 1:numel(QM)
        QSD.(QM{fni1}) = [Q.(QM{fni1}), QS.(QM{fni1})];
        QSD.MAE{fni1}  = Q.(QM{fni1}) - QS.(QM{fni1});
        QSD.RMSE(fni1) = mean(QSD.MAE{fni1}.^2).^.5;
      end
      cat_plot_boxplot( QSD.MAE,struct('ygrid',0,'style',4,'names',{cat_io_strrep(QM,'res_','')},'datasymbol','o', ...
              'usescatter',1,'groupcolor',cl));
      set(gca,'ygrid','on','ylim',[-3 3],'XTickLabelRotation',0);
      title(sprintf('MAE (%s)',segment{si})), xlabel('measures'); ylabel('error')
      if ~exist(fullfile(printdir,'boxplot'),'dir'), mkdir(fullfile(printdir,'boxplot')); end
      print(fh,fullfile(printdir,'boxplot',['MRART_boxplot_MAE_' qaversions{ qai } '_' segment{si} '.png']),'-dpng',pres)
      %%   
      bh = bar(QSD.RMSE); bh.CData = cl; bh.FaceColor = 'flat';
      ylim([0 3]); xlim([.4 6.6]); xticklabels({'NCR','ICR','RES','ECR','FEC','SIQR'}); 
      h = gca; h.YGrid = 'on';  h.XTickLabelRotation = 0; 
      title(sprintf('RMSE %s',segment{si})); xlabel('measures'); ylabel('error')
      for fi = 1:numel(QSD.RMSE)
        dt(fi) = text(fi-.45, QSD.RMSE(fi) + .1, sprintf('%0.3f',QSD.RMSE(fi)),'FontSize',8,'Color',cl(fi,:));
      end  
      set(gca,'ygrid','on','ylim',[0 3],'XTickLabelRotation',0);
      if ~exist(fullfile(printdir,'boxplot'),'dir'), mkdir(fullfile(printdir,'boxplot')); end
      print(fh,fullfile(printdir,'boxplot',['MRART_boxplot_RMSE_' qaversions{ qai } '_' segment{si} '.png']),'-dpng',pres)
    end
  
  
    %% find alignment xml and tsv
    mytable.ROC   = repmat({{'measure','AUC','ACC'}},1,4); 
    mytable.ANOVA = {'measure','F-value','p-value'}; 
    mark2rps = @(mark) min(100,max(0,105 - mark*10));
    gradlim  = [40 100]; % gradlim = [0.5 9.5]; 
    gradtick = 40:10:100; % 45:10:95; 
    for fni1 = 1:size(QFN,1)
      data = cell(1,3); 
      for ti = 2:size(tsv,1)
          name  = tsv{ti,1};
          group = tsv{ti,2};
    
          sid  = find( strcmp(Q.name,name),1);
          if ~isempty(sid) 
            if strcmp(QFN{fni1,1},'qualityratings')
              data{ group }(end+1)   = mark2rps( Q.(QFN{fni1,6})( sid )); 
            else
              data{ group }(end+1)   = Q.(QFN{fni1,6})( sid ); 
            end
            Q.group( sid , 1)      = group; 
          end
      end
      
      %if all(cell2mat(cellfun(@isnan,data,'UniformOutput',false)))
      %  continue
      %end
      
  
      %%
      if QFN{fni1,7}
        % ANOVA
        if fni1 == 1
          try
            % call stat toolbox licence
            [an_p(fni1), an_tab{fni1}, an_obj{fni1}] = anova1( mark2rps(Q.(QFN{fni1,6})) , uint8(Q.group) , 'off' ); 
          end
          pause(1); 
        end
        if strcmp(QFN{fni1,1},'qualityratings')
          [an_p(fni1), an_tab{fni1}, an_obj{fni1}] = anova1( mark2rps(Q.(QFN{fni1,6})) , uint8(Q.group) , 'off' ); 
        else
          [an_p(fni1), an_tab{fni1}, an_obj{fni1}] = anova1( Q.(QFN{fni1,6}) , uint8(Q.group) , 'off' ); 
        end
      
        % print boxplot figure
        if verb, fig = figure(40); else, fig = figure(); end
        fig.Visible = 'off';
        fig.Interruptible = 'off'; 
        fig.Position(3:4) = [130 200]; 
        fig.Name = sprintf('MR-ART - Boxplot - %s %s',qaversions{qai},strrep(QFN{fni1,6},'_','\_'));  
        if ~verb, fig.Visible = 'off'; else, fig.Visible = 'on'; end
  
        if strcmp(QFN{fni1,1},'qualityratings')
          cat_plot_boxplot(data,struct('ygrid',0,'style',4,'names',{{'no','light','strong'}},'datasymbol','o', ...
            'usescatter',1,'groupcolor',[0 0.8 0; 0.8 0.7 0; 0.8 0 0],'ylim',gradlim));
        else
          cat_plot_boxplot(data,struct('ygrid',0,'style',4,'names',{{'no','light','strong'}},'datasymbol','o', ...
            'usescatter',1,'groupcolor',[0 0.8 0; 0.8 0.7 0; 0.8 0 0]));
        end
        %%
        title(sprintf('MA groups %s', strrep(QFN{fni1,6},'_','\_'))); box on; 
        subtitle(sprintf('F=%0.1f, p=%0.1e', an_tab{fni1}{2,5}, an_p(fni1)))
        mytable.ANOVA = [mytable.ANOVA; 
          QFN(fni1,6), an_tab{fni1}(2,5), {an_p(fni1)}]; 
        if strcmp(QFN{fni1,1},'qualityratings')
          set(gca,'ygrid','on','ytick',gradtick,'XTickLabelRotation',0); 
        else
          set(gca,'ygrid','on','XTickLabelRotation',0); 
        end
        xlabel('groups'); ylabel(sprintf('%s (grades)',strrep(QFN{fni1,6},'_','\_'))); 
        if ~exist(fullfile(printdir,'boxplot'),'dir'), mkdir(fullfile(printdir,'boxplot')); end
        print(fig,fullfile(printdir,'boxplot',['MRART_boxplot_' strrep(QFN{fni1,6},'res_','') '_' qaversions{ qai } '_' segment{si} '.png']),'-dpng',pres)
        if ~verb, close(fig); end
      end
    
  
     
  
      %% age dependency
      fig = figure(12); fig.Position(3:4) = [300 200]; clf; 
      fig.Visible = 'off';
      fig.Interruptible = 'off'; 
      gcol = [0 0.8 0; 0.8 0.7 0; 0.8 0 0]; hold on; 
      gnam = {'no','light','strong'};
      if 0
        % select measure for tests
        fni2 = find(contains(QFN(:,2),'vol_rel_CGW')); if numel(fni2)>1, fni2 = fni2(1); end
      else
        fni2 = fni1; 
      end
      for gi = 1:max(Q.group)
        if strcmp(QFN{fni2,1},'qualityratings')
          sc = scatter(Q.age(Q.group==gi), mark2rps(Q.(QFN{fni2,6})(Q.group==gi)), 20); 
          [curve1{gi}, goodness, output] = fit( Q.age(Q.group==gi)', mark2rps(Q.(QFN{fni2,6})(Q.group==gi)),'poly1','robust','LAR');
        else
          if strcmp(QFN{fni2,6},'rWMV'), deg = 'poly2'; else, deg = 'poly1'; end
          sc = scatter(Q.age(Q.group==gi), Q.(QFN{fni2,6})(Q.group==gi), 20); 
          [curve1{gi}, goodness, output] = fit( Q.age(Q.group==gi & ~isnan(Q.(QFN{fni2,6}) ))', ...
            Q.(QFN{fni2,6})(Q.group==gi & ~isnan(Q.(QFN{fni2,6}))) , deg,'robust','LAR');
        end
        set(sc,'MarkerFaceColor',gcol(gi,:),'MarkerEdgeColor',gcol(gi,:), ...
          'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3); 
        mylegend{gi} = sprintf('%s (b=%0.3f/100a)',gnam{gi},curve1{gi}.p1 * 100); 
        Q.agefit.p1(gi) = curve1{gi}.p1; 
      end
      xlim([15 85]); 
      for gi = 1:max(Q.group)
        ph = plot(curve1{gi});
        set(ph,'Color',gcol(gi,:)); 
      end
      box on; grid on;
      title(sprintf('MA groups %s (diff(good/bad)=%0.3f)', strrep(QFN{fni2,6},'_','\_'),...
        abs(diff([mean(Q.(QFN{fni2,6})(Q.group==1)),mean(Q.(QFN{fni2,6})(Q.group==3))])))); 
      xlabel('age (years)'); 
      if strcmp(QFN{fni2,1},'qualityratings')
        ylabel(sprintf('%s (rps)',strrep(QFN{fni2,6},'_','\_'))); 
        ylim(gradlim - 10*(1-strcmp(QFN{fni2,2},'ICR'))); %set(gca,'YTick',45:10:95);
        legend(mylegend,'Location','SouthEast');
      else
        ylabel(sprintf('%s',strrep(QFN{fni2,6},'_','\_'))); 
        if mean(Q.(QFN{fni2,6})(Q.group==1)) > mean(Q.(QFN{fni2,6})(Q.group==3)) & ~strcmp(QFN{fni2,1},'subjectmeasures')
          legend(mylegend,'Location','SouthEast');
        else
          legend(mylegend,'Location','NorthEast');
        end
        
      end
      set(gca,'XTick',20:10:80);
      tdir = fullfile(printdir,'aging'); if ~exist(tdir,'dir'), mkdir(tdir); end
      print(fig, fullfile(tdir, sprintf('MRART_aging_%s_%s_%s.png', strrep(QFN{fni2,6},'res_','') , qaversions{qai} , segment{si})) , '-dpng',pres);
     
      
      
  
    
      
      %% outlier detection
      if fasttest, tss = 0.05; else tss = 0.01; end; ppos = [0 2 1 3];
      clear sens2 spec2 sensg2 specg2 cf
      if QFN{fni1,4}
        fig = figure; fig.Position(3:4) = [600 200]; 
        fig.Name = sprintf('MR-ART - ROC - %s',qaversions{qai});  
        if ~verb, fig.Visible = 'off'; end
  
        erths = [2.5 2.0 1.5]; % expert rating threshold for 3 groups (1=no MA, 2=light MAs, 3=severe MAs)
        for erthi = 1:numel(erths)
          erth = erths(erthi); 
          
          %  --------------------------------------------------------------------
          Q.train   = mod( round( 1/3:1/3:numel(Q.group)/3 )' , 2 ); % to get more or less all groups
          
          IQRfield  = QFN{fni1,6}; 
          model     = 10;
          cmodel    = 1; 
          rps       = 1; 
                
          if contains( QFN{fni1,2},  mriqcQFN)
            th      = min( Q.(QFN{fni1,2}) ):tss/10 * (max( Q.(QFN{fni1,2}) ) - min( Q.(QFN{fni1,2}) )):max( Q.(QFN{fni1,2}) );
            cf{erthi} = th; 
          else
            th      = 0.5:tss:6.5;   % global IQR threshold test range (school grades)
            cf{erthi} = -10:tss*10:30;     % protocoll-specific dIQR threshold test range (school grad range) 
          end
             
          %
          testsites     = 1; 
          sens2{erthi}  = {nan(numel(cf{erthi}),1),nan(numel(cf{erthi}),1)}; %#ok<*SAGROW> 
          spec2{erthi}  = sens2{erthi}; acc2{erthi} = sens2{erthi}; auc2{erthi} = sens2{erthi};
          Q.Nmn  = nan(size(Q.SIQR)); Q.Nsd = Q.Nmn;  Q.NXIQR = Q.Nmn;
          for i = 1:numel(cf{erthi}) % apply global IQR tresholds for ROC statistic 
            for ti = 1:2      % train vs. test
              M = Q.train == ti-1; nM = ~M; 
              switch erthi
                case 2, M(Q.group==2) = 0; nM(Q.group==2) = 0; 
              end
        
              if contains( QFN{fni1,2},  mriqcQFN)
                Q.NXIQR = Q.(IQRfield);
              else
                % get thresholds
                if rps
                  [ Q.NXIQR(M) , Q.Nmn(M) , Q.Nsd(M)] = cat_tst_qa_normer( mark2rps(  Q.(IQRfield)(M) ) ,  ...
                    struct('model',model,'figure',0,'cmodel',cmodel));
                else
                  [ Q.NXIQR(M) , Q.Nmn(M) , Q.Nsd(M)] = cat_tst_qa_normer(  Q.(IQRfield)(M) ,  ...
                    struct('model',model,'figure',0,'cmodel',cmodel));
                end
  
                % estimate for all values
                for gi = testsites
                  Q.Nmn(Q.site == gi) = cat_stat_nanmean( Q.Nmn(Q.site == gi) );
                  Q.Nsd(Q.site == gi) = cat_stat_nanmean( Q.Nsd(Q.site == gi) );
                end
                if cmodel == 1
                  if rps
                    Q.NXIQR = mark2rps(Q.(IQRfield)) - Q.Nmn;
                  else
                    Q.NXIQR = Q.(IQRfield) - Q.Nmn;
                  end
                else
                  if rps
                    Q.NXIQR = (mark2rps(Q.(IQRfield)) - Q.Nmn ) ./ Q.Nsd;
                  else
                    Q.NXIQR = (Q.(IQRfield) - Q.Nmn ) ./ Q.Nsd;
                  end
                end
              end
        
              TP = Q.group(M)> erth & Q.NXIQR(M) <  cf{erthi}(i); TPs = sum(TP);
              FP = Q.group(M)<=erth & Q.NXIQR(M) <  cf{erthi}(i); FPs = sum(FP);
              TN = Q.group(M)<=erth & Q.NXIQR(M) >= cf{erthi}(i); TNs = sum(TN);
              FN = Q.group(M)> erth & Q.NXIQR(M) >= cf{erthi}(i); FNs = sum(FN);
        
              sens2{erthi}{ti}(i) = TPs ./ max(1,TPs + FNs);
              spec2{erthi}{ti}(i) = TNs ./ max(1,TNs + FPs);
              acc2{erthi}{ti}(i)  = (TPs + TNs) / max(1,TPs + FNs + TNs + FPs);
        
              [~,~,~,auc2{erthi}{ti}(i)] = perfcurve( (Q.group(nM)>erth), cf{erthi}(i) - Q.NXIQR(nM), 'true');
            end
          end
          
          % get best value
          [~,mxcci1] = max( cat_stat_nanmean([ spec2{erthi}{1} , sens2{erthi}{1} ] , 2) ); % train
          [~,mxcci2] = max( cat_stat_nanmean([ spec2{erthi}{2} , sens2{erthi}{2} ] , 2) ); % test
          
  
          % figure
          % -----------------------------------------------------------------
          subplot('Position',[.06 + 0.33 * ppos(erthi) 0.12 0.27 0.77],'box','on'); cla; hold on; grid on;
          plot(1-spec2{erthi}{1} ,sens2{erthi}{1} ,'color',[0.8 0.0 0.6],'linewidth',0.5); 
          plot(1-spec2{erthi}{2} ,sens2{erthi}{2} ,'color',[0.0 0.4 0.8],'linewidth',0.5);
          plot(1-mean([spec2{erthi}{1},spec2{erthi}{2}],2) ,mean([sens2{erthi}{1},sens2{erthi}{2}],2) ,'color',[0.0 0.2 0.4],'linewidth',1);
          hold off; ylim([0.5 1.004]); xlim([-0.004 0.65]); ylim([0 1.0]); xlim([0 1])
          title(sprintf('ROC %s',cat_io_strrep(IQRfield,{'res_','_'},{'','\_'}))); 
          switch erthi
            case 1, subtitle(sprintf('no vs. light/severe MAs (n=%d)',sum(M))); 
            case 2, subtitle(sprintf('no/light vs. severe MAs (n=%d)',sum(M)));
            case 3, subtitle(sprintf('no vs. severe MAs (n=%d)',sum(M))); 
          end
          ax = gca; ax.FontSize = FS(2)*.8;
          xlabel('False positive rate (1-specificity)','FontSize',FS(2)*.9); 
          ylabel('True positive rate (sensitivity)','FontSize',FS(2)*.9); 
          lg = legend({... 
                  sprintf('run1: AUC=%0.3f\n (th=%0.2f rps: ACC=%0.3f)', auc2{erthi}{1}(mxcci2), cf{erthi}(mxcci2), acc2{erthi}{1}(mxcci2) ),...
                  sprintf('run2: AUC=%0.3f\n (th=%0.2f rps: ACC=%0.3f)', auc2{erthi}{2}(mxcci1), cf{erthi}(mxcci1), acc2{erthi}{2}(mxcci1) ),...
                  sprintf('avg.: AUC=%0.3f\n (th=%0.2f rps: ACC=%0.3f)',...
                    mean([auc2{erthi}{1}(mxcci2), auc2{erthi}{2}(mxcci1)]), ...
                    mean([cf{erthi}(mxcci2),      cf{erthi}(mxcci1)]), ...
                    mean([acc2{erthi}{1}(mxcci2), acc2{erthi}{2}(mxcci1)]))...
                 },'Location','southeast','Fontsize',FS(2)*.8); lg.Box = 'off';
          yticks(0:0.2:1);
          xticks(0:0.2:1);   
      
        end
       
       
        % print 
        tdir = fullfile(printdir,'ROC'); if ~exist(tdir,'dir'), mkdir(tdir); end
        print(fig, fullfile(tdir, sprintf('MRART_ROC_%s_%s_%s.png', strrep(IQRfield,'res_','') , qaversions{qai} , segment{si})) , '-dpng',pres);
        if ~verb, close(fig); end 
  
  
  
        % age dependency in normalized data
        fig = figure(12); fig.Position(3:4) = [300 200]; clf; 
        fig.Visible = 'off';
        fig.Interruptible = 'off'; 
        gcol = [0 0.8 0; 0.8 0.7 0; 0.8 0 0]; hold on; 
        gnam = {'no','light','strong'};
        for gi = 1:max(Q.group)
          sc = scatter(Q.age(Q.group==gi), Q.NXIQR(Q.group==gi), 20); 
          [curve1{gi}, goodness, output] = fit( Q.age(Q.group==gi)', Q.NXIQR(Q.group==gi) ,'poly1');
          set(sc,'MarkerFaceColor',gcol(gi,:),'MarkerEdgeColor',gcol(gi,:), ...
            'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3); 
          mylegend{gi} = sprintf('%s (b=%0.3f/100a)',gnam{gi},curve1{gi}.p1 * 100); 
          Q.agefit.p1(gi) = curve1{gi}.p1; 
        end
        xlim([15 85]); set(gca,'XTick',20:20:80); xlabel('age (years)'); 
        for gi = 1:max(Q.group)
          ph = plot(curve1{gi});
          set(ph,'Color',gcol(gi,:)); 
        end
        box on; grid on;
        title(sprintf('MA groups %s (diff(good/bad)=%0.3f)', strrep(QFN{fni1,6},'_','\_'),...
          abs(diff([mean(Q.(QFN{fni1,6})(Q.group==1)),mean(Q.(QFN{fni1,6})(Q.group==3))])))); 
        ylabel(sprintf('n%s (rps)',strrep(QFN{fni1,6},'_','\_'))); 
        legend(mylegend,'Location','SouthEast');
        
        tdir = fullfile(printdir,'aging_n'); if ~exist(tdir,'dir'), mkdir(tdir); end
        print(fig, fullfile(tdir, sprintf('MRART_aging_%s_%s_%s.png', strrep(QFN{fni1,6},'res_','') , qaversions{qai} , segment{si})) , '-dpng',pres);
     
  
  
  
  
        %% figure 2
        % -----------------------------------------------------------------
        fig = figure(); fig.Position(3:4) = [800 200];
        fig.Visible = 'off';
        fig.Interruptible = 'off'; 
        fig.Name = sprintf('MR-ART - ROC - %s',qaversions{qai});  
        if ~verb, fig.Visible = 'off'; end
  
        for i=1:2, spec2{4}{i} = mean([spec2{1}{i},spec2{2}{i},spec2{3}{i}],2); end
        for i=1:2, sens2{4}{i} = mean([sens2{1}{i},sens2{2}{i},sens2{3}{i}],2); end
        for i=1:2, auc2{4}{i}  = mean([auc2{1}{i},auc2{2}{i},auc2{3}{i}],2); end
        for i=1:2, acc2{4}{i}  = mean([acc2{1}{i},acc2{2}{i},acc2{3}{i}],2); end
        cf{4}   = mean([cf{1};cf{2};cf{3}],1); erths(4) = 3; 
        for erthi = 1:numel(erths)
          % get best value
          %[~,mxcci]  = max( cat_stat_nanmean([ spec2{erthi}{2} , sens2{erthi}{2} ] , 2) ); % train
          [~,mxcci1] = max( cat_stat_nanmean([ spec2{erthi}{1} , sens2{erthi}{1} ] , 2) ); % train
          [~,mxcci2] = max( cat_stat_nanmean([ spec2{erthi}{2} , sens2{erthi}{2} ] , 2) ); % train
          
  
          %erth = erths(erthi); % expert rating threshold for 3 groups (1=no MA, 2=light MAs, 3=severe MAs)
          subplot('Position',[.045 + 1/4 * ppos(erthi) 0.12 0.20 0.74],'box','on'); cla; hold on; grid on;
          plot(1-spec2{erthi}{1} ,sens2{erthi}{1} ,'color',min(1,[0.8 0.0 0.6]+.3),'linewidth',0.5); 
          plot(1-spec2{erthi}{2} ,sens2{erthi}{2} ,'color',min(1,[0.0 0.4 0.8]+.3),'linewidth',0.5);
          plot(1-mean([spec2{erthi}{1},spec2{erthi}{2}],2) ,mean([sens2{erthi}{1},sens2{erthi}{2}],2) ,'color',[0.0 0.2 0.4],'linewidth',1);
          hold off; ylim([0.5 1.004]); xlim([-0.004 0.65]); ylim([0 1.0]); xlim([0 1])
          title(sprintf('ROC %s',cat_io_strrep(IQRfield,{'res_','_'},{'','\_'}))); 
          switch erthi
            case 1, subtitle(sprintf('no vs. light/severe MAs (n=%d)',sum(M))); 
            case 2, subtitle(sprintf('no/light vs. severe MAs (n=%d)',sum(M)));
            case 3, subtitle(sprintf('no vs. severe MAs (n=%d)',sum(M))); 
            otherwise, subtitle(sprintf('average ROC of MAs',sum(M)));
          end
          ax = gca; ax.FontSize = FS(2)*.8;
          xlabel('False positive rate (1-specificity)','FontSize',FS(2)); 
          ylabel('True positive rate (sensitivity)','FontSize',FS(2)); 
          lg = legend({ ...
                  sprintf('run1: AUC=%0.3f\n (th=%0.2f rps: ACC=%0.3f)', auc2{erthi}{1}(mxcci2), cf{erthi}(mxcci2), acc2{erthi}{1}(mxcci2) ),...
                  sprintf('run2: AUC=%0.3f\n (th=%0.2f rps: ACC=%0.3f)', auc2{erthi}{2}(mxcci1), cf{erthi}(mxcci1), acc2{erthi}{2}(mxcci1) ),...
                  sprintf('avg.: AUC=%0.3f\n (th=%0.2f rps: ACC=%0.3f)', ...
                    mean([auc2{erthi}{1}(mxcci2), auc2{erthi}{2}(mxcci1)]), ...
                    mean([cf{erthi}(mxcci2),      cf{erthi}(mxcci1)]), ...
                    mean([acc2{erthi}{1}(mxcci2), acc2{erthi}{2}(mxcci1)]))...
                 },'Location','southeast','Fontsize',FS(2)*.8); lg.Box = 'off';
          yticks(0:0.2:1);
          xticks(0:0.2:1);      
        
          mytable.ROC{erthi} = [ mytable.ROC{erthi}; {IQRfield}, ...
            {mean([auc2{erthi}{1}(mxcci2),auc2{erthi}{2}(mxcci1)])}, {mean([acc2{erthi}{1}(mxcci2), acc2{erthi}{2}(mxcci1)])} ]; 
        end
        fprintf('  %s done.\n',IQRfield)
             
        
        
    
        % print 
        tdir = fullfile(printdir,'ROC4_only-for-fast-evaluation'); if ~exist(tdir,'dir'), mkdir(tdir); end
        print(fig, fullfile(tdir, sprintf('MRART_ROC4_%s_%s_%s.png', strrep(IQRfield,'res_','') , qaversions{qai} , segment{si})) , '-dpng',pres);
        close(fig); 
      end
  
  
  
    end
  
    
  
    
  
    %% scatter plot and correlation between dIQR and dGMV
    dlim = [0 4]; 
    if ~fasttest
      QM1 = find( contains( QFN(:,1) ,'subjectmeasures') );  
      %QM1 = find( contains( QFN(:,2) ,'SPM') ); 
      for fni1 = 1:numel(QM1)
      % * need to run for all QM compared to rGMV
      % * there is a lot of normal variance :(
    
        %if ~exist('fig4','var'), 
        if verb, fig = figure(41); else, fig = figure(); end; clf; 
        fig.Visible = 'off';
        fig.Interruptible = 'off'; 
        fig.Position(3:4) = [600 200];
        fig.Name = sprintf('MR-ART - rGMV changes - %s',qaversions{qai});  
        if ~verb, fig.Visible = 'off'; else, fig.Visible = 'on'; end
  
        QM2 = find( contains( QFN(:,1) ,'qualityratings') );
        for fni2 = 1:numel(QM2)
          dfield1     = ['d' QFN{QM1(fni1),2}];
          dfield2     = ['d' QFN{QM2(fni2),2}];
          dfield      = [dfield1 '_' dfield2]; 
          Q.(dfield)  = zeros(numel(Q.sub),1);
          Q.(dfield1) = zeros(numel(Q.sub),1);
          Q.(dfield2) = zeros(numel(Q.sub),1);
          for sxi = 1:numel(Q.sub)
            %%
            subids = find(contains(Q.sub,Q.sub{sxi}));
            G1     = find(Q.group0(subids)==1); 
            val1 = Q.(QFN{QM1(fni1),6})(subids(G1));
            val2 = Q.(QFN{QM2(fni2),6})(subids(G1));
          
            Q.(dfield1)(sxi,1) = Q.(QFN{QM1(fni1),6})(sxi) - val1;
            Q.(dfield2)(sxi,1) = Q.(QFN{QM2(fni2),6})(sxi) - val2;
          end
        
          %% robustfit
          subplot('Position',[0.07 + 0.90*(fni2-1)/numel(QM2) 0.15, 0.85/numel(QM2) .7]);
  
          plot(dlim,[0 0],'color',[0 0 0]); hold on; 
          hs = scatter( Q.(dfield2)( Q.(dfield2)>0 & Q.group==3) , Q.(dfield1)(Q.(dfield2)>0 & Q.group==3));
          hs.MarkerEdgeColor = [.8 0 .2]; hs.MarkerFaceColor =  hs.MarkerEdgeColor; 
          hs.MarkerFaceAlpha = .1; hs.MarkerEdgeAlpha = .15; hs.SizeData = 20; hs.Marker = 'v'; 
          hs = scatter( Q.(dfield2)( Q.(dfield2)>0 & Q.group==2) , Q.(dfield1)(Q.(dfield2)>0 & Q.group==2)); 
          hs.MarkerEdgeColor = [0.7 0.3 0]; hs.MarkerFaceColor =  hs.MarkerEdgeColor; 
          hs.MarkerFaceAlpha = .1; hs.MarkerEdgeAlpha = .15; hs.SizeData = 20;  hs.Marker = '^'; 
          ah = gca; ylim([-1,1] * max(ceil(abs(ah.YLim)*20)/20)); xlim(dlim); 
          ah.XTick = ah.XLim(1):round( diff([ah.XLim(1),ah.XLim(2)]) / 4  ,0):ah.XLim(2); 
          ah.YTick = ah.YLim(1):round( diff([ah.YLim(1),ah.YLim(2)]) / 10 ,2):ah.YLim(2); 
          
          if fni2>1, ylabel off; ah.YTickLabel = {}; end
          %Q.fit.(dfield) = fit( double(Q.(dfield2)(Q.(dfield2)>0 & Q.group==2)) , double(Q.(dfield1)(Q.(dfield2)>0 & Q.group==2)) ,'poly1'); hp = plot(Q.fit.(dfield)); hp.Color = [0 0.5 0];
          %Q.fit.(dfield) = fit( double(Q.(dfield2)(Q.(dfield2)>0 & Q.group==3)) , double(Q.(dfield1)(Q.(dfield2)>0 & Q.group==3)) ,'poly1'); hp = plot(Q.fit.(dfield)); hp.Color = [0.5 0 0];
          try
            Q.fit.(dfield) = fit( double(Q.(dfield2)(Q.(dfield2)>0 & ~isnan(Q.(dfield1)) & ~isnan(Q.(dfield2)))) , ...
              double(Q.(dfield1)( Q.(dfield2)>0 & ~isnan(Q.(dfield1)) & ~isnan(Q.(dfield2)))) ,'poly1'); 
            hp = plot(Q.fit.(dfield)); hp.Color = [0.6 0 0];
            legend off; grid on; box on; 
            if fni2==1, ylabel([QFN{QM1(fni1),6} ' change ']); else, ylabel(''); end 
            xlabel(strrep(dfield2,'_','\_'));
            subtitle(sprintf('%0.3f',Q.fit.(dfield).p1));
          end 
          title(sprintf('%s change',strrep(QFN{QM2(fni2),2},'_','\_')))
        end
   
        tdir = fullfile(printdir,'rGMVchanges'); if ~exist(tdir,'dir'), mkdir(tdir); end
        print(fig, fullfile(tdir, sprintf('MRART_%s_%s_%s.png', ['d' QFN{QM1(fni1),6}] , qaversions{qai}, segment{si} )) , '-dpng',pres);
        if ~verb, close(fig); end
      end
    end
    
    %%
    rig = {'no-lightSevere','noLight-severe','no-severe','average'}; 
    for ri = 1:4
      cat_io_csv( fullfile( printdir , sprintf('MRART_ROC%d%s_%s_%s_%s.csv', ri, rig{ri}, strrep(IQRfield,'res_','') , qaversions{qai} , segment{si})), mytable.ROC{ri} ) ; 
    end
    cat_io_csv( fullfile( printdir , sprintf('MRART_ANOVA_%s_%s_%s.csv', strrep(IQRfield,'res_','') , qaversions{qai}, segment{si})) ,  mytable.ANOVA ); 
  
    fprintf('%s done.\n',segment{si})
  end
end
fprintf('all done.\n')
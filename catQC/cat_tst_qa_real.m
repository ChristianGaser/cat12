%function cat_tst_qa_real6
% qa_group_test
% Test the automatic QA for 2-3 groups with expert quality rating.
% Problematic is the difference between image quality (the acutal
% rating) and that the final QA measure was scaled to match kappa.
% Therefore, it requires a rating for the final preprocessing quality 
% of each method and not only of the original image quality. 
% As compro1se noise rating can be used ...
% 
% Another problem is the scaling question of the marks ...

% TODO2022: 
%  * add machine learning outlier detection 
%  * add BIDS support for QC ?
%  * test the input variability of the methods? n-fold?
%
% Details 
%  * use two test scenarious (remove light, remove severe artefacts)?
%    - not fitting for all ratings! but also no useful effect
%
% Processing / Rating: 
%  * Rating CAMCAN
%

%#ok<*SAGROW,*AGROW>



% -- group test skript ----------------------------------------------
%  P          .. file-structure
%   .p0       .. segmentation files
%   .rating   .. expert ratings [0 - none, 0.25 ~ very very light; 0.5 ~ light; 0.75 ~ ; 1 - obvious; 1.5 - severe]
%   .csv      .. csv files with expert ratings and p0-filenames
%   .mdirs    .. ?
%   .xml      .. catxml files
%   .p0e      .. existing p0-files
%  X          .. XML-data structure
%  Q          .. XML-data (qc, global values)



  % the qafile is used to test different versions of the QC estimation [�IMPORTANT ]
  rerun     = 2; 
  fasttest  = 0; fast = {'full7','fast7'};
  qaversions = {
    'cat_vol_qa202412';
   ... 'cat_vol_qa202310'; 
   ... 'cat_vol_qa201901x';
  };
  qais  = 1:numel(qaversions);
  quai  = 1; 
%
segment = {'qcseg'}; si=1;
  for quai = qais
    qafile = qaversions{quai}; 

  
    % default parameter
    if ~exist('opt','var'), opt = struct(); end
    % printing options
    def.type      = '-dpng';  % dataype -depsc
    def.res       = '-r300';  % resolution
    def.dpi       = 72; 
    def.closefig  = 1;        % close figure after printing
    opt.maindir   = '/Volumes/SG5TB/MRData/202105_QA';          % main project directory [�IMPORTANT ]
    opt.resdir    = fullfile(opt.maindir,'+results',['private_' fast{fasttest+1}]);   % directory for results 
    opt.printdir  = fullfile(opt.maindir,'+figure' ,['private_' fast{fasttest+1}]);   % directory for figures
    opt           = cat_io_checkinopt(opt,def);
  
    % translate marks into percentage rating
    mark2rps    = @(mark) min(100,max(0,105 - mark*10));
    
    % create result directories
    if ~exist(opt.resdir,'dir'),   mkdir(opt.resdir); end
    if ~exist(opt.printdir,'dir'), mkdir(opt.printdir); end
    
    
  % ############
  % * add site variable or even sex and age
  % ############
    % files with ratings and additional variables 
    opt.rateCSV   = { % NIFTIs CSV column
      fullfile(opt.maindir,'ATLAS','ATLASR1.1a'),                   fullfile(opt.maindir,'+ratings','QC2022_Atlas.csv'),      'ATLASm',   5; % 2x avg
      fullfile(opt.maindir,'ATLAS','ATLASR1.1b'),                   fullfile(opt.maindir,'+ratings','QC2022_Atlas.csv'),      'ATLASo',   5; % 2x avg
      fullfile(opt.maindir,'private','Site J - Calgary preschool'), fullfile(opt.maindir,'+ratings','QC2022_Calgary.csv'),    'Calgary',  9; % 3x avg 
      fullfile(opt.maindir,'private','Site M - BEANstudy'),         fullfile(opt.maindir,'+ratings','QC2022_BEAN.csv'),       'BEAN',     5; % 1x 
      fullfile(opt.maindir,'private','Site O - ABIDE'),             fullfile(opt.maindir,'+ratings','QC2022_ABIDE2.csv'),     'ABIDE',    5; % 1x + checked
      fullfile(opt.maindir,'private','Site Q - NIH-HBD'),           fullfile(opt.maindir,'+ratings','QC2022_NIH.csv'),        'NIH',      5; % 1x  
      fullfile(opt.maindir,'IXI'),                                  fullfile(opt.maindir,'+ratings','QC2022_IXI.csv'),        'IXI',      5; % 1x  
      fullfile(opt.maindir,'private','fHKS'),                       fullfile(opt.maindir,'+ratings','QC2022_fHKS.csv'),       'fHKS',    5; % 1x avg
      ... PPMI ?    >   publication rules issue
      ... INDI ?    >   yes, it support a lot of cites
      ... NKI neu?? 
      ...   ds0000144 anxiety ?
      ...   ds0030097 AOMIC 1000 ? 
      ...   ds#       QTIM ? > yes, because I can use this later 
      ...   ds#       QTAB ? > yes, because I can use this later 
      };
    
    %  -- find p0-files --------------------------------------------------
    stime = cat_io_cmd('Find segmentation files:','','',1); fprintf('\n'); 
    
    
    % (1) QC ratings by groups as external ratings from function
    % ######## would be good to transfer this into the google table #########
    [good,bad]        = get_rating; 
    get_rating_sites  = { % fieldname , directory
     'NKI', 'Site G - NKI';
     'NYC', 'Site H - ADHD_NYC';
     'ORE', 'Site I - ADHD_ORE';
    };
    
    
    % (2) QC ratings from group-wise CSV files to good, to bad, and to P.p0{1} list
    P.p0{1}  = cell(0); stime2 = [];  % used different methods
    P.rating = zeros(0); 
    % search CSV files
    if size(opt.rateCSV,1)>0, cat_io_cmd('  CSV files:','g8','',1); stime=[]; fprintf('\n'); end
    for fi = 1:size(opt.rateCSV,1)
      stime = cat_io_cmd(sprintf('    %s:',opt.rateCSV{fi,3}),'g5','',1,stime); stime2 = stime; 
      
      % read table and remove empty name or rating fields 
      csv         = cat_io_csv(opt.rateCSV{fi,2}); 
      Pp0ff       = csv(2:end,1); 
      Pp0ffrating = csv(2:end,opt.rateCSV{fi,4});
      Pp0ffrating(cellfun('isempty',Pp0ff)) = [];  Pp0ff(cellfun('isempty',Pp0ff)) = [];   
      Pp0ff(cellfun('isempty',Pp0ffrating)) = [];  Pp0ffrating(cellfun('isempty',Pp0ffrating)) = [];   
      Pp0ffrating = cell2mat(Pp0ffrating); 
      P.csv.(opt.rateCSV{fi,3}) = csv; 
      
      % remove previous entries of the study on both list   
      if isfield( bad , opt.rateCSV{fi,3} ) 
        bad = rmfield( bad , opt.rateCSV{fi,3}); 
      else
        bad.(opt.rateCSV{fi,3}) = {}; 
      end
      if isfield( good , opt.rateCSV{fi,3} ) 
        good = rmfield( good , opt.rateCSV{fi,3}); 
      else
        good.(opt.rateCSV{fi,3}) = {}; 
      end
      
      % add the segmentation file and the ratings based on the csv list
      for pi = 1:size(Pp0ff,1) 
        % first get the file fromd different possible dirs .... 
        [pp,ff,ee] = spm_fileparts(Pp0ff{pi}); 
        if strcmp(ee,'.gz'), exnii = ''; else, exnii = '.nii'; end % catch .nii.gz cases
        
        % first try to find the default mri dir
        if exist(fullfile(opt.rateCSV{fi,1},pp,'mri'),'dir')
          file = fullfile( opt.rateCSV{fi,1} , pp , 'mri' , ['p0' ff exnii] ); 
        else % otherwise try the maindir
          file = fullfile( opt.rateCSV{fi,1} , pp , ['p0' ff exnii] );
        end
        
        % if nothing was found try to search for the files 
        if ~exist(file,'file')
          % search for the file (this is slow so we try to avoid it before)
          if exist(fullfile(opt.rateCSV{fi,1},pp,'mri'),'dir')
            file = cat_vol_findfiles(fullfile(opt.rateCSV{fi,1},pp,'mri'), ['p0*' ff exnii],struct('chararr',1));
          else
            file = cat_vol_findfiles(fullfile(opt.rateCSV{fi,1},pp), ['p0*' ff exnii],struct('chararr',1));
          end
          
          if ~exist(file,'file')
           fprintf('Miss: "%s"\n',Pp0ff{pi}); 
           file = ''; 
          end
        end
        
        % add the files to the list depending on the rating
        if ~isempty(file) 
          P.p0{1} = [P.p0{1}; file ];
          if Pp0ffrating(pi) > 0.5
            bad.( opt.rateCSV{fi,3} )  = [ bad.( opt.rateCSV{fi,3} );  file ]; 
          else
            good.( opt.rateCSV{fi,3} ) = [ good.( opt.rateCSV{fi,3} ); file ]; 
          end
          
          P.rating(end+1) = Pp0ffrating(pi); 
        end
      end
    end
    
    
    % add files from good/bad list
    catDBsites = setdiff( fieldnames(good) , opt.rateCSV(:,3) ) ;
    for sitei = 1:numel(catDBsites)
      sname = get_rating_sites{ strmatch(catDBsites{sitei},get_rating_sites(:,1)) , 1}; 
      sdir  = get_rating_sites{ strmatch(catDBsites{sitei},get_rating_sites(:,1)) , 2}; 
      stime2 = cat_io_cmd(sprintf('    %s:',sname),'g5','',1,stime2); 
  
      % add good images to file list
      for fi = 1:numel(good.(catDBsites{sitei}))
        [~,ff]  = fileparts(good.(catDBsites{sitei}){fi}); 
        file    = fullfile(fullfile(opt.maindir,'private'),sdir, 'mri',sprintf('p0%s.nii',ff));
        if ~exist(file,'file')
          file  = fullfile(fullfile(opt.maindir,'private'),sdir,sprintf('p0%s.nii',ff));
        end
        if ~exist(file,'file')
          file  = cat_vol_findfiles(fullfile(opt.maindir,'private',sdir),sprintf('p0%s*.nii',ff));
        end
        if ~isempty(file)
          P.p0{1} = [P.p0{1};file];
          P.rating(end+1) = 0; 
        end
      end
      
      % add bad images to file list
      for fi = 1:numel(bad.(catDBsites{sitei}))
        [~,ff]  = fileparts(bad.(catDBsites{sitei}){fi}); 
        file    = fullfile(fullfile(opt.maindir,'private'),sdir, 'mri',sprintf('p0%s.nii',ff));
        if ~exist(file,'file')
          file  = fullfile(fullfile(opt.maindir,'private'),sdir,sprintf('p0%s.nii',ff));
        end
        if ~exist(file,'file')
          file  = cat_vol_findfiles(fullfile(opt.maindir,'private',sdir),sprintf('p0%s*.nii',ff));
        end
        if ~isempty(file)
          P.p0{1} = [P.p0{1};file];
          P.rating(end+1) = 1; 
        end
      end
    end
    
    
    
    % search path of the QA-xml-files with good and bad
    stime2 = cat_io_cmd('  QC dirs:','g8','',1,stime2); 
    P.mdirs{1} = {
      'CAT12' fullfile(opt.maindir,'private','QA_good')  '<' [1.0 0.2 0.0] '-'; ...
      };
    P.mdirs{2} =  P.mdirs{1};
    for fi = 1:size(P.mdirs{1},1), P.mdirs{2}{fi,2} = strrep(P.mdirs{1}{fi,2},'QA_good','QA_bad'); end
    for ei = 1:2
      files           = cat_vol_findfiles(P.mdirs{ei}{1,2},'p0*.nii'); 
      P.p0{1}         = [P.p0{1}; files];
      P.rating(end+1:end+numel(files)) = (ei==2); % * ones(1,numel(files)); 
    end
    
    

    %% fasttest
    % ------------------------------------------------------------------------

    % get site directories
    if fasttest
      nlim        = 20;
      stime2 = cat_io_cmd(sprintf('  Limiting dataset to %d scans:',nlim),'g8','',1,stime2); 
      P.dirs      = spm_str_manip(P.p0{1},'hh'); 
      ses         = contains(P.dirs,'anat'); % get BIDS dirs
      P.dirs(ses) = spm_str_manip(P.dirs(ses),'hhh'); % rm sub BIDS dirs 
      P.udirs     = unique(P.dirs); % get sites
      for pi = 1:numel(P.udirs)
        Pdirs = contains(P.dirs,P.udirs{pi}); 
        if sum(Pdirs) > nlim 
          rm = Pdirs & cumsum(Pdirs)>nlim; 
          P.p0{1}(rm)   = [];
          P.dirs(rm)    = [];
          P.rating(rm)  = [];
        end
      end
    end

  
    %% estimate the qa and kappa values for the p0-files
    stime2 = cat_io_cmd('  XML files:','g8','',1,stime2); 
    if isfield(P,'xml'), P = rmfield(P,'xml'); end
    if isfield(P,'p0e'), P = rmfield(P,'p0e'); end
    % -- do QA ---------------------------------------------------------
    % xml-filenames 
    for pi = 1:numel(P.p0{1})
      [pp,ff] = fileparts(P.p0{1}{pi});
      [~,ppmri] = fileparts(pp); 
      if strcmp(ppmri,'mri')
        P.xml{1}{pi,1} = fullfile(strrep(pp,'mri','report'),[qafile '_' ff(3:end) '.xml']);
      else
        P.xml{1}{pi,1} = fullfile(pp,'report',[qafile '_' ff(3:end) '.xml']);
      end  
      P.p0e{1}(pi,1) = exist(P.xml{1}{pi},'file')>0;
    end
    


    
    
   %% (re)calcqa and load xml files
   % ------------------------------------------------------------------------
    stime2 = cat_io_cmd(sprintf('  Process QC of %d files:',sum(P.p0e{1})),'g8','',1,stime2); 
    P.qamat{1} = fullfile(opt.resdir,['group_' qafile '_NIR' P.mdirs{ei}{1,1} '.mat']);
    % run quality estimation 
    switch segment{si}
      case 'CAT'
        PX = P.p0{1}; 
        qaopt = struct('prefix',[qafile '_'],'write_csv',0,'mprefix','m','orgval',0,'rerun',rerun,'verb',2,'version',qafile);  % ##########################################
      case 'qcseg'
        PX = cat_io_strrep(P.p0{1},'/mri/p0','/'); 
        qaopt = struct('prefix',[qafile '_qcseg_'],'write_csv',0,'mprefix','m','orgval',0,'rerun',rerun,'verb',2,'version',qafile);  % ##########################################
    end
    if 1 %sum(P.p0e{1}==0)>0
      %eval(sprintf('%s(''p0'',P.p0{1}(P.p0e{1}==0),qaopt);',qafile)); 
     % cat_vol_qa('p0',P.p0{1}(P.p0e{1}==0),qaopt);
      cat_vol_qa('p0',PX(497:end),qaopt);
   %   cat_vol_qa('p0',P.p0{1},qaopt);
    end
    % if some files cannot process at all then have to remove them 
    for pi = sort( find( P.p0e{1}==0 ) , 'descend' )'
      if exist(P.xml{1}{pi},'file')==0
        cat_io_cprintf('err',sprintf('Error %s - missing measures. Remove entry\n', P.p0{1}{pi} ));
        P.p0{1}(pi)   = [];
        P.p0e{1}(pi)  = []; 
        P.xml{1}(pi)  = []; 
        P.rating(pi)  = []; 
      end
    end
    % load quality measures
    X.xml = cat_io_xml(P.xml{1});
    if 1
      %% bug correction
      reprocess = false(size(P.p0{1}));
      for xi=1:numel(X.xml)
        if isempty( X.xml(xi).filedata )
          reprocess(xi) = true; 
        end
      end
      if sum(reprocess)>0
        cat_vol_qa('p0',P.p0{1}(reprocess),qaopt);
        X.xml = cat_io_xml(P.xml{1});
      end
    end
    
    
    
    
  
    %%  -- prepare QA data ------------------------------------------------
    %   re-estimate new rating QR based on the existing measures QM
    Q = struct();
    stime = cat_io_cmd('  Prepare quality measurements:','g5','',1,stime2); 
    if isfield(P,'xml2'), P = rmfield(P,'xml2'); end
    if 1
      for  pi = 1:numel(P.xml{1})
        [pp,ff]       = spm_fileparts(P.p0{1}{pi}); 
        Q.fname{pi,1} = P.p0{1}{pi};
        Q.fpp{pi,1}   = pp;
        Q.fff{pi,1}   = ff;
        if 0
          try
            %  -- remarks ----------------------------------------------------- 
            X.xml2(pi) = cat_stat_marks('eval',0,X.xml(pi)); %,default); %,P.mdirs{1,1});
          catch 
            cat_io_cprintf('err','remarking failed for scan %d: %s',pi); 
          end
        end
      end
    end
    X.xml2 =  X.xml; 
    
    
    
    %  -- map QM data ----------------------------------------------------
    %  map values to the evaluation structure Q
    fieldxml = 'xml2'; % use original (xml) or updated (xml2) rating
    QM = {'NCR','ICR','res_RMS','res_ECR','res_ECRmm','FEC','contrastr','IQR','SIQR'};
    for q1 = 1:numel(QM)
      field   = QM{q1};
      fieldqm = 'qualityratings';
      for pi=1:numel(P.xml{1})
        try
          Q.(field)(pi,1)    = X.(fieldxml)(pi).(fieldqm).(QM{q1});
        catch
          Q.(field)(pi,1)    = nan; 
        end
      end
    end
    for pi=1:numel(P.xml{1})
      try
        Q.vx_vol(pi,:)  =  X.xml(pi).qualitymeasures.res_vx_vol; 
      catch
        Q.vx_vol        = nan(1,3);
      end
    end
    
    % remove low resolution entries
    if 1
      LR = any(Q.vx_vol>2,2); 
      
      P.p0{1}(LR)  = []; 
      P.xml{1}(LR) = []; 
      P.p0e{1}(LR) = [];
      P.rating(LR) = []; 
      
      X.xml(LR)    = [];
      X.xml2(LR)   = [];
      
      FN = fieldnames(Q); 
      for fni = 1:numel(FN)
        Q.(FN{fni})(LR) = [];
      end
    end

    % remove unclear entries
    if 0
      LR = P.rating>0.2 & P.rating<0.9; 
      
      P.p0{1}(LR)  = []; 
      P.xml{1}(LR) = []; 
      P.p0e{1}(LR) = [];
      P.rating(LR) = []; 
      
      X.xml(LR)    = [];
      X.xml2(LR)   = [];
      
      FN = fieldnames(Q); 
      for fni = 1:numel(FN)
        Q.(FN{fni})(LR) = [];
      end
    end
    
    
    
    % site alignment
    % -----------------------------------------------------------------------
    % Seperation by part of the path with a short name (siten) and an id (site)
    % -----------------------------------------------------------------------
    for pi = numel(X.xml):-1:1
      path = X.xml(pi).filedata.path; hhd = '';  
      for tpi = 1:6
        hhdo = hhd; 
        try
          [path,hhd,eed] = fileparts( path ); hhd = [hhd eed];
        catch
          hhd = 'nan'; 
        end
        switch hhd
          case 'Site A - sCMV',               Q.siten{pi} = 'CMV';    Q.site(pi) = 1; % ok
          case 'Site B - FHD-Boston',         Q.siten{pi} = 'FHD';    Q.site(pi) = 2; % ok
          case 'Site C - ADHD200_NYC',        Q.siten{pi} = 'AD2k';   Q.site(pi) = 3; % ok
          case 'Site D - fHKS',               Q.siten{pi} = 'fHKS';   Q.site(pi) = 4; % ok 
          case 'Site E - 99s_113c',           Q.siten{pi} = '99s';    Q.site(pi) = 5; % ok
          case 'Site F - FHD-Boston_read',    Q.siten{pi} = 'FHDB';   Q.site(pi) = 6; % ok
          case 'Site G - NKI',                Q.siten{pi} = 'NKI';    Q.site(pi) = 7; % ok
          case 'Site H - ADHD_NYC',           Q.siten{pi} = 'ADNYC';  Q.site(pi) = 8; % ok
          case 'Site I - ADHD_ORE',           Q.siten{pi} = 'ADORE';  Q.site(pi) = 9; % ok
          case 'IXI'  
            site        = double(P.csv.IXI{find(cellfun('isempty',strfind( P.csv.IXI(:,1) , X.xml(pi).filedata.file))==0,1,'first'),8}); 
            Q.siten{pi} = sprintf('IXI%d',site);    
            Q.site(pi)  = 9 + site; % 
          ... Calgary preschool: 
          ... - the protocoll is overall problematic and even the best images are extremly noisy
          ... - see also ONRC_2_part# (that uses a rescan concept for averaging)
          ... - you can try to keep only the best and worst?  |? remove only the worst (higher rating treshold?) 
          case 'Site J - Calgary preschool';  Q.siten{pi} = 'Calga';  Q.site(pi) = 13; % removed due to general low procoll qualy (but **highres**)
          case 'Site M - BEANstudy';          Q.siten{pi} = 'BEAN';  Q.site(pi) = 14; % R1 + checked, ok
  ... case 'Site N - twins';              Q.siten{pi} = 'Twins';  Q.site(pi) = 15; % need rating & need procesing
          case 'ABIDEII-BNI_1';               Q.siten{pi} = 'ABNI'; Q.site(pi) = 16; % R1 + checked, ok
          case 'ABIDEII-EMC_1';               Q.siten{pi} = 'AEMC'; Q.site(pi) = 17; % R1 + checked, ok
          case 'ABIDEII-ETH_1';               Q.siten{pi} = 'AETH'; Q.site(pi) = 18; % R1 + checked, ok
          case 'ABIDEII-GU_1';                Q.siten{pi} = 'AGU';  Q.site(pi) = 19; % R1 + checked, ok
          case 'ABIDEII-IP_1';                Q.siten{pi} = 'AIP';  Q.site(pi) = 20; % R1 + checked
          ... - the protocol is quite noisy (1.5 T?) but ok  
          case 'ABIDEII-IU_1';                Q.siten{pi} = 'AIU';  Q.site(pi) = 21; % R1 + checked, ok - highres
          case 'ABIDEII-KKI_1_29273_29322';   Q.siten{pi} = 'AKKI'; Q.site(pi) = 22; % R1 + checked
          case 'ABIDEII-KKI_1_29323_29372';   Q.siten{pi} = 'AKKI'; Q.site(pi) = 22; % R1 + checked
          case 'ABIDEII-KKI_1_29373_29423';   Q.siten{pi} = 'AKKI'; Q.site(pi) = 22; % R1 + checked
          case 'ABIDEII-KKI_1_29424_29485';   Q.siten{pi} = 'AKKI'; Q.site(pi) = 22; % R1 + checked
          case 'ABIDEII-KUL_3';               Q.siten{pi} = 'AKUL'; Q.site(pi) = 23; % R1 + chekced, ok, Nf=1
          case 'ABIDEII-NYU_1';               Q.siten{pi} = 'ANYU'; Q.site(pi) = 24; % R1 + checked, ok
          case 'ABIDEII-NYU_2';               Q.siten{pi} = 'ANYU'; Q.site(pi) = 24; % R1 + checked, ok
          case 'ABIDEII-OHSU_1';              Q.siten{pi} = 'AOHS'; Q.site(pi) = 25; % R1 + checked, ok, various resolutions
          ... ABIDE ORNC removed because: 
          ... - a fast/noisy protocol with many rescans was used 
          ... - in such noise data the variation by MAs is in range of "normal" SNR 
          ... - image averaging is an essential step here and the question could be if 
          ...   rescans with severe MAs could be detected before averaging 
          ...   well, I would argue that this should be part of the averaging 
          ...   routine that support a much more specific test-retest
  case 'ABIDEII-ONRC_2_part1';        Q.siten{pi} = 'AORN'; Q.site(pi) = 26; % R1 - **highres**, highnoise
  case 'ABIDEII-ONRC_2_part2';        Q.siten{pi} = 'AORN'; Q.site(pi) = 26; % R1
  case 'ABIDEII-ONRC_2_part3';        Q.siten{pi} = 'AORN'; Q.site(pi) = 26; % R1
  case 'ABIDEII-ONRC_2_part4';        Q.siten{pi} = 'AORN'; Q.site(pi) = 26; % R1
          ...
          case 'ABIDEII-SDSU_1';              Q.siten{pi} = 'ASDS'; Q.site(pi) = 27; % R1 + checked2, ok
          case 'ABIDEII-STANFORD';            Q.siten{pi} = 'ASTA'; Q.site(pi) = 28; % R1 + checked, ~ok
          case 'ABIDEII-TCD_1';               Q.siten{pi} = 'ATCD'; Q.site(pi) = 29; % R1 + checked2, ~ok
          case 'ABIDEII-UCD_1';               Q.siten{pi} = 'AUCD'; Q.site(pi) = 30; % R1 + checked2, ok, Nf<10 
          case 'ABIDEII-UCLA_1';              Q.siten{pi} = 'AUCL'; Q.site(pi) = 31; % R1 + checked,  ok 
          case 'ABIDEII-UCLA_Long';           Q.siten{pi} = 'AUCL'; Q.site(pi) = 31; % R1 + checked,  ok .. longitudinal
          case 'ABIDEII-UM';                  Q.siten{pi} = 'AUM'; Q.site(pi) = 32; % R1 + checked,  ok
          case 'ABIDEII-USM_1';               Q.siten{pi} = 'AUSM'; Q.site(pi) = 33; % R1 + checked,  ok, Nf=4, **highres**
          case 'Site Q - NIH-HBD'   
            try
              site        = P.csv.NIH{find(cellfun('isempty',strfind( P.csv.NIH(:,1) , X.xml(pi).filedata.file))==0,1,'first'),8}; 
              Q.siten{pi} = sprintf('DEV%d',site);    
              Q.site(pi)  = 33 + site; % R1
            catch
              Q.siten{pi} = sprintf('DEV1');
              Q.site(pi)  = 33 + 1;
            end
         case 'ATLASR1.1a'
           site        = str2double(hhdo(2:end)); 
           Q.siten{pi} = sprintf('ATM%d',site);    
           Q.site(pi)  = 40 + site; 
         case 'ATLASR1.1b'
           site        = str2double(hhdo(2:end)); 
           Q.siten{pi} = sprintf('ATO%d',site);    
           Q.site(pi)  = 51 + site; 
          ... case 'Site O - PPMI';               Q.siten{pi} = 'P';  Q.site(pi) = 34; % 
         % case 'IXI'
          otherwise
            if tpi == 5 && isempty(Q.siten{pi})
              %%
              fprintf('"%s" - %s\n',hhd,P.p0{1}{pi});
            end
            
            % #################
            % remove entries ? 
            % #################
        end
      end
    end
    for ti = unique(Q.site), Q.N(Q.site==ti) = sum(Q.site==ti); end
    %%
    
    
    
    if 0
    % -----------------------------------------------------------------------
    % This part is for manual setups of the variables and ratings. E.g., to 
    % check the expert ratings (that where done on a xyz-slice preview) for 
    % deeper problems and miss-alignment (i.e., wrong fields).
    % -----------------------------------------------------------------------
      
      
      % check resolution and other variables 
      fprintf('Site:      %4s %4s  %4s %4s %4s   %6s  >> %6s\n','pass','fail','x','y','z','min','RES'); 
      for si = unique(Q.site)
        xi = find( Q.site == si , 1, 'first'); 
        fprintf('Site %2d:   %4d %4d  %4.2fx%4.2fx%4.2f   %6.2f  >> %6.2f\n',...
          si, numel( Q.site == si & Q.group == 0), numel( Q.site == si & Q.group == 0), Q.vx_vol( xi , :), min(Q.vx_vol( xi , :)) , Q.res_RMS(xi) ); 
      end
      
      
      %% check for outlieres of the fast expert ratings
      %  --------------------------------------------------------------------
      %  I.e., passed images with very low IQR values and failed images with
      %  very high IQR values. In the fast rating procedure (based on a xyz-
      %  slice preview) local motion artifacts can be over-seen (often) or 
      %  over-interpretaded (rare).  
      %  --------------------------------------------------------------------
      clc
      tsite   = 42;  % select a site (see site definition above)
      tsort   = 0;  % sort the output with the most relant cases on top)
      % threshold between groups as the mean of the median of the lower
      % (passed) and upper (failed) median (the sides of the boxplot)
      trating = mean( [ median( Q.SIQR ( Q.site' == tsite & Q.group == 0 )) , ...
                        median( Q.SIQR ( Q.site' == tsite & Q.group == 1 )) ] ); 
      tgroup{1}   = Q.SIQR > trating & Q.site' == tsite & Q.group == 0 & Q.res_RMS < 3; % passed but low IQR
      tgroup{2}   = Q.SIQR < trating & Q.site' == tsite & Q.group == 1 & Q.res_RMS < 3; % failed but high IQR 
      texgroup{1} = P.p0{1}(tgroup{1});  
      texgroup{2} = P.p0{1}(tgroup{2}); 
      texgval{1}  = Q.SIQR(tgroup{1});  
      texgval{2}  = Q.SIQR(tgroup{2}); 
      id{1}       = find( tgroup{1} ); 
      id{2}       = find( tgroup{2} ); 
      tname       = {'Expert passed but low IQR (over-seen MAs)'; 
                     'Expert failed but relative good IQR (over-interpretated MAs)'};
      tcolor      = {[0 0.5 0]; [0.5 0 0] }; 
      torder      = {'descend'; 'ascend'}; 
      fprintf('Site %d "%s" (theshold=%0.2f):\n',tsite, Q.siten{find(Q.site==tsite,1,'first')},trating); 
      for gi = 1:2
        fprintf('%s:\n',tname{gi}); 
        if tsort, [~,sid] = sort(texgval{gi},torder{gi}); else, sid = 1:numel(texgval{gi}); end
        for fi = 1:numel(texgroup{gi})
          cat_io_cprintf(tcolor{gi},sprintf('  %4d - %0.2f - %s\n', id{gi}(sid(fi)), ...
            texgval{gi}(sid(fi)), spm_str_manip( texgroup{gi}{sid(fi)} ,'a70'))); 
        end
      end
      % general histogram plot of IQR values
      figure(933); clf(933);
      M       = Q.methods(:,1:size(Q.group,2))>0 & Q.group==0 & Q.site' == tsite; 
      [hg,hr] = hist(Q.(IQRfield)(M(:)),0.5:0.05:6.5);
      M       = Q.methods(:,1:size(Q.group,2))>0 & Q.group==1 & Q.site' == tsite; 
      [hb,hr] = hist(Q.(IQRfield)(M(:)),0.5:0.05:6.5); hold on
      fill([hr,hr(end)],[hb 0],[0.9 0.0 0.2],'facealpha',0.2); 
      fill([hr,hr(end)],[hg 0],[0.0 0.8 0.2],'facealpha',0.2);  
      plot(hr,hg,'color',[0.0 0.8 0.2],'Linestyle','-','linewidth',1.5); 
      plot(hr,hb,'color',[0.9 0.0 0.2],'Linestyle','-','linewidth',1.5); 
      xlim([0.5 6.5]); ylim([0 max([hg,hb])*1.2]); hold off; box on; grid on; 
      title('IQR histogram of expert rated scans','FontSize',FS(1),'FontWeight','bold'); 
      xlabel('IQR','FontSize',FS(2)); ylabel('number of scans','FontSize',FS(2)); 
      legend({'passed','failed'},'Location','northeast'); legend('boxoff'); 
      
      
      %% independent QC processing test (different xml-files!)
      %  ATLAS 346, 462, 436, 486 with incorrect segmentation (incorrect setup) 
      rerun = 2;     % FEC: 2,25,27,21 || ECR: 27,24, 20 ||�
      tside      = unique(Q.site); %1; % 1, 3, 5, 6, ... 18 23 30! 32 33 ... 41 42! 44 45 46  48 47 49 ... +11 
      slim       = 5; 
      qaversions = {'cat_vol_qa201901';  'cat_vol_qa202110'; 'cat_vol_qa'; 'cat_vol_qa202207b';'cat_vol_qa202210'}; 
      qaversions = {'cat_vol_qa202310'}; 
      qafield    = {'NCR','ICR','res_RMS','res_ECR','FEC','IQR', 'SIQR',}; clear testIQR; 
      for tsi = tside
        for qai = numel(qaversions)
          qaopt2 = struct('prefix','test_','write_csv',0,'mprefix','m','orgval',0,'rerun',rerun,'verb',2,'version',qaversions{qai});
          sideidg = find(Q.site == tside(tsi) & P.rating<0.5); if numel(sideidg)>slim, sideidg = sideidg(1:round(numel(sideidg)/slim):end); end %sideidg = []; 
          sideidb = find(Q.site == tside(tsi) & P.rating>=.5); if numel(sideidb)>slim, sideidb = sideidb(1:round(numel(sideidb)/slim):end); end
          sideidx = [sideidg , sideidb];
          if isempty(sideidx), continue; end
          fprintf('Run site %d - %s:\n',tside(tsi),Q.siten{sideidx(1)}); 
          qar = cat_vol_qa('p0',P.p0{1}( sideidx ),qaopt2); 
          %qar = cat_io_xml( P.xml{1}(sideidx) ); 
          for qaii = 1:numel(qafield), for qi = 1:numel(qar), testIQR{qaii}(qi) = qar(qi).qualityratings.(qafield{qaii}); end; end
          fhx=figure(939); fhx.Position(3:4) = [1000 200]; fhx.Name = sprintf('(%d) %s - %s',tside(tsi),Q.siten{sideidx(1)},qaversions{qai});
          for qaii = 1:numel(qafield)
            %nr = cat_tst_qa_normer( testIQR{qaii}' , struct('figure',0,'model',0,'cmodel',2));
            fhx=figure(939); subplot(1,numel(qafield),qaii); 
            cat_plot_boxplot({ testIQR{qaii}(P.rating( sideidx )<0.5); testIQR{qaii}(P.rating( sideidx )>=0.5) }, ...
              struct('usescatter',1,'style',4,'datasymbol','o','ygrid',0,'colormap',[0 .5 0; .5 0 0])); ylim([0 10]); grid on
            title(sprintf('%s',strrep(qafield{qaii},'_','\_'))); 
            
          end
          if ~exist(fullfile(opt.printdir,'fstsideboxplot'),'dir'), mkdir(fullfile(opt.printdir,'fstsideboxplot')); end
          print(fhx,fullfile(opt.printdir,'fstsideboxplot',sprintf('fig_perExpertGroups_%d-%s_%s_%s',tside(tsi),Q.siten{sideidx(1)},qaversions{qai},datestr(clock,'YYYYmm'))),opt.res,opt.type);
        
        end
      end
      %% rerun one site to update the CQ after sever changes  
      qaopt2 = struct('prefix',[qafile '_'],'write_csv',0,'mprefix','m','orgval',0,'rerun',rerun,'verb',2);
      eval(sprintf('qar = %s(''p0'',  P.p0{1}( Q.site == 32 ) ,qaopt2);',qafile)); %
      % - update Q?
      % - tmp-xml-name = filename 
       
      %%
      cat_tst_qa_normer( Q.(IQRfield)(Q.site' == 21 & Q.train == 1) , struct('figure',3,'model',0,'cmodel',2));
      
      %%
      for si=unique(Q.site) 
        cat_tst_qa_normer( Q.(IQRfield)(Q.site' == si), struct('figure',3,'model',4,'cmodel',2)); 
        title(sprintf('Site %d (n=%d)',si,sum(Q.site' == si))); pause(0.5); 
      end
    
    end
    
    
    % setup group and add grouped data
    allgood = {}; allbad = {};
    goodfields = fieldnames(good);
    for fi=1:numel(goodfields), allgood = [allgood; good.(goodfields{fi})]; end
    for fi=1:numel(goodfields), allbad  = [allbad;  bad.(goodfields{fi})];  end
    for pi = 1:numel(P.xml{1})
      [~,ff]        = fileparts(fileparts(X(1).xml(pi).filedata.path)); 
      Q.group(pi,1) = (~contains(ff,'QA_good') && contains(ff,'QA_bad') ) || ... 
                       ( ~strcmp(X(1).xml(pi).filedata.file,'anat') && ... 
                          any(cellfun('isempty',strfind( allbad , X(1).xml(pi).filedata.file ))==0 )) || ...
                       (  strcmp(X(1).xml(pi).filedata.file,'anat') && ...
                           any(cellfun('isempty',strfind( allbad , fullfile( ... 
                       spm_str_manip( X(1).xml(pi).filedata.path , 'hhht'), ...
                       spm_str_manip( X(1).xml(pi).filedata.path , 'hht'), ...
                       spm_str_manip( X(1).xml(pi).filedata.path , 'ht') ...
                       ) ))==0 ));
      Q.ff{pi,1}    = ff; 
  
      try
        Q.CMV(pi,1)   =  X(1).xml(pi).subjectmeasures.vol_rel_CGW(1);
        Q.GMV(pi,1)   =  X(1).xml(pi).subjectmeasures.vol_rel_CGW(2);
        Q.WMV(pi,1)   =  X(1).xml(pi).subjectmeasures.vol_rel_CGW(3);
      catch
        Q.CMV(pi,1)   = nan;
        Q.GMV(pi,1)   = nan; 
        Q.WMV(pi,1)   = nan; 
      end
    end
    % Q.GMVC = Q.GMV - repmat(mean(Q.GMV,1),size(Q.GMV,1),1) + mean(Q.GMV(:));  % should this not be group-site?   
    
    
  
  
    %  -- prepare group data ---------------------------------------------
    stime  = cat_io_cmd('  Prepare group data:','g5','',1,stime); 
    groups = '2'; groupsn = {'P';'F'};
    gcolor = lines(max(Q.site)); gcolor = repmat(gcolor,1,str2double(groups)); 
    gcolor = reshape(gcolor',3,numel(gcolor)/3)'; 
    Q.methods = 1; % ############################ REMOVE ??? ###############
    Q.methods = repmat(Q.methods,size(Q.NCR,1),1);
    %Q.train   = rand(size(Q.CMV))<0.3; % use random value
    Q.train   = mod(1:numel(Q.CMV),3)' > 0;        % split half, third, quad, ...
    % create groups
    FPgroups = cell(1,str2double(groups));
    for gi=1:str2double(groups)
      FPgroups{gi} = Q.group==(gi-1) & Q.methods(:,1:max(1));
    end
    
    
    
    
    
    
    %% -- create figure --------------------------------------------------
    
    IQRfields = {'NCR','SIQR','IQR','ICR','res_ECR','FEC','res_RMS'}; 
    %IQRfields = {'NCR','SIQR','IQR','res_ECR','FEC'}; 
  
    for IQRfieldi = 1:numel(IQRfields)
    
      IQRfield  = IQRfields{IQRfieldi};
      norm      = {'' 'n'};
      Pfields   = {IQRfield,['n' IQRfield]};
      FS        = [8 8 2] * (opt.dpi/100) * 2; 
      rmse      = @(x) mean( x.^2 ).^0.5;
    
      if all(isnan(Q.(IQRfield))), continue; end 

    % ##########  subdata   = 1; % test just a subset, run only good? 
      model     = 0; 
      cmodel    = 1; 
      onlygood  = 0; % print only passed values on the site specific second figure
      % mod-cmod:
      %  01 vs 02: the median model is surprisingly good  (cmod=2 is ok) 
      %  11 vs 12: quite good but similar to 0#           (cmod=2 is also not working)  
      %  21 vs 22: many positve outliers                  (cmod=2 is not working)
      %  31 vs 32: ok
      %  41 vs 42: ok
      %  51 vs 52: ok (cmod=2 is not working)
      % 6????

      testsites = unique(Q.site);
      testsites = setdiff( unique(Q.site) , [13 25 26 28 37:60]); 
      %testsites = [13 26 28 21 42 44 45 47 55 56 30 29 2    27 53    58 17 20];
      %testsites = [1 2 3 4 5 6 7 8 9              18 19 20    23 24 25    31 32 33]; 
      %testsites = [1 2 3 4 5 6 7 8 9              18 19 20    23 24 25    31 32 33]; 
      for ti = numel(testsites):-1:1
        if sum(Q.site==testsites(ti))<10 || sum(Q.site'==testsites(ti) & Q.group==0)<5 || sum(Q.site'==testsites(ti) & Q.group==1)<5, testsites(ti) = []; end 
      end
      
      
      % prepare normalized values with full dataset as far as this is only for visualization 
      for si = unique(Q.site)
        Mt = Q.site' == si; 
        Q.(['n' IQRfield])(Mt) = cat_tst_qa_normer( Q.(IQRfield)(Mt), ...
          struct('model',model,'figure',0,'cmodel',cmodel,'train',4));
      end
      stime = cat_io_cmd(sprintf('  Create figure (%s,%s,model=%d,cmodel=%d):',qafile,IQRfield,model,cmodel),'b','',1,stime); 
      
      for nci = numel(Pfields)
        % general histogram plot of N(S)IQR values ( 2 groups )
        fi  = 934 + model + 7*(cmodel-1);
        fdn = sprintf('%dsites_model%d_cmodel%d',numel(testsites),model,cmodel);
        fd  = figure(fi); clf(fi); fd.Name = fdn; fd.Position(3) = 600; fd.Position(4) = 400; 
        
        subplot(2,2,1);
        Mp       = Q.methods(:,1:size(Q.group,2))>0 & Q.group==0 & ...
          any(repmat(Q.site',1,numel(testsites)) == testsites,2); 
        Mf       = Q.methods(:,1:size(Q.group,2))>0 & Q.group==1 & ...
          any(repmat(Q.site',1,numel(testsites)) == testsites,2); 
        if nci == 1
          [hp,hr]  = hist(105 - 10*Q.([norm{nci} IQRfield])(Mp(:)),0:100);
          [hf,hr]  = hist(105 - 10*Q.([norm{nci} IQRfield])(Mf(:)),0:100); hold on
          myxlim   = [40 100]; 
        else
          [hp,hr]  = hist(0 - 10 * Q.([norm{nci} IQRfield])(Mp(:)),-80:0.5:40);
          [hf,hr]  = hist(0 - 10 * Q.([norm{nci} IQRfield])(Mf(:)),-80:0.5:40); hold on
          myxlim   = [-60 10]; 
        end
         
        th = cat_stat_nanmean( [cat_stat_nanmedian( Q.([norm{nci} IQRfield])( Mf(:) & Q.([norm{nci} IQRfield])(:)<cat_stat_nanmedian(Q.([norm{nci}  IQRfield])(Mf(:))) ) ) ... 
                                cat_stat_nanmedian( Q.([norm{nci} IQRfield])( Mp(:) & Q.([norm{nci} IQRfield])(:)>cat_stat_nanmedian(Q.([norm{nci}  IQRfield])(Mp(:))) ) ) ] ); 
        if nci == 1, thp = 105 - 10 * th; else, thp = -10 * th; end
       
        plot([thp thp],[0 max([hf,hp])*1.2],'Color',[0 0.5 1],'linewidth',0.5); 
        fill([hr,hr(1)],[hp 0],[0.0 0.8 0.2],'facealpha',0.2);  
        fill([hr,hr(1)],[hf 0],[0.9 0.0 0.2],'facealpha',0.1); 
        plot(hr,hp,'color',[0.0 0.8 0.2],'Linestyle','-','linewidth',1.5); 
        plot(hr,hf,'color',[0.9 0.0 0.2],'Linestyle','-','linewidth',1.5); 
        xlim(myxlim); ylim([0 max([hp,hf])*1.2]); hold off; box on; grid on; ax=gca; ax.XDir = 'reverse'; 
        title(sprintf('%s histogram of expert rated scans',['n' IQRfield]),'FontSize',FS(1),'FontWeight','bold'); 
        xlabel(['n' IQRfield],'FontSize',FS(2)); ylabel('number of scans','FontSize',FS(2)); 
    
    % ######## add DICE?
    
        legend(...
          {sprintf('threshold (%0.2f)',thp), ...
           sprintf('passed (%0.2f�%0.2f)',-10*cat_stat_nanmedian(Q.([norm{nci} IQRfield])(Mp(:))), 10*cat_stat_nanstd(Q.([norm{nci} IQRfield])(Mp(:)))), ...
           sprintf('failed (%0.2f�%0.2f)',-10*cat_stat_nanmedian(Q.([norm{nci} IQRfield])(Mf(:))), 10*cat_stat_nanstd(Q.([norm{nci} IQRfield])(Mf(:))))},...
           'Location','northeast'); legend('boxoff'); 
        subplot(2,2,2);
        if nci == 1
          cdata = { 105 - 10 * Q.(IQRfield)(Mp(:)) , 105 - 10 * Q.(IQRfield)(Mf(:)) };
        else
          cdata = { 0 - 10 * Q.(['n' IQRfield])(Mp(:)) , 0 - 10 * Q.(['n' IQRfield])(Mf(:)) }; 
        end
        cat_plot_boxplot( cdata ,struct('names',{{'passed','failed'}} ,'usescatter',1,... 
              'style',4,'datasymbol','o','groupcolor',[0 0.8 0; 0.8 0 0], ...
              'sort',0,'ylim',myxlim,'groupnum',1,'ygrid',1)); set(gca,'FontSize',FS(2)*0.8)
        hold on; plot([-0.5 2.5],[thp thp],'Color',[0 0.5 1],'linewidth',0.5); 
        title(sprintf('%s boxplot of expert rated groups',[norm{nci} IQRfield]),'FontSize',FS(1),'FontWeight','bold'); 
        ylabel(['n' strrep(IQRfield,'_','\_')],'FontSize',FS(2)); xlabel('expert groups','FontSize',FS(2));
           
        
        fprintf('\n'); 
        fprintf('    Passed: % 5.2f �% 5.2f\n', -10*cat_stat_nanmedian(Q.([norm{nci} IQRfield])(Mp(:))), 10*cat_stat_nanstd(Q.([norm{nci} IQRfield])(Mp(:))) );
        fprintf('    Failed: % 5.2f �% 5.2f\n', -10*cat_stat_nanmedian(Q.([norm{nci} IQRfield])(Mf(:))), 10*cat_stat_nanstd(Q.([norm{nci} IQRfield])(Mf(:))) ); 
        fprintf('    Incorr: % 5.2f  % 5.2f%% (MF: %0.2f%%, MP: %0.2f%%)\n', thp , ...
          cat_stat_nansum(Q.([norm{nci} IQRfield])(Mf(:)) < th) / sum(Mf(:)) + sum(Q.([norm{nci} IQRfield])(Mp(:)) > th) / sum(Mp(:)) * 100, ...
          cat_stat_nansum(Q.([norm{nci} IQRfield])(Mf(:)) < th) / sum(Mf(:)) * 100, sum(Q.([norm{nci}  IQRfield])(Mp(:)) > th) / sum(Mp(:)) * 100); 
        
        
        % general histogram plot of N(S)IQR values ( 3 groups ) 
        subplot(2,2,3);
        Mp       = Q.methods(:,1:size(Q.group,2))>0 & P.rating'<0.5 & ...
          any(repmat(Q.site',1,numel(testsites)) == testsites,2); 
        Mq       = Q.methods(:,1:size(Q.group,2))>0 & P.rating'>=0.5 & P.rating'<1.0 & ...
          any(repmat(Q.site',1,numel(testsites)) == testsites,2); 
        Mf       = Q.methods(:,1:size(Q.group,2))>0 & P.rating'>=1.0 & ...
          any(repmat(Q.site',1,numel(testsites)) == testsites,2); 
        if nci == 1
          [hp,hr]  = hist(105 - 10 * Q.([norm{nci} IQRfield])(Mp(:)),0:100);
          [hq,hr]  = hist(105 - 10 * Q.([norm{nci} IQRfield])(Mq(:)),0:100); hold on
          [hf,hr]  = hist(105 - 10 * Q.([norm{nci} IQRfield])(Mf(:)),0:100); hold on
        else
          [hp,hr]  = hist(0 - 10 * Q.([norm{nci} IQRfield])(Mp(:)),-80:0.5:40);
          [hq,hr]  = hist(0 - 10 * Q.([norm{nci} IQRfield])(Mq(:)),-80:0.5:40); hold on
          [hf,hr]  = hist(0 - 10 * Q.([norm{nci} IQRfield])(Mf(:)),-80:0.5:40); hold on
        end
        
        plot([thp thp],[0 max([hf,hp])*1.2],'Color',[0 0.5 1],'linewidth',0.5); 
        fill([hr,hr(1)],[hp 0],[0.0 0.8 0.2],'facealpha',0.2);  
        fill([hr,hr(1)],[hq 0],[0.9 0.7 0.0],'facealpha',0.3); 
        fill([hr,hr(1)],[hf 0],[0.9 0.0 0.2],'facealpha',0.1); 
        plot(hr,hp,'color',[0.0 0.8 0.2],'Linestyle','-','linewidth',1.5); 
        plot(hr,hq,'color',[0.9 0.7 0.0],'Linestyle','-','linewidth',1.5); 
        plot(hr,hf,'color',[0.9 0.0 0.2],'Linestyle','-','linewidth',1.5); 
        xlim(myxlim); ylim([0 max([hp,hf])*1.2]); hold off; box on; grid on; ax=gca; ax.XDir = 'reverse'; 
        title(sprintf('%s histogram of expert rated scans',[norm{nci} IQRfield]),'FontSize',FS(1),'FontWeight','bold'); 
        xlabel(['n' IQRfield],'FontSize',FS(2)); ylabel('number of scans','FontSize',FS(2)); 
        legend( ...
          {sprintf('threshold (%0.2f)',thp), ...
           sprintf('passed (%0.2f�%0.2f)'      , 10*cat_stat_nanmedian(Q.([norm{nci} IQRfield])(Mp(:))), 10*cat_stat_nanstd(Q.([norm{nci} IQRfield])(Mp(:)))), ...
           sprintf('questionable (%0.2f�%0.2f)', 10*cat_stat_nanmedian(Q.([norm{nci} IQRfield])(Mq(:))), 10*cat_stat_nanstd(Q.([norm{nci} IQRfield])(Mq(:)))), ...
           sprintf('failed (%0.2f�%0.2f)'      , 10*cat_stat_nanmedian(Q.([norm{nci} IQRfield])(Mf(:))), 10*cat_stat_nanstd(Q.([norm{nci} IQRfield])(Mf(:))))}, ...
           'Location','northeast'); legend('boxoff');
        %
        subplot(2,2,4);
        if nci == 1
          pdata = { 105 - 10 * Q.(IQRfield)(Mp(:)) , 105 - 10 * Q.(IQRfield)(Mq(:)) , 105 - 10 * Q.(IQRfield)(Mf(:)) };
        else
          pdata = { 0 - 10 * Q.(['n' IQRfield])(Mp(:)) , 0 - 10 * Q.(['n' IQRfield])(Mq(:)) , 0 - 10 * Q.(['n' IQRfield])(Mf(:)) };
        end
        cat_plot_boxplot( pdata ,struct('names',{{'passed','question.','failed'}},'usescatter',1,...'boxwidth',-1,... 
              'style',4,'datasymbol','o','groupcolor',[0 0.8 0; 0.8 0.7 0; 0.8 0 0], ...
              'sort',0,'ylim',myxlim,'groupnum',1,'ygrid',1)); set(gca,'FontSize',FS(2)*0.8); 
        hold on; plot([-0.5 3.5],[thp thp],'Color',[0 0.5 1],'linewidth',0.5); 
        title(sprintf('%s boxplot of expert rated scans',[norm{nci} IQRfield]),'FontSize',FS(1),'FontWeight','bold'); 
        ylabel(['n' strrep(IQRfield,'_','\_')],'FontSize',FS(2)); xlabel('expert groups','FontSize',FS(2));
      
        if 0
          fprintf('\n'); 
          fprintf('    Passed: % 5.2f �% 5.2f\n', 10*median(Q.([norm{nci} IQRfield])(Mp(:))), 10*std(Q.([norm{nci} IQRfield])(Mp(:))));
          fprintf('    Questo: % 5.2f �% 5.2f\n', 10*median(Q.([norm{nci} IQRfield])(Mq(:))), 10*std(Q.([norm{nci} IQRfield])(Mq(:)))); 
          fprintf('    Failed: % 5.2f �% 5.2f\n', 10*median(Q.([norm{nci} IQRfield])(Mf(:))), 10*std(Q.([norm{nci} IQRfield])(Mf(:)))); 
        end  
        
        % save figure
        ndir = fullfile(opt.printdir,'boxplot'); if ~exist(ndir,'dir'); mkdir(ndir); end
        print(fd,fullfile(opt.printdir,'boxplot',sprintf('fig_%s_perExpertGroups_%s_%s_%s',IQRfield,fdn,qafile,datestr(clock,'YYYYmm'))),opt.res,opt.type);
        if opt.closefig, close(fd); end
      end
      
      
      %%
      % =======================================================================
      % Figure Quality Values for Expert Groups 
      % =======================================================================
      if 0 
        for ai = 1:2
          for pfi = 1:numel(Pfields)
            for tfi = 1:2 % sort
              % create/use figure
              if exist('fh','var') && numel(fh)>=tfi && fh(tfi).isvalid
                clf(fh(tfi)); figure(fh(tfi))
              else
                fh(tfi) = figure('Name',sprintf('figure %d - %s %s values per expert group (mod=%d,cmod=%d)',tfi,Pfields{pfi},model,cmodel),'Position',...
                  [0 0 2000 400],'color',[1 1 1],'PaperPositionMode','auto');
              end
    
              % extract IQR values
              MMCR = cell(''); groupsn2 = cell(''); n = str2double(groups); gxi=1;
              for si = testsites %unique(Q.site)% (~isinf(Q.site))) %[1,3:7,9:max(Q.site)]
                for gi=1:n
                  if pfi == 1
                    MMCR{gxi} = mark2rps( Q.(Pfields{pfi})( FPgroups{gi} & repmat(Q.site'==si,1,size(FPgroups{gi},2)) ) ); 
                  else
                    MMCR{gxi} = 0 - 10 * Q.(Pfields{pfi})( FPgroups{gi} & repmat(Q.site'==si,1,size(FPgroups{gi},2)) ); 
                  end
                  if ai == 1
                    groupsn2{gxi} = sprintf('%s%s(%02d)', lower(groupsn{gi}), Q.siten{find(Q.site==si,1,'first')}, Q.site(find(Q.site==si,1,'first')) );              % identifier 
                  else
                    groupsn2{gxi} = sprintf('%s%02d', lower(groupsn{gi}), Q.site(find(Q.site==si,1,'first')) );  % just a shortcut 
                  end
                  gxi = gxi+1;
                end
                % sorting scores
                MMCRbadn(gxi - 2)     = min([numel(MMCR{gxi - 2}),numel(MMCR{gxi - 1})]); 
                MMCRbadn(gxi - 1)     = MMCRbadn(gxi - 2); 
                MMCRgoodsep(gxi - 2)  = (( median(MMCR{gxi - 2}( MMCR{gxi - 2} < median(MMCR{gxi - 2})) )) - ... 
                                         ( median(MMCR{gxi - 1}( MMCR{gxi - 1} > median(MMCR{gxi - 1})) ))) / ...
                                           median(MMCR{gxi - 2}); 
                MMCRgoodsep(gxi - 1)  = MMCRgoodsep(gxi - 2) - eps; 
              end
              % sort by group
              if tfi == 2 % 
                [~,MMCRmeansorti] = sort(MMCRgoodsep,'descend'); 
              else
                MMCRmeansorti = 1:numel(MMCR); 
              end
              boxMMCR     = MMCR(MMCRmeansorti); 
              boxgroupsn2 = groupsn2(MMCRmeansorti);
              boxMMCRbadn = MMCRbadn(MMCRmeansorti);
              % remove center with low number of cases
              mincases = 5;
              if mincases>0
                boxMMCR     = boxMMCR(boxMMCRbadn>=mincases); 
                boxgroupsn2 = boxgroupsn2(boxMMCRbadn>=mincases); 
                boxgcolor   = gcolor(MMCRbadn>=mincases,:); 
              end
    
              if onlygood
                boxMMCR     = boxMMCR(1:2:end); 
                boxgroupsn2 = boxgroupsn2(1:2:end); 
                boxgcolor   = gcolor(1:2:end,:); 
              end
    
              % create boxplot
              if pfi == 1
                [bp1,bp2] = cat_plot_boxplot(boxMMCR,struct('names',{boxgroupsn2},'boxwidth',-1,... 
                  'style',4,'groupcolor',boxgcolor, 'usescatter',1,...
                  'sort',0,'ylim',[20 100],'groupnum',1,'ygrid',1)); set(gca,'FontSize',FS(2)*0.8)
              else
                [bp1,bp2] = cat_plot_boxplot(boxMMCR,struct('names',{boxgroupsn2},'boxwidth',1.5 * (2 - onlygood),... 
                  'style',4,'groupcolor',boxgcolor, 'usescatter',1,...
                  'sort',0,'ylim',[-60 + 30*onlygood 10],'groupnum',1,'ygrid',1)); set(gca,'FontSize',FS(2)*0.8)
              end
              set(gca,'xTickLabelRotation',90)
              set(gca,'Position',[0.05 0.25 0.94 0.70]); 
              title(sprintf('%s by expert quality group',Pfields{pfi},numel(1)),...
                'FontSize',FS(1),'FontWeight','bold'); 
              xlabel(sprintf('qualitygroups p=pass (n=%d) and f=failed (n=%d) per site',...
                [sum(cellfun(@(x) numel(x),boxMMCR(1:2:end))),sum(cellfun(@(x) numel(x),boxMMCR(2:2:end)))]),...
                'FontSize',FS(2)); ylabel('rating score in percent','FontSize',FS(2));
    
              % print
              if strcmp(Pfields,'nSIQR')
                print(fh(tfi),fullfile(opt.printdir,sprintf('fig_ROC_anno%d_ordered%d_%s_model%d_cmodel%d_%s',...
                  ai>1,tfi>1,Pfields{pfi},model,cmodel)),opt.res,opt.type,datestr(clock,'YYYYmm'));
              else
                print(fh(tfi),fullfile(opt.printdir,sprintf('fig_ROC_anno%d_ordered%d_%s_%s',...
                  ai>1,tfi>1,Pfields{pfi})),opt.res,opt.type,datestr(clock,'YYYYmm'));
              end
           %   if opt.closefig, close(fd); end
    
    
              % average/std of each group
              if tfi == 1
                fprintf('\nMean and SD of the expert groups in rps: \n');
              else
                fprintf('\nSorted mean and SD of the expert groups in rps: \n');
              end  
              fprintf('Site:      %11s%11s','avg','SD'); for si=MMCRmeansorti(1:2:end), fprintf('%11s', groupsn2{si} ); end; fprintf('\n')
              fprintf('  mean good:'); fprintf('%10.2f ',[mean(cellfun(@(x) median(x) , MMCR(1:2:end))), ...
                std(cellfun(@(x) median(x) , MMCR(1:2:end))),  cellfun(@(x) median(x) , MMCR(1:2:end))]); fprintf('\n');
              fprintf('  std  good:'); fprintf('%10.2f ',[mean(cellfun(@(x) std(x)    , MMCR(1:2:end))), ...
                std(cellfun(@(x) std(x)    , MMCR(1:2:end))),  cellfun(@(x) std(x)    , MMCR(1:2:end))]); fprintf('\n');
              fprintf('  mean bad: '); fprintf('%10.2f ',[mean(cellfun(@(x) median(x) , MMCR(2:2:end))), ...
                std(cellfun(@(x) median(x) , MMCR(2:2:end))),  cellfun(@(x) median(x) , MMCR(2:2:end))]); fprintf('\n');
              fprintf('  std  bad: '); fprintf('%10.2f ',[mean(cellfun(@(x) std(x)    , MMCR(2:2:end))), ...
                std(cellfun(@(x) std(x)    , MMCR(2:2:end))),  cellfun(@(x) std(x)    , MMCR(2:2:end))]); fprintf('\n');
              fprintf('  qant diff:'); fprintf('%10.2f ',[ ...
                mean(cellfun(@(x,y) median(x(x<median(x)) - median(y(y>median(y)))) , MMCR(1:2:end), MMCR(2:2:end) )), ...
                std(cellfun(@(x,y)  median(x(x<median(x)) - median(y(y>median(y)))) , MMCR(1:2:end), MMCR(2:2:end) )), ...
                cellfun(@(x,y)      median(x(x<median(x)) - median(y(y>median(y)))) , MMCR(1:2:end), MMCR(2:2:end)) ]); fprintf('\n'); %#### better !!! ####
              fprintf('\n');
    
      % ################### print to file        
            end
          end
        end
      end
     
     
      
      
      if 0
        %% GM values for expert groups
        %  --------------------------------------------------------------------
        %  problem is here that the group are not balanced for age (and sex)
    %subplot('Position',[b2 pos{2,1}(2) 0.95-b2*1.1 pos{1,1}(4)]); 
        MMCV = cell(0); groupsn2 = cell(''); n = str2double(groups); gxi=1;
        MMCVM = cell(1,2); 
        %sitenx2 = upper('a_bcdef_ghi');
        for si=1:max(Q.site) % [1,3:7,9:max(Q.site)]
          for gi=1:n
            MMCV{gxi} = Q.GMV(M3{gi} & repmat(Q.site'==si,1,size(M3{gi},2))); 
            %sitenx   = Q.siten(Q.site'==si);
            %groupsn2{gxi} = sprintf('%s%s',lower(groupsn{gi}),sitenx2(si)); %sitenx{1});
            %xx = kmeans3D(MMC{gxi},3);
            gxi = gxi+1;
          end
          if 1
            sitemean    = mean(MMCV{gxi-2}); 
            MMCV{gxi-2} = MMCV{gxi-2} - sitemean;
            MMCV{gxi-1} = MMCV{gxi-1} - sitemean; 
            MMCVM{1}    = [MMCVM{1}; MMCV{gxi-2}];
            MMCVM{2}    = [MMCVM{2}; MMCV{gxi-1}]; 
          end
        end
        if 0 
          MMCRmeansorti
        end
        if 0
          %%
          cat_plot_boxplot(MMCVM,struct('boxwidth',-1,'names',{{'pass','failed'}},'usescatter',1,...
            'sort',0,'ylim',[-0.2 0.2],'groupnum',1,'ygrid',1)); set(gca,'FontSize',FS(2)*0.8)
        else      
          cat_plot_boxplot(MMCV,struct('boxwidth',-1,'names',{groupsn2},'usescatter',1,...
            'sort',0,'ylim',[-0.2 0.2],'groupnum',1,'ygrid',1,'groupcolor',gcolor)); set(gca,'FontSize',FS(2)*0.8)
          title(sprintf('GMV by expert quality group',numel(1)),...
            'FontSize',FS(1),'FontWeight','bold'); 
          set(gca,'xTickLabelRotation',90,'Position',[0.05 0.25 0.94 0.70]);
          xlabel(sprintf('qualitygroups p=pass (n=%d) and f=failed (n=%d) per site (A-I)',...
            [sum(cellfun(@(x) numel(x),MMCV(1:2:end))),sum(cellfun(@(x) numel(x),MMCV(2:2:end)))]),...
            'FontSize',FS(2)); ylabel('GMV','FontSize',FS(2));
    
        end
        
        fprintf('GMV diff: %0.2f\n\n',...
          mean( ((cellfun(@mean,MMCV(1:2:end)) - cellfun(@mean,MMCV(2:2:end))) ./ ...
                 (cellfun(@mean,MMCR(1:2:end)) - cellfun(@mean,MMCR(2:2:end)))).^2 ).^0.5 * 100); 
        
      end
         
        
      
      
      %% NCR, ICR, IQR, SIQR, NXIQR?
      IQRfield = IQRfields{IQRfieldi};
    
      fign      = 33 + model + 6*(cmodel-1); 
      fhtmp     = figure(fign); fhtmp.Position(3:4) = [1200 400]; clf(fign);
      if 0 % default
        th      = [0.5:0.02:1, 1:0.01:3, 3.02:0.02:4, 4.05:0.05:5, 5.1:0.1:6.5];   % global IQR threshold test range (school grades)
        cf      = [-1:0.1:-0.1, -0.08:0.02:0, 0:0.01:0.4, 0.42:0.02:0.8, 0.85:0.05:2, 2.1:0.1:3];      % protocoll-specific dIQR threshold test range (school grad range) 
      elseif 1 % default
        th      = [0.5:0.05:1.5, 1.52:0.02:2, 2.01:0.01:3, 3.02:0.02:3.5, 3.55:0.05:4, 4.1:0.1:6.5];   % global IQR threshold test range (school grades)
        cf      = [-1:0.1:-0.5, -0.45:0.05:-0.1, -0.08:0.02:0, 0:0.01:0.6, 0.62:0.02:1.0, 1.05:0.05:1.5, 1.6:0.1:2, 2.2:0.2:3];      % protocoll-specific dIQR threshold test range (school grad range) 
      else % fast test
        th      = [0.5:0.1:1, 1:0.05:3, 3.05:0.1:4, 4.1:0.2:6.5];   % global IQR threshold test range (school grades)
        cf      = [-1:0.2:-0.2, -0.18:0.05:0.5, 0.6:0.1:1, 1.2:0.2:3];   % protocoll-specific dIQR threshold test range (school grad range) 
      end
      erth      = 0.5;            % expert rating threshold
      nontest   = [0.15 0.95];
      usenontest = 1; 
      

      % testsites - a smaller sample of usefull sites is preferable  
      % 3 sites are inoptimal: 
      %  - Calgary (13) where all images have some kind of artifacts
      %  - ABIDE ... (26) with an high-res, high-noise rescans protocol for averaging 
      %  - ABIDE Standford (28) where nearly all scans have motion artifacts
%      testsites = unique(Q.site(1:4)); 
      testsites = setdiff( unique(Q.site) , [13 26 28]); 
     
      for ti = numel(testsites):-1:1
        if sum(Q.site==testsites(ti))<10 || sum(Q.site'==testsites(ti) & Q.group==0)<5 || sum(Q.site'==testsites(ti) & Q.group==1)<5, testsites(ti) = []; end 
      end
      fdn = sprintf('%dsites_model%d_cmodel%d',numel(testsites),model,cmodel);
     
      % (1) general site-unspecific threshold
      % _______________________________________________________________________
      %
      % Process all cases with the global treshold "th" in range of the quality 
      % measure. There is no further selection here - just process everything!
      % _______________________________________________________________________
      if 1 %~exist('sens','var')
        sens  = {nan(numel(th),1),nan(numel(th),1)}; spec = sens; acc = sens; auc = sens; 
        sensg = {nan(numel(th),max(Q.site)); nan(numel(th),max(Q.site))}; 
        specg = sensg; accg = sensg; aucg = sensg;
        for i = 1:numel(th) % apply global IQR tresholds for ROC statistic 
          for ti = 1:2      % train vs. test
            M  = any( repmat(Q.site',1,numel(testsites)) == repmat(testsites,size(Q.site,1)) , 2); 
            M  = M & Q.train == ti-1; 
            if usenontest, M  = M & (P.rating'<nontest(1) | P.rating'>nontest(2)); end
    
            TP = Q.group> erth & Q.(IQRfield) >  th(i); TPs = sum(TP(M));
            FP = Q.group<=erth & Q.(IQRfield) >  th(i); FPs = sum(FP(M));
            TN = Q.group<=erth & Q.(IQRfield) <= th(i); TNs = sum(TN(M));
            FN = Q.group> erth & Q.(IQRfield) <= th(i); FNs = sum(FN(M));
    
            sens{ti}(i) = TPs ./ max(1,TPs + FNs);
            spec{ti}(i) = TNs ./ max(1,TNs + FPs);
            acc{ti}(i)  = (TPs + TNs) / max(1,TPs + FNs + TNs + FPs);
    
            try
              [~,~,~,auc{ti}(i)] = perfcurve( Q.group(M) , Q.(IQRfield)(M) - th(i) , 'true');
            catch
              auc{ti}(i) = 1;  
            end
    
            % site-specific values
            for gi = testsites
              Mt = Q.site' == gi & Q.train == ti-1; 
              if usenontest, Mt = Mt & (P.rating'<nontest(1) | P.rating'>nontest(2)); end
    
              TP = Q.group> erth & Q.(IQRfield) >  th(i); TPs = sum(TP(Mt));
              FP = Q.group<=erth & Q.(IQRfield) >  th(i); FPs = sum(FP(Mt));
              TN = Q.group<=erth & Q.(IQRfield) <= th(i); TNs = sum(TN(Mt));
              FN = Q.group> erth & Q.(IQRfield) <= th(i); FNs = sum(FN(Mt));
    
              sensg{ti}(i,gi) = TPs ./ max(1,TPs + FNs);
              specg{ti}(i,gi) = TNs ./ max(1,TNs + FPs);
              accg{ti}(i,gi)  = (TPs + TNs) / max(1,TPs + FNs + TNs + FPs);
    
              try
                [~,~,~,aucg{ti}(i)] = perfcurve( Q.group(Mt) , Q.(IQRfield)(Mt) - th(i) , 'true');
              catch
                aucg{ti}(i) = 1; 
              end
            end  
          end
        end
      end
    
      
      sens2  = {nan(numel(cf),1),nan(numel(cf),1)}; spec2 = sens2; acc2 = sens2; auc2 = sens2;
      sensg2 = {nan(numel(cf),max(Q.site)); nan(numel(cf),max(Q.site))}; 
      specg2 = sensg2; accg2 = sensg2; aucg2 = sensg2;
      Q.Nmn  = nan(size(Q.IQR)); Q.Nsd = Q.Nmn;  Q.NXIQR = Q.Nmn;
      for i = 1:numel(cf) % apply global IQR tresholds for ROC statistic 
        for ti = 1:2      % train vs. test
          M  = any( repmat(Q.site',1,numel(testsites)) == repmat(testsites,size(Q.site,1)) , 2); 
          M  = M & Q.train == ti-1; 
          if usenontest, M  = M & (P.rating'<nontest(1) | P.rating'>nontest(2)); end
    
          if ti == 1
            [ Q.NXIQR(M) , Q.Nmn(M) , Q.Nsd(M)] = cat_tst_qa_normer( Q.(IQRfield)(M), ...
              struct('model',model,'figure',0,'cmodel',cmodel,'sites',Q.site(M)'));
          else
            for gi = testsites
              Q.Nmn(Q.site == gi) = cat_stat_nanmean( Q.Nmn(Q.site == gi) );
              Q.Nsd(Q.site == gi) = cat_stat_nanmean( Q.Nsd(Q.site == gi) );
            end
            if cmodel == 1
              Q.NXIQR = Q.(IQRfield) - Q.Nmn;
            else
              Q.NXIQR = (Q.(IQRfield) - Q.Nmn) ./ Q.Nsd;
            end
          end
    
          TP = Q.group> erth & Q.NXIQR >  cf(i); TPs = sum(TP(M));
          FP = Q.group<=erth & Q.NXIQR >  cf(i); FPs = sum(FP(M));
          TN = Q.group<=erth & Q.NXIQR <= cf(i); TNs = sum(TN(M));
          FN = Q.group> erth & Q.NXIQR <= cf(i); FNs = sum(FN(M));
    
          sens2{ti}(i) = TPs ./ max(1,TPs + FNs);
          spec2{ti}(i) = TNs ./ max(1,TNs + FPs);
          acc2{ti}(i)  = (TPs + TNs) / max(1,TPs + FNs + TNs + FPs);
    
          try
            [~,~,~,auc2{ti}(i)] = perfcurve( Q.group(M) , Q.NXIQR(M) - cf(i) , 'true');
          catch
            auc2{ti}(i) = 1; 
          end
        end
        % site-specific values
        for ti = 1:2
          for gi = testsites
            Mt = Q.site' == gi & Q.train == ti-1; 
            if usenontest, Mt = Mt & (P.rating'<nontest(1) | P.rating'>nontest(2)); end
  
            TP = Q.group> erth & Q.NXIQR >  cf(i); TPs = sum(TP(Mt));
            FP = Q.group<=erth & Q.NXIQR >  cf(i); FPs = sum(FP(Mt));
            TN = Q.group<=erth & Q.NXIQR <= cf(i); TNs = sum(TN(Mt));
            FN = Q.group> erth & Q.NXIQR <= cf(i); FNs = sum(FN(Mt));
    
            sensg2{ti}(i,gi) = TPs ./ max(1,TPs + FNs);
            specg2{ti}(i,gi) = TNs ./ max(1,TNs + FPs);
            accg2{ti}(i,gi)  = (TPs + TNs) / max(1,TPs + FNs + TNs + FPs);
    
            try
              [~,~,~,aucg2{ti}(i,gi)] = perfcurve( Q.group(Mt) , Q.NXIQR(Mt) - cf(i) , 'true');
            catch
              aucg2{ti}(i,gi) = 1; 
            end
          end  
        end
      end
      

      %% subfigure 1 with global value/threshold  
      [~,mxthi] = max( cat_stat_nanmean([ spec{2}  , sens{2}  ] , 2) ); % train
      [~,mxcci] = max( cat_stat_nanmean([ spec2{2} , sens2{2} ] , 2) ); % train
     
      subplot('Position',[0.05 0.1 0.28 0.85]); 
      hold on; box on;
      for gi = testsites
        plot(th,sensg{2}(:,gi),'color',min(1,max(0,[0.0 0.4 0.0]/2+0.6)),'linewidth',0.5);
        plot(th,specg{2}(:,gi),'color',min(1,max(0,[0.8 0.0 0.0]/2+0.6)),'linewidth',0.5);
      end
      plot(th,sens{2},'color',[0.0 0.5 0.0],'Linestyle','-','linewidth',1.5);  % [0.0 0.4 1.0]
      plot(th,spec{2},'color',[0.8 0.0 0.0],'Linestyle','-','linewidth',1.5);  % [0.9 0.0 0.2]
      plot( [th(mxthi) th(mxthi)] , [0 1.02] ,'color',[0.8 0.0 0.6],'Linestyle','-','linewidth',2);
      hold off; xlim([min(th),max(th)]); ylim([0 1.02]);
      title('sensitivity/specificity of site-unspecific threshold','FontSize',FS(1),'FontWeight','bold'); 
      xlabel([strrep(IQRfield,'_','\_') ' (site-unspecific)'],'FontSize',FS(2) ); ylabel('sensitivity/specificity','FontSize',FS(2)); 
      legend({'average  sensitivity','average  specificity','site  sensitivity ','site  specificity'},...
        'Location','Northeast'); legend('boxoff');
      
      
      % subfigure 2 with normalized value
      subplot('Position',[0.38 0.1 0.28 0.85]); 
      hold on; box on;
      for gi = testsites
        plot(cf,sensg2{2}(:,gi),'color',min(1,max(0,[0.0 0.4 0.0]/2+0.6)),'linewidth',0.5);
        plot(cf,specg2{2}(:,gi),'color',min(1,max(0,[0.8 0.0 0.0]/2+0.6)),'linewidth',0.5);
      end
      plot(cf,sens2{2},'color',[0.0 0.5 0.0],'Linestyle','-','linewidth',1.5); 
      plot(cf,spec2{2},'color',[0.8 0.0 0.0],'Linestyle','-','linewidth',1.5); 
      plot( [cf(mxcci) cf(mxcci)] , [0 1.02] ,'color',[0.0 0.4 0.8],'Linestyle','-','linewidth',2);
      hold off; xlim([min(cf),max(cf)]); ylim([0 1.02]);
      title(sprintf('sensitivity/specificity of site-specific threshold (N_{SITE}=%d,N_{SCANS}=%d)', ...
        numel(testsites),sum(~Q.train)),'FontSize',FS(1),'FontWeight','bold'); 
      xlabel(['N' strrep(IQRfield,'_','\_') ' (site-specific)'],'FontSize',FS(2)); ylabel('sensitivity/specificity','FontSize',FS(2)); 
      legend({'average  sensitivity','average  specificity','site  sensitivity ','site  specificity'},...
        'Location','Northeast'); legend('boxoff');
     

      % subfigure 3 with ROC
      subplot('Position',[0.70 0.1 0.28 0.85],'box','on'); cla; hold on; grid on;
      plot(1- mean(cell2mat(spec) ,2), mean(cell2mat(sens) ,2) ,'color',[0.8 0.0 0.6],'linewidth',1.5); 
      plot(1- mean(cell2mat(spec2),2), mean(cell2mat(sens2),2) ,'color',[0.0 0.4 0.8],'linewidth',1.5); 
      plot(1-spec{1}  ,sens{1}  ,'color',[0.8 0.0 0.6]/2+0.6,'linewidth',0.8); 
      plot(1-spec{2}  ,sens{2}  ,'color',[0.8 0.0 0.6]/2+0.6,'linewidth',0.8); 
      plot(1-spec2{1} ,sens2{1} ,'color',[0.0 0.4 0.8]/2+0.6,'linewidth',0.8); 
      plot(1-spec2{2} ,sens2{2} ,'color',[0.0 0.4 0.8]/2+0.6,'linewidth',0.8); 
      hold off; ylim([0.5 1.004]); xlim([-0.004 0.65]); ylim([0 1.0]); xlim([0 1])
      title('receiver operating characteristic (ROC)','FontSize',FS(1),'FontWeight','bold'); 
      xlabel('False positive rate (1-specificity)','FontSize',FS(2)); 
      ylabel('True positive rate (sensitivity)','FontSize',FS(2)); 
      legend({...sprintf('global threshold run1: AUC=%0.3f \n (th=%0.2f rps, ACC=%0.3f)',...
              ...  auc{1}(mxthi), mark2rps(th(mxthi)), acc{1}(mxthi) ),...
              ...sprintf('global threshold run2: AUC=%0.3f \n  (th=%0.2f rps, ACC=%0.3f)',...
              ...  auc{2}(mxthi), mark2rps(th(mxthi)), acc{2}(mxthi) ),...
              sprintf('%s (site-unspecific threshold): AUC=%0.3f \n (th=%0.2f rps, ACC=%0.3f)',...
                strrep(IQRfield,'_','\_'), ...
                mean([auc{1}(mxthi),auc{2}(mxthi)]), th(mxthi), mean([acc{1}(mxthi),acc{2}(mxthi)]) ),...
              ...sprintf('site specific threshold run1: AUC=%0.3f \n (cf=%0.2f, ACC=%0.3f)',... 
              ...  auc2{1}(mxcci), cf(mxcci), acc2{1}(mxcci) ),...
              ...sprintf('site specific threshold run2: AUC=%0.3f \n (cf=%0.2f, ACC=%0.3f)',... 
              ...  auc2{2}(mxcci), cf(mxcci), acc2{2}(mxcci) ),...
              sprintf('N%s (site-specific threshold): AUC=%0.3f \n (th=%0.2f rps, ACC=%0.3f)',...
                strrep(IQRfield,'_','\_'), ...
                mean([auc2{1}(mxcci),auc2{2}(mxcci)]), cf(mxcci), mean([acc2{1}(mxcci),acc2{2}(mxcci)]) ),...
              },'Location','southeast');
      
      %
      print(fhtmp, fullfile(opt.printdir,sprintf('fig_ROC_%s_%s_%s_%s', IQRfield, qafile, fdn, datestr(clock,'YYYYmm') )), opt.res, opt.type);
      %
      %if opt.closefig, close(fhtmp); end
    
    %%  print(fh1,fullfile(opt.printdir,sprintf('fig_ROC_%s',ff)),opt.res,opt.type);
             
      
     
      if 1
        %%
        opt.MarkColor  = cat_io_colormaps('marks+',10); 
        fprintf('\nResults per site (Model: %d-%d):',model,cmodel); 
        fprintf('\n  Site:%8s | ','All'); for ti = testsites, fprintf('%8s | ',sprintf('%d-%s',ti,Q.siten{find(Q.site==ti,1,'first')})); end
        fprintf('\n  N:   %8d | ',-1); for ti = testsites, fprintf('%8d | ',sum(Q.site==ti)); end
        fprintf('\n  AUC1:%8.3f | ',auc2{1}(mxcci)); for ti = testsites, cat_io_cprintf( opt.MarkColor( max(1,min(size(opt.MarkColor,1), round( 10 - cat_stat_nanmean( aucg2{1}(mxcci,ti) , 2 ) * size(opt.MarkColor,1)))) ,:), ...
                    sprintf('%8.2f | ',cat_stat_nanmean( aucg2{1}(mxcci,ti) , 2 ))); end
        fprintf('\n  AUC2:%8.3f | ',auc2{2}(mxcci)); for ti = testsites, cat_io_cprintf( opt.MarkColor( max(1,min(size(opt.MarkColor,1), round( 10 - cat_stat_nanmean( aucg2{2}(mxcci,ti) , 2 ) * size(opt.MarkColor,1)))) ,:), ...
                    sprintf('%8.2f | ',cat_stat_nanmean( aucg2{2}(mxcci,ti) , 2 ))); end
        fprintf('\n  ACC1:%8.3f | ',acc2{1}(mxcci)); for ti = testsites, cat_io_cprintf( opt.MarkColor( max(1,min(size(opt.MarkColor,1), round( 10 - cat_stat_nanmean( accg2{1}(mxcci,ti) , 2 ) * size(opt.MarkColor,1)))) ,:), ...
                    sprintf('%8.2f | ',cat_stat_nanmean( accg2{1}(mxcci,ti) , 2 ))); end
        fprintf('\n  ACC2:%8.3f | ',acc2{2}(mxcci)); for ti = testsites, cat_io_cprintf( opt.MarkColor( max(1,min(size(opt.MarkColor,1), round( 10 - cat_stat_nanmean( accg2{2}(mxcci,ti) , 2 ) * size(opt.MarkColor,1)))) ,:), ...
                    sprintf('%8.2f | ',cat_stat_nanmean( accg2{2}(mxcci,ti) , 2 ))); end
        fprintf('\n  MN:  %8.3f | ',mean([auc2{1}(mxcci),auc2{2}(mxcci),acc2{1}(mxcci),acc2{2}(mxcci)])); 
        for ti = testsites
          cat_io_cprintf( opt.MarkColor( max(1,min(size(opt.MarkColor,1), ...
            round( 10 - cat_stat_nanmean( [aucg2{1}(mxcci,ti), aucg2{2}(mxcci,ti), accg2{1}(mxcci,ti) , accg2{2}(mxcci,ti)] , 2 ) * size(opt.MarkColor,1)))) ,:), ...
                    sprintf('%8.2f | ',cat_stat_nanmean( [aucg2{1}(mxcci,ti), aucg2{2}(mxcci,ti), accg2{1}(mxcci,ti) , accg2{2}(mxcci,ti)] , 2 ))); 
        end
        fprintf('\n  MNG: %8.3f | ',mean([auc{1}(mxcci),auc{2}(mxcci),acc{1}(mxcci),acc{2}(mxcci)])); 
        for ti = testsites
          cat_io_cprintf( opt.MarkColor( max(1,min(size(opt.MarkColor,1), ...
            round( 10 - cat_stat_nanmean( [aucg{1}(mxcci,ti), aucg{2}(mxcci,ti), accg{1}(mxcci,ti) , accg{2}(mxcci,ti)] , 2 ) * size(opt.MarkColor,1)))) ,:), ...
                    sprintf('%8.2f | ',cat_stat_nanmean( [aucg{1}(mxcci,ti), aucg{2}(mxcci,ti), accg{1}(mxcci,ti) , accg{2}(mxcci,ti)] , 2 ))); 
        end
        fprintf('\n');
      end
     
      
    
      if 0
        %% -- correlations ---------------------------------------------------
        stime = cat_io_cmd('  Estimate correlations:','g5','',1,stime); 
      
        M  = true(size(Q.group,1),1); 
        M2 = repmat(M,1,size(Q.NCR,2)); M2 = M2(:);
        QMnames = {'NCR','ICR','RES','IQR','group','rCSFV','rGMV','rWMV'}; %,'method'
        QMcorrs = [ Q.NCR(M2) , Q.ICR(M2) , Q.res_RMS(M2) , Q.(IQRfield)(M2), ...
             Q.group(M2) Q.CMV(M2), Q.GMV(M2), Q.WMV(M2)];
         [C.QM.r, C.QM.p] = corr(QMcorrs,'type','Spearman');
         C.QM.rtab = [ [{''} QMnames]; [QMnames',num2cell(round((C.QM.r)*10000)/10000)] ];
         C.QM.ptab = [ [{''} QMnames]; [QMnames',num2cell(round((C.QM.p)*10000)/10000)] ];
      
        try C.QM = rmfield(C.QM,'txt'); end %#ok<*TRYNC> 
        T.ixiPPC = sprintf('\nPCC on IXI database:\n%s\n',...
          '------------------------------------------------------------------------------');
        for di=1:size(C.QM.rtab,2)-1, T.ixiPPC = sprintf('%s%8s   ',T.ixiPPC,C.QM.rtab{1,di}); end 
        for dj=3:size(C.QM.rtab,1)
          T.ixiPPC = sprintf('%s\n%8s   ',T.ixiPPC,C.QM.rtab{1,dj});
          for di=2:size(C.QM.rtab,2)
            if     C.QM.ptab{dj,di}<0.000001, star='***';
            elseif C.QM.ptab{dj,di}<0.000100, star='** ';
            elseif C.QM.ptab{dj,di}<0.010000, star='*  ';
            else                              star='   ';
            end
            if di<dj
              T.ixiPPC = sprintf('%s%8.4f%3s',T.ixiPPC,C.QM.rtab{dj,di},star);
            else
              T.ixiPPC = sprintf('%s%8s%3s',T.ixiPPC,'','');
            end
          end
        end
        C.QM.txt = T.ixiPPC;
        C.QM.txt = [C.QM.txt '\n\n\n'];
      
        % MANOVA
        [p,st] = anovan(Q.group(M2) ,{Q.NCR(M2) Q.ICR(M2) Q.res_RMS(M2) Q.(IQRfield)(M2) Q.methods(M2)},...
          'continuous',[1 2 3 4],'alpha',0.05,'display','off','varname',{'NCR','ICR','RES',IQRfield,'method'});
        %
        QMnames = {'NCR' 'ICR' 'res_RMS' 'methods' IQRfield,'GMV'};
        %C.QM.txt =  [C.QM.txt;C.QM.etab];
        clear px stx;
        for q1=1:numel(QMnames)
          [px{q1},stx{q1}] = anova1(Q.(QMnames{q1})(M2),Q.group(M2),'off');
          txt = sprintf('%10s: F=%8.3f p=%8.8f\n',QMnames{q1},stx{q1}{2,5},px{q1});
          C.QM.txt = [C.QM.txt txt]; 
          %C.QM.txt{size(C.QM.txt,1)+1:size(C.QM.txt,1)+1+size(stx{q1},1),1:size(stx{q1},2)} = stx{q1};
        end
        C.QM.txt = [C.QM.txt '\n\n'];
        %
        for i=1:size(st,1)
          for j=1:size(st,2)
            if isnumeric(st{i,j})
              C.QM.txt = sprintf('%s%12.7f  ',C.QM.txt,st{i,j});
            else
              C.QM.txt = sprintf('%s%12s  ',C.QM.txt,st{i,j});
            end
          end
          C.QM.txt = sprintf('%s\n',C.QM.txt);
        end
      
        
        %% -- save images ----------------------------------------------------
        stime = cat_io_cmd('  Print figures:','g5','',1,stime);
      
        fprintf('%s',C.QM.txt);
        
        ff = sprintf('qa_grouptest_%s_%s',qafile,datestr(clock,'YYYYmm'));
        f=fopen(fullfile(opt.printdir,sprintf('tst_%s.csv',ff)),'w'); if f~=-1, fprintf(f,'%s',C.QM.txt); fclose(f); end
        f=fopen(fullfile(opt.resdir  ,sprintf('tst_%s.csv',ff)),'w'); if f~=-1, fprintf(f,'%s',C.QM.txt); fclose(f); end
      
        print(fhtmp,fullfile(opt.printdir,sprintf('fig_%s_%s',ff)),opt.res,opt.type,datestr(clock,'YYYYmm'),'-dpng');
        print(fhtmp,fullfile(opt.resdir  ,sprintf('fig_%s_%s',ff)),opt.res,opt.type,datestr(clock,'YYYYmm'),'-dpng');
        
        %  if opt.closefig, close(fh1); end
        cat_io_cmd(' ','g5','',1,stime);
       
      end
    end
  end

function [good,bad] = get_rating
  %% further QA groups
  % == INDI AAa == ... many slice artifacts
  % == INDI AAb == movement and slice artifacts
  %{
  good.AAb = {
    'raw/INDI/HC/AAb/sub00306/INDI_HC_AAb_sub00306_T1_SD000000-RS00'
    'raw/INDI/HC/AAb/sub11043/INDI_HC_AAb_sub11043_T1_SD000000-RS00'
    'raw/INDI/HC/AAb/sub05580/INDI_HC_AAb_sub05580_T1_SD000000-RS00'
    'raw/INDI/HC/AAb/sub57025/INDI_HC_AAb_sub57025_T1_SD000000-RS00'
    'raw/INDI/HC/AAb/sub73168/INDI_HC_AAb_sub73168_T1_SD000000-RS00'
    'raw/INDI/HC/AAb/sub73812/INDI_HC_AAb_sub73812_T1_SD000000-RS00'
    'raw/INDI/HC/AAb/sub85257/INDI_HC_AAb_sub85257_T1_SD000000-RS00'
    'raw/INDI/HC/AAb/sub90950/INDI_HC_AAb_sub90950_T1_SD000000-RS00'
    }; 
  bad.AAb = {
    'raw/INDI/HC/AAb/sub07921/INDI_HC_AAb_sub07921_T1_SD000000-RS00'
    'raw/INDI/HC/AAb/sub39923/INDI_HC_AAb_sub39923_T1_SD000000-RS00'
    'raw/INDI/HC/AAb/sub53269/INDI_HC_AAb_sub53269_T1_SD000000-RS00'
    'raw/INDI/HC/AAb/sub57196/INDI_HC_AAb_sub57196_T1_SD000000-RS00'
    'raw/INDI/HC/AAb/sub59573/INDI_HC_AAb_sub59573_T1_SD000000-RS00'
    'raw/INDI/HC/AAb/sub64831/INDI_HC_AAb_sub64831_T1_SD000000-RS00'
    'raw/INDI/HC/AAb/sub78151/INDI_HC_AAb_sub78151_T1_SD000000-RS00'
    'raw/INDI/HC/AAb/sub97518/INDI_HC_AAb_sub97518_T1_SD000000-RS00'
    };
  %}
  % == NKI == movement artifacts
  good.NKI = {
      'raw/NKI/HC/NKI/1013090/NKI_HC_NKI_1013090_T1_SD000000-RS00'
      ...'raw/NKI/HC/NKI/1034049/NKI_HC_NKI_1034049_T1_SD000000-RS00' % light
      'raw/NKI/HC/NKI/1097782/NKI_HC_NKI_1097782_T1_SD000000-RS00'
      ...'raw/NKI/HC/NKI/1103722/NKI_HC_NKI_1103722_T1_SD000000-RS00' % light
      ...'raw/NKI/HC/NKI/1125244/NKI_HC_NKI_1125244_T1_SD000000-RS00' % light
      'raw/NKI/HC/NKI/1143655/NKI_HC_NKI_1143655_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/1192435/NKI_HC_NKI_1192435_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/1264721/NKI_HC_NKI_1264721_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/1271401/NKI_HC_NKI_1271401_T1_SD000000-RS00'
      ...'raw/NKI/HC/NKI/1285465/NKI_HC_NKI_1285465_T1_SD000000-RS00' % light
      'raw/NKI/HC/NKI/1288657/NKI_HC_NKI_1288657_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/1292527/NKI_HC_NKI_1292527_T1_SD000000-RS00'
      ... 1309 % light
      'raw/NKI/HC/NKI/1322220/NKI_HC_NKI_1322220_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/1334556/NKI_HC_NKI_1334556_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/1334913/NKI_HC_NKI_1334913_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/1351931/NKI_HC_NKI_1351931_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/1366839/NKI_HC_NKI_1366839_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/1382333/NKI_HC_NKI_1382333_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/1427581/NKI_HC_NKI_1427581_T1_SD000000-RS00'
      ... 'raw/NKI/HC/NKI/1476437/NKI_HC_NKI_1476437_T1_SD000000-RS00' % light
      ... 'raw/NKI/HC/NKI/1494230/NKI_HC_NKI_1494230_T1_SD000000-RS00' % light
      'raw/NKI/HC/NKI/1508861/NKI_HC_NKI_1508861_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/1522653/NKI_HC_NKI_1522653_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/1523112/NKI_HC_NKI_1523112_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/1591238/NKI_HC_NKI_1591238_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/1601547/NKI_HC_NKI_1601547_T1_SD000000-RS00'
      ... 1616 % light
      'raw/NKI/HC/NKI/1654606/NKI_HC_NKI_1654606_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/1709141/NKI_HC_NKI_1709141_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/1713513/NKI_HC_NKI_1713513_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/1734822/NKI_HC_NKI_1734822_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/1742775/NKI_HC_NKI_1742775_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/1757866/NKI_HC_NKI_1757866_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/1781836/NKI_HC_NKI_1781836_T1_SD000000-RS00'
      ... 'raw/NKI/HC/NKI/1834799/NKI_HC_NKI_1834799_T1_SD000000-RS00' % light
      'raw/NKI/HC/NKI/1915832/NKI_HC_NKI_1915832_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/1931386/NKI_HC_NKI_1931386_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/1935851/NKI_HC_NKI_1935851_T1_SD000000-RS00'
      ... 1940 % light
      'raw/NKI/HC/NKI/1943306/NKI_HC_NKI_1943306_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/1961098/NKI_HC_NKI_1961098_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/2021454/NKI_HC_NKI_2021454_T1_SD000000-RS00'
      ... 'raw/NKI/HC/NKI/2041109/NKI_HC_NKI_2041109_T1_SD000000-RS00' % light
      'raw/NKI/HC/NKI/2071045/NKI_HC_NKI_2071045_T1_SD000000-RS00'
      ... 'raw/NKI/HC/NKI/2113846/NKI_HC_NKI_2113846_T1_SD000000-RS00' % light
      'raw/NKI/HC/NKI/2139745/NKI_HC_NKI_2139745_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/2160826/NKI_HC_NKI_2160826_T1_SD000000-RS00' % is ok - WMHs?
      'raw/NKI/HC/NKI/2176073/NKI_HC_NKI_2176073_T1_SD000000-RS00'  
      'raw/NKI/HC/NKI/2221971/NKI_HC_NKI_2221971_T1_SD000000-RS00'
      ... 'raw/NKI/HC/NKI/2258252/NKI_HC_NKI_2258252_T1_SD000000-RS00' % light
      'raw/NKI/HC/NKI/2286053/NKI_HC_NKI_2286053_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/2286816/NKI_HC_NKI_2286816_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/2328270/NKI_HC_NKI_2328270_T1_SD000000-RS00' % light but more WMHs
      ... 'raw/NKI/HC/NKI/2357976/NKI_HC_NKI_2357976_T1_SD000000-RS00' % light
      'raw/NKI/HC/NKI/2362594/NKI_HC_NKI_2362594_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/2403029/NKI_HC_NKI_2403029_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/2436415/NKI_HC_NKI_2436415_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/2457957/NKI_HC_NKI_2457957_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/2475376/NKI_HC_NKI_2475376_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/2483617/NKI_HC_NKI_2483617_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/2513310/NKI_HC_NKI_2513310_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/2524225/NKI_HC_NKI_2524225_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/2531270/NKI_HC_NKI_2531270_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/2541244/NKI_HC_NKI_2541244_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/2631208/NKI_HC_NKI_2631208_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/2652676/NKI_HC_NKI_2652676_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/2678751/NKI_HC_NKI_2678751_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/2788776/NKI_HC_NKI_2788776_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/2842950/NKI_HC_NKI_2842950_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/2861923/NKI_HC_NKI_2861923_T1_SD000000-RS00'
      ... 'raw/NKI/HC/NKI/2864064/NKI_HC_NKI_2864064_T1_SD000000-RS00' % light
      ... 'raw/NKI/HC/NKI/2868630/NKI_HC_NKI_2868630_T1_SD000000-RS00' % light
      ... 'raw/NKI/HC/NKI/2882334/NKI_HC_NKI_2882334_T1_SD000000-RS00' % light
      ... 'raw/NKI/HC/NKI/2970212/NKI_HC_NKI_2970212_T1_SD000000-RS00' % light
      'raw/NKI/HC/NKI/2986140/NKI_HC_NKI_2986140_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3033715/NKI_HC_NKI_3033715_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3070385/NKI_HC_NKI_3070385_T1_SD000000-RS00'
      ... 'raw/NKI/HC/NKI/3081896/NKI_HC_NKI_3081896_T1_SD000000-RS00' % light
      'raw/NKI/HC/NKI/3094481/NKI_HC_NKI_3094481_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3106263/NKI_HC_NKI_3106263_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3166395/NKI_HC_NKI_3166395_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3209959/NKI_HC_NKI_3209959_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3261820/NKI_HC_NKI_3261820_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3304648/NKI_HC_NKI_3304648_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3329569/NKI_HC_NKI_3329569_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3346545/NKI_HC_NKI_3346545_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3362208/NKI_HC_NKI_3362208_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3366302/NKI_HC_NKI_3366302_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3374719/NKI_HC_NKI_3374719_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3409284/NKI_HC_NKI_3409284_T1_SD000000-RS00'
      ... 'raw/NKI/HC/NKI/3431295/NKI_HC_NKI_3431295_T1_SD000000-RS00' % med
      'raw/NKI/HC/NKI/3431940/NKI_HC_NKI_3431940_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3466763/NKI_HC_NKI_3466763_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3492484/NKI_HC_NKI_3492484_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3566919/NKI_HC_NKI_3566919_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3606220/NKI_HC_NKI_3606220_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3611212/NKI_HC_NKI_3611212_T1_SD000000-RS00'
      ...'raw/NKI/HC/NKI/3695138/NKI_HC_NKI_3695138_T1_SD000000-RS00' %light
      'raw/NKI/HC/NKI/3701005/NKI_HC_NKI_3701005_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3729976/NKI_HC_NKI_3729976_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3755111/NKI_HC_NKI_3755111_T1_SD000000-RS00' % very light
      'raw/NKI/HC/NKI/3808535/NKI_HC_NKI_3808535_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3811036/NKI_HC_NKI_3811036_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3815065/NKI_HC_NKI_3815065_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3827479/NKI_HC_NKI_3827479_T1_SD000000-RS00' % ok
      'raw/NKI/HC/NKI/3848143/NKI_HC_NKI_3848143_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3875444/NKI_HC_NKI_3875444_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3893245/NKI_HC_NKI_3893245_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3911767/NKI_HC_NKI_3911767_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3927656/NKI_HC_NKI_3927656_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3934325/NKI_HC_NKI_3934325_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3951328/NKI_HC_NKI_3951328_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3965194/NKI_HC_NKI_3965194_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/3989122/NKI_HC_NKI_3989122_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/4006955/NKI_HC_NKI_4006955_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/4015919/NKI_HC_NKI_4015919_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/4077433/NKI_HC_NKI_4077433_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/4087509/NKI_HC_NKI_4087509_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/4089896/NKI_HC_NKI_4089896_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/4097119/NKI_HC_NKI_4097119_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/4100790/NKI_HC_NKI_4100790_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/4119751/NKI_HC_NKI_4119751_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/4126393/NKI_HC_NKI_4126393_T1_SD000000-RS00'
      ... 'raw/NKI/HC/NKI/4143704/NKI_HC_NKI_4143704_T1_SD000000-RS00' % light
      'raw/NKI/HC/NKI/4154799/NKI_HC_NKI_4154799_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/4188346/NKI_HC_NKI_4188346_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/4277600/NKI_HC_NKI_4277600_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/4290056/NKI_HC_NKI_4290056_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/4323037/NKI_HC_NKI_4323037_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/4417007/NKI_HC_NKI_4417007_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/4487215/NKI_HC_NKI_4487215_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/4710111/NKI_HC_NKI_4710111_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/4791943/NKI_HC_NKI_4791943_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/5161117/NKI_HC_NKI_5161117_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/5279162/NKI_HC_NKI_5279162_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/5461463/NKI_HC_NKI_5461463_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/6471972/NKI_HC_NKI_6471972_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/6539040/NKI_HC_NKI_6539040_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/6692978/NKI_HC_NKI_6692978_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/6913939/NKI_HC_NKI_6913939_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/7055197/NKI_HC_NKI_7055197_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/8317157/NKI_HC_NKI_8317157_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/9006154/NKI_HC_NKI_9006154_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/9100911/NKI_HC_NKI_9100911_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/9536886/NKI_HC_NKI_9536886_T1_SD000000-RS00'
      ... 'raw/NKI/HC/NKI/9537916/NKI_HC_NKI_9537916_T1_SD000000-RS00' %light
      'raw/NKI/HC/NKI/9630905/NKI_HC_NKI_9630905_T1_SD000000-RS00'
      'raw/NKI/HC/NKI/9645370/NKI_HC_NKI_9645370_T1_SD000000-RS00'  
    }; 
  bad.NKI = {
    'raw/NKI/HC/NKI/1525492/NKI_HC_NKI_1525492_T1_SD000000-RS00'
    'raw/NKI/HC/NKI/1534265/NKI_HC_NKI_1534265_T1_SD000000-RS00'
    'raw/NKI/HC/NKI/1764982/NKI_HC_NKI_1764982_T1_SD000000-RS00' % MA
    'raw/NKI/HC/NKI/1875434/NKI_HC_NKI_1875434_T1_SD000000-RS00'
    'raw/NKI/HC/NKI/1898228/NKI_HC_NKI_1898228_T1_SD000000-RS00'
    'raw/NKI/HC/NKI/1903493/NKI_HC_NKI_1903493_T1_SD000000-RS00'
    'raw/NKI/HC/NKI/2005303/NKI_HC_NKI_2005303_T1_SD000000-RS00' % f
    'raw/NKI/HC/NKI/2465771/NKI_HC_NKI_2465771_T1_SD000000-RS00'
    'raw/NKI/HC/NKI/2494762/NKI_HC_NKI_2494762_T1_SD000000-RS00'
    'raw/NKI/HC/NKI/2505567/NKI_HC_NKI_2505567_T1_SD000000-RS00' % f
    'raw/NKI/HC/NKI/2674565/NKI_HC_NKI_2674565_T1_SD000000-RS00'
    'raw/NKI/HC/NKI/2675807/NKI_HC_NKI_2675807_T1_SD000000-RS00' % MAs
    'raw/NKI/HC/NKI/2784584/NKI_HC_NKI_2784584_T1_SD000000-RS00' % f  
    'raw/NKI/HC/NKI/2915821/NKI_HC_NKI_2915821_T1_SD000000-RS00' % f
    'raw/NKI/HC/NKI/3328556/NKI_HC_NKI_3328556_T1_SD000000-RS00'
    'raw/NKI/HC/NKI/3773269/NKI_HC_NKI_3773269_T1_SD000000-RS00'
    'raw/NKI/HC/NKI/3795193/NKI_HC_NKI_3795193_T1_SD000000-RS00'
    ... 'raw/NKI/HC/NKI/3895064/NKI_HC_NKI_3895064_T1_SD000000-RS00' % very light
    'raw/NKI/HC/NKI/4026825/NKI_HC_NKI_4026825_T1_SD000000-RS00' % f :(
    'raw/NKI/HC/NKI/4176156/NKI_HC_NKI_4176156_T1_SD000000-RS00'
    'raw/NKI/HC/NKI/4230470/NKI_HC_NKI_4230470_T1_SD000000-RS00'
    'raw/NKI/HC/NKI/5844518/NKI_HC_NKI_5844518_T1_SD000000-RS00' % f :(
    'raw/NKI/HC/NKI/8628499/NKI_HC_NKI_8628499_T1_SD000000-RS00' % f :(
    'raw/NKI/HC/NKI/9421819/NKI_HC_NKI_9421819_T1_SD000000-RS00' % f :(
    ... 'raw/NKI/HC/NKI/9780496/NKI_HC_NKI_9780496_T1_SD000000-RS00' % light
    };
  % == KKI == movement artifacts
  %{
  bad.KKI = {
      'raw/ADHD200/ADHD-COM/KKI/4362730/ADHD200_ADHD-COM_KKI_4362730_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/KKI/8337695/ADHD200_ADHD-COM_KKI_8337695_T1_SD000000-RS00'
      'raw/ADHD200/NG/KKI/0020008/ADHD200_NG_KKI_0020008_T1_SD000000-RS00'
    }; 
  %}
  good.NYC = {
      'raw/ADHD200/ADHD-COM/NYC/0010012/ADHD200_ADHD-COM_NYC_0010012_T1_SD000000-RS00'
      ... 'raw/ADHD200/ADHD-COM/NYC/0010013/ADHD200_ADHD-COM_NYC_0010013_T1_SD000000-RS00' % bad defaceing
      'raw/ADHD200/ADHD-COM/NYC/0010017/ADHD200_ADHD-COM_NYC_0010017_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010019/ADHD200_ADHD-COM_NYC_0010019_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010025/ADHD200_ADHD-COM_NYC_0010025_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010026/ADHD200_ADHD-COM_NYC_0010026_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010029/ADHD200_ADHD-COM_NYC_0010029_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010037/ADHD200_ADHD-COM_NYC_0010037_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010040/ADHD200_ADHD-COM_NYC_0010040_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010042/ADHD200_ADHD-COM_NYC_0010042_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010044/ADHD200_ADHD-COM_NYC_0010044_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010048/ADHD200_ADHD-COM_NYC_0010048_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010049/ADHD200_ADHD-COM_NYC_0010049_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010050/ADHD200_ADHD-COM_NYC_0010050_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010060/ADHD200_ADHD-COM_NYC_0010060_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010064/ADHD200_ADHD-COM_NYC_0010064_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010072/ADHD200_ADHD-COM_NYC_0010072_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010073/ADHD200_ADHD-COM_NYC_0010073_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010085/ADHD200_ADHD-COM_NYC_0010085_T1_SD000000-RS00' 
      'raw/ADHD200/ADHD-COM/NYC/0010086/ADHD200_ADHD-COM_NYC_0010086_T1_SD000000-RS00' 
      'raw/ADHD200/ADHD-COM/NYC/0010090/ADHD200_ADHD-COM_NYC_0010090_T1_SD000000-RS00' 
      ... 'raw/ADHD200/ADHD-COM/NYC/0010091/ADHD200_ADHD-COM_NYC_0010091_T1_SD000000-RS00' % light
      ... 'raw/ADHD200/ADHD-COM/NYC/0010095/ADHD200_ADHD-COM_NYC_0010095_T1_SD000000-RS00' % light
      'raw/ADHD200/ADHD-COM/NYC/0010101/ADHD200_ADHD-COM_NYC_0010101_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010103/ADHD200_ADHD-COM_NYC_0010103_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010109/ADHD200_ADHD-COM_NYC_0010109_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010118/ADHD200_ADHD-COM_NYC_0010118_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010126/ADHD200_ADHD-COM_NYC_0010126_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/1057962/ADHD200_ADHD-COM_NYC_1057962_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/2107638/ADHD200_ADHD-COM_NYC_2107638_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/1187766/ADHD200_ADHD-COM_NYC_1187766_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/1283494/ADHD200_ADHD-COM_NYC_1283494_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/1497055/ADHD200_ADHD-COM_NYC_1497055_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/1517240/ADHD200_ADHD-COM_NYC_1517240_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/1780174/ADHD200_ADHD-COM_NYC_1780174_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/1992284/ADHD200_ADHD-COM_NYC_1992284_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/2054438/ADHD200_ADHD-COM_NYC_2054438_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/2230510/ADHD200_ADHD-COM_NYC_2230510_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/2260910/ADHD200_ADHD-COM_NYC_2260910_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/2306976/ADHD200_ADHD-COM_NYC_2306976_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/2497695/ADHD200_ADHD-COM_NYC_2497695_T1_SD000000-RS00'
      ... 'raw/ADHD200/ADHD-COM/NYC/2682736/ADHD200_ADHD-COM_NYC_2682736_T1_SD000000-RS00' % light
      'raw/ADHD200/ADHD-COM/NYC/2741068/ADHD200_ADHD-COM_NYC_2741068_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/2950672/ADHD200_ADHD-COM_NYC_2950672_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/3235580/ADHD200_ADHD-COM_NYC_3235580_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/3457975/ADHD200_ADHD-COM_NYC_3457975_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/3601861/ADHD200_ADHD-COM_NYC_3601861_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/3679455/ADHD200_ADHD-COM_NYC_3679455_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/4116166/ADHD200_ADHD-COM_NYC_4116166_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/4154672/ADHD200_ADHD-COM_NYC_4154672_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/4187857/ADHD200_ADHD-COM_NYC_4187857_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/5164727/ADHD200_ADHD-COM_NYC_5164727_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/5971050/ADHD200_ADHD-COM_NYC_5971050_T1_SD000000-RS00'
      ... 'raw/ADHD200/ADHD-COM/NYC/9326955/ADHD200_ADHD-COM_NYC_9326955_T1_SD000000-RS00' % light
      'raw/ADHD200/ADHD-HIM/NYC/0010041/ADHD200_ADHD-HIM_NYC_0010041_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010001/ADHD200_ADHD-INA_NYC_0010001_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010002/ADHD200_ADHD-INA_NYC_0010002_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010007/ADHD200_ADHD-INA_NYC_0010007_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010011/ADHD200_ADHD-INA_NYC_0010011_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010028/ADHD200_ADHD-INA_NYC_0010028_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010033/ADHD200_ADHD-INA_NYC_0010033_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010035/ADHD200_ADHD-INA_NYC_0010035_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010051/ADHD200_ADHD-INA_NYC_0010051_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010056/ADHD200_ADHD-INA_NYC_0010056_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010071/ADHD200_ADHD-INA_NYC_0010071_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010075/ADHD200_ADHD-INA_NYC_0010075_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010081/ADHD200_ADHD-INA_NYC_0010081_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010087/ADHD200_ADHD-INA_NYC_0010087_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010096/ADHD200_ADHD-INA_NYC_0010096_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010104/ADHD200_ADHD-INA_NYC_0010104_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010105/ADHD200_ADHD-INA_NYC_0010105_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010106/ADHD200_ADHD-INA_NYC_0010106_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010107/ADHD200_ADHD-INA_NYC_0010107_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010108/ADHD200_ADHD-INA_NYC_0010108_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010115/ADHD200_ADHD-INA_NYC_0010115_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010116/ADHD200_ADHD-INA_NYC_0010116_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010129/ADHD200_ADHD-INA_NYC_0010129_T1_SD000000-RS00'
      ... 'raw/ADHD200/ADHD-INA/NYC/1023964/ADHD200_ADHD-INA_NYC_1023964_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/1208795/ADHD200_ADHD-INA_NYC_1208795_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/1471736/ADHD200_ADHD-INA_NYC_1471736_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/1511464/ADHD200_ADHD-INA_NYC_1511464_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/1740607/ADHD200_ADHD-INA_NYC_1740607_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/2773205/ADHD200_ADHD-INA_NYC_2773205_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/2821683/ADHD200_ADHD-INA_NYC_2821683_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/2983819/ADHD200_ADHD-INA_NYC_2983819_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/2996531/ADHD200_ADHD-INA_NYC_2996531_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/3441455/ADHD200_ADHD-INA_NYC_3441455_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/4060823/ADHD200_ADHD-INA_NYC_4060823_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/8009688/ADHD200_ADHD-INA_NYC_8009688_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/9907452/ADHD200_ADHD-INA_NYC_9907452_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010006/ADHD200_HC_NYC_0010006_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010008/ADHD200_HC_NYC_0010008_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010009/ADHD200_HC_NYC_0010009_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010010/ADHD200_HC_NYC_0010010_T1_SD000000-RS00'
      ...'raw/ADHD200/HC/NYC/0010014/ADHD200_HC_NYC_0010014_T1_SD000000-RS00' % light
      'raw/ADHD200/HC/NYC/0010021/ADHD200_HC_NYC_0010021_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010023/ADHD200_HC_NYC_0010023_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010024/ADHD200_HC_NYC_0010024_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010031/ADHD200_HC_NYC_0010031_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010034/ADHD200_HC_NYC_0010034_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010038/ADHD200_HC_NYC_0010038_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010039/ADHD200_HC_NYC_0010039_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010045/ADHD200_HC_NYC_0010045_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010052/ADHD200_HC_NYC_0010052_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010053/ADHD200_HC_NYC_0010053_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010054/ADHD200_HC_NYC_0010054_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010057/ADHD200_HC_NYC_0010057_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010058/ADHD200_HC_NYC_0010058_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010059/ADHD200_HC_NYC_0010059_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010063/ADHD200_HC_NYC_0010063_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010068/ADHD200_HC_NYC_0010068_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010069/ADHD200_HC_NYC_0010069_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010070/ADHD200_HC_NYC_0010070_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010079/ADHD200_HC_NYC_0010079_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010080/ADHD200_HC_NYC_0010080_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010082/ADHD200_HC_NYC_0010082_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010083/ADHD200_HC_NYC_0010083_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010084/ADHD200_HC_NYC_0010084_T1_SD000000-RS00' 
      'raw/ADHD200/HC/NYC/0010088/ADHD200_HC_NYC_0010088_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010089/ADHD200_HC_NYC_0010089_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010093/ADHD200_HC_NYC_0010093_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010097/ADHD200_HC_NYC_0010097_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010099/ADHD200_HC_NYC_0010099_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010100/ADHD200_HC_NYC_0010100_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010110/ADHD200_HC_NYC_0010110_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010112/ADHD200_HC_NYC_0010112_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010113/ADHD200_HC_NYC_0010113_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010114/ADHD200_HC_NYC_0010114_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010120/ADHD200_HC_NYC_0010120_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010121/ADHD200_HC_NYC_0010121_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010122/ADHD200_HC_NYC_0010122_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010124/ADHD200_HC_NYC_0010124_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010125/ADHD200_HC_NYC_0010125_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010128/ADHD200_HC_NYC_0010128_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/1000804/ADHD200_HC_NYC_1000804_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/1127915/ADHD200_HC_NYC_1127915_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/1320247/ADHD200_HC_NYC_1320247_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/1359325/ADHD200_HC_NYC_1359325_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/1435954/ADHD200_HC_NYC_1435954_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/1567356/ADHD200_HC_NYC_1567356_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/1700637/ADHD200_HC_NYC_1700637_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/1737393/ADHD200_HC_NYC_1737393_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/1854959/ADHD200_HC_NYC_1854959_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/1875084/ADHD200_HC_NYC_1875084_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/1884448/ADHD200_HC_NYC_1884448_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/1934623/ADHD200_HC_NYC_1934623_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/1995121/ADHD200_HC_NYC_1995121_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/2136051/ADHD200_HC_NYC_2136051_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/2730704/ADHD200_HC_NYC_2730704_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/2735617/ADHD200_HC_NYC_2735617_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/2907383/ADHD200_HC_NYC_2907383_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/2991307/ADHD200_HC_NYC_2991307_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/3011311/ADHD200_HC_NYC_3011311_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/3163200/ADHD200_HC_NYC_3163200_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/3349423/ADHD200_HC_NYC_3349423_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/3518345/ADHD200_HC_NYC_3518345_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/3650634/ADHD200_HC_NYC_3650634_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/3845761/ADHD200_HC_NYC_3845761_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/3999344/ADHD200_HC_NYC_3999344_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/4079254/ADHD200_HC_NYC_4079254_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/4084645/ADHD200_HC_NYC_4084645_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/4827048/ADHD200_HC_NYC_4827048_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/6206397/ADHD200_HC_NYC_6206397_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/6568351/ADHD200_HC_NYC_6568351_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/8415034/ADHD200_HC_NYC_8415034_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/8692452/ADHD200_HC_NYC_8692452_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/8834383/ADHD200_HC_NYC_8834383_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/9578663/ADHD200_HC_NYC_9578663_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/9750701/ADHD200_HC_NYC_9750701_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021002/ADHD200_NG_NYC_0021002_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021003/ADHD200_NG_NYC_0021003_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021005/ADHD200_NG_NYC_0021005_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021006/ADHD200_NG_NYC_0021006_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021007/ADHD200_NG_NYC_0021007_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021009/ADHD200_NG_NYC_0021009_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021010/ADHD200_NG_NYC_0021010_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021013/ADHD200_NG_NYC_0021013_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021014/ADHD200_NG_NYC_0021014_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021015/ADHD200_NG_NYC_0021015_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021016/ADHD200_NG_NYC_0021016_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021017/ADHD200_NG_NYC_0021017_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021018/ADHD200_NG_NYC_0021018_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021019/ADHD200_NG_NYC_0021019_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021020/ADHD200_NG_NYC_0021020_T1_SD000000-RS00'
      ... 'raw/ADHD200/NG/NYC/0021021/ADHD200_NG_NYC_0021021_T1_SD000000-RS00' % light
      'raw/ADHD200/NG/NYC/0021022/ADHD200_NG_NYC_0021022_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021023/ADHD200_NG_NYC_0021023_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021024/ADHD200_NG_NYC_0021024_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021025/ADHD200_NG_NYC_0021025_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021026/ADHD200_NG_NYC_0021026_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021028/ADHD200_NG_NYC_0021028_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021029/ADHD200_NG_NYC_0021029_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021030/ADHD200_NG_NYC_0021030_T1_SD000000-RS00'
      ... 'raw/ADHD200/NG/NYC/0021031/ADHD200_NG_NYC_0021031_T1_SD000000-RS00' % artifacts
      'raw/ADHD200/NG/NYC/0021032/ADHD200_NG_NYC_0021032_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021033/ADHD200_NG_NYC_0021033_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021034/ADHD200_NG_NYC_0021034_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021035/ADHD200_NG_NYC_0021035_T1_SD000000-RS00'
      ... 'raw/ADHD200/NG/NYC/0021036/ADHD200_NG_NYC_0021036_T1_SD000000-RS00' % light
      'raw/ADHD200/NG/NYC/0021037/ADHD200_NG_NYC_0021037_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021038/ADHD200_NG_NYC_0021038_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021039/ADHD200_NG_NYC_0021039_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021040/ADHD200_NG_NYC_0021040_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021041/ADHD200_NG_NYC_0021041_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021042/ADHD200_NG_NYC_0021042_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021043/ADHD200_NG_NYC_0021043_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021044/ADHD200_NG_NYC_0021044_T1_SD000000-RS00'
    };
  bad.NYC = {
      'raw/ADHD200/ADHD-COM/NYC/0010016/ADHD200_ADHD-COM_NYC_0010016_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010018/ADHD200_ADHD-COM_NYC_0010018_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010020/ADHD200_ADHD-COM_NYC_0010020_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010022/ADHD200_ADHD-COM_NYC_0010022_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010027/ADHD200_ADHD-COM_NYC_0010027_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010032/ADHD200_ADHD-COM_NYC_0010032_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010047/ADHD200_ADHD-COM_NYC_0010047_T1_SD000000-RS00' % light
      'raw/ADHD200/ADHD-COM/NYC/0010074/ADHD200_ADHD-COM_NYC_0010074_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010078/ADHD200_ADHD-COM_NYC_0010078_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/0010098/ADHD200_ADHD-COM_NYC_0010098_T1_SD000000-RS00' % light
      'raw/ADHD200/ADHD-COM/NYC/0010127/ADHD200_ADHD-COM_NYC_0010127_T1_SD000000-RS00' % light
      'raw/ADHD200/ADHD-COM/NYC/1099481/ADHD200_ADHD-COM_NYC_1099481_T1_SD000000-RS00' % light
      'raw/ADHD200/ADHD-COM/NYC/1918630/ADHD200_ADHD-COM_NYC_1918630_T1_SD000000-RS00' % f 
      'raw/ADHD200/ADHD-COM/NYC/2030383/ADHD200_ADHD-COM_NYC_2030383_T1_SD000000-RS00' % light
      ... 'raw/ADHD200/ADHD-COM/NYC/2570769/ADHD200_ADHD-COM_NYC_2570769_T1_SD000000-RS00' % light
      'raw/ADHD200/ADHD-COM/NYC/3174224/ADHD200_ADHD-COM_NYC_3174224_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/3433846/ADHD200_ADHD-COM_NYC_3433846_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/3542588/ADHD200_ADHD-COM_NYC_3542588_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/3619797/ADHD200_ADHD-COM_NYC_3619797_T1_SD000000-RS00' % light
      'raw/ADHD200/ADHD-COM/NYC/3653737/ADHD200_ADHD-COM_NYC_3653737_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/4095229/ADHD200_ADHD-COM_NYC_4095229_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/8697774/ADHD200_ADHD-COM_NYC_8697774_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/NYC/8915162/ADHD200_ADHD-COM_NYC_8915162_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-HIM/NYC/0010005/ADHD200_ADHD-HIM_NYC_0010005_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010015/ADHD200_ADHD-INA_NYC_0010015_T1_SD000000-RS00'
      ... 'raw/ADHD200/ADHD-INA/NYC/0010030/ADHD200_ADHD-INA_NYC_0010030_T1_SD000000-RS00' % light
      'raw/ADHD200/ADHD-INA/NYC/0010061/ADHD200_ADHD-INA_NYC_0010061_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010062/ADHD200_ADHD-INA_NYC_0010062_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010067/ADHD200_ADHD-INA_NYC_0010067_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/0010119/ADHD200_ADHD-INA_NYC_0010119_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/NYC/2297413/ADHD200_ADHD-INA_NYC_2297413_T1_SD000000-RS00'
      ... 'raw/ADHD200/ADHD-INA/NYC/2854839/ADHD200_ADHD-INA_NYC_2854839_T1_SD000000-RS00' % bad defaceing
      'raw/ADHD200/ADHD-INA/NYC/3349205/ADHD200_ADHD-INA_NYC_3349205_T1_SD000000-RS00' % light
      'raw/ADHD200/HC/NYC/0010003/ADHD200_HC_NYC_0010003_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010004/ADHD200_HC_NYC_0010004_T1_SD000000-RS00'
      ... 'raw/ADHD200/HC/NYC/0010036/ADHD200_HC_NYC_0010036_T1_SD000000-RS00' % light
      'raw/ADHD200/HC/NYC/0010043/ADHD200_HC_NYC_0010043_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010046/ADHD200_HC_NYC_0010046_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010055/ADHD200_HC_NYC_0010055_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010065/ADHD200_HC_NYC_0010065_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010066/ADHD200_HC_NYC_0010066_T1_SD000000-RS00'
      ... 'raw/ADHD200/HC/NYC/0010076/ADHD200_HC_NYC_0010076_T1_SD000000-RS00' light
      'raw/ADHD200/HC/NYC/0010077/ADHD200_HC_NYC_0010077_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010092/ADHD200_HC_NYC_0010092_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010094/ADHD200_HC_NYC_0010094_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010102/ADHD200_HC_NYC_0010102_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010111/ADHD200_HC_NYC_0010111_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010117/ADHD200_HC_NYC_0010117_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/0010123/ADHD200_HC_NYC_0010123_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/3243657/ADHD200_HC_NYC_3243657_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/3662296/ADHD200_HC_NYC_3662296_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/4164316/ADHD200_HC_NYC_4164316_T1_SD000000-RS00'
      'raw/ADHD200/HC/NYC/4562206/ADHD200_HC_NYC_4562206_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021008/ADHD200_NG_NYC_0021008_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021027/ADHD200_NG_NYC_0021027_T1_SD000000-RS00'
      'raw/ADHD200/NG/NYC/0021046/ADHD200_NG_NYC_0021046_T1_SD000000-RS00' % light
   }; 
  good.ORE = {
     'raw/ADHD200/ADHD-COM/ORE/1084283/ADHD200_ADHD-COM_ORE_1084283_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/ORE/1108916/ADHD200_ADHD-COM_ORE_1108916_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/ORE/1340333/ADHD200_ADHD-COM_ORE_1340333_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/ORE/1411223/ADHD200_ADHD-COM_ORE_1411223_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/ORE/2054310/ADHD200_ADHD-COM_ORE_2054310_T1_SD000000-RS00' % light
      'raw/ADHD200/ADHD-COM/ORE/2071989/ADHD200_ADHD-COM_ORE_2071989_T1_SD000000-RS00' % light
      'raw/ADHD200/ADHD-COM/ORE/2288903/ADHD200_ADHD-COM_ORE_2288903_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/ORE/2571197/ADHD200_ADHD-COM_ORE_2571197_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/ORE/3052540/ADHD200_ADHD-COM_ORE_3052540_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/ORE/3358877/ADHD200_ADHD-COM_ORE_3358877_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/ORE/4016887/ADHD200_ADHD-COM_ORE_4016887_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/ORE/6953386/ADHD200_ADHD-COM_ORE_6953386_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/ORE/1206380/ADHD200_ADHD-INA_ORE_1206380_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/ORE/1552181/ADHD200_ADHD-INA_ORE_1552181_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/ORE/2415970/ADHD200_ADHD-INA_ORE_2415970_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/ORE/3677724/ADHD200_ADHD-INA_ORE_3677724_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/ORE/3684229/ADHD200_ADHD-INA_ORE_3684229_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/ORE/4072305/ADHD200_ADHD-INA_ORE_4072305_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/ORE/5302451/ADHD200_ADHD-INA_ORE_5302451_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/1084884/ADHD200_HC_ORE_1084884_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/1418396/ADHD200_HC_ORE_1418396_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/1421489/ADHD200_HC_ORE_1421489_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/1502229/ADHD200_HC_ORE_1502229_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/1548937/ADHD200_HC_ORE_1548937_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/1647968/ADHD200_HC_ORE_1647968_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/2054998/ADHD200_HC_ORE_2054998_T1_SD000000-RS00' % light
      'raw/ADHD200/HC/ORE/2124248/ADHD200_HC_ORE_2124248_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/2155356/ADHD200_HC_ORE_2155356_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/2232413/ADHD200_HC_ORE_2232413_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/2426523/ADHD200_HC_ORE_2426523_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/2427434/ADHD200_HC_ORE_2427434_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/2578455/ADHD200_HC_ORE_2578455_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/2947936/ADHD200_HC_ORE_2947936_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/2959809/ADHD200_HC_ORE_2959809_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/3048401/ADHD200_HC_ORE_3048401_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/3051944/ADHD200_HC_ORE_3051944_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/3206978/ADHD200_HC_ORE_3206978_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/3212875/ADHD200_HC_ORE_3212875_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/3244985/ADHD200_HC_ORE_3244985_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/3302025/ADHD200_HC_ORE_3302025_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/3812101/ADHD200_HC_ORE_3812101_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/3848511/ADHD200_HC_ORE_3848511_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/3869075/ADHD200_HC_ORE_3869075_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/3899622/ADHD200_HC_ORE_3899622_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/4046678/ADHD200_HC_ORE_4046678_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/4103874/ADHD200_HC_ORE_4103874_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/4529116/ADHD200_HC_ORE_4529116_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/6592761/ADHD200_HC_ORE_6592761_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/8218392/ADHD200_HC_ORE_8218392_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/9499804/ADHD200_HC_ORE_9499804_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023000/ADHD200_NG_ORE_0023000_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023001/ADHD200_NG_ORE_0023001_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023003/ADHD200_NG_ORE_0023003_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023004/ADHD200_NG_ORE_0023004_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023005/ADHD200_NG_ORE_0023005_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023006/ADHD200_NG_ORE_0023006_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023007/ADHD200_NG_ORE_0023007_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023008/ADHD200_NG_ORE_0023008_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023010/ADHD200_NG_ORE_0023010_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023012/ADHD200_NG_ORE_0023012_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023013/ADHD200_NG_ORE_0023013_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023016/ADHD200_NG_ORE_0023016_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023017/ADHD200_NG_ORE_0023017_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023018/ADHD200_NG_ORE_0023018_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023019/ADHD200_NG_ORE_0023019_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023020/ADHD200_NG_ORE_0023020_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023024/ADHD200_NG_ORE_0023024_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023025/ADHD200_NG_ORE_0023025_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023026/ADHD200_NG_ORE_0023026_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023027/ADHD200_NG_ORE_0023027_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023028/ADHD200_NG_ORE_0023028_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023030/ADHD200_NG_ORE_0023030_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023033/ADHD200_NG_ORE_0023033_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023035/ADHD200_NG_ORE_0023035_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023036/ADHD200_NG_ORE_0023036_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023037/ADHD200_NG_ORE_0023037_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023038/ADHD200_NG_ORE_0023038_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023039/ADHD200_NG_ORE_0023039_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023040/ADHD200_NG_ORE_0023040_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023041/ADHD200_NG_ORE_0023041_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023042/ADHD200_NG_ORE_0023042_T1_SD000000-RS00'
  };
  bad.ORE = {
      'raw/ADHD200/ADHD-COM/ORE/1536593/ADHD200_ADHD-COM_ORE_1536593_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/ORE/1743472/ADHD200_ADHD-COM_ORE_1743472_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/ORE/2292940/ADHD200_ADHD-COM_ORE_2292940_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/ORE/2455205/ADHD200_ADHD-COM_ORE_2455205_T1_SD000000-RS00' % light-med
      'raw/ADHD200/ADHD-COM/ORE/2535204/ADHD200_ADHD-COM_ORE_2535204_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/ORE/2561174/ADHD200_ADHD-COM_ORE_2561174_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/ORE/2620872/ADHD200_ADHD-COM_ORE_2620872_T1_SD000000-RS00' % f
      'raw/ADHD200/ADHD-COM/ORE/3466651/ADHD200_ADHD-COM_ORE_3466651_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/ORE/3470141/ADHD200_ADHD-COM_ORE_3470141_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/ORE/3560456/ADHD200_ADHD-COM_ORE_3560456_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-COM/ORE/7333005/ADHD200_ADHD-COM_ORE_7333005_T1_SD000000-RS00' % light-med
      'raw/ADHD200/ADHD-HIM/ORE/2929195/ADHD200_ADHD-HIM_ORE_2929195_T1_SD000000-RS00' % f
      'raw/ADHD200/ADHD-HIM/ORE/3286474/ADHD200_ADHD-HIM_ORE_3286474_T1_SD000000-RS00' % f
      'raw/ADHD200/ADHD-INA/ORE/2559559/ADHD200_ADHD-INA_ORE_2559559_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/ORE/2790141/ADHD200_ADHD-INA_ORE_2790141_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/ORE/2845989/ADHD200_ADHD-INA_ORE_2845989_T1_SD000000-RS00'
      'raw/ADHD200/ADHD-INA/ORE/3652932/ADHD200_ADHD-INA_ORE_3652932_T1_SD000000-RS00' % f
      'raw/ADHD200/ADHD-INA/ORE/8720244/ADHD200_ADHD-INA_ORE_8720244_T1_SD000000-RS00' % f
      'raw/ADHD200/HC/ORE/1386056/ADHD200_HC_ORE_1386056_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/1481430/ADHD200_HC_ORE_1481430_T1_SD000000-RS00' % light
      'raw/ADHD200/HC/ORE/1664335/ADHD200_HC_ORE_1664335_T1_SD000000-RS00' % light
      'raw/ADHD200/HC/ORE/1679142/ADHD200_HC_ORE_1679142_T1_SD000000-RS00' % light
      'raw/ADHD200/HC/ORE/1696588/ADHD200_HC_ORE_1696588_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/2232376/ADHD200_HC_ORE_2232376_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/2409220/ADHD200_HC_ORE_2409220_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/2920716/ADHD200_HC_ORE_2920716_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/3162671/ADHD200_HC_ORE_3162671_T1_SD000000-RS00' % f light
      'raw/ADHD200/HC/ORE/4219416/ADHD200_HC_ORE_4219416_T1_SD000000-RS00'
      'raw/ADHD200/HC/ORE/8064456/ADHD200_HC_ORE_8064456_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023002/ADHD200_NG_ORE_0023002_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023011/ADHD200_NG_ORE_0023011_T1_SD000000-RS00'
      'raw/ADHD200/NG/ORE/0023031/ADHD200_NG_ORE_0023031_T1_SD000000-RS00'
  };
end


 

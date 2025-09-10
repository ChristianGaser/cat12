function cat_tst_qa_IXI( datadir0 , qaversions , segment, fasttest, rerun, dataset ) 
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
  if ~exist('dataset','var'), dataset = 'IXI'; end

  % ### datadir ###
  if ~exist( 'datadir0' , 'var' )
    switch dataset 
      case 'IXI',        datadir   = '/Volumes/SG5TB/MRData/202503_QA/IXI'; 
      case 'ADHD200',    datadir   = '/Volumes/WDE18TB/MRDataPP/202503_QA/ADHD200'; 
      case 'ABIDE',      datadir   = '/Volumes/WDE18TB/MRDataPP/202105_QA/private/Site O - ABIDE'; 
      case 'NKIrs',      datadir   = '/Volumes/WDE18TB/MRDataPP/202105_QA/private/Site G - NKI';
      case 'NKI',        datadir   = '/Volumes/WDE18TB/MRDataPP/202105_QA/NKI';
      case 'ADHD200NYC', datadir   = '/Volumes/WDE18TB/MRDataPP/202105_QA/private/Site H - ADHD_NYC';
      case 'ADHD200ORE', datadir   = '/Volumes/WDE18TB/MRDataPP/202105_QA/private/Site H - ADHD_ORE';
      case 'ds000256',   datadir   = '/Volumes/SG5TB/MRData/202503_QA/ds000256';
    end
  else
    switch dataset
      case 'IXI',        datadir   = fullfile(datadir0,'IXI-T1');
      case 'ADHD200',    datadir   = fullfile(datadir0,'ADHD200'); 
      case 'ABIDE',      datadir   = '/Volumes/WDE18TB/MRDataPP/202105_QA/private/Site O - ABIDE'; 
      case 'NKIrs',      datadir   = '/Volumes/WDE18TB/MRDataPP/202105_QA/private/Site G - NKI';
      case 'NKI',        datadir   = '/Volumes/WDE18TB/MRDataPP/202105_QA/NKI';
      case 'ADHD200NYC', datadir   = '/Volumes/WDE18TB/MRDataPP/202105_QA/private/Site H - ADHD_NYC';
      case 'ADHD200ORE', datadir   = '/Volumes/WDE18TB/MRDataPP/202105_QA/private/Site I - ADHD_ORE';
      case 'ds000256',   datadir   = fullfile(datadir0,'ds000256-ChildrenHeadMotionN24'); 
    end
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

  runPP      = 0; % run CAT/SPM preprocessing 
  useratings = 1; 
  printall   = 1; 

  resdir  = fullfile(fileparts(datadir), '+results',[dataset '_' fast{fasttest+1} '_202508']); %' datestr(clock,'YYYYmm')]); 
  if ~exist(resdir,'dir'), mkdir(resdir); end
  
  qias = 1:numel(qaversions);
  
  if runPP
    for si = 1:numel(segment)
      clear matlabbatch; 
      switch dataset
        case 'IXI',        IXIfiles = cat_vol_findfiles( datadir , 'IXI*.nii',struct('depth',0));
        case 'ADHD200',    IXIfiles = [
                              cat_vol_findfiles( fullfile( datadir , 'train_raw'), '*.nii.gz',struct('depth',1));
                              ...cat_vol_findfiles( fullfile( datadir , 'test_raw'),  '*.nii.gz',struct('depth',1));
                            ];
        case 'ABIDE',      IXIfiles = cat_vol_findfiles( datadir , 'anat.nii',struct('depth',5));
        case 'NKI',        IXIfiles = cat_vol_findfiles( datadir , 'NKI*.nii',struct('depth',0));
        case {'ADHD200NYC','ADHD200ORE'}
                           IXIfiles = cat_vol_findfiles( datadir , 'ADHD*.nii',struct('depth',0));
        case 'ds000256',   IXIfiles = cat_vol_findfiles( datadir , 'sub*T1w.nii.gz',struct('depth',3)); 
      end
      switch segment{si}
        case 'CAT'
          CATpreprocessing4qc;
          IXIfilesCAT = IXIfiles;
          switch dataset
            case {'IXI'}
              IXIfilesCAT( cellfun(@(x) exist(x,'file'),...
                spm_file(IXIfilesCAT,'prefix',['mri' filesep 'p0']))>0 ) = [];
            case {'ABIDE','ADHD200NYC','ADHD200ORE','NKI','NKIrs'}
              IXIfilesCAT( cellfun(@(x) exist(x,'file'),...
                spm_file(IXIfilesCAT,'prefix',['mri' filesep 'p0']))>0 ) = [];
            case 'ds000256'
              matlabbatch{1}.spm.tools.cat.estwrite.output.BIDS.BIDSyes.BIDSdir = ...
                fullfile('..','derivatives','CAT12.9_2874');
              IXIfilesCAT2 = strrep( IXIfilesCAT , datadir , fullfile( datadir, 'derivatives','CAT12.9_2874') ); 
              IXIfilesCAT( cellfun(@(x) exist(x,'file'), ...
                spm_file(IXIfilesCAT2, 'prefix','p0','ext','') )>0 ) = [];   % remove gz       
          end
          if ~isempty( IXIfilesCAT )
            matlabbatch{1}.spm.tools.cat.estwrite.data = IXIfilesCAT; 
            matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.lazy = 1; 
            spm_jobman('run',matlabbatch);
          end
        case 'SPM'
          switch dataset
            case 'IXI'
              SPMpreprocessing4qc;
              IXIfilesSPM = IXIfiles; 
              IXIfilesSPM( cellfun(@(x) exist(x,'file'),spm_file(IXIfilesSPM,'prefix','c1'))>0 ) = [];
              if ~isempty( IXIfilesSPM )
                matlabbatch{1}.spm.spatial.preproc.channel.vols = IXIfilesSPM;
                spm_jobman('run',matlabbatch);
              end
            otherwise
              error('%s is not prepared',dataset); 
          end
        case 'synthseg'
          error('synthseg is not prepared in the public script ... use MATLAB help')
        case 'qcseg'
          fprintf('No preprocessing required.\n\n');
      end
    end
  end
  
  
  
  
  %% (re)process QC values
  fprintf('Prepare %s: \n',dataset)
  for si = 1:numel(segment)
    for qai = qias
      switch segment{si}
        case 'CAT'
          switch dataset
            case {'IXI','ADHD200NYC','ADHD200ORE','NKI','NKIrs'}
              Pp0{si}{qai} = cat_vol_findfiles( fullfile( datadir , 'mri') , 'p0*',struct('depth',0));
            case 'ADHD200'
              Pp0{si}{qai} = cat_vol_findfiles( datadir , 'p0*',struct('depth',3));
            case 'ABIDE'
              Pp0{si}{qai} = cat_vol_findfiles( datadir , 'p0anat.nii',struct('depth',6));
            case 'ds000256'
              Pp0{si}{qai} = cat_vol_findfiles( datadir , 'p0*T1w.nii',struct('depth',5));
          end  
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
  fprintf('Prepare %s done. \n',dataset)
  
  %%
  fprintf('Process %s: \n',dataset)
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
  fprintf('Process %s done. \n',dataset)
  %%
  fprintf('Print %s: \n',dataset)
  switch dataset
    case 'IXI'
      Pcsv = fullfile( fileparts(datadir) , 'IXI_2019main.csv' ); 
      csv  = cat_io_csv(Pcsv,'','',struct('delimiter',';','csv','.')); 
    case 'ABIDE1'
      Pcsv = fullfile( datadir , '+META', 'Phenotypic_V1_0b.csv');
      csv  = cat_io_csv(Pcsv,'','',struct('delimiter',',','csv','.')); 
      csv(end,:) = []; 
    case 'ABIDE'
      Pcsv = fullfile( datadir , '+META', 'ABIDEII_Composite_Phenotypic2.csv' );
      csv  = cat_io_csv(Pcsv,'','',struct('delimiter',',','csv','.')); 
      csv(end,:) = []; 
    case 'ADHD200'
      Pcsv = fullfile( datadir , 'phenotypic', 'tables', 'ADHD200c.csv' );
      csv  = cat_io_csv(Pcsv,'','',struct('delimiter',',','csv','.')); 
    case 'ds000256'
      Pcsv = fullfile( datadir , 'participantsQC.tsv' ); 
      csv  = cat_io_csv(Pcsv,'','',struct('delimiter','\t','csv','.')); 
    case 'NKI'
      Pcsv = fullfile( datadir , 'participants.tsv' ); 
      csv  = cat_io_csv(Pcsv,'','',struct('delimiter','\t','csv','.')); 
    case {'ADHD200NYC','ADHD200ORE','NKIrs'}
      clear csv; 
      Pcsv = spm_file(cat_io_strrep(Pp0{1}{1},fullfile('mri','p0'),'vbmdb_'),'ext','.csv'); 
      for csvi = 1:numel(Pcsv)
        csvtmp = cat_io_csv(Pcsv{csvi},'','',struct('delimiter',',','csv','.')); 
        if ~isempty(csvtmp)
          if ~exist('csv','var'), csv = csvtmp; else, csv = [csv; csvtmp(2,:)]; end
        else
          if ~exist('csv','var'), csv = csvtmp; else, csv = [csv; repmat({''},1,size(csv,2))]; end
        end
      end
  end
  
  %% create figure
  fh = figure(39);
  fh.Visible = 'off';
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
              Pxml{pi,1} = fullfile(strrep(pp,[filesep 'mri'],[filesep 'report']),[qaversions{qai} '_' ff(3:end) '.xml']);
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
    
      clear NCR IQR ICR SIQR ECR ECRmm FEC CON TIV rGMV ID csvID age sex qc dx; ni = -1; 
      for xi = 1:numel(xml)
        %
        try
          CON(xi)   = xml(xi).qualityratings.contrastr;
          TIV(xi)   = xml(xi).subjectmeasures.vol_TIV; 
          rGMV(xi)  = xml(xi).subjectmeasures.vol_rel_CGW(2);  
        catch
          CON(xi)   = nan; 
          NCR(xi)   = nan; 
          ICR(xi)   = nan; 
          IQR(xi)   = nan;
          FEC(xi)   = nan; 
          SIQR(xi)  = nan; 
          TIV(xi)   = nan; 
          rGMV(xi)  = nan; 
          res_RMS(xi) = nan;
          siten{xi} = ''; 
          age(xi)   = 0; 
          sex(xi)   = 0; 
          continue
        end
        NCR(xi)   = xml(xi).qualityratings.NCR;
        ICR(xi)   = xml(xi).qualityratings.ICR;
        IQR(xi)   = xml(xi).qualityratings.IQR;
        res_RMS(xi)   = xml(xi).qualityratings.res_RMS;
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
        % get ID of the XML entry
      
        
        % get age of the subject with the ID 
        fname{xi} = xml(xi).filedata.file; %#ok<*SAGROW>
        siten{xi} = ''; 
        try
          switch dataset
            case 'IXI'
              seps      = find( fname{xi} == '-' );
              ID(xi)    = str2double( fname{xi}(4:6));
              csvID(xi) = find( cell2mat(csv(2:end,1)) == ID(xi) , 1 , 'first' ) + 1; % not all exist
              age(xi)   = csv{ csvID(xi) , end }; 
              sex(xi)   = csv{ csvID(xi) , 2 }==1;
              siten{xi} = fname{xi}(seps(1)+1:seps(2)-1);
            case 'ABIDE1'
              ID{xi}    = spm_file( fileparts( fileparts( fileparts(xml(xi).filedata.fname) )) , 'basename');
              csvID(xi) = find( contains( cellstr(num2str( ([csv{2:end,2}])')), ID{xi})==1 , 1 , 'first' ) + 1;  %([csv{2:end,2}])'
              age(xi)   = csv{ csvID(xi) , 5 }; 
              sex(xi)   = csv{ csvID(xi) , 6 }==1;
              siten{xi} = csv{ csvID(xi) , 1 } ; 
            case 'ABIDE' % ABIDE2
              ID{xi}    = spm_file( fileparts( fileparts( fileparts(xml(xi).filedata.fname) )) , 'basename');
              csvID(xi) = find( contains( cellstr(num2str( ([csv{2:end,2}])')), ID{xi})==1 , 1 , 'first' ) + 1;  %([csv{2:end,2}])'
              age(xi)   = csv{ csvID(xi) , 7 }; 
              sex(xi)   = csv{ csvID(xi) , 8 }==1;
              siten{xi} = csv{ csvID(xi) , 1 } ; 
            case 'ADHD200'
              ID{xi}    = spm_file( spm_file( xml(xi).filedata.fname, 'basename'),'ext','');
              csvID(xi) = find( cat_io_contains( cellfun(@num2str,csv(2:end,1),'UniformOutput',false) , ...
                num2str(str2double(ID{xi})) ) ==1 , 1 , 'first' ) + 1; 
              if ~isempty(csvID(xi))
                age(xi)   = csv{ csvID(xi) , 4 }; 
                sex(xi)   = csv{ csvID(xi) , 3 }==1;
                siten{xi} = num2str(csv{ csvID(xi) , 2 }); 
              end
              dx(xi)    = isempty(csv{ csvID(xi) , 6 }) | (csv{ csvID(xi), 6 }>0); 
              if isnumeric(csv{ csvID(xi) , 22 }) 
                qc(xi)  = csv{ csvID(xi) , 22 };
              else
                qc(xi)  = csv{ csvID(xi) , 18 };
              end
            case 'NKI'
              ID{xi}    = cat_io_strrep(fname{xi},{'_T1w','.nii','.gz','sub-'},{'','','',''});
              csvID(xi) = find( cat_io_contains( csv(2:end,1), ID{xi}(1:9)) ==1 , 1 , 'first' ) + 1; 
              age(xi)   = csv{ csvID(xi) , 2 }; 
              sex(xi)   = strcmp( csv{ csvID(xi) , 3 }, 'MALE');
            case {'ADHD200NYC','ADHD200ORE','NKIrs'}
              ID{xi}    = cat_io_strrep(fname{xi},{'_T1w','.nii','.gz'},{'','',''});
              csvID(xi) = find( cat_io_contains( spm_file(csv(2:end,1),'path','','ext',''), ID(xi)) , 1 , 'first' ) + 1; 
              age(xi)   = csv{ csvID(xi) , 6 }; 
              sex(xi)   = csv{ csvID(xi) , 11 }=='M';
            case 'ds000256' 
              ID{xi}    = cat_io_strrep(fname{xi},{'_T1w','.nii','.gz'},{'','',''});
              csvID(xi) = find( contains( csv(2:end,1), ID{xi} )==1, 1 , 'first' ) + 1; 
              age(xi)   = csv{ csvID(xi) , 2 }; 
              sex(xi)   = csv{ csvID(xi) , 3 }=='M';
              qc(xi)    = csv{ csvID(xi) , 5 };
          end
        catch
          csvID(xi) = ni; 
          ni        = ni - 1; 
          age(xi)   = 0;
          sex(xi)   = 0; 
          siten{xi} = 'XXX';
        end
      end

      % remove posible disease cases
      if exist('dx','var')
        dx = dx | age<=0 | contains(siten{xi},'XXX'); 
      else
        dx = age<=0;
      end

      for mii = 1:numel(measures)
        if exist(measures{mii},'var')
          eval(sprintf('%s(dx) = [];',measures{mii}));  
        end
      end
      res_RMS(dx) = []; 
      csvID( dx ) = []; 
      age( dx )   = []; 
      sex( dx )   = []; 
      siten( dx ) = []; 
      fname( dx ) = []; 
      xml( dx )   = []; 
      if exist('qc','var')
        qc( dx )  = []; 
      end
      %Pxml
      dx( dx )    = []; 

      % ignore cases without age
      age(age<=0) = 0; 
      for mii=1:numel(measures)
        if exist(measures{mii},'var')
          eval(sprintf('%s(%s<=0) = nan;',measures{mii},measures{mii}));  
        end
      end

      [siteu,~,site] = unique(siten);
      % remove undefined sites
      siteu(cellfun(@isempty,siteu)) = [];
      siteu(contains(siteu,'XXX')) = [];
           

      if useratings
        mark2rps = @(mark) min(100,max(0,105 - mark*10));
        ratings  = {'NCR','IQR','ICR','SIQR','ECR','FEC','CON','res_RMS'};
        for ri = 1:numel(ratings)
          if exist(ratings{ri},'var')
            eval(sprintf('%s = mark2rps(%s);',ratings{ri},ratings{ri})); 
          end
        end
        lgpos = 'South';
      else
        lgpos = 'North'; 
      end
  
      %% basic overview over dataset
      datasettable = cell(10,numel(siteu)+1);
      datasettable(1:5,1:2) = ...
        {'Datset:', dataset; 
         'N:'     , numel(age);
         'Sites'  , numel(siteu);
         'Age:'   , sprintf('%0.2f%s%0.2f',cat_stat_nanmean(age),char(177),cat_stat_nanstd(age)); 
         'Males:' , sum(sex==1)/numel(sex)};

      datasettable(7:7+3,1) = {'Site:';'N:';'Age';'Males'};  
      for sii = 1:numel(siteu)
        datasettable{7,sii+1}  = sprintf('%d) %s',sii,siteu{sii});
        datasettable{8,sii+1}  = sum(site==sii);
        datasettable{9,sii+1}  = sprintf('%0.2f%s%0.2f',cat_stat_nanmean(age(site==sii)),char(177),cat_stat_nanstd(age(site==sii))); 
        datasettable{10,sii+1} = sum(sex(site==sii)==1) / sum(site==sii);
        for ri = 1:numel(ratings)
          if exist(ratings{ri},'var')
            datasettable{10+ri*2-1,1} = sprintf('mean(%s):',ratings{ri});
            eval(sprintf(' datasettable{10+ri*2-1,sii+1} = cat_stat_nanmean(%s(site==sii));',ratings{ri}));
            datasettable{10+ri*2,1} = sprintf('std(%s):',ratings{ri});
            eval(sprintf(' datasettable{10+ri*2,sii+1} = cat_stat_nanstd(%s(site==sii));',ratings{ri}));
          end
        end
      end
      pfname = fullfile(resdir,sprintf('%s_dataset.csv',dataset)); 
      cat_io_csv(pfname,datasettable); 
      cat_io_cprintf('blue','    Print %s\n',pfname); 


%%
      usenormmain = numel(siteu)>1; 
      for usenorm = 0:usenormmain
        for ri = find(contains(measures,ratings)) 
          if usenorm
            scaling{ri} = [-15 45];
          else
            scaling{ri} = [40 100];
          end
        end

        if usenorm>0
          switch dataset
            case {'ABIDE', 'ADHD200', 'IXI'}
              for sdi = 1:numel(siteu)
                M = contains(siten,siteu{sdi});   
                if sum(M)>0
                  for fni=1:numel(ratings) 
                    if ~strcmp( ratings{fni} , 'res_RMS')
                      eval(sprintf(['%s(M) = -cat_tst_qa_normer( %s(M), ' ...
                        'struct(''figure'',0,''sites'',res_RMS(M),''siterf'',%0.0f,''cmodel'',1,''model'',0));'], ...
                        ratings{fni}, ratings{fni}, 10 + 100*contains('ABIDE',dataset) )); 
                    end
                  end
                end
              end
            otherwise 
              usenorm = 0; 
          end
        end      
  
    %
        avar  = {'age','TIV'}; % analysis variables (x-axis)
        avart = {'in years; '; 'in ml; '};
        switch dataset
          case 'IXI'  
            avars = {[20 90];[900 1900]};
            acols = [[0.0 0.4 0];   [0.2 0.6 0.2]; [0.4 .8 0.4];]; % aging 
            gcols = {[0.5 0 0];     [0 0 0.5 ];    [.5 .5 0]}; % group
            scols = [[0.0 0.4 0.8]; [1 0.5 0.0];   [0.8 .0 0.1];];
          otherwise
            avars = ...
              {[floor(min(age(age>0 & age<prctile(age,98)))/2)*2, ceil(max(age(age>0 & age<prctile(age,98)))/2)*2]; ...
               [floor(min(TIV(TIV>0 & TIV<prctile(TIV,98)))/100)*100, ceil(max(TIV(TIV>0 & TIV<prctile(TIV,98)))/100)*100]};
            acols = [[0.0 0.4 0];   [0.2 0.6 0.2]; [0.4 .8 0.4];]; % aging 
            gcols = {[0.5 0 0];     [0 0 0.5 ];    [.5 .5 0]}; % group
            scols = lines(numel(siteu)); %cat_io_colormaps('set1',numel(siteu)); %sites
        end
     
        %
        tables.agecorr   = {'measure','r1','r2','r3','ra','mean(ra)', 'p1','p2','p3','pa','mean(pa)', 'a1','a2','a3','aa','mean(aa)' }; 
        tables.tivcorr   = {'measure','r1','r2','r3','ra','mean(ra)', 'p1','p2','p3','pa','mean(pa)', 'a1','a2','a3','aa','mean(aa)' }; 
        tables.ageanova  = {'measure','F','Prob>F'}; 
        tables.siteanova = {'measure','F','Prob>F'}; 
        tables.tivanova  = {'measure','F','Prob>F'}; 
        tables.sexanova  = {'measure','F','Prob>F'}; 
        tables.sexmwwt   = {'measure','p','Z','U'}; 
        tables.sidestd   = {'measure','std(site1)','std(site2)','std(site3)'};
       %%
        for mi = 1:numel(measures)
          if usenorm && any( contains({'rGMV','rWMV','rCSFV','res_RMS'},measures{mi}) )
            continue; 
          end

          if exist(measures{mi},'var') && eval(sprintf('any(~isnan(%s))',measures{mi})) 
            for ai = 1:numel(avar)
      
              if ~strcmp( avar{ai} , measures{mi} )
                %% first subfigure with scatter plot of all values colored by density 
                clf
                marker = {'o'  's'  'd'  'v'  '^'  '>'  '<'  'p'  'h'  }; 
                subplot('position',[0.08 0.16 .25 .75]);
                hold on;
                switch dataset
                  case 'IXI'
                    eval(sprintf('scatter( %s(site==1)'' , %s(site==1)'' ,20, ''o'')', avar{ai},measures{mi}));
                    eval(sprintf('scatter( %s(site==2)'' , %s(site==2)'' ,20, ''^'')', avar{ai},measures{mi}));
                    eval(sprintf('scatter( %s(site==3)'' , %s(site==3)'' ,20, ''v'')', avar{ai},measures{mi}));
                  case {'ADHD200','ABIDE'}
                    for sitei = 1:numel(siteu)
                      %%
                      evalc( sprintf('sch = scatter( %s(site==%d & ~isnan(%s''))'' , %s(site==%d & ~isnan(%s''))'' ,20, ''o'')', ...
                        avar{ai}, sitei,  avar{ai}, measures{mi}, sitei, avar{ai}) );
                      sch.Marker = marker{mod(sitei,numel(marker)-1)+1};
                    end
                  otherwise
                    if exist('qc','var')
                      eval(sprintf('scatter( %s(qc==0)'' , %s(qc==0)'' ,20, ''o'')', avar{ai},measures{mi}));
                      eval(sprintf('scatter( %s(qc==1)'' , %s(qc==1)'' ,20, ''s'')', avar{ai},measures{mi}));
                      eval(sprintf('scatter( %s(qc==2)'' , %s(qc==2)'' ,20, ''s'')', avar{ai},measures{mi}));
                    else
                      eval(sprintf('scatter( %s'' , %s'' ,20, ''o'')', avar{ai},measures{mi}));
                    end
                end
                ax  = gca; 
                switch dataset
                  case {'ADHD200','ABIDE'}
                    ax.Position = [0.08 0.16 .27 .75];
                  otherwise
                    ax.Position = [0.08 0.16 .40 .75];
                end
                ax.FontSize = 10; 
                if usenorm
                  ax.YDir = 'reverse';
                end
                %
                hsd = findobj( get(ax,'Children') ,'Type','Scatter'); 
                for hsdi = 1:numel(hsd)
                  switch dataset
                    case 'ds000256'
                      scols = [0 0.8 0; 0.8 0.7 0; 0.8 0 0]; 
                      hsd(hsdi).MarkerEdgeColor = scols(hsdi,:);
                    otherwise
                      hsd(hsdi).MarkerEdgeColor = scols(numel(hsd) + 1 - hsdi,:);
                  end
                  hsd(hsdi).MarkerFaceColor = hsd(hsdi).MarkerEdgeColor; 
                  hsd(hsdi).MarkerFaceAlpha = 0.2; 
                  hsd(hsdi).MarkerEdgeAlpha = 0.4; 
                end
                box on; 
                title(sprintf('%s density in aging (%s, N=%d)',dataset,measures{mi},numel(NCR))); 
                myylim = scaling{mi}; % round( ax.YLim .* [0.9 1.1] * 4 ) / 4 ;
                xlim(avars{ai}); 
                % fitting biased by outiers ... this is pearson ...
                switch dataset
                  case 'IXI'
                    eval(sprintf('[R1,P1] = corrcoef( %s(site==1) , %s(site==1));',avar{ai},measures{mi}));  
                    eval(sprintf('[R2,P2] = corrcoef( %s(site==2) , %s(site==2));',avar{ai},measures{mi}));  
                    eval(sprintf('[R3,P3] = corrcoef( %s(site==3) , %s(site==3));',avar{ai},measures{mi}));  
                    eval(sprintf('[R4,P4] = corrcoef( %s , %s);',avar{ai},measures{mi}));  
                    P = cat_stat_nanmean([P1(2),P2(2),P3(2)]); 
                    R = cat_stat_nanmean([R1(2),R2(2),R3(2)]); 
                  otherwise
                    eval(sprintf('[R,P] = corrcoef( %s , %s);',avar{ai},measures{mi}));  
                end
                warning off; 
                if mi==1 && ai==1  % licence call
                  try
                    eval(sprintf('[curve1, goodness, output] = fit( %s(site==1)'', double(%s(site==1))'',''poly1'');', ...
                      avar{ai}, measures{mi})); cplot(1) = plot(curve1); cplot(1).Color = gcols{1}; 
                  end
                  pause(1); 
                end
                warning on
                
                
                % correlation 
                warning off; 
                switch dataset
                  case 'IXI'
                    %%
                    eval(sprintf('[curve1, goodness, output] = fit( %s(site(:)==1 & ~isnan(%s(:)))'', double(%s(site(:)==1 & ~isnan(%s(:))))'',''poly1'');', ...
                      avar{ai}, measures{mi}, measures{mi}, measures{mi})); cplot(1) = plot(curve1); cplot(1).Color = scols(1,:); 
                    eval(sprintf('[curve2, goodness, output] = fit( %s(site(:)==2 & ~isnan(%s(:)))'', double(%s(site(:)==2 & ~isnan(%s(:))))'',''poly1'');', ...
                      avar{ai}, measures{mi}, measures{mi}, measures{mi})); cplot(2) = plot(curve2); cplot(2).Color = scols(2,:);
                    eval(sprintf('[curve3, goodness, output] = fit( %s(site(:)==3 & ~isnan(%s(:)))'', double(%s(site(:)==3 & ~isnan(%s(:))))'',''poly1'');', ...
                      avar{ai}, measures{mi}, measures{mi}, measures{mi}));  cplot(3) = plot(curve3); cplot(3).Color = scols(3,:);
                    eval(sprintf('[curve4, goodness, output] = fit( %s(~isnan(%s(:)))'', double(%s(~isnan(%s(:))))'',''poly1'');', ...
                      avar{ai}, measures{mi}, measures{mi}, measures{mi}));  cplot(4) = plot(curve4); cplot(4).Color = [ 0.2 0.2 0.2 ]; cplot(4).LineStyle = '--'; 
                    pan = cat_stat_nanmean([curve1.p1, curve2.p1, curve3.p1]); 
                  otherwise
                    pan = 0; 
                    for sii = 1:numel(siteu)
                      try
                        eval(sprintf('[curves(sii), goodness, output] = fit( %s(site(:)==sii & ~isnan(%s(:)))'', double(%s(site(:)==sii & ~isnan(%s(:))))'',''poly1'');', ...
                          avar{ai}, measures{mi}, measures{mi}, measures{mi})); cplot(sii) = plot(curves(sii)); cplot(sii).Color = scols(sii,:); 
                        pan = pan + curves(sii).p1; 
                      catch
                      
                      end
                    end
                end
                legend off; 
                warning on;
                
                % anova
                try
                  tab = cell2table([num2cell(age'),num2cell(site),num2cell(TIV')],'Variablenames',{'age','site','TIV'},'rownames',cellstr(num2str(csvID','IXI%03d')));
                  eval(sprintf('[an0] = anova(tab,%s'');',measures{mi})); 
                  cat_io_csv(fullfile(resdir,'anovas',sprintf('IXI_anova_%s-age-site_%s_%s.csv', measures{mi}, qaversions{qai}, segment{si})),...
                    [{''} an0.stats.Properties.VariableNames; an0.stats.Properties.RowNames table2cell(an0.stats)]); 
                end
                
  
                % man-wikney  
                switch dataset
                  case 'IXI'
                    tables.sidestd = [ tables.sidestd ;
                      [ measures(mi) , num2cell( eval(sprintf('[std(%s(site==1)),std(%s(site==2)),std(%s(site==3))];',measures{mi},measures{mi},measures{mi})))] ]; 
                    if ai == 1 
                      tables.agecorr = [ tables.agecorr ; 
                        [ measures(mi) num2cell( [R1(2),R2(2),R3(2),R4(2),R,  P1(2),P2(2),P3(2),P4(2),P,  curve1.p1,curve2.p1,curve3.p1,curve4.p1,pan ]) ] ]; 
                    else
                      tables.tivcorr = [ tables.tivcorr ; 
                        [ measures(mi) num2cell( [R1(2),R2(2),R3(2),R4(2),R,  P1(2),P2(2),P3(2),P4(2),P,  curve1.p1,curve2.p1,curve3.p1,curve3.p1,pan ]) ] ]; 
                    end
                  otherwise
                end  
                ylim( myylim ); xlim(avars{ai}); 
                set(gca,'Ygrid',1)
                if useratings && any(contains(measures{mi},ratings)) && ~usenorm
                  set(gca,'Ytick',scaling{mi}(1)+5 : 10 : scaling{mi}(2));
                end
                if strcmp(measures{mi},'rGMV'), set(gca,'YTickLabel',num2str(get(gca,'YTick')','%0.2f') ); end
                if isnan( R(min(numel(R),2)) )
                  xlabel(sprintf('%s (%s)',avar{ai}, avart{ai})); 
                else
                  xlabel(sprintf('%s (%s b=%0.3f, \\it{r}\\rm{}=%0.3f, \\it{p}\\rm{}=%0.03f)',avar{ai}, avart{ai}, pan, R(min(numel(R),2)), P(min(numel(P),2)))); 
                end
                ylabel(measures{mi}); grid on; hold on; 
                switch dataset
                  case 'IXI'
                    legend({
                      sprintf('Guys (b=%0.3f, \\it{r}\\rm{}=%0.3f, \\it{p}\\rm{}=%0.3f)',curve1.p1, R1(2), P1(2) ), ...
                      sprintf('HH (b=%0.3f, \\it{r}\\rm{}=%0.3f, \\it{p}\\rm{}=%0.3f)'  ,curve2.p1, R2(2), P2(2)) , ...
                      sprintf('IOP (b=%0.3f, \\it{r}\\rm{}=%0.3f, \\it{p}\\rm{}=%0.3f)' ,curve3.p1, R3(2), P3(2)) , ...
                      },'Location',lgpos)
                  otherwise
                end
                if usenorm, ylabel(['N' get(get(gca,'Ylabel'),'String')]); end
                
      
                % second figure with poxplot of age groups
                switch dataset
                  case {'ABIDE','ADHD200'}
                    subplot('position',[0.425 0.16 .16 .75]);
                  otherwise
                    subplot('position',[0.575 0.16 .16 .75]);
                end
                if ai == 1
                  switch dataset
                    case 'IXI'
                      eval(sprintf('[an0,anova1age{mi}] = anova1(%s,1 +(age>40) + (age>60),''off'');',measures{mi})); 
                      bxtxt   = '{ %s(age<40); %s(age>40 & age<60); %s(age>60) }'; 
                      bxnames = {'20-40','40-60','60-90'}; 
                    otherwise
                      km = cat_stat_kmeans(age(age>prctile(age,10) & age<prctile(age,80)),3); 
                      km2 = round([min(age(age>0)), mean(km(1:2)), mean(km(2:3)), max(age)]);
                      bxtxt   = sprintf('{ %%s(age<%0.2f); %%s(age>%0.2f & age<%0.2f); %%s(age>%0.2f) }',km2(2),km(2:3),km(3)); 
                      eval(sprintf('[an0,anova1age{mi}] = anova1(%s,1 +(age>%d) + (age>%d),''off'');',measures{mi},km(2:3))); 
                      bxnames = {sprintf('<%d',km2(2)),sprintf('%d-%d',km2(2:3)),sprintf('>%d',km2(3))}; 
                  end
                  tables.ageanova = [ tables.siteanova; [ measures(mi) anova1age{mi}(2,4) anova1age{mi}(2,5) ] ]; 
                  cat_plot_boxplot( eval(sprintf(bxtxt,...
                    measures{mi},measures{mi},measures{mi})) , ...
                    struct('names',{bxnames},'style',4,'ylim',myylim,'ygrid',0, ...
                    'datasymbol','o','usescatter',1,'groupcolor',acols) );
                  title('age groups');  
                  xlabel('range (years)'); ylabel(measures{mi});
                else
                  if max(age)>60
                    eval(sprintf('[an0,anova1tiv{mi}] = anova1(%s,1 +(age>1300) + (age>1500),''off'');',measures{mi})); 
                    bxtxt   = '{ %s(TIV<1300); %s(TIV>1300 & TIV<1500); %s(TIV>1500) }'; 
                    bxnames = {'<1.3','~1.4','>1.5'}; 
                  else
                    km  = cat_stat_kmeans(TIV(TIV>prctile(TIV,10) & TIV<prctile(TIV,80)),3); 
                    km2 = round([prctile(TIV(TIV>0),5), mean(km(1:2)), mean(km(2:3)), prctile(TIV,95)]/100)*100;
                    bxtxt   = sprintf('{ %%s(TIV<%0.1f); %%s(TIV>%0.1f & TIV<%0.1f); %%s(TIV>%0.1f) }',km2(2),km(2:3),km(3)); 
                    eval(sprintf('[an0,anova1tiv{mi}] = anova1(%s,1 +(TIV>%d) + (TIV>%d),''off'');',measures{mi},km(2:3)));  
                    bxnames = {sprintf('<%0.1f',km2(2)/1000),sprintf('%0.1f-%0.1f',km2(2:3)/1000),sprintf('>%0.1f',km2(3)/1000)}; 
                  end
                  tables.tivanova = [ tables.siteanova; [ measures(mi) anova1tiv{mi}(2,4) anova1tiv{mi}(2,5) ] ]; 
                  cat_plot_boxplot( eval(sprintf(bxtxt,...
                    measures{mi},measures{mi},measures{mi})) , ...
                    struct('names',{bxnames},'style',4,'ylim',myylim,'ygrid',0, ...
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
                if usenorm, set(gca,'YDir','reverse');  ylabel(['N' get(get(gca,'Ylabel'),'String')]);  end
          
                % second figure with poxplot of centers
                 switch dataset
                   case {'ABIDE','ADHD200'}
                     subplot('position',[0.66 0.16 .32 .75]);
                   otherwise
                     subplot('position',[0.83 0.16 .16 .75]);
                 end
                if ai == 1 
                  eval(sprintf('[an0,anova1site{mi}] = anova1(%s,site,''off'');',measures{mi})); 
                  tables.siteanova = [ tables.siteanova; [ measures(mi) anova1site{mi}(2,4)' anova1site{mi}(2,5)' ] ]; 
                  switch dataset
                    case 'IXI'
                      cat_plot_boxplot( eval(sprintf('{ %s(site==1) ; %s(site==2); %s(site==3) }',...
                        measures{mi},measures{mi},measures{mi})) , ...
                        struct('names',{{'Guys','HH','IOP'}},'style',4,'ylim',myylim ,'ygrid',0, ... siteu
                        'datasymbol','o','usescatter',1,'groupcolor',scols) ); 
                    otherwise
                      if numel(siteu)>0
                        sitestr = ''; sitenam = {};
                        for sxi=1:numel(siteu)
                          sitestr = [sitestr sprintf('%s(site==%d);',measures{mi},sxi)]; %#ok<*AGROW>
                          sitenam = [sitenam siteu(sxi)];
                        end
                        cat_plot_boxplot( eval( sprintf('{ %s}',sitestr)) , ...
                          struct('style',4,'ylim',myylim ,'ygrid',0,'names',{siteu}, ...
                          'datasymbol','o','usescatter',1,'groupcolor',scols) ); 
                      else
                        cat_plot_boxplot( eval(sprintf('{ %s(qc==0) ; %s(qc==1); %s(qc==2) }',...
                          measures{mi},measures{mi},measures{mi})) , ...
                          struct('names',{{'no MA','lMA','sMA'}},'style',4,'ylim',myylim ,'ygrid',0, ... siteu
                          'datasymbol','o','usescatter',1,'groupcolor',scols) ); 
                      end
                  end  
                  title('scan sites');  
                  xlabel('center'); ylabel(measures{mi});
                else
                  if ~exist('mwwtest','file')
                    warning('Miss mwwtest function (https://github.com/dnafinder/mwwtest/blob/master/mwwtest.m) use anova1.')
                    eval(sprintf('[an0,anova1sex{mi}] = anova1(%s,sex,''off'');',measures{mi})); 
                    tables.sexanova = [ tables.sexanova; [ measures(mi) anova1sex{mi}(2,5:6) ] ]; 
                  else 
                    % https://github.com/dnafinder/mwwtest/blob/master/mwwtest.m
                    txt = evalc( sprintf('mwwt = mwwtest( %s(sex==0 & ~isnan(%s)) , %s(sex==1 & ~isnan(%s)) );',...
                      measures{mi}, measures{mi}, measures{mi}, measures{mi} ));
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
                if usenorm, set(gca,'YDir','reverse');  ylabel(['N' get(get(gca,'Ylabel'),'String')]);  end
                ax  = gca;  ax.FontSize = 10; hold on; 
                if 0 %strcmp(measures{mi},'IQR') || strcmp(measures{mi},'SIQR')
                  ph = plot([0 20],[1.9 1.9]); ph.Color = [1 0 0 ]; ph.LineStyle = '--'; ph.LineWidth = 1.5;     
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
                pfname = fullfile(pdir,sprintf('%s_%s_%s_%s_%s_norm%d.jpg',dataset,avar{ai}, measures{mi}, qaversions{qai}, segment{si}, usenorm)); 
                print(fh,pfname,'-r1200','-djpeg');
                cat_io_cprintf('blue','    Print %s\n',pfname); 
              end
            end
          end
        end
      end
      
      % print final table
      tf = fieldnames(tables); 
      for tfi = 1:numel(tf)
        pfname = fullfile(resdir,sprintf('%s_%s_%s_%s.csv', dataset, tf{tfi}, qaversions{qai}, segment{si}));
        cat_io_csv(pfname,tables.(tf{tfi})); 
        cat_io_cprintf('blue','    Save %s\n',pfname); 
      end
    end
  end
  fprintf('Print %s done. \n',dataset)

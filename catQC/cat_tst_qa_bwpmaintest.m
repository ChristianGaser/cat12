function cat_tst_qa_bwpmaintest( datadir, qaversions, segment, fasttest, rerun ) 
% -- BWP test script --------------------------------------------------
%
%  Requirements: 
%   1. Matlab with statistics and machine learning toolbox (robustfit).
%   2. Download and install SPM and CAT
%   3. Download IXI T1 data from: 
%      http://biomedic.doc.ic.ac.uk/brain-development/downloads/IXI/IXI-T1.tar
%
%   4. Specify in this script: 
%      1) the data directory "datadir" 
%      2) the QC version you would like to test (the file has to exist in the cat directory) 
%      3) the segmentation you would like to use
%
%  See also cat_tst_qa_main.
%  ------------------------------------------------------------------------


%#ok<*SAGROW,*AGROW,*UNRCH>

% full/faster test design with directory name 

cat_io_cprintf([0 0.5 0],'\n\n== Run cat_tst_qa_bwpmaintest ==\n') 

if ~license('test', 'Statistics_Toolbox')
  error('This function requires the "Statistics and Machine Learning Toolbox" of MATLAB.\n')
end

% ### datadir ###
if ~exist( 'datadir' , 'var' )
  opt.maindir = pwd; 
else
  opt.maindir = datadir; 
end

if ~exist( fullfile( opt.maindir, 'BWP') , 'dir')
  error('Cannot find the required "BWP" directory in "%s".', opt.maindir)
end  
if ~exist( fullfile( opt.maindir, 'BWPr') , 'dir')
  error('Cannot find the required "BWP" directory in "%s".', opt.maindir)
end
if ~exist( fullfile( opt.maindir, 'BWPgt') , 'dir')
  error('Cannot find the required "BWPgt" directory in "%s".', opt.maindir)
end


%% ### segmention ###
if ~exist( 'segment' , 'var')
  segment = {'CAT'}; % {'SPM','CAT','qcseg'}; % qcseg requires cat_vol_qa2024012
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
fast            = {'full','fast'}; 
corrtype        = 'Spearman'; 
opt.type        = '-depsc';
opt.res         = '-r300'; 
opt.dpi         = 90; 
opt.closefig    = 0; 
%opt.resdir      = fullfile(opt.maindir,'+results', ...
%  sprintf('%s_%s_%s', 'BWPmain', fast{fasttest+1}, datestr(clock,'YYYYmm')) );
opt.resdir      = fullfile(opt.maindir,'+results', ...
  sprintf('%s_%s_%s', 'BWPmain', fast{fasttest+1}, '202508' ));
recalc.qa       = 2 * rerun; 
recalc.kappa    = 0*2 * rerun;
recalc.runPP    = 1; 
if ~exist(opt.resdir,'dir'),   mkdir(opt.resdir); end

% preprocessing 
if recalc.runPP
  BWPfiles = [
    cat_vol_findfiles( fullfile( opt.maindir, 'BWP'),   'BWP*.nii'  ,struct('depth',0));
    cat_vol_findfiles( fullfile( opt.maindir, 'BWPr' ), 'rBWP*.nii' ,struct('depth',0));
    cat_vol_findfiles( fullfile( opt.maindir, 'BWPr' ), 'irBWP*.nii',struct('depth',0))];

  for si = 1:numel(segment)
    clear matlabbatch; 
    
    switch segment{si}
      case 'CAT'
        CATpreprocessing4qc;
        BWPfilesCAT = BWPfiles; 
        BWPfilesCAT( cellfun(@(x) exist(x,'file'),spm_file(BWPfilesCAT,'prefix',['mri' filesep 'p0']))>0 ) = [];
        if ~isempty( BWPfilesCAT )
          matlabbatch{1}.spm.tools.cat.estwrite.data = BWPfiles;
          matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.lazy = 1; 
          spm_jobman('run',matlabbatch);
        end
      case 'SPM'
        SPMpreprocessing4qc;
        BWPfilesSPM = BWPfiles; 
        BWPfilesSPM( cellfun(@(x) exist(x,'file'),spm_file(BWPfilesSPM,'prefix','c1'))>0 ) = [];
        if ~isempty( BWPfilesSPM )
          matlabbatch{1}.spm.spatial.preproc.channel.vols = BWPfilesSPM;
          spm_jobman('run',matlabbatch);
        end
      case 'synthseg'
        error('synthseg is not prepared in the public script ... use MATLAB help')
      case 'qcseg'
        fprintf('No preprocessing required.\n\n');
    end
  end
end

%
% * prepare figure for SPM comparison under simple conditions
%   - plot noise estimation (1-9%, 20-40 ABC bias, full res) for CAT vs SPM vs. GT
%   - MR-ART for SPM!

qais = 1:numel(qaversions);
qai  = 1; %#ok<NASGU> 


for qai = qais 
  % prepare things and print
  qafile = qaversions{qai}; 
  cat_io_cprintf('b',sprintf('\n%s:\n\n',qafile)); 
  if ~exist(opt.resdir,'dir'), mkdir(opt.resdir); end
  

  %  -- directories ----------------------------------------------------
  % search path of the QA-xml-files
  P.p0gt = fullfile(opt.maindir,'BWPgt','p0phantom_1.0mm_normal_cor.nii');
  P.bwpx  = {
    ...'GT'    fullfile(opt.maindir,'BWP')      '.' [0.0 0.5 0.0] '-'; ...
    'CAT12'    fullfile(opt.maindir,'BWP')      '<' [1.0 0.2 0.0] '-'; ... % UPDATE
    'SPM12'    fullfile(opt.maindir,'BWP')      '^' [0.0 0.8 0.0] '-'; ... % UPDATE
    ...'synthseg' fullfile(opt.maindir,'BWP')      '.' [0.0 0.5 0.0] '-'; ...
    'qcseg'    fullfile(opt.maindir,'BWP')      '.' [0.0 0.5 0.5] '-'; ...
    };

  methodx = [ any(contains(segment,'CAT')) any(contains(segment,'SPM')) any(contains(segment,'qcseg')) ]; 
  P.bwpx  = P.bwpx(methodx,:); 
  method  = 1; 
  
  for si = 1:size(P.bwpx,1)
    P.bwp = P.bwpx(si,:); 
  
    %  -- find p0-files --------------------------------------------------
    rsdir = 'r';
    % find specific BWP segmentations
    if strcmp(P.bwp(1,1),'qcseg') && ~strcmp( qaversions{qai} , 'cat_vol_qa202412' )
      cat_io_cprintf('err','ERROR: qcseg only supported for cat_vol_qa202412!\n\n');
      continue
    end
    for mi = 1 %method
      % find p0 images 
      if strcmp(P.bwp(mi,1),'CAT12')
        P.p0m{mi} = [cat_vol_findfiles(fullfile( P.bwp{mi,2}       ,'mri'),'p0B*0p*00x*00x*00.nii',struct('depth',0));
                     cat_vol_findfiles(fullfile([P.bwp{mi,2} rsdir],'mri'),'p0rB*0p*00x*00x*00.nii',struct('depth',0));
                     cat_vol_findfiles(fullfile([P.bwp{mi,2} rsdir],'mri'),'p0irB*0p*00x*00x*00.nii',struct('depth',0))];
      elseif strcmp(P.bwp(mi,1),'GT')
        P.p0m{mi} = [cat_vol_findfiles(fullfile( P.bwp{mi,2}       ,'mri'),'p0B*0p*00x*00x*00.nii',struct('depth',0));
                     cat_vol_findfiles(fullfile([P.bwp{mi,2} rsdir],'mri'),'p0rB*0p*00x*00x*00.nii',struct('depth',0));
                     cat_vol_findfiles(fullfile([P.bwp{mi,2} rsdir],'mri'),'p0irB*0p*00x*00x*00.nii',struct('depth',0))];
        P.p0m{mi} = repmat( P.p0gt, size(Pp0m{mi},1), size(Pp0m{mi},2)); 
      elseif strcmp(P.bwp(mi,1),'synthseg')
        P.p0m{mi} = [cat_vol_findfiles( P.bwp{mi,2}       ,'synthseg_p0B*0p*00x*00x*00.nii',struct('depth',0));
                     cat_vol_findfiles([P.bwp{mi,2} rsdir],'synthseg_p0rB*0p*00x*00x*00.nii',struct('depth',0));
                     cat_vol_findfiles([P.bwp{mi,2} rsdir],'synthseg_p0irB*0p*00x*00x*00.nii',struct('depth',0))];
      elseif strcmp(P.bwp(mi,1),'qcseg')
        P.p0m{mi} = [cat_vol_findfiles( P.bwp{mi,2}       ,'B*0p*00x*00x*00.nii',struct('depth',0));
                     cat_vol_findfiles([P.bwp{mi,2} rsdir],'rB*0p*00x*00x*00.nii',struct('depth',0));
                     cat_vol_findfiles([P.bwp{mi,2} rsdir],'irB*0p*00x*00x*00.nii',struct('depth',0))];
      else
        % SPM12
        P.c1{mi}  = [cat_vol_findfiles(P.bwp{mi,2},'c1B*0p*00x*00x*00.nii',struct('depth',0));
             cat_vol_findfiles([P.bwp{mi,2} rsdir],'c1rB*0p*00x*00x*00.nii',struct('depth',0));
             cat_vol_findfiles([P.bwp{mi,2} rsdir],'c1irB*0p*00x*00x*00.nii',struct('depth',0))]; 
        P.c2{mi}  = [cat_vol_findfiles(P.bwp{mi,2},'c2B*0p*00x*00x*00.nii',struct('depth',0)); 
             cat_vol_findfiles([P.bwp{mi,2} rsdir],'c2rB*0p*00x*00x*00.nii',struct('depth',0));
             cat_vol_findfiles([P.bwp{mi,2} rsdir],'c2irB*0p*00x*00x*00.nii',struct('depth',0))];
        P.c3{mi}  = [cat_vol_findfiles(P.bwp{mi,2},'c3B*0p*00x*00x*00.nii',struct('depth',0)); 
             cat_vol_findfiles([P.bwp{mi,2} rsdir],'c3rB*0p*00x*00x*00.nii',struct('depth',0));
             cat_vol_findfiles([P.bwp{mi,2} rsdir],'c3irB*0p*00x*00x*00.nii',struct('depth',0))];
        P.p0m{mi} = cat_io_strrep(P.c1{mi},'c1','p0'); 
        Pp0e      = false(1,numel(P.p0m{mi}));
        for fi = 1:numel(P.p0m{mi})
          Pp0e(fi) = exist(P.p0m{mi}{fi},'file') && ~cat_io_rerun(P.c1{mi}{fi},P.p0m{mi}{fi},0);
        end
        if any(Pp0e~=1)
          cat_io_cgw2seg(P.c3{mi}(Pp0e~=1),P.c1{mi}(Pp0e~=1),P.c2{mi}(Pp0e~=1));
        end
      end
      
      % 
      P.p0m{mi} = P.p0m{mi}(cellfun('isempty',strfind(P.p0m{mi},'_pn0_')));
      for fi = 1:numel(P.p0m{mi})
        [PPP{mi}{fi},FFF{mi}{fi}] = fileparts(P.p0m{mi}{fi});     
      end
    end
    % only data that was processed by all methods
    PP{1} = PPP{method(1)};
    FF{1} = FFF{method(1)};
    if numel(method)>1
      for mi = setdiff( method , 1)
        [FF{1},fii,fim]  = intersect(FF{1},FFF{mi});
        PP{mi}           = PPP{mi}(fim);
        for mii = setdiff( method(1:mi) , 1)
          FF{mii}        = FFF{mii}(fii); 
          PP{mii}        = PPP{mii}(fii);
        end  
      end
    end
    % remove old p0 files and create a new list that only includes test cases that are available for all methods 
    %P = rmfield(P,'p0');
    for mi = method
      for fi = 1:numel(FF{1})
        P.p0{mi}{fi,1} = fullfile(PP{mi}{fi},[FF{1}{fi} '.nii']);
      end
    end
    
    % remove cases for faster preliminary tests
    if fasttest
      % only field A, only extrem and average cases 
      % (1%, 5%, 9% noise; 20%, 60%, 100% inhomogeneity, 1x1x1, 1x1x2, 2x2x2 mm resolution) 
      removeEntries = {'pB_','pC_','pn3','pn7','rf040','rf080',... 
        'vx100x200x100','vx200x100x100','vx100x200x200','vx200x100x200','vx200x200x100'};
    else
      % remove these cases as they play no practical role
      removeEntries = {'vx100x200x200','vx200x100x200','vx200x200x100'};
    end
    if ~isempty(removeEntries)
      for pi = numel(P.p0{1}):-1:1
        [pp,ff,ee] = spm_fileparts( P.p0{1}{pi} );
        if contains(ff,removeEntries)
          for mi = method, P.p0{mi}(pi) = []; P.p0m{mi}(pi) = []; end
        end
      end
    end
    fprintf('  Found %d files (%d methods).  ',numel(P.p0{method(1)}),numel(method));
    
  
  
  
    % estimate the qa and kappa values for the p0-files
    Q = struct(); P.xml = {}; 
    Q.kappa = zeros(numel(P.p0{mi}),numel(method));
    for mi = method
      % -- do QA ---------------------------------------------------------
      % xml-filenames  
      P.p0a{mi} = P.p0{mi};
      for pi=1:numel(P.p0{mi})
        [pp,ff] = fileparts(P.p0{mi}{pi});
        if strcmp(P.bwp(mi,1),'CAT12')
          P.xml{mi}{pi,1} = fullfile(fileparts(pp),'report',[qafile '_' ff(3:end) '.xml']);
        elseif strcmp(P.bwp(mi,1),'SPM12')
          P.xml{mi}{pi,1} = fullfile(pp,[qafile '_spm_' ff(3:end) '.xml']); %'report',
          P.p0a{mi}{pi,1} = fullfile(pp,['c1' ff(3:end) '.nii']); % for kappa
        elseif strcmp(P.bwp(mi,1),'synthseg')
          P.xml{mi}{pi,1} = fullfile(pp,'report',[qafile '_' strrep(ff,'synthseg_p0','synthseg_') '.xml']);
        elseif strcmp(P.bwp(mi,1),'qcseg')
          P.xml{mi}{pi,1} = fullfile(pp,'report',[qafile '_qcseg_' ff '.xml']);
          P.p0a{mi}{pi,1} = fullfile(pp,['p0_qcseg_' ff '.nii']); % for kappa
        end
        P.p0e{mi}(pi) = exist(P.xml{mi}{pi},'file') | strcmp(P.bwp(mi,1),'qcseg') && (recalc.qa(1)<2);
        % check data
        if P.p0e{mi}(pi)
          if ~exist(P.xml{mi}{pi},'file')
            P.p0e{mi}(pi) = 0; 
          else
            xmlt = cat_io_xml(P.xml{mi}{pi});
            if ~isfield(xmlt,'filedata') || isempty(xmlt.filedata)
              P.p0e{mi}(pi) = 0; 
            end
          end
        end
      end
      
      %% (re)calcqa and load xml files
      P.qamat{mi} = fullfile(opt.resdir,['bwp_' qafile '_NIR' P.bwp{mi,1} '.mat']);
      if ~exist(P.qamat{mi},'file') || recalc.qa(1)>0 || any(P.p0e{mi}==0)
        qaopt = struct('prefix',[qafile '_'],'write_csv',0,'mprefix','m','orgval',0,'rerun',recalc.qa(1)); 
        if strcmp(P.bwp(mi,1),'synthseg')
          qav = [qaversions{qai} '_synthseg_']; 
        elseif strcmp(P.bwp(mi,1),'qcseg')
          qav = [qaversions{qai} '_qcseg_']; 
        elseif strcmp(P.bwp(mi,1),'SPM12')
          qav = [qaversions{qai} '_spm_']; 
        else
          qav = [qaversions{qai} '_']; 
        end
        if recalc.qa(1)>1
          cat_vol_qa('p0',P.p0a{mi},struct('prefix',qav,'version',qafile,'rerun',2));
        else
          cat_vol_qa('p0',P.p0a{mi}(P.p0e{mi}==0),struct('prefix',qav,'version',qafile,'rerun',0));
        end
        xml = cat_io_xml(P.xml{mi});
        save(P.qamat{mi},'xml');
      else
        load(P.qamat{mi},'xml');
      end
      X(mi).xml = xml;
      
  
  
  
      % -- estimate kappa --------------------------------------------------
      for ki=1:numel(P.p0{mi})
        if strcmp(P.bwp(mi,1),'qcseg')
          P.p0r(ki) = ~isempty(strfind(P.p0{mi}{ki},'rBWP')); 
        else
          P.p0r(ki) = ~isempty(strfind(P.p0{mi}{ki},'p0r')); 
        end
      end
      if fasttest
        P.kappamat{mi} = fullfile(opt.resdir,sprintf('bwp_kappa_NIR%d_%s_fast.mat',numel(P.p0{mi}),P.bwp{mi,1}));
      else
        P.kappamat{mi} = fullfile(opt.resdir,sprintf('bwp_kappa_NIR%d_%s.mat',numel(P.p0{mi}),P.bwp{mi,1}));
      end
      if (qai == 1 && ( ~exist(P.kappamat{mi},'file')) || recalc.kappa(1))
        %%
        clear val; 
        [~,valo] = eva_vol_calcKappa(P.p0{mi}(~P.p0r),{P.p0gt},struct('recalc',recalc.kappa(1),'realign',0,'realignres',1));
        [~,vali] = eva_vol_calcKappa(P.p0{mi}(P.p0r) ,{P.p0gt},struct('recalc',recalc.kappa(1),'realign',2,'realignres',2));
        val(1,:)                = valo(1,:);
        val(find(P.p0r==0)+1,:) = valo(2:end-2,:);
        val(find(P.p0r==1)+1,:) = vali(2:end-2,:); 
        val(end+1:end+2,:)      = vali(end-1:end,:); 
        save(P.kappamat{mi},'val')
      else
        load(P.kappamat{mi})
      end
  
      try
        if mi==1
          Q.kappa = mean(cell2mat(val(2:end-2,2:4)),2);
        else
          Q.kappa(:,mi) = mean(cell2mat(val(2:end-2,2:4)),2);
        end  
      catch
        %eva_vol_calcKappa(P.p0{mi}(P.p0r),{P.p0gt},struct('recalc',1,'realign',2,'realignres',1));
        %[~,val] = eva_vol_calcKappa(P.p0{mi},P.p0gt,struct('recalc',0,'realign',0,'realignres',1));
        %Q.kappa(:,mi) = mean(cell2mat(val(2:end-2,2:4)),2);
        %save(P.kappamat{mi},'val')
      end
      clear val
    end
    
    
    
    
    %% -- prepare QA data ------------------------------------------------
    cat_io_cmd('  Prepare quality measurements:','g5','',1); 
    nx = 4; 
  
    % create a new, simpler structure based on the xml-files
    mi = find(method>0,1,'first');
    Q.method  = repmat(1:5,size(Q.kappa,1),1);
    Q.noise   = nan(numel(P.xml{mi}),1); 
    Q.bias    = nan(numel(P.xml{mi}),1); 
    Q.field   = nan(numel(P.xml{mi}),1); 
    Q.vx_voli = nan(numel(P.xml{mi}),3); 
    Q.vx_vol  = nan(numel(P.xml{mi}),3); vx_vol_digit = 3;
    Q.vx_int  = nan(numel(P.xml{mi}),1); 
    Q.ECR     = nan(numel(P.xml{mi}),1); 
    Q.ECRgt   = nan(numel(P.xml{mi}),1); 
    Q.ECRmm   = nan(numel(P.xml{mi}),1); 
    Q.ECRmmgt = nan(numel(P.xml{mi}),1); 
    if fasttest
      Q.train = ones(numel(P.xml{mi}),1); % train/test subset (split half) 
    else
      %Q.train = rand(numel(Q.bias),2)' > 5; % train/test subset (split half) 
      Q.train = mod(1:numel(Q.bias),8)' < 4; % train/test subset (split half) 
    end
    for pi=1:numel(P.xml{mi})
      [~,ff] = spm_fileparts(P.xml{mi}{pi});
      pni  = strfind(ff,'pn');
      rfi  = strfind(ff,'rf');
      vxi  = strfind(ff,'vx');
  
      Q.noise(pi)    = str2double(ff(pni+2));
      Q.bias(pi)     = sign(0.5-(ff(rfi+5)=='n')).*str2double( ff(rfi+2:rfi+4));
      Q.field(pi)    = char(ff(rfi+6));
      Q.vx_int(pi,1) = contains(ff,'irBWP');
      Q.SIQRgt   = repmat(((Q.bias-20)/160*4)+1,1,numel(method));
      Q.vx_vol(pi,:) = [str2double(ff(vxi+2                : vxi+1+vx_vol_digit   )), ...
                        str2double(ff(vxi+3+vx_vol_digit   : vxi+2+vx_vol_digit*2 )), ...
                        str2double(ff(vxi+4+vx_vol_digit*2 : vxi+3+vx_vol_digit*3 ))]/100;
    end
    Q.NCRgt   = repmat(((Q.noise-1)/8*4)+1,1,numel(method));
    Q.ICRgt   = repmat(((Q.bias-20)/160*4)+1,1,numel(method));
    Q.RESgt   = repmat((mean(Q.vx_vol.^2,2).^0.5)*2,1,numel(method));
    Q.ECRgt   = Q.RESgt - 1; 
    Q.FECgt   = Q.NCRgt; 
    Q.IQRgt   = mean([Q.NCRgt, Q.ICRgt, Q.RESgt].^nx,2).^(1/nx);
    Q.SIQRgt  = mean([Q.NCRgt, Q.ICRgt, Q.RESgt, Q.ECRgt, Q.FECgt].^nx,2).^(1/nx);
    
    
    %  -- map QM data ----------------------------------------------------
    fieldxml = 'xml';
    QM = {'NCR','ICR','res_RMS','res_ECR','res_ECRmm','FEC','contrastr','IQR','SIQR'}; % ################## reprocess *IQR*
    for qmi = 1:numel(QM)
      % create a version with original values
      if any(strcmp(QM{qmi},{'NCR','ICR','res_ECR','res_ECRmm','FEC'}))
        field   = [QM{qmi} 'o']; 
        fieldqm = 'qualitymeasures';
      else
        field   = QM{qmi};
        fieldqm = 'qualityratings';
      end
      Q.(field) = zeros(numel(Q.noise),numel(method));
      for pi=1:numel(P.(fieldxml){mi})
        for mi = method
          try 
            Q.(field)(pi,mi) = X(mi).(fieldxml)(pi).(fieldqm).(QM{qmi}); 
          catch
            Q.(field)(pi,mi) = nan; 
          end
        end
      end
    end
    for pi=1:numel(P.(fieldxml){mi})
      try
        Q.vx_voli(pi,:)  =  X(mi).(fieldxml)(pi).qualitymeasures.res_vx_vol; 
      end
    end
    Q.CMV = zeros(numel(Q.noise),numel(method));
    Q.GMV = zeros(numel(Q.noise),numel(method));
    Q.WMV = zeros(numel(Q.noise),numel(method));
    for mi = method
      for pi=1:numel(P.xml{mi})
        try
          Q.CMV(pi,mi)   =  X(mi).xml(pi).subjectmeasures.vol_rel_CGW(1);
          Q.GMV(pi,mi)   =  X(mi).xml(pi).subjectmeasures.vol_rel_CGW(2);
          Q.WMV(pi,mi)   =  X(mi).xml(pi).subjectmeasures.vol_rel_CGW(3);
        catch
          fprintf('Failed: X(%d).xml(%d)\n',mi,pi)
          Q.CMV(pi,mi)   =  nan;
          Q.GMV(pi,mi)   =  nan;
          Q.WMV(pi,mi)   =  nan;
        end
      end
    end
    
    
    %  -- fit marks ------------------------------------------------------ 
    % * 0% noise and 0% bias are excluded because they are instable in SPM pp 
    %   and also smaller datasets (no fields) 
    % * more problematic is the detection of using low resolution data also
    %   for the general evaluation because the data is (i) no real standard,
    %   (ii) not representing the majority of data, and (iii) less accurate/robust
  for tstseti = 4  
    default = cat_stat_marks('default');
  
    opt.bwptestset = tstseti;  % ######## focus on 1 mm ! 
    switch opt.bwptestset
      case 1, M = (fasttest | Q.train) & Q.noise>0 & Q.bias>0 & Q.vx_int==0 & all(Q.vx_vol==1,2); % 1 mm focus for better results (in general the images have this resolution)
      case 2, M = (fasttest | Q.train) & Q.noise>0 & Q.bias>0 & Q.vx_int==0;                      % only non interpolated cases
      case 3, M = (fasttest | Q.train) & Q.noise>0 & Q.bias>0 & all(Q.vx_voli==1,2);              % only isotropic cases
      case 4, M = (fasttest | Q.train) & Q.noise>0 & Q.bias>0;                                    % all cases 
    end
    MECR  = (fasttest |  Q.train) & Q.noise>0 & Q.bias>0 & (all(Q.vx_vol==1,2) | all(Q.vx_vol==2,2) ); % here we need also the other resolutions 
    MECRT = (fasttest | ~Q.train) & Q.noise>0 & Q.bias>0; % here we need also the other resolutions 
    
  
    %  -- fit for rating -------------------------------------------------- 
    rmse = @(a,b) cat_stat_nanmean( (a-b).^2 ).^(1/.5);
    % NCR
    cat_io_cmd(sprintf('  Mark range of %s (N_train=%d/%d):',qafile,sum(M),sum(Q.train)),'n','',1); 
    [Q.fit.noiseNCR, Q.fit.noiseNCRstat] = robustfit(Q.NCRgt(M,1),Q.NCRo(M,1));
    default.QS{contains(default.QS(:,2),'NCR'),4} = ...
      [Q.fit.noiseNCR(1) + Q.fit.noiseNCR(2),Q.fit.noiseNCR(1) + Q.fit.noiseNCR(2) * 6];  %#ok<*FNDSB> 
    
    % ICR
    [Q.fit.biasICR, Q.fit.biasICRstat] = robustfit(Q.ICRgt(M,1),Q.ICRo(M,1));
    default.QS{contains(default.QS(:,2),'ICR'),4} = ...
      [Q.fit.biasICR(1) + Q.fit.biasICR(2),Q.fit.biasICR(1) + Q.fit.biasICR(2) * 6]; 
    
    % ECR
    ECRpos   = find(cellfun('isempty',strfind(default.QS(:,2),'ECR'  ))==0,1,'first'); %#ok<*STRCL1> 
    [Q.fit.resECR, Q.fit.resECRstat] = robustfit(Q.ECRgt(MECR,1),Q.res_ECRo(MECR,1));
    default.QS{ECRpos,4} = [Q.fit.resECR(1) + Q.fit.resECR(2),Q.fit.resECR(1) + Q.fit.resECR(2) * 6]; % we start with mark 2
    save(fullfile(opt.resdir,['qadefault_' qafile '.mat']),'default');
    fprintf('\n');
    
    % FEC
    FECpos = find(cellfun('isempty',strfind(default.QS(:,2),'FEC'))==0);
    try
      warning off; 
      [Q.fit.FEC, Q.fit.FECstat]  = robustfit(Q.FECgt(M,1),Q.FECo(M,1));
      warning on; 
      if ~isempty(FECpos)
        default.QS{FECpos,4} = round([Q.fit.FEC(1) + Q.fit.FEC(2), Q.fit.FEC(1) + Q.fit.FEC(2) * 6], -1);
      end
    catch
      Q.fit.FEC = [nan nan]; Q.fit.FECstat = struct('coeffcorr',nan(2,2),'p',nan(2,2)); 
      default.QS{find(cellfun('isempty',strfind(default.QS(:,2),'FEC'))==0),4} = [100 850];
    end

   
    %  -- remarks -------------------------------------------------------- 
    cat_io_cmd('  Remark:','g5','',1); 
    for mi = method
      for pi=1:numel(P.xml{mi})
        try
          X(mi).xml2(pi)   = cat_stat_marks('eval',0,X(mi).xml(pi),default);
          Q.NCR(pi,mi)     = X(mi).xml2(pi).qualityratings.NCR;
          Q.ICR(pi,mi)     = X(mi).xml2(pi).qualityratings.ICR;
          if isfield(  X(mi).xml2(pi).qualityratings,'res_ECR')
            Q.ECR(pi,mi)   = X(mi).xml2(pi).qualityratings.res_ECR;
          else
            Q.ECR(pi,mi)   = nan; 
          end
          Q.IQR(pi,mi)     = X(mi).xml2(pi).qualityratings.IQR;
          try 
            Q.SIQR(pi,mi)  = X(mi).xml2(pi).qualityratings.SIQR;
          catch
            Q.SIQR(pi,mi)  = X(mi).xml2(pi).qualityratings.IQR;
          end
          if isfield(  X(mi).xml2(pi).qualityratings,'FEC')
            Q.FEC(pi,mi)   = X(mi).xml2(pi).qualityratings.FEC;
          else
            Q.FEC(pi,mi)   = nan; 
          end
        catch
          Q.NCR(pi,mi)     = nan;
          Q.ICR(pi,mi)     = nan;
          Q.ECR(pi,mi)     = nan;
          Q.FEC(pi,mi)     = nan;
          Q.IQR(pi,mi)     = nan;
          Q.SIQR(pi,mi)    = nan;
        end
      end
    end    
    if all( isnan( Q.SIQR(:) ) )
      Q.SIQR = Q.IQR;
    end
  
    fprintf('\n    NCR:   [%8.4f, %8.4f] (r=%0.4f, p=%0.1e, RMSE=%0.4f)',...
      default.QS{contains(default.QS(:,2),'NCR'),4}, Q.fit.noiseNCRstat.coeffcorr(2), ...
      Q.fit.noiseNCRstat.p(2), rmse( Q.NCRgt(M,1), Q.NCR(M,1)) );
    fprintf('\n    ICR:   [%8.4f, %8.4f] (r=%0.4f, p=%0.1e, RMSE=%0.4f)',... - keep in mind that the BWP bias is scaled to 3 (i.e., mult. by 2 times!)
      default.QS{contains(default.QS(:,2),'ICR'),4}, Q.fit.biasICRstat.coeffcorr(2), ...
      Q.fit.biasICRstat.p(2), rmse( Q.ICRgt(M,1), Q.ICR(M,1)));
    fprintf('\n    ECR:   [%8.4f, %8.4f] (r=%0.4f, p=%0.1e, RMSE=%0.4f)',...
      default.QS{ECRpos,4}, Q.fit.resECRstat.coeffcorr(2), ...
      Q.fit.resECRstat.p(2), rmse( Q.ECRgt(M,1), Q.ECR(M,1)));
    fprintf('\n    FEC:   [%8.4f, %8.4f] (r=%-7.4f, p=%0.1e, RMSE=%0.4f)\n',...
      default.QS{FECpos,4}, Q.fit.FECstat.coeffcorr(2), ...
      Q.fit.FECstat.p(2), rmse( Q.FECgt(M,1), Q.FEC(M,1)) );
  
   
  
    
    % -- results -------------------------------------------------------- 
    % table for IQR values by method ...
    def.bstm    = 1;          % best mark
    def.wstm    = 6;          % worst mark
    def.wstmn   = 8.624;      % worst mark to get a 4 for values with std 
    def.bstl    = 0.5+eps*2;  % highest rating ... 0.5 because of rounding values
    def.wstl    = 10.5-eps*2; % lowest rating  ... to have only values with 1 digit .. but % scaling...
    
    setnan      = [1 nan];
    evallinear  = @(x,bst,wst)  setnan(isnan(x)+1) .* ... 
      (min(def.wstl,max(def.bstl,(sign(wst-bst)*x - sign(wst-bst)*bst) ./ abs(diff([wst ,bst])) .* abs(diff([def.bstm,def.wstm])) + def.bstm)));
    rms         = @(a,fact)   max(0,cat_stat_nanmean(a.^fact).^(1/fact));
    mark2rps    = @(mark) min(100,max(0,105 - mark*10));
     
     
    
    %% -- create figure --------------------------------------------------
  
    cat_io_cmd('  Create figure:','g5','',1); 
    testmethods = [1 0]; testmethods(numel(method)+1:end) = [];
  
    % BWP main figure
    sizex=9;
    FS=[7 6 2]*(opt.dpi/120)*2; b2 = 0.04/(9/9);  subx = 5; suby = 4; pos = cell(suby,subx); dimx=sizex*opt.dpi*1.2; 
    for x=1:subx
      for y=1:suby
        pos{y,x}=[(x-1)/subx + b2*1.2, (suby - y)/suby + b2*1.5 , 1/subx-b2*1.3 , 1/suby - b2*2];
      end
    end
  
    if exist('fh1','var') && ishghandle(fh1) 
      clf(fh1); %figure(fh1); 
    else
      fh1 = figure('Name','figure 1 - bwp main','Visible','off','Interruptible', 'off','Position', ...
        [0 0 dimx round(dimx*(suby/subx))],'color',[1 1 1],'PaperPositionMode','auto');
    end
    
  % ----------------------------------------------------------------------
  % A) BWP test variables
  %    - measures should be more accurate in high quality data
  %    - the resolution datasets for 1 and 2 mm are 3x smaller than the 
  %      1x1x2 and 1x2x2 datasets because of the permutations 
  % ----------------------------------------------------------------------
    rps  = 1; % ########## use rps is here a bit inoptimal 
    rpsnam = {'grade','%'};
    if rps, nlim = [40 100]; else, nlim = [0.5 6.5]; end 
    if rps, mlim = mark2rps(6:-1:1); else, mlim = 1:6; end 
    
    QM = {'NCR','ICR','res_RMS','res_ECR','FEC','SIQR'};   
    cl = [.8 0 .2;   0.8 0.6 0;  0.2 0.5 0;  .0 .5 .9; 0. 0 0.9;  0 0 0 ]; 
  
    % noise
    % --------------------------------------------------------------------
    subplot('Position',pos{1,1}); clear MMC
    M3=repmat(testmethods,size(M,1),1);
    M2=repmat(M,numel(method),1) & M3(:);
    [nln,~,nlid2]=unique(repmat(Q.noise(M),numel(method),1)); MM=repmat(Q.NCR(M2),numel(method),1);   %.*sum(M2>0,1)/numel(M2)
    for nlni=1:numel(nln); if rps, MMC{nlni} = mark2rps(MM(nlid2==nlni)); else, MMC{nlni} = MM(nlid2==nlni); end; end
    cat_plot_boxplot(MMC,struct('names',...
        [num2str(unique(Q.noise(M))) repmat('%',numel(num2str(unique(Q.noise(M)))),1)], ...
        'sort',0,'ylim',nlim,'groupnum',0,'ygrid',0,'groupcolor',...
        min(1,flip([1./(2:-1/(numel(MMC)/2-1):1) 1:1/((numel(MMC)-1)/2):2])' * cl(1,:) * .6), ...
        'style',4,'datasymbol','o','usescatter',1)); set(gca,'FontSize',FS(2)); 
    set(gca,'ylim',nlim,'ytick',mlim,'FontSize',FS(2),'YGrid','on'); ylabel('NCR'); 
    title(sprintf('noise (n=%d)', numel(M)),'FontSize',FS(1),'FontWeight','bold'); 
    xlabel('BWP noise levels'); ylabel(sprintf('NCR (%s)',rpsnam{rps+1})); 
    if ~rps, set(gca,'YDir','reverse'); end
  
    
    % bias
    % --------------------------------------------------------------------
    subplot('Position',pos{1,2}); clear MMC
    M3=repmat(testmethods,size(M,1),1);
    M2=repmat(M,numel(method),1) & M3(:);  
    [nln,~,nlid2]=unique(repmat(Q.bias(M),numel(method),1)); MM=repmat(Q.ICR(M2),numel(method),1);   
    for nlni=1:numel(nln); if rps, MMC{nlni} = mark2rps(MM(nlid2==nlni)); else, MMC{nlni} = MM(nlid2==nlni); end; end
    cat_plot_boxplot(MMC,struct('names',[num2str(nln) repmat('%',size(num2str(nln),1),1)], ...
      'sort',0,'ylim',nlim,'groupnum',0,'ygrid',0,'groupcolor',...
        min(1,flip([1./(2:-1/(numel(MMC)/2-1):1) 1:1/((numel(MMC)-1)/2):2])' * cl(2,:) * .6), ...
        'style',4,'datasymbol','o','usescatter',1)); set(gca,'FontSize',FS(2)); hold on
    set(gca,'ylim',nlim,'ytick',mlim,'FontSize',FS(2),'YGrid','on','XTickLabelRotation',0);
    title('inhomogeneity','FontSize',FS(1),'FontWeight','bold'); 
    xlabel('BWP inhomogeneity levels'); ylabel(sprintf('ICR (%s)',rpsnam{rps+1})); 
    if ~rps, set(gca,'YDir','reverse'); end
  
  
    % resolution
    % --------------------------------------------------------------------
    subplot('Position',pos{1,3},'box','on'); hold on; clear res
    res.sf   = 2; 
    res.rms  = @(a) (mean(a.^res.sf)).^(1/res.sf);
    res.defres = 1;
    res.res  = 0.25:0.1:3.25; 
    res.resp = res.defres:0.5:3.25;
    res.resa = res.defres:0.5:3.25;
    res.iso  = [res.res' res.res' res.res'];
    res.ani  = [res.defres*ones(numel(res.res),2) res.res'];
    res.ani2 = [res.defres*ones(numel(res.res),1) res.res' res.res'];
    res.isop = [res.resp' res.resp' res.resp'];
    res.anip = [res.defres*ones(numel(res.resa),2) res.resa'];
    res.anip2= [res.defres*ones(numel(res.resp),1) res.resp' res.resp'];
    plot(res.res,res.rms(res.iso'),'r-'); plot(res.res,res.rms(res.ani'),'b-'); plot(res.res,res.rms(res.ani2'),'g-');
    plot(res.resp,res.rms(res.isop'),'r+'); plot(res.resa,res.rms(res.anip'),'b+'); plot(res.resp,res.rms(res.anip2'),'g+');
    legend({'isotr. [R,R,R]',...
            sprintf('sliceres. [%0.0f %0.0f R]',res.defres,res.defres),...
            sprintf('sliceth. [R R %0.0f]',res.defres)},...
            'Location','Southwest','Fontsize',FS(2)*0.8);
    xlim([0.25 res.res(end)+eps]); ylim([0.25 res.res(end)]); 
    set(gca,'xtick',0.5:0.5:res.resp(end)+eps,'ytick',0:0.5:res.resp(end),'FontSize',FS(2));
    set(gca,'ylim',[0.25 3.25],'ytick',0.5:0.5:3,'yticklabel',num2str((1:6)'),'FontSize',FS(2));
    title('voxel resolution','FontSize',FS(1),'FontWeight','bold'); 
    xlabel('BWP resolution (R in mm)'); ylabel(sprintf('RES (%s)',rpsnam{rps+1})); 
    grid on; 
    if ~rps, set(gca,'YDir','reverse'); end
  
  
    % ECR
    % --------------------------------------------------------------------
    if isfield(  X(mi).xml2(pi).qualityratings,'res_ECR')
      subplot('Position',pos{1,4}); clear MMC
      M3=repmat(testmethods,size(MECRT,1),1);
      M2=repmat(MECRT,numel(method),1) & M3(:);  
      [nln,~,nlid2]=unique(repmat(Q.ECRgt(MECRT),numel(method),1)); MM=repmat(Q.ECR(M2),numel(method),1);   
      for nlni=1:numel(nln); if rps, MMC{nlni} = mark2rps(MM(nlid2==nlni)); else, MMC{nlni} = MM(nlid2==nlni); end; end
      switch size(MMC,2)
        case 2, MMCnames = {'1x1x1','2x2x2'}; 
        case 3, MMCnames = {'1x1x1','1x1x2','2x2x2'}; 
        case 4, MMCnames = {'1x1x1','1x1x2','1x2x2','2x2x2'}; 
      end
      cat_plot_boxplot(MMC,struct('names',{MMCnames}, ... [num2str(nln,'%0.2f')], ...
        'sort',0,'ylim',nlim,'groupnum',0,'ygrid',0, 'groupcolor',...
        min(1,flip([1./(2:-1/(numel(MMC)/2-1):1) 1:1/((numel(MMC)-1)/2):2])' * cl(4,:) * .6), ...
        'style',4,'datasymbol','o','usescatter',1)); set(gca,'FontSize',FS(2)); hold on
      set(gca,'ylim',nlim,'ytick',mlim,'FontSize',FS(2),'YGrid','on');
      title('edge resolution','FontSize',FS(1),'FontWeight','bold'); 
      xlabel('BWP resolution levels'); ylabel(sprintf('ECR (%s)',rpsnam{rps+1})); 
    end
    if ~rps, set(gca,'YDir','reverse'); end
        
  
    % FEC
    % --------------------------------------------------------------------
    subplot('Position',pos{1,5}); clear MMC
    M3=repmat(testmethods,size(M,1),1);
    M2=repmat(M,numel(method),1) & M3(:);
    [nln,~,nlid2]=unique(repmat(Q.noise(M),numel(method),1)); MM=repmat(Q.FEC(M2),numel(method),1);  
    for nlni=1:numel(nln); if rps, MMC{nlni} = mark2rps(MM(nlid2==nlni)); else, MMC{nlni} = MM(nlid2==nlni); end; end
    cat_plot_boxplot(MMC,struct('names',...
      [num2str(unique(Q.noise(M))) repmat('%',numel(num2str(unique(Q.noise(M)))),1)], ...
      'sort',0,'ylim',nlim,'groupnum',0,'ygrid',0,'groupcolor',...
       min(1,flip([1./(2:-1/(numel(MMC)/2-1):1) 1:1/((numel(MMC)-1)/2):2])' * cl(5,:) * .6), ...
      'style',4,'datasymbol','o','usescatter',1)); set(gca,'FontSize',FS(2))
    set(gca,'ylim',nlim,'ytick',mlim,'FontSize',FS(2),'YGrid','on'); ylabel('NCR');  
    title(sprintf('Fast Euler Characteristic', opt.bwptestset,numel(M)),'FontSize',FS(1),'FontWeight','bold'); 
    xlabel('BWP noise level'); ylabel(sprintf('FEC (%s)',rpsnam{rps+1})); 
    if ~rps, set(gca,'YDir','reverse'); end
  
    
    % RMSEs 
    % --------------------------------------------------------------------
     if rps
       NCRgtdiff   = mark2rps(Q.NCR(M,1))      - mark2rps(Q.NCRgt(M));
       ICRgtdiff   = mark2rps(Q.ICR(M,1))      - mark2rps(Q.ICRgt(M));
       RESgtdiff   = mark2rps(Q.res_RMS(M,1))  - mark2rps(Q.ECRgt(M));
       ECRgtdiff   = mark2rps(Q.ECR(M,1))      - mark2rps(Q.ECRgt(M));
       FECgtdiff   = mark2rps(Q.FEC(M,1))      - mark2rps(Q.FECgt(M));
       CONgtdiff   = mark2rps(Q.contrastr(M,1)) - mark2rps( median(Q.contrastr(M)) );
       IQRgtdiff   = mark2rps(Q.IQR(M,1))      - mark2rps(Q.IQRgt(M));
       SIQRgtdiff  = mark2rps(Q.SIQR(M,1))     - mark2rps(Q.SIQRgt(M));
       
       NCRgtdiffn  = (mark2rps(Q.NCR(M,1))     - mark2rps(Q.NCRgt(M)) )  ./  mark2rps(Q.NCRgt(M))  * 100;
       ICRgtdiffn  = (mark2rps(Q.ICR(M,1))     - mark2rps(Q.ICRgt(M)) )  ./  mark2rps(Q.ICRgt(M))  * 100;
       RESgtdiffn  = (mark2rps(Q.res_RMS(M,1)) - mark2rps(Q.ECRgt(M)) )  ./  mark2rps(Q.ECRgt(M))  * 100;
       ECRgtdiffn  = (mark2rps(Q.ECR(M,1))     - mark2rps(Q.ECRgt(M)) )  ./  mark2rps(Q.ECRgt(M))  * 100;
       FECgtdiffn  = (mark2rps(Q.FEC(M,1))     - mark2rps(Q.FECgt(M)) )  ./  mark2rps(Q.FECgt(M))  * 100;
       CONgtdiffn  = (mark2rps(Q.contrastr(M,1)) - mark2rps( median(Q.contrastr(M)) ) )  ./  mark2rps( median(Q.contrastr(M)) )  * 100;
       IQRgtdiffn  = (mark2rps(Q.IQR(M,1))     - mark2rps(Q.IQRgt(M)) )  ./  mark2rps(Q.IQRgt(M))  * 100;
       SIQRgtdiffn = (mark2rps(Q.SIQR(M,1))    - mark2rps(Q.SIQRgt(M)) ) ./  mark2rps(Q.SIQRgt(M)) * 100;
  
       rmseNCR     = rms( mark2rps(Q.NCR(M,1))      - mark2rps(Q.NCRgt(M)) , 2 );
       rmseICR     = rms( mark2rps(Q.ICR(M,1))      - mark2rps(Q.ICRgt(M)) , 2 );
       rmseRES     = rms( mark2rps(Q.res_RMS(M,1))  - mark2rps(Q.ECRgt(M)) , 2 );
       rmseECR     = rms( mark2rps(Q.ECR(M,1))      - mark2rps(Q.ECRgt(M)) , 2 );
       rmseFEC     = rms( mark2rps(Q.FEC(M,1))      - mark2rps(Q.FECgt(M)) , 2 );
       rmseCON     = rms( mark2rps(Q.contrastr(M,1)) - mark2rps( median(Q.contrastr(M)))  , 2 );
       rmseIQR     = rms( mark2rps(Q.IQR(M,1))      - mark2rps(Q.IQRgt(M)) , 2 );
       rmseSIQR    = rms( mark2rps(Q.SIQR(M,1))     - mark2rps(Q.SIQRgt(M)) , 2 );
       
       Bscaling    = { [-30 30] , [-30 30] }; 
     else
       NCRgtdiff   = (Q.NCR(M,1))      - (Q.NCRgt(M));
       ICRgtdiff   = (Q.ICR(M,1))      - (Q.ICRgt(M));
       RESgtdiff   = (Q.res_RMS(M,1))  - (Q.ECRgt(M));
       ECRgtdiff   = (Q.ECR(M,1))      - (Q.ECRgt(M));
       FECgtdiff   = (Q.FEC(M,1))      - (Q.FECgt(M));
       CONgtdiff   = (Q.contrastr(M,1)) - median(Q.contrastr(M));
       IQRgtdiff   = (Q.IQR(M,1))      - (Q.IQRgt(M));
       SIQRgtdiff  = (Q.SIQR(M,1))     - (Q.SIQRgt(M));
       
       NCRgtdiffn  = ((Q.NCR(M,1))     - (Q.NCRgt(M)) )  ./  (Q.NCRgt(M))  * 100;
       ICRgtdiffn  = ((Q.ICR(M,1))     - (Q.ICRgt(M)) )  ./  (Q.ICRgt(M))  * 100;
       RESgtdiffn  = ((Q.res_RMS(M,1)) - (Q.ECRgt(M)) )  ./  (Q.ECRgt(M))  * 100;
       ECRgtdiffn  = ((Q.ECR(M,1))     - (Q.ECRgt(M)) )  ./  (Q.ECRgt(M))  * 100;
       FECgtdiffn  = ((Q.FEC(M,1))     - (Q.FECgt(M)) )  ./  (Q.FECgt(M))  * 100;
       CONgtdiffn  = ((Q.contrastr(M,1)) - median(Q.contrastr(M)) )  ./  median(Q.contrastr(M))  * 100;
       IQRgtdiffn  = ((Q.IQR(M,1))     - (Q.IQRgt(M)) )  ./  (Q.IQRgt(M))  * 100;
       SIQRgtdiffn = ((Q.SIQR(M,1))    - (Q.SIQRgt(M)) ) ./  (Q.SIQRgt(M)) * 100;
  
       rmseNCR     = rms( (Q.NCR(M,1))      - (Q.NCRgt(M)) , 2 );
       rmseICR     = rms( (Q.ICR(M,1))      - (Q.ICRgt(M)) , 2 );
       rmseRES     = rms( (Q.res_RMS(M,1))  - (Q.ECRgt(M)) , 2 );
       rmseECR     = rms( (Q.ECR(M,1))      - (Q.ECRgt(M)) , 2 );
       rmseFEC     = rms( (Q.FEC(M,1))      - (Q.FECgt(M)) , 2 );
       rmseCON     = rms( (Q.contrastr(M,1)) - ( median(Q.contrastr(M)))  , 2 );
       rmseIQR     = rms( (Q.IQR(M,1))      - (Q.IQRgt(M)) , 2 );
       rmseSIQR    = rms( (Q.SIQR(M,1))     - (Q.SIQRgt(M)) , 2 );
            
       Bscaling    = { [-3 3]  , [-100 100] }; 
     end
    
     % absolute and RMS error
     subplot('Position',pos{3,4}); 
     cat_plot_boxplot({NCRgtdiff, ICRgtdiff, RESgtdiff, ECRgtdiff, FECgtdiff, SIQRgtdiff},...
       struct('names',{{'NCR','ICR','RES','ECR','FEC','SIQR'}},'sort',0,'ylim',Bscaling{1},...
       'groupcolor',cl,'style',4,'datasymbol','o','usescatter',1,'groupnum',0,'ygrid',0)); 
     h = gca; h.XTickLabelRotation = 0;  h.YGrid = 'on'; 
     ylabel(sprintf('rating (%s)',rpsnam{rps+1})); xlabel('(measure - goldstd)'); set(gca,'FontSize',FS(2));
     title('absolute error','FontSize',FS(1),'FontWeight','bold'); 
     %if ~rps, set(gca,'YDir','reverse'); end
  
     subplot('Position',pos{3,5}); 
     rmseval   = [ rmseNCR, rmseICR, rmseRES, rmseECR, rmseFEC, rmseSIQR ]; 
     bh = bar(rmseval(1:end)); bh.CData = cl; bh.FaceColor = 'flat';
     ylim(max(0,Bscaling{1})); xlim([.4 6.6]); xticklabels({'NCR','ICR','RES','ECR','FEC','SIQR'}); 
     h = gca; h.YGrid = 'on';  h.XTickLabelRotation = 0; 
     ylabel(sprintf('rating (%s)',rpsnam{rps+1})); xlabel('root(mean((measure - goldstd)^2))'); set(gca,'FontSize',FS(2));
     title('RMSE');
     for fi = 1:numel(rmseval)
       dt(fi) = text(fi-.45, rmseval(fi) + .1 + 1*(rps), sprintf('%0.3f',rmseval(fi)),'FontSize',8,'Color',cl(fi,:));
     end  
     fprintf('\nRMSE values from figure: \n')
     fprintf('Measure: %s\n', sprintf(repmat('%10s',1,6),'NCR','ICR','RES','ECR','FEC','SIQR') ); 
     fprintf('RMSE:    %s\n', sprintf('%10.4f',rmseval) );
     fprintf('\n')
     tabRMSE = [{'Measure:','NCR','ICR','RES','ECR','FEC','SIQR'}; 
          {'RMSE:'}, num2cell(rmseval) ]; 
     fname = fullfile(opt.resdir,sprintf('tst_bwpmaintest_%s_RMSEtable_%s.csv',qafile,segment{mi})); 
     cat_io_csv(fname, tabRMSE); 
     cat_io_cprintf('blue',sprintf('    Write %s\n',fname)); 
     %if ~rps, set(gca,'YDir','reverse'); end
  

  %
  % ----------------------------------------------------------------------
  % B) relation between different BWP aspects and the QM
  % ----------------------------------------------------------------------
     
     % The relative error is smaller compared to the absolute one because the
     % higher variance in low res data is weighted lower.
     % However, I am not sure if this is practically useful and clear. 
     % For development it was more relevant to know ...
  
     % NCR
     subplot('Position',pos{2,1}); clear MMC; 
     M2 = (fasttest | ~Q.train) & Q.noise==5 & Q.bias==20  & Q.vx_int==0;                       MMC{1}=repmat(Q.NCR(M2),numel(method),1);   
     M2 = (fasttest | ~Q.train) & Q.noise==5 & Q.bias==60  & Q.vx_int==0;                       MMC{2}=repmat(Q.NCR(M2),numel(method),1);   
     M2 = (fasttest | ~Q.train) & Q.noise==5 & Q.bias==100 & Q.vx_int==0;                       MMC{3}=repmat(Q.NCR(M2),numel(method),1);   
     M2 = (fasttest | ~Q.train) & Q.noise==5 & Q.bias>0    & Q.vx_int==0 &  all(Q.vx_vol==1,2); MMC{4}=repmat(Q.NCR(M2),numel(method),1);   
     M2 = (fasttest | ~Q.train) & Q.noise==5 & Q.bias>0    & Q.vx_int==0 &  all(Q.vx_vol==2,2); MMC{5}=repmat(Q.NCR(M2),numel(method),1);   
     M2 = (fasttest | ~Q.train) & Q.noise==5 & Q.bias>0    & Q.vx_int==1;                       MMC{6}=repmat(Q.NCR(M2),numel(method),1);  
     if rps, for ci=1:numel(MMC), MMC{ci} = mark2rps(MMC{ci}); end; end
     cat_plot_boxplot(MMC,struct('names',{{'B2','B6','B10','R','rR','irR'}}, 'subsets',[0 0 0 1 1 1],...
       'style',4,'datasymbol','o','usescatter',1,'sort',0,'ylim',nlim,'groupnum',0,'ygrid',0)); set(gca,'FontSize',FS(2)); hold on
     set(gca,'ylim',nlim,'ytick',mlim,'FontSize',FS(2),'ygrid',1);
     title('NCR at 5% noise (grad 3)','FontSize',FS(1),'FontWeight','bold'); 
     xlabel(sprintf('BWP bias (%d) / RES (%d-%d)',numel(MMC{1}),numel(MMC{4}),numel(MMC{5}))); ylabel(sprintf('NCR (%s)',rpsnam{rps+1}));
     if ~rps, set(gca,'YDir','reverse'); end
  
     % ICR
     subplot('Position',pos{2,2}); clear MMC; 
     M2 = (fasttest | ~Q.train) & Q.noise==1 & Q.bias==60  & Q.vx_int==0;                       MMC{1}=repmat(Q.ICR(M2),numel(method),1);   
     M2 = (fasttest | ~Q.train) & Q.noise==5 & Q.bias==60  & Q.vx_int==0;                       MMC{2}=repmat(Q.ICR(M2),numel(method),1);   
     M2 = (fasttest | ~Q.train) & Q.noise==9 & Q.bias==60  & Q.vx_int==0;                       MMC{3}=repmat(Q.ICR(M2),numel(method),1);   
     M2 = (fasttest | ~Q.train) & Q.noise>0  & Q.bias==60  & Q.vx_int==0 &  all(Q.vx_vol==1,2); MMC{4}=repmat(Q.ICR(M2),numel(method),1);   
     M2 = (fasttest | ~Q.train) & Q.noise>0  & Q.bias==60  & Q.vx_int==0 &  all(Q.vx_vol==2,2); MMC{5}=repmat(Q.ICR(M2),numel(method),1);   
     M2 = (fasttest | ~Q.train) & Q.noise>0  & Q.bias==60  & Q.vx_int==1;                       MMC{6}=repmat(Q.ICR(M2),numel(method),1);   
     if rps, for ci=1:numel(MMC), MMC{ci} = mark2rps(MMC{ci}); end; end
     cat_plot_boxplot(MMC,struct('names',{{'N1','N5','N9','R','rR','irR'}}, 'subsets',[0 0 0 1 1 1],...
       'style',4,'datasymbol','o','usescatter',1,'sort',0,'ylim',nlim,'groupnum',0,'ygrid',0)); set(gca,'FontSize',FS(2)); hold on
     set(gca,'ylim',nlim,'ytick',mlim,'FontSize',FS(2),'ygrid',1);
     title('60% bias (grad 2)','FontSize',FS(1),'FontWeight','bold'); 
     xlabel(sprintf('BWP noise (%d) / RES (%d-%d)',numel(MMC{1}),numel(MMC{4}),numel(MMC{5}))); ylabel(sprintf('ICR (%s)',rpsnam{rps+1}));
     if ~rps, set(gca,'YDir','reverse'); end
  
     % CON
     subplot('Position',pos{2,3}); clear MMC; 
     M2 = (fasttest | ~Q.train) & Q.noise>0  & Q.bias==60 & Q.vx_int==0 &  all(Q.vx_vol==1,2); MMC{1}=repmat(Q.contrastr(M2),numel(method),1);   
     M2 = (fasttest | ~Q.train) & Q.noise>0  & Q.bias==60 & Q.vx_int==0 &  all(Q.vx_vol==2,2); MMC{3}=repmat(Q.contrastr(M2),numel(method),1);   
     M2 = (fasttest | ~Q.train) & Q.noise>0  & Q.bias==60 & Q.vx_int==1;                       MMC{5}=repmat(Q.contrastr(M2),numel(method),1);   
     M2 = (fasttest | ~Q.train) & Q.noise==5 & Q.bias>0   & Q.vx_int==0 &  all(Q.vx_vol==1,2); MMC{2}=repmat(Q.contrastr(M2),numel(method),1);   
     M2 = (fasttest | ~Q.train) & Q.noise==5 & Q.bias>0   & Q.vx_int==0 &  all(Q.vx_vol==2,2); MMC{4}=repmat(Q.contrastr(M2),numel(method),1);   
     M2 = (fasttest | ~Q.train) & Q.noise==5 & Q.bias>0   & Q.vx_int==1;                       MMC{6}=repmat(Q.contrastr(M2),numel(method),1);   
     if rps, for ci=1:numel(MMC), MMC{ci} = mark2rps(MMC{ci}); end; end
     cat_plot_boxplot(MMC,struct('names',{{'N','B','rN','rB','irN','irB'}},'subsets',[0 0 1 1 0 0], ...
       'style',4,'datasymbol','o','usescatter',1,'sort',0,'ylim',nlim,'groupnum',0,'ygrid',0)); set(gca,'FontSize',FS(2)); hold on
     set(gca,'ylim',nlim,'ytick',mlim,'FontSize',FS(2),'ygrid',1);
     title('contrast','FontSize',FS(1),'FontWeight','bold'); 
     xlabel('BWP noise/bias/RES levels'); ylabel(sprintf('contrast (%s)',rpsnam{rps+1}));
     if ~rps, set(gca,'YDir','reverse'); end
  
     % ECR
     switch qafile
       case  'cat_vol_qa201901_202302'
         subplot('Position',pos{2,3}); clear MMC; 
         M2 = (fasttest | ~Q.train) & Q.noise>0  & Q.bias==60 & Q.vx_int==0 &  all(Q.vx_vol==2,2); MMC{1}=repmat(Q.ECR(M2),numel(method),1);   
         M2 = (fasttest | ~Q.train) & Q.noise>0  & Q.bias==60 & Q.vx_int==0 &  all(Q.vx_vol==1,2); MMC{3}=repmat(Q.ECR(M2),numel(method),1);   
         M2 = (fasttest | ~Q.train) & Q.noise>0  & Q.bias==60 & Q.vx_int==1 &  any(Q.vx_vol~=2,2); MMC{5}=repmat(Q.ECR(M2),numel(method),1);   
         M2 = (fasttest | ~Q.train) & Q.noise>0  & Q.bias==60 & Q.vx_int==1 &  all(Q.vx_vol==2,2); MMC{7}=repmat(Q.ECR(M2),numel(method),1);   
         M2 = (fasttest | ~Q.train) & Q.noise==5 & Q.bias>0   & Q.vx_int==0 &  all(Q.vx_vol==2,2); MMC{2}=repmat(Q.ECR(M2),numel(method),1);   
         M2 = (fasttest | ~Q.train) & Q.noise==5 & Q.bias>0   & Q.vx_int==0 &  all(Q.vx_vol==1,2); MMC{4}=repmat(Q.ECR(M2),numel(method),1);   
         M2 = (fasttest | ~Q.train) & Q.noise==5 & Q.bias>0   & Q.vx_int==1 &  any(Q.vx_vol~=2,2); MMC{6}=repmat(Q.ECR(M2),numel(method),1);   
         M2 = (fasttest | ~Q.train) & Q.noise==5 & Q.bias>0   & Q.vx_int==1 &  all(Q.vx_vol==2,2); MMC{8}=repmat(Q.ECR(M2),numel(method),1);   
        if rps, for ci=1:numel(MMC), MMC{ci} = mark2rps(MMC{ci}); end; end
        cat_plot_boxplot(MMC,struct('names',{{'r2N','r2B','r1N','r1B','ir1N','ir1B','ir2N','ir2B'}},'subsets',[0 0 1 1 0 0], ...
           'style',4,'datasymbol','o','usescatter',1,'sort',0,'ylim',nlim,'groupnum',0,'ygrid',0)); set(gca,'FontSize',FS(2)); hold on
         set(gca,'ylim',nlim,'ytick',mlim,'FontSize',FS(2),'ygrid',0);
         title('ECR at 5% noise','FontSize',FS(1),'FontWeight','bold'); 
         xlabel('BWP noise/bias/resolution levels'); ylabel(sprintf('ECR (%s)',rpsnam{rps+1}));
       otherwise
         subplot('Position',pos{2,4}); clear MMC; 
         M2 = (fasttest | ~Q.train) & Q.noise>0  & Q.bias==60 & Q.vx_int==0 &  all(Q.vx_vol==1,2); MMC{1}=repmat(Q.ECR(M2),numel(method),1);   
         M2 = (fasttest | ~Q.train) & Q.noise>0  & Q.bias==60 & Q.vx_int==0 &  all(Q.vx_vol==2,2); MMC{3}=repmat(Q.ECR(M2),numel(method),1);   
         M2 = (fasttest | ~Q.train) & Q.noise>0  & Q.bias==60 & Q.vx_int==1;                       MMC{5}=repmat(Q.ECR(M2),numel(method),1);   
         M2 = (fasttest | ~Q.train) & Q.noise==5 & Q.bias>0   & Q.vx_int==0 &  all(Q.vx_vol==1,2); MMC{2}=repmat(Q.ECR(M2),numel(method),1);   
         M2 = (fasttest | ~Q.train) & Q.noise==5 & Q.bias>0   & Q.vx_int==0 &  all(Q.vx_vol==2,2); MMC{4}=repmat(Q.ECR(M2),numel(method),1);   
         M2 = (fasttest | ~Q.train) & Q.noise==5 & Q.bias>0   & Q.vx_int==1;                       MMC{6}=repmat(Q.ECR(M2),numel(method),1);  
         if rps, for ci=1:numel(MMC), MMC{ci} = mark2rps(MMC{ci}); end; end
         cat_plot_boxplot(MMC,struct('names',{{'N','B','rN','rB','irN','irB'}},'subsets',[0 0 1 1 0 0], ...
           'style',4,'datasymbol','o','usescatter',1,'sort',0,'ylim',nlim,'groupnum',0,'ygrid',0)); set(gca,'FontSize',FS(2)); hold on
         set(gca,'ylim',nlim,'ytick',mlim,'FontSize',FS(2),'ygrid',1);
         title('ECR estimation','FontSize',FS(1),'FontWeight','bold'); 
         xlabel('BWP noise/bias/resolution levels'); ylabel(sprintf('ECR (%s)',rpsnam{rps+1}));
     end
     if ~rps, set(gca,'YDir','reverse'); end
  
     % FEC
     subplot('Position',pos{2,5}); clear MMC; 
     M2 = (fasttest | ~Q.train) & Q.noise==5 & Q.bias==20  & Q.vx_int==0;                       MMC{1}=repmat(Q.FEC(M2),numel(method),1);   
     M2 = (fasttest | ~Q.train) & Q.noise==5 & Q.bias==60  & Q.vx_int==0;                       MMC{2}=repmat(Q.FEC(M2),numel(method),1);   
     M2 = (fasttest | ~Q.train) & Q.noise==5 & Q.bias==100 & Q.vx_int==0;                       MMC{3}=repmat(Q.FEC(M2),numel(method),1);   
     M2 = (fasttest | ~Q.train) & Q.noise==5 & Q.bias>0    & Q.vx_int==0 &  all(Q.vx_vol==1,2); MMC{4}=repmat(Q.FEC(M2),numel(method),1);   
     M2 = (fasttest | ~Q.train) & Q.noise==5 & Q.bias>0    & Q.vx_int==0 &  all(Q.vx_vol==2,2); MMC{5}=repmat(Q.FEC(M2),numel(method),1);   
     M2 = (fasttest | ~Q.train) & Q.noise==5 & Q.bias>0    & Q.vx_int==1;                       MMC{6}=repmat(Q.FEC(M2),numel(method),1);  
     if rps, for ci=1:numel(MMC), MMC{ci} = mark2rps(MMC{ci}); end; end
     cat_plot_boxplot(MMC,struct('names',{{'B2','B6','B10','R','rR','irR'}}, 'subsets',[0 0 0 1 1 1],...
       'style',4,'datasymbol','o','usescatter',1,'sort',0,'ylim',nlim,'groupnum',0,'ygrid',0)); set(gca,'FontSize',FS(2)); hold on
     set(gca,'ylim',nlim,'ytick',mlim,'FontSize',FS(2),'ygrid',1);
     title('FEC at 5% noise (grad 3)','FontSize',FS(1),'FontWeight','bold'); 
     xlabel(sprintf('BWP bias (%d) /resolution levels (%d-%d)',numel(MMC{1}),numel(MMC{4}),numel(MMC{5}))); ylabel(sprintf('FEC (%s)',rpsnam{rps+1}));
     if ~rps, set(gca,'YDir','reverse'); end
    
  
   
  %%
  % ----------------------------------------------------------------------
  % C) relation between QM and kappa
  % ----------------------------------------------------------------------
    marklab = {'mark','percentage'};
    method = method; 
   
    % comments
    annotation('textbox',[0     0.965 0.04 0.04],'String','A','FontSize',FS(1)*1.5,'FontWeight','bold','EdgeColor','none')
    annotation('textbox',[0     0.715 0.04 0.04],'String','B','FontSize',FS(1)*1.5,'FontWeight','bold','EdgeColor','none')
    annotation('textbox',[0     0.465 0.04 0.04],'String','C','FontSize',FS(1)*1.5,'FontWeight','bold','EdgeColor','none')
    annotation('line'   ,[0 1],[0.76 0.76],'Color',repmat(0.2,1,3));
    annotation('line'   ,[0 1],[0.51 0.51],'Color',repmat(0.2,1,3));
  
    clear MX,
    MX{1} = M & Q.vx_int==0 &  all(Q.vx_vol==1,2); % 1mm not-interpolated 
    MX{2} = M & Q.vx_int==1;                       % 1mm interpolated
    MX{3} = M & Q.vx_int==0 & ~all(Q.vx_vol==1,2); % 2mm non-interpolated
    MXN    = {'1 mm','I(> 1 mm)','> 1 mm',};
    MXTL   = {'-' '-' '-','-'; 0.5 0.5 0.5 1};
    mrk    = {'o','^','v',''};
    kcol   = [0.2 0.5 0; 0 0.2 .8; .8 0 0]; 
  
    % (S)IQR
    if 1
      if strcmp(P.bwp(mi,1),'CAT12')
        ktx = {'kappa','SIQR',[0.725 0.975],[0 6]        ,[40 100];
               'kappa','GMV' ,[0.725 0.975],[0.425 0.485],[0.425 0.485];
               'GMV'  ,'SIQR',[0.425 0.485],[0 6]        ,[40 100]}; 
      elseif strcmp(P.bwp(mi,1),'SPM12')
        ktx = {'kappa','SIQR',[.55    0.95],[0 6]      ,[20 100];
               'kappa','GMV' ,[.55    0.95],[0.3 0.6],[0.30 0.6];
               'GMV'  ,'SIQR',[0.3  0.6],[0 6]      ,[20 100]}; 
      elseif strcmp(P.bwp(mi,1),'synthseg')
        ktx = {'kappa','SIQR',[.6    0.85],[0 6]      ,[40 100];
               'kappa','GMV' ,[.6    0.85],[0.42 0.53],[0.42 0.53];
               'GMV'  ,'SIQR',[0.42  0.53],[0 6]      ,[40 100]}; 
      end
    else % equal
      ktx = {'kappa','SIQR',[0.6 1],  [0 6]    ,[20 100];
             'kappa','GMV' ,[0.6 1],  [0.3 0.6],[0.3 0.6];
             'GMV'  ,'SIQR',[0.3 0.6],[0 6]    ,[20 100]}; 
    end
    for kti = 1:size(ktx,1)
      subplot('Position',pos{3,kti},'replace','box','on'); hold off; clear MMC
      for mi=method
        if method(mi)>0
          for mtxi = 1:numel(MX)
            if rps==0 || kti==2
              hs = scatter(Q.(ktx{kti,1})(MX{mtxi},mi), Q.(ktx{kti,2})(MX{mtxi},mi),14,mrk{mtxi}); hold on
            else
              hs = scatter(Q.(ktx{kti,1})(MX{mtxi},mi), mark2rps(Q.(ktx{kti,2})(MX{mtxi},mi)),14,mrk{mtxi}); hold on
            end
            set(hs,'markeredgecolor',kcol(mtxi,:),'markerfacecolor',kcol(mtxi,:),'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.4);
          end
        end
      end
      for mi=method
        if method(mi)>0
          for mtxi = 1:numel(MX)
            if rps==0 || kti==2
              Q.fit.methKapa{mi} = robustfit(Q.(ktx{kti,1})(MX{mtxi},mi),Q.(ktx{kti,2})(MX{mtxi},mi));
            else
              Q.fit.methKapa{mi} = robustfit(Q.(ktx{kti,1})(MX{mtxi},mi),mark2rps(Q.(ktx{kti,2})(MX{mtxi},mi)));
            end
            plot([0,1],[Q.fit.methKapa{mi}(1),Q.fit.methKapa{mi}(1) + Q.fit.methKapa{mi}(2) ],...
              MXTL{1,mtxi},'color',kcol(mtxi,:),'LineWidth',MXTL{2,mtxi});
          end
        end
      end
      hold off; 
      if (kti == 1 || kti == 3) && ~rps, set(gca,'YDir','reverse'); end
      switch kti
        case 1, set( gca, 'xdir', 'reverse' , 'YTick', 25:10:95, 'XTick', .55:.1:0.95, 'XTickLabelRotation',0);
        case 2, set( gca, 'xdir', 'reverse' , 'YTick', ktx{2,4}(1):round(diff(ktx{2,4})/6,2):ktx{2,4}(2), 'XTick', .55:.1:0.95, 'XTickLabelRotation',0);
        case 3, set( gca, 'xdir', 'reverse' , 'YTick', 25:10:95, 'XTick', ktx{3,4}(1):round(diff(ktx{2,4})/6,2):ktx{3,4}(2), 'XTickLabelRotation',0);
      end
      xlim(ktx{kti,3}); ylim(ktx{kti,4+rps}); 
      title(sprintf('%s vs. %s',ktx{kti,1},ktx{kti,2}),'FontSize',FS(1),'FontWeight','bold'); 
      xlabel(ktx{kti,1}); ylabel(sprintf('%s (%s)',ktx{kti,2},rpsnam{rps+1})); box on; grid on; legend off; 
      legend(MXN(1:3),'Location','Southwest','FontSize',6);
    end
  
    
    
    % -- correlations ---------------------------------------------------
    cat_io_cmd('  Estimate correlations:','g5','',1); 
    %M  = Q.noise>0 & Q.bias>0; 
    M2 = repmat(M,1,size(Q.kappa,2)) & repmat(method,size(Q.noise,1),1);
    QMnames = {'noise','bias','resRMS','NCR','ICR','RES','ECR','FEC','IQR','SIQR','Kappa','rCSFV','rGMV','rWMV'};
    QMcorrs = [ 
      repmat(Q.noise(M),numel(method),1), ...
      double(repmat(Q.bias(M),numel(method),1)) ,...
      Q.RESgt(M2), ...
      Q.NCR(M2) , Q.ICR(M2) , Q.res_RMS(M2) , Q.ECR(M2) , Q.FEC(M2), Q.IQR(M2), Q.SIQR(M2), ...
      Q.kappa(M2) , ...
      Q.CMV(M2),Q.GMV(M2),Q.WMV(M2)];
  
    [C.QM.r, C.QM.p] = corr(QMcorrs,'type',corrtype);
    C.QM.rtab = [ [{''} QMnames]; [QMnames',num2cell(round(C.QM.r*10000)/10000)] ];
    C.QM.ptab = [ [{''} QMnames]; [QMnames',num2cell(C.QM.p)] ];
    C.QM.etab = cell(1,size(C.QM.rtab,2)); C.QM.xtab =  C.QM.etab;
    C.QM.xtab{1} = sprintf('%s (r/p-value)',corrtype); 
    
    if isfield(C.QM,'txt'), C.QM = rmfield(C.QM,'txt'); end
    T.ixiPPC = sprintf('\n%s correlation on BWP:\n%s\n',corrtype,...
      '------------------------------------------------------------------------------');
    for di=2:size(C.QM.rtab,2)
      T.ixiPPC = sprintf('%s%8s   ',T.ixiPPC,C.QM.rtab{1,di}); 
      C.QM.xtab{1,di} = C.QM.rtab{1,di}; 
    end 
    for dj=2:size(C.QM.rtab,1)
      T.ixiPPC = sprintf('%s\n%8s   ',T.ixiPPC,C.QM.rtab{1,dj});
      C.QM.xtab{dj,1} = C.QM.rtab{dj,1}; 
      for di=2:size(C.QM.rtab,2)
        if     C.QM.ptab{dj,di}<0.000001, star='***';
        elseif C.QM.ptab{dj,di}<0.000100, star='** ';
        elseif C.QM.ptab{dj,di}<0.010000, star='*  ';
        else                              star='   ';
        end
        if di<dj
          T.ixiPPC = sprintf('%s%8.4f%3s',T.ixiPPC,C.QM.rtab{dj,di},star);
          C.QM.xtab{dj,di} = sprintf('%0.8e',C.QM.ptab{dj,di});
        else
          T.ixiPPC = sprintf('%s%8s%3s',T.ixiPPC,'','');
          C.QM.xtab{dj,di} = C.QM.rtab{dj,di};
        end
      end
    end
    fname = fullfile(opt.resdir,sprintf('tst_bwpmaintest_%s_rptable_%s.csv',qafile,segment{mi})); 
    cat_io_csv(fname, C.QM.xtab); 
    cat_io_cprintf('blue',sprintf('\n    Write %s\n',fname)); 
    

    C.QM.txt = T.ixiPPC;
    C.QM.txt = [C.QM.txt '\n\n\n'];
    fc = fullfile(opt.resdir,sprintf('tst_bwpmaintest_%s_%s.txt',qafile,segment{mi})); 
    f  = fopen(fc,'w'); if f~=-1, fprintf(f,C.QM.txt); fclose(f); end
    cat_io_cprintf('blue','    Write %s\n',fc); 
  
    
    
    
    
    % Kappa QM vs QMgt ... compare QR to the rated BWP values (method table)
    T.Kappa.rtab(1,:)  = {'method','NCR','ICR','RES','ECR','FEC','IQR','SIQR'};
    T.Kappa.ptab(1,:)  = {'method','NCR','ICR','RES','ECR','FEC','IQR','SIQR'};
    QMfields = {'NCR','ICR','res_RMS','res_ECR','FEC','IQR','SIQR'}; Q.res_RMSgt = Q.res_RMS;
    for mi=method
      T.Kappa.rtab{mi+1,1} = sprintf('Kappa %s',P.bwp{method(mi),1});
      T.Kappa.ptab{mi+1,1} = sprintf('Kappa %s',P.bwp{method(mi),1});
      QMcorrs = mark2rps([ Q.NCR(M2(:,1),method(mi)) , Q.ICR(M2(:,1),method(mi)) , ...
        Q.res_RMS(M2(:,1),method(mi)) , Q.ECR(M2(:,1),method(mi)) , Q.FEC(M2(:,1),method(mi)) , ...
        Q.IQR(M2(:,1),method(mi)) , Q.SIQR(M2(:,1),method(mi)) ]);
      [C.QM.r, C.QM.p] = corr(QMcorrs,Q.kappa(M2(:,1),method(mi)),'type',corrtype);
      for fi=1:numel(QMfields)
        T.Kappa.rtab{mi+1,fi+1} = C.QM.r(fi); 
        T.Kappa.ptab{mi+1,fi+1} = C.QM.p(fi);
      end
    end
    mi=mi+1;
    for fi=1:numel(QMfields)
      T.Kappa.rtab{mi+1,fi+1} = mean([T.Kappa.rtab{2:mi,fi+1}]); 
      T.Kappa.ptab{mi+1,fi+1} = 1; 
      T.Kappa.rtab{mi+1,1}    = 'Kappa mean';
      T.Kappa.rtab{mi+2,fi+1} = std([T.Kappa.rtab{2:mi,fi+1}]); 
      T.Kappa.ptab{mi+2,fi+1} = 1; 
      T.Kappa.rtab{mi+2,1}    = 'Kappa std';
    end
    mi=mi+3;
    T.Kappa.rtab{mi,1}    = 'Kappa all';
     QMcorrs = mark2rps([ Q.NCR(M2) , Q.ICR(M2) , Q.res_RMS(M2) , Q.ECR(M2) , Q.FEC(M2) , Q.IQR(M2) , Q.SIQR(M2) ]);
    [C.QM.r, C.QM.p] = corr(QMcorrs,Q.kappa(M2),'type',corrtype);
    for fi=1:numel(QMfields)
      T.Kappa.rtab{mi+1,fi+1} = C.QM.r(fi); 
      T.Kappa.ptab{mi+1,fi+1} = C.QM.p(fi);
    end
   
    
    if isfield(C.QM,'txt'), C.QM = rmfield(C.QM,'txt'); end
    T.ixiPPC = sprintf('\n%s on BWP:\n%s\n',corrtype,...
      '------------------------------------------------------------------------------');
    for di=1:size(T.Kappa.rtab,2), T.ixiPPC = sprintf('%s%8s   ',T.ixiPPC,T.Kappa.rtab{1,di}); end 
    for dj=2:size(T.Kappa.rtab,1)
      T.ixiPPC = sprintf('%s\n%12s ',T.ixiPPC,T.Kappa.rtab{dj,1});
      for di=2:size(T.Kappa.rtab,2)
        if     T.Kappa.ptab{dj,di}<0.000001, star='***';
        elseif T.Kappa.ptab{dj,di}<0.000100, star='** ';
        elseif T.Kappa.ptab{dj,di}<0.010000, star='*  ';
        else                                 star='   ';
        end
        T.ixiPPC = sprintf('%s%8.3f%3s',T.ixiPPC,T.Kappa.rtab{dj,di},star);
      end
    end
    C.QM.txt = T.ixiPPC;
    C.QM.txt = [C.QM.txt sprintf('\n\n\n')];
    mi = 1; 
    fc = fullfile(opt.resdir,sprintf('tst_bwpmaintest_%s_corr_%s.csv',qafile,segment{mi})); 
    f  = fopen(fc,'w'); if f~=-1, fprintf(f,C.QM.txt); fclose(f); end
    cat_io_cprintf('blue','    Write %s\n',fc); 
  
    
    
    
    % RMSE ... compare QR to the rated BWP values
    T.RMS       = {'method','NCR','ICR','ECR','RES'    ,'ECR','FEC','IQR','SIQR'}; 
    QMfields    = {         'NCR','ICR','ECR','res_RMS','ECR','FEC','IQR','SIQR'}; Q.res_RMSgt = Q.res_RMS; 
    MZ          = Q.noise>0 & Q.bias>0;
    for mi=1:numel(method)
      T.RMS{mi+1,1} = sprintf('RMSE %s',P.bwp{method(mi),1});
      for fi=1:numel(QMfields)
        T.RMS{mi+1,fi+1} = rms(mark2rps(Q.(QMfields{fi})(MZ,method(mi))) - ...
                               mark2rps(Q.([QMfields{fi} 'gt'])(MZ,1)),2); 
      end
    end
    if 0 %numel( methods ) > 1
      T.RMS{mi+numel(method),1} = sprintf('RMSE all');
      for fi=1:numel(QMfields)
        T.RMS{mi+numel(method),fi+1} = rms(reshape(mark2rps(Q.(QMfields{fi})(:,method)),1,numel(method)*size(Q.(QMfields{fi}),1)) - ...
                            reshape(mark2rps(repmat(Q.([QMfields{fi} 'gt'])(:,1),1,numel(method))),1,numel(method)*size(Q.(QMfields{fi}),1)),2); 
      end
    end
    % RMSE ... compare QR to the rated BWP values
    T.RMShq       = T.RMS(1,:);
    MH            = Q.noise>0 & Q.noise<5 & Q.bias>0; % & Q.bias<120; %& all(Q.vx_vol==1,2);
    for mi=1:numel(method)
      T.RMShq{mi+1,1} = sprintf('RMSE %s',P.bwp{method(mi),1});
      for fi=1:numel(QMfields)
        T.RMShq{mi+1,fi+1} = rms(mark2rps(Q.(QMfields{fi})(MH,method(mi))) - ...
                                 mark2rps(Q.([QMfields{fi} 'gt'])(MH,1)),2); 
      end
    end
    T.RMSlq       = T.RMS(1,:); 
    ML            = Q.noise>5 & Q.bias>=0; % & Q.field=='A';
    for mi=1:numel(method)
      T.RMSlq{mi+1,1} = sprintf('RMSE %s (%d)',P.bwp{method(mi),1});
      for fi=1:numel(QMfields)
        T.RMSlq{mi+1,fi+1} = rms(mark2rps(Q.(QMfields{fi})(ML,method(mi))) - ...
                                 mark2rps(Q.([QMfields{fi} 'gt'])(ML,1)),2); 
      end
    end
    if numel( method ) > 1
      T.RMS(mi+numel(method)+0,:) = [sprintf('RMSE mean (N=%d)',sum(MZ)),    num2cell(mean(cell2mat(T.RMS(2:end,2:end))))];
      T.RMS(mi+numel(method)+1,:) = [sprintf('RMSE std  (N=%d)',sum(MZ)),    num2cell(std( cell2mat(T.RMS(2:end,2:end))))];
    end
    T.RMS(mi+numel(method)+(numel(method)>1)*2+0,:) = [sprintf('RMSE HQ (N=%d)',sum(MH)), num2cell(mean(cell2mat(T.RMShq(2:end,2:end)),1))];
    T.RMS(mi+numel(method)+(numel(method)>1)*2+1,:) = [sprintf('RMSE LQ (N=%d)',sum(ML)), num2cell(mean(cell2mat(T.RMSlq(2:end,2:end)),1))];
    T.RMS(mi+numel(method)+(numel(method)>1)*2+2,:) = [sprintf('RMSE ALL(N=%d)',sum(MZ)), num2cell(mean(cell2mat(T.RMS(2:end,2:end)),1))];
    fc = fullfile(opt.resdir,sprintf('tst_bwpmaintest_%s_rms_%s.csv',qafile,segment{mi}));
    cat_io_csv(fc,T.RMS)
    cat_io_cprintf('blue','    Write %s\n',fc); 
    
  
    % -- save images ---------------------------------------------------  
    mi = 1; 
    cat_io_cmd('  Print figures:','g5','',1);
    fc = fullfile(opt.resdir,sprintf('tst_bwpmaintest_%s_%s.csv',qafile,segment{mi}));
    f  = fopen(fc,'w'); if f~=-1, fprintf(f,C.QM.txt); fclose(f); end
    cat_io_cprintf('blue','\n    Write %s\n',fc); 
  
    fc = fullfile(opt.resdir,sprintf('fig1_bwpmaintest_%s_%s',qafile,segment{mi}));
    print(fh1,fc,opt.res,'-dpng');
    cat_io_cprintf('blue','    Save  %s.png\n',fc); 
    if opt.closefig, close(fh1); end


    

    %% -- NCR vs. CNR  --------------------------------------------------- 
    % Based on a reviewer request here the comparison between the NCR vs. 
    % CNR definition. 
    % The NCR support a more linear scaling with lower/higher variance for 
    % good/bad ratings (lower values for good quality), whereas the CNR shows 
    % the oposit pattern, ie, higher/lower variance for good/bad values with
    % higher values for good quality.
    % The NCR focuses therefore more on finer separation of low quality data
    % whereas the CNR focuses more on further separation of high quality data. 
    % The relevance becomes more clear together with the Kappa ratings, 
    % where the CNR presents a clear linear relation.
    % However, the NCR is maybe more helpful in separating ultra-high
    % resolution data and is less depening on the RMS rating concept. 

    cat_io_cmd('  Print NCR vs. CNR figure:','g5','',1); fprintf('\n'); 
    
    % get raw measures for NCR and CNR
    for mi = 1
      for pi = 1:numel(X(mi).xml)
        Q2.NCR(pi,mi) = X(mi).xml(pi).qualityratings.NCR;
        Q2.CNR(pi,mi) = X(mi).xml(pi).qualitymeasures.contrast ./ ...
          (X(mi).xml(pi).qualitymeasures.NCR*X(mi).xml(pi).qualitymeasures.contrast);
   
        Q2.ICR(pi,mi) = X(mi).xml(pi).qualityratings.ICR;
        Q2.CIR(pi,mi) = X(mi).xml(pi).qualitymeasures.contrast ./ ...
          (X(mi).xml(pi).qualitymeasures.ICR*X(mi).xml(pi).qualitymeasures.contrast);
      end
    end
    
    corrNCR = [ corr(Q2.NCR,Q2.CNR,'Type','Spearman'), corr(Q2.NCR,Q2.CNR,'Type','Pearson') ];
    corrICR = [ corr(Q2.ICR,Q2.CIR,'Type','Spearman'), corr(Q2.ICR,Q2.CIR,'Type','Pearson') ];
    
    corrNCRkappa = [ corr(Q2.NCR,Q.kappa,'Type','Spearman'), corr(Q2.NCR,Q.kappa,'Type','Pearson') ]; 
    corrCNRkappa = [ corr(Q2.CNR,Q.kappa,'Type','Spearman'), corr(Q2.CNR,Q.kappa,'Type','Pearson') ]; 
    

    for figSize = 2%1:2
      fh1 = figure(333);
      if figSize==1, fh1.Position(3:4) = [1200 400]; else fh1.Position(3:4) = [750 250]; end
      clf(fh1); fh1.Visible = 'on'; 

      MM = Q.RESgt==2; % take only 1 mm to make it easier
      biases = unique(Q.bias); biascolors = cool(numel(biases)); 
      noises = unique(Q.noise); noisemarker = 'o<>^v'; lg = {}; 
      tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact'); 
      for ti = 1:3
        nexttile; hold on; 
              
        for nf = 1:numel(noises)
          for bf = 1:numel(biases)
            MMM = MM & Q.bias == biases(bf) & Q.noise == noises(nf); 
            switch ti
              case 1, sc(bf) = scatter(Q2.NCR(MMM),Q2.CNR(MMM),'filled');
              case 2, sc(bf) = scatter(Q2.NCR(MMM),Q.kappa(MMM),'filled');
              case 3, sc(bf) = scatter(Q2.CNR(MMM),Q.kappa(MMM),'filled');
            end
            sc(bf).Marker          = noisemarker(nf);
            sc(bf).MarkerFaceAlpha = .8; 
            sc(bf).MarkerEdgeAlpha = .0; 
            sc(bf).SizeData        = 50; 
            sc(bf).MarkerEdgeColor = biascolors(bf,:); 
            sc(bf).MarkerFaceColor = sc(bf).MarkerEdgeColor;
            lg = [lg; sprintf('nf=%0.0f%% bf=%0.0f%%',noises(nf),biases(bf))]; 
          end
        end
        box on; grid on; 
  
        switch ti
          case 1 
            title('CNR vs. NCR')
            ylabel('worse  << \bf NCR \rm >>  better')
            xlabel('better  << \bf CNR (mark) \rm >>  worse'); 
            if figSize==2, set(gca,'XTick',1:5); xlim([.5 5.5]); end 
            if figSize==1 || fasttest, legend(lg,'Location','Northeast','FontSize',7); end
          case 2 
            title('CNR vs. Kappa')
            ylabel('worse  << \bf Kappa \rm >>  better')
            xlabel('better  << \bf CNR (mark) \rm >>  worse'); 
            if figSize==2, set(gca,'XTick',1:5); xlim([.5 5.5]); end
            ylim([.84 .94]+0.01*(figSize==1)); 
          case 3 
            title('NCR vs. Kappa')
            ylabel('worse  << \bf Kappa \rm >>  better')
            xlabel('worse  << \bf NCR \rm >>  better'); 
            ylim([.84 .94]+0.01*(figSize==1));
            set(gca,'XTick',0:20:80);
        end
        %subtitle('(1 mm resolution with 3 bias fields per group)'); 
        %ax = gca; ax.YScale = 'log';
      end
      fc = fullfile(opt.resdir,sprintf('fig2_bwpmaintest_%s_rps%0.0f_CNRvsNCR%d_%s',qafile,rps,figSize,segment{mi})); 
      print(fh1,fc,opt.res,'-dpng');
      cat_io_cprintf('blue','    Save  %s.png\n',fc); 
      if opt.closefig, close(fh1); end
    end
    if opt.closefig, close(fh1); end

  end
    
  end  
end

fprintf('BWP done.\n')

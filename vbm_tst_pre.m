function vbm_tst_pre
% ______________________________________________________________________
% 
% - spezielle Testfälle zum auswählen???
% - methode ändert ihre nummer... problem???
% - bias:
%     original | method1 method2 ...
% ______________________________________________________________________
% Robert Dahnke 2013_01
% Structural Brain Mapping Group
% University Jena
%  
% $Id$
% ______________________________________________________________________
%#ok<*TRYNC>

  warning off; clear classes; warning on; 

    
  if ~exist('opt','var'), opt=struct(); end
  
% TEST PARAMETER
% ----------------------------------------------------------------------
  %def.results.printCase = 0; % [0=none,1=only vbm12,2=all]
  def.subStepSize = 1; % use every xth subject
  def.path        = '/Volumes/MyBook/MRData/vbm_tst';
  def.printspace  = [50 7 2]; % [name,columnwidth,fractonal digits] 
  def.SPM8path    = '/Users/dahnke/Neuroimaging/SPM8R4290_VBM8';
  def.VBM12ppath  = '/Users/dahnke/Neuroimaging/SPM8R4290_VBM12+';
  def.SPM12path   = '/Users/dahnke/Neuroimaging/spm12';
  def.RAWdir      = '+RAW';
  def.CSVdir      = '+CSV';
  def.Tdir        = 'T1';
  def.GTdir       = 'SEG';
  def.SRSdir      = 'SRS';
  def.initPATH    = [ ...
    'export FSLDIR=/Applications/FSL5.0; sh ${FSLDIR}/etc/fslconf/fsl.sh;' ...
    'export FREESURFER_HOME=/Applications/Freesurfer5.1; ${FREESURFER_HOME}/SetUpFreeSurfer.sh;' ...
    'MNI_DIR=${FREESURFER_HOME}/mni; ' ...
    'PATH=$PATH:$FSLDIR/bin:$MNI_DIR/bin; ' ...
    'source $FREESURFER_HOME/FreeSurferEnv.sh; ' ...
    'export FSLOUTPUTTYPE=NIFTI NIFTIexport=NIFTI; ' ...
    ];
  
  % set of measures for one talbe
  def.measure.BE    = 'dice'; % {'dice','jaccard','kappa', ...}
  def.measure.SEG   = 'kappa';
  def.measure.compress = 0 ;
  def.FailedTh      = 0.5;
  def.csv.delimiter = ',';
  def.csv.komma     = '.';
  
  def.stableMethods = {'SPM12nc','VBM12i'};
  def.stableRecalc  = 0;  
  opt.reprint       = 0;

  replaceSPMpath(def.SPM8path,def.SPM12path);
  
% TEST METHODS
% ----------------------------------------------------------------------
  def.methods   = { 
  % 'name'    NC segments t1corr recalc                       % new
  %'VBM12+'  1 {'pb' 'p0'} {'mc' 'mc'} {0.5  0.5}  0 % internal
%  'VBM12+'  1 {'p0'}      {'mg'}      {0.5}       0 % internal
%     'VBM12i'  0 {'p0'}      {'m'}       {1.8}       1 % internal noise correction
%     'SPMnc'   1 {'p0'}      {''}        {0.5}       0 
%     'SPM'     0 {'p0'}      {''}        {0.5}       0 
%    'SPM8'    0 {'p0'}      {''}        {0.5}       0 
%    'SPM12'   0 {'p0'}      {''}        {0.5}       1 
     'SPM8nc'  1 {'p0'}      {''}        {0.5}       0 
%     'SPM12nc' 1 {'p0'}      {''}        {0.5}       1 
%    'FSLncN3'   1 {'p0'}      {''}        {-inf}      1 
%    'FSL'     0 {'p0'}      {''}        {-inf}      1 
%    'VBM8'    0 {'p0'}      {'m'}       {0.5}       0 % internal noise correction
%     'VBM12'   0 {'p0'}      {'m'}       {1.5}       1 % internal noise correction
%    'N3'      0 {''}        {'m'}       {-inf}       0 
%    't1qa'   1 {'pa'}      {'mc'}      {0.80}       0
  };


% TEST DATASETS
% ----------------------------------------------------------------------
  % datasets with ground truth
  def.subdirs = {
%  'BWPC_noise'
%  'BWPC_bias'
%   'BWPC_resi'
%   'BWPC_resr'
%   'QA_good'
%   'QA_bad'
%   'BWPC_NIR'
%  'BWP_Collins'
%  'BWP_Collins_T2'
%  'BWP_Collins_PD'
%     'ADHD'
%     'ADNI'
%     'IXI'
%     'NIH'
%     'OASIS'
 %  'IBSRv1'
 %  'IBSRv2'
%   'SVE_LPBA40'
%      'BWP_20N'
%      'Tumorbase'
%      'private'
%     'BO'
%    'BO_inverse'     
%      'private_full'
%     'BWP_Collins'
 %   'SRS'
 %      'Apes'
   'Apes_uthscsa'
   };
 

  opt = checkinopt(opt,def);
  
  
  % =====
  
  % subscans ... sieh mal nach deinen anderen alten testskripten 
  
  % =====
  
 % ab hier wäre eine subfunction vielleicht hilfreich??? 
 for di=1:numel(opt.subdirs)
    fprintf('\n%s:\n',opt.subdirs{di})
   
   
    % find original files of each test subdir without NC files
    opt.RAW.dirs{di} = fullfile(opt.path,opt.RAWdir,def.subdirs{di});
    opt.RAW.T{di}    = vbm_findfiles(fullfile(opt.RAW.dirs{di},opt.Tdir),'*.nii');
    sanlmfiles = ~cellfun('isempty',strfind(opt.RAW.T{di} ,'sanlm'));
    opt.RAW.T{di}(sanlmfiles) = []; clear sanlmfiles;
    
    
    % ground truth (one file ore similar files like T) and 
    % STABLE ground truth (one file for each T file, but maybe incomplete)
    if ~exist(fullfile(opt.RAW.dirs{di},opt.GTdir),'dir')
      mkdir(fullfile(opt.RAW.dirs{di},opt.GTdir)); 
    end
    opt.RAW.pmT{di}  = vbm_findfiles(fullfile(opt.RAW.dirs{di},opt.GTdir),'pm*.nii');
    opt.RAW.p0T{di}  = vbm_findfiles(fullfile(opt.RAW.dirs{di},opt.GTdir),'p0*.nii');
    opt.RAW.m0T{di}  = vbm_findfiles(fullfile(opt.RAW.dirs{di},opt.GTdir),'m0*.nii');   
    for fi=1:opt.subStepSize:numel(opt.RAW.T{di})
      [pp,ff,ee]=fileparts(opt.RAW.T{di}{fi}); 
      opt.RAW.psT{di}{fi} = fullfile(opt.RAW.dirs{di},opt.GTdir,['ps' ff ee]);
      opt.RAW.msT{di}{fi} = fullfile(opt.RAW.dirs{di},opt.GTdir,['ms' ff ee]);
    end
    
    
    % noise correction of the original files, if it was not done before
    for fi=1:opt.subStepSize:numel(opt.RAW.T{di})
      [pp,ff,ee]=fileparts(opt.RAW.T{di}{fi}); 
      opt.RAW.nT{di}{fi} = fullfile(pp,['sanlm_' ff ee]);
      if ~exist(opt.RAW.nT{di}{fi},'file')
        fprintf('sanlm %s\n',opt.RAW.nT{di}{fi});
        vbm_vol_sanlm(struct('data',opt.RAW.T{di}{fi},'rician',0)); 
      end
    end
    
    
    % XML-file for qa values
    for fi=1:opt.subStepSize:numel(opt.RAW.T{di})
      [pp,ff]=fileparts(opt.RAW.T{di}{fi});
      opt.RAW.XML{di}{fi} = fullfile(pp,['vbm_' ff '.xml']);
    end
    
    
    
 

    
    
    % prepare test directories
    % ------------------------------------------------------------------
    % Now we get from the main RAW directory to the method test
    % directories. We have to copy the (noise-corrected) files.
    % ------------------------------------------------------------------
    for mi=1:size(opt.methods,1)
      
      % copy the files to a method specific test dir
      opt.method(mi).name     = opt.methods{mi,1};
      opt.method(mi).dirs{di} = fullfile(opt.path,opt.methods{mi,1},...
        strrep(opt.RAW.dirs{di},fullfile(opt.path,def.RAWdir),''));
      
      
      % test dir creation and total clean up
      if exist(opt.method(mi).dirs{di},'dir') && opt.methods{mi,6}==2
        rmdir(opt.method(mi).dirs{di},'s'); 
      end
      if ~exist(opt.method(mi).dirs{di},'dir')
        mkdir(opt.method(mi).dirs{di});
      end
      
      
      % copy files, if there is now file in the result dir 
      % for spm the corregistration will change the header, but I expect
      % no further interactions - otherwise a total clear up can remove
      % the whole directory. This is a better way to be sure that all 
      % files are estimated with the actual version of the method.
      for fi=1:opt.subStepSize:numel(opt.RAW.T{di})
        [pp,ff,ee] = fileparts(opt.RAW.T{di}{fi});
        ppn = fullfile(opt.path,opt.methods{mi,1},...
          strrep(strrep(pp,fullfile(opt.path,def.RAWdir),''),[filesep 'T1'],''));
        if ~exist(ppn,'dir'), mkdir(ppn); end
        opt.method(mi).T{di}{fi} = fullfile(ppn,[ff ee]);
        if ~exist(opt.method(mi).T{di}{fi},'file')
          if opt.methods{mi,2}
            linkfile(opt.RAW.nT{di}{fi},opt.method(mi).T{di}{fi},'f');
          else
            linkfile(opt.RAW.T{di}{fi} ,opt.method(mi).T{di}{fi},'f');
          end
        end
        % and the mat file
        opt.RAW.Tmat{di}{fi} = vbm_findfiles(pp,[ff '*.mat'],struct('chararr',1)); 
        opt.method(mi).Tmat{di}{fi} = vbm_findfiles(ppn,[ff '*.mat'],struct('chararr',1)); 
        if ~isempty(opt.method(mi).Tmat{di}{fi}) && exist(opt.method(mi).Tmat{di}{fi},'file')
          linkfile(opt.RAW.Tmat{di}{fi},opt.method(mi).T{di}{fi},'f');
        end
      end
  
      
      % Results files can vary a little bit, to allow substeps like for
      % VBM12 with different segmenations (1) paT, (2) pbT, and (3) p0T.
      for ii=1:numel(opt.methods{mi,3})
        for fi=1:opt.subStepSize:numel(opt.method(mi).T{di})
          [pp,ff,ee] = fileparts(opt.method(mi).T{di}{fi});
          opt.method(mi).([def.methods{mi,3}{ii} 'T']){di}{fi} = ...
            fullfile(pp,[def.methods{mi,3}{ii} ff ee]);
          opt.method(mi).([def.methods{mi,4}{ii} 'T']){di}{fi} = ...
            fullfile(pp,[def.methods{mi,4}{ii} ff ee]);
          opt.method(mi).XML{di}{fi} = fullfile(pp,['vbm12tst_' ff '.xml']);
        end
      end
      clear pp ff ee;
    end    
    
    
    % run preprocessing
    % ------------------------------------------------------------------
    % Run test if the last result entry is empty or if a recalc is set
    % for this method - the one we develope =) 
    % This is this the most costly part and may require grid computing
    % or use of the paralel processing toolbox.
    % ------------------------------------------------------------------
    for mi=1:size(opt.methods,1)  
      for fi=1:opt.subStepSize:numel(opt.method(mi).T{di})
        %fi=21; %17;%20;
        if opt.methods{mi,6}>0 || ...
            (~isempty(def.methods{mi,3}{1}) && ~exist(opt.method(mi).([def.methods{mi,3}{1} 'T']){di}{fi},'file')) || ...
            (~isempty(def.methods{mi,4}{1}) && ~exist(opt.method(mi).([def.methods{mi,4}{1} 'T']){di}{fi},'file'))
          
          starttime = clock; 
          fname = spm_str_manip(opt.method(mi).T{di}{fi},sprintf('a%d',opt.printspace(1)));
          fprintf(sprintf('%%%ds:',opt.printspace(1)),fname(1:min(numel(fname),opt.printspace(1))));
          
          try 
            switch opt.methods{mi,1}
              case 't1qa',
                try
                  vbm_vol_intscale(opt.method(mi).T{di}{fi},struct('verb',0,'write_maT',0));
                catch %#ok<CTCH>
                  [pp,ff,ee]=spm_fileparts(opt.method(mi).T{di}{fi});
                  createNullImage(fullfile(pp,[ff ee]),fullfile(pp,['pa' ff ee]));
                  createNullImage(fullfile(pp,[ff ee]),fullfile(pp,['mc' ff ee]));
                end
              case 'VBM12+', 
                try
                  cg_vbm12_test_rd(opt.method(mi).T{di}{fi},struct('debug',5,'verb',0));
                catch %#ok<CTCH>
                  [pp,ff,ee]=spm_fileparts(opt.method(mi).T{di}{fi});
                  if ~exist(fullfile(pp,['pb' ff ee]),'file'); 
                    createNullImage(fullfile(pp,[ff ee]),fullfile(pp,['pb' ff ee])); end
                  if ~exist(fullfile(pp,['p0' ff ee]),'file'); 
                    createNullImage(fullfile(pp,[ff ee]),fullfile(pp,['p0' ff ee])); end
%                   if ~exist(fullfile(pp,['mb' ff ee]),'file'); 
%                     createNullImage(fullfile(pp,[ff ee]),fullfile(pp,['mb' ff ee])); end
                  if ~exist(fullfile(pp,['mc' ff ee]),'file'); 
                    createNullImage(fullfile(pp,[ff ee]),fullfile(pp,['mc' ff ee])); end
                  sprintf('ERR\n');
                end
              case 'VBM8',              VBM8segment(opt.method(mi).T{di}{fi},opt.SPM8path,opt.SPM12path);
              case {'VBM12'},           VBM12segment(opt.method(mi).T{di}{fi},opt.SPM12path,opt.SPM12path,0);
              case {'VBM12i'},          VBM12segment(opt.method(mi).T{di}{fi},opt.SPM12path,opt.SPM12path,1);
              case {'SPM12','SPM12nc'}, SPM12segment(opt.method(mi).T{di}{fi},opt.SPM12path,opt.SPM12path);
              case {'SPM8','SPM8nc'},   SPM8newsegment(opt.method(mi).T{di}{fi},opt.SPM8path,opt.SPM12path);
              case {'SPM','SPMnc'},     SPM8segment(opt.method(mi).T{di}{fi},opt.SPM8path,opt.SPM12path);
              case {'FSL','FSLncN3'},   FSL(opt.method(mi).T{di}{fi},opt.initPATH);
              case 'FSLnc',             FSL(opt.method(mi).T{di}{fi},opt.initPATH,opt.RAW.psT{di}{fi});
              case {'N3','N3nc'},       N3(opt.method(mi).T{di}{fi},opt.initPATH);
              otherwise
            end
            fprintf(' %4.0fs\n',etime(clock,starttime));
          catch %#ok<CTCH>
            fprintf(' %4.0fs - ERROR\n',etime(clock,starttime));
          end
        end
      end
    end
    replaceSPMpath(def.SPM8path,def.SPM12path);  
    
    
      
    %% phantom generation 
    %  -----------------------------------------------------------------
    fprintf('Stable:          ');
    for fi=1:opt.subStepSize:numel(opt.method(1).T{di})
      stableFiles      = cell(size(opt.stableMethods)); 
      stableFilesExist = zeros(size(opt.stableMethods));

      [pp,ff] = spm_fileparts(opt.method(1).T{di}{fi});
      gtdir = fullfile(opt.RAW.dirs{di},opt.GTdir);
      if ~strcmp(spm_str_manip(pp,'t'),opt.subdirs{di})
        fnameparts = eval(['{'' ' strrep( opt.method(1).T{di}{fi},'/',''',''') '''}'])';
        fnameparts([1,end])=[]; 
        subdirpos = find(cellfun('isempty',strfind(fnameparts,opt.subdirs{di}))==0)+1;
        for sdp=subdirpos:numel(fnameparts), 
          gtdir = fullfile(gtdir,fnameparts{sdp});
        end
        if ~exist(gtdir,'dir'), mkdir(gtdir); end
      end
      
      opt.RAW.psT{di}{fi} = fullfile(gtdir,['ps' ff '.nii']);
      %opt.RAW.p0T{di}{fi} = fullfile(gtdir,['p0' ff '.nii']);
 
      for pfi=1:numel(opt.stableMethods)
        stableFiles{pfi}      = fullfile(def.path,opt.stableMethods{pfi},opt.subdirs{ii},['p0' ff '.nii']);
        stableFilesExist(pfi) = exist(stableFiles{pfi},'file');
      end    
      
      % run stable, if possible and necessary
      if (opt.stableRecalc || ~exist(opt.RAW.psT{di}{fi},'file')) && all(stableFilesExist)
        vbm_tst_staple_multilabels(char(stableFiles'),'',opt.RAW.psT{di}{fi});
        
        fprintf(sprintf('%s%%4d/%%4d',sprintf('%s',repmat('\b',1,9))),fi,numel(opt.RAW.T{di}));  
      end
      

    end
    fprintf(' .. done. \n');
    
    
    
    
    %% print images
    %  -----------------------------------------------------------------
    if opt.reprint 
      fprintf('Print:           ');
      for fi=1:opt.subStepSize:numel(opt.method(1).T{di})
        %try
          vbm_tst_cmpnii(opt.RAW.T{di}{fi});
        %end
        fprintf(sprintf('%s%%4d/%%4d',sprintf('%s',repmat('\b',1,9))),fi,numel(opt.RAW.T{di})); 
      end   
      fprintf(' .. done. \n');   
    end
    
    
    
    
    %% Estimate image QA parameter:
    %  -----------------------------------------------------------------
    % There are two phantom typs - phantoms with only one ground truth
    % for all datasets (BWP-Collins) and phantoms with a ground truth
    % dataset for each image (BWP-20N, SVE, ..).
    %   ds('l2','',[1 1 1],T/signal,p0T==3,T/signal,p0T/3,90)
%    fprintf('RAW-QA:          ');
%    for fi=1:opt.subStepSize:numel(opt.RAW.T{di})
%       if numel(opt.RAW.p0T{di})==1 && exist(opt.RAW.p0T{di}{1},'file')
%         [opt.RAW.qa(di,fi).p0T,opt.RAW.qam(di,fi).p0T] = vbm_tst_t1qa(...
%           opt.RAW.T{di}{fi},opt.RAW.p0T{di}{1},'',struct('verb',0));
%       elseif numel(opt.RAW.p0T{di})>1 && exist(opt.RAW.p0T{di}{fi},'file')
%         [opt.RAW.qa(di,fi).p0T,opt.RAW.qam(di,fi).p0T]  = vbm_tst_t1qa(...
%           opt.RAW.T{di}{fi},opt.RAW.p0T{di}{fi},'',struct('verb',0));   
%       elseif exist(opt.RAW.psT{di}{fi},'file')
%         [opt.RAW.qa(di,fi).psT,opt.RAW.qam(di,fi).psT] = vbm_tst_t1qa(...
%           opt.RAW.T{di}{fi},opt.RAW.psT{di}{fi},'',struct('verb',0));
%       end
%         
%         fprintf(sprintf('%s%%4d/%%4d',sprintf('%s',repmat('\b',1,9))),fi,numel(opt.RAW.T{di})); 
%       else
%       
%         S = vbm_io_xml(opt.RAW.XML{di}{fi});
%         if (numel(opt.RAW.p0T{di})==1 && exist(opt.RAW.p0T{di}{1},'file')) || ...
%            (numel(opt.RAW.p0T{di})>1 && exist(opt.RAW.p0T{di}{fi},'file'))
%           opt.RAW.qa(di,fi).p0T  = S.qa;
%           opt.RAW.qam(di,fi).p0T = S.qam;
%         else
%           opt.RAW.qa(di,fi).psT  = S.qa;
%           opt.RAW.qam(di,fi).psT = S.qam;
%         end
%       end
%    end
%    fprintf(' .. done. \n');
    clear T p0T;
   
    
    
    
    
    
    %% Evaluation 
    %{
    %  -----------------------------------------------------------------
    for mi=1:size(opt.methods,1)
      fprintf('%s-EV:          ',opt.methods{mi});
      for fi=1:opt.subStepSize:numel(opt.method(mi).T{di})   
        if def.methods{mi,6}>0 || ~exist(opt.method(mi).XML{di}{fi},'file') 
          % evaluation ...
          for ii=1:numel(opt.methods{mi,3})
            % real ground truth data
            if numel(opt.RAW.pmT{di})==1 && ...
              exist(opt.method(mi).([opt.methods{mi,3}{ii} 'T']){di}{fi},'file') && ...
              exist(opt.RAW.pmT{di}{1},'file') 
              [txt,tab,val] = vbm_tst_calc_kappa( ...
                opt.method(mi).([opt.methods{mi,3}{ii} 'T']){di}{fi}, ...
                opt.RAW.pmT{di}{1},opt.methods{mi,1},0);
            elseif numel(opt.RAW.p0T{di})==1 && ...
              exist(opt.method(mi).([opt.methods{mi,3}{ii} 'T']){di}{fi},'file') && ...
              exist(opt.RAW.p0T{di}{1},'file') 
              [txt,tab,val] = vbm_tst_calc_kappa( ...
                opt.method(mi).([opt.methods{mi,3}{ii} 'T']){di}{fi}, ...
                opt.RAW.p0T{di}{1},opt.methods{mi,1},0);
            
            elseif numel(opt.RAW.p0T{di})>1 && ...
              exist(opt.method(mi).([opt.methods{mi,3}{ii} 'T']){di}{fi},'file') && ...
              exist(opt.RAW.p0T{di}{fi},'file') 
            
              [txt,tab,val] = vbm_tst_calc_kappa( ...
                opt.method(mi).([opt.methods{mi,3}{ii} 'T']){di}{fi}, ...
                opt.RAW.p0T{di}{fi},opt.methods{mi,1},0);
            end

            % stable ground truth
            if exist(opt.method(mi).([opt.methods{mi,3}{ii} 'T']){di}{fi},'file') && ...
              exist(opt.RAW.psT{di}{fi},'file') 
              
              [txts,tabs,vals] = vbm_tst_calc_kappa( ...
                opt.method(mi).([opt.methods{mi,3}{ii} 'T']){di}{fi}, ...
                opt.RAW.psT{di}{fi},opt.methods{mi,1},0);
             
              %???????????
              % Biaskorrektur-test
              opt.RAW.qa(di,fi).psT = 0; %vbm_tst_t1qa( ...
%                 opt.method(mi).T{di}{fi}, ...
%                 opt.method(mi).([opt.methods{mi,3}{ii} 'T']){di}{fi}, ...
%                 opt.method(mi).([opt.methods{mi,4}{ii} 'T']){di}{fi}, struct('verb',0));
              
              vbm_io_xml(opt.RAW.XML{di}{fi},struct('qa',opt.RAW.qa(di,fi)),'write+');
            end
              
            try    res.method(mi).qa(di,fi).([opt.methods{mi,3}{ii} 'Tp0GT']) = val; 
            catch, res.method(mi).qa(di,fi).([opt.methods{mi,3}{ii} 'Tp0GT']) = struct();  %#ok<CTCH>
            end
            
            try    res.method(mi).qa(di,fi).([opt.methods{mi,3}{ii} 'TpsGT']) = vals; 
            catch, res.method(mi).qa(di,fi).([opt.methods{mi,3}{ii} 'TpsGT']) = struct();  %#ok<CTCH>
            end
           
            % save results
            vbm_io_xml(opt.method(mi).XML{di}{fi},struct('qa',res.method(mi).qa(di,fi)),'write+');
            
          end
          
          fprintf(sprintf('%s%%4d/%%4d',sprintf('%s',repmat('\b',1,9))),fi,numel(opt.RAW.T{di}));
        else
          % load times
          S = vbm_io_xml(opt.method(mi).XML{di}{fi}); 
          if isfield(S,'qa')
            fn=fieldnames(S.qa);
            for fni=1:numel(fn)
              res.method(mi).qa(di,fi).(fn{fni}) = S.qa.(fn{fni});
            end
          end
        end
        
      end
      fprintf(' .. done. \n');
    end
    %}
  end
  
  
  
  % write results 
  % --------------------------------------------------------------------
  %{

  
  result = cell(1,numel(opt.subdirs)); resultNF = result;
  resultqa = result; QAT=cell(numel(opt.subdirs),1); QATT=QAT;
  for di=1:numel(opt.subdirs)
    qan=5;
    % result{di}  = nan(numel(1:opt.subStepSize:numel(opt.method(mi).([opt.methods{mi,3}{ii} 'T']){di})),1);
    resultqa{di} = nan(numel(1:opt.subStepSize:numel(opt.method(mi).([opt.methods{mi,3}{ii} 'T']){di})),qan);
    for fi=1:opt.subStepSize:numel(opt.method(mi).([opt.methods{mi,3}{ii} 'T']){di})
   
      % Generation of the strings for the methodnames and its substep variables
      % ----------------------------------------------------------------
      ix=1;
      for mi=1:size(opt.methods,1)
        for ii=1:numel(opt.methods{mi,3})
          try
            %result(ix) = mean(res.method(mi).([opt.methods{mi,3}{ii} 'T']){di}{fi}.SEG.kappa(2)); %#ok<AGROW>
            result{di}(fi,ix) = mean(res.method(mi).qa(di,fi).([opt.methods{mi,3}{ii} 'Tp0GT']).SEG.kappa(2));
          catch %#ok<CTCH>
            try
              result{di}(fi,ix) = mean(res.method(mi).qa(di,fi).([opt.methods{mi,3}{ii} 'TpsGT']).SEG.kappa(2));
            catch %#ok<CTCH>
              result{di}(fi,ix) = nan;
            end
          end
          
          methodname{ix}=opt.methods{mi,1}; %#ok<AGROW>
          methodstep{ix}=opt.methods{mi,3}{ii}; %#ok<AGROW>
          ix=ix+1;
        end
      end

   
      
      % Generation of the header-lines for the command line output
      % ----------------------------------------------------------------
      if fi==1
        methodnamestr = ''; methodstepstr = ''; methodresstr = ''; methodresstrd='';
        for mi=1:numel(methodname)
          methodnamestr = [methodnamestr ...
            sprintf(sprintf('%%%ds ',opt.printspace(2)),methodname{mi})]; %#ok<AGROW>
          methodstepstr = [methodstepstr ...
            sprintf(sprintf('%%%ds ',opt.printspace(2)),methodstep{mi})]; %#ok<AGROW>
          methodresstr = [methodresstr ...
            sprintf('%%%d.%df ',opt.printspace(2),opt.printspace(3))]; %#ok<AGROW>
          methodresstrd = [methodresstrd ...
            sprintf('%%%dd ',opt.printspace(2))]; %#ok<AGROW>
        end
        Phdr.method = sprintf(...
          sprintf('%%%ds: %%%ds %%%ds %%%ds %%%ds %%%ds | %%s\n', ...
            opt.printspace(1),repmat(opt.printspace(2),1,qan)), ...
          sprintf('scans of subdir %s',opt.subdirs{di}),...
            'noise','bias','vol','isotr','BWC',methodnamestr);
        Phdr.version =  sprintf(...
          sprintf('%%%ds  %%%ds %%%ds %%%ds %%%ds %%%ds | %%s\n', ...
            opt.printspace(1),repmat(opt.printspace(2),1,qan)), ...
          '','','','','','',methodstepstr);
        Phdr.lines = sprintf(...
          sprintf('%%%ds: %%%ds %%%ds %%%ds %%%ds %%%ds | %%s\n', ...
            opt.printspace(1),repmat(opt.printspace(2),1,qan)), ...
          '','','','','','',methodstepstr);

        QAT{di} = sprintf('\n%s%s%s\n',Phdr.method,Phdr.version,repmat('-',size(Phdr.method)));
        QATT{di} = [{['scans of subdir ' opt.subdirs{di}],'noise','bias','vol','isotr','BWC'},methodname; ...
                    {''                  ,''     ,''    ,''   ,''     ,''   },methodstep];
        
        if di==1
          QATF    = sprintf('\n%s%s%s\n',Phdr.method,Phdr.version,repmat('-',size(Phdr.method)));
        end
      end
      
      
      
      
      
      % qa results of the ground truth or stable ground truth
      % ----------------------------------------------------------------
      try 
        resultqa{di}(fi,:) = [opt.RAW.qa(di,fi).p0T.noise,opt.RAW.qa(di,fi).p0T.bias_WMstd,...
                      opt.RAW.qa(di,fi).p0T.res_vol, opt.RAW.qa(di,fi).p0T.res_isotropy, ...
                      opt.RAW.qa(di,fi).p0T.contrast]; 
      catch %#ok<CTCH>
        try
         resultqa{di}(fi,:) = [opt.RAW.qa(di,fi).psT.noise,opt.RAW.qa(di,fi).psT.bias_WMstd,...
                      opt.RAW.qa(di,fi).psT.res_vol, opt.RAW.qa(di,fi).psT.res_isotropy, ...
                      opt.RAW.qa(di,fi).psT.contrast]; 
        end
      end

      
      
     % print brain extraction results (p0*)
     % print segmentation results (p0)
     % print bias correction (m0)
      [pp,ff]=fileparts(opt.method(1).T{di}{fi}); [pp,hh]=fileparts(pp); fn=[hh filesep ff]; clear pp ff;

      QATT{di} = [QATT{di}; fn num2cell([resultqa{di}(fi,:),result{di}(fi,:)])];
      QAline = sprintf(sprintf('%%%ds: %%%d.%df %%%d.%df %%%d.%df %%%d.%df %%%d.%df | %s\n', ...
            opt.printspace(1),[repmat(opt.printspace(2),1,qan);repmat(opt.printspace(3),1,qan)],methodresstr), ...
            spm_str_manip(fn,'a50'),resultqa{di}(fi,:),result{di}(fi,:));
      QAT{di} = [QAT{di} QAline]; 
      
      if any(result{di}(fi,:)<cell2mat([opt.methods{:,5}]))
        QATF = [QATF, QAline]; %#ok<AGROW>
      end
    end
    resultNF{di} = result{di}; resultNF{di}(result{di}<opt.FailedTh) = nan;
    
    
    QATT{di} = [QATT{di}; 
     sprintf('mean (>%0.2f)',opt.FailedTh),num2cell([vbm_stat_nanmean(resultqa{di}), vbm_stat_nanmean(resultNF{di})]); 
     sprintf('std  (>%0.2f)',opt.FailedTh),num2cell([vbm_stat_nanstd(resultqa{di}),  vbm_stat_nanstd(resultNF{di})]); 
     sprintf('min  (>%0.2f)',opt.FailedTh),num2cell([min(resultqa{di}),  min(resultNF{di})]);
     {sprintf('failed scans (of %d)',size(cell2mat(result'),1)),'','','','','',} ...
      num2cell( sum(result{di}<opt.FailedTh));
     {'failed %'    ,'','','','','',}  num2cell(100* (sum(cell2mat(result')<opt.FailedTh) / size(result{di},1))) ];
    
    QAT{di} = [QAT{di}, ...
      sprintf('%s\n',repmat('-',size(Phdr.method))), ...
      sprintf(sprintf('%%%ds: %%%d.%df %%%d.%df %%%d.%df %%%d.%df %%%d.%df | %s\n', ...
          opt.printspace(1),[repmat(opt.printspace(2),1,qan);repmat(opt.printspace(3),1,qan)],methodresstr), ...
          'mean',vbm_stat_nanmean(resultqa{di}), vbm_stat_nanmean(result{di})), ... 
      sprintf(sprintf('%%%ds: %%%d.%df %%%d.%df %%%d.%df %%%d.%df %%%d.%df | %s\n', ...
          opt.printspace(1),[repmat(opt.printspace(2),1,qan);repmat(opt.printspace(3),1,qan)],methodresstr), ...
          'std',vbm_stat_nanstd(resultqa{di}), vbm_stat_nanstd(result{di})), ... 
      ]; 
    
    vbm_io_csv(fullfile(opt.path,def.CSVdir,['QAT_' opt.subdirs{di} '.csv']),QATT{di},'','',opt.csv);
    
    
    
    
    
    %fprintf(QAT{di});
  end
  QATF = [QATF, ...
      sprintf('%s\n',repmat('-',size(Phdr.method))), ...
      sprintf(sprintf('%%%ds: %%%d.%df %%%d.%df %%%d.%df %%%d.%df %%%d.%df | %s\n', ...
          opt.printspace(1),[repmat(opt.printspace(2),1,qan);repmat(opt.printspace(3),1,qan)],methodresstr), ...
          sprintf('mean (>%0.2f)',opt.FailedTh),vbm_stat_nanmean(resultqa{di}), vbm_stat_nanmean(cell2mat(resultNF'))), ...
      sprintf(sprintf('%%%ds: %%%d.%df %%%d.%df %%%d.%df %%%d.%df %%%d.%df | %s\n', ...
          opt.printspace(1),[repmat(opt.printspace(2),1,qan);repmat(opt.printspace(3),1,qan)],methodresstr), ...
          sprintf('std  (>%0.2f)',opt.FailedTh),vbm_stat_nanstd(resultqa{di}), vbm_stat_nanstd(cell2mat(resultNF'))), ...
      sprintf('%s\n',repmat('-',size(Phdr.method))),...
      sprintf(sprintf('%%%ds: %%%ds %%%ds %%%ds %%%ds %%%ds | %s\n', ...
          opt.printspace(1),repmat(opt.printspace(2),1,qan),methodresstrd), ...
          sprintf('failed (of %d)',size(cell2mat(result'),1)),'','','','','', sum(cell2mat(result')<opt.FailedTh)), ... 
      sprintf(sprintf('%%%ds: %%%ds %%%ds %%%ds %%%ds %%%ds | %s\n', ...
          opt.printspace(1),repmat(opt.printspace(2),1,qan),methodresstr), ...
          'failed percent','','','','','', (100 * sum(cell2mat(result')<opt.FailedTh) ./ size(cell2mat(result'),1))), ... 
      sprintf('%s\n\n',repmat('-',size(Phdr.method))), ...
      ]; 
  
   fprintf(1,QATF);
   QATFtxt = fopen(fullfile(opt.path,def.CSVdir,['QATF_' datestr(clock,'YYYYmmdd-HHMMSS') '.txt']),'w');
   fprintf(QATFtxt,QATF);
   fclose(QATFtxt); clear QATFtxt
%%
 %}
  
  % 4) statistic
  
  % 5) result figures
  
end

function linkfile(source,destination,f)
  if isunix || ismac
    if exist('f','var') && f=='f'
      system(sprintf('ln -s -f %s %s',source,destination));
    else
      system(sprintf('ln -s %s %s',source,destination));
    end
  else
    if exist('f','var') && f=='f'
      copyfile(source,destination,f);
    else
      copyfile(source,destination);
    end
  end
end
function N3(file,initPATH) 

  wkd=pwd;
  [pp,ff,ee]=spm_fileparts(file); cd(pp);
  try 

    calcN3=[ ...
      ... filenames
      'pp=$(dirname $T); ff=$(basename $T .nii); ' ...
      'oT=$pp/$ff; mT=$pp/m$ff;' ...
      ... create mask and convert to mnc
...      'mri_convert ${oT}.nii  ${oT}.mnc; ' ...
      ... N3-nu-correction
...      'nu_correct -quiet -clobber -stop 0.01 -distance 40 -iterations 100 -fwhm 0.05 ${oT}.mnc ${mT}.mnc; cp -f ${mT}.mnc ${oT}.mnc ;' ...
...      'nu_correct -quiet -clobber -stop 0.01 -distance 20 -iterations 100 -fwhm 0.05 ${oT}.mnc ${mT}.mnc; cp -f ${mT}.mnc ${oT}.mnc;' ...
...      'nu_correct -quiet -clobber -stop 0.01 -distance 10 -iterations 100 -fwhm 0.05 ${oT}.mnc ${mT}.mnc; cp -f ${mT}.mnc ${oT}.mnc;' ...
...      ... convert back to nifti and delete files
...      'mri_convert ${mT}.mnc ${mT}.nii;' ...
...      'rm ${oT}.mnc ${mT}.mnc ${mT}.imp; ' ...
     'cp -f ${oT}.nii ${mT}.nii;'
      ];
   
    [SS,SR]=system(sprintf('%s T="%s"; %s',initPATH,file,calcN3)); 

  catch  %#ok<CTCH>
    createNullImage(fullfile(pp,[ff ee]),fullfile(pp,['m' ff ee]));
  end
      
  cd(wkd);

end
function FSL(file,initPATH,p0T) 
  
  if ~exist('p0T','var') || ~isempty(p0T);
    FSLtype = 'original'; p0T = '';
  elseif ~exist('p0T','var') || ~exist(p0T,'file')   
    return
  else
    FSLtype = 'optimal';
  end
  
  [pp,ff,ee] = spm_fileparts(file);
  Pp0 = fullfile(pp,sprintf('p0%s%s',ff,ee));
  Pn  = fullfile(pp,sprintf('n%s%s',ff,ee));
  Pm  = fullfile(pp,sprintf('n%s%s',ff,ee));
  
  if strfind(file,'Apes')
     
    % set atlas maps
    if ~isempty([strfind(lower(ff),'chimpanzee') strfind(lower(ff),'orangutan') strfind(lower(ff),'gorilla')])
      atlas = '/Users/dahnke/Neuroimaging/spm12/toolbox/vbm12/templates_1.50mm/TPM_gapes0_T1.nii';
    elseif ~isempty([strfind(lower(ff),'rhesus') strfind(lower(ff),'capuchin') strfind(lower(ff),'gibbon')])
      atlas = '/Users/dahnke/Neuroimaging/spm12/toolbox/vbm12/templates_1.50mm/TPM_lapes0_T1.nii'; %return
    elseif ~isempty([strfind(lower(ff),'squirrel')])
      atlas = '/Users/dahnke/Neuroimaging/spm12/toolbox/vbm12/templates_1.50mm/TPM_monkeys0_T1.nii'; return
    else
      atlas = '/Users/dahnke/Neuroimaging/spm12/toolbox/vbm12/templates_1.50mm/TPM_gapes1_T1.nii'; % macaqure, baboon, ...
    end
    
    % call N3 inhomogeneity correction
    if 1 %~exist(Pp0,'file') || ~exist(Pn,'file') 
      FSLtype = 'apes'; 
      N3(file,initPATH);
      copyfile(fullfile(pp,['m' ff ee]),fullfile(pp,['n' ff ee]));
    end
    
    betape = '';
  else
    betape = '';
  end
  
  wkd=pwd;
  [pp,ff,ee]=spm_fileparts(file); cd(pp);
  
  if      strfind(lower(ff),'t2'), modality='2';
  elseif  strfind(lower(ff),'pd'), modality='3';
  else                             modality='1';  
  end
  

  try
    calcFSL.init  = [ ...
      'ff=$(basename $T .nii); pp=$(dirname $T);' ...
      'oT=$pp/$ff; nT=$pp/n$ff; paT=$pp/pa$ff; psT=$pp/ps$ff; p0T=$pp/p0$ff; mT=$pp/m$ff; ebT=$pp/eb$ff;' ...
      'omat=$pp/${ff}_flirt.mat;' ...
      ];
    calcFSL.susan     = 'susan ${oT}.nii -1 1 3 1 0 ${nT}.nii;'; % noise correction 
    calcFSL.BET_FAST  = [ ...
      'BET2 ${nT}.nii ${paT} ' betape ';' ...% skull-stripping -S -R
      'fast -B -o ${oT} -t ' modality ' ${paT}.nii;']; % -B: output Ym;  -b: output: Yeb; 
    calcFSL.STABLE_FAST = [...
      'BET2 ${nT}.nii ${paT};' ...
      'cp %{oT}.nii ${psT}.nii;' ...
      'fslmath ${psT}.nii -mas ' p0T ';' ...
      'fast -b -B -t ' modality '-o ${oT} ${paT}.nii;'];
    calcFSL.cleanup = [ ...
      'mv -f ${oT}_restore.nii ${mT}.nii;' ...
      ... 'mv -f ${oT}_bias.nii ${ebT}.nii;' ...
      'rm -f ${oT}_pveseg.nii ${oT}_seg.nii ${oT}_mixeltype.nii;' ...
      ...'rm -f ${nT}.nii ${paT}.nii' ...
      'find ${pp} -name "*.mnc" -delete' ...
      ];
    calcFSL.flirt = [ ... linear registration
      'flirt -dof 12              -ref ${p0T}  -in ${Aa}.nii    -omat mat0;' ... init
      'convert_xfm -omat mat0 -inverse mat0;' ...
       'flirt -init mat0 -applyxfm -ref ${A} -in ${p0T}.nii -out ${p0T}.a1.nii;  ' ...
      'flirt -dof 12  -ref ${A}  -in ${p0T}.a1.nii     -omat ${omat};' ... fine
      'flirt -init ${omat} -applyxfm -ref ${A} -in ${oT}.nii -out ${mT}.affine.nii;  ' ... % ... auf n ändern
      ...'flirt -interp sinc -init ${omat} -applyxfm -ref ${A} -in ${oT}_pve_0.nii -out ${oT}_pve_0.affine.nii; ' ...
      ...'flirt -interp sinc -init ${omat} -applyxfm -ref ${A} -in ${oT}_pve_1.nii -out ${oT}_pve_1.affine.nii; ' ...
      ...'flirt -interp sinc -init ${omat} -applyxfm -ref ${A} -in ${oT}_pve_2.nii -out ${oT}_pve_2.affine.nii; ' ...
    ];
  %%
    switch FSLtype
      case 'original'
        % susan filtering, bet for skull-stripping
        calcFSLs = [calcFSL.init calcFSL.susan calcFSL.BET_FAST calcFSL.cleanup];
      case 'optimal'
        % sanlm-filtered image and ground truth segmentation for fast
        calcFSLs = [calcFSL.init calcFSL.susan calcFSL.STABLE_FAST calcFSL.cleanup];
      case 'apes'
        % sanlm-filtered image and ground truth segmentation for fast
        calcFSLs = [calcFSL.init calcFSL.BET_FAST calcFSL.cleanup];
    end
    
    
    if ~exist(Pp0,'file') || ~exist(Pm,'file') 
      [SS,SR] = system(sprintf('%s T="%s"; %s',initPATH,file,calcFSLs)); 


      switch modality
        case '1'
          vbm_io_cgw2seg(fullfile(pp,[ff '_pve_0.nii']), ...
                         fullfile(pp,[ff '_pve_1.nii']), ...
                         fullfile(pp,[ff '_pve_2.nii']),'FSL',1); %~strcmp(FSLtype,'apes'));
        case '2'                 
         vbm_io_cgw2seg(fullfile(pp,[ff '_pve_0.nii']), ...
                        fullfile(pp,[ff '_pve_1.nii']), ...
                        fullfile(pp,[ff '_pve_2.nii']),'FSL',1); %~strcmp(FSLtype,'apes'));
        case '3'                 
         vbm_io_cgw2seg(fullfile(pp,[ff '_pve_0.nii']), ...
                        fullfile(pp,[ff '_pve_2.nii']), ...
                        fullfile(pp,[ff '_pve_1.nii']),'FSL',1); %~strcmp(FSLtype,'apes'));
      end
    end
    
    % flirt realignment ... does not work very well...
    if 0 %strcmp(FSLtype,'apes')
      %%
      [SS,SR] = system(sprintf('%s T="%s"; A="%s"; As="%s"; Ass="%s"; Aa="%s"; %s',...
        initPATH,file,atlas,atlass,atlasss,atlasa,[calcFSL.init calcFSL.flirt])); 
    
      if 0 
        switch modality
          case '1'
           vbm_io_cgw2seg(fullfile(pp,[ff '_pve_0.affine.nii']), ...
                           fullfile(pp,[ff '_pve_1.affine.nii']), ...
                           fullfile(pp,[ff '_pve_2.affine.nii']),'FSL',1);
          case '2'                 
            vbm_io_cgw2seg(fullfile(pp,[ff '_pve_0.affine.nii']), ...
                           fullfile(pp,[ff '_pve_1.affine.nii']), ...
                           fullfile(pp,[ff '_pve_2.affine.nii']),'FSL',1);
          case '3'                 
            vbm_io_cgw2seg(fullfile(pp,[ff '_pve_0.affine.nii']), ...
                           fullfile(pp,[ff '_pve_2.affine.nii']), ...
                           fullfile(pp,[ff '_pve_1.affine.nii']),'FSL',1);
        end  
        movefile(fullfile(pp,['p0' ff '_pve_0..nii']),fullfile(pp,['p0' ff '.affine.nii']));   
      end
    end
  catch  %#ok<CTCH>
    %createNullImage(fullfile(pp,[ff ee]),fullfile(pp,['p0' ff ee]));
    %createNullImage(fullfile(pp,[ff ee]),fullfile(pp,['m' ff ee]));
  end
                 
  cd(wkd);

end
function SPM8segment(file,SPM8dir,SPMwkd)
%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev$)
%-----------------------------------------------------------------------
%   matlabbatch{1}.spm.spatial.coreg.estimate.ref     = cellstr(fullfile(spm('Dir'),'templates','T1.nii'));
%   matlabbatch{1}.spm.spatial.coreg.estimate.source  = cellstr(file);
%   matlabbatch{1}.spm.spatial.coreg.estimate.other   = {''};
%   matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
%   matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep  = [4 2];
%   matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol  = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
%   matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
%  

  matlabbatch{1}.spm.spatial.preproc.data           = cellstr(file);
  matlabbatch{1}.spm.spatial.preproc.output.GM      = [0 0 1];
  matlabbatch{1}.spm.spatial.preproc.output.WM      = [0 0 1];
  matlabbatch{1}.spm.spatial.preproc.output.CSF     = [0 0 1];
  matlabbatch{1}.spm.spatial.preproc.output.biascor = 1;
  matlabbatch{1}.spm.spatial.preproc.output.cleanup = 1;
  matlabbatch{1}.spm.spatial.preproc.opts.tpm = cellstr(char(...
    fullfile(SPM8dir,'tpm','grey.nii'),...
    fullfile(SPM8dir,'tpm','white.nii'),...
    fullfile(SPM8dir,'tpm','csf.nii')));   % Prior probability maps
  matlabbatch{1}.spm.spatial.preproc.opts.ngaus     = [2 2 2 4];
  matlabbatch{1}.spm.spatial.preproc.opts.regtype   = 'mni';
  matlabbatch{1}.spm.spatial.preproc.opts.warpreg   = 1;
  matlabbatch{1}.spm.spatial.preproc.opts.warpco    = 25;
  matlabbatch{1}.spm.spatial.preproc.opts.biasreg   = 0.0001;
  matlabbatch{1}.spm.spatial.preproc.opts.biasfwhm  = 60;
  matlabbatch{1}.spm.spatial.preproc.opts.samp      = 3;
  matlabbatch{1}.spm.spatial.preproc.opts.msk       = {''};

  pathchanges = replaceSPMpath(SPMwkd,SPM8dir);
  
  [pp,ff,ee]=fileparts(file);
  warning off; %#ok<*WNOFF>
  try
    clear global classes;
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
                 
    delete(fullfile(pp,[ff '_seg_sn.mat']));
    delete(fullfile(pp,[ff '_seg_inv_sn.mat']));
  catch e %#ok<NASGU>
    createNullImage(fullfile(pp,[ff ee]),fullfile(pp,['p0' ff ee]));
  end
  warning on; %#ok<*WNON>
  
  restoreSPMpath(pathchanges) 
  
  vbm_io_cgw2seg(fullfile(pp,['c3' ff '.nii']), ...
                 fullfile(pp,['c1' ff '.nii']), ...
                 fullfile(pp,['c2' ff '.nii']),'SPM');
                  
end
function SPM8newsegment(file,SPM8dir,SPMwkd)
%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev$)
%-----------------------------------------------------------------------
  matlabbatch{1}.spm.tools.preproc8.channel.vols      = cellstr(file);
  matlabbatch{1}.spm.tools.preproc8.channel.biasreg   = 0.0001;
  matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm  = 60;
  matlabbatch{1}.spm.tools.preproc8.channel.write     = [0 1];
  matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm     = {fullfile(spm('dir'),'TPM','TPM.nii,1')};
  matlabbatch{1}.spm.tools.preproc8.tissue(1).ngaus   = 2;
  matlabbatch{1}.spm.tools.preproc8.tissue(1).native  = [1 0];
  matlabbatch{1}.spm.tools.preproc8.tissue(1).warped  = [0 0];
  matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm     = {fullfile(spm('dir'),'TPM','TPM.nii,2')};
  matlabbatch{1}.spm.tools.preproc8.tissue(2).ngaus   = 2;
  matlabbatch{1}.spm.tools.preproc8.tissue(2).native  = [1 0];
  matlabbatch{1}.spm.tools.preproc8.tissue(2).warped  = [0 0];
  matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm     = {fullfile(spm('dir'),'TPM','TPM.nii,3')};
  matlabbatch{1}.spm.tools.preproc8.tissue(3).ngaus   = 2;
  matlabbatch{1}.spm.tools.preproc8.tissue(3).native  = [1 0];
  matlabbatch{1}.spm.tools.preproc8.tissue(3).warped  = [0 0];
  matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm     = {fullfile(spm('dir'),'TPM','TPM.nii,4')};
  matlabbatch{1}.spm.tools.preproc8.tissue(4).ngaus   = 3;
  matlabbatch{1}.spm.tools.preproc8.tissue(4).native  = [0 0];
  matlabbatch{1}.spm.tools.preproc8.tissue(4).warped  = [0 0];
  matlabbatch{1}.spm.tools.preproc8.tissue(5).tpm     = {fullfile(spm('dir'),'TPM','TPM.nii,5')};
  matlabbatch{1}.spm.tools.preproc8.tissue(5).ngaus   = 4;
  matlabbatch{1}.spm.tools.preproc8.tissue(5).native  = [0 0];
  matlabbatch{1}.spm.tools.preproc8.tissue(5).warped  = [0 0];
  matlabbatch{1}.spm.tools.preproc8.tissue(6).tpm     = {fullfile(spm('dir'),'TPM','TPM.nii,6')};
  matlabbatch{1}.spm.tools.preproc8.tissue(6).ngaus   = 2;
  matlabbatch{1}.spm.tools.preproc8.tissue(6).native  = [0 0];
  matlabbatch{1}.spm.tools.preproc8.tissue(6).warped  = [0 0];
  matlabbatch{1}.spm.tools.preproc8.warp.reg          = 4;
  matlabbatch{1}.spm.tools.preproc8.warp.affreg       = 'mni';
  matlabbatch{1}.spm.tools.preproc8.warp.samp         = 3;
  matlabbatch{1}.spm.tools.preproc8.warp.write        = [0 0];

  pathchanges = replaceSPMpath(SPMwkd,SPM8dir);
  
  [pp,ff,ee]=spm_fileparts(file);
  try
    clear global;
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);
                 
    if matlabbatch{mb}.spm.spatial.preproc.channel.write(2)            
      move(fullfile(pp,['BiasField_' ff '.nii']),fullfile(pp,['eb' ff '.nii']));
    end
             
  catch e  %#ok<NASGU>
    createNullImage(fullfile(pp,[ff ee]),fullfile(pp,['p0' ff ee]));
  end

  restoreSPMpath(pathchanges)
  
  vbm_io_cgw2seg(fullfile(pp,['c3' ff '.nii']), ...
                 fullfile(pp,['c1' ff '.nii']), ...
                 fullfile(pp,['c2' ff '.nii']),'SPM');
  
  if exist(fullfile(pp,[ff '_seg8.mat']),'file')
    delete(fullfile(pp,[ff '_seg8.mat']));
  end
end
function SPM12segment(file,SPM12dir,SPMwkd)
%   SPM12dirs    = [SPM12dir vbm_findfiles(SPM12dir,'*',struct('dirs','1')) ...
%                   fullfile(SPM12path,'toolbox','VBM12+')]';
%   SPM12baddirs = ~cellfun('isempty',strfind(SPM12dirs,'.svn')) | ... 
%                 ~cellfun('isempty',strfind(SPM12dirs,'private'))| ... 
%                 ~cellfun('isempty',strfind(SPM12dirs,'@'));
%   SPM12dirs(SPM12baddirs) = [];
%   
%   oldpath    = textscan(path,'%s','bufsize',2^16,'Delimiter',':'); oldpath=oldpath{1};
%   vbm12pathi = ~cellfun('isempty',strfind(oldpath,SPM12path));
%   
%   olddir = pwd; cd(SPM12dir);
%   [pp,ff,ee]=spm_fileparts(file);
% 
%   for pi=1:numel(vbm12pathi), if vbm12pathi(pi), rmpath(oldpath{pi}); end; end
%   for pi=numel(SPM12dirs):-1:1, addpath(SPM12dirs{pi},'-begin'); end

%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev$)
%-----------------------------------------------------------------------
  mb=1;
  matlabbatch{mb}.spm.spatial.preproc.channel.vols      = cellstr(file);
  matlabbatch{mb}.spm.spatial.preproc.channel.biasreg   = 0.0001;
  matlabbatch{mb}.spm.spatial.preproc.channel.biasfwhm  = 60;
  matlabbatch{mb}.spm.spatial.preproc.channel.write     = [0 1];
  matlabbatch{mb}.spm.spatial.preproc.tissue(1).tpm     = {'/Users/dahnke/Neuroimaging/spm12/tpm/TPM.nii,1'};
  matlabbatch{mb}.spm.spatial.preproc.tissue(1).ngaus   = 1;
  matlabbatch{mb}.spm.spatial.preproc.tissue(1).native  = [1 0];
  matlabbatch{mb}.spm.spatial.preproc.tissue(1).warped  = [0 0];
  matlabbatch{mb}.spm.spatial.preproc.tissue(2).tpm     = {'/Users/dahnke/Neuroimaging/spm12/tpm/TPM.nii,2'};
  matlabbatch{mb}.spm.spatial.preproc.tissue(2).ngaus   = 1;
  matlabbatch{mb}.spm.spatial.preproc.tissue(2).native  = [1 0];
  matlabbatch{mb}.spm.spatial.preproc.tissue(2).warped  = [0 0];
  matlabbatch{mb}.spm.spatial.preproc.tissue(3).tpm     = {'/Users/dahnke/Neuroimaging/spm12/tpm/TPM.nii,3'};
  matlabbatch{mb}.spm.spatial.preproc.tissue(3).ngaus   = 2;
  matlabbatch{mb}.spm.spatial.preproc.tissue(3).native  = [1 0];
  matlabbatch{mb}.spm.spatial.preproc.tissue(3).warped  = [0 0];
  matlabbatch{mb}.spm.spatial.preproc.tissue(4).tpm     = {'/Users/dahnke/Neuroimaging/spm12/tpm/TPM.nii,4'};
  matlabbatch{mb}.spm.spatial.preproc.tissue(4).ngaus   = 3;
  matlabbatch{mb}.spm.spatial.preproc.tissue(4).native  = [0 0];
  matlabbatch{mb}.spm.spatial.preproc.tissue(4).warped  = [0 0];
  matlabbatch{mb}.spm.spatial.preproc.tissue(5).tpm     = {'/Users/dahnke/Neuroimaging/spm12/tpm/TPM.nii,5'};
  matlabbatch{mb}.spm.spatial.preproc.tissue(5).ngaus   = 4;
  matlabbatch{mb}.spm.spatial.preproc.tissue(5).native  = [0 0];
  matlabbatch{mb}.spm.spatial.preproc.tissue(5).warped  = [0 0];
  matlabbatch{mb}.spm.spatial.preproc.tissue(6).tpm     = {'/Users/dahnke/Neuroimaging/spm12/tpm/TPM.nii,6'};
  matlabbatch{mb}.spm.spatial.preproc.tissue(6).ngaus   = 2;
  matlabbatch{mb}.spm.spatial.preproc.tissue(6).native  = [0 0];
  matlabbatch{mb}.spm.spatial.preproc.tissue(6).warped  = [0 0];
  matlabbatch{mb}.spm.spatial.preproc.warp.mrf          = 1;
  matlabbatch{mb}.spm.spatial.preproc.warp.cleanup      = 1;
  matlabbatch{mb}.spm.spatial.preproc.warp.reg          = [0 0.001 0.5 0.025 0.1];
  matlabbatch{mb}.spm.spatial.preproc.warp.affreg       = 'mni';
  matlabbatch{mb}.spm.spatial.preproc.warp.fwhm         = 0;
  matlabbatch{mb}.spm.spatial.preproc.warp.samp         = 3;
  matlabbatch{mb}.spm.spatial.preproc.warp.write        = [0 0];
  
  pathchanges = replaceSPMpath(SPMwkd,SPM12dir);
  
  warning off;
  [pp,ff,ee]=spm_fileparts(file);
  try
    clear global;
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);  
  
    vbm_io_cgw2seg(fullfile(pp,['c3' ff '.nii']), ...
                   fullfile(pp,['c1' ff '.nii']), ...
                   fullfile(pp,['c2' ff '.nii']),'SPM');
     
    try
      biasfile = fullfile(pp,['BiasField_' ff '.nii']); 
      if matlabbatch{mb}.spm.spatial.preproc.channel.write(2) && exist(biasfile,'file');
        movefile(biasfile,fullfile(pp,['eb' ff '.nii']));
      end
    end
    
    delete(fullfile(pp,['p' ff '_seg8.txt']));
    delete(fullfile(pp,[ff '_seg8.mat']));
  catch e %#ok<NASGU>
    createNullImage(fullfile(pp,[ff ee]),fullfile(pp,['p0' ff ee]));
    createNullImage(fullfile(pp,[ff ee]),fullfile(pp,['m' ff ee]));
  end
  warning on;

  restoreSPMpath(pathchanges)  
  
  if exist(fullfile(pp,[ff '_seg8.mat']),'file')
    delete(fullfile(pp,[ff '_seg8.mat']));
  end
  

end
function VBM12segment(file,SPM12dir,SPMwkd,LAS)
% get old matlab-paths, seperate VBM12 entries and replace them by VBM8
% entries

%-----------------------------------------------------------------------
% Job saved on 26-Mar-2013 17:27:18 by cfg_util (rev $Rev$)
% spm SPM - SPM12b (5298)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
  matlabbatch{1}.spm.tools.vbm.estwrite.data              = cellstr(file);
  matlabbatch{1}.spm.tools.vbm.estwrite.opts.tpm          = ...
    {'/Users/dahnke/Neuroimaging/spm12/tpm/TPM.nii'};
    %{'/Users/dahnke/Neuroimaging/spm12/toolbox/vbm12/templates_1.50mm/TPM.nii'};
  matlabbatch{1}.spm.tools.vbm.estwrite.opts.ngaus        = [1 1 2 3 4 2];  % [2 2 2 3 4 2] 
  matlabbatch{1}.spm.tools.vbm.estwrite.opts.biasreg      = 0.001; % 0.0001 - stronger correction
  matlabbatch{1}.spm.tools.vbm.estwrite.opts.biasfwhm     = 60;
  matlabbatch{1}.spm.tools.vbm.estwrite.opts.affreg       = 'mni';
  matlabbatch{1}.spm.tools.vbm.estwrite.opts.warpreg      = [0 0.001 0.5 0.025 0.1];
  matlabbatch{1}.spm.tools.vbm.estwrite.opts.samp         = 3;
  
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.dartelwarp.normhigh.darteltpm = ...
    {'/Users/dahnke/Neuroimaging/spm12/toolbox/vbm12/templates_1.50mm/Template_1_IXI550_MNI152.nii'};
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.sanlm     = 2;
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.cleanup   = 0;
%   matlabbatch{1}.spm.tools.vbm.estwrite.extopts.gcutstr   = 0.5;
%   matlabbatch{1}.spm.tools.vbm.estwrite.extopts.gcutCSF   = 0;  
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.vox       = 1.5;
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.bb        = [[-90 -126 -72];[90 90 108]];
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.print     = 1;
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.LAS       = LAS;
%   matlabbatch{1}.spm.tools.vbm.estwrite.extopts.WMHC      = 0;  
%   matlabbatch{1}.spm.tools.vbm.estwrite.extopts.WMHCstr   = 0.5; 
%   matlabbatch{1}.spm.tools.vbm.estwrite.extopts.BVC       = 1;
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.ROI       = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.extopts.surface   = 0;
  
  %%
  matlabbatch{1}.spm.tools.vbm.estwrite.output.GM.native        = 1;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.GM.warped        = 1;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.GM.modulated     = 1;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.GM.dartel        = 1;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.WM.native        = 1;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.WM.warped        = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.WM.modulated     = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.WM.dartel        = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.CSF.native       = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.CSF.warped       = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.CSF.modulated    = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.CSF.dartel       = 0;
  
  matlabbatch{1}.spm.tools.vbm.estwrite.output.jacobian.warped  = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.warps            = [0 0];
  
  matlabbatch{1}.spm.tools.vbm.estwrite.output.bias.native      = 1;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.bias.warped      = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.bias.affine      = 0;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.label.native     = 1;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.label.warped     = 1;
  matlabbatch{1}.spm.tools.vbm.estwrite.output.label.dartel     = 1;

%   matlabbatch{1}.spm.tools.vbm.estwrite.output.atlas.native     = 1;
%   matlabbatch{1}.spm.tools.vbm.estwrite.output.atlas.warped     = 0;
%   matlabbatch{1}.spm.tools.vbm.estwrite.output.atlas.dartel     = 0;
  
%   matlabbatch{1}.spm.tools.vbm.estwrite.output.pc.native        = 1;
%   matlabbatch{1}.spm.tools.vbm.estwrite.output.pc.warped        = 0;
%   matlabbatch{1}.spm.tools.vbm.estwrite.output.pc.modulated     = 0;
%   matlabbatch{1}.spm.tools.vbm.estwrite.output.pc.dartel        = 0;
%   matlabbatch{1}.spm.tools.vbm.estwrite.output.te.native        = 1;
%   matlabbatch{1}.spm.tools.vbm.estwrite.output.te.warped        = 0;
%   matlabbatch{1}.spm.tools.vbm.estwrite.output.te.modulated     = 0;
%   matlabbatch{1}.spm.tools.vbm.estwrite.output.te.dartel        = 0;

  

  %%
  pathchanges = replaceSPMpath(SPMwkd,SPM12dir);
   
  warning off;
  try
    clear global;
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);  
%cg_vbm_batch(namefile,writeonly,vbm_defaults)
    %delete(fullfile(pp,[ff '_seg8.mat']));
  catch %#ok<CTCH>
    createNullImage(fullfile(pp,[ff ee]),fullfile(pp,['p0' ff ee]));
    createNullImage(fullfile(pp,[ff ee]),fullfile(pp,['mn' ff ee]));
    createNullImage(fullfile(pp,[ff ee]),fullfile(pp,['mg' ff ee]));
  end
  warning on;


  restoreSPMpath(pathchanges)
  
end

function VBM8segment(file,VBM8dir,SPMwkd)
% % get old matlab-paths, seperate VBM12 entries and replace them by VBM8
% % entries
%   VBM8dirs    = [VBM8dir vbm_findfiles(VBM8dir,'*',struct('dirs','1'))]';
%   VBM8baddirs = ~cellfun('isempty',strfind(VBM8dirs,'.svn')) | ... 
%                 ~cellfun('isempty',strfind(VBM8dirs,'private'))| ... 
%                 ~cellfun('isempty',strfind(VBM8dirs,'@'));
%   VBM8dirs(VBM8baddirs) = [];
%   
%   oldpath    = textscan(path,'%s','bufsize',2^16,'Delimiter',':'); oldpath=oldpath{1};
%   vbm12pathi = ~cellfun('isempty',strfind(oldpath,VBM12path));
%   
%   olddir = pwd; cd(VBM8dir);
%   [pp,ff,ee]=spm_fileparts(file);
% 
%   for pi=1:numel(vbm12pathi), if vbm12pathi(pi), rmpath(oldpath{pi}); end; end
%   for pi=numel(VBM8dirs):-1:1, addpath(VBM8dirs{pi},'-begin'); end

%---------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev$)
%---------------------------------------------------------------------

  pathchanges = replaceSPMpath(SPMwkd,VBM8dir);

  matlabbatch{1}.spm.tools.vbm8.estwrite.data                   = cellstr(file);
  matlabbatch{1}.spm.tools.vbm8.estwrite.opts.tpm               = cellstr(fullfile(spm('Dir'),'toolbox','Seg','TPM.nii'));
  matlabbatch{1}.spm.tools.vbm8.estwrite.opts.ngaus             = [2 2 2 3 4 2];
  matlabbatch{1}.spm.tools.vbm8.estwrite.opts.biasreg           = 0.0001;
  matlabbatch{1}.spm.tools.vbm8.estwrite.opts.biasfwhm          = 30;
  matlabbatch{1}.spm.tools.vbm8.estwrite.opts.affreg            = 'mni';
  matlabbatch{1}.spm.tools.vbm8.estwrite.opts.warpreg           = 4;
  matlabbatch{1}.spm.tools.vbm8.estwrite.opts.samp              = 3;
  
  matlabbatch{1}.spm.tools.vbm8.estwrite.extopts.dartelwarp.normlow = struct([]);
  matlabbatch{1}.spm.tools.vbm8.estwrite.extopts.sanlm          = 2;
  matlabbatch{1}.spm.tools.vbm8.estwrite.extopts.mrf            = 0.15;
  matlabbatch{1}.spm.tools.vbm8.estwrite.extopts.cleanup        = 1;
  matlabbatch{1}.spm.tools.vbm8.estwrite.extopts.print          = 1;
  
  matlabbatch{1}.spm.tools.vbm8.estwrite.output.GM.native       = 0;
  matlabbatch{1}.spm.tools.vbm8.estwrite.output.GM.warped       = 0;
  matlabbatch{1}.spm.tools.vbm8.estwrite.output.GM.modulated    = 0;
  matlabbatch{1}.spm.tools.vbm8.estwrite.output.GM.dartel       = 0;
  matlabbatch{1}.spm.tools.vbm8.estwrite.output.WM.native       = 0;
  matlabbatch{1}.spm.tools.vbm8.estwrite.output.WM.warped       = 0;
  matlabbatch{1}.spm.tools.vbm8.estwrite.output.WM.modulated    = 0;
  matlabbatch{1}.spm.tools.vbm8.estwrite.output.WM.dartel       = 0;
  matlabbatch{1}.spm.tools.vbm8.estwrite.output.CSF.native      = 0;
  matlabbatch{1}.spm.tools.vbm8.estwrite.output.CSF.warped      = 0;
  matlabbatch{1}.spm.tools.vbm8.estwrite.output.CSF.modulated   = 0;
  matlabbatch{1}.spm.tools.vbm8.estwrite.output.CSF.dartel      = 0;
  matlabbatch{1}.spm.tools.vbm8.estwrite.output.bias.native     = 1;
  matlabbatch{1}.spm.tools.vbm8.estwrite.output.bias.warped     = 0;
  matlabbatch{1}.spm.tools.vbm8.estwrite.output.bias.affine     = 0;
  matlabbatch{1}.spm.tools.vbm8.estwrite.output.label.native    = 1;
  matlabbatch{1}.spm.tools.vbm8.estwrite.output.label.warped    = 0;
  matlabbatch{1}.spm.tools.vbm8.estwrite.output.label.dartel    = 0;
  matlabbatch{1}.spm.tools.vbm8.estwrite.output.jacobian.warped = 0;
  matlabbatch{1}.spm.tools.vbm8.estwrite.output.warps           = [0 0];

  
  [pp,ff,ee]=spm_fileparts(file);
  warning off;
  try
    clear global;
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);  

    delete(fullfile(pp,['p' ff '_seg8.txt']));
    delete(fullfile(pp,[ff '_seg8.mat']));
  catch e %#ok<NASGU>
    createNullImage(fullfile(pp,[ff ee]),fullfile(pp,['p0' ff ee]));
    createNullImage(fullfile(pp,[ff ee]),fullfile(pp,['m' ff ee]));
  end
  warning on;

  restoreSPMpath(pathchanges)
    
end
function createNullImage(FO,FN)
  if ~exist(FN,'file')
    V = spm_vol(FO);
    V.fname = FN;
    spm_write_vol(V,zeros(V.dim(1:3),'uint8'));
  end
end
function pathchanges = replaceSPMpath(oldSPMpath,newSPMpath)
% remove old SPM-directories from path and add the new SPM directories
  newSPMpath = cellstr(newSPMpath); 
  
  pathchanges.add = {};
  pathchanges.rem = {};
  pathchanges.olddir = pwd; 
 
  if ~strcmp(oldSPMpath,newSPMpath{1})
    newSPMpath=[newSPMpath;vbm_findfiles(newSPMpath,'*',struct('maxdepth',2,'dirs',1))];
    
    % find all directories in the old SPM path
    oldpath = textscan(path,'%s','bufsize',2^16,'Delimiter',':'); oldpath=oldpath{1};
    oldSPMpathi = ~cellfun('isempty',strfind(oldpath,oldSPMpath));
    for pi=1:numel(oldSPMpathi)
      if oldSPMpathi(pi)
        rmpath(oldpath{pi}); 
        pathchanges.rem{end+1} = oldpath{pi}; 
      end 
    end

    % add new paths
    oldpath2 = textscan(path,'%s','bufsize',2^16,'Delimiter',':'); oldpath2=oldpath2{1};
    for pi=numel(newSPMpath):-1:1
      %if all(cellfun('isempty',strfind(oldpath2,newSPMpath{pi})))
      warning off
        addpath(newSPMpath{pi},'-begin');
        pathchanges.add{end+1} = newSPMpath{pi};
      warning on
      %end
    end
  end
  cd(newSPMpath{1});
   
  % call SPM and complete paths
  %spm('defaults','fMRI'); %spm_check_installation('basic'); 
end
function restoreSPMpath(pathchanges)
% restore path bevor using 'replaceSPMpath' function
  
  for pi=1:numel(pathchanges.rem)
    addpath(pathchanges.rem{pi}); 
  end
  
  for pi=1:numel(pathchanges.add)
    warning off
    rmpath(pathchanges.add{pi}); 
    warning on  
  end
  
  cd(pathchanges.olddir );

end
function varargout = cat_io_dcm2bids_protocols(action,varargin)
%cat_io_dcm2bids_protocols. Sites, protocols, and overview tables of cat_io_dcm2bids.
%  The site dictionary (getSites, setupSites) and the protocol definitions 
%  (getProtocols, JSON files with the relevant MR parameters) are read once.
%  Each scan is tested against the protocols (testProtocols) and added to 
%  the overview tables of the database (updateDBtables, catDCM2BIDSdb-tables
%  with subjects, studies, and protocols, including shortened protocol JSON
%  files that can be used to define new protocols). cleanupVjson removes 
%  the patient-specific fields of a sidecar or keeps only the protocol (or 
%  BIDS) relevant fields.
%
%  varargout = cat_io_dcm2bids_protocols(action,varargin)
%
%  Actions (see the help of the local functions):
%    sites = cat_io_dcm2bids_protocols('getSites',Pcenterdict)
%    protocols = cat_io_dcm2bids_protocols('getProtocols',Pprodictdirs)
%    site = cat_io_dcm2bids_protocols('setupSites',Vjson,sites)
%    [match,pname,mismatchstr,Pnfailedid,Pfname] = cat_io_dcm2bids_protocols('testProtocols',Pjson,Vjson,Ptbldir,protocols,tol,job,sites)
%    Vjson = cat_io_dcm2bids_protocols('cleanupVjson',Vjson,type)
%    cat_io_dcm2bids_protocols('updateDBtables',Vjson,Pnii,Pdbdir,Ptbldir,fnameparts,sites,opts)
%
%  See also cat_io_dcm2bids.

  switch action
    case {'getSites','getProtocols','setupSites','testProtocols','cleanupVjson', ...
          'updateDBtables'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_protocols:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function sites = getSites(Pcenterdict)
  sites = {}; 
  if ~isempty(Pcenterdict) 
    for di = 1 % to support more you would have to match the fields 
      if ~isempty(Pcenterdict{di}) 
        if ~exist(Pcenterdict{di},'file')
          error('cat_io_dcm2bids:Pcenterdict','Center dictionary file "%s" is not existing.', Pcenterdict{di});
        end
        sites = [sites; struct2cell(cat_io_json(Pcenterdict{1}))']; %#ok<AGROW>
      end
    end
  end
end
% =========================================================================
function protocols = getProtocols(Pprodictdirs)
  protocols = cell(numel(Pprodictdirs),3); pdi = 0; 
  if isempty(Pprodictdirs) || isempty(Pprodictdirs{1}), return; end
  for di = 1:numel(Pprodictdirs)
    if isempty(Pprodictdirs{di}), continue; end  
    if ~exist(Pprodictdirs{di},'dir')
      error('cat_io_dcm2bids:protocoldir','Protocol directory %d  "%s" does not exist.\n',di,Pprodictdirs{di}); 
    end
    Pprodictsubdirs = cat_vol_findfiles(Pprodictdirs{di},'*.json');
    Pprodictsubdirs(cat_io_contains(Pprodictsubdirs,[filesep 'qc'])) = []; 
    if isempty(Pprodictsubdirs) % 'cat_io_dcm2bids:noProtocols',
      cat_io_cprintf('warn','No protocols found in:\n  %s\n',Pprodictdirs{di});
    end
    for fi = 1:numel(Pprodictsubdirs)
      pdi = pdi + 1; 
      protocols{pdi,1} = spm_file(Pprodictdirs{di},'basename'); 
      protocols{pdi,2} = Pprodictsubdirs{fi}; 
      protocols{pdi,3} = cat_io_json(Pprodictsubdirs{fi});
    end
  end
  protocols = protocols(1:max(1,pdi),:); % remove empty rows of directories without protocols
  if isempty(protocols{1}) % 'cat_io_dcm2bids:noProtocols',
    cat_io_cprintf('warn','No protocols found in any protocol directory:\n  %s\n', ...
      strjoin(Pprodictdirs(:)',sprintf('\n  ')));
  end
end
% =========================================================================
function site = setupSites(Vjson,sites)
  if ~isempty(sites)
    DevNam = cat_io_contains( sites(:,1) , Vjson.DeviceSerialNumber );
  else
    DevNam = 0; 
  end
  if any(DevNam)
    site = sites{find(DevNam,1,'first'),2};
  else
    site = Vjson.DeviceSerialNumber; 
  end
  site = cat_io_dcm2bids_bids('bidsLabel',site,1); 
end
% =========================================================================
function [match,pname,mismatchstr,Pnfailedid,Pfname] = testProtocols(Pjson,Vjson,Ptbldir,protocols,tol,job,sites)
% test protocols

  if isempty(protocols) || isempty(protocols{1,1})
  % no given protocols
    pmatchs         = inf; 
    pmatch          = 0; 
    mismatchstr{1}  = cell(0,3); 
    mismatchcnt(1)  = 0; 
    
  else
    pmatch      = ones(1,size(protocols,1)); 
    pmatchs     = zeros(1,size(protocols,1)); 
    mismatchstr = cell(1,size(protocols,1)); 
    mismatchcnt = zeros(1,size(protocols,1)); 
    for pri = 1:size(protocols,1)
      mismatchstr{pri} = cell(0,3); 
      FNpi = fieldnames(protocols{pri,3});
     % FNpi(cat_io_contains(FNpi,{'ConsistencyInfo','PulseSequenceDetails', ...
     %   'ImageComments','SequenceName','ProtocolName','SeriesDescription'})) = []; 
      
      pmatchs(pri) = numel(FNpi); 
      for fni = 1:numel(FNpi)
        if isfield( Vjson, FNpi{fni} ) 
          if strncmp(FNpi{fni},'x_',2), continue; end % comment
          
          if ( islogical( Vjson.(FNpi{fni})) || isnumeric( Vjson.(FNpi{fni})) ) && ...
             ( islogical( protocols{pri,3}.(FNpi{fni})) || isnumeric( protocols{pri,3}.(FNpi{fni})) )
            if numel( Vjson.(FNpi{fni}) ) == numel( protocols{pri,3}.(FNpi{fni}))
              matstr =  all( Vjson.(FNpi{fni}) >= protocols{pri,3}.(FNpi{fni})*(1-tol/100) & ...
                             Vjson.(FNpi{fni}) <= protocols{pri,3}.(FNpi{fni})/(1-tol/100)); 
              %matstr = max(0,min(1, max(0,abs( Vjson.(FNpi{fni}) - protocols{pri,3}.(FNpi{fni}) ) - (protocols{pri,3}.(FNpi{fni})*tol/100) ))); 
              mismatchcnt(pri) = mismatchcnt(pri) + (1-matstr) * .5; 
              pmatchn = matstr;
            else
              mismatchcnt(pri) = mismatchcnt(pri) + 0.5;
              pmatchn = 0; 
            end
          elseif ischar( Vjson.(FNpi{fni}) ) && ischar( protocols{pri,3}.(FNpi{fni}) )
            matstr = strcmp( Vjson.(FNpi{fni}) , protocols{pri,3}.(FNpi{fni}) );
            mismatchcnt(pri) = mismatchcnt(pri) + (1-matstr);
            pmatchn = matstr; 
          elseif iscell( Vjson.(FNpi{fni}) ) && iscell( protocols{pri,3}.(FNpi{fni}) )
            %pmatchn = strcmp( char(sort(protocols{pri,3}.(FNpi{fni}(:)))) , char(sort(Vjson.(FNpi{fni}(:)))) ); 
            matstr  = 0; 
            for pmi = 1:numel( Vjson.(FNpi{fni}) )
               matstr = matstr + ...
                (1 - max([0;cat_io_contains(protocols{pri,3}.(FNpi{fni}),Vjson.(FNpi{fni})(pmi))])); 
            end
            for pmi = 1:numel( protocols{pri,3}.(FNpi{fni}) )
              matstr = matstr + ...
                (1 - max([0;cat_io_contains(Vjson.(FNpi{fni}),protocols{pri,3}.(FNpi{fni})(pmi))])); 
            end
            mismatchcnt(pri) = mismatchcnt(pri) + max(0,matstr); % - tol/2);
            pmatchn = mismatchcnt(pri) < 1;
          else
            mismatchcnt(pri) = mismatchcnt(pri) + 1;
            pmatchn = 0; 
            % need refinement ! ... image type             
          end
          pmatch(pri) = pmatch(pri) & pmatchn; 
  
          if ~pmatchn
            if iscell( Vjson.(FNpi{fni}) ) 
              tstr1 = char(join(sort(Vjson.(FNpi{fni}))));
              tstr2 = char(join(sort(protocols{pri,3}.(FNpi{fni})))); 
            else
              if isscalar(Vjson.(FNpi{fni}))
                tstr1 = Vjson.(FNpi{fni});
              elseif ischar(Vjson.(FNpi{fni})) || isstring(Vjson.(FNpi{fni}))
                tstr1 = join(Vjson.(FNpi{fni}),1); 
              else
                tstr1 = sprintf('%dx%d %s', size(Vjson.(FNpi{fni}),1), ...
                  size(Vjson.(FNpi{fni}),2), class(Vjson.(FNpi{fni}))); 
              end
              if isscalar(protocols{pri,3}.(FNpi{fni}))
                tstr2 = protocols{pri,3}.(FNpi{fni});
              elseif ischar(protocols{pri,3}.(FNpi{fni})) || isstring(protocols{pri,3}.(FNpi{fni}))
                tstr2 = join(protocols{pri,3}.(FNpi{fni}),1); 
              else
                tstr2 = sprintf('%dx%d %s', size(protocols{pri,3}.(FNpi{fni}),1), ...
                  size(protocols{pri,3}.(FNpi{fni}),2), class(protocols{pri,3}.(FNpi{fni}))); 
              end
            end
            mismatchstr{pri} = [ mismatchstr{pri}; {FNpi{fni} tstr1 tstr2 }]; 
          end

          pmatch(pri) = pmatch(pri) & pmatchn; 
        end
      end
    end
  end
%% always print session 

  % display progress
  % =======================================================================
  fname1     = sprintf('%3d) %s', Vjson.SeriesNumber, spm_str_manip(strrep(sprintf('%s_%s_%s', ...
                strrep(Vjson.PatientID, strrep(char(Vjson.ScanDate),'-',''),''), ...
                strrep(char(Vjson.ScanDate),'-',''), Vjson.ProtocolName),'__','_'),'l55'));
  if isempty(protocols) || isempty(protocols{1,1})
  % no given protocols: no match but use the DICOM protocol name
    match      = 0;
    Pnfailedid = [];
    pname      = Vjson.ProtocolName;
    Pfname     = '';
    fprintf('%60s : ',fname1);
    datatype = cat_io_dcm2bids_bids('setupDatatype',pname);
    cat_io_cprintf([0 0 0.5],sprintf('%-50s%10s ', [datatype filesep pname], ''));
    return
  elseif all( mismatchcnt > 0.05 )
  % unknown protocol
    Pnfailed   = mismatchcnt; %cellfun(@(x) size(x,1),mismatchstr);
    Pnfailedid = find(Pnfailed == min(Pnfailed));% & Pnfailed < 6);
    pname      = Vjson.ProtocolName;
    Pfname     = '';
    fprintf('%60s : ',fname1);
    pname0 = spm_file(protocols{Pnfailedid(1),2},'basename');
    if min(Pnfailed) < 1
      % very close protocol
      Pfname     = protocols{Pnfailedid(1),2}; 
      cat_io_cprintf([.5 .5 0],sprintf('%-50s%10s ',pname0,'~'));
    elseif min(Pnfailed) < 6 
    % close protocol
      Pfname     = protocols{Pnfailedid(1),2}; 
      cat_io_cprintf([1 .5 0],sprintf('%-50s%10s ', ...
        pname0, sprintf('%2.0f/%2.0f',min(Pnfailed), numel(Pnfailed))));
    else
      datatype = cat_io_dcm2bids_bids('setupDatatype',pname);
      cat_io_cprintf([.7 0 0],sprintf('%-50s%10s ', ...
        sprintf('Unknown %s protocol',datatype), ...
        sprintf('%2.0f/%2.0f',min(Pnfailed), numel(Pnfailed))));
    end
  else
    Pnfailedid = {}; 
    %pname  = spm_file(protocols{find(pmatch==1 & max(pmatchs.*pmatch)==pmatchs,1,'first'),2},'basename'); 
    %Pfname = protocols{find(pmatch==1 & max(pmatchs.*pmatch)==pmatchs,1,'first'),2}; 
    pname  = spm_file(protocols{find(mismatchcnt == min(mismatchcnt),1,'first'),2},'basename'); 
    Pfname = protocols{find(mismatchcnt == min(mismatchcnt),1,'first'),2}; 
    pname0 = spm_str_manip( pname, 'l50');
    fprintf('%60s : ',fname1); 
    cat_io_cprintf([0 .5 0],sprintf('%-50s%10s ',pname0,''));
  end
  % match = 0 for no fitting protocol, otherwise the index of the fitting protocol
  [mincnt,matchid] = min(mismatchcnt);
  match = (mincnt <= 0.05) * matchid;

  % save non-fitting protocols
  % =======================================================================
  % This should support central integration of protocols. 
  % So we create a directory with the full and a shortened version of the protocol. 
  % The shortened version includes all fields used so far in existing protocols. 
  % In case of close protocols we copy and adapt these and create a difference table.
  % Moreover, we create a list of all cases to see how often the protocol is used.
  % =======================================================================
  
  %% Pjson,Vjson,protocols,tol,opts
  if isempty(protocols{1,1})
    proSubdir = 'undefined'; 
  elseif max(0,1 - min(mismatchcnt)) > .95
    proSubdir = 'conform';
  elseif max(0,1 - min(mismatchcnt)) > .5
    proSubdir = 'accepted';
  elseif min(Pnfailed) < 6
    proSubdir = 'semiconform';
  else    
    proSubdir = 'nonconform'; 
  end

  % evaluate protocols
  proname = strrep(strrep(Vjson.ProtocolName,'_','-'),' ','-'); 
  datatype = cat_io_dcm2bids_bids('setupDatatype',proname); 

  % study


  site     = setupSites(Vjson,sites);
  Pprodir  = fullfile(Ptbldir,'protocols-details',datatype,proSubdir,proname);
  if ~exist(Pprodir,'dir'), mkdir(Pprodir); end

  Vjsonc   = cleanupVjson(Vjson); 
  if ~match
    if min(Pnfailed) < 6
    % semi-match 
      for fi = 1:numel(Pnfailedid)
        Pprosubdir = fullfile(Pprodir,spm_file(protocols{Pnfailedid(fi),2},'basename')); 
        if ~exist(Pprosubdir,'dir'), mkdir(Pprosubdir); end
        Pprofilel = fullfile(Pprosubdir, ...
          sprintf('%s-%s-long.json', spm_file(protocols{Pnfailedid(fi),2},'basename'),site));
        cat_io_json(Pprofilel,Vjsonc);

        Pprofile = fullfile(Pprosubdir, ...
          sprintf('%s-%s.json', spm_file(protocols{Pnfailedid(fi),2},'basename'),site));
        % In this case the most similar protocol is used as template
        FNpi = intersect(fieldnames(protocols{Pnfailedid(fi),3}),fieldnames(Vjson));
        for fni = 1:numel(FNpi)
          Vjson3.(FNpi{fni}) = Vjson.(FNpi{fni});
        end
        cat_io_json(Pprofile,Vjson3); clear Vjson3;

        Pprofiled = fullfile(Pprosubdir,sprintf('%s-%s-diff.csv', ...
          spm_file(protocols{Pnfailedid(fi),2},'basename'),site));
        cat_io_csv(Pprofiled,mismatchstr{Pnfailedid(fi)});

        Pqc = spm_file(protocols{Pnfailedid(fi),2},'prefix','qc');
        if exist(Pqc,'file')
          Pprofileqc = fullfile(Pprosubdir,sprintf('qc%s-%s.json', ...
            spm_file(protocols{Pnfailedid(fi),2},'basename'),site));
          copyfile(Pqc,Pprofileqc);
        end
      end
    end
  end


  % create a list of cases to see how often the protocol is used
  Pprofilec = fullfile(Pprodir,sprintf('%s-%s-list.csv',proname,site));
  if exist(Pprofilec,'file')
    Yc = cat_io_csv(Pprofilec,'','',struct('convert2double',0));
    Yc{end+1,1} = Pjson;
    Yc = unique(Yc);
  else
    Yc{1,1} = Pjson; 
  end
  cat_io_csv(Pprofilec,Yc);
  
  % save the full unknown protocol
  Pprofilel = fullfile(Pprodir,sprintf('%s-%s-long.json',proname,site));
  cat_io_json(Pprofilel,Vjsonc);

  % save a shortened version as starting point   
  Pprofile  = fullfile(fileparts(Pjson),sprintf('%s-%s.json',proname,site));
  Pprofile2 = fullfile(Pprodir,sprintf('%s-%s.json',proname,site));
  FN = cellfun(@(x) fieldnames(x), protocols(:,3) , 'UniformOutput', false );
  FNM = {}; for fni=1:numel(FN), FNM = [FNM; FN{fni}]; end; FNM = unique(FNM); %#ok<AGROW>
  FNM = intersect(FNM,fieldnames(Vjson));
  
  if isempty(FNM), return; end

  for fni = 1:numel(FNM)
    Vjson4.(FNM{fni}) = Vjson.(FNM{fni});
  end
  Vjson4 = cleanupVjson(Vjson4); 
  cat_io_json(Pprofile, Vjson4);
  cat_io_json(Pprofile2,Vjson4);

end
% =========================================================================
function Vjson = cleanupVjson(Vjson,type)
% remove Patient specific entries. type={'basic'*|'protocol'}.

  if ~exist('type','var'), type = 'basic'; end

  FN    = fieldnames(Vjson); 
  Vjson = rmfield(Vjson, FN(cat_io_contains(FN,'Patient')));  
  
  % maybe 2-3 levels of data reduction ?
  % (0) all DCM2NII, (1) light cleanup, (2) severe cleanup

  switch type
    case 'basic'
      % negative list, i.e. remove these entries
      RFn = {
        'SeriesInstanceUID'; 'StudyInstanceUID'; 'StudyID'; 
        'InstitutionName'; 'InstitutionalDepartmentName'; 'InstitutionAddress'; 
        'StationName'; 
        ... 'ProcedureStepDescription'; 
        'BodyPartExamined';
        'AcquisitionTime'; 
        ... 'AcquisitionDateTime'
        'ShimSetting'; 'TxRefAmp'; 
        }; 
    case {'protocol','bids'}
      % positive list, i.e. keep only these entries
      % this should be for an overview of various files
      RFn = {
        'Modality'; 
        'MagneticFieldStrength'; 
        'Manufacturer';
        'PatientPosition';
        'MRAcquisitionType';
        'ScanningSequence';
        'SequenceVariant';
        'ScanOptions';
        'ImageType';
        'NonlinearGradientCorrection';
        'SliceThickness';
        'EchoTime';
        'RepetitionTime';
        'SpoilingState';
        'InversionTime';
        'FlipAngle';
        'PartialFourier';
        'BaseResolution';
        'DiffusionScheme';
        'PhaseResolution';
        'CoilString';
        'RefLinesPE';
        'CoilCombinationMethod';
        'MatrixCoilMode';
        'PercentPhaseFOV';
        'PercentSampling';
        'PhaseEncodingSteps';
        'AcquisitionMatrixPE';
        'DwellTime';
        'ReconMatrixPE';
        'InPlanePhaseEncodingDirectionDICOM';
        'DerivedVendorReportedEchoSpacing';
        ... 
        ... 'PulseSequenceName'; % not always and operator specific
        'ProcedureStepDescription';
        ...
        'PhaseEncodingDirection';
        ...'SliceTiming';
        'InPlanePhaseEncodingDirectionDICOM';
        'PixelBandwidth';
        'TotalReadoutTime';
        'EffectiveEchoSpacing';
        'DerivedVendorReportedEchoSpacing';
        'ParallelReductionFactorInPlane';
        'BandwidthPerPixelPhaseEncode';
        'EchoTrainLength';
        'MultibandAccelerationFactor';
        }; 
      if strcmp(type,'bids')
        % further technical (not identifying) fields of the BIDS sidecars, 
        % e.g. required by BIDS (TaskName for func, EchoTime1/2 for phasediff,
        % ASL and qMRI parameters)
        RFn = [RFn; {
          'TaskName'; 'SliceTiming'; 'SliceEncodingDirection'; 'SpacingBetweenSlices'; 
          'EchoTime1'; 'EchoTime2'; 'RepetitionTimeExcitation'; 'RepetitionTimePreparation'; 
          'NumberOfVolumesDiscardedByScanner'; 'NumberOfVolumesDiscardedByUser'; 'DelayTime'; 
          'ManufacturersModelName'; 'SoftwareVersions'; 'ReceiveCoilName'; 'PulseSequenceType'; 
          'SequenceName'; 'B0FieldIdentifier'; 'B0FieldSource'; 'MTState'; 'NumberShots'; 
          'ArterialSpinLabelingType'; 'PostLabelingDelay'; 'BackgroundSuppression'; 'M0Type'; 
          'TotalAcquiredPairs'; 'LabelingDuration'; 'PCASLType'; 'VascularCrushing'; 
          'AcquisitionVoxelSize'; 'BolusCutOffFlag'; 'LookLocker'; 'Units'; 
          'IntendedFor'}]; % remapped to the new file names (see remapIntendedFor in cat_io_dcm2bids_bids)
      end
      RFn = setdiff( fieldnames(Vjson), RFn );   
  end
  RFn   = intersect( RFn , fieldnames(Vjson) ); 
  Vjson = rmfield(Vjson,RFn); 
end
% =========================================================================
function updateDBtables(Vjson, Pnii, Pdbdir, Ptbldir, fnameparts, sites, opts) %#ok<INUSL> Pdbdir for the disabled scan/session lists
%updateDBtables. Update the overview tables of the database (subjects, 
%  studies, protocols) for a scan of the database directory Pdbdir 
%  (sub-*/ses-*/snr-*) with its JSON sidecar information Vjson and its 
%  image Pnii (if available, e.g. not for JSON-only outputs). 

%datetime( Vjson.AcquisitionDateTime , 'Format','yyyyMMdd-hhmmss'))

  % DB lists 
  Pstudytable    = spm_file(fullfile(Ptbldir,'studies'),  'ext', opts.tableformat); 
  Psubjecttable  = spm_file(fullfile(Ptbldir,'subjects'), 'ext', opts.tableformat); 
  Pprotocoltable = spm_file(fullfile(Ptbldir,'protocols'),'ext', opts.tableformat); 

  if isfield(Vjson,'ProcedureStepDescription')
    ProcedureStepDescription = getFileString(Vjson.ProcedureStepDescription);
  else
    ProcedureStepDescription = getFileString(Vjson.SeriesDescription);
  end
  % SCAN-LIST:
  %  - to avoid double entries ... but what to do ...
  %  - a flag (keep old / take new)
%%%%%%%  * this grows too strong >> internal variable + mat     
  if 0
    Thdr    = {'SeriesInstanceUID','DBpath'}; %,'IMPORTpath'}; 
    Tnewrow = {Vjson.SeriesInstanceUID, Pdbdir}; %, Pjson}; 
    cat_io_dcm2bids_helper('updateTable',Pscantable,Thdr,Tnewrow,1,opts.rerun);
  end

  % SESSION-LIST:
  %  - to avoid double entries ... but what to do ...
  %  - a flag (keep old / take new)
%%%%%%%  * this one is useless      
  if 0
    Thdr    = {'StudyInstanceUID','DBpath','ProcedureStepDescription'}; %,'IMPORTpath'}; 
    Tnewrow = {Vjson.StudyInstanceUID, spm_fileparts(Pdbdir), Vjson.ProcedureStepDescription}; %, Pjson}; 
    cat_io_dcm2bids_helper('updateTable',Psessiontable,Thdr,Tnewrow,1,opts.rerun);
  end


  % SUBJECT-LIST ("constant" data):
%%%%%%%  * this one is nice but a session counter would be nice 
%%%%%%%  * maybe the scan period 
%%%%%%%  * import path (not robust :/)
  site    = setupSites(Vjson,sites);
  Thdr    = {'PatientID','PatientName', 'PatientSex','PatientBirthDate','SiteID'}; 
  Tnewrow = {Vjson.PatientID, Vjson.PatientName, ...
             Vjson.PatientSex, Vjson.PatientBirthDate, ...
             site}; 
  cat_io_dcm2bids_helper('updateTable',Psubjecttable,Thdr,Tnewrow,1,opts.rerun);
  

  % STUDY-LIST (defined by ProcedureStepDescription with additional list of subjects):
  % - ProcedureStepDescription with subjects
  Pstudysubtable = spm_file(fullfile(Ptbldir,'study_subject_lists',ProcedureStepDescription),'ext',opts.tableformat); 
  Thdr    = {'Subjects'}; 
  Tnewrow = {Vjson.PatientID}; 
  nsub    = cat_io_dcm2bids_helper('updateTable',Pstudysubtable,Thdr,Tnewrow,1,opts.rerun);
  % - study-list with number of subjects
  Thdr    = {'ProcedureStepDescription','nSubjects'}; 
  Tnewrow = {ProcedureStepDescription,nsub}; 
  cat_io_dcm2bids_helper('updateTable',Pstudytable,Thdr,Tnewrow,1,opts.rerun);
% with counters for anat, func, fmap, dwi ?      
  

  % PROTOCOLS-LIST:
  datatype0         = cat_io_dcm2bids_bids('setupDatatype',Vjson.ProtocolName); 
  ProtocolName      = getFileString(strrep(fnameparts{2},'_','-'));
  %Pprotocolsubtable = spm_file(fullfile(Ptbldir,ProcedureStepDescription,datatype0,ProtocolName),'ext',opts.tableformat); 
  Pprotocolsubtable = spm_file(fullfile(Ptbldir,'protocols',datatype0,ProtocolName),'ext',opts.tableformat); 

  % save shortened protocol as json to be used as filter, with the image 
  % information (Dimensions [x y z volumes], VoxelSize [x y z]) if the
  % image is available, which is not used to separate protocols 
  Vjson0 = cleanupVjson(Vjson,'protocol'); 
  IFN    = {'Dimensions','VoxelSize'}; 
  for fni = 1:numel(IFN)
    if isfield(Vjson,IFN{fni}), Vjson0.(IFN{fni}) = Vjson.(IFN{fni}); end
  end
  Pjson0 = spm_file( Pprotocolsubtable,'ext','.json'); 
  ProtocolName0 = ProtocolName; 
  BFN = {'ProcedureStepDescription'};
  Vjson0 = rmfield(Vjson0,intersect(fieldnames(Vjson0),BFN)); 
  if ~exist(Pjson0,'file')
    cat_io_json( Pjson0, Vjson0); 
  else
    %Pjsons = cat_vol_findfiles( fullfile(Ptbldir,ProcedureStepDescription,datatype0)  , '*.json');
    Pjsons = cat_vol_findfiles( fullfile(Ptbldir,'protocols',datatype0)  , '*.json');
    Pfit = false(size(Pjsons)); 
    for pi = 1:numel(Pjsons)
      Vjson1   = cat_io_json( Pjsons{pi} ); 
      Vjson1   = rmfield(Vjson1,intersect(fieldnames(Vjson1),BFN)); 
      Pfit(pi) = cat_io_structEqual( rmfield(Vjson0,intersect(fieldnames(Vjson0),IFN)), ...
                              rmfield(Vjson1,intersect(fieldnames(Vjson1),IFN)) ); 
      % add the image information if it was not available before (JSON-only import)
      if Pfit(pi) && isfield(Vjson0,'VoxelSize') && ~isfield(Vjson1,'VoxelSize')
        Vjson1 = cat_io_json( Pjsons{pi} ); 
        for fni = 1:numel(IFN), Vjson1.(IFN{fni}) = Vjson0.(IFN{fni}); end
        cat_io_json( Pjsons{pi}, Vjson1 ); 
      end
    end
    if any( Pfit )
      % name of the matching protocol (variant), preferably of this protocol
      fiti = find( Pfit & startsWith(spm_file(Pjsons,'basename'),ProtocolName), 1); 
      if isempty(fiti), fiti = find(Pfit,1); end
      ProtocolName0 = spm_file(Pjsons{fiti},'basename'); 
    else
      ci = 0;
      while exist(Pjson0,'file')
        ci = ci+1;
        Pjson0 = spm_file(Pprotocolsubtable,'suffix',sprintf('_%d',ci),'ext','.json'); 
        ProtocolName0 = spm_file(ProtocolName,'suffix',sprintf('_%d',ci)); 
      end
      cat_io_json( Pjson0, Vjson0); 
    end
  end
  Thdr    = {'Subjects'}; 
  Tnewrow = {Vjson.PatientID}; 
  nsub    = cat_io_dcm2bids_helper('updateTable',Pprotocolsubtable,Thdr,Tnewrow,1,opts.rerun);

  % - study-list with number of subjects
  Thdr    = {'ProtocolName', 'nSubjects', ...
             'MagneticFieldStrength', 'Manufacturer', 'SoftwareVersions' ...
             'ProcedureStepDescription', ... % quite variant but may be useful to name the studies in future
             ...'PulseSequenceName', ... % quite variant
             'MRAcquisitionType', 'ScanningSequence', 'SequenceVariant', ...
             'ScanOptions', ...
             'ImageType', ...
             'NonlinearGradientCorrection',...
             'RepetitionTime', 'InversionTime', 'EchoTime', ...
             'FlipAngle', ...
             'BaseResolution', 'CoilString', 'CoilCombinationMethod', ...
             'MatrixCoilMode', 'ParallelReductionFactorInPlane', ...
             'PercentPhaseFOV','PercentSampling', ...
             'AcquisitionMatrixPE', 'ReconMatrixPE', ...
             'PixelBandwidth', 'DwellTime', ...
             'PhaseEncodingDirection', 'EchoTrainLength' ...
             ... 'BandwidthPerPixelPhaseEncode', 'EffectiveEchoSpacing'; 
             }; 
  Vjson0 = Vjson; 
  for fni = 1:numel(Thdr)
    if ~isfield(Vjson0,Thdr{fni}), Vjson0.(Thdr{fni}) = ''; end
    if iscell(Vjson0.(Thdr{fni}))
      Vjson0.(Thdr{fni}) = strjoin(Vjson0.(Thdr{fni})); 
    end
  end
  
  % add image dimensions (x, y, z, volumes) and voxel size (x, y, z) as 
  % separate columns (if the image is available)
  if exist(Pnii,'file')
    [vx,dim] = cat_io_dcm2bids_helper('niiHeaderInfo',Pnii);
  else
    vx = nan(1,3); dim = nan(1,4); 
  end
  for di = 1:4, Thdr{end+1} = sprintf('dim%d',di); Vjson0.(Thdr{end}) = dim(di); end %#ok<AGROW>
  for di = 1:3, Thdr{end+1} = sprintf('vx%d',di);  Vjson0.(Thdr{end}) = round(vx(di),4); end %#ok<AGROW>

  % new protocols, or with rerun replaced, or the image information of 
  % JSON-only imports added later 
  Tnewrow = {sprintf('%s_%s',datatype0,ProtocolName0),nsub}; 
  for fni = 3:numel(Thdr), Tnewrow{fni} = Vjson0.(Thdr{fni}); end
  if opts.rerun, reimport = 1; else, reimport = 2; end
  cat_io_dcm2bids_helper('updateTable',Pprotocoltable,Thdr,Tnewrow,1,reimport);


end
% =========================================================================
function str = getFileString(str)
  str = cat_io_strrep(str, ...
    {'ä','ü','ö',' '}, {'ae','ue','oe','_'});
  strd = double(str);
  str(strd<48 & strd~=45) = '_';
  str(strd>57 & strd<65)  = '_';
  str(strd>90 & strd<97 & strd~=95)  = '_';
  str(strd>122) = 'X';
end

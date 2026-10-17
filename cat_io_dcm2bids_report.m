function varargout = cat_io_dcm2bids_report(action,varargin)
%cat_io_dcm2bids_report. Report tables of cat_io_dcm2bids.
%  The report tables are created from the BIDS directory and its
%  derivatives (i.e., they can be created again at any time and fit to the
%  exported data) and contain only the BIDS labels of the subjects, i.e.,
%  they are for sharing and the evaluation by other software (e.g. R):
%    scans.<fmt>       one row per scan (BIDS image or sidecar) with the
%                      subject, session, site, age, sex, BIDS type, protocol
%                      test (closest protocol, status, PFM, PFR), quality
%                      measures and ratings (QC, coverage, phase-encoding
%                      reversal), orientation of the session, and main
%                      acquisition parameters
%    scans_<datatype>-<suffix>.<fmt>   the scans of one typical scan type
%                      (opts.scantypes) with type-specific results (e.g. the
%                      tissue volumes of the segmentation of anatomical scans)
%    sessions.<fmt>    one row per session (for the PI) with the combined
%                      ratings (protocol: worst PFR, quality: RMS of the best
%                      SQR of each scan type, orientation: SOR, worst coverage
%                      and reversal rating), the number of scans and ratings
%                      of the typical scan types, the tissue volumes of the
%                      best T1w scan (n/a if not available) with their 
%                      z-scores of the age/sex model of the sample and the 
%                      sample homogeneity of its mwc1 map (quartic mean 
%                      Z-score, see cat_io_dcm2bids_stats)
%    project.<fmt>     one row per dataset and site and for all sites of a
%                      dataset with the number of subjects, sessions, and
%                      scans, the demographics, median ratings, and the
%                      color of the site (site dictionary, see siteColors
%                      of cat_io_dcm2bids_protocols)
%  in <outdir>/<subdir>/<BIDSdir>-report/tables/<BIDSsubdir>/ with a JSON
%  file of each table with the description of its columns, where <fmt> is
%  opts.tableformat (csv or tsv) and missing values are n/a. The report 
%  figure (see cat_io_dcm2bids_reportfig) shows the tables as diagrams.
%  The regional results of each session (see writeRegionTables) are in
%  <outdir>/<subdir>/<BIDSdir>-report/regions/<BIDSsubdir>/:
%    <sub>_<ses>_anat.<fmt>  tissue volumes of the neuromorphometrics regions
%                            (segmentation of each anatomical scan)
%    <sub>_<ses>_dwi.<fmt>   FA and MD of the Mori white-matter regions
%    <sub>_<ses>_func.<fmt>  within-network connectivity of the 7 networks
%  with one row per scan and region and the JSON files anat.json, 
%  dwi.json, and func.json with the description of the columns.
%
%  varargout = cat_io_dcm2bids_report(action,varargin)
%
%  Actions (see the help of the local functions):
%    cat_io_dcm2bids_report('writeReportTables',Poutdir,job)
%
%  See also cat_io_dcm2bids, cat_io_dcm2bids_tables, cat_io_dcm2bids_stats,
%  cat_io_dcm2bids_reportfig.

  switch action
    case {'writeReportTables'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_report:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function writeReportTables(Poutdir,job)
%writeReportTables. Report tables of the BIDS directory of the job (with
%  its protocol subdirectories, see the help of this file). The tables are
%  created again in each run.
  PBIDS   = fullfile(Poutdir,job.BIDSdir,job.BIDSsubdir);
  Preport = fullfile(Poutdir,[job.BIDSdir '-report'],'tables',job.BIDSsubdir);
  Pold    = fullfile(Poutdir,[job.BIDSdir '-report'],sprintf('report_%s.%s',job.BIDSsubdir,job.opts.tableformat));
  if exist(Pold,'file'), delete(Pold); end % report table of older versions
  if exist(Preport,'dir'), rmdir(Preport,'s'); end
  S = collectScans(PBIDS);
  if isempty(S), return; end
  fmt = job.opts.tableformat;

  % scans (all and the typical scan types)
  [cols,desc] = scanColumns('');
  writeReport(fullfile(Preport,['scans.' fmt]),S,cols,desc,fmt);
  [types,tnam] = scanTypes(job.opts);
  for ti = 1:numel(types)
    St = S(strcmp({S.type},types{ti}));
    if isempty(St), continue; end
    [cols,desc] = scanColumns(types{ti});
    writeReport(fullfile(Preport,sprintf('scans_%s.%s',strrep(types{ti},'/','-'),fmt)),St,cols,desc,fmt);
  end

  % sessions and project
  E = sessionRecords(S,types,tnam);
  E = sessionStats(E);
  [cols,desc] = sessionColumns(types,tnam);
  writeReport(fullfile(Preport,['sessions.' fmt]),E,cols,desc,fmt);
  Q = projectRecords(S,E,job.subdir,types,tnam);
  Pdict = {}; if isfield(job,'dicts') && isfield(job.dicts,'Pcenterdict'), Pdict = job.dicts.Pcenterdict; end
  col = cat_io_dcm2bids_protocols('siteColors',Pdict,fullfile(Poutdir,'catDCM2BIDSdb'),{Q.site}); % (site dictionary)
  [Q.color] = col{:};
  [cols,desc] = projectColumns(types,tnam);
  writeReport(fullfile(Preport,['project.' fmt]),Q,cols,desc,fmt);
  fprintf('  Report tables: %s\n',Preport);

  % regional results of each session
  writeRegionTables(fullfile(Poutdir,[job.BIDSdir '-report'],'regions',job.BIDSsubdir),S,fmt);
end
% =========================================================================
function S = collectScans(PBIDS)
%collectScans. Records of all raw scans (sidecars of the MRI datatypes) of
%  the BIDS directory PBIDS (with its protocol subdirectories) with the
%  data of the participants and sessions files, the sidecar, and the
%  derivatives (QC, segmentation, and session registration).
  S = struct([]);
  if ~exist(PBIDS,'dir'), return; end
  P = cat_vol_findfiles(PBIDS,'sub-*.json');
  fs = regexptranslate('escape',filesep);
  P = P( cellfun('isempty',regexp(P,[fs '(derivatives|matchfiles)' fs],'once')) );
  P = P( ~cellfun('isempty',regexp(P,[fs '(anat|func|dwi|fmap|perf)' fs 'sub-[^' fs ']+\.json$'],'once')) );
  P = P( cellfun('isempty',regexp(spm_file(P,'basename'),'_(events|physio|stim|beh|scans|sessions)$','once')) );
  cache = containers.Map('KeyType','char','ValueType','any');
  for fi = 1:numel(P)
    [pp,name] = fileparts(P{fi});
    si    = strfind(P{fi},[filesep 'sub-']);
    root  = P{fi}(1:si(1)-1);
    parts = strsplit(name,'_');
    sub   = parts{1};
    ses   = regexp(name,'_(ses-[^_]+)','tokens','once'); if isempty(ses), ses = ''; else, ses = ses{1}; end
    ent   = regexp(name,'_((?:dir|echo|inv|part|flip|mt|rec|ce)-[^_]+)','tokens');
    ent   = strjoin(cellfun(@(x) ['_' x{1}],ent,'UniformOutput',false),'');
    sx    = parts{end};
    dt    = spm_file(pp,'basename');
    tok   = regexp(name,'_(part-(?:phase|real|imag)|inv-\d+)(?=_)','tokens');
    type  = [dt '/' sx strjoin(cellfun(@(x) ['_' x{1}],tok,'UniformOutput',false),'')];
    acq   = regexp(name,'_acq-(\d+)','tokens','once'); if isempty(acq), acq = ''; else, acq = acq{1}; end
    droot = fullfile(root,'derivatives',sub,ses,dt);

    % participants and sessions (cached)
    P0 = readTSVcached(cache,fullfile(root,'participants.tsv'));
    Ps = readTSVcached(cache,fullfile(root,sub,[sub '_sessions.tsv']));
    J  = readJSON(P{fi});
    Q  = readJSON(fullfile(droot,['catDCM2BIDSqc_' name '.json']));
    G  = readJSON(fullfile(droot,['catDCM2BIDSsegus_' name '.json']));
    if isempty(ses), sesname = sub; else, sesname = [sub '_' ses]; end
    A  = readJSONcached(cache,fullfile(root,'derivatives',sub,ses,['catDCM2BIDSaffine_' sesname '.json']));

    r = struct();
    r.dataset = spm_file(root,'basename');
    r.participant_id = sub; r.session_id = ses;
    r.site  = tsvValue(Ps,'session_id',ses,'site');
    r.age   = tsvNumber(Ps,'session_id',ses,'age');
    if isnan(r.age), r.age = tsvNumber(P0,'participant_id',sub,'age'); end
    r.sex   = sexCode(tsvValue(P0,'participant_id',sub,'sex'));
    r.datatype = dt; r.suffix = sx; r.entities = regexprep(ent,'^_',''); r.type = type;
    r.filename = strrep(fullfile(sub,ses,dt,name),filesep,'/'); 
    r.series   = jsonNumber(Q,{'protocol','SeriesNumber'}); % series of the scan (acq: first series of the acquisition)
    if isnan(r.series), r.series = str2double(acq); end
    r.protocol_dicom  = jsonPath(Q,{'protocol','ProtocolName'},'');
    r.protocol        = jsonPath(Q,{'protocol','protocol'},'');
    r.protocol_status = jsonPath(Q,{'protocol','status'},'');
    dev = jsonPath(Q,{'protocol','deviations'},{}); 
    if ischar(dev) || iscellstr(dev), dev = cellstr(dev); else, dev = {}; end %#ok<ISCLSTR>
    r.protocol_deviations = strjoin(dev,'; '); 
    if isempty(dev) && ~isempty(r.protocol), r.protocol_deviations = 'none'; end
    r.PFM = jsonNumber(Q,{'qualitymeasures','PFM'}); r.PFR = jsonNumber(Q,{'qualityrating','PFR'});
    r.SQR = jsonNumber(Q,{'qualityrating','SQR'});
    for fn = {'BSM','WSM','ISR','NSR','RES'}
      r.(fn{1}) = jsonNumber(Q,{'qualityrating',fn{1}});
      r.([fn{1} '_value']) = jsonNumber(Q,{'qualitymeasures',fn{1}});
    end
    r.BBM = jsonNumber(Q,{'qualitymeasures','BBM'}); r.BBR = jsonNumber(Q,{'qualityrating','BBR'});
    r.PRM = jsonNumber(Q,{'qualitymeasures','PRM'}); r.PRR = jsonNumber(Q,{'qualityrating','PRR'});
    for fn = {'OTM','ORM'}, r.(fn{1}) = jsonNumber(A,{'qualitymeasures',fn{1}}); end
    for fn = {'OTR','ORR','SOR'}, r.(fn{1}) = jsonNumber(A,{'qualityrating',fn{1}}); end
    r.manufacturer   = jsonPath(J,{'Manufacturer'},'');
    r.model          = jsonPath(J,{'ManufacturersModelName'},'');
    r.field_strength = jsonNumber(J,{'MagneticFieldStrength'});
    r.software       = jsonPath(J,{'SoftwareVersions'},'');
    r.TR = jsonNumber(J,{'RepetitionTime'}); r.TE = jsonNumber(J,{'EchoTime'});
    r.TI = jsonNumber(J,{'InversionTime'});  r.FA = jsonNumber(J,{'FlipAngle'});
    r.MB = jsonNumber(J,{'MultibandAccelerationFactor'});
    dim = [reshape(jsonPath(J,{'Dimensions'},[]),1,[]) nan(1,4)];
    vx  = [reshape(jsonPath(J,{'VoxelSize'},[]),1,[]) nan(1,3)];
    r.dim_x = dim(1); r.dim_y = dim(2); r.dim_z = dim(3); r.volumes = dim(4);
    if ~isnan(dim(1)) && isnan(dim(4)), r.volumes = 1; end
    r.vox_x = vx(1); r.vox_y = vx(2); r.vox_z = vx(3);
    % segmentation (anatomical scans)
    r.TIV  = jsonNumber(G,{'TIV'});   r.GMV  = jsonNumber(G,{'aGMV'});
    r.WMV  = jsonNumber(G,{'aWMV'});  r.CSFV = jsonNumber(G,{'aCSFV'});
    r.rGMV = jsonNumber(G,{'rGMV'});  r.rWMV = jsonNumber(G,{'rWMV'}); r.rCSFV = jsonNumber(G,{'rCSFV'});
    r.NCR  = jsonNumber(G,{'qc','NCR'}); r.TPMll = jsonNumber(G,{'qc','TPMll'});
    % DTI and rs-fMRI measures (see cat_io_dcm2bids_pp_dwi and cat_io_dcm2bids_pp_rsfmri)
    Jd = readJSON(fullfile(droot,['catDCM2BIDSdwi_' name '.json']));
    for fn = {'FA_WM','MD_WM','AD_WM','RD_WM','MD_WM_peakwidth','MD_brain'}, r.(fn{1}) = jsonNumber(Jd,{'global',fn{1}}); end
    r.DWI_outliers = jsonNumber(Jd,{'qc','outlier_volumes'});
    Jr = readJSON(fullfile(droot,['catDCM2BIDSrsfmri_' name '.json']));
    for fn = {'FD_mean','FD_max','FD_pct05','DVARS','tSNR'}, r.(fn{1}) = jsonNumber(Jr,{'qc',fn{1}}); end
    r.GCOR = jsonNumber(Jr,{'connectivity','GCOR'}); r.segregation = jsonNumber(Jr,{'connectivity','segregation'});
    r.FC_within = jsonNumber(Jr,{'connectivity','within'}); r.FC_between = jsonNumber(Jr,{'connectivity','between'});
    r.pplevel = max([jsonNumber(Jd,{'method','level'}) jsonNumber(Jr,{'method','level'})]); % (NaN without)
    % regional results (region tables, not columns)
    r.regional_anat = jsonPath(G,{'regional'},struct());
    r.regional_dwi  = jsonPath(Jd,{'regional'},struct());
    r.regional_func = jsonPath(Jr,{'connectivity','network'},struct());
    % normalized GM map, normalized FA map, and connectivity matrix (sample
    % homogeneity, not columns)
    r.Pmwc1 = ''; r.PwFA = ''; r.Pconn = '';
    for ext = {'.nii.gz','.nii'}
      if isempty(r.Pmwc1) && exist(fullfile(droot,['mwc1' name ext{1}]),'file'), r.Pmwc1 = fullfile(droot,['mwc1' name ext{1}]); end
      if isempty(r.PwFA) && exist(fullfile(droot,['catDCM2BIDSdwi_wFA_' name ext{1}]),'file'), r.PwFA = fullfile(droot,['catDCM2BIDSdwi_wFA_' name ext{1}]); end
    end
    if exist(fullfile(droot,['catDCM2BIDSrsfmri_conn_' name '.tsv']),'file'), r.Pconn = fullfile(droot,['catDCM2BIDSrsfmri_conn_' name '.tsv']); end
    if isempty(S), S = r; else, S(end+1) = r; end %#ok<AGROW>
  end
end
% =========================================================================
function writeRegionTables(Preg,S,fmt)
%writeRegionTables. Tables of the regional results of the scans S of each
%  session in the directory Preg (created again in each run, one 
%  subdirectory per dataset for several datasets): tissue volumes of the
%  neuromorphometrics regions (anat), FA and MD of the Mori white-matter 
%  regions (dwi), and the within-network connectivity of the 7 networks 
%  (func), with one row per scan and region (see the JSON files of the
%  columns, i.e. anat.json, dwi.json, and func.json). 
  if exist(Preg,'dir'), rmdir(Preg,'s'); end
  types = {
    'anat', 'regional_anat', {'scan','region','name','Vgm','Vwm','Vcsf'}, ...
      {'scan of the session (BIDS name without subject and session)', ...
       'region of the neuromorphometrics atlas of CAT (l/r for left/right)','name of the region', ...
       'gray matter volume in ml (native space, SPM segmentation, n/a for regions without this tissue)', ...
       'white matter volume in ml','cerebrospinal fluid volume in ml'}
    'dwi',  'regional_dwi',  {'scan','region','name','FA','MD','voxels'}, ...
      {'scan of the session (BIDS name without subject and session)', ...
       'white-matter region of the Mori atlas of CAT (l/r for left/right)','name of the region', ...
       'median fractional anisotropy in the white matter of the region (T1w segmentation)', ...
       'median mean diffusivity in mm2/s','number of voxels'}
    'func', 'regional_func', {'scan','network','FC_within'}, ...
      {'scan of the session (BIDS name without subject and session)', ...
       'network of Yeo et al. 2011 (7 networks of the Schaefer2018 parcels of CAT)', ...
       'mean within-network connectivity of the parcels (Fisher z)'}
    };
  mori = atlasNames('mori'); neuro = atlasNames('neuromorphometrics');
  multi = numel(unique({S.dataset})) > 1; 
  key = strcat({S.dataset},'|',{S.participant_id},'|',{S.session_id}); 
  for k = unique(key)
    Ss = S(strcmp(key,k{1}));
    if isempty(Ss(1).session_id), name = Ss(1).participant_id; else, name = [Ss(1).participant_id '_' Ss(1).session_id]; end
    if multi, Pd = fullfile(Preg,Ss(1).dataset); else, Pd = Preg; end
    for ti = 1:size(types,1)
      T = types{ti,3};
      for si = 1:numel(Ss)
        R = Ss(si).(types{ti,2}); 
        if ~isstruct(R) || isempty(fieldnames(R)), continue; end
        scan = regexprep(spm_file(Ss(si).filename,'basename'),'^sub-[^_]+_(ses-[^_]+_)?','');
        for fn = fieldnames(R)'
          v = R.(fn{1}); 
          switch types{ti,1}
            case 'anat', T(end+1,:) = {scan,fn{1},regionName(neuro,v,'ROIid'),field(v,'Vgm'),field(v,'Vwm'),field(v,'Vcsf')}; %#ok<AGROW>
            case 'dwi',  T(end+1,:) = {scan,fn{1},regionName(mori,fn{1},'ROIabbr'),field(v,'FA'),field(v,'MD'),field(v,'voxels')}; %#ok<AGROW>
            case 'func', T(end+1,:) = {scan,fn{1},v}; %#ok<AGROW>
          end
        end
      end
      if size(T,1) > 1
        cat_io_dcm2bids_helper('writeTable',fullfile(Pd,sprintf('%s_%s.%s',name,types{ti,1},fmt)),T,fmt);
        if ~exist(fullfile(Pd,[types{ti,1} '.json']),'file')
          cat_io_dcm2bids_helper('writeJSON',fullfile(Pd,[types{ti,1} '.json']), ...
            cell2struct(cellfun(@(d) struct('Description',d),types{ti,4},'UniformOutput',false),types{ti,3},2)); 
        end
      end
    end
  end
  if exist(Preg,'dir'), fprintf('  Region tables: %s\n',Preg); end
end
% =========================================================================
function A = atlasNames(name)
%atlasNames. ROIid, ROIabbr (if available), and ROIname of a CAT atlas.
  P = fullfile(fileparts(mfilename('fullpath')),'templates_MNI152NLin2009cAsym',[name '.csv']);
  C = regexp(fileread(P),'\r?\n','split'); C = C(~cellfun('isempty',strtrim(C)));
  C = cellfun(@(l) strsplit(l,';'),C,'UniformOutput',false); hdr = C{1}; C = vertcat(C{2:end});
  A.ROIid   = str2double(C(:,strcmp(hdr,'ROIid')));
  A.ROIname = strtrim(C(:,strcmp(hdr,'ROIname')));
  if any(strcmp(hdr,'ROIabbr')), A.ROIabbr = strtrim(C(:,strcmp(hdr,'ROIabbr'))); end
end
% =========================================================================
function n = regionName(A,v,by)
%regionName. Name of a region of the atlas names A by its ROIid (the field
%  ROIid of the record v) or its ROIabbr (v).
  n = ''; i = [];
  if strcmp(by,'ROIid') && isstruct(v) && isfield(v,'ROIid')
    i = find(A.ROIid == v.ROIid,1);
  elseif strcmp(by,'ROIabbr') && isfield(A,'ROIabbr')
    i = find(strcmp(A.ROIabbr,v),1);
  end
  if ~isempty(i), n = A.ROIname{i}; end
end
% =========================================================================
function x = field(v,fn)
%field. Numeric field fn of the record v (NaN if not available).
  if isstruct(v) && isfield(v,fn) && isnumeric(v.(fn)) && ~isempty(v.(fn)), x = double(v.(fn)); else, x = nan; end
end
% =========================================================================
function E = sessionRecords(S,types,tnam)
%sessionRecords. Records of the sessions (PI table) of the scan records S
%  with the combined ratings and the typical scan types.
  key = strcat({S.dataset},'|',{S.participant_id},'|',{S.session_id});
  [ukey,~,ki] = unique(key,'stable');
  E = struct([]);
  for ui = 1:numel(ukey)
    Sk = S(ki==ui);
    e  = struct('dataset',Sk(1).dataset,'participant_id',Sk(1).participant_id, ...
      'session_id',Sk(1).session_id,'site',Sk(1).site,'age',Sk(1).age,'sex',Sk(1).sex, ...
      'scans',numel(Sk));
    % combined ratings: worst protocol fit, RMS of the best quality of each
    % scan type, session orientation, worst coverage and reversal
    e.PFR = maxValue([Sk.PFR]);
    [ut,~,ti] = unique({Sk.type}); best = nan(1,numel(ut));
    for uj = 1:numel(ut), best(uj) = minValue([Sk(ti==uj).SQR]); end
    best = best(~isnan(best));
    if isempty(best), e.SQR = nan; else, e.SQR = sqrt(mean(best.^2)); end
    e.SOR = minValue([Sk.SOR]);
    e.BBR = maxValue([Sk.BBR]);
    e.PRR = maxValue([Sk.PRR]);
    for tj = 1:numel(types)
      St = Sk(strcmp({Sk.type},types{tj}));
      e.([tnam{tj} '_n'])   = numel(St);
      e.([tnam{tj} '_SQR']) = minValue([St.SQR]);
      e.([tnam{tj} '_PFR']) = maxValue([St.PFR]);
    end
    % tissue volumes of the best T1w scan with segmentation
    St = Sk(strcmp({Sk.type},'anat/T1w') & ~isnan([Sk.TIV]));
    if ~isempty(St)
      sqr = [St.SQR]; sqr(isnan(sqr)) = inf; [~,bi] = min(sqr); St = St(bi);
    end
    for fn = {'TIV','GMV','WMV','CSFV','rGMV','rWMV','rCSFV'}
      if isempty(St), e.(['T1w_' fn{1}]) = nan; else, e.(['T1w_' fn{1}]) = St.(fn{1}); end
    end
    if isempty(St), e.Pmwc1 = ''; else, e.Pmwc1 = St.Pmwc1; end
    % DTI and rs-fMRI measures of the best DWI and BOLD scan (with these 
    % measures) with the normalized FA map and the connectivity matrix
    e.PwFA = ''; e.Pconn = '';
    for tf = {{'dwi/dwi','dwi',{'FA_WM','MD_WM','MD_WM_peakwidth'},'PwFA'},{'func/bold','bold',{'FD_mean','tSNR','GCOR','segregation'},'Pconn'}}
      St = Sk(strcmp({Sk.type},tf{1}{1}) & ~isnan([Sk.(tf{1}{3}{1})]));
      if ~isempty(St)
        sqr = [St.SQR]; sqr(isnan(sqr)) = inf; [~,bi] = min(sqr); St = St(bi); e.(tf{1}{4}) = St.(tf{1}{4});
      end
      for fn = tf{1}{3}
        if isempty(St), e.([tf{1}{2} '_' fn{1}]) = nan; else, e.([tf{1}{2} '_' fn{1}]) = St.(fn{1}); end
      end
    end
    if isempty(E), E = e; else, E(end+1) = e; end %#ok<AGROW>
  end
end
% =========================================================================
function E = sessionStats(E)
%sessionStats. Statistics of the sessions (see cat_io_dcm2bids_stats): 
%  z-scores of the robust age/sex model of the sample of the TIV and the 
%  relative tissue volumes of the best T1w scan, of the white matter FA of
%  the best DWI scan, and of the system segregation of the best BOLD scan,
%  and the sample homogeneity (quartic mean Z-score, corrected for age and
%  sex) of the mwc1 map, of the normalized FA map (wFA), and of the 
%  connectivity matrix (Fisher z of the connections).
  fn = {'T1w_TIV','T1w_rGMV','T1w_rWMV','T1w_rCSFV','dwi_FA_WM','bold_segregation'};
  Y  = cell2mat(cellfun(@(f) reshape([E.(f)],[],1),fn,'UniformOutput',false));
  z  = cat_io_dcm2bids_stats('volumeModel',Y,[E.age],{E.sex});
  qmz = cat_io_dcm2bids_stats('homogeneity',{E.Pmwc1},[E.T1w_TIV],[E.age],{E.sex});
  qfa = cat_io_dcm2bids_stats('homogeneity',{E.PwFA},[],[E.age],{E.sex});
  qfc = cat_io_dcm2bids_stats('vectorHomogeneity',connVectors({E.Pconn}),[E.age],{E.sex});
  for ei = 1:numel(E)
    for fi = 1:numel(fn), E(ei).([fn{fi} '_z']) = z(ei,fi); end
    E(ei).T1w_QMZ = qmz(ei); E(ei).dwi_QMZ = qfa(ei); E(ei).bold_QMZ = qfc(ei);
  end
end
% =========================================================================
function X = connVectors(P)
%connVectors. Connections (upper triangle, Fisher z) of the connectivity
%  matrices P (cell of tables of cat_io_dcm2bids_pp_rsfmri, '' if not 
%  available) as rows of X (NaN for missing matrices or other sizes than 
%  the first one).
  X = []; m = 0;
  for i = 1:numel(P)
    x = [];
    if ~isempty(P{i}) && exist(P{i},'file')
      try
        C = readcell(P{i},'FileType','text','Delimiter','\t'); C = C(2:end,2:end);
        C(~cellfun(@(v) isnumeric(v) && isscalar(v),C)) = {nan}; R = cell2mat(C);
        if size(R,1) == size(R,2), x = R(triu(true(size(R)),1))'; end
      end
    end
    if m == 0 && ~isempty(x), m = numel(x); X = nan(numel(P),m); end
    if m > 0 && numel(x) == m, X(i,:) = x; end
  end
  if m == 0, X = nan(numel(P),1); end
end
% =========================================================================
function Q = projectRecords(S,E,project,types,tnam)
%projectRecords. Records of the project table: one row per dataset and
%  site and one for all sites of each dataset.
  Q = struct([]);
  [uds,~,di] = unique({E.dataset},'stable');
  for dj = 1:numel(uds)
    Ed = E(di==dj); Sd = S(strcmp({S.dataset},uds{dj}));
    sites = unique({Ed.site},'stable');
    if numel(sites) > 1, sites = [sites {'all'}]; end %#ok<AGROW>
    for sj = 1:numel(sites)
      if strcmp(sites{sj},'all'), Es = Ed; Ss = Sd;
      else, Es = Ed(strcmp({Ed.site},sites{sj})); Ss = Sd(strcmp({Sd.site},sites{sj}));
      end
      [usub,ui] = unique({Es.participant_id});
      q = struct('project',project,'dataset',uds{dj},'site',sites{sj}, ...
        'subjects',numel(usub),'sessions',numel(Es),'sessions_per_subject',numel(Es)/numel(usub), ...
        'scans',numel(Ss),'scans_per_session',numel(Ss)/numel(Es));
      for tj = 1:numel(types)
        q.([tnam{tj} '_per_subject']) = sum(strcmp({Ss.type},types{tj})) / numel(usub);
      end
      age = [Es(ui).age];
      sex = {Es(ui).sex}; sex = sex(~cellfun('isempty',sex)); % subjects with known sex
      if isempty(sex), q.males = nan; else, q.males = 100 * mean(strcmpi(sex,'M')); end
      q.age_mean = mean(age,'omitnan'); q.age_std = std(age,'omitnan');
      if sum(~isnan(age)) < 2, q.age_std = nan; end
      q.age_min  = minValue(age); q.age_max = maxValue(age);
      q.SQR = median([Es.SQR],'omitnan'); q.PFR = median([Es.PFR],'omitnan'); q.SOR = median([Es.SOR],'omitnan');
      st = {Ss.protocol_status}; st = st(~cellfun('isempty',st));
      if isempty(st), q.PCP = nan; else, q.PCP = 100 * mean(strcmp(st,'conform')); end
      if isempty(Q), Q = q; else, Q(end+1) = q; end %#ok<AGROW>
    end
  end
end
% =========================================================================
function [cols,desc] = scanColumns(type)
%scanColumns. Columns of the scan tables with their description (with the
%  tissue volumes of the segmentation for anat, DTI measures for dwi, and
%  rs-fMRI measures for func types, '' for the table of all scans).
  C = {
    'dataset'             'BIDS (protocol) directory of the scan'
    'participant_id'      'BIDS label of the subject'
    'session_id'          'BIDS label of the session'
    'site'                'site (abbreviation of the site dictionary)'
    'age'                 'age at the session in years'
    'sex'                 'sex of the participant (M, F, O)'
    'datatype'            'BIDS datatype (anat, func, dwi, fmap, perf)'
    'suffix'              'BIDS suffix'
    'entities'            'further BIDS entities besides sub, ses, task, and acq (e.g. dir-AP_part-phase)'
    'type'                'scan type (datatype/suffix with image part and inversion, as the render groups)'
    'filename'            'BIDS file of the scan (relative to the dataset, without extension)'
    'series'              'DICOM series number (otherwise of the acq entity, i.e. of the first series of the acquisition)'
    'protocol_dicom'      'DICOM protocol name'
    'protocol'            'closest protocol of the protocol test'
    'protocol_status'     'protocol status: conform, semiconform, or nonconform'
    'PFM'                 'protocol fit measure: weighted number of deviations from the closest protocol'
    'PFR'                 'protocol fit rating (marks 0.5-10.5, conform < 1.5, nonconform > 4.5)'
    'protocol_deviations' 'deviations from the closest protocol: field: scan value (protocol value)'
    'SQR'                 'overall quality rating of the scan (marks 0.5-10.5, RMS of the QC ratings)'
    'BSM'                 'rating of the between scan motion (marks)'
    'WSM'                 'rating of the within scan (slice) motion (marks)'
    'ISR'                 'rating of the inhomogeneity to signal ratio (marks)'
    'NSR'                 'rating of the noise to signal ratio (marks)'
    'RES'                 'rating of the resolution (marks)'
    'BSM_value'           'between scan motion (average correction of the realignment)'
    'WSM_value'           'within scan (slice) motion'
    'ISR_value'           'inhomogeneity to signal ratio'
    'NSR_value'           'noise to signal ratio'
    'RES_value'           'RMS of the voxel size in mm'
    'BBM'                 'missing brain coverage in percent of the expected coverage'
    'BBR'                 'rating of the brain coverage (marks)'
    'PRM'                 'correlation with the scan of opposite phase-encoding direction (high values: not reversed)'
    'PRR'                 'rating of the phase-encoding reversal (marks)'
    'OTM'                 'session orientation: RMSE of the origin offset to MNI in mm'
    'ORM'                 'session orientation: RMSE of the rotation to MNI in degree'
    'OTR'                 'rating of the origin offset of the session (marks)'
    'ORR'                 'rating of the rotation of the session (marks)'
    'SOR'                 'session orientation rating (marks)'
    'manufacturer'        'scanner manufacturer'
    'model'               'scanner model'
    'field_strength'      'magnetic field strength in Tesla'
    'software'            'scanner software version'
    'TR'                  'repetition time in seconds'
    'TE'                  'echo time in seconds'
    'TI'                  'inversion time in seconds'
    'FA'                  'flip angle in degree'
    'MB'                  'multiband acceleration factor'
    'dim_x'               'image matrix (x)'
    'dim_y'               'image matrix (y)'
    'dim_z'               'image matrix (z)'
    'volumes'             'number of volumes'
    'vox_x'               'voxel size (x) in mm'
    'vox_y'               'voxel size (y) in mm'
    'vox_z'               'voxel size (z) in mm'
    };
  if strncmp(type,'anat/',5)
    C = [C; {
      'TIV'    'total intracranial volume in ml (SPM segmentation)'
      'GMV'    'gray matter volume in ml'
      'WMV'    'white matter volume in ml'
      'CSFV'   'cerebrospinal fluid volume in ml'
      'rGMV'   'relative gray matter volume (GMV/TIV)'
      'rWMV'   'relative white matter volume (WMV/TIV)'
      'rCSFV'  'relative cerebrospinal fluid volume (CSFV/TIV)'
      'NCR'    'noise to contrast ratio of the segmentation'
      'TPMll'  'log-likelihood of the segmentation'
      }];
  end
  if strncmp(type,'dwi/',4)
    C = [C; {
      'FA_WM'           'median fractional anisotropy in the white matter (DTI, see cat_io_dcm2bids_pp_dwi)'
      'MD_WM'           'median mean diffusivity in the white matter in mm2/s'
      'AD_WM'           'median axial diffusivity in the white matter in mm2/s'
      'RD_WM'           'median radial diffusivity in the white matter in mm2/s'
      'MD_WM_peakwidth' 'peak width of the MD in the white matter (95th - 5th percentile, similar to PSMD) in mm2/s'
      'MD_brain'        'median mean diffusivity of the brain in mm2/s'
      'DWI_outliers'    'number of outlier volumes of the tensor fit (residual above 3 SD)'
      'pplevel'         'preprocessing level (see levelParameters in cat_io_dcm2bids_pp_anat)'
      }];
  end
  if strncmp(type,'func/',5)
    C = [C; {
      'FD_mean'     'mean framewise displacement in mm (Power et al. 2012, see cat_io_dcm2bids_pp_rsfmri)'
      'FD_max'      'maximum framewise displacement in mm'
      'FD_pct05'    'percentage of volumes with a framewise displacement above 0.5 mm'
      'DVARS'       'median signal change between volumes (RMS) in percent of the mean signal'
      'tSNR'        'temporal signal-to-noise ratio (median of the gray matter)'
      'GCOR'        'global correlation of the gray matter (Saad et al. 2013)'
      'FC_within'   'mean within-network connectivity (Fisher z, 7 networks of the Schaefer2018 100 parcels)'
      'FC_between'  'mean between-network connectivity (Fisher z)'
      'segregation' 'system segregation: (within - between) / within (Chan et al. 2014)'
      'pplevel'     'preprocessing level (see levelParameters in cat_io_dcm2bids_pp_anat)'
      }];
  end
  cols = C(:,1)'; desc = C(:,2)';
end
% =========================================================================
function [cols,desc] = sessionColumns(types,tnam)
%sessionColumns. Columns of the session table with their description.
  C = {
    'dataset'        'BIDS (protocol) directory'
    'participant_id' 'BIDS label of the subject'
    'session_id'     'BIDS label of the session'
    'site'           'site (abbreviation of the site dictionary)'
    'age'            'age at the session in years'
    'sex'            'sex of the participant (M, F, O)'
    'scans'          'number of scans of the session'
    'PFR'            'combined protocol fit rating: worst PFR of the scans (marks 0.5-10.5)'
    'SQR'            'combined quality rating: RMS of the best SQR of each scan type (marks)'
    'SOR'            'session orientation rating (marks)'
    'BBR'            'worst brain coverage rating of the scans (marks)'
    'PRR'            'worst phase-encoding reversal rating of the scans (marks)'
    };
  for ti = 1:numel(types)
    C = [C; {
      [tnam{ti} '_n']   sprintf('number of %s scans',types{ti})
      [tnam{ti} '_SQR'] sprintf('best quality rating SQR of the %s scans (marks)',types{ti})
      [tnam{ti} '_PFR'] sprintf('worst protocol fit rating PFR of the %s scans (marks)',types{ti})
      }]; %#ok<AGROW>
  end
  C = [C; {
    'T1w_TIV'   'total intracranial volume in ml (segmentation of the best T1w scan)'
    'T1w_GMV'   'gray matter volume in ml (best T1w scan)'
    'T1w_WMV'   'white matter volume in ml (best T1w scan)'
    'T1w_CSFV'  'cerebrospinal fluid volume in ml (best T1w scan)'
    'T1w_rGMV'  'relative gray matter volume (best T1w scan)'
    'T1w_rWMV'  'relative white matter volume (best T1w scan)'
    'T1w_rCSFV' 'relative cerebrospinal fluid volume (best T1w scan)'
    'T1w_TIV_z'   'z-score of the TIV of the robust age/sex model of the sample (see cat_io_dcm2bids_stats)'
    'T1w_rGMV_z'  'z-score of the relative gray matter volume of the robust age/sex model of the sample'
    'T1w_rWMV_z'  'z-score of the relative white matter volume of the robust age/sex model of the sample'
    'T1w_rCSFV_z' 'z-score of the relative cerebrospinal fluid volume of the robust age/sex model of the sample'
    'T1w_QMZ'     'sample homogeneity of the mwc1 map: quartic mean Z-score (high values: atypical map)'
    'dwi_FA_WM'           'median white matter FA of the best DWI scan (see cat_io_dcm2bids_pp_dwi)'
    'dwi_MD_WM'           'median white matter MD of the best DWI scan in mm2/s'
    'dwi_MD_WM_peakwidth' 'peak width of the white matter MD of the best DWI scan in mm2/s'
    'dwi_FA_WM_z'         'z-score of the white matter FA of the robust age/sex model of the sample'
    'dwi_QMZ'             'sample homogeneity of the normalized FA map (wFA): quartic mean Z-score (high values: atypical map)'
    'bold_FD_mean'        'mean framewise displacement of the best BOLD scan in mm (see cat_io_dcm2bids_pp_rsfmri)'
    'bold_tSNR'           'temporal signal-to-noise ratio of the best BOLD scan'
    'bold_GCOR'           'global correlation of the best BOLD scan'
    'bold_segregation'    'system segregation of the best BOLD scan'
    'bold_segregation_z'  'z-score of the system segregation of the robust age/sex model of the sample'
    'bold_QMZ'            'sample homogeneity of the connectivity matrix: quartic mean Z-score of the connections (high values: atypical connectivity)'
    }];
  cols = C(:,1)'; desc = C(:,2)';
end
% =========================================================================
function [cols,desc] = projectColumns(types,tnam)
%projectColumns. Columns of the project table with their description.
  C = {
    'project'              'project (output subdirectory)'
    'dataset'              'BIDS (protocol) directory'
    'site'                 'site (abbreviation of the site dictionary, all - all sites of the dataset)'
    'subjects'             'number of subjects'
    'sessions'             'number of sessions'
    'sessions_per_subject' 'average number of sessions per subject'
    'scans'                'number of scans'
    'scans_per_session'    'average number of scans per session'
    };
  for ti = 1:numel(types)
    C = [C; {[tnam{ti} '_per_subject'] sprintf('average number of %s scans per subject',types{ti})}]; %#ok<AGROW>
  end
  C = [C; {
    'males'    'percentage of male subjects (of the subjects with known sex)'
    'age_mean' 'mean age in years (first session of each subject)'
    'age_std'  'standard deviation of the age in years'
    'age_min'  'minimum age in years'
    'age_max'  'maximum age in years'
    'SQR'      'median combined quality rating of the sessions (marks)'
    'PFR'      'median combined protocol fit rating of the sessions (marks)'
    'SOR'      'median session orientation rating (marks)'
    'PCP'      'protocol conformity percentage (PCP): percentage of protocol conform scans'
    'color'    ['color of the site (RGB hex) of the site dictionary (Color of the scanner or of its city, scanners ' ...
                'with the same color differ in brightness, n/a if not defined), see siteColors of cat_io_dcm2bids_protocols']
    }];
  cols = C(:,1)'; desc = C(:,2)';
end
% =========================================================================
function writeReport(Ptable,R,cols,desc,fmt)
%writeReport. Write the records R with the columns cols as table and the
%  description of the columns as JSON file (BIDS style, with the rating
%  scale of the marks).
  T = cols;
  for ri = 1:numel(R)
    T(end+1,:) = cellfun(@(c) R(ri).(c),cols,'UniformOutput',false); %#ok<AGROW>
  end
  cat_io_dcm2bids_helper('writeTable',Ptable,T,fmt);
  D = struct();
  for ci = 1:numel(cols)
    D.(cols{ci}) = struct('Description',desc{ci});
  end
  D.RatingScale = ['Ratings are school marks from 0.5 (best) to 10.5 (worst), i.e. 1 very good, ' ...
    '2 good, 3 satisfactory, 4 sufficient, 5 critical, 6 unacceptable. Missing values are n/a.'];
  cat_io_dcm2bids_helper('writeJSON',spm_file(Ptable,'ext','.json'),D);
end
% =========================================================================
function [types,names] = scanTypes(opts)
%scanTypes. Typical scan types (datatype/suffix, e.g. anat/T1w, see
%  opts.scantypes) and their column names (the suffix, or the datatype and
%  suffix if the suffix is not unique, as in cat_io_dcm2bids_tables).
  types = cellstr(opts.scantypes);
  names = regexprep(types,'^[^/]*/','');
  dup   = cellfun(@(n) sum(strcmp(names,n)) > 1, names);
  names(dup) = strrep(types(dup),'/','_');
  names = regexprep(names,'[^a-zA-Z0-9_]','_');
end
% =========================================================================
function s = sexCode(s)
%sexCode. Sex as M, F, or O (e.g. of 'male' or 'f' of imported BIDS 
%  datasets), otherwise empty.
  switch lower(strtrim(s))
    case {'m','male'},   s = 'M';
    case {'f','female'}, s = 'F';
    case {'o','other'},  s = 'O';
    otherwise,           s = '';
  end
end
% =========================================================================
function T = readTSVcached(cache,P)
%readTSVcached. BIDS TSV file (header and rows as text) with cache.
  if isKey(cache,P), T = cache(P); return; end
  T = {};
  if exist(P,'file')
    txt = strtrim(cat_io_dcm2bids_helper('readText',P));
    L   = regexp(txt,'\r?\n','split');
    T   = cellfun(@(l) strsplit(l,sprintf('\t'),'CollapseDelimiters',false),L,'UniformOutput',false);
    n   = numel(T{1});
    T   = T(cellfun(@numel,T) == n);
    T   = vertcat(T{:});
  end
  cache(P) = T; %#ok<NASGU> handle object
end
% =========================================================================
function v = tsvValue(T,keycol,key,col)
%tsvValue. Value of the column col in the row with key in keycol ('' if not
%  available or n/a).
  v = '';
  if isempty(T), return; end
  kc = find(strcmp(T(1,:),keycol),1); vc = find(strcmp(T(1,:),col),1);
  if isempty(kc) || isempty(vc), return; end
  ri = find(strcmp(T(2:end,kc),key),1);
  if ~isempty(ri) && ~strcmpi(T{ri+1,vc},'n/a'), v = T{ri+1,vc}; end
end
% =========================================================================
function x = tsvNumber(T,keycol,key,col)
%tsvNumber. Numeric value of tsvValue (NaN if not available).
  x = str2double(tsvValue(T,keycol,key,col));
end
% =========================================================================
function J = readJSON(P)
%readJSON. JSON file as structure (empty structure if not available).
  J = struct();
  if exist(P,'file')
    try
      J = cat_io_json(P);
    catch
      J = struct(); % e.g. incomplete files
    end
  end
end
% =========================================================================
function J = readJSONcached(cache,P)
%readJSONcached. readJSON with cache (e.g. for the session registration).
  if isKey(cache,P), J = cache(P); return; end
  J = readJSON(P);
  cache(P) = J; %#ok<NASGU> handle object
end
% =========================================================================
function v = jsonPath(J,path,default)
%jsonPath. Value of the field path (e.g. {'protocol','status'}) of the
%  structure J (default if not available).
  v = J;
  for pi = 1:numel(path)
    if isstruct(v) && isfield(v,path{pi}), v = v.(path{pi}); else, v = default; return; end
  end
  if iscell(v) && ~iscellstr(v), v = default; end %#ok<ISCLSTR>
end
% =========================================================================
function x = jsonNumber(J,path)
%jsonNumber. Numeric scalar of the field path of J (NaN if not available).
  x = jsonPath(J,path,nan);
  if ~(isnumeric(x) || islogical(x)) || isempty(x), x = nan; else, x = double(x(1)); end
end
% =========================================================================
function x = minValue(v)
%minValue. Minimum of the finite values (NaN if none).
  v = v(isfinite(v)); if isempty(v), x = nan; else, x = min(v); end
end
% =========================================================================
function x = maxValue(v)
%maxValue. Maximum of the finite values (NaN if none).
  v = v(isfinite(v)); if isempty(v), x = nan; else, x = max(v); end
end

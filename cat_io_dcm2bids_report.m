function varargout = cat_io_dcm2bids_report(action,varargin)
%cat_io_dcm2bids_report. Report tables of cat_io_dcm2bids.
%  The report table summarizes each BIDS (protocol) directory with the 
%  number of subjects, sessions, and scans and the sex and age of the 
%  participants (writeReportTable). The subject and scan report tables are
%  not implemented yet (writeSubReportTSV, writeScanReportTSV).
%
%  Planned reports: 
%   basic fields: 
%    - sub, age-range, age(mn+sd), #MRI (session/subject)
%   quality-ratings (5): 
%    - anat/func/dwi/fmap/total
%   subject-measures: 
%    - anat:     TIV, rGMV, rWMV, rCSFV, GMT,
%    - func-rs:  network connectivity
%    - dwi:      AD,FD, ...
%
%  varargout = cat_io_dcm2bids_report(action,varargin)
%
%  Actions (see the help of the local functions):
%    cat_io_dcm2bids_report('writeReportTable',Poutdir,Pdirs,job)
%    cat_io_dcm2bids_report('writeSubReportTSV',sub,Poutdir,BIDSsubdir,subdir)
%    cat_io_dcm2bids_report('writeScanReportTSV',Vjson,sub,site,Poutdir,BIDSsubdir,subdir)
%
%  See also cat_io_dcm2bids.

  switch action
    case {'writeReportTable','writeSubReportTSV','writeScanReportTSV'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_report:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function writeReportTable(Poutdir,Pdirs,job)
%writeReportTable. Report table with one row per BIDS (protocol) directory 
%  of Pdirs with the number of subjects, sessions per subject, and scans per
%  subject (anat, dwi, func) and the sex and age of the participants.tsv: 
%    <outdir>/<subdir>/<BIDSdir>-report/report_<BIDSsubdir>.<tableformat>
%  The table is created again in each run. 
%  
%  Planned: 
%  study||protocol||#subjects|%male|minage|maxage||#anat/sub|#dwi/sub|#func/sub||aQR|dQR|fQR|SQR||Vgm|Vwm|Vcsf|TIV|mnFA|mn...
  Psubjecttable = fullfile(Poutdir,[job.BIDSdir '-report'],sprintf('report_%s.%s', job.BIDSsubdir, job.opts.tableformat));
  if exist(Psubjecttable,'file'), delete(Psubjecttable); end
  for pdi = 1:numel(Pdirs)
    Protocol = spm_file(Pdirs{pdi},'basename'); 
    if job.opts.protocolsubdirs
      Pparticipants = fullfile(Poutdir,job.BIDSdir,job.BIDSsubdir,Protocol,'participants.tsv');
    else
      Pparticipants = fullfile(Poutdir,job.BIDSdir,job.BIDSsubdir,'participants.tsv');
    end
    Thdr = {'project', 'protocol','#subjects', '#sessions/subject', ...
      '#anat/subjects', '#dwi/subjects', '#func/subjects', ...
      '%males', 'mean(age)', 'std(age)', 'min(age)', 'max(age)' }; 
    if ~exist(Pparticipants,'file')
      Tparticipants = Thdr;
    else
      Tparticipants = cat_io_csv(Pparticipants,'','',struct('delimiter','\t','convert2double',-1));
    end
    
    % look for available data
    if ~exist(Pdirs{pdi},'dir'), continue; end
    subjects = cat_vol_findfiles( Pdirs{pdi}, 'sub-*', struct('depth',1,'dirs',1));
    sessions = cat_vol_findfiles( Pdirs{pdi}, 'ses-*', struct('depth',2,'dirs',1));
    anat     = cat_vol_findfiles( Pdirs{pdi}, 'anat',  struct('depth',3,'dirs',1));
    func     = cat_vol_findfiles( Pdirs{pdi}, 'func',  struct('depth',3,'dirs',1));
    dwi      = cat_vol_findfiles( Pdirs{pdi}, 'dwi',   struct('depth',3,'dirs',1));

    % add row (only with participants, e.g. not if no scan was exported)
    if size(Tparticipants,1)>1
      Tnewrow = {job.subdir, Protocol, numel(subjects), numel(sessions)/numel(subjects), ...
        numel(anat)/numel(subjects), numel(dwi)/numel(subjects), numel(func)/numel(subjects), ...
        mean(cellfun(@(x) strcmp(x,'M'), Tparticipants(2:end,2))), mean(cell2mat(Tparticipants(2:end,3))), ...
        std(cell2mat(Tparticipants(2:end,3))), min(cell2mat(Tparticipants(2:end,3))), max(cell2mat(Tparticipants(2:end,3))), ...
        }; 
      cat_io_dcm2bids_helper('updateTable',Psubjecttable,Thdr,Tnewrow,pdi,0);
    end
  % create/extend BIDS csv files ???

    
  end
end
% =========================================================================
function writeSubReportTSV(sub,Poutdir,BIDSsubdir,subdir)
% Create report file with one subjects per row and simplified scan data
% this one could run on the scanreport-tsv files

return

  Preport = fullfile(Poutdir, subdir, ...
    sprintf('subjectreport_site-%s_study-%s_protocoldir-%s.tsv', Vjson.center, BIDSsubdir));
  if exist(Preport,'file')
    Treport = cat_io_csv(Preport, '','', struct('delimiter','\t','convert2double',0)); 
    pidpa = find(matches( Treport(2:end,1) , sub )) + 1;
    if isempty(pidpa) || pidpa<=0, pidpa = size(Treport,1) + 1; end
  else
    Treport = {
      'participant_id','participant_sex','participant_age',... 
      'center_id','study_id','study_date','study_timepoints', ...
      ... scanner_model_field-strength , ...
      'mr_anat_t1w','mr_func','mr_dwi','mr_fmap',''}; 
      pidpa = 2;
  end
  Treport(pidpa,:) = {
    sub, Vjson.PatientSex, round(Vjson.PatientAge), ...
    Vjson.center, Vjson.StudyID, Vjson.StudyID, ...     
    '',BIDSsubdir,};
  Treport(2:end,:) = sortrows(Treport(2:end,:)); 
  
  cat_io_csv(Preport,Treport,'','',struct('delimiter','\t')); 

  % subreports
end
% =========================================================================
function writeScanReportTSV(Vjson,sub,site,Poutdir,BIDSsubdir,subdir)
% Create report file with scans per row for further evaluation!
% - for all files
% - for each protocol set (extraction of the main report)
return

  Preport = fullfile(Poutdir, subdir, ...
    sprintf('scanreport_site-%s_study-%s_protocoldir-%s.tsv', site, BIDSsubdir));
  if exist(Preport,'file')
    Treport = cat_io_csv(Preport, '','', struct('delimiter','\t','convert2double',0)); 
    pidpa = find(matches( Treport(2:end,1) , sub )) + 1;
    if isempty(pidpa) || pidpa<=0, pidpa = size(Treport,1) + 1; end
    
  else
    Treport = {
      'participant_id','participant_sex','participant_age',... 
      ... 'center_id','study_id','study_date', ...
      ... scanner_model_field-strength , ...
      'mr_name','mr_para'}; 
    pidpa = 2;
  end
  Treport(pidpa,:) = {
    sub, Vjson.PatientSex, round(Vjson.PatientAge), ...
    ... site, Vjson.StudyID, Vjson.StudyID, ...     
    Vjson.ProtocolName, BIDSsubdir};
  Treport(2:end,:) = sortrows(Treport(2:end,:)); 
  
  cat_io_csv(Preport,Treport,'','',struct('delimiter','\t')); 

  % subreports

end

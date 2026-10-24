function varargout = cat_io_dcm2bids_export(action,varargin)
%cat_io_dcm2bids_export. Anonymized DICOM export of cat_io_dcm2bids.
%  The DICOM files of the scans of a run are exported with cleaned headers
%  and defaced images, e.g. to transfer the data to a central site that
%  runs the batch again (opts.dicomexport: 1 - scans with conform or
%  semiconform protocol, i.e. all scans without protocols, 2 - all scans).
%  The export uses the internal database, i.e. the scans of the input
%  selection that were imported (no localizers, scouts, derived or 2D
%  images, and no BIDS inputs, as they have no DICOM files).
%
%  THE EXPORT IS NOT AN ARCHIVE OR BACKUP OF THE ORIGINAL DATA, AS THE
%  HEADERS ARE CLEANED AND THE IMAGES ARE DEFACED. KEEP THE ORIGINAL DATA!
%
%  The export directory catDCM2BIDSdb-dicom next to the database contains
%    <yyyyMMdd-HHmmss>/  one directory per run with new scans and for each
%                        session the zip file sub-*_ses-*.zip with the
%                        DICOM files in the layout of the database, i.e.
%                        sub-*/ses-*/snr-*/<series>_<instance>.dcm (the
%                        original names are not used, as they can contain
%                        names), and its content table sub-*_ses-*.<fmt>
%                        with the images of the series (the converted 
%                        NIfTIs as subpath of the database without the
%                        DICOM directory of the file name, for the 
%                        anonymization level 2 the BIDS files), the
%                        directory of the series in the zip file, the 
%                        number of DICOM files, MB, and the acquisition
%                        duration (s, see acquisitionDuration in 
%                        cat_io_dcm2bids_protocols)
%    scans.<fmt>         exported scans of all runs (without personal data)
%    sessions.<fmt>      exported zip files of all runs
%    exports.mat         lists of the tables: a scan (DICOM series) is
%                        exported only once for an anonymization level,
%                        i.e. the run directories can be removed after the
%                        upload (remove a scan from the list to export it
%                        again)
%    readme.txt
%  where <fmt> is opts.tableformat.
%
%  Header: the changes are made in the byte stream of the files (explicit
%  VR little endian), i.e. all other elements stay unchanged, including the
%  private elements (e.g. the Siemens CSA headers used by dcm2niix).
%  - cleaned: name of the patient (replaced by the PatientID), names of
%    physicians and operators, addresses, phone numbers, other patient IDs,
%    accession and procedure IDs, and comments (see cleanedElements)
%  - birth date: the 15th of the month for opts.anonymize > 0 (for level 2
%    of the shifted date, i.e. the age with an accuracy of about a month; 
%    the age estimate of the import uses the 15th, i.e. the original and 
%    exported data give the same age, see assurePatientDCMfields in 
%    cat_io_dcm2bids_db)
%  - anonymization level 2 (opts.anonymize = 2, see anonymizeFile and
%    cat_io_dcm2bids_anonymize): the BIDS subject label as PatientID and
%    PatientName and the session label as StudyID, the dates of a session 
%    shifted to the date of its label and its times so that it starts at 
%    the center of its 3-hour interval of the day (01:30, 04:30, ..., 22:30,
%    i.e. the timing of scans and slices is kept, also dates and times in 
%    private elements, the birth date with the offset of the subject's 
%    first session), UIDs replaced by UIDs of the same length
%    (also in private elements), the scanner and institution by the site 
%    label (siteID, manufacturer, and field strength, without location, see
%    setupSites in cat_io_dcm2bids_protocols), no address, station, model,
%    coils, department, and study and procedure descriptions (see 
%    anonymizedElements; the SoftwareVersions are kept for the conversion),
%    parts of the patient name, the PatientID, the serial number, and the 
%    institution, city, model, and coil names in further elements (e.g. the
%    Siemens CSA headers) replaced by X, and the names
%    of the zip files and directories with the BIDS labels of the subject 
%    and session
%  - kept (level 1): PatientID (that has to be a pseudonym), sex, age, 
%    size, weight, dates and times, UIDs, scanner, institution, protocol and
%    sequence fields, image comments (e.g. the multiband factor), private 
%    elements
%  - check: the parts of the patient name are searched in the header of
%    each file (warning with the elements, e.g. a name in a description)
%  Image: the face is removed by the plane of spm_deface (registration of
%  the image of the series in the database or of its mean image for 4D
%  scans) as by the anonymization of the BIDS export, i.e. for anatomical
%  scans (opts.anonymize = 1) or all scans (opts.anonymize > 1). Pixels in
%  front of the plane are set to 0 (single-frame, multi-frame, and Siemens
%  mosaic images). Scans that cannot be defaced are not exported.
%  Limitations: only files with explicit VR little endian are exported (no
%  compressed or implicit transfer syntaxes), and the mosaic images of UIH
%  (e.g. DTI saved by slice, with the slice positions in a sequence instead
%  of the image position) are not supported by the defacing (not exported 
%  if they have to be defaced, e.g. all scans for opts.anonymize = 2; 
%  tested with the dcm_qa_uih data of github.com/neurolabusc, to be added 
%  with own UIH data).
%
%  varargout = cat_io_dcm2bids_export(action,varargin)
%
%  Actions (see the help of the local functions):
%    cat_io_dcm2bids_export('exportDICOM',Pjson,datatype,status,site,sub,ses,BIDSfile,Pmdbdir,opts)
%    [b,info] = cat_io_dcm2bids_export('anonymizeFile',b,plane,opts,X)
%    tags = cat_io_dcm2bids_export('cleanedElements')
%    E = cat_io_dcm2bids_export('dicomElements',b)
%
%  See also cat_io_dcm2bids, cat_io_dcm2bids_db, spm_deface.

  switch action
    case {'exportDICOM','anonymizeFile','cleanedElements','dicomElements'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_export:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function exportDICOM(Pjson,datatype,status,site,sub,ses,BIDSfile,Pmdbdir,opts)
%exportDICOM. Export of the DICOM series of the database scans Pjson
%  (sidecars) with their datatype, protocol status (see testProtocols in
%  cat_io_dcm2bids_protocols), site, subject, and session (BIDS labels), 
%  and BIDS file into the directory <Pmdbdir>-dicom, i.e. of the series that were not 
%  exported before (for this anonymization level). The lists of the tables
%  are saved after each session. The anonymization level 2 uses the BIDS 
%  labels and the date offsets of the registry of the database (see 
%  cat_io_dcm2bids_anonymize).
  Pexp = [Pmdbdir '-dicom'];
  Plst = fullfile(Pexp,'exports.mat');
  R = struct([]); Z = struct([]);
  if exist(Plst,'file'), load(Plst,'R','Z'); end % (records of older versions with other fields, see catRecords)
  A = [];
  if opts.anonymize > 1, A = cat_io_dcm2bids_anonymize('registry',Pmdbdir); end

  fprintf('\nDCM2BIDS - anonymized DICOM export\n');
  U = selectSeries(Pjson,datatype,status,site,sub,ses,BIDSfile,opts,A);
  if ~isempty(R) && ~isempty(U)
    key = @(uid,an) strcat(uid,'|',arrayfun(@(a) sprintf('%d',a),an,'UniformOutput',false));
    U   = U(~ismember(key({U.key},repmat(opts.anonymize,1,numel(U))),key({R.SeriesInstanceUID},[R.Anonymize])));
  end
  if isempty(U), fprintf('  No new scans.\n'); return; end

  run    = char(datetime('now','Format','yyyyMMdd-HHmmss'));
  Prun   = fullfile(Pexp,run);
  sesdir = cellfun(@fileparts,{U.snrdir},'UniformOutput',false);
  sess   = unique(sesdir,'stable');
  fprintf('  %d new scans of %d session(s) to %s\n',numel(U),numel(sess),Prun);
  tpm = []; % tissue probability map of spm_deface (loaded if required)
  for si = 1:numel(sess)
    Us = U(strcmp(sesdir,sess{si}));
    for ui = find([Us.deface] & cellfun('isempty',{Us.msg}))
      try
        if isempty(tpm), tpm = spm_load_priors8(spm_vol(fullfile(spm('Dir'),'tpm','TPM.nii'))); end
        Us(ui).plane = facePlane(Us(ui).Pjson,Pmdbdir,opts,tpm);
      catch e
        Us(ui).msg = sprintf('not exported, defacing failed (%s)',e.message);
      end
    end
    [Rs,Zs] = exportSession(Us,Prun,run,Pmdbdir,opts);
    R = cat_io_dcm2bids_helper('catRecords',R,Rs); Z = cat_io_dcm2bids_helper('catRecords',Z,Zs);
    if ~exist(Pexp,'dir'), mkdir(Pexp); end
    save(Plst,'R','Z');
    writeExportTables(Pexp,R,Z,opts);
  end
  if exist(Prun,'dir') && isempty(dir(fullfile(Prun,'*.zip'))), rmdir(Prun,'s'); end
end
% =========================================================================
function U = selectSeries(Pjson,datatype,status,site,sub,ses,BIDSfile,opts,A)
%selectSeries. Export units of the scans Pjson (database sidecars), i.e.
%  one unit per DICOM series of the scans with conform or semiconform
%  protocol (all scans without protocols) for opts.dicomexport = 1 or of
%  all scans for opts.dicomexport = 2, without BIDS inputs (no DICOM
%  files). A series is defaced if one of its scans is defaced by the
%  anonymization (anatomical scans or all scans for opts.anonymize > 1,
%  see anonymize in cat_io_dcm2bids_pp_anat). Units that cannot be
%  exported get a message. The scans of a series are listed with their 
%  sidecar (Pjson), BIDS file (bids), and acquisition duration (see 
%  acquisitionDuration in cat_io_dcm2bids_protocols, with the image 
%  dimensions). The names of the subject and session are the
%  names of the database, or for the anonymization level 2 the BIDS labels
%  sub and ses (with the anonymization X of the files, i.e. the labels and
%  the date offsets of the session and subject, and the hashed 
%  SeriesInstanceUID as key).
  U = struct('uid',{},'key',{},'Pjson',{},'bids',{},'duration',{},'scan',{},'snrdir',{},'subname',{},'sesname',{},'X',{}, ...
    'dcmdir',{},'deface',{},'datatype',{},'status',{},'series',{},'protocol',{},'plane',{},'msg',{});
  for si = 1:numel(Pjson)
    if isempty(Pjson{si}) || (opts.dicomexport == 1 && ~any(strcmp(status{si},{'conform','semiconform',''})))
      continue
    end
    Pimp = spm_file(Pjson{si},'prefix','catDCM2BIDSimport_');
    if ~exist(Pimp,'file'), continue; end
    I = cat_io_dcm2bids_db('readImportJSON',Pimp);
    if ~isfield(I,'DICOMdirectory') || isfield(I,'SourceFile'), continue; end % BIDS input
    J = cat_io_json(Pjson{si});
    if ~isfield(J,'SeriesInstanceUID') || isempty(J.SeriesInstanceUID), continue; end
    bids = ''; if numel(BIDSfile) >= si && ~isempty(BIDSfile{si}), bids = spm_file(BIDSfile{si},'ext','.nii.gz'); end
    Pn = spm_file(Pjson{si},'ext','.nii.gz'); if ~exist(Pn,'file'), Pn = spm_file(Pjson{si},'ext','.nii'); end
    if exist(Pn,'file'), [~,J.Dimensions] = cat_io_dcm2bids_helper('niiHeaderInfo',Pn); end
    dur = cat_io_dcm2bids_protocols('acquisitionDuration',J);
    deface = opts.anonymize > 1 || (opts.anonymize == 1 && strcmp(datatype{si},'anat'));
    ui = find(strcmp({U.uid},J.SeriesInstanceUID),1);
    if isempty(ui)
      snrdir = fileparts(Pjson{si}); [~,scan] = fileparts(snrdir);
      fnameparts = cat_io_dcm2bids_db('nameParts',spm_file(Pjson{si},'basename'));
      msg = '';
      if ~isfield(J,'PatientID') || isempty(J.PatientID)
        msg = 'not exported, no PatientID (the database subject would be named by the patient name)';
      elseif ~exist(I.DICOMdirectory,'dir')
        msg = 'not exported, DICOM directory not available';
      end
      [Psub,sesname] = fileparts(fileparts(snrdir)); [~,subname] = fileparts(Psub);
      key = J.SeriesInstanceUID; X = [];
      if ~isempty(A) % anonymization level 2
        [~,boffset] = cat_io_dcm2bids_anonymize('subject',A,subname);
        t = regexp(fnameparts{min(3,end)},'^(\d{8})(\d{6})?','tokens','once'); % study date and time
        offset = boffset; date = ''; time = 0; start = 0; % (start of the session at 00:00:00 without study time)
        if ~isempty(t) && ~isempty(t{2})
          time  = [3600 60 1] * str2double({t{2}(1:2); t{2}(3:4); t{2}(5:6)}); 
          start = cat_io_dcm2bids_anonymize('dayInterval',time); % (center of the 3-hour interval of the day)
        end
        try % (without valid study date only the shift of the subject)
          [~,offset] = cat_io_dcm2bids_anonymize('session',t{1},boffset,opts.sesunit); date = t{1};
        catch
          offset = boffset; date = '';
        end
        addr = strtrim(strsplit(fieldText(J,'InstitutionAddress'),',')); % (location, model, and coils replaced by X)
        wrd  = {fieldText(J,'InstitutionName'),addr{min(2,end)},fieldText(J,'ManufacturersModelName'), ...
          fieldText(J,'ReceiveCoilName'),fieldText(J,'CoilString')};
        X    = struct('code',sub{si}(5:end),'ses',regexprep(ses{si},'^ses-',''),'offset',offset, ...
          'boffset',boffset,'date',date,'time',time,'start',start,'site',site{si},'words',{wrd(cellfun(@numel,wrd) >= 4)},'salt',A.salt);
        subname = sub{si}; sesname = ses{si};
        key  = cat_io_dcm2bids_anonymize('hashUID',key,A.salt);
      end
      U(end+1) = struct('uid',J.SeriesInstanceUID,'key',key,'Pjson',{Pjson(si)},'bids',{{bids}},'duration',dur, ...
        'scan',scan,'snrdir',snrdir, ...
        'subname',subname,'sesname',sesname,'X',X,'dcmdir',I.DICOMdirectory,'deface',deface, ...
        'datatype',datatype{si},'status',status{si},'series',J.SeriesNumber,'protocol',fnameparts{min(2,end)}, ...
        'plane',[],'msg',msg); %#ok<AGROW>
    else
      U(ui).Pjson{end+1} = Pjson{si};
      U(ui).bids{end+1}  = bids;
      U(ui).duration(end+1) = dur;
      U(ui).deface       = U(ui).deface || deface;
    end
  end
end
% =========================================================================
function w = facePlane(Pjson,Pmdbdir,opts,tpm)
%facePlane. Face plane w of spm_deface in world coordinates (RAS, mm) of
%  the first magnitude image of a series in the database (its mean image
%  for 4D scans, as the anonymization of the BIDS export, see anonymize in
%  cat_io_dcm2bids_pp_anat), i.e. w*[x y z 1]' < 0 in front of the face.
%  Missing images are converted (e.g. for the JSON-only outputs).
  Pj = Pjson{1};
  for ji = 1:numel(Pjson)
    J = cat_io_json(Pjson{ji});
    if ~isfield(J,'ImageType') || ~any(ismember(upper(cellstr(J.ImageType)),{'P','PHASE'}))
      Pj = Pjson{ji}; break
    end
  end
  if opts.gzipi, niiext = '.nii.gz'; else, niiext = '.nii'; end
  Pn = spm_file(Pj,'ext',niiext);
  if ~exist(Pn,'file')
    [~,tn] = fileparts(tempname); Ptmp = fullfile(Pmdbdir,['catDCM2BIDStmp_' tn]); %#ok<NASGU> (evalc)
    c1 = onCleanup(@() removeDir(fullfile(Pmdbdir,['catDCM2BIDStmp_' tn])));
    evalc('cat_io_dcm2bids_importDCM(''convertImages'',{fileparts(Pj)},Pmdbdir,Ptmp,opts);');
    if ~exist(Pn,'file'), error('cat_io_dcm2bids_export:noImage','no image'); end
  end
  T  = tempname; mkdir(T);
  c2 = onCleanup(@() removeDir(T));
  if endsWith(Pn,'.gz'), Pu = gunzip(Pn,T); Pu = Pu{1}; else, Pu = Pn; end
  V = spm_vol(Pu);
  if numel(V) > 1
    Y = 0; for vi = 1:numel(V), Y = Y + spm_read_vols(V(vi)); end
    Pu = fullfile(T,'mean.nii');
    spm_write_vol(struct('fname',Pu,'dim',V(1).dim,'mat',V(1).mat,'dt',[16 0],'pinfo',[1 0 0]'),Y / numel(V));
  end
  M = spm_maff8(Pu,4,20,tpm,[],'mni');
  w = [0 -1.1 0.98 100] * M; % plane of spm_deface (in MNI space)
end
% =========================================================================
function [Rs,Zs] = exportSession(Us,Prun,run,Pmdbdir,opts)
%exportSession. Zip file of the series Us of a session with the anonymized
%  DICOM files (see exportFiles) and its content table (one row per image
%  of the database Pmdbdir, see the help of this file), and the records of
%  the export tables (Rs for the scans and Zs for the zip file).
  stime = tic; Rs = struct([]); Zs = struct([]);
  sub    = Us(1).subname; ses = Us(1).sesname;
  name   = [sub '_' ses];
  Pstage = fullfile(Prun,['.' name]); % files of the zip (removed at the end, also after errors)
  c = onCleanup(@() removeDir(Pstage));
  Us = exportFiles(Us,fullfile(Pstage,sub,ses),opts);
  ok = cellfun('isempty',{Us.msg});
  if any(ok)
    Pzip = fullfile(Prun,[name '.zip']);
    zipDirectory(Pstage,sub,Pzip);
    C = {'Image','Series','DICOMfiles','MB','AcquisitionDuration'};
    for ui = find(ok)
      for k = 1:numel(Us(ui).Pjson)
        if opts.anonymize > 1, img = Us(ui).bids{k}; else, img = dbImage(Us(ui).Pjson{k},Pmdbdir); end
        if isempty(img), img = 'n/a'; end % (no BIDS file, e.g. output 0)
        C(end+1,:) = {img, strjoin({sub,ses,Us(ui).scan},'/'), numel(Us(ui).files), ...
          round(Us(ui).bytes / 1e6,1), Us(ui).duration(k)}; %#ok<AGROW>
      end
    end
    cat_io_dcm2bids_helper('writeTable',fullfile(Prun,[name '.' opts.tableformat]),C,opts.tableformat);
    Fz = dir(Pzip); nf = sum(cellfun(@numel,{Us(ok).files})); mb = sum([Us(ok).bytes]) / 1e6;
    for ui = find(ok)
      Rs = [Rs struct('Export',run,'Zip',[name '.zip'],'Subject',sub,'Session',ses,'Scan',Us(ui).scan, ...
        'SeriesNumber',Us(ui).series,'ProtocolName',Us(ui).protocol,'Datatype',Us(ui).datatype, ...
        'Status',Us(ui).status,'Files',numel(Us(ui).files),'MB',round(Us(ui).bytes / 1e6,1), ...
        'Defaced',double(~isempty(Us(ui).plane)),'Anonymize',opts.anonymize, ...
        'AcquisitionDuration',Us(ui).duration(1), ...
        'NameCheck',strjoin(Us(ui).names,' '),'SeriesInstanceUID',Us(ui).key)]; %#ok<AGROW>
    end
    Zs = struct('Export',run,'Zip',[name '.zip'],'Subject',sub,'Session',ses,'Scans',nnz(ok), ...
      'Files',nf,'MB',round(mb,1),'ZipMB',round(Fz.bytes / 1e6,1),'Duration',round(toc(stime)));
    fprintf('  %s/%s: %d scans, %d files (%0.2f GB) > %s (%0.2f GB), %0.0f s\n',sub,ses,nnz(ok),nf, ...
      mb / 1000,[name '.zip'],Fz.bytes / 1e9,toc(stime));
  else
    fprintf('  %s/%s: no scans exported\n',sub,ses);
  end
  for ui = find(~ok), cat_io_cprintf('err','    %s: %s\n',Us(ui).scan,Us(ui).msg); end
  if opts.anonymize > 1, msg = 'replaced by X'; else, msg = 'kept'; end
  for ui = find(ok & ~cellfun('isempty',{Us.names}))
    cat_io_cprintf('err','    %s: parts of the patient name in %s (%s) - check before the upload!\n', ...
      Us(ui).scan,strjoin(Us(ui).names,' '),msg);
  end
end
% =========================================================================
function Us = exportFiles(Us,Pses,opts)
%exportFiles. Anonymized DICOM files of the series Us in the directory Pses
%  (scan subdirectories). The files of each DICOM directory are read once
%  and assigned to the series by their SeriesInstanceUID. A series with a
%  file that cannot be exported is not exported.
  [Us.files] = deal({}); [Us.bytes] = deal(0); [Us.names] = deal({});
  for d = unique({Us.dcmdir})
    ud = find(strcmp({Us.dcmdir},d{1}) & cellfun('isempty',{Us.msg}));
    if isempty(ud), continue; end
    F = dir(d{1}); F = F(~[F.isdir] & ~startsWith({F.name},'.'));
    for fi = 1:numel(F)
      [b,uid] = readSeriesFile(fullfile(d{1},F(fi).name),{Us(ud).uid});
      ui = ud(strcmp({Us(ud).uid},uid));
      if isempty(ui) || ~isempty(Us(ui).msg), continue; end
      try
        [b,info] = anonymizeFile(b,Us(ui).plane,opts,Us(ui).X);
      catch e
        Us(ui).msg = sprintf('not exported (%s)',e.message);
        cellfun(@delete,Us(ui).files); Us(ui).files = {};
        continue
      end
      Pd = fullfile(Pses,Us(ui).scan);
      if ~exist(Pd,'dir'), mkdir(Pd); end
      if isnan(info.instance), info.instance = numel(Us(ui).files) + 1; end
      Pf = fullfile(Pd,sprintf('%04d_%06d.dcm',info.series,info.instance)); k = 1;
      while exist(Pf,'file'), k = k + 1; Pf = fullfile(Pd,sprintf('%04d_%06d_%d.dcm',info.series,info.instance,k)); end
      fid = fopen(Pf,'w'); fwrite(fid,b,'uint8'); fclose(fid);
      Us(ui).files{end+1} = Pf;
      Us(ui).bytes        = Us(ui).bytes + numel(b);
      Us(ui).names        = unique([Us(ui).names info.names]);
    end
  end
  for ui = find(cellfun('isempty',{Us.msg}) & cellfun('isempty',{Us.files}))
    Us(ui).msg = 'not exported, no readable DICOM files (explicit VR little endian)';
  end
end
% =========================================================================
function [b,uid] = readSeriesFile(P,uids)
%readSeriesFile. Bytes b of the DICOM file P if its SeriesInstanceUID uid
%  (read from the first 64 kB) is one of uids, otherwise empty.
  b = uint8([]); uid = '';
  fid = fopen(P,'r'); if fid < 0, return; end
  c = onCleanup(@() fclose(fid));
  b = fread(fid,65536,'*uint8')';
  for full = [false true]
    if full, b = [b fread(fid,inf,'*uint8')']; end %#ok<AGROW>
    try
      E   = dicomElements(b,hex2dec('0020000E'));
      uid = elementValue(b,E,find(E.depth == 0 & E.tag == hex2dec('0020000E'),1));
    catch
      uid = '';
    end
    if ~isempty(uid) || numel(b) < 65536, break; end
  end
  if ~any(strcmp(uids,uid)), b = uint8([]); return; end
  b = [b fread(fid,inf,'*uint8')']; % rest of the file
end
% =========================================================================
function [b,info] = anonymizeFile(b,plane,opts,X)
%anonymizeFile. Cleaned header and, for a face plane (see facePlane), the
%  defaced image of the DICOM file b (bytes), with the series and instance
%  number and the elements with parts of the patient name (see nameCheck)
%  as info. Nested elements are blanked (same length), top-level elements
%  get new values (see cleanedElements). For the anonymization level 2 (X 
%  with the subject and session label, the date offsets of the session and
%  subject, site label, institution, city, model, and coil names, and the
%  secret salt, see 
%  cat_io_dcm2bids_anonymize), the PatientID and PatientName are the 
%  subject label, the StudyID is the session label, the scanner and 
%  institution are the site label (without location), the
%  address, station, model, coils, department, and descriptions of the 
%  study and procedure are removed (see anonymizedElements), and in place
%  all dates and times are
%  shifted (DA, DT, and TM to the session start at the center of its 
%  3-hour interval of the day, the birth date with the offset of the 
%  subject, and dates in private elements, see
%  shiftSession in cat_io_dcm2bids_anonymize), all 
%  UIDs (UI and in private elements, e.g. the Siemens CSA header) are 
%  replaced by UIDs of the same length, and parts of the patient name, the
%  PatientID, the serial number, and the institution, city, model, and coil
%  names (X.words) in further elements are replaced by X.
  if ~exist('X','var'), X = []; end
  E   = dicomElements(b);
  top = E.depth == 0;
  el  = @(tag) find(top & E.tag == hex2dec(tag),1);
  ts  = elementValue(b,E,el('00020010'));
  if ~strcmp(ts,'1.2.840.10008.1.2.1')
    error('cat_io_dcm2bids_export:syntax','transfer syntax %s is not supported',ts);
  end
  info.series   = str2double(elementValue(b,E,el('00200011')));
  info.instance = str2double(elementValue(b,E,el('00200013')));
  pid = elementValue(b,E,el('00100020'));
  if isempty(X), new = cell(0,2); ntag = []; else, new = anonymizedElements(X); ntag = hex2dec(new(:,1)); end

  % cleaned elements of all depths, the PatientName, nested PatientIDs
  % (e.g. of other patient IDs), and the elements of the level 2
  cln = ismember(E.tag,[hex2dec([cleanedElements; {'00100010'}]); ntag]) | ...
    (E.depth > 0 & E.tag == hex2dec('00100020'));
  cln = cln & ~strcmp(cellstr(E.vr),'SQ') & ~isnan(E.len) & E.len > 0;
  parts = unique(strsplit(lower(elementValue(b,E,el('00100010'))),{'^',' ','-','_',',','.'}));
  parts = parts(cellfun(@numel,parts) >= 3 & ~cellfun(@(p) contains(lower(pid),p),parts));
  info.names = nameCheck(b,E,parts,cln);

  if ~isempty(plane), b = defacePixels(b,E,plane); end
  for ei = find(cln & ~top)'
    b(E.voff(ei) + (0:E.len(ei) - 1)) = 32;
  end

  % anonymization level 2: dates, UIDs, and names and IDs (in place)
  if ~isempty(X)
    vr = cellstr(E.vr); ok = ~cln & ~isnan(E.len) & E.len > 0;
    sn = elementValue(b,E,el('00181000')); dates = {};
    for ei = find(ok & ismember(vr,{'DA','DT','TM'}))'
      v = E.voff(ei) + (0:E.len(ei) - 1);
      if ~strcmp(vr{ei},'TM'), dates = [dates regexp(char(b(v)),'(?<!\d)\d{8}','match')]; end %#ok<AGROW>
      if E.tag(ei) == hex2dec('00100030') % birth date of the subject
        b(v) = uint8(cat_io_dcm2bids_anonymize('shiftDate',char(b(v)),X.boffset));
      else % dates and times of the session (start at the center of its 3-hour interval)
        b(v) = uint8(cat_io_dcm2bids_anonymize('shiftSession',char(b(v)),vr{ei},X));
      end
    end
    for ei = find(ok & strcmp(vr,'UI'))'
      v = E.voff(ei) + (0:E.len(ei) - 1); u = char(b(v));
      h = strjoin(cellfun(@(x) cat_io_dcm2bids_anonymize('hashUID',x,X.salt),strsplit(strtrim(u(u ~= char(0))),'\'), ...
        'UniformOutput',false),'\');
      b(v) = uint8([h zeros(1,numel(v) - numel(h))]);
    end
    for ei = find(ok & mod(floor(E.tag / 65536),2) == 1 & ~ismember(vr,{'SQ','DA','DT','TM','UI'}))' % further private elements
      v = E.voff(ei) + (0:E.len(ei) - 1);
      b(v) = cat_io_dcm2bids_anonymize('hashUIDsInBytes',b(v),X.salt);
      b(v) = cat_io_dcm2bids_anonymize('shiftDatesInBytes',b(v),X,unique(dates));
      if numel(sn) >= 4 % serial number also in parts of UIDs (e.g. in compressed data)
        [i0,i1] = regexp(char(b(v)),['(?<!\d)' regexptranslate('escape',sn) '(?!\d)'],'start','end');
        for k = 1:numel(i0), b(v(i0(k):i1(k))) = 'X'; end
      end
    end
    ids = lower([{pid sn} X.words]);
    [~,pos] = nameCheck(b,E,[parts ids(cellfun(@numel,ids) >= 4)],cln);
    for k = 1:size(pos,1), b(pos(k,1):pos(k,2)) = 'X'; end
    pid = X.code;
  end

  % top-level elements with new values
  chg = cell(0,2);
  for ei = find(cln & top)'
    k = find(ntag == E.tag(ei),1);
    if E.tag(ei) == hex2dec('00100010'), v = pid; elseif isempty(k), v = ''; else, v = new{k,2}; end
    chg(end+1,:) = {ei,v}; %#ok<AGROW>
  end
  bd = el('00100030'); d = elementValue(b,E,bd);
  if opts.anonymize > 0 && numel(d) == 8 && all(isstrprop(d,'digit')) % (level 2: shifted date)
    chg(end+1,:) = {bd,[d(1:6) '15']};
  end
  b = spliceValues(b,E,chg);
end
% =========================================================================
function v = fieldText(S,fn)
%fieldText. Text of the field fn of the structure S ('' if not available).
  if isfield(S,fn) && ~isempty(S.(fn)) && ~isnumeric(S.(fn)), v = char(S.(fn)); else, v = ''; end
end
% =========================================================================
function new = anonymizedElements(X)
%anonymizedElements. Elements (group and element) with their new values of
%  the anonymization level 2 (X with the subject and session label and the
%  site label without location, see anonymizeFile). Nested elements are 
%  blanked. The SoftwareVersions are kept, as dcm2niix uses them for the 
%  conversion (e.g. of Siemens XA data).
  new = {
    '00100020', X.code              % PatientID
    '00181000', X.site              % DeviceSerialNumber
    '00080080', X.site              % InstitutionName
    '00080081', ''                  % InstitutionAddress
    '00081090', ''                  % ManufacturerModelName
    '00181250', ''                  % ReceiveCoilName
    '00181251', ''                  % TransmitCoilName
    '00081040', ''                  % InstitutionalDepartmentName
    '00081010', ''                  % StationName
    '00200010', X.ses               % StudyID (session label, e.g. w0052)
    '00081030', ''                  % StudyDescription
    '00400254', ''                  % PerformedProcedureStepDescription
    '00321060', ''                  % RequestedProcedureDescription
    '00400007', ''                  % ScheduledProcedureStepDescription
    };
end
% =========================================================================
function tags = cleanedElements
%cleanedElements. Elements (group and element) that are cleaned besides
%  the PatientName (replaced by the PatientID) and the PatientBirthDate
%  (15th of the month for opts.anonymize > 0).
  tags = {
    '00080050' % AccessionNumber
    '00080090' % ReferringPhysicianName
    '00080092' % ReferringPhysicianAddress
    '00080094' % ReferringPhysicianTelephoneNumbers
    '00081048' % PhysiciansOfRecord
    '00081050' % PerformingPhysicianName
    '00081060' % NameOfPhysiciansReadingStudy
    '00081070' % OperatorsName
    '00101000' % OtherPatientIDs
    '00101001' % OtherPatientNames
    '00101005' % PatientBirthName
    '00101040' % PatientAddress
    '00101060' % PatientMotherBirthName
    '00101090' % MedicalRecordLocator
    '00102000' % MedicalAlerts
    '00102110' % Allergies
    '00102154' % PatientTelephoneNumbers
    '00102180' % Occupation
    '001021B0' % AdditionalPatientHistory
    '00104000' % PatientComments
    '00321032' % RequestingPhysician
    '00324000' % StudyComments
    '00380010' % AdmissionID
    '00400006' % ScheduledPerformingPhysicianName
    '00400009' % ScheduledProcedureStepID
    '00400253' % PerformedProcedureStepID
    '00400280' % CommentsOnThePerformedProcedureStep
    '00401001' % RequestedProcedureID
    '00401400' % RequestedProcedureComments
    '00402400' % ImagingServiceRequestComments
    };
end
% =========================================================================
function [tags,pos] = nameCheck(b,E,words,cln)
%nameCheck. Elements of the header (before the pixel data, without the
%  cleaned elements cln) that contain one of the words (lower case, e.g. 
%  the parts of the patient name of at least 3 characters that are not 
%  part of the PatientID, in a description) and the first and last byte of
%  each occurrence (pos). A word is not part of another word or number.
  tags = {}; pos = zeros(0,2);
  if isempty(words), return; end
  pd = find(E.depth == 0 & E.tag == hex2dec('7FE00010'),1);
  if isempty(pd), h = lower(char(b)); else, h = lower(char(b(1:E.hoff(pd) - 1))); end
  for qi = 1:numel(words)
    num = all(isstrprop(words{qi},'digit'));
    inword = @(c) isstrprop(c,'alphanum') || (num && c == '.'); % numbers also with decimals
    for k = strfind(h,words{qi})
      e = k + numel(words{qi});
      if (k > 1 && inword(h(k - 1))) || (e <= numel(h) && inword(h(e))), continue; end
      ei = find(E.voff <= k & E.voff + E.len > k,1,'last');
      if isempty(ei) || cln(ei), continue; end
      tags{end+1} = sprintf('(%04X,%04X)',floor(E.tag(ei) / 65536),mod(E.tag(ei),65536)); %#ok<AGROW>
      pos(end+1,:) = [k e - 1]; %#ok<AGROW>
    end
  end
  tags = unique(tags);
end
% =========================================================================
function b = defacePixels(b,E,w)
%defacePixels. Set the pixels in front of the face plane w (world
%  coordinates RAS, see facePlane) to 0 for single-frame, multi-frame
%  (enhanced), and Siemens mosaic images (8 or 16 bit). The position of
%  the pixel in row r and column c (from 0) of a frame is p + c*sc*u +
%  r*sr*v in the LPS coordinates of DICOM, with the position p of its first
%  pixel, the directions u (along a row) and v (along a column) of the
%  image orientation, and the pixel spacing [sr sc]. The tiles of a mosaic
%  are the slices in the direction of the slice normal of the CSA header.
  top  = E.depth == 0;
  el   = @(tag) find(top & E.tag == hex2dec(tag),1);
  rows = elementNumbers(b,E,el('00280010')); cols = elementNumbers(b,E,el('00280011'));
  bits = elementNumbers(b,E,el('00280100')); spp  = elementNumbers(b,E,el('00280002'));
  pd   = el('7FE00010');
  if isempty(pd) || isnan(E.len(pd)) || isempty(rows) || isempty(cols) || isempty(bits) || ...
      (~isempty(spp) && spp ~= 1)
    error('cat_io_dcm2bids_export:pixelData','no supported pixel data');
  end
  switch bits
    case 16, cls = 'uint16'; nb = 2;
    case 8,  cls = 'uint8';  nb = 1;
    otherwise, error('cat_io_dcm2bids_export:bits','%d bit images are not supported',bits);
  end
  nf  = floor(E.len(pd) / (rows * cols * nb));
  pix = E.voff(pd) + (0:rows * cols * nf * nb - 1);
  Y   = reshape(typecast(b(pix),cls),cols,rows,nf);

  % position (nf x 3), orientation (nf x 6), and spacing (nf x 2) of the frames
  fg = el('52009230'); % per-frame functional groups of multi-frame images
  if isempty(fg)
    P = elementNumbers(b,E,el('00200032')); O = elementNumbers(b,E,el('00200037')); S = elementNumbers(b,E,el('00280030'));
  else
    sg = el('52009229'); % shared functional groups
    P  = frameValues(b,E,'00200032',fg,sg,nf); O = frameValues(b,E,'00200037',fg,sg,nf); S = frameValues(b,E,'00280030',fg,sg,nf);
  end
  if size(P,1) ~= nf || size(P,2) ~= 3 || size(O,2) ~= 6 || size(S,2) ~= 2
    error('cat_io_dcm2bids_export:geometry','no image geometry');
  end
  [c,r] = ndgrid(0:cols - 1,0:rows - 1);
  k = 0; nv = [0 0 0]; dz = 0;
  if contains(elementValue(b,E,el('00080008')),'MOSAIC')
    N  = str2double(csaItems(b,E,'NumberOfImagesInMosaic'));
    nv = str2double(csaItems(b,E,'SliceNormalVector'));
    dz = elementNumbers(b,E,el('00180088')); if isempty(dz), dz = elementNumbers(b,E,el('00180050')); end
    if isempty(N) || isnan(N(1)) || numel(nv) ~= 3 || any(isnan(nv)) || isempty(dz) || nf ~= 1
      error('cat_io_dcm2bids_export:mosaic','no mosaic geometry');
    end
    m  = ceil(sqrt(N(1))); tr = rows / m; tc = cols / m;
    k  = floor(r / tr) * m + floor(c / tc); r = mod(r,tr); c = mod(c,tc); % slice and pixel of the tile
    P  = P + O(4:6) * S(1) * (rows - tr) / 2 + O(1:3) * S(2) * (cols - tc) / 2; % first pixel of the first slice
  end
  a = w(1:3) .* [-1 -1 1]; % plane in LPS coordinates
  for fi = 1:nf
    d = a * P(fi,:)' + w(4) + a * O(fi,1:3)' * S(fi,2) * c + a * O(fi,4:6)' * S(fi,1) * r + a * nv(:) * dz(1) * k;
    f = Y(:,:,fi); f(d < 0) = 0; Y(:,:,fi) = f;
  end
  b(pix) = typecast(Y(:)','uint8');
end
% =========================================================================
function X = frameValues(b,E,tag,fg,sg,nf)
%frameValues. Values of the element tag of each frame (rows) of the
%  per-frame functional groups fg (by item, the first element of an item,
%  i.e. of the standard functional groups before private sequences, e.g. 
%  the Philips copy of the position) or of the shared functional groups sg
%  (the same values for all nf frames).
  ei = find(E.tag == hex2dec(tag) & E.top == fg);
  if isempty(ei) && ~isempty(sg), ei = find(E.tag == hex2dec(tag) & E.top == sg,1); end
  [~,ui] = unique(E.item(ei),'first'); ei = ei(ui);
  X = cell2mat(arrayfun(@(e) elementNumbers(b,E,e),ei(:),'UniformOutput',false));
  if size(X,1) == 1, X = repmat(X,nf,1); end
end
% =========================================================================
function v = csaItems(b,E,name)
%csaItems. Items (strings) of the element name of the Siemens CSA image
%  header (0029,1010) in the format SV10, e.g. NumberOfImagesInMosaic.
  v  = {};
  ei = find(E.depth == 0 & E.tag == hex2dec('00291010'),1);
  if isempty(ei) || isnan(E.len(ei)) || E.len(ei) < 16, return; end
  c = b(E.voff(ei) + (0:E.len(ei) - 1));
  if ~isequal(char(c(1:4)),'SV10'), return; end
  i32 = @(p) double(typecast(c(p:p+3),'int32'));
  p = 17;
  for ti = 1:i32(9)
    tname  = char(c(p:p+63)); tname = tname(1:find([tname char(0)] == char(0),1) - 1);
    nitems = i32(p + 76); p = p + 84;
    items  = cell(1,nitems);
    for ii = 1:nitems
      L = i32(p + 4); p = p + 16;
      s = char(c(p:p + L - 1)); s(s == char(0)) = []; items{ii} = strtrim(s);
      p = p + 4 * ceil(L / 4);
    end
    if strcmp(tname,name), v = items(~cellfun('isempty',items)); return; end
  end
end
% =========================================================================
function E = dicomElements(b,stoptag)
%dicomElements. Elements of the DICOM file b (bytes of a part 10 file with
%  explicit VR little endian) including the elements of sequences:
%    E.tag   .. group*65536 + element
%    E.vr    .. value representation (2 characters per row)
%    E.hoff  .. 1-based offset of the element
%    E.voff  .. 1-based offset of its value
%    E.len   .. value length (NaN if undefined)
%    E.lenb  .. bytes of the length field (2 or 4)
%    E.depth .. sequence depth (0 for top-level elements)
%    E.top   .. index of the top-level element
%    E.item  .. item of the top-level sequence (e.g. the frame of the
%               per-frame functional groups)
%  With stoptag, only the top-level elements up to this tag are read (e.g.
%  of the first bytes of a file). Undefined lengths are supported for
%  sequences, also for sequences with the VR UN (implicit VR), but not for
%  encapsulated (compressed) data.
  if ~exist('stoptag','var'), stoptag = inf; end
  n = numel(b);
  if n < 132 || ~isequal(char(b(129:132)),'DICM')
    error('cat_io_dcm2bids_export:noDICOM','No DICOM part 10 file.');
  end
  LONG = {'OB','OW','OF','SQ','UT','UN','OD','OL','OV','SV','UC','UR','UV'};
  cap = 1024; k = 0;
  [tag,hoff,voff,len,lenb,depth,top,item] = deal(zeros(cap,1)); vr = repmat(' ',cap,2);
  pos = 133; itemnr = 0;
  walkElements(n + 1,0,0,false);
  E = struct('tag',tag(1:k),'vr',vr(1:k,:),'hoff',hoff(1:k),'voff',voff(1:k),'len',len(1:k), ...
    'lenb',lenb(1:k),'depth',depth(1:k),'top',top(1:k),'item',item(1:k));

  function walkElements(stop,d,t,implicit)
    while pos < stop
      g = double(b(pos)) + 256 * double(b(pos + 1));
      if g == 65534, return; end % item or sequence delimiter of the caller
      e = double(b(pos + 2)) + 256 * double(b(pos + 3));
      if d == 0 && g * 65536 + e > stoptag, pos = stop; return; end
      if implicit
        v = 'UN'; L = u32(pos + 4); vo = pos + 8; lb = 4;
      else
        v = char(b(pos + 4:pos + 5));
        if any(strcmp(v,LONG))
          L = u32(pos + 8); vo = pos + 12; lb = 4;
        else
          L = double(b(pos + 6)) + 256 * double(b(pos + 7)); vo = pos + 8; lb = 2;
        end
      end
      k = k + 1;
      if k > numel(tag) % grow
        tag(2*k) = 0; hoff(2*k) = 0; voff(2*k) = 0; len(2*k) = 0; lenb(2*k) = 0;
        depth(2*k) = 0; top(2*k) = 0; item(2*k) = 0; vr(2*k,:) = ' ';
      end
      tag(k) = g * 65536 + e; vr(k,:) = v; hoff(k) = pos; voff(k) = vo; lenb(k) = lb; depth(k) = d;
      if d == 0, top(k) = k; t1 = k; itemnr = 0; else, top(k) = t; t1 = t; item(k) = itemnr; end
      if L == 4294967295
        if ~any(strcmp(v,{'SQ','UN'}))
          error('cat_io_dcm2bids_export:encapsulated','Encapsulated (compressed) data are not supported.');
        end
        len(k) = NaN; pos = vo;
        walkItems(inf,d + 1,t1,implicit || strcmp(v,'UN'));
        pos = pos + 8; % sequence delimiter
      elseif strcmp(v,'SQ')
        len(k) = L; pos = vo;
        walkItems(vo + L,d + 1,t1,implicit);
        pos = vo + L;
      else
        len(k) = L; pos = vo + L;
      end
    end
  end
  function walkItems(stop,d,t,implicit)
    while pos < stop
      e = double(b(pos + 2)) + 256 * double(b(pos + 3));
      if e == 57565, return; end % sequence delimiter (FFFE,E0DD)
      L = u32(pos + 4); pos = pos + 8;
      if d == 1, itemnr = itemnr + 1; end % items of top-level sequences
      if L == 4294967295
        walkElements(inf,d,t,implicit);
        pos = pos + 8; % item delimiter
      else
        walkElements(pos + L,d,t,implicit);
      end
    end
  end
  function x = u32(p)
    x = double(b(p)) + 256 * double(b(p + 1)) + 65536 * double(b(p + 2)) + 16777216 * double(b(p + 3));
  end
end
% =========================================================================
function v = elementValue(b,E,ei)
%elementValue. Text value of the element ei (without padding).
  if isempty(ei) || isnan(E.len(ei)) || E.len(ei) == 0, v = ''; return; end
  v = char(b(E.voff(ei) + (0:E.len(ei) - 1)));
  v = strtrim(v(v ~= char(0)));
end
% =========================================================================
function x = elementNumbers(b,E,ei)
%elementNumbers. Numeric values of the element ei (binary US, UL, SS, SL,
%  FL, and FD, otherwise text, e.g. DS and IS, with values separated by
%  backslashes).
  x = [];
  if isempty(ei) || isnan(E.len(ei)) || E.len(ei) == 0, return; end
  c = b(E.voff(ei) + (0:E.len(ei) - 1));
  switch E.vr(ei,:)
    case 'US', x = double(typecast(c,'uint16'));
    case 'UL', x = double(typecast(c,'uint32'));
    case 'SS', x = double(typecast(c,'int16'));
    case 'SL', x = double(typecast(c,'int32'));
    case 'FL', x = double(typecast(c,'single'));
    case 'FD', x = typecast(c,'double');
    otherwise, x = sscanf(strrep(char(c),'\',' '),'%f')';
  end
end
% =========================================================================
function b = spliceValues(b,E,chg)
%spliceValues. New values chg (element index and text) of top-level
%  elements (padded with a space to an even length).
  if isempty(chg), return; end
  [~,o] = sort([chg{:,1}]); chg = chg(o,:);
  seg = cell(1,3 * size(chg,1) + 1); p = 1;
  for ci = 1:size(chg,1)
    ei = chg{ci,1}; v = uint8(chg{ci,2}); if mod(numel(v),2), v(end+1) = 32; end %#ok<AGROW>
    h  = b(E.hoff(ei):E.voff(ei) - 1);
    if E.lenb(ei) == 2, h(end-1:end) = typecast(uint16(numel(v)),'uint8');
    else,               h(end-3:end) = typecast(uint32(numel(v)),'uint8'); end
    seg(3*ci-2:3*ci) = {b(p:E.hoff(ei) - 1), h, v};
    p = E.voff(ei) + E.len(ei);
  end
  seg{end} = b(p:end);
  b = [seg{:}];
end
% =========================================================================
function zipDirectory(P,sub,Pzip)
%zipDirectory. Zip file Pzip of the subdirectory sub of P (relative paths)
%  by the system zip with fast compression and without extra attributes
%  (macOS/Linux), otherwise (or if it fails) by the zip of MATLAB.
  if isunix
    q  = @(s) ['''' strrep(s,'''','''\''''') '''']; % single-quoted shell argument
    st = system(sprintf('cd %s && zip -r -q -X -1 %s %s',q(P),q(Pzip),q(sub)));
    if st == 0, return; end
    if exist(Pzip,'file'), delete(Pzip); end
  end
  zip(Pzip,sub,P);
end
% =========================================================================
function writeExportTables(Pexp,R,Z,opts)
%writeExportTables. Tables of the exported scans (scans.<fmt>) and zip
%  files (sessions.<fmt>) of all runs, and the readme file (see 
%  cat_io_dcm2bids_readme).
  T = {'scans',R; 'sessions',Z};
  for ti = 1:size(T,1)
    S = T{ti,2}; if isempty(S), continue; end
    C = struct2cell(S(:)); C = reshape(C,size(C,1),[])';
    cat_io_dcm2bids_helper('writeTable',fullfile(Pexp,[T{ti,1} '.' opts.tableformat]), ...
      [fieldnames(S)'; C],opts.tableformat);
  end
  cat_io_dcm2bids_readme('writeReadme','dicom',Pexp,struct('opts',opts));
end
% =========================================================================
function P = dbImage(Pjson,Pmdbdir)
%dbImage. Image of the sidecar Pjson (NIfTI, otherwise the sidecar of a 
%  JSON-only import) as subpath of the database Pmdbdir (with /), without
%  the DICOM directory of the file name (folder=protocol=time=series), as 
%  the original names can contain names.
  P = spm_file(Pjson,'ext','.nii.gz'); 
  if ~exist(P,'file'), P = spm_file(Pjson,'ext','.nii'); end
  if ~exist(P,'file'), P = Pjson; end
  rel = regexp(P,['^' regexptranslate('escape',Pmdbdir) '[/\\](.*)$'],'tokens','once');
  if isempty(rel), rel = regexp(P,'catDCM2BIDSdb[/\\](.*)$','tokens','once'); end
  if ~isempty(rel), P = rel{1}; end
  [pth,nam,ext] = fileparts(P); 
  P = strrep(fullfile(pth,[regexprep(nam,'^[^=]*=','') ext]),filesep,'/');
end
% =========================================================================
function removeDir(P)
%removeDir. Remove the directory P (if it exists).
  if exist(P,'dir'), rmdir(P,'s'); end
end

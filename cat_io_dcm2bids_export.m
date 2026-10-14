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
%  - birth date: the 15th of the month for opts.anonymize > 0 (the age
%    estimate of the import uses the 15th, i.e. the original and exported
%    data give the same age, see assurePatientDCMfields in cat_io_dcm2bids_db)
%  - kept: PatientID (that has to be a pseudonym), sex, age, size, weight,
%    dates and times, UIDs, scanner, institution, protocol and sequence
%    fields, image comments (e.g. the multiband factor), private elements
%  - check: the parts of the patient name are searched in the header of
%    each file (warning with the elements, e.g. a name in a description)
%  Image: the face is removed by the plane of spm_deface (registration of
%  the image of the series in the database or of its mean image for 4D
%  scans) as by the anonymization of the BIDS export, i.e. for anatomical
%  scans (opts.anonymize = 1) or all scans (opts.anonymize > 1). Pixels in
%  front of the plane are set to 0 (single-frame, multi-frame, and Siemens
%  mosaic images). Scans that cannot be defaced are not exported.
%
%  varargout = cat_io_dcm2bids_export(action,varargin)
%
%  Actions (see the help of the local functions):
%    cat_io_dcm2bids_export('exportDICOM',Pjson,datatype,status,Pmdbdir,opts)
%    [b,info] = cat_io_dcm2bids_export('anonymizeFile',b,plane,opts)
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
function exportDICOM(Pjson,datatype,status,Pmdbdir,opts)
%exportDICOM. Export of the DICOM series of the database scans Pjson
%  (sidecars) with their datatype and protocol status (see testProtocols
%  in cat_io_dcm2bids_protocols) into the directory <Pmdbdir>-dicom, i.e.
%  of the series that were not exported before (for this anonymization
%  level). The lists of the tables are saved after each session.
  Pexp = [Pmdbdir '-dicom'];
  Plst = fullfile(Pexp,'exports.mat');
  R = struct([]); Z = struct([]);
  if exist(Plst,'file'), load(Plst,'R','Z'); end

  fprintf('\nDCM2BIDS - anonymized DICOM export\n');
  U = selectSeries(Pjson,datatype,status,opts);
  if ~isempty(R) && ~isempty(U)
    key = @(uid,an) strcat(uid,'|',arrayfun(@(a) sprintf('%d',a),an,'UniformOutput',false));
    U   = U(~ismember(key({U.uid},repmat(opts.anonymize,1,numel(U))),key({R.SeriesInstanceUID},[R.Anonymize])));
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
    [Rs,Zs] = exportSession(Us,Prun,run,opts);
    R = [R Rs]; Z = [Z Zs]; %#ok<AGROW>
    if ~exist(Pexp,'dir'), mkdir(Pexp); end
    save(Plst,'R','Z');
    writeExportTables(Pexp,R,Z,opts);
  end
  if exist(Prun,'dir') && isempty(dir(fullfile(Prun,'*.zip'))), rmdir(Prun,'s'); end
end
% =========================================================================
function U = selectSeries(Pjson,datatype,status,opts)
%selectSeries. Export units of the scans Pjson (database sidecars), i.e.
%  one unit per DICOM series of the scans with conform or semiconform
%  protocol (all scans without protocols) for opts.dicomexport = 1 or of
%  all scans for opts.dicomexport = 2, without BIDS inputs (no DICOM
%  files). A series is defaced if one of its scans is defaced by the
%  anonymization (anatomical scans or all scans for opts.anonymize > 1,
%  see anonymize in cat_io_dcm2bids_pp_anat). Units that cannot be
%  exported get a message.
  U = struct('uid',{},'Pjson',{},'scan',{},'snrdir',{},'dcmdir',{},'deface',{},'datatype',{}, ...
    'status',{},'series',{},'protocol',{},'plane',{},'msg',{});
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
    deface = opts.anonymize > 1 || (opts.anonymize == 1 && strcmp(datatype{si},'anat'));
    ui = find(strcmp({U.uid},J.SeriesInstanceUID),1);
    if isempty(ui)
      snrdir = fileparts(Pjson{si}); [~,scan] = fileparts(snrdir);
      fnameparts = strsplit(spm_file(Pjson{si},'basename'),'=');
      msg = '';
      if ~isfield(J,'PatientID') || isempty(J.PatientID)
        msg = 'not exported, no PatientID (the database subject would be named by the patient name)';
      elseif ~exist(I.DICOMdirectory,'dir')
        msg = 'not exported, DICOM directory not available';
      end
      U(end+1) = struct('uid',J.SeriesInstanceUID,'Pjson',{Pjson(si)},'scan',scan,'snrdir',snrdir, ...
        'dcmdir',I.DICOMdirectory,'deface',deface,'datatype',datatype{si},'status',status{si}, ...
        'series',J.SeriesNumber,'protocol',fnameparts{min(2,end)},'plane',[],'msg',msg); %#ok<AGROW>
    else
      U(ui).Pjson{end+1} = Pjson{si};
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
function [Rs,Zs] = exportSession(Us,Prun,run,opts)
%exportSession. Zip file of the series Us of a session with the anonymized
%  DICOM files (see exportFiles) and its content table, and the records of
%  the export tables (Rs for the scans and Zs for the zip file).
  stime = tic; Rs = struct([]); Zs = struct([]);
  [Psub,ses] = fileparts(fileparts(Us(1).snrdir)); [~,sub] = fileparts(Psub);
  name   = [sub '_' ses];
  Pstage = fullfile(Prun,['.' name]); % files of the zip (removed at the end, also after errors)
  c = onCleanup(@() removeDir(Pstage));
  Us = exportFiles(Us,fullfile(Pstage,sub,ses),opts);
  ok = cellfun('isempty',{Us.msg});
  if any(ok)
    Pzip = fullfile(Prun,[name '.zip']);
    zipDirectory(Pstage,sub,Pzip);
    F   = dir(fullfile(Pstage,'**','*.dcm'));
    rel = strrep(strrep(fullfile({F.folder},{F.name}),[Pstage filesep],''),filesep,'/');
    cat_io_dcm2bids_helper('writeTable',fullfile(Prun,[name '.' opts.tableformat]), ...
      [{'File','Bytes'}; rel(:) {F.bytes}'],opts.tableformat);
    Fz = dir(Pzip); nf = sum(cellfun(@numel,{Us(ok).files})); mb = sum([Us(ok).bytes]) / 1e6;
    for ui = find(ok)
      Rs = [Rs struct('Export',run,'Zip',[name '.zip'],'Subject',sub,'Session',ses,'Scan',Us(ui).scan, ...
        'SeriesNumber',Us(ui).series,'ProtocolName',Us(ui).protocol,'Datatype',Us(ui).datatype, ...
        'Status',Us(ui).status,'Files',numel(Us(ui).files),'MB',round(Us(ui).bytes / 1e6,1), ...
        'Defaced',double(~isempty(Us(ui).plane)),'Anonymize',opts.anonymize, ...
        'NameCheck',strjoin(Us(ui).names,' '),'SeriesInstanceUID',Us(ui).uid)]; %#ok<AGROW>
    end
    Zs = struct('Export',run,'Zip',[name '.zip'],'Subject',sub,'Session',ses,'Scans',nnz(ok), ...
      'Files',nf,'MB',round(mb,1),'ZipMB',round(Fz.bytes / 1e6,1),'Duration',round(toc(stime)));
    fprintf('  %s/%s: %d scans, %d files (%0.2f GB) > %s (%0.2f GB), %0.0f s\n',sub,ses,nnz(ok),nf, ...
      mb / 1000,[name '.zip'],Fz.bytes / 1e9,toc(stime));
  else
    fprintf('  %s/%s: no scans exported\n',sub,ses);
  end
  for ui = find(~ok), cat_io_cprintf('err','    %s: %s\n',Us(ui).scan,Us(ui).msg); end
  for ui = find(ok & ~cellfun('isempty',{Us.names}))
    cat_io_cprintf('err','    %s: parts of the patient name in %s - check before the upload!\n', ...
      Us(ui).scan,strjoin(Us(ui).names,' '));
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
        [b,info] = anonymizeFile(b,Us(ui).plane,opts);
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
function [b,info] = anonymizeFile(b,plane,opts)
%anonymizeFile. Cleaned header and, for a face plane (see facePlane), the
%  defaced image of the DICOM file b (bytes), with the series and instance
%  number and the elements with parts of the patient name (see nameCheck)
%  as info. Nested elements are blanked (same length), top-level elements
%  get new values (see cleanedElements).
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

  % cleaned elements of all depths, the PatientName, and nested PatientIDs
  % (e.g. of other patient IDs)
  cln = ismember(E.tag,hex2dec(cleanedElements)) | E.tag == hex2dec('00100010') | ...
    (E.depth > 0 & E.tag == hex2dec('00100020'));
  cln = cln & ~strcmp(cellstr(E.vr),'SQ') & ~isnan(E.len) & E.len > 0;
  info.names = nameCheck(b,E,elementValue(b,E,el('00100010')),pid,cln);

  if ~isempty(plane), b = defacePixels(b,E,plane); end
  for ei = find(cln & ~top)'
    b(E.voff(ei) + (0:E.len(ei) - 1)) = 32;
  end
  chg = cell(0,2);
  for ei = find(cln & top)'
    if E.tag(ei) == hex2dec('00100010'), chg(end+1,:) = {ei,pid}; else, chg(end+1,:) = {ei,''}; end %#ok<AGROW>
  end
  bd = el('00100030'); d = elementValue(b,E,bd);
  if opts.anonymize > 0 && numel(d) == 8 && all(isstrprop(d,'digit'))
    chg(end+1,:) = {bd,[d(1:6) '15']};
  end
  b = spliceValues(b,E,chg);
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
function tags = nameCheck(b,E,pn,pid,cln)
%nameCheck. Elements of the header (before the pixel data, without the
%  cleaned elements cln) that contain a part of the patient name pn as
%  word (at least 3 characters and not part of the PatientID pid), e.g. a
%  name in a description.
  tags  = {};
  parts = unique(strsplit(lower(pn),{'^',' ','-','_',',','.'}));
  parts = parts(cellfun(@numel,parts) >= 3 & ~cellfun(@(p) contains(lower(pid),p),parts));
  if isempty(parts), return; end
  pd = find(E.depth == 0 & E.tag == hex2dec('7FE00010'),1);
  if isempty(pd), h = lower(char(b)); else, h = lower(char(b(1:E.hoff(pd) - 1))); end
  for qi = 1:numel(parts)
    for k = strfind(h,parts{qi})
      e = k + numel(parts{qi});
      if (k > 1 && isletter(h(k - 1))) || (e <= numel(h) && isletter(h(e))), continue; end % no word
      ei = find(E.voff <= k & E.voff + E.len > k,1,'last');
      if isempty(ei) || cln(ei), continue; end
      tags{end+1} = sprintf('(%04X,%04X)',floor(E.tag(ei) / 65536),mod(E.tag(ei),65536)); %#ok<AGROW>
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
%  per-frame functional groups fg (by item) or of the shared functional
%  groups sg (the same values for all nf frames).
  ei = find(E.tag == hex2dec(tag) & E.top == fg);
  if isempty(ei) && ~isempty(sg), ei = find(E.tag == hex2dec(tag) & E.top == sg,1); end
  [~,o] = sort(E.item(ei)); ei = ei(o);
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
%  files (sessions.<fmt>) of all runs, and the readme file.
  T = {'scans',R; 'sessions',Z};
  for ti = 1:size(T,1)
    S = T{ti,2}; if isempty(S), continue; end
    C = struct2cell(S(:)); C = reshape(C,size(C,1),[])';
    cat_io_dcm2bids_helper('writeTable',fullfile(Pexp,[T{ti,1} '.' opts.tableformat]), ...
      [fieldnames(S)'; C],opts.tableformat);
  end
  txt = {
    'catDCM2BIDS anonymized DICOM export (see cat_io_dcm2bids_export), e.g. to transfer the imported data to a central site that runs the batch again.'
    ''
    'THIS IS NOT AN ARCHIVE OR BACKUP OF THE ORIGINAL DATA, AS THE HEADERS ARE CLEANED AND THE IMAGES ARE DEFACED. KEEP THE ORIGINAL DATA!'
    ''
    '  <yyyyMMdd-HHmmss>/      one directory per run with the new scans'
    '    sub-*_ses-*.zip       DICOM files of a session in the layout of the database (sub-*/ses-*/snr-*/<series>_<instance>.dcm)'
    sprintf('    sub-*_ses-*.%-7s   content of the zip file',opts.tableformat)
    sprintf('  scans.%-14s   exported scans of all runs',opts.tableformat)
    sprintf('  sessions.%-11s   exported zip files of all runs',opts.tableformat)
    '  exports.mat             lists of the tables: each scan is exported only once (for an anonymization level),'
    '                          i.e. the run directories can be removed after the upload'
    ''
    'Header: names of the patient (replaced by the PatientID), physicians, and operators, addresses, phone numbers,'
    'other IDs, and comments are cleaned, the birth date is set to the 15th of the month (anonymization), whereas the'
    'fields used by catDCM2BIDS are kept (e.g. PatientID, sex, age, dates, UIDs, scanner, protocol, private fields).'
    'Images: defaced (anatomical scans, all scans for the extensive anonymization). The PatientID has to be a pseudonym!'
    };
  cat_io_dcm2bids_helper('writeText',fullfile(Pexp,'readme.txt'),[strjoin(txt,newline) newline]);
end
% =========================================================================
function removeDir(P)
%removeDir. Remove the directory P (if it exists).
  if exist(P,'dir'), rmdir(P,'s'); end
end

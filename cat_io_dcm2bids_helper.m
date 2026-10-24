function varargout = cat_io_dcm2bids_helper(action,varargin)
%cat_io_dcm2bids_helper. Shared helper functions of cat_io_dcm2bids.
%  Text and JSON files with UTF-8 encoding (readText, writeText, writeJSON),
%  paths (absPath, normPath, canonPath), NIfTI header information, single
%  volumes, the gzip handling of NIfTIs, and their export (niiHeaderInfo, 
%  niiVolume, prepNii, prepNiigz, gzipFile, exportNii), the index of the
%  derivatives in the database (registerDerivative, derivativeIndex), the
%  realigned 
%  volumes of 4D data resliced in memory (realignedVolumes, 
%  realignedVolume), and tables (updateTable, writeTable, and lists of 
%  records with catRecords).
%
%  varargout = cat_io_dcm2bids_helper(action,varargin)
%
%  Actions (see the help of the local functions):
%    txt = cat_io_dcm2bids_helper('readText',P)
%    cat_io_dcm2bids_helper('writeText',P,txt)
%    cat_io_dcm2bids_helper('writeJSON',P,S)
%    P = cat_io_dcm2bids_helper('absPath',P)
%    P = cat_io_dcm2bids_helper('normPath',P)
%    P = cat_io_dcm2bids_helper('canonPath',P)
%    P = cat_io_dcm2bids_helper('derivativesDir',Pbids)
%    [vx,dim,mat] = cat_io_dcm2bids_helper('niiHeaderInfo',P)
%    Y = cat_io_dcm2bids_helper('niiVolume',P,vol)
%    Po = cat_io_dcm2bids_helper('prepNii',Pi,opts,del)
%    Po = cat_io_dcm2bids_helper('prepNiigz',Pi,opts)
%    cat_io_dcm2bids_helper('gzipFile',P)
%    cat_io_dcm2bids_helper('exportNii',P,Pout,opts)
%    cat_io_dcm2bids_helper('registerDerivative',P,Pout)
%    D = cat_io_dcm2bids_helper('derivativeIndex',Pmdbdir)
%    cat_io_dcm2bids_helper('removeAcquisitionTime',P)
%    [W,geom,Pmerged] = cat_io_dcm2bids_helper('conversionInfo',P)
%    V = cat_io_dcm2bids_helper('realignedVolumes',P,Prp,opts)
%    Y = cat_io_dcm2bids_helper('realignedVolume',V,t)
%    entries = cat_io_dcm2bids_helper('updateTable',Ptable,Thdr,Tnewrow,id,reimport)
%    cat_io_dcm2bids_helper('writeTable',Ptable,T,format)
%    R = cat_io_dcm2bids_helper('catRecords',R,S)
%
%  See also cat_io_dcm2bids.

  switch action
    case {'readText','writeText','writeJSON','absPath','normPath','canonPath','derivativesDir', ...
          'niiHeaderInfo','niiVolume','prepNii','prepNiigz','gzipFile','exportNii','removeAcquisitionTime', ...
          'registerDerivative','derivativeIndex', ...
          'conversionInfo','realignedVolumes','realignedVolume', ...
          'updateTable','writeTable','catRecords'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_helper:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function [W,geom,Pmerged] = conversionInfo(P)
%conversionInfo. Warnings of dcm2niix W (cell) in the sidecar of the 
%  database image P (also with the prefixes of the anonymization and
%  denoising, see ConversionWarnings in cat_io_dcm2bids_importDCM), if one
%  of them describes a wrong geometry (geom, missing slices or varying 
%  slice distances), and the image of the merged echoes Pmerged of such
%  echo images (ConversionMergedEchoes, '' if not available), which is only
%  used for the QC and the alignment (mixed echo times). 
  W = {}; geom = false; Pmerged = ''; 
  [pth,nam] = fileparts(regexprep(P,'\.gz$','')); 
  Pj = fullfile(pth,[regexprep(nam,'^(sanlm_)?(anon_)?','') '.json']); 
  if ~exist(Pj,'file'), return; end
  try
    V = jsondecode(readText(Pj)); 
  catch
    return
  end
  if isfield(V,'ConversionWarnings') && ~isempty(V.ConversionWarnings)
    W = cellstr(V.ConversionWarnings); 
    geom = any(cat_io_contains(W,{'Interslice distance varies','Missing images'})); 
  end
  if isfield(V,'ConversionMergedEchoes') && ischar(V.ConversionMergedEchoes)
    Pm = fullfile(pth,V.ConversionMergedEchoes); 
    if ~exist(Pm,'file'), Pm = regexprep(Pm,'\.gz$',''); end % (gzip status)
    if ~exist(Pm,'file'), Pm = [Pm '.gz']; end
    if exist(Pm,'file'), Pmerged = Pm; end
  end
end
% =========================================================================
function txt = readText(P)
%readText. Read a text file (UTF-8).
  fid = fopen(P,'r','n','UTF-8'); txt = fread(fid,'*char')'; fclose(fid); 
end
% =========================================================================
function writeText(P,txt)
%writeText. Write a text file (UTF-8).
  fid = fopen(P,'w','n','UTF-8'); fprintf(fid,'%s',txt); fclose(fid); 
end
% =========================================================================
function writeJSON(P,S)
%writeJSON. Write a JSON file (UTF-8, readable format).
  fid = fopen(P,'w','n','UTF-8'); fprintf(fid,'%s',jsonencode(S,'PrettyPrint',true)); fclose(fid); 
end
% =========================================================================
function P = derivativesDir(Pbids)
%derivativesDir. Derivatives directory of catDCM2BIDS of the BIDS dataset
%  Pbids (BIDS derivatives dataset, i.e. Pbids/derivatives/catDCM2BIDS).
  P = fullfile(Pbids,'derivatives','catDCM2BIDS');
end
% =========================================================================
function P = absPath(P)
%absPath. Absolute path (relative to the current MATLAB directory). 
  if ~( strncmp(P,'/',1) || strncmp(P,'\\',2) || ~isempty(regexp(P,'^[A-Za-z]:','once')) )
    P = fullfile(pwd,P); 
  end
end
% =========================================================================
function P = normPath(P)
%normPath. Absolute path without "." and ".." components, repeated and 
%  final file separators (without resolving symbolic links). 
  P = absPath(P); 
  % keep the root (/, drive letter, or UNC server) and normalize the rest
  if ispc, sep = '\\/'; else, sep = '/'; end
  root  = regexp(P,['^(/|\\\\[^' sep ']+|[A-Za-z]:)'],'match','once'); 
  parts = regexp(P(numel(root)+1:end),['[^' sep ']+'],'match'); 
  pnorm = {}; 
  for ppi = 1:numel(parts)
    if strcmp(parts{ppi},'..')
      pnorm = pnorm(1:end-1); 
    elseif ~strcmp(parts{ppi},'.')
      pnorm{end+1} = parts{ppi}; %#ok<AGROW>
    end
  end
  if strcmp(root,'/'), P = ['/' strjoin(pnorm,'/')]; else, P = strjoin([{root} pnorm],filesep); end
  P = regexprep(P,['[' regexptranslate('escape',filesep) ']+$'],''); 
end
% =========================================================================
function P = canonPath(P)
%canonPath. Absolute path (see normPath) with resolved symbolic links if 
%  Java is available, e.g. to compare directories. 
  P = normPath(P); 
  if usejava('jvm')
    try %#ok<TRYNC>
      P = char(java.io.File(P).getCanonicalPath()); 
    end
  end
  P = regexprep(P,['[' regexptranslate('escape',filesep) ']+$'],''); 
end
% =========================================================================
function [vx,dim,mat] = niiHeaderInfo(P)
%niiHeaderInfo. Voxel size (x, y, z), image dimensions (x, y, z, volumes),
%  and the voxel-to-world matrix mat (as spm_vol: sform if available, 
%  otherwise qform, for voxel indices starting at 1) from the NIfTI header
%  without reading the image data, i.e., for .nii.gz files only the header
%  is decompressed (whereas spm_vol decompresses the whole file, e.g. for
%  the field of view of long time series).
  vx = nan(1,3); dim = nan(1,4); mat = nan(4); 
  try
    if numel(P) > 3 && strcmpi(P(end-2:end),'.gz')
      gis = java.util.zip.GZIPInputStream(java.io.FileInputStream(P)); 
      hdr = zeros(1,540,'uint8'); 
      for i = 1:540
        v = gis.read(); if v < 0, break; end
        hdr(i) = v; 
      end
      gis.close(); 
    else
      fid = fopen(P,'r'); hdr = fread(fid,540,'*uint8')'; fclose(fid); 
    end
    hsz  = typecast(hdr(1:4),'int32');           % 348 (NIfTI-1) or 540 (NIfTI-2)
    swap = ~any(hsz == [348 540]); 
    if swap, hsz = swapbytes(hsz); end
    if hsz == 348
      dm = typecast(hdr(41:56),'int16');         % dim at byte offset 40
      pd = typecast(hdr(77:108),'single');       % pixdim at byte offset 76
    else
      dm = typecast(hdr(17:80),'int64');         % dim at byte offset 16
      pd = typecast(hdr(105:168),'double');      % pixdim at byte offset 104
    end
    if swap, dm = swapbytes(dm); pd = swapbytes(pd); end
    vx  = double(abs(pd(2:4))); 
    dim = double(dm(2:5)); dim(dim<1 | (1:4) > dm(1)) = 1; 
    if nargout > 2, mat = niiMat(hdr,hsz,swap,double(pd)); end
  catch
    % fallback that reads the image 
    try
      evalc('V = spm_vol(P);'); 
      vx  = sqrt(sum(V(1).mat(1:3,1:3).^2)); 
      dim = [V(1).dim ones(1,3-numel(V(1).dim)) numel(V)]; 
      mat = V(1).mat; 
    end
  end
  dim = reshape(dim,1,[]); vx = reshape(vx,1,[]); 
end
% =========================================================================
function M = niiMat(hdr,hsz,swap,pd)
%niiMat. Voxel-to-world matrix of the NIfTI header bytes hdr (header size
%  hsz, byte swap, pixdim pd) as the mat of spm_vol (see @nifti/subsref): 
%  the sform for sform_code > 0, otherwise the qform (quaternion, qfac in
%  pixdim(1)), otherwise the voxel size around the image center, scaled to 
%  mm by the spatial units, for voxel indices starting at 1.
  if hsz == 348 % NIfTI-1 (byte offsets of the codes, quaternion, and sform)
    code = double(val(253:256,'int16')); q = double(val(257:280,'single')); s = double(val(281:328,'single')); un = hdr(124);
    dm = double(val(41:56,'int16'));
  else          % NIfTI-2
    code = double(val(345:352,'int32')); q = double(val(353:400,'double')); s = double(val(401:496,'double')); un = hdr(501);
    dm = double(val(17:80,'int64'));
  end
  if code(2) > 0
    M = [reshape(s,4,3)'; 0 0 0 1];
  elseif code(1) > 0
    b = q(1); c = q(2); d = q(3); a = sqrt(max(0,1 - b^2 - c^2 - d^2));
    if a < 1e-7, n = sqrt(b^2 + c^2 + d^2); b = b/n; c = c/n; d = d/n; a = 0; end
    R = [a^2+b^2-c^2-d^2  2*(b*c-a*d)      2*(b*d+a*c)
         2*(b*c+a*d)      a^2+c^2-b^2-d^2  2*(c*d-a*b)
         2*(b*d-a*c)      2*(c*d+a*b)      a^2+d^2-b^2-c^2];
    Z = pd(2:4); Z(Z < 0) = 1; if pd(1) < 0, Z(3) = -Z(3); end
    M = [R * diag(Z) q(4:6)'; 0 0 0 1];
  else
    vx = pd(2:4); vx(vx <= 0) = 1; o = (dm(2:4) + 1) / 2; % (no orientation: center of the image)
    M = [diag(vx) (-vx(:) .* (o(:) - 1)); 0 0 0 1];
  end
  M = M * [eye(4,3) [-1 -1 -1 1]'];
  sc = [1000 1 0.001]; u = double(bitand(un,7));
  if u >= 1 && u <= 3, M = diag([sc(u)*[1 1 1] 1]) * M; end
  function v = val(i,type)
    v = typecast(hdr(i),type); if swap, v = swapbytes(v); end
  end
end
% =========================================================================
function Y = niiVolume(P,vol)
%niiVolume. One volume (default the first) of a NIfTI-1 image (.nii or 
%  .nii.gz, scaled values) without reading or unzipping the whole file, 
%  e.g. for the comparison of the first volumes of long time series. 
  if ~exist('vol','var') || isempty(vol), vol = 1; end
  gz = numel(P) > 3 && strcmpi(P(end-2:end),'.gz'); 
  if gz
    gis = java.util.zip.GZIPInputStream(java.io.FileInputStream(P),65536); 
    ch  = java.nio.channels.Channels.newChannel(gis); 
    hdr = readBytes(ch,348); 
  else
    fid = fopen(P,'r'); 
    hdr = fread(fid,348,'*uint8')'; 
  end
  try
    if typecast(hdr(1:4),'int32') ~= 348 && swapbytes(typecast(hdr(1:4),'int32')) ~= 348
      error('cat_io_dcm2bids_helper:niiVolume','Only NIfTI-1 images are supported.'); 
    end
    swap = typecast(hdr(1:4),'int32') ~= 348; 
    hval = @(i,type) typecast(hdr(i),type); 
    dm   = double(hval(41:56,'int16')); 
    dt   = double(hval(71:72,'int16')); 
    bp   = double(hval(73:74,'int16')); 
    vo   = double(hval(109:112,'single')); 
    sl   = double(hval(113:116,'single')); 
    si   = double(hval(117:120,'single')); 
    if swap, dm = double(swapbytes(int16(dm))); dt = double(swapbytes(int16(dt))); bp = double(swapbytes(int16(bp))); 
      vo = double(swapbytes(single(vo))); sl = double(swapbytes(single(sl))); si = double(swapbytes(single(si))); end
    dim  = [dm(2:4) 1 1 1]; dim = dim(1:3); dim(dim < 1) = 1; 
    nb   = prod(dim) * bp/8; 
    off  = vo + (min(vol,max(1,prod(dm(5:min(end,dm(1)+1))))) - 1) * nb; 
    if gz
      for skip = [repmat(2^24,1,floor((off - 348)/2^24)) mod(off - 348,2^24)]
        readBytes(ch,skip); 
      end
      data = readBytes(ch,nb); 
      gis.close(); 
    else
      fseek(fid,off,'bof'); 
      data = fread(fid,nb,'*uint8')'; 
      fclose(fid); 
    end
    types = {2,'uint8'; 4,'int16'; 8,'int32'; 16,'single'; 64,'double'; 256,'int8'; 512,'uint16'; 768,'uint32'}; 
    Y = typecast(data,types{[types{:,1}] == dt,2}); 
    if swap, Y = swapbytes(Y); end
    Y = reshape(double(Y),dim); 
    if sl ~= 0 && isfinite(sl), Y = Y * sl + si; end
  catch e
    if gz, gis.close(); elseif ~isempty(fopen(fid)), fclose(fid); end
    rethrow(e); 
  end
end
% =========================================================================
function b = readBytes(ch,n)
%readBytes. Next n bytes of a (java) channel, e.g. of a gzip stream.
  b = zeros(1,0,'uint8'); 
  if n < 1, return; end
  bb = java.nio.ByteBuffer.allocate(n); 
  while bb.hasRemaining() && ch.read(bb) >= 0, end
  b  = reshape(typecast(bb.array(),'uint8'),1,[]); 
end
% =========================================================================
function Po = prepNii(Pi,opts,del)
  if opts.gzipi
    waschar = 0; 
    if ischar(Pi)
      waschar = 1; 
      Pi = cellstr(Pi); 
    end
    Po = Pi; 
    for fi = 1:numel(Pi)
      if strcmp( spm_file(Pi{fi},'ext'),'gz') 
        Po{fi} = spm_file(Pi{fi},'ext',''); 
      end
      if ~exist(Po{fi},'file') && exist(Pi{fi},'file')
        try
          gunzip(Pi{fi});
        catch
          % corrupted image (the sidecar stays, e.g. to convert it again)
          if exist(Pi{fi},'file'), delete(Pi{fi}); end
        end
        if exist('del','var') && del
          delete(Pi{fi}); 
        end
      end
    end
    if waschar
      Po = char(Po); 
    end
  else
    Po = Pi; 
  end 
end
% =========================================================================
function Po = prepNiigz(Pi,opts)
%prepNiigz. (Re)zip NIfTIs in case of opts.gzipi and remove the unzipped 
%  version. An existing .nii.gz is only updated if the .nii is newer (as 
%  gunzip keeps the file time, this is only the case for modified files). 
  if ~opts.gzipi, Po = Pi; return; end

  waschar = ischar(Pi); 
  if waschar, Pi = cellstr(Pi); end
  Po = Pi; 
  for fi = 1:numel(Pi)
    Pnii   = regexprep(Pi{fi},'\.gz$',''); 
    Po{fi} = [Pnii '.gz']; 
    if exist(Pnii,'file')
      Dnii = dir(Pnii); Dgz = dir(Po{fi}); 
      if isempty(Dgz) || Dnii.datenum > Dgz.datenum, gzipFile(Pnii); end
      delete(Pnii); 
    end
  end
  if waschar, Po = char(Po); end
end
% =========================================================================
function gzipFile(P)
%gzipFile. Compression of the file P to P.gz (P is kept), by the system 
%  gzip with the fastest compression (macOS/Linux, about 4 times faster than
%  the gzip of MATLAB for about 1% larger files, e.g. 16 instead of 65 s 
%  for a 0.9 GB rs-fMRI), otherwise (or if it fails) by the gzip of MATLAB.
  if isunix
    q  = @(s) ['''' strrep(s,'''','''\''''') '''']; % single-quoted shell argument
    st = system(sprintf('gzip -1 -c %s > %s',q(P),q([P '.gz']))); 
    if st == 0, return; end
  end
  gzip(P); 
end
% =========================================================================
function exportNii(P,Pout,opts)
%exportNii. Copy the NIfTI P of the database (zipped for opts.gzipi) to the
%  output file Pout (without extension) as .nii.gz (opts.gzipe) or .nii.
%  For the anonymization (opts.anonymize > 0), the acquisition time of the
%  description of the header is removed in P before (once, see 
%  removeAcquisitionTime).
  if isfield(opts,'anonymize') && opts.anonymize > 0, removeAcquisitionTime(P); end
  if opts.gzipe
    if opts.gzipi
      copyfile(P,[Pout '.nii.gz']); 
    else
      gzip(P); movefile([P '.gz'],[Pout '.nii.gz']); 
    end
  else
    if opts.gzipi
      gunzip(P); movefile(P(1:end-3),[Pout '.nii']); 
    else
      copyfile(P,[Pout '.nii']); 
    end
  end
end
% =========================================================================
function registerDerivative(P,Pout)
%registerDerivative. Index of the derivatives of the database: the result
%  P of the database (e.g. the label map c0 or a DTI map of a scan) as the
%  BIDS derivative Pout (without extension, i.e. with the extension of P, 
%  as for exportNii), also if it is not exported (opts.preprocessing = 1).
%  The renderings and the report tables use the database file for such 
%  derivatives (see derivativeIndex), i.e. the database is not part of the
%  BIDS directories. The index catDCM2BIDSderivatives.json is stored in the
%  session directory of the database (one process per subject, see 
%  cat_io_dcm2bids_parallel) with paths relative to the session directory
%  (P) and to the output directory (Pout) of the database.
  Pses    = fileparts(fileparts(P));                  % <db>/sub-*/ses-*/snr-*/P
  Poutdir = fileparts(fileparts(fileparts(Pses)));
  if ~strncmp(Pout,[Poutdir filesep],numel(Poutdir)+1), return; end
  ext  = regexp(P,'\.[^.\\/]+(\.gz)?$','match','once');
  Pidx = fullfile(Pses,'catDCM2BIDSderivatives.json');
  E    = struct('derivative',strrep([Pout(numel(Poutdir)+2:end) ext],filesep,'/'), ...
               'source',strrep(P(numel(Pses)+2:end),filesep,'/'));
  I = struct('Description','Database files of the BIDS derivatives (see registerDerivative in cat_io_dcm2bids_helper)', ...
    'derivatives',E([]));
  if exist(Pidx,'file')
    try
      J = jsondecode(readText(Pidx)); 
      if isfield(J,'derivatives') && isstruct(J.derivatives) && all(isfield(J.derivatives,{'derivative','source'}))
        I.derivatives = struct('derivative',{J.derivatives.derivative},'source',{J.derivatives.source})';
      end
    end
  end
  I.derivatives(strcmp({I.derivatives.derivative},E.derivative)) = [];
  I.derivatives(end+1) = E; 
  writeJSON(Pidx,I);
end
% =========================================================================
function D = derivativeIndex(Pmdbdir)
%derivativeIndex. Index of the derivatives of the database Pmdbdir (see
%  registerDerivative) as map of the BIDS derivative (full path) to its 
%  existing database file.
  D = containers.Map('KeyType','char','ValueType','char');
  if isempty(Pmdbdir) || ~exist(Pmdbdir,'dir'), return; end
  Poutdir = fileparts(regexprep(Pmdbdir,'[\\/]+$',''));
  Pidx    = dir(fullfile(Pmdbdir,'sub-*','ses-*','catDCM2BIDSderivatives.json'));
  for ii = 1:numel(Pidx)
    try
      J = jsondecode(readText(fullfile(Pidx(ii).folder,Pidx(ii).name)));
      for ei = 1:numel(J.derivatives)
        P = fullfile(Pidx(ii).folder,strrep(J.derivatives(ei).source,'/',filesep));
        if exist(P,'file'), D(fullfile(Poutdir,strrep(J.derivatives(ei).derivative,'/',filesep))) = P; end
      end
    end
  end
end
% =========================================================================
function writeTable(Ptable,T,format)
%writeTable. Write the cell table T (header and rows) as CSV or TSV file 
%  (format 'csv' or 'tsv', e.g. opts.tableformat) with UTF-8 encoding: 
%  numbers with up to 6 significant digits, numeric vectors separated by 
%  spaces, empty values and NaN as 'n/a' (as BIDS), CSV values with commas,
%  quotes, or line breaks in quotes (RFC 4180), and TSV values without 
%  tabs and line breaks. 
  if strcmpi(format,'tsv'), dl = sprintf('\t'); else, dl = ','; end
  L = cell(size(T,1),1); 
  for ri = 1:size(T,1)
    row = cell(1,size(T,2)); 
    for ci = 1:size(T,2)
      v = T{ri,ci}; 
      if isempty(v) || (isnumeric(v) && all(isnan(v(:))))
        v = 'n/a'; 
      elseif isnumeric(v) || islogical(v)
        v = strjoin(arrayfun(@(x) sprintf('%g',x),double(v(:))','UniformOutput',false),' '); 
      elseif iscell(v)
        v = strjoin(cellfun(@char,v(:)','UniformOutput',false),' '); 
      else
        v = char(v); 
      end
      if strcmp(dl,',')
        if any(ismember(v,[',"' newline char(13)])), v = ['"' strrep(v,'"','""') '"']; end
      else
        v = regexprep(v,'[\t\r\n]+',' '); 
      end
      row{ci} = v; 
    end
    L{ri} = strjoin(row,dl); 
  end
  if ~exist(fileparts(Ptable),'dir'), mkdir(fileparts(Ptable)); end
  writeText(Ptable,[strjoin(L,newline) newline]); 
end
function removeAcquisitionTime(P)
%removeAcquisitionTime. Remove the acquisition time of the description of
%  the header of the NIfTI P (.nii or .nii.gz, NIfTI-1 or NIfTI-2), e.g. 
%  "Time=133520.302;" of dcm2niix ("TE=37;Time=133520.302;phase=1;mb=8").
%  The header is only rewritten if the description contains a time. 
  if endsWith(P,'.gz')
    gz = java.util.zip.GZIPInputStream(java.io.FileInputStream(P)); 
    h  = zeros(1,540,'uint8'); n = 0; 
    while n < 540, x = gz.read(); if x < 0, break; end; n = n + 1; h(n) = x; end
    gz.close(); 
  else
    fid = fopen(P,'r'); h = fread(fid,540,'*uint8')'; fclose(fid); 
  end
  if numel(h) < 348, return; end
  switch double(typecast(h(1:4),'int32'))
    case 348, off = 148; % NIfTI-1
    case 540, off = 240; % NIfTI-2
    otherwise, return
  end
  d = char(h(off + (1:80))); d = d(1:find([d char(0)] == char(0),1) - 1); 
  if isempty(regexp(d,'Time=','once')), return; end
  d = regexprep(d,'Time=[0-9.]*;?',''); 
  Pu = regexprep(P,'\.gz$',''); 
  if endsWith(P,'.gz'), gunzip(P); end
  fid = fopen(Pu,'r+'); fseek(fid,off,'bof'); fwrite(fid,uint8([d zeros(1,80 - numel(d))]),'uint8'); fclose(fid); 
  if endsWith(P,'.gz'), gzipFile(Pu); delete(Pu); end
end
% =========================================================================
function entries = updateTable(Ptable,Thdr,Tnewrow,id,reimport)
%updateTable. Add a row to a table (csv) with the key in column id. 
%  An existing row with the same key is kept (reimport=0), replaced 
%  (reimport=1), or only its missing values (empty/NaN) are filled 
%  (reimport=2). A table with another header is converted to the current 
%  header (by the column names). entries is the number of rows. 
  if ~exist(Ptable,'file')
    if ~exist(fileparts(Ptable),'dir'), mkdir(fileparts(Ptable)); end
    cat_io_csv(Ptable,[Thdr;Tnewrow]);
    entries = 1; 
    return
  end
  Tfiles  = cat_io_csv(Ptable,'','',struct('convert2double',0));
  changed = 0; 

  % convert to the current header 
  if size(Tfiles,2) ~= numel(Thdr) || ~isequal(Tfiles(1,:),Thdr)
    T = repmat({''},size(Tfiles,1),numel(Thdr)); T(1,:) = Thdr; 
    [isc,ci] = ismember(Thdr,Tfiles(1,:)); 
    T(2:end,isc) = Tfiles(2:end,ci(isc)); 
    Tfiles = T; changed = 1; 
  end

  % row with the same key 
  keys = Tfiles(2:end,id); key = Tnewrow{id}; 
  if ~isempty(keys) && isnumeric(keys{1})
    if ischar(key), key = str2double(key); end
    hit = find(cellfun(@(k) isnumeric(k) && isequal(k,key), keys),1) + 1; 
  else
    hit = find(cellfun(@(k) (ischar(k) || isstring(k)) && strcmp(k,char(string(key))), keys),1) + 1; 
  end

  % convert fields to numbers if required
  if size(Tfiles,1) > 1
    for ci = 1:size(Tnewrow,2)
      if ischar(Tnewrow{ci}) && isnumeric(Tfiles{2,ci})
        Tnewrow{ci} = str2double(Tnewrow{ci}); 
      end
    end
  end

  missing = @(v) isempty(v) || (isnumeric(v) && all(isnan(v(:)))) || ...
    ((ischar(v) || isstring(v)) && any(strcmpi(char(v),{'nan','n/a'}))); 
  if isempty(hit)
    Tfiles(end+1,:) = Tnewrow; changed = 1; 
  elseif reimport == 1
    Tfiles(hit,:) = Tnewrow; changed = 1; 
  elseif reimport == 2
    for ci = 1:numel(Tnewrow)
      if missing(Tfiles{hit,ci}) && ~missing(Tnewrow{ci})
        Tfiles{hit,ci} = Tnewrow{ci}; changed = 1; 
      end
    end
  end
  if changed
    Tfiles(2:end,:) = sortrows(Tfiles(2:end,:),1);
    cat_io_csv(Ptable,Tfiles);
  end
  entries = size(Tfiles,1)-1;
end
% =========================================================================
function R = catRecords(R,S)
%catRecords. Concatenation of the record lists (structure arrays) R and S
%  with the fields of both, i.e. missing fields are NaN (e.g. records of a 
%  previous version of the batch without a new field).
  if isempty(R), R = S; return; end
  if isempty(S), return; end
  for fn = setdiff(fieldnames(S),fieldnames(R))', [R.(fn{1})] = deal(nan); end
  for fn = setdiff(fieldnames(R),fieldnames(S))', [S.(fn{1})] = deal(nan); end
  R = [R(:)' reshape(orderfields(S,R),1,[])];
end
% =========================================================================
function V = realignedVolumes(P,Prp,opts)
%realignedVolumes. Volumes of the (unzipped) 4D image P with the 
%  realignment of the QC (Prp, rp_*.txt with the rigid parameters of each 
%  volume relative to the first one, see cat_io_dcm2bids_qc), i.e. the
%  voxel-to-world mapping of volume t is spm_matrix(rp(t,:)) * V(1).mat, 
%  as the first volume (of each set) stays at its position (SPM realign), 
%  for the reslicing in memory (see realignedVolume) instead of a resliced
%  copy of the data. The first volume gets the mask of the voxels with 
%  data in all volumes (rmask, as the mask of spm_reslice) and the 
%  interpolation (rinterp: trilinear, B-spline of degree 4 for 
%  opts.hrrealign). Without realignment parameters, the volumes of P.
  V = spm_vol(P);
  if isempty(Prp) || ~exist(Prp,'file') || numel(V) < 2, return; end
  rp = load(Prp);
  if size(rp,1) ~= numel(V) || size(rp,2) < 6
    error('cat_io_dcm2bids_helper:realignedVolumes', ...
      'The realignment parameters (%d) do not fit the volumes (%d).',size(rp,1),numel(V));
  end
  M1 = V(1).mat; tiny = 5e-2; % (as spm_reslice)
  [x1,x2,x3] = ndgrid(1:V(1).dim(1),1:V(1).dim(2),1:V(1).dim(3));
  msk = true(V(1).dim);
  for t = 1:numel(V)
    V(t).mat = spm_matrix(rp(t,1:6)) * M1;
    M = V(t).mat \ M1;
    for d = 1:3
      x = M(d,1)*x1 + M(d,2)*x2 + M(d,3)*x3 + M(d,4);
      msk = msk & x >= 1 - tiny & x <= V(t).dim(d) + tiny;
    end
  end
  V(1).rmask = msk; V(1).rinterp = 1;
  if isfield(opts,'hrrealign') && opts.hrrealign, V(1).rinterp = 4; end
end
% =========================================================================
function Y = realignedVolume(V,t)
%realignedVolume. Volume t of the volumes V resliced into the grid of the
%  first volume with the realignment and mask of realignedVolumes (zeros 
%  outside the field of view of any volume), otherwise read (e.g. for 3D 
%  data or without realignment).
  if ~isfield(V,'rmask') || isempty(V(1).rmask) || t == 1 && isequal(V(1).mat,V(t).mat) && all(V(1).rmask(:))
    Y = spm_read_vols(V(t)); return
  end
  M = V(t).mat \ V(1).mat; d = V(1).dim;
  [x1,x2,x3] = ndgrid(1:d(1),1:d(2),1:d(3));
  X = {M(1,1)*x1 + M(1,2)*x2 + M(1,3)*x3 + M(1,4), M(2,1)*x1 + M(2,2)*x2 + M(2,3)*x3 + M(2,4), ...
       M(3,1)*x1 + M(3,2)*x2 + M(3,3)*x3 + M(3,4)};
  Vt = rmfield(V(t),intersect(fieldnames(V),{'rmask','rinterp'}));
  if V(1).rinterp > 1
    q = [V(1).rinterp*[1 1 1] 0 0 0];
    Y = spm_bsplins(spm_bsplinc(Vt,q),X{1},X{2},X{3},q);
  else
    Y = spm_sample_vol(Vt,X{1},X{2},X{3},1);
  end
  Y = reshape(Y,d); Y(~V(1).rmask) = 0;
end

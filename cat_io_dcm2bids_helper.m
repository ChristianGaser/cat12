function varargout = cat_io_dcm2bids_helper(action,varargin)
%cat_io_dcm2bids_helper. Shared helper functions of cat_io_dcm2bids.
%  Text and JSON files with UTF-8 encoding (readText, writeText, writeJSON),
%  paths (absPath, normPath, canonPath), NIfTI header information, single
%  volumes, the gzip handling of NIfTIs, and their export (niiHeaderInfo, 
%  niiVolume, prepNii, prepNiigz, gzipFile, exportNii), and tables 
%  (updateTable, writeTable, and lists of records with catRecords).
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
%    [vx,dim] = cat_io_dcm2bids_helper('niiHeaderInfo',P)
%    Y = cat_io_dcm2bids_helper('niiVolume',P,vol)
%    Po = cat_io_dcm2bids_helper('prepNii',Pi,opts,del)
%    Po = cat_io_dcm2bids_helper('prepNiigz',Pi,opts)
%    cat_io_dcm2bids_helper('gzipFile',P)
%    cat_io_dcm2bids_helper('exportNii',P,Pout,opts)
%    cat_io_dcm2bids_helper('removeAcquisitionTime',P)
%    entries = cat_io_dcm2bids_helper('updateTable',Ptable,Thdr,Tnewrow,id,reimport)
%    cat_io_dcm2bids_helper('writeTable',Ptable,T,format)
%    R = cat_io_dcm2bids_helper('catRecords',R,S)
%
%  See also cat_io_dcm2bids.

  switch action
    case {'readText','writeText','writeJSON','absPath','normPath','canonPath', ...
          'niiHeaderInfo','niiVolume','prepNii','prepNiigz','gzipFile','exportNii','removeAcquisitionTime', ...
          'updateTable','writeTable','catRecords'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_helper:unknownAction','Unknown action "%s".',action);
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
function [vx,dim] = niiHeaderInfo(P)
%niiHeaderInfo. Voxel size (x, y, z) and image dimensions (x, y, z, volumes) 
%  from the NIfTI header without reading the image data, i.e., for .nii.gz 
%  files only the header is decompressed. 
  vx = nan(1,3); dim = nan(1,4); 
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
  catch
    % fallback that reads the image 
    try
      evalc('V = spm_vol(P);'); 
      vx  = sqrt(sum(V(1).mat(1:3,1:3).^2)); 
      dim = [V(1).dim ones(1,3-numel(V(1).dim)) numel(V)]; 
    end
  end
  dim = reshape(dim,1,[]); vx = reshape(vx,1,[]); 
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


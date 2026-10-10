function varargout = cat_io_dcm2bids_helper(action,varargin)
%cat_io_dcm2bids_helper. Shared helper functions of cat_io_dcm2bids.
%  Text and JSON files with UTF-8 encoding (readText, writeText, writeJSON),
%  paths (absPath, normPath, canonPath), NIfTI header information and the 
%  gzip handling of NIfTIs (niiHeaderInfo, prepNii, prepNiigz), and tables
%  (updateTable).
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
%    Po = cat_io_dcm2bids_helper('prepNii',Pi,opts,del)
%    Po = cat_io_dcm2bids_helper('prepNiigz',Pi,opts)
%    entries = cat_io_dcm2bids_helper('updateTable',Ptable,Thdr,Tnewrow,id,reimport)
%
%  See also cat_io_dcm2bids.

  switch action
    case {'readText','writeText','writeJSON','absPath','normPath','canonPath', ...
          'niiHeaderInfo','prepNii','prepNiigz','updateTable'}
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
      if isempty(Dgz) || Dnii.datenum > Dgz.datenum, gzip(Pnii); end
      delete(Pnii); 
    end
  end
  if waschar, Po = char(Po); end
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

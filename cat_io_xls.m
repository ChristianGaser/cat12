function varargout = cat_io_xls(P,varargin)
%cat_io_xls. Write and read Excel files (xlsx) without Excel.
%  Writes cell tables as sheets of an xlsx file with formats and colors of
%  the cells (e.g. of QC ratings) and reads the sheets of xlsx files (e.g.
%  of Excel, LibreOffice, or Google Sheets). The files (Office Open XML)
%  are written and read directly, i.e. on all platforms and without Excel.
%  The action is given by the input: a table writes, otherwise it reads.
%
%  cat_io_xls(P,C[,opt])
%  [C,names] = cat_io_xls(P[,sheet])
%
%  P      xlsx file
%  C      cell table (rows x columns) of numbers, logicals, text, and empty
%         values (NaN and empty values are empty cells, vectors and cells
%         are joined to text separated by spaces), or a cell array of
%         tables {C1,C2,...} for several sheets
%  opt    options (one struct for all sheets or one per sheet)
%   .sheet       sheet name (default 'Sheet<i>')
%   .header      first row is the header (1, default: bold, frozen, and
%                with filter) or not (0)
%   .colorscale  colors of numeric columns (struct array):
%     .cols        columns (names of the header or numbers, columns that do
%                  not exist are ignored)
%     .map         colormap: 'marks' (default, CAT ratings 0.5-10.5 with
%                  the colors of the CAT reports), name of cat_io_colormaps
%                  or of a MATLAB colormap, or RGB matrix (n x 3, 0-1)
%     .lim         values of the first and last color ([min max], default
%                  range of the values, not used for 'marks')
%   .colormode   'fixed' cell colors (default) or 'scale' color scale
%                (conditional formatting with the colors at the min,
%                center, and max, i.e. Excel, LibreOffice, or Google Sheets
%                update the colors of changed values)
%
%  Text formats: text can start with TeX commands (as text of MATLAB
%  figures and \cellcolor of LaTeX tables) that define the format of its
%  cell, where a text that is a number after the commands is written as
%  number (e.g. '\cellcolor[marks]{2.3}2.3'):
%    \color{c}      text color
%    \cellcolor{c}  cell color
%    \bf, \it, \rm  bold, italic, and normal (followed by a space or a
%                   further command)
%  with the colors c as name (red, green, blue, cyan, magenta, yellow,
%  black, white, gray, darkGreen, orange, lightBlue), [rgb]{r g b} (0-1),
%  [RGB]{r g b} (0-255), [HTML]{RRGGBB}, or [marks]{m} (color of the CAT
%  rating m). Text on dark cell colors is white.
%
%  Reading gives the cells of the sheet (name or number, default 1) as
%  numbers (double), text (char), logicals, and empty cells ([]), where
%  numbers with a date format are given as text (yyyy-mm-dd, with HH:MM:SS
%  for times), and the names of all sheets. Formulas are given by their
%  last value. Formats and colors are not read.
%
%  Examples:
%    C = {'subject','status','SQR'; 'sub-01','\bf\color{red}nonconform',2.3; 'sub-02','conform',1.4};
%    cat_io_xls('qc.xlsx',C,struct('colorscale',struct('cols',{{'SQR'}},'map','marks')));
%    C = cat_io_xls('qc.xlsx');
%
%  See also cat_io_csv, cat_io_json, cat_io_colormaps.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  if nargin >= 2 && iscell(varargin{1})
    writeXLS(P,varargin{:});
  else
    [varargout{1:max(1,nargout)}] = readXLS(P,varargin{:});
  end
end
% =========================================================================
function writeXLS(P,C,opt)
%writeXLS. Write the cell table C (or the tables {C1,C2,...}) as sheets of
%  the xlsx file P with the options opt (see the help of this file).
  if nargin < 3 || isempty(opt), opt = struct(); end
  if isempty(C) || ~all(cellfun(@iscell,C(:))), C = {C}; end
  [pp,ff,ee] = fileparts(P); if isempty(ee), P = fullfile(pp,[ff '.xlsx']); end
  ns = numel(C);
  if isscalar(opt), opt = repmat(opt,1,ns); end
  if numel(opt) ~= ns, error('cat_io_xls:opt','One option struct or one per sheet (%d) expected.',ns); end

  % options of the sheets with unique sheet names
  def   = struct('sheet','','header',1,'colorscale',[],'colormode','fixed');
  O     = repmat(def,1,ns);
  names = cell(1,ns);
  for si = 1:ns
    for fn = fieldnames(opt(si))'
      if ~isfield(def,fn{1}), error('cat_io_xls:opt','Unknown option "%s".',fn{1}); end
      O(si).(fn{1}) = opt(si).(fn{1});
    end
    if ~ismember(O(si).colormode,{'fixed','scale'})
      error('cat_io_xls:colormode','Unknown colormode "%s" (fixed or scale).',O(si).colormode);
    end
    name = regexprep(char(O(si).sheet),'[\[\]:*?/\\]','_');
    if isempty(name), name = sprintf('Sheet%d',si); end
    name = name(1:min(31,end)); base = name; k = 1;
    while any(strcmpi(names(1:si-1),name)), k = k + 1; s = sprintf(' (%d)',k); name = [base(1:min(end,31-numel(s))) s]; end
    names{si} = name;
  end

  % cells of the sheets with their styles and the shared strings
  M = struct('font',containers.Map(),'fill',containers.Map(),'xf',containers.Map(), ...
    'fmt',containers.Map('KeyType','double','ValueType','any'));
  mapIndex(M.font,'<font><sz val="11"/><name val="Calibri"/><family val="2"/></font>');
  mapIndex(M.fill,'<fill><patternFill patternType="none"/></fill>');
  mapIndex(M.fill,'<fill><patternFill patternType="gray125"/></fill>');
  styleIndex(M,cellFormat());
  S = cell(1,ns); txt = cell(1,ns);
  for si = 1:ns
    [S{si},txt{si}] = sheetCells(C{si},O(si),M);
  end
  alltxt = cat(1,txt{:});
  alltxt = cellfun(@(s) s(1:min(end,32767)),alltxt,'UniformOutput',false);
  [sst,~,j] = unique(alltxt); n = 0;
  for si = 1:ns
    S{si}.sid(S{si}.kind == 2) = j(n + (1:numel(txt{si}))) - 1; n = n + numel(txt{si});
  end

  % package parts
  D = tempname; mkdir(D); cln = onCleanup(@() rmdir(D,'s'));
  ns0 = 'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"';
  rel = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships';
  hdr = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>';
  F   = {};
  ct  = '';
  for si = 1:ns
    F(end+1,:) = {sprintf('xl/worksheets/sheet%d.xml',si), [hdr sheetXML(S{si},si,ns0,rel)]}; %#ok<AGROW>
    ct = [ct sprintf(['<Override PartName="/xl/worksheets/sheet%d.xml" ContentType="application/' ...
      'vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'],si)]; %#ok<AGROW>
  end
  F(end+1,:) = {'[Content_Types].xml', [hdr '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">' ...
    '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>' ...
    '<Default Extension="xml" ContentType="application/xml"/>' ...
    '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>' ...
    '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>' ...
    '<Override PartName="/xl/sharedStrings.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sharedStrings+xml"/>' ...
    ct '</Types>']};
  F(end+1,:) = {'_rels/.rels', [hdr '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">' ...
    '<Relationship Id="rId1" Type="' rel '/officeDocument" Target="xl/workbook.xml"/></Relationships>']};
  wb = ''; dn = ''; wr = '';
  for si = 1:ns
    wb = [wb sprintf('<sheet name="%s" sheetId="%d" r:id="rId%d"/>',xmlText(names{si}),si,si)]; %#ok<AGROW>
    wr = [wr sprintf('<Relationship Id="rId%d" Type="%s/worksheet" Target="worksheets/sheet%d.xml"/>',si,rel,si)]; %#ok<AGROW>
    if ~isempty(S{si}.filter)
      dn = [dn sprintf('<definedName name="_xlnm._FilterDatabase" localSheetId="%d" hidden="1">''%s''!%s</definedName>', ...
        si-1,xmlText(strrep(names{si},'''','''''')),regexprep(S{si}.filter,'([A-Z]+)(\d+)','\$$1\$$2'))]; %#ok<AGROW>
    end
  end
  if ~isempty(dn), dn = ['<definedNames>' dn '</definedNames>']; end
  F(end+1,:) = {'xl/workbook.xml', [hdr '<workbook ' ns0 ' xmlns:r="' rel '"><bookViews><workbookView activeTab="0"/>' ...
    '</bookViews><sheets>' wb '</sheets>' dn '</workbook>']};
  F(end+1,:) = {'xl/_rels/workbook.xml.rels', [hdr '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">' ...
    wr sprintf('<Relationship Id="rId%d" Type="%s/styles" Target="styles.xml"/>',ns+1,rel) ...
    sprintf('<Relationship Id="rId%d" Type="%s/sharedStrings" Target="sharedStrings.xml"/>',ns+2,rel) '</Relationships>']};
  F(end+1,:) = {'xl/styles.xml', [hdr stylesXML(M,ns0)]};
  sst = xmlText(sst); si = '';
  if ~isempty(sst), si = sprintf('<si><t xml:space="preserve">%s</t></si>',sst{:}); end
  F(end+1,:) = {'xl/sharedStrings.xml', [hdr sprintf('<sst %s count="%d" uniqueCount="%d">',ns0,numel(alltxt),numel(sst)) ...
    si '</sst>']};
  for fi = 1:size(F,1)
    Pf = fullfile(D,F{fi,1}); if ~exist(fileparts(Pf),'dir'), mkdir(fileparts(Pf)); end
    fid = fopen(Pf,'w','n','UTF-8'); fprintf(fid,'%s',F{fi,2}); fclose(fid);
  end
  zip([P '.zip'],F([end-5 end-4 end-3 end-2 end-1 end 1:end-6],1),D); % ([Content_Types].xml first)
  movefile([P '.zip'],P,'f');
  clear cln
end
% =========================================================================
function [S,txt] = sheetCells(C,O,M)
%sheetCells. Cells of the table C with the options O of the sheet: kind
%  (0 empty, 1 number, 2 text), numbers num, styles xf (see styleIndex),
%  the conditional formats cf (colormode 'scale'), the filter range, and
%  the texts txt (in column order).
  [nr,nc] = size(C);
  if nr > 1048576 || nc > 16384, error('cat_io_xls:size','Table with %dx%d cells exceeds the size of a sheet.',nr,nc); end
  kind = zeros(nr,nc,'uint8'); num = nan(nr,nc); T = cell(nr,nc); xf = zeros(nr,nc);
  n   = cellfun('prodofsize',C);
  isn = cellfun('isclass',C,'double') & n == 1;
  if any(isn(:)), num(isn) = [C{isn}]; kind(isn & ~isnan(num)) = 1; end
  isc = cellfun('isclass',C,'char');
  if any(isc(:)), T(isc) = C(isc); kind(isc & n > 0) = 2; end
  for i = reshape(find(~isn & ~isc & n > 0),1,[])
    v = C{i};
    if (isnumeric(v) || islogical(v)) && isscalar(v)
      num(i) = double(v); kind(i) = ~isnan(num(i));
    elseif isnumeric(v) || islogical(v)
      T{i} = strjoin(arrayfun(@(x) sprintf('%g',x),double(v(:))','UniformOutput',false),' '); kind(i) = 2;
    elseif iscell(v)
      T{i} = strjoin(cellfun(@(x) char(string(x)),v(:)','UniformOutput',false),' '); kind(i) = 2;
    else
      T{i} = char(string(v)); kind(i) = 2;
    end
  end
  ii = find(kind == 1 & isinf(num)); kind(ii) = 2; T(ii) = arrayfun(@num2str,num(ii),'UniformOutput',false);
  if any(kind(:) == 2), T(kind == 2) = cellfun(@(s) reshape(s,1,[]),T(kind == 2),'UniformOutput',false); end

  % TeX commands of text cells (numbers if the text is a number)
  itex = find(kind == 2); itex = itex(strncmp(T(itex),'\',1)); ftex = cell(size(itex));
  for j = 1:numel(itex)
    [T{itex(j)},ftex{j},x] = texFormat(T{itex(j)});
    if ~isnan(x), num(itex(j)) = x; kind(itex(j)) = 1; elseif isempty(T{itex(j)}), kind(itex(j)) = 0; end
  end

  % header (bold)
  hdr = O.header && nr > 0 && nc > 0; names = {};
  if hdr
    xf(1,:) = styleIndex(M,cellFormat('bf',1));
    names = T(1,:); names(kind(1,:) ~= 2) = {''};
  end

  % colors of the numeric columns
  cf = {}; rows = (1+hdr:nr)'; L = colNames(nc);
  R  = O.colorscale; if isempty(R), R = struct('cols',{}); end
  for ri = 1:numel(R)
    if ~isfield(R(ri),'cols'), continue; end
    if ~isfield(R(ri),'map') || isempty(R(ri).map), R(ri).map = 'marks'; end
    if ~isfield(R(ri),'lim'), R(ri).lim = []; end
    if iscellstr(R(ri).cols) || ischar(R(ri).cols) %#ok<ISCLSTR>
      cols = find(ismember(names,cellstr(R(ri).cols)));
    else
      cols = R(ri).cols(R(ri).cols >= 1 & R(ri).cols <= nc);
    end
    if isempty(cols) || isempty(rows), continue; end
    v = num(rows,cols); v = v(kind(rows,cols) == 1);
    [cmap,bin,pts,nf] = ruleMap(R(ri),v);
    for k = reshape(cols,1,[])
      ok = rows(kind(rows,k) == 1);
      if isempty(ok), continue; end
      if strcmp(O.colormode,'scale')
        xf(ok,k) = styleIndex(M,cellFormat('nf',nf));
      else
        [ub,~,bi] = unique(bin(num(ok,k)));
        for u = 1:numel(ub)
          xf(ok(bi == u),k) = styleIndex(M,cellFormat('bc',rgbHex(cmap(ub(u),:)),'nf',nf));
        end
      end
    end
    if strcmp(O.colormode,'scale')
      sq = strjoin(cellfun(@(l) sprintf('%s%d:%s%d',l,rows(1),l,nr),L(cols),'UniformOutput',false),' ');
      cf{end+1} = sprintf(['<conditionalFormatting sqref="%s"><cfRule type="colorScale" priority="%d"><colorScale>' ...
        '%s%s</colorScale></cfRule></conditionalFormatting>'],sq,numel(cf)+1, ...
        sprintf('<cfvo type="num" val="%.15g"/>',pts), ...
        sprintf('<color rgb="FF%s"/>',rgbHex(cmap(bin(pts(1)),:)),rgbHex(cmap(bin(pts(2)),:)),rgbHex(cmap(bin(pts(3)),:)))); %#ok<AGROW>
    end
  end

  % formats of the TeX commands (on the formats of the header and colors)
  for j = 1:numel(itex)
    f = ftex{j};
    if ~any(structfun(@(x) ~isempty(x),f)), continue; end
    xf(itex(j)) = styleIndex(M,mergeFormat(M.fmt(xf(itex(j))),f));
  end

  % column widths (characters)
  w = zeros(1,nc);
  for k = 1:nc
    t  = cellfun('prodofsize',T(kind(:,k) == 2,k));
    x  = num(kind(:,k) == 1,k);
    if ~isempty(x), x = numel(sprintf('%.6g',max(abs(x)))) + any(x < 0) + 1; end
    w(k) = max([6; t(:)*1.1; x(:)]);
    if hdr && kind(1,k) == 2, w(k) = max(w(k),numel(T{1,k})*1.2 + 3); end
  end
  w = min(60,w + 1);

  txt = reshape(T(kind == 2),[],1);
  S = struct('kind',kind,'num',num,'sid',zeros(nr,nc),'xf',xf,'width',w,'hdr',hdr,'cf',{cf},'filter','');
  if hdr, S.filter = sprintf('A1:%s%d',L{max(1,nc)},nr); end
end
% =========================================================================
function x = sheetXML(S,si,ns0,rel)
%sheetXML. XML of the worksheet of the cells S (see sheetCells).
  [nr,nc] = size(S.kind); L = colNames(max(1,nc));
  X = cell(nr,nc);
  for k = 1:nc
    for kd = 0:2
      if kd == 0, r = find(S.kind(:,k) == 0 & S.xf(:,k) > 0); else, r = find(S.kind(:,k) == kd); end
      if isempty(r), continue; end
      switch kd
        case 0, s = sprintf(['<c r="' L{k} '%d" s="%d"/>\n'],[r'; S.xf(r,k)']);
        case 1, s = sprintf(['<c r="' L{k} '%d" s="%d"><v>%.15g</v></c>\n'],[r'; S.xf(r,k)'; S.num(r,k)']);
        case 2, s = sprintf(['<c r="' L{k} '%d" s="%d" t="s"><v>%d</v></c>\n'],[r'; S.xf(r,k)'; S.sid(r,k)']);
      end
      s = regexp(s(1:end-1),'\n','split');
      X(r,k) = s(:);
    end
  end
  R = cell(nr,1);
  for r = 1:nr
    c = [X{r,:}];
    if ~isempty(c), R{r} = sprintf('<row r="%d">%s</row>',r,c); end
  end
  if nr*nc > 0, dim = sprintf('A1:%s%d',L{nc},nr); else, dim = 'A1'; end
  sv = sprintf('<sheetView workbookViewId="0"%s>',repmat(' tabSelected="1"',1,si == 1));
  if S.hdr && nr > 1
    sv = [sv '<pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/><selection pane="bottomLeft"/>'];
  end
  cols = sprintf('<col min="%d" max="%d" width="%0.1f" customWidth="1"/>',[1:nc; 1:nc; S.width]);
  if ~isempty(cols), cols = ['<cols>' cols '</cols>']; end
  flt = ''; if ~isempty(S.filter), flt = sprintf('<autoFilter ref="%s"/>',S.filter); end
  x = ['<worksheet ' ns0 ' xmlns:r="' rel '"><dimension ref="' dim '"/><sheetViews>' sv '</sheetView></sheetViews>' ...
    '<sheetFormatPr defaultRowHeight="15"/>' cols '<sheetData>' [R{:}] '</sheetData>' flt [S.cf{:}] ...
    '<pageMargins left="0.7" right="0.7" top="0.75" bottom="0.75" header="0.3" footer="0.3"/></worksheet>'];
end
% =========================================================================
function x = stylesXML(M,ns0)
%stylesXML. XML of the styles (fonts, fills, and cell formats of M).
  part = @(m) mapList(m);
  fo = part(M.font); fi = part(M.fill); xf = part(M.xf);
  x = sprintf(['<styleSheet %s><fonts count="%d">%s</fonts><fills count="%d">%s</fills>' ...
    '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>' ...
    '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>' ...
    '<cellXfs count="%d">%s</cellXfs><cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>' ...
    '<dxfs count="0"/><tableStyles count="0"/></styleSheet>'], ...
    ns0,numel(fo),[fo{:}],numel(fi),[fi{:}],numel(xf),[xf{:}]);
end
% =========================================================================
function L = mapList(m)
%mapList. Keys of the map m (key -> index 0..n-1) in the order of the index.
  k = keys(m); v = cell2mat(values(m));
  L = cell(1,numel(k)); L(v+1) = k;
end
% =========================================================================
function id = mapIndex(m,key)
%mapIndex. Index of the key in the map m (added if new, 0-based).
  if isKey(m,key), id = m(key); else, id = double(m.Count); m(key) = id; end %#ok<NASGU> (handle)
end
% =========================================================================
function f = cellFormat(varargin)
%cellFormat. Format of a cell: bold bf, italic it, text color fc, cell color
%  bc (RGB hex), and number format nf (0 general, 2 0.00) with the
%  defaults of the normal cell and the fields given as name-value pairs.
  f = struct('bf',0,'it',0,'fc','','bc','','nf',0);
  for i = 1:2:numel(varargin), f.(varargin{i}) = varargin{i+1}; end
end
% =========================================================================
function f = mergeFormat(f0,f)
%mergeFormat. Format f with the undefined (empty) fields of the format f0.
  for fn = fieldnames(f)'
    if isempty(f.(fn{1})), f.(fn{1}) = f0.(fn{1}); end
  end
end
% =========================================================================
function id = styleIndex(M,f)
%styleIndex. Index of the cell format f (see cellFormat) in the styles M
%  (fonts, fills, and cell formats, added if new), where the text on dark
%  cell colors is white.
  fc = f.fc;
  if isempty(fc) && ~isempty(f.bc)
    c = sscanf(f.bc,'%2x')' / 255; c(c <= 0.04045) = c(c <= 0.04045) / 12.92;
    c(c > 0.04045) = ((c(c > 0.04045) + 0.055) / 1.055).^2.4;
    if [0.2126 0.7152 0.0722] * c(:) < 0.179, fc = 'FFFFFF'; end % (relative luminance, WCAG)
  end
  if isempty(fc), col = ''; else, col = sprintf('<color rgb="FF%s"/>',fc); end
  fo = mapIndex(M.font,sprintf('<font>%s%s<sz val="11"/>%s<name val="Calibri"/><family val="2"/></font>', ...
    repmat('<b/>',1,f.bf > 0),repmat('<i/>',1,f.it > 0),col));
  fi = 0;
  if ~isempty(f.bc)
    fi = mapIndex(M.fill,sprintf(['<fill><patternFill patternType="solid"><fgColor rgb="FF%s"/>' ...
      '<bgColor indexed="64"/></patternFill></fill>'],f.bc));
  end
  id = mapIndex(M.xf,sprintf('<xf numFmtId="%d" fontId="%d" fillId="%d" borderId="0" xfId="0"%s%s%s/>',f.nf,fo,fi, ...
    repmat(' applyNumberFormat="1"',1,f.nf > 0),repmat(' applyFont="1"',1,fo > 0),repmat(' applyFill="1"',1,fi > 0)));
  if ~isKey(M.fmt,id), M.fmt(id) = f; end
end
% =========================================================================
function [s,f,x] = texFormat(s)
%texFormat. Leading TeX commands of the text s (see the help of this file)
%  as cell format f (fields of cellFormat, empty if not defined), the text s
%  without the commands, and its number x (NaN if it is no number).
  f = struct('bf',[],'it',[],'fc',[],'bc',[],'nf',[]); x = nan; s0 = s;
  while ~isempty(s) && s(1) == '\'
    [t,e] = regexp(s,'^\\(?<cmd>cellcolor|color)(?<model>\[\w+\])?\{(?<val>[^}]*)\} ?','names','end','once');
    if ~isempty(t)
      c = texColor(regexprep(char(t.model),'[\[\]]',''),t.val);
      if strcmp(t.cmd,'color'), f.fc = c; else, f.bc = c; end
      s = s(e+1:end); continue
    end
    [t,e] = regexp(s,'^\\(bf|it|rm)(?![a-zA-Z]) ?','tokens','end','once');
    if isempty(t), break; end
    switch t{1}
      case 'bf', f.bf = 1;
      case 'it', f.it = 1;
      case 'rm', f.bf = 0; f.it = 0;
    end
    s = s(e+1:end);
  end
  if ~strcmp(s,s0) && ~isempty(regexp(s,'^\s*[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?\s*$','once'))
    x = str2double(s);
  end
end
% =========================================================================
function c = texColor(model,v)
%texColor. RGB hex code of the TeX color v of the color model ('' names,
%  rgb, RGB, HTML, or marks, see the help of this file).
  if isempty(model), model = 'name'; end
  switch model
    case 'name'
      N = {'black','000000'; 'white','FFFFFF'; 'red','FF0000'; 'green','00FF00'; 'blue','0000FF'; ...
           'cyan','00FFFF'; 'magenta','FF00FF'; 'yellow','FFFF00'; 'gray','808080'; ...
           'darkGreen','008000'; 'orange','FF8000'; 'lightBlue','ADD8E6'};
      i = find(strcmpi(N(:,1),strtrim(v)),1);
      if isempty(i), error('cat_io_xls:color','Unknown color "%s".',v); end
      c = N{i,2};
    case 'rgb',   c = rgbHex(sscanf(strrep(v,',',' '),'%f')');
    case 'RGB',   c = rgbHex(sscanf(strrep(v,',',' '),'%f')' / 255);
    case 'HTML'
      c = upper(regexprep(strtrim(v),'^#',''));
      if isempty(regexp(c,'^[0-9A-F]{6}$','once')), error('cat_io_xls:color','Unknown HTML color "%s".',v); end
    case 'marks', c = rgbHex(markColor(str2double(v)));
    otherwise,    error('cat_io_xls:color','Unknown color model "%s".',model);
  end
end
% =========================================================================
function [cmap,bin,pts,nf] = ruleMap(R,v)
%ruleMap. Colormap cmap of the colorscale rule R with the function bin of
%  the colormap index of values, the values pts of the colors of the color
%  scale (min, center, max), and the number format nf (2 for ratings),
%  where the limits are the range of the values v by default.
  if ischar(R.map) && strcmp(R.map,'marks')
    cmap = markColor([]); n = size(cmap,1);
    bin  = @(x) min(n-3,max(1,floor(x/9.5*n)));     % as cat_io_dcm2bids_reportfig
    pts  = [1 3.5 6]; nf = 2;
    return
  end
  if ischar(R.map)
    try
      cmap = cat_io_colormaps(R.map,64);
    catch
      cmap = feval(R.map,64);
    end
  else
    cmap = R.map;
  end
  lim = R.lim; v = v(isfinite(v));
  if isempty(lim) && isempty(v), lim = [0 1]; elseif isempty(lim), lim = [min(v) max(v)]; end
  if lim(2) <= lim(1), lim(2) = lim(1) + 1; end
  n   = size(cmap,1);
  bin = @(x) 1 + round((min(max(x,lim(1)),lim(2)) - lim(1)) / (lim(2) - lim(1)) * (n-1));
  pts = [lim(1) mean(lim) lim(2)]; nf = 0;
end
% =========================================================================
function c = markColor(m)
%markColor. Colors of the CAT ratings m (marks) as in the CAT reports (see
%  cat_io_dcm2bids_reportfig), the colormap for empty m.
  c = min(1,max(0,cat_io_colormaps('marks+',40)));
  if ~isempty(m), c = c(min(size(c,1)-3,max(1,floor(m(:)/9.5*size(c,1)))),:); end
end
% =========================================================================
function h = rgbHex(c)
%rgbHex. RGB hex code (RRGGBB) of the color c (0-1).
  if numel(c) < 3 || any(~isfinite(c(1:3))), error('cat_io_xls:color','Unknown RGB color.'); end
  h = sprintf('%02X%02X%02X',round(255*min(1,max(0,c(1:3)))));
end
% =========================================================================
function L = colNames(n)
%colNames. Column names A, B, ..., Z, AA, ... of the columns 1:n.
  L = cell(1,n);
  for k = 1:n
    s = ''; m = k;
    while m > 0, d = mod(m-1,26); s = [char(65+d) s]; m = (m-1-d)/26; end %#ok<AGROW>
    L{k} = s;
  end
end
% =========================================================================
function s = xmlText(s)
%xmlText. Text s (or cellstr) as XML text (escaped, without the control
%  characters that are not allowed in XML).
  if iscell(s), s = cellfun(@xmlText,s,'UniformOutput',false); return; end
  s(s < 32 & s ~= 9 & s ~= 10 & s ~= 13) = [];
  s = strrep(strrep(strrep(strrep(s,'&','&amp;'),'<','&lt;'),'>','&gt;'),'"','&quot;');
end
% =========================================================================
function [C,names] = readXLS(P,sheet)
%readXLS. Cells of the sheet (name or number, default 1) of the xlsx file P
%  and the names of all sheets (see the help of this file).
  if nargin < 2 || isempty(sheet), sheet = 1; end
  if ~exist(P,'file'), error('cat_io_xls:noFile','File "%s" does not exist.',P); end
  D = tempname; mkdir(D); cln = onCleanup(@() rmdir(D,'s'));
  try
    unzip(P,D);
  catch
    error('cat_io_xls:noXLSX','File "%s" is no xlsx file.',P);
  end

  % sheets of the workbook and their files
  wb    = readPart(fullfile(D,'xl','workbook.xml'));
  sh    = regexp(wb,'<sheet\s[^>]*>','match');
  names = xmlUnescape(regexp(sh,'(?<=\sname=")[^"]*','match','once'));
  rid   = regexp(sh,'\s\w+:id="([^"]*)"','tokens','once'); rid = cellfun(@(t) char([t{:}]),rid,'UniformOutput',false);
  if ischar(sheet), si = find(strcmp(names,sheet),1); else, si = sheet; end
  if isempty(si) || si < 1 || si > numel(sh)
    error('cat_io_xls:noSheet','Sheet "%s" does not exist in "%s".',num2str(sheet),P);
  end
  rels = regexp(readPart(fullfile(D,'xl','_rels','workbook.xml.rels')),'<Relationship\s[^>]*>','match');
  part = @(attr,val) partFile(D,rels,attr,val);
  Psh  = part('Id',rid{si});
  Psst = part('Type','/sharedStrings');
  Psty = part('Type','/styles');
  d1904 = ~isempty(regexp(wb,'date1904="(1|true)"','once'));

  % shared strings (text of all runs, without phonetic runs)
  sst = {};
  if ~isempty(Psst)
    t   = regexp(readPart(Psst),'<si[^>]*>(.*?)</si>','tokens');
    sst = cellfun(@(t) char([t{:}]),t,'UniformOutput',false);
    sst = xmlUnescape(regexprep(sst,{'<rPh[\s>].*?</rPh>','<rPr>.*?</rPr>','<[^>]*>'},{'','',''}));
  end

  % date formats of the cell formats
  isdate = false(0);
  if ~isempty(Psty)
    st  = readPart(Psty);
    fmt = regexp(st,'<numFmt\s[^>]*>','match');
    fid = str2double(regexp(fmt,'(?<=numFmtId=")\d+','match','once'));
    fcd = xmlUnescape(regexp(fmt,'(?<=formatCode=")[^"]*','match','once'));
    fcd = regexprep(fcd,{'\[[^\]]*\]','"[^"]*"','\\.'},{'','',''});
    dfm = fid(~cellfun('isempty',regexpi(fcd,'[dmyhs]','once')));
    xfs = regexp(regexp(st,'<cellXfs[^>]*>(.*?)</cellXfs>','tokens','once'),'<xf\s[^>]*>','match');
    if ~isempty(xfs)
      xid    = str2double(regexp(xfs{1},'(?<=numFmtId=")\d+','match','once'));
      isdate = ismember(xid,[14:22 27:36 45:47 50:58 dfm]);
    end
  end

  % cells of the sheet
  x    = readPart(Psh);
  rows = regexp(x,'<row(?<a>[^>]*)>(?<b>.*?)</row>','names');
  rnum = str2double(regexp({rows.a},'(?<=\sr=")\d+','match','once'));
  for ri = find(isnan(rnum)), if ri == 1, rnum(ri) = 1; else, rnum(ri) = rnum(ri-1) + 1; end, end
  cells = regexp({rows.b},'<c(?<a>[^>]*)>(?<b>.*?)</c>','names');
  nc  = cellfun(@numel,cells);
  if isempty(cells) || sum(nc) == 0, C = cell(0,0); return; end
  A   = [cells{nc > 0}];
  row = repelem(rnum(:),nc(:));
  ref = regexp({A.a},'(?<=\sr=")[A-Z]+','match','once');
  col = zeros(numel(A),1); [uref,~,ui] = unique(ref);
  for u = 1:numel(uref)
    if ~isempty(uref{u}), col(ui == u) = (26.^(numel(uref{u})-1:-1:0)) * (double(uref{u}(:)) - 64); end
  end
  if any(col == 0) % (cells without reference: next column of the row)
    n0 = 0;
    for i = 1:numel(A)
      if i == 1 || row(i) ~= row(i-1), n0 = 0; end
      if col(i) == 0, col(i) = n0 + 1; end
      n0 = col(i);
    end
  end
  typ = regexp({A.a},'(?<=\st=")\w+','match','once');
  sty = str2double(regexp({A.a},'(?<=\ss=")\d+','match','once')); sty(isnan(sty)) = 0;
  val = regexp({A.b},'<v>([^<]*)</v>','tokens','once');
  val(cellfun('isempty',val)) = {{''}}; val = cellfun(@(t) t{1},val,'UniformOutput',false);

  C   = cell(max(row),max(col));
  idx = sub2ind(size(C),row(:),col(:));
  num = str2double(val);
  ts  = strcmp(typ,'s') & ~cellfun('isempty',val);
  C(idx(ts)) = sst(num(ts) + 1);
  tn  = (cellfun('isempty',typ) | strcmp(typ,'n')) & ~isnan(num);
  dt  = tn & sty < numel(isdate); dt(dt) = isdate(sty(dt) + 1);
  C(idx(tn & ~dt)) = num2cell(num(tn & ~dt));
  for i = find(dt), C{idx(i)} = dateText(num(i),d1904); end
  for i = find(~ts & ~tn & ~cellfun('isempty',typ))
    switch typ{i}
      case 'inlineStr'
        C{idx(i)} = xmlUnescape(regexprep(A(i).b,{'<rPh[\s>].*?</rPh>','<rPr>.*?</rPr>','<[^>]*>'},{'','',''}));
      case {'str','e','d'}
        C{idx(i)} = xmlUnescape(val{i});
      case 'b'
        if ~isempty(val{i}), C{idx(i)} = num(i) ~= 0; end
    end
  end
  clear cln
end
% =========================================================================
function P = partFile(D,rels,attr,val)
%partFile. File of the relationship of the workbook with the attribute
%  attr (Id: equal, Type: ending with val), empty if it does not exist.
  if strcmp(attr,'Id')
    i = find(~cellfun('isempty',regexp(rels,['\sId="' regexptranslate('escape',val) '"'],'once')),1);
  else
    i = find(~cellfun('isempty',regexp(rels,['\sType="[^"]*' regexptranslate('escape',val) '"'],'once')),1);
  end
  P = '';
  if isempty(i), return; end
  t = regexp(rels{i},'(?<=\sTarget=")[^"]*','match','once');
  if strncmp(t,'/',1), P = fullfile(D,t(2:end)); else, P = fullfile(D,'xl',t); end
  if ~exist(P,'file'), P = ''; end
end
% =========================================================================
function x = readPart(P)
%readPart. XML of the part P (UTF-8) without the namespace prefixes of the
%  elements (e.g. <x:c> of some software) and with the empty rows, cells,
%  and strings as start and end tag (e.g. <c r="A1"/> as <c r="A1"></c>).
  fid = fopen(P,'r','n','UTF-8'); x = fread(fid,[1 inf],'*char'); fclose(fid);
  x = regexprep(x,'<(/?)\w+:','<$1');
  x = regexprep(x,'<(row|c|si)(\s[^>]*?|)/>','<$1$2></$1>');
end
% =========================================================================
function s = xmlUnescape(s)
%xmlUnescape. Text of the XML text s (or cellstr) with the entities and the
%  escaped characters of Excel (_xHHHH_).
  s = regexprep(s,{'&lt;','&gt;','&quot;','&apos;'},{'<','>','"',''''});
  s = regexprep(s,'&#x([0-9a-fA-F]+);','${char(hex2dec($1))}');
  s = regexprep(s,'&#(\d+);','${char(str2double($1))}');
  s = regexprep(s,'_x([0-9a-fA-F]{4})_','${char(hex2dec($1))}');
  s = strrep(s,'&amp;','&');
end
% =========================================================================
function s = dateText(x,d1904)
%dateText. Date text of the Excel date number x (yyyy-mm-dd, with HH:MM:SS
%  for times, only the time for x < 1).
  x = round(x * 86400) / 86400;
  if d1904, d = x + 695422; else, d = x + 693960; end
  if x < 1 && x > 0,   s = datestr(x,'HH:MM:SS'); %#ok<DATST>
  elseif mod(x,1) > 0, s = datestr(d,'yyyy-mm-dd HH:MM:SS'); %#ok<DATST>
  else,                s = datestr(d,'yyyy-mm-dd'); %#ok<DATST>
  end
end

function [filesfound,numberfound] = vbm_findfiles(varargin)
% ______________________________________________________________________
% vbm_findfiles  - Linux/UNIX-like find command/function
%
% returns a list in which all files are listed that match one of the
% patterns. This function was copyed form the 
%
% FORMAT:
% [files,number] = vbm_findfiles(startfolder, patterns, opts, ....)
%
% - startfolder     <char array> folder where to start search
% - patterns        <cell array <char array>> file patterns (shell-like)
%   (alt.)          <char array>  : single pattern
% - opt             <struct array>: optional parameters in the form of
%   .cellstr        <double array>: if set and !0 -> return as cellstr
%   .chararr        <double array>: if set and !0 -> return as char array
%   .depth          <double array>: sets both minimum and maximum depth
%   .dirs           <double array>: if set and !0 -> dirs instead of files
%   .maxdepth       <double array>: maximum depth where to search
%   .mindepth       <double array>: minimum depth where to search
%   .maxage         <double array>: seconds file must have changed in
%   .minage         <double array>: seconds file must not have changed in
%   .oneperdir      <double array>: if set and !0 -> only first match
%   .relative       <... array>   : if set -> !0 -> './' // char array
%
% when used as a command, multiple opts can be given as multiple arguments,
% seperated by spaces (' '):
%
% vbm_findfiles /search/folder *mypattern*.txt depth=3 oneperdir=1 relative=./
%
% when used in functional context, a second return value, the number
% of matching files, can be obtained:
%
% [files,number] = vbm_findfiles('/where','*.txt');
%
% NOTE: the minage/maxage feature only fully works when the system
% returns English-style month in calls to dir. i.e. under Linux, set the
% LANG environmental setting to 'en_US' before starting up MatLab

%
% TODO: make the returning object (cell) a persistent value to increase
%       speed by not having to return it from recursive calls !!!
%
% ______________________________________________________________________
% This function was part of the ... toolbox that I cannot find anymore.
% 
% ______________________________________________________________________
% $Id$ 

%#ok<*EFIND>

% - sanity checks: startfolder
  fsep = filesep;
  if nargin < 2
      dispdebug('vbm_findfiles: at least startfolder and patterns must be provided!',4);
      if nargout < 1, help(mfilename); return; end
      filesfound=cell(0);
      return;
  end
  startfolder=varargin{1};
  if ischar(startfolder) && ~isempty(startfolder)
      ispattern=find(startfolder=='?'|startfolder=='*');
      if ~isempty(ispattern), [filesfound,numberfound] = vbm_findfiles({varargin{1}},varargin{2:end}); return; end
  elseif iscell(startfolder) && ~isempty(startfolder)
      nstartfolder=cell(0);
      for nelem=1:prod(size(startfolder))
          if ~ischar(startfolder{nelem})
              error('ERROR: parameter startfolder must be of type char or cell of chars');
          end
          ispattern=find(startfolder{nelem}=='?'|startfolder{nelem}=='*');
          if isempty(ispattern)
              nstartfolder{end+1}=startfolder{nelem};
          else
              [pathparts,cparts]=splittocell(startfolder{nelem},'/\',1);
              for cpart=1:cparts
                  if ~isempty(find(pathparts{cpart}=='?'|pathparts{cpart}=='*')), break; end
              end
              if cpart==1
                  [pfolders,npfolders]=vbm_findfiles('.',pathparts{1},struct('dirs',1,'depth',1));
              else
                  [pfolders,npfolders]=vbm_findfiles(gluetostring({pathparts{1:(cpart-1)}},fsep),pathparts{cpart},struct('dirs',1,'depth',1));
              end
              if cpart < cparts
                  for ppart=1:npfolders
                      pfolders{ppart} = [pfolders{ppart} fsep gluetostring({pathparts{(cpart+1):end}},fsep)];
                  end
              end
              nstartfolder = [nstartfolder,pfolders];
          end
      end
      filesfound=cell(0);
      for nelem=1:length(nstartfolder)
          if ~isempty(find(nstartfolder{nelem}=='?'|nstartfolder{nelem}=='*'))
              filesfound = [filesfound(1:end) , vbm_findfiles({nstartfolder{nelem}},varargin{2:end})];
          elseif exist(nstartfolder{nelem},'dir') == 7
              filesfound = [filesfound(1:end) , vbm_findfiles(nstartfolder{nelem},varargin{2:end})];
          end
      end
      numberfound=length(filesfound);
      return;
  else
      error('ERROR: parameter startfolder must be of type char or cell of chars!');
  end
  if startfolder(size(startfolder,2)) ~= fsep, startfolder = [startfolder, fsep]; end
  if exist(startfolder,'dir') ~= 7, error(['ERROR: ' startfolder ' not found or no folder!']), end

  % - sanity checks: patterns
  patterns=varargin{2};
  if ~iscell(patterns)
      if ischar(patterns)
          patterns = { patterns };
      else
          error('ERROR: patterns must either be a single string or a cell array!');
      end
  end
  for count=1:size(patterns,2),
      if ~ischar(patterns{count})
          error('ERROR: all patterns cell contents must be char arrays!');
      end
  end

  % - option argument parsing
  if nargin < 3
      opt.dirs = 0;
      opt.maxdepth = 0;
      opt.mindepth = 0;
      opt.maxage = -1;
      opt.minage = -1;
      opt.oneperdir = 0;
      opt.relative = 0;
      opt.return = 'cellarr';
      opt.rfolder = startfolder;
  else
      opt=varargin{3};
      if ~isstruct(opt)
          clear opt;
          opt.dirs = 0;
          opt.maxdepth = 0;
          opt.mindepth = 0;
          opt.maxage = -1;
          opt.minage = -1;
          opt.oneperdir = 0;
          opt.relative = 0;
          opt.return = 'cellarr';
          opt.rfolder = startfolder;
          for acount=3:nargin
              if ischar(varargin{acount})
                  argnv=splittocell(varargin{acount},'=',1);
                  oname=argnv{1};
                  if size(argnv,2) > 1, oval=argnv{2}; else, oval=''; end
                  switch lower(oname)
                  case 'cellstr'
                      oval=str2double(oval);
                      if oval ~= 0, opt.return = 'cellstr'; end
                  case 'chararr'
                      oval=str2double(oval);
                      if oval ~= 0, opt.return = 'chararr'; end
                  case 'depth'
                      if str2double(oval) >= 0
                          opt.maxdepth=str2double(oval);
                          opt.mindepth=str2double(oval);
                      else
                          opt.maxdepth=0;
                          opt.mindepth=0;
                      end
                  case 'dirs'
                      oval=str2double(oval);
                      if oval == 0, opt.dirs = 0; else, opt.dirs = 1; end
                  case 'maxdepth'
                      if str2double(oval) >= 0, opt.maxdepth=str2double(oval); else, opt.maxdepth=0; end
                  case 'mindepth'
                      if str2double(oval) >= 0, opt.mindepth=str2double(oval); else, opt.mindepth=0; end
                  case 'maxage'
                      if str2double(oval) >= 0, opt.maxage=str2double(oval); else, opt.maxage=-1; end
                  case 'minage'
                      if str2double(oval) >= 0, opt.minage=str2double(oval); else, opt.minage=-1; end
                  case 'oneperdir'
                      oval=str2double(oval);
                      if oval == 0, opt.oneperdir=0; else, opt.oneperdir=1; end
                  case 'relative'
                      noval=str2double(oval);
                      if ~isnan(noval)
                          opt.relative=noval;
                          if noval < 1, opt.rfolder=startfolder; else, opt.rfolder=['.' fsep]; end
                      else
                          opt.relative=1;
                          opt.rfolder=oval;
                      end;
                  end
              end
          end
      else
          if ~isfield(opt,'dirs'), opt.dirs = 0; end
          if ~isfield(opt,'maxdepth'), if isfield(opt,'depth'), opt.maxdepth = opt.depth; else, opt.maxdepth = 0; end, end
          if ~isfield(opt,'mindepth'), if isfield(opt,'depth'), opt.mindepth = opt.depth; else, opt.mindepth = 0; end, end
          if ~isfield(opt,'maxage'), opt.maxage = -1; end
          if ~isfield(opt,'minage'), opt.minage = -1; end
          if ~isfield(opt,'oneperdir'), opt.oneperdir = 0; end
          if isfield(opt,'rfolder'), opt=rmfield(opt,'rfolder'); end
          if ~isfield(opt,'relative')
              opt.relative = 0;
              opt.rfolder = startfolder;
          else
              if ischar(opt.relative)
                  opt.rfolder=opt.relative;
                  opt.relative=1;
              else
                  if double(opt.relative) >= 1
                      opt.rfolder = ['.' fsep]; opt.relative = 1;
                  else
                      opt.rfolder = startfolder; opt.relative = 0;
                  end
              end
          end
          if ~isfield(opt,'return'), opt.return = 'cellarr'; end
      end
  end
  if isfield(opt,'cellstr') && opt.cellstr > 0, opt.return = 'cellstr'; end
  if isfield(opt,'chararr') && opt.chararr > 0, opt.return = 'chararr'; end
  if opt.dirs ~= 0, opt.dirs = 1; end
  if ~isa(opt.maxdepth,'double'), opt.maxdepth = 0; end
  if ~isa(opt.mindepth,'double'), opt.mindepth = 0; end
  if ~isa(opt.maxage,'double'), opt.maxage = -1; end
  if ~isa(opt.minage,'double'), opt.minage = -1; end
  if opt.oneperdir ~= 0, opt.oneperdir = 1; end
  if opt.relative ~= 0, opt.relative = 1; else opt.rfolder = startfolder; end
  opt.maxage=opt.maxage/86400; if opt.maxage<0, opt.maxage=-1; end
  opt.minage=opt.minage/86400; if opt.minage<0, opt.minage=-1; end

  % - make real call on behalf of files/dirs
  if opt.dirs == 0
      filesfound = findsubfiles(startfolder,patterns,1,opt.mindepth,opt.maxdepth,opt.minage,opt.maxage,opt.oneperdir,opt.rfolder);
  else
      filesfound = findsubdirs(startfolder,patterns,1,opt.mindepth,opt.maxdepth,opt.minage,opt.maxage,opt.oneperdir,opt.rfolder);
  end

  % - return the correct number of values
  if nargout > 1, numberfound=size(filesfound,2); end
  switch lower(opt.return)
  case {'chararr','cellstr'}
      filesfoundcstr=char([]);
      for fcount=1:size(filesfound,2), filesfoundcstr=strvcat(filesfoundcstr,filesfound{fcount}); end
      filesfound=filesfoundcstr;
  end
  if strcmp(lower(opt.return),'cellstr'), filesfound=cellstr(filesfound); end
end
% - end of vbm_findfiles(...)


% %%%%internal functions%%%%


function found = findsubfiles(path,patterns,adepth,sdepth,mdepth,mnage,mxage,operdir,relative)

  nfound = 0;
  found = { };
  mfilesep = filesep;

  % - first, recursively handle all subfolders, if depth is still valid
  if mdepth == 0 | adepth < mdepth
      ilist = dir(path);
      slist = size(ilist,1);
      [ilistd(1:slist)]=[ilist(:).isdir];
      ilistd=find(ilistd>0);
      slistd=size(ilistd,2);
      for count=1:slistd,
          filestoadd = { };
          if strcmp(ilist(ilistd(count)).name,'.') | strcmp(ilist(ilistd(count)).name,'..'), continue, end;
          filestoadd = findsubfiles([path ilist(ilistd(count)).name mfilesep],patterns,adepth+1,sdepth,mdepth,mnage,mxage,operdir,[relative ilist(ilistd(count)).name mfilesep]);
          sfound=size(filestoadd,2);
          if sfound > 0
              nfoundfrm = nfound + 1;
              nfoundnew = nfound + sfound;
              [found{nfoundfrm:nfoundnew}] = deal(filestoadd{:});
              nfound = nfoundnew;
          end
      end
  end

  % - then, if depth is valid, add files to the output
  if sdepth == 0 | sdepth <= adepth
      if any([mnage,mxage]>=0), rnow=now; end;
      spatt = size(patterns,2);
      for pcount=1:spatt,
          ilist = dir([path patterns{pcount}]);
          slist = size(ilist,1);
          if operdir == 1
              count=1;
              while count<=slist && (ilist(count).isdir>0 | (mnage>=0 && (rnow-datenum(ilist(count).date))<mnage) || (mxage>=0 && (rnow-datenum(ilist(count).date))>mxage)), count = count + 1; end
              if count <= slist
                  nfound = nfound + 1;
                  found{nfound} = [relative ilist(count).name];
              end
          else
              for count=1:slist,
                  if (ilist(count).isdir>0 | (mnage>=0 && (rnow-datenum(ilist(count).date))<mnage) || (mxage>=0 && (rnow-datenum(ilist(count).date))>mxage)), continue; end
                  nfound = nfound + 1;
                  found{nfound} = [relative ilist(count).name];
              end
          end
      end
  end

end
%   end of findsubfiles

function found = findsubdirs(path,patterns,adepth,sdepth,mdepth,mnage,mxage,operdir,relative)

  nfound = 0;
  found = { };
  mfilesep = filesep;

  % - first, recursively handle all subfolders, if depth is still valid
  if mdepth == 0 | adepth < mdepth
      ilist = dir(path);
      slist = size(ilist,1);
      [ilistd(1:slist)]=[ilist(:).isdir];
      ilistd=find(ilistd>0);
      slistd=size(ilistd,2);
      for count=1:slistd,
          filestoadd = { };
          if strcmp(ilist(ilistd(count)).name,'.') | strcmp(ilist(ilistd(count)).name,'..'), continue, end
          filestoadd = findsubdirs([path ilist(ilistd(count)).name mfilesep],patterns,adepth+1,sdepth,mdepth,mnage,mxage,operdir,[relative ilist(ilistd(count)).name mfilesep]);
          sfound=size(filestoadd,2);
          if sfound > 0
              nfoundfrm = nfound + 1;
              nfoundnew = nfound + sfound;
              [found{nfoundfrm:nfoundnew}] = deal(filestoadd{:});
              nfound = nfoundnew;
          end
      end
  end

  % - then, if depth is valid, add folders to the output
  if sdepth == 0 | sdepth <= adepth
      if any([mnage,mxage]>=0), rnow=now; end;
      spatt = size(patterns,2);
      for pcount=1:spatt,
          ilist = dir([path patterns{pcount}]);
          slist = size(ilist,1);
          if operdir == 1
              count=1;
              while count<=slist && (ilist(count).isdir==0 | (mnage>=0 && (rnow-datenum(ilist(count).date))<mnage) || (mxage>=0 && (rnow-datenum(ilist(count).date))>mxage)), count = count + 1; end
              while count <= slist
                  if strcmp(ilist(count).name,'.') | strcmp(ilist(count).name,'..'), count = count + 1; continue, end
                  nfound = nfound + 1;
                  found{nfound} = [relative ilist(count).name];
              end
          else
              for count=1:slist,
                  if (ilist(count).isdir==0 | (mnage>=0 && (rnow-datenum(ilist(count).date))<mnage) || (mxage>=0 && (rnow-datenum(ilist(count).date))>mxage)), continue; end
                  if strcmp(ilist(count).name,'.') | strcmp(ilist(count).name,'..'), continue, end
                  nfound = nfound + 1;
                  found{nfound} = [relative ilist(count).name];
              end
          end
      end
  end

end
%   end of findsubdirs
function [linetocell,cellcount] = splittocell(varargin)
  % splittocell  - split a delimited string into a cell array
  %
  % usage is straight forward:
  %
  % FORMAT:         [outcell,count] = splittocell(string[,delimiters,multi])
  %
  % Input fields:
  %    string       string to split
  %    delimiters   char array containing one or more delimiters
  %                 if left empty -> char(9) == <TAB>
  %    multi        must be '1' (numeric) to be effective, if set
  %                 multiple delimiters will be treated as one
  %
  % Output fields:
  %    outcell      cell array containing the tokens after split
  %    count        number of tokens in result

  % no arguments -> help me!
  if nargin == 0, help(mfilename); return; end

  % initialize return values and varargin{3}
  linetocell=cell(0);
  cellcount =0;
  multidelim=0;



  % do we have useful input ?
  if ~ischar(varargin{1}) | length(varargin{1})==0, return; end
  line=varargin{1};
  if size(line,2) ~= prod(size(line))
      dispdebug('splittocell: input must be a 1xN shaped char array!',4);
      return;
  end

  % are any other arguments specified
  if nargin < 2 | ~ischar(varargin{2})
      delimiter = char(9);
  else
      delimiter = reshape(varargin{2},1,prod(size(varargin{2})));
      if nargin > 2 & isnumeric(varargin{3}) & varargin{3} ~= 0, multidelim = 1; end
  end

  % multi-delimitting requested ?
  if multidelim == 0

      % set initial parameters
      ldelim=size(delimiter,2);
      lline =size(line,2);

      % find occurences of delimiter
      if ldelim==1
          cpos=[(1-ldelim),find(line==delimiter)];
      else
          cpos=[(1-ldelim),findstr(line,delimiter)];
      end
      lcpos =size(cpos,2);

      % any delimiter found at all ?
      if lcpos==1, cellcount=1; linetocell={line}; return; end

      % large array?
      if lcpos < 4096

          % line doesn't end with delimiter ?
          if cpos(lcpos) <= (lline-ldelim)
              % then make it look like it was...
              cpos =[cpos lline+1];
              lcpos=lcpos+1;
          end

          % extract substrings
          for dpos=1:(lcpos-1)
              linetocell{end+1} = line(cpos(dpos)+ldelim:cpos(dpos+1)-1);
          end

      else

          % get good rate
          crate = min(384,floor(lcpos^0.666));

          % iterate over parts
          linetocell={};
          for cmpos = 1:crate:(lcpos-crate)
              linetocell = [linetocell,splittocell(line(cpos(cmpos)+ldelim:cpos(cmpos+crate)-1),delimiter,multidelim)];
          end
          linetocell = [linetocell,splittocell(line(cpos(cmpos+crate)+ldelim:cpos(end)-1),delimiter,multidelim)];

      end

  else

      % set initial parameters
      ldelim=size(delimiter,2);
      lline =size(line,2);

      % find occurences of delimiter
      pdelim = [0];
      for cdelim=1:ldelim
          pdelim = union(pdelim,find(line==delimiter(cdelim)));
      end
      if pdelim(end) ~= lline, pdelim(end+1)=lline+1; end
      lpdel = size(pdelim,2);

      % extract substrings
      if pdelim(2)==1, linetocell{end+1} = ''; end
      for ppdel=1:(lpdel-1)
          if (pdelim(ppdel+1)-1) ~= pdelim(ppdel)
              linetocell{end+1} = line(pdelim(ppdel)+1:pdelim(ppdel+1)-1);
          end
      end

  end

  cellcount=length(linetocell);

end
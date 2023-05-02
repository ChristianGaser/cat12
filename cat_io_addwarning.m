function varargout = cat_io_addwarning(id,mess,level,nline,data,usebox,verb)
%cat_io_addwarning. Collect preprocessing warnings in CAT12
% Uses the global struture cat_err_res to collect warnings in the CAT12
% preprocessing function cat_run_job, cat_run_main, etc. 
% Structure will be exported in cat_tst
% 
% See also ../cat12/html/cat_methods_warnings.html
%
%   cat_io_addwarning(id,mess,level,nline,data,usebox)
%   warning_structure = cat_io_addwarning(level);
%
%   id    .. identifier
%   mess  .. message
%            The message will be printed as with the word WARNING at the 
%            beginning where \\n will not break the line and add some
%            spaces, e.g. 'message line 1 \\nmessage line 2':
%            >>WARNING: message line 1
%            >>         message line 2
%   level .. warning level
%            0 - note    - only relevant for experts/developer
%            1 - caution - uncitical aspectes that could be checked
%            2 - alert   - severe problems that should be checked
%            3 - warning - real matlab warning (full report)
%            4 - error   - real matlab error   (stops processing)
%   nline .. new line [before after] warning or by the following codes
%            1 - add new line in command line output before message
%            2 - add new line in command line output also after meassage
%            3 - add also some space to processing time at the end
%   data  .. structure with fields that are related to the warning,
%            e.g., parameter or test results that cause the warning
%            (in development)
%   usebox .. add a box around the message for command line output
%   verb   .. be verbose
%
% Examples: 
%   cat_io_addwarning('reset') 
%
%   cat_io_addwarning('err99','Processing failed',2,[0 1]);
%   cat_io_addwarning('warn0815','Processing bad',1,[0 1]);
%   cat_io_addwarning('note3','Comment',0,[0 1]);
%
%   cat_io_addwarning    % get all warnings & notes
%   cat_io_addwarning(2) % only get warnings
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  global cat_err_res
  
  if ~isfield(cat_err_res,'cat_warnings')
     cat_err_res.cat_warnings = struct('identifier',{},'message',{},'level',[],'data',{});
  elseif (nargin==1 && strcmp(id,'reset') )
     cat_err_res = rmfield(cat_err_res,'cat_warnings'); 
     cat_err_res.cat_warnings = struct('identifier',{},'message',{},'level',[],'data',{});
     if nargout>0
       varargout{1} = struct('identifier',{},'message',{},'level',[],'data',{});
     end
     return; 
  end

  if nargin > 1 && ischar( id )
    % variables
    if ~exist('nline','var'),   nline = [0 1]; end
    if ~exist('level','var'),   level = 1;  end
    if ~exist('data','var'),    data  = {}; end
    if ~exist('usebox','var'),  usebox = 1; end
    if ~exist('verb','var'),    verb   = 1; end
    % expert = cat_get_defaults('extopts.expertgui');  % not sure if the expert setting work
    
    if numel(nline) == 1
      switch nline
        case 0,     nline2 = [0 0];
        case 1,     nline2 = [1 0]; 
        case 2,     nline2 = [1 1];
        case 3,     nline2 = [1 2];
        otherwise,  nline2 = [1 1];
      end
    elseif numel(nline) == 2 
      nline2 = nline; 
    else
      nline2 = [0 1]; 
    end
    
    if nargin<2
      mess = id; 
    end
    if ~ischar(id)
      error('cat_io_addwarning:idstr','Identifier must be a char'); 
    end
    if ~ischar(mess)
      error('cat_io_addwarning:messstr','Message must be a char'); 
    end
    if ~isnumeric(level)
      error('cat_io_addwarning:levelnum','Level must be numeric'); 
    end
    level = max(0,min(2,level));
    
    cat_err_res.cat_warnings(end+1) = struct('identifier',id,'message',mess,'level',level,'data',{data});
    
    % Define a box that is open at the right side because it is not so
    % easys to replace each linebreak with the right number of spaces.
    usebox = usebox + 1; 
    messi  = [0,strfind(mess,'\\n'),numel(mess)];  
    bsize  = max([ 68 + 4 , diff(messi) + 10 + 4]); % normal + oversize
    box(1) = struct('s','','i','','e','');
    box(2) = struct(...
      's',['    '   repmat('-',1,bsize) '\n'], ... char(9559)
      'i',['    '   ' '], ...
      'e',['\n    ' repmat('-',1,bsize) '']); % char(9595)
    box(3) = struct(...
      's',['    '   repmat('=',1,bsize) '\n'], ... char(9559)
      'i',['    '   ' '], ...
      'e',['\n    ' repmat('=',1,bsize) '']); % char(9595)
    if strcmpi(spm_check_version,'octave') 
      % octave: "warning: range error for conversion to character value"
      box(4) = box(3); 
    else
      box(4) = struct(... % the chars are not vissible in the log file 
        's',['    '   char(9556) repmat(char(9552),1,bsize) '\n'], ... char(9559)
        'i',['    '   char(9553) ' '], ...
        'e',['\n    ' char(9562) repmat(char(9552),1,bsize) '']); % char(9595)
    end
    
    % print output
    if verb
      if nline2(1)>0, fprintf('\n'); end
      warnstr = strrep(mess,'\\n',['\n' box(usebox).i '             ']); 
      if level==0
        cat_io_cmd(sprintf([box(usebox).s box(usebox).i 'NOTE %02d:     ' id '\n' box(usebox).i '             ' warnstr box(usebox).e ],numel(cat_io_addwarning(0))),'note');
      elseif level==1
        cat_io_cmd(sprintf([box(usebox).s box(usebox).i 'WARNING %02d:  ' id '\n' box(usebox).i '             ' warnstr box(usebox).e ],numel(cat_io_addwarning(1))),'warning');
      elseif level==2
        cat_io_cmd(sprintf([box(usebox+1).s box(usebox+1).i 'ALERT %02d:    ' id '\n' box(usebox+1).i '             ' warnstr box(usebox+1).e ],numel(cat_io_addwarning(2))),'error');
      elseif level==3
        warning(id,warnstr)
      else
        error(id,warnstr)
      end
      if nline2(2) == 1 
        fprintf('\n'); 
      elseif nline2(2) == 2
        fprintf('\n'); cat_io_cmd(' ','');
      end
    end
  end
  
  if nargout || nargin==0 || nargin==1 
    if nargin == 0
      varargout{1} = cat_err_res.cat_warnings; 
    else
      if isnumeric( id )
        if numel(cat_err_res.cat_warnings)>0
          varargout{1} = cat_err_res.cat_warnings( [ cat_err_res.cat_warnings(:).level ] == id  ); 
        else
          varargout{1} = struct('identifier',{},'message',{},'level',[],'data',{});
        end
      elseif strcmp(id,'reset') && nargout>0 % required for octave
        cat_err_res = rmfield(cat_err_res,'cat_warnings'); 
        cat_err_res.cat_warnings = struct('identifier',{},'message',{},'level',[],'data',{});
        varargout{1} = struct('identifier',{},'message',{},'level',[],'data',{});
      else
        error('ERROR:cat_io_addwarning:incorrectInput',...
              ['Using only one input is limited to the keyword "reset" and \n' ...
               'to output collected warnings by level (0,1,2,3 or 4): \n' ...
               '  "cat_io_addwarning(''reset''); \n' ...
               '  "warning_structure = cat_io_addwarning(level);"\n']);
      end
    end
  end

return
function varargout = cat_io_addwarning(id,mess,level,nline,data)
%cat_io_addwarning. Collect preprocessing warnings in CAT12
% Uses the global struture cat_err_res to collect warnings in the CAT12
% preprocessing function cat_run_job, cat_run_main, etc. 
% Structure will be exported in cat_tst
% 
% See also ../cat12/html/cat_methods_warnings.html
%
%   cat_io_addwarning(id,mess,level,nline,data)
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
%   nline .. new line [before after] warning or by the following codes
%            1 - add new line in command line output before message
%            2 - add new line in command line output also after meassage
%            3 - add also some space to processing time at the end
%   data  .. structure with fields that are related to the warning,
%            e.g., parameter or test results that cause the warning
%            (in development)
%
% Examples: 
%   cat_io_addwarning('reset') 
%
%   cat_io_addwarning('err99','Precessing failed',2,[0 1]);
%   cat_io_addwarning('warn0815','Precessing bad',1,[0 1]);
%   cat_io_addwarning('note3','Comment',0,[0 1]);
%
%   cat_io_addwarning    % get all warnings & notes
%   cat_io_addwarning(2) % only get warnings
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id$

  global cat_err_res
  
  if ~isfield(cat_err_res,'cat_warnings')
     cat_err_res.cat_warnings = struct('identifier',{},'message',{},'level',[],'data',{});
  elseif (nargin==1 && strcmp(id,'reset') )
     cat_err_res = rmfield(cat_err_res,'cat_warnings'); 
     cat_err_res.cat_warnings = struct('identifier',{},'message',{},'level',[],'data',{});
     return; 
  end

  if nargin > 1 && ischar( id )
    if nargin == 2
      nline = [0 0]; 
      level = 1; 
    elseif nargin == 3 % old defintion
      nline = level; 
      level = 1;
    end 
    if ~exist('data','var'),  data = {}; end
    if ~exist('nline','var'), nline = [0 0]; end
    
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
      nline2 = [0 0]; 
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
    
    cat_err_res.cat_warnings(end+1) = struct('identifier',id,'message',mess,'level',level,'data',{data});
    
    if nline2(1)>0, fprintf('\n'); end
    warnstr = strrep(mess,'\\n','\n         '); 
    if level==0
      cat_io_cmd(sprintf(['NOTE:    ' warnstr]),'note');
    elseif level==1
      cat_io_cmd(sprintf(['CAUTION: ' warnstr]),'caution');
    else
      cat_io_cmd(sprintf(['ALERT:   ' warnstr]),'error');
      
    end
    if nline2(2) == 1 
      fprintf('\n'); 
    elseif nline2(2) == 2
      fprintf('\n'); cat_io_cmd(' ','')
    end
  end
  
  if nargout || nargin==0 || nargin==1 
    if nargin == 0
      varargout{1} = cat_err_res.cat_warnings; 
    else
      if numel(cat_err_res.cat_warnings)>0
        varargout{1} = cat_err_res.cat_warnings( [ cat_err_res.cat_warnings(:).level ] == id  ); 
      else
        varargout{1} = struct('identifier',{},'message',{},'level',[],'data',{});
      end
    end
  end

return
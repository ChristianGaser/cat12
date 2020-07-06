function varargout = cat_io_addwarning(id,mess,cls,nline,data)
%cat_io_addwarning. Collect preprocessing warnings in CAT12
%  Uses the global struture cat_err_res to collect warnings in the CAT12
%  preprocessing function cat_run_job, cat_run_main, etc. 
% 
%   cat_io_addwarning(id,mess,cls,nline,data)
%   warning_structure = cat_io_addwarning(cls);
%
%   id    .. identifier
%   mess  .. message
%            The message will be printed as with the word WARNING at the 
%            beginning where \\n will not break the line and add some
%            spaces, e.g. 'message line 1 \\nmessage line 2':
%            >>WARNING: message line 1
%            >>         message line 2
%   cls   .. warning classes
%            0 - note     (only relevant for experts/developer
%            1 - warning  (uncitical aspectes that could be checked)
%            2 - critical (severe problems that should be checked)
%   nline .. new line [before after] warning or by the following codes
%            1 - add new line in command line output before message
%            2 - add new line in command line output also after meassage
%            3 - add also some space to processing time at the end
%   data  .. structure with fields that are related to the warning,
%            e.g., parameter or test results that cause the warning
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
     cat_err_res.cat_warnings = struct('identifier',{},'message',{},'cls',[],'data',{});
  elseif (nargin==1 && strcmp(id,'reset') )
     cat_err_res = rmfield(cat_err_res,cat_warnings); 
     cat_err_res.cat_warnings = struct('identifier',{},'message',{},'cls',[],'data',{});
  end

  if nargin > 1 && ischar( id )
    if nargin == 3 % old defintion
      nline = cls; 
      cls   = 1;
    end 
    if ~exist('data','var'), data = {}; end
    
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
    
    
    cat_err_res.cat_warnings(end+1) = struct('identifier',id,'message',mess,'cls',cls,'data',{data});
    
    if nline2(1)>0, fprintf('\n'); end
    warnstr = strrep(mess,'\\n','\n         '); 
    if cls==0
      cat_io_cmd(sprintf(['WARNING: ' warnstr]),'comment');
    elseif cls==1
      cat_io_cmd(sprintf(['WARNING: ' warnstr]),'warning');
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
      varargout{1} = cat_err_res.cat_warnings; 
    end
  end

return
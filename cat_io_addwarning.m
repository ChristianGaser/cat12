function varargout = cat_io_addwarning(id,mess,nline)
%cat_io_addwarning. Collect preprocessing warnings in CAT12
%  Uses the global struture cat_err_res to collect warnings in the CAT12
%  preprocessing function cat_run_job, cat_run_main, etc. 
% 
%   cat_io_addwarning(id,mess,nline)
%
%   id    .. identifier
%   mess  .. message
%   nline .. new line 
%            1 - add new line in command line output before message
%            2 - add new line in command line output also after meassage
%            3 - add also some space to processing time at the end
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id$

  global cat_err_res
  
  if isfield(cat_err_res,'cat_warnings')
     cat_err_res.cat_warnings = struct('identifier',{},'message',{});
  end

  if nargin
    cat_err_res.cat_warnings(end+1) = struct('identifier',id,'message',mess);
    warnstr = strrep(mess,'\\n','\n'); 
    warnstr = strrep(warnstr,'\n','\n         '); 
    if exist('nline','var') && nline>0, fprintf('\n'); end
    cat_io_cmd(sprintf(['WARNING: ' warnstr]),'warning');
    if exist('nline','var') && nline>0
      if nline==2
        fprintf('\n'); 
      elseif nline == 3
        fprintf('\n'); cat_io_cmd(' ','')
      end
    end
  end
  
  if nargout
    varargout{1} = cat_err_res.cat_warnings; 
  end

return
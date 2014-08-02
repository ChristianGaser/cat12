function stime = vbm_io_cmd(str,style,strlength,verb,stime)
% ______________________________________________________________________
% Writes a string str with a specific style or color followed by a set 
% of blanks to fit a length of strlenght character. 
% The style variable can be a color that should be a 1x3 matrix with RGB
% value from 0 to 1. Otherwise, it can be a str 'error' (red), 'warning'
% (orange), or 'comment' (blue). If no style is given it is black. 
%
%   stime = vbm_io_cmd(str[,style,strlength,verb,stime])
%
% Example: 
%   stime = vbm_io_cmd('Testfunction','comment',63); 
%   pause(3);
%   fprintf('%3.0fs\n',etime(clock,stime));
%
%   Testfunction                                                      3s
%
%
%   stime = clock; 
%   pause(2.34);
%   stime = vbm_io_cmd('  Prepare segments','g5','',1,stime);
%   Testfunction2                                                     3s
%
% see also: vbm_io_cprintf for colored command line output.
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of neurology
%   University Jena
% ______________________________________________________________________
% $Id$ %

  if ~exist('verb','var') || isempty(verb), verb=1; end
  if ~exist('strlength','var') || isempty(strlength), strlength=66; end
  if verb
    if exist('stime','var') && ~isempty(stime)
      fprintf('% 4.0fs\n',etime(clock,stime));
    end
      
    if ~isempty(str) 
      if exist('style','var')
        vbm_io_cprintf(style,'%s%s',str,repmat(' ',1,1+strlength-length(str))); 
      else
        fprintf('%s:%s',str,repmat(' ',1,strlength-length(str))); 
      end
    end
  end
  stime = clock;  
end
function out = cat_stat_nansum(in, dim)
% ----------------------------------------------------------------------
% Average, not considering NaN values. Similar usage like sum() or 
% MATLAB nansum of the statistic toolbox. Process input as double
% due to errors in large single arrays and set data class of "out" 
% to the data class of "in" at the end of the processing.
% Use dim==0 to evaluate in(:) in case of dimension selection 
% (e.g., in(:,:,:,2) ).
%
% out = cat_stat_nansum(in,dim)
%
% Example 1: 
%   a = rand(4,6,3); 
%   a(rand(size(a))>0.5)=nan; 
%   av = cat_stat_nansum(a,3); 
%   am = nansum(a,3); % of the statistical toolbox ...
%   fprintf('%0.4f %0.4f\n',([av(:),am(:)])');
%
% Example 2 - special test call of example 1:
%   cat_stat_nansum('test')
%
% See also cat_stat_nanmedian, cat_stat_nanstd, cat_stat_nanmean.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  if nargin < 1
    help cat_stat_nansum;
    return;
  end;
  
  if ischar(in) && strcmp(in,'test')
    a = rand(4,6,3); 
    a(rand(size(a))>0.5)=nan; 
    av = cat_stat_nansum(a,3); 
    am = nansum(a,3); % of the statistical toolbox ...
    fprintf('%0.4f %0.4f\n',([av(:),am(:)])');
    out = nanmean(av(:) - am(:)); 
    return; 
  end
  
  if nargin < 2
    if size(in,1) ~= 1
      dim = 1;
    elseif size(in,2) ~= 1
      dim = 2;
    else 
      dim = 3; 
    end;
  end;
  
  if dim == 0 
    in  = in(:); 
    dim = 1; 
  end
  
  if isempty(in), out = 0; return; end
  
  % estimate mean
  %tp    = class(in); 
  tmpin = double(in); % single failed in large arrays
  tmpin(isnan(in(:))) = 0;
  out = sum(tmpin, dim);
  
  %eval(sprintf('out = %s(out);',tp)); % haha ... 
end
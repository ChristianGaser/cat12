function TF = cat_io_contains(str,pat,~,icTF)
%Support contains function for older matlabs. 
%
%   TF = cat_io_contains(str,pat[,'ignoreCase',TRUE|FALSE])
%
% See also contains or use cat_io_contains(1) to run the unit test. 
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id: 2558 2024-02-28 $

  
  switch nargin
    case 0
      help cat_io_contains;
    case 1
      unittest; 
    case {2,4}
      if ischar(str), str = cellstr(str); end
      if ischar(pat), pat = cellstr(pat); end
      
      if ~iscellstr(str)
        error('error:cat_io_contains','First argument must be text.')
      end
      if ~iscellstr(pat)
        error('error:cat_io_contains','Search term must be a text or pattern array.')
      end

      if nargin == 4 && icTF
        str = lower(str);
        pat = lower(pat);
      end

      TF = ~cellfun('isempty', strfind( str , pat(1))); %#ok<STRCL1> 
      if numel(pat) > 1 
        for str2i = 2:numel(pat)
          TF = TF | ~cellfun('isempty', strfind( str , pat{str2i}));  
        end
      end
      
    otherwise
      error('error:cat_io_contains:badInput','Wrong number of input elements!'); 
  end
end
% ======================================================================
function unittest
%unittest with cases from the MATLAB contains help.

  strpats = {
    { {'Mary Ann Jones','Paul Jay Burns','John Paul Smith'} {'Paul'} };
    { {'Mary Ann Jones','Christopher Matthew Burns','John Paul Smith'} {'Ann','Paul'} };
    { {'Mary Ann Jones','Christopher Matthew Burns','John Paul Smith'} 'Ann' };
    { {'Mary Ann Jones','Christopher Matthew Burns','John Paul Smith'} 'ann' };
    { 'peppers, onions, and mushrooms' 'onion'}; 
    { 'peppers, onions, and mushrooms' 'nothing'}; 
    { 1:2 'test'}; 
    { 'test' 1};
  };

  for spi = 1:numel(strpats)
    cat_io_cprintf('blue','\nTestcase %d:\n',spi)

    fprintf('Input1:          '); disp(strpats{spi}{1})
    fprintf('Input2:          '); disp(strpats{spi}{2})
    
    for cs = 0:1

      fprintf('contains:        '); 
      try
        if cs 
          disp(contains(strpats{spi}{1},strpats{spi}{2},'IgnoreCase',1))
        else
          disp(contains(strpats{spi}{1},strpats{spi}{2}))
        end
      catch e
        cat_io_cprintf('err',[e.message '\n']);
      end
  
      fprintf('cat_io_contains: '); 
      try
        if cs
          disp(cat_io_contains(strpats{spi}{1},strpats{spi}{2},'IgnoreCase',1))
        else
          disp(cat_io_contains(strpats{spi}{1},strpats{spi}{2}))
        end
      catch e
        cat_io_cprintf('err',[e.message '\n']);
      end
    end
  end
end
function Y = cat_vol_ctype(Y,type)
% ______________________________________________________________________
% Y = cat_vol_ctype(Y[,type])
%
% Convert datatype with checking of min/max, nan, and rounding for 
% [u]int[8|16], single, double, and char. Default round type is 'uint8'. 
% Y has to be a matrix or cell. 
%
% This function is only written for our private use, mostly to convert 
% single to uint8. We did not check for special behavior, for extremly 
% high values or special rounding issues, or converting to larger 
% classes etc.!
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$
% ______________________________________________________________________
 
  if nargin==0, help cat_vol_ctype; return; end
  if nargin==1 && ischar(Y) && strcmp(Y,'test'), testctype; return; end

  types = {'int8','int16','int32','int64','single','float32','float64'...
           'uint8','uint16','uint32','uint64','double'};

  if ~exist('type','var')
    type = 'uint8';
  else
    type  = lower(type);
    % check for additional information such as -le and cut name 
    ind = strfind(type,'-');
    if ~isempty(ind)
      type = type(1:ind-1);
    end
    if ~any(contains(types,type))  
      error('MATLAB:SPM:CAT:cat_vol_ctype:UnknownType', ...
            ['ERROR: cat_vol_ctype: unknown data type ''%s'' ' ...
             '(only [u]int[8|16], single, and double).'],type);
    end
  end
  % use single for logical arrays to be compatible 
  type = cat_io_strrep(type, ...
    {'float32', 'float64'}, ...
    {'single',  'double'});

  
  if iscell(Y)
    % recall function 
    for yi = 1:numel(Y)
      Y{yi} = cat_vol_ctype(Y{yi}, type);
    end
  else
    type  = types{contains(types, type)};
 
    % prepare conversion
    if contains('int', type) 
      switch class(Y)
        case {'single','double'}
          % replace nan
          Y = single(Y);
          Y(isnan(Y)) = 0; 
          Y = round( min( single(intmax(type)), max(single(intmin(type)), Y )));
        case {'uint8','uint16'}
          Y = min( uint16(intmax(type)), max(uint16(intmin(type)), uint16(Y) ));
        case {'uint32','uint64'}
          Y = min( uint32(intmax(type)), max(uint32(intmin(type)), uint32(Y) ));
        case {'int8','int16'}
          Y = min(  int16(intmax(type)), max( int16(intmin(type)), uint16(Y) ));
        case {'int32','int64'}
          Y = min(  int64(intmax(type)), max( int64(intmin(type)), uint64(Y) ));
        otherwise
          % this is not working for very old matlab versions
          eval(sprintf('Y = min( double(intmax(''%s'')), max( double(intmin(''%s'')), double(Y) ))',type,type));
      end
    elseif contains(type,'single')
      Y = min( single(realmax(type)), max(-single(realmax(type)), single(Y) ));
    elseif contains(type,'double')
      Y = min( double(realmax(type)), max(-double(realmax(type)), double(Y) ));
    end
    
    % convert
    eval(sprintf('Y = %s(Y);', type));
  end
  
end
function testctype
%% unit test with two major cases: 
%   cell B: double to single/double/(u)int(8,16,32,64) 
%   cell c: uint16 to single/double/(u)int(8,16,32,64)

  ncases = {'(uint8)', 'single', 'double', ...
            'uint8', 'uint16', 'uint32', 'uint64', ...
            'int8' , 'int16',  'int32',  'int64'};
  tval   = 512; 
  A      = randn(10,10,10) * tval; 
  C      = cell(1,numel(ncases)); B = C; 

  % default uint8
  C{1}   = cat_vol_ctype( A );
  B{1}   = cat_vol_ctype( int16( round(A) ) );
  % other cases
  for ci = 2:numel(ncases) 
    C{ci} = cat_vol_ctype( A                 , ncases{ci});
    B{ci} = cat_vol_ctype( int16( round(A) ) , ncases{ci});
  end     


  %% plot results
  fh = figure(38478); clf
  fh.Name  = 'cat_vol_ctype unit test';
  fh.Color = [1 1 1];
  
  % plot C
  for ci = 1:numel(ncases) 
    subplot(4,6,ci); 
    imagesc( C{ci}(:,:,round(size(A,3)/2)) );
    title(ncases{ci}); caxis([-tval,tval]); 
    axis equal off; 
  end
  
  % plot B
  for ci = 1:numel(ncases) 
    subplot(4,6,ci + 1 +numel(ncases)); 
    imagesc( B{ci}(:,:,round(size(A,3)/2)) );
    title(ncases{ci}); caxis([-tval,tval])
    axis equal off; 
  end

  fprintf('Test cases single: uint8, char, single, double, uint[8,16,32,64], uint[8,16,32,64]:\n'); 
  disp(C)

  fprintf('Test cases int16:  uint8, char, single, double, uint[8,16,32,64], uint[8,16,32,64]:\n'); 
  disp(B)

  fprintf('Convert Cell:\n')
  cat_vol_ctype( { A, int16( round(A) ) } )
  
  fprintf('Test error:\n')
  cat_vol_ctype( A ,'char');
  

end
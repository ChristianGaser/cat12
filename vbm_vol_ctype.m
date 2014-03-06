function Y=vbm_vol_ctype(Y,type)
% ______________________________________________________________________
% Y=vbm_vol_conv(Y[,type])
%
% Convert datatype with checking of min/max, nan, and rounding for 
% [u]int[8|16], single, double, and char. Default round type is 'uint8'. 
% Y has to be a matrix or cell. 
%
% This function is only writen for our private use, mostly to convert 
% single to uint8. I did not check for special behavior, for extremly 
% high values or special rounding issues, or converting to larger 
% classes etc.!
% ______________________________________________________________________
%
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
%
% $Id$
% ______________________________________________________________________
  
  types = {'int8','int16','int32','int64','single',...
           'uint8','uint16','uint32','uint64','double'};

  if ~exist('type','var');
    type = 'uint8';
  else
    type  = lower(type);
    if all(cellfun('isempty',strfind(types,type)))
      error('MATLAB:SPM:VBM:vbm_vol_ctype:UnknownType', ...
            ['ERROR: vbm_vol_ctype: unknown data type ''%s'' ' ...
             '(only [u]int[8|16], single, and double).'],type);
    end
  end

  
  if iscell(Y)
    % recall function 
    for yi=1:numel(Y)
      Y{yi}=vbm_vol_ctype(Y{yi},type);
    end
  else
    type = types{find(~cellfun('isempty',strfind(types,type)),1,'first')};
 
    % prepare convertation
    if ~isempty(strfind('int',type)) || ~isempty(strfind('char',type))
      switch class(Y)
        case {'single','double'}
          Y = single(Y);
          Y(isnan(Y)) = 0;
          Y = round(min(single(intmax(type)),...
                    max(single(intmin(type)),Y)));
        otherwise
          Y = int64(Y);
          Y = round(min(int64(intmax(type)),...
                    max(int64(intmin(type)),Y)));
      end
    elseif ~isempty(strfind(type,'single'))
      Y = min(single(realmax(type)),max(single(realmin(type)),Y));
    elseif ~isempty(strfind(type,'double'))
      Y = min(double(realmax(type)),max(double(realmin(type)),Y));
    end
    
    % convert
    eval(sprintf('Y = %s(Y);',type));
  end
  
end
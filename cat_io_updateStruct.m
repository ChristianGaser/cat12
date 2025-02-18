function S=cat_io_updateStruct(S,SN,RepByEmpty,ind)
% _________________________________________________________________________
% Add and/or updates entries of the structure 'SN' to the structure 'SI'.
% If 'RepByEmpty=0', fields of SI where only overwriten by nonemtpy fields 
% of 'SN'. If the structure has multiple entries, 'ind' can be used to set 
% a specific field, otherwise all entries will be replaced. 
%
%   S = cat_io_updateStruct(SI,SN[,RepByEmpty,ind])
%   
% WARNING:
%   This function is still in developent! Be careful by using it, due to
%   unexpected behaviour. Updating of structures is a complex topic with 
%   many subcases and here only a simple alignment is used!
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  % check input
  if ~exist('RepByEmpty','var'), RepByEmpty = 0; end
  if ~exist('ind','var') 
    if numel(SN) <= 1, ind = 1; else ind = 1:numel(SN); end
  end
  ind = min(single(intmax),max(1,ind));
  
  
  if numel(SN) > 1
    % multiple element update 
    for sni=ind
      S = cat_io_updateStruct(S,SN(sni),RepByEmpty,sni);
    end
  else  
    % single element update
    fnS = fieldnames(SN);
    for fnSi=1:numel(fnS)
      if isfield(S,fnS{fnSi}) 
        if (ischar(SN.(fnS{fnSi})) && RepByEmpty) || ~isempty(SN.(fnS{fnSi}))
          if isstruct(SN.(fnS{fnSi})) 
            % if the field is a structure too, cat_io_updateStruct has
            % to be used recursive
            if ~isstruct(S(1).(fnS{fnSi}))
              S = rmfield(S,fnS{fnSi});
            end
            if numel(S)<ind
              % if the field does not exist yet, we use the first
              % element for initialization 
              try
                S(ind).(fnS{fnSi}) = cat_io_updateStruct(S(1).(fnS{fnSi}),SN.(fnS{fnSi}),0);
              end
            else
              try
                S(ind).(fnS{fnSi}) = cat_io_updateStruct(S(ind).(fnS{fnSi}),SN.(fnS{fnSi}),RepByEmpty);
              end
            end
          else
            S(ind).(fnS{fnSi}) = SN.(fnS{fnSi});
          end
        end
      else
        S(ind).(fnS{fnSi}) = SN.(fnS{fnSi});
      end
    end
  end
end
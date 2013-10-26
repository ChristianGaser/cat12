function S=vbm_io_updateStruct(S,SN,RepByEmpty)
% _________________________________________________________________________
% Add and/or updates entries of the structure SN to the structure SI. If 
% RepByEmpty=0, fields of SI where only overwriten by nonemtpy fields of 
% SN.
%
%   SO=vbm_io_updateStruct(SI,SN,RepByEmpty)
%   
% _________________________________________________________________________
% $Id$
  if ~exist('empty','var'), RepByEmpty=0; end

  fnS = fieldnames(SN);
  for fnSi=1:numel(fnS)
    if isfield(S,fnS{fnSi}) 
      if RepByEmpty || ~isempty(SN.(fnS{fnSi}))
        if isstruct(SN.(fnS{fnSi})) 
          S.(fnS{fnSi}) = vbm_io_updateStruct(S.(fnS{fnSi}),SN.(fnS{fnSi}),RepByEmpty);
        else
          S.(fnS{fnSi}) = SN.(fnS{fnSi});
        end
      end
    else
      S.(fnS{fnSi}) = SN.(fnS{fnSi});
    end
  end
end
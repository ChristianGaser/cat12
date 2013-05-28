function S=vbm_io_updateStruct(S,SN,RepByEmpty)
  if ~exist('empty','var'), RepByEmpty=0; end

  fnS = fieldnames(SN);
  for fnSi=1:numel(fnS)
    if isfield(S,fnS{fnSi}) 
      if RepByEmpty || ~isempty(SN.(fnS{fnSi}))
        if isstruct(SN.(fnS{fnSi})) 
          S.(fnS{fnSi}) = vbm_io_updateStruct(S.(fnS{fnSi}),SN.(fnS{fnSi}));
        else
          S.(fnS{fnSi}) = SN.(fnS{fnSi});
        end
      end
    else
      S.(fnS{fnSi}) = SN.(fnS{fnSi});
    end
  end
end
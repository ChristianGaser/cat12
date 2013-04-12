function S=vbm_io_updateStruct(S,SN)
  fnS = fieldnames(SN);
  for fnSi=1:numel(fnS)
    if isfield(S,fnS{fnSi}) && isstruct(SN.(fnS{fnSi}))
      S.(fnS{fnSi}) = vbm_io_updateStruct(S.(fnS{fnSi}),SN.(fnS{fnSi}));
    else
      S.(fnS{fnSi}) = SN.(fnS{fnSi});
    end
  end
end
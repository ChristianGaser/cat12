function [TH2,T] = cat_io_struct2table(S,F,verb)
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  if ~exist('F','var') || isempty(F) || (ischar(F) && F=='*') || ...
    (iscell(F) && ~all(cellfun('isempty',strfind(F,'*'))))
    F = getFN(S,1024);
  end

  TH=cell(1,numel(F)); 
  for fi=1:numel(F)
    TH{fi} = strrep(F{fi},'.','_');
  end

  if ~exist('verb','var'), verb = 1; end

  spm_clf('Interactive'); 
  spm_progress_bar('Init',numel(S),'struct2table','entry');


  T=cell(numel(S),numel(F)); TH2=TH;
  for si=1:numel(S)
    % be verbose with progressbar
    if numel(S)>10, end

    fields=''; ef=0;
    for fi=1:numel(F)
      try
        eval(sprintf('iF=numel(S(si).%s);',F{fi}));
        eval(sprintf('cF=isnumeric(S(si).%s);',F{fi}));
        if cF && iF>1 && iF<10
          for iFi=1:iF
            TH2{fi+ef}=sprintf('%s_%d',TH{fi},iFi);
            eval(sprintf('T{si,fi+ef} = S(si).%s(iFi);',F{fi}));
            if iFi<iF, ef=ef+1; end   
          end
          if ef>100; return; end
        else
          TH2{fi+ef}=TH{fi};
          eval(sprintf('T{si,fi+ef} = S(si).%s;',F{fi}));
        end
      catch
        fields=sprintf('%s,%s',fields,F{fi});  
        if verb
          fprintf('Miss field %d - ''%s'' in %d!\n',fi,F{fi},si);
        end
        T{si,fi+ef} = [];
        TH2{fi+ef}  = TH{fi};
      end
    end
    if ~isempty(fields)
      if verb
        fprintf('Miss field [%s] in %d!\n',fields(2:end),si);
      end
    end

    spm_progress_bar('Set',si);
  end

  % clear progressbar
  spm_progress_bar('Clear');
end
% =========================================================================
function FNS = getFN(SS,dimlim)
%getFN(S). Recursive extraction of structure elements as string to eval. 

  if ~exist('dimlim','var'), dimlim = 10; end

  if isempty(SS)
    FNS = SS;
  else
    S   = SS(1);
    FN  = fieldnames(S);
    FNS = {};
    for fni = 1:numel(FN)
      % need this for useful order of fields
      acc = num2str( 1 + round( log10( numel( S.(FN{fni}) ))) );

      if isstruct( S.(FN{fni}) )
        % recursive call in case of structures
        FNI = getFN(S.(FN{fni}),dimlim); 
        if numel(S.(FN{fni})) == 1
          for fnii = 1:numel(FNI)
            FNI{fnii} = [FN{fni} '.' FNI{fnii}]; 
          end
        else
          FNI = {};
          for fnii = 1:numel(FNI)
            for sii = 1:numel(S.(FN{fni}))
              FNI = [FNI; sprintf(['%s(%0' acc 'd).%s'], FN{fni}, sii, FNI{fnii})]; %#ok<AGROW> 
            end
          end
        end
      elseif ischar( S.(FN{fni}) ) 
        FNI{1} = sprintf('%s', FN{fni} ); 
      elseif iscellstr( S.(FN{fni}) ) %#ok<ISCLSTR> 
        FNI = {};
        for fnii = 1:min(dimlim,numel( S.(FN{fni}) ))
          FNI = [FNI; sprintf(['%s{%0' acc 'd}'],FN{fni},fnii) ]; %#ok<AGROW> 
        end
      else
        if numel( S.(FN{fni}) ) == 1
          FNI{1} = sprintf('%s',FN{fni});
        else
          % just extract a limited number of elements
          FNI = {};
          for fnii = 1:min(dimlim,numel( S.(FN{fni}) ))
            FNI = [FNI; sprintf(['%s(%0' acc 'd)'],FN{fni},fnii) ]; %#ok<AGROW> 
          end
        end
      end
      FNS = [FNS; FNI]; %#ok<AGROW> 
    end
    FNS = unique(FNS);
  end
end
function S=cat_io_mergeStruct(S,SN,ri,id)
% _________________________________________________________________________
% Merge to structures 'S' and 'SN'. 
%
%   S = cat_io_mergeStruct(S,SN[,[],id])
% 
%   where id>0 updates elements S(id)  
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
% $Id: 2558 2024-02-28 $

  % check input
  maxri = 20; 
  if ~exist('id','var') || isempty(id), id=0; end
  if ~exist('ri','var') || isempty(ri), ri=0; end
  
  SFN  = fieldnames(S);
  SNFN = fieldnames(SN);

  %% add new empty fields in additional element in S
  NSFN = setdiff(SNFN,SFN);
  if ~isempty(NSFN)
    ne   = numel(S)+1; 
    for ni = 1:numel(NSFN)
      if isnumeric(SN(1).(NSFN{ni})) 
        S(ne).(NSFN{ni}) = [];
      elseif islogical(SN(1).(NSFN{ni}))
        S(ne).(NSFN{ni}) = false; 
      elseif ischar(SN(1).(NSFN{ni}))
        S(ne).(NSFN{ni}) = '';
      elseif isstruct(SN(1).(NSFN{ni}))
        if ri<maxri
          Stmp = cat_io_mergeStruct(struct(),SN(1).(NSFN{ni})(1),ri+1);
        else
          Stmp = struct(); 
        end
        S(ne).(NSFN{ni}) = Stmp(1);
      elseif iscell(SN(1).(NSFN{ni}))
        S(ne).(NSFN{ni}) = {};
      end
    end
    S(ne) = []; 
  end

  %%
  NSSFN = setdiff(SFN,SNFN);
  if ~isempty(NSSFN)
    sne   = numel(SN) + 1; 
    for ni = 1:numel(NSSFN)
      if isnumeric(S(1).(NSSFN{ni})) 
        SN(sne).(NSSFN{ni}) = [];
      elseif islogical(S(1).(NSSFN{ni}))
        SN(sne).(NSSFN{ni}) = false(0);
      elseif ischar(S(1).(NSSFN{ni}))
        SN(sne).(NSSFN{ni}) = '';
      elseif isstruct(S(1).(NSSFN{ni}))
        if ri<maxri
          Stmp = cat_io_mergeStruct(struct(),S(1).(NSSFN{ni})(1),ri+1);
        else
          Stmp = struct(); 
        end
        SN(sne).(NSSFN{ni}) = Stmp(1);
      elseif iscell(S(1).(NSSFN{ni}))
        SN(sne).(NSSFN{ni}) = {};
      end
    end
    SN(sne) = [];
  end

  %%
  S   = orderfields(S);
  SN  = orderfields(SN);
  
  ns  = numel(S); 
  for sni=1:numel(SN)
    if id==0
      S(ns + sni) = SN(sni);
    else
      S(id + sni - 1) = SN(sni);
    end
  end
end
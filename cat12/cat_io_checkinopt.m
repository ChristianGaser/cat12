function res = cat_io_checkinopt(opt, def, cond)
% format: res = cat_io_checkinopt(opt,def,cond)
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id: 2558 2024-02-28 $


  if ~exist('def','var'),  def  = []; end
  if ~exist('cond','var'), cond = []; end

  res = def; 
  if ~isfield(res,'do'),   res.do   = 1; end   
  if ~isfield(res,'verb'), res.verb = 0; end
  if numel(opt)>1, error('ERROR:checkinopt:optsize','ERROR: the size of the parameter struct ''opt'' should be 1!'); end
  if numel(def)>1, error('ERROR:checkinopt:optsize','ERROR: the size of the parameter struct ''def'' should be 1!'); end
  
  % only elements of def will be in res... do not check for subfields!
  %fields = intersect(fieldnames(opt),fieldnames(def)); 
  res = cat_io_updateStruct(def,opt,1);
  
  for r=1:numel(cond)
    str=cond{r}; str=strrep(str,'opt.','res.');str=strrep(str,'def.','res.');
    if ~eval(str)
      error('Condition ''%s'' do not fit: %s',str,evalc('res'));
    end
  end
return
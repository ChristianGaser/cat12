function res = checkinopt(opt, def, cond)
% format: res = checkin(opt,def,cond)
% _________________________________________________________________________
% Robert Dahnke 2011_01
% Center of Neuroimaging 
% University Jena
% $Id$
% _________________________________________________________________________  

  if ~exist('def','var'),  def=[]; end
  if ~exist('cond','var'), cond=[]; end

  res = def; 
  %res.opt = opt; res.def = def; res.cond = cond;
  if ~isfield(res,'do'),   res.do   = 1; end   
  if ~isfield(res,'verb'), res.verb = 0; end
  if numel(opt)>1, error('ERROR:checkinopt:optsize','ERROR: the size of the parameter struct ''opt'' should be 1!'); end
  if numel(def)>1, error('ERROR:checkinopt:optsize','ERROR: the size of the parameter struct ''def'' should be 1!'); end
  
  % only elments of def will be in res... do not check for subfields!
  %fields = intersect(fieldnames(opt),fieldnames(def)); 
  res = vbm_io_updateStruct(def,opt,1);
  
  for r=1:numel(cond)
    str=cond{r}; str=strrep(str,'opt.','res.');str=strrep(str,'def.','res.');
    if ~eval(str),
      error('Condition ''%s'' do not fit: %s',str,evalc('res'));
    end
  end
return
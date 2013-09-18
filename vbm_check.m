function varargout = vbm_check(action,varargin)
% Check input files
% ______________________________________________________________________
% Blabla
%
% [INtype,F,V,T] = vbm_check_files('checkinfiles',IN);
% OUT            = vbm_check_files('checkinfiles',INtype,F,V,T);
% res            = vbm_check_files('checkinopt',opt,def,cond);
% 
% INtype = [1=F(file) | 2=V(hdr) | 3=V(hdr)+T(volume) | 4=T(volume)]
% ______________________________________________________________________
% Robert Dahnke 2012_10
% Structural Brain Mapping Group
% University Jena
%  
% $Id$
% ______________________________________________________________________

  if nargin==0
    error('MATLAB:vbm_check','Need some intput!\n');
  end
  
  warnstat = warning('query');
  warning('off','all'); 
  % 'QFORM0 representation has been rounded.');
  switch action
    case 'checkinfiles',  
      switch num2str(nargout)
        case '1',  varargout{1} = vbm_checkinfiles(varargin);
        case '2', [varargout{1},varargout{2}] = vbm_checkinfiles(varargin);
        case '3', [varargout{1},varargout{2},varargout{3}] = vbm_checkinfiles(varargin);
        case '4', [varargout{1},varargout{2},varargout{3},varargout{4}] = vbm_checkinfiles(varargin);
      end      
    case 'checkoutfiles', varargout{1} = vbm_checkoutfiles(varargin);
    case 'checkinopt',    varargout{1} = vbm_checkinopt(varargin);
    otherwise, error('MATLAB:vbm_check','Unknown action ''%s''',action);       
  end
  warning(warnstat(1).state,'all');
end
function varargout = vbm_checkinfiles(varargin)
  if numel(varargin{1})>0, IN       = varargin{1}{1}; end
  if numel(varargin{1})>1, readvols = varargin{1}{2}; else readvols = 0; end;


  % 1=F (file), 2=V (hdr), 3=V (hdr) + T (volume), 4= T (volume)
  if ischar(IN) || iscell(IN) && all(cellfun('isclass',IN,'char')) 
    INtype = 1;
  elseif isstruct(IN) && isfield(IN,'fname') && ~isfield(IN,'dat')
    INtype = 2;
  elseif isstruct(IN) && isfield(IN,'fname') && isfield(IN,'dat')
    INtype = 3;    
  elseif isnumeric(IN)
    INtype = 4; 
  else
    error('MATLAB:vbm_check_files','Unknown input.\n');
  end
  
  % create output
  varargout{1} = INtype;
  varargout{2} = {''};
  varargout{3} = struct();
  varargout{4} = [];
  
  switch num2str(INtype)
    case '1'
      % have to check for no image files
      IN=cellstr(IN);
      for i=numel(IN):-1:1
        if ~isempty(IN{i})
          [pp,ff,ee]=spm_fileparts(IN{i});
          if ~exist(fullfile(pp,[ff ee]),'file')
            fprintf('Cannot find ''%s''!\n',IN{i});
            IN(i)=[]; %:end-1)=IN(i+1:end); % remove file
          else
            [pp,ff,ee]=spm_fileparts(IN{i}); dd=dir(fullfile(pp,[ff,ee]));
            if ~any(strcmp(ee,{'.nii','.img'})) || dd.bytes<2^16
              IN(i)=[]; %IN(i:end-1)=IN(i+1:end); % remove file
            end
          end
        end
      end
    % ########## Error management if file not exist  
      if nargout>1,             varargout{2} = IN; end
      if nargout>2,             varargout{3} = spm_vol(char(IN)); end
      if nargout>3 && readvols, varargout{4} = spm_read_vols(varargout{3}); end
    case '2',
      if nargout>1,             varargout{2} = cellstr([IN.fname]); end
      if nargout>2,             varargout{3} = IN; end
      if nargout>3 && readvols, varargout{4} = spm_read_vols(varargout{3}); end
    case '3',
      if nargout>3 && readvols
        if numel(IN)==1,        varargout{4} = IN.dat; 
        else                    varargout{4} = {IN.dat};
        end
      end
      if nargout>2,             varargout{3} = IN; end
      if nargout>1,             varargout{2} = cellstr([IN.fname]); end
    case '4', 
      if nargout>1,             varargout{2} = {''}; end
      if nargout>2,             varargout{3} = struct(); end
      if nargout>3,             varargout{4} = IN; end
  end
end
function varargout = vbm_checkoutfiles(varargin)
  if numel(varargin{1})>0, INtype = varargin{1}{1}; end
  if numel(varargin{1})>1, F      = [varargin{1}{2}(:).fname]; end
  if numel(varargin{1})>1, V      = varargin{1}{2}; end
 %if numel(varargin{2})>3, T      = varargin{1}{4}; end; ... later to save memory

  varargout{1} = 'ERROR';
  switch num2str(INtype)
    case '1', 
      if nargout>0 && exist('F','var')
        varargout{1} = F; 
      else 
        error('MATLAB:vbm_check','Need second input for filename output!\n');
      end
    case '2', 
      if nargout>0 && exist('V','var')
        varargout{1} = V; 
      else 
        error('MATLAB:vbm_check','Need third input for header output!\n');
      end
    case '3', 
      if nargout>0 && exist('T','var')
        varargout{1} = V; 
        varargout{1}.dat = varargin{1}{4}; 
      else
        error('MATLAB:vbm_check','Need fourth input for volume output!!\n');
      end
    case '4', 
      if nargout>0, 
        varargout{1} = varargin{1}{4};
      else
        error('MATLAB:vbm_check','Need fourth input for volume output!!\n');
      end
    otherwise, error('MATLAB:vbm_check','Unkown INtype ''%d''!\n',INtype);
  end
end
function varargout = vbm_checkinopt(varargin)
  if numel(varargin{1})>0, opt  = varargin{1}{1}; 
  else  error('MATLAB:vbm_check_files:checkinopt','Need at least one intput!\n');
  end
  if numel(varargin{1})>1, def  = varargin{1}{2}; else def=[]; end
  if numel(varargin{1})>2, cond = varargin{1}{3}; else cond=[]; end

  res = def; 
  %res.opt = opt; res.def = def; res.cond = cond;
  if ~isfield(res,'do'),   res.do   = 1; end   
  if ~isfield(res,'verb'), res.verb = 0; end
  
  % only elments of def will be in res... do not check for subfields!
  %fields = intersect(fieldnames(opt),fieldnames(def)); 
  fields = fieldnames(opt); 
  for fn = 1:numel(fields), res.(fields{fn}) = opt.(fields{fn}); end
  
  for r=1:numel(cond)
    str=cond{r}; str=strrep(str,'opt.','res.');str=strrep(str,'def.','res.');
    if ~eval(str),
      error('Condition ''%s'' do not fit: %s',str,evalc('res'));
    end
  end
  
  % set output filenames if possible
  if isfield(res,'pre') && isfield(res,'fname')
    fields = fieldnames(res.pre);
    for fn = 1:numel(fields)
      opt.fnames = vbm_io_handle_pre(res,'fname',res.pre.(fields{fn}));
    end
  end
  
  % set output variable
  varargout{1}=res;
end
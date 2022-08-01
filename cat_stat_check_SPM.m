function job = cat_stat_check_SPM(job)
% cat_stat_check_SPM to check Z-score across sample using design matrix in SPM.mat
% and check design matrix for orthogonality 
%
% Calls cat_stat_homogeneity and splits data into different samples
% according to the defined block (for cross-sectional data) or 
% subject effects (for longitudinal data).
% Furthermore, the design matrix is used to adjust the data and
% to test for orthogonality of the parameters.
%
% Required job fields are:
%
% .SPM           - SPM.mat filename
% .use_unsmoothed_data - use unsmoothed data if found
% .adjust_data       - adjust data using design matrix
% .check_zscore      - check sample homogeneity using Z-score
% .check_ortho       - check for design orthogonality
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

if ~nargin
  error('No argument given.\n');
else
  spmmat = fullfile(char(job.dir),'SPM.mat');
  if exist(spmmat,'file'), load(spmmat);
  else error('File %s not found.',spmmat); end
  
  if isfield(job.check_SPM,'check_SPM_ortho')
    check_ortho  = job.check_SPM.check_SPM_ortho;
  else
    check_ortho = false;
  end
  check_zscore = isfield(job.check_SPM.check_SPM_zscore,'do_check_zscore');
  
  if check_zscore
    use_unsmoothed_data = job.check_SPM.check_SPM_zscore.do_check_zscore.use_unsmoothed_data;
    adjust_data         = job.check_SPM.check_SPM_zscore.do_check_zscore.adjust_data;
  else
    use_unsmoothed_data = false;
    adjust_data         = false;
  end
end

if check_zscore
  fprintf('\n-------------------------------------------\n');
  fprintf('   Check sample homogeneity\n');
  fprintf('-------------------------------------------\n');

  % get some parameters from SPM
  xX     = SPM.xX;
  VY     = SPM.xY.VY;    

  % check whether for the first file unsmoothed data exists
  % and find begin of (unsmoothed) filename
  [pth,nam,ext] = spm_fileparts(VY(1).fname);
  unsmoothed_found = 0;
  ind_str = 2;
  if  (nam(1) == 's')
    while ind_str < 6
      newname = fullfile(pth,[nam(ind_str:end) ext]);
      
      % check for unsmmothed file
      if exist(newname,'file')
        unsmoothed_found = 1;
        break
      end
      
      ind_str = ind_str + 1;
    end
  end
  
  % check whether all other data also have related unsmoothed data
  if unsmoothed_found
    VY_unsmoothed = VY;
    for i=1:numel(VY)
      [pth,nam,ext] = spm_fileparts(VY(i).fname);
      VY_unsmoothed(i).fname = fullfile(pth,[nam(ind_str:end) ext]);
      
      % if no unsmoothed data could be found disable use
      % of unsmoothed data
      if ~exist(newname,'file')
        unsmoothed_found = 0;
      end
      
    end
  end
  
  % allow to use unsmoothed data
  if unsmoothed_found
    if use_unsmoothed_data
      fprintf('\nUse unsmoothed data\n');
      VY = VY_unsmoothed;
    end
  else
    if use_unsmoothed_data
      fprintf('\nNo unsmoothed data found. Use original data from design matrix.\n');
    end
  end
  
  if spm_mesh_detect(VY)
    mesh_detected = 1;
  else
    mesh_detected = 0;
  end
  
  % sometimes xX.iB and xX.iH are not correct and cannot be used to reliably recognize the design
  xX = correct_xX(xX);
  
  % check for longitudinal designs (paired t-test, flexible factorial)
  repeated_anova = ~isempty(xX.iB);
  
  if repeated_anova
    % coding of subject factor should be hopefully always 2nd column of xX.I
    cl = 2;
    n_samples = max(xX.I(:,cl));
    
    % if number of samples gets too large we just use one sample
    if n_samples > 50
      n_samples = 1;
      cl = 1;
    end
  else
    if ~isempty(xX.iH)
      n_samples = length(xX.iH);
      [rw,cl] = find(xX.I == length(xX.iH)); 
    elseif ~isempty(xX.iC)
      n_samples = length(xX.iC);
      [rw,cl] = find(xX.I == length(xX.iC)); 
    else
      error('Design cannot be recognized.')
    end
  end
  
  % always use last found column
  cl = max(cl);
  
  for i=1:numel(VY)
    if ~exist(char(VY(i).fname),'file')
      fprintf('Error: File %s could not be found.\nPlease check that data or analysis have not moved.\n',char(VY(i).fname));
      return
    end
  end
  
  % select data for each sample
  job_check_zscore.data = cell(n_samples,1);
  for i=1:n_samples
    ind = find(xX.I(:,cl)==i);
    job_check_zscore.data{i} = char(VY(ind).fname);
  end
  
  % don't use parameter files for quality measures
  job_check_zscore.data_xml = '';
  
  % adjust data using whole design matrix
  if adjust_data && repeated_anova % Don't adjust data for long. designs
    fprintf('Disable adjustment for longitudinal data.\n');
    job_check_zscore.c = [];
  elseif adjust_data % adjust using design matrix
    job_check_zscore.c{1} = SPM.xX.X;
    fprintf('Data are adjusted using design matrix.\n');
  else % Don't adjust data
    job_check_zscore.c = [];
  end
  
  % check for global scaling
  if ~all(SPM.xGX.gSF==1)
    job_check_zscore.gSF = SPM.xGX.gSF;
  end
  
  % batch mode 
  if nargin > 0 && isfield(job.check_SPM.check_SPM_zscore.do_check_zscore,'save')
    job_check_zscore.fname  = job.check_SPM.check_SPM_zscore.do_check_zscore.fname; 
    job_check_zscore.outdir = job.check_SPM.check_SPM_zscore.do_check_zscore.outdir; 
    job_check_zscore.save   = job.check_SPM.check_SPM_zscore.do_check_zscore.save; 
  end

  % use external defined window
  job_check_zscore.new_fig = true;

  % also rescue fecatorial design
  job_check_zscore.factorial_design = job;
  
  cat_stat_homogeneity(job_check_zscore);
end

if check_ortho
  fprintf('\n-------------------------------------------\n');
  fprintf('    Check design orthogonality\n');
  fprintf('-------------------------------------------\n');
  h = check_orthogonality(SPM.xX);
  
  if nargin > 0 && isfield(job.check_SPM.check_SPM_zscore,'do_check_zscore') && isfield(job.check_SPM.check_SPM_zscore.do_check_zscore,'save')
    %%
    if ~isempty(job.check_SPM.check_SPM_zscore.do_check_zscore.fname)
    dpi = cat_get_defaults('print.dpi'); 
    if isempty(dpi), dpi = 150; end

    if isempty(job.check_SPM.check_SPM_zscore.do_check_zscore.outdir{1}), job.check_SPM.check_SPM_zscore.do_check_zscore.outdir{1} = pwd; end

    % save
    warning('OFF','MATLAB:print:UIControlsScaled');
    fname = fullfile(job.check_SPM.check_SPM_zscore.do_check_zscore.outdir{1},[job.check_SPM.check_SPM_zscore.do_check_zscore.fname 'DesignOrthogonality.png']);
    try, print(h, '-dpng', '-opengl', sprintf('-r%d',dpi), fname); end
    warning('ON','MATLAB:print:UIControlsScaled');
    end

    
    if job.check_SPM.check_SPM_zscore.do_check_zscore.save>1
      close(h)
    end
  end
end

%---------------------------------------------------------------

function h = check_orthogonality(varargin)
% modified function desorth from spm_DesRep.m for checking design orthogonality

if ~isstruct(varargin{1})
  error('design matrix structure required')
else
  xX = varargin{1};
end

%-Locate DesMtx (X), scaled DesMtx (nX) & get parameter names (Xnames)
%--------------------------------------------------------------------------
if isfield(xX,'xKXs') && ...
    ~isempty(xX.xKXs) && isstruct(xX.xKXs)
  X = xX.xKXs.X;
elseif isfield(xX,'X') && ~isempty(xX.X)
  X = xX.X;
else
  error('Can''t find DesMtx in this structure!')
end

[nScan,nPar] = size(X);

if isfield(xX,'nKX') && ~isempty(xX.nKX)
  inX = 1; else inX = 0; end

if isfield(xX,'name') && ~isempty(xX.name)
  Xnames = xX.name; else Xnames = {}; end

%-Compute design orthogonality matrix
%--------------------------------------------------------------------------
% columns with covariates and nuisance parameters must be mean corrected
% to reliably estimate orthogonality 
X(:,[xX.iC xX.iG]) = X(:,[xX.iC xX.iG]) - repmat(mean(X(:,[xX.iC xX.iG])), nScan, 1);

tmp = sqrt(sum(X.^2));
O = X'*X./kron(tmp',tmp);
tmp = sum(X);
tmp = abs(tmp)<eps*1e5;
bC  = kron(tmp',tmp);


%-Display
%==========================================================================
% create figure
ws = spm('Winsize','Graphics');
FS = cat_get_defaults('extopts.fontsize');

h = figure(20);
clf(h);

set(h,'MenuBar','none','Position',[10 10 0.85*ws(3) 0.85*ws(4)],'NumberTitle','off',...
  'Color',[1 1 1]);

%-Title
%--------------------------------------------------------------------------
hTax = axes('Position',[0.03,0,0.94,1],...
  'DefaultTextFontSize',FS,...
  'XLim',[0,1],'YLim',[0,1],...
  'Visible','off');

str='Statistical analysis: Design orthogonality';
text(0.5,0.95,str,'Fontsize',FS+5,'Fontweight','Bold',...
  'HorizontalAlignment','center');

line('Parent',hTax,...
  'XData',[0.3 0.7],'YData',[0.92 0.92],'LineWidth',3,'Color','r');


%-Display design matrix
%--------------------------------------------------------------------------
hDesMtx = axes('Position',[.07 .4 .6 .4]);
if inX    %-Got a scaled DesMtx
  hDesMtxIm = image((xX.nKX + 1)*32);
else
  hDesMtxIm = image((spm_DesMtx('sca',X,     Xnames) + 1)*32);
end

STick = spm_DesRep('ScanTick',nScan,32);
PTick = spm_DesRep('ScanTick',nPar,32);

set(hDesMtx,'TickDir','out',...
  'XTick',PTick,'XTickLabel','',...
  'YTick',STick,'YTickLabel','')
  ylabel('design matrix')

%-Parameter names
if ~isempty(Xnames)
  axes('Position',[.07 .8 .6 .1],'Visible','off',...
    'DefaultTextFontSize',FS-1,'DefaultTextInterpreter','TeX',...
    'XLim',[0,nPar]+0.5)
  for i=PTick, text(i,.05,Xnames{i},'Rotation',90), end
end

%-Setup callbacks to allow interrogation of design matrix
%--------------------------------------------------------------------------
set(hDesMtxIm,'UserData',...
  struct('X',X,'Xnames',{Xnames}))
set(hDesMtxIm,'ButtonDownFcn','spm_DesRep(''SurfDesMtx_CB'')')

%-Design orthogonality
%----------------------------------------------------------------------
hDesO = axes('Position',[.07 .18 .6 .2]);
tmp = 1-abs(O); tmp(logical(tril(ones(nPar),-1))) = 1;
hDesOIm = image(tmp*64);

% check for longitudinal designs (paired t-test, flexible factorial)
repeated_anova = ~isempty(xX.iB);

% orthogonality is difficult to judge for longitudinal designs
if ~repeated_anova
  % one-sample t-test or regression (just one sample effect)
  if numel(xX.iH) == 1
    ind_ortho = [xX.iC xX.iG];
  % interaction designs (have "@sF" in their name)
  elseif ~isempty(strfind(xX.name,'@sF'))
    ind_nointeract = [];
    iCG = [xX.iC xX.iG];
    for i=1:numel(iCG)
      if isempty(strfind(xX.name{iCG(i)},'@sF'))
        ind_nointeract = [ind_nointeract iCG(i)];
      end
    end
    ind_ortho = [xX.iH ind_nointeract];
  else
    ind_ortho = [xX.iH xX.iC xX.iG];
  end

  if any(any(bC(ind_ortho,ind_ortho)))
    fprintf('\nRed frames indicate combinations between covariates/nuisance variables');
    if numel(xX.iH) > 1
      fprintf(' and group factors.\n');
    else fprintf('\n'); end
    fprintf('Orthogonality between nuisance parameters and parameters of interest should be carefully checked for high (absolute) values that point to co-linearity (correlation) between these variables.\n');
    fprintf('In case of such a high co-linearity nuisance parameters should be preferably used with global scaling (optionally with AnCova instead of proportional scaling).\n');
    fprintf('For more information please check the manual or the online help.\n\n');
  end
  
  for i=1:size(ind_ortho,2)
    for j=1:size(ind_ortho,2)
      if bC(ind_ortho(i),ind_ortho(i)) & i>j
        rectangle('Position', [ind_ortho(i)-0.5 ind_ortho(j)-0.5 1 1],...
          'EdgeColor','r');
        fprintf('Orthogonality between %s and %s:\t %g\n',xX.name{ind_ortho(i)},...
          xX.name{ind_ortho(j)},O(ind_ortho(i),ind_ortho(j)));
      end
    end
  end
end

set(hDesO,'Box','off','TickDir','out',...
  'XaxisLocation','top','XTick',PTick,'XTickLabel','',...
  'YaxisLocation','right','YTick',PTick,'YTickLabel','',...
  'YDir','reverse')
tmp = [1,1]'*[[0:nPar]+0.5];
line('Xdata',tmp(1:end-1)','Ydata',tmp(2:end)')

xlabel('design orthogonality')
set(get(hDesO,'Xlabel'),'Position',[0.5,nPar,0],...
  'HorizontalAlignment','left',...
  'VerticalAlignment','top')
set(hDesOIm,...
  'UserData',struct('O',O,'bC',bC,'Xnames',{Xnames}),...
  'ButtonDownFcn','spm_DesRep(''SurfDesO_CB'')')

if ~isempty(Xnames)
  axes('Position',[.69 .18 0.01 .2],'Visible','off',...
    'DefaultTextFontSize',FS+1,...
    'DefaultTextInterpreter','TeX',...
    'YDir','reverse','YLim',[0,nPar]+0.5)
  for i=PTick
    text(0,i,Xnames{i},'HorizontalAlignment','left')
  end
end


%-Design descriptions
%--------------------------------------------------------------------------
str = '';
%line('Parent',hTax,...
%  'XData',[0.3 0.7],'YData',[0.16 0.16],'LineWidth',3,'Color','r')
hAx = axes('Position',[0.03,0.08,0.94,0.08],'Visible','off');
xs = struct('Measure',  ['abs. value of cosine of angle between ',...
       'columns of design matrix'],...
    'Scale',  {{  'black - colinear (cos=+1/-1)';...
        'white - orthogonal (cos=0)';...
        'gray  - not orthogonal or colinear'}});

set(hAx,'Units','Pixel');
AxPos = get(hAx,'Position');
set(hAx,'YLim',[0,AxPos(4)])

dy = FS; y0 = floor(AxPos(4)) -dy; y = y0;

text(0.3,y,str,...
  'HorizontalAlignment','Center',...
  'FontWeight','Bold','FontSize',FS+2)
y=y-2*dy;

for sf = fieldnames(xs)'
  text(0.3,y,[strrep(sf{1},'_',' '),' :'],...
    'HorizontalAlignment','Right','FontWeight','Bold',...
    'FontSize',FS)
  s = xs.(sf{1});
  if ~iscellstr(s), s={s}; end
  for i=1:numel(s)
    text(0.31,y,s{i},'FontSize',FS)
    y=y-dy;
  end
end

colormap(gray)


%---------------------------------------------------------------

function xX = correct_xX(xX)

% vector of covariates and nuisance variables
iCG = [xX.iC xX.iG];
iHB = [xX.iH xX.iB];

% set columns with covariates and nuisance variables to zero
X = xX.X;
X(:,iCG) = 0;

ncol = size(X,2);

% calculate sum of columns
% The idea behind this is that for each factor the sum of all of its columns should be "1".
Xsum = zeros(size(X));
for i=1:ncol
  % only sum up columns without covariates and nuisance variables
  if isempty(find(iCG==i))
  Xsum(:,i) = sum(X(:,1:i),2);
  end
end

% find columns where all entries are constant except zeros entries
% that indicate columns with covariates and nuisance variables
ind = find(any(diff(Xsum))==0 & sum(Xsum)>0);

% no more than 2 factors expected
if length(ind) > 2
  error('Weird design was found that cannot be analyzed correctly.');
end

% correction is only necessary if 2 factors (iH/iB) were found
if length(ind) > 1
  iF = cell(length(ind),1);

  j = 1;
  % skip columns with covariates and nuisance variables
  while find(iCG==j),  j = j + 1; end

  for i=j:length(ind)
  iF{i} = [j:ind(i)];
  
  j = ind(i)+1;
  % skip columns with covariates and nuisance variables
  while find(iCG==j), j = j + 1; end
  end
  
  % not sure whether this will always work but usually iB (subject effects) should be larger than iH (time effects)
%  if length(iF{1}) > length(iF{2})
if 0 % will be probably not always correct 
  xX.iB = iF{1};
  xX.iH = iF{2};
  else
  xX.iB = iF{2};
  xX.iH = iF{1};
  end
end

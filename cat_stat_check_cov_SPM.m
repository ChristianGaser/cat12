function cat_stat_check_cov_SPM(SPM)
%cat_stat_check_cov_SPM to check covariance across sample using design matrix in SPM.mat
%
% Calls cat_stat_check_cov and splits data into different samples
% according to the defined block (for cross-sectional data) or 
% subject effects (for longitudinal data).
% Furthermore, the design matrix is used to adjust the data.
% If a contrast is defined the columns of this contrast are 
% excluded from the nuisance parameters. 
% Then, only the remaining columns of the design matrix are used as 
% nuisance parameters and not the parameters of interest.
%
%_______________________________________________________________________
% Christian Gaser
% $Id: cat_stat_check_cov_SPM.m 1267 2018-02-02 14:47:02Z gaser $

if nargin == 0
  P = spm_select(1,'SPM.mat','Select SPM.mat');
  load(P)
end

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

% check whether all other data also have respective unsmmothed data
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
  if spm_input('Use unsmoothed data?',1,'yes|no',[1,0],1)
    fprintf('\nUse unsmoothed data\n');
    VY = VY_unsmoothed;
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
  n_samples = length(xX.iB);
  [rw,cl] = find(xX.I == length(xX.iB)); % find column which codes subject factor (length(xX.iB) -> n_samples)  
else
  if ~isempty(xX.iH)
    n_samples = length(xX.iH);
    [rw,cl] = find(xX.I == length(xX.iH)); 
  else
    error('Weird design found that cannot be recognized.')
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
if mesh_detected
  job.data_surf = cell(n_samples,1);
  for i=1:n_samples
    ind = find(xX.I(:,cl)==i);
    job.data_surf{i} = char(VY(ind).fname);
  end
else
  job.data_vol = cell(n_samples,1);
  for i=1:n_samples
    ind = find(xX.I(:,cl)==i);
    job.data_vol{i} = char(VY(ind).fname);
  end
  job.gap = 3;
end

% don't use parameter files for quality measures
job.data_xml = '';

% exclude columns with parameters of interest based on given contrast
Ic = spm_input('Select contrast for adjusting data','+1','m',...
                {SPM.xCon.name,'Adjust with whole design matrix','Do not adjust data'});


switch Ic
case numel(SPM.xCon) + 2 % Don't adjust data
  job.c = [];
case numel(SPM.xCon) + 1 % Duse whole design matrix
  job.c{1} = xX.X;
  fprintf('Data are adjusted with whole design matrix.\n');
otherwise
  % initially use whole design matrix as counfound
  job.c{1} = xX.X;
  c = SPM.xCon(Ic).c;
  [indi, indj] = find(c~=0);
  job.c{1}(:,indi) = [];
  fprintf('Data are adjusted using parameters of interest defined in contrast.\n');
end

cat_stat_check_cov(job);

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

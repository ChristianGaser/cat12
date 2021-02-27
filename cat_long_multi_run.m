function out = cat_long_multi_run(job)
% Call cat_long_main for multiple subjects
%
% Christian Gaser
% $Id$

warning off;

if isdeployed, job.nproc = 0; end

if isfield(job.extopts,'vox') && job.extopts.vox ~= 1.5
  job.extopts.vox = 1.5;
  fprintf('\n---------------------------------------------------------------------------\n');
  fprintf('Warning: No change of output voxel size possible for longitudinal pipeline!\n');
  fprintf('---------------------------------------------------------------------------\n');
end

% use some options from GUI or default file
opts        = job.opts;
extopts     = job.extopts;
output      = job.output;
modulate    = job.modulate;
dartel      = job.dartel;
ROImenu     = job.ROImenu;
longmodel   = job.longmodel;
surfaces    = job.output.surface;

setappdata(0,'job',job);

if isfield(job,'delete_temp')  
  delete_temp = job.delete_temp;
else
  delete_temp = 1;
end

% modify job.subj w.r.t. different selection options
if isfield(job,'datalong') 
  if isfield(job.datalong,'subjects')
    job.data  = {};
  
    for si = 1:numel(job.datalong.subjects)
      for ti = 1:numel(job.datalong.subjects{si})
        job.subj(si).mov{ti,1} = job.datalong.subjects{si}{ti}; 
      end
      job.data = [job.data; job.datalong.subjects{si}];
    end
  else
    job.data  = {};
  
    n_ti = numel(job.datalong.timepoints{1});
    for si = 1:numel(job.datalong.timepoints)
      for ti = 1:numel(job.datalong.timepoints{si})
        job.subj(ti).mov(si) = job.datalong.timepoints{si}(ti); 
      end
      % check that number of time points does not differ
      if numel(job.datalong.timepoints{si}) ~= n_ti
        error('Number of time points differ between the subjects. Please take care to select the same number of time points for all subjects!');
      end
      job.data = [job.data job.datalong.timepoints{si}];
    end
  end
  
  % remove datalong field to prevent that modification of job.subjs is called again
  job = rmfield(job,'datalong');
end

job_name = fullfile(spm('dir'),'toolbox','cat12','cat_long_main.txt');

% we have to copy the original txt-file to a matlab file because for deployed versions
% matlab files will be always pre-compiled, but we need the original matlab file untouched
m_job_name = strrep(job_name,'.txt','.m');
if isdeployed
  [status, mesg] = copyfile(job_name,m_job_name,'f');
  if ~status
    error(mesg);
  end
end

% mirror jobs for all subjects
jobs = repmat({m_job_name}, 1, numel(job.subj));

inputs = cell(1, numel(job.subj));

if cat_get_defaults('extopts.subfolders')
  mrifolder    = 'mri';
  surffolder   = 'surf';
  labelfolder  = 'label';
  reportfolder = 'report';
else
  mrifolder    = '';
  surffolder   = ''; 
  labelfolder  = '';
  reportfolder = ''; 
end

out.surf = cell(''); out.thick = cell(''); out.mwp1 = cell('');
out.catreport = cell(''); out.catroi = cell('');
for i=1:numel(job.subj)
  out.sess(i).warps = cell(1,1);
  [pth,nam,ext,num] = spm_fileparts(job.subj(i).mov{1});
  out.sess(i).warps{1} = fullfile(pth,mrifolder,['avg_y_', nam, ext, num]);

  out.sess(i).files = cell(numel(job.subj(i).mov),1);
  m = numel(job.subj(i).mov);
  data = cell(m,1);
  for j=1:m
    [pth,nam,ext,num] = spm_fileparts(job.subj(i).mov{j});
    switch modulate
    case 0
      out.sess(i).files{j} = fullfile(pth,mrifolder,['wp1r', nam, ext, num]);
    case 1
      out.sess(i).files{j} = fullfile(pth,mrifolder,['mwp1r', nam, ext, num]);
    case 2
      out.sess(i).files{j} = fullfile(pth,mrifolder,['m0wp1r', nam, ext, num]);
    end
    data{j} = job.subj(i).mov{j};

    out.mwp1      = [out.mwp1       fullfile(pth,mrifolder   ,['mwp1r'          nam ext num])]; 
    out.surf      = [out.surf       fullfile(pth,surffolder  ,['lh.central.r'   nam '.gii'])]; 
    out.thick     = [out.thick      fullfile(pth,surffolder  ,['lh.thickness.r' nam])]; 
    out.catreport = [out.catreport  fullfile(pth,reportfolder,['cat_r'          nam '.xml'])]; 
    out.catroi    = [out.catroi     fullfile(pth,labelfolder ,['catROI_r'       nam '.xml'])]; 
      
  end
  inputs{1,i} = data;
    
end

% split job and data into separate processes to save computation time
if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index')) && numel(job.subj) > 1
  if nargout==1
    varargout{1} = cat_parallelize(job,mfilename,'subj');
  else
    cat_parallelize(job,mfilename,'subj');
  end
  return
else
  spm_jobman('run',jobs,inputs{:}); 
end
warning on;

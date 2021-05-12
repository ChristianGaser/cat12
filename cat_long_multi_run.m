function out = cat_long_multi_run(job)
% Call cat_long_main for multiple subjects
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
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
useprior    = job.enablepriors;
surfaces    = job.output.surface;
longTPM     = job.longTPM;

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
  txt_fileid = fopen(job_name,'r');
  txt_contents = fread(txt_fileid);
  fclose(txt_fileid);

  m_fileid = fopen(m_job_name,'r');
  m_contents = fread(m_fileid);
  fclose(m_fileid);
  
  % check whether length of txt- and m-file differs or content differs and only then the txt-file will be copied
  % this allows to pre-install the m-file on systems where this file is read-only
  if (length(txt_contents) == length(m_contents) && any(txt_contents ~= m_contents)) || (length(txt_contents) ~= length(m_contents))
    [status, mesg] = copyfile(job_name,m_job_name,'f');
    if ~status
      fprintf(mesg);
      fprintf('\nIf you do not have write permissions, the administrator should copy the %s file to %s after installing the precompiled version. This prevents overwriting the read-only file.\n',job_name,m_job_name);
      return
    end
  end
  
end

% mirror jobs for all subjects
jobs = repmat({m_job_name}, 1, numel(job.subj));

inputs = cell(1, numel(job.subj));

out.surf = cell(''); out.thick = cell(''); out.mwp1 = cell('');
out.catreport = cell(''); out.catroi = cell('');

for i=1:numel(job.subj)
  [mrifolder, reportfolder, surffolder, labelfolder] = cat_io_subfolders(job.subj(i).mov{1},job);
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

function out = cat_long_multi_run(job)
% Call cat_long_main for multiple subjects
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
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
printlong   = job.printlong; 
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

  is_copied = zeros(numel(job.data),1);

  % create BIDS structure and unzip or copy file to BIDS-folder
  if isfield(job.output,'BIDS')
    if isfield(job.output.BIDS,'BIDSyes')
    
      c = 1;
      is_copied = ones(numel(job.data),1);
      for ti = 1:numel(job.subj)

        BIDSfolder = job.output.BIDS.BIDSyes.BIDSfolder;
        
        % get path of first data set and find "sub-" BIDS part
        name1 = spm_file(job.data{c},'fpath');
        ind = min(strfind(name1,'sub-'));

        if ~isempty(ind)
          % remove leading ".." for real BIDS structure
          BIDSfolder = strrep(BIDSfolder,['..' filesep],'');

          length_name = length(name1);

          % Shorten path until "sub-" indicator is found and add additional
          % relative paths to get BIDSfolder relative to "sub-" directories.
          % This is necessary because there might be additional session 
          % folders and more
          while length_name > ind
            name1 = spm_file(name1,'fpath');
            length_name = length(name1);
            BIDSfolder = ['..' filesep BIDSfolder];
          end
        end

        % we need this in job.extopts for cat_io_subfolders
        job.extopts.BIDSfolder = BIDSfolder;
      
        for si = 1:numel(job.subj(ti).mov)
          [pth,nam,ext] = spm_fileparts(job.subj(ti).mov{si});
          rootfolder = cat_io_subfolders(job.subj(ti).mov{si},job);
          % remove additional mri subfolder because registration first will 
          % work in the upper folder
          ind = strfind(rootfolder,[filesep 'mri']);
          if ~isempty(ind), rootfolder(ind:end) = []; end
          
          name = fullfile(pth,[nam ext]);
          newdir = fullfile(pth,rootfolder);
          % uncompress nii.gz files and change file name for job
          if strcmp(ext,'.gz')
            fname = gunzip(name,newdir);
            fprintf('Uncompress %s to %s\n',name,newdir);
            job.subj(ti).mov{si} = char(fname);
            job.data{c} = char(fname);
          else
            if ~exist(newdir), mkdir(newdir); end
            is_copied(c) = 2;
            s = copyfile(name,newdir);
            if ~s
              error('Could not write %s to %s',name,newdir);
            else
              fprintf('Copy %s to %s\n',name,newdir);
              job.subj(ti).mov{si} = fullfile(newdir,[nam ext]);
              job.data{c} = fullfile(newdir,[nam ext]);
            end
          end
          c = c + 1;
        end
      end
    
    else
      job.output  = rmfield(job.output,'BIDS');
      if isfield(job.extopts,'BIDSfolder'), job.extopts = rmfield(job.extopts,'BIDSfolder'); end
      output = job.output;
      extopts = job.extopts;

    end
  end
  
  % also uncompress gz-files for non-BIDS structure
  if ~isfield(job.extopts,'BIDSfolder') || isempty(job.extopts.BIDSfolder)
    c = 1;
    for ti = 1:numel(job.subj)
      for si = 1:numel(job.subj(ti).mov)
        [pth,nam,ext] = spm_fileparts(job.subj(ti).mov{si});
        % uncompress nii.gz files and change file name for job
        if strcmp(ext,'.gz')
          is_copied(c) = 1;
          fname = gunzip(job.subj(ti).mov{si});
          job.subj(ti).mov{si} = char(fname);
          job.data{c} = char(fname);
          fprintf('Uncompress %s\n',job.subj(ti).mov{si});
        end
        c = c + 1;
      end
    end
  else
    % remove BIDS fields because files are now already copied to BIDS-folder
    job.output  = rmfield(job.output,'BIDS');
    job.extopts = rmfield(job.extopts,'BIDSfolder');
    output = job.output;
    extopts = job.extopts;
  end
  
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
  
  out.sess(i).warps    = cell(1,1);
  [pth,nam,ext,num]    = spm_fileparts(job.subj(i).mov{1});
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
    
  % save XML Parameter
  if ~(isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index')) && numel(job.subj) > 1 )
    [pp,ff]    = spm_fileparts(job.subj(i).mov{1});
    longxml    = fullfile( pp , reportfolder , ['catlong_' ff '.xml'] ); 
    jobsx      = rmfield(job,{'data','subj'});
    jobsx.subj = job.subj(i); 
    jobsx.out  = out; 
    jobsx.dirs = struct('mrifolder',mrifolder, 'reportfolder', reportfolder, ...
      'surffolder', surffolder, 'labelfolder', labelfolder, 'pp1', pp, 'ff1', ff);
    cat_io_xml(longxml,struct('parameter',jobsx));
  end
  
  
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

for i=1:numel(job.data)
  [pth,nam,ext] = spm_fileparts(job.data{i});
  if exist('is_copied','var') && is_copied(i)
    spm_unlink(fullfile(pth,[nam ext]));
    if is_copied(i) == 2
      fprintf('Remove copied file %s\n',fullfile(pth,[nam ext]));
    else
      fprintf('Remove unzipped file %s\n',fullfile(pth,[nam ext]));
    end
  end
end

warning on;
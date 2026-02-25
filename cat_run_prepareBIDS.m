function [job, BIDSfolder, logdir] = cat_run_prepareBIDS(job)


  % extract BIDS directory from GUI parameters
  if isfield(job.output,'BIDS') 
    if isfield(job.output.BIDS,'BIDSyes') 
      BIDSfolder = job.output.BIDS.BIDSyes.BIDSfolder;
      resdircase = 0; 
    elseif isfield(job.output.BIDS,'BIDSrel')
      BIDSfolder = job.output.BIDS.BIDSrel.BIDSfolder;
      resdircase = 1;
    elseif isfield(job.output.BIDS,'relative')  
      BIDSfolder = job.output.BIDS.relative.BIDSfolder;
      resdircase = 2;
    else
      BIDSfolder = ''; 
      resdircase = 0; 
    end
  else
    BIDSfolder = ''; 
    resdircase = 0; 
  end
  % also for longitudinal calls
  job.extopts.BIDSfolder = BIDSfolder; 
  job.extopts.resdircase = resdircase; 
  

  %% check for previous CAT BIDS processing 
  if ~isfield(job,'useprior') % long issue
    previousCATderivatives = cat_io_contains(job.data, fullfile('derivatives','CAT')) | ...
      cat_io_contains(job.data, [filesep 'mri' filesep] ); 
    if any(previousCATderivatives)
      [~,C] = spm_str_manip( job.data(previousCATderivatives ) ,'C');
      cat_io_cprintf('warn', sprintf( ...
          ['\nWARNING: Remove files from previous CAT derivatives and/or mri directories: \n' ...
           '\n         %s ' ...
           '%s \n\n\n'], C.s, char( strcat( '\n           .', C.m ))' ));
      
      job.data(previousCATderivatives) = []; 
    end
  end
  
  %% get BIDS
  [sfiles,sfilesBIDS,BIDSsub,devdir,rdevdir,logdir,mdevdir] = cat_io_checkBIDS(job.data, BIDSfolder, resdircase);
  
  %% evaluate input 
  if BIDSfolder
    if numel(unique(sfilesBIDS)) > 1  
      cat_io_cprintf('warn', sprintf( ...
        ['\nWARNING: Please note that you selected %d subjects with BIDS directory structure and %d without. ' ...
         '\n'], sum(sfilesBIDS), sum(~sfilesBIDS) ));
    end
    if numel(unique(mdevdir)) > 1
      [umdevdir,imdevdir] = unique(mdevdir);
      BIDStst = {'\n           .','\n [BIDS]    .'}; 
      [~,C] = spm_str_manip( umdevdir ,'C');
      cat_io_cprintf('warn',  sprintf( ...
        ['\nWARNING: The results will be writen into multiple relative BIDS directories under: \n' ...
         '\n         %s ' ...
         '%s \n\n\n'], C.s, char( strcat( [BIDStst( 1+sfilesBIDS(imdevdir)) ]', C.m ))' ));
    end
  end
  
  %% added path and filesnames - maybe better as separate structure
  job.filedata.help        = [ 'Structure directory and file names. \n' ...
                               ' logdir      .. path for log file \n' ...
                               ' rawdir      .. origin directory of the RAW data \n' ...
                               ' BIDSfolder  .. relative path to the main result directory \n' ...
                               ' BIDSsubs    .. BIDS subject name \n' ...
                               ' BIDSdirs    .. absolution (full) path to the result directories \n' ...
                               ' BIDSrdirs   .. relative (full) path to the result directories \n' ];
  job.filedata.rawdir      = sfiles; 
  job.filedata.logdir      = logdir;
  job.filedata.isBIDS      = sfilesBIDS;
  job.filedata.BIDSsubs    = BIDSsub;
  job.filedata.BIDSdirs    = devdir;
  job.filedata.BIDSrdirs   = rdevdir;
  job.filedata.BIDSfolder  = BIDSfolder; 
end
  


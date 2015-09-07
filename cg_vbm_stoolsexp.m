function stoolsexp = cg_vbm_stoolsexp
%_______________________________________________________________________
% wrapper for calling VBM surface utilities
%_______________________________________________________________________
% Robert Dahnke and Christian Gaser
% $Id$
%_______________________________________________________________________



%% average surface mesh
%-----------------------------------------------------------------------

  avg.data         = cfg_files;
  avg.data.tag     = 'data';
  avg.data.name    = 'Sample';
  avg.data.filter  = 'gifti';
  avg.data.ufilter = '^[rl]h.central';
  avg.data.num     = [1 Inf];
  avg.data.help    = {
    'Select surfaces.'
    };

  avg.meshsmooth         = cfg_entry;
  avg.meshsmooth.tag     = 'meshsmooth';
  avg.meshsmooth.name    = 'Surface Smoothing Iterations';
  avg.meshsmooth.strtype = 'r';
  avg.meshsmooth.num     = [1 Inf];
  avg.meshsmooth.val     = {[0 2 32]};
  avg.meshsmooth.help    = {
    'Smoothing of the average surface. '
    ''
    };

  avg.surfside         = cfg_menu;
  avg.surfside.tag     = 'surfside';
  avg.surfside.name    = 'Side Handling';
  avg.surfside.labels  = {'separate','mirror'};
  avg.surfside.values  = {1,2};
  avg.surfside.val     = {1};
  avg.surfside.help    = {
    'Handling of the cortical hemispheres.'
    ''
    };
 
  avg.surfname         = cfg_entry;
  avg.surfname.tag     = 'surfname';
  avg.surfname.name    = 'Surface Filename';
  avg.surfname.strtype = 's';
  avg.surfname.num     = [1 Inf];
  avg.surfname.val     = {'average'};
  avg.surfname.help    = {'Name of the surface.'};

  avg.outdir         = cfg_files;
  avg.outdir.tag     = 'outdir';
  avg.outdir.name    = 'Output Directory';
  avg.outdir.filter  = 'dir';
  avg.outdir.ufilter = '.*';
  avg.outdir.num     = [0 1];
  avg.outdir.val{1}  = {''};
  avg.outdir.dir     = fullfile(spm('dir'),'toolbox','vbm12');
  avg.outdir.help    = {'Select a directory where files are written.'};

  avg.main          = cfg_exbranch;
  avg.main.tag      = 'avg_surf';
  avg.main.name     = 'Average Surface Mesh';
  avg.main.val      = {
    avg.data ...
    avg.meshsmooth ...
    avg.surfside ...
    avg.surfname ...
    avg.outdir ...
    };
  avg.main.vfiles   = @vfiles_avg;  
  avg.main.prog     = @vbm_surf_avg;
  avg.main.help     = {
    'Averaging of cortical surfaces.'
    ''
    };

  %}
  
%% data smoothing
%-----------------------------------------------------------------------
  smooth.data         = cfg_files;
  smooth.data.tag     = 'data';
  smooth.data.name    = 'Sample';
  smooth.data.filter  = 'any';
  smooth.data.ufilter = '^[rl]h.(?!cent|sphe|defe).*';
  smooth.data.num     = [1 Inf];
  smooth.data.help    = {'Select surface data (texture) files for smoothing.'};
  
  smooth.fwhm         = cfg_entry;
  smooth.fwhm.tag     = 'fwhm';
  smooth.fwhm.name    = 'Smoothing filter size in fwhm';
  smooth.fwhm.strtype = 'r';
  smooth.fwhm.num     = [1 1];
  smooth.fwhm.val     = {15};
  smooth.fwhm.help    = {
    'Select filter size for smoothing. For cortical thickness a good starting value is 15mm, while other surface parameters based on cortex folding (e.g. gyrification, cortical complexity) need a larger filter size of about 25mm.'};
 
  smooth.main      = cfg_exbranch;
  smooth.main.tag  = 'datasmooth';
  smooth.main.name = 'Smooth Surface Data';
  smooth.main.val  = {
    smooth.data ...
    smooth.fwhm ...
  };
  smooth.main.vfiles = @vfiles_smooth;
  smooth.main.prog = @vbm_surf_smooth;
  smooth.main.help = {
    'Gaussian smoothing of surface data (texture).'
    ''
  }; 



%% resample surface (mesh and data)
%-----------------------------------------------------------------------
  resample.data         = cfg_files;
  resample.data.tag     = 'data';
  resample.data.name    = 'Surface Data Files';
  resample.data.filter  = 'any';
  resample.data.ufilter = '^[rl]h.(?!cent|sphe|defe).*';
  resample.data.num     = [1 Inf];
  resample.data.help    = {'Select surfaces files for resampling to template space.'};

  resample.main         = cfg_exbranch;
  resample.main.tag     = 'surfresamp';
  resample.main.name    = 'Resample Surface Data';
  resample.main.val     = {
    resample.data
    };
  resample.main.vfiles  = @vfiles_resample;
  resample.main.prog    = @vbm_surf_resample;
  resample.main.help    = {
    'In order to analyze surface parameters all data have to be resampled into template space and the resampled data have to be finally smoothed. Resampling is done using the warped coordinates of the resp. sphere.'};





%% Toolset
%-----------------------------------------------------------------------
  
  stoolsexp = cfg_choice;
  stoolsexp.name   = 'Surface Expert Tools';
  stoolsexp.tag    = 'stoolsexp';
  stoolsexp.values = {...
    resample.main, ...
    smooth.main, ...
    avg.main, ...
    };

return



%% Result files
%_______________________________________________________________________
function vf = vfiles_smooth(job)
  vf    = job.data; 
  sinfo = vbm_surf_info(job.data);
  for i=1:numel(vf)
    vf(i) = vbm_surf_rename(sinfo(i),'dataname',sprintf('s%d%s',job.fwhm,sinfo(i).dataname));
  end
return;
function vf = vfiles_resample(job)
  vf    = job.data; 
  sinfo = vbm_surf_info(job.data);
  for i=1:numel(vf)
    sinfo(i).resampled = 1; 
    vf(i) = vbm_surf_rename(sinfo(i));
  end
return;
function vf = vfiles_avg(job)
  if isempty(job.outdir{1})
    outdir = spm_fileparts(job.data{1});
  else
    outdir = job.outdir{1};
  end
  
  if job.surfside==2
    side = {'lh','rh'}; 
  else
    side = {};
    sinfo = vbm_surf_info(job.data);
    if ~isempty(strfind([sinfo.side],'lh')), side{end+1}= 'lh'; end
    if ~isempty(strfind([sinfo.side],'rh')), side{end+1}= 'rh'; end
  end
  vf    = cell(numel(side),numel(job.meshsmooth));
  for si = 1:numel(side)
    for smi=1:numel(job.meshsmooth)
      if job.meshsmooth(smi)>0
        vf{si,smi} = fullfile(outdir,sprintf('%s.%s_%dmm.gii',side{si},job.surfname,job.meshsmooth(smi)));
      else
        vf{si,smi} = fullfile(outdir,sprintf('%s.%s.gii',side{si},job.surfname));
      end   
    end
  end
 
return;

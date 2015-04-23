function stoolsexp = cg_vbm_stoolsexp
%_______________________________________________________________________
% wrapper for calling VBM surface utilities
%_______________________________________________________________________
% Robert Dahnke and Christian Gaser
% $Id: cg_vbm_stools.m 689 2015-04-19 15:37:49Z dahnke $
%_______________________________________________________________________



%% average surface mesh
%-----------------------------------------------------------------------

  avg.data_surf         = cfg_files;
  avg.data_surf.tag     = 'data_surf';
  avg.data_surf.name    = 'Sample';
  avg.data_surf.filter  = 'gifti';
  avg.data_surf.ufilter = 'resampled';
  avg.data_surf.num     = [1 Inf];
  avg.data_surf.help    = {
    'Select surfaces.'
    };

  avg.surfsmooth         = cfg_entry;
  avg.surfsmooth.tag     = 'surfsmooth';
  avg.surfsmooth.name    = 'Surface smoothing iterations';
  avg.surfsmooth.strtype = 'r';
  avg.surfsmooth.num     = [1 Inf];
  avg.surfsmooth.val     = {[0 2 32]};
  avg.surfsmooth.help    = {
    'Smoothing of the average surface. '
    ''
    };

  avg.surfside         = cfg_menu;
  avg.surfside.tag     = 'surfside';
  avg.surfside.name    = 'Side handling';
  avg.surfside.labels  = {'separate','mirror'};
  avg.surfside.values  = {1,2};
  avg.surfside.val     = {1};
  avg.surfside.help    = {
    'Handling of the cortical hemispheres.'
    ''
    };
 
  avg.surfname         = cfg_entry;
  avg.surfname.tag     = 'surfsurfname';
  avg.surfname.name    = 'Surface Filename';
  avg.surfname.strtype = 's';
  avg.surfname.num     = [1 Inf];
  avg.surfname.val     = {'average'};
  avg.surfname.help    = {'Name of the surface.'};

  avg.outdir         = cfg_files;
  avg.outdir.tag     = 'outdir';
  avg.outdir.name    = 'Output directory';
  avg.outdir.filter  = 'dir';
  avg.outdir.ufilter = '.*';
  avg.outdir.num     = [1 1];
  avg.outdir.dir     = fullfile(spm('dir'),'toolbox','vbm12');
  avg.outdir.help    = {'Select a directory where files are written.'};

  avg.avg_surf      = cfg_exbranch;
  avg.avg_surf.tag  = 'avg_surf';
  avg.avg_surf.name = 'Average surface mesh';
  avg.avg_surf.val  = {
    avg.data_surf ...
    avg.surfsmooth ...
    avg.surfside ...
    avg.surfname ...
    avg.outdir ...
    };
  avg.avg_surf.prog = @vbm_surf_display; %@vbm_surf_avg;
  avg.avg_surf.help = {
    'Averaging of cortical surfaces.'
    ''
    };

  
  
%% data smoothing
%-----------------------------------------------------------------------
  data_smooth         = cfg_files;
  data_smooth.tag     = 'data_smooth';
  data_smooth.name    = 'Sample';
  data_smooth.filter  = 'any';
  data_smooth.ufilter = '[rl]h.(?!cent|sphe|defe).*';
  data_smooth.num     = [1 Inf];
  data_smooth.help    = {'Select surface data (texture) files for smoothing.'};
  
  fwhm_smooth         = cfg_entry;
  fwhm_smooth.tag     = 'fwhm';
  fwhm_smooth.name    = 'Smoothing filter size in fwhm';
  fwhm_smooth.strtype = 'r';
  fwhm_smooth.num     = [1 1];
  fwhm_smooth.val     = {15};
  fwhm_smooth.help    = {
    'Select filter size for smoothing. For cortical thickness a good starting value is 15mm, while other surface parameters based on cortex folding (e.g. gyrification, cortical complexity) need a larger filter size of about 25mm.'};
 
  datasmooth      = cfg_exbranch;
  datasmooth.tag  = 'datasmooth';
  datasmooth.name = 'Smooth surface data';
  datasmooth.val  = {
    data_smooth ...
    fwhm_smooth ...
  };
  datasmooth.vfiles = @vfiles_datasmooth;
  datasmooth.prog = @vbm_surf_smooth;
  datasmooth.help = {
    'Gaussian smoothing of surface data (texture).'
    ''
  }; 




%% resample surface (mesh and data)
%-----------------------------------------------------------------------
  data_surfdata         = cfg_files;
  data_surfdata.tag     = 'data_surf';
  data_surfdata.name    = 'Surfaces parameters';
  data_surfdata.filter  = 'any';
  data_surfdata.ufilter = '^[lr]h.';
  data_surfdata.num     = [1 Inf];
  data_surfdata.help    = {'Select surfaces parameter files for resampling to template space.'};

  resample_data      = cfg_exbranch;
  resample_data.tag  = 'surfresamp';
  resample_data.name = 'Resample surface parameters';
  resample_data.val  = {data_surfdata};
  resample_data.prog = @vbm_surf_display; %@vbm_surf_resample;
  resample_data.help = {
    'In order to analyze surface parameters all data have to be rsampled into template space and the rsampled data have to be finally smoothed. Resampling is done using the warped coordinates of the resp. sphere.'};

  



%% Toolset
%-----------------------------------------------------------------------
  
  stoolsexp = cfg_choice;
  stoolsexp.name   = 'Surface Expert Tools';
  stoolsexp.tag    = 'stoolsexp';
  stoolsexp.values = {...
    resample_data, ...
    datasmooth, ...
    avg.avg_surf, ...
    };

return



%% Result files
%_______________________________________________________________________
function vf = vfiles_datasmooth(job)
  vf = job.data_smooth;
  for i=1:size(job.data_smooth,1),
      [pth,nam,ext] = spm_fileparts(job.data_smooth{i});
      vf{i} = fullfile(pth,sprintf('s%d.%s%s%s',job.fhwm,nam,ext));
  end;
return;
function vf = vfiles_resample_data(job)
  vf = job.data_smooth;
  for i=1:size(job.data_smooth,1),
      [pth,nam,ext] = spm_fileparts(job.data_smooth{i});
      vf{i} = fullfile(pth,sprintf('s%d.%s%s%s',job.fhwm,nam,ext));
  end;
return;

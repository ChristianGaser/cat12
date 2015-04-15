function stools = cg_vbm_stools
% wrapper for calling VBM utilities
%
%_______________________________________________________________________
% Christian Gaser
% $Id$

%_______________________________________________________________________




%% common fields
%-----------------------------------------------------------------------  
  outdir         = cfg_files;
  outdir.tag     = 'outdir';
  outdir.name    = 'Output directory';
  outdir.filter  = 'dir';
  outdir.ufilter = '.*';
  outdir.num     = [1 1];
  outdir.help    = {'Select a directory where files are written.'};

  surfname         = cfg_entry;
  surfname.tag     = 'surfname';
  surfname.name    = 'Surface name';
  surfname.strtype = 's';
  surfname.num     = [1 Inf];
  surfname.val     = {'average'};
  surfname.help    = {'Name of the surface.'};
  
  datafieldname         = cfg_entry;
  datafieldname.tag     = 'datafieldname';
  datafieldname.name    = 'Datafield name';
  datafieldname.strtype = 's';
  datafieldname.num     = [1 Inf];
  datafieldname.val     = {'intensity'};
  datafieldname.help    = {'Name of the extracted data.'};
 
  fwhm         = cfg_entry;
  fwhm.tag     = 'fwhm';
  fwhm.name    = 'Smoothing filter size in fwhm';
  fwhm.strtype = 'r';
  fwhm.num     = [1 1];
  fwhm.val     = {15};
  fwhm.help    = {
    'Select filter size for smoothing. For cortical thickness a good starting value is 15mm, while other surface parameters based on cortex folding (e.g. gyrification, cortical complexity) need a larger filter size of about 25mm.'};
 

  
  
%% Surface covariance and quality assurance
%-----------------------------------------------------------------------
  data_surf_cov         = cfg_files;
  data_surf_cov.tag     = 'data_surf';
  data_surf_cov.name    = 'Sample';
  data_surf_cov.filter  = 'gifti';
  data_surf_cov.ufilter = 'resampled';
  data_surf_cov.num     = [1 Inf];
  data_surf_cov.help    = {'Select resampled surfaces parameter files.'};

  sample_cov         = cfg_repeat;
  sample_cov.tag     = 'sample';
  sample_cov.name    = 'Data';
  sample_cov.values  = {data_surf_cov};
  sample_cov.num     = [1 Inf];
  sample_cov.help = {...
  'Specify data for each sample. If you specify different samples the mean correlation is displayed in seperate boxplots for each sample.'};

  c         = cfg_entry;
  c.tag     = 'c';
  c.name    = 'Vector';
  c.help    = {'Vector of nuisance values'};
  c.strtype = 'r';
  c.num     = [Inf 1];

  transform         = cfg_repeat;
  transform.tag     = 'transform';
  transform.name    = 'Nuisance variable';
  transform.values  = {c};
  transform.num     = [0 Inf];
  transform.help    = {'This option allows for the specification of nuisance effects to be removed from the data. A potential nuisance parameter can be age. In this case the variance explained by age will be removed prior to the calculation of the correlation.'};

  data_xml = cfg_files;
  data_xml.name = 'XML files';
  data_xml.tag  = 'data_xml';
  data_xml.filter = 'xml';
  data_xml.num  = [1 Inf];
  data_xml.help   = {
  'These are the xml-files that are saved during segmentation. Please note, that the order of the xml-files must be the same as the other data files..'};

  qam         = cfg_repeat;
  qam.tag     = 'qam';
  qam.name    = 'Load quality measures';
  qam.values  = {data_xml};
  qam.num     = [0 Inf];
  qam.help    = {'This option allows to also load the quality measures that are saved in the xml-files. Please note, that the order of the xml-files must be the same as the other data files.'};

  check_mesh_cov      = cfg_exbranch;
  check_mesh_cov.tag  = 'check_mesh_cov';
  check_mesh_cov.name = 'Check sample homogeneity of surfaces';
  check_mesh_cov.val  = {sample_cov,qam,transform};
  check_mesh_cov.prog = @cg_check_cov;
  check_mesh_cov.help = {
  'If you have a reasonable sample size artefacts are easily overseen. In order to identify surfaces with poor image quality or even artefacts you can use this function. Surfaces have to be rsampled to the template space (e.g. normalized images). The idea of this tool is to check the correlation of all files across the sample.'
  ''
  'The correlation is calculated between all images and the mean for each image is plotted using a boxplot and the indicated filenames. The smaller the mean correlation the more deviant is this surface from the sample mean. In the plot outliers from the sample are usually isolated from the majority of images which are clustered around the sample mean. The mean correlation is plotted at the y-axis and the x-axis reflects the image order. Images are plotted from left to right which is helpful if you have selected the images in the order of different sub-groups.'};

    
  
  
%% surface measures
%-----------------------------------------------------------------------  
  data_surf_extract         = cfg_files;
  data_surf_extract.tag     = 'data_surf';
  data_surf_extract.name    = 'Sample';
  data_surf_extract.filter  = 'gifti';
  data_surf_extract.ufilter = '^[lr]h.central';
  data_surf_extract.num     = [1 Inf];
  data_surf_extract.help    = {'Select surfaces to extract values.'};
  
  GI        = cfg_menu;
  GI.name   = 'Gyrification index';
  GI.tag    = 'GI';
  GI.labels = {'none','yes'};
  GI.values = {0,1};
  GI.val    = {1};
  GI.help   = {
    'Extract gyrification index (GI) based on absolute mean curvature. The method is described in Luders et al. NeuroImage, 29: 1224-1230, 2006.'
  };

  FD        = cfg_menu;
  FD.name   = 'Cortical complexity (fractal dimension)';
  FD.tag    = 'FD';
  FD.labels = {'none','yes'};
  FD.values = {0,1};
  FD.val    = {0};
  FD.help   = {
    'Extract Cortical complexity (fractal dimension) which is described in Yotter et al. Neuroimage, 56(3): 961-973, 2011.'
    ''
    'Warning: Estimation of cortical complexity is very slow!'
    ''
  };

  SD        = cfg_menu;
  SD.name   = 'Sulcus depth';
  SD.tag    = 'SD';
  SD.labels = {'none','yes'};
  SD.values = {0,1};
  SD.val    = {1};
  SD.help   = {
    'Extract log10-transformed sulcus depth based on the euclidian distance between the central surface and its convex hull.'
    ''
    'Log-transformation is used to render the data more normally distributed.'
    ''
  };

  SA        = cfg_menu;
  SA.name   = 'Surface area';
  SA.tag    = 'SA';
  SA.labels = {'none','yes'};
  SA.values = {0,1};
  SA.val    = {1};
  SA.help   = {
    'Extract log10-transformed local surface area using re-parameterized tetrahedral surface. The method is described in Winkler et al. NeuroImage, 61: 1428–1443, 2012.'
    ''
    'Log-transformation is used to render the data more normally distributed.'
    ''
  };

  surfextract      = cfg_exbranch;
  surfextract.tag  = 'surfextract';
  surfextract.name = 'Extract surface parameters';
  surfextract.val  = {data_surf_extract,GI,FD,SD};
  surfextract.prog = @vbm_surf_parameters;
  surfextract.help = {'Using this option several surface parameters can be extracted that can be further analyzed.'};

  
  

%% data smoothing
%-----------------------------------------------------------------------
  data_smooth         = cfg_files;
  data_smooth.tag     = 'data_smooth';
  data_smooth.name    = 'Sample';
  data_smooth.filter  = 'any';
  data_smooth.ufilter = '(?=thickness|gyri|frac|logs)';
  data_smooth.num     = [1 Inf];
  data_smooth.help    = {'Select surface data files for smoothing.'};
  
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
  datasmooth.prog = @vbm_surf_display; %@vbm_surf_smooth;
  datasmooth.help = {
    'Gaussian smoothing of surface data surfaces.'
    ''
  }; 


% extract volumetric data
%-----------------------------------------------------------------------  
  v2s.res         = cfg_entry;
  v2s.res.tag     = 'res';
  v2s.res.name    = 'Sampling Size';
  v2s.res.strtype = 'r';
  v2s.res.num     = [1 1];
  v2s.res.val     = {0.5};
  v2s.res.help    = {
    'Resolution of grid along normals [mm]. '
    ''
  };

  v2s.origin         = cfg_entry;
  v2s.origin.tag     = 'origin';
  v2s.origin.name    = 'Origin';
  v2s.origin.strtype = 'r';
  v2s.origin.num     = [1 1];
  v2s.origin.val     = {-1};
  v2s.origin.help    = {
    'Origin (start point) of grid along normals [mm]. Give negative values for origin outside the surface.'
    ''
  };

  v2s.length         = cfg_entry;
  v2s.length.tag     = 'length';
  v2s.length.name    = 'Length';
  v2s.length.strtype = 'r';
  v2s.length.num     = [1 1];
  v2s.length.val     = {15};
  v2s.length.help    = {
    'Lenght of grid along normals [mm].'
    ''
  };

  v2s.interp         = cfg_menu;
  v2s.interp.tag     = 'interp';
  v2s.interp.name    = 'Interpolation Type';
  v2s.interp.labels  = {'nearest','linear','cubic'};
  v2s.interp.values  = {{'nearest_neighbour'},{'linear'},{'cubic'}};
  v2s.interp.val     = {{'nearest_neighbour'}};
  v2s.interp.help    = {
    'Volume extration interpolation type. '
    ' -linear:            Use linear interpolation.'
    ' -nearest_neighbour: Use nearest neighbour interpolation (Default).'
    ' -cubic:             Use cubic interpolation.'
    ''
  };

  v2s.mapping         = cfg_menu;
  v2s.mapping.tag     = 'func';
  v2s.mapping.name    = 'Mapping function';
  v2s.mapping.labels  = {'mean','min','max','exp'}; %,'range'
  v2s.mapping.values  = {{'average'},{'min'},{'max'},{'exp'}}; %,'range'
  v2s.mapping.val     = {{'average'}};
  v2s.mapping.help    = {
    'Volume extration interpolation type. '
    ' -average:  Use average for mapping along normals.'
    ...' -range:    Count number of values in range for mapping along normals.'
    ...'            If any value is out of range values will be counted only until this point'
    ...'            Default value: 3.40282e+38 3.40282e+38'
    ' -max:      Use maximum value for mapping along normals (Default). '
    ...'            Optionally a 2nd volume can be defined to output its value at the maximum value of the 1st volume.'
    ' -min:      Use minimum value for mapping along normals. '
    ...'            Optionally a 2nd volume can be defined to output its value at the minimum value of the 1st volume.'
    ' -exp:      Use exponential average of values for mapping along normals.'
    ...'            The argument defines the distance in mm where values are decayed to 50% '...
    ...'            (recommended value is 10mm).'
    ...'            Default value: 3.40282e+38'
    ...' -sum:               Use sum of values for mapping along normals.'
    '' ...
  };

  % range +2 value, min/max +maskvol, exp +1 value
  % further filter option by thickness and sulcus/gyrus witdh
  

%% extract volumetric data in individual space
%-----------------------------------------------------------------------  
  data_surf_sub         = cfg_files;
  data_surf_sub.tag     = 'data_mesh';
  data_surf_sub.name    = 'Sample';
  data_surf_sub.filter  = 'gifti';
  data_surf_sub.ufilter = '^[rl]h.central.*';
  data_surf_sub.num     = [1 Inf];
  data_surf_sub.help    = {'Select subject surface files.'};
  
  data_sub         = cfg_files; 
  data_sub.tag     = 'data_vol';
  data_sub.name    = 'Volumes';
  data_sub.filter  = 'image';
  data_sub.ufilter = '^(?!wmr|wp|w0rp|wc).*'; % no normalized images
  data_sub.num     = [1 Inf];
  data_sub.help    = {
    'Select subject space volumes.'
  };

  vol2surf      = cfg_exbranch;
  vol2surf.tag  = 'vol2surf';
  vol2surf.name = 'Extract Volume Data by Individual Surface';
  vol2surf.val  = {
    data_surf_sub ...
    data_sub ...
    datafieldname ...
  };
  vol2surf.prog = @vbm_surf_display; %@vbm_surf_vol2surf;
  vol2surf.help = {
    'Extract volumetric data from individual space.'
    ''
  };

  vol2surfexp      = cfg_exbranch;
  vol2surfexp.tag  = 'vol2surfexp';
  vol2surfexp.name = 'Extract Volume Data by Individual Surface (expert)';
  vol2surfexp.val  = {
    data_surf_sub ...
    data_sub ...
    datafieldname ...
    v2s.res v2s.origin v2s.length ...
    v2s.interp ...
    v2s.mapping ...
  };
  vol2surfexp.prog = @vbm_surf_display; %@vbm_surf_vol2surf;
  vol2surfexp.help = {
    'Extract volumetric data from individual space.'
    ''
  };


%% extract volumetric data in template space
%-----------------------------------------------------------------------  
  data_surf_avg         = cfg_files; 
  data_surf_avg.tag     = 'data_mesh';
  data_surf_avg.name    = 'Sample';
  data_surf_avg.filter  = 'gifti';
  data_surf_avg.ufilter = '^[rl]h.*';
  data_surf_avg.num     = [1 2];
  data_surf_avg.dir     = fullfile(spm('dir'),'toolbox','vbm12');
  data_surf_avg.help    = {'Select template surface files.'};

  data_norm         = cfg_files; 
  data_norm.tag     = 'data_vol';
  data_norm.name    = 'Volumes';
  data_norm.filter  = 'image';
  data_norm.ufilter = '^(?=wm|wp|w0rp|wc).*'; % only normalized images
  data_norm.num     = [1 Inf];
  data_norm.help    = {
    'Select normalized space volumes.'
  };

  vol2tempsurf      = cfg_exbranch;
  vol2tempsurf.tag  = 'vol2surf';
  vol2tempsurf.name = 'Extract Volume Data by Template Surface';
  vol2tempsurf.val  = {
    data_surf_avg ...
    data_norm ...
    datafieldname ...
  };
  vol2tempsurf.prog = @vbm_surf_display; %@vbm_surf_vol2surf;
  vol2tempsurf.help = {
    'Extract volumetric data from template/normalized space.'
    ''
  };

  vol2tempsurfexp      = cfg_exbranch;
  vol2tempsurfexp.tag  = 'vol2tempsurfexp';
  vol2tempsurfexp.name = 'Extract Volume Data by Template Surface (expert)';
  vol2tempsurfexp.val  = {
    data_surf_avg ...
    data_norm ...
    datafieldname ...
    v2s.res v2s.origin v2s.length ...
    v2s.interp ...
    v2s.mapping ...
  };
  vol2tempsurfexp.prog = @vbm_surf_display; %@vbm_surf_vol2surf;
  vol2tempsurfexp.help = {
    'Extract volumetric data from template space.'
    'The template surface was generated by VBM12 surface processing [1] of the average wmr*.nii images of VBM12.   '
    ''
    '[1] Dahnke, R., Yotter, R. A., and Gaser, C. 2012.'
    'Cortical thickness and central surface estimation.'
    'Neuroimage, 65C:336?348.'
    ''
    '  WARNING: ONLY FOR VISUAL PURPOSE! '
    '           ALTHOUGH THERE IS A GOOD OVERLAPE TO THE TEMPLATE,'
    '           THERE IS ONLY A SMALL OVERLAPE TO THE INDIVIDUAL ANATOMY!'
    '           DO NOT USE THE RESULTS FOR FURTHER STATISTIC!'
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

  



%% average surface mesh
%-----------------------------------------------------------------------
  data_surf_avg         = cfg_files;
  data_surf_avg.tag     = 'data_surf';
  data_surf_avg.name    = 'Sample';
  data_surf_avg.filter  = 'gifti';
  data_surf_avg.ufilter = 'resampled';
  data_surf_avg.num     = [1 Inf];
  data_surf_avg.help    = {'Select surfaces.'};


  surfsmooth         = cfg_entry;
  surfsmooth.tag     = 'surfsmooth';
  surfsmooth.name    = 'Surface smoothing iterations';
  surfsmooth.strtype = 'r';
  surfsmooth.num     = [1 Inf];
  surfsmooth.val     = {[0 2 32]};
  surfsmooth.help    = {
    'Smoothing of the average surface. '
    ''
    };

  surfside         = cfg_menu;
  surfside.tag     = 'surfside';
  surfside.name    = 'Side handling';
  surfside.labels  = {'separate','mirror'};
  surfside.values  = {1,2};
  surfside.val     = {1};
  surfside.help    = {
    'Handling of the cortical hemispheres.'
    ''
    };

  outdir_fsavg         = cfg_files;
  outdir_fsavg.tag     = 'outdir';
  outdir_fsavg.name    = 'Output directory';
  outdir_fsavg.filter  = 'dir';
  outdir_fsavg.ufilter = '.*';
  outdir_fsavg.num     = [1 1];
  outdir_fsavg.dir     = fullfile(spm('dir'),'toolbox','vbm12');
  outdir_fsavg.help    = {'Select a directory where files are written.'};
  
  avg_surf      = cfg_exbranch;
  avg_surf.tag  = 'avg_surf';
  avg_surf.name = 'Average surface mesh';
  avg_surf.val  = {
    data_surf_avg ...
    surfsmooth ...
    surfside ...
    surfname ...
    outdir_fsavg ...
  };
  avg_surf.prog = @vbm_surf_display; %@vbm_surf_avg;
  avg_surf.help = {
    'Averaging of cortical surfaces.'
    ''
  };



% surface calculations 
% ----------------------------------------------------------------------
% estimation per person (individual and group sampling):
% g groups with i datafiles and i result datafile
 
  data_surf         = cfg_files;
  data_surf.tag     = 'data_surf';
  data_surf.name    = 'Sample';
  data_surf.filter  = 'gifti';
  data_surf.ufilter = 'resampled';
  data_surf.num     = [1 Inf];
  data_surf.help    = {'Select surfaces.'};
 
  
  sample_calc         = cfg_repeat;
  sample_calc.tag     = 'sample';
  sample_calc.name    = 'Data';
  sample_calc.values  = {data_surf};
  sample_calc.num     = [1 Inf];
  sample_calc.help = {...
  'Specify data for each sample. All sample must have the same size and same entry order.'};

  
  surffunction         = cfg_entry;
  surffunction.tag     = 'surffunction';
  surffunction.name    = 'Function';
  surffunction.strtype = 's';
  surffunction.num     = [1 Inf];
  surffunction.val     = {'mean(S)'};
  surffunction.help    = {
    'Function to evaluate the surfaces.'
    's1 + s2'
    'mean(S)'
  };

  surfcalc      = cfg_exbranch;
  surfcalc.tag  = 'surfcalc';
  surfcalc.name = 'Surface value calculations (multi-mesh).';
  surfcalc.val  = {
    sample_calc ...
    surffunction ...
    surfside ...
    surfname ...
    outdir ...
  };
  surfcalc.prog = @vbm_surf_display; %@vbm_surf_calc;
  surfcalc.help = {
    'Mathematical operations for a set of surface data.'
  };




% surface calculations
% ----------------------------------------------------------------------
% estimation per person (only group sampling):
% i datafiles and 1 datafiles
 
  data_surf         = cfg_files;
  data_surf.tag     = 'data_surf';
  data_surf.name    = 'Sample';
  data_surf.filter  = 'gifti';
  data_surf.ufilter = 'resampled';
  data_surf.num     = [1 Inf];
  data_surf.help    = {'Select surfaces.'};
 
  surffunction         = cfg_entry;
  surffunction.tag     = 'surffunction';
  surffunction.name    = 'Function';
  surffunction.strtype = 's';
  surffunction.num     = [1 Inf];
  surffunction.val     = {'mean(S)'};
  surffunction.help    = {
    'Function to evaluate the surfaces.'
    's1 + s2'
    'mean(S)'
  };

  surfcalcavg      = cfg_exbranch;
  surfcalcavg.tag  = 'vbm_surf_calcavg';
  surfcalcavg.name = 'Surface value calculations (single-mesh).';
  surfcalcavg.val  = {
    data_surf ...
    surffunction ...
    surfside ...
    surfname ...
    outdir ...
  };
  surfcalcavg.prog = @vbm_surf_display; %@vbm_surf_calcavg;
  surfcalcavg.help = {
    'Mathematical operations for surface data.'
  };





%% Resample surfaces 
%-----------------------------------------------------------------------
  data_surf         = cfg_files;
  data_surf.tag     = 'data_surf';
  data_surf.name    = 'Surfaces parameters';
  data_surf.filter  = 'any';
  data_surf.ufilter = '^[lr]h.[tgfl][hyro][irag][cias]';
  data_surf.num     = [1 Inf];
  data_surf.help    = {'Select surfaces parameter files for resampling to template space.'};

  surfresamp      = cfg_exbranch;
  surfresamp.tag  = 'surfresamp';
  surfresamp.name = 'Resample and smooth surface parameters';
  surfresamp.val  = {data_surf,fwhm};
  surfresamp.prog = @vbm_surf_resamp;
  surfresamp.help = {
  'In order to analyze surface parameters all data have to be rsampled into template space and the rsampled data have to be finally smoothed. Resampling is done using the warped coordinates of the resp. sphere.'};




%% Resample surfaces FreeSurfer
%-----------------------------------------------------------------------
  data_fs         = cfg_files;
  data_fs.tag     = 'data_fs';
  data_fs.name    = 'Freesurfer subject directories';
  data_fs.filter  = 'dir';
  data_fs.ufilter = '.*';
  data_fs.num     = [1 Inf];
  data_fs.help    = {'Select subject folders of freesurfer data to rsample thickness data.'};

  surfresamp_fs      = cfg_exbranch;
  surfresamp_fs.tag  = 'surfresamp_fs';
  surfresamp_fs.name = 'Resample and smooth existing freesurfer thickness data';
  surfresamp_fs.val  = {data_fs,fwhm,outdir};
  surfresamp_fs.prog = @vbm_surf_resamp_freesurfer;
  surfresamp_fs.help = {
  'If you have existing freesurfer thickness data this function can be used to resample these data, smooth the resampled data, and convert freesurfer data to gifti format.'};






%% Toolset
%-----------------------------------------------------------------------
  
  stools = cfg_choice;
  stools.name   = 'Surface Tools';
  stools.tag    = 'stools';
  stools.values = {check_mesh_cov,surfextract,surfresamp,surfresamp_fs,...
    resample_data,datasmooth,avg_surf,...
    vol2surf,vol2surfexp,vol2tempsurf,vol2tempsurfexp,surfcalc,surfcalcavg};

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

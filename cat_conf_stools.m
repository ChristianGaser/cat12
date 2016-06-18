function stools = cat_conf_stools(expert)
%_______________________________________________________________________
% wrapper for calling CAT surface utilities
%_______________________________________________________________________
% Robert Dahnke and Christian Gaser
% $Id$
%_______________________________________________________________________

  % try to estimate number of processor cores
  try
    numcores = max(feature('numcores'),1);
  catch
    numcores = 1;
  end
  if ~exist('expert','var')
    expert = 0; % switch to de/activate further GUI options
  end

  % parallelize
  % ____________________________________________________________________
  nproc         = cfg_entry;
  nproc.tag     = 'nproc';
  nproc.name    = 'Split job into separate processes';
  nproc.strtype = 'w';
  nproc.val     = {numcores};
  nproc.num     = [1 1];
  nproc.help    = {
    'In order to use multi-threading the CAT12 segmentation job with multiple subjects can be split into separate processes that run in the background. You can even close Matlab, which will not affect the processes that will run in the background without GUI. If you do not want to run processes in the background then set this value to 0.'
      ''
      'Keep in mind that each process needs about 1.5..2GB of RAM, which should be considered to choose the right number of processes.'
      ''
      'Please further note that no additional modules in the batch can be run except CAT12 segmentation. Any dependencies will be broken for subsequent modules.'
    };
 
  
  % do not process, if result allready exist
  % ____________________________________________________________________
  lazy         = cfg_menu;
  lazy.tag     = 'lazy';
  lazy.name    = 'Lazy processing';
  lazy.labels  = {'yes','no'};
  lazy.values  = {1,0};
  lazy.val     = {0};
  lazy.help    = {
    'Do not process data if the result exist. '
  };


%% Surface correlation and quality check
%-----------------------------------------------------------------------
  data_surf_cov         = cfg_files;
  data_surf_cov.tag     = 'data_surf';
  data_surf_cov.name    = 'Sample';
  data_surf_cov.filter  = 'gifti';
  data_surf_cov.ufilter = 'resampled';
  data_surf_cov.num     = [1 Inf];
  data_surf_cov.help    = {'Select resampled surfaces parameter files.'};
  
  data_xml = cfg_files;
  data_xml.name = 'XML files';
  data_xml.tag  = 'data_xml';
  data_xml.filter = '^cat_.*xml';
  data_xml.num  = [1 Inf];
  data_xml.help   = {
  'These are the xml-files that are saved during segmentation in the report folder. Please note, that the order of the xml-files must be the same as the other data files.'};

  sample_cov         = cfg_repeat;
  sample_cov.tag     = 'sample';
  sample_cov.name    = 'Data';
  sample_cov.values  = {data_surf_cov};
  sample_cov.num     = [1 Inf];
  sample_cov.help = {...
  'Specify data for each sample. If you specify different samples the mean correlation is displayed in separate boxplots for each sample.'};

  qam         = cfg_repeat;
  qam.tag     = 'qam';
  qam.name    = 'Load quality measures';
  qam.values  = {data_xml};
  qam.num     = [0 Inf];
  qam.help    = {'This option allows to also load the quality measures that are saved in the xml-files. Please note, that the order of the xml-files must be the same as the other data files.'};

  c         = cfg_entry;
  c.tag     = 'c';
  c.name    = 'Vector';
  c.help    = {'Vector of nuisance values'};
  c.strtype = 'r';
  c.num     = [Inf 1];

  nuisance         = cfg_repeat;
  nuisance.tag     = 'nuisance';
  nuisance.name    = 'Nuisance variable';
  nuisance.values  = {c};
  nuisance.num     = [0 Inf];
  nuisance.help    = {...
  'This option allows for the specification of nuisance effects to be removed from the data. A potential nuisance parameter can be age. In this case the variance explained by age will be removed prior to the calculation of the correlation.'};

  check_mesh_cov      = cfg_exbranch;
  check_mesh_cov.tag  = 'check_mesh_cov';
  check_mesh_cov.name = 'Check sample homogeneity of surfaces';
  check_mesh_cov.val  = {sample_cov,qam,nuisance};
  check_mesh_cov.prog = @cat_stat_check_cov;
  check_mesh_cov.help = {
  'If you have a reasonable sample size artefacts are easily overseen. In order to identify surfaces with poor image quality or even artefacts you can use this function. Surfaces measures have to be resampled to the template space (e.g. normalized data). The idea of this tool is to check the correlation of all files across the sample.'
  ''
  'The correlation is calculated between all surfaces measures and the mean for each surface measures is plotted using a boxplot and the indicated filenames. The smaller the mean correlation the more deviant is this surface measures from the sample mean. In the plot outliers from the sample are usually isolated from the majority of images which are clustered around the sample mean. The mean correlation is plotted at the y-axis and the x-axis reflects the image order. Images are plotted from left to right which is helpful if you have selected the images in the order of different sub-groups.'};

    
  
  
%% surface measures
%-----------------------------------------------------------------------  
  data_surf_extract         = cfg_files;
  data_surf_extract.tag     = 'data_surf';
  data_surf_extract.name    = 'Central Surfaces';
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


  if expert
    area        = cfg_menu;
    area.name   = 'surface area';
    area.tag    = 'area';
    area.labels = {'none','yes'};
    area.values = {0,1};
    area.val    = {1};
    area.help   = {
      'WARNING: This GI measures is still in development and not varified yet!\n\n It extract the average gyrification index (AGI) as local area relation between the individual central and group average surface (FS average).'
    };
  
    % Two cases ... one with the stardard average surface and one were
    % another average surface metric file "area" of a group can be choosen.
    GIA        = cfg_menu;
    GIA.name   = 'Average gyrification index (only resampled output)';
    GIA.tag    = 'GIA';
    GIA.labels = {'none','yes'};
    GIA.values = {0,1};
    GIA.val    = {1};
    GIA.help   = {
      'WARNING: This GI measures is still in development and not varified yet!\n\n It extract the average gyrification index (AGI) as local area relation between the individual central and group average surface (FS average).'
    };
  
    GII        = cfg_menu;
    GII.name   = 'Inflating gyrification index ';
    GII.tag    = 'GII';
    GII.labels = {'none','yes'};
    GII.values = {0,1};
    GII.val    = {1};
    GII.help   = {
      'WARNING: This GI measures is still in development and not varified yet!\n\n It extract a inflating gyrification index (IGI) as local area relation between the individual central and an inflated version of it.'
    };
  
    GIS        = cfg_menu;
    GIS.name   = 'Spherical gyrification index (only resampled output)';
    GIS.tag    = 'GIS';
    GIS.labels = {'none','yes'};
    GIS.values = {0,1};
    GIS.val    = {1};
    GIS.help   = {
      'WARNING: This GI measures is still in development and not varified yet!\n\n It extract a spherical mapping based gyrification index (SGI) as local area relation between the individual central and the hull surface.'
    };
    
    GIL        = cfg_menu;
    GIL.name   = 'Laplacian gyrification index';
    GIL.tag    = 'GIL';
    GIL.labels = {'none','yes'};
    GIL.values = {0,1};
    GIL.val    = {1};
    GIL.help   = {
      'WARNING: This GI measures is still in development and not varified yet!\n\n Extract Laplacian gyrification index (LGI) as local area relation between the individual central and the hull surface [Dahnke:2010].'
    };
  end

  
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
    'Extract sqrt-transformed sulcus depth based on the euclidean distance between the central surface and its convex hull.'
    ''
    'Transformation with sqrt-function is used to render the data more normally distributed.'
    ''
  };

  SA        = cfg_menu;
  SA.name   = 'Surface area';
  SA.tag    = 'SA';
  SA.labels = {'none','yes'};
  SA.values = {0,1};
  SA.val    = {1};
  SA.help   = {
    'Extract log10-transformed local surface area using re-parameterized tetrahedral surface. The method is described in Winkler et al. NeuroImage, 61: 1428-1443, 2012.'
    ''
    'Log-transformation is used to render the data more normally distributed.'
    ''
  };

  surfextract      = cfg_exbranch;
  surfextract.tag  = 'surfextract';
  surfextract.name = 'Extract additional surface parameters';
  if expert > 1
    surfextract.val  = {data_surf_extract,GI,GIA,GII,GIL,GIS,FD,SD,nproc,lazy};
  else
    surfextract.val  = {data_surf_extract,GI,FD,SD,nproc};
  end
  surfextract.prog = @cat_surf_parameters;
  surfextract.help = {'Additional surface parameters can be extracted that can be used for statistical analysis.'};




%% map volumetric data
%-----------------------------------------------------------------------  
  v2s.datafieldname         = cfg_entry;
  v2s.datafieldname.tag     = 'datafieldname';
  v2s.datafieldname.name    = 'Texture Name';
  v2s.datafieldname.strtype = 's';
  v2s.datafieldname.num     = [1 Inf];
  v2s.datafieldname.val     = {'intensity'};
  v2s.datafieldname.help    = {
    'Name of the texture as part of the filename.'
    ''
    '  [rh|lh].TEXTURENAME[.resampled].subjectname[.gii]' 
    };
 
  v2s.interp         = cfg_menu;
  v2s.interp.tag     = 'interp';
  v2s.interp.name    = 'Interpolation Type';
  v2s.interp.labels  = {'Nearest neighbour','Linear','Cubic'};
  v2s.interp.values  = {{'nearest_neighbour'},{'linear'},{'cubic'}};
  v2s.interp.val     = {{'linear'}};
  v2s.interp.help    = {
    'Volume extration interpolation type. '
    ' -linear:            Use linear interpolation (default).'
    ' -nearest_neighbour: Use nearest neighbour interpolation.'
    ' -cubic:             Use cubic interpolation.'
    ''
  };
  
  % sample function 
  v2s.sample         = cfg_menu;
  v2s.sample.tag     = 'sample';
  v2s.sample.name    = 'Sample Function';
  v2s.sample.labels  = {'Mean','Maximum','Minimum','Absolute maximum'};
  v2s.sample.values  = {{'avg'},{'max'},{'min'},{'maxabs'}};
  v2s.sample.val     = {{'avg'}};
  v2s.sample.help    = {
    'Sample function to combine the values of the grid along the surface normals.'
  };

  %% -- sampling points and average function
 
  % startpoint
  v2s.abs_startpoint         = cfg_entry;
  v2s.abs_startpoint.tag     = 'startpoint';
  v2s.abs_startpoint.name    = 'Startpoint';
  v2s.abs_startpoint.strtype = 'r';
  v2s.abs_startpoint.val     = {-0.5};
  v2s.abs_startpoint.num     = [1 1];
  v2s.abs_startpoint.help    = {
    'Absolut position of the startpoint of the grid along the surface normals in mm. Give negative value for a startpoint outside the surface (CSF direction). '
  };
  v2s.rel_startpoint = v2s.abs_startpoint;
  v2s.rel_startpoint.val     = {0};
  v2s.rel_startpoint.help    = {
    'Relative position of the startpoint of the grid along the surface normals from the center of a tissue class. Give negative value for a startpoint outside the surface (CSF direction). '
  };
  
  % stepsize
  v2s.abs_stepsize         = cfg_entry;
  v2s.abs_stepsize.tag     = 'stepsize';
  v2s.abs_stepsize.name    = 'Stepsize';
  v2s.abs_stepsize.strtype = 'r';
  v2s.abs_stepsize.val     = {0.5};
  v2s.abs_stepsize.num     = [1 1];
  v2s.abs_stepsize.help    = {
    'Absolute stepsize in mm of the grid along the surface beginning from the startpoint. '
  };
  v2s.rel_stepsize = v2s.abs_stepsize; 
  v2s.rel_stepsize.help    = {
    'Relative stepsize based on the thickness/depth of the tissue class of the grid along the surface beginning from the startpoint. '
  };

  % endpoint
  v2s.abs_endpoint         = cfg_entry;
  v2s.abs_endpoint.tag     = 'endpoint';
  v2s.abs_endpoint.name    = 'Endpoint';
  v2s.abs_endpoint.strtype = 'r';
  v2s.abs_endpoint.val     = {+0.5};
  v2s.abs_endpoint.num     = [1 1];
  v2s.abs_endpoint.help    = {
    'Absolut position of the endpoint of the grid along the surface normals in mm. Give negative value for a startpoint outside the surface (CSF direction). '
  };
  v2s.rel_endpoint = v2s.abs_endpoint;
  v2s.rel_endpoint.val     = {1};
  v2s.rel_endpoint.help    = {
    'Relative position of the endpoint of the grid along the surface normals from the center of a tissue class. Give negative value for a startpoint outside the surface (CSF direction). '
  };

  % tissue class
  v2s.abs_class         = cfg_menu;
  v2s.abs_class.tag     = 'class';
  v2s.abs_class.name    = 'Tissue Class';
  v2s.abs_class.labels  = {'GM'};
  v2s.abs_class.values  = {'GM'};
  v2s.abs_class.val     = {1};
  v2s.abs_class.help    = {
    'Tissue class for which the relative positions are estimated.'
  };
  v2s.rel_class = v2s.abs_class; 
  if expert==2
    v2s.abs_class.labels  = {'GM','WM','CSF'};
    v2s.abs_class.values  = {'GM','WM','CSF'};
  end
  
  % absolute position
  v2s.abs_mapping         = cfg_branch;
  v2s.abs_mapping.tag     = 'abs_mapping';
  v2s.abs_mapping.name    = 'Absolution Position From a Tissue Boundary';
  v2s.abs_mapping.val   = {
    v2s.abs_class ...
    v2s.abs_startpoint ...
    v2s.abs_stepsize ...
    v2s.abs_endpoint ...
  }; 
  v2s.tissue.help    = {
    'Map volumetric data from abolute position(s) from a tissue boundary.'
  };
  
  %% relative mapping
  v2s.rel_mapping         = cfg_branch;
  v2s.rel_mapping.tag     = 'rel_mapping';
  v2s.rel_mapping.name    = 'Relative Position Within a Tissue Class';
  v2s.rel_mapping.val   = {
    v2s.rel_class ...
    v2s.rel_startpoint ...
    v2s.rel_stepsize ...
    v2s.rel_endpoint ...
  };
  v2s.rel_mapping.help    = {
    'Map volumetric data from relative positions from the center of a tissue class.'
  };

  %% -- Mapping function

  v2s.mapping         = cfg_choice;
  v2s.mapping.tag     = 'mapping';
  v2s.mapping.name    = 'Mapping Function';
  v2s.mapping.values  = {
    v2s.abs_mapping ...
    v2s.rel_mapping ...
  }; 
  v2s.mapping.help    = {
    'Volume extration type. '
    '  Absolution Position From a Tissue Boundary:'
    '    Extract a set of values around a tissue boundary with a specified absolute sample '
    '    distance and combine these values.'
    '  Relative Position Within a Tissue Class:' 
    '    Extract a set of values around the center of a tissue class with a specified relative sample'
    '    distance and combine these values.'
    '' 
  };
  v2s.mapping.val = {v2s.rel_mapping};



% extract volumetric data in individual space
%-----------------------------------------------------------------------  
  v2s.data_surf_sub_lh         = cfg_files;
  v2s.data_surf_sub_lh.tag     = 'data_mesh_lh';
  v2s.data_surf_sub_lh.name    = '(Left) Individual Surfaces';
  v2s.data_surf_sub_lh.filter  = 'gifti';
  v2s.data_surf_sub_lh.ufilter = '^lh.central.(?!nofix).*';
  v2s.data_surf_sub_lh.num     = [1 Inf];
  v2s.data_surf_sub_lh.help    = {
    'Select left subject surface files (do not select the *.nofix.* surface).'
    'Right side will automatically processed.'
    };
   
  v2s.data_sub         = cfg_files; 
  v2s.data_sub.tag     = 'data_vol';
  v2s.data_sub.name    = 'Volumes in Native Space';
  v2s.data_sub.filter  = 'image';
  v2s.data_sub.ufilter = '^(?!wm|wp|m0wp|mwp|wc).*'; % no normalized images
  v2s.data_sub.num     = [1 Inf];
  v2s.data_sub.help    = {
    'Select volumes in native (subject) space.'
  };

  v2s.vol2surf      = cfg_exbranch;
  v2s.vol2surf.tag  = 'vol2surf';
  v2s.vol2surf.name = 'Map Volume (Native Space) to Individual Surface';
  if expert
    v2s.vol2surf.val = {
      v2s.data_sub ...
      v2s.data_surf_sub_lh ...
      v2s.sample ...
      v2s.interp ...
      v2s.datafieldname ...
      v2s.mapping ...
      };
  else
    v2s.vol2surf.val = {
    v2s.data_sub ...
    v2s.data_surf_sub_lh ...
    v2s.sample ...
    v2s.datafieldname ...
    v2s.mapping ...
    };
  end
  v2s.vol2surf.prog = @cat_surf_vol2surf;
  v2s.vol2surf.help = {
    'Map volume (native space) to individual surface.'
    ''
  };



%% extract volumetric data in template space
%-----------------------------------------------------------------------  
  v2s.data_surf_avg_lh         = cfg_files; 
  v2s.data_surf_avg_lh.tag     = 'data_mesh_lh';
  v2s.data_surf_avg_lh.name    = '(Left) Template Hemisphere';
  v2s.data_surf_avg_lh.filter  = 'gifti';
  v2s.data_surf_avg_lh.ufilter = '^lh.*';
  v2s.data_surf_avg_lh.num     = [1 1];
  v2s.data_surf_avg_lh.val{1}  = {fullfile(spm('dir'),'toolbox','cat12','templates_surfaces','lh.central.Template_T1_IXI555_MNI152.gii')};
  v2s.data_surf_avg_lh.dir     = fullfile(spm('dir'),'toolbox','cat12');
  v2s.data_surf_avg_lh.help    = {
    'Select left template surface file. '
    'Right hemisphere will automatically processed.'
    };
    
  v2s.data_norm         = cfg_files; 
  v2s.data_norm.tag     = 'data_vol';
  v2s.data_norm.name    = 'Spatially Normalized Volumes';
  v2s.data_norm.filter  = 'image';
  v2s.data_norm.ufilter = '.*';
  v2s.data_norm.num     = [1 Inf];
  v2s.data_norm.help    = {
    'Select spatially normalized volumes (in template space). The output file will be named according to the volume names'
    ''
    '  [rh|lh].central.Template_T1_IXI555_MNI152.VOLUMENAME.gii' 
  };

  v2s.vol2tempsurf      = cfg_exbranch;
  v2s.vol2tempsurf.tag  = 'vol2surftemp';
  v2s.vol2tempsurf.name = 'Map Normalized Volume to Template Surface';
  v2s.vol2tempsurf.val  = {
    v2s.data_norm ...
    v2s.data_surf_avg_lh ...
    v2s.sample ...
    v2s.interp ...
    v2s.datafieldname ...
    v2s.mapping ...
  };
  v2s.vol2tempsurf.prog = @cat_surf_vol2surf;
  v2s.vol2tempsurf.help = {
    'Map spatially normalized data (in template space) to template surface.'
    'The template surface was generated by CAT12 surface processing [1] of the average of 555 Dartel-normalized images of the IXI database that were also used to create the IXI Dartel template.   '
    ''
    '  [1] Dahnke, R., Yotter, R. A., and Gaser, C. 2012.'
    '  Cortical thickness and central surface estimation.'
    '  Neuroimage, 65C:336?348.'
    ''
    '  WARNING: This function is primarily thought in order to project statistical results '
    '           to the template surface. Do not use the output of this function'
    '           for any statistical analysis. Statistical analysis should only use'
    '           resampled data of individual volume data (in native space) that were mapped'
    '           to individual surfaces.'
    ''
  };


%% surface to ROI (in template space)
%  ---------------------------------------------------------------------
%  * CxN  cdata files [ thickness , curvature , ... any other file ] in template space [ resampled ]
%  * M    atlas files [ choose files ]
%  * average measures [ mean , std , median , max , min , mean95p ]
% 
%  - csv export (multiple files) > catROIs_ATLAS_SUBJECT
%  - xml export (one file)       > catROIs_ATLAS_SUBJECT
%  ---------------------------------------------------------------------
%  surface ROIs have sidewise index?!
%  ---------------------------------------------------------------------

% set of cdata files
  if expert && 0 % this is not ready now
    s2r.cdata         = cfg_files;
    s2r.cdata.tag     = 'cdata';
    s2r.cdata.name    = '(Left) Surface Data Files';
    s2r.cdata.filter  = 'any';
    s2r.cdata.ufilter = 'lh.(?!cent|sphe|defe).*';
    s2r.cdata.num     = [1 Inf];
    s2r.cdata.help    = {'Surface data sample. Both sides will be processed'};
  else % only smoothed/resampled
    s2r.cdata         = cfg_files;
    s2r.cdata.tag     = 'cdata';
    s2r.cdata.name    = '(Left) Surface Data Files';
    s2r.cdata.filter  = 'any';
    s2r.cdata.ufilter = '^lh.(?!cent|sphe|defe).*';
    s2r.cdata.num     = [1 Inf];
    s2r.cdata.help    = {'Surface data sample. Both sides will be processed'};
  end
  
  s2r.cdata_sample         = cfg_repeat;
  s2r.cdata_sample.tag     = 'cdata_sub.';
  s2r.cdata_sample.name    = 'Surface Data Sample';
  s2r.cdata_sample.values  = {s2r.cdata};
  s2r.cdata_sample.num     = [1 Inf];
  s2r.cdata_sample.help = {[...
    'Specify data for each sample (i.e. thickness, gyrification, ...). ' ...
    'All samples must have the same size and same order. ' ...
    ''
  ]};

% ROI files
  s2r.ROIs         = cfg_files;
  s2r.ROIs.tag     = 'rdata';
  s2r.ROIs.name    = '(Left) ROI atlas files';
  s2r.ROIs.filter  = 'any';
  if expert 
    s2r.ROIs.ufilter = 'lh.aparc.*';
  else
    s2r.ROIs.ufilter = 'lh.aparc_[a2009s|DKT40JT].*'; % not yet working for all atlases
  end
  s2r.ROIs.dir     = fullfile(spm('dir'),'toolbox','cat12','atlases_surfaces'); 
  s2r.ROIs.num     = [1 Inf];
  s2r.ROIs.help    = {'These are the ROI atlas files. Both sides will be processed.'};

% ROI area
  s2r.area         = cfg_menu;
  s2r.area.tag     = 'area';
  s2r.area.name    = 'Estimate ROI Area';
  s2r.area.labels  = {'no','yes'};
  s2r.area.values  = {0,1};
  s2r.area.val     = {1}; 
  s2r.area.help    = {'Estimate area of each ROI.'};
 
% ROI area
  s2r.vernum         = cfg_menu;
  s2r.vernum.tag     = 'vernum';
  s2r.vernum.name    = 'Count ROI Vertices';
  s2r.vernum.labels  = {'no','yes'};
  s2r.vernum.values  = {0,1};
  s2r.vernum.val     = {1}; 
  s2r.vernum.help    = {'Count vertices of each ROI.'};
  
% average mode within a ROI  
  % mean
  s2r.avg.mean         = cfg_menu;
  s2r.avg.mean.tag     = 'mean';
  s2r.avg.mean.name    = 'Mean Estimation';
  s2r.avg.mean.labels  = {'no','yes'};
  s2r.avg.mean.values  = {0,1};
  s2r.avg.mean.val     = {1}; 
  s2r.avg.mean.help    = {'Set mean value estimation per ROI.'}; 
  % std
  s2r.avg.std         = cfg_menu;
  s2r.avg.std.tag     = 'std';
  s2r.avg.std.name    = 'STD Estimation';
  s2r.avg.std.labels  = {'no','yes'};
  s2r.avg.std.values  = {0,1};
  s2r.avg.std.val     = {1}; 
  s2r.avg.std.help    = {'Set standard deviation estimation per ROI.'}; 
  % min
  s2r.avg.min         = cfg_menu;
  s2r.avg.min.tag     = 'min';
  s2r.avg.min.name    = 'Minimum Estimation';
  s2r.avg.min.labels  = {'no','yes'};
  s2r.avg.min.values  = {0,1};
  s2r.avg.min.val     = {0}; 
  s2r.avg.min.help    = {'Set minimum estimation per ROI.'};   
  % max
  s2r.avg.max         = cfg_menu;
  s2r.avg.max.tag     = 'max';
  s2r.avg.max.name    = 'Maximum Estimation';
  s2r.avg.max.labels  = {'no','yes'};
  s2r.avg.max.values  = {0,1};
  s2r.avg.max.val     = {0}; 
  s2r.avg.max.help    = {'Set maximum estimation per ROI.'};   
  % median
  s2r.avg.median         = cfg_menu;
  s2r.avg.median.tag     = 'median';
  s2r.avg.median.name    = 'Median Estimation';
  s2r.avg.median.labels  = {'no','yes'};
  s2r.avg.median.values  = {0,1};
  s2r.avg.median.val     = {0}; 
  s2r.avg.median.help    = {'Set median estimation per ROI.'};   
  % all functions
  s2r.avg.main         = cfg_branch;
  s2r.avg.main.tag     = 'avg';
  s2r.avg.main.name    = 'ROI Average Functions';
  s2r.avg.main.val     = {
    s2r.avg.mean ...
    s2r.avg.std ...
    s2r.avg.min ...
    s2r.avg.max ...
    s2r.avg.median ...
  };

%% main function
  surf2roi      = cfg_exbranch;
  surf2roi.tag  = 'surf2roi';
  surf2roi.name = 'Extract ROI-based surface values';
  switch expert
  case 2
    surf2roi.val  = {
      s2r.cdata_sample ...
      s2r.ROIs ...
      nproc ... 
      ... s2r.vernum ... 
      ... s2r.area ... does not work yet
      s2r.avg.main};
  case 1
    surf2roi.val  = {
      s2r.cdata_sample ...
      s2r.ROIs};
  case 0
    surf2roi.val  = {s2r.cdata_sample};
  end
  surf2roi.prog = @cat_surf_surf2roi;
  surf2roi.help = {
    'While ROI-based values for VBM (volume) data are automatically saved in the label folder as XML file it is necessary to additionally extract these values for surface data. This has to be done after preprocessing the data and creating cortical surfaces. '
    ''
    'You can extract ROI-based values for cortical thickness but also for any other surface parameter that was extracted using the ''Extract Additional Surface Parameters'' function.'
    ''
    'Please note that these values are extracted from data in native space without any smoothing. As default the mean inside a ROI is calculated and saved as XML file in the label folder.'
  };

%% roi to surface  
%  ---------------------------------------------------------------------

% ROI files
  r2s.ROIs         = cfg_files;
  r2s.ROIs.tag     = 'rdata';
  r2s.ROIs.name    = 'ROI atlas files';
  r2s.ROIs.filter  = 'any';
  r2s.ROIs.ufilter = '.*';
  r2s.ROIs.dir     = fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm'); 
  r2s.ROIs.num     = [1 Inf];
  r2s.ROIs.help    = {'These are the indivudal ROI atlas files from the label directory. Choose CSV or XML files.'};

 % atlas used for extraction .. if xml 
  r2s.atlas         = cfg_entry;
  r2s.atlas.tag     = 'atlas';
  r2s.atlas.name    = 'Atlas maps';
  r2s.atlas.help    = {'Enter the name of the atlas maps, e.g. "LPBA40", or "all".'};
  r2s.atlas.strtype = 's';
  r2s.atlas.val     = {'all'};
  r2s.atlas.num     = [1 Inf];
  
% datafield in the CSV or XML file eg the GM thickness
  r2s.data         = cfg_entry;
  r2s.data.tag     = 'fields';
  r2s.data.name    = 'Datafields';
  r2s.data.help    = {'Enter the name of the data fields, e.g. "Vgm", "Igm", "Tgm", "mySurfData" or "all"'};
  r2s.data.strtype = 's';
  r2s.data.val     = {'all'};
  r2s.data.num     = [1 Inf];
  
% surface 
  r2s.surf        = cfg_menu;
  r2s.surf.name   = 'Surface type for mapping';
  r2s.surf.tag    = 'surf';
  r2s.surf.labels = {'FS average','Dartel average','subject'};
  r2s.surf.values = {'freesurfer','dartel','subject'};
  r2s.surf.val    = {'freesurfer'};
  r2s.surf.help   = {'Surface type for value projection. '};
 
% average function? - not required 
 
%% main function
  roi2surf      = cfg_exbranch;
  roi2surf.tag  = 'roi2surf';
  roi2surf.name = 'Map ROI data to the surface';
  roi2surf.val  = {
    r2s.ROIs ...
    r2s.atlas ...
    r2s.data ...
    r2s.surf ...
    };
  roi2surf.prog = @cat_roi_roi2surf;
  roi2surf.help = {
    ''
  };
%% surface calculations 
%  ---------------------------------------------------------------------
%  estimation per subject (individual and group sampling):
%  g groups with i datafiles and i result datafile
 
  sc.cdata         = cfg_files;
  sc.cdata.tag     = 'cdata';
  sc.cdata.name    = 'Surface Data Files';
  sc.cdata.filter  = 'any';
  sc.cdata.ufilter = '[rl]h.*';
  sc.cdata.num     = [1 Inf];
  sc.cdata.help    = {'These are the surface data files that are used by the calculator.  They are referred to as s1, s2, s3, etc in the order they are specified.'};
  
  sc.cdata_sub         = cfg_files;
  sc.cdata_sub.tag     = 'cdata';
  sc.cdata_sub.name    = 'Surface Data Files';
  sc.cdata_sub.filter  = 'any';
  sc.cdata_sub.ufilter = '[rl]h.(?!cent|sphe|defe).*gii';
  sc.cdata_sub.num     = [1 Inf];
  sc.cdata_sub.help    = {'These are the surface data files that are used by the calculator.  They are referred to as s1, s2, s3, etc in the order they are specified.'};
   
  sc.cdata_sample         = cfg_repeat;
  sc.cdata_sample.tag     = 'cdata_sub.';
  sc.cdata_sample.name    = 'Surface Data Sample';
  sc.cdata_sample.values  = {sc.cdata_sub};
  sc.cdata_sample.num     = [1 Inf];
  sc.cdata_sample.help = {...
    'Specify data for each sample. All samples must have the same size and same order.'};
 
  sc.outdir         = cfg_files;
  sc.outdir.tag     = 'outdir';
  sc.outdir.name    = 'Output Directory';
  sc.outdir.filter  = 'dir';
  sc.outdir.ufilter = '.*';
  sc.outdir.num     = [0 1];
  sc.outdir.val{1}  = {''};
  sc.outdir.help    = {
    'Files produced by this function will be written into this output directory.  If no directory is given, images will be written  to current working directory.  If both output filename and output directory contain a directory, then output filename takes precedence.'
  };
  
  sc.surfname         = cfg_entry;
  sc.surfname.tag     = 'dataname';
  sc.surfname.name    = 'Output Filename';
  sc.surfname.strtype = 's';
  sc.surfname.num     = [1 Inf];
  sc.surfname.val     = {'output'};
  sc.surfname.help    = {'The output surface data file is written to current working directory unless a valid full pathname is given.  If a path name is given here, the output directory setting will be ignored.'};
 
  sc.dataname         = cfg_entry;
  sc.dataname.tag     = 'dataname';
  sc.dataname.name    = 'Texture Name';
  sc.dataname.strtype = 's';
  sc.dataname.num     = [1 Inf];
  sc.dataname.val     = {'output'};
  sc.dataname.help    = {
    'Name of the texture as part of the filename.'
    ''
    '  [rh|lh].TEXTURENAME[.resampled].subjectname[.gii]' 
  };

  sc.expression         = cfg_entry;
  sc.expression.tag     = 'expression';
  sc.expression.name    = 'Expression';
  sc.expression.strtype = 's';
  sc.expression.num     = [1 Inf];
  sc.expression.val     = {'s1'};
  sc.expression.help    = {
    'Example expressions (f):'
    '  * Mean of six surface textures (select six texture files)'
    '    f = ''(s1+s2+s3+s4+s5+s6)/6'''
    '  * Make a binaray mask texture at threshold of 100'
    '    f = ''(s1>100)'''
    '  * Make a mask from one texture and apply to another'
    '    f = ''s2.*(s1>100)'''
    '        - here the first texture is used to make the mask, which is applied to the second texture'
    '  * Sum of n texures'
    '    f = ''s1 + s2 + s3 + s4 + s5 + ...'''
    '  * Sum of n texures (when reading data into a data-matrix - use dmtx arg)'
    '    f = mean(S)'
    ''
  };

  sc.dmtx         = cfg_menu;
  sc.dmtx.tag     = 'dmtx';
  sc.dmtx.name    = 'Data Matrix';
  sc.dmtx.labels  = {
    'No - don''t read images into data matrix',...
    'Yes -  read images into data matrix'
  };
  sc.dmtx.values  = {0,1};
  sc.dmtx.val     = {0};
  sc.dmtx.help    = {
    'If the dmtx flag is set, then textures are read into a data matrix S (rather than into separate variables s1, s2, s3,...). The data matrix should be referred to as S, and contains textures in rows. Computation is vertex by vertex, S is a NxK matrix, where N is the number of input textures, and K is the number of vertices per plane.'
  };


% ----------------------------------------------------------------------

  surfcalc      = cfg_exbranch;
  surfcalc.tag  = 'surfcalc';
  surfcalc.name = 'Surface Calculator';
  surfcalc.val  = {
    sc.cdata ...
    sc.surfname ...
    sc.outdir ...
    sc.expression ...
    sc.dmtx ...
  };
  surfcalc.prog = @cat_surf_calc;
  surfcalc.help = {
    'Mathematical operations for surface data (textures).'
    'It works similar to ''spm_imcalc''.  The input surface data must have the same number of entries (e.g. data of the same hemisphere of a subject or resampled data).'
  };


  surfcalcsub      = cfg_exbranch;
  surfcalcsub.tag  = 'surfcalcsub';
  surfcalcsub.name = 'Surface Calculator (subject-wise)';
  surfcalcsub.val  = {
    sc.cdata_sample ...
    sc.dataname ...
    sc.outdir ...
    sc.expression ...
    sc.dmtx ...
  };
  surfcalcsub.prog = @cat_surf_calc;
  surfcalcsub.help = {
    'Mathematical operations for surface data sets (textures).'
    'In contrast to the ''Surface Calculator'' it allows to apply the same expression to multiple subjects. Please note that a fixed name structure is expected: [rh|lh].TEXTURENAME[.resampled].subjectname[.gii]. Here TEXTURENAME will be replaced by the output name.'
  };




%% Resample and smooth surfaces 
%-----------------------------------------------------------------------
  data_surf         = cfg_files;
  data_surf.tag     = 'data_surf';
  data_surf.name    = 'Surfaces Data';
  data_surf.filter  = 'any';
  if expert > 1
    data_surf.ufilter = '^[lr]h.';
  else
    data_surf.ufilter = '[lr]h.(?!cent|sphe|defe).*';
  end
  data_surf.num     = [1 Inf];
  data_surf.help    = {'Select Surfaces Data Files for Resampling to Template Space.'};

  fwhm         = cfg_entry;
  fwhm.tag     = 'fwhm';
  fwhm.name    = 'Smoothing Filter Size in FWHM';
  fwhm.strtype = 'r';
  fwhm.num     = [1 1];
  fwhm.val     = {15};
  fwhm.help    = {
    'Select filter size for smoothing. For cortical thickness a good starting value is 15mm, while other surface parameters based on cortex folding (e.g. gyrification, cortical complexity) need a larger filter size of about 25mm. For no filtering use a value of 0.'};

  surfresamp      = cfg_exbranch;
  surfresamp.tag  = 'surfresamp';
  surfresamp.name = 'Resample and Smooth Surface Data';
  if expert > 1
    surfresamp.val  = {data_surf,fwhm,nproc,lazy};
  else
    surfresamp.val  = {data_surf,fwhm,nproc};
  end
  surfresamp.prog = @cat_surf_resamp;
  surfresamp.help = {
  'In order to analyze surface parameters all data have to be resampled into template space and the resampled data have to be finally smoothed. Resampling is done using the warped coordinates of the resp. sphere.'};




%% Resample and smooth FreeSurfer surfaces
%-----------------------------------------------------------------------
  data_fs         = cfg_files;
  data_fs.tag     = 'data_fs';
  data_fs.name    = 'Freesurfer Subject Directories';
  data_fs.filter  = 'dir';
  data_fs.ufilter = '.*';
  data_fs.num     = [1 Inf];
  data_fs.help    = {'Select subject folders of freesurfer data to resample thickness data.'};

  outdir         = cfg_files;
  outdir.tag     = 'outdir';
  outdir.name    = 'Output Directory';
  outdir.filter  = 'dir';
  outdir.ufilter = '.*';
  outdir.num     = [0 1];
  outdir.help    = {'Select a directory where files are written.'};

  surfresamp_fs      = cfg_exbranch;
  surfresamp_fs.tag  = 'surfresamp_fs';
  surfresamp_fs.name = 'Resample and Smooth Existing FreeSurfer Thickness Data';
  surfresamp_fs.val  = {data_fs,fwhm,outdir};
  surfresamp_fs.prog = @cat_surf_resamp_freesurfer;
  surfresamp_fs.help = {
  'If you have existing freesurfer thickness data this function can be used to resample these data, smooth the resampled data, and convert freesurfer data to gifti format.'};

%% Flipsides
%-----------------------------------------------------------------------
  flip.cdata         = cfg_files;
  flip.cdata.tag     = 'cdata';
  flip.cdata.name    = 'Surface Data Files';
  flip.cdata.filter  = 'any';
  flip.cdata.ufilter = '^s.*mm\.lh.*';
  flip.cdata.num     = [1 Inf];
  flip.cdata.help    = {'Texture maps that should be flipped/mirrored from right to left.'};
  

  flipsides      = cfg_exbranch;
  flipsides.tag  = 'flipsides';
  flipsides.name = 'Flip right to left hemisphere';
  flipsides.val  = {flip.cdata};
  flipsides.prog = @cat_surf_flipsides;
  flipsides.help = {
  'This function flip the right hemisphere to the left side, to allow side'
  ''
  ''};


%% Toolset
%-----------------------------------------------------------------------
  
  stools = cfg_choice;
  stools.name   = 'Surface Tools';
  stools.tag    = 'stools';
  if expert==2
    stools.values = { ...
      check_mesh_cov, ...
      surfextract, ...
      surfresamp, ...
      surfresamp_fs,...
      v2s.vol2surf, ...
      v2s.vol2tempsurf, ...
      surfcalc, ...
      surfcalcsub, ...
      surf2roi, ...
      roi2surf, ...
      flipsides, ...
      ... roicalc, ...
      };    
  elseif expert==1
    stools.values = { ...
      check_mesh_cov, ...
      surfextract, ...
      surfresamp, ...
      surfresamp_fs,...
      v2s.vol2surf, ...
      v2s.vol2tempsurf, ...
      surfcalc, ...
      surfcalcsub, ...
      surf2roi, ...
      ... roi2surf, ...
      flipsides, ...
      ... roicalc, ...
      };
  else
    stools.values = { ...
      check_mesh_cov, ...
      surfextract, ...
      surfresamp, ...
      surfresamp_fs,...
      surf2roi, ...
      v2s.vol2surf, ...
      v2s.vol2tempsurf, ...
      surfcalc, ...
      surfcalcsub, ...
      };
  end
return


%% Result files
%_______________________________________________________________________

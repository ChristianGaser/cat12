function stools = cg_vbm_stools(expert)
%_______________________________________________________________________
% wrapper for calling VBM surface utilities
%_______________________________________________________________________
% Robert Dahnke and Christian Gaser
% $Id$
%_______________________________________________________________________



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
  data_xml.filter = 'xml';
  data_xml.num  = [1 Inf];
  data_xml.help   = {
  'These are the xml-files that are saved during segmentation. Please note, that the order of the xml-files must be the same as the other data files..'};

  sample_cov         = cfg_repeat;
  sample_cov.tag     = 'sample';
  sample_cov.name    = 'Data';
  sample_cov.values  = {data_surf_cov};
  sample_cov.num     = [1 Inf];
  sample_cov.help = {...
  'Specify data for each sample. If you specify different samples the mean correlation is displayed in seperate boxplots for each sample.'};

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

  transform         = cfg_repeat;
  transform.tag     = 'transform';
  transform.name    = 'Nuisance variable';
  transform.values  = {c};
  transform.num     = [0 Inf];
  transform.help    = {...
  'This option allows for the specification of nuisance effects to be removed from the data. A potential nuisance parameter can be age. In this case the variance ',...
  'explained by age will be removed prior to the calculation of the correlation.'};

  check_mesh_cov      = cfg_exbranch;
  check_mesh_cov.tag  = 'check_mesh_cov';
  check_mesh_cov.name = 'Check sample homogeneity of surfaces';
  check_mesh_cov.val  = {sample_cov,qam,transform};
  check_mesh_cov.prog = @cg_check_cov;
  check_mesh_cov.help = {
  'If you have a reasonable sample size artefacts are easily overseen. In order to identify surfaces with poor image quality or even artefacts you can use this function. Surfaces measures have to be resampled to the template space (e.g. normalized data). The idea of this tool is to check the correlation of all files across the sample.'
  ''
  'The correlation is calculated between all surfaces measures and the mean for each surface measures is plotted using a boxplot and the indicated filenames. The smaller the mean correlation the more deviant is this surface measures from the sample mean. In the plot outliers from the sample are usually isolated from the majority of images which are clustered around the sample mean. The mean correlation is plotted at the y-axis and the x-axis reflects the image order. Images are plotted from left to right which is helpful if you have selected the images in the order of different sub-groups.'};

    
  
  
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
  surfextract.name = 'Extract additional surface parameters';
  surfextract.val  = {data_surf_extract,GI,FD,SD};
  surfextract.prog = @vbm_surf_parameters;
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
    'Name of the surface data part of the filename.'
    ''
    '  [s*mm.][rh|lh].DATANAME[.resampled|].subjectname[.gii]' 
    };
 
  v2s.interp         = cfg_menu;
  v2s.interp.tag     = 'interp';
  v2s.interp.name    = 'Interpolation Type';
  v2s.interp.labels  = {'nearest','linear','cubic'};
  v2s.interp.values  = {{'nearest_neighbour'},{'linear'},{'cubic'}};
  v2s.interp.val     = {{'linear'}};
  v2s.interp.help    = {
    'Volume extration interpolation type. '
    ' -linear:            Use linear interpolation.'
    ' -nearest_neighbour: Use nearest neighbour interpolation (Default).'
    ' -cubic:             Use cubic interpolation.'
    ''
  };
  
  
  %% -- absolute position from a boundary --
  
  v2s.boundary_class         = cfg_menu;
  v2s.boundary_class.tag     = 'class';
  v2s.boundary_class.name    = 'Surface';
  v2s.boundary_class.labels  = {'central','inner','outer'}; % hull?
  v2s.boundary_class.values  = {1 2 3};
  % "inner" and "outer" are not yet prepared...
  v2s.boundary_class.labels  = {'central'}; % hull?
  v2s.boundary_class.values  = {1};
  v2s.boundary_class.val     = {1};
  v2s.boundary_class.help    = {
    'Surface used for distance description.'
  };

  v2s.boundary_pos         = cfg_entry;
  v2s.boundary_pos.tag     = 'pos';
  v2s.boundary_pos.name    = 'Absolute Position';
  v2s.boundary_pos.strtype = 'r';
  v2s.boundary_pos.val     = {0};
  v2s.boundary_pos.num     = [1 1];
  v2s.boundary_pos.help    = {
    'Absolute position from surface. Use negative values for deeper positions pointing inwards.'
    'All values are limited by the maximum possible distance within a cortical structure such as gyri or sulci.'
  };
  
  
  v2s.boundary         = cfg_exbranch;
  v2s.boundary.tag     = 'boundary';
  v2s.boundary.name    = 'Absolute Position From a Surface Boundary';
  v2s.boundary.val     = {
    v2s.boundary_class ...
    v2s.boundary_pos ...
    };
  v2s.boundary.help    = {
    'Map volumetric data at an absolute position from a surface.'
    'A value of -1 from the central surface will map GM values at a position of 1 mm inwards to the central surface. '
    'A value of 1 from the central surface will map GM values at a position of 1 mm outwards to the central surface. '
  };

  %% -- relative position within a tissue class
   
  v2s.tissue_class         = cfg_menu;
  v2s.tissue_class.tag     = 'class';
  v2s.tissue_class.name    = 'Tissue Class';
  v2s.tissue_class.labels  = {'GM','WM'};
  v2s.tissue_class.values  = {1 2};
  v2s.tissue_class.val     = {1};
  v2s.tissue_class.help    = {
    'Tissue class for which the relative positions are estimated.'
  };

  v2s.tissue_pos         = cfg_entry;
  v2s.tissue_pos.tag     = 'pos';
  v2s.tissue_pos.name    = 'Relative Position';
  v2s.tissue_pos.strtype = 'r';
  v2s.tissue_pos.val     = {0};
  v2s.tissue_pos.num     = [1 1];
  v2s.tissue_pos.help    = {
    'Relative position within the tissue class, where 1 describe the deepest position and 0 the outer position.'
    'For GM a value of 0 describes the GM/CSF interface and a value of 1 describes the GM/WM interface.'
    'For WM a value of 0 describes the WM/GM interface and a value of 1 describes the WM centerline or skeleton as a thinned version of the WM.'
  };
  
  v2s.tissue         = cfg_branch;
  v2s.tissue.tag     = 'tissue';
  v2s.tissue.name    = 'Relative Position Within a Tissue Class';
  v2s.tissue.val     = {
    v2s.tissue_class ...
    v2s.tissue_pos ...
    };
  v2s.tissue.help    = {
    'Map volumetric data from a relative position within a tissue class.'
  };

  %% -- tissue sample --
  
  v2s.tissuerange_stepsize         = cfg_entry;
  v2s.tissuerange_stepsize.tag     = 'stepsize';
  v2s.tissuerange_stepsize.name    = 'Stepsize';
  v2s.tissuerange_stepsize.strtype = 'r';
  v2s.tissuerange_stepsize.val     = {0.1};
  v2s.tissuerange_stepsize.num     = [1 1];
  v2s.tissuerange_stepsize.help    = {
    'Relative stepsize of the sample points centered around 0.5.  For example a stepsize of 0.3 will result in 3 sample points at a relative position of 0.2, 0.5 and 0.8, '
    'while a stepsize of 0.1 will result in sample points 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0.' 
  };

  v2s.tissuerange_sample         = cfg_menu;
  v2s.tissuerange_sample.tag     = 'sample';
  v2s.tissuerange_sample.name    = 'Sample Function';
  v2s.tissuerange_sample.labels  = {'mean','max','min',}; %,'median','std'
  v2s.tissuerange_sample.values  = {{'average'},{'max'},{'min'}}; %,{'median'},{'std'}
  v2s.tissuerange_sample.val     = {{'average'}};
  v2s.tissuerange_sample.help    = {
    'Tissue boundary used for distance description.'
  };
  
  v2s.tissuerange         = cfg_branch;
  v2s.tissuerange.tag     = 'tissuerange';
  v2s.tissuerange.name    = 'Sample Across Several Relative Positions Within a Tissue Class';
  v2s.tissuerange.val     = {
    v2s.tissue_class ...
    v2s.tissuerange_stepsize ...
    v2s.tissuerange_sample ...
    };
  v2s.boundary.help    = {
    'Extract a set of values within a tissue class with a specified relative sample distance and average these values by mean, median, min, max or standard deviation'
  };

  %% -- boundary sample? 
  
  v2s.boundaryrange_class         = cfg_menu;
  v2s.boundaryrange_class.tag     = 'class';
  v2s.boundaryrange_class.name    = 'Tissue Class';
  v2s.boundaryrange_class.labels  = {'GM','WM'}; % hull
  v2s.boundaryrange_class.values  = {1 2};
  v2s.boundaryrange_class.val     = {1};
  v2s.boundaryrange_class.help    = {
    'Tissue boundary used for distance description.'
  };

  v2s.boundaryrange_stepsize         = cfg_entry;
  v2s.boundaryrange_stepsize.tag     = 'stepsize';
  v2s.boundaryrange_stepsize.name    = 'Stepsize';
  v2s.boundaryrange_stepsize.strtype = 'r';
  v2s.boundaryrange_stepsize.val     = {0.1};
  v2s.boundaryrange_stepsize.num     = [1 1];
  v2s.boundaryrange_stepsize.help    = {
    'Absolute stepsize of the sample points centered around the tissue boundary. '
    'The function use 5 sample points.  The maximum value is 1 (to avoid problems with cortical structures).'
    'This means a stepsize of 0.5 will generate 5 sample points with absolute distance to the choosen boundary of -1, -0.5, 0.0, 0.5, 1.0 mm.'
  };
 
  v2s.boundaryrange_sample         = cfg_menu;
  v2s.boundaryrange_sample.tag     = 'sample';
  v2s.boundaryrange_sample.name    = 'Sample Function';
  v2s.boundaryrange_sample.labels  = {'mean','max','min',}; %,'median','std'
  v2s.boundaryrange_sample.values  = {{'average'},{'max'},{'min'}}; %,{'median'},{'std'}
  v2s.boundaryrange_sample.val     = {{'average'}};
  v2s.boundaryrange_sample.help    = {
    'Tissue boundary used for distance description.'
  };
  
  v2s.boundaryrange         = cfg_branch;
  v2s.boundaryrange.tag     = 'boundaryrange';
  v2s.boundaryrange.name    = 'Absolute Position From a Tissue Boundary (Sample)';
  v2s.boundaryrange.val     = {
    v2s.boundaryrange_class ...
    v2s.boundaryrange_stepsize ...
    v2s.boundaryrange_sample ...
    };
  v2s.boundary.help    = {
    'Extract a set of values around a tissue boundary with a specified absolute sample distance and average these values by mean, median, minimum, maximum or standard deviation'
  };


  %% -- Mapping function

  v2s.mapping         = cfg_choice;
  v2s.mapping.tag     = 'mapping';
  v2s.mapping.name    = 'Mapping Function';
  v2s.mapping.values  = {
...      v2s.tissuerange ...
...      v2s.tissue ...
...      v2s.boundaryrange ...
    v2s.boundary ...
  }; 
  v2s.mapping.help    = {
    'Volume extration type. '
    '  tissue-range:'
    '    extract a set of values within a tissue class with a specified relative sample '
    '    distance and average these values by mean, median, minimum, maximum or standard deviation'
    '  tissue-based:' 
    '    extract one value with a specified relative position within a tissue class'
    '  boundary-range:' 
    '    extract a set of values within a tissue class with a specified relative sample'
    '    distance and average these values by mean, median, minimum, maximum or standard deviation'
    '  boundary-based: '
    '    extract one value from a specified absolute distance from a tissue interface'
    '' 
  };
  v2s.mapping.val     = {v2s.boundary};
  


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
  v2s.data_sub.name    = 'Volumes in native space';
  v2s.data_sub.filter  = 'image';
  v2s.data_sub.ufilter = '^(?!wm|wp|w0rp|wc).*'; % no normalized images
  v2s.data_sub.num     = [1 Inf];
  v2s.data_sub.help    = {
    'Select volumes in native (subject) space.'
  };

  v2s.vol2surf      = cfg_exbranch;
  v2s.vol2surf.tag  = 'vol2surf';
  v2s.vol2surf.name = 'Map Volume (Native Space) to Individual Surface';
  v2s.vol2surf.val = {
    v2s.data_sub ...
    v2s.data_surf_sub_lh ...
    v2s.datafieldname ...
    v2s.interp ...
    v2s.mapping ...
    };
% not yet working because of coordination issues
%  v2s.vol2surf.prog = @vbm_surf_vol2surf;
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
  v2s.data_surf_avg_lh.val{1}  = {fullfile(spm('dir'),'toolbox','vbm12','templates_surfaces','lh.central.Template_T1_IXI555_MNI152.gii')};
  v2s.data_surf_avg_lh.dir     = fullfile(spm('dir'),'toolbox','vbm12');
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
    'Select spatially normalized volumes (in template space).'
  };

  v2s.vol2tempsurf      = cfg_exbranch;
  v2s.vol2tempsurf.tag  = 'vol2surftemp';
  v2s.vol2tempsurf.name = 'Map Normalized Volume to Template Surface';
  v2s.vol2tempsurf.val  = {
    v2s.data_norm ...
    v2s.data_surf_avg_lh ...
    v2s.datafieldname ...
    v2s.interp ...
    v2s.mapping ...
  };
  v2s.vol2tempsurf.prog = @vbm_surf_vol2surf;
  v2s.vol2tempsurf.help = {
    'Map spatially normalized data (in template space) to template surface.'
    'The template surface was generated by VBM12 surface processing [1] of the average of 555 Dartel-normalized images of the IXI database that were also used to create the IXI Dartel template.   '
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


% surface calculations 
% ----------------------------------------------------------------------
% estimation per subject (individual and group sampling):
% g groups with i datafiles and i result datafile
 
  sc.cdata         = cfg_files;
  sc.cdata.tag     = 'cdata';
  sc.cdata.name    = 'Surface Data Files';
  sc.cdata.filter  = 'any';
  sc.cdata.ufilter = '^s.mm.*';
  sc.cdata.num     = [1 Inf];
  sc.cdata.help    = {'These are the surface data files that are used by the calculator.  They are referred to as s1, s2, s3, etc in the order that they are specified.'};
  
  sc.cdata_sub         = cfg_files;
  sc.cdata_sub.tag     = 'cdata';
  sc.cdata_sub.name    = 'Surface Data Files';
  sc.cdata_sub.filter  = 'any';
  sc.cdata_sub.ufilter = '[rl]h.(?!cent|sphe|defe).*';
  sc.cdata_sub.num     = [1 Inf];
  sc.cdata_sub.help    = {'These are the surface data files that are used by the calculator.  They are referred to as s1, s2, s3, etc in the order that they are specified.'};
   
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
    'Name of the surface data part of the filename.'
    ''
    '  [s*mm.][rh|lh].DATANAME[.resampled|].subjectname[.gii]' 
  };

  sc.expression         = cfg_entry;
  sc.expression.tag     = 'expression';
  sc.expression.name    = 'Expression';
  sc.expression.strtype = 's';
  sc.expression.num     = [1 Inf];
  sc.expression.val     = {' '};
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
  surfcalc.prog = @vbm_surf_calc;
  surfcalc.help = {
    'Mathematical operations for surface data (textures).'
    'It works similar to ''spm_imcalc''.  The input surface data must have the same number of entries.  This means that the must came from same hemisphere of a subject, or the have to be resampled.'
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
  surfcalcsub.prog = @vbm_surf_calc;
  surfcalcsub.help = {
    'Mathematical operations for surface data sets (textures).'
    'In contrast to the ''Surface Calculator'' it allows the definition of datasets sets ''si'' for multiple subjects.  Therefore, each sample requires surface data of the same subjects to evalute the expression for each subject.'
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
    data_surf.ufilter = '[lr]h.(?!cent|sphe|defe|gyrus|sulcuswidth).*';
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
    'Select filter size for smoothing. For cortical thickness a good starting value is 15mm, while other surface parameters based on cortex folding (e.g. gyrification, cortical complexity) need a larger filter size of about 25mm.'};
  
  surfresamp      = cfg_exbranch;
  surfresamp.tag  = 'surfresamp';
  surfresamp.name = 'Resample and Smooth Surface Data';
  surfresamp.val  = {data_surf,fwhm};
  surfresamp.prog = @vbm_surf_resamp;
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
  surfresamp_fs.prog = @vbm_surf_resamp_freesurfer;
  surfresamp_fs.help = {
  'If you have existing freesurfer thickness data this function can be used to resample these data, smooth the resampled data, and convert freesurfer data to gifti format.'};





%% Toolset
%-----------------------------------------------------------------------
  
  stools = cfg_choice;
  stools.name   = 'Surface Tools';
  stools.tag    = 'stools';
  stools.values = { ...
    check_mesh_cov, ...
    surfextract, ...
    surfresamp, ...
    surfresamp_fs,...
    v2s.vol2surf, ...
    v2s.vol2tempsurf, ...
    surfcalc, ...
    surfcalcsub ...
    };

return


%% Result files
%_______________________________________________________________________

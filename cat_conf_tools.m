function tools = cat_conf_tools(expert)
% wrapper for calling CAT utilities
% 
% tools = cat_conf_tools(expert)
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  if ~exist('expert','var'), expert = 1; end
  
% multi use fields
% -----------------------------------------------------------------------

  % just used once 
  data_xml          = cfg_files;
  data_xml.name     = 'Quality measures';
  data_xml.tag      = 'data_xml';
  data_xml.filter   = 'xml';
  data_xml.ufilter  = '^cat_.*\.xml$';
  data_xml.val      = {{''}};
  data_xml.num      = [0 Inf];
  data_xml.help     = {
    'Select the quality measures that are saved during segmentation as xml-files in the report folder.'
    'Please note, that the order of the xml-files should be the same as the other data files.'
    };

  outdir            = cfg_files;
  outdir.tag        = 'outdir';
  outdir.name       = 'Output directory';
  outdir.filter     = 'dir';
  outdir.ufilter    = '.*';
  outdir.num        = [0 1];
  outdir.help       = {'Select a directory where files are written.'};
  outdir.val{1}     = {''};

  subdir            = cfg_entry;
  subdir.tag        = 'subdir';
  subdir.name       = 'Additional output directory';
  subdir.strtype    = 's';
  subdir.num        = [0 Inf];
  subdir.val        = {'EVA_volmod'};
  subdir.help       = {
    'The directory is created within the choosen output directory. If no name is given no subdirecty is created. ' ''};

  data              = cfg_files; 
  data.tag          = 'data';
  data.name         = 'Volumes';
  data.filter       = 'image';
  data.ufilter      = '.*';
  data.num          = [1 Inf];
  data.help         = {''};

  
  % Do not process, if result already exists and is younger than the
  % original image, i.e., if the original was changed then it will be 
  % processed again. The function is quite helpful to develop and test
  % SPM batches and avoid reprocessing of slow steps.
  lazy         = cfg_menu;
  lazy.tag     = 'lazy';
  lazy.name    = 'Lazy processing';
  lazy.labels  = {'Yes','No'};
  lazy.values  = {1,0};
  lazy.val     = {0};
  lazy.help    = {
    'Do not process data if the result already exists. '
  };


% CHECK USE
  data_vol          = cfg_files;
  data_vol.name     = 'Sample data';
  data_vol.tag      = 'data_vol';
  data_vol.filter   = 'image';
  data_vol.num      = [1 Inf];
  data_vol.help     = {'These are the (spatially registered) data. They must all have the same image dimensions, orientation, voxel size etc. Furthermore, it is recommended to use unsmoothed files.'};
  
  % also this is a separate function that is used for the results
  spm_type          = cfg_menu; %
  spm_type.tag      = 'spm_type';
  spm_type.name     = 'Data type of the output images';
  spm_type.labels   = {'same','uint8','int8','uint16','int16','single'};
  spm_type.values   = {0 2 256 512 4 16};
  spm_type.val      = {16};
  spm_type.help     = {
    'SPM data type of the output image. Single precision is recommended, but  uint16 also provides good results. Internal scaling supports a relative high accuracy for the limited number of bits, special values such as NAN and INF (e.g. in the background) will be lost and NAN is converted to 0, -INF to the minimum, and INF to the maximum value. '
    ''
  };

  % also this limit is a separate function that is used for the noise filter
  % and therefore included here
  intlim            = cfg_entry;
  intlim.tag        = 'intlim';
  intlim.name       = 'Global intensity limitation';
  intlim.strtype    = 'r';
  intlim.num        = [1 inf];
  intlim.val        = {100};
  intlim.help       = {
    'General intensity limitation to remove strong outliers by using 100%% of the original histogram values. '
    'You can also specify the lower and upper boundary seperatlym, e.g. [80 99], what will keep only 80%% of the (many) low (background values) but 99%% of the (rew) high intensity (skull) values. ' 
    ''
  };

  prefix            = cfg_entry;
  prefix.tag        = 'prefix';
  prefix.name       = 'Filename prefix';
  prefix.strtype    = 's';
  prefix.num        = [0 Inf];
  prefix.val        = {''};
  prefix.help       = {''};

  suffix            = cfg_entry;
  suffix.tag        = 'suffix';
  suffix.name       = 'Filename suffix';
  suffix.strtype    = 's';
  suffix.num        = [0 Inf];
  suffix.val        = {''};
  suffix.help       = {''};

  fname             = prefix; 
  fname.name        = 'Filename';
  fname.tag         = 'fname';
  fname.val         = {'CATcheckdesign_'}; 
  fname.help        = {'Basic filename to save figures.'};
   
  save              = cfg_menu;
  save.name         = 'Save & close windows';
  save.tag          = 'save';
  save.labels       = {'Save & close','Save only','No'};
  save.values       = {2,1,0};
  save.val          = {0};
  save.help         = {'Save and close figures for batch processing.'};
  
  verb                                = cfg_menu;
  verb.tag                            = 'verb';
  verb.name                           = 'Verbose output';
  verb.labels                         = {'No' 'Yes' 'Yes (Details)'};
  verb.values                         = {0 1 2};
  verb.val                            = {1};
  verb.hidden                         = expert<1;
  verb.help                           = {
    'Be more or less verbose. '
    ''
    };
% get subbatches
% -------------------------------------------------------------------------
  [T2x,T2x_surf,F2x,F2x_surf] = conf_T2x;
  [check_cov, check_cov2]     = conf_check_cov(data_xml,outdir,fname,save,expert);
  quality_measures            = conf_quality_measures;
  [defs,defs2]                = conf_vol_defs;
  nonlin_coreg                = cat_conf_nonlin_coreg;
  createTPM                   = conf_createTPM(data_vol,expert,suffix,outdir); 
  createTPMlong               = conf_createTPMlong(data_vol,expert);
  long_report                 = conf_long_report(data_vol,data_xml,expert);
  headtrimming                = conf_vol_headtrimming(intlim,spm_type,prefix,suffix,verb,lazy,expert);
  check_SPM                   = conf_stat_check_SPM(outdir,fname,save,expert); 
  showslice                   = conf_stat_showslice_all(data_vol);
  maskimg                     = conf_vol_maskimage(data,prefix);
  calcvol                     = conf_stat_TIV;
  spmtype                     = conf_io_volctype(data,intlim,spm_type,prefix,suffix,verb,expert,lazy);
  calcroi                     = conf_roi_fun(outdir);
  [ROI,sROI,ROIsum]           = cat_conf_ROI(expert);
  resize                      = conf_vol_resize(data,prefix,expert,outdir);
  avg_img                     = conf_vol_average(data,outdir);
  realign                     = conf_vol_series_align(data);
  shootlong                   = conf_shoot(expert); 
  [sanlm,sanlm2]              = conf_vol_sanlm(data,intlim,spm_type,prefix,suffix,lazy,expert);
  biascorrlong                = conf_longBiasCorr(data,expert,prefix);
  data2mat                    = conf_io_data2mat(data,outdir);
  boxplot                     = conf_io_boxplot(outdir,subdir,prefix,expert);
  getCSVXML                   = cat_cfg_getCSVXML(outdir,expert);
  file_move                   = conf_io_file_move; 
  %urqio                       = conf_vol_urqio; % this cause problems
  iqr                         = conf_stat_IQR(data_xml);
  qa                          = conf_vol_qa(expert,outdir);
  
  
% create main batch 
% -------------------------------------------------------------------------
  tools = cfg_choice;
  tools.name   = 'Tools';
  tools.tag    = 'tools';
  tools.values = { ...
    showslice, ...                        cat.stat.pre 
    ... qa, ...                           cat.stat.pre
    qa, ...                               
    check_cov, ...                        cat.stat.pre
    check_cov2, ...                       cat.stat.pre
    quality_measures, ...                     cat.stat.pre
    check_SPM, ...                        cat.stat.pre
    ...
    calcvol, ...                          cat.stat.pre
    calcroi, ...                          cat.stat.pre
      ROIsum, ...
    iqr, ...                              cat.stat.pre
    ...
    T2x, F2x, T2x_surf, F2x_surf, ...     cat.stat.models?
    ...
    ... SPLIT THIS FILE ?!
    ...
    sanlm, ...                            cat.pre.vtools.
    sanlm2, ...
    maskimg, ...                          cat.pre.vtools.
    spmtype, ...                          cat.pre.vtools.
    headtrimming, ...                     cat.pre.vtools.
    resize, ...
    ...
    realign, ...                          cat.pre.long.?          % hidden
    shootlong,...                         cat.pre.long.?          % hidden
    biascorrlong,...                      cat.pre.long.?          % hidden
    createTPMlong, ...                    cat.pre.long.createTPM  % hidden
    long_report, ...                  cat.pre.long.report     % hidden
    ...
    createTPM, ...                        
    nonlin_coreg, ...                     cat.pre.vtools.
    defs, ...                             cat.pre.vtools.
    defs2, ...                            cat.pre.vtools.
    avg_img, ...                          cat.pre.vtoolsexp.
    data2mat, ...                         cat.pre.vtools.
    ...
    boxplot, ...                          cat.stat.eval ... print of XML data by the boxplot function and saving as images  
    getCSVXML, ...                        cat.stat.eval ... read out of XML/CSV data and export as batch dependency 
    file_move, ...
    ...                                   
    };
  
  %RD202005: the cause problems at Christians installation ... check it 
  %if expert 
  %  tools.values = [tools.values,{urqio}]; 
  %end
return
%_______________________________________________________________________
function file_move = conf_io_file_move

% ---------------------------------------------------------------------
% file_move Move/Delete Files
% ---------------------------------------------------------------------


% ---------------------------------------------------------------------
% files Files to move/copy/delete
% ---------------------------------------------------------------------
files         = cfg_files;
files.tag     = 'files';
files.name    = 'Files to move/copy/delete';
files.help    = {'These files will be moved, copied or deleted.'};
files.filter = {'any'};
files.ufilter = '.*';
files.num     = [0 Inf];
% ---------------------------------------------------------------------
% moveto Move to
% ---------------------------------------------------------------------
moveto         = cfg_files;
moveto.tag     = 'moveto';
moveto.name    = 'Move to';
moveto.help    = {'Files will be moved to the specified directory.'};
moveto.filter = {'dir'};
moveto.ufilter = '.*';
moveto.num     = [1 1];
% ---------------------------------------------------------------------
% copyto Copy to
% ---------------------------------------------------------------------
copyto         = cfg_files;
copyto.tag     = 'copyto';
copyto.name    = 'Copy to';
copyto.help    = {'Files will be moved to the specified directory.'};
copyto.filter = {'dir'};
copyto.ufilter = '.*';
copyto.num     = [1 1];
% ---------------------------------------------------------------------
% moveto Move to
% ---------------------------------------------------------------------
moveto1         = cfg_files;
moveto1.tag     = 'moveto';
moveto1.name    = 'Move to';
moveto1.help    = {'Files will be moved to the specified directory.'};
moveto1.filter = {'dir'};
moveto1.ufilter = '.*';
moveto1.num     = [1 1];
% ---------------------------------------------------------------------
% pattern Pattern
% ---------------------------------------------------------------------
pattern         = cfg_entry;
pattern.tag     = 'pattern';
pattern.name    = 'Pattern';
pattern.help    = {'The regular expression pattern to look for.'};
pattern.strtype = 's';
pattern.num     = [1  Inf];
% ---------------------------------------------------------------------
% repl Replacement
% ---------------------------------------------------------------------
repl         = cfg_entry;
repl.tag     = 'repl';
repl.name    = 'Replacement';
repl.help    = {'This string (or pattern) will be inserted instead.'};
repl.strtype = 's';
repl.num     = [1  Inf];
% ---------------------------------------------------------------------
% patrep Pattern/Replacement Pair
% ---------------------------------------------------------------------
patrep         = cfg_branch;
patrep.tag     = 'patrep';
patrep.name    = 'Pattern/Replacement Pair';
patrep.val     = {pattern repl };
% ---------------------------------------------------------------------
% patreplist Pattern/Replacement List
% ---------------------------------------------------------------------
patreplist         = cfg_repeat;
patreplist.tag     = 'patreplist';
patreplist.name    = 'Pattern/Replacement List';
patreplist.help    = {'Regexprep supports a list of multiple patterns and corresponding replacements. These will be applied to the filename portion (without path, without extension) one after another. E.g., if your filename is ''testimage(.nii)'', and you replace ''test'' with ''xyz'' and ''xyzim'' with ''newtestim'', the final filename will be ''newtestimage.nii''.'};
patreplist.values  = {patrep };
patreplist.num     = [1 Inf];
% ---------------------------------------------------------------------
% unique Unique Filenames
% ---------------------------------------------------------------------
unique         = cfg_menu;
unique.tag     = 'unique';
unique.name    = 'Unique Filenames';
unique.help    = {
                  'If the regexprep operation results in identical output filenames for two or more input files, these can not be written/renamed to their new location without loosing data. If you are sure that your regexprep patterns produce unique filenames, you do not need to care about this.'
                  'If you choose to append a running number, it will be zero-padded to make sure alphabetical sort of filenames returns them in the same order as the input files are.'
                  }';
unique.labels = {
                 'Don''t Care'
                 'Append Index Number'
                 }';
unique.values = {
                 false
                 true
                 }';
% ---------------------------------------------------------------------
% moveren Move and Rename
% ---------------------------------------------------------------------
moveren         = cfg_branch;
moveren.tag     = 'moveren';
moveren.name    = 'Move and Rename';
moveren.val     = {moveto1 patreplist unique };
moveren.help    = {'The input files will be moved to the specified target folder. In addition, their filenames (not paths, not extensions) will be changed by replacing regular expression patterns using MATLABs regexprep function. Please consult MATLAB help and HTML documentation for how to specify regular expressions.'};

ren         = cfg_branch;
ren.tag     = 'ren';
ren.name    = 'Rename';
ren.val     = {patreplist unique };
ren.help    = {'The input files will be moved to the specified target folder. In addition, their filenames (not paths, not extensions) will be changed by replacing regular expression patterns using MATLABs regexprep function. Please consult MATLAB help and HTML documentation for how to specify regular expressions.'};


% ---------------------------------------------------------------------
% copyto Copy to
% ---------------------------------------------------------------------
copyto1         = cfg_files;
copyto1.tag     = 'copyto';
copyto1.name    = 'Copy to';
copyto1.help    = {'Files will be moved to the specified directory.'};
copyto1.filter = {'dir'};
copyto1.ufilter = '.*';
copyto1.num     = [1 1];
% ---------------------------------------------------------------------
% pattern Pattern
% ---------------------------------------------------------------------
pattern         = cfg_entry;
pattern.tag     = 'pattern';
pattern.name    = 'Pattern';
pattern.help    = {'The regular expression pattern to look for.'};
pattern.strtype = 's';
pattern.num     = [1  Inf];
% ---------------------------------------------------------------------
% repl Replacement
% ---------------------------------------------------------------------
repl         = cfg_entry;
repl.tag     = 'repl';
repl.name    = 'Replacement';
repl.help    = {'This string (or pattern) will be inserted instead.'};
repl.strtype = 's';
repl.num     = [1  Inf];
% ---------------------------------------------------------------------
% patrep Pattern/Replacement Pair
% ---------------------------------------------------------------------
patrep         = cfg_branch;
patrep.tag     = 'patrep';
patrep.name    = 'Pattern/Replacement Pair';
patrep.val     = {pattern repl };
% ---------------------------------------------------------------------
% patreplist Pattern/Replacement List
% ---------------------------------------------------------------------
patreplist         = cfg_repeat;
patreplist.tag     = 'patreplist';
patreplist.name    = 'Pattern/Replacement List';
patreplist.help    = {'Regexprep supports a list of multiple patterns and corresponding replacements. These will be applied to the filename portion (without path, without extension) one after another. E.g., if your filename is ''testimage(.nii)'', and you replace ''test'' with ''xyz'' and ''xyzim'' with ''newtestim'', the final filename will be ''newtestimage.nii''.'};
patreplist.values  = {patrep };
patreplist.num     = [1 Inf];
% ---------------------------------------------------------------------
% unique Unique Filenames
% ---------------------------------------------------------------------
unique         = cfg_menu;
unique.tag     = 'unique';
unique.name    = 'Unique Filenames';
unique.help    = {
                  'If the regexprep operation results in identical output filenames for two or more input files, these can not be written/renamed to their new location without loosing data. If you are sure that your regexprep patterns produce unique filenames, you do not need to care about this.'
                  'If you choose to append a running number, it will be zero-padded to make sure alphabetical sort of filenames returns them in the same order as the input files are.'
                  }';
unique.labels = {
                 'Don''t Care'
                 'Append Index Number'
                 }';
unique.values = {
                 false
                 true
                 }';
unique.val    = {false};                
% ---------------------------------------------------------------------
% copyren Copy and Rename
% ---------------------------------------------------------------------
copyren         = cfg_branch;
copyren.tag     = 'copyren';
copyren.name    = 'Copy and Rename';
copyren.val     = {copyto1 patreplist unique };
copyren.help    = {'The input files will be copied to the specified target folder. In addition, their filenames (not paths, not extensions) will be changed by replacing regular expression patterns using MATLABs regexprep function. Please consult MATLAB help and HTML documentation for how to specify regular expressions.'};
% ---------------------------------------------------------------------
% delete Delete
% ---------------------------------------------------------------------
delete         = cfg_const;
delete.tag     = 'delete';
delete.name    = 'Delete';
delete.val     = {false};
delete.help    = {'The selected files will be deleted.'};
% ---------------------------------------------------------------------
% action Action
% ---------------------------------------------------------------------
action         = cfg_choice;
action.tag     = 'action';
action.name    = 'Action';
action.values  = {moveto copyto moveren copyren ren delete };


file_move         = cfg_exbranch;
file_move.tag     = 'file_move';
file_move.name    = 'Move/Rename/Delete Files';
file_move.val     = {files action };
file_move.help    = {'Move, rename or delete files.'};
file_move.prog    = @cat_io_file_move;
file_move.vout    = @vout_file_move;

return
%_______________________________________________________________________
function long_report = conf_long_report(data_vol,data_xml,expert)
% -------------------------------------------------------------------------
% Batch to create a final report of the processing of a set of files of one
% (or multiple) subject(s).
% 
% RD202201: start of development for fast visualisation of longitudinal and
%           test-retest data
% -------------------------------------------------------------------------
  data_vol.name         = 'Volume Data Files';
  data_vol.num          = [0 Inf];
  data_vol.val{1}       = {''};
  
  data_surf             = cfg_files;
  data_surf.tag         = 'data_surf';
  data_surf.name        = '(Left) Surface Data Files';
  data_surf.filter      = 'any';
  data_surf.ufilter     = 'lh.(?!cent|pial|white|sphe|defe|hull|pbt).*';
  data_surf.num         = [0 Inf];
  data_surf.help        = {'Surface data files. Both sides will be processed'};
  data_surf.val{1}      = {''};
   
  avg_vol               = data_vol; 
  avg_vol.tag           = 'avg_vol';
  avg_vol.name          = 'Volume Average Data File (In Development)'; % ###### not implemented yet ######
  avg_vol.num           = [0 1];
  avg_vol.help          = {'Segmentation of an average volume T1 to estimate further measures.' ''}; 
  avg_vol.val{1}        = {''};
  avg_vol.hidden        = expert<2;

  avg_surf              = data_surf; 
  avg_surf.tag          = 'avg_vol';
  avg_surf.name         = '(Left) Surface Average Data File (In Development)'; % ###### not implemented yet ######
  avg_surf.num          = [0 1];
  avg_surf.help         = {'Surface/thickness of an average volume T1 to estimate further measures.' ''}; 
  avg_surf.val{1}       = {''}; 
  avg_surf.hidden       = expert<2;
  
  % selected automatically ... need further controlling routines for covariance analysis
  xmls                  = data_xml; 
  xmls.name             = 'XML Data Files (In Development)'; % ###### not implemented yet ######
  xmls.hidden           = expert<2; 
  
  timepoints            = cfg_entry;
  timepoints.tag        = 'timepoints';
  timepoints.name       = 'Timepoints (In Development)'; % ###### not implemented yet ######
  timepoints.help       = {'Define difference between timepoints in years. '}; 
  timepoints.strtype    = 'r';
  timepoints.num        = [0 inf];
  timepoints.val        = {[]}; 
  timepoints.hidden     = expert<2;
  
  
  % == options ==
  smoothvol             = cfg_entry;
  smoothvol.tag         = 'smoothvol';
  smoothvol.name        = 'Volumetric Smoothing';
  smoothvol.help        = {'FWHM of volumetric smoothing in mm.'}; 
  smoothvol.strtype     = 'r';
  smoothvol.num         = [1 1];
  smoothvol.val         = {3}; 
  
  smoothsurf            = cfg_entry;
  smoothsurf.tag        = 'smoothsurf';
  smoothsurf.name       = 'Thickness Smoothing';
  smoothsurf.help       = {'Amount of surface-based smoothing in mm'}; 
  smoothsurf.strtype    = 'r';
  smoothsurf.num        = [1 1];
  smoothsurf.val        = {12}; 
  
  midpoint              = cfg_menu;
  midpoint.tag          = 'midpoint';
  midpoint.name         = 'Scaling (In Development)'; % ###### not implemented yet ######
  midpoint.labels       = {
    'first image'
    'mean'
    };
  midpoint.values       = {0;1};
  midpoint.val          = {0};
  midpoint.help         = {'Data scaling by first image or by mean value. ' ''}; 
  midpoint.hidden       = expert<2;
  
  boxplot               = cfg_menu;
  boxplot.tag           = 'boxplot';
  boxplot.name          = 'Boxplot (In Development)'; % ###### not implemented yet ######
  boxplot.labels        = {
    'no'
    'yes'
    };
  boxplot.values        = {0;1};
  boxplot.val           = {0};
  boxplot.help          = {'Use boxplots.' ''};
  boxplot.hidden        = expert < 2;
  
  plotGMWM              = cfg_menu;
  plotGMWM.tag          = 'plotGMWM';
  plotGMWM.name         = 'Plot WM and GM in one figure'; % ###### not implemented yet ######
  plotGMWM.labels       = {
    'no'
    'yes'
    };
  plotGMWM.values       = {0;1};
  plotGMWM.val          = {1};
  plotGMWM.help         = {'Plot WM and GM in one figure. ' ''};
  plotGMWM.hidden       = expert < 2;
  
  opts                  = cfg_exbranch;
  opts.tag              = 'opts';
  opts.name             = 'Options';
  opts.val              = {smoothvol smoothsurf midpoint plotGMWM}; 
  opts.help             = {'Specify some processing options.' ''};
  opts.hidden           = expert<1;
  
  % == output ==
  vols                  = cfg_menu;
  vols.tag              = 'vols';
  vols.name             = 'Difference Maps';
  vols.labels           = {'No';'Yes'};
  vols.values           = {0,1};
  vols.val              = {0};
  vols.help             = {'Write difference volume maps.' ''};
  
  surfs                 = cfg_menu;
  surfs.tag             = 'surfs';
  surfs.name            = 'Difference Surfaces Data Files';
  surfs.labels          = {'No';'Yes'};
  surfs.values          = {0,1};
  surfs.val             = {0};
  surfs.help            = {'Write difference surface data files.' ''};
  
  xml                   = cfg_menu;
  xml.tag               = 'xml';
  xml.name              = 'XML';
  xml.labels            = {'No';'Yes'};
  xml.values            = {0,1};
  xml.val               = {1};
  xml.help              = {'Write combined XML file.' ''};
  
  output                = cfg_exbranch;
  output.tag            = 'output';
  output.name           = 'Write Output Data';
  output.val            = {vols surfs xml}; 
  output.help           = {'Specify output data.' ''};
  output.hidden         = expert<1;
  
  % == main ==
  long_report           = cfg_exbranch;
  long_report.tag       = 'long_report';
  long_report.name      = 'Longitudinal Report';
  if expert
    long_report.val     = {data_vol avg_vol data_surf avg_surf xmls timepoints opts output};
  else
    long_report.val     = {data_vol data_surf};
  end  
  long_report.prog      = @cat_long_report;
  %long_report.vout      = @vout_long_report; 
  long_report.hidden    = expert<1;
  long_report.help      = {
    };
return
%_______________________________________________________________________
function getCSVXML = cat_cfg_getCSVXML(outdir,expert)
% -------------------------------------------------------------------------
% Batch to read out of XML/CSV data and export as batch dependency/file.
% 
% RD202104
% -------------------------------------------------------------------------

  % n-files, e.g. XML for direct extraction or nii/gii as selector
  files               = cfg_files;
  files.num           = [1 Inf];
  files.tag           = 'files';
  files.name          = 'Subjects';
  files.filter        = 'any';
  files.help          = {'Select XML/NIFTI/GIFTI files of subjects those XML/CSV data should be extracted. '};
  
  % 0..1-file ... maybe n later
  csvfile             = cfg_files;
  csvfile.num         = [0 1];
  csvfile.tag         = 'csvfile';
  csvfile.name        = 'CSV file';
  csvfile.filter      = 'any';
  csvfile.ufilter     = '.*csv';
  csvfile.val         = {{''}};
  csvfile.help        = {
   ['Select one CSV file that contains further information, e.g. age or sex.  The first line has to be the header with the name of the variables.  ' ...
    'The first row has to include an unique identifier for the selected subjects files give above, e.g. the subject ID, the filename, or path if the filename is not unique. ' ...
    'For instance, a file IXI_IOP_493 can be identified by the subject ID 493 given in the IXI CSV table. ' ...
    'However, filenames in BIDS are not suited for identification and you has to specify the "Path/filename selector" to select the directory entry that include the ID. ']
    ''
    };
 
 
  % set of variables names for extraction ... preselection TIV IQR ...
  % the variables were extracted and a depency for each created
  % The CSV selection is a bit more tricky. 
  fields              = cfg_entry;
  fields.tag          = 'fields';
  fields.name         = 'XML and CSV fieldnames';
  fields.strtype      = 's+';
  fields.num          = [0 inf]; 
  fields.val          = {{'ALLCSV'}};
  fields.help         = {
   ['Enter the fieldnames (XML) or columns names (CSV) you want to get here. ' ...
    'The fieldnames where used to create the depency object and will be converted to variables. ' ...
    'The CSV columns should be useable as variable otherwise you have to rename them in the file. ' ...
    'You can use ALLCSV to create variables of all CSV header fields. ' ] 
    ''
   ['E.g. the catxml contain the TIV in the subfield "subjectmeasures.vol_TIV" that will create the depency variable "subjectmeasures_vol_TIV". ' ...
    'To extract one value of a matrix or cell field use the matlab specification, e.g., to extract the GM value from "subjectmeasures.vol_CGW" use "subjectmeasures.vol_CGW(2)"' ...
    'For the CSV file the rows has to be select by using the same name, but will be converted into a variable a field "sex(1=m,2=f)" will result in "sex_1m_2f". ']
    ''
    'Enter one field per row (what creates a cellstr), e.g.:'
    '  AGE'                           
    '  SEX_ID_(1=m,2=f)'             
    '  subjectmeasures.vol_TIV'       
    '  subjectmeasures.vol_abs_CGW(2)'
    '  ALLCSV'
    ''
    };
  
    
  % quality measures (expert)
  QMfield               = cfg_menu;
  QMfield.tag           = 'xmlfieldsQualityMeasures';
  QMfield.name          = 'Image quality';
  QMfield.labels        = {
    'Noise Contrast Ratio (NCR)'
    'Inhomogeny Contrast Ratio (ICR)'
    'Resolution RMSE (resRMS)'
    'Minimum tissue contrast'
    };
  QMfield.values        = {
    'qualitymeasures.NCR'
    'qualitymeasures.ICR'
    'qualitymeasures.res_RMS'
    'qualitymeasures.contrast'
    };
  QMfield.val           = {'qualitymeasures.NCR'};
  QMfield.help          = {'CAT preprocessing image quality measures (not normalized).' ''};
  
  
  % quality ratings 
  QRfield               = cfg_menu;
  QRfield.tag           = 'xmlfieldsQualityRating';
  QRfield.name          = 'Image quality ratings';
  QRfield.labels        = {
    'Image Quality Rating (IQR)'
    'Noise Contrast Ratio (NCR)'
    'Inhomogeny Contrast Ratio (ICR)'
    'Resolution RMSE (resRMS)'
    'Minimum tissue contrast'
    };
  QRfield.values        = {
    'qualityratings.IQR'
    'qualityratings.NCR'
    'qualityratings.ICR'
    'qualityratings.res_RMS'
    'qualityratings.contrast'
    };
  QRfield.val           = {'qualityratings.IQR'};
  QRfield.help          = {'CAT preprocessing image quality ratings (normalized marks).' ''};
  
  
  % surface measures
  SMfield               = cfg_menu;
  SMfield.tag           = 'xmlfieldsSurfMeasure';
  SMfield.name          = 'Surface quality';
  SMfield.labels        = {
    ... 'Surface Euler number'
    'Surface defect area'
    'Surface defect number'
    'Surface intensity RMSE'
    'Surface position RMSE'
    'Surface self-intersections'
    };
  SMfield.values        = {
    ... 'qualitymeasures.SurfaceEulerNumber'
    'qualitymeasures.SurfaceDefectArea'
    'qualitymeasures.SurfaceDefectNumber'
    'qualitymeasures.SurfaceIntensityRMSE'
    'qualitymeasures.SurfacePositionRMSE'
    'qualitymeasures.SurfaceSelfIntersections'
    };
  SMfield.val           = {'qualitymeasures.SurfaceDefectArea'};
  SMfield.help          = {'CAT preprocessing surface quality measures (not normalized). ' ''};
  
  
  % segmenation measures
  USMfield               = cfg_menu;
  USMfield.tag           = 'xmlfieldsSPMmeasures';
  USMfield.name          = 'Unified segmentation validation measures';
  USMfield.labels        = {
    'SPM log-likelyhood'
    'SPM tissue peak 1 (def. GM)'
    'SPM tissue peak 2 (def. WM)'
    'SPM tissue peak 3 (def. CSF1)'
    'SPM tissue peak 4 (def. CSF2)'
    'SPM tissue volume 1 (GM)'
    'SPM tissue volume 2 (WM)'
    'SPM tissue volume 3 (CSF)'
    'SPM tissue volume 4 (HD1)'
    'SPM tissue volume 5 (HD2)'
    'SPM tissue volume 6 (BG)'
    ...'CAT skull-stripping parameter'
    ...'CAT high BG parameter'
    };
  USMfield.values        = {
    'SPMpreprocessing.ll'
    'SPMpreprocessing.mn(1)'
    'SPMpreprocessing.mn(2)'
    'SPMpreprocessing.mn(3)'
    'SPMpreprocessing.mn(4)'
    'ppe.SPMvols0(1)'
    'ppe.SPMvols0(2)'
    'ppe.SPMvols0(3)'
    'ppe.SPMvols0(4)'
    'ppe.SPMvols0(5)'
    'ppe.SPMvols0(6)'
    ...'ppe.skullstrippedpara'
    ...'ppe.highBGpara'
    ...reg.ll
    ...reg.dt, rmsdt
    };
  USMfield.val           = {'SPMpreprocessing.ll'};
  USMfield.help          = {'SPM preprocessing measures for evaluation of the preprocessing. The tissue peaks depend on the defined number of SPM peaks within a class (default=[1 1 2 3 4 2]). The volumes depend on the TPM that are by default GM, MW, CSF, HD1 (hard tissue), HD2 (soft tissue), backgroun (BG). ' ''};
  USMfield.hidden        = expert<2; 
  
  
  % individual measures
  IMfield               = cfg_menu;
  IMfield.tag           = 'xmlfieldsMorphMeasures';
  IMfield.name          = 'Morphometric measures';
  IMfield.labels        = {
    'Total Intracranial Volume (TIV)'
    'Total Surface Area (TSA)'
    'Mean cortical thickness'
    'Cortical thickness standard deviation'
    'Relative CSF volume'
    'Relative GM  volume'
    'Relative WM  volume'
    'Relative WMH volume'
    'Absolute CSF volume'
    'Absolute GM  volume'
    'Absolute WM  volume'
    'Absolute WMH volume'
    };
  IMfield.values        = {
    'subjectmeasures.vol_TIV'
    'subjectmeasures.surf_TSA'
    'subjectmeasures.dist_thickness{1}(1)'
    'subjectmeasures.dist_thickness{1}(2)'
    'subjectmeasures.vol_rel_CGW(1)'
    'subjectmeasures.vol_rel_CGW(2)'
    'subjectmeasures.vol_rel_CGW(3)'
    'subjectmeasures.vol_rel_CGW(4)'
    'subjectmeasures.vol_abs_CGW(1)'
    'subjectmeasures.vol_abs_CGW(2)'
    'subjectmeasures.vol_abs_CGW(3)'
    'subjectmeasures.vol_abs_CGW(4)'
    ...'ppe.reg.rmsdt'
    ...'ppe.reg.rmsdtc'
    };
  IMfield.val           = {'subjectmeasures.vol_TIV'};
  IMfield.help          = {'Global morphometric measures. ' ''};
  
  
  % individual measures
  PDfield               = cfg_menu;
  PDfield.tag           = 'xmlfieldsQualityMeasures';
  PDfield.name          = 'Predefined XML fields';
  PDfield.labels        = {
    'Total Intracranial Volume (TIV)'
    'Total Surface Area (TSA)'
    ...
    'Image Quality Rating (IQR)'
    'Noise Contrast Ratio (NCR)'
    'Inhomogeny Contrast Ratio (ICR)'
    'Resolution RMSE (resRMS)'
    'Minimum tissue contrast'
    ...
    'Surface defect area'
    'Surface defect number'
    };
  PDfield.values        = {
    'subjectmeasures.vol_TIV'
    'subjectmeasures.surf_TSA'
    ...
    'qualitratings.IQR'
    'qualitratings.NCR'
    'qualitratings.ICR'
    'qualitratings.res_RMS'
    'qualitratings.contrast'
    ...
    'qualitymeasures.SurfaceDefectArea'
    'qualitymeasures.SurfaceDefectNumber'
    };
  PDfield.val           = {'subjectmeasures.vol_TIV'};
  PDfield.help          = {'Predefined XML fields. ' ''};
  
  % - groupname     ... %  opt.names       = [];            % array of group names
  setname               = cfg_entry;
  setname.tag           = 'setname';
  setname.name          = 'Name';
  setname.help          = {'Name of the dataset that replaces the number of the set. ' ''}; 
  setname.strtype       = 's';
  setname.num           = [0 Inf];
  setname.val           = {''};
 
  % - title 
  ftitle              = setname; 
  ftitle.name         = 'title';
  ftitle.tag          = 'Plot title';
  ftitle.help         = {'Name of figure' ''};
  ftitle.val          = {''};
  % - yname (measure/scala)
  fname               = setname; 
  fname.name          = 'name';
  fname.tag           = 'Measure name';
  fname.help          = {'Name of the measure ploted at the y-axis. ' ''};
  %
  fspec               = setname; 
  fspec.name          = 'name';
  fspec.tag           = 'Measure name';
  fspec.help          = {'Name of the measure ploted at the y-axis. ' ''};
  %  opt.ylim        = [-inf inf];    % y-axis scaling
  ylim                = cfg_entry;
  ylim.tag            = 'ylim';
  ylim.name           = 'y-axis limits';
  ylim.help           = {'Limitation of x-axis. '}; 
  ylim.strtype        = 'r';
  ylim.num            = [1 2];
  ylim.val            = {[-inf inf]}; 
  %  opt.subsets     = false(1,numel(data)); 
  
  xmlfield0           = cfg_exbranch;
  xmlfield0.tag       = 'xmlfieldsFull';
  xmlfield0.name      = 'Data field (full)';
  xmlfield0.val       = { ftitle  , fname , fspec , ylim}; 
  xmlfield0.help      = {'Specify set properties such as name or color' ''};
  
  
  xmlfield            = cfg_entry;
  xmlfield.tag        = 'xmlfieldsSimple';
  xmlfield.name       = 'Data field';
  xmlfield.help       = { 
   ['Specify field for data extraction that result in one value per file, e.g., ' ...
    'measures.vol_rel_CGW(1) to extract the first (CSF) volume value. '] ''};
  xmlfield.strtype    = 's';
  xmlfield.num        = [1 Inf];
  xmlfield.def        = @(val) 'subjectmeasures.vol_TIV';
  
  xmlfields           = cfg_repeat;
  xmlfields.tag       = 'xmlfields';
  xmlfields.name      = 'XML fields';
  if expert>1
    xmlfields.values  = {xmlfield,xmlfield0,PDfield,USMfield,QMfield,QRfield,SMfield,IMfield};
  else
    xmlfields.values  = {xmlfield,xmlfield0,PDfield};
  end
  xmlfields.val       = {};
  xmlfields.num       = [0 Inf];
  xmlfields.forcestruct = 1;
  xmlfields.help      = {'Specify manually grouped XML files.'};
  
  % ------
  
  csvdelkom           = cfg_menu;
  csvdelkom.tag       = 'seg';
  csvdelkom.name      = 'CSV delimiter and komma';
  csvdelkom.labels    = {',.',';,',';.',' ,',' .'}; % ... space/tab? ' ,',' .' 
  csvdelkom.values    = {',.',';,',';.',' ,',' .'};
  csvdelkom.val       = {',.'}; 
  csvdelkom.help      = {'Delimiter and komma in the CSV file. '};

  write               = cfg_entry;
  write.tag           = 'fname';
  write.name          = 'Write outputs file name';
  write.strtype       = 's';
  write.val           = {''};
  write.num           = [0 inf];
  write.help          = {
   ['Write outputs in multiple TXT files and add the number of subjects and the name of the variable, ' ...
    'e.g. "IXI" with "555" subjects and "subjectmeasures.vol_TIV" will result in "IXI555_subjectmeasures_vol_TIV". ' ...
    'Do not write anything when empty. ']
    ''
    };

  
  csvid               = cfg_entry;
  csvid.tag           = 'csvIDfd';
  csvid.name          = 'Subpath selector';
  csvid.strtype       = 'w';
  csvid.num           = [0 inf]; 
  csvid.val           = {[]}; 
  csvid.help          = {
   ['Because the filename (=0 or []) does not allways defines the subject ID you can select another directory of the file path. ' ...
    'E.g., for the file ".../myProject/GROUP/SUB01/TP01/T1w/001.nii" you have to define the 3rd ancestor (=-3), ' ...
    'whereas ".../myProject/GROUP/SUB01/TP01/T1w/report/catxml_001.xml" would require the 4th ancestor (=-4). ']
    ''
    };
  
  filesel             = cfg_entry;
  filesel.tag         = 'filesel';
  filesel.name        = 'Filepart selector';
  filesel.help        = {'Limitation of x-axis. '}; 
  filesel.strtype     = 'w';
  filesel.num         = [0 inf];
  filesel.val         = {[]}; 
  filesel.help        = {'Specify a part of the filename, e.g. by 1 to select "IXI002" from "IXI002-Guys-0815-T1.nii". No intput uses the full filename. Two inputs can ' };
  
  fileseps            = cfg_entry;
  fileseps.tag        = 'fileseps';
  fileseps.name       = 'Filename seperators';
  fileseps.strtype    = 's';
  fileseps.val        = {'_-.'};
  fileseps.num        = [0 inf];
  fileseps.help       = {
   'Seperators used within the filename. E.g. to select "IXI002" from "IXI002-Guys-0815-T1.nii" by defining also the ID filename selector with "1". '
    };
  
  % path filename varialbes ? 
  % pathsel + filesel + name
  % filenameselector    
  name            = cfg_entry;
  name.tag        = 'name';
  name.name       = 'Field / variable name';
  name.strtype    = 's';
  name.val        = {};
  name.num        = [1 inf];
  name.help       = {
   'Select a unique name for the variable / field.'};
  
  fnamefield          = cfg_exbranch;
  fnamefield.tag      = 'fnamefields';
  fnamefield.name     = 'Filename field';
  fnamefield.val      = {csvid filesel fileseps name};
  fnamefield.help     = {'' ''};
  
  fnamefields         = cfg_repeat;
  fnamefields.tag     = 'fnamefields';
  fnamefields.name    = 'Filename fields';
  fnamefields.values  = {fnamefield}; 
  fnamefields.val     = {};
  fnamefields.num     = [0 Inf];
  %fnamefields.forcestruct = 0;
  fnamefields.help    = {'Selectors to define the subject ID by a given path/filename, e.g., the IXI filename also include a site ID and weighting: "IXI002-Guys-0815-T1.nii'};
  
  % ############ improve help
  idselector          = cfg_exbranch;
  idselector.tag      = 'idselector';
  idselector.name     = 'ID filename selector';
  idselector.val      = {csvid filesel fileseps};
  idselector.help     = {'Selectors to define the subject ID by a given path/filename, e.g., the IXI filename also include a site ID and weighting: "IXI002-Guys-0815-T1.nii'};
    
  verb                                = cfg_menu;
  verb.tag                            = 'verb';
  verb.name                           = 'Verbose output';
  verb.labels                         = {'No' 'Yes' 'Yes (Details)'};
  verb.values                         = {0 1 2};
  verb.val                            = {1};
  %verb.hidden                         = expert<1;
  verb.help                           = {
    'Be more or less verbose. '
    ''
    };
  
  getCSVXML           = cfg_exbranch;
  getCSVXML.tag       = 'getCSVXML';
  getCSVXML.name      = 'XML/CSV readout';
  getCSVXML.val       = {files csvfile csvdelkom xmlfields fields fnamefields idselector outdir write verb};
  getCSVXML.prog      = @cat_stat_getCSVXMLfield;
  getCSVXML.vout      = @vout_stat_getCSVXML;
  getCSVXML.hidden    = expert<1;
  getCSVXML.help      = {
    'This batch allows to extract XML and CSV entries and filename-parts for a given list of a subset of (processed) files to use the in statistical models.  ' 
    '' % XML block
   ['The XML extraction of the CAT*.xml allows for instance to extract informations from the CAT prepcrocessing, such as global volumes, image quality ratings or preprocessing setting, ' ...
    'You need to select the processed files that define the used subjects, e.g. the "catxml_*.xml" or the "lh.thickness.*.gii" files or the original images. ' ...
    'Moreover, specific atlas regions can be extracted from CAT atlas XML files. '] % is this useful ? Use a predefined batch for it, eg. by automaticly analyse the CSV atlas files
    '' % CSV block
   ['Many databases use CSV files to store information such as "IXI_ID; SEX_ID; HEIGHT; WEIGHT; ...; AGE" for the IXI dataset.  ' ...
    'You also need to select the delimiter type of the CSV file and specify how to distinguish the part of the filename that contains the subject ID (this must be the first column in the CSV file) and the fieldnames (e.g. "SEX_ID" or "AGE").']
    '' % FILE selector block ?
    'You can write the results into a text file or you can use the DEPENDENCY function of the SPM batches by choosing the column-wise output vectors. '
    ''
    'Problems can occure if the ID is not fully unique, i.e. if a ID (e.g. 1) is part of another ID (e.g. 101), or if an ID is used multiple times. '
    };
return

%_______________________________________________________________________
function resize = conf_vol_resize(data,prefix,expert,outdir)
% -------------------------------------------------------------------------
% Simple function to resize and scale images. 
% 
% RD202005
% -------------------------------------------------------------------------

           
  % developer with matrix values
  res             = cfg_entry;
  res.tag         = 'res';
  res.name        = 'Resolution';
  res.strtype     = 'r';
  res.num         = [1 inf]; 
  res.val         = {1}; 
  res.help        = {
   ['Voxel resolution in mm.  For isotropocic resolution you can enter 1 value (e.g., 1.2 for 1.2x1.2x1.2 mm' char(179) ...
   '), otherwise you have to specify all 3 dimensions. '] 
  };

  % this is a special expert case
  scale           = cfg_entry;
  scale.tag       = 'scale';
  scale.name      = 'Scaling';
  scale.strtype   = 'r';
  scale.num       = [1 inf]; 
  scale.val       = {1}; 
  scale.hidden    = expert<2; 
  scale.help      = {
   'Scaling to resize the object by changing the voxel size. A value of 0.5 for instance will half the xyz scale of the object.  '
  };

  % trim data (increase or decrease boundary box) 
  trim           = cfg_entry;
  trim.tag       = 'trim';
  trim.name      = 'Trimming';
  trim.strtype   = 'r';
  trim.num       = [1 6]; 
  trim.val       = {[0 0 0 0 0 0]}; 
  trim.hidden    = expert<2; 
  trim.help      = {
   'Change boundary box by adding (positive values) or removing (negative values) voxel on each side (-x,+x,-y,+y,-z,+z). '
  };


  % use header to resample to 
  Pref            = data; 
  Pref.tag        = 'Pref';
  Pref.name       = 'Alternative image space';
  Pref.num        = [1 1]; 
  Pref.val        = {''}; % this is not working
  Pref.help       = {[
    'Alternative output space to resample to another image. ' ...
    'Leave empty to use space of input image. ' ...
    'Set voxel resolution to 0 to use the resolution of this image. ']};
  
  % main setting
  restype         = cfg_choice; 
  restype.tag     = 'restype';
  restype.name    = 'Operation';
  restype.values  = {res,Pref,scale,trim};
  restype.val     = {res}; 
  restype.help    = {'The images can be resize to (i) a specific resolution and (ii) to the space of another images (like in ImCalc). '};
  if ~scale.hidden
    restype.help  = [ restype.help; {'Moreover, the data of the volume can be rescaled, e.g., to adopt template data for other species. The images are not resliced in this case. '}];
  end
    
  % imcalc interpolation field
  imcalc            = spm_cfg_imcalc;
  method            = imcalc.val{6}.val{3}; 
  if expert>1 
  % extended version with additional filtering filtering
  % there are different levels available (FWHM size) but I want to keep it simple
    method.labels{9}  = 'Trilinear (with smooth downsampling)';
    method.values{9}  = -2001; 
    method.labels{10} = '5th Degree Sinc (with light smooth downsampling)';
    method.values{10} = -2005; 
    method.help      = [ method.help'; { '    Trilinear / 5th Degree Sinc (with smooth downsampling)';  '    - If image dimensions are downsampled, prior Gaussian filtering allows denoising and simulation of the partial volume effect.  The FWHM can be defined as the ratio of the new to the original voxel size:   vx_vol_org ./ vx_vol_org - 1.  E.g. an image of 0.2x0.2x0.5 mm downsampled to 0.5x0.5x0.5 mm supports smoothing with FWHM=[3 3 0], which reduces noise along the downsampled axis. '; ''}];  
  end
  clear imcalc
  
  prefix.val      = {'r'};
  prefix.help     = {
    'Use "auto" to add resolution automatically, e.g., "r0.8_*.nii" for final resolution or "rx0.5_*.nii" for the scaling parameter. '
    'If you want the original resolution use 0 in the resolution setting (autoprefix "rorg_") or 1 in the scaling setting. ' 
    };
  resize          = cfg_exbranch;
  resize.tag      = 'resize';
  resize.name     = 'Resize images';
  resize.val      = {data,restype,method,prefix,outdir};
  resize.prog     = @cat_vol_resize;
  resize.vfiles   = @vout_resize;
  resize.vout     = @vout_resize;
  resize.help     = {'Interpolation of images.' ''};
return

%_______________________________________________________________________
function shootlong = conf_shoot(expert)
% -------------------------------------------------------------------------
% This is slightly modified version of the original Shooting that allows to 
% specify another default file. It is required to remove slight movements of
% between scans in longidudinal data.
% Although this is a general batch it is here defined as a private function 
% called only from the CAT longidudinal batch.
% 
% RD202005
% -------------------------------------------------------------------------

  % get shooting toolbox definition
  shoot = tbx_cfg_shoot; 
  
  % find the create template batch
  FN = cell(1,numel(shoot.values)); for fni=1:numel( shoot.values ), FN{fni} = shoot.values{fni}.name; end
  fi = find(cellfun('isempty',strfind(FN,'Run Shooting (create Templates)'))==0,1);
  
  % field to select an m-file with similar to shooting defaults.
  dfile              = cfg_files; 
  dfile.tag          = 'dfile';
  dfile.name         = 'Shooting default file';
  dfile.filter       = 'm';
  dfile.ufilter      = '.*';
  dfile.num          = [0 1];
  dfile.val          = {{''}};
  dfile.help         = {'Select one Shooting default matlab m-file.  If empty Shooting "spm_shoot_defaults" is used. '};
  
  % creat new version 
  shootlong          = shoot.values{fi}; 
  shootlong.prog     = @cat_spm_shoot_template; 
  shootlong.hidden   = expert<2; 
  shootlong.val      = [shootlong.val {dfile}]; 
  
return

%_______________________________________________________________________
function createTPM = conf_createTPM(data,expert,name,outdir)
% -------------------------------------------------------------------------
% Batch to create own templates based on a Shooting template or a CAT pre-
% processing
% 
% RD202005
% -------------------------------------------------------------------------

  % Shooting template input files
  tfiles                 = data; 
  tfiles.tag             = 'tfiles';
  tfiles.name            = 'First spatial registration template';
  tfiles.help            = {'Select all Dartel/Shooting template volumes, i.e., you have to select in general 6 Dartel or 5 Shooting template volumes. ' ''};
  tfiles.ufilter         = '^Template.*';
  tfiles.num             = [1 Inf];
  
  % local intensity normalized T1 input images
  mfiles                 = data; 
  mfiles.tag             = 'mfiles';
  mfiles.name            = 'Normalized tissue maps';
  mfiles.help            = {[ ...
    'Select the averaged normalized tissue volumes.  ' ...
    'Use "CAT Apply Deformation" or "SPM Deformation" tools to apply the mapping of Dartel/Shooting ' ...
    'to map all tissue maps from the individual to the template space.  ' ...
    'Next, use cat_vol_avg to average each tissue class.  ' ...
    'Use CAT expert mode to write all tissue classes p1 to p6 to native space (use the TPM output field for p4-6) and correct WMHs to WM (WMHC=3).  ' ...
    'If this field is empty the template tissues are use to simulate a T1 image.  ' '']};
  mfiles.ufilter         = '.*';
  mfiles.val             = {''}; 
  mfiles.num             = [0 Inf];
  
  % normalized tissue maps
  pfiles                 = data; 
  pfiles.tag             = 'pfiles';
  pfiles.name            = 'Normalized intensity maps';
  pfiles.help            = {[
    'Select all normalized tissue volumes.  ' ...
    'Use "CAT Apply Deformation" or "SPM Deformation" tools to apply the mapping of Dartel/Shooting ' ...
    'to map all local intensity normalized T1 maps from the individual to the template space.  ' ...
    'Next use cat_vol_avg to average all images.  ' ...
    'Use CAT expert mode to write the local intensity normalized maps in native space.  ' ...
    'If this field is empty the template tissues are use instead. ' '']};
  pfiles.ufilter         = '.*';
  pfiles.val             = {''}; 
  pfiles.num             = [0 Inf];
  
  % atlas maps
  afiles                 = data; 
  afiles.tag             = 'afiles';
  afiles.name            = 'Atlas maps';
  afiles.help            = {[
    'Select all atlas maps in native space usually written in the label directory. ' ...
    'Use "CAT Apply Deformation" or "SPM Deformation" tools to apply the mapping of Dartel/Shooting ' ...
    'to map all individual atlas maps from the individual to the template space.  ' ...
    'Next use cat_vol_avg to average all images with discrete interpolation. ' ...
    'Use CAT expert mode to write out selected volume atlases maps into native space.  ' ... 
    'If this field is empty no atlas files are generated and no Template is generated.  ']
    };
  afiles.ufilter         = '.*';
  afiles.val             = {''}; 
  afiles.num             = [0 Inf];
  
  logfile                = data; 
  logfile.tag            = 'logfile';
  logfile.name           = 'Log/report file';
  logfile.help           = {'Select a file where the report is added at the end. '};
  logfile.ufilter         = '.*';
  logfile.val             = {''}; 
  logfile.num             = [0 1];
  
  % input files
  files               = cfg_branch;
  files.tag           = 'files';
  files.name          = 'Files';
  files.val           = {tfiles pfiles mfiles afiles logfile};
  files.help          = {'Define input files. All images has to be in the same space having the same resolution. ' ''}; 
  
  
  
  % options
  fstrength            = cfg_menu;
  fstrength.tag        = 'fstrength';
  fstrength.name       = 'Filter strength';
  fstrength.labels     = {'small' 'medium' 'strong'};
  fstrength.values     = {2 3 4};
  fstrength.val        = {2};
  fstrength.help       = {'Main filter control parameter with 3 levels. ' ''};

  
  % trimming?   > main batch (= boundary box optimization)
  % resolution? > main batch
  
  
  % name 
  name.tag             = 'name';
  name.name            = 'Template name'; 
  name.val             = {'MyTemplate'};
  
  % verb
  verb                 = cfg_menu;
  verb.tag             = 'verb';
  verb.name            = 'Verbose output';
  verb.labels          = {'No' 'Yes'};
  verb.values          = {0 1};
  verb.val             = {1};
  verb.help            = {'Be verbose.' ''};

  % input files
  opt                  = cfg_branch;
  opt.tag              = 'opt';
  opt.name             = 'Options';
  opt.val              = {fstrength,name,outdir,verb,};
  opt.help             = {'Main options.' ''}; 
  
  % verb
  verb                 = cfg_menu;
  verb.tag             = 'verb';
  verb.name            = 'Verbose output';
  verb.labels          = {'No' 'Yes'};
  verb.values          = {0 1};
  verb.val             = {1};
  verb.help            = {'Be verbose.' ''};

  
  %{
  def.write.name        = 'MyTemplate'; % template name
  def.write.outdir      = '';           % main output directory
  def.write.subdir      = '';           % create sub directory
  def.write.TPM         = 1;            % write TPM 
  def.write.TPMc        = 1;            % write seperate TPM classes
  def.write.TPM4        = 1;            % write 4 class TPM  
  def.write.TPM4c       = 1;            % write seperate 4 class TPM
  def.write.T1          = 1;            % write T1  
  def.write.T2          = 1;            % write T2 
  def.write.GS          = 1;            % write Shooting template
  def.write.DT          = 1;            % create and write Dartel template
  job.write.brainmask   = 1;            % write brainmask
  %}
  
  % input files
  write               = cfg_branch;
  write.tag           = 'write';
  write.name          = 'Output';
  write.val           = {};
  write.help          = {'' ''}; 
  
  % main
  createTPM        = cfg_exbranch;
  createTPM.tag    = 'createTPMlong';
  createTPM.name   = 'TPM creation';
  createTPM.val    = {files, opt, write};
  createTPM.prog   = @cat_vol_createTPM;
  createTPM.vfiles = @vout_createTPM;
  createTPM.vout   = @vout_createTPM;
  createTPM.hidden = expert<2;
  createTPM.help   = {
    'Create individual TPMs for preprocessing by using Dartel/Shooting templates. ' 
   ['SPM uses TPMs with 6 tissue classes (GM,WM,CSF,HD1,HD2,BG), whereas the head classes (HD) can be empty.   ' ...
    'However, the affine normalized or a soft non-linear normalized space is expected to obtain best result (see options in cat_main_registration).  ' ...
    'A resolution of 1.5 mm seems to be quite optimal as far as we have to smooth anyway.  ' ...
    'The images will be filtered in different ways to allow soft meanderings of anatomical structures.  ' ...
    'WMHs should probably be corrected to WM (WMHC=2) in the average preprocessing.' ]
      ''};
return

%_______________________________________________________________________
function createTPMlong = conf_createTPMlong(data,expert)
% -------------------------------------------------------------------------
% This is a special version of the cat_vol_createTPM batch only for the
% longitudinal preprocessing without further GUI interaction and well
% defined input.
% 
% RD202005
% -------------------------------------------------------------------------

  % update input
  data.tag             = 'files';
  data.name            = 'GM segments';
  data.help            = {'Select GM segments.  The other tissue classes (2-6) will be selected automaticelly. '};
  
  fstrength            = cfg_menu;
  fstrength.tag        = 'fstrength';
  fstrength.name       = 'Filtermodel';
  fstrength.labels     = {
    'very small (plasticity)'
    'small (plasticty/aging)'
    'medium (aging/development)'
    'strong (development)'};
  fstrength.values     = {1 2 3 4};
  fstrength.val        = {1};
  fstrength.help       = {
    ['Main filter control parameter with 4 settings, (1) very small for variations in plasticity, ' ...
     '(2) small for changes in pasticity/short time aging, (3) medium for changes in long-time aging '...
     'and short-time development, and (3) strong for large variations in long-time development. '] 
     ''
  };

  verb                 = cfg_menu;
  verb.tag             = 'verb';
  verb.name            = 'Verbose output';
  verb.labels          = {'No' 'Yes'};
  verb.values          = {0 1};
  verb.val             = {1};
  verb.help            = {
    'Be verbose.'
    ''
    };

  writeBM                 = cfg_menu;
  writeBM.tag             = 'writeBM';
  writeBM.name            = 'Write brainmask';
  writeBM.labels          = {'No' 'Yes'};
  writeBM.values          = {0 1};
  writeBM.val             = {1};
  writeBM.help            = {
    'Save brainmask image.'
    ''
    };
  
  % main
  createTPMlong        = cfg_exbranch;
  createTPMlong.tag    = 'createTPMlong';
  createTPMlong.name   = 'Longitudinal TPM creation';
  createTPMlong.val    = {data,fstrength,writeBM,verb};
  createTPMlong.prog   = @cat_long_createTPM;
  createTPMlong.vfiles = @vout_createTPMlong;
  createTPMlong.vout   = @vout_createTPMlong;
  createTPMlong.hidden = expert<1; 
  createTPMlong.help   = {
    'Create individual TPMs for longitudinal preprocessing. This is a special version of the cat_vol_createTPM batch only for the longitudinal preprocessing without further GUI interaction and well defined input. '
   ['There has to be 6 tissue classes images (GM,WM,CSF,HD1,HD2,BG) that can be in the native space, the affine or a non-linear normalized space.  ' ...
    'However, the affine normalized or a soft non-linear normalized space is expected to give the best result (see options in cat_main_registration).  ' ...
    'A resolution of 1.5 mm seems to be quite optimal as far as we have to smooth anyway.  ' ...
    'The images will be filtered in different ways to allow soft meanderings of anatomical structures.  ' ...
    'WMHs should probably be corrected to WM (WMHC=2) in the average preprocessing.' ]
      ''};
return

%_______________________________________________________________________
function iqr = conf_stat_IQR(data_xml)
%  ------------------------------------------------------------------------
  iqr_name         = cfg_entry;
  iqr_name.tag     = 'iqr_name';
  iqr_name.name    = 'Output file';
  iqr_name.strtype = 's';
  iqr_name.num     = [1 Inf];
  iqr_name.val     = {'IQR.txt'};
  iqr_name.help    = {'The output file is written to current working directory unless a valid full pathname is given'};

  iqr       = cfg_exbranch;
  iqr.tag   = 'iqr';
  iqr.name  = 'Get Weighted Overall Image Quality';
  iqr.val   = {data_xml,iqr_name};
  iqr.prog  = @cat_stat_IQR;
  iqr.help  = {'This function reads weighted overall image quality from saved xml-files.' ''};
return

%_______________________________________________________________________
function longBiasCorr = conf_longBiasCorr(data,expert,prefix)
% -------------------------------------------------------------------------
% Longitudinal bias correction by using the average segmentation.
% See cat_long_biascorr.
%
% RD202010: First tests showed clear improvements of the timepoints but the
%           whole pipeline seems to be less affected.
%           Hence, corrections are maybe more relevant for plasticity
%           studies or in case of artifacts.
% -------------------------------------------------------------------------

  images        = data;
  images.tag    = 'images';
  images.name   = 'Realigned images of one subject';
  images.num    = [1 inf];
  
  segment       = data; 
  segment.tag   = 'segment';
  segment.num   = [1 1];
  segment.name  = 'Average tissue segmentation of one subject';

  bstr                 = cfg_menu;
  bstr.tag             = 'str';
  bstr.name            = 'Strength of correction';
  bstr.labels          = {'no correction','small','medium','strong','very strong'};
  bstr.values          = {0,0.25,0.5,0.75,1.0};
  bstr.val             = {0.5};
  bstr.help            = {
    'Strength of bias correction.'
    ''
    };

  longBiasCorr        = cfg_exbranch;
  longBiasCorr.tag    = 'longBiasCorr';
  longBiasCorr.name   = 'Longitudinal Bias Correction';
  longBiasCorr.val    = {images,segment,bstr,prefix};
  longBiasCorr.prog   = @cat_long_biascorr;
  longBiasCorr.vout   = @vout_conf_longBiasCorr;
  longBiasCorr.hidden = expert<0; 
  longBiasCorr.help   = {'Bias correction based on the segmentation of the average map.' ''};
return

%_______________________________________________________________________
function qa = conf_vol_qa(expert,outdir) 
% Batch for estimation of image quality by a given input segmentation. 
% There was the idea of a relative common batch that allows to use a wide
% set of maps to allow personal adaptions, e.g., to measure in background
% regions or to use atlas maps for region-specific results. However, this 
% becomes quite complex and would focus on experts that have so specific 
% knowledge that they better write there own code. 
% So I try to keep it simple here to support our image quality measures 
% also for other tissue segmentation, e.g. by SPM, FSL, FreeSurfer.  

  % update input
  data            = cfg_files;
  data.tag        = 'images';
  data.name       = 'Images';
  data.help       = {'Select images that should be evaluated.'};
  data.filter     = 'image';
  data.ufilter    = '.*';
  data.num        = [1 Inf];  
  
  catlab          = data; 
  catlab.ufilter  = '^p0.*';
  catlab.tag      = 'catp0'; 
  catlab.name     = 'Default with CAT label map';
  catlab.help     = {['Select CAT label map with brain tissues (p0*.nii).  Also label maps created by other tissue segmentations can be used, ' ...
    'as long the following labeling is used: CSF=1, GM=2, and WM=3 with intermediate PVE values (e.g., 2.32 for 68% GM and 32% WM.  ']};
  
  catsegp         = data; 
  catsegp.ufilter = '^p1.*';
  catsegp.tag     = 'catp1'; 
  catsegp.name    = 'Default with CAT segment maps';
  catsegp.help    = {'Select corresponing CAT GM tissue segments of the selected images above (p1*.nii).  The WM and CSF maps were selected automatically.  ' ''};
  
  spmsegc         = data; 
  spmsegc.ufilter = '^c1.*';
  spmsegc.tag     = 'spmc1'; 
  spmsegc.name    = 'Default with SPM segment maps';
  spmsegc.help    = {'Select corresponing individual SPM GM tissue segments of the selected images above (c1*.nii).  The WM and CSF maps were selected automatically. ' ''};

  gm              = data; 
  gm.tag          = 'gm'; 
  gm.name         = 'GM segment';
  gm.help         = {'Select images with GM segmentation. ' ''};
  wm              = data; 
  wm.tag          = 'wm'; 
  wm.name         = 'WM segment';
  wm.help         = {'Select images with WM segmentation. ' ''};
  cm              = data; 
  cm.tag          = 'cm'; 
  cm.name         = 'CSF segment';
  cm.help         = {'Select images with CSF segmentation. ' ''};
  seg             = cfg_exbranch; 
  seg.tag         = 'seg';
  seg.name        = 'Brain tissue segmentation';
  seg.help        = {'Select tissue segments of other segmentations' ''}; 
  
% FSL segment maps  

% FS label map

  model           = cfg_choice; 
  model.tag       = 'model';
  model.name      = 'Segmentation';
  model.values    = {catlab,catsegp,spmsegc,seg}; 
  model.val       = {catlab}; 
  model.help      = {[ ...
    'Select a input segmentation for the estimation of the quality measures. ' ...
    'The default model is developed for typcial structural T1/T2/PD-based images with a given brain tissue classification. ']}; 
   
  
  % main options
  % -----------------------------------------------------------------------
  prefix          = cfg_entry;
  prefix.tag      = 'prefix';
  prefix.name     = 'Filename prefix';
  prefix.strtype  = 's';
  prefix.num      = [0 Inf];
  prefix.val      = {'qc_'};
  prefix.help     = {'Specify the string to be prepended to the filenames of the XML file(s). ' ''};

  
% Definition of own XML subfield to extend the CAT-XML file 
%{
  fdname          = cfg_entry; 
  fdname.tag      = 'fdname';
  fdname.name     = 'XML-fieldname';
  fdname.strtype  = 's';
  fdname.num      = [0 Inf];
  fdname.val      = {'qc'};
  fdname.help     = {
    'Specify the field name in the "quality_measure" subfield that will include the quality measurements and ratings. ' ''}; 
  
  update          = cfg_menu;
  update.tag      = 'fdupdate';
  update.name     = 'Update result';
  update.labels   = {'No' 'Yes'};
  update.values   = {0 1};
  update.val      = {1}; 
  update.hidden   = ~expert; 
  update.help     = {['Update (replace) an existing datafield. ' ...
    'Otherwise, the data and time of this estimation process are added to the new fieldname. '];''};
%}
  
  verb            = cfg_menu;
  verb.tag        = 'verb';
  verb.name       = 'Print results';
  verb.labels     = {'0' '1'};
  verb.values     = {0 1};
  verb.val        = {1}; 
  verb.help       = {'Print progress and results. ';''};

  outdir.val{1}   = {'report'}; 
  
  opts            = cfg_branch;
  opts.tag        = 'opts';
  opts.name       = 'Options';
  opts.val        = {outdir, prefix, verb }; % fdname, update,
  opts.help       = {'Basic options. ' ''};

  % main
  qa              = cfg_exbranch;
  qa.tag          = 'iqe';
  qa.name         = 'Image quality estimation';
  qa.val          = {data, model, opts};
  qa.prog         = @cat_vol_qa; 
  qa.vfiles       = @vout_qa; % XML files + values
  qa.hidden       = expert<2;
  qa.help         = {'Image quality estimation based on a set of images and a given set of input segmentation defined by different models. '};
return
  
%_______________________________________________________________________
function [sanlm,sanlm2] = conf_vol_sanlm(data,intlim,spm_type,prefix,suffix,lazy,expert)

  % --- update input variables ---
  data.help         = {'Select images for filtering.'};
  
  prefix.val        = {'sanlm_'};
  prefix.help       = {
    'Specify the string to be prepended to the filenames of the filtered image file(s). Default prefix is "samlm_". Use the keyword "PARA" to add the name of the filter, e.g., "classic" or "optimized-medium".'
    ''
  };
  if expert>1
    prefix.help       = {
      'Specify the string to be prepended to the filenames of the filtered image file(s). Default prefix is "samlm_". Use the keyword "PARA" to add the strength of filtering, e.g. "sanlm_PARA" result in "sanlm_NC#_*.nii".'
      ''
    };
    suffix.val        = {''};
    suffix.help       = {
      'Specify the string to be appended to the filenames of the filtered image file(s). Default suffix is ''''.  Use the keyword "PARA" to add the name of the filter, e.g., "classic" or "optimized-medium".'
      ''
    };
  end
  

  % --- new fields ---
  rician            = cfg_menu;
  rician.tag        = 'rician';
  rician.name       = 'Rician noise';
  rician.labels     = {'Yes' 'No'};
  rician.values     = {1 0};
  rician.val        = {0};
  rician.help       = {
    'MRIs can have Gaussian or Rician distributed noise with uniform or nonuniform variance across the image. If SNR is high enough (>3) noise can be well approximated by Gaussian noise in the foreground. However, for SENSE reconstruction or DTI data a Rician distribution is expected. Please note that the Rician noise estimation is sensitive for large signals in the neighbourhood and can lead to artefacts, e.g. cortex can be affected by very high values in the scalp or in blood vessels.'
    ''
  };

  % remove artifacts
  outlier           = cfg_entry;
  outlier.tag       = 'outlier';
  outlier.name      = 'Strength of outlier correction';
  outlier.strtype   = 'r';
  outlier.num       = [1 1];
  outlier.val       = {1};
  outlier.help      = {
    'Remove strong outliers (salt and pepper noise) with more than n times of the average local correction strength. Larger values will result in stronger corrections, whereas lower values result in less corrections. Changes will be more visible in high quality areas/images.' 
  };

  % developer with matrix values
  NCstr           = cfg_entry;
  NCstr.tag       = 'NCstr';
  NCstr.name      = 'Strength of noise corrections';
  NCstr.strtype   = 'r';
  NCstr.num       = [1 1]; %inf]; % this case did not work with yet
  NCstr.def       = @(val) cat_get_defaults('extopts.NCstr', val{:});
  NCstr.hidden    = expert<1;
  NCstr.help      = {
   ['Strength of the spatial adaptive (sub-resolution) non-local means (SANLM) noise correction. Please note that the filter strength is automatically estimated. Change this parameter only for specific conditions. ' ...
    'Typical values are: none (0), classic (1), light (2), medium (3|-inf), strong (4), heavy (5). The "classic" option use the ordinal SANLM filter without further adaptations. The "light" option uses the half filter strength of "medium" cases. The "strong" option use 8-times of the "medium" filter strength. Sub-resolution filtering is only used in case of high image resolution below 0.8 mm or in case of the "heavy" option. ' ...
    'For the global modified scheme use smaller values (>0) for less denoising, higher values (<=1) for stronger denoising, and "inf" for an automatic estimated threshold. Negative values control the local adaptive scheme, with the default "-inf"|"-1", that successfully tested on a variety of scans. Use higher values (>-1,<0) for less filtering and lower values "<-1" for stronger filtering. The value 0 will turn off any noise correction.']
    ''
  };

  % noise correction level
  NCstrm            = cfg_menu;
  NCstrm.tag        = 'NCstr';
  NCstrm.name       = 'Strength of Noise Corrections';
  NCstrm.def        = @(val) cat_get_defaults('extopts.NCstr', val{:});
  NCstrm.help       = {
    ['Strength of the (sub-resolution) spatial adaptive non local means (SANLM) noise correction. Please note that the filter strength is automatically estimated. Change this parameter only for specific conditions. ' ...
     'The "light" option applies half of the filter strength of the adaptive "medium" cases, whereas the "strong" option uses the full filter strength, force sub-resolution filtering and applies an additional iteration. Sub-resolution filtering is only used in case of high image resolution below 0.8 mm or in case of the "strong" option.']
    ['If you have scans with low amount of noise then use the "light" option. If you have data that was resampled or interpolated in some way (i.e., even within scanning/reconstruction) then the noise is often blurred over multiple voxels and has to be handled on lower resolutions available for the try the "strong" or "heavy" filter setting.' ...
     'If you have multiple scans that should be averaged than you should use the ".. for average" filter settings.']
     'The filter will always leave some low amount of noise in the data that is assumed by preprocessing routines such as the tissue classification with its Gaussian fitting.'
     ''
  };
  if expert
    NCstrm.values   = {2 -inf 4 5 12 14};
    NCstrm.labels   = {'light (adapted half strength; 2)','medium (adapted; default; -1|3|-inf)','strong (low-resolution filtering; 4)','heavy (low-resolution filtering with 2 iterations; 5)','light for averaging (adapted half strength; 12)','strong for averaging (low-resolution filtering; 14)',};
  else
    NCstrm.values   = {2 -inf 4 5};
    NCstrm.labels   = {'light','medium (default)','strong','heavy'};
  end

  addnoise          = cfg_entry;
  addnoise.tag      = 'addnoise';
  addnoise.name     = 'Strength of additional noise in noise-free regions';
  addnoise.strtype  = 'r';
  addnoise.val      = {0.5}; 
  addnoise.num      = [1 1];
  addnoise.help     = {
    'Add minimal amount of noise in regions without any noise to avoid image segmentation problems. This parameter defines the strength of additional noise as percentage of the average signal intensity. '
    ''
  };

  replaceNANandINF         = cfg_menu;
  replaceNANandINF.tag     = 'replaceNANandINF';
  replaceNANandINF.name    = 'Replace NAN and INF';
  replaceNANandINF.labels  = {'Yes' 'No'};
  replaceNANandINF.values  = {1 0};
  replaceNANandINF.val     = {1};
  replaceNANandINF.help    = {
    'Replace NAN by 0, -INF by the minimum and INF by the maximum of the image.'
    ''
    };

  % relative value vs. on/off
  if expert
    relativeFilterStengthLimit          = cfg_entry;
    relativeFilterStengthLimit.tag      = 'relativeFilterStengthLimit';
    relativeFilterStengthLimit.name     = 'Factor of relative filter strength limit';
    relativeFilterStengthLimit.strtype  = 'r';
    relativeFilterStengthLimit.num      = [1 1];
    relativeFilterStengthLimit.val      = {1};
    relativeFilterStengthLimit.hidden   = expert<2;
    relativeFilterStengthLimit.help     = {
      'Limit the relative noise correction to avoid over-filtering of low intensity areas. Low values will lead to less filtering in low intensity areas, whereas high values will be closer to the original filter. INF deactivates the filter. '
      ''
    };
  else
    relativeFilterStengthLimit          = cfg_menu;
    relativeFilterStengthLimit.tag      = 'relativeFilterStengthLimit';
    relativeFilterStengthLimit.name     = 'Use relative filter strength';
    relativeFilterStengthLimit.labels   = {'Yes' 'No'};
    relativeFilterStengthLimit.values   = {1 0};
    relativeFilterStengthLimit.val      = {1};
    relativeFilterStengthLimit.hidden   = expert<2;
    relativeFilterStengthLimit.help     = {
      'Limit the relative noise correction to avoid over-filtering of low intensities areas.'
      ''
      };
  end

  relativeIntensityAdaption             = cfg_entry;
  relativeIntensityAdaption.tag         = 'relativeIntensityAdaption';
  relativeIntensityAdaption.name        = 'Strength of relative intensity adaptation';
  relativeIntensityAdaption.strtype     = 'r';
  relativeIntensityAdaption.num         = [1 1];
  relativeIntensityAdaption.val         = {1};
  relativeIntensityAdaption.hidden      = expert<2;
  relativeIntensityAdaption.help        = {
    'Strength of relative intensity adaptation, with 0 for no adaptation and 1 for full adaptation. The SANLM filter is often very successful in the background and removed nearly all noise. However, routines such as the SPM Unified Segmentation expect Gaussian distribution in all regions and is troubled by regions with too low variance. Hence, a relative limitation of SANLM correction is added here that is based on the bias reduced image intensity. '
    ''
  };

  % very special parameter ...
  % -----------------------------------------------------------------------
  iter         = cfg_entry;
  iter.tag     = 'iter';
  iter.name    = 'Number of additional sub-resolution iterations';
  iter.strtype = 'r';
  iter.num     = [1 1];
  iter.val     = {0};
  iter.hidden  = expert<1;
  iter.help    = {
    'Choose number of additional iterations that can further reduce sub-resolution noise but also anatomical information, e.g. larger blood vessel or small gyri/sulci.'
    ''
  };

  iterm         = cfg_entry;
  iterm.tag     = 'iterm';
  iterm.name    = 'Number of additional iterations';
  iterm.strtype = 'r';
  iterm.num     = [1 1];
  iterm.val     = {0};
  iterm.hidden  = expert<1;
  iterm.help    = {
    'Choose number of additional iterations that can further reduce noise but also anatomical information, e.g. smaller blood-vessels.'
    ''
  };

  relativeIntensityAdaptionTH         = cfg_entry;
  relativeIntensityAdaptionTH.tag     = 'relativeIntensityAdaptionTH';
  relativeIntensityAdaptionTH.name    = 'Strength of smoothing of the relative filter strength limit';
  relativeIntensityAdaptionTH.strtype = 'r';
  relativeIntensityAdaptionTH.num     = [1 1];
  relativeIntensityAdaptionTH.val     = {2};
  relativeIntensityAdaptionTH.hidden  = expert<2;
  relativeIntensityAdaptionTH.help    = {
    'Smoothing of the relative filter strength limitation.'
    ''
  };

  resolutionDependency                = cfg_menu;
  resolutionDependency.tag            = 'resolutionDependency';
  resolutionDependency.name           = 'Resolution depended filtering';
  resolutionDependency.labels         = {'Yes' 'No'};
  resolutionDependency.values         = {1 0};
  resolutionDependency.val            = {0};
  resolutionDependency.hidden         = expert<2;
  resolutionDependency.help           = {
    'Resolution depending filtering with reduced filter strength in data with low spatial resolution defined by the "Range of resolution dependency".'
    ''
    };

  resolutionDependencyRange           = cfg_entry;
  resolutionDependencyRange.tag       = 'resolutionDependencyRange';
  resolutionDependencyRange.name      = 'Range of resolution dependency';
  resolutionDependencyRange.strtype   = 'r';
  resolutionDependencyRange.num       = [1 2];
  resolutionDependencyRange.val       = {[1 2.5]};
  resolutionDependencyRange.hidden    = expert<1;
  resolutionDependencyRange.help      = {
    'Definition of the spatial resolution for "full filtering" (first value) and "no filtering" (second value), with [1 2.5] for typical structural data of humans. '
    ''
  };

  resolutionReduction                 = cfg_menu;
  resolutionReduction.tag             = 'red';
  resolutionReduction.name            = 'Low resolution filtering';
  resolutionReduction.labels          = {'Yes (allways)' 'Yes (only highres <0.8 mm)' 'No'};
  resolutionReduction.values          = {11 1 0};
  resolutionReduction.val             = {0};
  %resolutionReduction.hidden          = expert<1;
  resolutionReduction.help            = {
    'Some MR images were interpolated or use a limited frequency spectrum to support higher spatial resolution with acceptable scan-times (e.g., 0.5x0.5x1.5 mm on a 1.5 Tesla scanner). However, this can result in "low-frequency" noise that can not be handled by the standard NLM-filter. Hence, an additional filtering step is used on a reduces resolution. As far as filtering of low resolution data will also remove anatomical information the filter use by default maximal one reduction with a resolution limit of 1.6 mm. I.e. a 0.5x0.5x1.5 mm image is reduced to 1.0x1.0x1.5 mm, whereas a 0.8x0.8x0.4 mm images is reduced to 0.8x0.8x0.8 mm and a 1x1x1 mm dataset is not reduced at all. '
    ''
    };

  verb                                = cfg_menu;
  verb.tag                            = 'verb';
  verb.name                           = 'Verbose output';
  verb.labels                         = {'No' 'Yes'};
  verb.values                         = {0 1};
  verb.val                            = {1};
  verb.hidden                         = expert<1;
  verb.help                           = {
    'Be verbose.'
    ''
    };
  
  sharpening         = cfg_entry;
  sharpening.tag     = 'sharpening';
  sharpening.name    = 'Sharpening';
  sharpening.strtype = 'r';
  sharpening.num     = [1 1];
  sharpening.val     = {1};
  sharpening.hidden  = expert<2;
  sharpening.help    = {
    'By smoothing heavily noisy areas, fine structures and local contrasts such as cerebellar sublobuli can disappear.  The effect is similar to a real photo of a meadow at night and short exposure time (e.g., with a lot of noise), where the denoising filter merges everything into one large smooth area.  The sharpening tries to preserve local contrasts. '
    'Sharpening is only applied in case of the optimized filters.'
    ''
  };
  
  % -----------------------------------------------------------------------


  
  
  nlm_default           = cfg_branch;
  nlm_default.tag       = 'classic';
  nlm_default.name      = 'Classic SANLM filter';
  nlm_default.val       = {};
  nlm_default.help      = {
      'Classical SANLM filter without further adaptations, i.e. strong filtering on the full resolution.' 
  }; 

  nlm_optimized         = cfg_branch;
  nlm_optimized.tag     = 'optimized';
  nlm_optimized.name    = 'Optimized filter';
  nlm_optimized.val     = {NCstrm};
  nlm_optimized.help    = {
      'Optimized SANLM filter with predefined parameter settings.' 
  }; 

  nlm_expert          = cfg_branch;
  nlm_expert.tag      = 'expert';
  nlm_expert.name     = 'Optimized filter (expert options)';
  nlm_expert.hidden   = expert<0;
  nlm_expert.val      = {NCstr iter iterm outlier addnoise relativeIntensityAdaption relativeIntensityAdaptionTH relativeFilterStengthLimit resolutionDependency resolutionDependencyRange resolutionReduction lazy};
  nlm_expert.help     = {
      'Optimized SANLM filter with all parameters.' 
  }; 

  nlmfilter             = cfg_choice;
  nlmfilter.tag         = 'nlmfilter';
  nlmfilter.name        = 'Filter type';
  if expert
    nlmfilter.values    = {nlm_default nlm_optimized nlm_expert};
  else
    nlmfilter.values    = {nlm_default nlm_optimized};
  end
  nlmfilter.val         = {nlm_optimized}; 
  if expert
    nlmfilter.help      = {
      'Selection between the classical SANLM filter and an optimized SANLM filter with predefined settings or detailed parameterization. The classic filter is often too strong in normal data that was not interpolated or resampled and the default CAT12 preprocessing uses the medium optimized version. '
      ''
    }; 
  else
    nlmfilter.help      = {
      'Selection between the classical SANLM filter and an optimized SANLM filter. The classic filter is often too strong in normal data that was not interpolated or resampled and the default CAT12 preprocessing uses the medium optimized version. ' 
      ''
    };
  end
  
  % V1
  sanlm                 = cfg_exbranch;
  sanlm.tag             = 'sanlm';
  sanlm.name            = 'Spatially adaptive non-local means (SANLM) denoising filter';
  
  intlim.hidden         = expert<2;
  
  sanlm.val             = {data spm_type prefix suffix intlim rician replaceNANandINF nlmfilter};
  sanlm.prog            = @cat_vol_sanlm;
  sanlm.vout            = @vout_sanlm;
  sanlm.help            = {
    'This function applies an spatial adaptive (sub-resolution) non-local means denoising filter to the data. This filter will remove noise while preserving edges. The filter strength is automatically estimated based on the standard deviation of the noise. '
    ''
    'This filter is internally used in the segmentation procedure anyway. Thus, it is not necessary (and not recommended) to apply the filter before segmentation.'
    ''
  };

  % V2
  sanlm2                = sanlm; 
  sanlm2.tag            = 'sanlm2';
  sanlm2.name           = 'Spatially adaptive non-local means (SANLM) denoising filter V2';
  sanlm2.val            = {data spm_type prefix suffix intlim rician sharpening replaceNANandINF nlmfilter};
  sanlm2.prog           = @cat_vol_sanlm2;
  sanlm2.hidden         = expert<2;

return

%_______________________________________________________________________
function spmtype = conf_io_volctype(data,  intlim,  spm_type,prefix,suffix,verb,expert,lazy)
  % update variables 
  data.help           = {'Select images for data type conversion';''};
  
  intlim.tag          = 'range';
  intlim.num          = [1 inf];
  
  prefix.val          = {'PARA'};
  prefix.help         = {
    'Specify the string to be prepended to the filenames of the converted image file(s). Default prefix is "PARA" that is replaced by the chosen datatype.'
    ''
  };
  
  suffix.hidden       = expert<2;
  suffix.help         = {
    'Specify the string to be prepended to the filenames of the converted image file(s). Default prefix is ''''. Use "PARA" to add the datatype to the filename.'
    ''
  };

  spm_type.labels(1)  = []; % remove native case
  spm_type.values(1)  = []; % remove native case
  spm_type.tag        = 'ctype';
  
  lazy.hidden         = expert<1; 
  verb.hidden         = expert<1; 
  
  intscale            = cfg_menu;
  intscale.name       = 'Intensity scaling';
  intscale.tag        = 'intscale';
  if expert>1
    intscale.labels   = {'No (round in case of integer)','Yes (0:1)','Yes (-1:1)','Yes (0:max)','Yes (min:max)', 'Yes (0:256)'};
    intscale.values   = {0,1,-1,inf,-inf,2};
  else
    intscale.labels   = {'No (round in case of integer)','Yes (0:1)','Yes (-1:1)','Yes (0:max)','Yes (min:max)'};
    intscale.values   = {0,1,-1,inf,-inf};
  end
  intscale.val        = {1};
  intscale.help       = {'Normalize image intensities in range between (i) 0 and 1 (0:1) or (ii) between -1 and 1 (-1:1) balanced around zero, i.e., the absoluted values are ranged betweed 0 and 1. '};

  % new
  spmtype             = cfg_exbranch;
  spmtype.tag         = 'spmtype';
  spmtype.name        = 'Image data type converter'; 
  spmtype.val         = {data prefix suffix intlim spm_type intscale verb lazy};
  spmtype.prog        = @cat_io_volctype;
  spmtype.vout        = @vout_volctype;
  spmtype.help        = {
    'Convert the image data type to reduce disk-space.'
    'Uses 99.99% of the main intensity histogram to avoid problems due to outliers. Although the internal scaling supports a relative high accuracy for the limited number of bits, special values such as NAN and INF will be lost!'
    ''
  };
return

%_______________________________________________________________________
function headtrimming = conf_vol_headtrimming(intlim,spm_type,prefix,suffix,verb,lazy,expert)

  suffix.hidden         = expert<1; 
  intlim.hidden         = expert<1; 
  lazy.hidden           = expert<1;
  intlim.num            = [1 inf]; 
  
  % update input variables
  intlim1               = intlim;
  intlim1.tag           = 'intlim1';
  intlim1.name          = 'Global intensity limitation for masking';
  intlim1.val           = {90};
  intlim1.hidden        = expert<1; 
  intlim1.num           = [1 Inf];
  intlim1.help          = {'General intensity limitation to remove strong outliers by using 90% of the original histogram values. Too high values will include background noise and do not allow trimming, whereas to low values will cut objects with low values (e.g. by image inhomogeneities). ' ''};

  prefix.val            = {'trimmed_'};

  intscale              = cfg_menu;
  intscale.name         = 'Intensity scaling';
  intscale.tag          = 'intscale';
  intscale.labels       = {'No','Yes','Yes (force positive values)'};
  intscale.values       = {0,1,3};
  intscale.val          = {0};
  intscale.help         = {'Normalize image intensities in range between 0 and 1 (unsigned integer or forced scaling option) or -1 to 1 (signed integer / float). '};

  
  % many subjects
  simages               = cfg_files;
  simages.tag           = 'simages';
  simages.name          = 'Source images';
  simages.help          = {'Select images for trimming (e.g. T1 images).' ''};
  simages.filter        = 'image';
  simages.ufilter       = '.*';
  simages.num           = [1 Inf];

  images1               = cfg_files;
  images1.tag           = 'oimages';
  images1.name          = 'Images';
  images1.help          = {'Select other images that should be trimmed similar to the source images (e.g. coregistrated images).' ''};
  images1.filter        = 'image';
  images1.ufilter       = '.*';
  images1.num           = [1 Inf];

  oimages               = cfg_repeat;
  oimages.tag           = 'oimages';
  oimages.name          = 'Other images';
  oimages.help          = {'Select other images that should be trimmed similar to the source images. For example, the source images are a set of T1 images, whereas the second set may be a set of coregistered images of the same subjects with the same image dimensions.' ''};
  oimages.values        = {images1};
  oimages.val           = {};
  oimages.num           = [0 Inf];

  manysubjects          = cfg_branch;
  manysubjects.tag      = 'manysubjects';
  manysubjects.name     = 'Many subjects';
  manysubjects.val      = {simages oimages};
  manysubjects.help     = {
    'Create stacks of images of one class that include the same number of many subjects:'
    '  { {S1T1, S2T1,...} {S1T2, S2T2, ...} ... }'
    ''
  };

  % manyimages 
  subjectimages         = cfg_files;
  subjectimages.tag     = 'subjectimages';
  subjectimages.name    = 'Subject';
  subjectimages.help    = {
    'Select all images of one subject that are in the same space and should be trimmed together.'
  };
  if expert
    subjectimages.help  = [ subjectimages.help; {
    'The first image is used to estimate the trimming. ' 
    ''
    }];
  else
    subjectimages.help  = [ subjectimages.help; {
      'In general the first image is used to estimate the trimming (see "Average images" option). ' 
      ''
    }];
  end
  subjectimages.filter  = 'image';
  subjectimages.ufilter = '.*';
  subjectimages.num     = [1 Inf];

  manyimages            = cfg_repeat;
  manyimages.tag        = 'manyimages';
  manyimages.name       = 'Many images';
  manyimages.help       = {
    'Collect images of each subject that should be trimmed together and are in the same space.' 
    '  { {S1T1, S1T2,...} {S2T1, S2T2, ...}  {S2T1, S2T2 } ... }'
    ''};
  manyimages.values     = {subjectimages};
  manyimages.val        = {};
  manyimages.num        = [1 Inf];

  % image selection type
  timages               = cfg_choice;
  timages.tag           = 'image_selector';
  timages.name          = 'Select type of image selection'; 
  timages.values        = {manyimages manysubjects}; 
  timages.val           = {manyimages};
  timages.help          = {
    'Select "many images" if you have a small number of subjects with a VARYING number of images.'
    'Select "many subjects" if you have a large number of subject with the SAME number of images.'
  };


  pth                   = cfg_entry;
  pth.tag               = 'pth';
  pth.name              = 'Percentual trimming threshold';
  pth.strtype           = 'r';
  pth.num               = [1 1];
  pth.val               = {0.4};
  pth.help              = {'Percentual treshold for trimming. Lower values will result in a wider mask, ie. more air, whereas higher values will remove more air but maybe also brain regions with very low intensity.' ''};

  open                  = cfg_entry;
  open.tag              = 'open';
  open.name             = 'Size of morphological opening of the mask';
  open.strtype          = 'n';
  open.num              = [1 1];
  open.val              = {2};
  open.hidden           = expert<1; 
  open.help             = {'The morphological opening of the mask allows to avoid problems due to noise in the background. However, too large opening will also remove the skull or parts of the brain.' ''};

  addvox                = cfg_entry;
  addvox.tag            = 'addvox';
  addvox.name           = 'Add voxels around mask';
  addvox.strtype        = 'w';
  addvox.num            = [1 1];
  addvox.val            = {2};
  addvox.hidden         = expert<1; 
  addvox.help           = {'Add # voxels around the original mask to avoid to hard masking.' ''};

  mask                  = cfg_menu;
  mask.name             = 'Final masking with source image';
  mask.tag              = 'mask';
  mask.labels           = {'Yes','No'};
  mask.values           = {1,0};
  mask.val              = {1};
  mask.hidden           = false; % expert < 1; % the field is important if multiple images are used and the first one is skull-stripped
                                 % but this stripping/masking should not applied to the other images
  mask.help             = {'Use source image for trimming and final masking (e.g. for skull-stripping in longitudinal pipeline).'};


  % don't change data type
  spm_type.val         = {0};
  spm_type.tag         = 'ctype';
  % --- main ---
  headtrimming         = cfg_exbranch;
  headtrimming.tag     = 'datatrimming';
  headtrimming.name    = 'Image data trimming'; 
  headtrimming.val     = {timages prefix suffix mask intlim1 pth open addvox intlim spm_type intscale verb lazy};
  headtrimming.prog    = @cat_vol_headtrimming;
  headtrimming.vout    = @vout_headtrimming;
  headtrimming.help    = {
    'Remove air around the head and convert the image data type to save disk-space but also to reduce memory-space and load/save times. Corresponding images have to have the same image dimenions. '
    'Uses 99.99% of the main intensity histogram to avoid problems due to outliers. Although the internal scaling supports a relative high accuracy for the limited number of bits, special values such as NAN and INF will be lost!'
    ''
  };
return

%_______________________________________________________________________
function maskimg = conf_vol_maskimage(data,prefix)
 
  % update input variables
  data.name       = 'Select images';
  data.help       = {'Select images for lesion or brain masking';''};
 
  prefix.val      = {'msk_'};
  prefix.help     = {
    'Specify the string to be prepended to the filenames of the masked image file(s).'
    ''
  };

  % lesion mask
  mask = data; 
  mask.tag        = 'mask';
  mask.name       = 'Select lesion mask images';
  mask.help       = {'Select (additional) lesion mask images that describe the regions that should be set to zero.';''};
  mask.num        = [0 Inf];
  
  % brain mask
  bmask = data; 
  bmask.tag       = 'bmask';
  bmask.name      = 'Optionally select additional brain mask images';
  bmask.help      = {'Select (additional) brain mask images that describe the regions that should remain in the image.';''};
  bmask.num       = [0 Inf];
  bmask.val       = {{''}}; 
  
  % recalc
  recalc          = cfg_menu;
  recalc.tag      = 'recalc';
  recalc.name     = 'Reprocess';
  recalc.help     = {'If an output image already exist then use this image rather than the original input image for additional masking. This allows you to add lesions from other lesion images.'};
  recalc.labels   = {'Yes' 'No'};
  recalc.values   = {1 0};
  recalc.val      = {1};
  
  % main
  maskimg         = cfg_exbranch;
  maskimg.tag     = 'maskimg';
  maskimg.name    = 'Manual image (lesion) masking'; 
  maskimg.val     = {data mask bmask recalc prefix};
  maskimg.prog    = @cat_vol_maskimage;
  maskimg.vout    = @vout_maskimg;
  maskimg.help    = {
    'Mask images to avoid segmentation and registration errors in brain lesion. The number of mask images has to be equal to the number of the original images. Voxels inside the lesion mask(s) and outside the brainmask(s) will be set to zero. '
    'If you have multiple lesion masks than add them with the original images, eg. "images = {sub01.nii; sub02.nii; sub01.nii}" and "mask = {sub01_lesion1.nii; sub02_lesion1.nii; sub01_lesion2.nii}". Alternatively, you can choose only one original image and a various number of mask files.'
    ''
  };

return

%_______________________________________________________________________
function [defs,defs2] = conf_vol_defs()

  field           = cfg_files;
  field.tag       = 'field';
  field.name      = 'Deformation Fields';
  field.filter    = 'image';
  field.ufilter   = '^(i)?y_.*\.nii$';
  field.num       = [1 Inf];
  field.help      = {[
    'Select deformation fields for all subjects. ' ...
    'Use the "y_*.nii" to project data from subject to template space, and the "iy_*.nii" to map data from template to individual space. ' ...
    'Both deformation maps can be created in the CAT preprocessing by setting the "Deformation Field" flag to forward or inverse.' ... 
  ]};

  field1          = cfg_files;
  field1.tag      = 'field1';
  field1.name     = 'Deformation Field';
  field1.filter   = 'image';
  field1.ufilter  = '^(i)?y_.*\.nii$';
  field1.num      = [1 1];
  field1.help     = {[
    'Select the deformation field of one subject.' ...
    'Use the "y_*.nii" to project data from subject to template space, and the "iy_*.nii" to map data from template to individual space.' ...
    'Both deformation maps can be created in the CAT preprocessing by setting the "Deformation Field" flag to forward or inverse.' ... 
  ]};

  images1         = cfg_files;
  images1.tag     = 'images';
  images1.name    = 'Images';
  images1.help    = {'Select images to be warped. Note that there should be the same number of images as there are deformation fields, such that each flow field warps one image.'};
  images1.filter  = 'image';
  images1.ufilter = '.*';
  images1.num     = [1 Inf];

  images          = cfg_repeat;
  images.tag      = 'images';
  images.name     = 'Images';
  images.help     = {'The flow field deformations can be applied to multiple images. At this point, you choose how many images each flow field should be applied to.'};
  images.values   = {images1};
  images.num      = [1 Inf];

  interp          = cfg_menu;
  interp.name     = 'Interpolation';
  interp.tag      = 'interp';
  interp.labels   = {
    'Nearest neighbour','Trilinear','2nd Degree B-spline',...
    '3rd Degree B-Spline ','4th Degree B-Spline ','5th Degree B-Spline',...
    '6th Degree B-Spline','7th Degree B-Spline','Categorical'};
  interp.values   = {0,1,2,3,4,5,6,7,-1};
  interp.val      = {1};
  interp.help     = {
    'The method by which the images are sampled when being written in a different space.'
    '    Nearest Neighbour:     - Fastest, but not normally recommended.'
    '    Bilinear Interpolation:     - OK for PET, or realigned fMRI.'
    '    B-spline Interpolation:     - Better quality (but slower) interpolation/* \cite{thevenaz00a}*/, especially with higher degree splines.  Can produce values outside the original range (e.g. small negative values from an originally all positive image). Do not use B-splines when there is any region of NaN or Inf in the images. '
    '    Categorical Interpolation:  - Slow (particularly when there are lots of categories). This is intended to warp categorical images such as label maps.'
  }';

  modulate        = cfg_menu;
  modulate.tag    = 'modulate';
  modulate.name   = 'Modulate image (preserve volume)';
  modulate.labels = {'No','Affine + non-linear (SPM12 default)','Non-linear only'};
  modulate.values = {0 1 2};
  modulate.val    = {0};
  modulate.help   = {
    '"Modulation" is to compensate for the effect of spatial normalisation. Spatial normalisation causes volume changes due to affine transformation (global scaling) and non-linear warping (local volume change). The SPM default is to adjust spatially normalised grey matter (or other tissue class) by using both terms and the resulting modulated images are preserved for the total amount of grey matter. Thus, modulated images reflect the grey matter volumes before spatial normalisation. However, the user is often interested in removing the confound of different brain sizes and there are many ways to apply this correction. We can use the total amount of GM, GM+WM, GM+WM+CSF, or manual estimated total intracranial volume (TIV). Theses parameters can be modeled as nuisance parameters (additive effects) in an AnCova model or used to globally scale the data (multiplicative effects): '
    ''
    '% Correction   Interpretation'
    '% ----------   --------------'
    '% nothing      absolute volume'
    '% globals      relative volume after correcting for total GM or TIV (multiplicative effects)'
    '% AnCova       relative volume that can not be explained by total GM or TIV (additive effects)'
    ''
    'Modulated images can be optionally saved by correcting for non-linear warping only. Volume changes due to affine normalisation will be not considered and this equals the use of default modulation and globally scaling data according to the inverse scaling factor due to affine normalisation. I recommend this option if your hypothesis is about effects of relative volumes which are corrected for different brain sizes. This is a widely used hypothesis and should fit to most data. The idea behind this option is that scaling of affine normalisation is indeed a multiplicative (gain) effect and we rather apply this correction to our data and not to our statistical model. These modulated images are indicated by "m0" instead of "m". '
    ''
  };

  bb             = cfg_entry;
  bb.tag         = 'bb';
  bb.name        = 'Bounding box';
  bb.help        = {'The bounding box (in mm) of the volume which is to be written (relative to the anterior commissure).'};
  bb.strtype     = 'r';
  bb.num         = [2 3];
  bb.val         = {[NaN NaN NaN; NaN NaN NaN]};
  
  vox             = cfg_entry;
  vox.tag         = 'vox';
  vox.name        = 'Voxel sizes';
  vox.help        = {'The voxel sizes (x, y & z, in mm) of the written normalised images.'};
  vox.strtype     = 'r';
  vox.num         = [1 3];
  vox.val         = {[NaN NaN NaN]};

  images1.help    = {'Select images to be warped for this subject.'};
  defs            = cfg_exbranch;
  defs.tag        = 'defs';
  defs.name       = 'Apply deformations (many images)';
  defs.val        = {field1,images1,bb,vox,interp,modulate};
  defs.prog       = @cat_vol_defs;
  defs.vfiles     = @vout_defs;
  defs.help       = {'This is an utility for applying a deformation field of one subject to many images.'};

  defs2           = cfg_exbranch;
  defs2.tag       = 'defs2';
  defs2.name      = 'Apply deformations (many subjects)';
  defs2.val       = {field,images,bb,vox,interp,modulate};
  defs2.prog      = @cat_vol_defs;
  defs2.vfiles    = @vout_defs2;
  defs2.help      = {'This is an utility for applying deformation fields of many subjects to images.'};
return

%_______________________________________________________________________
function realign  = conf_vol_series_align(data)
  
  data.help       = {
  'Select all images for this subject'};

  tim             = cfg_entry;
  tim.tag         = 'times';
  tim.name        = 'Times';
  tim.strtype     = 'e';
  tim.val         = {NaN};
  tim.num         = [1 Inf];
  tim.help        = {'Specify the times of the scans in years. If you leave the default NaN value the standard warping regularization will be used for all scans.'};

  bparam          = cfg_entry;
  bparam.tag      = 'bparam';
  bparam.name     = 'Bias Regularisation';
  bparam.help     = {
    'MR images are usually corrupted by a smooth, spatially varying artifact that modulates the intensity of the image (bias). These artifacts, although not usually a problem for visual inspection, can impede automated processing of the images.'
    ''
    'An important issue relates to the distinction between variations in the difference between the images that arise because of the differential bias artifact due to the physics of MR scanning, and those that arise due to shape differences.  The objective is to model the latter by deformations, while modelling the former with a bias field. We know a priori that intensity variations due to MR physics tend to be spatially smooth. A more accurate estimate of a bias field can be obtained by including prior knowledge about the distribution of the fields likely to be encountered by the correction algorithm. For example, if it is known that there is little or no intensity non-uniformity, then it would be wise to penalise large estimates of the intensity non-uniformity.'
    'Knowing what works best should be a matter of empirical exploration, as it depends on the scans themselves.  For example, if your data has very little of the artifact, then the bias regularisation should be increased.  This effectively tells the algorithm that there is very little bias in your data, so it does not try to model it.'
  }';
  bparam.strtype  = 'e';
  bparam.num      = [1 1];
  bparam.val      = {1e7};

  setCOM        = cfg_menu;
  setCOM.tag    = 'setCOM';
  setCOM.name   = 'Use center-of-mass to set origin';
  setCOM.help   = { ...
      ''
      'Use center-of-mass to roughly correct for differences in the position between image and template. This will internally correct the origin. '
      ''
      'If affine registration fails you can try to disable this option and/or set the origin manually. '
    };
  setCOM.def    = @(val) cat_get_defaults('extopts.setCOM', val{:});
  setCOM.labels = {'No','Yes'};
  setCOM.values = {0 1};


  wparam          = cfg_entry;   
  wparam.tag      = 'wparam';
  wparam.name     = 'Warping Regularisation';
  wparam.help     = {
    'Registration involves simultaneously minimising two terms.  One of these is a measure of similarity between the images (mean-squared difference in the current situation), whereas the other is a measure of the roughness of the deformations.  This measure of roughness involves the sum of the following terms:',...
    '* Absolute displacements need to be penalised by a tiny amount.  The first element encodes the amount of penalty on these.  Ideally, absolute displacements should not be penalised, but it is often necessary for technical reasons.',...
    '* The `membrane energy'' of the deformation is penalised (2nd element), usually by a relatively small amount. This penalises the sum of squares of the derivatives of the velocity field (ie the sum of squares of the elements of the Jacobian tensors).',...
    '* The `bending energy'' is penalised (3rd element). This penalises the sum of squares of the 2nd derivatives of the velocity.',...
    '* Linear elasticity regularisation is also included (4th and 5th elements).  The first parameter (mu) is similar to that for linear elasticity, except it penalises the sum of squares of the Jacobian tensors after they have been made symmetric (by averaging with the transpose).  This term essentially penalises length changes, without penalising rotations.',...
    '* The final term also relates to linear elasticity, and is the weight that denotes how much to penalise changes to the divergence of the velocities (lambda).  This divergence is a measure of the rate of volumetric expansion or contraction.',...
    'Note that regularisation is specified based on what is believed to be appropriate for a year of growth.  The specified values are divided by the number of years time difference.' 
  };
  wparam.strtype  = 'e';
  wparam.num      = [1 5];
  wparam.val      = {[0 0 100 25 100]};
  % Change to (eg): wparam.val     = {[0 0 100 25 12]};

  write_rimg          = cfg_menu;
  write_rimg.tag      = 'write_rimg';
  write_rimg.name     = 'Save rigidly registered images';
  write_rimg.help     = {'Do you want to save the rigidly registered images? The resliced images are named the same as the originals, except that they are prefixed by ''r''.'};
  write_rimg.labels   = {'Save','Dont save'};
  write_rimg.values   = { 1 0 };
  write_rimg.val      = {1};

  write_avg           = cfg_menu;
  write_avg.tag       = 'write_avg';
  write_avg.name      = 'Save Mid-point average';
  write_avg.help      = {'Do you want to save the mid-point average template image? This is likely to be useful for groupwise alignment, and is prefixed by ''avg_'' and written out in the same directory of the first time point data. Please note that with rigid registration a weighted median/mean is stored instead of the average image. In areas with low stdev the mean is used and in areas with larger stdev the median is more weighted.'};
  write_avg.labels    = {'Save','Dont save'};
  write_avg.values    = { 1 0 };
  write_avg.val       = {1};

  write_jac           = cfg_menu;
  write_jac.tag       = 'write_jac';
  write_jac.name      = 'Save Jacobians';
  write_jac.help      = {'Do you want to save a map of the Jacobian determinants?  Some consider these useful for morphometrics (although the divergences of the initial velocities may be preferable). Each map of Jacobians encodes the relative volume (at each spatial location) between the scan and the median time-point average. Values less than one indicate contraction (over time), whereas values greater than one indicate expansion.  These files are prefixed by ``j_'''' and written out in the same directory of the first time point data.'};
  write_jac.labels    = {'Save','Dont save'};
  write_jac.values    = { 1 0 };
  write_jac.val       = {1};

  write_def           = cfg_menu;
  write_def.tag       = 'write_def';
  write_def.name      = 'Deformation Fields';
  write_def.help      = {'Deformation fields can be saved to disk, and used by the Deformations Utility. Deformations are saved as y_*.nii files, which contain three volumes to encode the x, y and z coordinates.  They are written in the same directory as the corresponding image.'};
  write_def.labels    = {'Save','Dont save'};
  write_def.values    = { 1 0 };
  write_def.val       = {0};

  use_brainmask         = cfg_menu;
  use_brainmask.name    = 'Use Brainmask';
  use_brainmask.tag     = 'use_brainmask';
  use_brainmask.labels  = {'Yes','No'};
  use_brainmask.values  = {1,0};
  use_brainmask.val     = {1};
  use_brainmask.help    = {'Use brainmask at last level of rigid body registration to obtain better registration by considering brain regions only.'};

  reduce                = cfg_menu;
  reduce.name           = 'Reduce Bounding Box';
  reduce.tag            = 'reduce';
  reduce.labels         = {'Yes','No'};
  reduce.values         = {1,0};
  reduce.val            = {1};
  reduce.help           = {
    'Reduce bounding box at final resolution level because usually there is a lot of air around the head after registration of multiple scans. This helps to save memory and time for later use of these registered images.'
    ''
    'Please note that this option can only be used for rigid registration and will be disabled for non-linear registration.'
  };

  nonlin                = cfg_branch;
  nonlin.tag            = 'nonlin';
  nonlin.name           = 'Non-linear registration';
  nonlin.val            = {tim wparam write_jac write_def};
  nonlin.help           = {''};

  rigid                 = cfg_const;
  rigid.tag             = 'rigid';
  rigid.name            = 'Rigid body registration';
  rigid.val             = {1};
  rigid.help            = {'Rigid registration only'};

  reg                   = cfg_choice;
  reg.name              = 'Registration Method';
  reg.tag               = 'reg';
  reg.values            = {rigid nonlin};
  reg.val               = {rigid};
  reg.help              = {'Choose between rigid body and non-linear registration. The non-linear registration is using the methods of the Longitudinal Toolbox and is thought for data over longer periods, where the deformations can then be used to calculate local volume changes, which are then multiplied (modulated) by the segmented mean image. Rigid body registration can be used to detect more subtle effects over shorter periods of time (e.g. brain plasticity or training effects after a few weeks or even shorter times).'};

  noise                 = cfg_entry;
  noise.tag             = 'noise';
  noise.name            = 'Noise Estimate';
  noise.help            = {'.'};
  noise.strtype         = 'e';
  noise.num             = [Inf Inf];
  noise.val             = {NaN};
  noise.help            = {'Specify the standard deviation of the noise in the images.  If a scalar is entered, all images will be assumed to have the same level of noise.  For any non-finite values, the algorithm will try to estimate the noise from fitting a mixture of two Rician distributions to the intensity histogram of each of the images, and assuming that the Rician with the smaller overall intensity models the intensity distribution of air in the background. This works reasonably well for simple MRI scans, but less well for derived images (such as averages) and it fails badly for scans that are skull-stripped.  The assumption used by the registration is that the residuals, after fitting the model, are i.i.d. Gaussian. The assumed standard deviation of the residuals is derived from the estimated Rician distribution of the air.'
  };

  realign               = cfg_exbranch;
  realign.tag           = 'series';
  realign.name          = 'Longitudinal Registration';
  realign.val           = {data noise setCOM bparam use_brainmask reduce reg write_rimg write_avg};
  realign.help          = {
    'Longitudinal registration of series of anatomical MRI scans for a single subject.  It is based on inverse-consistent alignment among each of the subject''s scans, and incorporates a bias field correction.  Prior to running the registration, the scans should already be in very rough alignment, although because the model incorporates a rigid-body transform, this need not be extremely precise.  Note that there are a bunch of hyper-parameters to be specified.  If you are unsure what values to take, then the defaults should be a reasonable guess of what works.  Note that changes to these hyper-parameters will impact the results obtained.'
    ''
    'The alignment assumes that all scans have similar resolutions and dimensions, and were collected on the same (or very similar) MR scanner using the same pulse sequence.  If these assumption are not correct, then the approach will not work as well.  There are a number of settings (noise estimate, regularisation etc). Default settings often work well, but it can be very helpful to try some different values, as these can have a large effect on the results.'
    ''
    'The resliced images are named the same as the originals, except that they are prefixed by ''r''.'
  };
  realign.prog          = @cat_vol_series_align;
  realign.vout          = @vout_realign;

return

%_______________________________________________________________________
function [T2x,T2x_surf,F2x,F2x_surf] = conf_T2x

  data_T2x          = cfg_files;
  data_T2x.tag      = 'data_T2x';
  data_T2x.name     = 'Data';
  data_T2x.filter   = {'image'};
  data_T2x.ufilter  = '^spmT.*';
  data_T2x.num      = [1 Inf];
  data_T2x.help     = {'Select spmT-data to transform or convert.'};

  sel               = cfg_menu;
  sel.name          = 'Convert t value to';
  sel.tag           = 'sel';
  sel.labels        = {'p','-log(p)','correlation coefficient cc','standard Normal (z-score) distribution','apply thresholds without conversion'};
  sel.values        = {1,2,3,6,5};
  sel.val           = {2};
  sel.help          = {'Select conversion of t-value'};

  thresh05          = cfg_entry;
  thresh05.tag      = 'thresh05';
  thresh05.name     = 'Threshold';
  thresh05.help     = {''};
  thresh05.strtype  = 'r';
  thresh05.num      = [1 1];
  thresh05.val      = {0.05};

  thresh001         = cfg_entry;
  thresh001.tag     = 'thresh001';
  thresh001.name    = 'Threshold';
  thresh001.help    = {''};
  thresh001.strtype = 'r';
  thresh001.num     = [1 1];
  thresh001.val     = {0.001};

  kthresh           = cfg_entry;
  kthresh.tag       = 'kthresh';
  kthresh.name      = 'Extent (voxels)';
  kthresh.help      = {'Enter the extent threshold in voxels'};
  kthresh.strtype   = 'r';
  kthresh.val       = {0};
  kthresh.num       = [1 1];

  noniso            = cfg_menu;
  noniso.name       = 'Correct for non-isotropic smoothness';
  noniso.tag        = 'noniso';
  noniso.labels     = {'Yes','No'};
  noniso.values     = {1,0};
  noniso.val        = {1};
  noniso.help       = {'Correct for non-isotropic smoothness for cluster extent thresholds.'};

  none              = cfg_const;
  none.tag          = 'none';
  none.name         = 'None';
  none.val          = {1};
  none.help         = {'No threshold'};

  k                 = cfg_branch;
  k.tag             = 'k';
  k.name            = 'k-value';
  k.val             = {kthresh, noniso };
  k.help            = {''};

  fwe               = cfg_branch;
  fwe.tag           = 'fwe';
  fwe.name          = 'FWE';
  fwe.val           = {thresh05 };
  fwe.help          = {''};

  fdr               = cfg_branch;
  fdr.tag           = 'fdr';
  fdr.name          = 'FDR';
  fdr.val           = {thresh05 };
  fdr.help          = {''};

  fwe2              = cfg_branch;
  fwe2.tag          = 'fwe2';
  fwe2.name         = 'FWE';
  fwe2.val          = {thresh05, noniso };
  fwe2.help         = {''};

  uncorr            = cfg_branch;
  uncorr.tag        = 'uncorr';
  uncorr.name       = 'uncorrected';
  uncorr.val        = {thresh001 };
  uncorr.help       = {''};

  kuncorr           = cfg_branch;
  kuncorr.tag       = 'kuncorr';
  kuncorr.name      = 'uncorrected';
  kuncorr.val       = {thresh05, noniso };
  kuncorr.help      = {''};

  En                = cfg_branch;
  En.tag            = 'En';
  En.name           = 'Expected voxels per cluster';
  En.val            = {noniso };
  En.help           = {''};

  inverse           = cfg_menu;
  inverse.name      = 'Show also inverse effects (e.g. neg. values)';
  inverse.tag       = 'inverse';
  inverse.labels    = {'Yes','No'};
  inverse.values    = {1,0};
  inverse.val       = {0};
  inverse.help      = {'Show also inverse effects (e.g. neg. values). This is not valid if you convert to (log) p-values.'};

  threshdesc        = cfg_choice;
  threshdesc.name   = 'Threshold type peak-level';
  threshdesc.tag    = 'threshdesc';
  threshdesc.values = {none uncorr fdr fwe};
  threshdesc.val    = {uncorr};
  threshdesc.help   = {'Select method for voxel threshold'};

  cluster           = cfg_choice;
  cluster.name      = 'Cluster extent threshold';
  cluster.tag       = 'cluster';
  cluster.values    = {none k En kuncorr fwe2};
  cluster.val       = {none};
  cluster.help      = {'Select method for extent threshold'};

  conversion        = cfg_branch;
  conversion.tag    = 'conversion';
  conversion.name   = 'Conversion';
  conversion.val    = {sel threshdesc inverse cluster};
  conversion.help   = {''};

  atlas             = cfg_menu;
  atlas.name        = 'Atlas Labeling';
  atlas.tag         = 'atlas';
  atlas.labels{1}   = 'None';
  atlas.values{1}   = 'None';
  list              = spm_atlas('List','installed');
  for i=1:numel(list)
    atlas.labels{i+1} = list(i).name;
    atlas.values{i+1} = list(i).name;
  end
  atlas.val         = {'None'};
  atlas.help        = {
    'Select atlas for labeling. The prepending atlas name ''dartel_'' indicates that this atlas was created using Dartel spatial normalization with the Dartel IXI template as default.'
    ''
    'Please note, that you can install additional atlases for CAT12 using the command ''cat_install_atlases''. '
  };



  % T2x volumes
  % -----------------------------------------------------------------------
  T2x      = cfg_exbranch;
  T2x.tag  = 'T2x';
  T2x.name = 'Threshold and transform spmT images';
  T2x.val  = {data_T2x,conversion,atlas};
  T2x.prog = @cat_stat_spm2x;
  T2x.vout = @vout_stat_spm2x;
  T2x.help = {
    'This function transforms t-maps to P, -log(P), r, or z-maps.'
    'The following formulas are used:'
    '--------------------------------'
    'correlation coefficient:'
    '            t'
    '  r = ------------------'
    '        sqrt(t^2 + df)'
    'p-value:'
    '  p = 1-spm_Tcdf'
    'log p-value:'
    '  -log10(1-P) = -log(1-spm_Tcdf)'
    'For the latter case of log transformation this means that a p-value of p=0.99 (0.01) is transformed to a value of 2.'
    'Examples:'
    'p-value  -log10(1-P)'
    '0.1      1'
    '0.05     1.30103 (-log10(0.05))'
    '0.01     2'
    '0.001    3'
    '0.0001   4'
    'All maps can be thresholded using height and extent thresholds and you can also apply corrections for multiple comparisons based on family-wise error (FWE) or false discovery rate (FDR). You can easily threshold and/or transform a large number of spmT-maps using the same thresholds.'
    'Naming convention of the transformed files:'
    '   Type_Contrast_Pheight_Pextent_K_Neg'
    '   Type:      P    - p-value'
    '              logP - log p-value'
    '              R    - correlation coefficient'
    '              T    - t-value'
    '              Z    - z-value'
    '   Contrast:  name used in the contrast manager with replaced none valid'
    '              strings'
    '   Pheight:   p    - uncorrected p-value in % (p<0.05 will coded with "p5")'
    '              pFWE - p-value with FWE correction in %'
    '              pFDR - p-value with FDR correction in %'
    '   Pextent:   pk    - uncorr. extent p-value in % (p<0.05 coded with "p5")'
    '              pkFWE - extent p-value with FWE correction in %'
    '   K:         extent threshold in voxels'
    '   Neg:       image also shows thresholded inverse effects (e.g. neg. '
    '              values) '
  }';


  % T2x surfaces
  % -----------------------------------------------------------------------
  
  % Do not use 3D atlases for surfaces
  data_T2x.filter = {'gifti'};
  
  % surfaces
  T2x_surf        = T2x;
  T2x_surf.val    = {data_T2x,conversion};
  T2x_surf.tag    = 'T2x_surf';
  T2x_surf.name   = 'Threshold and transform spmT surfaces';
  T2x_surf.vout   = @vout_stat_spm2x_surf;
  
  
  % F2x volumes
  % -----------------------------------------------------------------------

  data_F2x          = cfg_files;
  data_F2x.tag      = 'data_F2x';
  data_F2x.name     = 'Data';
  data_F2x.filter   = {'image'};
  data_F2x.ufilter  = '^spmF.*';
  data_F2x.num      = [1 Inf];
  data_F2x.help     = {'Select spmF-data to select.'};

  sel               = cfg_menu;
  sel.name          = 'Convert F value to';
  sel.tag           = 'sel';
  sel.labels        = {'p','-log(p)','coefficient of determination R^2','apply thresholds without conversion'};
  sel.values        = {1,2,3,4};
  sel.val           = {2};
  sel.help          = {'Select conversion of F-value'};

  none              = cfg_const;
  none.tag          = 'none';
  none.name         = 'None';
  none.val          = {1};
  none.help         = {'No threshold'};

  cluster           = cfg_choice;
  cluster.name      = 'Cluster extent threshold';
  cluster.tag       = 'cluster';
  cluster.values    = {none k En kuncorr fwe2};
  cluster.val       = {none};
  cluster.help      = {'Select method for extent threshold'};

  conversion        = cfg_branch;
  conversion.tag    = 'conversion';
  conversion.name   = 'Conversion';
  conversion.val    = {sel threshdesc cluster};
  conversion.help   = {''};

  F2x               = cfg_exbranch;
  F2x.tag           = 'F2x';
  F2x.name          = 'Threshold and transform spmF images';
  F2x.val           = {data_F2x,conversion,atlas};
  F2x.prog          = @cat_stat_spm2x;
  F2x.vout          = @vout_stat_spm2x;
  F2x.help          = {
    'This function transforms F-maps to P, -log(P), or R2-maps.'
    'The following formulas are used:'
    '--------------------------------'
    'coefficient of determination R2:'
    '               1'
    '  R2 = ------------------'
    '        1 + F*(p-1)/n-p)'
    'p-value:'
    '  p = 1-spm_Fcdf'
    'log p-value:'
    '  -log10(1-P) = -log(1-spm_Fcdf)'
    'For the last case of log transformation this means that a p-value of p=0.99 (0.01) is transformed to a value of 2.'
    'Examples:'
    'p-value  -log10(1-P)'
    '0.1      1'
    '0.05     1.30103 (-log10(0.05))'
    '0.01     2'
    '0.001    3'
    '0.0001   4'
    'All maps can be thresholded using height and extent thresholds and you can also apply corrections for multiple comparisons based on family-wise error (FWE) or false discovery rate (FDR). You can easily threshold and/or transform a large number of spmT-maps using the same thresholds.'
    'Naming convention of the transformed files:'
    '   Type_Contrast_Pheight_K'
    '   Type:      P    - p-value'
    '              logP - log p-value'
    '              R2   - coefficient of determination'
    '   Contrast:  name used in the contrast manager with replaced none valid'
    '              strings'
    '   Pheight:   p    - uncorrected p-value in % (p<0.05 will coded with "p5")'
    '              pFWE - p-value with FWE correction in %'
    '              pFDR - p-value with FDR correction in %'
    '   K:         extent threshold in voxels'
  }';


  % F2x surfaces
  % -----------------------------------------------------------------------

  % Do not use 3D atlases for surfaces
  data_F2x.filter  = {'gifti'};
  
  F2x_surf         = F2x;
  F2x_surf.val     = {data_F2x,conversion};
  F2x_surf.tag     = 'F2x_surf';
  F2x_surf.name    = 'Threshold and transform spmF surfaces';
  F2x_surf.vout    = @vout_stat_spm2x_surf;

return

%_______________________________________________________________________
function showslice = conf_stat_showslice_all(data_vol)
  data_vol.help = {'Select all images. Images have to be in the same orientation with same voxel size and dimension (e.g. normalized images)'};

  scale           = cfg_menu;
  scale.tag       = 'scale';
  scale.name      = 'Proportional scaling?';
  scale.labels    = {'No','Yes'};
  scale.values    = {0 1};
  scale.val       = {0};
  scale.help      = {'This option should be only used if image intensity is not scaled (e.g. T1 images) or if images have to be scaled during statistical analysis (e.g. modulated images).'};

  orient          = cfg_menu;
  orient.tag      = 'orient';
  orient.name     = 'Spatial orientation';
  orient.labels   = {'axial','coronal','sagittal'};
  orient.values   = {3 2 1};
  orient.val      = {3};
  orient.help     = {'Spatial orientation of slice.'};

  slice           = cfg_entry;
  slice.tag       = 'slice';
  slice.name      = 'Selected slice (in mm)?';
  slice.strtype   = 'r';
  slice.num       = [1 1];
  slice.val       = {0};
  slice.help      = {'Choose slice in mm.'};

  showslice       = cfg_exbranch;
  showslice.tag   = 'showslice';
  showslice.name  = 'Display one slice for all images';
  showslice.val   = {data_vol,scale,orient,slice};
  showslice.prog  = @cat_stat_showslice_all;
  showslice.help  = {'This function displays a selected slice for all images and indicates the respective filenames which is useful to check image quality for a large number of files in a circumscribed region (slice).'};

%_______________________________________________________________________
function quality_measures = conf_quality_measures
  
  data          = cfg_files;
  data.name     = 'Sample data';
  data.tag      = 'data';
  data.filter   = {'image','mesh'};
  data.num      = [1 Inf];
  data.help     = {'These are the (spatially registered or resampled) data. They must all have the same data dimension, orientation, voxel or mesh size etc. Furthermore, it is recommended to use unsmoothed files.'};

  globals        = cfg_menu;
  globals.tag    = 'globals';
  globals.name   = 'Global scaling with TIV';
  globals.labels = {'Yes', 'No'};
  globals.values = {1 0};
  globals.val    = {0};
  globals.help    = {
    'This option is to correct mean z-scores for TIV by global scaling. It is only meaningful for VBM data.'
    ''
  };

  csv_name         = cfg_entry;
  csv_name.tag     = 'csv_name';
  csv_name.name    = 'Output csv file';
  csv_name.strtype = 's';
  csv_name.num     = [1 Inf];
  csv_name.val     = {'Quality_measures.csv'};
  csv_name.help    = {
    'The output file is written to current working directory unless a valid full pathname is given. The following parameters are saved:'
    '  Mean z-score - low values indicate more similarity/homogeneity to sample'
    '  Weighted overall image quality (IQR) - low values mean better image quality before preprocessing'
    '  Normalized product of IQR and Mean z-score - low values point to good image quality before preprocessing and large homogeneity to sample after preprocessing'
    '  Euler Number (for surfaces only) - lower numbers point to better quality of surface extraction'
    '  Size of topology defects (for surfaces only) - smaller size points to better quality of surface extraction'
    ''
    };

  quality_measures         = cfg_exbranch;
  quality_measures.tag     = 'quality_measures';
  quality_measures.name    = 'Check sample homogeneity for very large samples using mean z-score';
  quality_measures.val     = {data,globals,csv_name};
  quality_measures.prog    = @cat_stat_quality_measures;
  quality_measures.help    = {
    'In order to identify data with poor image quality or even artefacts you can use this function. In contrast to the Check Homogeneity tool this function can be also applied to very large samples, but provides no graphical output.'
    'The saved quality parameters in the csv-file can be then used with external analysis tools. The following parameters are saved:'
    '  Mean z-score - low values indicate more similarity/homogeneity to sample'
    '  Weighted overall image quality (IQR) - low values mean better image quality before preprocessing'
    '  Normalized product of IQR and Mean z-score - low values point to good image quality before preprocessing and large homogeneity to sample after preprocessing'
    '  Euler Number (for surfaces only) - lower numbers point to better quality of surface extraction'
    '  Size of topology defects (for surfaces only) - smaller size points to better quality of surface extraction'
    ''
  };

%_______________________________________________________________________
function [check_cov, check_cov2] = conf_check_cov(data_xml,outdir,fname,save,expert) 
 
  % --- update input data ---
  data_xml.name     = 'Quality measures (leave emtpy for autom. search)';
  data_xml.help     = {
    'Select optional the quality measures that are saved during segmentation as xml-files in the report folder. This allows to additionally analyze image quality parameters such as noise, bias, and weighted overall image quality.'
    'Please note, that the order of the xml-files should be the same as the other data files.'
    'Leave empty for automatically search for these xml-files.'
    };
  
  % --- further data ---
  c                 = cfg_entry;
  c.tag             = 'c';
  c.name            = 'Vector/Matrix';
  c.help            = {'Vector or matrix of nuisance values'};
  c.strtype         = 'r';
  c.num             = [Inf Inf];

  nuisance          = cfg_repeat;
  nuisance.tag      = 'nuisance';
  nuisance.name     = 'Nuisance variable';
  nuisance.values   = {c};
  nuisance.num      = [0 Inf];
  nuisance.help     = {'This option allows for the specification of nuisance effects to be removed from the data. A potential nuisance parameter can be TIV if you check segmented data with the default modulation. In this case the variance explained by TIV will be removed prior to the calculation of the correlation. Another meaningful nuisance effect is age. This parameter should be defined for all samples as one variable and may also contain several columns.'};

  gap               = cfg_entry;
  gap.tag           = 'gap';
  gap.name          = 'Separation';
  gap.strtype       = 'n';
  gap.num           = [1 1];
  gap.val           = {3};
  gap.hidden        = expert<2;
  gap.help          = {'To speed up calculations you can define that covariance is estimated only every x voxel. Smaller values give slightly more accurate covariance, but will be much slower.'};

  data_vol          = cfg_files;
  data_vol.name     = 'Sample data';
  data_vol.tag      = 'data_vol';
  data_vol.filter   = {'image','resampled.*\.gii$'};
  data_vol.num      = [1 Inf];
  data_vol.help     = {'These are the (spatially registered or resampled) data. Volumes must all have the same image dimensions, orientation, voxel size and surfaces should be resampled with the same parameters. Furthermore, it is recommended to use unsmoothed files (i.e. for volumes).'};

  sample            = cfg_repeat;
  sample.tag        = 'sample';
  sample.name       = 'Data';
  sample.values     = {data_vol};
  sample.num        = [1 Inf];
  sample.help       = {'Specify data for each sample. If you specify different samples the mean correlation is displayed in separate boxplots (or violin plots) for each sample.'};


  check_cov         = cfg_exbranch;
  check_cov.tag     = 'check_cov';
  check_cov.name    = 'Check Sample Homogeneity';
  if expert>1
    check_cov.val     = {sample,data_xml,gap,nuisance,outdir,fname,save};
  else
    check_cov.val     = {sample,data_xml,gap,nuisance};
  end
  check_cov.prog    = @cat_stat_check_cov;
  check_cov.help    = {
    'In order to identify data with poor data quality or even artefacts you can use this function. 3D images have to be in the same orientation with same voxel size and dimension (e.g. normalized images without smoothing) while surfaces have to be resampled and smoothed using the same parameters. The idea of this tool is to check the correlation of all data across the sample.'
    ''
    'The correlation is calculated between all data and the mean for each data is plotted using a boxplot and the indicated filenames. The smaller the mean correlation the more deviant is this data from the sample mean. In the plot, outliers from the sample are usually isolated from the majority of data which are clustered around the sample mean. The mean correlation is plotted at the y-axis and the x-axis reflects the data order.'
    'If you have loaded quality measures, you can also display the ratio between weighted overall image quality (IQR) and mean correlation. These two are the most important measures for assessing data quality.'
  };


  % --- main ---
  check_cov2        = check_cov; 
  check_cov2.tag    = 'check_cov2';
  check_cov2.name   = 'Check sample homogeneity of 3D data (new exp. version)';
  check_cov2.val    = {sample,gap,nuisance};
  check_cov2.prog   = @cat_stat_check_cov2;

%_______________________________________________________________________
function check_SPM = conf_stat_check_SPM(outdir,fname,save,expert) 

  outdir.hidden               = expert<2;
  fname.hidden                = expert<2; 
  save.hidden                 = expert<2;

  spmmat                      = cfg_files;
  spmmat.tag                  = 'spmmat';
  spmmat.name                 = 'Select SPM.mat';
  spmmat.filter               = {'mat'};
  spmmat.ufilter              = '^SPM\.mat$';
  spmmat.num                  = [1 1];
  spmmat.help                 = {'Select the SPM.mat file that contains the design specification.'};

  % check_SPM_cov
  use_unsmoothed_data         = cfg_menu;
  use_unsmoothed_data.name    = 'Use unsmoothed data if found';
  use_unsmoothed_data.tag     = 'use_unsmoothed_data';
  use_unsmoothed_data.labels  = {'Yes','No'};
  use_unsmoothed_data.values  = {1,0};
  use_unsmoothed_data.val     = {1};
  use_unsmoothed_data.help    = {'Check for sample homogeneity results in more reliable values if unsmoothed data are used. Unsmoothed data contain more detailed information about differences and similarities between the data.'};

  adjust_data                 = cfg_menu;
  adjust_data.name            = 'Adjust data using design matrix';
  adjust_data.tag             = 'adjust_data';
  adjust_data.labels          = {'Yes','No'};
  adjust_data.values          = {1,0};
  adjust_data.val             = {1};
  adjust_data.help            = {'This option allows to use nuisance and group parameters from the design matrix to obtain adjusted data. In this case the variance explained by these parameters will be removed prior to the calculation of the correlation. Furthermore, global scaling (if defined) is also applied to the data.'};

  do_check_cov                = cfg_branch;
  do_check_cov.tag            = 'do_check_cov';
  do_check_cov.name           = 'Yes';
  do_check_cov.val            = {use_unsmoothed_data adjust_data ,outdir,fname,save};
  do_check_cov.help           = {''};

  none                        = cfg_const;
  none.tag                    = 'none';
  none.name                   = 'No';
  none.val                    = {1};
  none.help                   = {''};

  check_SPM_cov               = cfg_choice;
  check_SPM_cov.name          = 'Check for sample homogeneity';
  check_SPM_cov.tag           = 'check_SPM_cov';
  check_SPM_cov.values        = {none do_check_cov};
  check_SPM_cov.val           = {do_check_cov};
  check_SPM_cov.help          = {
    'In order to identify images with poor image quality or even artefacts you can use this function. The idea of this tool is to check the correlation of all files across the sample using the files that are already defined in SPM.mat.'
    ''
    'The correlation is calculated between all images and the mean for each image is plotted using a boxplot (or violin plot) and the indicated filenames. The smaller the mean correlation the more deviant is this image from the sample mean. In the plot outliers from the sample are usually isolated from the majority of images which are clustered around the sample mean. The mean correlation is plotted at the y-axis and the x-axis reflects the image order'
  };

  check_SPM_ortho             = cfg_menu;
  check_SPM_ortho.name        = 'Check for design orthogonality';
  check_SPM_ortho.tag         = 'check_SPM_ortho';
  check_SPM_ortho.labels      = {'Yes','No'};
  check_SPM_ortho.values      = {1,0};
  check_SPM_ortho.val         = {1};
  check_SPM_ortho.help        = {'Review Design Orthogonality.'};

  check_SPM                   = cfg_exbranch;
  check_SPM.tag               = 'check_SPM';
  check_SPM.name              = 'Check design orthogonality and homogeneity';
  check_SPM.val               = {spmmat,check_SPM_cov,check_SPM_ortho};
  check_SPM.prog              = @cat_stat_check_SPM;
  check_SPM.help              = {'Use design matrix saved in SPM.mat to check for sample homogeneity of the used data and for orthogonality of parameters.'};

%_______________________________________________________________________
function calcvol = conf_stat_TIV
  calcvol_name         = cfg_entry;
  calcvol_name.tag     = 'calcvol_name';
  calcvol_name.name    = 'Output file';
  calcvol_name.strtype = 's';
  calcvol_name.num     = [1 Inf];
  calcvol_name.val     = {'TIV.txt'};
  calcvol_name.help    = {
  'The output file is written to current working directory unless a valid full pathname is given.'};

  calcvol_TIV         = cfg_menu;
  calcvol_TIV.tag     = 'calcvol_TIV';
  calcvol_TIV.name    = 'Save values';
  calcvol_TIV.labels  = {'TIV only' 'TIV & GM/WM/CSF/WMH'};
  calcvol_TIV.values  = {1 0};
  calcvol_TIV.val     = {1};
  calcvol_TIV.help    = {'You can save either the total intracranial volume (TIV) only or additionally also save the global volumes for GM, WM, CSF, and WM hyperintensities.'
  ''
  };

  calcvol_savenames         = cfg_menu;
  calcvol_savenames.tag     = 'calcvol_savenames';
  calcvol_savenames.name    = 'Add filenames';
  calcvol_savenames.labels  = {'Values only' 'Add file names' 'Add folders and file names'};
  calcvol_savenames.values  = {0 1 2};
  calcvol_savenames.val     = {0};
  calcvol_savenames.help    = {'You can either save only the values (that can be easily read with spm_load) or also add file names (and folders) to 1st column.'
  ''
  };

  clear data_xml
  data_xml = cfg_files;
  data_xml.name = 'XML files';
  data_xml.tag  = 'data_xml';
  data_xml.filter = 'xml';
  data_xml.ufilter = '^cat_.*\.xml$';
  data_xml.num  = [1 Inf];
  data_xml.help   = {...
  'Select xml-files that are saved during segmentation in the report folder.'};

  calcvol       = cfg_exbranch;
  calcvol.tag   = 'calcvol';
  calcvol.name  = 'Estimate TIV and global tissue volumes';
  calcvol.val   = {data_xml,calcvol_TIV,calcvol_savenames,calcvol_name};
  calcvol.prog  = @cat_stat_TIV;
  calcvol.vout  = @vout_stat_TIV;
  calcvol.help  = {
  'This function reads raw volumes for TIV/GM/WM/CSF/WM hyperintensities (WMH) and saves values in a txt-file. These values can be read with the matlab command: vol = spm_load. If you choose to save all values the entries for TIV/GM/WM/CSF/WMH are now saved in vol(:,1) vol(:,2) vol(:,3), vol(:,4), and vol(:,5) respectively.'
  ''
  'You can use TIV either as nuisance in an AnCova model or as user-specified globals with the "global calculation" option depending on your hypothesis. The use of TIV as nuisance or globals is recommended for modulated data where both the affine transformation and the non-linear warping of the registration are corrected for. '
  ''
  };

%_______________________________________________________________________
function calcroi = conf_roi_fun(outdir)   
  roi_xml               = cfg_files;
  roi_xml.name          = 'XML files';
  roi_xml.tag           = 'roi_xml';
  roi_xml.filter        = 'xml';
  roi_xml.ufilter       = '^catROI.*\.xml$';
  roi_xml.num           = [1 Inf];
  roi_xml.help          = {'These are the xml-files that are saved in the label folder after CAT12 segmentation.'};

  % NOT USED
  %{
  usefolder             = cfg_menu;
  usefolder.tag         = 'folder';
  usefolder.name        = 'Use foldername';
  usefolder.labels      = {'Yes' 'No'};
  usefolder.values      = {1 0};
  usefolder.val         = {0};
  usefolder.help        = {'Use foldername to describe the subject.'};
  %}
  
  point                 = cfg_menu;
  point.tag             = 'point';
  point.name            = 'Decimal point';
  point.labels          = {',','.'};
  point.values          = {',','.'};
  point.val             = {'.'};
  point.help            = {'Decimal point.'}; % that has to be unequal to the column delimiter.'};

  % tab "\t" does not work and so we automatically switch in case of decimal 
  % point "," to delimiter ";".
  %{
  delimiter             = cfg_menu;
  delimiter.tag         = 'delimiter';
  delimiter.name        = 'column delimiter';
  delimiter.labels      = {',',';',' '};
  delimiter.values      = {',',';',' '};
  delimiter.val         = {','};
  delimiter.help        = {'Delimiter between columns.'};
  %}

  calcroi_name          = cfg_entry;
  calcroi_name.tag      = 'calcroi_name';
  calcroi_name.name     = 'Output file';
  calcroi_name.strtype  = 's';
  calcroi_name.num      = [1 Inf];
  calcroi_name.val      = {'ROI'};
  calcroi_name.help     = {'The mean volume values in mL (e.g. GM volume) or the mean surface values (e.g. thickness) are written to the current working directory unless a valid full pathname is given. The output file will also include the name of the atlas and the measure (e.g. Vgm). The file is using tabstops to separate values in order to easily import the file into Excel or SPSS or any other software for subsequent analysis.'};

  calcroi               = cfg_exbranch;
  calcroi.tag           = 'calcroi';
  calcroi.name          = 'Estimate mean/volume inside ROI';
  calcroi.val           = {roi_xml,point,outdir,calcroi_name}; 
  %calcroi.val   = {roi_xml,usefolder,point,outdir,calcroi_name}; % usefolder is never used
  calcroi.prog          = @(job)cat_roi_fun('exportSample',job);
  calcroi.help          = {
    'This function reads values inside a ROI from different atlases (that were selected for CAT12 segmentation) and saves either the mean volume values in mL (e.g. GM volume) or the mean surface values (e.g. thickness) for all data in a csv-file. '
    'Missed values are replaced by NaN.'
  };

%_______________________________________________________________________
function urqio = conf_vol_urqio
%  ------------------------------------------------------------------------
%  Ultra-High Resolution Quantitative Image Optimization
%  ------------------------------------------------------------------------

  % -- Data ---
  r1          = cfg_files; 
  r1.tag      = 'r1';
  r1.name     = 'R1-Volumes';
  r1.filter   = 'image';
  r1.ufilter  = '.*';
  r1.num      = [1 Inf];
  r1.help     = {'Select R1 weighted images.'};

  pd          = cfg_files; 
  pd.tag      = 'pd';
  pd.name     = 'PD-Volumes';
  pd.filter   = 'image';
  pd.ufilter  = '.*';
  pd.num      = [1 Inf];
  pd.help     = {'Select PD weighted images.'};

  r2s         = cfg_files; 
  r2s.tag     = 'r2s';
  r2s.name    = 'R2s-Volumes';
  r2s.filter  = 'image';
  r2s.ufilter = '.*';
  r2s.num     = [1 Inf];
  r2s.help    = {'Select R2s weighted images.'};

  data        = cfg_branch;
  data.tag    = 'data';
  data.name   = 'Input data';
  data.val    = {r1 pd r2s}; 
  data.help   = {
    'Input Images.'
  };

  
  % --- Parameter ---
  spm         = cfg_menu;
  spm.tag     = 'spm';
  spm.name    = 'Use SPM Preprocessing';
  spm.labels  = {'No','Yes'};
  spm.values  = {0 1};
  spm.val     = {1};
  spm.help    = {
    'Use SPM preprocessing if the data is not skull-stripped.'
  };
  
  bc          = cfg_menu;
  bc.tag      = 'bc';
  bc.name     = 'Bias Correction';
  bc.labels   = {'No','light','medium','strong'};
  bc.values   = {0 0.5 1 2};
  bc.val      = {1};
  bc.help     = {
    'Additional bias correction that is important for detection and correction of blood vessels.'
    ''
    'The correction uses a simple tissue classification and local filter approaches to estimate the local signal intensity in the WM and GM segment, e.g. a minimum/maximum filter in the WM for PD/T1 images.  Next, unclassified voxels were approximated and smoothed depending on the defined strength.  '
    ''
  };

  in          = cfg_menu;
  in.tag      = 'in';
  in.name     = 'Intensity Normalization';
  in.labels   = {'No','Yes'};
  in.values   = {0 1};
  in.val      = {1};
  in.help     = {
    'Additional global intensity normalization that is also important for detection and correction of blood vessels.'
    ''
  };

  bvc         = cfg_menu;
  bvc.tag     = 'bvc';
  bvc.name    = 'Blood Vessel Correction';
  bvc.labels  = {'No','Yes'};
  bvc.values  = {0 1};
  bvc.val     = {1};
  bvc.help    = {
    'Correction of blood vessels with high intensity in T1/R1/R2s and low intensity in PD images by CSF-like intensities. '
    ''
  };

  ss          = cfg_menu;
  ss.tag      = 'ss';
  ss.name     = 'Apply Skull-Stripping';
  ss.labels   = {'No','Yes'};
  ss.values   = {0 1};
  ss.val      = {1};
  ss.help     = {
    'Write skull-stripped images. '
    ''
  };

  nc          = cfg_menu;
  nc.tag      = 'nc';
  nc.name     = 'Noise Correction';
  nc.labels   = {'No','Yes'};
  nc.values   = {0 1};
  nc.val      = {1};
  nc.help     = {
    'Noise corrections of the final images.'
    ''
  };

  prefix         = cfg_entry;
  prefix.tag     = 'prefix';
  prefix.name    = 'Filename prefix';
  prefix.strtype = 's';
  prefix.num     = [0 Inf];
  prefix.val     = {'catsyn_'};
  prefix.help    = {
    'Prefix of output files.'};


  opts        = cfg_branch;
  opts.tag    = 'opts';
  opts.name   = 'Parameter';
  opts.val    = {spm bc in bvc ss nc prefix}; 
  opts.help   = {
    'Parameter settings for image correction.'
  };


  % --- Output ---
  pdo         = cfg_menu;
  pdo.tag     = 'pd';
  pdo.name    = 'PD Output';
  pdo.labels  = {'No','Yes'};
  pdo.values  = {0 1};
  pdo.val     = {1}; 
  pdo.help    = {
    'Write PD output images.'
  };

  t1o         = cfg_menu;
  t1o.tag     = 't1';
  t1o.name    = 'T1 Output';
  t1o.labels  = {'No','Yes'};
  t1o.values  = {0 1};
  t1o.val     = {1}; 
  t1o.help    = {
    'Write synthesized T1 output images based on the PD image.'
  };

  r1o         = cfg_menu;
  r1o.tag     = 'r1';
  r1o.name    = 'R1 Output';
  r1o.labels  = {'No','Yes'};
  r1o.values  = {0 1};
  r1o.val     = {1}; 
  r1o.help    = {
    'Write R1 output images.'
  };

  r2so        = cfg_menu;
  r2so.tag    = 'r2s';
  r2so.name   = 'R2s Output';
  r2so.labels = {'No','Yes'};
  r2so.values = {0 1};
  r2so.val    = {1}; 
  r2so.help   = {
    'Write R2s output images.'
  };

  bvco        = cfg_menu;
  bvco.tag    = 'bv';
  bvco.name   = 'Blood Vessel Output';
  bvco.labels = {'No','Yes'};
  bvco.values = {0 1};
  bvco.val    = {0}; 
  bvco.help   = {
    'Write map of blood vessels.'
  };
    
  output      = cfg_branch;
  output.tag  = 'output';
  output.name = 'Output';
  output.val  = {r1o r2so pdo t1o bvco}; 
  output.help = {
    'Output images.'
  };

  
  % --- main ---
  % batch mode - output is undefined!
  urqio       = cfg_exbranch;
  urqio.tag   = 'urqio';
  urqio.name  = 'Ultra-High Resolution Quantitative Image Optimization';
  urqio.val   = {data opts output};
  urqio.prog  = @cat_vol_urqio;
  %urqio.vout  = @vout_urqio;
  urqio.help  = {
    'Additional correction of high resolution PD, R1, and R2s weighted images that includes another bias correction, intensity normalization, and blood vessel correction step. '
    ''
    'WARNING: This tool is in development and was just tested on a small set of subjects!'
  };
function boxplot = conf_io_boxplot(outdir,subdir,name,expert)

  files                 = cfg_files;
  files.tag             = 'files';
  files.name            = 'Files';
  files.help            = {'files.' ''};
  files.filter          = 'any';
  files.ufilter         = '.*.xml';
  files.num             = [3 Inf];
  
  % - groupname     ... %  opt.names       = [];            % array of group names
  setname               = cfg_entry;
  setname.tag           = 'setname';
  setname.name          = 'Name';
  setname.help          = {'Name of the dataset that replaces the number of the set. ' ''}; 
  setname.strtype       = 's';
  setname.num           = [0 Inf];
  setname.val           = {''};
 
  % - groupcolor    ... %  has to be generated later
  setcolor              = cfg_entry;
  setcolor.tag          = 'setcolor';
  setcolor.name         = 'Color';
  setcolor.strtype      = 'r';
  setcolor.num          = [0 Inf];
  setcolor.help         = {
   ['Color of the dataset that replaces the default color definied by the color table below. ' ...
    'The color has to be defined as 3x1 RGB value between 0 and 1, e.g. [0 0.5 0] for dark green. '] ''}; 
  setcolor.val          = {''};
  
  % - groupcolor    ... %  has to be generated later
  color                 = cfg_menu;
  color.tag             = 'setcolor';
  color.name            = 'Color';
  color.labels          = { ...
    'colormap'
    'red (light)';    'red';    'red (dark)'; 
    'orange (light)'; 'orange'; 'orange (dark)'; 
    'yellow (light)'; 'yellow'; 'yellow (dark)'; 
    'green (light)';  'green';  'green (dark)'; 
    'cyan (light)';   'cyan';   'cyan (dark)'; 
    'blue (light)';   'blue';   'blue (dark)';
    'violet (light)'; 'violet'; 'violet (dark)';
    'gray (light)';   'gray';   'gray (dark)';
    };
  color.values          = {
    ''; 
    [1   1/3 1/3]; [1   0   0  ]; [2/3 0   0  ]; % red
    [1   2/3 1/3]; [1   1/2 0  ]; [2/3 1/3 0  ]; % orange
    [1   1   1/3]; [1   1   0  ]; [2/3 2/3 0  ]; % yellow
    [1/3 1   1/3]; [0   1   0  ]; [0   2/3 0  ]; % green
    [1/3 1   1  ]; [0   1   1  ]; [0   2/3 2/3]; % cyan
    [1/3 1/3 1  ]; [0   0   1  ]; [0   0   2/3]; % blue
    [1   1/3 1  ]; [1   0   1  ]; [2/3 0   2/3]; % violet
    [3/4 3/4 3/4]; [1/2 1/2 1/2]; [1/4 1/4 1/4]; % gray
    };
  color.val             = {''};
  color.help            = {
    'Color of the dataset that replaces the default color definied by the colormap below. ' ''}; 
  
  % subset .. was not working
  %{
  subset                = cfg_menu;
  subset.tag            = 'subset';
  subset.name           = 'Subset';
  subset.labels         = {'W','G'};
  subset.values         = {0 1};
  subset.def            = @(val) 0;
  subset.help           = {'Subset G with gray background' ''};
  %}
  
  % as structure with subfields
  datasetxml            = cfg_exbranch;
  datasetxml.tag        = 'data';
  datasetxml.name       = 'Dataset';
  datasetxml.val        = { files , setname, color};
  datasetxml.help       = {'Specify major properties of a dataset with predfined colors.' ''};
  
  datasetxml2            = cfg_exbranch;
  datasetxml2.tag        = 'data';
  datasetxml2.name       = 'Dataset (color definition)';
  datasetxml2.val        = { files , setname, setcolor}; 
  datasetxml2.help       = {'Specify major properties of a dataset with own color definition.' ''};
  
  datasets              = cfg_repeat;
  datasets.tag          = 'data';
  datasets.name         = 'Datasets';
  datasets.values       = {datasetxml datasetxml2};
  datasets.val          = {};
  datasets.num          = [1 Inf];
  datasets.help         = {'Specify manually grouped XML files their name and color. '};

  
  % or as xmlparagroup where all xml-files are selected and then internally differentiated 
  % - xmlfiles 
  % - XML grouping parameter(s) (eg. software.version, parameter.extopts.collcorr ) ... 
  % ... the idea is nice but what if multiple parameter change? 
  % ... it is more easy, clearer and saver to force manual grouping or may focus on some elements 
  %     Computer + SPM + CAT revision
  
  
  
  % ---
  
  % quality measures (expert)
  QMfield               = cfg_menu;
  QMfield.tag           = 'xmlfields';
  QMfield.name          = 'Image quality';
  QMfield.labels        = {
    'Noise Contrast Ratio (NCR)'
    'Inhomogeny Contrast Ratio (ICR)'
    'Resolution RMSE (resRMS)'
    'Minimum tissue contrast'
    };
  QMfield.values        = {
    'qualitymeasures.NCR'
    'qualitymeasures.ICR'
    'qualitymeasures.res_RMS'
    'qualitymeasures.contrast'
    };
  QMfield.val           = {'qualitymeasures.NCR'};
  QMfield.help          = {'CAT preprocessing image quality measures (not normalized).' ''};
  
  
  % quality ratings 
  QRfield               = cfg_menu;
  QRfield.tag           = 'xmlfields';
  QRfield.name          = 'Image quality ratings';
  QRfield.labels        = {
    'Noise Contrast Ratio (NCR)'
    'Inhomogeny Contrast Ratio (ICR)'
    'Resolution RMSE (resRMS)'
    'Minimum tissue contrast'
    };
  QRfield.values        = {
    'qualitratings.NCR'
    'qualitratings.ICR'
    'qualitratings.res_RMS'
    'qualitratings.contrast'
    };
  QRfield.val           = {'qualitratings.NCR'};
  QRfield.help          = {'CAT preprocessing image quality ratings (normalized marks).' ''};
  
  
  % surface measures
  SMfield               = cfg_menu;
  SMfield.tag           = 'xmlfields';
  SMfield.name          = 'Surface quality';
  SMfield.labels        = {
    ... 'Surface Euler number'
    'Surface defect area'
    'Surface defect number'
    'Surface intensity RMSE'
    'Surface position RMSE'
    'Surface self-intersections'
    };
  SMfield.values        = {
    ... 'qualitymeasures.SurfaceEulerNumber'
    'qualitymeasures.SurfaceDefectArea'
    'qualitymeasures.SurfaceDefectNumber'
    'qualitymeasures.SurfaceIntensityRMSE'
    'qualitymeasures.SurfacePositionRMSE'
    'qualitymeasures.SurfaceSelfIntersections'
    };
  SMfield.val           = {'qualitymeasures.SurfaceDefectArea'};
  SMfield.help          = {'CAT preprocessing surface quality measures (not normalized). ' ''};
  
  
  % segmenation measures
  USMfield               = cfg_menu;
  USMfield.tag           = 'xmlfields';
  USMfield.name          = 'Unified segmentation validation measures';
  USMfield.labels        = {
    'SPM log-likelyhood'
    'SPM tissue peak 1 (def. GM)'
    'SPM tissue peak 2 (def. WM)'
    'SPM tissue peak 3 (def. CSF1)'
    'SPM tissue peak 4 (def. CSF2)'
    'SPM tissue volume 1 (GM)'
    'SPM tissue volume 2 (WM)'
    'SPM tissue volume 3 (CSF)'
    'SPM tissue volume 4 (HD1)'
    'SPM tissue volume 5 (HD2)'
    'SPM tissue volume 6 (BG)'
    ...'CAT skull-stripping parameter'
    ...'CAT high BG parameter'
    };
  USMfield.values        = {
    'SPMpreprocessing.ll'
    'SPMpreprocessing.mn(1)'
    'SPMpreprocessing.mn(2)'
    'SPMpreprocessing.mn(3)'
    'SPMpreprocessing.mn(4)'
    'ppe.SPMvols0(1)'
    'ppe.SPMvols0(2)'
    'ppe.SPMvols0(3)'
    'ppe.SPMvols0(4)'
    'ppe.SPMvols0(5)'
    'ppe.SPMvols0(6)'
    ...'ppe.skullstrippedpara'
    ...'ppe.highBGpara'
    ...reg.ll
    ...reg.dt, rmsdt
    };
  USMfield.val           = {'SPMpreprocessing.ll'};
  USMfield.help          = {'SPM preprocessing measures for evaluation of the preprocessing. The tissue peaks depend on the defined number of SPM peaks within a class (default=[1 1 2 3 4 2]). The volumes depend on the TPM that are by default GM, MW, CSF, HD1 (hard tissue), HD2 (soft tissue), backgroun (BG). ' ''};
  USMfield.hidden        = expert<2; 
  
  
  % individual measures
  IMfield               = cfg_menu;
  IMfield.tag           = 'xmlfields';
  IMfield.name          = 'Morphometric measures';
  IMfield.labels        = {
    'Total Intracranial Volume (TIV)'
    'Total Surface Area (TSA)'
    'Mean cortical thickness'
    'Cortical thickness standard deviation'
    'Relative CSF volume'
    'Relative GM  volume'
    'Relative WM  volume'
    'Relative WMH volume'
    'Absolute CSF volume'
    'Absolute GM  volume'
    'Absolute WM  volume'
    'Absolute WMH volume'
    };
  IMfield.values        = {
    'subjectmeasures.vol_TIV'
    'subjectmeasures.surf_TSA'
    'subjectmeasures.dist_thickness{1}(1)'
    'subjectmeasures.dist_thickness{1}(2)'
    'subjectmeasures.vol_rel_CGW(1)'
    'subjectmeasures.vol_rel_CGW(2)'
    'subjectmeasures.vol_rel_CGW(3)'
    'subjectmeasures.vol_rel_CGW(4)'
    'subjectmeasures.vol_abs_CGW(1)'
    'subjectmeasures.vol_abs_CGW(2)'
    'subjectmeasures.vol_abs_CGW(3)'
    'subjectmeasures.vol_abs_CGW(4)'
    ...'ppe.reg.rmsdt'
    ...'ppe.reg.rmsdtc'
    };
  IMfield.val           = {'subjectmeasures.vol_TIV'};
  IMfield.help          = {'Global morphometric measures. ' ''};
  
  
  % - title 
  ftitle              = setname; 
  ftitle.name         = 'title';
  ftitle.tag          = 'Plot title';
  ftitle.help         = {'Name of figure' ''};
  ftitle.val          = {''};
  % - yname (measure/scala)
  fname             = setname; 
  fname.name          = 'name';
  fname.tag           = 'Measure name';
  fname.help          = {'Name of the measure ploted at the y-axis. ' ''};
  %
  fspec               = setname; 
  fspec.name          = 'name';
  fspec.tag           = 'Measure name';
  fspec.help          = {'Name of the measure ploted at the y-axis. ' ''};
  %  opt.ylim        = [-inf inf];    % y-axis scaling
  ylim                = cfg_entry;
  ylim.tag            = 'ylim';
  ylim.name           = 'y-axis limits';
  ylim.help           = {'Limitation of x-axis. '}; 
  ylim.strtype        = 'r';
  ylim.num            = [1 2];
  ylim.val            = {[-inf inf]}; 
  %  opt.subsets     = false(1,numel(data)); 
  
  xmlfield0           = cfg_exbranch;
  xmlfield0.tag       = 'xmlfields';
  xmlfield0.name      = 'Data field (complex)';
  xmlfield0.val       = { ftitle  , fname , fspec , ylim}; 
  xmlfield0.help      = {'Specify set properties such as name or color' ''};
  
  
  xmlfield            = cfg_entry;
  xmlfield.tag        = 'xmlfields';
  xmlfield.name       = 'Data field (simple)';
  xmlfield.help       = { 
   ['Specify field for data extraction that result in one value per file, e.g., ' ...
    'measures.vol_rel_CGW(1) to extract the first (CSF) volume value. '] ''};
  xmlfield.strtype    = 's';
  xmlfield.num        = [1 Inf];
  xmlfield.def        = @(val) 'subjectmeasures.vol_TIV';
  
  xmlfields           = cfg_repeat;
  xmlfields.tag       = 'xmlfields';
  xmlfields.name      = 'XML-fields';
  if expert
    xmlfields.values  = {xmlfield,xmlfield0,USMfield,QMfield,QRfield,SMfield,IMfield};
  else
    xmlfields.values  = {xmlfield,xmlfield0,QRfield,SMfield,IMfield};
  end
  xmlfields.val       = {};
  xmlfields.num       = [1 Inf];
  xmlfields.forcestruct;
  xmlfields.help      = {'Specify manually grouped XML files.'};
  
  % ------
  
  
  % main figure title and xlabel (ylabel is defined by the fieldselector)
  title               = cfg_entry;
  title.tag           = 'title';
  title.name          = 'Title';
  title.help          = {'Title of the plot append to the field specific title definition.'
    'If the first char is a + than the title is added to automatic generated titles. ' 
    ''}; 
  title.strtype       = 's';
  title.num           = [0 inf];
  title.val           = {''};
  
  xlabel              = name;
  xlabel.tag          = 'xlabel';
  xlabel.name         = 'Xlabel';
  xlabel.val          = {'groups'};
  xlabel.help         = {'General group name, e.g., methods, versions, parameter A. ' ''};
  
            
  %  opt.style       = 0;             % violin-plot: 0 - box plot; 1 - violin plot; 2 - violin + thin box plot 
  %  opt.violin      = 0;             % violin-plot: 0 - box plot; 1 - violin plot; 2 - violin + thin box plot 
  style               = cfg_menu;
  style.tag           = 'style';
  style.name          = 'Plotting Style';
  style.labels        = {'Box-plot','Violin-lot','Violin-box-lot','Density-plot'};
  style.values        = {0 1 2 3};
  style.def           = @(val) 0;
  style.help          = {'Type of data plot.' ''};
  
  % colorset    - menu? jet, hsv, ...
  colorset            = cfg_menu;
  colorset.tag        = 'colorset';
  colorset.name       = 'Colorset';
  colorset.labels     = {'Parula','Jet','HSV','Hot','Cool','Spring','Summer','Autumn','Winter','Lines','Prism'};
  colorset.values     = {'parula','jet','hsv','hot','cool','spring','summer','autumn','winter','lines','prism'};
  colorset.def        = @(val) 'jet';
  colorset.help       = {'Default colors of the boxes.' ''};
  
  
  % figformat   - 
  fsize               = cfg_entry;
  fsize.tag           = 'fsize';
  fsize.name          = 'Figure size in cm';
  fsize.help          = {
    'Define height/size of the figure. Using only one value only defines the heigh. '
    'Default matlab figure is [9.8778 7.4083] cm. ' 
    }; 
  fsize.strtype       = 'r';
  fsize.num           = [1 2];
  fsize.val           = {[4.5 3.6]}; 
  
  %  opt.notched     = 0;             % thinner at median [0 1] with 1=0.5
  notched             = cfg_menu;
  notched.tag         = 'notched';
  notched.name        = 'Notched boxes';
  notched.labels      = {'No','Yes'};
  notched.values      = {0 1};
  notched.def         = @(val) 0;
  notched.help        = {'Notched boxes that are thinner at the median.' ''};

  %  opt.symbol      = '+o';          % outlier symbols
  symbol              = cfg_entry;
  symbol.tag          = 'symbol';
  symbol.name         = 'Outlier symbols';
  symbol.help         = {'Notched boxes that are thinner at the median.' ''};
  symbol.strtype      = 's';
  symbol.num          = [1 2];
  symbol.val          = {'+o'}; 
  
  %  opt.maxwhisker  = 1.5;           % 
  maxwhisker          = cfg_entry;
  maxwhisker.tag      = 'maxwhisker';
  maxwhisker.name     = 'Maximum whisker';
  maxwhisker.help     = {' '}; 
  maxwhisker.strtype  = 'r';
  maxwhisker.num      = [1 1];
  maxwhisker.val      = {1.5}; 
  
%  opt.sort        = 0;             % no sorting
%                  = 1;             % sort groups (ascending)
%                  = 2;             % sort groups (descending)[inactive]
%                  = [index];       % or by a index matrix
  sortdata            = cfg_menu;
  sortdata.tag        = 'sort';
  sortdata.name       = 'Sort groups';
  sortdata.labels     = {'No','Ascending','Descending'};
  sortdata.values     = {0 1 2};
  sortdata.def        = @(val) 0;
  sortdata.help       = {'' ''};

%  opt.fill        = 1;             % filling of boxes: 0 - no filling; 0.5 - half-filled boxes; 1 - filled boxes
  fill                = cfg_menu;
  fill.tag            = 'fill';
  fill.name           = 'Plot box';
  fill.labels         = {'No filling','half-filling','full-filling'};
  fill.values         = {0 0.5 1};
  fill.def            = @(val) 1;
  fill.help           = {'Filling of boxes: 0 - no filling; 0.5 - half-filled boxes; 1 - filled boxes' ''};
  fill.hidden         = expert<1; 
  
  menubar             = cfg_menu;
  menubar.tag         = 'menubar';
  menubar.name        = 'Display MATLAB menu bars';
  menubar.labels      = {'No','Yes'};
  menubar.values      = {0 1};
  menubar.def         = @(val) 0;
  menubar.help        = {'Print number of elements of each group in the groupname' ''};

  close               = cfg_menu;
  close.tag           = 'close';
  close.name          = 'Close figures';
  close.labels        = {'No','Yes'};
  close.values        = {0 1};
  close.def           = @(val) 0;
  close.help          = {'Close figures after export. ' ''};

  
%  opt.groupnum    = 1;             % add number of elements
  groupnum            = cfg_menu;
  groupnum.tag        = 'fill';
  groupnum.name       = 'Number group elements';
  groupnum.labels     = {'No','Yes'};
  groupnum.values     = {0 1};
  groupnum.def        = @(val) 1;
  groupnum.help       = {'Print number of elements of each group in the groupname' ''};
  groupnum.hidden     = expert<1; 

% [opt.groupmin    = 5;]            % minimum number of non-nan-elements in a group [inactive]

%  opt.ygrid       = 1;             % activate y-grid-lines
  ygrid               = cfg_menu;
  ygrid.tag           = 'ygrid';
  ygrid.name          = 'Plot y-grid lines';
  ygrid.labels        = {'No','Yes'};
  ygrid.values        = {0 1};
  ygrid.def           = @(val) 1;
  ygrid.help          = {'Plot y-grid-lines' ''};
  ygrid.hidden        = expert<1; 
  
%  opt.gridline    = '-';           % grid line-style


%  opt.box         = 1;             % plot box
  box                 = cfg_menu;
  box.tag             = 'box';
  box.name            = 'Plot box';
  box.labels          = {'No','Yes'};
  box.values          = {0 1};
  box.def             = @(val) 1;
  box.help            = {'' ''};

%  opt.outliers    = 1;             % plot outliers
  outliers            = cfg_menu;
  outliers.tag        = 'outliers';
  outliers.name       = 'Plot outliers';
  outliers.labels     = {'No','Yes'};
  outliers.values     = {0 1};
  outliers.def        = @(val) 1;
  outliers.help       = {'' ''};
  outliers.hidden     = expert<1; 
  
%  opt.boxwidth    = 0.8;           % width of box
  boxwidth            = cfg_entry;
  boxwidth.tag        = 'boxwidth';
  boxwidth.name       = 'Boxwidth'; 
  boxwidth.help       = {'Main width of the boxes. ' ''}; 
  boxwidth.strtype    = 'r';
  boxwidth.num        = [1 1];
  boxwidth.val        = {0.8}; 
  boxwidth.hidden     = expert<1; 
  
%  opt.groupcolor  = [R G B];       % matrix with (group)-bar-color(s) 
%                                     use jet(numel(data)) 
%                                     or other color functions
  groupcolormap               = cfg_menu;
  groupcolormap.tag           = 'groupcolormap';
  groupcolormap.name          = 'Colormap';
  groupcolormap.labels        = {'jet','hsv','warm','cold'};
  groupcolormap.values        = {'jet','hsv','warm','cold'};
  groupcolormap.def           = @(val) 1;
  groupcolormap.help          = {'' ''};

  % colormapdef
  colormapdef            = cfg_entry;
  colormapdef.tag        = 'colormapdef';
  colormapdef.name       = 'Own colormap'; 
  colormapdef.help       = {'Definition of a own colormap.' ''}; 
  colormapdef.strtype    = 'r';
  colormapdef.num        = [inf 3];
  colormapdef.val        = {jet(10)};
  
%  opt.symbolcolor = 'r';           % color of symbols
  symbolcolor               = cfg_menu;
  symbolcolor.tag           = 'symbolcolor';
  symbolcolor.name          = 'Outlier color';
  symbolcolor.labels        = {'red','green','blue','black','magenta','cyan','yellow'};
  symbolcolor.values        = {'r','g','b','n','m','c','y'};
  symbolcolor.val           = {'r'};
  symbolcolor.help          = {'Color of outlier symbols.' ''};


%  opt.showdata    = 0;             % show data points: 0 - no; 1 - as points; 2 - as short lines (barcode plot)
  showdata               = cfg_menu;
  showdata.tag           = 'showdata';
  showdata.name          = 'Show datapoints';
  showdata.labels        = {'No','Points','Bars'};
  showdata.values        = {0 1 2};
  showdata.def           = @(val) 0;
  showdata.help          = {'Show data points as points or as short lines (barcode plot).' ''};

%  opt.median      = 2;             % show median: 0 - no; 1 - line; 2 - with different fill colors 
  median               = cfg_menu;
  median.tag           = 'median';
  median.name          = 'Median style';
  median.labels        = {'No','Line','Brightend'};
  median.values        = {0 1 2};
  median.def           = @(val) 2;
  median.help          = {'Show median as line or with different fill colors. ' ''};
  
%  opt.edgecolor   = 'none';        % edge color of box 
  boxedgecolor         = cfg_menu;
  boxedgecolor.tag     = 'edgecolor';
  boxedgecolor.name    = 'Use edgecolor for boxes';
  boxedgecolor.labels  = {'No','Yes'};
  boxedgecolor.values  = {'none','-n'};
  boxedgecolor.val     = {'none'};
  boxedgecolor.help    = {'Show median as line or with different fill colors. ' ''};

%  opt.trans       = 0.25;          % transparency of the boxes

%  opt.sat         = 0.50;          % saturation of the boxes

  % opt.hflip = 0; 
  hflip               = cfg_menu;
  hflip.tag           = 'hflip';
  hflip.name          = 'Flip data';
  hflip.labels        = {'No','Yes'};
  hflip.values        = {0 1};
  hflip.def           = @(val) 0;
  hflip.help          = {'Flip data orientation.' ''};
  
  % opt.vertical    = 1;  % boxplot orientation 
  vertical            = cfg_menu;
  vertical.tag        = 'vertical';
  vertical.name       = 'Box orientation';
  vertical.labels     = {'Horizontal','Vertical'};
  vertical.values     = {0 1};
  vertical.def        = @(val) 1;
  vertical.help       = {'Boxplot orientation. ' ''};

  % fontsize
  fontsize            = cfg_entry;
  fontsize.tag        = 'FS';
  fontsize.name       = 'Fontsize';
  fontsize.help       = {'Main font size of the figure. ' ''}; 
  fontsize.strtype    = 'r';
  fontsize.num        = [1 1];
  fontsize.val        = {10}; 
  
  % main parameter structure passed to cat_plot_boxplot
  opts           = cfg_exbranch;
  opts.tag       = 'opts';
  opts.name      = 'Options';
  opts.val       = { title, xlabel , ...
    ygrid, style , colorset , hflip , vertical , fsize, fontsize}; 
  % boxedgecolor
  opts.help      = {'Specify the thickness of specific ROIs.' ''};
 
  extopts          = opts;
  extopts.tag      = 'extopts';
  extopts.name     = 'Extended options';
  extopts.val      = {
    showdata, sortdata, ...
    maxwhisker, symbol, symbolcolor, ...
    median, notched, boxwidth , boxedgecolor, fill, ...
    menubar, ...
  };
  
%     * figure size
%     * plot table 
% E
  
  % ------
  type                = cfg_menu;
  type.tag            = 'type';
  type.name           = 'Export data type';
  type.labels         = {'none','fig','png','epsc','fig + png','fig + png + epsc'};
  type.values         = {'','fig','png','epsc','fig png','fig png epsc'};
  type.def            = @(val) 'fig png epsc';
  type.help           = {'Type of output files.' ''};

  name.tag            = 'name';
  name.name           = 'Prefix';
  name.help           = {'Additional prefix'}; 
  name.val            = {''};
  
  subdir.val          = {'CAT_boxplot'};

  % dpi
  
  output              = cfg_exbranch;
  output.tag          = 'output';
  output.name         = 'Output'; 
  output.val          = {outdir,subdir,name,type,close}; 
  output.help         = {''};
  
  % -----
  boxplot             = cfg_exbranch;
  boxplot.tag         = 'boxplot';
  boxplot.name        = 'XML boxplot';
  if expert
    boxplot.val       = {datasets,xmlfields,opts,extopts,output};
  else
    boxplot.val       = {datasets,xmlfields,opts,output};
  end
  boxplot.prog        = @cat_plot_boxplot;
  boxplot.hidden      = expert<1; 
  %boxplot.vout        = @vout_io_boxplot;
  boxplot.help        = {''};
  
return
 
%_______________________________________________________________________
function avg_img = conf_vol_average(data,outdir)
% image average
% -------------------------------------------------------------------------

  % update input functions
  data.name       = 'Select images';
  data.help       = {'Select images for calculating average.'};

  outdir.help     = {'Select a directory where files are written otherwise the path of the first image will be used.'};

  % filename
  output          = cfg_entry;
  output.tag      = 'output';
  output.name     = 'Output Filename';
  output.help     = {
    'The output image is written to current working directory unless a valid full pathname is given. If a path name is given here, the output directory setting will be ignored.'
    'If the field is left empty, i.e. set to '''', then the name of the 1st input image, preprended with ''i'', is used (change this letter in the spm_defaults if necessary).'
  };
  output.strtype  = 's';
  output.num      = [0 Inf];
  output.val      = {'avg.nii'};

  % main
  avg_img         = cfg_exbranch;
  avg_img.tag     = 'avg_img';
  avg_img.name    = 'Image Average';
  avg_img.val     = {data output outdir};
  avg_img.help    = {'This function is for calculating the average of a set of images, which should be of same dimension and voxel size (i.e. after spatial registration).'};
  avg_img.prog    = @cat_vol_avg;
  avg_img.vout    = @vout_avg;

%_______________________________________________________________________
function data2mat = conf_io_data2mat(data,outdir)
% -------------------------------------------------------------------------
% Batch to save Matlab mat files of surface or resampled volume data for use
% with machine learning tools

  resolution         = cfg_entry;
  resolution.tag     = 'resolution';
  resolution.name    = 'Spatial resolution for resampling';
  resolution.strtype = 'r';
  resolution.num     = [1 1];
  resolution.val     = {4};
  resolution.help    = {
    'Volume data can be saved with a lower spation resolution which is especially helpful with further use with machine learning tools such as relevance/support vector approaches or Gaussian Process models. Spatial structure of the data is not considered.'
    'Recommended resampling values are 3-8mm. For BrainAGE we obtained bes prediction accuracy with values of 4 or 8mm.'
  };

  c                  = cfg_entry;
  c.tag              = 'c';
  c.name             = 'Vector/Matrix';
  c.help             = {'Vector or matrix of nuisance values'};
  c.strtype          = 'r';
  c.num              = [Inf Inf];

  nuisance           = cfg_repeat;
  nuisance.tag       = 'nuisance';
  nuisance.name      = 'Nuisance variable';
  nuisance.values    = {c};
  nuisance.num       = [0 Inf];
  nuisance.help      = {'This option allows for the specification of nuisance effects to be removed from the data. A potential nuisance parameter can be TIV if you check segmented data with the default modulation. In this case the variance explained by TIV will be removed from the data. Another meaningful nuisance effect is age. This parameter should be defined for all samples as one variable and may also contain several columns.'};

  mask = data; 
  mask.tag           = 'mask';
  mask.name          = 'Select brain mask image';
  mask.def           = @(val) cat_get_defaults('extopts.brainmask', val{:});
  mask.help          = {'Select additionally mask image to exclude non-brain areas.';''};
  mask.num           = [0 1];

  data.name          = 'Sample volume data';
  data.tag           = 'data';
  data.filter        = 'image';
  data.num           = [1 Inf];
  data.help          = {'Select spatially registered data. They must all have the same image dimensions, orientation, voxel size etc. Furthermore, it is recommended to use smoothed files with further use with machine learning tools.'};

  sample             = cfg_repeat;
  sample.tag         = 'sample';
  sample.name        = 'Data';
  sample.values      = {data};
  sample.num         = [1 Inf];
  sample.help        = {'Specify data for each sample. If you specify different samples a label variable will be also saved that decodes the samples.'};

  vol_data           = cfg_exbranch;
  vol_data.tag       = 'vol_data';
  vol_data.name      = 'Volume data';
  vol_data.val       = {sample,mask,resolution};
  vol_data.help      = {''};
  
  data.name          = 'Sample surface data';
  data.filter        = 'mesh';
  data.ufilter       = 'resampled';
  data.help          = {'Select resampled and smoothed surface data. They must all have the same mesh size (32k or 164k).'};

  sample.values      = {data};
  surf_data          = cfg_exbranch;
  surf_data.tag      = 'surf_data';
  surf_data.name     = 'Surface data';
  surf_data.val      = {sample};
  surf_data.help     = {''};
  
  data_type          = cfg_choice;
  data_type.tag      = 'data_type';
  data_type.name     = 'Select data type';
  data_type.values   = {vol_data surf_data};
  data_type.val      = {vol_data};
  data_type.help     = {'Choose between volume and surface data.'};
  
  fname              = cfg_entry; 
  fname.name         = 'Filename';
  fname.tag          = 'fname';
  fname.val          = {'Data.mat'}; 
  fname.help         = {'Filename to save data matrix.'};

  data2mat           = cfg_exbranch;
  data2mat.tag       = 'data2mat';
  data2mat.name      = 'Save volume or surface data as mat-file';
  data2mat.val       = {data_type,nuisance,fname,outdir};
  data2mat.prog      = @cat_io_data2mat;
  data2mat.vout      = @vout_io_data2mat;
  data2mat.help      = {
    'Save spatially registered volume or resampled surface data as Matlab data matrix for further use with machine learning tools.'
    'Volume data can be optionally masked to remove non-brain areas.'
    'A mat-file will be saved with the following parameters:'
    '  Y     - data matrix with size number of subjects x number of voxels/vertices'
    '  label - label of samples'
    '  ind   - index for volume or surface data inside mask'
    '  dim   - dimension of original data'
  };

%_______________________________________________________________________
function vf = vout_defs(job)

PU = job.field1;
PI = job.images;

vf = cell(numel(PI),1);
for i=1:numel(PU), % ########### RD202201: i is not used  
    for m=1:numel(PI),
        [pth,nam,ext,num] = spm_fileparts(PI{m});
        
        switch job.modulate
        case 2
            filename = fullfile(pth,['m0w' nam ext num]);
        case 1
            filename = fullfile(pth,['mw' nam ext num]);
        case 0
            filename = fullfile(pth,['w' nam ext num]);
        otherwise 
            error('incorrect - DEP')
        end;
        vf{m} = filename;
    end
end

return;
%_______________________________________________________________________
function vf = vout_defs2(job)

  PU = job.field;
  PI = job.images;

  vf = cell(numel(PU),numel(PI));
  for i=1:numel(PU),
      for m=1:numel(PI),
          [pth,nam,ext,num] = spm_fileparts(PI{m}{i});

          switch job.modulate
          case 2
              filename = fullfile(pth,['m0w' nam ext num]);
          case 1
              filename = fullfile(pth,['mw' nam ext num]);
          case 0
              filename = fullfile(pth,['w' nam ext num]);
          otherwise 
            error('incorrect - DEP')
          end;
          vf{i,m} = filename;
      end
  end

return;
%_______________________________________________________________________
function cdep = vout_urqio(job)
  %%
  cdep = cfg_dep;
  if job.output.r1
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'R1 Images';
    cdep(end).src_output = substruct('.','data','()',{1},'.',[job.opts.prefix 'r1_'],'()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
  if job.output.pd
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'PD Images';
    cdep(end).src_output = substruct('.','data','()',{1},'.',[job.opts.prefix 'pd_'],'()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
  if job.output.t1
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'T1 Images';
    cdep(end).src_output = substruct('.','data','()',{1},'.',[job.opts.prefix 't1_'],'()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
  if job.output.r2s==1 || job.output.r2s==3
    cdep(end+1)           = cfg_dep;
    cdep(end).sname      = 'R2s nobc Images';
    cdep(end).src_output = substruct('.','data','()',{1},'.',[job.opts.prefix 'nobc_r2s_'],'()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end 
  if job.output.r2s==2 || job.output.r2s==3
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'R2s bc Images';
    cdep(end).src_output = substruct('.','data','()',{1},'.',[job.opts.prefix 'bc_r2s_'],'()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end 
  if job.output.bv
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Blood Vessel Images';
    cdep(end).src_output = substruct('.','data','()',{1},'.',[job.opts.prefix 'bv_'],'()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
  if numel(cdep)>1
    cdep(1)=[];
  end
%%
return;
%_______________________________________________________________________
function dep = vout_sanlm(varargin)
  %job.returnOnlyFilename = 1; 
  %vf = cat_vol_sanlm(job); 
  
  dep(1)            = cfg_dep;
  dep(1).sname      = 'SANLM Images';
  dep(1).src_output = substruct('.','files');
  dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
return;
%_______________________________________________________________________
function dep = vout_maskimg(varargin)
  %job.returnOnlyFilename = 1; 
  %vf = cat_vol_maskimage(job); 
  
  dep            = cfg_dep;
  dep.sname      = 'Masked Images';
  dep.src_output = substruct('.','files');
  dep.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
return;
%_______________________________________________________________________
function cdep = vout_headtrimming(job)
  job.returnOnlyFilename = 1; 
  %vf = cat_vol_headtrimming(job); 
  vf = job; 
  cdep = cfg_dep;

  if isfield(vf.image_selector,'manysubjects')
    cdep(end).sname      = 'source images';
    cdep(end).src_output = substruct('.','image_selector','.','manysubjects','.','simages'); 
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    if isfield(vf.image_selector.manysubjects,'oimages')
      for i=1:numel(vf.image_selector.manysubjects.oimages)
        cdep(end+1)          = cfg_dep;%#ok<AGROW>
        cdep(end).sname      = sprintf('other images %d',i);
        cdep(end).src_output = substruct('.','image_selector','.','manysubjects','.','oimages','{}',{i}); 
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
      end
    end
  elseif isfield(vf.image_selector,'subjectimages') 
    % image-wise
    % - first image
    cdep(end).sname      = sprintf('first images of all subjects');
    cdep(end).src_output = substruct('.','image_selector','.','firstimages'); 
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    % - other images
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = sprintf('other images of all subjects');
    cdep(end).src_output = substruct('.','image_selector','.','otherimages'); 
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    % subject-wise ... 
    % the substruct seems to be correct, but it does not work, 
    % probably because each cdep entry has to be unique  
    % RD201810
    %{
    for si=1:numel(vf.image_selector.subjectimages)
      cdep(end+1)          = cfg_dep;%#ok<AGROW>
      cdep(end).sname      = sprintf('all imgages of subject %d',si);
      cdep(end).src_output = substruct('.','image_selector','.','subjectimages','{}',{si}); 
      cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    % single scans
    for si=1:numel(vf.image_selector.subjectimages)
      for fi=1:numel(vf.image_selector.subjectimages{si})
        cdep(end+1)          = cfg_dep;%#ok<AGROW>
        cdep(end).sname      = sprintf('subject %d image %d',si,fi);
        cdep(end).src_output = substruct('.','image_selector','.','subjectimages','{}',{si},'()',{fi}); 
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
      end
    end
    %}
  else
    for i=1:numel(vf.images)
      cdep(i)            = cfg_dep;
      cdep(i).sname      = sprintf('image %d',i);
      cdep(i).src_output = substruct('.','images','{}',{i}); 
      cdep(i).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
  end

return;

%_______________________________________________________________________
function dep = vout_volctype(varargin)
 % job.returnOnlyFilename = 1; 
 % vf = cat_io_volctype(job);
    
  dep            = cfg_dep;
  dep.sname      = 'Converted Images';
  dep.src_output = substruct('.','files','()',{':'});
  dep.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
return;

%_______________________________________________________________________
function dep = vout_createTPM(varargin)
  dep(1)              = cfg_dep;
  dep(1).sname        = 'TPM';
  dep(1).src_output   = substruct('.','tpm','()',{':'});
  dep(1).tgt_spec     = cfg_findspec({{'filter','image','strtype','e'}});
  
  dep(end+1)          = cfg_dep;
  dep(end).sname      = 'T1';
  dep(end).src_output = substruct('.','t1','()',{':'});
  dep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});

  dep(end+1)          = cfg_dep;
  dep(end).sname      = 'atlases';
  dep(end).src_output = substruct('.','atlas','()',{':'});
  dep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});

  dep(end+1)          = cfg_dep;
  dep(end).sname      = 'brainmask';
  dep(end).src_output = substruct('.','brainmask','()',{':'});
  dep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
return;

%_______________________________________________________________________
function dep = vout_createTPMlong(varargin)
  dep            = cfg_dep;
  dep.sname      = 'Longitudinal TPMs';
  dep.src_output = substruct('.','tpm','()',{':'});
  dep.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});

  dep            = cfg_dep;
  dep.sname      = 'Longitudinal TPMs Tissues';
  dep.src_output = substruct('.','tpmtiss','()',{':'},'()',{':'});
  dep.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});  
return;

%_______________________________________________________________________
function dep = vout_conf_longBiasCorr(varargin)
  dep            = cfg_dep;
  dep.sname      = 'Longitudinal Bias Corrected';
  dep.src_output = substruct('.','bc','()',{':'});
  dep.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
return

%_______________________________________________________________________
function dep = vout_resize(varargin)
  dep            = cfg_dep;
  dep.sname      = 'Resized';
  dep.src_output = substruct('.','res','()',{':'});
  dep.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
return;
%_______________________________________________________________________
function vf = vout_qa(job)
  if isfield(job,'data')
    s  = cellstr(char(job.data)); vf = s; 
  elseif isfield(job,'images')
    s  = cellstr(char(job.images)); vf = s; 
  else  
    s  = {}; 
    vf = {}; 
  end
  for i=1:numel(s)
      [pth,nam,ext,num] = spm_fileparts(s{i});
      if isfield(job,'prefix') % old 
        vf{i} = fullfile(pth,[job.prefix,nam,ext,num]);
      elseif isfield(job,'opts') && isfield(job.opts,'prefix') 
        vf{i} = fullfile(pth,[job.opts.prefix,nam,ext,num]);
      else
        vf{i} = fullfile(pth,[nam,ext,num]);
      end  
  end
return;

%------------------------------------------------------------------------
function dep = vout_avg(job)
  dep            = cfg_dep;
  if ~ischar(job.output) || strcmp(job.output, '<UNDEFINED>')
      dep.sname  = 'Average Image';
  else
      dep.sname  = sprintf('Average Image %s', job.output);
  end
  dep.src_output = substruct('.','files');
  dep.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
return

%------------------------------------------------------------------------
function dep = vout_stat_TIV(varargin)
  dep            = cfg_dep;
  dep.sname      = 'TIV';
  dep.src_output = substruct('.','calcvol','()',{':'});
  dep.tgt_spec   = cfg_findspec({{'strtype','e','strtype','r'}});
return

%_______________________________________________________________________
function dep = vout_io_data2mat(varargin)
    
  dep            = cfg_dep;
  dep.sname      = 'Saved mat-file';
  dep.src_output = substruct('.','fname','()',{':'});
  dep.tgt_spec   = cfg_findspec({{'strtype','e','strtype','r'}});
return;

%------------------------------------------------------------------------
function cdep = vout_realign(job)
  ind  = 1;
  if job.write_avg
      cdep(ind)            = cfg_dep;
      cdep(ind).sname      = 'Midpoint Average';
      cdep(ind).src_output = substruct('.','avg','()',{':'});
      cdep(ind).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
      ind = ind + 1;
  end
  if job.write_rimg
      cdep(ind)            = cfg_dep;
      cdep(ind).sname      = 'Realigned images';
      cdep(ind).src_output = substruct('.','rimg','()',{':'});
      cdep(ind).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
      ind = ind + 1;
  end
  if isfield(job.reg,'nonlin') && job.reg.nonlin.write_jac
      cdep(ind)            = cfg_dep;
      cdep(ind).sname      = 'Jacobian Diff';
      cdep(ind).src_output = substruct('.','jac','()',{':'});
      cdep(ind).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
      ind = ind + 1;
  end
  if isfield(job.reg,'nonlin') && job.reg.nonlin.write_def
      cdep(ind)            = cfg_dep;
      cdep(ind).sname      = 'Deformation (1)';
      cdep(ind).src_output = substruct('.','def1','()',{':'});
      cdep(ind).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
      ind = ind + 1;

      cdep(ind)            = cfg_dep;
      cdep(ind).sname      = 'Deformation (2)';
      cdep(ind).src_output = substruct('.','def2','()',{':'});
      cdep(ind).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end

return
 
%------------------------------------------------------------------------
function dep = vout_stat_spm2x(job)
  dep            = cfg_dep;
  dep.sname      = 'Transform & Threshold spm volumes';
  dep.src_output = substruct('.','Pname');
  dep.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
return

%------------------------------------------------------------------------
function dep = vout_stat_spm2x_surf(job)
  dep            = cfg_dep;
  dep.sname      = 'Transform & Threshold spm surfaces';
  dep.src_output = substruct('.','Pname');
  dep.tgt_spec   = cfg_findspec({{'filter','gifti','strtype','e'}});
return

%------------------------------------------------------------------------
function dep = vout_stat_getCSVXML(job)
  dep = cfg_dep;

  job.dep = 1; 
  if isempty( job.files )
    out = cat_stat_getCSVXMLfield(job); 

    FN = fieldnames(out);
    if iscell(job.fields)
      for fni = 1:numel(FN)
        dep(end + (fni>1))  = cfg_dep; 
        dep(end).sname      = sprintf('%s',FN{fni});
        dep(end).src_output = substruct('.',FN{fni},'()',{':'});
        dep(end).tgt_spec   = cfg_findspec({}); %{{'filter',,'strtype','e'}});
      end
    end
  end
return
function dep = vout_file_move(job)

% Define virtual output for cfg_run_move_file. Output can be passed on to
% either a cfg_files or an evaluated cfg_entry.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id$

rev = '$Rev$'; %#ok

if ~isfield(job.action,'delete')
    dep = cfg_dep;
    dep.sname = 'Moved/Rename/Copied Files';
    dep.src_output = substruct('.','files');
    dep.tgt_spec   = cfg_findspec({{'class','cfg_files','strtype','e'}});
else
    dep = [];
end
return
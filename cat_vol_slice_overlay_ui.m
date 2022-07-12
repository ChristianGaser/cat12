function cat_vol_slice_overlay_ui
% Wrapper to cat_vol_slice_overlay
% Call help for slice_overlay for any additional help
% 
% Additional fields to slice_overlay:
% OV.name       - char array of filenames for overlay that can be interactively
%                 selected
% OV.slices_str - char array of slice values (e.g. '-32:2:20')
%                 use empty string for automatically estimating slices with
%                 local maxima
% OV.xy         - define number of columns and rows
%                 comment this out for interactive selection or set the values
%                 to [Inf 1] for using one row and automatically estimate number
%                 of columns or use [1 Inf] for using one column
% OV.atlas      - define atlas for labeling
%                 comment this out for interactive selection
%                 or use 'none' for no atlas information
% OV.save       - save result as png/jpg/pdf/tif
%                 comment this out for interactive selection or use '' for not 
%                 saving any file or use just file extension (png/jpg/pdf/tif) to 
%                 automatically estimate filename to save
% OV.FS         - normalized font size (default 0.08)
% OV.name_subfolder
%               - if result is saved as image use up to 2 subfolders to add their 
%                 names to the filename (default 1)
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

% use default T1 from Shooting
%OV.reference_image = fullfile(spm('dir'),'toolbox','cat12','templates_MNI152NLin2009cAsym','Template_T1.nii');
% or its masked version
OV.reference_image = fullfile(spm('dir'),'toolbox','cat12','templates_MNI152NLin2009cAsym','Template_T1_masked.nii');

OV.reference_range = [0.2 1.0];                        % intensity range for reference image
OV.opacity = Inf;                                      % transparency value for overlay (<1)
OV.cmap    = jet;                                      % colormap for overlay

% char array of file names
OV.name = char(fullfile(cat_get_defaults('extopts.pth_templates'),'Template_4_GS.nii,1'),...
               fullfile(cat_get_defaults('extopts.pth_templates'),'cobra.nii'));
                
% range for each file
% Use range 0..0 if you want to autoscale range.
% If you are using log. scaling, check the highest p-value in the table
% and approximate the range; e.g. for a max. p-value of p=1e-7 and a
% threshold of p<0.001 use a range of [3 7]. Check cat_stat_spmT2x.m for details.
% If you are unsure, simply use the autorange option by using a range of [0 0].
% The log-scaled values are calculated by -log10(1-p):
% p-value       -log10(1-P)
%  0.1           1
%  0.05          1.30103 (-log10(0.05))
%  0.01          2
%  0.001         3
%  0.0001        4

% Number of fields in range should be the same as number of files (see above)
% or define one field, which is valid for all.
% Be careful: intensities below the lower range are not shown!
OV.range   =[[0.5 1]; [126 139]];

% OV.func can be used to set the image to defined values (e.g. NaN) for the given range
%OV.func = 'i1(i1>log10(0.05) & i1<-log10(0.05))=NaN;';

% selection of slices and orientations
% if OV.slices_str is an empty string then slices with local maxima are estimated automatically
OV.slices_str = char('','0:2:36','-40:5:-5');
OV.transform  = char('axial','sagittal','coronal');

% define output format of slices
OV.labels.format = '%3.1f';

% define number of columns and rows
% comment this out for interactive selection
%OV.xy = [3 5];

% or use Inf to automatically estimate the number of necessray rows or columns
OV.xy = [Inf 1]; % use one row and automatically estimate number of columns

% save result as png/jpg/pdf/tif
% comment this out for interactive selection or use 'none' for not 
% saving any file or use just file extension (png/jpg/pdf/tif) to automatically
% estimate filename to save
OV.save = 'png';

% if result is saved as image use up to 2 subfolders to add their names to the filename (default 1)
OV.name_subfolder = 2;

% Comment this out if you wish slice overview
OV.overview = [];

% Comment this out if you wish slice labels
OV.labels = [];

% Comment this out if you wish colorbar
OV.cbar = [];

% Normalized font size
OV.FS = 0.08;

% define atlas for labeling
% comment this out for interactive selection
% or use 'none' for skipping atlas information
OV.atlas = 'cat12_neuromorphometrics';

% call slice overlay with that settings
cat_vol_slice_overlay(OV)

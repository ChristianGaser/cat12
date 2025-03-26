% Computational Anatomy Toolbox
% Version 2664 (CAT12.9) 2025-03-10
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$
% ==========================================================================
% Description
% ==========================================================================
% This toolbox is a collection of extensions to the segmentation algorithm 
% of SPM12 (Wellcome Department of Cognitive Neurology) to provide computational
% morphometry. It is developed by Christian Gaser and Robert Dahnke (Jena
% University Hospital, Departments of Psychiatry and Neurology) and is available 
% to the scientific community under the terms of the GNU General Public License as 
% published by the Free Software Foundation; either version 2 of the License, 
% or (at your option) any later version.
%
% If you use any CAT12 code for commercial application, please email 
% christian.gaser@uni-jena.de.
%
% General files
%   INSTALL.txt                      - installation instructions
%   CHANGES.txt                      - changes in revisions
%   Contents.m                       - this file
%
% Core functions
%   cat12.m
%   cat_amap.m                       - compilation wrapper for cat_amap.c
%   cat_main.m
%   cat_run.m                        - runtime funtion for CAT12
%   cat_run_job.m
%   cat_run_newcatch.m
%   cat_run_oldcatch.m
%   cat_defaults.m                   - sets the defaults for CAT12
%   cat_get_defaults.m               - defaults for CAT12
%   spm_cat12.m                      - toolbox wrapper to call CAT12
%   
% Utilities
%   cat_sanlm.m                      - Spatial Adaptive Non Local Means Denoising Filter
%   cat_update.m                     - check for new updates
%   cat_debug.m                      - print debug information for SPM12 and CAT12
%   cat_ornlm.m
%   cat_plot_boxplot.m
%   slice_overlay.m                  - overlay tool
%   sliderPanel.m
%   
% Input & Output
%   cat_io_3Dto4D.m
%   cat_io_FreeSurfer.m
%   cat_io_cgw2seg.m
%   cat_io_checkinopt.m
%   cat_io_cmd.m
%   cat_io_colormaps.m
%   cat_io_cprintf.m
%   cat_io_csv.m
%   cat_io_handle_pre.m
%   cat_io_img2nii.m
%   cat_io_matlabversion.m
%   cat_io_remat.m
%   cat_io_seg2cgw.m
%   cat_io_struct2table.m
%   cat_io_updateStruct.m
%   cat_io_writenii.m
%   cat_io_xml.m
%   cat_io_xml2csv.m
%   
% Longitudinal batch
%   cat_long_main.m                  - longitudinal batch mode
%   cat_long_multi_run.m             - call cat_long_main for multiple subjects
%   
% Statistics
%   cat_stat_calc_stc.m
%   cat_stat_homogeneity.m           - check sample homogeneity across sample
%   cat_stat_marks.m
%   cat_stat_nanmean.m
%   cat_stat_nanmedian.m
%   cat_stat_nanstat1d.m
%   cat_stat_nanstd.m
%   cat_stat_nansum.m
%   cat_stat_showslice_all.m         - show 1 slice of all images
%   cat_stat_spm.m
%   cat_stat_spmF2x.m                - transformation of F-maps to P, -log(P), R2 maps
%   cat_stat_spmT2x.m                - transformation of t-maps to P, -log(P), r or d-maps
%   cat_stat_TIV.m                   - read total intracranial volume (TIV) from xml-files
%   
% Surface functions
%   cat_surf_avg.m
%   cat_surf_calc.m
%   cat_surf_createCS.m
%   cat_surf_display.m
%   cat_surf_info.m
%   cat_surf_parameters.m
%   cat_surf_rename.m
%   cat_surf_render.m
%   cat_surf_resamp.m
%   cat_surf_resamp_freesurfer.m
%   cat_surf_resample.m
%   cat_surf_smooth.m
%   cat_surf_vol2surf.m
%   
% Test and experimental functions
%   cat_tst_BWPsliceartifact.m
%   cat_tst_CJV.m
%   cat_tst_calc_kappa.m
%   cat_tst_qa.m
%   cat_tst_staple_multilabels.m
%   
% Volume functions
%   cat_vol_approx.m
%   cat_vol_atlas.m
%   cat_vol_average.m
%   cat_vol_calc_roi.m
%   cat_vol_correct_slice_scaling.m
%   cat_vol_ctype.m
%   cat_vol_defs.m                   - apply deformations to images
%   cat_vol_findfiles.m
%   cat_vol_groupwise_ls.m
%   cat_vol_imcalc.m
%   cat_vol_iscale.m
%   cat_vol_morph.m                  - morphological operations to 3D data
%   cat_vol_nanmean3.m
%   cat_vol_partvol.m
%   cat_vol_pbt.m
%   cat_vol_resize.m
%   cat_vol_sanlm.m                  - GUI for cat_sanlm
%   cat_vol_series_align.m
%   cat_vol_set_com.m
%   cat_vol_slice_overlay.m          - wrapper for overlay tool slice_overlay
%   cat_vol_slice_overlay_ui.m       - example for user interface for overlay wrapper cat_slice_overlay.m
%   cat_vol_smooth3X.m
%      
% Batch mode
%   cat_batch_long.m                 - batch mode wrapper for spm_jobman for longitudinal pipeline
%   cat_batch_spm.m                  - batch mode wrapper for spm_jobman for SPM12
%   cat_batch_vbm.m                  - batch mode wrapper for spm_jobman for CAT12
%   
% Check input and files
%   cat_check.m
%   cat_check_system_output.m
%   
% Configuration
%   cat_conf_extopts.m
%   cat_conf_long.m
%   cat_conf_opts.m
%   cat_conf_stools.m
%   cat_conf_stoolsexp.m
%   cat_conf_tools.m                 - wrapper for calling CAT12 utilities
%   tbx_cfg_cat.m
%
% Templates/Atlases volumes
%   Template_?_GS.nii                - Geodesic Shooting template of 555 subjects from IXI database
%                                      in MNI152NLin2009cAsym space provided for 6 different iteration steps
%   Template_T1.nii                  - average of 555 T1 images of IXI database in MNI152 space after 
%                                      Geodesic Shooting
%   aal3.*                           - AAL3 atlas
%   anatomy2.*                       - Anatomy2 atlas
%   cobra.*                          - CoBra atlas
%   hammers.*                        - Hammers atlas
%   ibsr.*                           - IBSR atlas
%   julichbrain.*                    _ JulichBrain atlas
%   mori.*                           - Mori atlas
%   thalamus.*                       - Thalamic Nuclei atlas
%   neuromorphometrics.*             - Neuromorphometrics atlas
%   lpba40.nii                       - LPBA40 atlas
%   brainmask.nii                    - brainmask
%   cat12.nii                        - partitions for hemispheres, subcortical structures and vessels
%
% Templates/Atlases surfaces
%   ?h.central.Template_T1.gii
%                                    - central surface of Dartel average brain for result mapping
%   ?h.central.freesurfer.gii        - central surface of freesurfer fsaverage
%   ?h.inflated.freesurfer.gii       - inflated surface of freesurfer fsaverage
%   ?h.sphere.freesurfer.gii         - spherical surface of freesurfer fsvaverage




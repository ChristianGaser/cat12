function tools = cg_vbm8_tools
% wrapper for calling VBM utilities
%
%_______________________________________________________________________
% Christian Gaser
% $Id$

rev = '$Rev$';

%_______________________________________________________________________

data = cfg_files;
data.tag  = 'data';
data.name = 'Volumes';
data.filter = 'image';
data.ufilter = '.*';
data.num     = [1 Inf];
data.help = {[...
'Select raw data (e.g. T1 images) for processing. ',...
'This assumes that there is one scan for each subject. ',...
'Note that multi-spectral (when there are two or more registered ',...
'images of different contrasts) processing is not yet implemented ',...
'for this method.']};

%------------------------------------------------------------------------

data_T2x = cfg_files;
data_T2x.tag  = 'data';
data_T2x.name = 'Volumes';
data_T2x.filter = 'image';
data_T2x.ufilter = '^spmT.*\.[in][im][gi]$';
data_T2x.num     = [1 Inf];
data_T2x.help = {'Select spmT-images to convert.'};

T2x = cfg_exbranch;
T2x.tag = 'T2x';
T2x.name = 'Threshold and transform spmT-maps';
T2x.val = {data_T2x};
T2x.prog   = @cg_spmT2x;

p0 = '';
p1 = 'This function transforms t-maps to P, -log(P), r or d-maps.';
p2 = 'The following formulas are used:';
p3 = '--------------------------------';
p4 = 'correlation coefficient:';
p5 = '          sign(t)';
p6 = 'r = ------------------';
p7 = '           df';
p8 = '    sqrt(------ + 1)';
p9 = '          t*t';
p10='effect-size';
p11='           2r';
p12='d = ----------------';
p13='    sqrt(1-sqr(r))';
p14='p-value';
p15='p = 1-spm_Tcdf';
p16='log p-value';
p17='-log10(1-P) = -log(1-spm_Tcdf)';
p18=['For the last case of log transformation this means that a p-value of p=0.99 (0.01) is ',...
'transformed to a value of 2.'];
p19='Examples:';
p20='p-value  -log10(1-P)';
p21='0.1      1';
p22='0.05     1.3';
p23='0.01     2';
p24='0.001    3';
p25='0.0001   4';
p26=['All maps can be thresholded using height and extent thresholds and you can also apply corrections ',...
'for multiple comparisons based on family-wise error (FWE) or false discovery rate (FDR). You can easily ',...
'threshold and/or transform a large number of spmT-maps using the same thresholds.'];
p27='Naming convention of the transformed files:';
p28='   Type_Contrast_Pheight_Pextent_K_Neg';
p29='   Type:      P    - p-value';
p30='              logP - log p-value';
p31='              R    - correlation coefficient';
p32='              D    - effect size';
p33='              T    - t-value';
p34='   Contrast:  name used in the contrast manager with replaced none valid';
p35='              strings';
p36='   Pheight:   p    - uncorrected p-value in % (p<0.05 will coded with "p5")';
p37='              pFWE - p-value with FWE correction in %';
p38='              pFDR - p-value with FDR correction in %';
p39='   Pextent:   pk    - uncorr. extent p-value in % (p<0.05 coded with "p5")';
p40='              pkFWE - extent p-value with FWE correction in %';
p41='   K:         extent threshold in voxels';
p42='   Neg:       image also shows thresholded inverse effects (e.g. neg. ';
p43='              values) ';

T2x.help = {p1,p0,p2,p0,p3,p4,p3,p0,p5,p6,p7,p8,p9,p0,p3,p10,p3,p11,p12,p13,p0,p3,...
	p14,p3,p15,p0,p3,p16,p3,p17,p0,p18,p0,p19,p20,p21,p22,p23,p24,p25,p0,p26,p0,p27,p28,p0,...
	p29,p30,p31,p32,p33,p0,p34,p35,p0,p36,p37,p38,p0,p39,p40,p0,p41,p0,p42,p43};
%------------------------------------------------------------------------

data.help = {[...
'Select all images. Images have to be in the same orientation with same voxel size and dimension ',...
'(e.g. normalized images)']};

scale = cfg_menu;
scale.tag = 'scale';
scale.name = 'Proportional scaling?';
scale.labels = {'no','yes'};
scale.values = {0 1};
scale.val = {0};
scale.help = {[...
'This option should be only used if image intensity is not scaled (e.g. T1 images) ',...
'or if images have to be scaled during statistical analysis (e.g. modulated images).']};

check_sd_sdname = cfg_entry;
check_sd_sdname.tag = 'sd_name';
check_sd_sdname.name = 'Output standard deviation file';
check_sd_sdname.strtype = 's';
check_sd_sdname.num = [1 Inf];
check_sd_sdname.val  = {'SD.nii'};
check_sd_sdname.help  = {[...
'The output file is written to current working directory ',...
'unless a valid full pathname is given. If you do not want to write this file ',...
'leave name empty.']};

check_sd_meanname = cfg_entry;
check_sd_meanname.tag = 'mean_name';
check_sd_meanname.name = 'Output mean file';
check_sd_meanname.strtype = 's';
check_sd_meanname.num = [1 Inf];
check_sd_meanname.val  = {'Mean.nii'};
check_sd_meanname.help  = {[...
'The output file is written to current working directory ',...
'unless a valid full pathname is given. If you do not want to write this file ',...
'leave name empty.']};

check_sd = cfg_exbranch;
check_sd.tag = 'check_sd';
check_sd.name = 'Check sample homogeneity across sample';
check_sd.val = {data,scale,check_sd_meanname,check_sd_sdname};
check_sd.prog   = @cg_check_sample_sd;
check_sd.help = {[...
'If you have a reasonable sample size artefacts are easily overseen. In order to identify images with poor image quality ',...
'or even artefacts you can use this function. Images have to be in the same orientation with same voxel size and dimension ',...
'(e.g. normalized images). The idea of this tool is to check the standard deviation across the sample.'],...
'',[...
'Standard deviation is caclulated by the sum of the squared distance of each image from the sample mean. Hence, the squared ',...
'distance of one image from the sample mean represents the amount to which this images deviates from the sample mean.'],...
'',[...
'The squared distance to mean is calculated for each image and plotted using a boxplot and the indicated filenames. The larger ',...
'the squared distance the more deviant is this image from the sample mean. In the Òsquared distance to mean plotÓ outliers from ',...
'the sample are usually isolated from the majority of images which are clustered around the sample mean. The squared distance to ',...
'mean is plotted at the y-axis and the x-axis reflects the image order. Images are plotted from left to right which is helpful if ',...
'you have selected the images in the order of different sub-groups. Furthermore this is also useful for fMRI images which can be ',...
'also used with this tool. The proportional scaling option should be only used if image intensity is not scaled (e.g. T1 images) ',...
'or if images have to be scaled during statistical analysis (e.g. modulated images).']};

%------------------------------------------------------------------------

data.help = {[...
'Select all images. Images have to be in the same orientation with same voxel size and dimension ',...
'(e.g. normalized images)']};

showslice_slice = cfg_entry;
showslice_slice.tag = 'slice';
showslice_slice.name = 'Slice (in mm)?';
showslice_slice.strtype = 'e';
showslice_slice.num = [1 1];
showslice_slice.val  = {0};
showslice_slice.help = {[...
'Choose slice in mm.']};

showslice = cfg_exbranch;
showslice.tag = 'showslice';
showslice.name = 'Display one slice for all images';
showslice.val = {data,scale,showslice_slice};
showslice.prog   = @cg_showslice_all;
showslice.help = {[...
'This function displays a selected slice for all images and indicates the respective filenames which is useful to check image quality ',...
'for a large number of files in a circumscribed region (slice).']};

%------------------------------------------------------------------------

calcvol_files = cfg_files;
calcvol_files.tag  = 'data';
calcvol_files.name = 'Volumes';
calcvol_files.filter = '*';
calcvol_files.ufilter = 'seg8.*\.txt$';
calcvol_files.num     = [1 Inf];
calcvol_files.help = {[...
'Select all *_seg8.txt files containing raw volumes, which were saved by VBM8 toolbox.']};

calcvol_name = cfg_entry;
calcvol_name.tag = 'calcvol_name';
calcvol_name.name = 'Output file';
calcvol_name.strtype = 's';
calcvol_name.num = [1 Inf];
calcvol_name.val  = {'raw_volumes.txt'};
calcvol_name.help  = {[...
'The output file is written to current working directory ',...
'unless a valid full pathname is given']};

calcvol = cfg_exbranch;
calcvol.tag = 'calcvol';
calcvol.name = 'Read raw volumes (GM/WM/CSF/Total)';
calcvol.val = {calcvol_files,calcvol_name};
calcvol.prog   = @execute_calcvol;
calcvol.help = {[...
'This function reads raw volumes for GM/WM/CSF/Total and saves values in a txt-file. ',...
'These values can be read with the matlab command: vol = spm_load. The values for GM/WM/CSF/TOTAL ',...
'are now saved in vol(:,1) vol(:,2) vol(:,3) and vol(:,4).'],...
'',[...
'You can use these variables either as nuisance in an AnCova model or as user-specified globals with ',...
'the "global calculation" option. Depending on your hypothesis and/or your data you can just use gray ',...
'matter ("gm") or calculate the sum of gray/white matter with "gm+wm". The use of raw volumes as ',...
'nuisance or globals is only recommended for modulated data. These data are corrected for size changes ',... 
'due to spatial  normalization and are thought to be in raw (un-normalized) space. In contrast, un-modulated ',...
'data are yet corrected for differences in size due to spatial normalization to a ',...
'reference brain and there is no need to correct for these differences again.']};

%------------------------------------------------------------------------

tools = cfg_choice;
tools.name = 'Tools';
tools.tag  = 'tools';
tools.values = {check_sd,showslice,calcvol,T2x};

return

%------------------------------------------------------------------------
function execute_calcvol(p)
%
% calculate raw volumes all tissue classes
%

fprintf('%35s\t%5s\t%5s\t%5s\t%5s\n','Name','GM','WM','CSF','Total');
fid = fopen(p.calcvol_name,'w');
for i=1:length(p.data)
	tmp = load(deblank(p.data{i}));
    [pth,nam]     = spm_fileparts(p.data{i});
	if numel(tmp)==3
		fprintf(fid,'%3.2f\t%3.2f\t%3.2f\t%3.2f\n',tmp(1),tmp(2),...
			tmp(3),sum(tmp));
		fprintf('%35s\t%5.2f\t%5.2f\t%5.2f\t%5.2f\n',nam(1:end-4),tmp(1),tmp(2),...
			tmp(3),sum(tmp));
	else
		error(['Wrong format in ' p.data{i}]);
	end
end
if fclose(fid)==0
	fprintf('\nValues saved in %s.\n',p.calcvol_name);
end

return
%------------------------------------------------------------------------

function tools = cg_vbm8_tools
% wrapper for calling VBM utilities
%
%_______________________________________________________________________
% $Id$

entry = inline(['struct(''type'',''entry'',''name'',name,',...
        '''tag'',tag,''strtype'',strtype,''num'',num,''help'',{{}})'],...
        'name','tag','strtype','num');

files = inline(['struct(''type'',''files'',''name'',name,',...
        '''tag'',tag,''filter'',fltr,''num'',num,''help'',{{}})'],...
        'name','tag','fltr','num');

mnu = inline(['struct(''type'',''menu'',''name'',name,',...
        '''tag'',tag,''labels'',{labels},''values'',{values},''help'',{{}})'],...
        'name','tag','labels','values');

branch = inline(['struct(''type'',''branch'',''name'',name,',...
        '''tag'',tag,''val'',{val},''help'',{{}})'],...
        'name','tag','val');

%_______________________________________________________________________

data = files('Data','data','image',[1 Inf]);
data.help = {[...
'Select raw data (e.g. T1 images) for processing. ',...
'This assumes that there is one scan for each subject. ',...
'Note that multi-spectral (when there are two or more registered ',...
'images of different contrasts) processing is not yet implemented ',...
'for this method.']};

%------------------------------------------------------------------------
%------------------------------------------------------------------------

data_label = files('Data','data','[cp]3.*\.[in][im][gi]$',[1 Inf]);
data_label.help = {[...
'Select segmented csf image. The segmented images are used for labeling where CSF is coded with "1", ',...
'GM with "2", and WM with "3". This numeration differs from that used in SPM, but is in ',...
'accordance to the commonly used labeling.'],...
'',[...
'To compute labeling you have to write all segmented images (csf, gm, and wm) before.']};

th_label    = entry('Threshold for label','th_label','e',[1 1]);
th_label.val = {0};
th_label.help = {[...
'Choose threshold for the sum of probabilities of GM+WM+CSF to mask label.']};

label      = branch('Label segmentations','label',{data_label,th_label});
label.prog   = @execute_label;
label.help = {[...
'Use segmented images to compute a labeled image where the tissue class with the maximal value ',...
'is used for labeling. The resulting file is indicated by the index 0.']};

%------------------------------------------------------------------------

data_T2x = files('Data','data','^spmT.*\.[in][im][gi]$',[1 Inf]);
data_T2x.help = {'Select spmT-images to convert.'};

T2x      = branch('Threshold and transform spmT-maps','T2x',{data_T2x});
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

check_sd_files = files('Data','data','image',[1 Inf]);
check_sd_files.help = {[...
'Select all images. Images have to be in the same orientation with same voxel size and dimension ',...
'(e.g. normalized images)']};

check_sd_scale    = mnu('Proportional scaling?','scale',{'no','yes'},{0,1});
check_sd_scale.val = {0};
check_sd_scale.help = {[...
'This option should be only used if image intensity is not scaled (e.g. T1 images) ',...
'or if images have to be scaled during statistical analysis (e.g. modulated images).']};

check_sd_sdname      = entry('Output standard deviation file','sd_name','s',[1 Inf]);
check_sd_sdname.val  = {'SD.nii'};
check_sd_sdname.help  = {[...
'The output file is written to current working directory ',...
'unless a valid full pathname is given. If you do not want to write this file ',...
'leave name empty.']};

check_sd_meanname      = entry('Output mean file','mean_name','s',[1 Inf]);
check_sd_meanname.val  = {'Mean.nii'};
check_sd_meanname.help  = {[...
'The output file is written to current working directory ',...
'unless a valid full pathname is given. If you do not want to write this file ',...
'leave name empty.']};

check_sd      = branch('Check sample homogeneity across sample','check_sd',{check_sd_files,check_sd_scale,check_sd_meanname,check_sd_sdname});
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

showslice_files = files('Data','data','image',[1 Inf]);
showslice_files.help = {[...
'Select all images. Images have to be in the same orientation with same voxel size and dimension ',...
'(e.g. normalized images)']};

showslice_scale    = mnu('Proportional scaling?','scale',{'no','yes'},{0,1});
showslice_scale.val = {0};
showslice_scale.help = {[...
'This option should be only used if image intensity is not scaled (e.g. T1 images) ',...
'or if images have to be scaled during statistical analysis (e.g. modulated images).']};

showslice_slice    = entry('Slice (in mm)?','slice','e',[1 1]);
showslice_slice.val = {0};
showslice_slice.help = {[...
'Choose slice in mm.']};

showslice      = branch('Display one slice for all images','showslice',{showslice_files,showslice_scale,showslice_slice});
showslice.prog   = @cg_showslice_all;
showslice.help = {[...
'This function displays a selected slice for all images and indicates the respective filenames which is useful to check image quality ',...
'for a large number of files in a circumscribed region (slice).']};

%------------------------------------------------------------------------

calcvol_files = files('Data','data','seg8.*\.txt$',[1 Inf]);
calcvol_files.help = {[...
'Select all *_seg8.txt files containing raw volumes, which were saved by VBM8 toolbox.']};

calcvol_name      = entry('Output Filename','calcvol_name','s',[1 Inf]);
calcvol_name.val  = {'raw_volumes.txt'};
calcvol_name.help  = {[...
'The output file is written to current working directory ',...
'unless a valid full pathname is given']};

calcvol      = branch('Read raw volumes (GM/WM/CSF/Total)','calcvol',{calcvol_files,calcvol_name});
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

tools.type = 'repeat';
tools.name = 'Tools';
tools.tag  = 'tools';
tools.values = {check_sd,showslice,calcvol,T2x,label};

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
%------------------------------------------------------------------------
function execute_label(p)
%
% calculate label of all tissue classes using probalities
%

spm_progress_bar('init',length(p.data),'Create labeling','completed');
for i=1:length(p.data)
    % image containing "p3" is csf
    csf = char(char(p.data{i}));
    [pth,nam,ext] = spm_fileparts(csf);
    % find "p3" or "c3" string
    ind = strfind(nam,'p3');
    if isempty(ind)
	    ind = strfind(nam,'c3');
    end
    % and replace it for gm and wm
    wm = fullfile(pth,[nam(1:ind(1)) '2' nam(ind(1)+2:end) ext]);
    gm = fullfile(pth,[nam(1:ind(1)) '1' nam(ind(1)+2:end) ext]);
    
    if ~exist(wm)
        error(['File ' wm ' not found']);
    end
    if ~exist(gm)
        error(['File ' gm ' not found']);
    end
    % use csf-gm-wm for correct labeling
    V = spm_vol(str2mat(csf,gm,wm));

    gwc = spm_read_vols(V);
    
    % label segmentations
    % label is only defined if sum of all tissues is > 0.1
    mask_gwc = (sum(gwc,4) > 0);

    label = zeros(V(1).dim(1:3));
    label(find((gwc(:,:,:,1) >= gwc(:,:,:,3)) & (gwc(:,:,:,1) >= gwc(:,:,:,2)))) = 1;
    label(find((gwc(:,:,:,2) >= gwc(:,:,:,3)) & (gwc(:,:,:,2) >= gwc(:,:,:,1)))) = 2;
    label(find((gwc(:,:,:,3) >= gwc(:,:,:,1)) & (gwc(:,:,:,3) >= gwc(:,:,:,2)))) = 3;

    clear gwc
    label = label.*mask_gwc;
    
    % replace p3 or c3 with pl/cl
    V(1).fname = fullfile(pth,[nam(1:ind(1)) '0' nam(ind(1)+2:end) ext]);
    V(1).pinfo = [0 0 1]';
    spm_write_vol(V(1),label);
    spm_progress_bar('set',i);
end
spm_progress_bar('clear');

return
%------------------------------------------------------------------------


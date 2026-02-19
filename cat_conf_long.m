function varargout = cat_conf_long(varargin)
% Configuration file for longitudinal data
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

newapproach = 0; 

if newapproach && nargin>0 
  [dep,varargout{1},varargout{2}] = vout_long(varargin{1});
  return
end

try
  expert = cat_get_defaults('extopts.expertgui');
catch %#ok<CTCH>
  expert = 0; 
end

% try to estimate number of processor cores
try
  numcores = cat_get_defaults('extopts.nproc');
  % because of poor memory management use only half of the cores for windows
  if ispc
    numcores = round(numcores/2);
  end
  numcores = max(numcores,1);
catch
  numcores = 0;
end

% force running in the foreground if only one processor was found or for compiled version
% or for Octave
if numcores == 1 || isdeployed || strcmpi(spm_check_version,'octave'), numcores = 0; end

%------------------------------------------------------------------------
nproc         = cfg_entry;
nproc.tag     = 'nproc';
nproc.name    = 'Split job into separate processes';
nproc.strtype = 'w';
nproc.val     = {numcores};
nproc.num     = [1 1];
nproc.hidden  = numcores <= 1 || isdeployed;
nproc.help    = {
    'In order to use multi-threading the CAT segmentation job with multiple subjects can be split into separate processes that run in the background. You can even close Matlab, which will not affect the processes that will run in the background without GUI. If you do not want to run processes in the background then set this value to 0.'
    ''
    'Keep in mind that each process needs about 1.5..2GB of RAM, which should be considered to choose the appropriate number of processes.'
    ''
    'Please further note that no additional modules in the batch can be run except CAT segmentation. Any dependencies will be broken for subsequent modules.'
  };
%------------------------------------------------------------------------
% files long with two different selection schemes
% - timepoints-subjects
% - subjects-timepoints
data              = cfg_files;
data.tag          = 'data';
data.name         = 'Volumes';
data.filter       = {'image','.*\.(nii.gz)$'};
data.ufilter      = '.*';
data.num          = [1 Inf];
data.help         = {'Select the same number and order of subjects for each time point. '};

timepoint         = data; 
timepoint.tag     = 'timepoints';
timepoint.name    = 'Timepoint';

timepoints        = cfg_repeat;
timepoints.tag    = 'timepoints';
timepoints.name   = 'Timepoints';
timepoints.values = {timepoint};
timepoints.num    = [2 Inf];
timepoints.help   = {'Specify time points. '};

subjlong          = data; 
subjlong.num      = [2 Inf];
subjlong.tag      = 'subjects';
subjlong.name     = 'Subject';
subjlong.help     = {'Select all longitudinal T1 images for this subject. '};

subjects          = cfg_repeat;
subjects.tag      = 'subjects';
subjects.name     = 'Subjects';
subjects.values   = {subjlong};
subjects.num      = [1 Inf];
subjects.help     = {'Specify subjects. '};

datalong          = cfg_choice;
datalong.tag      = 'datalong';
datalong.name     = 'Longitudinal data';
datalong.values   = {timepoints subjects};
datalong.val      = {subjects};
datalong.help     = {
 ['Select mode of longitudinal data selection for time points or subjects. ' ...
  'In case of "timepoints" you can create multiple time points where each time point has to contain the same number and order of subjects. ' ...
  'If you have a varying number of time points for each subject you have to use the "subjects" mode where you have to define the files of each subject separately. ']
}; 

%------------------------------------------------------------------------
longmodel        = cfg_menu;
longmodel.tag    = 'longmodel';
longmodel.name   = 'Longitudinal Model';
longmodel.labels = {
  'Optimized for detecting small changes (i.e. plasticity/learning effects)', ...
  'Optimized for detecting large changes (i.e. aging effects)', ...
  'Optimized for detecting large changes with brain/head growth (i.e. developmental effects)', ...
  'Save both plasticity and aging models'};
longmodel.values = {1 2 0 3};
longmodel.val    = {3};
if expert 
  % Add the internal values and the special plasticity & aging model for 
  % developer only because it is not fully working now (RD20220317).
  longmodel.labels{1} = [longmodel.labels{1}(1:end-1) '; 1)']; 
  longmodel.labels{2} = [longmodel.labels{2}(1:end-1) '; 2)']; 
  longmodel.labels{3} = [longmodel.labels{3}(1:end-1) '; 3)']; 
  longmodel.labels{4} = [longmodel.labels{3}(1:end-1) '; 0)']; 
  if expert > 1
    longmodel.labels{5} = [longmodel.labels{3}(1:end-1) ' V2; 4)']; 
    longmodel.values{5} = 4;
  end
end
longmodel.help = {
'The longitudinal pre-processing in CAT has been developed and optimized to detect subtle effects over shorter periods of time (e.g. brain plasticity or training effects after a few weeks or even shorter periods of time) and is less sensitive to detect larger changes over longer periods of time (e.g. ageing or developmental effects). To detect larger effects, we also offer a model that additionally takes into account deformations between time points. The use of deformations between the time points makes it possible to estimate and detect larger changes, while subtle effects over shorter periods of time in the range of weeks or a few months can be better detected with the model for small changes.'
''
'Unlike the plasticity and ageing models, the developmental pipeline must include a time point-independent registration to adjust the growth of the brain/head. '
''
'Please note that due to the additional warping and modulation steps, the resulting files are saved with "mwmwp1r" for gray matter instead of "mwp1r"'
''
};
if expert % add some further detail for the combined model that is only available for experts
  longmodel.help{end-3} = [longmodel.help{end-1} ...
    ' It therefore also requires independent processing and cannot be processed together with the other two models.'];
end

% The heavy option is at the limit and the images starts to look artifical.
% However, this could be relevant of strong artifacts and plasticity studies.  
bstr                 = cfg_menu;
bstr.tag             = 'bstr';
bstr.name            = 'Strength of final longitudinal bias correction (IN DEVELOPMENT)';
bstr.labels          = {'no correction','light','medium','strong'}; %,'heavy'};
bstr.values          = {0,0.25,0.5,0.75}; %,1.0
bstr.val             = {0};
bstr.hidden          = expert<2; 
bstr.help            = {
  'Strength of final longitudinal bias correction that utilize the average segmentation for further subtile corrections. Test also higher SPM bias correction that also incooperates the information from the average by using the longTPM. Use weaker corrections if the points in time are far apart or if the imgages are less affected by inhomogeneities. Only use stong corrections in case of severe inhomogeneities or artefacts and check the results! '
  'This correction was introduced in CAT12.7 (2020/10) and is still under test! So use it carefully! '
  ''
};

prepavg              = cfg_menu;
prepavg.tag          = 'prepavg';
prepavg.name         = 'Optimize orignal data before averaging (IN DEVELOPMENT)';
prepavg.labels       = {'no preparation','denoising','denosing+trimming'};
prepavg.values       = {0,1,2};
prepavg.val          = {2};
prepavg.hidden       = expert<2; 
prepavg.help         = {
  'Denoising, trimming and intensity scaling of the original timepoint data before creating the average. '
  ''
};

avgLASWMHC              = cfg_menu;
avgLASWMHC.tag          = 'avgLASWMHC';
avgLASWMHC.name         = 'Handling of LAS and WMHC on the average (IN DEVELOPMENT)';
avgLASWMHC.labels       = {'classic (AVG=TP)','reduced LAS','reduce LAS & extra WMH class','reduce LAS in TPs & extra WMH class'};
avgLASWMHC.values       = {0,1,2,3};
avgLASWMHC.val          = {0};
avgLASWMHC.hidden       = expert<2; 

avgLASWMHC.help         = {
 ['The use of the LAS for the creation of the indiviudal TPM as well as in the time point specific processing can result in a overestimation of subcortical GM. ' ...
  'A lower correction in the average or the time point is therefore maybe better suited. '];
 ['The correction of WMHs to the WM is most similar to the original SPM TPM in principle. ' ...
  'However, the many GM-like values of large WMHs seams to bias the WM peak resulting in WM over- and GM underestimation. ' ...
  'Without or with temporar correction (WMHC=0, WMHC=1) the values of WMHs are maybe corrected to GM what can causes problems in the WMHC of the time points. ' ...
  'Using an extra WMH class (WMHC=3) may reduce this bias. '];
  ''
};

enablepriors          = cfg_menu;
enablepriors.tag      = 'enablepriors';
enablepriors.name     = 'Use priors for longitudinal data';
enablepriors.labels   = {'No','Yes'};
enablepriors.values   = {0 1};
enablepriors.val      = {1};
enablepriors.hidden   = expert<1;
enablepriors.help     = {
  'The average image is used as a first estimate for affine transformation, segmentation and surface extraction. The idea is that by initializing with the average image we can reduce random variations and improve the robustness and sensitivity of the entire longitudinal pipeline. Furthermore, it significantly increases the speed of the surface extraction.'
  ''
};

%------------------------------------------------------------------------
delete_temp        = cfg_menu;
delete_temp.tag    = 'delete_temp';
delete_temp.name   = 'Delete temporary files';
delete_temp.labels = {'No','Yes'};
delete_temp.values = {0 1};
delete_temp.val    = {1};
delete_temp.help = {
'Temporary files such as the native segmentations, deformation fields or any processed data from the average image are usually removed after preprocessing. However, if you like to keep these files (for debugging) you can enable this option.'
''
};

longTPM        = cfg_menu;
longTPM.tag    = 'longTPM';
longTPM.name   = 'Use longitudinal TPM from avg image';
longTPM.labels = {'No','Yes'};
longTPM.values = {0 1};
longTPM.val    = {1};
longTPM.hidden = expert<1; 
longTPM.help = {
'Use longitudinal TPM from average image.'
};

printlong         = cfg_menu;
printlong.tag     = 'printlong';
printlong.name    = 'Create CAT long report';
printlong.labels  = {'No','Yes (volume only)','Yes (volume and surfaces)'};
printlong.values  = {0 1 2};
printlong.hidden  = expert < 1;
printlong.def     = @(val)cat_get_defaults('extopts.print', val{:});
printlong.help    = {
  'Create final longitudinal CAT report that requires Java.'
};


%------------------------------------------------------------------------

% boundary box
bb          = cfg_entry;
bb.strtype  = 'r';
bb.num      = [inf inf];
bb.tag      = 'bb';
bb.name     = 'Bounding box';
bb.val      = {12}; 
bb.hidden   = expert < 1;
bb.help     = {
  'The bounding box describes the dimensions of the volume to be written starting from the anterior commissure in mm.  It should include the entire brain (or head in the case of the Boundary Box of the SPM TPM) and additional space for smoothing the image.  The MNI 9-mm boundary box is optimized for CATs MNI152NLin2009cAsym template and supports filter cores up to 10 mm.  Although this box support 12 mm filter sizes theoretically, slight interference could occur at the edges and larger boxes are recommended for safety. '
  'Additionally, it is possible to use the boundary box of the TPM or the template for special (animal) templates with strongly different boundary boxes. '
  ''
  'The boundary box or its id (BBid see table below) has to be entered. '
  ''
  '  NAME         BBID        BOUNDARY BOX                          SIZE ?             FILESIZE $   '
  '  TMP BB            0        boundary box of the template (maybe too small for smoothing!)         '
  '  TPM BB            1        boundary box of the TPM                                               ' 
  '  MNI SPM          16      [ -90  -126  -72;  90  90  108 ]      [121x145x121]      4.2 MB (100%)'
  '  MNI CAT          12      [ -84  -120  -72;  84  84    96 ]      [113x139x113]      3.8 MB ( 84%)'
  '  ? - for 1.5 mm; $ - for 1.5 mm uint8'
  ''
};


%------------------------------------------------------------------------
extopts = cat_conf_extopts(expert); 
if ~expert
  extopts.val = [ extopts.val {bb} ];
end
opts    = cat_conf_opts(expert);
output  = cat_conf_output(expert); 
%------------------------------------------------------------------------

% RD202007: Allow only lh+rh surface processing in long mode although this
%           this does not help to update default settings via function handle.         
clear FN; for vi = 1:numel(output.val), FN{vi} = output.val{vi}.tag; end
surf = find(cellfun('isempty',strfind(FN,'surface'))==0); 
output.val{surf}.labels = {'No','Yes'};
output.val{surf}.values = {0 1}; 


long = cfg_exbranch;
long.name = 'CAT: Segment longitudinal data';
long.tag  = 'long';
if newapproach % new way - not working
  
  % remove major output fields
  clear FN; for vi = 1:numel(output.val), FN{vi} = output.val{vi}.tag; end
  removefields = {'warps','jacobianwarped'};
  for vi = 1:numel(removefields)
    output.val(find(cellfun('isempty',strfind(FN,removefields{vi}))==0)) = []; 
  end
  
  % remove subfields
  removefields = {'native'};
  for vim = 1:numel(output.val)
    clear FN; for vi = 1:numel(output.val{vim}.val), FN{vi} = output.val{vim}.val{vi}.tag; end
    if numel(output.val{vim}.val)
      for vi = 1:numel(removefields)
        output.val{vim}.val(find(cellfun('isempty',strfind(FN,removefields{vi}))==0)) = []; 
      end
    end
  end

  if expert
    output.val = [output.val, delete_temp]; 
  end
  long.val  = {datalong,longmodel,prepavg,bstr,avgLASWMHC,nproc,opts,extopts,output}; 
  long.vout = @vout_long2;
else
  % old appraoch
  %------------------------------------------------------------------------
  modulate        = cfg_menu;
  modulate.tag    = 'modulate';
  modulate.name   = 'Modulated GM/WM segmentations';
  modulate.labels = {'No','Yes'};
  modulate.values = {0 1};
  modulate.val    = {1};
  modulate.help = {
  '"Modulation" is to compensate for the effect of spatial normalisation. Spatial normalisation causes volume changes due to affine transformation (global scaling) and non-linear warping (local volume change). After modulation the resulting modulated images are preserved for the total amount of grey matter signal in the normalised partitions. Thus, modulated images reflect the tissue volumes before spatial normalisation. However, the user is almost always interested in removing the confound of different brain sizes and there are many ways to apply this correction. In contrast to previous VBM versions I now recommend to use total intracranial volume (TIV) as nuisance parameter in an AnCova model. '
  ''
  'Please note that I do not use the SPM modulation where the original voxels are projected into their new location in the warped images because this method introduces aliasing artifacts. Here, I use the scaling by the Jacobian determinants to generate "modulated" data. '
  ''
  'For longitudinal data the modulation is actually not necessary because normalization estimates for one subject are the same for all time points and thus modulation will be also the same for all time points. However, modulation might be useful if you want to compare the baseline images in a cross-sectional design in order to test whether there are any differences between the groups at the beginning of the longitudinal study. '
  ''
  };


  %------------------------------------------------------------------------
  dartel        = cfg_menu;
  dartel.tag    = 'dartel';
  dartel.name   = 'DARTEL export';
  if expert
    dartel.labels = {'No','Rigid (SPM12 default)','Affine','Both'};
    dartel.values = {0 1 2 3};
  else
    dartel.labels = {'No','Rigid (SPM12 default)','Affine'};
    dartel.values = {0 1 2};
  end
  dartel.val    = {0};
  dartel.help   = {
  'This option is to export data into a form that can be used with DARTEL. The SPM default is to only apply rigid body transformation. However, a more appropriate option is to apply affine transformation, because the additional scaling of the images requires less deformations to non-linearly register brains to the template.'
  ''
  'Please note, that this option is only useful if you intend to create a customized DARTEl template for your longittudinal data. The DARTEL exported segmentations is saved for the average image of all time points for one subject (and also for the data of all time points) and can be used in order to create a customized template with the DARTEL toolbox. The resulting flow fields can be finally applied to the respective native segmentations (e.g. p1/p2 images) to obtain normalized segmentations according to the newly created DARTEL template.'
  ''
  };
  
  % extract only the surface, ROI, sROI, and BIDS menu
  FN      = cell(1,numel(output.val));  for fni=1:numel(output.val), FN{fni} = output.val{fni}.tag; end
  ROI     = output.val{ setdiff( find(cellfun('isempty',strfind(FN,'ROImenu'))==0) , ...
                     find(cellfun('isempty',strfind(FN,'sROImenu'))==0,1) ) }; 
  BIDS    = output.val{find(cellfun('isempty',strfind(FN,'BIDS'))==0)};
  surface = output.val{find(cellfun('isempty',strfind(FN,'surface'))==0)};

  output.val = {BIDS,surface};
  
  delete_temp.hidden = expert<1;
  
  long.val  = {datalong,longmodel,enablepriors,prepavg,bstr,avgLASWMHC,nproc,opts,extopts,output,ROI,longTPM,modulate,dartel,printlong,delete_temp};
  
% does not yet work! 
% long.vout = @vout_long;
end
long.prog = @cat_long_multi_run;

long.help = {
'This option offers customized processing of longitudinal data. Please note that two different longitudinal models are offered. The first longitudinal model is optimized for processing and capturing smaller changes over time in response to short-term plasticity effects (e.g. from learning and training). This model will probably not be as sensitive for larger longitudinal changes where large parts of the brain change over time (e.g. atrophy due to Alzheimers disease or aging). This is due to the effect of estimating the parameters of spatial registration  from the deformations of all time points and then applying them to all time points. If a large atrophy occurs between the time points, this can lead to a displacement of tissue boundaries and might result in areas with reduced volumes over time, which are surrounded by areas with increased volume due to these displacement problems. For data with larger volume changes over time, you should choose the longitudinal model, which is optimized to detect larger changes. This model also takes into account the deformations between time points and prevents the problems mentioned above.'
''
'Please note that surface-based preprocessing and ROI estimates are not affected by the selected longitudinal model, as the realigned images are used independently to create cortical surfaces, thickness, or ROI estimates.'
''
'Furthermore, we use an idea that was introduced by Reuter et al. for the Freesurfer software. Here the processing of the individual time points is initialized by the processed results from the (unbiased) average image. This reduces random variations in the processing procedure and improves the robustness and sensitivity of the overall longitudinal analysis.'
''
};

%------------------------------------------------------------------------
varargout{1} = long; 
return;
%------------------------------------------------------------------------
 

%------------------------------------------------------------------------
function dep = vout_long(job)

% this is not yet working!
if isfield(job.datalong,'subjects')
  job.subj = job.datalong.subjects;
else
  job.subj = job.datalong.timepoints;
end

for k=1:numel(job.subj)
    cdep            = cfg_dep;
    cdep.sname      = sprintf('Segmented longitudinal data (Subj %d)',k);
    cdep.src_output = substruct('.','sess','()',{k},'.','files');
    cdep.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    if k == 1
        dep = cdep;
    else
        dep = [dep cdep];
    end
end;

% add all surface/thickness files! of all time points and all subjects
if job.output.surface
    for k=1:numel(job.subj)
        dep(end+1)          = cfg_dep;
        dep(end).sname      = 'mwp1 Images';
        dep(end).src_output = substruct('.','mwp1','()',{':'});
        dep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    for k=1:numel(job.subj)
        dep(end+1)          = cfg_dep;
        dep(end).sname      = 'Left Central Surfaces';
        dep(end).src_output = substruct('.','surf','()',{':'});
        dep(end).tgt_spec   = cfg_findspec({{'filter','gifti','strtype','e'}});
    end
    for k=1:numel(job.subj)
        dep(end+1)          = cfg_dep;
        dep(end).sname      = 'Left Thickness';
        dep(end).src_output = substruct('.','thick','()',{':'});
        dep(end).tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}});
    end
    for k=1:numel(job.subj)
        dep(end+1)          = cfg_dep;
        dep(end).sname      = 'CAT Report';
        dep(end).src_output = substruct('.','catreport','()',{':'});
        dep(end).tgt_spec   = cfg_findspec({{'filter','xml','strtype','e'}});
    end
    for k=1:numel(job.subj)
        dep(end+1)          = cfg_dep;
        dep(end).sname      = 'ROI XML File';
        dep(end).src_output = substruct('.','catroi','()',{':'});
        dep(end).tgt_spec   = cfg_findspec({{'filter','xml','strtype','e'}});
    end
end


%------------------------------------------------------------------------

%------------------------------------------------------------------------
function [dep,out,inputs] = vout_long2(job)
    inputs = cell(1, numel(job.subj));

    [mrifolder, reportfolder, surffolder] = cat_io_subfolders(job.subj(1).mov,job);
    
    for i=1:numel(job.subj),
      %%
        out.subj(i).warps = cell(1,1);
        if iscell(job.subj(i).mov)
            [pth,nam,ext,num] = spm_fileparts(job.subj(i).mov{1});
        else
            [pth,nam,ext,num] = spm_fileparts(job.subj(i).mov);
        end
        out.subj(i).warps{1} = fullfile(pth,mrifolder,['avg_y_', nam, ext, num]);

        out.subj(i).files = cell(numel(cellstr(job.subj(i).mov)),1);
        m = numel(cellstr(job.subj(i).mov)); % number of scans of this subject
%%
        data = cell(m,1);
        for j=1:m
            if iscell(job.subj(i).mov)
              [pth,nam,ext,num] = spm_fileparts(job.subj(i).mov{j});
            else
              [pth,nam,ext,num] = spm_fileparts(job.subj(i).mov);
            end
            % for output (DEP)
            volumes = {
              'GM'  1; 
              'WM'  2; 
              'CSF' 3; 
              'WMH' 7;
              };

            for txi = 1:size(volumes,1)
              %%
              tissue = { 
                volumes{txi,1}  'warped' 1      sprintf('wp%d',volumes{txi,2})   sprintf('wp%dr',volumes{txi,2})    ''; 
                volumes{txi,1}  'mod'    [1 3]  sprintf('mwp%d',volumes{txi,2})  sprintf('mwp%dr',volumes{txi,2})   ''; 
                volumes{txi,1}  'mod'    [2 3]  sprintf('m0wp%d',volumes{txi,2}) sprintf('m0wp%dr',volumes{txi,2})  ''; 
                volumes{txi,1}  'dartel' [1 3]  sprintf('rp%da',volumes{txi,2})  sprintf('rp%dr',volumes{txi,2})    '_affine'; 
                volumes{txi,1}  'dartel' [2 3]  sprintf('rp%dr',volumes{txi,2})  sprintf('rp%dr',volumes{txi,2})    '_affine'; 
              };
              for ti = 1:size(tissue,1)
                if isfield(job.output,tissue{ti,1}) && isfield(job.output.(tissue{ti,1}),tissue{ti,2}) && ...
                  any( job.output.(tissue{ti,1}).(tissue{ti,2}) == tissue{ti,3} )
                  out.subj(i).(tissue{ti,4}){j,1} = fullfile(pth,mrifolder,[(tissue{ti,5}), nam, (tissue{ti,6}), ext, num]);
                end
              end
            end
            %%
            if isfield(job.output,'labelnative') && job.output.labelnative
              out.subj(i).p0{j,1}   = fullfile(pth,mrifolder,['p0r', nam, ext, num]);
            end
            if isfield(job.output.bias,'warped') && job.output.bias.warped
              out.subj(i).wm{j,1}   = fullfile(pth,mrifolder,['wmr', nam, ext, num]);
            end
            if job.output.surface
              out.subj(i).surface{j,1}   = fullfile(pth,surffolder,['lh.central.'  , nam, ext, num]);
              out.subj(i).thickness{j,1} = fullfile(pth,surffolder,['lh.thickness.', nam, ext, num]);
            end

          % for input
          if iscell(job.subj(i).mov)
            data{j} = job.subj(i).mov{j};
          else
            data{j} = cellstr(job.subj(i).mov);
          end
        end
        
        inputs{1,i} = data;
    end
    
    %%
    maps = {...
      ... 'files','warps', ... internal 
      'wp1','wp2','wp3','wp7',...             % unmodulated (warped==1)
      'mwp1','mwp2','mwp3','mwp7',...         % modulated (mod==1 | mod==3)
      'm0wp1','m0wp2','m0wp3','m0wp7', ...    % modulated (mod==2 | mod==3)
      'rp1a','rp2a','rp3a','rp7a',...         % dartel affine (dartel==1 |?dartel==3)
      'rp1r','rp2r','rp3r','rp7r', ...        % dartel rigid  (dartel==2 |?dartel==3)
      'p0','wm', ...                          % 
      'surface','thickness' ...               % surface
      };
    
    for mi=1:numel(maps)
        for k=1:numel(job.subj)
          if isfield(out.subj(k),maps{mi})
            cdep            = cfg_dep;
            cdep.sname      = sprintf('%s files of Subject %d',maps{mi},k);
            cdep.src_output = substruct('.','subj','()',{k},'.',maps{mi});
            cdep.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
            if ~exist('dep','var'); 
                dep = cdep;
            else
                dep = [dep cdep];
            end
          end
        end
    end
    %%
return
%------------------------------------------------------------------------

function Atlas_JulichBrain2CAT12
% Convert JulichBrain atlas to CAT12 atlas
%
% https://search.kg.ebrains.eu/instances/f1fe19e8-99bd-44bc-9616-a52850680777
%_______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

%Vlh = spm_vol(spm_select(1,'image','Select lh image',{},pwd,'JulichBrain_MPMAtlas_l_N10_nlin2Stdicbm152asym2009c'));
%Vrh = spm_vol(spm_select(1,'image','Select rh image',{},pwd,'JulichBrain_MPMAtlas_r_N10_nlin2Stdicbm152asym2009c'));
Vlh = spm_vol(fullfile(spm('dir'),'toolbox','cat12','development','JulichBrainAtlas_3.1_207areas_MPM_lh_MNI152.nii'));
Vrh = spm_vol(fullfile(spm('dir'),'toolbox','cat12','development','JulichBrainAtlas_3.1_207areas_MPM_rh_MNI152.nii'));
lh = spm_read_vols(Vlh);
rh = spm_read_vols(Vrh);

% set doubled defined areas at midline to 0
lhgt0 = lh>0;
rhgt0 = rh>0;
ind = find(lhgt0 & rhgt0);
lh(ind) = 0;
rh(ind) = 0;

% left/right coding
atlas = zeros(size(lh));
atlas(lhgt0) = 2*lh(lhgt0) - 1; 
atlas(rhgt0) = 2*rh(rhgt0); 
 
% replace remaining holes with median value
holes = atlas > 0;
holes = (cat_vol_morph(holes,'c') - holes) > 0;
tmp = atlas; %cat_vol_median3c(single(atlas));
atlas(holes) = double(tmp(holes));

Vlh.fname = fullfile(spm('dir'),'toolbox','cat12','development','julichbrain_refined.nii');
Vlh.descrip = 'JulichBrain v3.1';
Vlh.dt(1) = 4; 
spm_write_vol(Vlh,atlas);

%% roi information
%Slh = cat_io_xml(spm_select(1,'xml','Select lh xml',{},pwd,'JulichBrain_MPMAtlas_l_N10_nlin2Stdicbm152asym2009c'));
%Srh = cat_io_xml(spm_select(1,'xml','Select rh xml',{},pwd,'JulichBrain_MPMAtlas_r_N10_nlin2Stdicbm152asym2009c'));
Slh = cat_io_xml(fullfile(spm('dir'),'toolbox','cat12','development','JulichBrainAtlas_3.1_207areas_MPM_lh_MNI152.xml'));
Srh = cat_io_xml(fullfile(spm('dir'),'toolbox','cat12','development','JulichBrainAtlas_3.1_207areas_MPM_rh_MNI152.xml'));

fid = fopen( fullfile(spm('dir'),'toolbox','cat12','templates_MNI152NLin2009cAsym','julichbrain.csv'),'w');
fprintf(fid,'ROIid;ROIabbr;ROIname\n');

Xlh = Slh.Structures.Structure;
Xrh = Srh.Structures.Structure;
n_lh = numel(Xlh);
n_rh = numel(Xrh);
n = n_lh + n_rh;
Name = cell(n,1);
ID = zeros(n,1);

for i=1:n_lh
  Name{(2*i)-1} = Xlh(i).CONTENT;
  ID((2*i)-1) = 2*Xlh(i).ATTRIBUTE.grayvalue-1;
end

for i=1:n_rh
  Name{2*i} = Xrh(i).CONTENT;
  ID(2*i) = 2*Xrh(i).ATTRIBUTE.grayvalue;
end

% write csv file for CAT12
for i=1:n
  % create short name without not needed information
  ShortName = strrep(Name{i},'Area','');
  ind = strfind(ShortName,'(');
  if ~isempty(ind)
    ShortName = ShortName(1:ind-1);
  end
  ShortName = strrep(ShortName,' ','');
  
  % write left/right information
  if rem(i,2)
    fprintf(fid,'%d;l%s;Left %s\n',ID(i),ShortName,Name{i});
  else
    fprintf(fid,'%d;r%s;Right %s\n',ID(i),ShortName,Name{i});
  end
end

fclose(fid);

% use imcalc to get to template space
matlabbatch{1}.spm.util.imcalc.input = {
                                        fullfile(cat_get_defaults('extopts.pth_templates'),'aal3.nii')
                                        Vlh.fname
                                        };
matlabbatch{1}.spm.util.imcalc.output = 'julichbrain3.nii';
matlabbatch{1}.spm.util.imcalc.outdir = {fullfile(spm('dir'),'toolbox','cat12','templates_MNI152NLin2009cAsym')};
matlabbatch{1}.spm.util.imcalc.expression = 'i2';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 0;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

spm_jobman('run',matlabbatch)

V = spm_vol('julichbrain3.nii');
vol = round(spm_read_vols(V));
V.pinfo(1) = 1;
spm_write_vol(V,vol);


%% Surfaces

% Jubrain surfaces
Alh = gifti(fullfile(spm('dir'),'toolbox','cat12','development','lh.JulichBrainAtlas_3.1.label.gii')); 
Arh = gifti(fullfile(spm('dir'),'toolbox','cat12','development','rh.JulichBrainAtlas_3.1.label.gii')); 

% template surfaces
Slh = gifti(fullfile(spm('dir'),'toolbox','cat12','templates_surfaces','lh.sphere.freesurfer.gii')); 
Srh = gifti(fullfile(spm('dir'),'toolbox','cat12','templates_surfaces','rh.sphere.freesurfer.gii')); 

% template surfaces
S32lh = gifti(fullfile(spm('dir'),'toolbox','cat12','templates_surfaces_32k','lh.sphere.freesurfer.gii')); 
S32rh = gifti(fullfile(spm('dir'),'toolbox','cat12','templates_surfaces_32k','rh.sphere.freesurfer.gii')); 

% simple mapping via sphere
A32lh = cat_surf_fun('cdatamapping', S32lh, Slh, Alh.cdata);
A32rh = cat_surf_fun('cdatamapping', S32rh, Srh, Arh.cdata);

% prepare values for FreeSurfer annot standard (correct for undefined entry)
rgbalh = int64(round(Alh.labels.rgba(2:end,:)*255)); rgbalh(:,4) = 0; 
rgbarh = int64(round(Arh.labels.rgba(2:end,:)*255)); rgbarh(:,4) = 0; 

% save maps
fname = {
  fullfile(spm('dir'),'toolbox','cat12','atlases_surfaces','lh.JulichBrainAtlas_3.1.freesurfer.annot'); 
  fullfile(spm('dir'),'toolbox','cat12','atlases_surfaces','rh.JulichBrainAtlas_3.1.freesurfer.annot');
  }; 
Alhfs = Alh.cdata; Arhfs = Arh.cdata;
for i=1:max(Alh.cdata)
  Alhfs(Alh.cdata==i) = rgbalh(i,1) + rgbalh(i,2)*2^8 + rgbalh(i,3)*2^16 + rgbalh(i,4)*2^24; 
end
for i=1:max(A32rh)
  Arhfs(Arh.cdata==i) = rgbarh(i,1) + rgbarh(i,2)*2^8 + rgbarh(i,3)*2^16 + rgbarh(i,4)*2^24; 
end
cat_io_FreeSurfer('write_annotation', fname{1}, ...
  1:numel(Alh.cdata), Alhfs, struct('numEntries',max(Alh.cdata), 'orig_tab','JulichBrainAtlas_3.1.label', ...
    'struct_names', {Alh.labels.name'},'table',rgbalh)); 
cat_io_FreeSurfer('write_annotation', fname{2}, ...
  1:numel(Arh.cdata), Arhfs, struct('numEntries',max(Arh.cdata), 'orig_tab','JulichBrainAtlas_3.1.label', ...
    'struct_names', {Alh.labels.name'},'table',rgbarh)); 
cat_io_cprintf('blue','  %s\n',fname{1})


% save 32k data
fname32 = {
  fullfile(spm('dir'),'toolbox','cat12','atlases_surfaces_32k','lh.JulichBrainAtlas_3.1.freesurfer.annot'); 
  fullfile(spm('dir'),'toolbox','cat12','atlases_surfaces_32k','rh.JulichBrainAtlas_3.1.freesurfer.annot');
  }; 
A32lhfs = A32lh; A32rhfs = A32rh;
for i=1:max(A32lh)
  A32lhfs(A32lh==i) = rgbalh(i,1) + rgbalh(i,2)*2^8 + rgbalh(i,3)*2^16 + rgbalh(i,4)*2^24; 
end
for i=1:max(A32rh)
  A32rhfs(A32lh==i) = rgbarh(i,1) + rgbarh(i,2)*2^8 + rgbarh(i,3)*2^16 + rgbarh(i,4)*2^24; 
end
cat_io_FreeSurfer('write_annotation', fname32{1}, ...
  1:numel(A32lh), A32lhfs, struct('numEntries',max(A32lh), 'orig_tab','JulichBrainAtlas_3.1.label', ...
    'struct_names', {Alh.labels.name'},'table',rgbalh));  
cat_io_FreeSurfer('write_annotation', fname32{2}, ...
  1:numel(A32rh), A32rhfs, struct('numEntries',max(A32rh), 'orig_tab','JulichBrainAtlas_3.1.label', ...
    'struct_names', {Alh.labels.name'},'table',rgbarh));
cat_io_cprintf('blue','  %s\n',fname32{1})















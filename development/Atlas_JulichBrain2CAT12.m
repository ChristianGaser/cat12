function Atlas_JulichBrain2CAT12
% Convert JulichBrain atlas to CAT12 atlas
%_______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

Vlh = spm_vol(spm_select(1,'image','Select lh image',{},pwd,'JulichBrain_MPMAtlas_l_N10_nlin2Stdicbm152asym2009c'));
Vrh = spm_vol(spm_select(1,'image','Select rh image',{},pwd,'JulichBrain_MPMAtlas_r_N10_nlin2Stdicbm152asym2009c'));
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
tmp = cat_vol_median3c(single(atlas));
atlas(holes) = double(tmp(holes));

Vlh.fname = 'julichbrain_refined.nii';
Vlh.descrip = 'JulichBrain v2.2';
spm_write_vol(Vlh,atlas);

Slh = cat_io_xml(spm_select(1,'xml','Select lh xml',{},pwd,'JulichBrain_MPMAtlas_l_N10_nlin2Stdicbm152asym2009c'));
Srh = cat_io_xml(spm_select(1,'xml','Select rh xml',{},pwd,'JulichBrain_MPMAtlas_r_N10_nlin2Stdicbm152asym2009c'));

fid = fopen('julichbrain.csv','w');
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

matlabbatch{1}.spm.util.imcalc.input = {
                                        fullfile(cat_get_defaults('extopts.pth_templates'),'aal3.nii')
                                        'julichbrain_refined.nii,1'
                                        };
matlabbatch{1}.spm.util.imcalc.output = 'julichbrain.nii';
matlabbatch{1}.spm.util.imcalc.outdir = {''};
matlabbatch{1}.spm.util.imcalc.expression = 'i2';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 0;
matlabbatch{1}.spm.util.imcalc.options.dtype = 2;

spm_jobman('run',matlabbatch)

V = spm_vol('julichbrain.nii');
vol = round(spm_read_vols(V));
V.pinfo(1) = 1;
spm_write_vol(V,vol);

function correct_HCP_annot
% replace right hemipshere names and labels with that from left hemipshere to be compatible to other
% annot-files
%_______________________________________________________________________
% Christian Gaser
% $Id$

[vl, fl, cl] = cat_io_FreeSurfer('read_annotation','lh.HCP-MMP1.annot');
[vr, fr, cr] = cat_io_FreeSurfer('read_annotation','rh.HCP-MMP1.annot');

n = cl.numEntries;

for i=1:n
  struct_names = cl.struct_names(i);
  
  % don't change first "???" entry
  if i>1
    % remove L_/R_ and _ROI parts of the name
    cl.struct_names(i) = {struct_names{1}(3:end-4)};
    cr.struct_names(i) = {struct_names{1}(3:end-4)};
  end
  
  % replace labels
  fr(find(fr==cr.table(i,5))) = cl.table(i,5);
  
  % correct colortable
  cr.table(i,:) = cl.table(i,:);
end

cat_io_FreeSurfer('write_annotation', 'atlases_surfaces/lh.aparc_HCP_MMP1.freesurfer.annot', vl, fl, cl);
cat_io_FreeSurfer('write_annotation', 'atlases_surfaces/rh.aparc_HCP_MMP1.freesurfer.annot', vr, fr, cr);

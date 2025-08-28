function cat_stat_coord_diff(P)
% compute coordinate differences and eucldiean distance between two 
% or more surface meshes. The difference and the euclidean distance is always 
% estimated to the first surface: surfX - surf1
%
% output name will be diff_{name_surfX}
% and euclidean_{name_surfX}
%
% FORMAT cat_stat_coord_diff(P)
% P     - filenames for image 1..X
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________

if nargin == 0
  % select images for each subject
  don = 0;
  for i = 1:1000,
    Ps = spm_select([0 Inf],'^lh.(central|pial|white)',['Select all surfaces for time point ' num2str(i)]);
    if size(Ps,2) < 1, don = 1; ; break; end;
    P{i} = Ps;
  end
end

coord = {'x','y','z'};
hemi = {'lh.','rh.'};

avg_vertices = cell(1,numel(P));

% go through all time points
for side = hemi
  for i = 1:numel(P)
    sname = char(strrep(cellstr(P{i}),'lh.',side));
    
    G = gifti(sname);
    avg_vertices{i} = zeros(size(G(1).vertices));
    
    % go through all subjects for that time point and average coordinates
    for j = 1:numel(G)
      avg_vertices{i} = avg_vertices{i} + G(j).vertices;
    end
    avg_vertices{i} = avg_vertices{i}/numel(G);
    
    % calc difference to baseline surface
    sinfo = cat_surf_info(sname(1,:));
    if i > 1
      fprintf('Calculate s%d-s1\n',i);
      diff_coord = cell(1,3);
      for k=1:3
        diff_coord{k} = avg_vertices{i}(:,k) - avg_vertices{1}(:,k);
        outname = cat_surf_rename(sname(1,:),'ee','','dataname',['diff' coord{k} '_' sinfo.dataname])
        cat_io_FreeSurfer('write_surf_data',outname{1},diff_coord{k});
      end  
      outname = cat_surf_rename(sname(1,:),'ee','','dataname',['euclidean_' sinfo.dataname]);
      e = sqrt(diff_coord{1}.^2+diff_coord{2}.^2+diff_coord{3}.^2);
      cat_io_FreeSurfer('write_surf_data',outname{1},e);
    else
      % write zeros for the 1st image
      outname = cat_surf_rename(sname(1,:),'ee','','dataname',['euclidean_' sinfo.dataname]);
      e = zeros(size(avg_vertices{1}(:,1)));
      cat_io_FreeSurfer('write_surf_data',outname{1},e);
    end
  end
end

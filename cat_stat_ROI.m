function cat_stat_ROI(p)
%cat_stat_ROI to save mean values inside ROI for many subjects
%
%_______________________________________________________________________
% Christian Gaser
% $Id: cat_stat_ROI.m 834 2016-01-18 13:26:38Z gaser $

% ROI measures to search for
ROI_measures = char('Vgm','Vwm','Vcsf','mean_thickness');
n_ROI_measures = size(ROI_measures,1);
n_data = length(p.roi_xml);

[path, roi_name, ext] = fileparts(p.calcroi_name);

for i=1:n_data        
  xml = convert(xmltree(deblank(p.roi_xml{i})));

  if ~isfield(xml,'ROI')
    error('XML file contains no ROI information.');
  end

  % remove leading catROI*_ part from name
  [path2, ID] = fileparts(p.roi_xml{i});
  ind = strfind(ID,'_');
  ID = ID(ind(1)+1:end);

  atlases = fieldnames(xml.ROI);
  n_atlases = numel(atlases);
  
  for j=1:n_atlases
    if ~isfield(xml.ROI.(atlases{j}),'tr')
      error('Missing mandatory tr-field in XML file.');
    end
  
    n_ROIs = numel(xml.ROI.(atlases{j}).tr) - 1; % ignore header
    hdr = xml.ROI.(atlases{j}).tr{1}.td;
    
    for k=1:numel(hdr)
      for l=1:n_ROI_measures

        % check for field with ROI names
        if strcmp(hdr{k},'ROIappr') || strcmp(hdr{k},'ROIabbr') || strcmp(hdr{k},'lROIname') || strcmp(hdr{k},'rROIname')
          name_index = k;  
        end

        % look for pre-defined ROI measures
        if strcmp(hdr{k},deblank(ROI_measures(l,:)))
        
          % create filename with information about atlas and measures and print ROI name
          if (i==1) 
            out_name = fullfile(path,[ roi_name '_' deblank(atlases{j}) '_' hdr{k} '.csv']);
            fid{j,k} = fopen(out_name,'w');
            fprintf('Save values in %s\n',out_name);

            fprintf(fid{j,k},'Name\t');
            for m=1:n_ROIs
              fprintf(fid{j,k},'%s\t',char(xml.ROI.(atlases{j}).tr{m+1}.td(name_index)));
            end
          end

          % print ROI values
          fprintf(fid{j,k},'\n%s\t',ID);
          for m=1:n_ROIs
            fprintf(fid{j,k},'%s\t',char(xml.ROI.(atlases{j}).tr{m+1}.td(k)));
          end
          
          % close files after last dataset
          if (i==n_data)
            fclose(fid{j,k});
          end
                              
        end        
      end
    end
  end
end

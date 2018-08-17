function check_pipeline_ROIs(P)
%check_pipeline_ROIs to check distribution of ROI-values for different
% cat12 releases using check_pipeline.sh
%
% _________________________________________________________________________
% $Id$

if nargin == 0
  P = spm_select(Inf,'csv','Select csv files');
end

for i=1:size(P,1)
  csv_name = deblank(P(i,:));
  C = cat_io_csv(csv_name);
  [pth, name] = spm_fileparts(csv_name);
  
  [m,n] = size(C);
  for j=1:m
    for k=1:n
      if strcmp(C{j,k},'NaN'), C{j,k} = 0; end
    end
  end

  rev = num2str(cell2mat(C(2:end,1)));
  
  roi_names = C(1,:);

  roi_values = cell2mat(C(2:end,2:end));
  median_values = median(roi_values);
  roi_values = roi_values - median_values;

  figure
  cat_plot_boxplot(roi_values,struct('names',rev,'violin',2,'showdata',1,'outliers',0,'ylim',[-3 3]));
  pos = get(gcf,'Position');
  pos(3:4) = [1000 500];
  set(gcf,'MenuBar','none','Name',name,'Position',pos)
  title('Neuromorphometric atlas: Difference to median [ml]')
  saveas(gcf,[name '.png']);
  
end
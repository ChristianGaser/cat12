function cat_stat_IQR(p)
%cat_stat_IQR to read weighted overall image quality (IQR) from xml-files
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id: 2558 2024-02-28 $

fid = fopen(p.iqr_name,'w');

if fid < 0
  error('No write access: check file permissions or disk space.');
end

cat_progress_bar('Init',length(p.data_xml),'Load xml-files','subjects completed')
for i=1:length(p.data_xml)
    xml = cat_io_xml(deblank(p.data_xml{i})); 
    try
      iqr = xml.qualityratings.IQR;
    catch % also try to use old versions
      try
        iqr = xml.QAM.QM.rms;
      catch % give up
        iqr = nan; 
      end
    end

    [pth,nam]     = spm_fileparts(p.data_xml{i});
    fprintf(fid,'%s\n',iqr);
    fprintf('%s\n',iqr);
    cat_progress_bar('Set',i);  
end
cat_progress_bar('Clear');


if fclose(fid)==0
  fprintf('\nValues saved in %s.\n',p.iqr_name);
end

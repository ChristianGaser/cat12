function cat_stat_IQR(p)
%cat_stat_IQR to read weighted overall image quality (SIQR) from xml-files
% For xml-files of older CAT versions without SIQR the former image quality
% rating (IQR) is used instead.
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

fid = fopen(p.iqr_name,'w');

if fid < 0
  error('No write access: check file permissions or disk space.');
end

nIQR = 0; % number of xml-files without SIQR
cat_progress_bar('Init',length(p.data_xml),'Load xml-files','subjects completed')
for i=1:length(p.data_xml)
    xml = cat_io_xml(deblank(p.data_xml{i}));
    try
      % SIQR replaced IQR, but older xml-files only contain IQR
      if isfield(xml.qualityratings,'SIQR')
        iqr = xml.qualityratings.SIQR;
      else
        iqr = xml.qualityratings.IQR; nIQR = nIQR + 1;
      end
    catch % also try to use old versions
      try
        iqr = xml.QAM.QM.rms; nIQR = nIQR + 1;
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

if nIQR
  cat_io_cprintf('warn',sprintf(['\nWarning: %d of %d xml-files are from older CAT versions without SIQR. ' ...
    'The former IQR is saved for these files, which is not directly comparable to SIQR.\n'],nIQR,length(p.data_xml)));
end

if fclose(fid)==0
  fprintf('\nValues saved in %s.\n',p.iqr_name);
end

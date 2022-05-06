function varargout = cat_stat_TIV(p)
%cat_stat_TIV to read total intracranial volume (TIV) from xml-files
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

if ~isfield(p,'calcvol_savenames')
  p.calcvol_savenames = 0;
end

if ~p.calcvol_TIV
  fprintf('%60s\t%7s\t%7s\t%7s\t%7s\t%7s\n','Name','Total','GM','WM','CSF','WMH');
end
fid = fopen(p.calcvol_name,'w');

if fid < 0
  error('No write access: check file permissions or disk space.');
end

if p.calcvol_TIV
  calcvol = zeros(length(p.data_xml),1);
else
  calcvol = zeros(length(p.data_xml),5);
end

spm_progress_bar('Init',length(p.data_xml),'Load xml-files','subjects completed')
for i=1:length(p.data_xml)
  xml = cat_io_xml(deblank(p.data_xml{i})); 
  try
    tmp  = xml.subjectmeasures.vol_abs_CGW; 
  catch % use nan
    if p.calcvol_TIV
      tmp = nan; 
    else
      tmp  = nan(1,5);
    end
  end
    
  if isfield(xml,'filedata')
    filename = spm_str_manip(deblank(fullfile(xml.filedata.path,xml.filedata.file)),'a60');
  else
    [pth,file] = fileparts(deblank(p.data_xml{i}));
    filename = spm_str_manip(deblank(p.data_xml{i}),'a60');
  end
  
  % only save TIV
  if p.calcvol_TIV
    switch p.calcvol_savenames
      case 0
        fprintf(fid,'%7.2f\n',sum(tmp));
      case 1
        fprintf(fid,'%s\t%7.2f\n',file,sum(tmp));
      case 2
        fprintf(fid,'%s\t%7.2f\n',fullfile(pth,file),sum(tmp));
    end
    calcvol(i) = sum(tmp);
    fprintf('%60s\t%7.2f\n',filename,sum(tmp));
  else % also save GM/WM/CSF 
    switch p.calcvol_savenames
      case 0
        fprintf(fid,'%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\n',sum(tmp),tmp(2),tmp(3),tmp(1),tmp(4));
      case 1
        fprintf(fid,'%s\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\n',file,sum(tmp),tmp(2),tmp(3),tmp(1),tmp(4));
      case 2
        fprintf(fid,'%s\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\n',fullfile(pth,file),sum(tmp),tmp(2),tmp(3),tmp(1),tmp(4));
    end
    calcvol(i,:) = [sum(tmp),tmp(2),tmp(3),tmp(1),tmp(4)];
    fprintf('%60s\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\n',filename,sum(tmp),tmp(2),tmp(3),tmp(1),tmp(4));
  end
  spm_progress_bar('Set',i);  
end
spm_progress_bar('Clear');

if fclose(fid)==0
  fprintf('\nValues saved in %s.\n',p.calcvol_name);
    if nargout == 1
      varargout{1}.calcvol = calcvol;
    end
end

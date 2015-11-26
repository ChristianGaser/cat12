function cat_stat_TIV(p)
%cat_stat_TIV to read total intracranial volume (TIV) from xml-files
%
%_______________________________________________________________________
% Christian Gaser
% $Id$

if ~p.calcvol_TIV
  fprintf('%35s\t%5s\t%5s\t%5s\t%5s\t%5s\n','Name','Total','GM','WM','CSF','WMH');
end
fid = fopen(p.calcvol_name,'w');

if fid < 0
	error('No write access: check file permissions or disk space.');
end

spm_progress_bar('Init',length(p.data_xml),'Load xml-files','subjects completed')
for i=1:length(p.data_xml)
    xml = convert(xmltree(deblank(p.data_xml{i})));
    tmp = eval(xml.QAS.SM.vol_abs_CGW);

    [pth,nam]     = spm_fileparts(p.data_xml{i});
    % only save TIV
    if p.calcvol_TIV
        fprintf(fid,'%3.2f\n',sum(tmp));
        fprintf('%35s\t%5.2f\n',nam(5:end),sum(tmp));
    else % also save GM/WM/CSF 
        fprintf(fid,'%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\n',sum(tmp),tmp(3),tmp(1),...
			tmp(2),tmp(4));
        fprintf('%35s\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\n',nam(5:end),sum(tmp),tmp(3),tmp(1),...
			tmp(2),tmp(4));
    end
    spm_progress_bar('Set',i);  
end
spm_progress_bar('Clear');


if fclose(fid)==0
	fprintf('\nValues saved in %s.\n',p.calcvol_name);
end

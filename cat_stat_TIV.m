function cat_stat_TIV
%cat_stat_TIV to read total intracranial volume (TIV) from xml-files
%
%_______________________________________________________________________
% Christian Gaser
% $Id: cat_stat_TIV.m 769 2015-11-17 23:06:11Z gaser $

P = spm_select(Inf,'xml','select xml-files');

spm_progress_bar('Init',size(P,1),'Load xml-files','subjects completed')
for i=1:size(P,1)
    xml = convert(xmltree(deblank(P(i,:))));

    fprintf('%s\n',xml.QAS.SM.vol_TIV);
    spm_progress_bar('Set',i);  
end
spm_progress_bar('Clear');
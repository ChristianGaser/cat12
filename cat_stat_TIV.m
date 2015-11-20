function TIV = cat_stat_TIV
%cat_stat_TIV to read total intracranial volume (TIV) from xml-files
%
%_______________________________________________________________________
% Christian Gaser
% $Id$

P = spm_select(Inf,'xml','select xml-files');

n = size(P,1);

TIV = zeros(n,1);
spm_progress_bar('Init',n,'Load xml-files','subjects completed')
for i=1:n
    xml = convert(xmltree(deblank(P(i,:))));
    TIV(i) = str2double(xml.QAS.SM.vol_TIV);
    spm_progress_bar('Set',i);  
end
spm_progress_bar('Clear');
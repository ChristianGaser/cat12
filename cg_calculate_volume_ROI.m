function cg_calculate_volume_ROI

Vatlas = spm_vol('Dartel_v17.img');
atlas = round(spm_read_vols(Vatlas));

fid = fopen('Dartel_v17.txt','r');
S = textscan(fid,'%d%s','Delimiter','\t');
fclose(fid)

for i=1:length(S{1})
    fprintf('%2d\t%s\n',S{1}(i),char(S{2}(i,:)));
end

ROIs = spm_input('Define ROI IDs',1,'e','',Inf);

P = spm_select(Inf,'^m.*','Select modulated images');
V = spm_vol(P);
n = length(V);

ind_left = cell(length(ROIs),1);
for j = 1:length(ROIs)
    ind_left{j} = find(atlas == ROIs(j)+1);
end

ind_right = cell(length(ROIs),1);
for j = 1:length(ROIs)
    ind_right{j} = find(atlas == ROIs(j));
end

fprintf('\n%10s\tleft/right volume\n','Name');

[tmp fname] = spm_str_manip(char(V.fname),'C');

scale = abs(det(V(1).mat(1:3,1:3)));
for i=1:n
    vol = spm_read_vols(V(i));
    sum_left  = 0;
    sum_right = 0;
    if ~isempty(fname)
        fprintf('%10s\t',char(fname.m(i,:)));
    else
        fprintf('%10s\t',char(V(i).fname));
    end
    for j=1:length(ROIs)
        left  = scale*sum(vol(ind_left{j}));
        right = scale*sum(vol(ind_right{j}));
        fprintf('%5.2f\t%5.2f\t',left,right);
        sum_left  = sum_left + left;
        sum_right = sum_right + right;
    end
   fprintf('%5.2f\t%5.2f\n',sum_left,sum_right);
end
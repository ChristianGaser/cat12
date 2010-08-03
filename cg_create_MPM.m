function se_create_MPM

global defaults
global st
spm_defaults

spm_figure('GetWin','Interactive');
Vb = spm_vol(spm_select(1,'image','Select Background',[],spm('Dir','se_anatomy')));
AnatMask = spm_vol(spm_select(1,'image','select AnatMask'));

MAPname = spm_input('Title of anatomical map',1,'s','');
smooth = spm_input('Do smoothed PMaps exist ?','+1','y/n',[1,0],2);
if ~smooth; kernel = spm_input('smoothing {FWHM in mm}','+1','i',6); end
MAP = struct('name',{},'GV',{},'ref',{},'smoothed',{},'VOL',{});

P = spm_select(Inf,'image','select PMaps');
V = spm_vol(P);
nrParts = length(V);

[tmp fname] = spm_str_manip(char(V.fname),'C');

for i=1:nrParts

    MAP(i).ref = V(i).fname;
    MAP(i).name = char(fname.m(i,:));
    MAP(i).smoothed = [MAP(i).ref(1:end-size(spm_str_manip(MAP(i).ref,'t'),2)) 's' spm_str_manip(MAP(i).ref,'t')];
end

if ~smooth
    spm('Pointer','Watch');
    spm('FigName','Smooth: working');
    spm_progress_bar('Init',nrParts,'Smoothing','Volumes Complete');
    for i = 1:nrParts
        Q = deblank(MAP(i).ref);
        [pth,nm,xt,vr] = fileparts(deblank(Q));
        U = fullfile(pth,['s' nm xt vr]);
        spm_smooth(Q,U,kernel);
        spm_progress_bar('Set',i);
    end
    spm_progress_bar('Clear');
    spm('Pointer');
end


count = 1;
for i=1:nrParts
    Vi(i) = spm_vol(MAP(i).ref);
    Vs(i) = spm_vol(MAP(i).smoothed);
    MAP(i).GV = count;
    count = count + 2;
    fprintf('%d\t%s\n',MAP(i).GV, MAP(i).name);
end

errorMap = se_check_PMap(Vi);
if errorMap
    warning('No original PBmaps!');
end


switch spm('FnBanner');
    case 'SPM2'
        Vo    = struct('fname',	[MAPname '.img'],...
            'dim',		[Vb.dim(1:3),spm_type('uint8')],...
            'mat',		Vb.mat,...
            'pinfo',	[1 0 0]',...
            'descrip',	'Maximum probability map');
    otherwise
        Vo    = struct('fname',	[MAPname '.img'],...
            'dim',		Vb.dim(1:3),...
            'dt',      [spm_type('uint8') spm_platform('bigend')],...
            'mat',		Vb.mat,...
            'pinfo',	[1 0 0]',...
            'descrip',	'Maximum probability map');
end

dmtx = 1; hold = 0;
n   = prod(size(Vi));

dat   = spm_read_vols(Vb);

Yp = zeros(size(dat));

spm_progress_bar('Init',Vo.dim(3),'MPM calculation');
save all
for p = 1:Vo.dim(3),
    if any(any(dat(:,:,p)))
        M = inv(spm_matrix([0 0 -p 0 0 0 1 1 1])*inv(Vo.mat)*Vi(i).mat);
        for i = 1:n
            X(:,:,i) = spm_slice_vol(Vi(i),M,Vo.dim(1:2),0);
        end
        if any(any(sum(X,3)>75))
            [h,j] = find(max(X,[],3)>75 | sum(X,3)>125);
            for px = 1:size(h,1)
                 if prod(size(find(X(h(px),j(px),:) == max(X(h(px),j(px),:))))) == 1
                    Yp(h(px),j(px),p) = MAP(find(X(h(px),j(px),:) == max(X(h(px),j(px),:)))).GV;
                else
                    infrage = find(X(h(px),j(px),:) == max(X(h(px),j(px),:)));
                    msk = [h(px)-1 j(px)-1 p-1; h(px)+0 j(px)-1 p-1; h(px)+1 j(px)-1 p-1;
                        h(px)-1 j(px)+0 p-1; h(px)+0 j(px)+0 p-1; h(px)+1 j(px)+0 p-1;
                        h(px)-1 j(px)+1 p-1; h(px)+0 j(px)+1 p-1; h(px)+1 j(px)+1 p-1;
                        h(px)-1 j(px)-1 p+0; h(px)+0 j(px)-1 p+0; h(px)+1 j(px)-1 p+0;
                        h(px)-1 j(px)+0 p+0; h(px)+0 j(px)+0 p+0; h(px)+1 j(px)+0 p+0;
                        h(px)-1 j(px)+1 p+0; h(px)+0 j(px)+1 p+0; h(px)+1 j(px)+1 p+0;
                        h(px)-1 j(px)-1 p+1; h(px)+0 j(px)-1 p+1; h(px)+1 j(px)-1 p+1;
                        h(px)-1 j(px)+0 p+1; h(px)+0 j(px)+0 p+1; h(px)+1 j(px)+0 p+1;
                        h(px)-1 j(px)+1 p+1; h(px)+0 j(px)+1 p+1; h(px)+1 j(px)+1 p+1;]';
                    surround = [];
                    for area = 1:size(infrage,1)
                        sx = spm_sample_vol(Vi(infrage(area)),msk(1,:),msk(2,:),msk(3,:),0);
                        sx(sx==0) = NaN;
                        surround(area) = mean(sx(~isnan(sx)));
                    end
                    if prod(size(find(surround == max(surround)))) == 1
                        Yp(h(px),j(px),p) = MAP(infrage(find(surround == max(surround)))).GV;
                    else
                        infrage =  infrage(find(surround == max(surround)));
                        sX = []; for i=1:prod(size(infrage)); sX(i) = spm_sample_vol(Vs(infrage(i)),h(px),j(px),p,0); end
                        if prod(size(infrage(find(sX == max(sX))))) == 1
                            Yp(h(px),j(px),p) = MAP(infrage(find(sX == max(sX)))).GV;
                        else
                            infrage = infrage(find(sX == max(sX)));
                            msk = [h(px)-1 j(px)-1 p-1; h(px)+0 j(px)-1 p-1; h(px)+1 j(px)-1 p-1;
                                h(px)-1 j(px)+0 p-1; h(px)+0 j(px)+0 p-1; h(px)+1 j(px)+0 p-1;
                                h(px)-1 j(px)+1 p-1; h(px)+0 j(px)+1 p-1; h(px)+1 j(px)+1 p-1;
                                h(px)-1 j(px)-1 p+0; h(px)+0 j(px)-1 p+0; h(px)+1 j(px)-1 p+0;
                                h(px)-1 j(px)+0 p+0; h(px)+0 j(px)+0 p+0; h(px)+1 j(px)+0 p+0;
                                h(px)-1 j(px)+1 p+0; h(px)+0 j(px)+1 p+0; h(px)+1 j(px)+1 p+0;
                                h(px)-1 j(px)-1 p+1; h(px)+0 j(px)-1 p+1; h(px)+1 j(px)-1 p+1;
                                h(px)-1 j(px)+0 p+1; h(px)+0 j(px)+0 p+1; h(px)+1 j(px)+0 p+1;
                                h(px)-1 j(px)+1 p+1; h(px)+0 j(px)+1 p+1; h(px)+1 j(px)+1 p+1;]';
                            surround = [];
                            for area = 1:size(infrage,1)
                                sx = spm_sample_vol(Vs(infrage(area)),msk(1,:),msk(2,:),msk(3,:),0);
                                sx(sx==0) = NaN;
                                surround(area) =  mean(sx(~isnan(sx)));
                            end
                            if prod(size(find(surround == max(surround)))) == 1
                                Yp(h(px),j(px),p) = MAP(infrage(find(surround == max(surround)))).GV;
                            else
                            whos 
                            infrage(1)
                                Yp(h(px),j(px),p) = MAP(infrage(1)).GV;
                            end
                        end
                    end
                end
            end
        end
    end
    spm_progress_bar('Set',p);
end
spm_progress_bar('Clear')


out = Yp;
clear X;
spm_progress_bar('Init',Vo.dim(3),'Filter 1');
for z=2:Vb.dim(3)-1,
    spm_progress_bar('Set',z);
    [h,j] = find(Yp(:,:,z)<100 & Yp(:,:,z)>0);
    for px = 1:size(h,1)
        x=h(px); y=j(px); p = z;
        surround = Yp(x-1:x+1,y-1:y+1,z-1:z+1);
        if sum(surround > 99)>13
            sux = []; for i=unique(surround(surround>0))'; sux(find([MAP.GV] == i)) = sum(surround == i);  end
            if (any(sux>13)),  out(x,y,z) = MAP(find(sux == max(sux))).GV;
            else
                X = [];
                for i=1:n; X(i) = spm_sample_vol(Vi(i),x,y,z,0);  end
                if sum(X == max(X)) == 1, out(x,y,z) = MAP(find(X == max(X))).GV;
                else
                    infrage = find(X == max(X)); surround = [];
                    msk = [h(px)-1 j(px)-1 p-1; h(px)+0 j(px)-1 p-1; h(px)+1 j(px)-1 p-1; h(px)-1 j(px)+0 p-1; h(px)+0 j(px)+0 p-1; h(px)+1 j(px)+0 p-1;
                        h(px)-1 j(px)+1 p-1; h(px)+0 j(px)+1 p-1; h(px)+1 j(px)+1 p-1; h(px)-1 j(px)-1 p+0; h(px)+0 j(px)-1 p+0; h(px)+1 j(px)-1 p+0;
                        h(px)-1 j(px)+0 p+0; h(px)+0 j(px)+0 p+0; h(px)+1 j(px)+0 p+0; h(px)-1 j(px)+1 p+0; h(px)+0 j(px)+1 p+0; h(px)+1 j(px)+1 p+0;
                        h(px)-1 j(px)-1 p+1; h(px)+0 j(px)-1 p+1; h(px)+1 j(px)-1 p+1; h(px)-1 j(px)+0 p+1; h(px)+0 j(px)+0 p+1; h(px)+1 j(px)+0 p+1;
                        h(px)-1 j(px)+1 p+1; h(px)+0 j(px)+1 p+1; h(px)+1 j(px)+1 p+1;]';
                    for area = 1:size(infrage,1)
                        sx = spm_sample_vol(Vi(infrage(area)),msk(1,:),msk(2,:),msk(3,:),0);
                        sx(sx==0) = NaN; surround(area) =  mean(sx(~isnan(sx)));
                    end
                    if sum(surround == max(surround)) == 1
                        out(x,y,z) = MAP(infrage(find(surround == max(surround)))).GV;
                    else, sX = [];
                        infrage = infrage(find(surround == max(surround)));
                        for i=1:numel(infrage); sX(i) = spm_sample_vol(Vs(infrage(i)),h(px),j(px),p,0); end
                        if prod(size(infrage(find(sX == max(sX))))) == 1
                            out(x,y,z) = MAP(infrage(find(sX == max(sX)))).GV;
                        else
                            infrage = infrage(find(sX == max(sX)));
                            msk = [h(px)-1 j(px)-1 p-1; h(px)+0 j(px)-1 p-1; h(px)+1 j(px)-1 p-1; h(px)-1 j(px)+0 p-1; h(px)+0 j(px)+0 p-1; h(px)+1 j(px)+0 p-1;
                                h(px)-1 j(px)+1 p-1; h(px)+0 j(px)+1 p-1; h(px)+1 j(px)+1 p-1; h(px)-1 j(px)-1 p+0; h(px)+0 j(px)-1 p+0; h(px)+1 j(px)-1 p+0;
                                h(px)-1 j(px)+0 p+0; h(px)+0 j(px)+0 p+0; h(px)+1 j(px)+0 p+0; h(px)-1 j(px)+1 p+0; h(px)+0 j(px)+1 p+0; h(px)+1 j(px)+1 p+0;
                                h(px)-1 j(px)-1 p+1; h(px)+0 j(px)-1 p+1; h(px)+1 j(px)-1 p+1; h(px)-1 j(px)+0 p+1; h(px)+0 j(px)+0 p+1; h(px)+1 j(px)+0 p+1;
                                h(px)-1 j(px)+1 p+1; h(px)+0 j(px)+1 p+1; h(px)+1 j(px)+1 p+1;]'; surround = [];
                            for area = 1:size(infrage,1)
                                sx = spm_sample_vol(Vs(infrage(area)),msk(1,:),msk(2,:),msk(3,:),0);
                                sx(sx==0) = NaN; surround(area) =  mean(sx(~isnan(sx)));
                            end
                            if prod(size(infrage(find(surround == max(surround))))) == 1
                                out(x,y,z) = MAP(infrage(find(surround == max(surround)))).GV;
                            end
                        end
                    end
                end
            end
        end
    end
end
spm_progress_bar('Clear')


Yp = out;

clear X;
spm_progress_bar('Init',Vo.dim(3),'Filter 2');
for z=2:Vb.dim(3)-1,
    spm_progress_bar('Set',z);
    [h,j] = find(Yp(:,:,z)>100);
    for px = 1:size(h,1)
        x=h(px); y=j(px); p = z;
        surround = Yp(x-1:x+1,y-1:y+1,z-1:z+1);
        if prod(size(find(surround > 99)))>16
            sux = []; for i=unique(surround(surround>0))'; sux(find([MAP.GV] == i)) = prod(size(find(surround == i)));  end
            if (any(sux>16)),
                out(x,y,z) = MAP(find(sux == max(sux))).GV;
            end
        end
    end
end
spm_progress_bar('Clear')


hemi = spm_read_vols(AnatMask);
ind_hemi = find(hemi>0);
% add 1 for left hemisphere
out(ind_hemi) = out(ind_hemi) + round(hemi(ind_hemi)) - 1;

mat = Vo.mat; M = Vo.mat;
Vo = spm_create_vol(Vo); for p=1:Vo.dim(3), Vo = spm_write_plane(Vo,out(:,:,p),p); end;  save([MAPname '.mat'],'mat','M');
switch spm('FnBanner');
    case 'SPM2'
        Vo= spm_close_vol(Vo);
    otherwise
end


MAP(1).MaxMap = Vo;

for locNr=1:size(MAP,2)
    MAP(locNr).XYZ = [];
    for z=1:Vb.dim(3)
        if any(any(out(:,:,z) == MAP(locNr).GV))
            [x, y] = find(out(:,:,z) == MAP(locNr).GV);
            MAP(locNr).XYZ = [MAP(locNr).XYZ [x y z*ones(size(x,1),1)]'];
        end
    end
    MAP(locNr).XYZmm = [];
    MAP(locNr).XYZmm(1,:) = MAP(locNr).XYZ(1,:)+MAP(1).MaxMap.mat(1,4);
    MAP(locNr).XYZmm(2,:) = MAP(locNr).XYZ(2,:)+MAP(1).MaxMap.mat(2,4);
    MAP(locNr).XYZmm(3,:) = MAP(locNr).XYZ(3,:)+MAP(1).MaxMap.mat(3,4);
end

MAP(1).orient = 1;

spm_progress_bar('Init',size(MAP,2),'PMap');
for p = 1:size(MAP,2)
    MAP(p).Z = spm_sample_vol(spm_vol(MAP(p).ref),MAP(p).XYZ(1,:),MAP(p).XYZ(2,:),MAP(p).XYZ(3,:),0);
    spm_progress_bar('Set',p);
end
spm_progress_bar('clear');

%AnatMask = spm_vol([spm('Dir','se_anatomy') filesep 'AnatMask.img']);
for i=1:size(MAP,2)
    tmp = spm_sample_vol(AnatMask,MAP(i).XYZ(1,:),MAP(i).XYZ(2,:),MAP(i).XYZ(3,:),0);
    MAP(i).LR = zeros(size(MAP(i).Z));
    MAP(i).LR(tmp == 2) = -1;
    MAP(i).LR(tmp == 1) = 1;
    MAP(i).VOL = [size(find(tmp == 2),2) size(find(tmp == 1),2)];
end



for i=1:size(MAP,2); MAP(i).PMap = spm_vol([spm('Dir','se_anatomy') filesep 'PMaps' filesep spm_str_manip(MAP(i).ref,'t')]); end

spm_figure('GetWin','Interactive');
spm_progress_bar('Init',size(MAP,2),'Preparing data');
for m=1:size(MAP,2)
    MAP(m).allXYZ  = [];
    MAP(m).allZ    = [];
    dat = spm_read_vols(MAP(m).PMap);
    for p = 1:MAP(m).PMap.dim(3)
        [i,j,z] = find(dat(:,:,p));
        if any(i), MAP(m).allXYZ = [MAP(m).allXYZ [i'; j'; p*ones(1,length(i))]]; MAP(m).allZ = [MAP(m).allZ (z')/25]; end
    end
    MAP(m).allLR = spm_get_data(AnatMask,MAP(m).allXYZ);
    MAP(m).allLR(MAP(m).allLR  == 2) = -1;
    spm_progress_bar('Set',m);
end
spm_figure('Clear','Interactive');

MAP = rmfield(MAP,'PMap');

try
    load(fullfile(spm('Dir','se_anatomy'),'cmap2.mat'));
    MAP(1).cmap = cmap;
end

try
    load(fullfile(spm('Dir','spm'),'rend','render_single_subj.mat'));
    if (exist('rend') ~= 1), % Assume old format...
        rend = cell(size(Matrixes,1),1);
        for i=1:size(Matrixes,1),
            rend{i}=struct('M',eval(Matrixes(i,:)),...
                'ren',eval(Rens(i,:)),...
                'dep',eval(Depths(i,:)));
            rend{i}.ren = rend{i}.ren/max(max(rend{i}.ren));
        end;
    end;

    for i=1:length(rend),
        rend{i}.max=0;
        rend{i}.data = cell(1,1);
        if issparse(rend{i}.ren),
            d = size(rend{i}.ren);
            B1 = spm_dctmtx(d(1),d(1));
            B2 = spm_dctmtx(d(2),d(2));
            rend{i}.ren = B1*rend{i}.ren*B2';
            rend{i}.dep = exp(B1*rend{i}.dep*B2')-1;
        end;
        msk = find(rend{i}.ren>1);rend{i}.ren(msk)=1;
        msk = find(rend{i}.ren<0);rend{i}.ren(msk)=0;
    end;
    MAP(1).rend = rend;
end

save([MAPname '_MPM.mat'],'MAP')

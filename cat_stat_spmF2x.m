function cat_stat_spmF2x(vargin)
%CG_SPMF2X transformation of F-maps to P, -log(P), R2 maps
%
% The following formulas are used:
%
% --------------------------------
% coefficient of determination R2
% --------------------------------
%
%           F*(n-1)
% R2 = ------------------
%        n-p + F*(n-1)
%
% --------------------------------
% p-value
% --------------------------------
%
% p = 1-spm_Fcdf
%
% --------------------------------
% log p-value
% --------------------------------
%
% -log10(1-P) = -log(1-spm_Fcdf)
%
% For the last case of log transformation this means that a p-value of p=0.99 (0.01)
% is transformed to a value of 2
%
% Examples:
% p-value   -log10(1-P)
% 0.1       1
% 0.05      1.3
% 0.01      2
% 0.001     3
% 0.0001    4
%
% All maps can be thresholded using height and extent thresholds and you can 
% also apply corrections for multiple comparisons based on family-wise error 
% (FWE) or false discovery rate (FDR). You can easily threshold and/or 
% transform a large number of spmF-maps using the same thresholds.
%
% Naming convention of the transformed files:
%   Type_Contrast_Pheight_K
%
%   Type:      P    - p-value
%              logP - log p-value
%              R2   - coefficient of determination
%              F    - F-value
%
%   Contrast:  name used in the contrast manager with replaced none valid 
%              strings
%    
%   Pheight:   p    - uncorrected p-value in % (p<0.05 will coded with "p5")
%              pFWE - p-value with FWE correction in %
%              pFDR - p-value with FDR correction in %
%              
%   Pextent:   pk    - uncorr. extent p-value in % (p<0.05 coded with "p5")
%              pkFWE - extent p-value with FWE correction in %
%
%   K:         extent threshold in voxels
%
%_______________________________________________________________________
% Christian Gaser
% $Id$

rev = '$Rev$';

if nargin == 1
    P = char(vargin.data_F2x);

    sel = vargin.conversion.sel;

    if isfield(vargin.conversion.threshdesc,'fwe')
        adjustment = 1;
        u0  = vargin.conversion.threshdesc.fwe.thresh05;
    elseif isfield(vargin.conversion.threshdesc,'fdr')
        adjustment = 2;
        u0  = vargin.conversion.threshdesc.fdr.thresh05;
    elseif isfield(vargin.conversion.threshdesc,'uncorr')
        adjustment = 0;
        u0  = vargin.conversion.threshdesc.uncorr.thresh001;
    elseif isfield(vargin.conversion.threshdesc,'none')
        adjustment = -1;
        u0  = -Inf;
    end
    
    if isfield(vargin.conversion.cluster,'fwe2')
        extent_FWE = 1;
        pk  = vargin.conversion.cluster.fwe2.thresh05;
        noniso = vargin.conversion.cluster.fwe2.noniso;
    elseif isfield(vargin.conversion.cluster,'kuncorr')
        extent_FWE = 0;
        pk  = vargin.conversion.cluster.kuncorr.thresh05;
        noniso = vargin.conversion.cluster.kuncorr.noniso;
    elseif isfield(vargin.conversion.cluster,'k')
        extent_FWE = 0;
        pk  = vargin.conversion.cluster.k.kthresh;
        noniso = vargin.conversion.cluster.k.noniso;
    elseif isfield(vargin.conversion.cluster,'En')
        extent_FWE = 0;
        pk = -1;
        noniso = vargin.conversion.cluster.En.noniso;
    else
        extent_FWE = 0;
        pk = 0;
        noniso = 0;
    end
end

if nargin < 1
    P = spm_select(Inf,'^spmF.*(img|nii|gii)','Select F-images');

    sel = spm_input('Convert F value to?',1,'m',...
    '1-p|-log(1-p)|coefficient of determination R^2',1:3, 2);

    %-Get height threshold
    %-------------------------------------------------------------------
    str = 'FWE|FDR|uncorr|none';
    adjustment = spm_input('p value adjustment to control','+1','b',str,[1 2 0 -1],1);
    
    switch adjustment
    case 1 % family-wise false positive rate
        %---------------------------------------------------------------
        u0  = spm_input('p value (family-wise error)','+0','r',0.05,1,[0,1]);
    case 2 % False discovery rate
        %---------------------------------------------------------------    
        u0  = spm_input('p value (false discovery rate)','+0','r',0.05,1,[0,1]);
    case 0  %-NB: no adjustment
        % p for conjunctions is p of the conjunction SPM
        %---------------------------------------------------------------
        u0  = spm_input('threshold {F or p value}','+0','r',0.001,1);
    otherwise  %-NB: no threshold
        % p for conjunctions is p of the conjunction SPM
        %---------------------------------------------------------------
        u0  = -Inf;
    end

    if adjustment > -1 
        pk     = spm_input('extent threshold {k or p-value}','+1','r',0,1);
    else
        pk = 0;
    end
    if (pk < 1) && (pk > 0)
        extent_FWE = spm_input('p value (extent)','+1','b','uncorrected|FWE corrected',[0 1],1);
    end

    if pk ~= 0
        noniso = spm_input('Correct for non-isotropic smoothness?','+1','b','no|yes',[0 1],2);
    else
        noniso = 0;
    end
end


switch adjustment
case 1 % family-wise false positive rate
    p_height_str = '_pFWE';
case 2 % False discovery rate
    p_height_str = '_pFDR';
case 0 %-NB: no adjustment
    p_height_str = '_p';
otherwise  %-NB: no threshold
    p_height_str = '';
end

for i=1:size(P,1)
    [pth,nm] = spm_fileparts(deblank(P(i,:)));

    SPM_name = fullfile(pth, 'SPM.mat');
    
    % SPM.mat exist?
    if ~exist(SPM_name,'file')
       error('SPM.mat not found')
    end

    if strcmp(nm(1:6),'spmF_0') 
        Ic = str2double(nm(length(nm)-2:length(nm)));
    else
        error('Only spmF_0* files can be used');
    end

    load(SPM_name);
    xCon = SPM.xCon;
    df   = [xCon(Ic).eidf SPM.xX.erdf];
    STAT = xCon(Ic).STAT;
    R    = SPM.xVol.R;                  %-search Volume {resels}
    R    = R(1:find(R~=0,1,'last'));    % eliminate null resel counts
    S    = SPM.xVol.S;                  %-search Volume {voxels}
    XYZ  = SPM.xVol.XYZ;                %-XYZ coordinates
    FWHM = SPM.xVol.FWHM;
    v2r  = 1/prod(FWHM(~isinf(FWHM)));  %-voxels to resels

    Vspm   = cat(1,xCon(Ic).Vspm);
    Vspm.fname = fullfile(pth,Vspm.fname);

    if ~isfield(SPM.xVol,'VRpv')
        noniso = 0;
    end

    if noniso
        SPM.xVol.VRpv.fname = fullfile(pth,SPM.xVol.VRpv.fname);
    end

    switch adjustment
    case 1 % family-wise false positive rate
    %---------------------------------------------------------------
       u  = spm_uc(u0,df,STAT,R,1,S);

    case 2 % False discovery rate
    %---------------------------------------------------------------
       u  = spm_uc_FDR(u0,df,STAT,1,Vspm,0);

    otherwise  %-NB: no adjustment
    % p for conjunctions is p of the conjunction SPM
    %---------------------------------------------------------------
       if (u0 <= 1) && (u0 > 0)
           u = spm_u(u0,df,STAT);
       else
           u = u0;
       end
    end

    F = spm_data_read(Vspm.fname,'xyz',XYZ);

    %-Calculate height threshold filtering
    %-------------------------------------------------------------------    
    Q      = find(F > u);

    %-Apply height threshold
    %-------------------------------------------------------------------
    F      = F(:,Q);
    XYZ    = XYZ(:,Q);
    if isempty(Q)
        fprintf('No voxels survive height threshold u=%0.2g\n',u);
    end


    %-Extent threshold
    %-----------------------------------------------------------------------
    if ~isempty(XYZ)
        
       if (pk < 1) && (pk > 0)
            if extent_FWE
                Pk = 1;
                k = 0;
                while (Pk >= pk && k<S)
                    k = k + 1;
                    [Pk, Pn] = spm_P(1,k*v2r,u,df,STAT,R,1,S);
                end
                p_extent_str = ['_pkFWE' num2str(pk*100)];
            else
                Pn = 1;
                k = 0;
                while (Pn >= pk && k<S)
                    k = k + 1;
                    [Pk, Pn] = spm_P(1,k*v2r,u,df,STAT,R,1,S);
                end
                p_extent_str = ['_pk' num2str(pk*100)];
            end
        elseif (pk < 0)
            k = 0;
            [P2, Pn2, Em, En] = spm_P(1,k,u,df,STAT,R,1,S);
            k = ceil(En/v2r);
            p_extent_str = '_En';
        else
            k = pk;
            p_extent_str = '';
        end

        %-Calculate extent threshold filtering
        %-------------------------------------------------------------------
        if  ~spm_mesh_detect(xCon(Ic(1)).Vspm)
            A = spm_clusters(XYZ);
        else
            T = false(SPM.xVol.DIM');
            T(XYZ(1,:)) = true;
            G = export(gifti(SPM.xVol.G),'patch');
            A = spm_mesh_clusters(G,T)';
            A = A(XYZ(1,:));
        end

        if noniso
            fprintf('Use local RPV values to correct for non-stationary of smoothness.\n');

            Q     = [];
            if isfield(SPM.xVol,'G') % mesh detected?
                [N2,F2,XYZ2,A2,L2]  = spm_mesh_max(F,XYZ,gifti(SPM.xVol.G));
            else
                [N2,F2,XYZ2,A2,L2]  = spm_max(F,XYZ);
            end

            % sometimes max of A and A2 differ, thus we have to use the smaller value
            for i2 = 1:min([max(A) max(A2)])
                %-Get LKC for voxels in i-th region
                %----------------------------------------------------------
                LKC = spm_data_read(SPM.xVol.VRpv.fname,'xyz',L2{i2});
                
                %-Compute average of valid LKC measures for i-th region
                %----------------------------------------------------------
                valid = ~isnan(LKC);
                if any(valid)
                    LKC = sum(LKC(valid)) / sum(valid);
                else
                    LKC = v2r; % fall back to whole-brain resel density
                end
                
                %-Intrinsic volume (with surface correction)
                %----------------------------------------------------------
                IV   = spm_resels([1 1 1],L2{i2},'V');
                IV   = IV*[1/2 2/3 2/3 1]';
                k_noniso = IV*LKC/v2r;

                % find corresponding cluster in spm_clusters if cluster exceeds threshold
                if k_noniso >= k
                    ind2 = find(A2==i2);
                    for i = 1:min([max(A) max(A2)])
                        j = find(A == i);
                        if length(j)==N2(ind2)
                            if any(ismember(XYZ2(:,ind2)',XYZ(:,j)','rows'))
                                Q = [Q j];
                                break
                            end
                        end
                    end
                end
            end
        else
            Q     = [];
            for i = 1:min(max(A))
                j = find(A == i);
                if length(j) >= k; Q = [Q j]; end
            end
        end


        % ...eliminate voxels
        %-------------------------------------------------------------------
        F     = F(:,Q);
        XYZ   = XYZ(:,Q);
        if isempty(Q)
            fprintf('No voxels survived extent threshold k=%3.1f\n',k);
        end

    else

        k = 0;

    end % (if ~isempty(XYZ))

    if ~isempty(Q)
    
       switch sel
       case 1
          F2x = 1-spm_Fcdf(F,df);
          F2x_name = 'P_';
       case 2
          F2x = -log10(max(eps,1-spm_Fcdf(F,df)));
          F2x_name = 'logP_';
       case 3
    	  	F2x = (df(2)-1)*F./(df(2) - df(1)+F*(df(2) -1));
		      F2x_name = 'R2_';
       end

       str_num = deblank(xCon(Ic).name);

       % replace spaces with "_" and characters like "<" or ">" with "gt" or "lt"
       str_num(strfind(str_num,' ')) = '_';
       strpos = strfind(str_num,' > ');
       if ~isempty(strpos), str_num = [str_num(1:strpos-1) '_gt_' str_num(strpos+1:end)]; end
       strpos = strfind(str_num,' < ');
       if ~isempty(strpos), str_num = [str_num(1:strpos-1) '_lt_' str_num(strpos+1:end)]; end
       strpos = strfind(str_num,'>');
       if ~isempty(strpos), str_num = [str_num(1:strpos-1) 'gt' str_num(strpos+1:end)]; end
       strpos = strfind(str_num,'<');
       if ~isempty(strpos), str_num = [str_num(1:strpos-1) 'lt' str_num(strpos+1:end)]; end
       str_num = spm_str_manip(str_num,'v');
           
       if isfield(SPM.xVol,'G')
            ext = '.gii';
       else ext = '.nii'; end

       if u0 > -Inf
           name = [F2x_name str_num p_height_str num2str(u0*100) p_extent_str '_k' num2str(k) ext];
       else
           name = [F2x_name str_num ext];
       end
       fprintf('Save %s\n', name);
    
       out = deblank(fullfile(pth,name));

       %-Reconstruct (filtered) image from XYZ & F pointlist
       %-----------------------------------------------------------------------
       Y      = zeros(Vspm.dim);
       OFF    = XYZ(1,:) + Vspm.dim(1)*(XYZ(2,:)-1 + Vspm.dim(2)*(XYZ(3,:)-1));
       Y(OFF) = F2x;

       VO = Vspm;
       VO.fname = out;
       VO.dt = [spm_type('float32') spm_platform('bigend')];
       VO = spm_data_hdr_write(VO);
       spm_data_write(VO,Y);
    
    end
end

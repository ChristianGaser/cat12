function cg_spmT2x(vargin)
%CG_SPMT2X transformation of t-maps to P, -log(P), r or d-maps
%
% The following formulas are used:
%
% --------------------------------
% correlation coefficient:
% --------------------------------
%          sign(t)
% r = ------------------
%            df
%     sqrt(------ + 1)
%           t*t
%
% --------------------------------
% effect-size
% --------------------------------
%            2t
% d = ----------------
%         sqrt(df)
%
% --------------------------------
% p-value
% --------------------------------
%
% p = 1-spm_Tcdf
%
% --------------------------------
% log p-value
% --------------------------------
%
% -log10(1-P) = -log(1-spm_Tcdf)
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
% transform a large number of spmT-maps using the same thresholds.
%
% Naming convention of the transformed files:
%   Type_Contrast_Pheight_Pextent_K_Neg
%
%   Type:      P    - p-value
%              logP - log p-value
%              R    - correlation coefficient
%              D    - effect size
%              T    - t-value
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
%   Neg:       image also shows thresholded inverse effects (e.g. neg. 
%              values) 
%_______________________________________________________________________
% Christian Gaser
% $Id$

rev = '$Rev$';

if nargin == 1
    P = [];
    for i=1:numel(vargin.data)
        P = strvcat(P,deblank(vargin.data{i}));
    end
end

spm2 = 0;
if strcmp(spm('ver'),'SPM2'), spm2 = 1; end

if nargin < 1
    if spm2
        P = spm_get(Inf,'spmT*.img','Select normalized files');
    else
        P = spm_select(Inf,'^spmT.*\.img$','Select images');
    end
end

sel = spm_input('Convert t value to?',1,'m',...
    '1-p|-log(1-p)|correlation coefficient cc|effect size d|apply thresholds without conversion',1:5, 2);

%-Get height threshold
%-------------------------------------------------------------------
str = 'FWE|FDR|none';
adjustment = spm_input('p value adjustment to control','+1','b',str,[],1);

switch adjustment
case 'FWE' % family-wise false positive rate
    %---------------------------------------------------------------
    u0  = spm_input('p value (family-wise error)','+0','r',0.05,1,[0,1]);
    p_height_str = '_pFWE';

case 'FDR' % False discovery rate
    %---------------------------------------------------------------    
    u0  = spm_input('p value (false discovery rate)','+0','r',0.05,1,[0,1]);
    p_height_str = '_pFDR';

otherwise  %-NB: no adjustment
    % p for conjunctions is p of the conjunction SPM
    %---------------------------------------------------------------
    u0  = spm_input(['threshold {T or p value}'],'+0','r',0.001,1);
    p_height_str = '_p';

end

pk     = spm_input('extent threshold {k or p-value}','+1','r',0,1,[0,Inf]);
if (pk < 1) & (pk > 0)
    extent_FWE = spm_input('p value (extent)','+1','b','uncorrected|FWE corrected',[0 1],1);
end

if sel > 2
    neg_results = spm_input('Show also inverse effects (e.g. neg. values)','+1','b','yes|no',[1 0],2);
else
    neg_results = 0;
end

for i=1:size(P,1)
    spmT = deblank(P(i,:));
    Vspm = spm_vol(spmT);   
    [pth,nm,xt,vr] = fileparts(spmT);

    SPM_name = fullfile(pth, ['SPM.mat' vr]);
    
    % SPM.mat exist?
    if ~exist(SPM_name)
       error('SPM.mat not found')
    end

    if strcmp(nm(1:6),'spmT_0') 
        Ic = str2num(nm(length(nm)-2:length(nm)));
    else
        error('Only spmT_0* files can be used');
    end

    load(SPM_name);
    xCon = SPM.xCon;
    df   = [xCon(Ic).eidf SPM.xX.erdf];
    STAT = xCon(Ic).STAT;
    R    = SPM.xVol.R;          %-search Volume {resels}
    S    = SPM.xVol.S;          %-search Volume {voxels}
    XYZ  = SPM.xVol.XYZ;            %-XYZ coordinates
    FWHM = SPM.xVol.FWHM;
    v2r  = 1/prod(FWHM(~isinf(FWHM)));          %-voxels to resels

    if i==1 
        if isfield(SPM.xVol,'VRpv')
            noniso = spm_input('Correct for non-isotropic smoothness?','+1','b','no|yes',[0 1],2);
        else
            noniso = 0;
        end
    end

    if noniso
        SPM.xVol.VRpv = spm_vol(fullfile(pth,SPM.xVol.VRpv.fname));
    end

    switch adjustment
    case 'FWE' % family-wise false positive rate
    %---------------------------------------------------------------
       u  = spm_uc(u0,df,STAT,R,1,S);

    case 'FDR' % False discovery rate
    %---------------------------------------------------------------
       u  = spm_uc_FDR(u0,df,STAT,1,Vspm,0);

    otherwise  %-NB: no adjustment
    % p for conjunctions is p of the conjunction SPM
    %---------------------------------------------------------------
       if u0 <= 1
           u = spm_u(u0,df,STAT);
       else
           u = u0;
       end
    end

    Z = spm_get_data(Vspm,XYZ);

    %-Calculate height threshold filtering
    %-------------------------------------------------------------------    
    if neg_results
        Q      = find((Z > u) | (Z < -u));
    else
        Q      = find(Z > u);
    end

    %-Apply height threshold
    %-------------------------------------------------------------------
    Z      = Z(:,Q);
    XYZ    = XYZ(:,Q);
    if isempty(Q)
        fprintf('No voxels survive height threshold u=%0.2g\n',u);
    end


    %-Extent threshold
    %-----------------------------------------------------------------------
    if ~isempty(XYZ)

        if (pk < 1) & (pk > 0)
            if extent_FWE
                Pk = 1;
                k = 0;
                while (Pk >= pk & k<S)
                    k = k + 1;
                    [Pk Pn] = spm_P(1,k*v2r,u,df,STAT,R,1,S);
                end
                p_extent_str = ['_pkFWE' num2str(pk*100)];
            else
                Pn = 1;
                k = 0;
                while (Pn >= pk & k<S)
                    k = k + 1;
                    [Pk Pn] = spm_P(1,k*v2r,u,df,STAT,R,1,S);
                end
                p_extent_str = ['_pk' num2str(pk*100)];
            end
        else
            k = pk;
            p_extent_str = '';
        end
        
        %-Calculate extent threshold filtering
        %-------------------------------------------------------------------
        if noniso
            fprintf('Use local RPV values to correct for non-stationary of smoothness.\n');
            A     = spm_clusters(XYZ);
            [N Z2 XYZ2 A2 V2R] = spm_max_nS(Z,XYZ,SPM.xVol.VRpv);
            N = V2R/v2r;
            Q     = [];
            for i = 1:max(A)
                j = find(A == i);
                if N(min(find(A2 == i))) >= k; Q = [Q j]; end
            end
        else
            A     = spm_clusters(XYZ);
            Q     = [];
            for i = 1:max(A)
                j = find(A == i);
                if length(j) >= k; Q = [Q j]; end
            end
        end


        % ...eliminate voxels
        %-------------------------------------------------------------------
        Z     = Z(:,Q);
        XYZ   = XYZ(:,Q);
        if isempty(Q)
        fprintf('No voxels survived extent threshold k=%0.2g\n',k);
        end

    else

        k = 0;

    end % (if ~isempty(XYZ))

    if ~isempty(Q)
    
       switch sel
       case 1
          t2x = 1-spm_Tcdf(Z,df(2));
          t2x_name = 'P_';
       case 2
          t2x = -log10(max(eps,1-spm_Tcdf(Z,df(2))));
          t2x_name = 'logP_';
       case 3
          t2x = sign(Z).*(1./((df(2)./((Z.*Z)+eps))+1)).^0.5;
          t2x_name = 'R_';
       case 4
          t2x = 2*Z/sqrt(df(2));
          t2x_name = 'D_';
       case 5
          t2x = Z;
          t2x_name = 'T_';
       end

       str_num = deblank(xCon(Ic).name);

       % replace spaces with "_" and characters like "<" or ">" with "gt" or "lt"
       str_num(findstr(str_num,' ')) = '_';
       strpos = findstr(str_num,' > ');
       if ~isempty(strpos), str_num = [str_num(1:strpos-1) '_gt_' str_num(strpos+1:end)]; end
       strpos = findstr(str_num,' < ');
       if ~isempty(strpos), str_num = [str_num(1:strpos-1) '_lt_' str_num(strpos+1:end)]; end
       strpos = findstr(str_num,'>');
       if ~isempty(strpos), str_num = [str_num(1:strpos-1) 'gt' str_num(strpos+1:end)]; end
       strpos = findstr(str_num,'<');
       if ~isempty(strpos), str_num = [str_num(1:strpos-1) 'lt' str_num(strpos+1:end)]; end
       str_num = spm_str_manip(str_num,'v');
    
       if neg_results
            neg_str = '_bi'; 
       else
            neg_str = '';
       end
       
       name = [t2x_name str_num p_height_str num2str(u0*100) p_extent_str '_k' num2str(k) neg_str '.nii'];
       fprintf('Save %s\n', name);
    
       out = deblank(fullfile(pth,name));

       %-Reconstruct (filtered) image from XYZ & Z pointlist
       %-----------------------------------------------------------------------
       Y      = zeros(Vspm.dim(1:3));
       OFF    = XYZ(1,:) + Vspm.dim(1)*(XYZ(2,:)-1 + Vspm.dim(2)*(XYZ(3,:)-1));
       Y(OFF) = t2x;

       VO = Vspm;
       VO.fname = out;
       if spm2
            VO.dim(4) = spm_type('int16');
       else
            VO.dt = [spm_type('int16') spm_platform('bigend')];
       end
       spm_write_vol(VO,Y);
    
    end
end

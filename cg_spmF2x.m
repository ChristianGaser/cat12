function cg_spmF2x(vargin)
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
%   K:         extent threshold in voxels
%
%_______________________________________________________________________
% Christian Gaser
% $Id: cg_spmF2x.m 119 2009-04-24 13:10:00Z gaser $

rev = '$Rev: 119 $';

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
        P = spm_get(Inf,'spmF*.img','Select normalized files');
    else
        P = spm_select(Inf,'^spmF.*\.img$','Select images');
    end
end

sel = spm_input('Convert F value to?',1,'m',...
	'1-p|-log(1-p)|coefficient of determination R^2',[1 2 3], 2);

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
    u0  = spm_input(['threshold {F or p value}'],'+0','r',0.001,1);
    p_height_str = '_p';

end

k     = spm_input('extent threshold {k}','+1','r',0,1,[0,Inf]);

for i=1:size(P,1)
    spmF = deblank(P(i,:));
    Vspm = spm_vol(spmF);   
    [pth,nm,xt,vr] = fileparts(spmF);

    SPM_name = fullfile(pth, ['SPM.mat' vr]);
    
    % SPM.mat exist?
    if ~exist(SPM_name)
       error('SPM.mat not found')
    end

    if strcmp(nm(1:6),'spmF_0') 
        Ic = str2num(nm(length(nm)-2:length(nm)));
    else
        error('Only spmF_0* files can be used');
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

    F = spm_get_data(Vspm,XYZ);

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
        
        %-Calculate extent threshold filtering
        %-------------------------------------------------------------------
        A     = spm_clusters(XYZ);
        Q     = [];
        for i = 1:max(A)
            j = find(A == i);
            if length(j) >= k; Q = [Q j]; end
        end


        % ...eliminate voxels
        %-------------------------------------------------------------------
        F     = F(:,Q);
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
          F2x = 1-spm_Fcdf(F,df);
          F2x_name = 'P_';
       case 2
          F2x = -log10(max(eps,1-spm_Fcdf(F,df)));
          F2x_name = 'logP_';
       case 3
    	  	f2x = (df(2)-1)*F./(df(2) - df(1)+F*(df(2) -1));
		      F2x_name = 'R2';
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
           
       name = [F2x_name str_num p_height_str num2str(u0*100) '_k' num2str(k) '.nii'];
       fprintf('Save %s\n', name);
    
       out = deblank(fullfile(pth,name));

       %-Reconstruct (filtered) image from XYZ & F pointlist
       %-----------------------------------------------------------------------
       Y      = zeros(Vspm.dim(1:3));
       OFF    = XYZ(1,:) + Vspm.dim(1)*(XYZ(2,:)-1 + Vspm.dim(2)*(XYZ(3,:)-1));
       Y(OFF) = F2x;

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

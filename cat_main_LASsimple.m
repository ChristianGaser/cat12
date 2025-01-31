function Yml = cat_main_LASsimple(Ysrc,Ycls,Tth,LASstr)
%cat_main_LASsimple. Local Intensity Normalization.
% ______________________________________________________________________
%
% Highly simplified version of the Local Adapative Segmenation (LAS) 
% functions cat_main_LAS that only did minimal filtering of the 
% classification to avoid strong outliers such as blood vessels. 
%
% It is important to avoid high intensity blood vessels in the process, 
% because they will push down local WM and GM intensity. 
%
% Based on this values a intensity transformation is used. Compared to 
% the global correciton this has to be done for each voxel. To save time
% only a rough linear transformation is used.
% ______________________________________________________________________
%
%   Yml = cat_main_LASsimple(Ysrc,Ycls[,T3th,LASstr])
% 
%   Yml     .. local intensity correct image (T1w: 0-1 = BG-WM)
%   Ysrc    .. (bias corrected) T1 image
%   Ycls    .. SPM tissue class map
%   Tth     .. structure with tissue thresholds of CSF, GM, and WM in Ysrc
%   LASstr  .. parameter to control the strenght of the correction. 
%              (0 - slight, 1 - strong, default = 0.5) 
%
% This is an exclusive subfunction of cat_main. See also cat_main_LAS
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$ 

  % in case of given threshold resort them for Ycls order (GM,WM,CSF)
  if ~exist('Tth','var'), T3th = nan(1,3); else, T3th = Tth([2,3,1]); end 
  if ~exist('LASstr','var'), LASstr = .5; end
  
  %% Estimation of local tissue intensity for GM, WM, and CSF.
  %  In constrast to cat_main_LAS(s), the classification is only minimal 
  %  refined to avoid outliers from blood vessels.
  Ylab    = cell(1,3); 
  minYsrc = min(Ysrc(:))+1; 
  Ysrc    = Ysrc + minYsrc; 
  for ci = [1,2,3] % classes
    % tissue values
    Yi       = Ysrc .* (Ycls{ci}>128); 
    
    % estimate global threshold if not given
    if isnan(T3th(ci)) 
      T3th(ci) = cat_stat_nanmedian(Ysrc(Ycls{ci}>128));
    end

    % remove outliers
    Yi(Yi>T3th(ci) + 4*cat_stat_nanstd(Yi(Yi(:)>0))) = 0; 
    Yi(Yi>T3th(ci) + 4*cat_stat_nanstd(Yi(Yi(:)>0))) = 0; 
    
    % first approximation for local outlier removal 
    Yw       = cat_vol_approx(Yi); 
    Yi(Yi>Yw*1.2 | Yi<Yw/1.2) = 0; 
    
    % final approximation 
    Ylab{ci} = cat_vol_approx(Yi); 
    
    % smoothness defined by LASstr
    Ylab{ci} = cat_vol_smooth3X(Ylab{ci},4 * max(0,1-LASstr) ) - minYsrc; 
    
    % scaling 
    Ylab{ci} = Ylab{ci} / cat_stat_nanmedian(Ylab{ci}(Ycls{ci}>128)) * T3th(ci);
  end
  Ysrc = Ysrc - minYsrc; 

  %% order in case of T2; 
  [~,si] = sort(T3th); 
  Ylab   = Ylab(si([2,3,1])); % reorder 
 
  %% scaling and final correction of Yml similar to global map 
  Yml = zeros(size(Ysrc)); 
  Yml = Yml + ( (Ysrc>=Ylab{2}                ) .* (3 + (Ysrc - Ylab{2}) ./ max(eps,Ylab{2} - Ylab{3})) ); % scale highest tissue (WM in T1w)
  Yml = Yml + ( (Ysrc>=Ylab{1} & Ysrc<Ylab{2} ) .* (2 + (Ysrc - Ylab{1}) ./ max(eps,Ylab{2} - Ylab{1})) ); % scale second highest tissue (GM in T1w)
  Yml = Yml + ( (Ysrc>=Ylab{3} & Ysrc<Ylab{1} ) .* (1 + (Ysrc - Ylab{3}) ./ max(eps,Ylab{1} - Ylab{3})) ); % scale third highest tissue (CSF in T1w)
  Yml = Yml + ( (Ysrc< Ylab{3}                ) .* (    (Ysrc - minYsrc) ./ max(eps,Ylab{3} - minYsrc)) ); % scale background

  % rescale and general limits
  Yml = Yml / 3; 
  Yml(isnan(Yml) | Yml<0) = 0; Yml(Yml>2) = 2;
end 
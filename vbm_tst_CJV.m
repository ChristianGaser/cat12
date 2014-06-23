function CJV = vbm_tst_CJV(P,Pp0)
% ______________________________________________________________________
% Function to estimate the CJV (coefficent of joint variation) in
% images.
% 
%  CJV = vbm_tst_CJV(P,Pp0)
%  
%  P    ... set of n images (cellstr or char)
%  Pp0  ... set of 1 (ground thruth) or n images (cellstr of char)
%  CJV  ... matrix of n CJV values of each image
%
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
%
% $Id$
% ______________________________________________________________________

  if iscell(P)   && size(P,1)  <size(P,2), end% P=char(P);     end
  if iscell(Pp0) && size(Pp0,1)<size(P,2), end%Pp0=char(Pp0); end

  
  % get headers
  V   = spm_vol(char(P));
  Vp0 = spm_vol(char(Pp0));
  
  % get group gt-image
  if numel(Vp0)==1
    Yp0 = spm_read_vols(Vp0);
  end
  
  % get images
  CJV = zeros(size(V));
  if ~isempty(P)
    for i=1:numel(V)
      try
        Y = spm_read_vols(V(i));
        if numel(Vp0)>1
          Yp0 = spm_read_vols(Vp0(i));
        end

        CJV(i) = ( vbm_stat_nanstd(Y(Yp0>2.5))/2 + ...
                   vbm_stat_nanstd(Y(Yp0>1.5 & Yp0<2.5))/2 )  / ...
                 ( vbm_stat_nanmean(Y(Yp0>2.5)) - ...
                   vbm_stat_nanmean(Y(Yp0>1.5 & Yp0<2.5)) );
      catch
        vbm_io_cprintf('err',sprintf('Error: %s\n',P(i,:)));
        CJV(i) = nan;
      end
    end
  else
    vbm_io_cprintf('err',sprintf('Error: no images\n'));
  end
end
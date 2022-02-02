function [img,range] = cat_plot_glassbrain(P)
%cat_plot_glassbrain. Creat glassbrain images
%
% [img,range] = cat_plot_glassbrain(P)
%
% P     .. filename or image header
% img   .. glassbrain rendering (cell with 3 views + colormap)
% range .. minimum and maximum intensities
%
% Used in cat_long_report.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  if isstruct(P)
    V = P; 
    if isfield(V,'dat')
      Y = P.dat(:,:,:); 
    else 
      Y = spm_read_vols(V); 
    end    
  else
    V = spm_vol(P);
    Y = spm_read_vols(V); 
  end
 
  if 0
  % load futher atlas maps, eg. to outline the ventricles or the cerebellum
    LAB = cat_get_defaults('extopts.LAB');
    VT1 = spm_vol(char(cat_get_defaults('extopts.T1')));
    YT1 = spm_read_vols(VT1);
    YT1 = cat_vol_sample(VT1,V,YT1,1); 
    VA  = spm_vol(char(cat_get_defaults('extopts.cat12atlas')));
    YA  = spm_read_vols(VA); 

    YCT  = smooth3(YA==LAB.CT | YA==(LAB.CT+1)) & YT1>0.1;
    YV   = smooth3(YA==LAB.VT | YA==(LAB.VT+1));
  end
  
  %% main 
  img{1} = flip(rot90(cat_stat_nansum(Y,3),1),2); 
  img{2} = rot90(cat_stat_nansum(shiftdim(Y,2),3),2);
  img{3} = rot90(cat_stat_nansum(flip(shiftdim(Y,1),2),3),-1); 
  imin   = floor(min([img{1}(:);img{2}(:);img{3}(:)])); 
  imax   = ceil( max([img{1}(:);img{2}(:);img{3}(:)] ));
  img{4} = flip( (imin:(imax-imin)/128:imax)' , 1); 
  range  = [imin imax];
  
  if 0 
    %% test area 
    %d=1; figure; imagesc(log10(max(1,img{d}-2)) + img{d}), colormap pink; axis equal off; 
    %d=1; figure; imagesc(img{d}), colormap pink; axis equal off; 
    d=1; 
    figure(9939); clf(9939);  hold on; 
    image(img{d});
    contour(log(img{d}),[0.2 0.2],'color',[0 0 0]);
    if 0
      contour(log(imgct{d}),[0.5 0.5],'color',[0 0 0]);
      contour(log(imgv{d}),[0.2 0.2],'color',[0 0 1]);
    end
    colormap(flip(hot)); 
    axis equal off; hold off; 
  end
  

end

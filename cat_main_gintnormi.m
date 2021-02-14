function Ysrc = cat_main_gintnormi(Ym,Tth)
% cat_main_gintnormi
% ______________________________________________________________________
% 
% Inverse function of cat_main_gintnorm that restore the original 
% intensity levels of an intensity normalized image Ym. 
%
%   Ysrc = cat_main_gintnormi(Ym,Tth)
% 
%   Ysrc .. image with original intensity levels
%   Ym   .. intensity normalized image
%   Tth  .. data structure with the intensiy transformation
%
% Example:
%  % old values in the original image Ysrc
%  Tth.T3thx = [BGth CSFth GMth WMth WMth+diff([GMth,WMth]) ]; 
%  % new values of the output images of the tresholds Tth.T3thx
%  Tth.T3th  = 0:1/3:4/3;                                      
%  % div by 3 due to old definitions
%  Ym = cat_main_gintnormi(Ysrc/3,Tth); 
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id$


  % RD202102: There should be no inversion but in the test-retest of
  %           Buchert sanon-102748-00003-00001-1 it happend.
  Tth.T3th = sort(Tth.T3th); 
    
  if 1
    % RD2020: old version without interpolation is a bit inoptimal for the histogram 
    T3th  = Tth.T3thx; 
    T3thx = Tth.T3th; 
  else
    % use interpolation to avoid steps in the histogram
    % makima is better but not available in oder MATLABs :/
    %try 
    %  T3th  = interp1(Tth.T3thx,1:0.1:numel(Tth.T3thx),'makima');
    %  T3thx = interp1(Tth.T3th ,1:0.1:numel(Tth.T3th) ,'makima');
    %catch
      T3th  = interp1(Tth.T3thx,1:0.1:numel(Tth.T3thx),'spline');
      T3thx = interp1(Tth.T3th ,1:0.1:numel(Tth.T3th) ,'spline');
    %end 
  end
  
  if all(T3th==T3thx), Ysrc = Ym; return; end
  
  [T3th,si] = sort(T3th);
  T3thx = T3thx(si);
  
  isc=1;
  Ym = Ym*3;
  Ysrc = Ym; 
  for i=2:numel(T3th)
    M = Ym>T3th(i-1) & Ym<=T3th(i);
    Ysrc(M(:)) = T3thx(i-1) + (Ym(M(:)) - T3th(i-1))/diff(T3th(i-1:i))*diff(T3thx(i-1:i));
  end
  M  = Ym>=T3th(end); 
  Ysrc(M(:)) = numel(T3th)/isc/6 + (Ym(M(:)) - T3th(i))/diff(T3th(end-1:end))*diff(T3thx(i-1:i));    
return
function out = upfirdn2dMex(x, h, p, q)
%UPFIRDN2DMEX  Upsample, apply a specified FIR filter, and downsample a signal.
%   UPFIRDN2DMEX(X,H,P,Q) is a cascade of three systems applied to input signal X:
%         (1) Upsampling by P (zero insertion).  P defaults to 1 if not 
%             specified.
%         (2) FIR filtering with the filter specified by the impulse response 
%             given in H.
%         (3) Downsampling by Q (throwing away samples).  Q defaults to 1 if not 
%             specified.
%
% X must be a 2D matrix and h a vector. For that special case it is compatible to the
% Matlab function upfirdn.m from the Signal Processing Toolbox.
%
% Christian Gaser
% $Id: upfirdn2dMex.m 224 2009-12-02 23:39:15Z gaser $

rev = '$Rev: 224 $';

disp('Compiling upfirdn2dMex.c')

pth = fileparts(which(mfilename));
p_path = pwd;
cd(pth);
mex -O upfirdn2dMex.c 
cd(p_path);

out = upfirdn2dMex(x, h, p, q);

return

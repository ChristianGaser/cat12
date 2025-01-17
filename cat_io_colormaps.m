
function [C,XML] = cat_io_colormaps(Cname,ncolors)
% _________________________________________________________________________
% Create CAT Colormaps.
%
%   [C,XML] = cat_io_colormaps(Cname,ncolors)
%
%   Cname - colormap
%   colormaps internally used in CAT12
%   'marks'
%   'marks+'
%   'turbo'
%   'BCGWHw'
%   'BGWHn'
%   'magentadk'
%   'magenta'
%   'orange'
%   'blue'
%   'BCGWHw'
%   'BCGWHwov'
%   'BCGWHn'
%   'BCGWHn2'
%   'BCGWHgov'
%   'BCGWHnov'
%   'BCGWHcheckcov'
%   'curvature'
%   'hotinv'
%   'cold'
%   'coldinv'
%   'BWR'
%   categorical colormaps from RColorBrewer
%   https://rdrr.io/cran/RColorBrewer/man/ColorBrewer.html
%   'accent'
%   'dark2'
%   'paired'
%   'set1'
%   'set2'
%   'set3'
%   categorical colormaps from ggsci
%   https://nanx.me/ggsci/articles/ggsci.html
%   'nejm'
%   'jco'
%   'jama'
%   'd3'
% 
%   ncolors - number of colors
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$
  
  % number of colors:
  if ~exist('Cname','var')
    Cname = 'marks+';
  end
  if ~exist('ncolors','var')
    ncolors = [];
  else 
    if ncolors<1   
      error('MATLAB:cat_io_colormaps','Need at least one Color');
    elseif ncolors>2^12
      error('MATLAB:cat_io_colormaps', ...
       'Wow, why do you want to create a rainbow???. Please use less colors or change me.\n');
    end
  end
  
  % expect continuous colormaps by default
  cmap_categorical = false;
  
  % load basic colormap
  switch Cname
    case 'marks+', 
      C = [ 
            0.0000    0.4000    0.0000  % 0 - excellent (dark green)
            0.0000    0.8000    0.0000  % 1 - excellent (light green)
            0.4000    0.6000    0.1000  % 2 - good      (yellow-green)
            1.0000    0.6000    0.4000  % 3 - ok        (yellow-orange)
            1.0000    0.3000    0.0000  % 4 - bad       (red-orange)
            0.8000    0.2000    0.0000  % 5 - very bad  (red)
            0.7000    0.0000    0.0000  % 6 - unusable  (dark red)
            0.6000    0.0000    0.0000  % 7 - unusable  (dark red)
            0.5000    0.0000    0.0000  % 8 - unusable  (dark red)
            0.4000    0.0000    0.0000  % 9 - unusable  (dark red)
          ];
    case 'marks'
      C = [ %JET
            0.0000    0.0000    0.5625  % 4 - -3              (dark blue)
            0.0000    0.0000    1.0000  % 3 - -2              (blue)
            0.0000    1.0000    1.0000  % 2 - -1              (cyan)
            0.0000    1.0000    0.0000  % 1 -  0 normal case  (green)
            1.0000    1.0000    0.0000  % 2 - +1              (yellow)
            1.0000    0.0000    0.0000  % 3 - +2              (red)
            0.5104    0.0000    0.0000  % 4 - +3              (dark red)
          ];
    % vbm-output
    % GMT output
    % ...
    case 'trafficlight'
      C = [
        0.000    0.500   0.000 % dk green
        0.100    0.800   0.000 % green
        0.900    0.700   0.000 % yellow
        1.000    0.000   0.000 % red
        0.500    0.000   0.000 % dk red
      ];
    case 'trafficlight2'
      C = [
        0.000    0.000   0.900 % blue
        0.000    0.900   0.000 % green
        0.900    0.900   0.000 % orange
        1.000    0.000   0.000 % red
        0.700    0.000   0.700 % pink
      ];
    case 'accent'
      C = accent; 
      cmap_categorical = true;
    case 'dark2'
      C = dark2; 
      cmap_categorical = true;
    case 'paired'
      C = paired; 
      cmap_categorical = true;
    case 'set1'
      C = set1; 
      cmap_categorical = true;
    case 'set2'
      C = set2; 
      cmap_categorical = true;
    case 'set3'
      C = set3; 
      cmap_categorical = true;
    case 'nejm'
      C = nejm; 
      cmap_categorical = true;
    case 'jco'
      C = jco; 
      cmap_categorical = true;
    case 'jama'
      C = jama; 
      cmap_categorical = true;
    case 'd3'
      C = d3; 
      cmap_categorical = true;
    case 'turbo'
      C = turbo; 
    case 'magentadk'
      C = [0.95 0.95 0.95; 0.7 0.2 0.7];
    case 'magenta'
      C = [0.95 0.95 0.95; 1.0 0.4 1.0];
    case 'orange'
      C = [0.95 0.95 0.95; 0.8 0.4 0.6];
    case 'blue'
      C = blue;
    case 'BCGWHw'
      C = BCGWHw;
    case 'BCGWHwov'
      C = BCGWHwov;
    case 'BCGWHn'
      C = BCGWHn;
    case 'BCGWHn2';
      C = BCGWHnov;
    case 'BCGWHgov'
      C = BCGWHgov;
    case 'BCGWHnov'
      C = BCGWHnov;
    case 'BCGWHcheckcov'
      C = BCGWHcheckcov;
    case 'curvature';
      C = [ 
            0.9900    0.9900    0.9900 
            0.9500    0.9000    0.8000 
            0.9700    0.8500    0.6000 
            1.0000    0.8000    0.3000 
            1.0000    0.6000    0.0000 
            1.0000    0.3000    0.0000 
            1.0000    0.0000    0.0000  
            0.5000    0.0000    0.0000  
            0.0000    0.0000    0.0000  
          ];
    case 'hotinv';
      C = hotinv;
    case 'hot';
      C = hotinv; C = C(end:-1:1,:);
    case 'cold';
      C = hotinv; C = C(end:-1:1,:); C = [C(:,3),C(:,2),C(:,1)];
    case 'coldinv';
      C = hotinv; C = [C(:,3),C(:,2),C(:,1)];
    case 'BWR';
      CR = hotinv; 
      CB = [CR(:,3),CR(:,2),CR(:,1)]; CB = CB(end:-1:1,:);
      C  = [CB;CR(2:end,:,:,:)];
    otherwise, error('MATLAB:cat_io_colormaps','Unknown Colormap ''%s''\n',Cname);
  end
  if isempty(ncolors), ncolors = size(C,1); end
  
  % interpolate colormap, if colormap is categorical only use interpolation for larger number of colors
  if (size(C,1)<ncolors && ~cmap_categorical) || (size(C,1)<ncolors && cmap_categorical)
    ss    = (size(C,1)-1) / (ncolors-1);
    [X,Y] = meshgrid(1:ss:size(C,1),1:3);
    C     = interp2(1:size(C,1),1:3,C',X,Y)'; 
    XML   = cellstr([ dec2hex(round(min(255,max(0,C(:,1)*255)))), ...
             dec2hex(round(min(255,max(0,C(:,2)*255)))), ...
             dec2hex(round(min(255,max(0,C(:,3)*255)))) ]);
  elseif (size(C,1)>ncolors && ~cmap_categorical)
    ss    = (size(C,1)+1) / (ncolors);
    [X,Y] = meshgrid(1:ss:size(C,1)+1,1:3);  X = min(size(C,1),X);
    C     = interp2(1:size(C,1),1:3,C',X,Y)'; 
    XML   = cellstr([ dec2hex(round(min(255,max(0,C(:,1)*255)))), ...
             dec2hex(round(min(255,max(0,C(:,2)*255)))), ...
             dec2hex(round(min(255,max(0,C(:,3)*255)))) ]);
  % limit categorical colormaps
  elseif size(C,1)>ncolors && cmap_categorical
    C = C(1:ncolors,:);
  end
 
  
end


% from RColorBrewer
% https://rdrr.io/cran/RColorBrewer/man/ColorBrewer.html
function C=accent
  C = [
    127 201 127
    190 174 212
    253 192 134
    255 255 153
    56  108 176
    240 2   127
    191 91  23
    102 102 102
  ]/255;
end

% from RColorBrewer
% https://rdrr.io/cran/RColorBrewer/man/ColorBrewer.html
function C=dark2
  C = [
    27  158 119
    217 95  2
    117 112 179
    231 41  138
    102 155 30
    230 171 2
    166 118 29
    102 102 102 
  ]/255;
end

% from RColorBrewer
% https://rdrr.io/cran/RColorBrewer/man/ColorBrewer.html
function C=paired
  C = [
    166 206 227
    31  120 180
    178 223 138
    51  160 44
    251 154 153
    227 26  28
    253 191 111
    255 127 0
    202 178 214
    106 61  154
    255 255 153
    177 89  40
  ]/255;
end

% from RColorBrewer
% https://rdrr.io/cran/RColorBrewer/man/ColorBrewer.html
function C=set1
  C = [
    228 26  28
    55  126 184
    77  175 74
    152 78 163
    255 127 0
    255 255 51
    166 86  40
    247 129 191
    153 153 153
  ]/255;
end

% from RColorBrewer
% https://rdrr.io/cran/RColorBrewer/man/ColorBrewer.html
function C=set2
  C = [
    102 194 165
    252 141 98
    141 160 203
    231 138 195
    166 216 84
    255 217 47
    229 196 148
    179 179 179
  ]/255;
end

% from RColorBrewer
% https://rdrr.io/cran/RColorBrewer/man/ColorBrewer.html
function C=set3
  C = [
    141 211 199
    255 255 179
    190 186 218
    251 128 114
    128 177 211
    253 180 98
    179 222 105
    252 205 229
    217 217 217
    188 128 189
    204 235 197
    255 237 111
  ]/255;
end

% from ggsci
% https://nanx.me/ggsci/articles/ggsci.html
% extended by additional colors
function C=nejm
  C = [
    '#BC3C29'
    '#0072B5'
    '#E18727'
    '#20854E'
    '#7876B1'
    '#6F99AD'
    '#FFDC91'
    '#EE4C97'
    '#8C564B'
    '#BCBD22'
    '#00A1D5'
    '#374E55'
    '#003C67'
    '#8F7700'
    '#7F7F7F'    
    '#353535'    
  ];
  C = hex2rgb(C);
end

% from ggsci
% https://nanx.me/ggsci/articles/ggsci.html
function C=jama
  C = [
    '#374E55'
    '#DF8F44'
    '#00A1D5'
    '#B24745'
    '#79AF97'
    '#6A6599'
    '#80796B'
  ];
  C = hex2rgb(C);
end

% from ggsci
% https://nanx.me/ggsci/articles/ggsci.html
function C=jco
  C = [
    '#0073C2'
    '#EFC000'
    '#868686'
    '#CD534C'
    '#7AA6DC'
    '#003C67'
    '#8F7700'
    '#3B3B3B'
    '#A73030'
    '#4A6990'
  ];
  C = hex2rgb(C);
end

% from ggsci
% https://nanx.me/ggsci/articles/ggsci.html
function C=d3
  C = [
    '#1F77B4'
    '#FF7F0E'
    '#2CA02C'
    '#D62728'
    '#9467BD'
    '#8C564B'
    '#E377C2'
    '#7F7F7F'
    '#BCBD22'
    '#17BECF'
  ];
  C = hex2rgb(C);
end

function C=hotinv
  C = [ 
    0.9900    0.9900    0.9900 
    0.9500    0.9000    0.6000 
    1.0000    0.8000    0.3000 
    1.0000    0.6000    0.0000 
    1.0000    0.3000    0.0000 
    1.0000    0.0000    0.0000  
    0.5000    0.0000    0.0000  
    0.0000    0.0000    0.0000  
  ]; 
end

function C=BCGWHcheckcov
  C = [ 
         0         0         0
    0.0131    0.0281    0.0915
    0.0261    0.0562    0.1830
    0.0392    0.0843    0.2745
    0.0523    0.1124    0.3660
    0.0654    0.1405    0.4575
    0.0784    0.1686    0.5490
    0.1221    0.2437    0.6134
    0.1658    0.3188    0.6779
    0.2095    0.3938    0.7423
    0.2532    0.4689    0.8067
    0.2969    0.5440    0.8711
    0.3406    0.6190    0.9356
    0.3843    0.6941    1.0000
    0.3494    0.7219    0.9091
    0.3144    0.7497    0.8182
    0.2795    0.7775    0.7273
    0.2446    0.8053    0.6364
    0.2096    0.8332    0.5455
    0.1747    0.8610    0.4545
    0.1398    0.8888    0.3636
    0.1048    0.9166    0.2727
    0.0699    0.9444    0.1818
    0.0349    0.9722    0.0909
         0    1.0000         0
    0.1667    1.0000         0
    0.3333    1.0000         0
    0.5000    1.0000         0
    0.6667    1.0000         0
    0.8333    1.0000         0
    1.0000    1.0000         0
    1.0000    0.8333         0
    1.0000    0.6667         0
    1.0000    0.5000         0
    1.0000    0.3333         0
    1.0000    0.1667         0
    1.0000         0         0
    1.0000    0.0621    0.0719
    1.0000    0.1242    0.1438
    1.0000    0.1863    0.2157
    1.0000    0.2484    0.2876
    1.0000    0.3105    0.3595
    1.0000    0.3725    0.4314
    1.0000    0.4346    0.5033
    1.0000    0.4967    0.5752
    1.0000    0.5588    0.6471
    1.0000    0.6209    0.7190
    1.0000    0.6830    0.7908
    1.0000    0.7451    0.8627
    0.9663    0.7077    0.8424
    0.9325    0.6703    0.8220
    0.8988    0.6329    0.8016
    0.8651    0.5956    0.7812
    0.8314    0.5582    0.7608
    0.7976    0.5208    0.7404
    0.7639    0.4834    0.7200
    0.7302    0.4460    0.6996
    0.6965    0.4086    0.6792
    0.6627    0.3712    0.6588
    0.6290    0.3339    0.6384
    0.5953    0.2965    0.6180
    0.5616    0.2591    0.5976
    0.5278    0.2217    0.5773
    0.4941    0.1843    0.5569
  ]; 
end

function C=BCGWHgov
C = [
         0.95    0.95      0.95
         0.5    0.5        0.95
         0    0.5         1
         0    1         0.5
         0.5    1.0000         0
    0.4000    0.4000         0
    0.8000         0         0
    0.9000    0.4314    0.4627
    1.0000    0.8627    0.9255
    1.0000    0.4314    0.9627
    1.0000         0    1.0000
    ...0.7882         0    1.0000
    1              1    1.0000
    0.7882         0    1.0000
  ];
end

function C=BCGWHnov
C = [
         0         0         0
    ...0.0174    0.4980    0.7403
    ...0.8084    0.9216    1.0000
    ...0.6784    0.9216    1.0000
         0.0    0.05        .5
         0.0    0.4         1  % CSF
         0.0    0.7        0.1 %
         0    0.9500         0 % GM
    1.0000    1.0000         0 % 
    0.8000         0         0 % WM
    0.9000    0.4314    0.4627
    1.0000    0.8627    0.9255
    1.0000    0.4314    0.9627
    1.0000         1    1.0000
    0.7882         1    1.0000
    1              1    1.0000
  ];
end

function C=BCGWHn
  C = [
    0.0392    0.1412    0.4157
    0.0349    0.2366    0.4806
    0.0305    0.3320    0.5455
    0.0261    0.4275    0.6105
    0.0218    0.5229    0.6754
    0.0174    0.6183    0.7403
    0.0131    0.7137    0.8052
    0.0087    0.8092    0.8702
    0.0044    0.9046    0.9351
         0    1.0000    1.0000
         0    0.9163    0.8333
         0    0.8327    0.6667
         0    0.7490    0.5000
         0    0.6653    0.3333
         0    0.5817    0.1667
         0    0.4980         0
         0    0.5984         0
         0    0.6988         0
         0    0.7992         0
         0    0.8996         0
         0    1.0000         0
    0.3333    1.0000         0
    0.6667    1.0000         0
    1.0000    1.0000         0
    1.0000    0.8902         0
    1.0000    0.7804         0
    1.0000    0.6706         0
    1.0000    0.5608         0
    1.0000    0.4510         0
    0.9333    0.3007         0
    0.8667    0.1503         0
    0.8000         0         0
    0.8154    0.0462    0.0603
    0.8308    0.0923    0.1207
    0.8462    0.1385    0.1810
    0.8615    0.1846    0.2413
    0.8769    0.2308    0.3017
    0.8923    0.2769    0.3620
    0.9077    0.3231    0.4223
    0.9231    0.3692    0.4827
    0.9385    0.4154    0.5430
    0.9538    0.4615    0.6033
    0.9692    0.5077    0.6637
    0.9846    0.5538    0.7240
    1.0000    0.6000    0.7843
    0.9974    0.6214    0.7935
    0.9948    0.6429    0.8026
    0.9922    0.6643    0.8118
    0.9895    0.6858    0.8209
    0.9869    0.7072    0.8301
    0.9843    0.7286    0.8392
    0.9817    0.7501    0.8484
    0.9791    0.7715    0.8575
    0.9765    0.7929    0.8667
    0.9739    0.8144    0.8758
    0.9712    0.8358    0.8850
    0.9686    0.8573    0.8941
    0.9660    0.8787    0.9033
    0.9634    0.9001    0.9124
    0.9608    0.9216    0.9216
  ]; 
end

function C=BCGWHnov_old
  C = [
    0.0392    0.1412    0.4157
    0.0349    0.2366    0.4806
    0.0305    0.3320    0.5455
    0.0261    0.4275    0.6105
    0.0218    0.4980    0.6754
    0.0174    0.4980    0.7403
    0.0131    0.4980    0.8052
    0.0087    0.4980    0.8702
    0.0044    0.4980    0.9351
         0    0.4980    1.0000
         0    0.4980    0.8333
         0    0.4980    0.6667
         0    0.4980    0.5000
         0    0.4980    0.3333
         0    0.4980    0.1667
         0    0.4980         0
         0    0.6117         0
         0    0.7255         0
         0    0.8392         0
         0    0.9529         0
         0    1.0000         0
    0.2405    0.9608    0.0013
    0.4810    0.9216    0.0026
    0.7216    0.8824    0.0039
    0.8608    0.9412    0.0020
    1.0000    1.0000         0
    1.0000    0.8588         0
    1.0000    0.6549         0
    1.0000    0.4510         0
    0.9000    0.2255         0
    0.8000         0         0
    0.8200    0.0600    0.0784
    0.8400    0.1200    0.1569
    0.8600    0.1800    0.2353
    0.8800    0.2400    0.3137
    0.9000    0.3000    0.3922
    0.9200    0.3600    0.4706
    0.9400    0.4200    0.5490
    0.9600    0.4800    0.6274
    0.9800    0.5400    0.7059
    1.0000    0.6000    0.7843
    0.9749    0.5400    0.7808
    0.9498    0.4800    0.7772
    0.9247    0.4200    0.7737
    0.8996    0.3600    0.7702
    0.8745    0.3000    0.7667
    0.8494    0.2400    0.7631
    0.8243    0.1800    0.7596
    0.7992    0.1200    0.7561
    0.7741    0.0600    0.7525
    0.7490         0    0.7490
    0.7102         0    0.7102
    0.6714         0    0.6714
    0.6327         0    0.6327
    0.5939         0    0.5939
    0.5551         0    0.5551
    0.5163         0    0.5163
    0.4776         0    0.4776
    0.4388         0    0.4388
    0.4000         0    0.4000
         0         0         0
   ];
end

function C=BCGWHw
  C = [
    1.0000    1.0000    1.0000
    0.9741    0.9843    0.9961
    0.9482    0.9686    0.9922
    0.9224    0.9530    0.9882
    0.8965    0.9373    0.9843
    0.8706    0.9216    0.9804
    0.8322    0.9216    0.9843
    0.7937    0.9216    0.9882
    0.7553    0.9216    0.9922
    0.7168    0.9216    0.9961
    0.6784    0.9216    1.0000
    0.5686    0.8470    0.8882
    0.4588    0.7725    0.7765
    0.3059    0.6810    0.5177
    0.1529    0.5895    0.2588
         0    0.4980         0
         0    0.5984         0
         0    0.6988         0
         0    0.7992         0
         0    0.8996         0
         0    1.0000         0
    0.2500    1.0000         0
    0.5000    1.0000         0
    0.7500    1.0000         0
    1.0000    1.0000         0
    1.0000    0.8627         0
    1.0000    0.7255         0
    1.0000    0.5882         0
    1.0000    0.4510         0
    0.9333    0.3007         0
    0.8667    0.1503         0
    0.8000         0         0
    0.8500    0.1000    0.2000
    0.9000    0.2000    0.4000
    0.9500    0.3000    0.6000
    1.0000    0.4000    0.8000
    1.0000    0.4672    0.8286
    1.0000    0.5345    0.8571
    1.0000    0.6017    0.8857
    1.0000    0.6689    0.9143
    1.0000    0.7361    0.9429
    1.0000    0.8034    0.9714
    1.0000    0.8706    1.0000
    0.9500    0.7868    1.0000
    0.9000    0.7029    1.0000
    0.8500    0.6191    1.0000
    0.8000    0.5353    1.0000
    0.7500    0.4515    1.0000
    0.7000    0.3676    1.0000
    0.6500    0.2838    1.0000
    0.6000    0.2000    1.0000
    0.5250    0.1750    1.0000
    0.4500    0.1500    1.0000
    0.3750    0.1250    1.0000
    0.3000    0.1000    1.0000
    0.2250    0.0750    1.0000
    0.1500    0.0500    1.0000
    0.0750    0.0250    1.0000
         0         0    1.0000
         0         0         0
  ]; 
end

function C=BCGWHwov
  C = [
    1.0000    1.0000    1.0000
    0.9741    0.9843    0.9961
    0.9482    0.9686    0.9922
    0.9224    0.9530    0.9882
    0.8965    0.9373    0.9843
    0.8706    0.9216    0.9804
    0.8225    0.9216    0.9853
    0.7745    0.9216    0.9902
    0.7264    0.9216    0.9951
    0.6784    0.9216    1.0000
    0.6052    0.8719    0.9255
    0.5320    0.8222    0.8510
    0.4588    0.7725    0.7765
    0.2294    0.6352    0.3882
         0    0.4980         0
         0    0.6117         0
         0    0.7255         0
         0    0.8392         0
         0    0.9529         0
         0    1.0000         0
    0.2405    0.9608    0.0013
    0.4810    0.9216    0.0026
    0.7216    0.8824    0.0039
    0.8608    0.9412    0.0020
    1.0000    1.0000         0
    1.0000    0.8588         0
    1.0000    0.6549         0
    1.0000    0.4510         0
    0.9000    0.2255         0
    0.8000         0         0
    0.8200    0.0600    0.0784
    0.8400    0.1200    0.1569
    0.8600    0.1800    0.2353
    0.8800    0.2400    0.3137
    0.9000    0.3000    0.3922
    0.9200    0.3600    0.4706
    0.9400    0.4200    0.5490
    0.9600    0.4800    0.6274
    0.9800    0.5400    0.7059
    1.0000    0.6000    0.7843
    0.9749    0.5400    0.7808
    0.9498    0.4800    0.7772
    0.9247    0.4200    0.7737
    0.8996    0.3600    0.7702
    0.8745    0.3000    0.7667
    0.8494    0.2400    0.7631
    0.8243    0.1800    0.7596
    0.7992    0.1200    0.7561
    0.7741    0.0600    0.7525
    0.7490         0    0.7490
    0.7102         0    0.7102
    0.6714         0    0.6714
    0.6327         0    0.6327
    0.5939         0    0.5939
    0.5551         0    0.5551
    0.5163         0    0.5163
    0.4776         0    0.4776
    0.4388         0    0.4388
    0.4000         0    0.4000
         0         0         0
   ];
end

function C = blue
  C = [
    0.02  0.25  0.50
    0.20  0.45  0.75
    0.40  0.80  0.95
    0.80  0.95  0.98
    0.94  0.96  0.99
    ];
end

function C = turbo
%TURBO   Turbo colormap.
%   TURBO(M) returns an M-by-3 matrix containing the turbo colormap, a
%   variant of the jet colormap that is more perceptually uniform.
C = [
     0.18995, 0.07176, 0.23217;
     0.19483, 0.08339, 0.26149;
     0.19956, 0.09498, 0.29024;
     0.20415, 0.10652, 0.31844;
     0.20860, 0.11802, 0.34607;
     0.21291, 0.12947, 0.37314;
     0.21708, 0.14087, 0.39964;
     0.22111, 0.15223, 0.42558;
     0.22500, 0.16354, 0.45096;
     0.22875, 0.17481, 0.47578;
     0.23236, 0.18603, 0.50004;
     0.23582, 0.19720, 0.52373;
     0.23915, 0.20833, 0.54686;
     0.24234, 0.21941, 0.56942;
     0.24539, 0.23044, 0.59142;
     0.24830, 0.24143, 0.61286;
     0.25107, 0.25237, 0.63374;
     0.25369, 0.26327, 0.65406;
     0.25618, 0.27412, 0.67381;
     0.25853, 0.28492, 0.69300;
     0.26074, 0.29568, 0.71162;
     0.26280, 0.30639, 0.72968;
     0.26473, 0.31706, 0.74718;
     0.26652, 0.32768, 0.76412;
     0.26816, 0.33825, 0.78050;
     0.26967, 0.34878, 0.79631;
     0.27103, 0.35926, 0.81156;
     0.27226, 0.36970, 0.82624;
     0.27334, 0.38008, 0.84037;
     0.27429, 0.39043, 0.85393;
     0.27509, 0.40072, 0.86692;
     0.27576, 0.41097, 0.87936;
     0.27628, 0.42118, 0.89123;
     0.27667, 0.43134, 0.90254;
     0.27691, 0.44145, 0.91328;
     0.27701, 0.45152, 0.92347;
     0.27698, 0.46153, 0.93309;
     0.27680, 0.47151, 0.94214;
     0.27648, 0.48144, 0.95064;
     0.27603, 0.49132, 0.95857;
     0.27543, 0.50115, 0.96594;
     0.27469, 0.51094, 0.97275;
     0.27381, 0.52069, 0.97899;
     0.27273, 0.53040, 0.98461;
     0.27106, 0.54015, 0.98930;
     0.26878, 0.54995, 0.99303;
     0.26592, 0.55979, 0.99583;
     0.26252, 0.56967, 0.99773;
     0.25862, 0.57958, 0.99876;
     0.25425, 0.58950, 0.99896;
     0.24946, 0.59943, 0.99835;
     0.24427, 0.60937, 0.99697;
     0.23874, 0.61931, 0.99485;
     0.23288, 0.62923, 0.99202;
     0.22676, 0.63913, 0.98851;
     0.22039, 0.64901, 0.98436;
     0.21382, 0.65886, 0.97959;
     0.20708, 0.66866, 0.97423;
     0.20021, 0.67842, 0.96833;
     0.19326, 0.68812, 0.96190;
     0.18625, 0.69775, 0.95498;
     0.17923, 0.70732, 0.94761;
     0.17223, 0.71680, 0.93981;
     0.16529, 0.72620, 0.93161;
     0.15844, 0.73551, 0.92305;
     0.15173, 0.74472, 0.91416;
     0.14519, 0.75381, 0.90496;
     0.13886, 0.76279, 0.89550;
     0.13278, 0.77165, 0.88580;
     0.12698, 0.78037, 0.87590;
     0.12151, 0.78896, 0.86581;
     0.11639, 0.79740, 0.85559;
     0.11167, 0.80569, 0.84525;
     0.10738, 0.81381, 0.83484;
     0.10357, 0.82177, 0.82437;
     0.10026, 0.82955, 0.81389;
     0.09750, 0.83714, 0.80342;
     0.09532, 0.84455, 0.79299;
     0.09377, 0.85175, 0.78264;
     0.09287, 0.85875, 0.77240;
     0.09267, 0.86554, 0.76230;
     0.09320, 0.87211, 0.75237;
     0.09451, 0.87844, 0.74265;
     0.09662, 0.88454, 0.73316;
     0.09958, 0.89040, 0.72393;
     0.10342, 0.89600, 0.71500;
     0.10815, 0.90142, 0.70599;
     0.11374, 0.90673, 0.69651;
     0.12014, 0.91193, 0.68660;
     0.12733, 0.91701, 0.67627;
     0.13526, 0.92197, 0.66556;
     0.14391, 0.92680, 0.65448;
     0.15323, 0.93151, 0.64308;
     0.16319, 0.93609, 0.63137;
     0.17377, 0.94053, 0.61938;
     0.18491, 0.94484, 0.60713;
     0.19659, 0.94901, 0.59466;
     0.20877, 0.95304, 0.58199;
     0.22142, 0.95692, 0.56914;
     0.23449, 0.96065, 0.55614;
     0.24797, 0.96423, 0.54303;
     0.26180, 0.96765, 0.52981;
     0.27597, 0.97092, 0.51653;
     0.29042, 0.97403, 0.50321;
     0.30513, 0.97697, 0.48987;
     0.32006, 0.97974, 0.47654;
     0.33517, 0.98234, 0.46325;
     0.35043, 0.98477, 0.45002;
     0.36581, 0.98702, 0.43688;
     0.38127, 0.98909, 0.42386;
     0.39678, 0.99098, 0.41098;
     0.41229, 0.99268, 0.39826;
     0.42778, 0.99419, 0.38575;
     0.44321, 0.99551, 0.37345;
     0.45854, 0.99663, 0.36140;
     0.47375, 0.99755, 0.34963;
     0.48879, 0.99828, 0.33816;
     0.50362, 0.99879, 0.32701;
     0.51822, 0.99910, 0.31622;
     0.53255, 0.99919, 0.30581;
     0.54658, 0.99907, 0.29581;
     0.56026, 0.99873, 0.28623;
     0.57357, 0.99817, 0.27712;
     0.58646, 0.99739, 0.26849;
     0.59891, 0.99638, 0.26038;
     0.61088, 0.99514, 0.25280;
     0.62233, 0.99366, 0.24579;
     0.63323, 0.99195, 0.23937;
     0.64362, 0.98999, 0.23356;
     0.65394, 0.98775, 0.22835;
     0.66428, 0.98524, 0.22370;
     0.67462, 0.98246, 0.21960;
     0.68494, 0.97941, 0.21602;
     0.69525, 0.97610, 0.21294;
     0.70553, 0.97255, 0.21032;
     0.71577, 0.96875, 0.20815;
     0.72596, 0.96470, 0.20640;
     0.73610, 0.96043, 0.20504;
     0.74617, 0.95593, 0.20406;
     0.75617, 0.95121, 0.20343;
     0.76608, 0.94627, 0.20311;
     0.77591, 0.94113, 0.20310;
     0.78563, 0.93579, 0.20336;
     0.79524, 0.93025, 0.20386;
     0.80473, 0.92452, 0.20459;
     0.81410, 0.91861, 0.20552;
     0.82333, 0.91253, 0.20663;
     0.83241, 0.90627, 0.20788;
     0.84133, 0.89986, 0.20926;
     0.85010, 0.89328, 0.21074;
     0.85868, 0.88655, 0.21230;
     0.86709, 0.87968, 0.21391;
     0.87530, 0.87267, 0.21555;
     0.88331, 0.86553, 0.21719;
     0.89112, 0.85826, 0.21880;
     0.89870, 0.85087, 0.22038;
     0.90605, 0.84337, 0.22188;
     0.91317, 0.83576, 0.22328;
     0.92004, 0.82806, 0.22456;
     0.92666, 0.82025, 0.22570;
     0.93301, 0.81236, 0.22667;
     0.93909, 0.80439, 0.22744;
     0.94489, 0.79634, 0.22800;
     0.95039, 0.78823, 0.22831;
     0.95560, 0.78005, 0.22836;
     0.96049, 0.77181, 0.22811;
     0.96507, 0.76352, 0.22754;
     0.96931, 0.75519, 0.22663;
     0.97323, 0.74682, 0.22536;
     0.97679, 0.73842, 0.22369;
     0.98000, 0.73000, 0.22161;
     0.98289, 0.72140, 0.21918;
     0.98549, 0.71250, 0.21650;
     0.98781, 0.70330, 0.21358;
     0.98986, 0.69382, 0.21043;
     0.99163, 0.68408, 0.20706;
     0.99314, 0.67408, 0.20348;
     0.99438, 0.66386, 0.19971;
     0.99535, 0.65341, 0.19577;
     0.99607, 0.64277, 0.19165;
     0.99654, 0.63193, 0.18738;
     0.99675, 0.62093, 0.18297;
     0.99672, 0.60977, 0.17842;
     0.99644, 0.59846, 0.17376;
     0.99593, 0.58703, 0.16899;
     0.99517, 0.57549, 0.16412;
     0.99419, 0.56386, 0.15918;
     0.99297, 0.55214, 0.15417;
     0.99153, 0.54036, 0.14910;
     0.98987, 0.52854, 0.14398;
     0.98799, 0.51667, 0.13883;
     0.98590, 0.50479, 0.13367;
     0.98360, 0.49291, 0.12849;
     0.98108, 0.48104, 0.12332;
     0.97837, 0.46920, 0.11817;
     0.97545, 0.45740, 0.11305;
     0.97234, 0.44565, 0.10797;
     0.96904, 0.43399, 0.10294;
     0.96555, 0.42241, 0.09798;
     0.96187, 0.41093, 0.09310;
     0.95801, 0.39958, 0.08831;
     0.95398, 0.38836, 0.08362;
     0.94977, 0.37729, 0.07905;
     0.94538, 0.36638, 0.07461;
     0.94084, 0.35566, 0.07031;
     0.93612, 0.34513, 0.06616;
     0.93125, 0.33482, 0.06218;
     0.92623, 0.32473, 0.05837;
     0.92105, 0.31489, 0.05475;
     0.91572, 0.30530, 0.05134;
     0.91024, 0.29599, 0.04814;
     0.90463, 0.28696, 0.04516;
     0.89888, 0.27824, 0.04243;
     0.89298, 0.26981, 0.03993;
     0.88691, 0.26152, 0.03753;
     0.88066, 0.25334, 0.03521;
     0.87422, 0.24526, 0.03297;
     0.86760, 0.23730, 0.03082;
     0.86079, 0.22945, 0.02875;
     0.85380, 0.22170, 0.02677;
     0.84662, 0.21407, 0.02487;
     0.83926, 0.20654, 0.02305;
     0.83172, 0.19912, 0.02131;
     0.82399, 0.19182, 0.01966;
     0.81608, 0.18462, 0.01809;
     0.80799, 0.17753, 0.01660;
     0.79971, 0.17055, 0.01520;
     0.79125, 0.16368, 0.01387;
     0.78260, 0.15693, 0.01264;
     0.77377, 0.15028, 0.01148;
     0.76476, 0.14374, 0.01041;
     0.75556, 0.13731, 0.00942;
     0.74617, 0.13098, 0.00851;
     0.73661, 0.12477, 0.00769;
     0.72686, 0.11867, 0.00695;
     0.71692, 0.11268, 0.00629;
     0.70680, 0.10680, 0.00571;
     0.69650, 0.10102, 0.00522;
     0.68602, 0.09536, 0.00481;
     0.67535, 0.08980, 0.00449;
     0.66449, 0.08436, 0.00424;
     0.65345, 0.07902, 0.00408;
     0.64223, 0.07380, 0.00401;
     0.63082, 0.06868, 0.00401;
     0.61923, 0.06367, 0.00410;
     0.60746, 0.05878, 0.00427;
     0.59550, 0.05399, 0.00453;
     0.58336, 0.04931, 0.00486;
     0.57103, 0.04474, 0.00529;
     0.55852, 0.04028, 0.00579;
     0.54583, 0.03593, 0.00638;
     0.53295, 0.03169, 0.00705;
     0.51989, 0.02756, 0.00780;
     0.50664, 0.02354, 0.00863;
     0.49321, 0.01963, 0.00955;
     0.47960, 0.01583, 0.01055];
end

function rgb = hex2rgb(hex)
  rgb = reshape(sscanf(hex(:,2:end)','%2x'),3,[]).'/255;
end
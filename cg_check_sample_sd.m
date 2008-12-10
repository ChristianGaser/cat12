function cg_check_sample_sd(vargin)
%cg_check_sample_sd	to check standard deviation across sample
%
% Images have to be in the same orientation with same voxel size
% and dimension (e.g. normalized images)
% An output image will be save with SD at each voxel.
%_______________________________________________________________________
% Christian Gaser
% $Id$

if nargin == 1
	P = [];
	for i=1:numel(vargin.data)
		P = strvcat(P,deblank(vargin.data{i}));
	end
	norm = vargin.scale;
	sd_name = vargin.sd_name;
	mean_name = vargin.mean_name;
end

if nargin < 1
	P = spm_select(Inf,'image','Select images');
end

V = spm_vol(P);
n = size(P,1);

if length(V)>1 & any(any(diff(cat(1,V.dim),1,1),1))
	error('images don''t all have same dimensions'), end
if any(any(any(diff(cat(3,V.mat),1,3),3)))
	error('images don''t all have same orientation & voxel size'), end

if nargin < 1
	norm = spm_input('Prop. scaling (e.g. for T1- or modulated images)?',1,'yes|no',[1 0],2);
	sd_name = spm_input('Output SD filename',1,'s','SD.img');
	mean_name = spm_input('Output Mean filename',1,'s','Mean.img');
end

if ~isempty(sd_name)
	Vo1 = V(1);
    Vo1.dt = [spm_type('int16') spm_platform('bigend')];
	Vo1.fname = sd_name;

	Vo1.descrip = 'std';
end

if ~isempty(mean_name)
	Vo2 = V(1);
	Vo2.dt = [spm_type('int16') spm_platform('bigend')];
	Vo2.fname = mean_name;

	Vo.descrip = 'mean';
end


% global scaling
if norm
	gm=zeros(size(V,1),1);
	disp('Calculating globals...');
	for i=1:size(V,1), gm(i) = spm_global(V(i)); end
	gm_all = mean(gm);
	for i=1:n
		V(i).pinfo(1:2,:) = gm_all*V(i).pinfo(1:2,:)/gm(i);
	end
end


%-Start progress plot
%-----------------------------------------------------------------------
mean_data = zeros(V(1).dim(1:3));
std_data  = zeros(V(1).dim(1:3));
squared_distance = zeros(n,1);

spm_progress_bar('Init',2*n,'volumes completed');
for i=1:n
    vol = spm_read_vols(V(i));
    mean_data = mean_data + vol/n;
	spm_progress_bar('Set',i);    
end

mask = (mean_data ~= 0);
sz_mask = sum(mask(:));

for i=1:n
    vol = spm_read_vols(V(i));
    dist2 = (mean_data(mask) - vol(mask)).^2;
    squared_distance(i) = sum(dist2(~isnan(dist2)));
    std_data(mask) = std_data(mask) + dist2;
	spm_progress_bar('Set',i+n);    
end

std_data = (std_data/(n-1)).^0.5;

spm_progress_bar('Clear');

% sort files
fprintf('Squared distance to mean:\n');
[squared_distance_sorted, ind] = sort(squared_distance);
for i=1:n
    fprintf('%s: %3.2f\n',V(ind(i)).fname,squared_distance_sorted(i));
end

Fgraph = spm_figure('GetWin','Graphics');
spm_figure('Clear',Fgraph);
FS      = spm('FontSizes');

[tmp fname] = spm_str_manip(char(V.fname),'C');
fprintf('Compressed filenames: %s    \n',tmp);

xpos = 2*(0:n-1)/(n-1);
for i=1:n
    text(xpos(i),squared_distance(i),fname.m{i},'FontSize',FS(8),'HorizontalAlignment','center')
end

hold on
cg_boxplot({squared_distance});
set(gca,'XTick',[],'XLim',[-.25 2.25]);
if max(squared_distance) > min(squared_distance)
  set(gca,'YLim',[0.9*min(squared_distance) 1.1*max(squared_distance)]);
end
title(sprintf('Boxplot: squared distance to mean    \nCommon filename: %s*%s',fname.s,fname.e),'FontSize',FS(12),'FontWeight','Bold');
ylabel('<----- low (good quality) --- squared distance from mean --- large (poor quality)------>    ','FontSize',FS(10),'FontWeight','Bold');
xlabel('<----- first --- file order --- last ------>    ','FontSize',FS(10),'FontWeight','Bold');
hold off

% check for replicates
for i=1:n
  for j=1:n
	if (i>j) & (squared_distance(i) == squared_distance(j))
		[s,differ] = unix(['diff ' V(i).fname ' ' V(j).fname]);
		if (s==0), fprintf(['\nWarning: ' V(i).fname ' and ' V(j).fname ' are same files?\n']); end
	end
  end
end

%-Write output image (uses spm_write_vol - which calls spm_write_plane)
%-----------------------------------------------------------------------
if ~isempty(sd_name)
	Vo1 = spm_write_vol(Vo1,std_data);
	disp(['Writing ' Vo1.fname]);
end

if ~isempty(mean_name)
	Vo2 = spm_write_vol(Vo2,mean_data);
	disp(['Writing ' Vo2.fname]);
end

%-End
%-----------------------------------------------------------------------
spm_progress_bar('Clear')

show = spm_input('Show files with CheckReg (beginning from bad)?','+1','yes|no',[1 0],2);
if show
	number = min([n 15]);
	number = spm_input('How many files ?','+1','e',number);
	
	list = str2mat(V(ind(n:-1:1)).fname);
	list2 = list(1:number,:);
	spm_check_registration(list2)
end
return

%-----------------------------------------------------------------------
function s = cg_boxplot (data,notched,symbol,vertical,maxwhisker)
%-----------------------------------------------------------------------
% usage: s = cg_boxplot (data,notched,symbol,vertical,maxwhisker);
%
% The box plot is a graphical display that simultaneously describes several 
% important features of a data set, such as center, spread, departure from 
% symmetry, and identification of observations that lie unusually far from
% the bulk of the data.
%
% data is a matrix with one column for each dataset, or data is a cell
% vector with one cell for each dataset.
% notched = 1 produces a notched-box plot. Notches represent a robust 
% estimate of the uncertainty about the median.
% notched = 0 (default) produces a rectangular box plot. 
% notched in (0,1) produces a notch of the specified depth.
% notched values outside [0,1] are amusing if not exactly practical.
% symbol sets the symbol for the outlier values, default symbol for
% points that lie outside 3 times the interquartile range is 'o',
% default symbol for points between 1.5 and 3 times the interquartile
% range is '+'. 
% Examples
% symbol = '.' points between 1.5 and 3 times the IQR is marked with
% '.' and points outside 3 times IQR with 'o'.
% symbol = ['x','*'] points between 1.5 and 3 times the IQR is marked with
% 'x' and points outside 3 times IQR with '*'.
% vertical = 0 makes the boxes horizontal, by default vertical = 1.
% maxwhisker defines the length of the whiskers as a function of the IQR
% (default = 1.5). If maxwhisker = 0 then boxplot displays all data  
% values outside the box using the plotting symbol for points that lie
% outside 3 times the IQR.   
%
% The returned matrix s has one column for each dataset as follows:
%
%    1  minimum
%    2  1st quartile
%    3  2nd quartile (median)
%    4  3rd quartile
%    5  maximum
%    6  lower confidence limit for median
%    7  upper confidence limit for median
%
% Example
%
%   title("Grade 3 heights");
%   tics("x",1:2,["girls";"boys"]);
%   axis([0,3]);
%   boxplot({randn(10,1)*5+140, randn(13,1)*8+135});
%

% Author: Alberto Terruzzi <t-albert@libero.it>
% Version: 1.4
% Created: 6 January 2002
% Copyright (C) 2002 Alberto Terruzzi
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% modified by Christian Gaser (christian.gaser@uni-jena.de)
% original version was written for octave by Alberto Terruzzi

% assign parameter defaults
if nargin < 1 || nargin > 5
   error('s = boxplot (data,notch,symbol,vertical,maxwhisker)')
end
if nargin < 5, maxwhisker = 1.5; end
if nargin < 4, vertical = 1; end
if nargin < 3, symbol = ['+','o']; end
if nargin < 2, notched = 0; end

if length(symbol)==1, symbol(2)=symbol(1); end

if notched==1, notched=0.25; end
a = 1-notched;

% figure out how many data sets we have
if iscell(data), 
  nc = length(data);
else
  if isvector(data), data = data(:); end
  nc = columns(data);
end

% compute statistics
% s will contain
%    1,5    min and max
%    2,3,4  1st, 2nd and 3rd quartile
%    6,7    lower and upper confidence intervals for median
s = zeros(7,nc);
box = zeros(1,nc);
whisker_x = ones(2,1)*[1:nc,1:nc];
whisker_y = zeros(2,2*nc);
outliers_x = [];
outliers_y = [];
outliers2_x = [];
outliers2_y = [];

for i=1:nc
  % Get the next data set from the array or cell array
  if iscell(data)
    col = data{i}(:);
  else
    col = data(:,i);
  end
  % Skip missing data
  col(isnan(col)) = [];
  % Remember the data length
  nd = length(col);
  box(i) = nd;
  if (nd > 1)
    % min,max and quartiles
    s(1:5,i) = [min(col) prctile(col,[25 50 75]) max(col)]';
    % confidence interval for the median
    est = 1.57*(s(4,i)-s(2,i))/sqrt(nd);
    s(6,i) = max([s(3,i)-est, s(2,i)]);
    s(7,i) = min([s(3,i)+est, s(4,i)]);
    % whiskers out to the last point within the desired inter-quartile range
    IQR = maxwhisker*(s(4,i)-s(2,i));
    whisker_y(:,i) = [min(col(col >= s(2,i)-IQR)); s(2,i)];
    whisker_y(:,nc+i) = [max(col(col <= s(4,i)+IQR)); s(4,i)];
    % outliers beyond 1 and 2 inter-quartile ranges
    outliers = col((col < s(2,i)-IQR & col >= s(2,i)-2*IQR) | (col > s(4,i)+IQR & col <= s(4,i)+2*IQR));
    outliers2 = col(col < s(2,i)-2*IQR | col > s(4,i)+2*IQR);
    outliers_x = [outliers_x; i*ones(size(outliers))];
    outliers_y = [outliers_y; outliers];
    outliers2_x = [outliers2_x; i*ones(size(outliers2))];
    outliers2_y = [outliers2_y; outliers2];
  elseif (nd == 1)
    % all statistics collapse to the value of the point
    s(:,i) = col;
    % single point data sets are plotted as outliers.
    outliers_x = [outliers_x; i];
    outliers_y = [outliers_y; col];
  else
    % no statistics if no points
    s(:,i) = NaN;
  end
end

% Note which boxes don't have enough stats
chop = find(box <= 1);
    
% Draw a box around the quartiles, with width proportional to the number of
% items in the box. Draw notches if desired.
box = box*0.3/max(box);
quartile_x = ones(11,1)*[1:nc] + [-a;-1;-1;1;1;a;1;1;-1;-1;-a]*box;
quartile_y = s([3,7,4,4,7,3,6,2,2,6,3],:);

% Draw a line through the median
median_x = ones(2,1)*[1:nc] + [-a;+a]*box;
median_y = s([3,3],:);

% Chop all boxes which don't have enough stats
quartile_x(:,chop) = [];
quartile_y(:,chop) = [];
whisker_x(:,[chop,chop+nc]) = [];
whisker_y(:,[chop,chop+nc]) = [];
median_x(:,chop) = [];
median_y(:,chop) = [];

% Add caps to the remaining whiskers
cap_x = whisker_x;
cap_x(1,:) = cap_x(1,:) - 0.05;
cap_x(2,:) = cap_x(2,:) + 0.05;
cap_y = whisker_y([1,1],:);

% Do the plot
if vertical
	plot(quartile_x, quartile_y, 'b-')
	hold on
	plot(whisker_x, whisker_y, 'b-')
	plot(cap_x, cap_y, 'b-')
	plot(median_x, median_y, 'r-')
	plot(outliers_x, outliers_y, [symbol(1),'r'])
        plot(outliers2_x, outliers2_y, [symbol(2),'r']);
else
	plot(quartile_y, quartile_x, 'b-')
	hold on
	plot(whisker_y, whisker_x, 'b-')
	plot(cap_y, cap_x, 'b-')
	plot(median_y, median_x, 'r-')
	plot(outliers_y, outliers_x, [symbol(1),'r'])
        plot(outliers2_y, outliers2_x, [symbol(2),'r']);
end

hold off
return

%-----------------------------------------------------------------------
function y = prctile(x,p);
%-----------------------------------------------------------------------
%PRCTILE gives the percentiles of the sample in X.
%   Y = PRCTILE(X,P) returns a value that is greater than P percent
%   of the values in X. For example, if P = 50  Y is the median of X. 
%
%   P may be either a scalar or a vector. For scalar P, Y is a row   
%   vector containing Pth percentile of each column of X. For vector P,
%   the ith row of Y is the P(i) percentile of each column of X.

%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision$  $Date$

[prows pcols] = size(p);
if prows ~= 1 & pcols ~= 1
    error('P must be a scalar or a vector.');
end
if any(p > 100) | any(p < 0)
    error('P must take values between 0 and 100');
end

if (~any(isnan(x)))
   y = prctilecol(x,p);
else                    % if there are NaNs, process each column
   if (size(x,1) == 1)
      x = x';
   end
   c = size(x,2);
   np = length(p);
   y = zeros(np,c);
   for j=1:c
      xx = x(:,j);
      xx = xx(~isnan(xx));
      y(:,j) = prctilecol(xx,p)';
   end
   if (min(size(x)) == 1)
      y = y';
   end
end

return
      
%-----------------------------------------------------------------------
function y = prctilecol(x,p);
%-----------------------------------------------------------------------
xx = sort(x);
[m,n] = size(x);

if m==1 | n==1
    m = max(m,n);
	if m == 1,
	   y = x*ones(length(p),1);
	   return;
	end
    n = 1;
    q = 100*(0.5:m - 0.5)./m;
    xx = [min(x); xx(:); max(x)];
else
    q = 100*(0.5:m - 0.5)./m;
    xx = [min(x); xx; max(x)];
end

q = [0 q 100];
y = interp1(q,xx,p);

return

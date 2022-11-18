function cat_batch_spm(batchname,varargin)
% wrapper for using spm8 batch mode (see cat_batch_cat12.sh)
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

if nargin < 1
	fprintf('Syntax: cat_batch_spm(batchname)\n');
	exit
end
if nargin > 1
% initialize variables given by the shell script to adapt batches

	if mod(numel(varargin),2) ~= 0
		fprintf('Syntax: cat_batch_spm(batchname,''var1name'',var1,...)\n');
		exit
	end
	for vi = 1:2:numel(varargin)-1
		try
			eval(sprintf('%s = varargin{vi+1};', varargin{vi} )); 
		catch
			fprintf('Cannot evaluate variable %0.0f "%s = varargin{%0.0f+1}; " \n', (vi+1)/2, varargin{vi}, (vi+1)/2 );
			fprintf('Syntax: cat_batch_spm(batchname,''var1name'',var1,...)\n');
			exit
		end	
	end
	clear vi; 
	
	% show result
	whos	
end


% set up SPM enviroment (the batch may need SPM variables/functions)
spm_get_defaults;
global defaults
spm_jobman('initcfg'); 

if ~exist(batchname,'file')
	fprintf('Batchfile %s not found\n',batchname);
	exit
end

eval(batchname)

if ~exist('matlabbatch','var')
	fprintf('Batchfile %s did not returned variable matlabbatch.\n', batchname);
	exit
end

warning off
try
	spm_jobman('run',matlabbatch);
catch
	fprintf('Batchfile %s was not running successfully.\n', batchname);

	exit
end	
warning off
exit

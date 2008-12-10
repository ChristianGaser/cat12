function cg_vbm8_debug
%cg_vbm8_debug	print debug information for SPM8 and VBM8
%
% FORMAT cg_vbm8_debug
%
%__________________________________________________________________________
% Christian Gaser
% $Id$

% print last error
fprintf('\nLast error message:\n');
fprintf('-------------------------------------------------------------------------------------\n');
fprintf('-------------------------------------------------------------------------------------\n');
try
	er = lasterror;
	fprintf('%s\n',er.message);
	if isfield(er,'stack')
		for i=1:length(er.stack)
			fprintf('%s at line %g\n',char(er.stack(i).file),er.stack(i).line);
		end
	end
catch
	fprintf('%s\n',lasterr);
end

fprintf('-------------------------------------------------------------------------------------\n');
fprintf('-------------------------------------------------------------------------------------\n');

% find release Id of VBM8
mext = {'.m','.c','.h','.man'};  %- MEX
vbmdir = spm('Dir','cg_config_vbm8.m');
d = dir(fullfile(vbmdir,'*'));
f = {d(~[d.isdir]).name};
d = {d([d.isdir]).name}; d = {d{~ismember(d,{'.' '..'})}};
L = max(cellfun('length',f));
pat = '\$Id: (\S+) (\d+) ([0-9-]+) ([0-9:]+Z) (\w+) \$';
Id = [];
for i=1:length(f)
    [pathstr, name, ext] = fileparts(f{i});
    if ismember(ext,mext)
  	  fid = fopen(fullfile(vbmdir,f{i}),'r');
  	  if fid == -1, continue; end
  	  V = 'none';
  	  while 1
  		  tline = fgetl(fid);
  		  if ~ischar(tline), break, end
  		  tok = regexp(tline, pat, 'tokens');
  		  if ~isempty(tok), V = tok{1}{2}; break; end
  	  end
  	  fclose(fid);
  	  if ~strcmp(V,'none')
  	  	  Id = [Id str2num(V)];
	  end
    end
end

vbm_ver = max(Id);

fprintf('\nVersion information:\n');
fprintf('-------------------------------------------------------------------------------------\n');
fprintf('VBM8:\Revision %d\n',vbm_ver);

ver


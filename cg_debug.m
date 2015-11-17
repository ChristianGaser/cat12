function cg_debug
%cg_debug	print debug information for SPM12 and CAT12
%
% FORMAT cg_debug
%
%__________________________________________________________________________
% Christian Gaser
% $Id$

rev = '$Rev$';

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

fprintf('\nVersion information:\n');
fprintf('-------------------------------------------------------------------------------------\n');

ver


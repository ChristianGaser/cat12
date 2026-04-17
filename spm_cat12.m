function spm_cat12(varargin)
% Legacy wrapper — redirects to spm_CAT after cat12 to CAT migration.
%
% After updating from the old cat12 layout the old cat_update.m calls
% spm_cat12.  This file lives only in the new CAT folder, so MATLAB
% finds it here.  We delegate to spm_CAT which handles the migration
% (removing the old cat12 directory) and then launches CAT normally.
% ______________________________________________________________________

spm_CAT(varargin{:});

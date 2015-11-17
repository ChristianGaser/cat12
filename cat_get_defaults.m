function varargout = cat_get_defaults(defstr, varargin)
% Get/set the defaults values associated with an identifier
% FORMAT defval = cat_get_defaults(defstr)
% Return the defaults value associated with identifier "defstr". 
% Currently, this is a '.' subscript reference into the global  
% "defaults" variable defined in spm_defaults.m.
%
% FORMAT cat_get_defaults(defstr, defval)
% Sets the cat value associated with identifier "defstr". The new
% cat value applies immediately to:
% * new modules in batch jobs
% * modules in batch jobs that have not been saved yet
% This value will not be saved for future sessions of SPM. To make
% persistent changes, edit cat_defaults.m.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% based on Volkmar Glauches version of
% spm_get_defaults
% $Id$

global cat12;
if isempty(cat12)
    cat_defaults;
end

if nargin == 0
    varargout{1} = cat12;
    return
end

% construct subscript reference struct from dot delimited tag string
tags = textscan(defstr,'%s', 'delimiter','.');
subs = struct('type','.','subs',tags{1}');

if nargin == 1
    varargout{1} = subsref(cat12, subs);
else
    cat12 = subsasgn(cat12, subs, varargin{1});
end

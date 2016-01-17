function [CATrel, CATver]  = cat_version
% check for CAT revision
%
% FORMAT [CATrel, CATver] = cat_version
% 
% This function will retrieve the CAT release and version and is a
% modified version of spm_version.m
%_______________________________________________________________________
% Christian Gaser
% $Id$

persistent CAT_VER;
v = CAT_VER;
if isempty(CAT_VER)
    v = struct('Name','','Version','','Release','','Date','');
    % try Contents.txt and then Contents.m
    try
        vfile = fullfile(spm('Dir'),'toolbox','cat12','Contents.txt');
        if ~exist(vfile,'file')
          vfile = fullfile(spm('Dir'),'toolbox','cat12','Contents.m');
        end
        fid = fopen(vfile,'rt');
        if fid == -1, error('Can''t open %s.',vfile); end
        l1 = fgetl(fid); l2 = fgetl(fid);
        fclose(fid);
        l1 = strtrim(l1(2:end)); l2 = strtrim(l2(2:end));
        t  = textscan(l2,'%s','delimiter',' '); t = t{1};
        v.Name = l1; v.Date = t{4};
        v.Version = t{2}; v.Release = t{3}(2:end-1);
    catch
        error('Can''t obtain CAT Revision information.');
    end
    CAT_VER = v;
end
CATrel = v.Release;
CATver = v.Version;

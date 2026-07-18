function out = cat_conf_BIDSinput_run(job)
%cat_conf_BIDSinput_run. Batch prog of the BIDS input add-on.
%
% ADD-ON (not part of upstream CAT12). See cat_conf_BIDSinput.
%
%   out = cat_conf_BIDSinput_run(job)
%
%   out.cross         .. cross-sectional files (subjects with <=1 session)
%   out.longsubjects  .. 1xS cell, one cellstr of files per longitudinal
%                        subject (subjects with >=2 sessions)
%   out.long          .. all longitudinal files in one list
%   out.files         .. every matching file
%
% The module itself does not process anything - it only resolves the
% dataset and prints the report. The actual preprocessing is done by the
% CAT12 modules that consume these outputs as batch dependencies.
%
% See also cat_conf_BIDSinput, cat_io_BIDSquery
% ______________________________________________________________________
%
% BIDS input add-on for CAT12 (fork).
% ______________________________________________________________________

  res = cat_conf_BIDSinput_query(job, 0);

  out               = struct();
  out.files         = res.files;
  out.cross         = res.cross;
  out.long          = res.long;
  out.longsubjects  = res.longsubjects;
  out.longsublabels = res.longsublabels;
  out.subjects      = res.subjects;
  out.bidsdir       = res.bidsdir;

  if isempty(out.files)
    msg = ['  No input files were found. Check the dataset root and the ' ...
           'filters above.\n\n'];
    try
      cat_io_cprintf('err', msg);   % needs a Java command window
    catch
      fprintf(msg);
    end
  end
end

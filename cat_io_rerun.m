function run = cat_io_rerun(files,filedates,verb)
%cat_io_rerun(f1,fd). Test if a file f1 is newer than another file/date fd.  
% This function is used to estimated if a file is newer than another given 
% file or date. For instance file is the result of anther file that was 
% changed in the meantime, it has to be reprocessed. 
%
%  run = cat_io_rerun(files,filedates)
% 
%  run      .. logical vector with the number of given files
%              cell if directories or wildcards are used
%  files    .. filenames (cellstr or char)
%  filedat  .. filenames (cellstr or char) or datetimes or datenum
%
% Examples: 
%  1) Is the working directory younger than the SPM dir?
%     cat_io_rerun(pwd,spm('dir'); 
%
%  2) Is the working directory younger than one month?
%     cat_io_rerun(pwd,clock - [0 1 0 0 0 0]) 
%   
%  3) Is this function younger than one year?
%     cat_io_rerun(which('cat_io_rerun'),clock - [1 0 0 0 0 0]) 
%
% _________________________________________________________________________
% Robert Dahnke
% $Id$

  if ~exist('verb','var'), verb = 0; end
  files = cellstr(files);
  if iscellstr(filedates) || ischar(filedates)
    filedates = cellstr(filedates);
    if numel(filedates) == 1
      filedates = repmat(filedates,numel(files),1);
    else
      if ~isempty(filedates) && numel(files) ~= numel(filedates)
        error('ERROR:cat_io_rerun:inputsize','Number of files and filedates has to be equal.\n')
      end
    end  
  else 
    if size(filedates,1)
      filedates = repmat(filedates,numel(files),1);
    end
  end
  
  run = ones(size(files)); 
  for fi = 1:numel(files)
    [pp,ff,ee] = spm_fileparts(files{fi}); files{fi} = fullfile(pp,[ff ee]); % remove additional dimensions ",1" 
    if ~exist(files{fi},'file')
      run(fi) = 1; 
    else 
      fdata = dir(files{fi});
      if numel(fdata)>1
        run = num2cell(run); 
      end
      
      if exist('filedates','var') && iscellstr(filedates) 
        [pp,ff,ee] = spm_fileparts(filedates{fi}); filedates{fi} = fullfile(pp,[ff ee]); % remove additional dimensions ",1" 
        if exist(filedates{fi},'file')
          fdata2 = dir(filedates{fi});
          if verb
            fprintf('Input 1: %50s: %s\n',fdata.name,datestr(fdata.datenum)); 
            fprintf('Input 2: %50s: %s\n',fdata2.name,datestr(fdata2.datenum));
          end
          if numel(fdata)>1
            run{fi} = [fdata(:).datenum] < fdata2.datenum;
          else
            run(fi) = fdata.datenum < fdata2.datenum;
          end
        end
      elseif ~isempty(filedates) && isdatetime( filedates(fi,:) )
        if numel(fdata)>1
          run{fi} = [fdata(:).datenum] < datenum( filedates(fi,:) );
        else
          run(fi) = fdata.datenum < datenum( filedates(fi,:) );
        end
      else
        run(fi) = 1;
      end
    end
  end
end

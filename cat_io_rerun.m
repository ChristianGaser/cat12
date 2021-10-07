function run = cat_io_rerun(files,filedates,verb)
%cat_io_rerun(f1,fd). Test if a file f1 is newer than another file/date fd.  
% This function is used to estimated if a file is newer than another given 
% file or date. For instance file is the result of another file that was 
% changed in the meantime, it has to be reprocessed. 
%
%  run = cat_io_rerun(files,filedates,verb)
% 
%  run      .. logical vector with the number of given files
%              cell if directories or wildcards are used
%  files    .. filenames (cellstr or char)
%  filedat  .. filenames (cellstr or char) or datetimes or datenum
%
% Examples: 
%  1) Is the working directory younger/newer than the SPM dir?
%     cat_io_rerun(pwd,spm('dir'); 
%
%  2) Is the working directory younger/newer than one month?
%     cat_io_rerun(pwd,clock - [0 1 0 0 0 0]) 
%   
%  3) Is this function younger than one year?
%     cat_io_rerun(which('cat_io_rerun'),clock - [1 0 0 0 0 0]) 
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$
  
  if ~exist('verb','var'), verb = 0; end
  files = cellstr(files);

  % only use that function in developer mode because it's simply too dangerous if files
  % are not processed if already existing and parameter changed
  if cat_get_defaults('extopts.expertgui') < 2
    run = ones(size(files));
    return
  end
  
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
  exf = ones(size(files)); 
  for fi = 1:numel(files)
    [pp,ff,ee] = spm_fileparts(files{fi}); files{fi} = fullfile(pp,[ff ee]); % remove additional dimensions ",1" 
    if ~exist(files{fi},'file')
      if numel(files)>1
        run{fi} = 1; 
        if verb
          fprintf(' Input file 2-%02d does not exist: %40s\n',fi,spm_str_manip( files{fi}, 'a40'));
        end
      else
        run(fi) = 1; 
        if verb
          fprintf(' Input file 2 does not exist: %43s\n',spm_str_manip( files{fi}, 'a43'));
        end
      end
      exf = 0; 
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
            if numel(files)==1
              fprintf(' Input file 1: %50s: %s\n',spm_str_manip( fdata.name , 'a50'),datestr(fdata.datenum) ); 
              fprintf(' Input file 2: %50s: %s\n',spm_str_manip( fdata2.name, 'a50'),datestr(fdata2.datenum));
            else
              if fi == 1
                fprintf(' Input file 1:     %50s: %s\n',   spm_str_manip( fdata.name ,'a50'),datestr(fdata.datenum) ); 
              end
              fprintf(' Input file 2-%02d: %50s: %s\n',fi,spm_str_manip( fdata2.name,'a50'),datestr(fdata2.datenum));
            end
          end
          if numel(fdata)>1
            run{fi} = [fdata(:).datenum] > fdata2.datenum;
          else
            run(fi) = fdata.datenum > fdata2.datenum;
          end
        end
      elseif ~isempty(filedates) && isdatetime( filedates(fi,:) )
        if numel(fdata)>1
          run{fi} = [fdata(:).datenum] > datenum( filedates(fi,:) );
        else
          run(fi) = fdata.datenum > datenum( filedates(fi,:) );
        end
      else
        if numel(fdata)>1
          run{fi} = 1; 
        else
          run(fi) = 1;
        end
      end
    end
  end
  if verb 
    if (iscell(run) && any(cell2mat(run))) || ( ismatrix(run) && any(run) )
      if all(exf)
        cat_io_cprintf([0.5 0.0 0.0],' Reprocessing is required. \n'); 
      elseif all(exf==0) && numel(files)>1
        cat_io_cprintf([0.5 0.0 0.0],' (Re)processing is required. \n'); 
      else
        cat_io_cprintf([0.5 0.0 0.0],' Processing is required. \n');
      end
    else
      cat_io_cprintf([0.0 0.5 0.0],' Reprocessing is NOT required. \n'); 
    end
  end
end

function run = cat_io_lazy(files,filedates,verb,force)
%cat_io_lazy. Test if a file is newer than another file.  
% This function is used to estimated if a file is newer than another given 
% file or date. For instance file is the result of anther file that was 
% changed in the meantime, it has to be reprocessed. 
%
%  run = cat_io_lazy(files,filedates,verb)
% 
%  run      .. logical vector with the number of given files
%              cell if directories or wildcards are used
%  files    .. filenames (cellstr or char)
%  filedat  .. filenames (cellstr or char) or datetimes or datenum
%  verb     .. print details about the files and about the result 
%               (default = 0.5 only display if reprocessing is NOT reqired)
%  force    .. use also in non developer mode (default = 0)
%
% Examples: 
%  1) Is the working directory younger than the SPM dir?
%     cat_io_lazy(pwd,spm('dir'); 
%
%  2) Is the working directory younger than one month?
%     cat_io_lazy(pwd,clock - [0 1 0 0 0 0]) 
%   
%  3) Is this function younger than one year?
%     cat_io_lazy(which('cat_io_lazy'),clock - [1 0 0 0 0 0]) 
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  if ~exist('verb','var'),  verb  = 0.5; end
  if ~exist('force','var'), force = 1; end
  
  % only use that function in developer mode because it's simply too dangerous if files
  % are not processed if already existing and parameter changed
  if cat_get_defaults('extopts.expertgui') < 2 && ~force 
    run = zeros(size(files));
    return
  end
  
  files = cellstr(files);
  if iscellstr(filedates) || ischar(filedates)
    filedates = cellstr(filedates);
    if numel(filedates) == 1
      filedates = repmat(filedates,numel(files),1);
    else
      if ~isempty(filedates) && numel(files) ~= numel(filedates)
        error('ERROR:cat_io_lazy:inputsize','Number of files and filedates has to be equal.\n')
      end
    end  
  else 
    if size(filedates,1)
      filedates = repmat(filedates,numel(files),1);
    end
  end
  
  run = ones(size(files)); 
  for fi = 1:numel(files)
    if ~exist(files{fi},'file')
      run(fi) = 1; 
    else 
      fdata = dir(files{fi});
      if numel(fdata)>1
        run = num2cell(run); 
      end
      if exist('filedates','var') && iscellstr(filedates) && exist(filedates{fi},'file')
        fdata2 = dir(filedates{fi});
        if numel(fdata)>1
          run{fi} = [fdata(:).datenum] < fdata2.datenum;
        else
          run(fi) = fdata.datenum < fdata2.datenum;
        end
      elseif ~isempty(filedates) 
        if numel(fdata)>1
          run{fi} = [fdata(:).datenum] < datenum( filedates(fi,:) );
        else
          run(fi) = fdata.datenum < datenum( filedates(fi,:) );
        end
      end
      % be verbose only if verb>=1 or if no reprocessing is required  
      if verb >= 1 || (verb && ~( (iscell(run) && any(cell2mat(run))) || ( ismatrix(run) && any(run) ) ))
        fprintf('\n'); 
        if numel(files)==1
          fprintf('\n Input file 1: %50s: %s\n',spm_str_manip( fdata.name , 'a50'),datestr(fdata.datenum) ); 
          fprintf('\n Input file 2: %50s: %s\n',spm_str_manip( fdata2.name, 'a50'),datestr(fdata2.datenum));
        else
          if fi == 1
            fprintf('\n Input file 1:     %50s: %s\n',   spm_str_manip( fdata.name ,'a50'),datestr(fdata.datenum) ); 
          end
          fprintf('\n Input file 2-%02d: %50s: %s\n',fi,spm_str_manip( fdata2.name,'a50'),datestr(fdata2.datenum));
        end
      end
    end
  end
  % be verbose only if verb>=1 or if no reprocessing is required  
  if verb >= 1 || ~( (iscell(run) && any(cell2mat(run))) || ( ismatrix(run) && any(run) ) ) 
    if verb >= 1 && ( (iscell(run) && any(cell2mat(run))) || ( ismatrix(run) && any(run) ) )
      if all(exf)
        cat_io_cprintf([0.5 0.0 0.0],' Reprocessing is required. \n'); 
      elseif all(exf==0) && numel(files)>1
        cat_io_cprintf([0.5 0.0 0.0],' (Re)processing is required. \n'); 
      else
        cat_io_cprintf([0.5 0.0 0.0],' Processing is required. \n');
      end
    elseif verb 
      cat_io_cprintf([0.0 0.5 0.0],' Reprocessing is NOT required. \n'); 
    end
  end
end

function run = cat_io_rerun(files,filedates,verb,force)
%cat_io_rerun(f1,fd). Test if a file f1 is newer than another file/date fd.  
% This function is used to estimated if a file is newer than another given 
% file or date. For instance file is the result of another file that was 
% changed in the meantime, it has to be reprocessed. 
%
%  run = cat_io_rerun(files,filedates,verb,force)
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
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$
  
  if ~exist('verb','var'),  verb  = 0.5; end
  if ~exist('force','var'), force = 0; end
  files = cellstr(files);

  % only use that function in developer mode because it's simply too dangerous
  % if files are not processed if already existing and parameter changed
  if force>=0 && (cat_get_defaults('extopts.expertgui') < 2 || force~=0)
    if verb, cat_io_cprintf([0.5 0.0 0.0],' Reprocessing! \n'); end
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
        run(fi) = 1; 
        if verb
          fprintf(' Input file 2-%02d does not exist: %70s\n',fi,spm_str_manip( files{fi}, 'a70'));
        end
      else
        run(fi) = 1; 
        if verb
          fprintf(' Input file 2 does not exist: %73s\n',spm_str_manip( files{fi}, 'a73'));
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
          if numel(fdata)>1
            run(fi) = [fdata(:).datenum] >= fdata2.datenum;
          else
            run(fi) = fdata.datenum >= fdata2.datenum;
          end

          % be verbose only if verb>=1 or if no reprocessing is required  
          if verb >= 1 || (verb && ~( (iscell(run) && any(cell2mat(run))) || ( ismatrix(run) && any(run) ) ))
            if fi==1, fprintf('\n'); end
            if numel(files)==1 && numel(filedates)==1
              fprintf(' Input file 1: %80s: %s\n',spm_str_manip( fdata.name , 'a80'),datestr(fdata.datenum) ); 
              fprintf(' Input file 2: %80s: %s -',spm_str_manip( fdata2.name, 'a80'),datestr(fdata2.datenum));
            elseif numel(files) == numel(filedates)
              fprintf(' Input file %02d-1: %80s: %s\n',fi,spm_str_manip( fdata.name ,'a80'),datestr(fdata.datenum) ); 
              fprintf(' Input file %02d-2: %80s: %s -',fi,spm_str_manip( fdata2.name,'a80'),datestr(fdata2.datenum));
            elseif numel(files) == 1
              if fi == 1
                fprintf(' Input file 1-%02d: %80s: %s\n',fi,spm_str_manip( fdata.name ,'a80'),datestr(fdata.datenum) ); 
              end
              fprintf(' Input file 1-%02d: %80s: %s -',fi,spm_str_manip( fdata2.name,'a80'),datestr(fdata2.datenum));
            else 
              error('Error:cat_io_rerun:inputerror','Size input1 does not match size of input2 ([n1,n2]: [1,1], [1,n], or [n,n]).')
            end
            if run(fi), cat_io_cprintf([0.5 0.0 0.0],' reprocess\n'); else, cat_io_cprintf([0.0 0.5 0.0],' do not process\n'); end
          end
        elseif verb > 1
          if numel(files)==1 && numel(filedates)==1
            cat_io_cprintf([0.5 0.0 0.0],' Input file 2: %80s: %s\n',spm_str_manip( filedates{fi}, 'a80'),'missing');          
          elseif numel(files) == numel(filedates)
            cat_io_cprintf([0.5 0.0 0.0],' Input file %02d-1: %80s: %s\n',fi,spm_str_manip( files{fi}    ,'a80'),datestr(fdata.datenum) ); 
            cat_io_cprintf([0.5 0.0 0.0],' Input file %02d-2: %80s: %s\n',fi,spm_str_manip( filedates{fi},'a80'),'missing');
          else
            cat_io_cprintf([0.5 0.0 0.0],' Input file 2-%02d: %80s: %s\n',fi,spm_str_manip( filedates{fi},'a80'),'missing');
          end
        end
      elseif ~isempty(filedates) && isdatetime( filedates(fi,:) )
        if numel(fdata)>1
          run{fi} = [fdata(:).datenum] >= datenum( filedates(fi,:) );
        else
          run(fi) = fdata.datenum >= datenum( filedates(fi,:) );
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

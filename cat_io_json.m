function S = cat_io_json(varargin)
%cat_io_json. Read json file.   

  S = struct();

  if nargin == 1
    if isstruct(varargin{1})
      % job case
      files = job.files; 
    else
      % just a file
      files = cellstr(varargin{1});
    end
  end

  for fi = 1:numel(files)
    if ~exist(files{fi},'file')
      cat_io_cprintf('err','ERROR: Miss "%s"\n',files{fi})
    else
      fid = fopen(files{fi}); 
      raw = fread(fid,inf); 
      str = char(raw'); 
      fclose(fid); 
      val = jsondecode(str);

      if fi == 1
        S = val; 
      else
        S = cat_io_mergeStruct(S,val); 
      end
    end
  end
end
  
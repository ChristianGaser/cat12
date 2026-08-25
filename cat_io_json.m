function S = cat_io_json(varargin)
%cat_io_json. Write and read json files.   
% Json did not allow separation between row/column matrix's and reading 
% always result in row matrix's! 
 
  if nargin == 1 && ischar(varargin{1}) && strcmp( varargin{1} , 'unittest')
    unittest;
  elseif nargin == 1
    S = readjson(varargin{:});
  elseif nargin == 2
    writejson(varargin{:})
  else
    error('Incorrect number of inputs.\n'); 
  end

end
% =========================================================================
function S = readjson(varargin)
%readjson. Read json files.
  S = struct();
  
  if isstruct(varargin{1})
    % job case
    files = job.files; 
  else
    % just a file
    files = cellstr(varargin{1});
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
% =========================================================================
function writejson(varargin)
%writejson. write json files.

  if iscell( varargin{1} ) && numel( varargin{1} ) == numel( varargin{2} )
  % write multiple files
  
    for i = 1:numel( varargin{1} )
      if isstruct( varargin{2}(i) )
        cat_io_json(varargin{1}{i}, varargin{2}(i));
      elseif iscell( varargin{2} )
        cat_io_json(varargin{1}{i}, varargin{2}{i});
      end
    end
  
    return
  end
  
  % check input properties
  if ~ischar( varargin{1} ) && isfile( varargin{1} )
    error('First argument has to be char and valide file name.')
  end
  if ~strcmp( spm_file( varargin{1} , 'ext') ,'json')
    error('File extension has to be .json');
  end

  if ~isstruct( varargin{2} ) 
    error('Second argument has to be a structure.')
  end

  if ~exist( spm_fileparts( varargin{1}), 'dir') 
    mkdir( spm_fileparts( varargin{1}) );
  end

  % write data
  f0  = fopen( varargin{1} ,'w');
  txt = jsonencode(varargin{2}); %,'PrettyPrint', true);
  fwrite( f0, txt); 
  fclose( f0 );

end
% =========================================================================
function unittest
%unittest. Quick visual test.

  %% prepare data
  data(1,1) = struct('numbers',[1 2 3],'char','text1','cellstr',{{'Hello','World'}}); 
  data(2,1) = struct('numbers',[4 5 6],'char','text2','cellstr',{{'Goodbye','World'}}); 
  fname{1} = fullfile(tmffolder,'testfile1.json'); 
  fname{2} = fullfile(tmffolder,'testfile2.json'); 
  
  % write and read single record into one file
  cat_io_json( fname{1}, data(1)); 
  data1 = cat_io_json( fname{1} ); 
  fprintf('\nTest1 -  single record - one file: \n'), disp(data(1)); disp(data1);

  % write and read multiple records into one file
  cat_io_json( fname{1}, data); 
  data2 = cat_io_json( fname{1} ); 
  fprintf('\nTest2 -  multiple record - one file: \n'), disp(data); disp(data2);

  % write and read multiple records into multiple files
  cat_io_json( fname, data); 
  data3 = cat_io_json( fname ); 
  fprintf('\nTest3 -  multiple record - multiple file: \n'), disp(data); disp(data3);

  % cleanup
  for i=1:numel(fname), try delete(fname{i}); end; end
end
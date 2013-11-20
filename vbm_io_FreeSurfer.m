function varargout=vbm_io_FreeSurfer(action,varargin)

  switch action
    case 'write_surf'
      write_surf(varargin{1}, varargin{2}.vertices, varargin{2}.faces);
    case 'read_surf'
      [varargout{1}.vertices,varargout{1}.faces] = read_surf(varargin{1}); 
      varargout{1}.faces = varargout{1}.faces+1;   
    case 'write_surf_data'
      write_wfile(varargin{1},varargin{2});
    case 'read_surf_data'
      [varargout{1},varargout{2}] = read_wfile(varargin{1});
    otherwise
  end

end

function write_surf(fname, vert, face)
% write_surf - FreeSurfer I/O function to write a surface file
% 
% write_surf(fname, vert, face)
% 
% writes a surface triangulation into a binary file
% fname - name of file to write
% vert  - Nx3 matrix of vertex coordinates
% face  - Mx3 matrix of face triangulation indices
% 
% The face matrix here must be matlab compatible
% (no zero indices).  It is converted to FreeSurfer
% indices that start at zero.
% 
% See also freesurfer_read_surf, freesurfer_write_curv, freesurfer_write_wfile

  if(nargin ~= 3)
    fprintf('USAGE: freesurfer_write_surf(fname, vert, face)\n');
    return;
  end

  if size(vert,2) ~= 3,
      error('vert must be Nx3 matrix');
  end

  if size(face,2) ~= 3,
      error('face must be Mx3 matrix');
  end

  %fprintf('...subtracting 1 from face indices for FreeSurfer compatibility.\n');
  face = face - 1;

  % open it as a big-endian file
  fid = fopen(fname, 'wb', 'b') ;

  TRIANGLE_FILE_MAGIC_NUMBER = 16777214 ;
  fwrite3(fid, TRIANGLE_FILE_MAGIC_NUMBER);

  vnum = size(vert,1) ;  % number of vertices
  fnum = size(face,1) ;  % number of faces

  % Ouput a couple of text lines with creation date
  fprintf(fid,'created by %s on %s\n\n',getenv('USER'),datestr(now)); % creation date 

  fwrite(fid, vnum,'int32');
  fwrite(fid, fnum,'int32');

  % reshape vert into column array and write
  vert = reshape(vert',size(vert,1)*size(vert,2),1);
  fwrite(fid, vert,'float32');

  % reshape face into column array and write
  face = reshape(face',size(face,1)*size(face,2),1);
  fwrite(fid, face,'int32');

  fclose(fid) ;

end

function [vertex_coords, faces] = read_surf(fname)
  %
  % [vertex_coords, faces] = read_surf(fname)
  % reads a the vertex coordinates and face lists from a surface file
  % note that reading the faces from a quad file can take a very long
  % time due to the goofy format that they are stored in. If the faces
  % output variable is not specified, they will not be read so it 
  % should execute pretty quickly.
  %


  %
  % read_surf.m
  %
  % Original Author: Bruce Fischl
  % CVS Revision Info:
  %    $Author$
  %    $Date$
  %    $Revision$
  %
  % Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
  %
  % Terms and conditions for use, reproduction, distribution and contribution
  % are found in the 'FreeSurfer Software License Agreement' contained
  % in the file 'LICENSE' found in the FreeSurfer distribution, and here:
  %
  % https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
  %
  % Reporting: freesurfer@nmr.mgh.harvard.edu
  %


  %fid = fopen(fname, 'r') ;
  %nvertices = fscanf(fid, '%d', 1);
  %all = fscanf(fid, '%d %f %f %f %f\n', [5, nvertices]) ;
  %curv = all(5, :)' ;

  % open it as a big-endian file


  %QUAD_FILE_MAGIC_NUMBER =  (-1 & 0x00ffffff) ;
  %NEW_QUAD_FILE_MAGIC_NUMBER =  (-3 & 0x00ffffff) ;

  TRIANGLE_FILE_MAGIC_NUMBER =  16777214 ;
  QUAD_FILE_MAGIC_NUMBER =  16777215 ;

  fid = fopen(fname, 'rb', 'b') ;
  if (fid < 0)
    error('MATLAB:vbm_io_FreeSurfer:read_surf','could not open curvature file %s.', fname) ;
  end
  magic = fread3(fid) ;

  if(magic == QUAD_FILE_MAGIC_NUMBER)
    vnum = fread3(fid) ;
    fnum = fread3(fid) ;
    vertex_coords = fread(fid, vnum*3, 'int16') ./ 100 ; 
    faces = zeros(fnum,4,'single');
    if (nargout > 1)
      for i=1:fnum
        for n=1:4
          faces(i,n) = fread3(fid) ;
        end
      end
    end
  elseif (magic == TRIANGLE_FILE_MAGIC_NUMBER)
    fgets(fid) ;
    fgets(fid) ;
    vnum = fread(fid, 1, 'int32') ;
    fnum = fread(fid, 1, 'int32') ;
    vertex_coords = fread(fid, vnum*3, 'float32') ; 
    faces = fread(fid, fnum*3, 'int32') ;
    faces = reshape(faces, 3, fnum)' ;
  end

  vertex_coords = reshape(vertex_coords, 3, vnum)' ;
  fclose(fid) ;
end

function fwrite3(fid, val)
%
% fwrite3.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author$
%    $Date$
%    $Revision$
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

% write a 3 byte integer out of a file
%fwrite(fid, val, '3*uchar') ;
  
  b1 = bitand(bitshift(val, -16), 255) ;
  b2 = bitand(bitshift(val, -8), 255) ;
  b3 = bitand(val, 255) ; 
  fwrite(fid, b1, 'uchar') ;
  fwrite(fid, b2, 'uchar') ;
  fwrite(fid, b3, 'uchar') ;
end

function [retval] = fread3(fid)
  % [retval] = fd3(fid)
  % read a 3 byte integer out of a file


  %
  % fread3.m
  %
  % Original Author: Bruce Fischl
  % CVS Revision Info:
  %    $Author$
  %    $Date$
  %    $Revision$
  %
  % Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
  %
  % Terms and conditions for use, reproduction, distribution and contribution
  % are found in the 'FreeSurfer Software License Agreement' contained
  % in the file 'LICENSE' found in the FreeSurfer distribution, and here:
  %
  % https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
  %
  % Reporting: freesurfer@nmr.mgh.harvard.edu
  %

  b1 = fread(fid, 1, 'uchar') ;
  b2 = fread(fid, 1, 'uchar') ;
  b3 = fread(fid, 1, 'uchar') ;
  retval = bitshift(b1, 16) + bitshift(b2,8) + b3 ;

end

function err = write_wfile(fname, w, v)
% err = write_wfile(fname, w, <v>)
% 
% writes a vector into a binary 'w' file
%  fname - name of file to write to
%  w     - vector of values to be written
%  v     - 0-based vertex numbers 
%          (assumes 0 to N-1 if not present or empty).
%
% See also read_wfile.
%


%
% write_wfile.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author$
%    $Date$
%    $Revision$
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%


  err = 1;

  if(nargin ~= 2 && nargin ~= 3)
    fprintf('USAGE: err = write_wfile(fname, w, <v>) \n');
    return;
  end

  vnum = length(w) ;

  % Handle when v not given or is empty %
  if (exist('v','var') ~= 1), v = []; end
  if (isempty(v)), v = 0:vnum-1; end

  % open it as a big-endian file
  fid = fopen(fname, 'wb', 'b') ;
  if(fid == -1)
    fprintf('ERROR: could not open %s\n',fname);
    return;
  end

  fwrite(fid, 0, 'int16') ;
  fwrite3(fid, vnum) ;
  for i=1:vnum
    fwrite3(fid, v(i)) ;          % vertex number (0-based)
    fwrite(fid,  w(i), 'float') ; % vertex value
  end

  fclose(fid) ;

  err = 0;
end

function [w,v] = read_wfile(fname)
%
% [w,v] = read_wfile(fname)
% reads a vector into a binary 'w' file
% fname - name of file to write to
% w     - vector of values to be written
% v     - vector of vertex indices (0-based)
%
% See also write_wfile.
%


%
% read_wfile.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author$
%    $Date$
%    $Revision$
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%


  w = [];
  v = [];

  if(nargin ~= 1)
    fprintf('[w,v] = read_wfile(fname)\n');
    return;
  end

  % open it as a big-endian file
  fid = fopen(fname, 'rb', 'b') ;
  if (fid < 0)
    error('MATLAB:vbm_io_FreeSurfer:read_wfile','ERROR: Could not open w file %s.', fname) ;
  end

  fread(fid, 1, 'int16') ;
  vnum = fread3(fid) ;
  w = zeros(vnum,1) ;
  v = zeros(vnum,1) ;
  for i=1:vnum
    v(i) = fread3(fid) ; % vertex number (0-based)
    w(i) = fread(fid, 1, 'float') ; % vertex value
  end

  fclose(fid) ;
end
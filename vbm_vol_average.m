function VO = vbm_vol_average(P,fname,PT,dt,nr)
% ______________________________________________________________________
% Creates median images of a set of files P or volumes V with the same
% image properties.
%
%   VO = vbm_vol_median(V[,fname,PT,tpye,nr])
%   VO = vbm_vol_median(P[,fname,PT,dt,nr])
%
%   P     = char or cell array with filenames
%   PT    = template image with another resolution (for interpolation)
%   V     = SPM volume structure
%   VO    = output volume
%   dt    = for [median,mean,std] images with spm nii datatype
%           (0=no image, 2=uint8, 4=int16, 16=single, ...)
%   nr    = use number in fname
%   fname = name of the outputfile, add median/mean/std automaticly
%            mypath/myname.nii > mypath/median003_myname.nii ...
% ______________________________________________________________________
% Robert Dahnke 2013_08
% Structural Brain Mapping Group
% University Jena
%  
% $Id$
% ______________________________________________________________________
  
  if isstruct(P)
    V = P; clear P;
    for fi = 1:numel(V), P{fi} = V(fi).fname; end
  end
  P = char(P);
  if isempty(P), VO=struct(); return; end
  
  if ~exist('nr','var') || isemtpy(nr)
    nr = 0;
  end
  if nr==0
    nstr = '';
  else
    nstr=num2str(size(P,1),'%3.0f');
  end
  
  if ~exist('fname','var') || isempty(fname)
    [pp,ff,ee] = spm_fileparts(P(1,:));
    fname2{1} = fullfile(pp,['median' nstr '_'  ff ee]);
    fname2{2} = fullfile(pp,['mean'   nstr '_'  ff ee]);
    fname2{3} = fullfile(pp,['std'    nstr '_'  ff ee]);
  else
    if iscell(fname)
      fname2 = fname; 
    else
      [pp,ff,ee] = spm_fileparts(fname);
      fname2{1} = fullfile(pp,['median' nstr '_' ff ee]);
      fname2{2} = fullfile(pp,['mean'   nstr '_' ff ee]);
      fname2{3} = fullfile(pp,['std'    nstr '_' ff ee]);      
    end
  end

  if ~exist('dt','var') || isempty(dt)
    % median, mean, sd
    dt = [16 16 16];
  end
  if size(P,1)<2,
    warging('MATLAB:vbm_vol_average:input','WARNING:vbm_vol_average:to less input (n=1)!\n');
    VO=struct(); return; 
  end
  if size(P,1)<3,
    warging('MATLAB:vbm_vol_average:input','WARNING:vbm_vol_average:to less input for median and std (n=2)!\n');
    dt(1)=0; dt(1)=0;
  end
  
  
  if numel(dt)==1; dt=repmat(dt,1,3); end
  
  if exist('PT','var') && ~isempty(PT) && exist(PT,'file')
  % reslicing
    if iscell(PT) && size(PT,1)<size(PT,2), PT=PT'; end
    PT = char(PT);
    if ~exist(PT,'file')
      error('ERROR:vbm_vol_median:PT','ERROR:vbm_vol_median:PT-file ''%s'' does not exist.',PT);
    end
    
    para.reslice.mask   = 0;
    para.reslice.mean   = 0;
    para.reslice.interp = 1;
    para.reslice.which  = 1;
    para.reslice.wrap   = [0 0 0];
    para.reslice.prefix = 'vbmvolmedian_r';
    
    spm_reslice([cellstr(PT);cellstr(P)],para.reslice);
    
    PR=cell(size(P,1),1);
    for fi = 1:size(P,1)
      [pp,ff,ee] = spm_fileparts(P(fi,:));
      PR{fi} = fullfile(pp,[para.reslice.prefix  ff ee]);
    end
    P = char(PR); clear PR;
    vbm_spm_reprint; 
  end
  
  PO = {};
  for i=1:3
    if dt(i)>0
      % median

      
      V  = spm_vol(P); if exist(fname2{i},'file'), delete(fname2{i}); end
      VO1 = V(1); VO1.fname = fname2{i}; VO1.dt(1) = 16;
      switch i
        case 1, VO1.descript = sprintf('median image of %s scans',size(P,1));
        case 2, VO2.descript = sprintf('mean image of %s scans',size(P,1));
        case 3, VO3.descript = sprintf('std image of %s scans',size(P,1));
      end
      
      VO1 = spm_create_vol(VO1);
      Y  = zeros([VO1.dim(1:2),1,size(P,1)],'single');
      for p=1:V(1).dim(3)
        for fi = 1:size(P,1), 
           Y(:,:,1,fi) = single(spm_slice_vol(V(fi),spm_matrix([0 0 p]),V(fi).dim(1:2),0));
        end
        switch i
          case 1, YO = nanmedian(Y,4);
          case 2, YO = nanmean(Y,4);
          case 3, YO = nanstd(Y,1,4);
        end
        VO1 = spm_write_plane(VO1,YO,p);
      end
      
      if exist('dt','var')
        Y = spm_read_vols(VO1);
        VO1.dt(1) = dt(i); 
        if dt(i)==2 || dt(i)==4
          VO1.pinfo(1) = round(max(Y(:))/(16^dt(i)-1));
        end
        VO1 = rmfield(VO1,'private');
        if exist(fname2{i},'file'), delete(fname2{i}); end
        
        spm_write_vol(VO1,Y);
      end  
      PO=[PO,fname2(i)]; %#ok<AGROW>
    end
  end
  if nargout==1, VO = spm_read_vols(PO); end
  
  if exist('PT','file')
    for fi = 1:numel(P)
      delete(P(fi,:))
    end
  end
end
function vbm_spm_reprint(str,lines)
  if ~exist('str','var'), str = ''; end
  if ~exist('lines','var'), lines=3; end
  if lines>0
    fprintf(sprintf('%s',repmat('\b',1,lines*73+1)));
  else
    fprintf(sprintf('%s',repmat('\b',1,-lines)));
  end
  fprintf(str);
end

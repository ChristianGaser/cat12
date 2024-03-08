function varagout = cat_vol_mimcalc(job)
% Just a small batch to run imcalc for multiple subject in the same,
% e.g., to apply some function or general correction. 
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id: cat_main_gintnormi.m 1834 2021-05-28 14:45:20Z dahnke $

  SVNid = '$Rev: 1901 $';
  
  def.var     = {};
  def.verb    = 1;
  def.cleanup = 0;
  job = cat_io_checkinopt(job,def);

  BIDSdirname = ['derivatives' filesep 'mimcalc'];
  varagout = {};
  
  % spm banner
  if job.verb
    spm('FnBanner',mfilename,SVNid);
  end

  for ri = 1:numel(job.images)
    for si = 1:numel(job.images{ri})
      [a,b,c,d] = checkBIDS(job.images{ri}(si),BIDSdirname);
      sfiles{si}{ri}=a; sfilesBIDS{si}(ri)=b; BIDSsub{si}{ri}=c; outdir{si}{ri} = d{1};
    end
  end
  %%
    
  % if isempty(job.outdir)
  %   outdir = spm_fileparts(job.images{si}{1}); % home of the first image type
  % else
  % end
   
  for si = 1:numel(job.images{1})    
    job2 = rmfield(job,{'prefix','images'}); 
    Pcleanup = {};
   
    % handle zipped BIDS 
    for ri = 1:numel(job.images)
      [~,ff,ee] = spm_fileparts(job.images{ri}{si});
      if strcmp(ee,'.gz')
        if ~exist(fullfile( outdir{si}{ri} , ff),'file')
          gunzip(job.images{ri}{si},outdir{si}{ri});
        end
        job2.input{ri} = fullfile( outdir{si}{ri} , ff);  
      else
        if ~strcmp( spm_fileparts(job.images{ri}{si}) , outdir{si}{ri} )
          copyfile(job.images{ri}{si},outdir{si}{ri});
        end
        job2.input{ri} = fullfile( outdir{si}{ri} , [ff ee]); 
      end
      Pcleanup{ri} = job2.input{ri};
    end
    
    [~,ff,ee]      = spm_fileparts(job2.input{1});
    job2.output    = [job.prefix,ff,ee];
    job2.outdir{1} = outdir{si}{ri}; 
    
    % call
    out = my_spm_imcalc(job2);
    
    % display progress file 
    vararout{1}.Pname{si} = out;
  end
  
  % remove temporary files in case of BIDS
  if job.cleanup
    for fi = 1:numel(Pcleanup)
      if exist(Pcleanup{fi},'file'), delete(Pcleanup{fi}); end
    end
  end

  vararout{1}.Pname(cellfun('isempty',vararout{1}.Pname)==1) = [];  
  
  
  % spm banner
  if job.verb
    fprintf('done.\n');
    spm_progress_bar('Clear')
  end
end
function out = my_spm_imcalc(job)
  [p,nam,ext] = spm_fileparts(job.output);
  if isempty(p)
      if isempty(job.outdir{1})
          p = pwd;
      else
          p = job.outdir{1};
      end
  end
  if isempty(nam)
      nam = [spm_get_defaults('imcalc.prefix') spm_file(job.input{1},'basename')];
      ext = ['.' spm_file(job.input{1},'ext')];
  end
  if isempty(ext)
      ext = spm_file_ext;
  end
  out.files = { fullfile(p,[nam ext]) };
  extra_vars = {};
  if numel(job.var)
      extra_vars = { job.var };
  end
  options = {};
  if isfield(job,'options')
    options = job.options; 
  end

  switch lower(job.expression)
    case 'approx'
      if numel(job.input)>1
        error('cat_vol_mimcalc:approx:filenumber','Only one file per subject allowed')
      end
      
      V  = spm_vol(char(job.input));
      Y  = spm_read_vols(V); 
      Y2 = cat_vol_approx(Y); 
      V2 = V; V2.fname = out.files{1}; 
      spm_write_vol(V2,Y2); 
     
      cmd = 'spm_image(''display'',''%s'')';
      fprintf('ImCalc-approx Image: %s\n',spm_file(out.files{1},'link',cmd));
    case {'msk-gm','msk-wm','msk-csf','msk-brain','mask-gm','mask-wm','mask-brain',...
        'brainmask','brainmsk','gm-mask','wm-mask','csf-msk'}
      if numel(job.input)>1
        error('cat_vol_mimcalc:msk:filenumber','Only one file per subject allowed')
      end

      if ~isempty(strfind(lower(job.expression),'gm') )
        Pmsk = cellstr([char(cat_get_defaults('extopts.shootingtpm')),',1']);
      elseif ~isempty(strfind(lower(job.expression),'wm') )
        Pmsk = cellstr([char(cat_get_defaults('extopts.shootingtpm')),',2']);
      elseif ~isempty(strfind(lower(job.expression),'brain') )
        Pmsk = cat_get_defaults('extopts.brainmask');
      else
        error('cat_vol_mimcalc:msk:case','Unknown masking case (use "msk-gm","msk-wm","msk-brain").')
      end

      cat_vol_imcalc(char( [ cellstr(job.input) ; Pmsk] ), out.files{1}, 'i1 .* (i2>0.1)');
     
      cmd = 'spm_image(''display'',''%s'')';
      fprintf('ImCalc-mask Image: %s\n',spm_file(out.files{1},'link',cmd));
      
    otherwise
      try
          cat_vol_imcalc(char(job.input), out.files{1}, job.expression, options, extra_vars{:});

          cmd = 'spm_image(''display'',''%s'')';
          fprintf('ImCalc Image: %s\n',spm_file(out.files{1},'link',cmd));
      catch
          cat_io_cprintf('err',sprintf('ImCalc Image: %s failed \n',out.files{1}));
          out.files{1} = '';
      end    
  end
end
function [sfiles,sfilesBIDS,BIDSsub,devdir] = checkBIDS(sfiles,BIDsdirname) 
  sfilesBIDS  = false(size(sfiles)); 
  BIDSsub     = ''; 
  devdir      = cell(size(sfiles));

  % if BIDS structure is detectected than use only the anat directory 
  for sfi = numel(sfiles):-1:1
    % detect BIDS directories 
    sdirs = strsplit(sfiles{sfi},filesep); 
    if strcmpi(sdirs{end-1}(1:min(4,numel(sdirs{end-1}))),'anat'), ana = 1; else, ana = 0; end
    if strcmpi(sdirs{end-2}(1:min(4,numel(sdirs{end-2}))),'ses-'), ses = 1; else, ses = 0; end
    if ses==0
      if strcmpi(sdirs{end-2}(1:min(4,numel(sdirs{end-2}))),'sub-'), sub = 1; else, sub = 0; end
    else
      if strcmpi(sdirs{end-3}(1:min(4,numel(sdirs{end-3}))),'sub-'), sub = 1; else, sub = 0; end
    end
    
    % differentiate between cross and long cases
    if  ses && sub && ~ana, sfiles(sfi) = []; end
    if  ses && sub, BIDSsub = sdirs{end-3}; devi = numel(sdirs)-3; end % long
    if ~ses && sub, BIDSsub = sdirs{end-2}; devi = numel(sdirs)-2; end % cross
    sfilesBIDS(sfi) = sub; 
    
    % setup result directory - without BIDS the default is used
    devdir{sfi}     = '';
    for di = 2:numel(sdirs)-1
      if sub && di == devi, devdir{sfi} = [devdir{sfi} filesep BIDsdirname]; end % add some directories inbetween
      devdir{sfi} = [devdir{sfi} filesep sdirs{di}]; 
    end
  end
end
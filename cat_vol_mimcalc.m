function out = cat_vol_mimcalc(job)
% Just a small batch to run imcalc for multiple subject in the same,
% e.g., to apply some function or general correction. 
% Extended by coregistrations and BIDS options.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id: cat_main_gintnormi.m 1834 2021-05-28 14:45:20Z dahnke $

  SVNid = '$Rev: 1901 $';
  
  def.var        = {};
  def.verb       = 1;
  def.cleanup    = 1;
  def.outdir     = {};
  def.BIDSdir    = 'derivatives/mimcalc';
  def.prefix     = ''; 
  def.suffix     = '';
  def.options.coreg = 0; 
  
  job = cat_io_checkinopt(job,def);

  out = struct();
  
  % spm banner
  if job.verb
    spm('FnBanner',mfilename,SVNid);
  end

  % prepare file and directory handling
  outdir = cell(size(job.images)); 
  for ri = 1:numel(job.images)
    for si = 1:numel(job.images{ri})
      [~,~,~,d] = checkBIDS(job.images{ri}(si),job.BIDSdir);
      outdir{si}{ri} = d{1};
      if ~isempty(job.outdir) && iscell( job.outdir ) && ~isempty( job.outdir{1} )
        fx = strfind(outdir{si}{ri},job.BIDSdir); 
        if ~isempty(fx) && cat_io_contains(job.BIDSdir,'derivatives') % combine out and BIDsdir
          outdir{si}{ri} = fullfile( job.outdir{1} , outdir{si}{ri}(fx(end):end) ); 
        elseif ~cat_io_contains(job.BIDSdir,'derivatives') % do not add subdirs
          outdir{si}{ri} = fullfile(job.outdir{1},job.BIDSdir); 
        else
          outdir{si}{ri} = job.outdir{1}; 
        end
      end
      if ~exist(outdir{si}{ri},'dir'), mkdir(outdir{si}{ri}); end
    end
  end
       
  for si = 1:numel(job.images{1})    
    job2 = rmfield(job,{'prefix','images','suffix','outdir'}); 
    Pcleanup = cell(size(job.images));
   
    % handle zipped BIDS 
    for ri = 1:numel(job.images)
      [pp,ff,ee] = spm_fileparts(job.images{ri}{si});
      job.images{ri}{si} = fullfile(pp,[ff ee]); 
      cleanup = 0;
      if strcmp(ee,'.gz')
        if ~exist(fullfile( outdir{si}{ri} , ff),'file')
          gunzip(job.images{ri}{si},outdir{si}{ri});
          cleanup = 1; 
        end
        job2.input{ri} = fullfile( outdir{si}{ri} , ff);  
      else
        if ~strcmp( spm_fileparts(job.images{ri}{si}) , outdir{si}{ri} ) && ...
          ~exist( spm_file( job.images{ri}{si},'path',outdir{si}{ri}) ,'file')
          copyfile(job.images{ri}{si},outdir{si}{ri}); 
        end
        job2.input{ri} = fullfile( outdir{si}{ri} , [ff ee]); 
      end
      if cleanup
        Pcleanup{ri} = job2.input{ri};
      end
    end
    
    % filenameparts
    [~,ff,ee]      = spm_fileparts(job2.input{1});

    % update prefix
    if ~isempty( strfind(job.prefix,'\f') )
      ff = ff(1+numel(strfind(job.prefix,'\f')):end); 
      prefix = strrep(job.prefix,'\f',''); 
    else
      prefix = job.prefix; 
    end
    if ~isempty( strfind(job.prefix,'\b') )
      ff = ff(1:end - numel(strfind(prefix,'\b'))); 
      prefix = strrep(prefix,'\b',''); 
    end
    
    % update suffix
    if ~isempty( strfind(job.suffix,'\b') )
      ff = ff(1:end - numel(strfind(job.suffix,'\b'))); 
      suffix = strrep(job.suffix,'\b',''); 
    else 
      suffix = job.suffix; 
    end
    if ~isempty( strfind(job.suffix,'\f') )
      ff = ff(1+numel(strfind(job.suffix,'\f')):end); 
      suffix = strrep(suffix,'\f',''); 
    end

    job2.output    = [prefix,ff,suffix,ee];
    job2.outdir{1} = outdir{si}{1}; 
    
    % call
    out1 = my_spm_imcalc(job2);
    
    % display progress file 
    out.Pname{si,1} = out1.files{1};
  
    % remove temporary files in case of BIDS
    if job.cleanup
      for ri = 1:numel(Pcleanup)
        if ~isempty(Pcleanup{ri}) && exist(Pcleanup{ri},'file'), delete(Pcleanup{ri}); end
      end
    end
  end

  out.Pname(cellfun('isempty',out.Pname)==1) = [];  
  
  % spm banner
  if job.verb
    spm_progress_bar('Clear')
  end
end
function out = my_spm_imcalc(job)
  [p,nam,ext] = spm_fileparts(job.output);
  if isempty(p)
      if isempty(job.outdir) || (iscell(job.outdir) && isempty(job.outdir{1}))
          p = pwd;
      elseif iscell(job.outdir)
          p = job.outdir{1};
      else
          p = job.outdir;
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

  if job.options.coreg
  % coregistration to support simple processing of T1w./T2w 
    V  = spm_vol(char(job.input)); Vr = V; 
    for vi = 2:numel(V) 
      evalc( sprintf(['cormats{vi} = spm_coreg( V(1) , V(vi) , ' ...
                 'struct( ''sep'' , [8 4 2 1] , ''fwhm'' , [7 7] , ' ...
                '''graphics'' , 0 ) );']) ); 
      Y = spm_read_vols(Vr(vi)); 
      Vr(vi).fname = spm_file( Vr(vi).fname , 'path', job.outdir,'prefix','r'); 
      Vr(vi).mat   = spm_matrix(cormats{vi}) \ eval(sprintf('V(vi).mat')); %#ok<USENS>
      spm_write_vol(Vr(vi), Y); 
      job.input{vi} = Vr(vi).fname;
    end
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

      if cat_io_contains( lower(job.expression),'gm')
        Pmsk = cellstr([char(cat_get_defaults('extopts.shootingtpm')),',1']);
      elseif cat_io_contains( lower(job.expression),'wm')
        Pmsk = cellstr([char(cat_get_defaults('extopts.shootingtpm')),',2']);
      elseif cat_io_contains( lower(job.expression),'brain')
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

  if job.options.coreg
    for vi = 2:numel(V) 
      if exist(Vr(vi).fname,'file')
        delete(Vr(vi).fname)
      end
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
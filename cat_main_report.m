function out = cat_main_report(job)
%cat_main_report. Retrospective creation of the CAT report
% 
%  out = cat_main_report(job)
% 
%  job
%   .files    .. xml-files for CAT report generation
%   .Pm       .. use other input for bias corrected files
%   .Pp0      .. use other input for segmentation files
%   .print    .. print setting
%   .outdir   .. use different output directory (if empty print 
%                at the original path)
% 
% Examples: 
%  cat_main_report
%  cat_main_report(struct('outdir',pwd));   
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  if ~exist('job','var') || ~isfield(job,'files') || isempty(job.files) || isempty(job.files{1})
    job.files = cellstr(spm_select([1 inf],'any','Select CAT XML files for Reports:',{},pwd,'^cat_.*\.xml$'));
  end
  if isempty(job.files) || ischar(job.files) || isempty(job.files{1}) , out = struct(); return; end

  def.files     = {};   % cat_ (or catlong_) XML input files
  def.Pm        = {};   % bias corrected images to replace default setting
  def.Pp0       = {};   % segmentation label maps to replace default setting
                        % surfaces? 
  def.print     = 2;    % render options for develpers?
  def.outdir    = {};
  
  job = cat_io_checkinopt(job,def);
  job.outdir = cellstr(job.outdir);

  % check size Fm
  if ~isempty(job.Pm) && iscell(job.Pm) && ~isempty(job.Pm{1}) && (numel(job.files) ~= numel(job.Pm))
    error(sprintf('%s:incorrectInputPm',mfilename), ...
      'The number of xml and bias corrected images has to be equal.');
  end

  % check size Fp0 
  if ~isempty(job.Pp0) && iscell(job.Pp0) && ~isempty(job.Pp0{1}) && (numel(job.files) ~= numel(job.Pp0))
    error(sprintf('%s:incorrectInputPp0',mfilename), ...
      'The number of xml and label maps has to be equal.');
  end

  out.files = {}; 
  out.filesj = {}; 
  for fi = 1:numel(job.files)

    % read xml; 
    xml = cat_io_xml(job.files{fi});

    if isfield(xml,'error') %|| ~isfield(xml,'SPMpreprocessing')
      cat_io_cprintf('err',sprintf('  Error "%s" was not processed correctly.\n', job.files{fi})); 
      continue; 
    end
    
    % load Ym
    if isfield( job, 'Pm' ) && iscell(job.Pm) && ~isempty( job.Pm ) && ~isempty( job.Pm{fi} )  && exist( job.Pm{fi} , 'file' )
      Vmi = spm_vol( job.Pm{fi} ); 
    elseif exist( xml.filedata.Fm , 'file' )
      Vmi = spm_vol( xml.filedata.Fm  ); 
    elseif exist( xml.filedata.F , 'file' )
      Vmi = spm_vol( xml.filedata.F );
    else
      [pp,ff] = spm_fileparts(job.files{fi});
      ppm = [pp(1:end-5) strrep(pp(end-5:end),'report','mri')];
      ffm = ['m' ff(5:end) '.nii'];
      Pmi = fullfile(ppm,ffm); 
      if exist(Pmi,'file')
        Vmi = spm_vol( Pmi );
      else
        ppm = [pp(1:end-6) strrep(pp(end-5:end),'report','')];
        ffm = [ff(5:end) '.nii'];
        Pmi = fullfile(ppm,ffm); 
        if exist(Pmi,'file')
          Vmi = spm_vol( Pmi ); 
        else
          Vmi = []; 
        end
      end
    end

    % load Yp0 
    if isfield( job, 'Pp0' ) && iscell(job.Pp0) && ~isempty( job.Pp0 ) && ~isempty( job.Pp0{fi} ) && exist( job.Pp0{fi} , 'file' )
      Vp0 = spm_vol( job.Pp0{fi} );
    elseif exist( xml.filedata.Fp0 , 'file' )
      Vp0 = spm_vol( xml.filedata.Fp0 );
    else
      [pp,ff] = spm_fileparts(job.files{fi});
      ppm = [pp(1:end-6) strrep(pp(end-5:end),'report','mri')];
      ffm = ['p0' ff(5:end) '.nii'];
      Pp0 = fullfile(ppm,ffm); 
      if exist(Pp0,'file')
        Vp0 = spm_vol( Pp0 );
      else
        Vp0 = []; 
      end
    end

    if ~isempty( Vmi )
      Ymi = spm_read_vols( Vmi );
    else
      Ymi = [];
    end
    if ~isempty( Vp0 )
      Yp0 = spm_read_vols( Vp0 );
    else
      % try to locate the segmentation file p0*.nii 
      if isfield(xml,'filedata')
        Pp0 = xml.filedata.Fp0;
      end
      if ~exist(Pp0,'file')
        Pp0a = strrep(Pp0,[filesep 'mri' filesep],filesep); 
        if exist(Pp0a,'file'), Pp0 = Pp0a; end
      end
      % if found load it otherwise keep it empty
      if exist( Pp0 ,'file')
        Vp0 = spm_vol( Pp0 ); 
        Yp0 = spm_read_vols( Vp0 );
      else
        Yp0 = []; 
      end
    end

    %% create if Yy is avaiable 
    Pc = spm_file( xml.filedata.F, 'prefix', 'cat_', 'path', spm_fileparts(xml.filedata.Fp0) );
    Py = spm_file( xml.filedata.F, 'prefix', 'y_'  , 'path', spm_fileparts(xml.filedata.Fp0) );
    if exist( Pc , 'file' )
      Vl1 = spm_vol( Pc );
      Yl1 = spm_read_vols( Vl1 );
    elseif 0 % exist( Py , 'file' ) ... to create the segmentation overlay not working yet
      Vy  = spm_vol( Py );
      for di = 1:3
        Yy(:,:,:,di) = single( Vy.private.dat(:,:,:,1,di) );
      end
      
      % this is not working (orientation incorrect) 
      Yyi = spm_diffeo('invdef', Yy , Vp0.dim , eye(4), eye(4)); %Vy.mat, Vp0.mat ); 
      Vli = spm_vol( xml.parameter.extopts.cat12atlas{1} );
      Yli = cat_vol_ctype( cat_vol_sample(Vp0,Vli,Yy,0) );
    else
      Yl1 = Yp0>0; 
    end

    if isfield(xml.subjectmeasures,'dist_thickness')
      %%
      clear Psurf; 
      surfdir = strrep( spm_fileparts( xml.filedata.Fp0 ) , [filesep 'mri'] , [filesep 'surf'] );
      side = {'lh','rh'};
      for si = 1:2
        Psurf(si) = struct( ...
          'Pcentral', spm_file( xml.filedata.F, 'prefix', sprintf('%s.%s.',side{si},'central')  , 'path', surfdir ,'ext','.gii'), ...
          'Ppial'   , spm_file( xml.filedata.F, 'prefix', sprintf('%s.%s.',side{si},'pial')     , 'path', surfdir ,'ext','.gii'), ...
          'Pwhite'  , spm_file( xml.filedata.F, 'prefix', sprintf('%s.%s.',side{si},'white')    , 'path', surfdir ,'ext','.gii'), ...
          'Pthick'  , spm_file( xml.filedata.F, 'prefix', sprintf('%s.%s.',side{si},'thickness'), 'path', surfdir ,'ext','') );
      end
      if ~exist(Psurf(1).Pcentral,'file') % this files are typically deleted 
        Psurf = '';
      end

    else
      Psurf = {};
    end
    if ~exist('Psurf','var'), Psurf = ''; end


    % prepare CAT report variables as good as possible
    xmljob          = xml.parameter;
    xmljob.files    = job.files(fi);  
    xmljob.extopts.print    = job.print;
    xmljob.output.surface   = ~isempty(Psurf);
    if ~isempty( job.outdir ) && ~isempty( job.outdir{1} )
      [~,ff] = spm_fileparts( job.files{fi} );
      ff = strrep(strrep(ff,'catlong_','catlongreport'),'cat_','catreport_'); 
      xmljob.imgprint.fname  = fullfile(job.outdir{1}, [ff '.pdf']);
      xmljob.imgprint.fnamej = fullfile(job.outdir{1}, [ff '.jpg']);
    end
 
    % QA
    xmlqa             = xml;   
    
    % SPM preprocessing parameters
    if isfield( xml , 'SPMpreprocessing' )
      xmlres          = xml.SPMpreprocessing; 
    else
      xmlres          = struct('setCOM',NaN);
    end
    if isfield(xml,'ppe')
      xmlres.ppe      = xml.ppe; 
    else
      xmlres.ppe      = struct(); 
    end
    xmlres.do_dartel  = 0; 
    xmlres.stime      = clock; 
    xmlres.image      = Vmi;
    xmlres.image0     = Vmi; 
    if isfield(xml.parameter,'vbm') && ~isfield(xml.parameter,'extopts')
      xml.parameter.extopts = xml.parameter.vbm.extopts;
    end
    xmlres.bb         = [-72 -108 -72; 72 72 72] + [-1 -1 -1; 1 1 1] .* xml.parameter.extopts.bb;

    % prepare str for parameters 
    str = cat_main_reportstr(xmljob, xmlres, xmlqa);

    % print reprot
    cat_main_reportfig(Ymi, Yp0, [], Psurf, xmljob, xmlqa, xmlres, str);

    % set dependency
    if ~isempty( job.outdir ) && ~isempty( job.outdir{1} )
      out.files{fi}  = xmljob.imgprint.fname; 
      out.filesj{fi} = xmljob.imgprint.fnamej;
    else
      out.files{fi} = spm_file( xml.filedata.F, 'prefix', 'catreport_' , ....
        'path', strrep( spm_fileparts( xml.filedata.Fp0 ), [filesep 'mri'] , [filesep 'report'] ), ...
        'ext','.pdf');
      out.filesj{fi} = spm_file( xml.filedata.F, 'prefix', 'catreportj_' , ....
        'path', strrep( spm_fileparts( xml.filedata.Fp0 ), [filesep 'mri'] , [filesep 'report'] ), ...
        'ext','.jpg');      
    end
  end

  % remove unprocessed cases
  out.files(  cellfun(@isempty,out.files  ) ) = [];
  out.filesj( cellfun(@isempty,out.filesj ) ) = [];
end




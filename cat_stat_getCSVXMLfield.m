function out = cat_stat_getCSVXMLfield(job)
%cat_stat_getCSVXMLfield. Extract file/line specific fields from CSV/XML
% This function extract data from given CSV and XML files. 
%
%  out = cat_stat_getCSVXMLfield(job)
%
%  job         .. matlabjob structure
%   .files     .. XML files (cellstr)
%   .csvfile   .. CSV file (cellstr)
%   .xmlsets   .. predefined XML fields (default = 'full'; 
%                  'default', 'expert', 'developer', 'full')
%   .fields    .. set of searchstring for field extraction
%   .csvIDcol  .. column that defines the fields in the CSV file (def. = 1)
%
%   .csvIDfd   .. field that defines the dataIDs (eg. subject)
%   .outdir    .. output directory (default = '', i.e., current directory)
%   .write     .. datatype to export, e.g. {'mat','csv','tsv','txt'}
%   .fname     .. output file name, default '<fn>_nf<nf>_nl<nl>', where 
%                 <fn>, <nf>, and <nl> are replaced by filename, number of 
%                 fields (columns) and lines.
%   .verb      .. be verbose  
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$ 

%
% TODO
% == XML == 
% * multi entry probem, fields with more than 1 element (in case of char lines)
%   - creation of subentries *_1 ... *_n in case of less elements 
%     >> parameter for limit? with default = 9?
%
% * check of known XML types?
%   - cat_* vs. catlong_* vs. catROI*
%   - mixing/merging of types?
%     > merging (adding fields) >> independend second batch
%     > mixing (adding lines) >> not useful fil there is no overlap
%     > would require to analyse all XML fields
%       >> create error while processing and exclude files
%
% * test of correct processing
%   - cat_        .. ok
%   - catlong_    .. ok but only parameters
%   - catROI(s)_  .. NOT OK only first entry, char converting problem!
%   - test of missing fields 
%   - missing information? >> missing entires? 
%
% == CSV ==
% * field search with columnsID ?
% 
% == CSV/XML == 
% * ... 
%
% == general
% * CODING of variables 
%
% * documentation
% * tests
% * cleanup
%    


%#ok<*ASGLU,*AGROW>
  
  def.files     = {''}; % n-files, e.g. XML for direct extraction or nii/gii as selector
  def.csvfile   = {''}; % 0..1-file ... maybe n later
  def.xmlfields = {};   % manuel predefinition of XML field 
  def.label     = {};   % manuel predefinition of CAT ROI labels ... currently no used
  def.labels    = {};   % manuel predefinition of CAT ROI labels ... currently no used
  def.xmlsets   = '';   %
  def.xmlmatlim = 9;    % limitation to split of matrix elements
  def.fnamefields         = struct('');
  def.idselector.csvIDfd  = [1 1];
  def.idselector.fileseps = '_-.';
  def.fields    = {};   % set of variables names for extraction ... preselection TIV IQR ...
                        % the variables were extracted and a depency for each created
  def.seg       = ',.'; % delimiter/komma setting {',.',';,',';.'} .. for export
  def.csvIDcol  = 1;
  def.csvIDfd   = 0;    % defines compontent of the pathname to match the subject
  def.outdir    = '';   % output directory
  def.write     = {'mat','csv','tsv'}; 
  def.fname     = '';   % write output files‚
  def.dep       = 0;    % internal variable to test dependencies
  def.verb      = 1; 
    
  job = cat_io_checkinopt(job,def); 
  out = struct(); 

  % simplification of some variables 
  job.csvfile = char(job.csvfile); 
  job.outdir  = char(job.outdir); 
  %if numel(job.csvIDfd)==1, job.csvIDfd(2) = job.csvIDfd; end % currently not use 
   
  if job.dep 
    % prepare dependencies in matlab batch mode and handle empty inputs 
    job.files = job.files(1); 
    job.verb  = 0; 
    if ... (isempty(job.xmlsets) && isempty(job.xmlfields)) || ... isempty(job.fields) && 
      isempty(job.files) || strcmp(char(job.files),'<UNDEFINED>') || strcmp(char(job.files),'<')
      return
    end
  end 
  

  
  %% try to find XML files
  %  ----------------------------------------------------------------------
  %  The original goal was to use also other nifti/gifti input generally
  %  used for statistical analysis.
  if ~isempty( job.files{1} )
    [pp,ff,ee] = spm_fileparts(job.files{1});  
    [~,fname]  = spm_fileparts(job.files{1}); 
    Pxml = job.files; 

    % get original file 
    for xi = 1:numel(Pxml)
      [pp,ff,ee] = spm_fileparts(Pxml{xi}); 
      Pnii{xi} = fullfile(pp,[ff(5:end) ee]);
    end
  else
    Pxml  = {};
    Pnii  = {};
    fname = ''; 
  end
    
  % check XML files
  if ~isempty(Pxml)
    for fi = 1:numel(Pxml)   
      if ~exist(Pxml{fi},'file') 
        if exist('sinfo','var')
          Pxml{fi} = fullfile( sinfo(fi).pp , [ 'cat_' sinfo(fi).name '.xml' ]); % same dir? 
        else
          [pp,ff,ee] = spm_fileparts(job.files{fi});  
          Pxml{fi} = fullfile( pp , [ 'cat_' sinfo(fi).name '.xml' ]); % same dir? 
        end
        if ~exist(Pxml{fi},'file')
          Pxml{fi} = fullfile( fileparts(sinfo(fi).pp) , [ 'cat_' sinfo(fi).name '.xml' ]); 
          if ~exist(Pxml{fi},'file')
            cat_io_cprintf('warn','  ERROR: Cannot find XML input %d: "%s".\n',fi,Pxml{fi})
          end
        end
      end
    end
  end

  Pxmlff = spm_str_manip(Pxml,'tr');
  if any( cat_io_contains(Pxmlff,'cat_') ) && any( cat_io_contains(Pxmlff,'catlong_') ) || ...
     any( cat_io_contains(Pxmlff,'cat_') ) && any( cat_io_contains(Pxmlff,'catROI') )   
    error(sprintf('%s:mixedCATxmls',mfilename),'Mixxed CAT XML files are not supported.'); 
  end

  % check volume label files if ROIs where selected
  %{
  if ~isempty( job.label )
    Plabel  = cellstr(spm_file(Pnii,'prefix',['label' filesep 'catROI_'] ,'ext','xml'));
    for fi = 1:numel(Plabel)   
      if ~exist(Plabel{fi},'file') 
        Plabel{fi} = cellstr(spm_file(Pnii,'prefix','catROI_' ,'ext','xml'));
        if ~exist(Plabel{fi},'file') 
          Plabel{fi} = ''; 
        end
      end
    end
  end
  %}



    
  
  % create useful fieldnames of the given (search) fields (both XML as CSV)
  replaceset = {
      {' ','\t',',','','&'}         '_' ;         
      {'(',')','{','}','[',']','='}  ''  ; 
    };
  if ~isempty( job.fields ) && ~isempty( job.fields{1} )
    fields = genvarname(cat_io_strrep(cat_io_strrep(cat_io_strrep( job.fields ,replaceset{1,1},replaceset{1,2}),replaceset{2,1},replaceset{2,2}),'__','_'));
    fields = cat_io_strrep(fields,{'x0x2E','0x2E'},{'.','.'});
    for si = numel(fields):-1:1
      if ~isnan( str2double(job.fields{si}) ), fields(si) = []; end
    end 
  else
    fields = {};
  end




  % (predefined) XML fields
  % -----------------------------------------------------------------------
  if ~isempty( job.xmlsets ) && ~isempty( job.files ) && ~isempty( job.files{1} )
    % field can be listed multiple times as we use unqiue at the end
    xmlfn0 = {
      ... subjectmeasures
       'subjectmeasures.vol_TIV'
       'subjectmeasures.surf_TSA'
       'subjectmeasures.vol_abs_WMH'
       ... quality rating
       'qualityratings.IQR'
       ... software 
       'software.version_cat'
       'software.revision_cat'
       'software.date'
       };
    xmlfn1 = {
        ... subjectmeasures
        'subjectmeasures.vol_TIV'
        'subjectmeasures.surf_TSA'
        'subjectmeasures.dist_thickness{1}(1)'
        'subjectmeasures.dist_thickness{1}(2)'
        'subjectmeasures.vol_rel_CGW(1)'
        'subjectmeasures.vol_rel_CGW(2)'
        'subjectmeasures.vol_rel_CGW(3)'
        'subjectmeasures.vol_rel_CGW(4)'
        'subjectmeasures.vol_abs_CGW(1)'
        'subjectmeasures.vol_abs_CGW(2)'
        'subjectmeasures.vol_abs_CGW(3)'
        'subjectmeasures.vol_abs_CGW(4)'
        ... quality ratings
        'qualityratings.IQR'
        'qualityratings.NCR'
        'qualityratings.ICR'
        'qualityratings.res_RMS'
        'qualityratings.contrast'
        ... quality measures
        'qualitymeasures.NCR'
        'qualitymeasures.ICR'
        'qualitymeasures.res_RMS'
        'qualitymeasures.contrast'
        ... quality surface measures 
        'qualitymeasures.SurfaceEulerNumber'
        'qualitymeasures.SurfaceDefectArea'
        'qualitymeasures.SurfaceDefectNumber'
          };
    xmlfn2 = {
        ... quality surface measures 
        'qualitymeasures.SurfaceEulerNumber'
        'qualitymeasures.SurfaceDefectArea'
        'qualitymeasures.SurfaceDefectNumber'
        'qualitymeasures.SurfaceIntensityRMSE'
        'qualitymeasures.SurfacePositionRMSE'
        'qualitymeasures.SurfaceSelfIntersections'    
        ... preprocessing measures SPM 
        'SPMpreprocessing.ll'
        'SPMpreprocessing.mn(1)'
        'SPMpreprocessing.mn(2)'
        'SPMpreprocessing.mn(3)'
        'SPMpreprocessing.mn(4)'
        'ppe.SPMvols0(1)'
        'ppe.SPMvols0(2)'
        'ppe.SPMvols0(3)'
        'ppe.SPMvols0(4)'
        'ppe.SPMvols0(5)'
        'ppe.SPMvols0(6)'
        ...'ppe.skullstrippedpara'
        ...'ppe.highBGpara'
        ...reg.ll
        ...reg.dt, rmsdt
      };

%#################    
% ========= multiple entry export ... ROI label fiels ??? ======= ???
% ok for tables but not for dependencies 
%#################    
    % read first file to get the XML fields to create the dependency output
    % assumes that all files are equal!
    xml     = cat_io_xml( Pxml{1} );
    [xmlfull,xmlcell] = cat_io_struct2cell( xml );

    % match/select fields
    if ~isempty( job.fields ) && ~isempty( job.fields{1} ) %~isempty( fields ) 
      xmlfull = xmlfull( cat_io_contains( xmlfull, fields ) );
    end

    % define xmlfields depending on the sets 
    [pp,ff,ee] = spm_fileparts(job.files{1}); 
    if numel(ff) > 4 && strcmp( ff(1:4) , 'cat_' )
      switch job.xmlsets 
        case 'default'
          xmlfields = unique( xmlfn0 );
        case 'expert'
          xmlfields = unique( [ xmlfn0 ; xmlfn1 ]);
        case 'developer'
          xmlfields = unique( [ xmlfn0 ; xmlfn1 ; xmlfn2 ]);
        case 'full'
          xmlfields = xmlfull;
      end
    elseif ( numel(ff) > 7 && strcmp( ff(1:7) , 'catROI_'  ) ) || ... 
           ( numel(ff) > 8 && strcmp( ff(1:7) , 'catROIs_' ) )
      % ROI handling 
      % ... select ATLAS ?
      % ... select ROI ? 
      % ... no dependencies only tables ? 
      %    side merge?

    elseif numel(ff) > 8 && strcmp( ff(1:8) , 'catlong_' )
      xmlfields = xmlfull;
    else
      xmlfields = xmlfull;
    end
  else
    xmlfields = {};  
  end

  % match/select
  if ~isempty( job.fields ) && ~isempty( job.fields{1} )
    xmlfields = xmlfields( cat_io_contains( lower(xmlfields), lower(fields) ) );
  end

  if ~isempty( job.xmlfields )
    for fxi = 1:numel(job.xmlfields)
      FN = fieldnames( job.xmlfields{fxi} ); 
      for fni = 1:numel(FN)
        xmlfields = [xmlfields; job.xmlfields{fxi}.(FN{fni}) ];
      end
    end
  end



  

  %% get CSV data
  %  ----------------------------------------------------------------------
  % also as output list
  out = rmfield( out ,'files');  
  
  if ~isempty(job.csvfile) && ~strcmp(job.csvfile,'<UNDEFINED>') 
  % rename fielnames
    if ~exist(job.csvfile,'file')
      error('cat_stat_getCSVXMLfield:noCSVfile','ERROR: No CSV file at this path "%s".\n',job.csvfile); 
    end
    
    [~,fname] = spm_fileparts(job.csvfile);  

    % get table
    if isfield(job,'dep') && job.dep
      csv = cat_io_csv(job.csvfile,1); % read only header
    else
      csv = cat_io_csv(job.csvfile);
    end
    csvvars = genvarname(cat_io_strrep(cat_io_strrep(cat_io_strrep( csv(1,1:end) , ...
      replaceset{1,1},replaceset{1,2}),replaceset{2,1},replaceset{2,2}),'__','_'))';
  
    % get original ids
    if ~( isfield(job,'dep') && job.dep ) 
      csvids  = csv(2:end,job.csvIDcol);
      try
        lx      = max(length(num2str(max(cell2mat(csvids))))); 
        if isnumeric( cell2mat(csvids(1)) )
          for si = 1:numel(csvids)
            csvids{si} = num2str(csvids{si},sprintf('%%0%dd',lx)); 
          end
        end
      end
    end

    % remove lines
    if isfield(job,'csvselcol') && job.csvselcol > 0 && job.csvselcol < size(csv,2)
      % remove empty entries
      csvrm = cellfun(@isempty, csv(2:end,job.csvselcol) ); 
      csv([false;csvrm],:)  = [];
      csvids(csvrm)         = [];

      % remove NaN and entries
      if isnumeric( csv(2:end,job.csvselcol) )
        csvrm = isnan(cell2mat(csv(2:end,job.csvselcol)));
      else
        csvrm = cat_io_contains(lower(csv(2:end,job.csvselcol)),'n/a');
      end
      csv([false;csvrm],:)  = [];
      csvids(csvrm)         = [];
      
      % remove zeros 
      if isnumeric( csv(2:end,job.csvselcol) )
        csvrm = cell2mat(csv(2:end,job.csvselcol)) == 0;
        csv([false;csvrm],:)  = [];
        csvids(csvrm)         = [];
      end
    elseif isfield(job,'csvselcol') && job.csvselcol >= size(csv,2)
      error(sprintf('%s:csvselcol',mfilename),'There is no column number %d. ',job.csvselcol)
    end


    % match/select
    if ~isempty( job.fields ) && ~isempty( job.fields{1} )
      if isempty( fields ) && all(cellfun(@str2num, job.fields) < 0)
        % if only negative values are given than we need to init them
        csvsvari = 1:numel(csvvars);
      else 
        csvsvari = [];
        for si = 1:numel( job.fields )
          if ~isnan( str2double(job.fields{si}))  && str2double(job.fields{si})  > 0
            csvsvari = [ csvsvari ; si ];
          end
        end
      end
      csvsvari = [csvsvari find( cat_io_contains( lower(csvvars), lower(fields) ) )];
      for si = 1:numel( job.fields )
        if  ~isnan( str2double(job.fields{si}))  && str2double(job.fields{si}) < 0
          csvsvari = setdiff( csvsvari , si );
        end
      end
      csvsvars = csvvars( csvsvari );  
      csvsvarj = [];
    else
      csvsvars = csvvars; 
      csvsvari = 1:numel(csvvars);
      csvsvarj = [];
    end

    if job.verb 
      fprintf('Found %d CSV fields:  \n',numel(csvsvars)); 
      for fni=1:numel(csvsvars), cat_io_cprintf('b','% 6d) %s\n',fni,csvsvars{fni}); end
      fprintf('\n');
      if isempty(csvvars)
        fprintf('The CSV has the following fields (If it looks strange, check the deliminator):\n  '); 
        for fni=1:numel(csvvars), cat_io_cprintf('b','%s  ',csvvars{fni}); end
        fprintf('\n');
      end
    end
  else
    csvvars  = {}; 
    [csvsvars,csvsvari,csvsvarj] = intersect( csvvars , fields );
  end
  
  %%
  nsubs = numel(Pxml); 
  if ~isempty(job.fnamefields) && isfield( job.fnamefields(fi) , 'csvIDfd' ) && ~isempty(job.fnamefields(fi).csvIDfd)
    %%
    for fi = 1:numel(job.fnamefields)
      for si = 1:nsubs
        if exist('Pxml','var') && ~isempty(Pxml)
          pp = Pxml{si}; ff=''; P{si} = Pxml{si}; 
        end
        for pi = 0:job.fnamefields(fi).csvIDfd(1)
          [pp,ff0] = spm_fileparts(pp);
          if job.fnamefields(fi).csvIDfd(1)==pi || ( numel(job.fnamefields(fi).csvIDfd)>1 && pi>job.fnamefields(fi).csvIDfd(2) )
            if diff(job.fnamefields(fi).csvIDfd)
              ff = [ff0 filesep ff];   
            else
              ff = ff0; 
            end
          end
        end
        
        [pp,ff,ee] = spm_fileparts(ff); 
        pp  = cat_io_strrep(ff,num2cell(job.fnamefields(fi).fileseps),repmat({filesep},size(job.fnamefields(fi).fileseps))); 
        pps = textscan(pp,'%s','Delimiter',filesep); 
        if numel(job.fnamefields(fi).filesel)>1
          ff  = pps{1}{job.fnamefields(fi).filesel(1):job.fnamefields(fi).filesel(2)};
        elseif isscalar(job.fnamefields(fi).filesel)       
          try
            ff  = pps{1}{job.fnamefields(fi).filesel(1)};
          catch
            ppsfn = ''; for ppsi=1:numel(pps{1}), ppsfn = [ ppsfn '  ' pps{1}{ppsi} ]; end
            error('cat_stat_getCSVXMfield:fnamefield.filesel_too_large', sprintf( ...
              'Filename fieldselector is too small/large. \n  %d Fields: %s \n',numel(pps{1}),ppsfn)); 
          end
        end

        if ischar(ff) 
          out.(job.fnamefields(fi).name){si} = ff; 
        else 
          out.(job.fnamefields(fi).name)(si) = ff; 
        end
      end
    end
  elseif exist('csv','var') && size(csv,1) > 1
    nsubs = size(csv,1) - 1; 
  else
    nsubs = 0; 
  end
  

  

  
  if isfield(job,'dep') && job.dep
  % -----------------------------------------------------------------------
  % prepare output for dependencies
  % -----------------------------------------------------------------------

% ========== what in case of (XML) fields with mutliple elements? ==============  
    [xmlvars,xmlvarsi] = setdiff( xmlfields , csvsvars );
    xmlvarsf = cat_io_strrep(xmlvars,{'.','(',')','{','}'},{'_','','','',''}); 
    for fni = 1:numel(xmlvarsf)
      eval(sprintf('out.%s = nan;',xmlvarsf{fni})); 
    end 
    for fni = 1:numel(csvsvars)
      eval(sprintf('out.%s = nan;',csvsvars{fni})); 
    end

    return
  end
  
  
  %% CSV read out
  warns = {}; 
  for fi = 1:numel(csvsvars)
    for si = 1:nsubs
      %%
      if ~isempty( Pnii )
        if 0 %exist('Pxml','var') && ~isempty(Pxml) % not working yet
          pp = Pxml{si}; ff=''; P{si} = Pxml{si}; 
        elseif exist('Pnii','var') && ~isempty(Pnii)
          pp = Pnii{si}; ff=''; P{si} = Pnii{si}; 
        else 
          pp = ''; ff = ''; 
        end
        
        %%
        if ~isempty( job.idselector.csvIDfd ) || job.idselector.csvIDfd~=0
          for pi = 0:job.idselector.csvIDfd(1)
            [pp,ff0] = spm_fileparts(pp);
            if  job.idselector.csvIDfd(1)==pi || pi>job.idselector.csvIDfd(2)
              if diff(job.idselector.csvIDfd)
                ff = [ff0 filesep ff];   
              else
                ff = ff0; 
              end
            end
          end
          [pp,ff,ee] = spm_fileparts(ff); 
  
          pp  = cat_io_strrep(ff,num2cell(job.idselector.fileseps),repmat({filesep},size(job.idselector.fileseps))); 
          pps = textscan(pp,'%s','Delimiter',filesep); 
          if numel(job.idselector.filesel)>1
            ff  = pps{1}{job.idselector.filesel(1):job.idselector.filesel(2)};
          elseif isscalar( job.idselector.filesel )      
            ff  = pps{1}{job.idselector.filesel(1)};
          end
        end
      
        %%
        ids = []; idsl = [];
        for sii = 1:numel(csvids)
          if ~isempty( strfind( csvids{sii} ,ff ) ) ||  ~isempty( strfind( ff , csvids{sii} ) )
            ids = [ids sii]; idsl = [idsl length(csvids{sii})];
          end
        end
      else
        ids  = si;
        idsl = si;
      end
     
      %% remove double entries by bad ids (eg. the id="1" can also be found in "101" 
      ids( idsl < max(idsl))  = [];
      % idsl( idsl < max(idsl)) = [];
      if isempty(ids)
        cat_io_cprintf('err',sprintf('Cannot find data for "%s". You selected "%s" of the filename.\n',P{si},ff));
        out.(csvsvars{fi}){si,1} = nan; 
      elseif numel(ids)>1
        if isscalar( unique( [csv{ids + 1,1}] ) )
          warn = sprintf('      Found multiple possible entries for subject "%s". Take the first entry. Check your csv file! \n',P{si}); 
          out.(csvsvars{fi}){si,1} = csv{ids + 1,csvsvari(fi)}; 
        else
          warn = sprintf('      Found multiple possible entries (index=[ %s]) for subject in "%s". Check file!\n',sprintf('%d ',ids),P{si}); 
          out.(csvsvars{fi}){si,1} = nan;
        end
        if cat_io_contains(warns,warn) 
          cat_io_cprintf('warn',warn);
          warns = [warns warn]; 
        end
      else
        out.(csvsvars{fi}){si,1} = csv{ids + 1,csvsvari(fi)}; 
      end
      if ~isempty(ids) && fi == 1, out.ids{si,1} = csv{ids + 1}; end
    end    
    if  isnumeric( out.(csvsvars{fi}){1} )
      clear temp
      try 
        %error
        temp = cell2mat(out.(csvsvars{fi})); 
        out = rmfield(out,csvsvars{fi});
        out.(csvsvars{fi}) = temp;
      catch
        temp = zeros(numel(out.(csvsvars{fi})),1); 
        for si = 1:numel(out.(csvsvars{fi}))
          try 
            temp(si) = out.(csvsvars{fi}){si}; 
          catch
            try 
              temp(si) = str2double( out.(csvsvars{fi}){si} ); 
            catch
              temp(si) = nan; 
            end
          end
        end
        out = rmfield(out,csvsvars{fi});
        out.(csvsvars{fi}) = temp;
      end
    end
  end
  
  
  if ~isempty(Pxml)
    if job.verb, fprintf('  Read XML files .. '); end
    XML = cat_io_xml(Pxml); 

    % remove CSV fields
    [xmlvars,xmlvarsi] = setdiff( xmlfields , csvsvars );
    xmlvars0  = xmlvars; %fields(xmlvarsi);
    xmlvars00 = fieldnames(XML); 
    xmlvarsf  = cat_io_strrep(xmlvars,{'.','(',')','{','}'},{'_','','','',''}); 

    % check if fields exist 
    for fi = numel(xmlvars0):-1:1
      try
        eval(sprintf('XML(1).%s;',xmlvars0{fi})); 
      catch
        xmlvars0(fi) = [];
        xmlvars(fi)  = [];
        xmlvarsf(fi) = []; 
      end
    end

    % print existing fields
    if job.verb 
      fprintf('Found %d XML fields:  \n',numel(xmlvars)); 
      for fni=1:numel(xmlvars), cat_io_cprintf('b','    %s \n',xmlvars{fni}); end
      fprintf('\n');
    end

    %%
    for si = 1:numel(XML)
      for fi = 1:numel(xmlvars)
% =========== SPLIT multiple values ? - if yes then how many? ===========        
        try
          clear val;
          eval(sprintf('val = XML(si).%s;',xmlvars{fi})); 
        catch 
          if isnumeric(out.(xmlvarsf{fi})(1))
            val = NaN;
          elseif iscell(out.(xmlvarsf{fi})(1))
            if isnumeric(out.(xmlvarsf{fi}){1})
              var = NaN; %#ok<NASGU>
            else
              var = 'NaN'; %#ok<NASGU>
            end
          else
            val = 'NaN'; 
          end
        end
        if iscell( val ) 
          if ~iscell( val{1} )
            out.(xmlvarsf{fi}){si,1} = val; 
          else
            out.(xmlvarsf{fi}){si,1} = 'NA'; 
          end
        elseif ischar( val )
          out.(xmlvarsf{fi}){si,1} = val; 
        else
          try
            out.(xmlvarsf{fi})(si,1) = val(1); 
          catch
            if si>2 & ischar(out.(xmlvarsf{fi}){si-1,1})
              out.(xmlvarsf{fi}){si,1} = 'NA'; 
            else
              out.(xmlvarsf{fi}){si,1} = NaN; 
            end
          end
        end
      end
    end

    
    % check if fields are missing
    if ~isempty(job.csvfile) && ~strcmp(job.csvfile,'<UNDEFINED>')
      missed = setdiff( job.fields , [ job.fields(csvsvarj); xmlvars0 ] ); 
    else
      missed = setdiff( job.fields , xmlvars0 ); 
    end
    if numel(missed)>0 && ~isempty(missed{1})
      if job.dep
        fprintf('getCSVXML setup: Cannot find %d fields:  ',numel(missed)); 
      else
        fprintf('  Missed %d fields:      ',numel(missed)); 
      end
      for fni=1:numel(missed), cat_io_cprintf('r','%s \n',missed{fni}); end
      fprintf('\n'); 

      fprintf('Further %d XML fields: \n',numel(xmlvars00)); 
      for fni=1:numel(xmlvars00), cat_io_cprintf('b','  %s... \n',xmlvars00{fni}); end
      fprintf('\n\n');

    end


    if job.dep
      eval(out)
      return
    end
  end
  if numel(out)==0
    warning('cat_stat_getCSVXMLfields:dissMatch', ...
      'No matching fields at all. Stop processing.'); 
    return
  end
    



  %% coding
  % -----------------------------------------------------------------------
  if 0 %job.coding
    for fni = 1:numel(FN)
      scsv{1,fni} = FN{fni};
      if iscell(out.(FN{fni})) 
        if ischar(out.(FN{fni}){1})
          %% convert ordinary values into integer and save the coding as well as the original values seperatly
          try
            [out2.([FN{fni} '_org']), out2.([FN{fni} '_code']), out2.([FN{fni} 'c']) ] = unique(out.(FN{fni}));
          catch
            out2.([FN{fni} '_org'])  = out.(FN{fni}); 
            out2.([FN{fni} '_code']) = out.(FN{fni}); 
            out2.([FN{fni} 'c'])     = out.(FN{fni}); 
          end
        end
      end
    end
  end
  


  %% prepare and write output 
  %  ----------------------------------------------------------------------
  FN   = fieldnames(out); 
  nSN  = numel(out.(FN{1}));

  % udpate filenames
  if exist('fname','var')   
    job.fname = cat_io_strrep(job.fname,{'<FNAME>'  , '<fname>'  , '<fn>', '<FN>'}, fname);
    job.fname = cat_io_strrep(job.fname,{'<NFIELDS>', '<nfields>', '<nf>', '<NF>'}, sprintf('%d',numel(FN)-1)); 
    job.fname = cat_io_strrep(job.fname,{'<NLINES>' , '<nlines>' , '<nl>', '<NL>'}, sprintf('%d',nSN));
  end

  % create output directory and prepare filenames 
  if isempty(job.outdir), job.outdir = pwd; end
  Pcsv = fullfile(job.outdir,sprintf('%s.csv',job.fname)); 
  Ptsv = fullfile(job.outdir,sprintf('%s.tsv',job.fname)); 
  Pmat = fullfile(job.outdir,sprintf('%s.mat',job.fname)); 
  Ptxt = cell(numel(FN),1); 
  scsv = cell(numel(job.files),numel(FN)); 


  % create table and write TXT output
  for fni = 1:numel(FN)
    if cat_io_contains('txt',job.write)
      Ptxt{fni} = fullfile(job.outdir,sprintf('%s_%s.txt',job.fname,FN{fni})); 
      h = fopen(Ptxt{fni},'w');
    else
      h = 0; 
    end
    scsv{1,fni} = FN{fni};
    if iscell(out.(FN{fni})) 
      for si = 1:numel(out.(FN{fni}))
        if iscell(out.(FN{fni}){si})
          try
            for cii = 1:numel(out.(FN{fni}{si}))
              if ischar(out.(FN{fni}){si})
                scsv{si+1,li} = [scsv{si+1,li} sprintf('%s ',out.(FN{fni}){si})];
                if h, fprintf(h ,'%s ',out.(FN{fni}){si}); end
              elseif out.(FN{fni}){si} == round(out.(FN{fni}){si})
                scsv{si+1,li} = [scsv{si+1,li} sprintf('%d ',out.(FN{fni}){si})];
                if h, fprintf(h ,'%d ',out.(FN{fni}){si}); end
              else
                scsv{si+1,li} = [scsv{si+1,li} sprintf('%f ',out.(FN{fni}){si})];
                if h, fprintf(h ,'%f ',out.(FN{fni}){si}); end
              end
            end
          end
          if h, fprintf(h ,'\n'); end
        elseif ischar(out.(FN{fni}){si})
          scsv{si+1,fni} = sprintf('"%s"',out.(FN{fni}){si,:});
          if h, fprintf(h ,'"%s"\n',out.(FN{fni}){si,:}); end
        elseif out.(FN{fni}){si} == round(out.(FN{fni}){si})
          scsv{si+1,fni} = sprintf('%d',out.(FN{fni}){si});
          if h, fprintf(h ,'%d\n',out.(FN{fni}){si}); end
        else            
          scsv{si+1,fni} = sprintf('%f',out.(FN{fni}){si});
          if h, fprintf(h ,'%f\n',out.(FN{fni}){si}); end
        end
      end
    elseif all( ischar(out.(FN{fni})) )
      for si = 1:numel(out.(FN{fni}))
        scsv{si+1,fni} = sprintf('"%s"',out.(FN{fni})(si,:));
        if h, fprintf(h ,'"%s"\n',out.(FN{fni})(si,:)); end
      end
    elseif isnumeric(out.(FN{fni})) && all( out.(FN{fni}) == round(out.(FN{fni})) )
      for si = 1:numel(out.(FN{fni}))
        scsv{si+1,fni} = sprintf('%d',out.(FN{fni})(si));
        if h, fprintf(h ,'%d\n',out.(FN{fni})(si,:)); end
      end
    elseif mean( isnumeric(out.(FN{fni})) > 0.5)
      for si = 1:numel(out.(FN{fni}))
        if isnumeric(out.(FN{fni})(si))
          scsv{si+1,fni} = sprintf('%f',out.(FN{fni})(si));
          if h, fprintf(h ,'%f\n',out.(FN{fni}){si}); end
        else
          try
            scsv{si+1,fni} = sprintf('%f',double(out.(FN{fni})(si)));
            if h, fprintf(h ,'%f\n',double(out.(FN{fni}){si})); end
          catch
            scsv{si+1,fni} = nan; 
            if h, fprintf(h ,'%f\n',nan); end
          end
        end
      end
    else
      if si>1 
        if isnumeric(out.(FN{fni})(si-1))
          try
            scsv{si+1,fni} = sprintf('%f',out.(FN{fni})(si));
            if h, fprintf(h ,'%f\n',out.(FN{fni}){si}); end
          catch
            scsv{si+1,fni} = nan; 
            if h, fprintf(h ,'%f\n',nan); end
          end
        elseif ischar(out.(FN{fni})(si-1))
          try
            scsv{si+1,fni} = sprintf('"%s"',out.(FN{fni})(si,:));
            if h, fprintf(h ,'"%s"\n',out.(FN{fni}){si,:}); end
          catch
            scsv{si+1,fni} = nan; 
            if h, fprintf(h ,'"NA"\n'); end
          end
        else
          scsv{si+1,fni} = nan; 
          if h, fprintf(h ,'"NA"\n'); end
        end
      else
        scsv{si+1,fni} = nan; 
        if h, fprintf(h ,'"NA"\n'); end
      end
    end
  end
  
  % just print table
  if job.verb > 1
    fprintf('\nDipslay first 10 rows and columns:\n'); 
    disp(scsv(1:10,1:10))
  end

  % write CSV file
  if cat_io_contains('csv',job.write)
    cat_io_csv(Pcsv,scsv,struct('delimiter',job.seg(1),'komma',job.seg(2)));  
    if job.verb
      fprintf('  %s\n',spm_file( Pcsv ,'link',sprintf('open(''%s'')',Pcsv)));    
    end
  end

  % write TSV file
  if cat_io_contains('tsv',job.write)
    cat_io_csv(Ptsv,scsv);  
    if job.verb
      fprintf('  %s\n',spm_file( Ptsv ,'link',sprintf('open(''%s'')',Ptsv)));    
    end
  end

  % write MAT file
  if cat_io_contains('mat',job.write)
    cat_io_csv(Pmat,'scsv');  
    if job.verb
      fprintf('  %s\n',spm_file( Pmat ,'link',sprintf('open(''%s'')',Pmat)));    
    end
  end
  

  if job.verb, fprintf('Done\n'); end  
end
  
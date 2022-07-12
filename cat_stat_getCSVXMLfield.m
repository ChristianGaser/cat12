function out = cat_stat_getCSVXMLfield(job)
%cat_stat_getCSVXMLfield(job). Extract subject specific fields from CSV/XML
%
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$ 


%#ok<*ASGLU,*AGROW>

  SVNid = '$Rev$';
  
  def.files     = {};   % n-files, e.g. XML for direct extraction or nii/gii as selector
  def.csvfile   = {''}; % 0..1-file ... maybe n later
  def.xmlfields = {};
  def.label     = {};
  def.labels    = {};
  def.fields    = {};   % set of variables names for extraction ... preselection TIV IQR ...
                        % the variables were extracted and a depency for each created
  def.seg       = ',.'; % delimiter/komma setting {',.',';,',';.'} 
  def.csvIDfd   = 0;    % defines compontent of the pathname to match the subject
  def.outdir    = '';   % output directory
  def.fname     = '';   % write output files
  def.dep       = 0;    % internal variable to test dependencies
  def.verb      = 1; 
  
  job = cat_io_checkinopt(job,def); 
  out = struct(); 
  
  job.csvfile = char(job.csvfile); 
  job.outdir  = char(job.outdir); 
  if numel(job.csvIDfd)==1, job.csvIDfd(2) = job.csvIDfd; end
 
  % merge XML and normal fields
  %job.fields = job.fields 
  
  if job.dep 
    job.files = job.files(1); 
    job.verb  = 0; 
    if isempty(job.fields) || isempty(job.files) || strcmp(char(job.files),'<UNDEFINED>')
      return
    end
  else
    % new banner
    if job.verb, spm('FnBanner',mfilename,SVNid); end
  end 
  
  
  
  %% try to find XML files
  %if job.verb, fprintf('  Select XML files\n'); end
  % analyse dir structure
  [pp,ff,ee] = spm_fileparts(job.files{1});  
  switch ee
    case '.xml'
      Pxml = job.files; 
      for xi = 1:numel(Pxml)
        [pp,ff,ee] = spm_fileparts(Pxml{xi}); 
        Pnii{xi} = fullfile(pp,[ff(5:end) ee]);
      end
    case '.gii'
      % easy
      sinfo = cat_surf_info( job.files ); 
      Pxml  = {sinfo(:).catxml}';
      Pnii  = {sinfo(:).catxml}';
      %Pnii  = cellstr( [char(spm_str_manip(Pxml,'hh')) repmat(filesep,numel(Pxml),1) char(spm_str_manip(Pxml,'t'))] );
    case '.nii'
      % ########### this may need an external function ##########  
      % cat_vol_info
      %for fi = 1:numel(job.files)   
      %  [ppi,ffi,eei] = spm_fileparts(job.files{fi});
      %end
      
      %[pp0,pp1] = spm_fileparts(pp); 
      %if ~strcmp(pp1,'surf') || ~strcmp(pp1,'mri') || ~strcmp(pp1,'label')
      %  pp0 = pp; pp1 = ''; 
      %else
      %  pp1 = 'report'; 
      %end
      Pnii = job.files; 
      Pxml = cellstr(spm_file(Pnii,'prefix',['report' filesep 'cat_'],'ext','xml'));

      %error('cat_stat_getCSVXMLfield:badFileType','ERROR: Unsupportet file type "%s"',ee);
    otherwise
      error('cat_stat_getCSVXMLfield:badFileType','ERROR: Unsupportet file type "%s".\n',ee);
  end
  
    % check XML files
  if exist('Pxml','var') && ~isempty(Pxml)
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
  elseif exist('Pnii','var') && ~isempty(Pnii)
    for fi = 1:numel(Pnii)   
      if ~exist(Pnii{fi},'file')
        cat_io_cprintf('warn','  ERROR: Cannot find XML input %d: "%s".\n',fi,Pnii{fi})
      end
    end
  end
  
  
  % check volume label files if ROIs where selected
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
  
  
  % check surface label files if ROIs where selected
  if ~isempty( job.labels )
    Plabels = cellstr(spm_file(Pnii,'prefix',['label' filesep 'catROIs_'],'ext','xml'));
    for fi = 1:numel(Plabels)   
      if ~exist(Plabels{fi},'file') 
        Plabels{fi} = cellstr(spm_file(Pnii,'prefix','catROIs_' ,'ext','xml'));
        if ~exist(Plabels{fi},'file') 
          Plabels{fi} = ''; 
        end
      end
    end
  end
  

    
  
  % create useful fieldnames
  replaceset = {
    {' ','\t',',','.','&'}  '_' ;         
    {'(',')','{','}','[',']','='}  ''; 
    };
  if ~isempty( job.fields )
    fields = genvarname(cat_io_strrep(cat_io_strrep(cat_io_strrep( job.fields ,replaceset{1,1},replaceset{1,2}),replaceset{2,1},replaceset{2,2}),'__','_'));
  else
    fields = {};
  end
  xmlfields = {};  
  if ~isempty( job.xmlfields )
    for fxi = 1:numel(job.xmlfields)
      FN = fieldnames( job.xmlfields{fxi} ); 
      for fni = 1:numel(FN)
        xmlfields = [xmlfields; job.xmlfields{fxi}.(FN{fni}) ];
      end
    end
  end
  
  
  %% get CSV data
  % also as output list
  clear out; 
  out.files = job.files;
  
  if ~isempty(job.csvfile) && ~strcmp(job.csvfile,'<UNDEFINED>')
  % rename fielnames
   
    if job.verb, fprintf('  Read CSV files .. '); end
    if ~exist(job.csvfile,'file')
      error('cat_stat_getCSVXMLfield:noCSVfile','ERROR: No CSV file at this path "%s".\n',job.csvfile); 
    end
    
    % get table
    csv = cat_io_csv(job.csvfile,'','',struct('delimiter',job.seg(1),'komma',job.seg(2)));
    % get field names
    csvvars = genvarname(cat_io_strrep(cat_io_strrep(cat_io_strrep( csv(1,1:end) ,replaceset{1,1},replaceset{1,2}),replaceset{2,1},replaceset{2,2}),'__','_'));
    %% get original ids
    csvids  = csv(2:end,1);
    lx      = max(length(num2str(max(cell2mat(csvids))))); 
    if isnumeric( cell2mat(csvids(1)) )
      for si = 1:numel(csvids)
        csvids{si} = num2str(csvids{si},sprintf('%%0%dd',lx)); 
      end
    end
    %% csvidsl = cellfun('length',csvids);
    if ~isempty( [ cell2mat(strfind(fields,'ALLCSV')) cell2mat(strfind(fields,'ALL')) ] )
      csvsvars = csvvars; 
      csvsvari = 1:numel(csvvars); 
      csvsvarj = unique( [ cell2mat(strfind(fields,'ALLCSV')) cell2mat(strfind(fields,'ALL')) ] ); 
    else
      [csvsvars,csvsvari,csvsvarj] = intersect( csvvars , fields );
    end
    if job.verb 
      fprintf('Found %d CSV fields:  ',numel(csvsvars)); 
      for fni=1:numel(csvsvars), cat_io_cprintf('b','%s  ',csvsvars{fni}); end
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
  
  
  nsubs = max([numel(Pxml),numel(Pnii)]); 
  if ~isempty(job.fnamefields)
    %%
    for fi = 1:numel(job.fnamefields)
      for si = 1:nsubs
        if exist('Pxml','var') && ~isempty(Pxml)
          pp = Pxml{si}; ff=''; P{si} = Pxml{si}; 
        elseif exist('Pnii','var') && ~isempty(Pnii)
          pp = Pnii{si}; ff=''; P{si} = Pnii{si}; 
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
        elseif numel(job.fnamefields(fi).filesel)==1        
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
  end
  
  
  
  
  %% CSV read out
  warns = {}; 
  for fi = 1:numel(csvsvars)
    for si = 1:nsubs
      %%
      if 0 %exist('Pxml','var') && ~isempty(Pxml) % not working yet
        pp = Pxml{si}; ff=''; P{si} = Pxml{si}; 
      elseif exist('Pnii','var') && ~isempty(Pnii)
        pp = Pnii{si}; ff=''; P{si} = Pnii{si}; 
      end
      %%
      if ~isempty( job.idselector.csvIDfd ) || job.idselector.csvIDfd~=0
        for pi = 0:job.idselector.csvIDfd(1)
          [pp,ff0] = spm_fileparts(pp);
          if  job.csvIDfd(1)==pi || pi>job.idselector.csvIDfd(2)
            if diff(job.idselector.csvIDfd)
              ff = [ff0 filesep ff];   
            else
              ff = ff0; 
            end
          end
        end
        [pp,ff,ee] = spm_fileparts(ff); 

%%
        pp  = cat_io_strrep(ff,num2cell(job.idselector.fileseps),repmat({filesep},size(job.idselector.fileseps))); 
        pps = textscan(pp,'%s','Delimiter',filesep); 
        if numel(job.idselector.filesel)>1
          ff  = pps{1}{job.idselector.filesel(1):job.idselector.filesel(2)};
        elseif numel(job.idselector.filesel)==1        
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
     
      %% remove double entries by bad ids (eg. the id="1" can also be found in "101" 
      ids( idsl < max(idsl))  = [];
      % idsl( idsl < max(idsl)) = [];
      if isempty(ids)
        cat_io_cprintf('err',sprintf('Cannot find data for "%s". You selected "%s" of the filename.\n',P{si},ff));
        out.(csvsvars{fi}){si,1} = nan; 
      elseif numel(ids)>1
        if numel( unique( [csv{ids + 1,1}] ) )==1
          warn = sprintf('      Found multiple possible entries for subject "%s". Take the first entry. Check your csv file! \n',P{si}); 
          out.(csvsvars{fi}){si,1} = csv{ids + 1,csvsvari(fi)}; 
        else
          warn = sprintf('      Found multiple possible entries (index=[ %s]) for subject in "%s". Check file!\n',sprintf('%d ',ids),P{si}); 
          out.(csvsvars{fi}){si,1} = nan;
        end
        if isempty( strfind(warns,warn) )
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
  
  
  
  
  %% extract XML data
  % read XML data
  if ~isempty( job.label ) && ~isempty(Plabel)
    % read data (verbose)
    if job.verb, fprintf('  Read volume label files .. '); end
    LAB   = cat_io_xml(Plabel); 
    ROI   = job.label;
    ROIf  = cat_io_strrep(ROI,{'.','(',')','{','}'},{'_','','','',''});  
    ROIid = zeros(size(ROI)); 
    for ri = 1:numel(ROI)
      f1 = find(ROI{ri}=='(',1); 
      fe = find(ROI{ri}==')',1); 
      if ~isempty(f1) && ~isempty(fe)
       ROIid(ri) = str2double( ROI{ri}(f1+1:fe-1));
      end
    end
    if job.verb, fprintf(' done. \n'); end
    
    %% check if fields exist 
    for fi = numel(ROI):-1:1
      try
        eval(sprintf('LAB(1).%s;',ROI{fi})); 
      catch
        ROI(fi)  = [];
        ROIf(fi) = []; 
      end
    end
    
    %% set up variables
    for si = 1:numel(LAB)
      for fi = 1:numel(ROI)
        try
          eval(sprintf('val = LAB(si).%s;',ROI{fi})); 
        catch 
          if isnumeric(outr.(ROIf{fi}){1})
            val = NaN;
          else
            val = 'NaN'; 
          end
        end
        outr.(ROIf{fi})(si) = val; 
      end
    end
    
    %% write result as file
    for fi = 1:numel(ROI)
      Plabelr{fi} = fullfile(job.outdir,sprintf('%s%d_label_%s_%s.txt',job.fname,numel(job.files),ROIf{fi},...
        LAB(1).neuromorphometrics.names{ ROIid(fi) })); 
      hr = fopen(Plabelr{fi},'w'); 
      fprintf('  Write volume label: "%s"\n',spm_file( Plabelr{fi} ,'link',sprintf('open(''%s'')',Plabelr{fi}))); 
      for si = 1:numel(outr.(ROIf{fi})) %spm_file( Ptxt{fni} ,'link',sprintf('open(''%s'')',Ptxt{fni})));    
        fprintf(hr,'%f\n', outr.(ROIf{fi})(si));
      end
      fclose(hr); 
    end
  end
  
  
  
  %%
  if ~isempty( job.labels )
    if ~isempty(Plabels)
      if job.verb, fprintf('  Read surface label files .. '); end
      LABs  = cat_io_xml(Plabels); 
      ROIs  = job.label;
      ROIsf = cat_io_strrep(ROIs,{'.'},{'_'}); 
    end
  end
  
  
  
  if ~isempty(Pxml)
    if job.verb, fprintf('  Read XML files .. '); end
    XML = cat_io_xml(Pxml); 

    % remove CSV fields
    [xmlvars,xmlvarsi] = setdiff( xmlfields , csvsvars );
    xmlvars0 = xmlvars; %fields(xmlvarsi);
    xmlvarsf = cat_io_strrep(xmlvars,{'.','(',')','{','}'},{'_','','','',''}); 

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
      fprintf('Found %d XML fields:  ',numel(xmlvars)); 
      for fni=1:numel(xmlvars), cat_io_cprintf('b','%s ',xmlvars{fni}); end
      fprintf('\n');
    end

    for si = 1:numel(XML)
      for fi = 1:numel(xmlvars)
        try
          eval(sprintf('val = XML(si).%s;',xmlvars{fi})); 
        catch 
          if isnumeric(out.(xmlvarsf{fi}){1})
            val = NaN;
          else
            val = 'NaN'; 
          end
        end
        if iscell( val )
          out.(xmlvarsf{fi}){si} = val; 
        else
          out.(xmlvarsf{fi})(si) = val; 
        end
      end
    end

    
    % check if fields are missing
    if ~isempty(job.csvfile) && ~strcmp(job.csvfile,'<UNDEFINED>')
      missed = setdiff( job.fields , [ job.fields(csvsvarj); xmlvars0 ] ); 
    else
      missed = setdiff( job.fields , xmlvars0 ); 
    end
    if ~isempty(missed) 
      if job.dep
        fprintf('getCSVXML setup: Cannot find %d fields:  ',numel(missed)); 
      else
        fprintf('  Missed %d fields:      ',numel(missed)); 
      end
      for fni=1:numel(missed), cat_io_cprintf('r','%s ',missed{fni}); end
      fprintf('\n');

      fprintf('Further %d XML fields: ',numel(xmlvars)); 
      for fni=1:numel(csvvars), cat_io_cprintf('b','%s  ',csvvars{fni}); end
      fprintf('\n');

    end


    if job.dep
      return
    end
  end
  
  
    
  
  
  %% write output files
  if ~isempty(job.fname) 
    FN   = fieldnames(out); 
    
    if isempty(job.outdir), job.outdir = pwd; end
    Pcsv = fullfile(job.outdir,sprintf('%s%d.csv',job.fname,numel(job.files))); 
    Ptxt = cell(numel(FN),1); 
    scsv  = cell(numel(job.files),numel(FN)); 
    
    for fni = 1:numel(FN)
      Ptxt{fni} = fullfile(job.outdir,sprintf('%s%d_%s.txt',job.fname,numel(job.files),FN{fni})); 
      
      h = fopen(Ptxt{fni},'w'); 
      scsv{1,fni} = FN{fni};
      if iscell(out.(FN{fni})) 
        if ischar(out.(FN{fni}){1})
          %% convert ordinary values into integer and save the coding as well as the original values seperatly
          %out.([FN{fni} 'o']) = out.(FN{fni});
          [out.([FN{fni} '_org']), out.([FN{fni} '_code']), out2.(FN{fni}) ] = unique(out.(FN{fni}));
          Ptxto{fni} = fullfile(job.outdir,sprintf('%s%d_%s_original.txt',job.fname,numel(job.files),FN{fni})); 
          Ptxtc{fni} = fullfile(job.outdir,sprintf('%s%d_%s_coding.txt'  ,job.fname,numel(job.files),FN{fni})); 
          ho = fopen(Ptxto{fni},'w'); 
          hc = fopen(Ptxtc{fni},'w'); 
          fprintf(hc,'Coding of "%s":\n',FN{fni}); 
          for cii = 1:numel(out.([FN{fni} '_org']))
            fprintf(hc,'%d: %s\n', out.([FN{fni} '_code'])(cii), out.([FN{fni} '_org']){cii});
          end
          fclose(hc); 
          
        end
        %%
        for si = 1:numel(out.(FN{fni}))
          if iscell(out.(FN{fni}){si})
            for cii = 1:numel(out.(FN{fni}{si}))
              if ischar(out.(FN{fni}){si})
                scsv{si+1,li} = [scsv{si+1,li} sprintf('%d,',out2.(FN{fni}){si})];
                fprintf(h ,'%d,',out.([FN{fni} '_code']){si});
                fprintf(ho,'%s,',out.([FN{fni} '_org']){si});
                scsv{si+1,li} = [scsv{si+1,li} sprintf('%s,',out.(FN{fni}){si})];
              elseif out.(FN{fni}){si} == round(out.(FN{fni}){si})
                fprintf(h ,'%d,',out.(FN{fni}){si});
                fprintf(ho,'%d,',out.(FN{fni}){si});
                scsv{si+1,li} = [scsv{si+1,li} sprintf('%d,',out.(FN{fni}){si})];
              else
                fprintf(h ,'%f,',out.(FN{fni}){si});
                fprintf(ho,'%f,',out.(FN{fni}){si});
                scsv{si+1,li} = [scsv{si+1,li} sprintf('%f,',out.(FN{fni}){si})];
              end
            end
            fprintf(h ,'\n');
            fprintf(ho,'\n');
          elseif ischar(out.(FN{fni}){si})
            fprintf(h ,'%d\n',out.([FN{fni} '_code']));
            fprintf(ho,'%s\n',out.(FN{fni}){si});
            scsv{si+1,fni} = sprintf('%s',out.(FN{fni}){si});
          elseif out.(FN{fni}){si} == round(out.(FN{fni}){si})
            fprintf(h ,'%d\n',out.(FN{fni}){si});
%            fprintf(ho,'%d\n',out.(FN{fni}){si});
            scsv{si+1,fni} = sprintf('%d',out.(FN{fni}){si});
          else            
            fprintf(h ,'%f\n',out.(FN{fni}){si});
            fprintf(ho,'%f\n',out.(FN{fni}){si});
            scsv{si+1,fni} = sprintf('%f',out.(FN{fni}){si});
          end
        end
      elseif ischar(out.(FN{fni}))
%        fprintf(h ,'%s\n',out.(FN{fni}));
        fprintf(h ,'%d\n',out.([FN{fni} '_code']));
        fprintf(ho,'%s\n',out.(FN{fni}));
        for si = 1:numel(out.(FN{fni}))
          scsv{si+1,fni} = sprintf('%s',out.(FN{fni})(si));
        end
      elseif out.(FN{fni}) == round(out.(FN{fni}))
        fprintf(h ,'%d\n',out.(FN{fni}));
        for si = 1:numel(out.(FN{fni}))
          scsv{si+1,fni} = sprintf('%d',out.(FN{fni})(si));
        end
      else
        fprintf(h ,'%f\n',out.(FN{fni}));
        for si = 1:numel(out.(FN{fni}))
          scsv{si+1,fni} = sprintf('%f',out.(FN{fni})(si));
        end
      end
      fclose(h); 
      try 
        if ~isempty(out.(FN{fni})) && iscell( out.(FN{fni}) ) && ischar(out.(FN{fni}){1})
          fclose(ho);
        end
      catch 
        disp; 
      end
%      out = rmfield(out,FN{fni}); out.(FN{fni}) = out2.(FN{fni}); 
    end
    
    if job.verb
      fprintf('  Write %d files to ',numel(FN)); cat_io_cprintf('b',sprintf('%s',job.outdir)); fprintf(':\n'); 
      for fni = 1:numel(Ptxt)
        fprintf('    %s\n',spm_file( Ptxt{fni} ,'link',sprintf('open(''%s'')',Ptxt{fni})));    
      end
    end
    
    if job.verb > 1
      fprintf('\n'); 
      disp(scsv)
    end
    
    % write CSV filef
    cat_io_csv(Pcsv,scsv,struct('delimiter',job.seg(1),'komma',job.seg(2)));  
 
  end
  
  if job.verb, fprintf('Done\n'); end  
end
  
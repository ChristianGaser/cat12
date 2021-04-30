function out = cat_stat_getCSVXMLfield(job)
%cat_stat_getCSVXMLfield(job). Extract subject specific fields from CSV/XML
%
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$ 


%#ok<*ASGLU,*AGROW>

  SVNid = '$Rev$';
  
  def.files     = {};   % n-files, e.g. XML for direct extraction or nii/gii as selector
  def.csvfile   = {''}; % 0..1-file ... maybe n later
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
  if job.verb, fprintf('  Select XML files\n'); end
  % analyse dir structure
  [pp,ff,ee] = spm_fileparts(job.files{1});  
  switch ee
    case '.xml'
      Pxml = job.files; 
    case '.gii'
      % easy
      sinfo = cat_surf_info( job.files ); 
      Pxml  = {sinfo(:).catxml}';
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


      error('cat_stat_getCSVXMLfield:badFileType','ERROR: Unsupportet file type "%s"',ee);
    otherwise
      error('cat_stat_getCSVXMLfield:badFileType','ERROR: Unsupportet file type "%s"',ee);
  end
  % check XML files
  for fi = 1:numel(Pxml)   
    if ~exist(Pxml{fi},'file')
      cat_io_cprintf('warn','  ERROR: Cannot find input %d: "%s".',fi,Pxml{fi})
    end
  end
  
  
  % create useful fieldnames
  replaceset = {
    {' ','\t',',','.','&'}  '_' ;         
    {'(',')','{','}','[',']','='}  ''; 
    };
  fields = genvarname(cat_io_strrep(cat_io_strrep(cat_io_strrep( job.fields ,replaceset{1,1},replaceset{1,2}),replaceset{2,1},replaceset{2,2}),'__','_'));

  
  %% get CSV data
  % also as output list
  clear out; 
  out.files = job.files;
  
  if ~isempty(job.csvfile) && ~strcmp(job.csvfile,'<UNDEFINED>')
  % rename fielnames
   
    if job.verb, fprintf('  Read CSV files .. '); end
    if ~exist(job.csvfile,'file')
      error('cat_stat_getCSVXMLfield:noCSVfile','ERROR: No CSV file at this path "%s"',job.csvfile); 
    end
    
    % get table
    csv = cat_io_csv(job.csvfile,'','',struct('delimiter',job.seg(1),'komma',job.seg(2)));
    % get field names
    csvvars = genvarname(cat_io_strrep(cat_io_strrep(cat_io_strrep( csv(1,1:end) ,replaceset{1,1},replaceset{1,2}),replaceset{2,1},replaceset{2,2}),'__','_'));
    % get original ids
    csvids  = csv(2:end,1);
    if isnumeric( cell2mat(csvids(1)) )
      for si = 1:numel(csvids)
        csvids{si} = num2str(csvids{si},'%d'); 
      end
    end
    %csvidsl = cellfun('length',csvids);
    [csvsvars,csvsvari,csvsvarj] = intersect( csvvars , fields );
    if job.verb 
      fprintf('Found %d CSV fields:  ',numel(csvsvars)); 
      for fni=1:numel(csvsvars), cat_io_cprintf('b','%s  ',csvsvars{fni}); end
      fprintf('\n');
    end
    for fi = 1:numel(csvsvars)
      for si = 1:numel(Pxml)
        %%
        pp = Pxml{si}; ff='';
        for pi = 0:job.csvIDfd(1)
          [pp,ff0] = spm_fileparts(pp);
          if  job.csvIDfd(1)==pi || pi>job.csvIDfd(2)
            if diff(job.csvIDfd)
              ff = [ff0 filesep ff];   
            else
              ff = ff0; 
            end
          end
        end
        [pp,ff,ee] = spm_fileparts(ff); 
        %%
        ids = []; idsl = [];
        for sii = 1:numel(csvids)
          if ~isempty( strfind(ff,csvids{sii}) )
            ids = [ids sii]; idsl = [idsl length(csvids{sii})];
          end
        end
        % remove double entries by bad ids (eg. the id="1" can also be found in "101" 
        ids( idsl < max(idsl))  = [];
        %idsl( idsl < max(idsl)) = [];
        if isempty(ids)
          cat_io_cprintf('err',sprintf('Cannot find data for "%s".\n',Pxml{si}));
          out.(csvsvars{fi}){si,1} = nan; 
        elseif numel(ids)>1
          %%
          if numel( unique( [csv{ids + 1,1}] ) )==1
            cat_io_cprintf('warn',sprintf('      Found multiple possible entries for subject "%s". Take the first entry. Check your csv file! \n',Pxml{si}));
            out.(csvsvars{fi}){si,1} = csv{ids + 1,csvsvari(fi)}; 
          else
            cat_io_cprintf('warn',sprintf('      Found multiple possible entries for subject "%s". Check file!\n',Pxml{si}));
            out.(csvsvars{fi}){si,1} = nan; 
          end
          %csvidslmin = abs(length(ff) - csvidsl);
        else
          out.(csvsvars{fi}){si,1} = csv{ids + 1,csvsvari(fi)}; 
        end
        if fi == 1, out.ids{si,1} = ids; end
        
      end    
      if  isnumeric( out.(csvsvars{fi}){1} )
        clear temp
        temp = cell2mat(out.(csvsvars{fi})); 
        out = rmfield(out,csvsvars{fi});
        out.(csvsvars{fi}) = temp; 
      end
    end
  else
    csvsvars = {};
  end                  
  
  
  
  
  %% extract XML data
  % read XML data
  if job.verb, fprintf('  Read XML files .. '); end
  XML = cat_io_xml(Pxml); 
  
  % remove CSV fields
  [xmlvars,xmlvarsi] = setdiff( fields , csvsvars );
  xmlvars0 = job.fields(xmlvarsi);
   
  % check if fields exist 
  for fi = numel(xmlvars0):-1:1
    try
      eval(sprintf('XML(1).%s;',xmlvars0{fi})); 
    catch
      xmlvars0(fi) = [];
      xmlvars(fi)  = []; 
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
        eval(sprintf('val = XML(si).%s;',xmlvars0{fi})); 
      catch 
        val = 'NaN';
      end
      if iscell( val )
        out.(xmlvars{fi}){si} = val; 
      else
        out.(xmlvars{fi})(si) = val; 
      end
    end
  end

  
  
  % check if fields are missing
  missed = setdiff( job.fields , [ job.fields(csvsvarj)  xmlvars0 ] ); 
  if ~isempty(missed) 
    if job.dep
      fprintf('getCSVXML setup: Cannot find %d fields:  ',numel(missed)); 
    else
      fprintf('  Missed %d fields:      ',numel(missed)); 
    end
    for fni=1:numel(missed), cat_io_cprintf('r','%s ',missed{fni}); end
    fprintf('\n');
  end
  
  
  if job.dep
    return
  end
  
  
  %% write output files
  if ~isempty(job.fname)
    FN   = fieldnames(out); 
    
    if isempty(job.outdir), job.outdir = pwd; end
    if job.verb, fprintf('  Write %d files to ',numel(FN)+1); cat_io_cprintf('b',sprintf('%s',job.outdir)); end
    
    Pcsv = fullfile(job.outdir,sprintf('%s%d.csv',job.fname,numel(job.files))); 
    csv  = cell(numel(job.files),numel(FN)); 
  
    for fni = 1:numel(FN)
      P = fullfile(job.outdir,sprintf('%s%d_%s.txt',job.fname,numel(job.files),FN{fni})); 
      
      h = fopen(P,'w'); 
      csv{1,fni} = FN{fni};
      if iscell(out.(FN{fni})) 
        for si = 1:numel(out.(FN{fni}))
          if iscell(out.(FN{fni}){si})
            for cii = 1:numel(out.(FN{fni}{si}))
              if ischar(out.(FN{fni}){si})
                fprintf(h,'%s,',out.(FN{fni}){si});
                csv{si+1,li} = [csv{si+1,li} sprintf('%s,',out.(FN{fni}){si})];
              elseif out.(FN{fni}){si} == round(out.(FN{fni}){si})
                fprintf(h,'%d,',out.(FN{fni}){si});
                csv{si+1,li} = [csv{si+1,li} sprintf('%d,',out.(FN{fni}){si})];
              else
                fprintf(h,'%f,',out.(FN{fni}){si});
                csv{si+1,li} = [csv{si+1,li} sprintf('%f,',out.(FN{fni}){si})];
              end
            end
            fprintf(h,'\n');
          elseif ischar(out.(FN{fni}){si})
            fprintf(h,'%s\n',out.(FN{fni}){si});
            csv{si+1,fni} = sprintf('%s',out.(FN{fni}){si});
          elseif out.(FN{fni}){si} == round(out.(FN{fni}){si})
            fprintf(h,'%d\n',out.(FN{fni}){si});
            csv{si+1,fni} = sprintf('%f',out.(FN{fni}){si});
          else            
            fprintf(h,'%f\n',out.(FN{fni}){si});
            csv{si+1,fni} = sprintf('%f',out.(FN{fni}){si});
          end
        end
      elseif ischar(out.(FN{fni}))
        fprintf(h,'%s\n',out.(FN{fni}));
        for si = 1:numel(out.(FN{fni}))
          csv{si+1,fni} = sprintf('%s',out.(FN{fni})(si));
        end
      elseif out.(FN{fni}) == round(out.(FN{fni}))
        fprintf(h,'%d\n',out.(FN{fni}));
        for si = 1:numel(out.(FN{fni}))
          csv{si+1,fni} = sprintf('%d',out.(FN{fni})(si));
        end
      else
        fprintf(h,'%f\n',out.(FN{fni}));
        for si = 1:numel(out.(FN{fni}))
          csv{si+1,fni} = sprintf('%f',out.(FN{fni})(si));
        end
      end
      fclose(h); 
    end
    
    % write CSV filef
    cat_io_csv(Pcsv,csv,struct('delimiter',job.seg(1),'komma',job.seg(2)));  
  end
  
  if job.verb, fprintf('\nDone\n'); end  
end
  
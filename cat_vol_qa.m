function varargout = cat_vol_qa(action,varargin)
% CAT Preprocessing T1 Quality Control
% ______________________________________________________________________
% 
% Estimation of image quality measures like noise, inhomogeneity,
% contrast, resolution, etc. and scaling for school marks. 
%
% [QAS,QAM] = cat_vol_qa(action,varargin)
% 
%  action ..
%  1) 'caterr'    .. short version without image analysis used in the main
%                    CAT preprocessing. 
%
%  2) 'cat12ver'  .. process with an older cat_vol_qa version
%       To use older versions, they have to be renamed (including all
%       internal calls) and added to the update_rating subfunction here. 
%
%  3) 'cat12'     .. CAT internal preprocessing interface 
%     (this is the processing case that is also called in all other cases)
%
%       [QAS,QAM] = cat_vol_qa('cat12',Yp0,Po,Ym,res[,opt])
% 
%       Pp0 .. segmentation files (p0*.nii)
%       Po  .. original files (*.nii)
%       Pm  .. modified files (m*.nii)
%       Yp0 .. segmentation image matrix
%       Ym  .. modified image matrix
%
%       opt            = parameter structure
%       opt.verb       = verbose level  [ 0=nothing | 1=points | 2*=times ]
%       opt.redres     = resolution in mm for intensity scaling [ 4* ];
%       opt.write_csv  = final cms-file
%       opt.write_xml  = images base xml-file
%       opt.sortQATm   = sort QATm output
%       opt.orgval     = original QAM results (no marks)
%       opt.recalc     =
%       opt.avgfactor  = 
%       opt.prefix     = prefix of xml output file (default cat_*.xml) 
%
%  4) 'p0' .. direct call
%     
%       qa = cat_vol_qa('p0',Pp0,opt);
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
%
% $Id$
% ______________________________________________________________________

%#ok<*ASGLU>

  % init output
  QAS = struct(); 
  QAR = struct(); 
  
  if nargin == 0, help cat_vol_qa; return; end

  if ischar(action) && strcmp(action,'getdef')
    varargout{1} = upate_rating(struct(),varargin{1},1);
    return
  end

  if isstruct(action)
    if isfield(action,'reportfolder') && isempty(action.reportfolder)
      mrifolder    = '';
      reportfolder = ''; 
    else
      mrifolder    = 'mri';
      reportfolder = 'report'; 
    end
  else
    try
      switch action
        case 'cat12err'
          [mrifolder, reportfolder] = cat_io_subfolders(varargin{1}.job.data{1},varargin{1}.job);
        case 'cat12'
          [mrifolder, reportfolder] = cat_io_subfolders(varargin{2},varargin{6}.job);
        case 'cat12ver'
          if isfield(varargin{6},'reportfolder') && isempty(varargin{6}.reportfolder)
            mrifolder    = '';
            reportfolder = ''; 
          else
            mrifolder    = 'mri';
            reportfolder = 'report'; 
          end
        otherwise
          [mrifolder, reportfolder] = cat_io_subfolders(varargin{4}.catlog,varargin{6}.job);
      end
    catch
      mrifolder    = 'mri';
      reportfolder = 'report'; 
    end
  end
  

  % handling of batch calls, where action is the batch job-structure from
  % SPM, and other actions used later
  action2 = action; 
  if isstruct(action)
    % SPM batch structure 
    if isfield(action,'model')
      if isfield(action.model,'catp0')
        Po  = action.images;
        Pp0 = action.model.catp0; 
        if numel(Po)~=numel(Pp0) && isscalar(Pp0)
          Pp0 = repmat(Pp0,numel(Po),1);
        end
        Pm  = action.images;
        action.data = Pp0;
      elseif isfield(action.model,'spmc0')
        Po  = action.images;
        Pp0 = action.model.spmc0; 
        if numel(Po)~=numel(Pp0) && isscalar(Pp0)
          Pp0 = repmat(Pp0,numel(Po),1);
        end
        Pm  = action.images;
        action.data = Pp0;  
      elseif isfield(action.model,'spmc1')
        %% SPM case where the spmc1 should point to all other segmentation files
       
        % prepare and check other input files
        spmp0ne = false(numel(action.model.spmc1),1); 
        for si = 1:numel(action.model.spmc1)
          [pp,ff,ee]               = spm_fileparts( action.model.spmc1{si} ); 
          action.model.spmp0{si,1} = fullfile(pp,['c0' ff(3:end) ee]);
          action.model.spmc2{si,1} = fullfile(pp,['c2' ff(3:end) ee]);
          action.model.spmc3{si,1} = fullfile(pp,['c3' ff(3:end) ee]);
          spmp0ne(si,1)            = ~exist(action.model.spmp0{si},'file');
        end        

        % prepare batch that produce label maps
        if sum(spmp0ne)>0 
          fprintf('Prepare CGW-label maps:\n')
          mjob.images{1}   = action.images(spmp0ne); 
          mjob.images{2}   = action.model.spmc3(spmp0ne); % CSF 
          mjob.images{3}   = action.model.spmc1(spmp0ne); % GM
          mjob.images{4}   = action.model.spmc2(spmp0ne); % WM
          mjob.expression  = 'i1*0 + i2*1 + i3*2 + i4*3';  % the first one is just for the name
          mjob.prefix      = 'c0'; 
          mjob.verb        = 0; 
          mjob.cleanup     = 0;
          mjob.ignoreBIDS  = 1; 
          cat_vol_mimcalc(mjob);
        end

        %% run QC
        action2 = rmfield(action,'model'); 
        action2.model.spmc0  = action.model.spmp0;
        action2.reportfolder = ''; 
       
        out = cat_vol_qa(action2,varargin); 

        varargout{1}.data = action2.images;
        varargout{2} = out; 
        for pi = 1:numel(action2.images)
          [pp,ff,ee] = spm_fileparts(action2.images{pi});
          varargout{1}.xmls{pi} = fullfile(pp,reportfolder,[action2.opts.prefix ff '.xml']);
        end
        return 
      
      elseif isfield(action.model,'seg')
        %% error('no implemented yet')
        fprintf('Prepare CGW-label maps:\n')
        mjob.images{1,1} = action.images; 
        mjob.images{2,1} = action.model.seg.cm;
        mjob.images{3,1} = action.model.seg.gm;
        mjob.images{4,1} = action.model.seg.wm;
        mjob.expression  = 'i1*0 + i2*1 + i3*2 + i4*3';  % the first one is just for the name
        mjob.prefix      = 'p0'; %
        mjob.verb        = 0; 
        cat_vol_mimcalc(mjob);

        action2 = rmfield(action,'model'); 
        action2.model.catp0 = spm_file(action.images,'prefix','p0');  

        varargout{:} = cat_vol_qa(action2,varargin); 
        return 

      end
    end
    if isfield(action,'data')
      Pp0 = action.data;
    end
    action = 'p0';
  end
  if nargin==3 && isstruct(varargin{2}) && isstruct(varargin{2})
    opt  = cat_check('checkinopt',varargin{2},defaults);
    nopt = 1; 
  elseif nargin>6 && isstruct(varargin{6}) && isstruct(varargin{6})
    opt  = cat_check('checkinopt',varargin{6},defaults);
    nopt = 1;   
  else
    if isstruct(action2)
      opt = cat_check('checkinopt',action2.opts,defaults);
    else
      opt = defaults;
    end
    nopt = 0;
    if isfield(opt,'recalc') && opt.recalc, opt.rerun = opt.recalc; end
  end
  if cat_io_contains(opt.prefix,'VERSION')
    opt.prefix = strrep( opt.prefix , 'VERSION', strrep( opt.version ,'_','')); 
  end
  if isfield(opt,'model') && isfield(opt.model,'spmc1')
    opt.reportfolder = ''; 
  else 
    opt.reportfolder = reportfolder;
  end
  

  % check input by action
  switch action
    case 'p0'
    % segment image cases
      if nargin<=3 && ( ~exist('Pp0','var') || isempty(Pp0) )
        if (nargin-nopt)<2  
          Pp0 = cellstr(spm_select(inf,'image',...
            'select p0-segment image',{},pwd,'^p0.*')); 
          if isempty(Pp0{1}), return; end
        else
          Pp0 = varargin{1};
        end
        if numel(action)==2
          Po = Pp0; Pm = Pp0;
          for fi=1:numel(Pp0)
            [pp,ff,ee] = spm_fileparts(Pp0{fi}); 
            if strcmp(ee,'.gz'), ee =[ff(end-3:end) ee]; ff = ff(1:end-4);  end
            [ppa,ppb] = spm_fileparts(pp); 
            if strcmp(ppb,'mri'), ppo = ppa; else, ppo = pp; end 

            deri = min([ strfind(ppo,[filesep 'derivatives' filesep 'CAT']), ...
                         strfind(ppo,[filesep 'derivatives' filesep 'SPM']), ...
                         strfind(ppo,[filesep 'derivatives' filesep 'T1P'])]); 
            if isempty( deri )
              if cat_io_contains(ff,'qcseg') 
                ff = strrep(ff,'p0_qcseg_',''); 
                Po{fi} = fullfile(ppo,[ff ee]);
              elseif strcmp(ff(1:2),'p0') || strcmp(ff(1:2),'c1')
                Po{fi} = fullfile(ppo,[ff(3:end) ee]); 
              elseif cat_io_contains(ff,'synthseg_p0')
                Po{fi} = fullfile(ppo,[strrep(ff,'synthseg_p0','') ee]); 
              end
            else
              BIDShome = fileparts(ppo(1:deri(1)));
              fsep     = strfind( ppo(deri(1) + 16:end) , filesep ) + deri(1) + 16;
              if isempty(fsep)
                Pofi = cat_vol_findfiles(BIDShome,[ff(3:end) '.nii.gz']); 
                if isempty(Pofi)
                  Pofi = cat_vol_findfiles(BIDShome,[ff(3:end) '.nii']); 
                end
                Po{fi} = Pofi{1}; 
              else
                Po{fi} = fullfile( BIDShome, ppo(fsep(1):end) , [ff(3:end) ee] );
              end
            end
            Pm{fi} = fullfile(pp,['m'  ff(3:end) ee]);
            if ~exist(Pm{fi},'file'), Pm{fi}=[Pm{fi} '.gz']; end
            if ~exist(Pm{fi},'file'), Pm{fi}=''; end
          end
        else
          Po = cellstr(spm_select(repmat(numel(Pp0),1,2),...
            'image','select original image(s)',{},pwd,'.*')); 
          Pm = cellstr(spm_select(repmat(numel(Pp0),1,2),...
            'image','select modified image(s)',{},pwd,'.*')); 
        end

      elseif nargin<=5 && ( ~exist('Pp0','var') || isempty(Pp0) )
        Pp0 = varargin{1};
        Po  = varargin{2};
        Pm  = varargin{3};

      elseif ( ~exist('Pp0','var') || isempty(Pp0) )
        error('MATLAB:cat_vol_qa:inputerror',...
          'Wrong number/structure of input elements!'); 
      end

    case {'cat12ver', 'cat12'}
      % nothing to do

    case 'cat12err'
      % again ?
      %opt  = cat_check('checkinopt',varargin{6},defaults);  

    otherwise

      error('MATLAB:cat_vol_qa:inputerror',...
        'Wrong number/structure of input elements!'); 
  end
  if ~exist('species','var'), species='human'; end
    
  
  %
  % --------------------------------------------------------------------
  [QA,QMAfn]  = cat_stat_marks('init'); 
  % remove res_ECR in case of given older versions 
  if any( strcmp(opt.version, opt.versions0) )
    QMAfn( cat_io_contains(QMAfn,'res_ECR'))   = [];
  end
  % remove FEC in case of given older versions 
  if any( strcmp(opt.version, opt.versions1) )
    QMAfn( cat_io_contains(QMAfn,'FEC'))   = [];
  end
  stime       = clock;
  


  
  % Print options
  % --------------------------------------------------------------------
  Cheader = {'scan'};
  Theader = sprintf(sprintf('%%%ds:',opt.snspace(1)-1),'scan');
  Tline   = sprintf('%%5d) %%%ds:',opt.snspace(1)-8);
  Tline2  = sprintf('%%5d) %%6s%%%ds:',opt.snspace(1)-14); 
  Tavg    = sprintf('%%%ds:',opt.snspace(1)-1);
  TlineE  = sprintf('%%5d) %%%ds: %%s',opt.snspace(1)-7);
  for fi=1:numel(QMAfn)
    Cheader = [Cheader QMAfn{fi}]; %#ok<AGROW>
    QAMfni  = strrep( QMAfn{fi} ,'_','');
    Theader = sprintf(sprintf('%%s%%%ds',opt.snspace(2)),Theader,...
                cat_io_strrep( QAMfni,{'contrastr';'resECR'},{'CON';'ECR'}) ); %(1:min(opt.snspace(2)-1,numel(QMAfn{fi}))) 
    Tline   = sprintf('%s%%%d.%df',Tline,opt.snspace(2),opt.snspace(3));
    Tline2  = sprintf('%s%%%d.%df',Tline2,opt.snspace(2),opt.snspace(3));
    Tavg    = sprintf('%s%%%d.%df',Tavg,opt.snspace(2),opt.snspace(3));
  end
  Cheader = [Cheader 'IQR']; 
  Theader = sprintf(sprintf('%%s%%%ds',opt.snspace(2)),Theader,'IQR');
  if ~any( strcmp(opt.version,opt.versions0) ) % ~any( cat_io_contains(opt.version,opt.versions0) )
    Tline   = sprintf('%s%%%d.%df',Tline,opt.snspace(2),opt.snspace(3));
    Tline2  = sprintf('%s%%%d.%df',Tline2,opt.snspace(2),opt.snspace(3));
    Cheader = [Cheader 'SIQR']; 
    Theader = sprintf(sprintf('%%s%%%ds',opt.snspace(2)),Theader,'SIQR');
    Tavg    = sprintf('%s%%%d.%df',Tavg,opt.snspace(2),opt.snspace(3));
  end
  Tline   = sprintf('%s%%%d.%df%%s\n',Tline,opt.snspace(2),opt.snspace(3));
  Tline2  = sprintf('%s%%%d.%df\n',Tline2,opt.snspace(2),opt.snspace(3));
  Tavg    = sprintf('%s%%%d.%df\n',Tavg,opt.snspace(2),opt.snspace(3));
  
  
  
  % estimation part    
  switch action
    case 'cat12'
    % Direct call of the specific QC version with input images given by the 
    % varargin structure used in the CAT pipeline (processing of one case)
%sprintf('[QAS,QAR] = %s(''cat12'',varargin{:});', opt.version)
if isstruct(varargin{end-1}), varargin{end-1}.write_xml = 0; end

      eval(sprintf('[QAS,QAR] = %s(''cat12'',varargin{:});', opt.version));
      QAR  = upate_rating(QAS,opt.version);

    case 'p0'    
    % Default case of multiple input files where we have to load the images
    % and will also show the resulting ratings (batch processing of many 
    % files.

      % return for empty input
      if isempty(Pp0) || (isempty(Pp0{1}) && numel(Pp0)<=1) 
        cat_io_cprintf('com','No images for QA!\n'); 
        return
      end
 
      % remove num entry from SPM GUI input
      Po  = cat_io_strrep(Po, {'.nii.gz,1','.nii,1'},{'.nii.gz','.nii'}); 
      Pm  = cat_io_strrep(Pm, {'.nii.gz,1','.nii,1'},{'.nii.gz','.nii'}); 
      Pp0 = cat_io_strrep(Pp0,{'.nii.gz,1','.nii,1'},{'.nii.gz','.nii'}); 

      % if files are zipped
      Poe = cellfun(@(x) exist(x,'file'), Po); 
      Po(Poe==0) = spm_file(Po(Poe==0),'ext','.nii.gz'); clear Poe 


      % name segmentation if possible
      if isempty(Pp0{1})
        Pp0 = Po; 
        [pp,ff,ee] = spm_fileparts(Po{1});
        segment = 'internal'; 
      else
        [pp,ff,ee] = spm_fileparts(Pp0{1});
        switch ff(1:2)
          case 'sy',  segment = 'synthseg'; 
          case 'c0',  segment = 'SPM'; 
          case 'c1',  segment = 'SPM'; 
          case 'p0',  segment = 'CAT'; 
          otherwise,  segment = 'internal'; 
        end
      end
      % print title
      if opt.verb>1
        fprintf('\nCAT Preprocessing T1 Quality Control (');
        [ver_cat, rev_cat] = cat_version;
        cat_io_cprintf('blue',' %s',opt.version );
        cat_io_cprintf('n',' ,'); 
        cat_io_cprintf('blue',' %s',segment );
        cat_io_cprintf('n',' , %s ):\n',sprintf('Rev: %s %s', ver_cat, rev_cat) );
        fprintf('\n%s\n%s\n',  Theader,repmat('-',size(Theader)));  
      end

      % preare variables
      qamat   = nan(numel(Po),numel(QMAfn));
      qamatm  = nan(numel(Po),numel(QMAfn));
      mqamatm = 10.5*ones(numel(Po),1 + ~any( strcmp(opt.version,opt.versions0) ));
      QAS     = struct(); 
      QAR     = struct(); 
      
      QAR.mark2rps = @(mark) min(100,max(0,105 - mark*10)) + isnan(mark).*mark;



      % loop for multiple files
      % -------------------------------------------------------------------
      for fi = 1:numel(Pp0)
        stime1 = clock; 
        
        % setup the XML file name
        [pp,ff,ee] = spm_fileparts( strrep(Pp0{fi},'.nii.gz','.nii') ); ff(1:2) = []; 
        [ppa,ppb] = spm_fileparts(pp); 
        if strcmp(ppb,'mri'), ppo = fullfile(ppa,reportfolder); else, ppo = pp; end 
        sfile   = fullfile(ppo,[opt.prefix ff '.xml']); 

        if ~exist( sfile ,'file') 
        % test if file may exist in raw directory 
          ppo    = spm_fileparts( strrep(Po{fi},'.nii.gz','.nii') );
          sfileo = spm_file(Po{fi},'path',ppo,'prefix',opt.prefix,'ext','.xml'); 
          if exist(sfileo,'file')
            sfile = sfileo; 
          else
            sfileo = spm_file(Po{fi},'path',fullfile(ppo,'report'),'prefix',opt.prefix,'ext','.xml'); 
            if exist(sfileo,'file')
              sfile = sfileo; 
            end
          end
        end
        
        if ~exist( sfile ,'file') 
        % not processed (eg. no additional comment)
          run   = 1; 
          rerun = '';
        elseif opt.rerun == 2  ||  strcmp(opt.prefix,'tmp')
        % forced processing
          run   = 1;
          if cat_get_defaults('extopts.expertgui') > 1
            rerun = ' forced';
          else
            rerun = ''; 
          end
        elseif (opt.rerun == 0 && cat_io_rerun(Po{fi}, sfile , 0 , 0)) || ... load if the QC file is available and newer than the input
               (opt.rerun == 1 && ( cat_io_rerun(Po{fi}, sfile , 0 , 0 ) || ...  load if the QC file is newer than the input and function
               (0 && cat_io_rerun(which(opt.version), sfile , 0 , 0)) ) ) 
          run   = 1;
          rerun = ' updated';
        else
          run   = 0; 
          rerun = ' loaded';
        end

        if ~run
          % If no reprocessing is required then just try to load the values 
          % from the QC XML but update the rating
          try 
            QASfi  = cat_io_xml(sfile); 
            QARfi  = upate_rating(QASfi,opt.version);
            cat_io_xml( sfile, QARfi, 'write+');

            %QASfix = cat_io_updateStruct(QASfi, QARfi);
            
            % try to update the QC structure
            [QAS, QAR, qamat, qamatm, mqamatm] = updateQAstructure(QAS, ...
              QASfi, QAR, QARfi, qamat, qamatm, mqamatm, QMAfn, fi);
            
            rerun   = ' loaded';
          catch e
          % in case of problems we have to reprocess
            run     = 1; 
            rerun   = ' reprocessed';
          end
        end

        clear e
        if run
        % try to run QC estimation 

          % load images
          try
            [Yp0,Ym,Vo] = getImages(Pp0,Po,Pm,fi); 
          catch e
            % setup default output
            opt2      = opt; 
            opt2.subj = fi;
            opt2.job  = cat_get_defaults;
            opt2.job.channel.vols{fi}     = [Po{fi} ',1'];
            opt2.job.data{fi}             = [Po{fi} ',1'];
            opt2.job.extopts.darteltpms   = {};
            opt2.job.extopts.shootingtpms = {};
            opt2.job.extopts.subfolder    = isfield(action,'model') && isfield(action.model,'spmc1');  
            opt2.caterr     = struct();
            opt2.caterrtxt  = ''; 
            
            [QASfi,QARfi] = cat12err(opt2,mrifolder,reportfolder);

            % try to update the QC structure
            [QAS, QAR, qamat, qamatm, mqamatm] = updateQAstructure(QAS, ...
              QASfi, QAR, QARfi, qamat, qamatm, mqamatm, QMAfn, fi);

            continue
          end

          % general function called from CAT
          if ~exist( Po{fi} ,'file')
            continue
          end

          evalc('res.image = spm_vol(Po{fi});'); 

          if ~isempty(Yp0)
            try
              [QASfi,QARfi] = cat_vol_qa('cat12ver',Yp0,Vo,Ym,res,species,opt,Pp0{fi});
            catch
              cat_io_cprintf('warn',sprintf('Failed ...run cat_vol_qa202412')); 
              opt2 = opt; opt2.version = 'cat_vol_qa202412';
              [QASfi,QARfi] = cat_vol_qa('cat12ver',Yp0,Vo,Ym,res,species,opt2,Pp0{fi});
            end
          else
            opt2 = opt; opt2.version = 'cat_vol_qa202412';
            [QASfi,QARfi] = cat_vol_qa('cat12ver',Yp0,Vo,Ym,res,species,opt2,Pp0{fi});
          end
          try
            % try to update the QC structure
            [QAS, QAR, qamat, qamatm, mqamatm] = updateQAstructure(QAS, ...
              QASfi, QAR, QARfi, qamat, qamatm, mqamatm, QMAfn, fi);
            
          catch
            %% this is not a good solution !
            opt2      = opt; 
            opt2.subj = fi;
            opt2.job  = cat_get_defaults;
            opt2.job.channel.vols{fi}     = [Po{fi} ',1'];
            opt2.job.data{fi}             = [Po{fi} ',1'];
            opt2.job.extopts.darteltpms   = {};
            opt2.job.extopts.shootingtpms = {};
            opt2.caterr     = struct();
            opt2.caterrtxt  = ''; 
            
            [QASfi,QARfi] = cat12err(opt2,mrifolder,reportfolder);

            %% try to update the QC structure
            [QAS, QAR, qamat, qamatm, mqamatm] = updateQAstructure(QAS, ...
              QASfi, QAR, QARfi, qamat, qamatm, mqamatm, QMAfn, fi);
          
          end
        end

        % print result
        if opt.verb>1 
          % update the rerun parameter
          rerun = sprintf('%s%3.0fs',rerun, etime(clock,stime1)); 

          if exist('e','var')
          % write short error-message in case of error
            switch e.identifier
              case {'cat_vol_qa:noYo','cat_vol_qa:noYm','cat_vol_qa:badSegmentation'}
                em = e.identifier;
              case 'cat_vol_qa:missingVersion'
                rethrow(e);
              case 'MATLAB:badsubscript'
              otherwise
                em = ['ERROR:\n' repmat(' ',1,10) e.message '\n'];
                for ei=1:numel(e.stack)
                  em = sprintf('%s%s%5d: %s\n',em,repmat(' ',1,10),...
                    e.stack(ei).line(end),e.stack(ei).name);
                end  
            end

            [pp,ff] = spm_fileparts(Po{fi});
            QAS(fi).filedata.fnames = [spm_str_manip(pp,sprintf('k%d',floor( (opt.snspace(1)-19) /3) - 1)),'/',...
                                 spm_str_manip(ff,sprintf('k%d',(opt.snspace(1)-19) - floor((opt.snspace(1)-14)/3)))];
            cat_io_cprintf(opt.MarkColor(end,:),sprintf(TlineE,fi,Pp0{fi},[em '\n']));
          else
          % write measurements and time
            if opt.orgval 
              cat_io_cprintf(opt.MarkColor(min(size(opt.MarkColor,1),max(1,floor( mqamatm(fi,end)/9.5 * ...
                size(opt.MarkColor,1)))),:),sprintf(Tline,fi,...
                spm_str_manip(QAS(fi).filedata.fname,['a' num2str(opt.snspace(1) - 14)]),...
                qamat(fi,:), max(1,min(9.5,mqamatm(fi,:))), rerun));
            else
              %%
              cat_io_cprintf(opt.MarkColor(min(size(opt.MarkColor,1),max(1,floor( mqamatm(fi,end)/9.5 * ...
                size(opt.MarkColor,1)))),:),sprintf(Tline,fi,...
                spm_str_manip(QAS(fi).filedata.fname,['a' num2str(opt.snspace(1) - 14)]),...
                qamatm(fi,:), max(1,min(9.5,mqamatm(fi,:))), rerun));
            end
          end
        end
      end
      

      % sort by mean mark
      % -------------------------------------------------------------------
      if opt.sortQATm && numel(Po)>1
        % sort matrix
        [smqamatm,smqamatmi] = sort(mqamatm(:,end),'ascend');
        sqamatm  = qamatm(smqamatmi,:);
        sqamat   = qamat(smqamatmi,:); 

        % print matrix
        if opt.verb>0
          fprintf('%s\n',repmat('-',size(Theader))); 
          for fi = 1:numel(QAS)
            if isfield( QAS(smqamatmi(fi)), 'filedata') && isfield( QAS(smqamatmi(fi)).filedata, 'fname')
              fname = spm_str_manip( QAS(smqamatmi(fi)).filedata.fname , ['a' num2str(opt.snspace(1) - opt.snspace(2) - 14)] );
            else
              fname = 'FILENAME ERROR';
            end
            if opt.orgval 
              cat_io_cprintf(opt.MarkColor(max(1,min(size(opt.MarkColor,1),...
                round( mqamatm(smqamatmi(fi),end)/9.5 * ...
                size(opt.MarkColor,1)))),:),sprintf(...
                Tline2,fi,sprintf('(%d)',smqamatmi(fi)),...
                fname, sqamat(fi,:),max(1,min(10.5,mqamatm(smqamatmi(fi),:)))));
            else
              cat_io_cprintf(opt.MarkColor(max(1,min(size(opt.MarkColor,1),...
                round( mqamatm(smqamatmi(fi),end)/9.5 * ...
                size(opt.MarkColor,1)))),:),sprintf(...
                Tline2,fi,sprintf('(%d)',smqamatmi(fi)),...
                fname, sqamatm(fi,:),mqamatm(smqamatmi(fi),:)));
            end
          end
        end
      end
      % print the results for each scan 
      if opt.verb>1 && numel(Pp0)>1
        fprintf('%s\n',repmat('-',size(Theader)));  
        if opt.orgval 
          fprintf(Tavg,'mean', cat_stat_nanmean(qamat,1), cat_stat_nanmean(mqamatm,1));   %#ok<CTPCT>
          fprintf(Tavg,'std' , cat_stat_nanstd(qamat,1),  cat_stat_nanstd(mqamatm,1));    %#ok<CTPCT>  
        else
          fprintf(Tavg,'mean', cat_stat_nanmean(qamatm,1), cat_stat_nanmean(mqamatm,1));   %#ok<CTPCT>
          fprintf(Tavg,'std' , cat_stat_nanstd(qamatm,1),  cat_stat_nanstd(mqamatm,1));    %#ok<CTPCT>  
        end 
      end
      if opt.verb>0, fprintf('\n'); end


      
      % result tables (cell structures)
      % ----------------------------------------------------------------
      if nargout>2 && opt.write_csv
        QAT   = [Cheader(1:end-1); ... there is no mean for the original measures
                 Po               , num2cell(qamat); ...
                 'mean'           , num2cell(cat_stat_nanmean(qamat,1)); ...
                 'std'            , num2cell( cat_stat_nanstd(qamat,1,1))];
        QATm  = [Cheader; ...
                 Po               , num2cell(qamatm)          , ...
                                    num2cell(cat_stat_nanmean(qamatm,2)); ...
                 'mean'           , num2cell(cat_stat_nanmean(qamatm,1))  , ...
                                    num2cell(cat_stat_nanmean(mqamatm,1)); ...
                 'std'            , num2cell( cat_stat_nanstd(qamatm,1,1)), ...
                                    num2cell( cat_stat_nanstd(mqamatm,1))];


        % write csv results
        % --------------------------------------------------------------
        if opt.write_csv
          pp = spm_fileparts(Pp0{1});
          cat_io_csv(fullfile(pp,reportfolder,[opt.prefix num2str(numel(Vo),'%04d') ...
            'cat_vol_qa_values.csv']),QAT);
          cat_io_csv(fullfile(pp,reportfolder,[opt.prefix num2str(numel(Vo),'%04d') ...
            'cat_vol_qa_marks.csv']),QATm);
        end
      end 
      
      if opt.verb>0
        fprintf('Quality Control for %d subject was done in %0.0fs\n', ...
          numel(Pp0),etime(clock,stime)); fprintf('\n');
      end



    case 'cat12err'
      opt = cat_check('checkinopt',varargin{1},defaults);
      QAS = cat12err(opt,mrifolder,reportfolder);


    case 'cat12ver'
      % main processing with subversions 
      [pp,ff,ee] = spm_fileparts(strrep(varargin{2}.fname,'.nii.gz','.nii'));   
      
      % Call of different versions of the QC:
      % -------------------------------------------------------------------
      % estimation of the measures for the single case by different versions. 
      % To use other older cat_vol_qa versions copy them into a path and 
      % rename the filename and the all internal use of the functionname.
      % Extend the default variable versions0 for older functions without
      % SIQR and res_ECR measure. 
      % Older versions may use different parameters - check similar
      % pepared verions. 
      % -------------------------------------------------------------------
      if isfield(opt,'version') 
        if ~exist(opt.version,'file')
          error('cat_vol_qa:missingVersion','Selected QC version "%s" is not available! ',opt.version); 
        elseif ~strcmp(opt.version,mfilename)
          switch opt.version
            % in older versions some parameters where defined different 
            % and we need to update them
            case {'cat_vol_qa201602'}
              % here the 
              vx_vol  = sqrt(sum(varargin{2}.mat(1:3,1:3).^2));
              Yp0toC  = @(Yp0,c) 1-min(1,abs(Yp0-c));
              qa.subjectmeasures.vol_TIV = sum(varargin{1}(:)>0) ./ prod(vx_vol) / 1000;
              for i = 1:3
                qa.subjectmeasures.vol_abs_CGW(i) = sum( Yp0toC(varargin{1}(:),i)) ./ prod(vx_vol) / 1000; 
                qa.subjectmeasures.vol_rel_CGW(i) = qa.subjectmeasures.vol_abs_CGW(i) ./ ...
                                                     qa.subjectmeasures.vol_TIV; 
              end
              varargin2 = varargin; 
              varargin2{6}.job.extopts.subfolders = ~isempty(reportfolder); 
              eval(sprintf('[QAS,QAR] = %s(''cat12'',varargin2{1:4},struct(),varargin2{5:end});',opt.version));
            otherwise
              varargin2 = varargin; 
              varargin2{6}.job.extopts.subfolders = ~isempty(reportfolder); 
              eval(sprintf('[QAS,QAR] = %s(''cat12'',varargin2{:});',opt.version));
          end
        else
        % setup the current/default version 
          varargin2 = varargin; 
          varargin2{6}.version = 'cat_vol_qa202310'; 
          varargin2{6}.job.extopts.subfolders = ~isempty(reportfolder); 
          eval(sprintf('[QAS,QAR] = %s(''cat12'',varargin2{:});', varargin2{6}.version));
        end
      end

      % redo the quality rating, ie. measure scaling
      QAR = upate_rating(QAS,opt.version);
      QAS.subjectratings = QAR.subjectratings; 
      QAS.qualityratings = QAR.qualityratings; 
      

      if nargin>6
        pp0     = spm_fileparts(varargin{7}); 
        [ppx,ffx] = spm_fileparts(pp0); 
        if strcmp(ffx,'mri')
          sfile   = fullfile(ppx,reportfolder,[opt.prefix ff '.xml']); 
          catfile = fullfile(ppx,reportfolder,['cat_' ff '.xml']); 
        else
          sfile   = fullfile(pp0,[opt.prefix ff '.xml']); 
          catfile = fullfile(pp0,['cat_' ff '.xml']); 
        end
      else
        sfile   = fullfile(pp,reportfolder,[opt.prefix ff '.xml']); 
        catfile = fullfile(pp,reportfolder,['cat_'  ff '.xml']); 
      end  
      if exist(sfile,'file')
        S  = cat_io_xml( sfile ); 
        SN = cat_io_mergeStruct( S , QAS , [], 1); 
      elseif exist(catfile,'file')
        S  = cat_io_xml( catfile ); 
        SN = cat_io_mergeStruct( S , QAS , [], 1); 
      else
        SN = QAS;
      end
      
      cat_io_xml( sfile, SN );

      
    otherwise
      % catched before
  end

  % export 
  if strcmp(action,'cat12') % exist('Pp0','var') && isscalar(Pp0) && opt.write_xml
    QAS.qualityratings = QAR.qualityratings;
    QAS.subjectratings = QAR.subjectratings;
    QAS.ratings_help   = QAR.help;

    [pp,ff] = spm_fileparts(QAS.filedata.fname);  
    cat_io_xml( fullfile(pp,reportfolder,[opt.prefix ff '.xml']) ,QAS,'write+'); 
    cat_io_xml( fullfile(pp,reportfolder,[opt.prefix ff '.xml']) ,QAR,'write+'); %struct('QAS',QAS,'QAM',QAM)
  end

  if (isempty(varargin) || isstruct(varargin{1}) || isstruct(action)) && exist('Pp0','var')
  % SPM batch output case
    varargout{1}.data = Pp0;
    for pi = 1:numel(Pp0)
      [pp,ff,ee] = spm_fileparts(Pp0{pi});
      varargout{1}.xmls{pi} = fullfile(pp,reportfolder,[opt.prefix ff '.xml']);
    end
  else
  % processing output case
    if nargout>1, varargout{2} = QAR; end
    if nargout>0, varargout{1} = QAS; end 
  end
end
%==========================================================================
function [QAS, QAR, qamat, qamatm, mqamatm] = ...
  updateQAstructure(QAS, QASfi, QAR, QARfi, qamat, qamatm, mqamatm, QMAfn, fi)
%updateQAstructure. Update the quality structure.

  try
    QAS = cat_io_updateStruct(QAS,QASfi,0,fi);
    QAR = cat_io_updateStruct(QAR,QARfi,0,fi);
  catch
    fprintf('ERROR-Struct');
  end
  
  % color for the differen mark cases 
  for fni = 1:numel(QMAfn)
    if isfield(QAS(fi).qualitymeasures,QMAfn{fni})
      try
        qamat(fi,fni)  = QAS(fi).qualitymeasures.(QMAfn{fni});
        qamatm(fi,fni) = QAR(fi).qualityratings.(QMAfn{fni});
      catch
        qamat(fi,fni)  = QASfi.qualitymeasures.(QMAfn{fni});
        qamatm(fi,fni) = QARfi.qualityratings.(QMAfn{fni});
      end
    end    
  end
  try
    mqamatm(fi,1) = QAR(fi).qualityratings.IQR;
  catch
    mqamatm(fi,1) = QASfi.qualityratings.IQR;
  end

  if size(mqamatm,2)==2
    try
      mqamatm(fi,2) = QAR(fi).qualityratings.SIQR;
    catch
      mqamatm(fi,2) = QARfi.qualityratings.SIQR;
    end
  end

  mqamatm(fi,:) = max(0,min(10.5, mqamatm(fi,:)));
          
end
%==========================================================================
function [Yp0,Ym,Vo,p0rmse] = getImages(Pp0,Po,Pm,fi)
%getImages. Load major images for QC processing.
%
%  [Yp0,Ym,Vo] = getImages(Pp0,Po,Pm,fi)
%
%  Yp0 .. (resampled) segmentation as labelmap (0-BG, 1-CSF, 2-GM, 3-WM)
%  Ym  .. bias-corrected image
%  Vo  .. original image

  if numel(Pp0{fi})>1 && Pp0{fi}(end-1)==',', Pp0{fi}(end-1:end) = []; end
  if numel(Po{fi})>1  && Po{fi}(end-1)==',',  Po{fi}(end-1:end) = []; end
  if numel(Pm{fi})>1  && Pm{fi}(end-1)==',',  Pm{fi}(end-1:end) = []; end

  if (isempty(Po{fi}) || ~exist(Po{fi},'file')) && ...
      ~isempty(Pm{fi}) && exist(Pm{fi},'file')
    Po{fi} = Pm{fi}; 
    if 0 %cat_get_defaults('extopts.expertgui') % RD202508: deactived as we only use the original now
      cat_io_cprintf('warn','Cannot find/open original image use bias corrected:   %80s\n',Pm{fi}) 
    end
  elseif (isempty(Pm{fi}) || ~exist(Pm{fi},'file')) && ...
         ~isempty(Po{fi}) && exist(Po{fi},'file')
    if exist( [Pm{fi} '.gz'] , 'file')
      Pm{fi} = [Pm{fi} '.gz']; 
    else
      Pm{fi} = Po{fi}; 
      if 0 %cat_get_defaults('extopts.expertgui') % RD202508: deactived as we only use the original now
        cat_io_cprintf('warn','Cannot find/open bias corrected image use original:   %80s\n',Po{fi}) 
      end
    end
  end
  if 0 %isempty(Pp0{fi}) || ~exist(Pp0{fi},'file')
    % RD202508: deactived as we run a simple segmentation + warning elsewhere  
    cat_io_cprintf('err','Cannot find/open segmentation: \n  %s\n',Pp0{fi}) 
  end
  

  % handle gzipped original data 
  [pp,ff,ee] = spm_fileparts(Po{fi});
  if exist(fullfile(pp,[ff ee]),'file')
    evalc('Vo  = spm_vol(Po{fi});');
  %elseif exist(fullfile(pp,[ff ee '.gz']),'file')
  %  gunzip(fullfile(pp,[ff ee '.gz']));
  %  Vo  = spm_vol(Po{fi});
  %  delete(fullfile(pp,[ff ee '.gz'])); 
  else
    error('cat_vol_qa:noYo','No original image.');
  end

  % load further images - bias corrected Ym and segmentation Yp0
  evalc('Vm  = spm_vol(Pm{fi})');
  [pp0,ff0,ee0] = spm_fileparts(Pp0{fi});
  if cat_io_contains(ff0,'qcseg')
    Yp0 = []; 
  else
    switch ff0(1:2) 
      case {'p0','sy','c0'} % cat and other label map (use c0 for SPM)
        evalc('Vp0 = spm_vol(Pp0{fi});');
        if ~isempty(Vm) && any(Vp0.dim ~= Vm.dim)
          [Vx,Yp0] = cat_vol_imcalc(Vp0,Vm,'i1',struct('interp',2,'verb',0)); % linear interp
        else
          evalc('Yp0 = single(spm_read_vols(Vp0))');
        end
      case 'c1'
        tval = [2 3 1];
        for ci = 1:3
          evalc('Vc = spm_vol(fullfile(pp0,sprintf(''c%d%s%s'',ci,ff0(3:end),ee0)));');
          if ci == 1
            evalc('Yp0 = tval(ci) * single(spm_read_vols(Vc))');
          else
            evalc('Yp0 = Yp0 + tval(ci) * single(spm_read_vols(Vc))');
          end
        end 
      otherwise
        Yp0 = []; 
    end
  end
  Yp0(isnan(Yp0) | isinf(Yp0)) = 0; 
   
  % Internal bias correction to handle the original images as the processed
  % bias corrected image in CAT is also intensity normalized and endoised.
  vx_vol  = sqrt(sum(Vo.mat(1:3,1:3).^2));
  evalc('Ym  = single(spm_read_vols(spm_vol(Po{fi})))');
  %WMth = cat_stat_nanmedian(Ym(Yp0>2.9)); 
  if ~isempty(Yp0)
    Ym(isnan(Yp0) | isinf(Yp0)) = 0; 
    %Yw = cat_vol_morph( Yp0>2.95 , 'e',0,vx_vol)  &  cat_vol_morph( Yp0>2.25 , 'e',1,vx_vol); RD20250907: erode to hard for SPM 
    Yw  = Yp0>2.9; % ### RD20250910: however this is a bit too simple in case of segmentation bugs and a check for local outliers would be good
    Yb  = cat_vol_approx( Ym .* Yw + Yw .* min(Ym(:)) ,'rec') - min(Ym(:)); 
    Ym  = Ym ./ max(eps,Yb);
  end
  
  % Detection of possible preprocessing issues in case the segmentation 
  % varies strongly from the original image. We intensity normalize here 
  % the Yp0 map to fit to the Ym. WMHs could cause problems but this is 
  % even intended.
  if ~isempty(Yp0)
    Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));
    T3th = zeros(1,3); for ci = 1:3, T3th(ci) = cat_stat_nanmedian(Ym(  Yp0toC(Yp0(:),ci) > 0.95 )); end 
    Yp0t = zeros(size(Yp0)); for ci = 1:3, Yp0t = Yp0t + T3th(ci) .* Yp0toC(Yp0,ci); end
    p0rmse = (cat_stat_nanmean(Ym(Yp0(:)>0) - Yp0t(Yp0(:)>0)).^2)^0.5; 
    if p0rmse>0.2
      cat_io_cprintf('warn', ['Segmentation is maybe not fitting to the image ' ...
        '(RMSE(Ym,Yp0)=%0.2f)?:\n  %s\n  %s\n'], p0rmse,Pm{fi}, Pp0{fi}); 
    end      
  end
end
%==========================================================================
function [QAS,QAR] = cat12err(opt,mrifolder,reportfolder)
%cat12err. Create short report in case of CAT preprocessing error. 
% This report contain basic parameters used for the CAT error report
% creation in cat_io_report.

  % file information
  % -----------------------------------------------------------------------
  [pp,ff,ee]          = spm_fileparts(opt.job.channel.vols{opt.subj});
  [QAS.filedata.path,QAS.filedata.file] = ...
                        spm_fileparts(opt.job.channel.vols{opt.subj});
  QAS.filedata.fname  = opt.job.data{opt.subj};
  QAS.filedata.F      = fullfile(pp,[ff ee]);
  QAS.filedata.Fm     = fullfile(pp,mrifolder,['m'  ff ee]);
  QAS.filedata.Fp0    = fullfile(pp,mrifolder,['p0' ff ee]);
  QAS.filedata.fnames = [ 
    spm_str_manip(pp,sprintf('k%d', ...
      floor( max(opt.snspace(1)-19-ff,opt.snspace(1)-19)/3) - 1)), '/',...
    spm_str_manip(ff,sprintf('k%d',...
     (opt.snspace(1)-19) - floor((opt.snspace(1)-14)/3))), ...
     ];
  
  % load header for resolution 
  if exist(QAS.filedata.F,'file')
    V = spm_vol(QAS.filedata.F); 
  elseif exist(QAS.filedata.Fm,'file')
    V = spm_vol(QAS.filedata.Fm); 
  elseif exist(QAS.filedata.Fp0,'file')
    V = spm_vol(QAS.filedata.Fp0); 
  else
    vx_vol = [0 0 0];
  end
  if ~exist('vx_vol','var')
    vx_vol = sqrt(sum(V.mat(1:3,1:3).^2));
  end

  % software, parameter and job information
  % -----------------------------------------------------------------------
  [ver_cat, rev_cat] = cat_version;
  ver_cat = ver_cat(4:end); % remove leading CAT
  [nam,rev_spm] = spm('Ver');
  QAS.software.version_spm = rev_spm;
  if strcmpi(spm_check_version,'octave')
    QAS.software.version_octave = version;  
  else
    A = ver;
    for i=1:length(A)
      if strcmp(A(i).Name,'MATLAB')
        QAS.software.version_matlab = A(i).Version; 
      end
    end
    clear A
  end

  % 1 line: Matlab, SPM, CAT version number and GUI and experimental mode 
  if ispc,      OSname = 'WIN';
  elseif ismac, OSname = 'MAC';
  else,         OSname = 'LINUX';
  end
  
  QAS.software.system       = OSname;
  QAS.software.version_cat  = ver_cat;
  if ~isfield(QAS.software,'version_segment')
    QAS.software.version_segment = rev_cat;
  end
  QAS.software.revision_cat = rev_cat;
  try
    QAS.hardware.numcores = max(cat_get_defaults('extopts.nproc'),1);
  catch
    QAS.hardware.numcores = 1;
  end
  
  % save important preprocessing parameters 
  QAS.parameter.opts        = opt.job.opts;
  QAS.parameter.extopts     = rmfield(opt.job.extopts,...
    {'LAB','atlas','satlas','darteltpms','shootingtpms','fontsize'});
  QAS.parameter.caterr      = opt.caterr; 
  QAS.error                 = opt.caterrtxt; 

  % redo the quality rating, ie. measure scaling
  QAS.qualitymeasures.NCR         = NaN; 
  QAS.qualitymeasures.ICR         = NaN; 
  QAS.qualitymeasures.contrastr   = NaN; 
  QAS.qualitymeasures.res_ECR     = NaN; 
  QAS.qualitymeasures.res_RMS     = cat_stat_nanmean(vx_vol.^2).^0.5;
  QAS.subjectmeasures.vol_rel_CGW = [NaN NaN NaN]; 
  QAS.subjectmeasures.SQR         = NaN; 
  
  QAR = upate_rating(QAS,opt.version);
      
  % export 
  if opt.write_xml
    cat_io_xml(fullfile(pp,reportfolder,[opt.prefix ff '.xml']),QAS,'write+');    
    cat_io_xml(fullfile(pp,reportfolder,[opt.prefix ff '.xml']),QAR,'write+');    
  end
end
%==========================================================================
function def = defaults
%default. cat_vol_qa default parameters. 

  def.verb       = 2;         % verbose level  [ 0=nothing | 1=points | 2*=results ]
  def.write_csv  = 2;         % final cms-file [ 0=do not write |1=write | 2=overwrite ] 
  def.write_xml  = 1;         % images base xml-file
  def.sortQATm   = 1;         % sort QATm output
  def.orgval     = 0;         % original QAM results (no marks)
  def.prefix     = 'cat_';    % prefix of QC variables
  def.method     = 'spm';     % used 
  def.snspace    = [100,7,3];
  %def.nogui      = exist('XT','var');
  def.rerun      = 1;         % 0-load if exist, 1-reprocess if "necessary", 2-reprocess 
  def.version    = 'cat_vol_qa201901x';
  def.MarkColor  = cat_io_colormaps('marks+',40); 
  def.versions0  = {'cat_vol_qa201602'};  % no ECR
  def.versions1  = {'cat_vol_qa201602','cat_vol_qa201901','cat_vol_qa202110','cat_vol_qa202205'};  % no FEC
end
%==========================================================================
function QARfi = upate_rating(QASfi,version,getdef)
% update rating of stable versions 

  if ~exist('getdef','var'), getdef = 0; end

  ndef =  cat_stat_marks('default');
%version = 'default'; 
  switch version
    case {'cat_vol_qa201602','cat_vol_qa201901'}
      % robust versions with minimal changes/differences 
      ndef.noise = [0.046797 0.397905]; 
      ndef.bias  = [0.338721 2.082731];
    case 'cat_vol_qa201901x'  
      % final refined version of robust version 201901 ###############
      ndef.noise = [  0.0183,   0.0868]; 
      ndef.bias  = [  0.2270,   1.3949];
      ndef.ECR   = [  0.0202,   0.1003]; 
      ndef.FEC   = [130.0000, 470.0000]; 
    case {'cat_vol_qa202110'}
      % changes between 2019/01 and 2021/10 result a bit different version 
      % partial better, partially worse (regular successor of version 201901) 
      ndef.noise = [0.054406 0.439647]; 
      ndef.bias  = [0.190741 1.209683];
    case {'cat_vol_qa202110x'}
      % revised 202110 version
      ndef.noise = [    0.0183 0.1362]; 
      ndef.bias  = [    0.1823 1.1144];
      ndef.ECR   = [   -0.0041 1.1918]; 
    case {'cat_vol_qa202205'}
      % latest regular QC version as successor of the 202110 (~202205)
      ndef.noise = [0.056985 0.460958];
      ndef.bias  = [0.187620 1.206548]; 
    case {'cat_vol_qa202207b'}
      % attempt to reorganize the QC and improve accuracy that was finally 
      % not successful 
      ndef.noise = [0.026350 0.203736]; 
      ndef.bias  = [0.120658 0.755107]; 
    case {'cat_vol_qa202301'} 
      % latest reworked QC version as successor of the standard qa with
      % 202207b as mid version that was not successful 
      ndef.noise = [    0.0255 0.1887]; 
      ndef.bias  = [    0.1996 1.2632]; 
      ndef.ECR   = [    0.0596 1.1890]; 
    case {'cat_vol_qa202310', 'cat_vol_qa'}
      % the 202310(dd) represents a new start ###############
      ndef.noise = [    0.0326 0.2417]; %0.0322 0.2418]; 
      ndef.bias  = [    0.1781 0.9393]; %0.1784 0.9356]; 
      ndef.ECR   = [  0.3597,   1.4316]; %0.3489 1.4117]; %-0.0364 1.1497];
      ndef.FEC   = [110.0000, 480.0000];
    case 'cat_vol_qa202412'  
      % final refined version of robust version 201901 ###############
      ndef.noise = [  0.0095,   0.0740]; %[  0.0172,   0.1234];
      ndef.bias  = [  0.1695,   0.9907]; %[  0.1668,   1.0234];
      ndef.ECR   = [  0.3859,   1.2101]; %[  0.4141,   1.5532]; %-0.0364 1.1497];
      ndef.FEC   = [160.0000, 420.0000]; %[100.0000, 630.0000];
    case 'default'
    otherwise
      warning('missing scaling definition'); 
      ndef.noise = [0.5  0.5]; 
      ndef.bias  = [0.33 2.0];
  end
 
  ndef.QS{find(cellfun('isempty',strfind(ndef.QS(:,2),'NCR'))==0,1),4} = ndef.noise; %#ok<*STRCL1> 
  ndef.QS{find(cellfun('isempty',strfind(ndef.QS(:,2),'ICR'))==0,1),4} = ndef.bias;
  if isfield(ndef,'FEC'),   ndef.QS{find(cellfun('isempty',strfind(ndef.QS(:,2),'FEC'))==0,1),4} = ndef.FEC; end
  if isfield(ndef,'ECR'),   ndef.QS{find(cellfun('isempty',strfind(ndef.QS(:,2),'ECR'))==0,1),4} = ndef.ECR; end
 
  if ~getdef
    QARfi = cat_stat_marks('eval',1,QASfi,ndef);
  else
    QARfi = ndef;
  end
end

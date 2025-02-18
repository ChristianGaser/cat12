function spm_cat12(varargin)
% ______________________________________________________________________
% CAT12 Toolbox wrapper to start CAT with different user modes or 
% default files.  Changing the user mode requires restarting of CAT and
% SPM.  The expert user mode allows to control further parameters and  
% semi-evaluated functions, whereas the developer mode contain parameter
% for internal tests and unsafe functions.
% 
%   cat12(action)
%   
%   CAT user modes:
%     action = ['default','expert','developer'] 
%
%   CAT default files for other species (in development):
%     action = ['oldwoldmonkeys'|'greaterapes']
%
%   CAT start with own default files:
%     action = 'select' 
%     action = 'mypath/cat_defaults_mydefaults'
%
% ______________________________________________________________________
% Christian Gaser, Robert Dahnke
% $Id$


rev = '$Rev$';
global deffile;
global cprintferror;  % temporary, because of JAVA errors in cat_io_cprintf ... 20160307
%try clearvars -global deffile;  end %#ok<TRYNC>

% start cat with different default file
catdir = fileparts(mfilename('fullpath')); 
catdef = fullfile(catdir,'cat_defaults.m');

% add path for octave mex-files
if strcmpi(spm_check_version,'octave')
  if ismac
    uts = uname;
    if strcmp(uts.machine,'x86_64')
      addpath(fullfile(catdir,'mexmaci'));
    else
      addpath(fullfile(catdir,'mexmaca'));
    end
  elseif isunix
    addpath(fullfile(catdir,'mexa64'));
  elseif ispc
    addpath(fullfile(catdir,'mexw64'));
  end
end

% check that CAT12 is installed in the correct folder
pth = fileparts(mfilename('fullpath'));
[pth2, nam] = fileparts(pth);
if ~strcmp(nam,'cat12')
  spm('alert!',sprintf('Please check that you do not have multiple CAT12 installations in your path!\nYour current CAT12 version is installed in %s but should be installed in %s',pth,fullfile(catdir)),'WARNING');
end

% check that mex-files on MAC are not blocked
try
  feval(@cat_sanlm,single(rand(6,6,6)),1,3);
catch
  if strcmpi(spm_check_version,'octave')
    oldpath = pwd;
    cd(pth)
    compile
    feval(@cat_sanlm,single(rand(6,6,6)),1,3);
    cd(oldpath)
    
    % check that patches and updates exist
    if ~exist('savexml') || ~exist('fcnchk')
      error('Please update and patch SPM12 first')
    end
  elseif ismac 
    CATDir = fullfile(catdir);
    web('https://en.wikibooks.org/wiki/SPM/Installation_on_64bit_Mac_OS_(Intel)#Troubleshooting');
    cat_io_cmd(sprintf('\nThe following commands will be executed as administrator to allow execution of CAT12 binaries and mex-files.\n Please now type admin password to call sudo\n'),'warning');
    cat_io_cmd(sprintf('You can also break that command here and run the commands that are listed on the open website under Troubleshooting manually.\n'),'warning');
    cmd = ['sudo xattr -r -d com.apple.quarantine ' CATDir];
    system(cmd); fprintf([cmd '\n']);
    cmd = ['sudo find ' CATDir ' -name *.mexmac* -exec spctl --add {} \;'];
    system(cmd); fprintf([cmd '\n']);
    cmd = ['sudo chmod a+x ' CATDir '/CAT.mac*/CAT*'];
    system(cmd); fprintf([cmd '\n']);
    cmd = ['sudo find ' CATDir ' -name *.mexmac* -exec xattr -d com.apple.quarantine {} \;'];
  end
end


% get expert level except for standalone installation
expert   = cat_get_defaults('extopts.expertgui'); 

if nargin==0 && (isempty(deffile) || strcmp(deffile,catdef))
  deffile = catdef; 
  if ~strcmp(cat_get_defaults('extopts.species'),'human') || expert>0 
    restartspm = 1;
  else
    restartspm = 0;
  end
elseif nargin==1 
  deffile = varargin{1}; 
  restartspm = 1;
elseif nargin==2
  deffile = varargin{1};
  catdef  = varargin{2};
  restartspm = 1;
else
  deffile = catdef; 
  restartspm = 1;
end


% choose filesspecies 
speciesdisp = '';
switch cat_get_defaults('extopts.species')
  case {'select','human','default','expert','developer'}
    % nothing to do
  otherwise
    % load default to remove animal settings
    try clearvars -global cat; end %#ok<TRYNC>
    [deffile_pp,deffile_ff] = fileparts(catdef);
    oldwkd = cd; 
    cd(deffile_pp);
    try clearvars -global cat; end %#ok<TRYNC>
    clear cat;
    eval(deffile_ff);
    eval('global cat;'); 
    cd(oldwkd);
end
switch lower(deffile) 
  case 'select'
    deffile = spm_select(1,'batch','Select CAT default file!','',catdir);
    if isempty(deffile) 
      return
    end
  case {'default','human'}
    mycat  = cat_get_defaults; 
    mycat.extopts.expertgui = 0;
    restartspm = 1;
    deffile = catdef; 
  case 'expert'
    mycat  = cat_get_defaults; 
    mycat.extopts.expertgui = 1;
    restartspm = 1;
    deffile = catdef; 
  case 'developer'
    mycat  = cat_get_defaults; 
    mycat.extopts.expertgui = 2;
    restartspm = 1;
    deffile = catdef; 
  case {'greaterapes','lesserapes','oldworldmonkeys','newworldmonkeys','mammals','chimpanzees','dogs',...
        'greaterape' ,'lesserape' ,'oldworldmonkey' ,'newworldmonkey', 'mammal', 'chimpanzee' ,'dog', ...
        'baboons', 'macaques', ...
        'baboon' ,'macaca', 'macaque'}
    switch lower(deffile)
      case {'greaterapes','greaterape'},          species = 'ape_greater';     speciesdisp = ' (greater apes)';
      %case {'lesserapes','lesserape'},            species = 'ape_lesser';      speciesdisp = ' (lesser apes)';
      case {'oldworldmonkeys','oldworldmonkey'},  species = 'monkey_oldworld'; speciesdisp = ' (oldworld monkeys)';
      %case {'newworldmonkeys','newworldmonkey'},  species = 'monkey_newworld'; speciesdisp = ' (newworld monkeys)';
      %case {'mammals','mammal'},                  species = 'mammal';          speciesdisp = ' (mammal)';
      case {'chimpanzees','chimpanzee'},          species = 'chimpanzee';      speciesdisp = ' (chimpanzee)';
      case {'macaque','macaques'},                species = 'macaque';         speciesdisp = ' (macaque)';
      case {'baboons','baboon'},                  species = 'baboon';          speciesdisp = ' (baboon)';
      case {'dogs','dog'},                        species = 'dog';             speciesdisp = ' (dogs)';
      otherwise
        error('CAT:unreadySpecies','Templates of species "%s" are not ready yet.\n',deffile);
    end  
    
    mycat                      = cat_get_defaults;
    % change TPM and user higher resolution and expect stronger bias
    mycat.opts.tpm             = {fullfile(catdir,'templates_animals',[species '_TPM.nii'])};
    mycat.opts.biasreg         = 0.001;                                  % less regularisation 
    mycat.opts.biasfwhm        = 30;                                     % stronger fields 
    mycat.opts.samp            = 1;                                      % smaller resampling
    mycat.opts.affreg          = 'none';                                 % no affine regularisation
    % use species specific templates, higher resolution, stronger corrections and less affine registration (by SPM) 
    mycat.extopts.species      = species;  
    %mycat.extopts.brainscale   = 200; % non-human brain volume in cm3 (from literature) or scaling in mm (check your data)
    mycat.extopts.darteltpm    = {fullfile(catdir,'templates_animals',[species '_Template_1.nii'])}; % Indicate first Dartel template
    mycat.extopts.shootingtpm  = {fullfile(catdir,'templates_animals',[species '_Template_0_GS.nii'])}; % Indicate first Shooting template
    mycat.extopts.cat12atlas   = {fullfile(catdir,'templates_animals',[species '_cat.nii'])};        % VBM atlas with major regions for VBM, SBM & ROIs
    mycat.extopts.brainmask    = {fullfile(catdir,'templates_animals',[species '_brainmask.nii'])};  % brainmask for affine registration
    mycat.extopts.T1           = {fullfile(catdir,'templates_animals',[species '_T1.nii'])};         % T1 for affine registration
    mycat.extopts.sanlm        = 3;                                     % ISARNLM for stronger corrections
    mycat.extopts.restype      = 'best';        
    mycat.extopts.resval       = [0.50 0.30];                           % higher internal resolution 
    %mycat.extopts.APP          = 1070;                                  % less affine registration, but full corrections (by SPM)
    mycat.extopts.vox          = 0.50;                                  % voxel size for normalized data 
    mycat.extopts.bb           = [[-inf -inf -inf];[inf inf inf]];      % template default
    mycat.extopts.WMHC         = 0;                                     % not in primates yet
    mycat.extopts.expertgui    = 2;                                     % set to expert later ...
    mycat.extopts.ignoreErrors = 1;  
    switch species
      case 'monkey_oldworld'
        mycat.extopts.atlas = { ... 
          fullfile(catdir,'templates_animals','monkey_oldworld_atlas_inia19NeuroMaps.nii') 1 {'csf','gm','wm'} 1; 
          };
      case 'chimpanzee'
        mycat.extopts.atlas = { ... 
          fullfile(catdir,'templates_animals','chimpanzee_atlas_davi.nii') 1 {'csf','gm','wm'} 1; 
          };
      otherwise
        mycat.extopts.atlas = {}; 
        mycat.output.ROI    = 0;
    end
    
    restartspm = 1;
    deffile    = catdef; 
  otherwise
    % lazy input - no extension 
    [deffile_pp,deffile_ff,deffile_ee] = fileparts(deffile);
    if isempty(deffile_ee)
      deffile_ee = '.m';
    end
    % lazy input - no directory
    if isempty(deffile_pp) 
      if exist(fullfile(pwd,deffile_ff,deffile_ee),'file') 
        deffile_pp = pwd; 
      else
        deffile_pp = fullfile(catdir); 
      end
    end
    deffile = fullfile(deffile_pp,[deffile_ff,deffile_ee]); 

    if isempty(deffile) || ~exist(deffile,'file')
      help spm_cat12;
      error('CAT:unknownDefaultFile','Unknown action or nonexisting default file "%s".\n',deffile);
    end
end


% The cat12 global variable is created and localy destroyed, because we want to call the cat12 function. 
if exist('mycat','var') 
  try clearvars -global cat; end %#ok<TRYNC>
  eval('global cat; cat = mycat;'); 
else
  % set other defaultfile
  oldwkd = cd; 
  cd(deffile_pp);
  try clearvars -global cat; end %#ok<TRYNC>
  clear cat;
  eval(deffile_ff);
  eval('global cat;'); 
  cd(oldwkd);
end

% initialize SPM 
eval('global defaults;'); 
% this is required to initialize the atlas variable for default users
if restartspm 
  clear defaults; 
  spm_jobman('initcfg');
end
clear cat;

% initialize atlas variable 
exatlas  = cat_get_defaults('extopts.atlas'); 
for ai = 1:size(exatlas,1)
  if exatlas{ai,2}<=expert && exist(exatlas{ai,1},'file')
    [pp,ff,ee]  = spm_fileparts(exatlas{ai,1}); 

    % if output.atlases.ff does not exist then set it by the default file value
    if isempty(cat_get_defaults(['output.atlases.' ff]))
      cat_get_defaults(['output.atlases.' ff], exatlas{ai,4})
    end
  end
end

exsatlas  = cat_get_defaults('extopts.satlas'); 
for ai = 1:size(exsatlas,1)
  if exsatlas{ai,3}<=expert && exist(exsatlas{ai,2},'file')
    name = exsatlas{ai,1};

    % if output.atlases.ff does not exist then set it by the default file value
    if isempty(cat_get_defaults(['output.satlases.' name]))
      cat_get_defaults(['output.satlases.' name], exsatlas{ai,4})
    end
  end
end

 

% temporary, because of JAVA errors in cat_io_cprintf ... 20160307
if expert<2
  cprintferror=1;
end

spm('FnBanner',mfilename,cat_version);
[Finter,Fgraph] = spm('FnUIsetup','CAT12.9');
url = fullfile(fileparts(mfilename('fullpath')),'doc','cat.html');

% open interactive help for newer version because display of html pages does not work anymore
if cat_io_matlabversion > 20212
  % SPM header image
  Pposter = fullfile( catdir, 'doc', 'images', 'CAT_Poster.jpg'); 
  F = spm_figure('GetWin'); 
  spm_figure('clear',F); 
  Fpos = get(F,'Position'); 
  h = imshow(imread(Pposter)); 
  set(get(h,'Parent'),'Position',[0 0 1 1]);
  set(F,'Position',Fpos);
  
  % CAT help 
  pos = Fgraph.Position(3:4);
  web(url,'-noaddressbox','-new')
else
  spm_help('!Disp',url,'',Fgraph,'Computational Anatomy Toolbox for SPM12');
end

% check that binaries for surface tools are running 
cat_system('CAT_3dVol2Surf');

%% add some directories 
spm_select('PrevDirs',{fullfile(catdir)});

%% command line output
cat_io_cprintf('silentreset');
switch expert
  case 0, expertguitext = '';
  case 1, expertguitext = ['Expert Mode' speciesdisp];
  case 2, expertguitext = ['Developer Mode' speciesdisp];
end
cat_io_cprintf([0.0 0.0 0.5],sprintf([ ...
    '\n' ...
    '   _______  ___  _______    \n' ...
    '  |  ____/ / _ \\\\ \\\\_   _/   ' expertguitext '\n' ...
    '  | |___  / /_\\\\ \\\\  | |     Computational Anatomy Toolbox\n' ...
    '  |____/ /_/   \\\\_\\\\ |_|     CAT12.9 - https://neuro-jena.github.io\n\n']));
cat_io_cprintf([0.0 0.0 0.5],' CAT default file:\n\t%s\n\n',deffile); 

% call GUI
cat12('fig'); 

% force use of PET modality in SPM12 to avoid problems of very low variance in spm_spm.m
warning off
spm('chmod','PET');
warning on

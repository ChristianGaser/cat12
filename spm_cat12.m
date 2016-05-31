function spm_cat12(varargin)
% ______________________________________________________________________
% CAT12 Toolbox wrapper to call CAT functions.
% 
%   spm_cat12 
%     .. start with CAT default parameter file
%   spm_cat12('gui')
%     .. start with default file of another species (in development)
%   spm_cat12(species) 
%     .. start with default file of another species (in development)
%        species = ['oldwoldmonkey'|'newwoldmonkey'|'greaterape'|'lesserape']
%   spm_cat12('mypath/cat_defaults_mydefaults') 
%     .. start CAT with another default parameter file
% ______________________________________________________________________
% Christian Gaser
% $Id$

% ______________________________________________________________________
% Development:
%   spm_cat12('mypath/cat_defaults_mydefaults',1) 
%     .. restart SPM for GUI updates
% ______________________________________________________________________

rev = '$Rev$';
global deffile;
global cprintferror;  % temporary, because of JAVA errors in cat_io_cprintf ... 20160307
%try clearvars -global deffile;  end %#ok<TRYNC>

% start cat with different default file
catdir = fullfile(spm('dir'),'toolbox','cat12'); 
catdef = fullfile(catdir,'cat_defaults.m');
if nargin==0 && (isempty(deffile) || strcmp(deffile,catdef))
  deffile = catdef; 
  restartspm = 0;
elseif nargin==1 
  deffile = varargin{1}; 
  restartspm = 1;
else
  deffile = catdef; 
  restartspm = 1;
end


% choose files
switch lower(deffile) 
  case {'select','choose'}
    deffile = spm_select(1,'batch','Select CAT default file!','',catdir);
    if isempty(deffile) 
      return
    end
  case 'default'
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
  %{
  % GUI for primates requires updates of the default files and some tests  
  case 'gui'
    deffile = spm_input('Species class',1,'human|ape|monkey',...
      {'human','ape','monkey'},1);
    deffile = deffile{1}; 
    
    switch lower(deffile)
      %case 'human'
      %  deffile = spm_input('Species class','+1','adult|child|neonate|fetus|other',...
      %    {'human_adult','human_child','human_neonate','human_fetus','human_other'},1); 
      %  deffile = deffile{1};
      case 'ape'
        deffile = spm_input('Species class','+1','greater|lesser|other',...
          {'ape_greater','ape_lesser','other'},1);
        deffile = deffile{1};
      case 'monkey'
        deffile = spm_input('Species class','+1','old world|new world|other',...
          {'monkey_oldworld','monkey_newworld','other'},1);
        deffile = deffile{1};
    end
    %}
end


switch lower(deffile)
  case 'human'
    deffile = catdef; 
  case {'monkey_oldworld','oldwoldmonkey','cat_defaults_monkey_oldworld','cat_defaults_monkey_oldworld.m'}
    deffile = fullfile(catdir,'templates_animals','cat_defaults_monkey_oldworld.m');
  case {'monkey_newworld','newworldmonkey','cat_defaults_monkey_newworld','cat_defaults_monkey_newworld.m'}
    deffile = fullfile(catdir,'templates_animals','cat_defaults_monkey_newworld.m');
  case {'ape_greater','greaterape','cat_defaults_ape_greater','cat_defaults_ape_greater.m'}
    deffile = fullfile(catdir,'templates_animals','cat_defaults_ape_greater.m');
  case {'ape_lesser','lesserape','cat_defaults_ape_lesser','cat_defaults_ape_lesser.m'}
    deffile = fullfile(catdir,'templates_animals','cat_defaults_ape_lesser.m');
end

if exist('mycat','var') 
  try clearvars -global cat; end %#ok<TRYNC>
  eval('global cat; cat = mycat;'); 
else
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
      deffile_pp = fullfile(spm('dir'),'toolbox','cat12'); 
    end
  end
  deffile = fullfile(deffile_pp,[deffile_ff,deffile_ee]); 

  % check if file exist
  if ~exist(deffile,'file')
    error('CAT:miss_cat_default_file','Can''t find CAT default file "%"','deffile'); 
  end

  % set other defaultfile
  % The cat12 global variable is created and localy destroyed, because we 
  % want to call the cat12 function. 
  %if 1 %nargin>0 %~strcmp(catdef,deffile) 
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
if isempty(defaults) || (nargin==2 && varargin{2}==1) || restartspm
  clear defaults; 
  spm_jobman('initcfg');
end
clear cat;

% temporary, because of JAVA errors in cat_io_cprintf ... 20160307
if cat_get_defaults('extopts.expertgui')<2
  cprintferror=1;
end


spm('FnBanner',mfilename,rev);
[Finter,Fgraph] = spm('FnUIsetup','CAT12');
url = fullfile(spm('Dir'),'toolbox','cat12','html','cat.html');
spm_help('!Disp',url,'',Fgraph,'Computational Anatomy Toolbox for SPM12');

[ST, RS] = cat_system('CAT_DumpCurv -h');
% because status will not give 0 for help output we have to check whether we can find the
% keyword "Usage" in output
if isempty(strfind(RS,'Usage'));
  if ispc
    [ST, RS] = system('systeminfo.exe');
  else
    [ST, RS] = system('uname -a');
  end
  cat_io_cmd(sprintf('\nWARNING: Surface processing will not work because CAT-binaries are not compatible to your system:\n%s\n',RS),'warning');
  fprintf('\n\nFor future support of your system please send this message to christian.gaser@uni-jena.de\n\n');
end

%% command line output
switch cat_get_defaults('extopts.expertgui')
  case 0, expertguitext = '';
  case 1, expertguitext = 'Expert Mode';
  case 2, expertguitext = 'Developer Mode';
end
cat_io_cprintf([0.0 0.0 0.5],sprintf([ ...
    '\n' ...
    '   _______  ___  _______    \n' ...
    '  |  ____/ / _ \\\\ \\\\_   _/   ' expertguitext '\n' ...
    '  | |___  / /_\\\\ \\\\  | |     Computational Anatomy Toolbox\n' ...
    '  |____/ /_/   \\\\_\\\\ |_|     CAT12 - http://www.neuro.uni-jena.de\n\n']));
cat_io_cprintf([0.0 0.0 0.5],' CAT default file:\n\t%s\n\n',deffile); 

% call GUI
cat12('fig'); 
  


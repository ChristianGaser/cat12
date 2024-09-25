function cat_progress_bar(action,varargin)
% Display progress and remaining time with the same syntax as for spm_progress_bar.
%
% FORMAT cat_progress_bar('Init', n_steps [, process_name, type])
% Initialises progress tool with additional specific name "process_name".  
% The "type" can be used to use a command line counter ('cmd') instead of  
% the popup window or to switch off this function ('off'). 
%
% FORMAT cat_progress_bar('Set',value)
% Sets iteration.
%
% FORMAT cat_progress_bar('Clear')
% Clears the progress window.
% 
% Example 1 - percentage progess bar window: 
%   cat_progress_bar('Init', 10, 'Bar', 'bar')
%   for i=1:10, cat_progress_bar('Set',i), pause(.25); end
%   cat_progress_bar('Clear')
%  
% Example 2 - command line progress:
%   cat_progress_bar('Init', 10, 'CMD', 'cmd')
%   for i=1:10, cat_progress_bar('Set',i), pause(.25); end
%   cat_progress_bar('Clear')
%
% Example 3 - percentage command line progress:
%   cat_progress_bar('Init', 10, '', 'cmd%')
%   for i=1:10, cat_progress_bar('Set',i), pause(.25); end
%   cat_progress_bar('Clear')
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

persistent sum_time time_old n_iterations Fwaitbar progress_step Fwaitbartype Fwaitbarname

%if ~nargin, action = 'Init'; end % this cannot work as init requires input
if ~nargin, help cat_progress_bar; return; end

if strcmpi(action,'init')
  if nargin > 1
    % catch possible errors of older calls with additional field, eg. 
    %  cat_progress_bar('Init', 10 ,'CAT-Preprocessing','Volumes Complete');
    switch varargin{2}
      case {'bar','cmd','cmd%'}
        bartype = varargin{2};
      otherwise
        bartype = 'bar'; 
    end
  else
    bartype = 'bar';  
  end  
elseif strcmpi(action,'off') || strcmpi(action,'silent') || strcmpi(action,'quite') || strcmpi(action,'')
  return
else
  if exist('Fwaitbartype','var')
    if ischar(Fwaitbartype)
      bartype = Fwaitbartype;
    else
      if strcmpi(action,'clear')
        return
      else
        error('cat_progress_bar:nobar','No progess bar!');
      end
    end
  else
   if strcmpi(action,'clear')
     return
   else
     error('cat_progress_bar:nobar','No progess bar!');
   end
  end
end
switch bartype
    case {'bar','cmd','cmd%'}
      switch lower(action)
          % Initialise
          %-------------------------------------------------------------------
          case 'init'
              n_iterations = varargin{1};

              % from older versions that were not printed
              %{
              if nargin > 2
                arg3 = varargin{2};
              else
                arg3 = '';
              end
              %}

              if nargin > 1 
                arg2 = varargin{2};
                if ~strcmp(bartype,'bar') && arg2(end)~=':'
                  arg2(end+1) = ':'; 
                end
              else
                if strcmp(bartype,'bar')
                  arg2 = 'Computing';
                else
                  arg2 = ''; 
                end
              end
              
              Fwaitbarname = arg2; 
              Fwaitbartype = bartype;
              sum_time = 0;
              time_old = clock;
              progress_step = 1;
      
              switch Fwaitbartype
                case 'bar'
                  try
                    Fwaitbar = waitbar(0,arg2,'Name',arg2);
                  catch
                    Fwaitbar = waitbar(0,arg2);
                  end
                  
                  % don't know whether this works, but I have used the parameters from spm_progress_bar
                  set(Fwaitbar,'IntegerHandle','off','InvertHardcopy','on','PaperPositionMode','auto','Tag','Interactive');
      
                case 'cmd'
                  fprintf('% 6d/% 6d',0,n_iterations); 
                case 'cmd%'
                  fprintf('%s% 6.2f ',arg2,0); 
              end
              
          % Set
          %-------------------------------------------------------------------
          case 'set'
              iter = varargin{1};
              
              % save time if we don't have to show every step
              if rem(iter, progress_step), return; end
              
              % estimate time for remaining iterations
              diff_time = etime(clock, time_old);
              sum_time = sum_time + diff_time;
              avg_time = sum_time/iter;
              remain_time = avg_time*(n_iterations-iter);
      
              % if execution time is much shorter than displaying the ui we skip
              % most of that displaying to save time
              if diff_time < 0.05; progress_step = progress_step*2; end
              
              % add 0.5s to remaining time to prevent that at the end 0s is
              % displayed for a longer time
              if iter/n_iterations*100 > 1
                str = sprintf('%.f%% (%s remaining)',iter/n_iterations*100,time2str(0.5+remain_time));
              else
                str = sprintf('%.f%%',iter/n_iterations*100);
              end
              switch Fwaitbartype
                case 'bar'
                  if ishandle(Fwaitbar), waitbar(iter/n_iterations,Fwaitbar,str); end
                case 'cmd'
                  fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b% 6d/% 6d',min(999999,iter),n_iterations);
                case 'cmd%'
                  fprintf('\b\b\b\b\b\b\b% 6.2f%%',min(999,iter/n_iterations*100));
              end

              % save old values
              time_old = clock;
          
          % Clear
          %-------------------------------------------------------------------
          case 'clear'
              switch Fwaitbartype
                case 'bar'
                  if ishandle(Fwaitbar), delete(Fwaitbar); end
                case 'cmd'
                  fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b             \b\b\b\b\b\b\b\b\b\b\b\b\b');
                case 'cmd%'
                  fprintf(sprintf('%s',repmat('\b',1,numel(Fwaitbarname)))); 
                  fprintf('\b\b\b\b\b\b\b       \b\b\b\b\b\b\b\b');
              end
      
          % Error
          %-------------------------------------------------------------------
          otherwise
              error('Unknown action string');
      end
  otherwise
    error('error:cat_progress_bar:bartype',sprintf('Unknown bartype "%s"',bartype)); 
end

return

function str = time2str(t)
minutes = t/60;
hours   = t/3600;
days    = hours/24;

if days > 2
  str = sprintf('%d days %02.1f h', floor(days),24*(days-floor(days)));
elseif days > 1
  str = sprintf('%d day %02.1f h', floor(days),24*(days-floor(days)));
elseif hours > 1
  str = sprintf('%d:%02.0f h', floor(hours),60*(hours-floor(hours)));
elseif minutes > 1
  str = sprintf('%d:%02.0f min',floor(minutes),60*(minutes-floor(minutes)));
else
  str = sprintf('%d s',round(t));
end

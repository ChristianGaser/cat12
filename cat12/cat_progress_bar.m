function cat_progress_bar(action,varargin)
% Display progress and remaining time with the same syntax as for spm_progress_bar
%
% FORMAT cat_progress_bar('Init', n_steps, process_name)
% Initialises progress tool
%
% FORMAT cat_progress_bar('Set',value)
% Sets iteration.
%
% FORMAT cat_progress_bar('Clear')
% Clears the progress window.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id: 2558 2024-02-28 $

persistent sum_time time_old n_iterations Fwaitbar progress_step

if ~nargin, action = 'Init'; end


switch lower(action)
    % Initialise
    %-------------------------------------------------------------------
    case 'init'
        n_iterations = varargin{1};
        if nargin > 2, arg2 = varargin{2}; else arg2 = 'Computing';  end
        
        sum_time = 0;
        time_old = clock;
        progress_step = 1;

        try
          Fwaitbar = waitbar(0,arg2,'Name',arg2);
        catch
          Fwaitbar = waitbar(0,arg2);
        end
        
        % don't know whether this works, but I have used the parameters from spm_progress_bar
        set(Fwaitbar,'IntegerHandle','off','InvertHardcopy','on','PaperPositionMode','auto','Tag','Interactive');

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
        if ishandle(Fwaitbar), waitbar(iter/n_iterations,Fwaitbar,str); end
          
        % save old values
        time_old = clock;
    
    % Clear
    %-------------------------------------------------------------------
    case 'clear'
        if ishandle(Fwaitbar), delete(Fwaitbar); end

    % Error
    %-------------------------------------------------------------------
    otherwise
        error('Unknown action string');
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

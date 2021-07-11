function cat_io_senderrormail(xml_file,log_file,err_txt)
% ______________________________________________________________________
% Function that prepares a mail with errors in CAT. 
% xml_file - xml-file of errornous data
% log_file - log-file of errornous data
% err_txt  - optional error text 
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (http://www.neuro.uni-jena.de)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Revision$  $Date$

%#ok<*TRYNC>
  
  % program version
  [nam,rev_cat] = cat_version;
  [nam,rev_spm] = spm('Ver');

  % read xml-file
  fid = fopen(xml_file,'r');
  if nargin > 2
    err_txt = [err_txt '\n\n'];
  else
    err_txt = '';
  end
  while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    err_txt = [err_txt tline '\n'];
  end   
  fclose(fid);

  % initial message  
  emailSubject = sprintf('CAT %s error',rev_cat); 
  emailSubject = strrep(emailSubject,' ','%20');  
  mainBody     = sprintf('Hi Christian,\\n\\nan error has occurred in CAT %s with SPM %s with MATLAB %s under %s.\\n%s',...
    rev_cat , rev_spm, version('-release') , computer); 
    
  mainBody = sprintf('%s\\nBest regards\\n\\n\\n',mainBody);
  mainBody = convert_utf([mainBody err_txt]);
  
  % final creation of the mail 
  recipients = 'vbmweb@gmail.com';  
  
  %% open xml-file if mail does not work
  status = web(sprintf('mailto:"%s"?subject="%s"&body="%s"',recipients,emailSubject,mainBody),'-browser');
  if status
    % try  with shorter text
    status = web(sprintf('mailto:"%s"?subject="%s"&body="%s"',recipients,emailSubject,mainBody(1:8000)),'-browser');
    if status
      web(sprintf('mailto:"%s"?subject="%s"',recipients,emailSubject),'-browser');
    end
  end

  % show information anyway because sometimes mail does not open even if status is zero
  alert_txt = sprintf('Please send mail to %s and copy content of %s and %s into mail if not already done or attach these files.',recipients,xml_file,log_file);
  spm('alert',alert_txt);
  edit(xml_file);
  edit(log_file);
  fprintf('\n\n%s\n',alert_txt)  
end

function txt = convert_utf(txt)

  txt = strrep(txt,':','%20');
  txt = strrep(txt,'"','''');
  txt = strrep(txt,' ','%20');
  txt = strrep(txt,'\n','%0A');
  txt = strrep(txt,'!','%21');
  txt = strrep(txt,'#','%23');
  txt = strrep(txt,'%%','%25');
  txt = strrep(txt,'*','%2A');
  txt = strrep(txt,'/','%2F');
  txt = strrep(txt,'<','%3C');
  txt = strrep(txt,'>','%3E');
  txt = strrep(txt,'?','%3F');
  txt = strrep(txt,'(','%28');
  txt = strrep(txt,')','%29');
  txt = strrep(txt,'\\\\\\\\','\');
  txt = strrep(txt,'\\\\','\');
  txt = strrep(txt,'\\','\');
end
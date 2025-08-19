function cat_main_reportcmd(job,res,qa)
% ______________________________________________________________________
% 
%   Display report of CAT preprocessing in the command window and 
%   cleanup some figures. 
%
%   cat_main_reportcom(job,res,qa)
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

  VT0 = res.image0(1);
  [pth,nam] = spm_fileparts(VT0.fname); 

  % command window output
  QMC       = cat_io_colormaps('marks+',17);
  GMC       = cat_io_colormaps('turbo',45);
  GMC       = GMC ./ repmat( max(1,sum(GMC,2)) , 1 , 3);  % make bright values darker 
  color     = @(QMC,m) QMC(max(1,min(size(QMC,1),isnan(m) + round(((m-1)*3)+1))),:);
  colorgmt  = @(GMC,m) GMC(max(1,min(size(GMC,1),isnan(m) + round(((m-0.5)*10)+1))),:);
  mark2rps  = @(mark) min(100,max(0,105 - mark*10)) + isnan(mark).*mark;
  grades    = {'A+','A','A-','B+','B','B-','C+','C','C-','D+','D','D-','E+','E','E-','F'};
  mark2grad = @(mark) grades{max(min(numel(grades),max(max(isnan(mark)*numel(grades),1),round((mark+2/3)*3-3))))};
  
  % definition of subfolders
  [mrifolder, reportfolder, surffolder, labelfolder] = cat_io_subfolders(VT0.fname,job);

  mrifolder = fullfile(pth,mrifolder);
  [stat, val] = fileattrib(mrifolder);
  if stat, mrifolder = val.Name; end

  reportfolder = fullfile(pth,reportfolder);
  [stat, val] = fileattrib(reportfolder);
  if stat, reportfolder = val.Name; end

  surffolder = fullfile(pth,surffolder);
  [stat, val] = fileattrib(surffolder);
  if stat, surffolder = val.Name; end

  labelfolder = fullfile(pth,labelfolder);
  [stat, val] = fileattrib(labelfolder);
  if stat, labelfolder = val.Name; end

  % output conclusion
  fprintf('\n%s',repmat('-',1,72));
  fprintf(1,'\nCAT preprocessing takes %0.0f minute(s) and %0.0f second(s).\n', ...
    floor(round(etime(clock,res.stime))/60),mod(round(etime(clock,res.stime)),60));
  
  % image quality (just use real to avoid some rare problems with irrational values that should not occur anymore)
  cat_io_cprintf(color(QMC,real(qa.qualityratings.SIQR)), sprintf('Structural Image Quality Rating (SIQR):  %5.2f%%%% (%s)\n',...
    mark2rps(real(qa.qualityratings.IQR)),mark2grad(real(qa.qualityratings.IQR))));
  
  % processing quality
  % coming soon
  
  % print GMV / GMT values (for developer in colors)
  if job.extopts.expertgui > 0
    if job.extopts.expertgui > 1 && isfield(qa.subjectmeasures,'dist_thickness')
      col = colorgmt(GMC,qa.subjectmeasures.dist_thickness{1}(1)); 
    elseif job.extopts.expertgui > 1 
      col = colorgmt(GMC, qa.subjectmeasures.vol_rel_CGW(2) * 5 ); % simple translation to thickness
    else
      col = [0 0 0];
    end
    cat_io_cprintf(col, sprintf('Relative gray matter volume (GMV/TIV):   %5.2f%%%% (%5.2f / %5.2f ml)\n',...
      qa.subjectmeasures.vol_rel_CGW(2)*100 , qa.subjectmeasures.vol_abs_CGW(2) , qa.subjectmeasures.vol_TIV ));
    if isfield(qa.subjectmeasures,'dist_thickness')
      cat_io_cprintf(col, sprintf('Gray matter thickness (GMT):           %5.2f %s %4.2f mm\n',...
         qa.subjectmeasures.dist_thickness{1}(1), native2unicode(177, 'latin1'), qa.subjectmeasures.dist_thickness{1}(2) ) );
    end
  end
  
  
  
  %% print subfolders
  if job.extopts.subfolders
    fprintf('Segmentations are saved in %s\n',mrifolder);
    fprintf('Reports are saved in %s\n',reportfolder);
    if job.output.ROI
      fprintf('Labels are saved in %s\n',labelfolder);
    end
    if job.output.surface && exist('Psurf','var') && ~isempty(Psurf)
      fprintf('Surface measurements are saved in %s\n',surffolder);
    end
  end

  fprintf('%s\n\n',repmat('-',1,72));
  
  % finish diary entry of "../report/cmdln_*.txt"
  % read diary and add the command-line output to the *.xml and *.mat file
  diary off; 
  try %#ok<TRYNC>
    fid  =fopen(res.catlog);
    txt = fread(fid,200000,'uint8=>char');
    fclose(fid); 
    txt2 = textscan(txt,'%s','Delimiter',''); 
    cat_io_xml(fullfile(reportfolder,['cat_' nam '.xml']),struct(...
      'catlog',txt2),'write+'); % here we have to use the write+!
  end    
  
  spm_progress_bar('Clear');
  cat_progress_bar('Clear');

end
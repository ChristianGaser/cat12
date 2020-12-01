function str = cat_main_reportstr(job,res,qa) 
% ______________________________________________________________________
% 
% Prepare text output for CAT report. This function may heavily change
% due to parameter changes. However, this only effects the output of the
% CAT report. Called from cat_main. 
%
%   str = cat_main_reportstr(job,res,qa,cat_warnings)
%
%   str          .. cellstrings with 3 elements with two strings
%   str{1}       .. full width table 
%   str{2:3}     .. half width table 
%
%   job          .. SPM/CAT parameter structure
%   res          .. SPM preprocessing structure
%   qa           .. CAT quality and subject information structure
%
%   See also cat_main_reportfig and cat_main_reportcmd.
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id$

%#ok<*AGROW,*NOCOM>

  QMC   = cat_io_colormaps('marks+',17);
  color = @(QMC,m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
   
  %mark2str2 = @(mark,s,val) sprintf(sprintf('\\\\bf\\\\color[rgb]{%%0.2f %%0.2f %%0.2f}%s',s),color(QMC,mark),val);
  marks2str   = @(mark,str) sprintf('\\bf\\color[rgb]{%0.2f %0.2f %0.2f}%s',color(QMC,mark),str);
  mark2rps    = @(mark) min(100,max(0,105 - mark*10)) + isnan(mark).*mark;
  grades      = {'A+','A','A-','B+','B','B-','C+','C','C-','D+','D','D-','E+','E','E-','F'};
  mark2grad   = @(mark) grades{min(numel(grades),max(max(isnan(mark)*numel(grades),1),round((mark+2/3)*3-3)))};

  catdef = cat_get_defaults; 
  
% blue for none defaults
% orange for warning changes 
  npara  = '\color[rgb]{0  0   0}'; 
  cpara  = '\color[rgb]{0  0   1}'; 
  wpara  = '\color[rgb]{1  0.3 0}'; 
  
  
  % CAT GUI parameter:
  % --------------------------------------------------------------------
  str{1} = [];

  % use red output if a beta version is used
  catv = qa.software.revision_cat; isbeta = strfind(lower(catv),'beta'); 
  if ~isempty(isbeta), catv = [catv(1:isbeta-1) '\color[rgb]{0.8 0 0}' catv(isbeta:isbeta+3 ) '\color[rgb]{0.8 0 0}' catv(isbeta+4:end)]; end

% red version in case of old version?
  % 1 line: Matlab, SPM12, CAT12 version number and GUI and experimental mode 
  str{1} = [str{1} struct('name', 'Version: OS / Matlab / SPM12 / CAT12:','value',...
    sprintf('%s / %s / %s / %s (%s)',qa.software.system,qa.software.version_matlab,...
    qa.software.version_spm,qa.software.version_cat,catv))];
  % add CAT segmentation version if not current
  if ~isempty(qa.software.version_segment)
    str{1}(end).name = [str{1}(end).name(1:end-1) ' / seg:']; 
    str{1}(end).value = [str{1}(end).value ' / \color[rgb]{0 0.2 1}' qa.software.version_segment '']; 
  end 
  % write GUI mode
  if     job.extopts.expertgui==1, str{1}(end).value = [str{1}(end).value '\bf\color[rgb]{0 0.2 1}e']; 
  elseif job.extopts.expertgui==2, str{1}(end).value = [str{1}(end).value '\bf\color[rgb]{0 0.2 1}d'];
  end  
  % write experimental flag
  if str{1}(end).value(end)==')'
    if job.extopts.experimental, str{1}(end).value = [str{1}(end).value(1:end-1) '\bf\color[rgb]{0 0.2 1}x\color[rgb]{0 0 0})']; end  
  else
    if job.extopts.experimental, str{1}(end).value = [str{1}(end).value '\bf\color[rgb]{0 0.2 1}x']; end  
  end
  if job.extopts.ignoreErrors > 2, str{1}(end).value = [str{1}(end).value '    \bf\color[rgb]{0.8 0 0}Ignore Errors!']; end  
  
  
  % 2 lines: TPM, Template, Normalization method with voxel size
  if strcmp(res.tpm(1).fname,catdef.opts.tpm{1}), cp{1} = npara; else, cp{1} = cpara; end
  str{1} = [str{1} struct('name', 'Tissue Probability Map:','value',[cp{1} strrep(spm_str_manip(res.tpm(1).fname,'k40d'),'_','\_')])];
  if res.do_dartel
    if job.extopts.regstr==0 % Dartel
      if strcmp(job.extopts.darteltpm{1},catdef.extopts.darteltpm{1}), cp{1} = npara; else, cp{1} = cpara; end
      str{1} = [str{1} struct('name', 'Dartel Registration to: ',...
                        'value',[cp{1} strrep(spm_str_manip(job.extopts.darteltpm{1},'k40d'),'_','\_')])];
    elseif job.extopts.regstr==4 % Dartel
      if strcmp(job.extopts.shootingtpm{1},catdef.extopts.shootingtpm{1}), cp{1} = npara; else, cp{1} = cpara; end
      str{1} = [str{1} struct('name', 'Shooting Registration to: ',...
                        'value',[cp{1} strrep(spm_str_manip(job.extopts.shootingtpm{1},'k40d'),'_','\_')])];
    else
      if strcmp(job.extopts.shootingtpm{1},catdef.extopts.shootingtpm{1}), cp{1} = npara; else, cp{1} = cpara; end
      if job.extopts.expertgui==0
        str{1} = [str{1} struct('name','Optimized Shooting Registration to:',...
                          'value',strrep(spm_str_manip(job.extopts.shootingtpm{1},'k40d'),'_','\_'))];
      else
        str{1} = [str{1} struct('name', sprintf('Optimized Shooting Registration (regstr:%s) to :',sprintf('%g ',job.extopts.regstr)),...
                          'value',[cp{1} strrep(spm_str_manip(job.extopts.shootingtpm{1},'k40d'),'_','\_')])];
      end
    end
  end

  % 1 line 1: Affreg
  if strcmp(job.opts.affreg,catdef.opts.affreg), cp{1} = npara; else, cp{1} = cpara; end
  if isfield(job,'useprior') && ~isempty(job.useprior) && exist(char(job.useprior),'file')
    affstr = 'AVGprior';
  else
    affstr = job.opts.affreg; 
  end
  str{1} = [str{1} struct('name', 'affreg:','value',sprintf('%s{%s}',cp{1},affstr))];
  % 1 line 2: APP

  if job.extopts.APP == catdef.extopts.APP, cp{1} = npara; else, cp{1} = cpara; end
  APPstr = {'none','light','full','','','animal'}; APPstr{1071} = 'default'; APPstr{1145} = 'default2'; 
  str{1}(end).name  = [str{1}(end).name(1:end-1) ' / APP ']; 
  str{1}(end).value = [str{1}(end).value sprintf(' / %s{%s}',cp{1},APPstr{job.extopts.APP+1})];

% #####
% one new super SPM parameter?
% #####
  % 1 line 3: biasstr / biasreg+biasfwhm
  if job.opts.biasacc>0
    if job.opts.biasacc == catdef.opts.biasstr, cp{1} = npara; else, cp{1} = cpara; end % yes, catdef.opts.biasstr!
    biasacc = {'ultralight','light','medium','strong','heavy'};
    str{1}(end).name  = [str{1}(end).name(1:end-1) ' / biasstr '];  
    str{1}(end).value = [str{1}(end).value sprintf(' / %s{%s}',cp{1},biasacc{round(job.opts.biasacc*4)+1})];
    if job.extopts.expertgui % add the value
      str{1}(end).value = [str{1}(end).value sprintf('(%0.2f;breg:%0.0e;bfwhm:%0.0f)',job.opts.biasacc,job.opts.biasreg,job.opts.biasfwhm)]; 
    end
  elseif job.opts.biasstr>0
    if job.opts.biasstr == catdef.opts.biasstr, cp{1} = npara; else, cp{1} = cpara; end
    biasstr = {'ultralight','light','medium','strong','heavy'};
    str{1}(end).name  = [str{1}(end).name(1:end-1) ' / biasstr '];  
    str{1}(end).value = [str{1}(end).value sprintf(' / %s{%s}',cp{1},biasstr{round(job.opts.biasstr*4)+1})];
    if job.extopts.expertgui % add the value
      str{1}(end).value = [str{1}(end).value sprintf('(%0.2f;breg:%0.0e;bfwhm:%0.0f)',job.opts.biasstr,job.opts.biasreg,job.opts.biasfwhm)]; 
    end
  else
    if job.opts.biasreg  == catdef.opts.biasreg,  cp{1} = npara; else, cp{1} = cpara; end
    if job.opts.biasfwhm == catdef.opts.biasfwhm, cp{2} = npara; else, cp{2} = cpara; end
    str{1}(end).name  = [str{1}(end).name(1:end-1) ' / biasreg / biasfwhm'];
    str{1}(end).value = [str{1}(end).value sprintf(' / %s{%0.0e} / %s{%0.2f}',cp{1},job.opts.biasreg,cp{2},job.opts.biasfwhm)]; 
  end
  % 1 line 3: SPM segmentation accuracy with samp and tol
  if isfield(job.opts,'acc') && job.opts.acc>0
    %str{1} = [str{1} struct('name', '','value','')];
    if job.opts.acc == catdef.opts.acc, cp{3} = npara; else, cp{3} = cpara; end
    if job.extopts.expertgui==0
      accstr = {'ultra low','low','std','high','ultra high'};
      str{1}(end).name  = [str{1}(end).name(1:end-1) ' / accuracy '];  
      str{1}(end).value = [str{1}(end).value sprintf(' / %s{%s}',co,accstr{round(job.opts.acc*4)+1})];
    else % add the value
      if job.opts.samp == catdef.opts.samp, cp{1} = npara; else, cp{1} = cpara; end
      if job.opts.tol  == catdef.opts.tol,  cp{2} = npara; else, cp{2} = cpara; end
      str{1}(end).name  = [str{1}(end).name(1:end-1) ' / acc (samp/tol) '];  
      str{1}(end).value = [str{1}(end).value sprintf('%s|%0.2f} (%s{%0.2f}/%s{%0.0e})',cp{3},job.opts.acc,cp{1},job.opts.samp,cp{2},job.opts.tol)]; 
    end
  else
    if job.extopts.expertgui
      %str{1} = [str{1} struct('name', '','value','')];
      if job.opts.samp == catdef.opts.samp, cp{1} = npara; else, cp{1} = cpara; end
      if job.opts.tol  == catdef.opts.tol,  cp{2} = npara; else, cp{2} = cpara; end
      str{1}(end).name  = [str{1}(end).name(1:end-1) '/ samp / tol '];
      str{1}(end).value = [str{1}(end).value sprintf('/ %s{%0.2f} / %s{%0.0e}',cp{1},job.opts.samp,cp{2},job.opts.tol)]; 
    end
  end


  % 1 line: adaptive noise parameter ( MRFstr + SANLM + NCstr )
  NCstr.labels = {'none','full','light','medium','strong','heavy'};
  NCstr.values = {0 1 2 -inf 4 5}; 
  defstr  = {'none','ultralight','light','medium','strong','heavy',... sanlm vs. isarnlm
             'ultralight+','ultralight+','light+','medium+','strong+','heavy+'};
  defstrm = @(x) defstr{ round(max(0,min(2,x))*4) + 1 + (x>0) + (x>1)};
  str{1} = [str{1} struct('name', 'Noise reduction:','value','')]; 
  if job.extopts.NCstr == catdef.extopts.NCstr, cp{1} = npara; else, cp{1} = cpara; end 
  if job.extopts.NCstr==0 
    if job.extopts.mrf==0
      str{1}(end).value = [cpara 'no noise correction'];
    else
      if job.extopts.expertgui==0
        str{1}(end).value = [cpara 'MRF']; 
      else
        str{1}(end).value = sprintf('%sMRF(%0.2f)',cpara,job.extopts.mrf); 
      end  
    end
  else
    str{1}(end).value = sprintf('SANLM(%s{%s})',cp{1},NCstr.labels{find(cell2mat(NCstr.values)==job.extopts.NCstr,1,'first')});
  end

  if job.extopts.NCstr~=0 && job.extopts.mrf
    if job.extopts.expertgui==0
      str{1}(end).value = '+MRF'; 
    else
      str{1}(end).value = sprintf('%s{+MRF(%0.2f)}',cpara,job.extopts.mrf); 
    end 
  end


  % 1 line(s): LASstr / GCUTstr / CLEANUPstr
  if job.extopts.LASstr     == catdef.extopts.LASstr,     cp{1} = npara; else, cp{1} = cpara; end 
  if job.extopts.gcutstr    == catdef.extopts.gcutstr,    cp{2} = npara; else, cp{2} = cpara; end 
  if job.extopts.cleanupstr == catdef.extopts.cleanupstr, cp{3} = npara; else, cp{3} = cpara; end 
  gcutstr  = {'none','SPM','GCUT','APRG','APRG2'};  
  if ~job.extopts.expertgui
    str{1}(end).name  = 'LAS strength / Skull-Stripping:';
    str{1}(end).value = sprintf('%s{%s} / %s{%s}',...
      cp{1},defstrm(job.extopts.LASstr),...
      cp{2},gcutstr{ceil(job.extopts.gcutstr+2)}); 
  else
    str{1}(end).name  = 'LASstr / GCUTstr / CLEANUPstr:';
    str{1}(end).value = sprintf('%s{%s(%0.2f)} / %s{%s(%0.2f)} / %s{%s(%0.2f)}',...
      cp{1},defstrm(job.extopts.LASstr),job.extopts.LASstr,...
      cp{2},gcutstr{ceil(job.extopts.gcutstr+2)},job.extopts.gcutstr,...
      cp{3},defstrm(job.extopts.cleanupstr),job.extopts.cleanupstr); 
  end
  if job.extopts.spm_kamap == catdef.extopts.spm_kamap,  cp{1} = npara; else, cp{1} = cpara; end 
  if job.extopts.WMHC      == catdef.extopts.WMHC,       cp{2} = npara; else, cp{2} = cpara; end 
  if job.extopts.SLC       == catdef.extopts.SLC,        cp{3} = npara; else, cp{3} = cpara; end 
  if job.extopts.collcorr  == catdef.extopts.collcorr,   cp{4} = npara; else, cp{4} = cpara; end 
  restype = char(fieldnames(job.extopts.restypes));
  if strcmp(restype, catdef.extopts.restype), cp{5} = npara; else, cp{5} = cpara; end 
% ############
% WMHC in case of SPM segmentation is SPM
% ############
  if job.extopts.expertgui
    str{1} = [str{1} struct('name', 'KAMAP / WMHC / SLC / collcorr / restype:','value',...
           sprintf('%s{%d} / %s{%d} / %s{%d} / %s{%d} / %s{%s}',...
          cp{1},job.extopts.spm_kamap, cp{2},job.extopts.WMHC, cp{3},job.extopts.SLC, cp{4},job.extopts.collcorr, cp{5},restype))];
  else
    kamapstr  = {'SPM US','','KAMAP'};  
    wmhcstr   = {'none (WMH=GM)','temporary (WMH=GM)','(WMH=WM)','own class'};
    str{1} = [str{1} struct('name', 'Initial Segmentation / WMH Correction / Int. Res.:',...
      'value',sprintf('%s{%s} / %s{%s} / %s{%s}',cp{1},kamapstr{job.extopts.spm_kamap+1}, cp{2},wmhcstr{job.extopts.WMHC+1}, cp{5},restype))];
  end
  if ~strcmp('native',restype) && ~strcmp('best',restype)
    if any(job.extopts.restypes.(restype) == catdef.extopts.resval),  cp{1} = npara; else, cp{1} = cpara; end 
    str{1}(end).value = [str{1}(end).value sprintf('(%s{%0.2f %0.2f})',cp{1}, job.extopts.restypes.(restype))];
  end 

  % line 8: resolution parameter
  % RD202007:  I am sure if a colorrating works here - it is just to much.
  v3 = sprintf('%s',native2unicode(179, 'latin1')); 
  if job.output.surface
    str{1} = [str{1} struct('name', 'Voxel resolution (original > internal > PBT; vox):',...
           'value',sprintf('%4.2fx%4.2fx%4.2f > %4.2fx%4.2fx%4.2f > %4.2f%s mm%s; %4.2f%s mm%s ', ... 
           qa.qualitymeasures.res_vx_vol,qa.qualitymeasures.res_vx_voli,job.extopts.pbtres,...
           v3,v3,job.extopts.vox(1),v3,v3))];
  else
    str{1} = [str{1} struct('name', 'Voxel resolution (original > intern; vox):',...
           'value',sprintf('%4.2fx%4.2fx%4.2f mm%s > %4.2fx%4.2fx%4.2f mm%s; %4.2f%s mm%s', ...
           qa.qualitymeasures.res_vx_vol,v3,qa.qualitymeasures.res_vx_voli,...
           v3,job.extopts.vox(1),v3,v3))];
  end       
  % str{1} = [str{1} struct('name', 'Norm. voxel size:','value',sprintf('%0.2f mm',job.extopts.vox))]; % does not work yet 

  % line 9: surface parameter
  
  
  % Image Quality measures:
  % --------------------------------------------------------------------
  str{2} =       struct('name', '\bfImage and Preprocessing Quality:','value',''); 
  str{2} = [str{2} struct('name',' Resolution:','value',marks2str(qa.qualityratings.res_RMS,...
    sprintf('%5.2f%% (%s)',mark2rps(qa.qualityratings.res_RMS),mark2grad(qa.qualityratings.res_RMS))))];
  str{2} = [str{2} struct('name',' Noise:','value',marks2str(qa.qualityratings.NCR,...
    sprintf('%5.2f%% (%s)',mark2rps(qa.qualityratings.NCR),mark2grad(qa.qualityratings.NCR))))];
  str{2} = [str{2} struct('name',' Bias:','value',marks2str(qa.qualityratings.ICR,...
    sprintf('%5.2f%% (%s)',mark2rps(qa.qualityratings.ICR),mark2grad(qa.qualityratings.ICR))))]; % not important and more confusing 
  str{2} = [str{2} struct('name','\bf Weighted average (IQR):','value',marks2str(qa.qualityratings.IQR,...
    sprintf('%5.2f%% (%s)',mark2rps(qa.qualityratings.IQR),mark2grad(qa.qualityratings.IQR))))];
  if isfield(qa.qualitymeasures,'SurfaceEulerNumber') && ~isempty(qa.qualitymeasures.SurfaceEulerNumber)
    if job.extopts.expertgui
      str{2} = [str{2} struct('name',' Mean surface Euler number:','value',marks2str(qa.qualityratings.SurfaceEulerNumber,...
                sprintf('%g', qa.qualitymeasures.SurfaceEulerNumber)))]; 
              
      if isfield(qa.qualitymeasures,'SurfaceDefectNumber') && ~isempty(qa.qualitymeasures.SurfaceDefectNumber) 
        str{2}(end).name  = [str{2}(end).name(1:end-8)  ' / defect number:'];
        str{2}(end).value = [str{2}(end).value ' / ' marks2str(qa.qualityratings.SurfaceDefectNumber,...
                sprintf('%0.2f', qa.qualitymeasures.SurfaceDefectNumber)) ];
      end
            
    else
      str{2} = [str{2} struct('name',' Mean surface Euler number:','value',sprintf('%g', qa.qualitymeasures.SurfaceEulerNumber))]; 
    end
  end
  
  if isfield(qa.qualitymeasures,'SurfaceDefectArea') && ~isempty(qa.qualitymeasures.SurfaceDefectArea) && ~any(job.output.surface == [5 6])
    if job.extopts.expertgui
      str{2} = [str{2} struct('name',' Mean topology defects size:','value',marks2str(qa.qualityratings.SurfaceDefectArea,...
                sprintf('%0.2f%%', qa.qualitymeasures.SurfaceDefectArea)))];

      if isfield(qa.qualitymeasures,'SurfaceSelfIntersections') && ~isempty(qa.qualitymeasures.SurfaceSelfIntersections)
        str{2}(end).name  = [str{2}(end).name(1:end-6)  ' / self-inters. size:'];
        str{2}(end).value = [str{2}(end).value ' / ' marks2str(qa.qualityratings.SurfaceSelfIntersections,...
                sprintf('%0.2f%%', qa.qualitymeasures.SurfaceSelfIntersections)) ];
      end
   
    else
      str{2} = [str{2} struct('name',' Mean size of topology defects:','value',sprintf('%0.2f%%', qa.qualitymeasures.SurfaceDefectArea))];
    end
    
  end
  
  if job.extopts.expertgui && isfield(qa.qualityratings,'SurfaceIntensityRMSE')
      str{2} = [str{2} struct('name',' Surface intensity / position RMSE:','value',[ marks2str( qa.qualityratings.SurfaceIntensityRMSE ,...
        sprintf('%0.2f', qa.qualitymeasures.SurfaceIntensityRMSE)) ' / ' ...
        marks2str( qa.qualityratings.SurfacePositionRMSE ,sprintf('%0.2f', qa.qualitymeasures.SurfacePositionRMSE) ) ] ) ];
  end

  % Subject Measures
  % --------------------------------------------------------------------
  % Volume measures

  % header
  str{3} = struct('name', '\bfVolumes:','value',sprintf('%5s %5s %5s ','CSF','GM','WM')); 
  if job.extopts.WMHC>2, str{3}(end).value = [str{3}(end).value sprintf('%5s ','WMH')]; end
  if job.extopts.SLC>0,  str{3}(end).value = [str{3}(end).value sprintf('%5s ','SL')];  end

  % absolute volumes
  if job.extopts.WMHC<=2 && isfield(qa,'subjectmeasure') && isfield(qa.subjectmeasures,'vol_rel_WMH') && ...
    ( (qa.subjectmeasures.vol_rel_WMH>0.01 || qa.subjectmeasures.vol_rel_WMH/qa.subjectmeasures.vol_rel_CGW(3)>0.02) )
    if job.extopts.WMHC == 2
      str{3} = [str{3} struct('name', ' Absolute volume:','value',...
        sprintf('%5.0f %5.0f {\\bf\\color[rgb]{1 0 1}%5.0f} ', qa.subjectmeasures.vol_abs_CGW(1:3)))];
    else
      str{3} = [str{3} struct('name', ' Absolute volume:','value',...
        sprintf('{%5.0f \\bf\\color[rgb]{1 0 1}%5.0f} %5.0f', qa.subjectmeasures.vol_abs_CGW(1:3)))];
    end      
  else
    str{3} = [str{3} struct('name', ' Absolute volume:','value',...
      sprintf('%5.0f %5.0f %5.0f ', qa.subjectmeasures.vol_abs_CGW(1:3)))];
  end
  if job.extopts.WMHC>2,  str{3}(end).value = [str{3}(end).value sprintf('%5.1f ',qa.subjectmeasures.vol_abs_CGW(4))]; end 
  if job.extopts.SLC>0,   str{3}(end).value = [str{3}(end).value sprintf('%5.1f ',qa.subjectmeasures.vol_abs_CGW(5))]; end
  str{3}(end).value = [str{3}(end).value 'cm' native2unicode(179, 'latin1')];

  % relative volumes
  if job.extopts.WMHC<=2 && isfield(qa,'subjectmeasure') && isfield(qa.subjectmeasures,'vol_rel_WMH') && ...
    ( (qa.subjectmeasures.vol_rel_WMH>0.01 || qa.subjectmeasures.vol_rel_WMH/qa.subjectmeasures.vol_rel_CGW(3)>0.02) )
    if job.extopts.WMHC == 2
      str{3} = [str{3} struct('name', ' Relative volume:','value',...
        sprintf('%5.1f %5.1f {\\bf\\color[rgb]{1 0 1}%5.1f} ', qa.subjectmeasures.vol_rel_CGW(1:3)*100))];
    else
      str{3} = [str{3} struct('name', ' Relative volume:','value',...
        sprintf('{%5.1f \\bf\\color[rgb]{1 0 1}%5.1f} %5.1f ', qa.subjectmeasures.vol_rel_CGW(1:3)*100))];
    end
  else
    str{3} = [str{3} struct('name', ' Relative volume:','value',...
      sprintf('%5.1f %5.1f %5.1f ', qa.subjectmeasures.vol_rel_CGW(1:3)*100))];
  end
  if job.extopts.WMHC>2,  str{3}(end).value = [str{3}(end).value sprintf('%5.1f ',qa.subjectmeasures.vol_rel_CGW(4)*100)]; end 
  if job.extopts.SLC>0,   str{3}(end).value = [str{3}(end).value sprintf('%5.1f ',qa.subjectmeasures.vol_rel_CGW(5)*100)]; end
  str{3}(end).value = [str{3}(end).value '%'];

  % warning if many WMH were found and not handled as extra class 
  if job.extopts.WMHC<=2 && isfield(qa,'subjectmeasure') && isfield(qa.subjectmeasures,'vol_rel_WMH') && ...
    ( (qa.subjectmeasures.vol_rel_WMH>0.01 || qa.subjectmeasures.vol_rel_WMH/qa.subjectmeasures.vol_rel_CGW(3)>0.02) )
    if job.extopts.WMHC == 2
      str{3}(end-1).value = [str{3}(end-1).value sprintf('\\color[rgb]{1 0 1} (WM inc. %0.0fcm%s WMHs)', ...
        qa.subjectmeasures.vol_abs_WMH,native2unicode(179, 'latin1'))];  
      str{3}(end).value   = [str{3}(end).value   sprintf('\\color[rgb]{1 0 1} (WM inc. %0.1f%% WMHs)', qa.subjectmeasures.vol_rel_WMH)];  
      %str{3}(end).value = [str{3}(end).value sprintf('\\bf\\color[rgb]{1 0 1} WMHs %0.1f%% > WM!', qa.subjectmeasures.vol_rel_WMH * 100)];  
    else
      str{3}(end-1).value = [str{3}(end-1).value sprintf('\\color[rgb]{1 0 1} (GM inc. %0.0fcm%s WMHs)', ...
        qa.subjectmeasures.vol_abs_WMH,native2unicode(179, 'latin1'))];   
      str{3}(end).value   = [str{3}(end).value sprintf('\\color[rgb]{1 0 1} (GM inc. %0.1f%% WMHs)', qa.subjectmeasures.vol_rel_WMH * 100)];  
      %str{3}(end).value = [str{3}(end).value sprintf('\\bf\\color[rgb]{1 0 1} WMHs %0.1f%% > GM!', qa.subjectmeasures.vol_rel_WMH * 100)];  
    end
  end
  
  str{3} = [str{3} struct('name', ' TIV:','value', sprintf(['%0.0f cm' native2unicode(179, 'latin1')],qa.subjectmeasures.vol_TIV))];  
  if isfield(qa.subjectmeasures,'surf_TSA') && job.extopts.expertgui>1
    str{3}(end).name  = [str{3}(end).name  ' / TSA:']; 
    str{3}(end).value = [str{3}(end).value sprintf(' / %0.0f cm%s' ,qa.subjectmeasures.surf_TSA,char(178))];  
  end

  % Surface measures - Thickness, (Curvature, Depth, ...)
  %if cellfun('isempty',strfind({Psurf(:).Pcentral},'ch.')), thstr = 'Cerebral Thickness'; else thstr = 'Thickness'; end
  thstr = 'Thickness';
  if isfield(qa.subjectmeasures,'dist_thickness') && ~isempty(qa.subjectmeasures.dist_thickness)
    str{3} = [str{3} struct('name', ['\bf' thstr ':'],'value',sprintf('%5.2f%s%5.2f mm', ...
           qa.subjectmeasures.dist_thickness{1}(1),native2unicode(177, 'latin1'),qa.subjectmeasures.dist_thickness{1}(2)))];
         
    % we warn only if WMHC is off ... without WMHC you have not thresholds!
    %if job.extopts.WMHC==0 && (qa.subjectmeasures.vol_rel_WMH>0.01*3 || ... % 3 times higher treshold 
    %     qa.subjectmeasures.vol_rel_WMH/qa.subjectmeasures.vol_rel_CGW(3)>0.02*3)
    %  str{3}(end).value   = [str{3}(end).value   sprintf('\\color[rgb]{1 0 1} (may biased by WMHs!)')]; 
    %end
         
    if isfield(qa.subjectmeasures,'dist_gyruswidth') && ~isnan(qa.subjectmeasures.dist_gyruswidth{1}(1))
      str{3} = [str{3} struct('name', '\bfGyruswidth:','value',sprintf('%5.2f%s%5.2f mm', ...
             qa.subjectmeasures.dist_gyruswidth{1}(1),native2unicode(177, 'latin1'),qa.subjectmeasures.dist_gyruswidth{1}(2)))];
    end
    if isfield(qa.subjectmeasures,'dist_sulcuswidth') && ~isnan(qa.subjectmeasures.dist_sulcuswidth{1}(1))
      str{3} = [str{3} struct('name', '\bfSulcuswidth:','value',sprintf('%5.2f%s%5.2f mm', ...
             qa.subjectmeasures.dist_sulcuswidth{1}(1),native2unicode(177, 'latin1'),qa.subjectmeasures.dist_sulcuswidth{1}(2)))];
    end
  end

  % Preprocessing Time
  str{2} = [str{2} struct('name','\bfProcessing time:','value',sprintf('%02.0f:%02.0f min', ...
    floor(round(etime(clock,res.stime))/60),mod(round(etime(clock,res.stime)),60)))]; 

%% Warnings
% key changes? 
%  - use backup pipeline (with #?)
%  - use AMAP/SPM/mixed ?
%  - add skull-stripping, background, ... settings? 
	str{2} = [str{2} struct('name', '','value','')]; % empty line
  wname  = {'note','warning','alert'}; 
  wcolor = [0 0 1; 1 0.5 0; 0.8 0 0]; 
  valnl  = 80; 
  for wi=2:-1:0
    warn = cat_io_addwarning(wi); 
    wn   = numel(warn);
    if wn
      if job.extopts.expertgui > -wi
        msg  =  strrep(warn(1).identifier,'_','\_') ; 
        for wmi=2:wn
          msg = [msg ', ' strrep( warn(wmi).identifier,'_','\_') ];
% linebreak may cause other problems ...          
          if numel(msg)>valnl, msg = [msg '\\n']; valnl = valnl + 80; end 
        end
        
        if wn~=1, wnamepl='s'; else, wnamepl=''; end
        str{2} = [str{2} struct('name', ...
          sprintf('\\color[rgb]{%0.1f %0.1f %0.1f}\\bf{%d %s%s} (%s)',...
            wcolor(wi+1,:), wn,wname{wi+1},wnamepl, msg),'value','')]; 
      end
    end
  end  


  %% adding one space for correct printing of bold fonts
  for ssi = 1:numel(str)
    for si = 1:numel(str{ssi})
      str{ssi}(si).name   = [str{ssi}(si).name  '  '];   
      str{ssi}(si).value  = [str{ssi}(si).value '  '];
    end
  end

end
function cstr = isdefault(job,field,str)
  
  eval(['catval = job.' field ';']); 
  catdef = cat_get_defaults(field);
  if strcmp(catval,catdef)
    cstr = str;
  else
    cstr = ['\color[rgb]{0 0 1}{' str '}']; 
  end
end
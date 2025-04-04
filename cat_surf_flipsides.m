function Prdata = cat_surf_flipsides(job)
% ______________________________________________________________________
% Mirror surface side. 
%
% There are to major ways (i) using the mapping between the templates, 
% and (ii) estimating an individual mapping. The first way only need a
% resampling of the data, whereas the second way requires a spherical 
% mapping (~5 minutes per hemisphere). 
%
% To keep things easy the first way used the resampled gifti files 
% s#mm.[rl]h.TEXTURE.resampled.FILENAME.gii as input and ...
% 
%   * template | subject
%
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$



  
% Todo:
% 
% * Hier muss ein allgemeinenes Konzept vorliegen, mit dem sich sowohl eine
%   automatische (z.B. fuer cat_surf_calc), als auch manuelle Nutzung (GUI)
%   moeglich ist. 
%
% * Ich fuerchte ich mach die dinge komplizierter als sie sind ...
%
% Facts: 
%   1. Die Daten muessen IMMER gesampled werden.
%   2. Fuer eine Analyse kommen nur GIFTIs in Frage. 
%   3. Fuer eine Analyse spielt die Seite keine Rolle. 
%   >> Als in/output waere die s#mm.rh.DATA.resampled.SUBJECT.gii wohl sinnvoll  
%   >> Als Name koennte man entweder eine dritte seite Lh/Rh nutzen (=mesh) oder 
%      einen weiteren Term 'flip' oder 'flipped' einfuegen 
%   >> Die 32k Meshes sind flip-ready!
%   >> Nutzt man den Subject-case waere die Orignaldaten rh.DATA.SUBJECT 
%      wohl besser (achtung hier koennten auch resampled reinrutschen)
%   >> Analysefertig ausgabe sinnvoll ... andere Faelle machens nur unnoetig komplex
%
% Frage: 
%   1. Spielt die Reihefolge des Remeshings eine Rolle?
%   >> Wahrscheinlich besser die Originaloberflaeche auf die andere Seite zu mappen im individuellen Fall.
%   >> Beim template-Fall spiel das wohl eher keine Rolle
%
% * Namensgebebung: 
%      rh.thickness.MaMu99.gii > Lh.thickness.MaMu99.gii
%      s#mm.rh.thickness.resampled.MaMu99.gii
%
% * Flip man nur texturen oder ganze Oberflaechen?
%   - wegen Analyse ganze Oberflaechen
%
% * GUI Anbindung
%   - job.cdata (lh > rh, rh > lh) ... am besten nur eine seite waehlbar, so bekommst du ein sichereres reslutat 
%   - job.type (template | subject) ...
%
% * automatische Ansbindung fuer surf_calc
%   - es wird kein bild geschrieben, sondern die daten werden als output
%     uebergeben

% Ablauf...
% * Input sollten die Originaldaten (Texturen) sein.
% Subject (aufwendiger - genauer)
% * Die Daten werden geflipt und ein Mapping zur anderen fsavg bestimmt.
% * Resampling
% Template (einfacher - schneller)
% * Die Daten werden gesampled 
% * Das vorbestimmte Mapping zwischen den fsavgs wird angewandt.


  if ~exist('job','var'), job=struct(); end

  def.verb           = 1; 
  def.cdata          = {}; 
  def.type           = 'template'; % template | subject
  def.recalctemplate = 0; 
  def.usetemplate    = 1; 
  def.debug          = 0; 
  job = cat_io_checkinopt(job,def);
  
  if isempty(job.cdata); job.cdata = cellstr(spm_select(inf,'gifti','Select surface','','','^s.*')); end
  job.cdata = cellstr(job.cdata);
  if isempty(job.cdata)
    fprintf('Nothing to do.\n'); 
    return;
  end

  
  %% fsavg sphere and central
  Psphere{1}   = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces','lh.sphere.freesurfer.gii');
  Psphere{2}   = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces','rh.sphere.freesurfer.gii');
  Pcentral{1}  = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces','lh.central.freesurfer.gii');
  Pcentral{2}  = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces','rh.central.freesurfer.gii');
  % fsavg flipping spherical mapping
  Preg{1}      = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces','lh.sphere.rhreg.freesurfer.gii');
  Preg{2}      = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces','rh.sphere.lhreg.freesurfer.gii');
  % fsavg flipping 
  Pflipsphere{1}  = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces','rh.sphere.lhflip.freesurfer.gii');
  Pflipsphere{2}  = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces','lh.sphere.rhflip.freesurfer.gii');
  Pflipcentral{1} = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces','rh.central.lhflip.freesurfer.gii');
  Pflipcentral{2} = fullfile(fileparts(mfilename('fullpath')),'templates_surfaces','lh.central.rhflip.freesurfer.gii');

  
  
  %% spherical registration of lh to rh template 
  if job.recalctemplate
    for si=1:2
      %% flipping
      flipsurf(Pcentral{si},Pflipcentral{si});
      flipsurf(Psphere{si},Pflipsphere{si});    
  
      % registration 
      cmd = sprintf('CAT_WarpSurf -type 0 -i "%s" -is "%s" -t "%s" -ts "%s" -ws "%s"',...
        Pflipcentral{si},Pflipsphere{si},Pcentral{~(si-1)+1},Psphere{~(si-1)+1},Preg{~(si-1)+1});
      cat_system(cmd,job.debug);
    end
  end
  
  
  %%
  sides  = {'lh','rh'};
  %sidesi = {'flh','frh'}; % further letters ... problem [rl]h.* 
  %sidesi = {'ih','ih'};   % other letter ... problem [rl]h.* > [rli]h, but ih for both left and right 
  sidesi = {'rh','lh'};    % simple flipped ... flip information in resampled field  
  Prdata = cell(numel(job.cdata),1);
  Prmesh = cell(numel(job.cdata),1);
  if job.usetemplate
    % this is the simple case 
    
    sinfo = cat_surf_info(job.cdata);
     
    for di = 1:numel(job.cdata)
      %%
      try
        switch sinfo(di).side
          case 'lh', si = 1; 
          case 'rh', si = 2; 
        end

        Prdata(di) = cat_surf_rename(sinfo(di).Pdata,'side',sidesi{si},...
          'dataname',sinfo(di).dataname,'templateresampled',sprintf('%sresampled',sides{si}));
        Prmesh(di) = cat_surf_rename(sinfo(di).Pmesh,'side',sidesi{si},...
          'dataname','central','templateresampled',sprintf('%sresampled',sides{si}));

        % flip gifti
        gii = gifti(sinfo(di).Pmesh); 
        [pp,ff] = spm_fileparts(sinfo(di).Pmesh);
        gii.vertices(:,1) = -gii.vertices(:,1);
        gii.faces = [gii.faces(:,2),gii.faces(:,1),gii.faces(:,3)];
        cat_io_FreeSurfer('write_surf',fullfile(pp,[ff 'tmp']),gii); 

        cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s" "%s" "%s"',...
          fullfile(pp,[ff 'tmp']),Preg{~(si-1)+1},Psphere{si},Prmesh{di},job.cdata{di},Prdata{di});
        err = cat_system(cmd,job.debug); 

        if err
          cat_io_cprintf('err','Case "%s" did not work.\n',job.cdata{di}); 
          continue;
        end

        gc = gifti(Prdata{di}); gm = gifti(Prmesh{di}); gm.cdata = gc.cdata; save(gm,Prdata{di});    

        % delete temporary files
        if exist(Prmesh{di},'file'), delete(Prmesh{di}); end
        if exist(fullfile(pp,[ff 'tmp']),'file'), delete(fullfile(pp,[ff 'tmp'])); end

        if job.verb
          fprintf('Display %s\n',spm_file(Prdata{di},'link','cat_surf_display(''%s'')'));
        end
      catch
        cat_io_cprintf('err','Case "%s" did not work.\n',job.cdata{di}); 
      end
    end
  else
    error('not ready now');
    
    %{
    for di = 1:numel(job.cdata)
      sinfo = cat_surf_info(job.cdata{di});

      switch sinfo.side
        case 'lh', si = 1; 
        case 'rh', si = 2; 
      end

      % flipping
      flipsurf(Pcentral{si},Pflipcentral{si});
      flipsurf(Psphere{si},Pflipsphere{si});    
  
      % registration 
      cmd = sprintf('CAT_WarpSurf -type 0 -i "%s" -is "%s" -t "%s" -ts "%s" -ws "%s"',...
        Pflipcentral{si},Pflipsphere{si},Pcentral{~(si-1)+1},Psphere{~(si-1)+1},Preg{~(si-1)+1});
      cat_system(cmd,opt.debug);
    
      end
    %}
    
  end
  
end
function flipsurf(P,Po)
  gii  = gifti(P);
  gii.vertices(:,1) = -gii.vertices(:,1);
  gii.faces = [gii.faces(:,2),gii.faces(:,1),gii.faces(:,3)];
  save(gii,Po);
end
  
  
  
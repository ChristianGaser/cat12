function varargout = cat_io_dcm2bids_reportfig(action,varargin)
%cat_io_dcm2bids_reportfig. Report figure of cat_io_dcm2bids.
%  Diagrams of the report tables (see cat_io_dcm2bids_report) for the PI
%  and project level as A4 pages: an overview page of the BIDS directory
%  and one page per site (for several sites) with
%    header        project, BIDS directory, site (scanner), number of
%                  subjects, sessions, and scans, and the key ratings
%                  (median quality SQR, protocol fit PFR, and orientation
%                  SOR of the sessions, and the conform scans)
%    sex and age   sex of the subjects (pie) and age pyramid (age at the
%                  first session of each subject, males left, females
%                  right, see opts.reportfig.agerange and agestep)
%    recruitment   sessions per month (from the date of the session label,
%                  i.e. not for coded session labels of anonymize > 1)
%    sites         site table (overview page of several sites)
%    completeness  protocol status of the best scan of each typical scan
%                  type (opts.scantypes) and session, or missing scans
%    ratings       session ratings (SQR, PFR, SOR, BBR, PRR) and image
%                  quality SQR of the best scan of each scan type and
%                  session as box plots (marks 0.5 to 6.5, worse ratings
%                  are clipped as in the CAT report)
%    volumes       relative tissue volumes of the best T1w scan by age with
%                  the robust age/sex model of the sample
%    homogeneity   sample homogeneity of the mwc1 maps (quartic mean
%                  Z-score) by the T1w quality
%    check         sessions with the worst ratings (opts.reportfig.nworst)
%  The values are given in or next to the diagrams. The pages are saved as
%    <outdir>/<subdir>/<BIDSdir>-report/report_<BIDSsubdir>.pdf
%  and the overview is shown in the SPM Graphics window in interactive
%  sessions. The statistics of the volumes and the homogeneity are columns
%  of the session table (see cat_io_dcm2bids_stats).
%
%  varargout = cat_io_dcm2bids_reportfig(action,varargin)
%
%  Actions (see the help of the local functions):
%    Ppdf = cat_io_dcm2bids_reportfig('writeReportFigure',Poutdir,job)
%
%  See also cat_io_dcm2bids_report, cat_io_dcm2bids_stats, cat_plot_boxplot.

  switch action
    case {'writeReportFigure'}
      [varargout{1:nargout}] = feval(action,varargin{:});
    otherwise
      error('cat_io_dcm2bids_reportfig:unknownAction','Unknown action "%s".',action);
  end
end
% =========================================================================
function Ppdf = writeReportFigure(Poutdir,job)
%writeReportFigure. Report figure of the report tables of the BIDS
%  directory of the job (see the help of this file).
  Prep = fullfile(Poutdir,[job.BIDSdir '-report']);
  Ptab = fullfile(Prep,'tables',job.BIDSsubdir);
  Ppdf = fullfile(Prep,['report_' job.BIDSsubdir '.pdf']);
  if exist(Ppdf,'file'), delete(Ppdf); end
  fmt  = job.opts.tableformat;
  S    = readTable(fullfile(Ptab,['scans.' fmt]));
  E    = readTable(fullfile(Ptab,['sessions.' fmt]));
  Q    = readTable(fullfile(Ptab,['project.' fmt]));
  if isempty(E), Ppdf = ''; return; end
  F    = setupFigure(job,E,Q);

  % overview page and one page per site (for several sites)
  pages = {''};
  if numel(F.sites) > 1, pages = [pages F.sites]; end
  fh = figure('Visible','off','Color','w','Units','centimeters','Position',[1 1 F.pw F.ph], ...
    'PaperUnits','centimeters','PaperSize',[F.pw F.ph],'PaperPosition',[0 0 F.pw F.ph], ...
    'MenuBar','none','ToolBar','none');
  fhclean = onCleanup(@() close(fh));
  ws = warning('off','MATLAB:print:ContentTypeImageSuggested'); % vector output of many objects
  for pgi = 1:numel(pages)
    clf(fh);
    if isempty(pages{pgi})
      si = true(numel(S),1); ei = true(numel(E),1);
    else
      si = strcmp(textCol(S,'site'),pages{pgi}); ei = strcmp(textCol(E,'site'),pages{pgi});
    end
    drawPage(fh,F,S(si),E(ei),pages{pgi},pgi,numel(pages));
    exportgraphics(fh,Ppdf,'ContentType','vector','Padding','figure','Append',pgi > 1);
    if pgi == 1, showSPM(fh); end
  end
  warning(ws);
  fprintf('  Report figure: %s\n',Ppdf);
end
% =========================================================================
function F = setupFigure(job,E,Q)
%setupFigure. Common settings of all pages: page size, fonts, colors, scan
%  types, sites, and the statistics of the whole sample (volume model and
%  homogeneity threshold), i.e., the site pages show the same model.
  F.pw = 21; F.ph = 29.7; F.fs = 7;                 % A4 portrait in cm and base font size
  F.project = job.subdir; F.dataset = job.BIDSsubdir;
  F.ropts   = job.opts.reportfig;
  F.marks   = [0.5 6.5];                            % rating range (worse ratings are clipped)
  MarkColor = cat_io_colormaps('marks+',40);        % quality colors (as command line)
  F.col2mark = @(val) min(1,max(0,MarkColor(min(size(MarkColor,1)-3, ...
    max(1,floor( val/9.5 * size(MarkColor,1)))),:)));
  F.light   = @(c,a) 1 - a*(1 - c);                 % lighter color (a = 0 white, 1 color)
  F.tcol    = [0.15 0.15 0.15];                     % text color
  [F.types,F.tnam] = scanTypes(job.opts);
  F.sites   = unique(textCol(E,'site'))'; F.sites(cellfun('isempty',F.sites)) = [];
  F.sitecol = cat_io_colormaps('nejm',max(1,numel(F.sites)));
  F.sexnam  = {'M','F','O',''}; F.sexlab = {'male','female','other','n/a'};
  F.sexcol  = [0.25 0.45 0.80; 0.85 0.35 0.35; 0.55 0.55 0.55; 0.80 0.80 0.80];
  F.statnam = {'conform','semiconform','nonconform','untested','missing'};
  F.statcol = [F.col2mark(1); F.col2mark(3); F.col2mark(5.5); 0.62 0.62 0.62; 0.90 0.90 0.90];
  F.tisnam  = {'GM','WM','CSF'};
  F.tiscol  = [0.20 0.60 0.30; 0.90 0.55 0.15; 0.25 0.45 0.85];
  F.Q       = Q;

  % volume model and homogeneity threshold of the whole sample (see
  % cat_io_dcm2bids_stats)
  Y = 100 * [numCol(E,'T1w_rGMV') numCol(E,'T1w_rWMV') numCol(E,'T1w_rCSFV')];
  [~,F.M] = cat_io_dcm2bids_stats('volumeModel',Y,numCol(E,'age'),textCol(E,'sex'));
  qmz = numCol(E,'T1w_QMZ'); qmz = qmz(isfinite(qmz));
  if numel(qmz) >= 3, F.qmzthr = mean(qmz) + 2*std(qmz); else, F.qmzthr = nan; end
end
% =========================================================================
function drawPage(fh,F,S,E,site,pgi,npg)
%drawPage. One page of the report figure with the scans S and sessions E
%  of all sites (site = '') or of one site.
  drawHeader(fh,F,S,E,site);
  drawSex(fh,F,E,[1.3 21.9 3.8 3.8]);
  drawAge(fh,F,E,[7.0 21.8 7.0 4.4]);
  if isempty(site) && numel(F.sites) > 1
    drawTimeline(fh,F,E,[1.6 17.9 10.2 2.5]);
    drawSiteTable(fh,F,[13.0 17.2 7.0 3.6]);
  else
    drawTimeline(fh,F,E,[1.6 17.9 18.4 2.5]);
  end
  drawCompleteness(fh,F,S,E,[2.4 12.9 7.0 3.3]);
  drawRatings(fh,F,E,{'SQR','PFR','SOR','BBR','PRR'},{'SQR','PFR','SOR','BBR','PRR'}, ...
    [12.0 12.9 8.0 3.3],'Session ratings (combined ratings of the scans)');
  drawRatings(fh,F,E,strcat(F.tnam,'_SQR'),F.tnam,[1.6 7.7 8.2 3.3], ...
    'Image quality SQR (best scan of each session)');
  drawVolumes(fh,F,E,[12.0 7.7 5.3 3.3]);
  drawHomogeneity(fh,F,E,[1.6 1.8 8.2 3.9]);
  drawWorst(fh,F,E,[11.0 1.1 9.0 5.1]);

  % footer
  ax = textAxes(fh,F,[1.0 0.3 19.0 0.4]);
  text(ax,0,0.5,sprintf('cat_io_dcm2bids report figure of %s/%s, %s',F.project,F.dataset, ...
    char(datetime('now','Format','yyyy-MM-dd HH:mm'))),'FontSize',F.fs-1,'Color',[.5 .5 .5],'Interpreter','none');
  text(ax,1,0.5,sprintf('page %d/%d',pgi,npg),'FontSize',F.fs-1,'Color',[.5 .5 .5], ...
    'HorizontalAlignment','right');
end
% =========================================================================
function drawHeader(fh,F,S,E,site)
%drawHeader. Title, numbers of subjects, sessions, and scans, the scanner
%  of site pages, and the key ratings as colored boxes.
  if isempty(site) && isscalar(F.sites), site = F.sites{1}; end % project of one site
  ax = textAxes(fh,F,[1.3 27.2 10.5 1.9]);
  text(ax,0,0.90,sprintf('DCM2BIDS report: %s',F.project),'FontSize',F.fs+7,'FontWeight','bold', ...
    'Interpreter','none','Color',F.tcol);
  if isempty(site)
    str = sprintf('BIDS directory %s: overview of %d site(s)',F.dataset,numel(F.sites));
  else
    str = sprintf('BIDS directory %s: site %s',F.dataset,site);
  end
  text(ax,0,0.56,str,'FontSize',F.fs+2,'Interpreter','none','Color',F.tcol);
  nsub = numel(unique(textCol(E,'participant_id')));
  text(ax,0,0.30,sprintf('%d subjects, %d sessions, %d scans',nsub,numel(E),numel(S)), ...
    'FontSize',F.fs,'Interpreter','none','Color',F.tcol);
  if ~isempty(site) % scanner of the site
    sc = strtrim(strcat(textCol(S,'manufacturer'),{' '},textCol(S,'model')));
    fs = numCol(S,'field_strength'); fs = fs(isfinite(fs));
    if ~isempty(fs), sc = strcat(sc,sprintf(' %gT',median(fs))); end
    text(ax,0,0.08,sprintf('%s, software %s',mostFrequent(sc),valueCounts(textCol(S,'software'))), ...
      'FontSize',F.fs-1,'Interpreter','none','Color',F.tcol);
  end

  % key ratings (medians of the sessions and conform scans)
  st = textCol(S,'protocol_status'); st = st(~cellfun('isempty',st));
  if isempty(st), conform = nan; else, conform = 100 * mean(strcmp(st,'conform')); end
  key = {'quality SQR',     median(numCol(E,'SQR'),'omitnan'), 0
         'protocol PFR',    median(numCol(E,'PFR'),'omitnan'), 0
         'orientation SOR', median(numCol(E,'SOR'),'omitnan'), 0
         'protocol conform',conform, 1};
  for ki = 1:size(key,1)
    ax = textAxes(fh,F,[12.15 + (ki-1)*2.0 27.35 1.85 1.35]);
    v  = key{ki,2};
    if key{ki,3}, mark = 1 + 5 * (1 - v/100); vstr = sprintf('%0.0f%%',v); % 100% as 1, 0% as 6
    else,         mark = v; vstr = sprintf('%0.2f',v);
    end
    if isfinite(v), col = F.light(F.col2mark(mark),0.6); else, col = [.92 .92 .92]; vstr = 'n/a'; end
    patch(ax,[0 1 1 0],[0 0 1 1],col,'EdgeColor',[.6 .6 .6]);
    text(ax,0.5,0.75,key{ki,1},'FontSize',F.fs-1,'HorizontalAlignment','center','Color',F.tcol);
    text(ax,0.5,0.32,vstr,'FontSize',F.fs+4,'FontWeight','bold','HorizontalAlignment','center','Color',F.tcol);
  end
  ax = textAxes(fh,F,[1.0 26.9 19.0 0.1]);
  line(ax,[0 1],[0.5 0.5],'Color',[.6 .6 .6],'LineWidth',0.5);
end
% =========================================================================
function drawSex(fh,F,E,pos)
%drawSex. Pie of the sex of the subjects with the number and percentage.
  [~,ui] = unique(textCol(E,'participant_id'),'stable');
  sex = textCol(E(ui),'sex');
  sex(~ismember(sex,F.sexnam(1:3))) = {''};
  cnt = cellfun(@(s) sum(strcmp(sex,s)),F.sexnam);
  ax  = newAxes(fh,F,pos);
  if sum(cnt) > 0
    nz = find(cnt > 0);
    h  = pie(ax,cnt(nz));
    for i = 1:numel(nz) % percentage in the slice (the numbers are in the age table)
      set(h(2*i-1),'FaceColor',F.sexcol(nz(i),:),'EdgeColor','w');
      p = 100 * cnt(nz(i)) / sum(cnt); str = '';
      if p >= 8, str = sprintf('%0.0f%%',p); end
      set(h(2*i),'String',str,'Position',0.55 * get(h(2*i),'Position'),'FontSize',F.fs, ...
        'FontWeight','bold','Color','w','HorizontalAlignment','center');
    end
  end
  axis(ax,'equal','off');
  title(ax,sprintf('Sex (%d subjects)',numel(sex)),'FontSize',F.fs+1,'Visible','on');
end
% =========================================================================
function drawAge(fh,F,E,pos)
%drawAge. Age pyramid of the subjects (age at the first session, males
%  left, females right) with the counts of each bin and the age statistics
%  of each sex next to the diagram.
  [~,ui] = unique(textCol(E,'participant_id'),'stable');
  age = numCol(E(ui),'age'); sex = textCol(E(ui),'sex');
  edges = F.ropts.agerange(1):F.ropts.agestep:F.ropts.agerange(2);
  cen   = edges(1:end-1) + diff(edges)/2;
  agec  = min(edges(end),max(edges(1),age));      % ages outside the range in the first/last bin
  cM = histcounts(agec(strcmp(sex,'M')),edges);
  cF = histcounts(agec(strcmp(sex,'F')),edges);
  ax = newAxes(fh,F,pos);
  title(ax,sprintf('Age at the first session (%d-year bins)',F.ropts.agestep),'FontSize',F.fs+1);
  hold(ax,'on');
  xm = max([cM cF 1]) * 1.3;
  barh(ax,cen,-cM,1,'FaceColor',F.sexcol(1,:),'EdgeColor','w','LineWidth',0.3);
  barh(ax,cen, cF,1,'FaceColor',F.sexcol(2,:),'EdgeColor','w','LineWidth',0.3);
  for bi = find(cM > 0), text(ax,-cM(bi) - xm*0.02,cen(bi),num2str(cM(bi)),'FontSize',F.fs-1.5, ...
      'HorizontalAlignment','right','VerticalAlignment','middle','Color',F.tcol); end
  for bi = find(cF > 0), text(ax, cF(bi) + xm*0.02,cen(bi),num2str(cF(bi)),'FontSize',F.fs-1.5, ...
      'HorizontalAlignment','left','VerticalAlignment','middle','Color',F.tcol); end
  plot(ax,[0 0],edges([1 end]),'Color',[.4 .4 .4],'LineWidth',0.5);
  xt = 0:niceStep(max([cM cF 1])/3):max([cM cF 1]);
  set(ax,'XLim',[-xm xm],'YLim',edges([1 end]),'YTick',edges(1):10:edges(end), ...
    'XTick',unique([-xt xt]),'XTickLabel',arrayfun(@(x) num2str(abs(x)),unique([-xt xt]),'UniformOutput',false));
  xlabel(ax,'subjects','FontSize',F.fs); ylabel(ax,'age in years','FontSize',F.fs);
  text(ax,-xm*0.95,edges(end) - diff(edges([1 end]))*0.04,F.sexlab{1},'FontSize',F.fs,'Color',F.sexcol(1,:), ...
    'FontWeight','bold','VerticalAlignment','top');
  text(ax, xm*0.95,edges(end) - diff(edges([1 end]))*0.04,F.sexlab{2},'FontSize',F.fs,'Color',F.sexcol(2,:), ...
    'FontWeight','bold','VerticalAlignment','top','HorizontalAlignment','right');

  % age statistics of each sex
  ax  = textAxes(fh,F,[pos(1)+pos(3)+0.5 pos(2) 20 - pos(1) - pos(3) - 0.5 pos(4)]);
  grp = {'all',true(size(sex)); 'male',strcmp(sex,'M'); 'female',strcmp(sex,'F'); ...
         'other/n/a',~strcmp(sex,'M') & ~strcmp(sex,'F')};
  xc  = [0 0.26 0.50 0.82];
  hdr = {'','n (%)','mean \pm SD','range'};
  for ci = 1:4, text(ax,xc(ci),0.92,hdr{ci},'FontSize',F.fs,'FontWeight','bold','Color',F.tcol); end
  for gi = 1:size(grp,1)
    a = age(grp{gi,2} & isfinite(age)); y = 0.92 - gi*0.16;
    if gi == 1, col = F.tcol; elseif gi == 4, col = F.sexcol(3,:); else, col = F.sexcol(gi-1,:); end
    text(ax,xc(1),y,grp{gi,1},'FontSize',F.fs,'Color',col);
    text(ax,xc(2),y,sprintf('%d (%0.0f%%)',sum(grp{gi,2}),100*sum(grp{gi,2})/max(1,numel(sex))),'FontSize',F.fs,'Color',F.tcol);
    if numel(a) > 1
      text(ax,xc(3),y,sprintf('%0.1f \\pm %0.1f',mean(a),std(a)),'FontSize',F.fs,'Color',F.tcol);
    elseif isscalar(a)
      text(ax,xc(3),y,sprintf('%0.1f',a),'FontSize',F.fs,'Color',F.tcol);
    end
    if ~isempty(a)
      text(ax,xc(4),y,sprintf('%0.0f-%0.0f',min(a),max(a)),'FontSize',F.fs,'Color',F.tcol);
    end
  end
  na = sum(~isfinite(age));
  if na > 0, text(ax,0,0.92 - 5*0.16,sprintf('%d subject(s) without age',na),'FontSize',F.fs-1,'Color',[.4 .4 .4]); end
end
% =========================================================================
function drawTimeline(fh,F,E,pos)
%drawTimeline. Sessions per month (stacked by site) of the session date of
%  the session labels (ses-YYYYMMDD) with the number of each month.
  sid = textCol(E,'session_id');
  tok = regexp(sid,'^ses-(\d{4})(\d{2})\d{2}','tokens','once');
  ym  = nan(numel(sid),2);
  for i = 1:numel(tok), if ~isempty(tok{i}), ym(i,:) = str2double(tok{i}); end, end
  ok  = ym(:,1) >= 1980 & ym(:,1) <= 2100 & ym(:,2) >= 1 & ym(:,2) <= 12;
  ax  = newAxes(fh,F,pos);
  if ~any(ok)
    title(ax,'Sessions per month','FontSize',F.fs+1); axis(ax,'off');
    text(ax,0.5,0.5,'no session dates (e.g. coded session labels of anonymize > 1)','FontSize',F.fs, ...
      'HorizontalAlignment','center','Color',[.4 .4 .4]);
    return
  end
  mi   = (ym(:,1) - 1) * 12 + ym(:,2);          % month index
  m0   = min(mi(ok)); nm = max(mi(ok)) - m0 + 1;
  site = textCol(E,'site');
  sites = F.sites; if isempty(sites), sites = {''}; end
  cnt = @(sel) accumarray(reshape(mi(sel),[],1) - m0 + 1,1,[nm 1]); % sessions per month
  C = zeros(nm,numel(sites));
  for si = 1:numel(sites)
    C(:,si) = cnt(ok & strcmp(site,sites{si}));
  end
  C(:,end+1) = cnt(ok & ~ismember(site,sites)); % sessions without site
  hb = bar(ax,1:nm+1,[C; zeros(1,size(C,2))],1,'stacked','EdgeColor','w','LineWidth',0.2); % (empty row to stack one month)
  for si = 1:numel(hb)
    if si <= numel(sites), set(hb(si),'FaceColor',F.sitecol(si,:)); else, set(hb(si),'FaceColor',[.6 .6 .6]); end
  end
  hold(ax,'on');
  tot = sum(C,2);
  if nm <= 48
    for i = find(tot' > 0)
      text(ax,i,tot(i),num2str(tot(i)),'FontSize',F.fs-2,'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom','Color',F.tcol);
    end
  end
  lab = arrayfun(@(m) sprintf('%04d-%02d',floor((m-1)/12)+1,mod(m-1,12)+1),m0:m0+nm-1,'UniformOutput',false);
  step = max(1,ceil(nm/max(2,round(pos(3)/1.6))));
  set(ax,'XLim',[0.4 nm+0.6],'YLim',[0 max(tot)*1.3 + 0.5],'XTick',1:step:nm,'XTickLabel',lab(1:step:nm));
  ylabel(ax,'sessions','FontSize',F.fs);
  str = sprintf('Sessions per month (%d sessions, %s to %s)',sum(ok),lab{1},lab{end});
  if any(~ok), str = sprintf('%s, %d without date',str,sum(~ok)); end
  title(ax,str,'FontSize',F.fs+1);
end
% =========================================================================
function drawSiteTable(fh,F,pos)
%drawSiteTable. Table of the sites of the project table with the number of
%  subjects, sessions, and scans, the median ratings (colored marks), and
%  the conform scans, i.e., the legend of the site colors.
  Q  = F.Q; qs = textCol(Q,'site');
  rows = [F.sites {'all'}];
  ax = textAxes(fh,F,pos);
  title(ax,'Sites','FontSize',F.fs+1,'Visible','on');
  col = {'site','subj.','sess.','scans','SQR','PFR','SOR','conf.'};
  fld = {'site','subjects','sessions','scans','SQR','PFR','SOR','conform'};
  xc  = [0.04 0.36 0.48 0.60 0.71 0.81 0.91 1.0];
  rh  = min(0.13,0.9/(numel(rows)+1));
  for ci = 1:numel(col)
    text(ax,xc(ci),0.95,col{ci},'FontSize',F.fs-1,'FontWeight','bold','Color',F.tcol, ...
      'HorizontalAlignment',align(ci));
  end
  for ri = 1:numel(rows)
    qi = find(strcmp(qs,rows{ri}),1); y = 0.95 - ri*rh;
    if isempty(qi), continue; end
    if ri <= numel(F.sites)
      patch(ax,[0 0.025 0.025 0],y + rh*[-0.3 -0.3 0.3 0.3],F.sitecol(ri,:),'EdgeColor','none');
    end
    for ci = 1:numel(col)
      v = nan; if isfield(Q,fld{ci}), v = Q(qi).(fld{ci}); end
      if ischar(v), str = v;
      elseif ci >= 5 && ci <= 7 && isfinite(v)
        str = sprintf('%0.2f',v);
        patch(ax,xc(ci) + [-0.085 0.005 0.005 -0.085],y + rh*[-0.42 -0.42 0.42 0.42], ...
          F.light(F.col2mark(v),0.6),'EdgeColor','none');
      elseif ci == 8 && isfinite(v), str = sprintf('%0.0f%%',v);
      elseif isfinite(v), str = sprintf('%d',round(v));
      else, str = 'n/a';
      end
      fw = 'normal'; if ri == numel(rows), fw = 'bold'; end % all sites
      text(ax,xc(ci),y,str,'FontSize',F.fs-1,'Color',F.tcol,'HorizontalAlignment',align(ci), ...
        'FontWeight',fw,'Interpreter','none');
    end
  end
  function a = align(ci)
    if ci == 1, a = 'left'; else, a = 'right'; end
  end
end
% =========================================================================
function drawCompleteness(fh,F,S,E,pos)
%drawCompleteness. Protocol status of the best scan of each typical scan
%  type and session (conform, semiconform, nonconform, without protocol
%  test, or missing scan) in percent of the sessions.
  ax  = newAxes(fh,F,pos);
  sk  = strcat(textCol(S,'participant_id'),'|',textCol(S,'session_id'));
  ek  = strcat(textCol(E,'participant_id'),'|',textCol(E,'session_id'));
  st  = textCol(S,'protocol_status'); ty = textCol(S,'type');
  rnk = 4 * ones(numel(S),1);                    % without protocol test
  for ri = 1:3, rnk(strcmp(st,F.statnam{ri})) = ri; end
  nt  = numel(F.types); P = zeros(nt,5); nses = numel(E);
  for ti = 1:nt
    for ei = 1:nses
      r = rnk(strcmp(sk,ek{ei}) & strcmp(ty,F.types{ti}));
      if isempty(r), P(ti,5) = P(ti,5) + 1; else, P(ti,min(r)) = P(ti,min(r)) + 1; end
    end
  end
  P  = 100 * P / max(1,nses);
  hb = barh(ax,1:nt+1,[P; zeros(1,5)],0.7,'stacked','EdgeColor','w','LineWidth',0.3); % (empty row to stack one type)
  for ci = 1:5, set(hb(ci),'FaceColor',F.statcol(ci,:)); end
  hold(ax,'on');
  for ti = 1:nt
    x0 = 0;
    for ci = 1:5
      if P(ti,ci) >= 9
        text(ax,x0 + P(ti,ci)/2,ti,sprintf('%0.0f%%',P(ti,ci)),'FontSize',F.fs-1.5, ...
          'HorizontalAlignment','center','Color',F.tcol);
      end
      x0 = x0 + P(ti,ci);
    end
    text(ax,102,ti,sprintf('%0.0f/%d',nses*(1 - P(ti,5)/100),nses),'FontSize',F.fs-1,'Color',F.tcol);
  end
  set(ax,'YDir','reverse','YLim',[0.4 nt+0.6],'YTick',1:nt,'YTickLabel',F.tnam,'XLim',[0 100], ...
    'TickLabelInterpreter','none');
  xlabel(ax,'sessions in %','FontSize',F.fs);
  title(ax,'Completeness and protocol status (sessions)','FontSize',F.fs+1);
  % legend below the axes
  lax = textAxes(fh,F,[pos(1)-0.8 pos(2)-1.15 pos(3)+1.6 0.35]);
  xl  = linspace(0,1,6);
  for ci = 1:5
    patch(lax,xl(ci) + [0 0.03 0.03 0],[0.15 0.15 0.85 0.85],F.statcol(ci,:),'EdgeColor',[.7 .7 .7]);
    text(lax,xl(ci) + 0.04,0.5,F.statnam{ci},'FontSize',F.fs-1.5,'Color',F.tcol);
  end
end
% =========================================================================
function drawRatings(fh,F,E,fields,names,pos,ttl)
%drawRatings. Box plots of the ratings of the fields of the sessions E
%  (marks 0.5 to 6.5 with the colors of the marks, worse ratings are
%  clipped) with the data points (colored by site), the median, and the
%  number of sessions of each rating.
  ax = newAxes(fh,F,pos);
  n  = numel(fields); data = cell(1,n); site = textCol(E,'site');
  for fi = 1:n, data{fi} = numCol(E,fields{fi}); end
  nv = cellfun(@(d) sum(isfinite(d)),data);
  if ~any(nv)
    title(ax,ttl,'FontSize',F.fs+1); axis(ax,'off');
    text(ax,0.5,0.5,'no ratings','FontSize',F.fs,'HorizontalAlignment','center','Color',[.4 .4 .4]);
    return
  end
  datac = cellfun(@(d) min(F.marks(2),max(F.marks(1),d(isfinite(d)))),data,'UniformOutput',false);
  datac(nv == 0) = {nan};
  set(0,'CurrentFigure',fh); set(fh,'CurrentAxes',ax);
  cat_plot_boxplot(datac,struct('names',{names},'ylim',F.marks,'ygrid',0,'fontsize',F.fs, ...
    'groupcolor',repmat([.55 .55 .62],n,1),'showdata',0,'outliers',0,'median',2,'trans',0.4,'sat',0.4));
  hold(ax,'on');
  % mark colors in the background (as the QC ratings)
  hp = gobjects(1,6);
  for mi = 1:6
    hp(mi) = patch(ax,[0.4 n+1 n+1 0.4],mi + [-0.5 -0.5 0.5 0.5],F.light(F.col2mark(mi),0.35),'EdgeColor','none');
  end
  uistack(hp,'bottom');
  % data points (by site, behind the boxes, clipped ratings as triangles),
  % median, and number of each rating
  hd = gobjects(0);
  for fi = find(nv > 0)
    d = data{fi}; ok = isfinite(d); d = d(ok); s = site(ok);
    [~,so] = sort(d); x = zeros(size(d));
    x(so) = fi + 0.32 * (mod((0:numel(d)-1)' * 0.618,1) - 0.5);
    dc = min(F.marks(2),max(F.marks(1),d));
    mk = repmat('o',size(d)); mk(d > F.marks(2)) = '^'; mk(d < F.marks(1)) = 'v';
    [us,~,ui] = unique(s);
    for si = 1:numel(us)
      c = siteColor(F,us{si});
      for m = 'o^v'
        g = ui == si & mk == m;
        if any(g)
          hd(end+1) = plot(ax,x(g),dc(g),m,'MarkerSize',2.5,'MarkerFaceColor',c, ...
            'MarkerEdgeColor',F.light(c,0.6),'LineWidth',0.2); %#ok<AGROW>
        end
      end
    end
    md = median(d);
    text(ax,fi + 0.40,min(F.marks(2)-0.2,max(F.marks(1)+0.2,md)),sprintf('%0.2f',md),'FontSize',F.fs-1.5, ...
      'Color',F.tcol,'FontWeight','bold');
  end
  if ~isempty(hd), uistack(hd,'bottom'); end
  uistack(hp,'bottom');
  set(ax,'XLim',[0.4 n+1],'YLim',F.marks,'YTick',1:6,'YDir','normal','Layer','top', ...
    'XTick',1:n,'TickLabelInterpreter','tex','XTickLabel',cellfun(@(nm,k) ...
    sprintf('%s\\newline{\\fontsize{%g}n=%d}',strrep(nm,'_','\_'),F.fs-2,k),names,num2cell(nv),'UniformOutput',false));
  ylabel(ax,'mark','FontSize',F.fs);
  title(ax,ttl,'FontSize',F.fs+1);
end
% =========================================================================
function drawVolumes(fh,F,E,pos)
%drawVolumes. Relative tissue volumes of the best T1w scan of each session
%  by age (marker by sex) with the robust age/sex model of the sample and
%  the volume statistics next to the diagram (z-scores of the session
%  table, see cat_io_dcm2bids_stats).
  ax  = newAxes(fh,F,pos);
  title(ax,'Tissue volumes (best T1w scan)','FontSize',F.fs+1);
  age = numCol(E,'age'); sex = textCol(E,'sex');
  fn  = {'T1w_rGMV','T1w_rWMV','T1w_rCSFV'};
  Y   = 100 * [numCol(E,fn{1}) numCol(E,fn{2}) numCol(E,fn{3})];
  Z   = [numCol(E,[fn{1} '_z']) numCol(E,[fn{2} '_z']) numCol(E,[fn{3} '_z'])];
  ok  = isfinite(age) & any(isfinite(Y),2);
  if ~any(ok)
    axis(ax,'off');
    text(ax,0.5,0.5,{'no T1w segmentation','(preprocessing)'},'FontSize',F.fs,'HorizontalAlignment','center','Color',[.4 .4 .4]);
    return
  end
  hold(ax,'on');
  mk = {'s','o','^','^'}; xr = [min(age(ok)) max(age(ok))] + [-1 1];
  for ti = 1:3
    for gi = 1:4
      g = ok & isfinite(Y(:,ti)) & strcmp(sex,F.sexnam{gi});
      if gi == 4, g = ok & isfinite(Y(:,ti)) & ~ismember(sex,F.sexnam(1:3)); end
      plot(ax,age(g),Y(g,ti),mk{gi},'MarkerSize',2.5,'MarkerFaceColor',F.tiscol(ti,:), ...
        'MarkerEdgeColor',F.light(F.tiscol(ti,:),0.5),'LineWidth',0.2);
    end
    o = ok & abs(Z(:,ti)) > 3; % outliers of the model
    plot(ax,age(o),Y(o,ti),'o','MarkerSize',5,'MarkerEdgeColor','k','LineWidth',0.5);
    M = F.M(ti); % model at the mean sex (median for the robust z-score)
    if M.n > 0
      xa = [0 0]; if isfinite(M.mage), xa = xr - M.mage; end
      plot(ax,xr,M.b(1) + M.b(2) * xa,'-','Color',F.tiscol(ti,:),'LineWidth',1);
    end
  end
  set(ax,'XLim',xr,'YLim',[0 max(60,ceil(max(Y(:))/10)*10)]);
  xlabel(ax,'age in years','FontSize',F.fs); ylabel(ax,'volume in % of TIV','FontSize',F.fs);

  % statistics of the volumes next to the diagram
  tax = textAxes(fh,F,[pos(1)+pos(3)+0.3 pos(2) 20.2 - pos(1) - pos(3) - 0.3 pos(4)+0.4]);
  tiv = numCol(E,'T1w_TIV'); y = 0.97; dy = 0.105;
  for ti = 1:3
    v = Y(isfinite(Y(:,ti)),ti); M = F.M(ti);
    text(tax,0,y,sprintf('%s %0.1f \\pm %0.1f%%',F.tisnam{ti},mean(v),std(v)),'FontSize',F.fs-1, ...
      'Color',F.tiscol(ti,:),'FontWeight','bold'); y = y - dy;
    if startsWith(M.method,'robust linear')
      text(tax,0,y,sprintf('  %+0.2f%%/10 years',10*M.b(2)),'FontSize',F.fs-1.5,'Color',F.tcol); y = y - dy;
    end
  end
  t = tiv(isfinite(tiv));
  if ~isempty(t)
    text(tax,0,y,sprintf('TIV %0.0f \\pm %0.0f ml',mean(t),std(t)),'FontSize',F.fs-1,'Color',F.tcol); y = y - dy;
  end
  text(tax,0,y,sprintf('n = %d sessions',sum(ok)),'FontSize',F.fs-1.5,'Color',F.tcol); y = y - dy;
  text(tax,0,y,sprintf('|z| > 3: %d (black ring)',sum(any(abs(Z(ok,:)) > 3,2))),'FontSize',F.fs-1.5,'Color',F.tcol); y = y - dy;
  plot(tax,0.04,y,'s','MarkerSize',3,'MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','none');
  text(tax,0.10,y,'male','FontSize',F.fs-1.5,'Color',F.tcol);
  plot(tax,0.50,y,'o','MarkerSize',3,'MarkerFaceColor',[.5 .5 .5],'MarkerEdgeColor','none');
  text(tax,0.56,y,'female','FontSize',F.fs-1.5,'Color',F.tcol);
end
% =========================================================================
function drawHomogeneity(fh,F,E,pos)
%drawHomogeneity. Sample homogeneity of the mwc1 map of the best T1w scan
%  (quartic mean Z-score, QMZ) by the quality rating of the T1w scans
%  (colored by site) with the threshold of atypical maps (mean + 2 SD of
%  the whole sample, as cat_stat_homogeneity) and the labels of the
%  sessions above the threshold.
  ax  = newAxes(fh,F,pos);
  ttl = 'Sample homogeneity of the mwc1 maps (T1w)';
  title(ax,ttl,'FontSize',F.fs+1);
  qmz = numCol(E,'T1w_QMZ'); sqr = numCol(E,'T1w_SQR'); site = textCol(E,'site');
  ok  = isfinite(qmz) & isfinite(sqr);
  if ~any(ok)
    axis(ax,'off');
    text(ax,0.5,0.5,{'no sample homogeneity','(at least 3 mwc1 maps of the preprocessing)'},'FontSize',F.fs, ...
      'HorizontalAlignment','center','Color',[.4 .4 .4]);
    return
  end
  hold(ax,'on');
  yr = [min(qmz(ok)) max([qmz(ok); F.qmzthr])]; yr = yr + [-0.1 0.15] * max(diff(yr),0.05);
  for mi = 1:6 % mark colors of the T1w quality
    patch(ax,mi + [-0.5 0.5 0.5 -0.5],yr([1 1 2 2]),F.light(F.col2mark(mi),0.35),'EdgeColor','none');
  end
  if isfinite(F.qmzthr)
    plot(ax,F.marks,[F.qmzthr F.qmzthr],'--','Color',[.6 0 0],'LineWidth',0.6);
    text(ax,F.marks(2) - 0.05,F.qmzthr,sprintf('mean + 2 SD = %0.3f',F.qmzthr),'FontSize',F.fs-1.5, ...
      'Color',[.6 0 0],'HorizontalAlignment','right','VerticalAlignment','bottom');
  end
  x = min(F.marks(2),max(F.marks(1),sqr));
  [us,~,ui] = unique(site);
  for si = 1:numel(us)
    g = ok & ui == si; c = siteColor(F,us{si});
    if any(g), plot(ax,x(g),qmz(g),'o','MarkerSize',3,'MarkerFaceColor',c,'MarkerEdgeColor',F.light(c,0.6),'LineWidth',0.2); end
  end
  % labels of the (at most 6) sessions above the threshold without overlap
  hi  = find(ok & qmz > F.qmzthr);
  [~,so] = sort(qmz(hi),'descend'); hl = hi(so(1:min(end,6)));
  yl  = qmz(hl); for i = 2:numel(yl), yl(i) = min(yl(i),yl(i-1) - 0.07*diff(yr)); end
  sub = textCol(E,'participant_id'); ses = textCol(E,'session_id');
  for i = 1:numel(hl)
    if x(hl(i)) > 4, xo = -0.12; ha = 'right'; else, xo = 0.12; ha = 'left'; end
    text(ax,x(hl(i)) + xo,yl(i),shortLabel([sub{hl(i)} ' ' ses{hl(i)}],34),'FontSize',F.fs-2.5,'Color',F.tcol, ...
      'Interpreter','none','HorizontalAlignment',ha);
  end
  set(ax,'XLim',F.marks,'YLim',yr,'XTick',1:6,'Layer','top');
  xlabel(ax,'T1w quality SQR (mark, best T1w scan)','FontSize',F.fs); ylabel(ax,'QMZ','FontSize',F.fs);
  str = sprintf('\\fontsize{%g}n = %d, QMZ %0.3f \\pm %0.3f, above threshold: %d',F.fs-1,sum(ok), ...
    mean(qmz(ok)),std(qmz(ok)),numel(hi));
  if any(isfinite(qmz) & ~ok), str = sprintf('%s (%d without T1w SQR)',str,sum(isfinite(qmz) & ~ok)); end
  title(ax,{ttl,str},'FontSize',F.fs+1);
end
% =========================================================================
function drawWorst(fh,F,E,pos)
%drawWorst. Sessions to check: sessions with the worst ratings (combined
%  and per scan type) and atypical mwc1 maps with the reasons, i.e., the
%  ratings of mark 3.5 or worse (at most 3), otherwise the worst rating.
  ax  = textAxes(fh,F,pos);
  title(ax,'Sessions to check (worst ratings)','FontSize',F.fs+1,'Visible','on');
  fld = [{'SQR'} strcat(F.tnam,'_SQR') {'PFR','SOR','BBR','PRR'}];
  lab = [{'SQR'} strcat(F.tnam,' SQR') {'PFR','SOR','BBR','PRR'}];
  R   = cell2mat(cellfun(@(f) numCol(E,f),fld,'UniformOutput',false));
  qmz = numCol(E,'T1w_QMZ');
  worst = max(R,[],2,'omitnan');
  score = worst; score(qmz > F.qmzthr) = max(score(qmz > F.qmzthr),4.5); % atypical maps at least as mark 4.5
  [~,so] = sort(score,'descend','MissingPlacement','last');
  so = so(isfinite(score(so))); so = so(1:min(end,F.ropts.nworst));
  if isempty(so)
    text(ax,0.5,0.5,'no ratings','FontSize',F.fs,'HorizontalAlignment','center','Color',[.4 .4 .4]);
    return
  end
  sub = textCol(E,'participant_id'); ses = textCol(E,'session_id'); site = textCol(E,'site');
  xc  = [0 0.42 0.58 0.72]; rh = 0.95 / (F.ropts.nworst + 1);
  hdr = {'subject / session','site','worst','reasons (mark 3.5 or worse)'};
  for ci = 1:4, text(ax,xc(ci),0.97,hdr{ci},'FontSize',F.fs-1,'FontWeight','bold','Color',F.tcol); end
  for ri = 1:numel(so)
    ei = so(ri); y = 0.97 - ri*rh;
    text(ax,xc(1),y + rh*0.18,shortLabel(sub{ei},30),'FontSize',F.fs-1.5,'Color',F.tcol,'Interpreter','none');
    text(ax,xc(1),y - rh*0.22,ses{ei},'FontSize',F.fs-1.5,'Color',[.4 .4 .4],'Interpreter','none');
    text(ax,xc(2),y,site{ei},'FontSize',F.fs-1.5,'Color',F.tcol,'Interpreter','none');
    if isfinite(worst(ei))
      patch(ax,xc(3) + [0 0.1 0.1 0],y + rh*[-0.35 -0.35 0.35 0.35],F.light(F.col2mark(worst(ei)),0.6),'EdgeColor','none');
      text(ax,xc(3) + 0.05,y,sprintf('%0.2f',worst(ei)),'FontSize',F.fs-1,'FontWeight','bold', ...
        'HorizontalAlignment','center','Color',F.tcol);
    end
    [v,vi] = sort(R(ei,:),'descend','MissingPlacement','last');
    sel = find(v >= 3.5 & isfinite(v),3); if isempty(sel) && isfinite(v(1)), sel = 1; end
    rs  = arrayfun(@(i) sprintf('%s %0.1f',lab{vi(i)},v(i)),sel,'UniformOutput',false);
    if qmz(ei) > F.qmzthr, rs{end+1} = sprintf('QMZ %0.3f',qmz(ei)); end %#ok<AGROW>
    text(ax,xc(4),y,strjoin(rs,', '),'FontSize',F.fs-1.5,'Color',F.tcol,'Interpreter','none');
  end
end
% =========================================================================
function showSPM(fh)
%showSPM. Copy of the page in the SPM Graphics window in interactive
%  sessions (fonts scaled to the window height).
  if ~usejava('desktop') || ~feature('ShowFigureWindows'), return; end
  try
    fg = spm_figure('GetWin','Graphics'); spm_figure('Clear',fg);
    copyobj(findobj(allchild(fh),'flat','Type','axes'),fg); % all diagrams and texts are axes
    u0 = get(fg,'Units'); set(fg,'Units','centimeters'); p = get(fg,'Position'); set(fg,'Units',u0);
    sc = p(4) / 29.7;
    h  = findall(fg,'-property','FontSize');
    for hi = 1:numel(h), set(h(hi),'FontSize',get(h(hi),'FontSize') * sc); end
  catch e
    cat_io_cprintf('warn','  Report figure not shown in the SPM Graphics window: %s\n',e.message);
  end
end
% =========================================================================
function ax = newAxes(fh,F,pos)
%newAxes. Axes of a diagram at the position pos ([x y w h] in cm).
  ax = axes('Parent',fh,'Units','normalized','Position',[pos(1)/F.pw pos(2)/F.ph pos(3)/F.pw pos(4)/F.ph], ...
    'FontSize',F.fs,'Box','on','Layer','top','Color','w','XColor',[.3 .3 .3],'YColor',[.3 .3 .3], ...
    'TickLength',[0.005 0.005],'LineWidth',0.5);
  ax.Toolbar = []; disableDefaultInteractivity(ax); % no interactive elements (faster)
end
% =========================================================================
function ax = textAxes(fh,F,pos)
%textAxes. Invisible axes (0 to 1) for text at the position pos (in cm).
  ax = axes('Parent',fh,'Units','normalized','Position',[pos(1)/F.pw pos(2)/F.ph pos(3)/F.pw pos(4)/F.ph], ...
    'XLim',[0 1],'YLim',[0 1],'Visible','off','HitTest','off');
  ax.Toolbar = []; disableDefaultInteractivity(ax);
  hold(ax,'on');
end
% =========================================================================
function st = niceStep(v)
%niceStep. Step of 1, 2, or 5 times a power of 10 (at least 1) of the value v.
  e  = 10^floor(log10(max(v,1))); f = max(v,1) / e;
  m  = [1 2 5 10]; st = e * m(find(m >= f - eps,1));
end
% =========================================================================
function c = siteColor(F,site)
%siteColor. Color of the site (gray without site).
  si = find(strcmp(F.sites,site),1);
  if isempty(si), c = [.5 .5 .5]; else, c = F.sitecol(si,:); end
end
% =========================================================================
function s = shortLabel(s,n)
%shortLabel. Label of at most n characters (shortened in the middle).
  if numel(s) > n, s = [s(1:ceil(n/2)-1) '...' s(end-floor(n/2)+3:end)]; end
end
% =========================================================================
function [types,names] = scanTypes(opts)
%scanTypes. Typical scan types (datatype/suffix, see opts.scantypes) and
%  their column names (as in cat_io_dcm2bids_report).
  types = cellstr(opts.scantypes);
  names = regexprep(types,'^[^/]*/','');
  dup   = cellfun(@(n) sum(strcmp(names,n)) > 1, names);
  names(dup) = strrep(types(dup),'/','_');
  names = regexprep(names,'[^a-zA-Z0-9_]','_');
end
% =========================================================================
function R = readTable(P)
%readTable. Report table (csv or tsv, see writeTable of
%  cat_io_dcm2bids_helper) as structure array (one element per row) with
%  numeric columns (n/a as NaN) and text columns (n/a as '').
  R = struct([]);
  if ~exist(P,'file'), return; end
  if strcmpi(spm_file(P,'ext'),'tsv'), dl = sprintf('\t'); else, dl = ','; end
  C = readcell(P,'FileType','text','Delimiter',dl,'DatetimeType','text','DurationType','text');
  if size(C,1) < 2, return; end
  hdr = matlab.lang.makeValidName(cellfun(@char,C(1,:),'UniformOutput',false));
  V   = C(2:end,:); n = size(V,1);
  isna = cellfun(@(x) isa(x,'missing') || (ischar(x) && strcmp(x,'n/a')),V);
  isnm = cellfun(@(x) isnumeric(x) || islogical(x),V);
  vals = cell(1,numel(hdr));
  for ci = 1:numel(hdr)
    if all(isna(:,ci) | isnm(:,ci))
      x = nan(n,1); x(isnm(:,ci)) = cellfun(@double,V(isnm(:,ci),ci));
      vals{ci} = num2cell(x);
    else
      v = V(:,ci); v(isna(:,ci)) = {''};
      v(isnm(:,ci)) = cellfun(@(x) sprintf('%.15g',x),v(isnm(:,ci)),'UniformOutput',false);
      vals{ci} = v;
    end
  end
  R = cell2struct([vals{:}],hdr,2);
end
% =========================================================================
function x = numCol(R,fn)
%numCol. Numeric column fn of the records R as column vector (NaN if not
%  available or not numeric).
  x = nan(numel(R),1);
  if isempty(R) || ~isfield(R,fn), return; end
  for i = 1:numel(R)
    v = R(i).(fn);
    if (isnumeric(v) || islogical(v)) && ~isempty(v), x(i) = double(v(1)); end
  end
end
% =========================================================================
function c = textCol(R,fn)
%textCol. Text column fn of the records R as column cell ('' if not
%  available), also of numeric values (e.g. a site named by the serial
%  number of the scanner).
  c = repmat({''},numel(R),1);
  if isempty(R) || ~isfield(R,fn), return; end
  for i = 1:numel(R)
    v = R(i).(fn);
    if ischar(v), c{i} = v;
    elseif isnumeric(v) && isscalar(v) && isfinite(v), c{i} = sprintf('%.15g',v);
    end
  end
end
% =========================================================================
function v = mostFrequent(c)
%mostFrequent. Most frequent (non-empty) text of the cell c.
  c = c(~cellfun('isempty',c));
  if isempty(c), v = ''; return; end
  [u,~,ui] = unique(c);
  [~,mi] = max(accumarray(ui(:),1));
  v = u{mi};
end
% =========================================================================
function v = valueCounts(c)
%valueCounts. Different (non-empty) texts of the cell c with their counts
%  (most frequent first).
  c = c(~cellfun('isempty',c));
  if isempty(c), v = 'n/a'; return; end
  [u,~,ui] = unique(c);
  n = accumarray(ui(:),1);
  [n,so] = sort(n,'descend'); u = u(so);
  v = strjoin(arrayfun(@(i) sprintf('%s (%d)',u{i},n(i)),1:numel(u),'UniformOutput',false),', ');
end

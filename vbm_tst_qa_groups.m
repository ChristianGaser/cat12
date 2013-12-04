function vbm_tst_qa_groups(files,names,qm,printgroup,name)
% ______________________________________________________________________
% 
% vbm_tst_qa_groups(files,names,qm,printgroup,name)
% 
% group(i).name  = ''
% group(i).files = {}
%
% Usergroup           #Scanner  #Subjects  #Protocols  #Groups
%  1) Scientist       1+        >20        1           1+          
%  2) Developer       1+        <10        1+          1
%
% ______________________________________________________________________
% Optimierung von QA Parameter:
% 1)  Einzeloptimierung: 
%     Entfernung von negativen Ausreißern mit starkem Rauschen, Inhomogenität, ...
% 2A) Intragruppenoptimierung: 
%     Verringern der Gruppenvarianz durch das Entfernen von Ausreißern 
% 2B) Intergruppenoptimierung: 
%     Verringern der Varianz durch Entfernen von Ausreißern
% ______________________________________________________________________
% $Id$


  %% 0) defaults:

  if ~exist('names','var'), names = cellstr(num2str((1:numel(files))')); end
  namesn = names; for fi=1:numel(files), namesn{fi} = sprintf('%s (%d)',names{fi},numel(files{fi})); end

  if numel(files) > numel(names), error('MATLAB:vbm_tst_qa_groups:groupnames','Error: To small number of group names.\n'); end

  if ~exist('qm','var')
    qm = {
      'noisec'         'low (good) <--- noisec (noise (1/CNR) after correction) ---> strong (bad)'
      'biasc'          'low (good) <--- biasc (inhomogeneity after correction) ---> strong (bad)'
      'NERR'           'low (good) <--- NERR (NoiseEdgeResolutionRelation - artifacts) ---> strong (bad)'
      };
  end

  if ~exist('printgroup','var'), printgroup = 0; end
  if ~exist('name','var'), name = ''; end
  
  gir  = 1:numel(files);
  %gir  = randperm(numel(files)); % char(64 + randperm(numel(files)));
  girn = cell(numel(files),1); for fi=1:numel(files), girn{fi} = sprintf('%02d: %s (%d)',gir(fi),names{fi},numel(files{fi})); end

  


  %% 1) get data from xml files 
  ci=0; data = struct(); tmp = cell(1,1);
  for gi=1:numel(files)
    for fi=1:numel(files{gi})
      ci = ci+1;
      data.id(ci,1)      = gi;
      data.group{ci,1}   = names{gi};
      data.groupn{ci,1}  = namesn{gi};
      data.gir{ci,1}     = gir(gi);
      data.girn{ci,1}    = girn(gi);

      tmp{gi}{fi} = vbm_io_xml(files{gi}{fi});

      for qmi=1:size(qm,1)
        try
          if isinf(tmp{gi}{fi}.qa.(qm{qmi,1})(1))
            data.(qm{qmi,1})(ci,1) = nan;
            data.full(ci,qmi)      = nan;
          else
            data.(qm{qmi,1})(ci,1) = tmp{gi}{fi}.qa.(qm{qmi,1})(1);
            data.full(ci,qmi)      = tmp{gi}{fi}.qa.(qm{qmi,1})(1);
          end
        catch
          data.(qm{qmi,1})(ci,1)   = nan;
          data.full(ci,qmi)        = nan;
        end
      end
    end
  end



  %% 2) estimate median (and other group values) to sort groups
  group = struct(); 
  for qmi=1:size(qm,1)
    for gi=1:numel(files)
      group.(qm{qmi,1}).data{gi} = data.(qm{qmi,1})(strcmp(data.group,names{gi}));
    end
  end
  
  %{
  group = struct(); 
  for qmi=1:size(qm,1)
    group.(qm{qmi,1}).nsortval = [];
    group.(qm{qmi,1}).sortval = [];
    for gi=1:numel(files)
      group.(qm{qmi,1}).data{gi} = data.(qm{qmi,1})(strcmp(data.group,names{gi}));
   
      group.(qm{qmi,1}).n(gi)    = sum(strcmp(data.group,names{gi})); 
      group.(qm{qmi,1}).med(gi)  = vbm_stat_nanmedian(group.(qm{qmi,1}).data{gi});
      group.(qm{qmi,1}).mean(gi) = vbm_stat_nanmean(group.(qm{qmi,1}).data{gi});
      group.(qm{qmi,1}).std(gi)  = vbm_stat_nanstd(group.(qm{qmi,1}).data{gi});
      group.(qm{qmi,1}).sortval  = [ group.(qm{qmi,1}).sortval ;
        repmat(group.(qm{qmi,1}).med(gi)',group.(qm{qmi,1}).n(gi),1)];
        %repmat(group.(qm{qmi,1}).mean(gi)' + (group.(qm{qmi,1}).std(gi)/4)',group.(qm{qmi,1}).n(gi),1)];
      group.(qm{qmi,1}).nsortval  = [ group.(qm{qmi,1}).nsortval ; group.(qm{qmi,1}).med(gi)];

    end
    [tmp,group.(qm{qmi,1}).sorti]  = sort(group.(qm{qmi,1}).sortval);
    [tmp,group.(qm{qmi,1}).sortin] = sort(group.(qm{qmi,1}).nsortval);
  end
%}
%{
  % sort groups 
  datas = struct();
  for qmi=1:size(qm,1)
    datas.(qm{qmi,1}).id           = group.(qm{qmi,1}).sorti;
    datas.(qm{qmi,1}).group        = data.group(group.(qm{qmi,1}).sorti);
    datas.(qm{qmi,1}).groupn       = data.groupn(group.(qm{qmi,1}).sorti);
    datas.(qm{qmi,1}).gir          = data.gir(group.(qm{qmi,1}).sorti);
    datas.(qm{qmi,1}).girn         = data.girn(group.(qm{qmi,1}).sorti);
    for qmsi=1:size(qm,1)
      datas.(qm{qmi,1}).(qm{qmsi}) = data.(qm{qmi,1})(group.(qm{qmi,1}).sorti);
    end
    for gi=1:numel(names)
      datas.(qm{qmi,1}).data{gi}   = group.(qm{qmi,1}).data{group.(qm{qmi,1}).sortin(gi)};
    end
    datas.(qm{qmi,1}).full         = data.full(group.(qm{qmi,1}).sorti,:);
  end
  %namess = names(group.(qm{qmi,1}).sortin);
%}


  %% 3) MANOVA
  [mano.d,mano.p,mano.stats] = manova1(data.full,data.group); 


  %% 4) Create result structures and tables 
  % gmdist = mahalanobis distance between groups
  % mmdist = mahalanobis distance within groups
  % scans  = mahalanobis distance for each scan to the group mean
  tab.gmdist = [{''}, names'; names, num2cell(mano.stats.gmdist)];
  for ni=1:numel(names), 
    tab.mmdist{ni,1} = names{ni};
    tab.mmdist{ni,2} = mean(mano.stats.mdist(strcmp(data.group,names{ni})));
    tab.mmdist{ni,4} = median(mano.stats.mdist(strcmp(data.group,names{ni})));
    tab.mmdist{ni,3} = std(mano.stats.mdist(strcmp(data.group,names{ni})));
    mano.stats.mmdist{ni} = mano.stats.mdist(strcmp(data.group,names{ni}));
%    tab.smdist.(namess{ni}) = [files{ni} num2cell(mano.stats.mdist(strcmp(data.group,names{ni}))')];
  end
  %gi=1; for fi=1:size(tab.smdist.(namess{gi}),1), fprintf('% 80s %8.4f\n',tab.smdist.(namess{gi}){fi,1},tab.smdist.(namess{gi}){fi,2}); end



  %% 5) Display results
  % Creating a nice plot depent strongly on the number of input groups and
  % the main goal of the analysis. If we have to identify each group, only a 
  % limited number of groups is possible in the standard SPM/VBM graphic 
  % window. But, if we only want to now how good one group to a set of other
  % groups is, it is possilbe to mark only this center and compare it to the 
  % othters. 
  % The results are printed only for the standard SPM graphics figure,
  % because it make less sense to create larger figures that can only read on
  % a computer and not printed in a paper. 

  grouplimits(1) = 120; % boxplot 
  grouplimits(2) = 60;  % groupnames/legend
  grouplimits(3) = 20;  % groupnames/legend
  fontsize(1)    = 11;
  fontsize(2)    = fontsize(1) * 0.9;
  fontsize(3)    = fontsize(1) * 0.9 - numel(names)/50;
  plotn          = size(qm,1)+4;
  if numel(names)<grouplimits(3)
    plotwidth    = 0.94;
    pnames       = names;
    lpos         = [0.92,1-7*(0.90/plotn)+numel(names)*0.0065,0,0];
  else 
    plotwidth    = 0.8;
    pnames       = gir;
    lpos         = [0.92,1-5*(0.90/plotn)+numel(names)*0.0065,0,0];
  end
  

  if numel(names)<grouplimits(1)
    hand.spmfig = spm_figure('GetWin','Graphics');
    spm_figure('Clear','Graphics');
   % annotation('textbox',[0.01 0.95 0.20 0.03],'String',name,'FontWeight','bold',...
   %   'FontSize',fontsize(1)*1.1,'EdgeColor','none','BackgroundColor',[1 1 1]);

    % BOXPLOTs for each QM
    for qmi=1:size(qm,1)
      hand.sp{1,qmi} = subplot(plotn,1,qmi);
      vbm_plot_boxplot(group.(qm{qmi,1}).data,struct('names',{pnames},'sort',1));
      set(hand.sp{1,qmi},'FontSize',fontsize(2),'Position',[0.06 - numel(names)/500,1.0-qmi*(0.90/plotn),plotwidth,0.65/plotn]); 
      title([qm{qmi,2} '        '],'FontSize',fontsize(1),'FontWeight','bold'); %,'Position',[40 0.2]);
    end
    
    % BOXPLOT for mahalanobis distance
    qmi=qmi+1;
    hand.sp{1,qmi} = subplot(plotn,1,qmi);
    set(hand.sp{1,qmi},'FontSize',fontsize(2),'Position',[0.06 - numel(names)/500,1.0-qmi*(0.90/plotn),plotwidth,0.65/plotn]);
    vbm_plot_boxplot(mano.stats.mmdist,struct('names',{pnames},'sort',1));
    title('in Mahalanobis distance within group    ','FontSize',fontsize(1),'FontWeight','bold'); 
    
    % BOXPLOT for mahalanobis distance
    qmi=qmi+1; hand.sp{1,qmi} = subplot(plotn,1,qmi);
    set(hand.sp{1,qmi},'FontSize',fontsize(2),'Position',[0.06 - numel(names)/500,1.0-qmi*(0.90/plotn),plotwidth,0.65/plotn]);
    gmdist=mano.stats.gmdist; gmdist(eye(size(gmdist,1))==1)=nan;
    vbm_plot_boxplot(mat2cell(gmdist,ones(1,size(gmdist,1)),size(gmdist,1)),struct('names',{pnames},'sort',1));
    title('mean Mahalanobis distance to other groups      ','FontSize',fontsize(1),'FontWeight','bold'); 
    
    % Groupmatix
    qmi=qmi+1; hand.sp{1,qmi} = subplot(plotn,2,plotn*2-1);
    set(hand.sp{1,qmi},'FontSize',fontsize(2),'Position',[0.06 - numel(names)/500,0.03,0.40,2*0.95/plotn]);
    [tmp,sorti]=sort(median(mano.stats.gmdist,2)); 
    imagesc(mano.stats.gmdist(sorti,sorti)); daspect([1  1  1]);
    set(gca,'XTick',1:numel(names),'XTickLabel',pnames(sorti),'YTick',1:numel(names),'YTickLabel',pnames(sorti));
    title('Mahalanobis distance between group means      ','FontSize',fontsize(1),'FontWeight','bold');  
    
    % MANOVA
    qmi=qmi+1; hand.sp{1,qmi} = subplot(plotn,2,plotn*2);
    set(hand.sp{1,qmi},'FontSize',fontsize(2),'Position',[0.50,0.03,0.45 - 0.2*(numel(names)>grouplimits(2)),2*0.95/plotn]);
    axissize=max(mano.stats.canon,[],1) - min(mano.stats.canon,[],1);
    if axissize(1)>axissize(2)
      hand.gscatter{1} = gscatter(mano.stats.canon(:,1),mano.stats.canon(:,2),data.groupn,[],'.xo+*',[],'off');
      xlabel('C1'); ylabel('C2'); %axis equal;
    else
      hand.gscatter{1} = gscatter(mano.stats.canon(:,2),mano.stats.canon(:,1),data.groupn,[],'.xo+*',[],'off');
      xlabel('C2'); ylabel('C1'); %axis equal;
    end
    %set(hand.gscatter{1},'MarkerSize',markersize,'LineWidth',1)
    %set(hand.sp{1,qmi},'FontSize',fontsize(2));
    qms=qm{1}; for qmi=2:size(qm,1), qms=[qms ', ' qm{qmi,1}]; end %#ok<AGROW>
    title(sprintf('Principle components of a MANOVA test     \n of %s    ',qms),'FontSize',fontsize(1),'FontWeight','bold');
    
    
    colormap jet;
    warning off;
    legend(girn,'Position',lpos);
    warning on;
    colormap(vbm_io_colormaps('marks',50)); 
  end

  %%

  %{
  if printgroup
    for gi=6 %1:numel(files)

      spm_figure('TurnPage',spm_figure('#page'),hand.spmfig);
      hand.spmfig = spm_figure('FindWin','Graphics');

      for qmi=1:size(qm,1)
        sgroup = datas.(qm{qmi,1}).group; ci=1;
        for ni=setxor(gi,1:numel(names))
          sgroup(~cellfun('isempty',strfind(datas.(qm{qmi,1}).group,names{ni}))) = ...
            cellstr(repmat( gir(ci) ,sum(~cellfun('isempty',strfind(datas.(qm{qmi,1}).group,names{ni}))),1)); ci=ci+1;
        end  
        hand.sp{gi+1,qmi}   = subplot(plotn,1,mod(qmi-1,plotn)+1);
        set(hand.sp{gi+1,qmi},'FontSize',fontsize(2));
        %boxplot(datas.(qm{qmi,1}).(qm{qmi,1}),sgroup,'notch','on'); % 'plotstyle','compact')

        title(qm{qmi,2},'FontSize',fontsize(1),'FontWeight','bold');
      end
      %ylim([0,round(6*max(datas.(qm{qmi,1}).(qm{qmi,1})))/4]);

      % 
      qmi = qmi + 1;
      if qmi>1 && mod(qmi,plotn)==0,
        spm_figure('NewPage',hand.spmfig); 
        spm_figure('TurnPage',spm_figure('#page'),hand.spmfig);
        %hand.spmfig = spm_figure('FindWin','Graphics');
        qmi = qmi + 1;
      end

      sgroup = data.groupn; ci=1;
      for ni=setxor(gi,1:numel(names))
        sgroup(~cellfun('isempty',strfind(data.group,names{ni}))) = ...
          cellstr(repmat( girn{ci} ,sum(~cellfun('isempty',strfind(data.group,names{ni}))),1)); ci=ci+1;
      end  
      %% groupcolor = repmat([0 0 1],numel(sgroup),1); groupcolor(gi,:)=[1 0 0];
      groupcolor = max(0,colormap(jet(round(size(files,1)*1.4)))-0.1);
      groupcolor(groupcolor(:,1)>0.8 & groupcolor(:,2)>0.8 & groupcolor(:,3)>0.8,:)=[];
      groupcolor(sum(groupcolor(:,1:2),2)>1.7,:)=max(0,groupcolor(sum(groupcolor(:,1:2),2)>1.7,:) - 0.3);
      groupcolor(size(files,1)+1:end,:)=[]; groupcolor(gi,:)=[0.9 0 0];
      groupmarker = repmat('x+o.sd^v><ph',1,numel(sgroup)); groupmarker(gi)='*';
      hand.sp{gi+1,qmi}   = subplot(plotn,1,mod(qmi-1,plotn)+1 : mod(qmi-1,plotn)+2); 
      %hand.gscatter{gi+1} = gscatter(mano.stats.canon(:,2),mano.stats.canon(:,1),sgroup,cold,'.o+*');
      hand.gscatter{gi+1} = gscatter(mano.stats.canon(:,2),mano.stats.canon(:,1),sgroup,groupcolor,groupmarker);
      set(hand.gscatter{gi+1},'MarkerSize',markersize,'LineWidth',1)
      set(hand.sp{gi+1,qmi},'FontSize',fontsize(2));
      title(sprintf('Principle Components of a %s \n of %s','Manova',qms),'FontSize',fontsize(1),'FontWeight','bold');

      % spm_figure('NewPage',hand.spmfig); 
    end

  end
  %}

end


%-----------------------------------------------------------------------
function varargout=vbm_plot_SDbar(mdata,sdata,names,sortdata)
% ______________________________________________________________________
% plot bars with std
%
% vbm_plot_SDbar(vbm_stat_nanmean(gmdist,1),vbm_stat_nanstd(gmdist,1),pnames,1); 
%    
% ______________________________________________________________________
% $Id$

  if exist('sortdata','var') && sortdata
    [mdata,sorti] = sort(mdata);
    sdata = sdata(sorti);
    names = names(sorti);
  end

  nc = numel(names);
  whisker_x = ones(2,1)*[1:nc,1:nc];
  whisker_y = zeros(2,2*nc);
  whisker_y(1,:) = reshape(repmat(mdata + sdata,1,2),1,nc*2);
  whisker_y(2,:) = reshape(repmat(mdata,1,2),1,nc*2);
  
  chop = find(isnan(sdata')==1); nc = numel(names);
  whisker_x(:,[chop,chop+nc]) = [];
  whisker_y(:,[chop,chop+nc]) = [];
  names(chop) = [];

  % Add caps to the remaining whiskers
  cap_x = whisker_x;
  cap_x(1,:) = cap_x(1,:) - 0.1;
  cap_x(2,:) = cap_x(2,:) + 0.1;
  cap_y = whisker_y([1,1],:);
  
  h{1}=bar(mdata,'b'); hold on;
  h{2}=plot(whisker_x, whisker_y, 'b-');
  h{2}=plot(cap_x, cap_y, 'b-');
  
  set(gca,'XTick',1:numel(names),'XTickLabel',names);
  
  if nargout>0
    varargout{1}=h;
  end
  
  hold off;
end
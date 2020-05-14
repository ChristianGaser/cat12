function hRes = cat_stat_results
% just a batch used to create cat_stat_spm_results_ui 
% delete 202006 or later



  %% call modified SPM statistic 
  if 0%exist('xSPM','var');
    xSPM.thresDesc = 'FWE';
    [hReg,xSPM] = cat_stat_spm_results_ui('Setup',xSPM); 
  else
    [hReg,xSPM] = cat_stat_spm_results_ui; 
  end
  
  % create SPM result table
  %spm_list('List',xSPM,hReg); 
  
  
  %% SPM figure
  hRes.Fgraph       = spm_figure('GetWin','Graphics');
  hRes.FgraphC      = get( hRes.Fgraph ,'children'); 
  hRes.Fmenu3d      = findobj( hRes.FgraphC ,'Type','uicontextmenu','Tag','');
  hRes.FmenuTab     = findobj( hRes.FgraphC ,'Type','uicontextmenu','Tag','TabDat');
  hRes.FgraphAx     = findobj( hRes.FgraphC,'Type','Axes');
  hRes.FgraphAxPos  = cell2mat(get( hRes.FgraphAx , 'Position'));
  hRes.Ftext        = findobj(hRes.Fgraph,'Type','Text');
  
  % find the SPM string within the surface axis 
  stext = get(hRes.Ftext,'String'); 
  hRes.Ftext3dspm   = findobj(hRes.Fgraph,'Type','Text','String', ...
    stext{ find(~cellfun('isempty',strfind(stext,'SPM\{'))) } );
  set(hRes.Ftext3dspm,'visible','off','HitTest','off');
  
  % find the SPM result text
  hRes.Ftext3dres = get( findobj(hRes.Fgraph,'Type','Text','Color',[0.7 0.7 0.7]),'parent'); 
  for axi=1:numel(hRes.Ftext3dres), set( hRes.Ftext3dres{axi},'HitTest','off'); end
  
  % fine red lines of the SPM result table
  hRes.Fline        = findobj(hRes.Fgraph,'Type','Line','Tag','');
  hRes.Flinewhite   = findobj(hRes.Fgraph,'Color',[1 1 1]); %findobj(hRes.Fgraph,'Type','Line','Color',[1 1 1]);
  hRes.Fsurf        = findobj(hRes.Fgraph,'Type','Line','Tag','CrossBar');
  set( hRes.Fsurf , 'LineWidth', 2, 'MarkerSize', 200, 'Color', [1 0 0])
  
  hRes.FlineAx      = get(hRes.Fline,'parent');
  hRes.Flabels      = [ hRes.FgraphAx( hRes.FgraphAxPos(:,1) == 0.65); hRes.FgraphAx( hRes.FgraphAxPos(:,1) == 0.02)]; 
  
  % make nice contrast box that is a bit larger than the orinal boxes
  hRes.Fcons        = hRes.FgraphAx( hRes.FgraphAxPos(:,1) == 0.65 & hRes.FgraphAxPos(:,2) > 0.6 ) ;  
  for axi = 1:numel( hRes.Fcons ), l = get( hRes.Fcons(axi) , 'ylim'); set(hRes.Fcons(axi) , 'box','on','ylim', round(l) + [-0.015 0.015]); end
  
  % remove non integer values
  hRes.Fdesm        = hRes.FgraphAx( hRes.FgraphAxPos(:,1) == 0.65 & hRes.FgraphAxPos(:,2) < 0.6 ) ;  
  xt = get(hRes.Fdesm,'xtick'); xt(round(xt)~=xt) = []; set(hRes.Fdesm,'xtick',xt); 
  
  
  %
  hRes.Fval         = hRes.FgraphAx( hRes.FgraphAxPos(:,1) > 10); 
  hRes.Fsurf        = hRes.FgraphAx( hRes.FgraphAxPos(:,1) == 0.05);
  %hRes.Flabels      = setdiff( hRes.FgraphAx , hRes.Fval ); 
  set(hRes.Fline,'HitTest','off')
  for axi = 1:numel( hRes.Flabels ), set(hRes.Flabels(axi),'HitTest','off'); end
  for axi = 1:numel( hRes.FlineAx ), set(hRes.FlineAx{axi},'HitTest','off'); end
 % rotate3d(hRes.FmenuTab,'off') 
 
 
 
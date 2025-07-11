function cat_tst_BWP(Pmethod,Presdir)
% 
% cat_tst_BWP(Pmethod,Presdir)
%
% See cat_tst_main.
%

% TODO: 
% * separate evaluation of 1mm case and all 5 resolution cases 
%   > overlap measures? stat tests?
%
% - addapt fontsize with more figures
% - single plot of each subfigure?
% - overlapping plot not really helpful for mor than 6 cases 
% - overlapping plot miss label
% - missed data -> missed method (add print)

%#ok<*SAGROW
Presdir = fullfile(Presdir,'BWP'); 
  
% find data
datas = {'thickness','pbt'};
tcase = {'100mm','mm'}; 
for datasi = 1:2
  Pdata = cell(size(Pmethod,1),1); 
  Pexist = zeros(1,size(Pmethod,1),1);
  for pi = 1:size(Pmethod,1)
    if exist( fullfile(Pmethod{pi,2},'01_Bestoff') ,'dir')
      Pdata{pi} = char( cat_vol_findfiles( fullfile(Pmethod{pi,2},'01_Bestoff') , sprintf( 'lh.%s.BWPT*100mm',datas{datasi}) ));
    elseif exist( Pmethod{pi,2} ,'dir')
      Pdata{pi} = char( cat_vol_findfiles( Pmethod{pi,2} , sprintf( 'lh.%s.BWPT*100mm',datas{datasi}) )) ;
    else
      error('ERROR: Miss startfolder: %s\n',Pmethod{pi,2} );
    end
    Pexist(pi) = ~isempty(Pdata{pi});
    if ~Pexist(pi), cat_io_cprintf('err',sprintf('Miss %s/%s\n',Pmethod{pi,2}, sprintf( 'lh.%s.BWPT*100mm',datas{datasi}))); end
  end
  Pdata(Pexist==0) = []; 
  Pmethod(Pexist==0,:) = []; 
  
  % plot histogram
  cph = cat_plot_histogram(char(Pdata), struct('color',cell2mat(Pmethod(:,3))));
  yrange = round( get(get(cph(1),'Parent'),'YLim')/2,2)*3; 
  %  
  for fi = 1:2
    fh = figure(394); clf(fh)
    if fi==1 
      subplots = size(Pmethod,1); 
      xt = ceil(size(Pmethod,1) / 4); 
      tiledlayout(xt, ceil(size(Pmethod,1) /  xt),   ...
        'TileSpacing', 'compact', 'Padding', 'compact'); 
    else
      subplots = 1; 
    end
  
    for si = 1:subplots
      if fi==1
        % single print on one figure
        ax    = nexttile; hold on; 
        spi   = si;
        copyobj( cph(si), ax); 
      else
        ax    = axes; 
        spi   = 1:size(Pmethod,1); 
        copyobj( allchild( get(cph(1),'parent') ), ax); 
      end
  
     
      %%
      xlim(ax,[1,3]); ylim(ax,yrange); 
      ax.XTick = 1:.5:3; 
      grid(ax,'on'); box(ax,'on'); 
      xlabel('thickness (mm)'); 
      if fi==1
        title(ax,sprintf('BWPT-%s',Pmethod{si,1}));
        md = zeros(1,3); std = md; pk = md;
      else
        title(ax,'Brain Web Phantom Thickness'); 
        md = zeros(size(Pmethod,1),3); std = md; pk = md;
      end
  
  
      %% estimate some numbers
      for pi = spi
        cdata{pi}  = [
          cat_io_FreeSurfer('read_surf_data',Pdata{pi}); 
          cat_io_FreeSurfer('read_surf_data',strrep(Pdata{pi},[filesep 'lh'],[filesep 'rh']))]; %#ok<AGROW>
      
        tval = 1.5:0.5:2.5; d=cell(1,3); xp = zeros(1,3);
        for di = 1:3
          range = [tval(di)-0.2 tval(di)+0.2];
          d{di} = smooth( cph(pi).YData( cph(pi).XData(:)>range(1) & cph(pi).XData(:)<range(2) ) ); 
          spm_smooth(d{di},d{di},5); % smooth also beyond borders
          [md(pi,:),std(pi,:)] = cat_stat_kmeans( cdata{pi}(cdata{pi}>0.5 & cdata{pi}<2.9) ,3); 
          [pk(pi,di), xp(di)] = max( d{di} ); 
          pkv(pi,di) = cph(pi).XData(xp(di) + nnz( cph(pi).XData(:)<range(1)) ); %#ok<AGROW>
          pkx(pi,di) = abs(pkv(pi,di) - tval(di)); %#ok<AGROW>
          px(pi,di)  = cph(pi).XData(xp(di) + nnz( cph(pi).XData(:)<range(1) )); %#ok<AGROW>
          if fi==1
            plot(ax, repmat( cph(pi).XData(xp(di) + nnz( cph(pi).XData(:)<range(1) )), 1, 2), ...
                          [0 cph(pi).YData(xp(di) + nnz( cph(pi).XData(:)<range(1) ))], ...
                          'Color', Pmethod{spi,3}, 'LineWidth', 1  ); 
          end
        end
        subtitle(sprintf('MPH=%0.4f, MAPE=%0.4f',mean(pk(pi,:)),mean(pkx(pi,:))))
      end
      pk = pk / max(pk(:)); 
    end
  
    ff = sprintf('cat_tst_BWP_plot%d_%s',fi,datas{datasi});
    if ~exist(Presdir,'dir'), mkdir(Presdir); end
    print(fh,fullfile(Presdir,ff),'-r300','-dpng');

    if fi==2
      title(sprintf('Rusak 2021 (N=%d)',numel(Pdata{pi})));
      subtitle('[MPH=mean peak height (higher=better), MAPE=mean absolute peak error (lower=better)]'); 
      try %######## not working
        legend(  strcat( Pmethod(:,1), ...
          repmat(' (MPH=',5,1), num2str(mean(pk,2),'%0.3f'), ...
          repmat(' MAPE=',5,1), num2str(mean(pkx,2),'%0.3f'), ...
          repmat(')',size(pk,1),1) ));
      end
      
      ff = sprintf('cat_tst_BWP_%s',datas{datasi});
      tab = [{'method',     'peaks','','',            'NMPH=normaized mean peak height (higher=better)','','','',   'MAPE=mean absolute peak error (lower=better)','','',''; 
              ''      ,     '1.5mm','2.0mm','2.5mm',  'mean','1.5mm','2.0mm','2.5mm',                    'mean','1.5mm','2.0mm','2.5mm'}
              Pmethod(:,1),  cellstr(num2str(px(:,1),'%0.4f')),  cellstr(num2str(px(:,2),'%0.4f')),  cellstr(num2str(px(:,3),'%0.4f')), cellstr(num2str(mean(pk,2),'%0.4f')), ...
                             cellstr(num2str(pk(:,1),'%0.4f')),  cellstr(num2str(pk(:,2),'%0.4f')),  cellstr(num2str(pk(:,3),'%0.4f')),   ...
                             cellstr(num2str(mean(pkx,2),'%0.4f')), cellstr(num2str(pkx(:,1),'%0.4f')), cellstr(num2str(pkx(:,2),'%0.4f')), cellstr(num2str(pkx(:,3),'%0.4f'))  ];
      cat_io_csv(fullfile(Presdir, sprintf('%s.csv',ff)),tab); 
  
      cell2table( tab )
    end
  end
  try %#ok<TRYNC>
    close(get(get(cph(1),'parent'),'parent')); 
  end
end

%%

close(fh); 


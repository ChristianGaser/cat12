function cat_tst_BWP(Pmethod,Presdir0,datas,fast)
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
% - overlapping plot not really helpful for mor than 6 cases > skip
% - overlapping plot miss labels
% - missed data -> missed method (add print)

%#ok<*SAGROW
Presdir = fullfile(Presdir0,'BWP'); 
if ~exist(Presdir,'dir'), mkdir(Presdir); end
if ~exist('datas','var'), datas = {'thickness','pbt'}; end
     
% find surface data (texture) for analysis
tcase = {'100mm','mm'}; fasti = 2-fast;  
for datasi = 1:numel(datas)
  % get files from the Bestoff or Phantom (multiple resolutions) directoy 
  Pdata = cell(size(Pmethod,1),1); 
  for pi = 1:size(Pmethod,1)
    if fasti == 2  &&  exist( fullfile(Pmethod{pi,2},'07_Phantoms') ,'dir') 
      Pdata{pi} = cat_vol_findfiles( fullfile(Pmethod{pi,2},'07_Phantoms') , sprintf( 'lh.%s.BWPT*%s' ,datas{datasi},tcase{fasti}) );
    elseif exist( fullfile(Pmethod{pi,2},'01_Bestoff') ,'dir') 
      Pdata{pi} = cat_vol_findfiles( fullfile(Pmethod{pi,2},'01_Bestoff') , sprintf( 'lh.%s.BWPT*%s' ,datas{datasi},tcase{fasti}) );
    elseif exist( Pmethod{pi,2} ,'dir')
      Pdata{pi} = cat_vol_findfiles( Pmethod{pi,2} , sprintf( 'lh.%s.BWPT*%s',datas{datasi},tcase{fasti}) );
    else
      error('ERROR: Miss startfolder: %s\n',Pmethod{pi,2} );
    end
  end
  

  %% plot histogram
  if fasti == 1
  % just one case
    Pdata = cellfun(@(x) char(x), Pdata, 'UniformOutput', false); 
    cph = cat_plot_histogram(char(Pdata), struct('color',cell2mat(Pmethod(:,3))));  
    yrange = round( get(get(cph(1),'Parent'),'YLim')*1.5 / 2,2)*2; 
  else
  % extend version for multiple intput datasets per method (resolution levels)

    % get main plot to set up range
    Pdatatmp = [Pdata{:}]; Pdatatmp = Pdatatmp(1:6:end); 
    cph = cat_plot_histogram(char(Pdatatmp), struct('color',cell2mat(Pmethod(:,3))));  
    ax0 = cph.Parent; fg0 = ax0(1).Parent; fg0.Visible = false;
    %yrange = round( get(get(cph(1),'Parent'),'YLim')/2,2); 
    yrange = [0. .03]; 
  
    for pi = 1:size(Pmethod,1)
      %% get basic plot
      cphpi{pi} = cat_plot_histogram(char(Pdata{pi}), struct('color',jet(numel(Pdata{pi}))));
      ax1 = cphpi{pi}.Parent; fg1 = ax1(1).Parent; fg1.Visible = false;

      % optimize plot
      orgax = get(cphpi{pi}(1),'parent'); 
      delete(findobj(orgax.Children,'DisplayName',''))
      orgax.FontSize = 13;
      hdata   = reshape( [cphpi{pi}(:).YData], numel( cphpi{pi}(1).YData ), numel(cphpi{pi}) )'; 
      hdatamn = mean( hdata ); 
      hold on
      pmn = plot( cphpi{pi}(1).XData , hdatamn,'-'); pmn.Color = [0 0 0]; pmn.LineWidth = 2;
      orgfig  = get(orgax(1),'parent'); 
      xlim(orgax,[0,4]); ylim(orgax,yrange); 
      orgax.XTick = .5:.5:3.5; 
      grid(orgax,'on'); box(orgax,'on'); 
      ylabel(datas{datasi}); 
      xlabel('thickness (mm)'); 
      title(orgax,Pmethod(pi,1));
      legend({'0.7 mm','1.00 mm','1.25 mm','1.50 mm','1.75 mm','2.00 mm','average'});
      ff = sprintf('cat_tst_BWP_multires%d_plot3_%s_%s',fasti-1,datas{datasi},Pmethod{pi,1});
      print(orgfig(1),fullfile(Presdir,ff),'-r300','-dpng');

      cph(pi).XData = cphpi{pi}(:).XData;
      cph(pi).YData = mean( hdata );
      cph(pi).Color = Pmethod{pi,3};

    end

  end  



  %% create two figures:
  %  fi==1 with subplots for all figures
  %  fi==2 with all methods in one figure (not so usefull)
  fi=1; si = 1; %#ok<NASGU>
  for fi = 1 % :2
    %%
    fh = figure(394); clf(fh); set(fh,'Visible',false); 
    if fi==1 
      subplots = size(Pmethod,1); 
      xt = ceil(size(Pmethod,1) .^ .5 * 2/3); 
      if size(Pmethod,1)>6 & size(Pmethod,1)<13, xt = 2; end
      tiledlayout(xt, ceil(size(Pmethod,1) /  xt),   ...
        'TileSpacing', 'compact', 'Padding', 'compact'); 
    else
      subplots = 1; 
    end
    fh.Position(3:4) = [200   300] .* [ceil(size(Pmethod,1) /  xt),xt]; 
    %%
    for si = 1:subplots
      if fi==1
        ax    = nexttile; hold on; 
        spi   = si;
        copyobj( cphpi{si}, ax); 
        lh = ax.Children; 
        mh = copyobj( cph(si), ax); 
        set(lh,'Color',min(1,.85 + cph(si).Color*0.15))
        set(mh,'LineWidth',1.5); 

        close(get(get(cphpi{si}(1),'parent'),'parent')); 
      else
        % single print on one figure
        ax    = axes; 
        spi   = 1:size(Pmethod,1); 
        orgax = allchild( get(cph(1),'parent') ); 
        copyobj( orgax(numel(orgax)/2+1:numel(orgax)) , ax); 
        ax.FontSize = 13;
        % add marker for highest value
        hold on; 
        axC = flip( ax.Children); % reverse order
        for axci = 1:numel(axC)
          [maxy,maxxi] = max(axC(axci).YData); 
          maxx =  axC(axci).XData(maxxi);
          plx = plot(maxx,maxy); 
          plx.Color  = Pmethod{axci,3}; 
          plx.Marker = Pmethod{axci,4}; 
          plx.MarkerEdgeColor = plx.Color;
          plx.MarkerFaceColor = plx.Color;
        end

      end
  
     
      %%
      xlim(ax,[1,3]); ylim(ax,yrange); 
      ax.XTick = 1:.5:3; 
      grid(ax,'on'); box(ax,'on'); 
      ylabel(datas{datasi}); 
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
        if fasti == 1
          cdata{pi}  = [
            cat_io_FreeSurfer('read_surf_data',Pdata{pi}); 
            cat_io_FreeSurfer('read_surf_data',strrep(Pdata{pi},[filesep 'lh'],[filesep 'rh']))]; %#ok<AGROW>
        else
          if 0
            cdata{pi} = cph(pi).YData;
          else
            cdata{pi} = []; 
            for xi = 1:numel(Pdata{pi})
              cdata{pi} = [ cdata{pi} ; 
                cat_io_FreeSurfer('read_surf_data',Pdata{pi}{xi}); 
                cat_io_FreeSurfer('read_surf_data',strrep(Pdata{pi}{xi},[filesep 'lh'],[filesep 'rh']))]; 
            end
          end
        end
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
  
   
    if fi==2
      title(sprintf('BWP (N=%d)', numel(Pdata) ));
      subtitle('[MPH=mean peak height (higher=better), MAPE=mean absolute peak error (lower=better)]'); 
      try 
        legend(  strcat( Pmethod(:,1), ...
          repmat(' (MD=',size(Pmethod,1),1), num2str(md(:,2),'%0.3f'), ...
          repmat(', MPH=',size(Pmethod,1),1), num2str(mean(pk,2),'%0.3f'), ...
          repmat(', MAPE=',size(Pmethod,1),1), num2str(mean(pkx,2),'%0.3f'), ...
          repmat(')',size(pk,1),1) ));
      end
      
      tab = [{'method',     'peaks','','',            'NMPH=normaized mean peak height (higher=better)','','','',   'MAPE=mean absolute peak error (lower=better)','','',''; 
              ''      ,     '1.5mm','2.0mm','2.5mm',  'mean','1.5mm','2.0mm','2.5mm',                    'mean','1.5mm','2.0mm','2.5mm'}
              Pmethod(:,1),  cellstr(num2str(px(:,1),'%0.4f')),  cellstr(num2str(px(:,2),'%0.4f')),  cellstr(num2str(px(:,3),'%0.4f')), cellstr(num2str(mean(pk,2),'%0.4f')), ...
                             cellstr(num2str(pk(:,1),'%0.4f')),  cellstr(num2str(pk(:,2),'%0.4f')),  cellstr(num2str(pk(:,3),'%0.4f')),   ...
                             cellstr(num2str(mean(pkx,2),'%0.4f')), cellstr(num2str(pkx(:,1),'%0.4f')), cellstr(num2str(pkx(:,2),'%0.4f')), cellstr(num2str(pkx(:,3),'%0.4f'))  ];
      ff = sprintf('cat_tst_BWP_multires%d_%s',fasti-1,datas{datasi});
      cat_io_csv(fullfile(Presdir, sprintf('%s.csv',ff)),tab); 
  
      % print table
      tab(1,:) = {'method','peak','peak','peak','NMPH','NMPH','NMPH','NMPH',   'MAPE','MAPE','MAPE','MAPE'}; 
      cell2table( tab(3:end,:) ,'VariableNames',strcat( strcat( tab(1,:)',tab(2,:)')') )
    end

    ff = sprintf('cat_tst_BWP_multires%d_plot%d_%s',fasti-1,fi,datas{datasi});
    print(fh,fullfile(Presdir,ff),'-r300','-dpng');

  end
  try %#ok<TRYNC>
    close(get(get(cph(1),'parent'),'parent')); 
  end
end

%%

close(fh); 


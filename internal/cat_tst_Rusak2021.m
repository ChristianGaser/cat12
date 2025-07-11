function cat_tst_Rusak2021(Pmethod0,Presdir0)
% 
% TODO:
% - add local change maps ?
% - add other measures (rGMV, intensity, ...)
% - add SPM regression analysis 
% - eval all values by round?

%#ok<*SAGROW

% data path definition
if exist('Pmethod0','var')
  Pmethod = Pmethod0;
  Presdir = Presdir0;
else
  Praw = '/Volumes/WDE18TB/MRDataPP/202505_CS4';
  Pmethod = {
   ...'C12-CS22'    fullfile(Praw,'derivatives','CAT12.9_2679_CS2')                     [.0 .8 .0] 'd'
  ... 'C12-CS24'    fullfile(Praw,'derivatives','CAT12.9_2679_CS24')                    [.3 .6 .0] 'p'
   'C12-CS24R5'  fullfile(Praw,'derivatives','CAT12.9_2679_CS24')                    [.3 .6 .0] 'p'
  ... 'C12-CS41'    fullfile(Praw,'derivatives','CAT12.9_2679_CS41')                    [.8 .0 .6] '<'
  ... 'C12-CS42'    fullfile(Praw,'derivatives','CAT12.9_2679_CS42')                    [.6 .0 .8] '>' % miss last subject
  ... 'C12-CS42p0'  fullfile(Praw,'derivatives','CAT12.9_2701_CS42p0')                  [.6 .0 .8] '>' % no big difference to Ym-Ybased
   'C12-CS42R5'  fullfile(Praw,'derivatives','CAT12.9_2701_CS42R5')                  [.8 .0 .9] '>' % no big difference to Ym-Ybased
   'S25-CS22'    fullfile(Praw,'derivatives','CAT12.9_2701_SPMorg_CS22','SPM25')     [.8 .0 .0] 'd'
  ... 'S25-CS24'    fullfile(Praw,'derivatives','CAT12.9_2701_SPMorg_CS24','SPM25')     [.8 .0 .0] 'd'
  ... 'S25-CS24A'   fullfile(Praw,'derivatives','CAT12.9_2701_SPMamap_CS24','SPM25')    [.4 .0 .0] 'p' % miss subdir
   ...'S25-CS42R3'  fullfile(Praw,'derivatives','CAT12.9_2701_SPM_CS42R3','SPM25')    [.4 .0 .0] 'p' % miss subdir Russak failed
   ...'T1Prep'      fullfile(Praw,'derivatives','T1Prep')                               [.0 .1 .9] '^'
   ...'T1PrepA'     fullfile(Praw,'derivatives','T1PrepAMAP')                           [.0 .5 .3] 'd' % not processed
   ...'T1PrepA1'    fullfile(Praw,'derivatives','T1PrepAMAPds1')                        [.0 .3 .6] 'v'
   ...'T1Prep'     '/Users/dahnke/Documents/T1Prep/test_folder'                                           [.0 .4 .8] '^'
 ...  'T1PrepA'     '/Users/dahnke/Documents/T1Prep/test_folder_amap'                                      [.0 .6 .4] 'v'
  };
  Presdir = '/Volumes/WDE18TB/MRDataPP/202505_CS4/derivatives/results'; 
end


% add Rusak subdir 
for pi = 1:size(Pmethod,1)
  if any(contains(Pmethod{pi,1},'T1Prep'))
    Pmethod(pi,2) = fullfile(Pmethod(pi,2), 'surf'); 
  else
    Pmethod(pi,2) = fullfile(Pmethod(pi,2), '16_Rusak2021n36'); 
  end
  if ~exist( Pmethod{pi,2} ,'dir'), error('Miss directory %s',Pmethod{pi,2}); end
end
Presdir = fullfile(Presdir,'Rusak2021'); 


% load data
fprintf('Load data: '); 
mnthick = nan(1,size(Pmethod,1));
mdthick = nan(1,size(Pmethod,1));
Psurf   = cell(1,size(Pmethod,1));
Pdata   = cell(1,size(Pmethod,1));
datas   = 'thickness'; 
for pi = 1:size(Pmethod,1)
  fprintf('\n  Method %2d/%2d: %16s', pi, size(Pmethod,1), Pmethod{pi} )
  if pi == 1 
    Psurf{pi}  = cat_vol_findfiles(Pmethod{pi,2},'lh.central.sub-ADNI*simGMatrophy*');
    Pdata{pi}  = cat_vol_findfiles(Pmethod{pi,2},sprintf('lh.%s.sub-ADNI*simGMatrophy*',datas));
    Pexist     = zeros(numel(Pdata{pi}), size(Pmethod,1)); 
  else 
    Psurf{pi}  = cat_io_strrep(Psurf{1},  Pmethod{1,2}, Pmethod{pi,2} );
    Pdata{pi}  = cat_io_strrep(Pdata{1}, Pmethod{1,2}, Pmethod{pi,2} );
  end
  % find missing files
  Pexist(:,pi) = cellfun(@(x) exist(x,'file'),Pdata{pi});
end


% remove missing files in all testsets
if any(Pexist(:)==0)
  cat_io_cprintf('err','\n\nMissing files:'); 
  for pi = 1:size(Pmethod,1)
    if all(Pexist(:,pi)==0)
      cat_io_cprintf('err','Miss all files in %s',Pmethod{pi,2})
    elseif ~all(Pexist(:,pi)==2)
      Psurf{pi}( Pexist(:,pi)==0 )
    end
  end
else
  cat_io_cprintf([0 0.5 0],'\n\nNo missing files!\n'); 
end
% removal not use
%for pi = 1:size(Pmethod,1)
%  Psurf{pi}( sum(Pexist==0,2)>0 ) = [];
%  Pfiles{pi}( sum(Pexist==0,2)>0 ) = [];
%end
fprintf('\nLoad data 2: '); 
S = cell(size(Pmethod,1),numel(Pdata{pi})); cdata = S; 
for pi = 1:size(Pmethod,1)
  fprintf('\n  Method %d/%d: %16s', pi, size(Pmethod,1), Pmethod{pi} )
  for si = 1:numel(Pdata{pi})
    txt = sprintf('    Load subject %3d/%3d', si,numel(Pdata{pi}));
    if si > 1, fprintf(repmat('\b',1,numel(txt))); end
    fprintf(txt); 

    try
      S{pi,si}       = gifti(Psurf{pi}{si}); 
      cdata{pi,si}   = cat_io_FreeSurfer('read_surf_data',Pdata{pi}{si}); 
    catch
      S{pi,si}       = struct(); 
      cdata{pi,si}   = nan; 
    end
    mnthick(si,pi) = cat_stat_nanmean(cdata{pi,si}); 
    mdthick(si,pi) = cat_stat_nanmedian(cdata{pi,si}); 
  end
end


% get basic information from filename
ff          = spm_str_manip(Pdata{pi},'t');
ffsubid     = strfind(ff{1},'sub-ADNI');
subject     = cellfun(@(ffc) ffc(ffsubid+8:ffsubid+8+7),ff,'UniformOutput',false); 
subid = zeros(size(subject)); sd = unique(subject); for si = 1:numel(sd), subid( contains(subject,sd{si}) ) = si; end
atrophy     = cellfun(@(ffc) ffc(end-5:end-2),ff,'UniformOutput',false); 
atrophynum  = cellfun( @str2num, atrophy); 
age         = 20 + 80*atrophynum; 
sex         = zeros(size(atrophynum)); 
HC          = ones(size(atrophynum)); 
tiv         = 1000 * ones(size(atrophynum,1),size(Pmethod,1)); 
% #### get real values from ADNI files or is there a CSV? #######
scan        = strcat(subject,'_',atrophy); 
fprintf('\nLoad data done. \n\n'); 


%% local mapping
%% ??? 
% - what with CC? masking?
fprintf('Map data: \n');
rmse = @(x) mean(x.^2).^.5;

for pi = 1:size(Pmethod,1)
  fprintf('  Method %d/%d: %16s\n', pi, Pmethod{pi}, size(Pmethod,1))

  for si = 1:numel(subject) 
    fprintf('    Load subject %d - %s\n', si, scan{si})
    if atrophynum(si) > 0
      cdatam{pi,si} = cat_surf_fun('cdatamapping', ...
        S{pi,contains(scan,[subject{si} '_0.00'])},  S{pi,si},  cdata{pi,si}); 
    else
      cdatam{pi,si} = cdata{pi,si}; 
    end
    try
      cdatadiff(si,pi) = rmse( cdatam{pi,si} - cat_stat_nanmean( [ ...
        (cdatam{pi,contains(scan,[subject{si} '_0.00'])} + 0.00), ...
        (cdatam{pi,contains(scan,[subject{si} '_0.02'])} + 0.02), ...
        (cdatam{pi,contains(scan,[subject{si} '_0.05'])} + 0.05) ] ) ); 
    catch
      cdatadiff(si,pi) = rmse( cdatam{pi,si} - cdatam{pi,contains(scan,[subject{si} '_0.00'])} ); 
    end
  end
end
%% this is not working
if 1
  clear cdatasub cdatamed0 cdatadiff
  for pi = 1:size(Pmethod,1)
    for si = 1:numel(subject) 
      cdatamed0{pi,si} = cdatam{pi,si}; % + atrophynum(si); 
    end
    for ssi = 1:max(subid(:))
      cdatasub{pi,ssi} = median( [ cdatamed0{pi,subid==ssi}] ,2 );
    end
    %
   % for si = 1:numel(subject) 
   %   cdatadiff(pi,si) = rmse( cdatam{pi,si} -  cdatasub{pi,subid(si)} ); 
   % end
  end
end

%% create statistics ... as (i) cross-sectional and (ii) as longitudinal analysis
localstat = 0; 
onlyAllreadySmoothed = 0; 
if localstat
  for subset = 4 %1:3 
    for pi=1:size(Pmethod,1)
      % VBM 
      %matlabbatch = cat_tst_agereg( Pgmv{pi}, Presdir, Pmethod{pi,1}, 'GMV', age, sex, cohort==1, tiv(:,pi), 1, onlyAllreadySmoothed);
      %if ~isempty(matlabbatch), spm_jobman('run',matlabbatch); end
    
      % SBM
      for sm = {'pbt'} %,'thickness'} ... thickness not working - similar results ...
        %%
        nsub = numel(unique(subid) );
        Pdatapi = strrep( Pdata{pi} ,sprintf('lh.%s.',datas), sprintf('lh.%s.',sm{1})); 
        switch subset
          case 1, ss = (1:numel(age))'; Presdirss = fullfile( Presdir, '6TP_1.0mm'); 
          case 2, ss = repmat( [1 0 0 0 1 1], 1, nsub)'>0; Presdirss = fullfile(Presdir, '3TP_1.0mm');  % three time points with 0 0.5 and 1 mm loss
          case 3, ss = repmat( [1 0 1 1 0 0], 1, nsub)'>0; Presdirss = fullfile(Presdir, '3TP_0.1mm');  % three time points with 0 0.05 and .1 mm loss
          case 4, ss = repmat( [1 0 0 1 0 0], 1, nsub)'>0; Presdirss = fullfile(Presdir, '2TP_0.5mm');  % two time points with .5 mm loss
        end
        matlabbatch = cat_tst_agereg( Pdatapi(ss), [Presdirss '_long'], Pmethod{pi,1}, sm{1}, age(ss), sex(ss), HC(ss), tiv(ss,pi), 1, onlyAllreadySmoothed, subid(ss));
        if ~isempty(matlabbatch) 
          try
            spm_jobman('run',matlabbatch); 
          catch 
            cat_io_cprintf('err',sprintf('Error running matlabbatch "%s"',Presdirss)); 
          end  
        end
        if subset == 1
          matlabbatch = cat_tst_agereg( Pdatapi(ss), [Presdirss '_cs'], Pmethod{pi,1}, sm{1}, age(ss), sex(ss), HC(ss), tiv(ss,pi), 1, onlyAllreadySmoothed);
          if ~isempty(matlabbatch) 
            try
              spm_jobman('run',matlabbatch); 
            catch 
              cat_io_cprintf('err',sprintf('Error running matlabbatch "%s"',Presdirss)); 
            end  
          end
        end
      end
    end
  end
end


%% create figure
mnthickloss = zeros(numel(subject), size(Pmethod,1)); 
mdthickloss = mnthickloss;
fh = figure(33);  % set(fh,'Visible',0,'Interruptible',0); 
for fi = 1:2
  %%
  fprintf('Analyse Rusak data: \n'); 
  acor   = zeros(1,size(Pmethod,1)); afit = cell(size(acor)); 
  afitp1 = zeros(1,size(Pmethod,1));
  clf, hold on; 
  fh.Position(3:4) = [600*fi 800]; 

  % print for legend
  if fi==1
    for pi = 1:size(Pmethod,1)
      plot(0,0,'Color',Pmethod{pi,3},'Marker', ...
        Pmethod{pi,4},'MarkerFaceColor',Pmethod{pi,3});
    end
  else
    xt = ceil(size(Pmethod,1) / 5); 
    tiledlayout(xt, ceil(size(Pmethod,1) /  xt),   ...
      'TileSpacing', 'compact', 'Padding', 'compact'); 
  end


  % main print
  for pi = 1:size(Pmethod,1)
    if fi==2, fhagex = nexttile; hold on; end
  
    for si = 1:numel(subject) 
      try
        if 0
          mnthickloss(si,pi) = mnthick(si,pi) - mnthick( contains(scan,[subject{si} '_0.00']), pi );
          mdthickloss(si,pi) = mdthick(si,pi) - mdthick( contains(scan,[subject{si} '_0.00']), pi );
        else
          mnthickloss(si,pi) = mnthick(si,pi) - mean(cdatasub{pi,subid(si)});
          mdthickloss(si,pi) = mdthick(si,pi) - median(cdatasub{pi,subid(si)});
        end
      catch
        % no 0.00 case
        mnthickloss(si,pi) = nan;
        mdthickloss(si,pi) = nan;
      end
    end
    subs = unique(subject); 
    if 1
      for si = 1:numel(subs)
        subi         = contains(scan,subs{si}); 
        ph           = plot( atrophynum(subi) , mnthickloss(subi, pi) ); 
        ph.Color     = max(0,min(1,Pmethod{pi,3}.^.5));
        ph.LineWidth = .25;
        ph.Marker    = Pmethod{pi,4};
        ph.LineStyle = '--';
        ph.MarkerFaceColor = ph.Color;
      end
    end
  
    % plot fit (only linear tissue loss!)
    nonan = ~isnan(atrophynum) & ~isnan(mnthickloss(:,pi)); 
    if any(nonan)
      [afit{pi},bfit{pi}]  = fit( atrophynum(nonan) , mnthickloss(nonan,pi), 'poly1');
      afitp1(pi) = afit{pi}.p1; 
      arsq(pi)   = bfit{pi}.rsquare; 
    else
      afit{pi}   = []; 
      afitp1(pi) = nan; 
      arsq(pi)   = nan; 
    end
    if ~isempty(afit{pi})
      acor(pi)   = corr( atrophynum(nonan) , mnthickloss(nonan,pi) );  
    else
      acor(pi)   = nan; 
    end

    % plot confidence interval ... ah not really happy with this one
    xint = linspace(0, 90,90);
    cover = 0.05:.2:.95;
    for ci = 1:numel(cover)
      if ~isempty(afit{pi})
        CIO = predint(afit{pi},xint,cover(ci),'obs');
        p = fill([xint'; flip(xint')],[CIO(:,1); flip(CIO(:,2))], Pmethod{pi,3}); 
        p.FaceAlpha = .1; p.EdgeAlpha = .0; p.EdgeColor = Pmethod{pi,3};
      end
    end
  
    if ~isempty(afit{pi})
      fph = plot( afit{pi});
      fph.Color     = Pmethod{pi,3};
      fph.LineWidth = 2;
      fph.MarkerFaceColor = fph.Color;
    end

    % format
    grid on; box on; 
    ylim([-1.2,0.1]);
    xlabel('simulated atrophy (mm)'); 
    ylabel('estimated atrophy (mm)');
    
    xlim([0,1]);
    ylim([-1.1 .1]); 
    plot([0 1],[0 -1],'Color',[0.5 .5 .5]); 
    set(gca,'FontSize',16,'XTick',0:.2:1); 
    if fi==2
      legend off; 
      title(sprintf('%s',Pmethod{pi,1}));
      st = subtitle( strcat( strcat('n=', num2str(nnz(Pexist(:,pi)==2),'%0.0f'), ...
        ', r=', num2str(acor(pi)','%0.2f'), ...
        ', p1=', num2str(afitp1(pi)','%0.2f'), ...
        ', r^2=', num2str(arsq(pi),'%0.2f') ))); 
      st.FontSize = st.FontSize-2; 
    end  
  end
  if fi==1
    title('Rusak2021 atrophy phantom'); 
    lgd = legend( strcat( strcat(Pmethod(:,1), ...
      ' (n=', num2str(sum(Pexist==2)','%0.0f'), ...
      ', r=', num2str(acor','%0.4f'), ...
      ', p1=', num2str(afitp1','%0.4f'), ...
      ', r^2=', num2str(arsq','%0.4f') , ')' ))); 
  end  

  % ... save ...
  ff = sprintf('cat_tst_Rusak2021_plot%d',fi);
  if ~exist(Presdir,'dir'), mkdir(Presdir); end
  print(fh,fullfile(Presdir,ff),'-r300','-dpng');
end
%close(fh);


%% write CSV table 
tab = [{'method'},{'rho'},{'p1'},{'r^2'},{'n'},{'failed'}; 
       Pmethod(:,1), num2cell(acor'), num2cell(afitp1'), ...
        num2cell(arsq'), num2cell(sum(Pexist==0)'), num2cell(sum(Pexist==2)') ]; 
cat_io_csv(fullfile(Presdir, sprintf('%s.csv',ff)),tab); 
 
fprintf('done.\n\n'); 


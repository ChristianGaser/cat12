function cat_tst_qa_ATLAS( datadir0, qaversions, segment, fasttest ) 
%% Evaluation of CATQC in ATLAS 
%  ------------------------------------------------------------------------
%  QC evaluation of the ATLAS lesion dataset for original and lesion-masked 
%  images. In the ideal case the measures are somewhat similar and the QC 
%  measures are not or less affected by the tissue changes in lesions.
% 
%  Requirements: 
%   0. Download and install SPM and CAT
%   1. Download ATLAS T1 data from: 
%
%   2. Specify in this script: 
%      1) the data directory "datadir" 
%      2) the QC version you would like to tests (the file has to exist in the cat directory) 
%      3) the segmentation you would like to use
%
%  ------------------------------------------------------------------------

cat_io_cprintf([0 0.5 0],'\n\n== Run cat_tst_qa_ATLAS ==\n') 
  
% ### datadir ###
if ~exist( 'datadir0' , 'var' )
  datadir  = '/Volumes/SG5TB/MRData/202503_QA/ATLAS_R1.1'; 
else
  datadir  = fullfile(datadir0,'ATLAS'); 
end

% ### QC version ### 
if ~exist( 'qaversions' , 'var')
  qaversions = {
    ...'cat_vol_qa201901';  % classic version (quite stable since 2016)
    'cat_vol_qa201901x'; % refined, debugged version of 201901 
    ...'cat_vol_qa202110';  % second classic version (successor of 201901)
    ...'cat_vol_qa202110x'; % refined, debugged version of 202110   
    ...'cat_vol_qa202205';  % last regular version before update (successor of 202110, stopped)
    ...'cat_vol_qa202310';  % redesigned version based on 201901 and 202110  * default *
    ...'cat_vol_qa202412';  % experimental version with internal segmentation >> qcseg
     };
end

% ### segmention ###
if ~exist( 'segment' , 'var')
  segment = {'CAT'}; % {'SPM','CAT','qcseg'}; % qcseg requires cat_vol_qa2024012
end

if ~exist( 'fasttest', 'var'), fasttest = 0; end 
fast = {'full','fast'};
  
resdir     = fullfile(fileparts(datadir), '+results', ['ATLAS_' fast{fasttest+1} '_' datestr(clock,'YYYYmm')]); 
if ~exist(resdir,'dir'), mkdir(resdir); end

% directories
runPP      = 1; 
rerun      = 1; 

% prepare data
% - gunzip 
% - imcalc 


if runPP
  for si = 1:numel(segment)
    clear matlabbatch; 
    ATLASfiles = cat_vol_findfiles( datadir , 'ATLAS*.nii',struct('depth',3));
    switch segment{si}
      case 'CAT'
        CATpreprocessing4qc;
        matlabbatch{1}.spm.tools.cat.estwrite.data = ATLASfiles;
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.lazy = 1; 
        spm_jobman('run',matlabbatch);
      case 'SPM'
        SPMpreprocessing4qc;
        matlabbatch{1}.spm.spatial.preproc.channel.vols = ATLASfiles;
        spm_jobman('run',matlabbatch);
      case 'synthseg'
        error('synthseg is not prepared in the public script ... use MATLAB help')
      case 'qcseg'
        fprintf('No preprocessing required.\n\n');
    end
  end
end



%
qais = 1:numel(qaversions);
% (re)process QC values
datadira = fullfile(datadir,'ATLASR1.1a'); % msk
datadirb = fullfile(datadir,'ATLASR1.1b'); % org
switch segment{si}
   case 'CAT'  
     Pp0b = cat_vol_findfiles( datadirb , 'p0ATLAS*.nii',struct('depth',3)); %spm_fileparts(datadir) 
   case 'SPM'
     Pp0b = cat_vol_findfiles( datadirb , 'c1ATLAS*.nii',struct('depth',2)); %spm_fileparts(datadir) 
   case 'synthseg'
   case 'qcseg'
     Pp0b = cat_vol_findfiles( datadirb , 'p0_qcseg_ATLAS*.nii',struct('depth',2)); %spm_fileparts(datadir) 
end
Pp0a = cat_io_strrep(Pp0b,{'ATLASR1.1b';'ATLASorg'},{'ATLASR1.1a';'ATLASmsk'}); %cat_vol_findfiles( datadira , 'p0*'); %spm_fileparts(datadir) 
if fasttest % 5 ~> 50 scans
  Pp0a = Pp0a(1:12:end);
  Pp0b = Pp0b(1:12:end);
end
Pp0 = [Pp0a, Pp0b]; Pp0 = Pp0'; Pp0 = Pp0(:);  % [msk, org]
for qai = qais
   switch segment{si}
     case 'CAT',    qcv = [qaversions{qai} '_']; 
     case 'SPM',    qcv = [qaversions{qai} '_spm_']; 
     case 'qcseg',  qcv = [qaversions{qai} '_qcseg_']; Pp0 = cat_io_strrep(Pp0,[filesep 'mri' filesep 'p0'],filesep); 
    end
  cat_vol_qa('p0',Pp0,struct('prefix',qcv,'version',qaversions{ qai },'rerun',rerun));
end
verb = 'on';

%
clear NCR ICR IQR ECR SIQR rGMV
for qai = qais

  %% find xml files 
  switch segment{si}
     case 'CAT'
        Pxml{1} = spm_file(cat_io_strrep(Pp0a,['mri' filesep 'p0'],['report' filesep]), ...
          'prefix',sprintf('%s_',qaversions{ qai }),'ext','.xml'); 
        Pxml{2} = cat_io_strrep(Pxml{1},{'ATLASR1.1a';'ATLASmsk'},{'ATLASR1.1b';'ATLASorg'});
     case 'SPM'
        Pxml{1} = spm_file(cat_io_strrep(Pp0a,['c1'],['report' filesep]), ...
          'prefix',sprintf('%s_spm_',qaversions{ qai }),'ext','.xml'); 
        Pxml{2} = cat_io_strrep(Pxml{1},{'ATLASR1.1a';'ATLASmsk'},{'ATLASR1.1b';'ATLASorg'});
     case 'qcseg'
        Pxml{1} = spm_file(cat_io_strrep(Pp0a,['p0_'],['report' filesep]), ...
          'prefix',sprintf('%s_',qaversions{ qai }),'ext','.xml'); 
        Pxml{2} = cat_io_strrep(Pxml{1},{'ATLASR1.1a';'ATLASmsk'},{'ATLASR1.1b';'ATLASorg'});
    end
  for xsi = numel(Pxml{1}):-1:1
    if ~exist(Pxml{1}{xsi},'file') || ~exist(Pxml{2}{xsi},'file')
      Pxml{1}(xsi) = [];
      Pxml{2}(xsi) = [];
    end
  end

  % load xml data
  xml = cell(1,2); for i=1:numel(Pxml), xml{i}  = cat_io_xml(Pxml{i}); end

  % extract measure
  for i=1:numel(Pxml)
    for si = 1:numel(xml{i})
      try
        fname{si,i} = xml{i}(si).filedata.file;
        TIV(si,i)   = xml{i}(si).subjectratings.vol_TIV;
        NCR(si,i)   = xml{i}(si).qualityratings.NCR;
        ICR(si,i)   = xml{i}(si).qualityratings.ICR;
        if isfield(xml{i}(si).qualityratings,'FEC')
          FEC(si,i) = xml{i}(si).qualityratings.FEC;
        else
          FEC(si,i) = nan; 
        end
        ECR(si,i)   = xml{i}(si).qualityratings.res_ECR;
        IQR(si,i)   = xml{i}(si).qualityratings.IQR;
        SIQR(si,i)  = xml{i}(si).qualityratings.SIQR;
        rGMV(si,i)  = xml{i}(si).subjectmeasures.vol_rel_CGW(2);
       
        % get ID of the XML entry
        seps        = find( fname{si} == '_' );
        siten{si,i} = fname{si}(seps(1)+2:seps(1)+5);    
      end
    end
  end
  dIQR = diff(IQR,1,2); % msk - org
  rmse = @(x) mean(x.^2).^.5;  

  % store for later
  QAVdiff{qai} = dIQR; 
  if 0
   
    %numel(Pp0b( abs(QAVdiff{qai}) > .5 )) / numel(Pp0b); 
    %%
    ss = find(abs(QAVdiff{qai}) > 1);
    %% ss = ss(11); 
    switch segment{si}
      case 'CAT',    qcv = [qaversions{qai} '_']; 
      case 'qcseg',  qcv = [qaversions{qai} '_qcseg_']; 
    end
    cat_vol_qa('p0',Pp0a( ss ) ,struct('prefix',qcv,'version',qaversions{ qai },'rerun',1));   
    cat_vol_qa('p0',Pp0b( ss ) ,struct('prefix',qcv,'version',qaversions{ qai },'rerun',1));   
  end

  
  %% create figure
  if 1
    for sm = 0:1
      if ~exist('fh','var') || ~isvalid(fh), fh = figure(39); end 
      if ~fasttest, fh.Visible = verb; end
      if sm, fh.Position(3:4) = [200 160]; else, fh.Position(3:4) = [300 250]; end 
      ss  = 0.05 * (sm + 1); clf(fh);  
      try
        hst = histogram(dIQR,round(-max(abs(dIQR))/ss)*ss - ss/2:ss:round(max(abs(dIQR))/ss)*ss);
      catch
        hst = histogram(dIQR);
      end
      ylabel('number');  grid on; 
      if sm 
        title('masked vs. original'); smstr = '_sm';
        xlabel(sprintf('SIQR (RMSE=%0.3f)', rmse(dIQR) )); hold on; 
      else
        title('SIQR difference between masked and original images'); smstr = '';
        xlabel(sprintf('SIQR difference (RMSE=%0.4f) %s', rmse(dIQR), strrep(qaversions{qai},'_','\_') )); hold on; 
        %ph = plot([0 0],ylim); ph.Color = [0 0 0.5]; ph.LineStyle = '--'; ph.LineWidth = 1;    
      end
      
      xlim([-3  3]); %xlim([-1  1] * ceil(max(abs(dIQR)))); 
      ylim([ 0  1] * ceil( max(hst.Values)*1.1 / 10)*10);
    
      % print figure
      print(fh,fullfile(resdir  ,sprintf('ATLAS_%s%s',qaversions{qai},smstr)),'-r600','-djpeg');
    end
  end

  %% create figure
  for sm = 0:1
    if ~exist('fh','var'), fh = figure(39); end  
    if ~fasttest, fh.Visible = verb; end
    if sm, fh.Position(3:4) = [160 160]; else, fh.Position(3:4) = [150 250]; end 
    clf(fh);  
    if sm 
      title('masked vs. original'); smstr = '_sm';
      xlabel(sprintf('SIQR (RMSE=%0.3f)', rmse(dIQR) )); hold on; 
    else
      title('SIQR difference between masked and original images'); smstr = '';
      xlabel(sprintf('SIQR difference (RMSE=%0.4f) %s', rmse(dIQR), strrep(qaversions{qai},'_','\_') )); hold on; 
      %ph = plot([0 0],ylim); ph.Color = [0 0 0.5]; ph.LineStyle = '--'; ph.LineWidth = 1;    
    end
   cat_plot_boxplot( QAVdiff(qai) , struct('names',{{sprintf('SIQR (RMSE=%0.3f)', rmse(dIQR) )}},'usescatter',1,...
      'style',4,'datasymbol','o','ylim', [-3 3],'ygrid',0 )); %[-.1 .1] * ceil(max(abs( QAVdiff{qai} ))) ) ); 
   title('masked vs. original'); set(gca,'YGrid',1); 

    % print figure
    print(fh,fullfile(resdir  ,sprintf('ATLAS_boxplot_%s%s',qaversions{qai},smstr)),'-r600','-djpeg');
  end
end




%% final figure
fh = figure(38); fh.Position(3:4) = [150 300];
cat_plot_boxplot( QAVdiff , struct('names',{qaversions},'style',4,...
  'usescatter',1,'datasymbol','o','ylim',[-3 3],'ygrid',0) ); 
title('SIQR difference between masked and original images'); 
xlabel('QA versions'); ylabel('SIQR difference');
ax  = gca;  ax.FontSize = 10; ax.YGrid = 1; hold on; 
%ph = plot([0.4 numel(QAVdiff)+.6],[0 0]); ph.Color = [0 0 0.5]; ph.LineStyle = '--'; ph.LineWidth = 1;      
ylim([-3 3]); 

% print figure
print(fh,fullfile(resdir  ,sprintf('ATLAS_boxplot_all')),'-r600','-djpeg');



%% estimate size of leason 
Po = spm_file(cat_io_strrep(Pp0a,['mri' filesep 'p0'],'')); 
for pi = 1:numel(Pp0a)
  Plesmat{pi} = spm_file(Pp0a{pi},'prefix','lesionvol_','ext','.mat'); 
  if exist(Plesmat{pi},'file')
    load(Plesmat{pi},'Vles');
  else
    try
      V = spm_vol(Pp0a{pi}); 
    catch
      V = spm_vol(Po{pi}); 
    end
    Y = spm_read_vols(V);
    TIVvx = sum(Y(:) > 0);

    V = spm_vol(Po{pi}); 
    Y = spm_read_vols(V); 
    Y = smooth3( Y ) > .5; 
  
    Vles = sum(Y(:)==0) ./ TIVvx; 
    save(Plesmat{pi},'Vles');
  end
  Vlesion(pi) = Vles; 
end


%
fh     = figure(37); clf; hold on
fh.Position(3:4) = [300 200]; 
cmap   = cat_io_colormaps('nejm',numel(qais)); 
marker = {'o'  's'  'd'  'v'  '^'  '>'  '<'  'p'  'h'  }; 
for qai = qais
  hs = scatter( Vlesion , QAVdiff{qai} ); hold on; 
  hs.MarkerEdgeColor = cmap(qai,:); hs.MarkerFaceColor =  hs.MarkerEdgeColor; 
  hs.MarkerFaceAlpha = .2; hs.MarkerEdgeAlpha = .2; hs.SizeData = 20;  hs.Marker = marker{qai}; 
end
ylim([-3,3]);
for qai = qais
  lesfit{qai} = fit( Vlesion(~isinf(Vlesion))' , QAVdiff{qai}(~isinf(Vlesion))  ,'poly1'); hp = plot(lesfit{qai}); hp.Color = cmap(qai,:); %#ok<*SAGROW> 
end
if numel(qaversions)>1, legend(strrep(qaversions,'_','\_')); else, legend off; end
set(gca,'XTick',0:10); 
grid on; box on; 
ylabel('SIQR error (masked - raw)');   
xlabel('lesion volume in % of the TIV');
title('ATLAS SQIR difference by lesion volume');

% print figure
print(fh,fullfile(resdir  ,sprintf('ATLAS_scatterplot_lesionsize')),'-r1200','-djpeg');

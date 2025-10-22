function cat_tst_qa_ATLAS2( datadir0, qaversions, segment, fasttest, recalcQC ) 
%% Evaluation of CATQC in ATLAS 
%  ------------------------------------------------------------------------
%  QC evaluation of the ATLAS lesion dataset for original and lesion-masked 
%  images. In the ideal case the measures are somewhat similar and the QC 
%  measures are not or less affected by the tissue changes in lesions.
% 
%  Requirements: 
%   1. Matlab with curve fitting toolbox (fit)
%   2. Download and install SPM and CAT
%   3. Download ATLAS T1 data from: 
%
%   4. Specify in this script: 
%      1) the data directory "datadir" 
%      2) the QC version you would like to tests (the file has to exist in the cat directory) 
%      3) the segmentation you would like to use
%
%  See also cat_tst_qa_main.
%  ------------------------------------------------------------------------

cat_io_cprintf([0 0.5 0],'\n\n== Run cat_tst_qa_ATLAS ==\n') 
 
if license('test', 'Curve_Fitting_Toolbox')
    error('This function requires the "Curve Fitting Toolbox" of MATLAB.\n')
end

% ### datadir ###
if ~exist( 'datadir0' , 'var' )
  datadir  = '/Volumes/WDE18TB/MRData/Dahnke2025_QC/ATLAS_2/Training'; 
else
  datadir  = fullfile(datadir0,'ATLAS_2','Training'); 
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
  
resdir     = fullfile( fileparts( fileparts(datadir) ), '+results', ['ATLAS_' fast{fasttest+1} '_202508']); %' datestr(clock,'YYYYmm')]); 
if ~exist(resdir,'dir'), mkdir(resdir); end

% directories
rerun      = recalcQC; 

% prepare data
% - gunzip 
% - imcalc 

ATLASfiles = cat_vol_findfiles( datadir , 'sub*.nii.gz',struct('depth',5));
if ~isempty(ATLASfiles)
  gunzip(ATLASfiles); 
  for fi = 1:numel( ATLASfiles), delete(ATLASfiles{fi}); end
end


%% mask atlas files 
ATLASfiles    = cat_vol_findfiles( datadir , 'sub*T1w.nii',struct('depth',5));
mskATLASfiles = cat_vol_findfiles( datadir , 'masked_sub*T1w.nii',struct('depth',5));
if numel(ATLASfiles) > numel(mskATLASfiles) 
  t1w = cat_vol_findfiles( datadir , 'sub*T1w.nii' ,struct('depth',5)); 
  msk = cat_vol_findfiles( datadir , 'sub*mask.nii',struct('depth',5));
  te  = cellfun(@exist, spm_file(t1w,'prefix','masked_'));
  matlabbatch{1}.spm.tools.cat.tools.maskimg.data   = t1w(~te);
  matlabbatch{1}.spm.tools.cat.tools.maskimg.mask   = msk(~te); 
  matlabbatch{1}.spm.tools.cat.tools.maskimg.bmask  = {''};
  matlabbatch{1}.spm.tools.cat.tools.maskimg.recalc = 1;
  matlabbatch{1}.spm.tools.cat.tools.maskimg.prefix = 'masked_';
  spm_jobman('run',matlabbatch); clear matlabbatch
  mskATLASfiles = cat_vol_findfiles( datadir , 'masked_sub*.nii',struct('depth',5));
end



%% preprocessing
for si = 1:numel(segment)
  clear matlabbatch; 
  switch segment{si}
    case 'CAT'
      CATpreprocessing4qc;
      ATLASfilesCAT = [ATLASfiles; mskATLASfiles];
      ATLASfilesCAT( cellfun(@(x) exist(x,'file'),spm_file(ATLASfilesCAT,'prefix',['mri' filesep 'p0']))>0 ) = [];
      if ~isempty( ATLASfilesCAT )
        matlabbatch{1}.spm.tools.cat.estwrite.data = ATLASfilesCAT; 
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.lazy = 1; 
        spm_jobman('run',matlabbatch);
      end
    case 'SPM'
      SPMpreprocessing4qc;
      ATLASfilesSPM = [ATLASfiles; mskATLASfiles];
      ATLASfilesSPM( cellfun(@(x) exist(x,'file'),spm_file(ATLASfilesSPM,'prefix','c1'))>0 ) = [];
      if ~isempty( ATLASfilesSPM )
        matlabbatch{1}.spm.spatial.preproc.channel.vols = ATLASfilesSPM;
        spm_jobman('run',matlabbatch);
      end
    case 'synthseg'
      error('synthseg is not prepared in the public script ... use MATLAB help')
    case 'qcseg'
      fprintf('No preprocessing required.\n\n');
  end
end


%%
for si = 1:numel(segment)
  qais = 1:numel(qaversions);
  % (re)process QC values
  switch segment{si}
     case 'CAT'  
       Pp0b = cat_vol_findfiles( datadir , 'p0masked_sub*.nii',struct('depth',6)); %spm_fileparts(datadir) 
     case 'SPM'
       Pp0b = cat_vol_findfiles( datadir , 'c1masked_sub*.nii',struct('depth',5)); %spm_fileparts(datadir) 
     case 'synthseg'
     case 'qcseg'
       Pp0b = cat_vol_findfiles( datadir , 'p0_qcseg_masked_sub*.nii',struct('depth',5)); %spm_fileparts(datadir) 
  end
  Pp0a = cat_io_strrep(Pp0b,{'masked_'},{''}); %cat_vol_findfiles( datadira , 'p0*'); %spm_fileparts(datadir) 
  if fasttest % 5 ~> 50 scans
    Pp0a = Pp0a(1:12:end);
    Pp0b = Pp0b(1:12:end);
  end
  Pp0 = [Pp0a, Pp0b]; Pp0 = Pp0'; Pp0 = Pp0(:);  % [msk, org]
  if rerun
    for qai = qais
       switch segment{si}
         case 'CAT',    qcv = [qaversions{qai} '_']; 
         case 'SPM',    qcv = [qaversions{qai} '_spm_']; 
         case 'qcseg',  qcv = [qaversions{qai} '_qcseg_']; Pp0 = cat_io_strrep(Pp0,[filesep 'mri' filesep 'p0'],filesep); 
        end
      cat_vol_qa('p0',Pp0,struct('prefix',qcv,'version',qaversions{ qai },'rerun',rerun));
    end
  end
  verb = 'on';



  %% estimate size of leason 
  fprintf('Get lesion size:            '); 
  switch segment{si}
    case 'CAT'
      Po = spm_file(cat_io_strrep(Pp0b,{['mri' filesep 'p0masked_'],'T1w'},{'','label-L_desc-T1lesion_mask'})); 
    case 'SPM'
      Po = spm_file(cat_io_strrep(Pp0b,{'c1masked_','T1w'},{'','label-L_desc-T1lesion_mask'})); 
  end
  for pi = 1:numel(Pp0b)
    Plesmat{pi} = spm_file(Pp0b{pi},'prefix','lesionvol_','ext','.mat'); 
    if exist(Plesmat{pi},'file')
      load(Plesmat{pi},'Vles');
    else
      fprintf('\b\b\b\b\b\b\b\b\b%4d/%4d',pi,numel(Pp0b)); 
      try
        V = spm_vol(Pp0a{pi}); 
      catch
        V = spm_vol(Po{pi}); 
      end
      Y = spm_read_vols(V);
      TIVvx = sum(Y(:) > 0.5);
  
      V = spm_vol(Po{pi}); 
      Y = spm_read_vols(V); 
    
      Vles = sum(Y(:) > 0.5) ./ TIVvx; 
      save(Plesmat{pi},'Vles');
    end
    Vlesion(pi) = Vles; 
  end


  %%
  clear NCR ICR IQR ECR SIQR rGMV
  for qai = qais
  
    %% find xml files 
    switch segment{si}
       case 'CAT'
          Pxml{1} = spm_file(cat_io_strrep(Pp0b,['mri' filesep 'p0'],['report' filesep]), ...
            'prefix',sprintf('%s_',qaversions{ qai }),'ext','.xml'); 
          Pxml{2} = cat_io_strrep(Pxml{1},{'masked_'},{''});
       case 'SPM'
          Pxml{1} = spm_file(cat_io_strrep(Pp0b,'c1',['report' filesep]), ...
            'prefix',sprintf('%s_spm_',qaversions{ qai }),'ext','.xml'); 
          Pxml{2} = cat_io_strrep(Pxml{1},{'masked_'},{''});
       case 'qcseg'
          Pxml{1} = spm_file(cat_io_strrep(Pp0b,'p0_',['report' filesep]), ...
            'prefix',sprintf('%s_',qaversions{ qai }),'ext','.xml'); 
          Pxml{2} = cat_io_strrep(Pxml{1},{'masked_'},{''});
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
      for ssi = 1:numel(xml{i})
        try
          fname{ssi,i} = xml{i}(ssi).filedata.file;
          TIV(ssi,i)   = xml{i}(ssi).subjectratings.vol_TIV;
          NCR(ssi,i)   = xml{i}(ssi).qualityratings.NCR;
          ICR(ssi,i)   = xml{i}(ssi).qualityratings.ICR;
          if isfield(xml{i}(ssi).qualityratings,'FEC')
            FEC(ssi,i) = xml{i}(ssi).qualityratings.FEC;
          else
            FEC(ssi,i) = nan; 
          end
          ECR(ssi,i)   = xml{i}(ssi).qualityratings.res_ECR;
          IQR(ssi,i)   = xml{i}(ssi).qualityratings.IQR;
          SIQR(ssi,i)  = xml{i}(ssi).qualityratings.SIQR;
          rGMV(ssi,i)  = xml{i}(ssi).subjectmeasures.vol_rel_CGW(2);
         
          % get ID of the XML entry
          seps        = find( fname{si} == '_' );
          siten{ssi,i} = fname{ssi}(seps(1)+2:seps(1)+5);    
        end
      end
    end
    dIQR = diff(SIQR,1,2); % msk - org
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
      for sm = 0
        if ~exist('fh','var') || ~isvalid(fh), fh = figure(39); end 
        fh.Interruptible = 'off'; fh.Visible = 'off'; 
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
        print(fh,fullfile(resdir  ,sprintf('ATLAS_%s%s_%s',qaversions{qai},smstr,segment{si})),'-r600','-djpeg');
      end
    end
  
    %% create figure
    for sm = 0%:1
      if ~exist('fh','var'), fh = figure(39); end  
      fh.Interruptible = 'off'; fh.Visible = 'off'; 
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
      print(fh,fullfile(resdir  ,sprintf('ATLAS_boxplot_%s%s_%s',qaversions{qai},smstr,segment{si})),'-r600','-djpeg');
    end
  
  
    %% lesion size figure
    fh     = figure(37); clf; hold on
    fh.Position(3:4) = [300 200];  fh.Interruptible = 'off'; fh.Visible = 'off'; 
    cmap   = cat_io_colormaps('nejm',numel(qais)); 
    marker = {'o'  's'  'd'  'v'  '^'  '>'  '<'  'p'  'h'  }; 
  
    hs = scatter( Vlesion , QAVdiff{qai} ); hold on; 
    hs.MarkerEdgeColor = cmap(1,:); hs.MarkerFaceColor =  hs.MarkerEdgeColor; 
    hs.MarkerFaceAlpha = .2; hs.MarkerEdgeAlpha = .2; hs.SizeData = 20;  hs.Marker = marker{1}; 
    ylim([-3,3]);
    lesfit{qai} = fit( Vlesion(~isinf(Vlesion))' , QAVdiff{qai}(~isinf(Vlesion))  ,'poly1'); hp = plot(lesfit{qai}); hp.Color = cmap(1,:); %#ok<*SAGROW> 
    legend(strrep(qaversions{qai},'_','\_')); 
    
    set(gca,'XTick',0:0.05:0.3,'XTickLabel',0:5:30); xlim([0 .3]); 
    grid on; box on; 
    ylabel('SIQR error (masked - raw)');   
    xlabel('lesion volume in % of the TIV');
    title('ATLAS SQIR difference by lesion volume');
    
    % print figure
    print(fh,fullfile(resdir  ,sprintf('ATLAS_lesionsize_%s_%s',qaversions{qai},segment{si})),'-r1200','-djpeg');
  
  
  end
  
end
fprintf('ATLAS done.\n\n');   


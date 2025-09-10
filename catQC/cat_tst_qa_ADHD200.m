function cat_tst_qa_ADHD200( datadir, qaversions, segment, fasttest)    
%cat_tst_qa_ADHD200.
%
%   1) Downlaod raw data: 
%      https://fcon_1000.projects.nitrc.org/indi/adhd200/index.html
%   2) Unzip T1 files into a train_raw and test_raw directory
%   3) There should be also a phenotypic directory with the
%      ADHD200.csv with age, sex and site information 
%


 %% setup
  dataset = 'ADHD200';
  datadir = fullfile(datadir,dataset); 
  CATver  = 'CAT12.9_2874'; 
  resdir  = fullfile(datadir,'+results','ADHD200');
  if fasttest
    files = cat_vol_findfiles( fullfile( datadir , 'test_raw') , '*.nii.gz',struct('depth',1));
  else
    files = [
      cat_vol_findfiles( fullfile( datadir , 'train_raw'), '*.nii.gz',struct('depth',1));
      cat_vol_findfiles( fullfile( datadir , 'test_raw'), '*.nii.gz',struct('depth',1));
      ];
  end
  if ~exist(resdir,'dir'), mkdir(resdir); end


%% run preprocessing
  for si = 1:numel(segment)
    clear matlabbatch; 
    switch segment{si}
      case 'CAT'
        CATpreprocessing4qc;
        filesCAT = files;
        matlabbatch{1}.spm.tools.cat.estwrite.output.BIDS.BIDSyes.BIDSdir = ...
          fullfile('..','derivatives',CATver);
        filesCAT2 = strrep( filesCAT , datadir , fullfile( datadir, 'derivatives',CATver) ); 
        filesCAT( cellfun(@(x) exist(x,'file'), ...
          spm_file(filesCAT2, 'prefix','p0') )>0 ) = [];          
        if ~isempty( filesCAT )
          matlabbatch{1}.spm.tools.cat.estwrite.data = filesCAT; 
          matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.lazy = 1; 
          spm_jobman('run',matlabbatch);
        end
      case 'SPM'
        SPMpreprocessing4qc;
        IXIfilesSPM = files; 
        IXIfilesSPM( cellfun(@(x) exist(x,'file'),spm_file(IXIfilesSPM,'prefix','c1'))>0 ) = [];
        if ~isempty( IXIfilesSPM )
          matlabbatch{1}.spm.spatial.preproc.channel.vols = IXIfilesSPM;
          spm_jobman('run',matlabbatch);
        end
      case 'synthseg'
        error('synthseg is not prepared in the public script ... use MATLAB help')
      case 'qcseg'
        fprintf('No preprocessing required.\n\n');
    end
  end


% run QC

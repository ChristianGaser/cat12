function cat_run_newcatch(job,estwrite,tpm,subj)
% ______________________________________________________________________
% This function contains the new matlab try-catch block.
% The new try-catch block has to be in a separate file to avoid an error.
%
% See also cat_run_newcatch.
% ______________________________________________________________________
% $Revision$  $Date$
  try
    cat_run_job(job,estwrite,tpm,subj); 
  catch caterr 
    % add further information for special errors
    if isempty(caterr.identifier)
      switch caterr.message
        case 'insufficient image overlap'
          adderr = MException('SPM:AlignmentError','There is not enough overlap in the images to obtain a solution.');
        otherwise
          adderr = MException('SPM:CAT:cat_main',strrep(caterr.message,'\','\\'));
      end
      caterr = addCause(caterr,adderr);
    end
    
    cat_io_cprintf('err',sprintf('\n%s\nCAT Preprocessing error: %s: %s \n%s\n%s\n%s\n', ...
      repmat('-',1,72),caterr.identifier,...
      spm_str_manip(job.channel(1).vols{subj},'a60'),...
      repmat('-',1,72),caterr.message,repmat('-',1,72)));  

    % write error report
    caterrtxt = cell(numel(caterr.stack),1);
    for si=1:numel(caterr.stack)
      cat_io_cprintf('err',sprintf('%5d - %s\n',caterr.stack(si).line,caterr.stack(si).name));  
      caterrtxt{si} = sprintf('%5d - %s\n',caterr.stack(si).line,caterr.stack(si).name); 
    end
    cat_io_cprintf('err',sprintf('%s\n',repmat('-',1,72)));  

    % delete template files 
    [pth,nam,ext] = spm_fileparts(job.channel(1).vols{subj}); 
    % delete noise corrected image
    if exist(fullfile(pth,['n' nam ext]),'file')
      try %#ok<TRYNC>
        delete(fullfile(pth,['n' nam ext]));
      end
    end

    % save cat xml file
    caterrstruct = struct();
    for si=1:numel(caterr.stack)
      caterrstruct(si).line = caterr.stack(si).line;
      caterrstruct(si).name = caterr.stack(si).name;  
      caterrstruct(si).file = caterr.stack(si).file;  
    end
    cat_tst_qa('cat12err',struct('write_csv',0,'write_xml',1,'caterrtxt',caterrtxt,'caterr',caterrstruct,'job',job));
    
    % rethrow error 
    if ~cat_get_defaults('extopts.ignoreErrors')
      rethrow(caterr); 
    end 
  end
end

%{

Special errors that may should get further descriptions
========================================================================        
- MATLAB:badfid_mx
- MATLAB:error:nonScalarInput
- CAT-ERROR: MATLAB:nonExistentField in ./Apes_BET/ncapuchin_diva: QMm fehlt in 226 - cat_stat_marks
- squirrel_2 > Bus Error!


- Warning: ERROR: cat_vol_morph - lab - no object! 
  * Ventricle detection                                              
  * LAS
  * gcut


% QA warning handling
========================================================================        
------------------------------------------------------------------------
CAT-ERROR: MATLAB:nonExistentField in ./Apes_BET/ncapuchin_diva:
Reference to non-existent field 'QMm'.
------------------------------------------------------------------------
  226 - cat_stat_marks
 1368 - cat_main
  372 - run_job
  100 - cat_run
   29 - cfg_run_cm
 1631 - local_runcj
  915 - cfg_util
  456 - fill_run_job
  227 - spm_jobman
  113 - cat_batch_vbm
------------------------------------------------------------------------
========================================================================        
        
MATLAB:cat_io_xml:write: Can't write XML-file '/Volumes/vbmDB/MRData/cat12tst/results/deffiles/cat_defaults/BWPC_full/cat_BWP_HC_Collins_t1_pn7_rf040pA.xml'!
1010:Too many outputs requested.  Most likely cause is missing [] around
left hand side that has a comma separated list expansion.
        
        
========================================================================        
        ROI estimation:                                                    
  ROI mapping data to group space                                 
------------------------------------------------------------------------
CAT Preprocessing error: MATLAB:badfid_mx: 
./BWPC_full/BWP_HC_Collins_t1_pn7_rf060pC.nii 
./BWPC_noise/BWP_HC_Collins_t1_pn1_rf040pA.nii 
------------------------------------------------------------------------
Invalid file identifier.  Use fopen to generate a valid file identifier.
------------------------------------------------------------------------
  109 - read_hdr_raw
   30 - read_hdr
   26 - nifti
   18 - spm_vol_nifti
  102 - spm_vol_hdr
   60 - spm_vol
 2697 - cat_vol_ROInorm
 1328 - cat_main
  247 - cat_run_job
   10 - cat_run_newcatch
  141 - run_job
  100 - cat_run
   29 - cfg_run_cm
 1631 - local_runcj
  915 - cfg_util
  456 - fill_run_job
  227 - spm_jobman
  113 - cat_batch_vbm
------------------------------------------------------------------------
========================================================================     
------------------------------------------------------------------------
CAT12 r0:    ./deffiles/cat_defaults/BO_inverse/IXI014-HH-1236-PD.nii
------------------------------------------------------------------------
CAT Preprocessing error: MATLAB:error:nonScalarInput: ./deffiles/cat_defaults/BO_inverse/IXI014-HH-1236-PD.nii 
------------------------------------------------------------------------
Formatted arguments cannot be non-scalar numeric matrices.
------------------------------------------------------------------------
   34 - cat_run_job
   10 - cat_run_newcatch
  141 - run_job
  100 - cat_run
   29 - cfg_run_cm
 1631 - local_runcj
  915 - cfg_util
  456 - fill_run_job
  227 - spm_jobman
  113 - cat_batch_vbm
------------------------------------------------------------------------
========================================================================     
------------------------------------------------------------------------
CAT Preprocessing error: MATLAB:extraneousStrucRefBlocks: ./SRS_SRS_Jena_DaRo81_T1_20140718-rs07_SI30_TrioTim_MPRAGE08mms009a1001.nii 
------------------------------------------------------------------------
Field reference for multiple structure elements that is followed by more reference blocks is an error.
------------------------------------------------------------------------
   31 - cat_run_job
   10 - cat_run_newcatch
  141 - run_job
  100 - cat_run
   29 - cfg_run_cm
 1631 - local_runcj
  915 - cfg_util
  456 - fill_run_job
  227 - spm_jobman
  113 - cat_batch_vbm
------------------------------------------------------------------------
>> multidimension error
========================================================================     
------------------------------------------------------------------------
CAT12 r0:    ./deffiles/cat_defaults/BO_inverse/IXI014-HH-1236-PD.nii
------------------------------------------------------------------------

------------------------------------------------------------------------
CAT Preprocessing error: MATLAB:error:nonScalarInput: ./deffiles/cat_defaults/BO_inverse/IXI014-HH-1236-PD.nii 
------------------------------------------------------------------------
Formatted arguments cannot be non-scalar numeric matrices.
------------------------------------------------------------------------
   34 - cat_run_job
   10 - cat_run_newcatch
  141 - run_job
  100 - cat_run
   29 - cfg_run_cm
 1631 - local_runcj
  915 - cfg_util
  456 - fill_run_job
  227 - spm_jobman
  113 - cat_batch_vbm
------------------------------------------------------------------------
========================================================================     


    %}

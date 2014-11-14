function cg_vbm_run_newcatch(job,estwrite,tpm,subj)
% ______________________________________________________________________
% This function contains the new matlab try-catch block.
% The new try-catch block has to be in a separate file to avoid an error.
%
% See also cg_vbm_run_newcatch.
% ______________________________________________________________________
% $Revision$  $Date$
  try
    cg_vbm_run_job(job,estwrite,tpm,subj);
  catch vbmerr 

    % add further information for special errors
    if isempty(vbmerr.identifier)
      switch vbmerr.message
        case 'insufficient image overlap'
          adderr = MException('SPM:AlignmentError','There is not enough overlap in the images to obtain a solution.');
        otherwise
          adderr = MException('SPM:VBM:cg_vbm_write',strrep(vbmerr.message,'\','\\'));
      end
      vbmerr = addCause(vbmerr,adderr);
    end
    
    vbm_io_cprintf('err',sprintf('\n%s\nVBM Preprocessing error: %s: %s \n%s\n%s\n%s\n', ...
      repmat('-',1,72),vbmerr.identifier,...
      spm_str_manip(job.channel(1).vols{subj},'a60'),...
      repmat('-',1,72),vbmerr.message,repmat('-',1,72)));  

    % write error report
    for si=1:numel(vbmerr.stack)
      vbm_io_cprintf('err',sprintf('%5d - %s\n',vbmerr.stack(si).line,vbmerr.stack(si).name));  
    end
    vbm_io_cprintf('err',sprintf('%s\n',repmat('-',1,72)));  

    % delete template files 
    [pth,nam,ext] = spm_fileparts(job.channel(1).vols{subj}); 
    % delete bias map
    if exist(fullfile(pth,['bf' nam(2:end) ext]),'file')
      delete(fullfile(pth,['bf' nam(2:end) ext]));
    end
    % delete noise corrected image
    if exist(fullfile(pth,['n' nam(2:end) ext]),'file')
      delete(fullfile(pth,['n' nam(2:end) ext]));
    end

    % rethrow error 
    if ~cg_vbm_get_defaults('extopts.ignoreErrors')
      rethrow(vbmerr); 
    end 
  end
end

%{

Special errors that may should get further descriptions
========================================================================        
- MATLAB:badfid_mx
- MATLAB:error:nonScalarInput
- VBM-ERROR: MATLAB:nonExistentField in ./Apes_BET/ncapuchin_diva: QMm fehlt in 226 - vbm_stat_marks
- squirrel_2 > Bus Error!


- Warning: ERROR: vbm_vol_morph - lab - no object! 
  * Ventricle detection                                              
  * LAS
  * gcut


% QA warning handling
========================================================================        
------------------------------------------------------------------------
VBM-ERROR: MATLAB:nonExistentField in ./Apes_BET/ncapuchin_diva:
Reference to non-existent field 'QMm'.
------------------------------------------------------------------------
  226 - vbm_stat_marks
 1368 - cg_vbm_write
  372 - run_job
  100 - cg_vbm_run
   29 - cfg_run_cm
 1631 - local_runcj
  915 - cfg_util
  456 - fill_run_job
  227 - spm_jobman
  113 - cg_vbm_batch
------------------------------------------------------------------------
========================================================================        
        
MATLAB:vbm_io_xml:write: Can't write XML-file '/Volumes/vbmDB/MRData/vbm12tst/results/deffiles/cg_vbm_defaults/BWPC_full/vbm_BWP_HC_Collins_t1_pn7_rf040pA.xml'!
1010:Too many outputs requested.  Most likely cause is missing [] around
left hand side that has a comma separated list expansion.
        
        
========================================================================        
        ROI estimation:                                                    
  ROI mapping data to group space                                 
------------------------------------------------------------------------
VBM Preprocessing error: MATLAB:badfid_mx: 
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
 2697 - vbm_vol_ROInorm
 1328 - cg_vbm_write
  247 - cg_vbm_run_job
   10 - cg_vbm_run_newcatch
  141 - run_job
  100 - cg_vbm_run
   29 - cfg_run_cm
 1631 - local_runcj
  915 - cfg_util
  456 - fill_run_job
  227 - spm_jobman
  113 - cg_vbm_batch
------------------------------------------------------------------------
========================================================================     
------------------------------------------------------------------------
VBM12 r0:    ./deffiles/cg_vbm_defaults/BO_inverse/IXI014-HH-1236-PD.nii
------------------------------------------------------------------------
VBM Preprocessing error: MATLAB:error:nonScalarInput: ./deffiles/cg_vbm_defaults/BO_inverse/IXI014-HH-1236-PD.nii 
------------------------------------------------------------------------
Formatted arguments cannot be non-scalar numeric matrices.
------------------------------------------------------------------------
   34 - cg_vbm_run_job
   10 - cg_vbm_run_newcatch
  141 - run_job
  100 - cg_vbm_run
   29 - cfg_run_cm
 1631 - local_runcj
  915 - cfg_util
  456 - fill_run_job
  227 - spm_jobman
  113 - cg_vbm_batch
------------------------------------------------------------------------
========================================================================     
------------------------------------------------------------------------
VBM Preprocessing error: MATLAB:extraneousStrucRefBlocks: ./SRS_SRS_Jena_DaRo81_T1_20140718-rs07_SI30_TrioTim_MPRAGE08mms009a1001.nii 
------------------------------------------------------------------------
Field reference for multiple structure elements that is followed by more reference blocks is an error.
------------------------------------------------------------------------
   31 - cg_vbm_run_job
   10 - cg_vbm_run_newcatch
  141 - run_job
  100 - cg_vbm_run
   29 - cfg_run_cm
 1631 - local_runcj
  915 - cfg_util
  456 - fill_run_job
  227 - spm_jobman
  113 - cg_vbm_batch
------------------------------------------------------------------------
>> multidimension error
========================================================================     
------------------------------------------------------------------------
VBM12 r0:    ./deffiles/cg_vbm_defaults/BO_inverse/IXI014-HH-1236-PD.nii
------------------------------------------------------------------------

------------------------------------------------------------------------
VBM Preprocessing error: MATLAB:error:nonScalarInput: ./deffiles/cg_vbm_defaults/BO_inverse/IXI014-HH-1236-PD.nii 
------------------------------------------------------------------------
Formatted arguments cannot be non-scalar numeric matrices.
------------------------------------------------------------------------
   34 - cg_vbm_run_job
   10 - cg_vbm_run_newcatch
  141 - run_job
  100 - cg_vbm_run
   29 - cfg_run_cm
 1631 - local_runcj
  915 - cfg_util
  456 - fill_run_job
  227 - spm_jobman
  113 - cg_vbm_batch
------------------------------------------------------------------------
========================================================================     


    %}

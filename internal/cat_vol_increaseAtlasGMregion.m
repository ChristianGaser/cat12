function cat_vol_increaseAtlasGMregion
%_______________________________________________________________________
% Simple internal function to improve atlas maps that were thresholded  
% by tissue probality such as the neuromorphometrics atlas. 
%_______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
%
% $Id$
% ______________________________________________________________________

  Pt = spm_select([1 1],'image','select TPM');
  Pa = cellstr(spm_select([1 Inf],'image','select atlas maps'));
  pth = 0.0;
  method = 'vbdist'; % vbdist, downcut 
  
  %% ROIs
  % region of changes
  DelRoi   = {'Cerebral White Matter'};
  % region that should not change
  noDilRoi = {
    'Ventricle'
    'Accumbens'
    'Amygdala'
    'Hippocampus'
    'Pallidum'
    'Putamen'
    'Thalamus'
    'Caudate'
    'Thalamus'
  };

  for ai=1:numel(Pa)
    %% prepare data
    [pp,ff,ee] = spm_fileparts(Pa{ai});
    
    % find csv-file
    csv = cat_io_csv(fullfile(pp,[ff '.csv']));
    id  = cell2mat(csv(2:end,1));
    roi = csv(2:end,3);
    
    DelRoiIds = []; 
    for dri=1:numel(DelRoi)
      DelRoiIds  = [DelRoiIds;id(~cellfun('isempty',strfind(roi,DelRoi{dri})))]; %#ok<AGROW>
    end
    DelRoiIds = unique(DelRoiIds);
    
    noDilRoiIds = [];
    for ndri=1:numel(noDilRoi)
      noDilRoiIds  = [noDilRoiIds;id(~cellfun('isempty',strfind(roi,noDilRoi{ndri})))]; %#ok<AGROW>
    end
    noDilRoiIds = unique(noDilRoiIds);
    
    %% load images
    Va = spm_vol(Pa{ai});
    Vt = spm_vol(Pt);
    
    Ya = spm_read_vols(Va);
    Yt = spm_read_vols(Vt);
    
    Ya2 = Ya;
    for dri=1:numel(DelRoiIds)
      Ya2(Ya2==DelRoiIds(dri))=0;
    end
    
    switch method 
      case 'vbdist'
        % simple distance
        [YD,YI] = cat_vbdist(single(Ya2>0),Yt>pth); Ya3=Ya(YI);

        Ya3(Ya3==0) = Ya(Ya3==0);
        for dri=1:numel(noDilRoiIds)
          Ya3(Ya3==noDilRoiIds(dri) & Ya~=noDilRoiIds(dri)) = ...
            Ya( Ya3==noDilRoiIds(dri) & Ya~=noDilRoiIds(dri));
        end
        Ya2 = Ya3;
      case 'downcut'
        % region-growing ... simpler is better! this is just for tests
        Ya2 = single(Ya2); Ya2(Yt<pth)=nan;
        Ya4 = cat_vol_downcut(Ya2,single(Yt),0.0); 

        Ya4(Ya4==0) = Ya(Ya4==0);
        for dri=1:numel(noDilRoiIds)
          Ya4(Ya4==noDilRoiIds(dri) & Ya~=noDilRoiIds(dri)) = ...
            Ya( Ya4==noDilRoiIds(dri) & Ya~=noDilRoiIds(dri));
        end
        Ya2 = Ya4;
    end
    % median-filter?
    % Ya2 = cat_vol_median3c(Ya2,Ya2>0);
    
    %% display things
    ds('l2','',1.5,Yt,single(Ya3)>0,single(Ya)/80,single(Ya2)/80,60)

    %% write result
    copyfile(fullfile(pp,[ff ee]),fullfile(pp,[ff '_orginal' ee]));
    spm_write_vol(Va,Ya2);
    
  end
end




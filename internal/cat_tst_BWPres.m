% This function use the BWP (Collins) to create resolution reduces 
% reinterpolated phantoms

BWPpaths = {
  '/Volumes/MyBook/MRData/vbm_tst/+RAW/BWP_Collins/T1/BWP_HC_Collins_t1_pn0_rf00000.nii'; 
};

resolutions  = single(1.25:0.25:3);
resolutions2 = single([1.25:0.25:3,5,7,9]);
res.iso      = [resolutions'  * ones(1,3);
                resolutions2.^(1/3)' * ones(1,3)];
res.aniso    = [(resolutions2' * zeros(1,2) + 1) resolutions2'];

clc; fprintf('BWP Resolutionupdate:\n');
for ph = 1:numel(BWPpaths)
  % load image
  hdr = spm_vol(BWPpaths{ph});
  img = single(spm_read_vols(hdr));
  
  resnames = fieldnames(res);
  for rni = 1:numel(resnames)
    for ri = 1:size(res.(resnames{rni}),1)
      % name
      if res.(resnames{rni})(ri,1) == res.(resnames{rni})(ri,2) && ...
         res.(resnames{rni})(ri,1) == res.(resnames{rni})(ri,3)
        iso = 'i';
      else
        iso = 'a';
      end
      
      hdri = hdr; hdrr = hdr; 
      [pp,ff,rr] = spm_fileparts(hdr.fname);
      hdri.fname = fullfile(pp,sprintf('%s_I%s%4.0fx%4.0fmy%s',ff,iso,... res.(resnames{rni})(ri,:).^3,...
        res.(resnames{rni})(ri,1:2:3)*1000,rr));
      hdrr.fname = fullfile(pp,sprintf('%s_R%s%4.0fx%4.0fmy%s',ff,iso,... res.(resnames{rni})(ri,:).^3,...
        res.(resnames{rni})(ri,1:2:3)*1000,rr));
      
      fprintf('%s\n',hdri.fname);
     
      
      
      % reduce resolution
      if 1
        [Rx,Ry,Rz] = meshgrid(1:res.(resnames{rni})(ri,2):size(img,2), ...
                              1:res.(resnames{rni})(ri,1):size(img,1), ...
                              1:res.(resnames{rni})(ri,3):size(img,3));                            
        imgr = vbm_vol_interp3f(img,Rx,Ry,Rz,'cubic');
        [Rx,Ry,Rz] = meshgrid(1:1/res.(resnames{rni})(ri,2):size(imgr,2), ...
                              1:1/res.(resnames{rni})(ri,1):size(imgr,1), ...
                              1:1/res.(resnames{rni})(ri,3):size(imgr,3));     
        imgit = vbm_vol_interp3f(imgr,Rx,Ry,Rz,'cubic');     
      else
        % the offset is necessary for even cases of the reduction 
        % (to have similar smoothing for all images)
        os = pi/10;
      
        [Rx,Ry,Rz] = meshgrid(0.5 + os:res.(resnames{rni})(ri,2):size(img,2) + 0.5 + os, ...
                              0.5 + os:res.(resnames{rni})(ri,1):size(img,1) + 0.5 + os, ...
                              0.5 + os:res.(resnames{rni})(ri,3):size(img,3) + 0.5 + os);                            
        imgr = vbm_vol_interp3f(img,Rx,Ry,Rz,'cubic');
        [Rx,Ry,Rz] = meshgrid(1-os/res.(resnames{rni})(ri,2):1/res.(resnames{rni})(ri,2):size(imgr,2), ...
                              1-os/res.(resnames{rni})(ri,1):1/res.(resnames{rni})(ri,1):size(imgr,1), ...
                              1-os/res.(resnames{rni})(ri,3):1/res.(resnames{rni})(ri,3):size(imgr,3));     
        imgit = vbm_vol_interp3f(imgr,Rx,Ry,Rz,'cubic');     
      end
      
      imgi  = zeros(size(img),'single');
      simgi = min(size(img),size(imgit));
      imgi(1:simgi(1),1:simgi(2),1:simgi(3)) = imgit(1:simgi(1),1:simgi(2),1:simgi(3));
      
      
      % save image 
      hdrr.dim = size(imgr); 
      hdrr.mat([1,6,11,13:15]) = (hdr.mat([1,6,11,13:15]) - 1*[zeros(1,3) ones(1,3).*(0.5+os)]) .* ...
        [res.(resnames{rni})(ri,:) 1./res.(resnames{rni})(ri,:)] + ...
        1*[zeros(1,3) ones(1,3).*res.(resnames{rni})(ri,:)/2];
      
      %spm_write_vol(hdrr,imgr);
     % vbm_io_writenii(hdrr,imgr)
      vbm_io_writenii(hdri,imgi)
     % spm_write_vol(hdri,imgi);
      clear imgr imgi imgit;
    end
    fprintf('\n');
  end
  
  
  
end



function [vbmSTC,Ystc] = vbm_qa_calcSTC(Yp0,VT0,trans,te,res)
% Subject Template  (TE)
% ----------------------------------------------------------------------
% This measure shoold describe the difference between our expectation
% from the mean group probability map and the subject. Strong variation
% can represent 
%   (1) strong anatomical variations of this subject, and 
%   (2) normalisation error (that are often caused by special anatomies
%       or be previous preprocessing errors)
% Stronger changes are expected in with growing distance from the core
% of the WM. 
% ----------------------------------------------------------------------

  opt.tpm = 1;
  if opt.tpm 
    %% VBM-Dartel template
    [pth1,nam1,ext1] = spm_fileparts(char(cg_vbm_get_defaults('extopts.darteltpm')));
    VclsA = spm_vol(fullfile(pth1,[strrep(nam1,'Template_1','Template_6'),ext1]));
    YclsA = cell(1,3);
    for i=1:2
      YclsA{i} = single(spm_sample_vol(VclsA(i), ...
                        double(trans.atlas.Yy(:,:,:,1)), ...
                        double(trans.atlas.Yy(:,:,:,2)), ...
                        double(trans.atlas.Yy(:,:,:,3)), 1));
      YclsA{i} = reshape(YclsA{i},size(Yp0));
    end

    % now we need to create a CSF probability map (for the next correction)
    Yclsb = vbm_vol_smooth3X(vbm_vol_morph((YclsA{1} + YclsA{2})>0.3,'lc',2),2)>0.5;
    for i=1:2, YclsA{i} = YclsA{i} .* smooth3(Yclsb); end
    YclsA{3}   = (Yclsb & smooth3(YclsA{1} + YclsA{2})<0.6) .* ...
                 ((Yclsb - single(YclsA{1} + YclsA{2}) ./ ...
                 median(YclsA{1}(Yclsb) + YclsA{2}(Yclsb)))); 
    % final correction for maximum probability of 1
    YclsAsum   = (YclsA{1} + YclsA{2} + YclsA{3}) .* Yclsb;
    for i=1:3, YclsA{i} = (YclsA{i}./max(eps,YclsAsum)) .* Yclsb; end
    Yp0A = YclsA{1}*2 + YclsA{2}*3 + YclsA{3} .* Yclsb;
    clear YclsA YclsAsum VclsA Yclsb;

  else
    %% SPM-Tissue template
    VclsB = spm_vol(res.tpm(1).fname);
    YclsB = cell(1,7);
    for i=1:3
      YclsB{i} = single(spm_sample_vol(VclsB(i), ...
                       double(trans.atlas.Yy(:,:,:,1)), ...
                       double(trans.atlas.Yy(:,:,:,2)), ...
                       double(trans.atlas.Yy(:,:,:,3)), 1));
      YclsB{i} = reshape(YclsB{i},d);
    end

    % now we need to create a CSF probability map (for the next correction)
    Yclsb = vbm_vol_smooth3X(vbm_vol_morph((YclsB{1} + YclsB{2})>0.3,'lc',2),2)>0.5;
    for i=1:3, YclsB{i} = YclsB{i} .* smooth3(Yclsb); end
    % final correction for maximum probability of 1
    YclsBsum   = (YclsB{1} + YclsB{2} + YclsB{3}) .* Yclsb;
    for i=1:3, YclsB{i} = (YclsB{i}./max(eps,YclsBsum)) .* Yclsb; end
    Yp0A = YclsB{1}*2 + YclsB{2}*3 + YclsB{3} .* Yclsb;
    clear YclsB YclsBsum Yclsb VclsB;
  end

  % Now we can estimate the difference maps for each intensity/Yp0b map.
  % But finally only our segment/Yp0b map is important, because other
  % non-intensity scaled images will have higher errors due to the
  % intensity scaling.
  Ystc = abs(max(1,Yp0A)-max(1,Yp0)); % we are not interessed in skull-stripping differences... maybe later ;-)
  spm_smooth(Ystc,Ystc,8);  % we are only interessed on larger changes

  vbm_io_writenii(VT0,Ystc,'te', ...
    'group expectation map (matching of template after normalization)', ...
    'uint8',[0,1/255],min([1 0 0 0],cell2mat(struct2cell(te)')),0,trans);
  vbm_io_writenii(VT0,Ystc,'te', ...
    'group expectation map (matching of template after normalization)', ...
    'uint8',[0,1/255],min([0 1 2 2],cell2mat(struct2cell(te)')),0,trans);
  vbmSTC = sum(Ystc(:))./sum(Yp0(:)>0);

end
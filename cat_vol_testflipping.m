function [flipped,flippedval,stime] = cat_vol_testflipping(varargin)
%cat_vol_testflipping. Test input image for LR flipping. 
% Side flipping of images can accidently occurring by converting images types
% (eg. mnc2nii in minc-tools) or by realigning images. Most brains of larger
% animals show a slight positive rotation (i.e. against the clock) that can 
% be seen quite good on the top view of both hemishpheres. The left hemi-
% sphere is a bit smaller in frontal regions but its occipital lob is a bit
% longer or moves a bit more to the right. Moreover the 
%
%  isflipped = cat_vol_testflipping(obj,Affine,method,stime,verb)
%  
% _________________________________________________________________________
% Robert Dahnke
% $Id$

  if nargin == 1; 
    %job    = varargin; 
  else
    obj    = varargin{1};
    Affine = varargin{2};
    method = varargin{3};
    if nargin>=4
      stime  = varargin{4};
    else
      stime  = clock; 
    end
    if nargin>=5
      verb   = varargin{5};
    else
      verb   = 0; 
    end
  end


  switch method
    case 1
      % simple fast test and warning if the x axis is positive that is
      % untypical for (structural) MRI.
      
      mati = spm_imatrix( Affine );  
      if sign( mati(7) ) == -1
        flipped    = 1; 
        flippedval = 1; 
        mid  = [mfilename 'cat_run_job:PossiblyFlippedInput'];
        msg  = sprintf(['Positive x coordinate detected. \\\\n' ...
          'Please check image for _possible_ side flipping. \\\\n'...
          'You can check for the expected surface rotation (against the clock). \\\\n' ... 
          'Use spm_image to correct flipping by changing the \\\\n' ...
          'sign of the x axis scaling from 1 to -1. ']); 
        cat_io_addwarning(mid,msg,2,[0 2]);
      else
        flipped    = 0;
        flippedval = 0; 
      end
      
    case 2
      %% test for flipping
      %  RD202009: IN DEVELOPMENT - not working in all cases
      %  The flipping test focuses on the following aspects
      %  (i)   better alignment to normal side-biased average template,
      %  (ii)  more neutral alignment of flipped data also after additional
      %        flipping, and  
      %  (iii) less shearing in unflipped data.
      % 
      %  Just draw and test the different cases on paper 
      %     _____       _____       _____
      %    /  \  \     /  /  \     /  |  \
      %   /   |   \   /   |   \   /   \   \
      %   \___\___/   \___/___/   \___|___/  
      %   unflipped    flipped     flipped
      %                            aligned

      
      % The test takes about 20 to 120 seconds and we need some command line output. 
      if nargin >= 4
        stime = cat_io_cmd('Test side flipping:','','',1,stime); 
      end
      
      % Apply affine registration to use a rigid transformation. 
      % Use only rounded parts of the initial affine registration to avoid
      % a too strong bias - at least reset the shearing component to 0 but
      % try to use the othere information to avoid miss-alignments and long
      % processing. Unclear how well this works
      affscale        = spm_imatrix(Affine); 
      affscale(1:3)   = fix(affscale(1:3) / 3  ) * 3;    % remove translation bias
      %affscale(4:9)   = fix(affscale(4:9) * 100) / 100;  % remove rotation bias 
      affscale(10:12) = 0;                               % remove shearing bias
      affscale        = spm_matrix(affscale);

      %obj.fwhm = 0; 
      obj.imageflipped(1)     = obj.image(1); 
      obj.imageflipped(1).mat = affscale * obj.image.mat;
      obj.image.mat           = affscale * obj.image.mat; 
      % flipped image
      obj.imageflipped(1).mat = [-1 0 0 0; 0 1 0 0 ; 0 0 1 0; 0 0 0 1] * obj.image(1).mat;
      
      % add artificial positive rotation that is worse for flipped images
      imat  = spm_imatrix(obj.image(1).mat);         imat(6) = imat(6) - 0.01; obj.image(1).mat        = spm_matrix(imat); 
      imat  = spm_imatrix(obj.imageflipped(1).mat);  imat(6) = imat(6) - 0.01; obj.imageflipped(1).mat = spm_matrix(imat); 

      % Estimate ideal rigid registration for both flipped and unflipped
      % It is expected that the unflipped images can achieve a higher
      % log-likelyhood (llu vs. llf)
      % This test is expecting a allready aligned images from cat_run_job
      % and use lower number of interation to increase speed.
      [Affineu,llu ]  = spm_maff8(obj.image(1)       ,obj.samp,obj.fwhm,obj.tpm,eye(4) ,'rigid',10);
      [Affinef,llf ]  = spm_maff8(obj.imageflipped(1),obj.samp,obj.fwhm,obj.tpm,eye(4) ,'rigid',10);
      
      % Use only one interation to just estimate the log-likelyhood.
      % For the additional flipped image it is expected that the alignment 
      % is better for a flipped image because its alignment is already a
      % compromise, whereas a unflipped images with match better to the 
      % template before but not after flipping. 
      [Affinebu,llbu] = spm_maff8(obj.imageflipped(1),obj.samp,obj.fwhm,obj.tpm,Affineu,'rigid',1); %#ok<ASGLU>
      [Affinebf,llbf] = spm_maff8(obj.image(1)       ,obj.samp,obj.fwhm,obj.tpm,Affinef,'rigid',1); %#ok<ASGLU>

      % registration compontents 
      mati            = spm_imatrix( Affineu );  
      matiflipped     = spm_imatrix( Affinef );  

      % main test variables
      flippedval(1)   = (sum(abs( mati(10:11) )) / sum(abs( matiflipped(10:11) )) - 1); % less shearing
      flippedval(2)   = ((llf  ./ llu ) - 1) * 100; % the correct image fits better to the TPM 
      flippedval(3)   = ((llbu ./ llbf) - 1) * 100; % if we flip the perfect align image the overlap of the correct image is _lower_
      
      % main variable with weighing and limitations
      flipped         = min(0.95,max(-.95, flippedval(1) / 3 +  ... 
                                           flippedval(2) * 2 +  ... 
                                           flippedval(3) / 2 ));
                                      
      % just print some values
      if verb
        disp([mati llu llbu 0 0 0;matiflipped llf llbf flippedval])
      end
        
      %% warning/error message
      %  Create the final message as note/warning/alert, depending on the 
      %  estimated probability of flipping. 
      if flipped>0.3 
        fprintf('\n'); 
        mstr = {'mostLikelyUNflippedInput','possiblyFlippedInput','probablyFlippedInput','mostLikelyFlippedInput'};
        fstr = {'most likely NOT','possibly','probably','most likely'};
        mid  = [mfilename 'cat_run_job:%s',mstr{ round(flipped * (numel(fstr) - 1) + 1 ) }];
        msg  = sprintf(['The image is %s flipped (probability %0.0f%%%%%%%%)! \\\\n'...
          'You can check the surface rotation (against the clock). \\\\n' ... 
          'Use spm_image to correct flipping by changing the \\\\n' ...
          'scaling of the x axis scaling from 1 to -1. '], ...
          fstr{ round(flipped * (numel(fstr) - 1) + 1 ) }, flipped * 100 ); 

        cat_io_addwarning(mid,msg,round(flipped * 3 - 1),[0 2],flippedval);
      end
    otherwise 
      error('Unknown method %s.',method)
  end
end
function [trans,reg,Affine] = cat_main_registration(job,dres,Ycls,Yy,Ylesion)
% ______________________________________________________________________
%  Spatial registration function of cat_main preprocessing that include
%  the SPM DARTEL and (optimized) SHOOTING registration approaches. 
%
%  There are 4 image spaces with image different dimensions and resolution: 
%    (1) individual/original  (IR)  - properties of the (interpolated) anatomical image 
%    (2) template             (TR)  - properties of the registration template image
%    (3) registration         (RR)  - properties for the registration 
%    (4) output/analyse       (AR)  - final resolution of the normalized images
%
%  DARTEL runs depending on the output resolution and RR = AR, whereas 
%  the original DARTEL runs always on the TR. It takes typically about  
%  3 Minutes/Subject for RR = 1.5 mm. 
%  SHOOTING needs about 10 Minutes/Subject for RR = 1.5 mm and much more 
%  for higher RR and especially the low smooth deformation are important. 
%  The optimized version used therefore reduced resolution level to speedup
%  the iteration and allow more smooth deformations and reduce the final 
%  costs. The changes between the deformation and the matching compared to 
%  the templates are used as iteration criteria.
%  In 2021 a boundary box parameter was added that allows to change the size
%  of the final processed images.
%  Moreover, a backup routine is available that use the given nonlinear 
%  registration from the input (normally from the unified segmentation).
%
%  Main control parameter:
%    job.extopts.regstr ..
%      * Main cases:
%        0  .. DARTEL, 
%        eps - 1  .. Optimized Shooting with low (eps; fast) to high
%                    to high quality (1; slow) with 0.5 as default
%        2  .. Optimized Shooting with fixed resolutions (3:(3-TR)/4:TR)
%        3  .. Optimized Shooting with fixed resolutions (TR/2:TR/4:TR)
%        4  .. Default Shooting
% 
%      * Fixed resolution level with interpolation:
%        11 .. hard deformations
%        12 .. medium deformations
%        13 .. soft deformations
%
%      * Fixed resolution level without interpolation (max res = TR):
%        21 .. hard deformations
%        22 .. medium deformations
%        23 .. soft deformations
%
%  Structure:
%    [trans,reg] = cat_main_registration(job,res,Ycls,Yy,Ylesion)
%   
%    trans .. output variable that is used in cat_main and cat_io_writenii
%    reg   .. output variable that include the summarize of the deformation
% 
%    job   .. SPM job structure with cat_main parameter
%    res   .. SPM parameter structure with further fields from CAT processing
%    Ycls  .. tissue classification that is used for deformation
%    Yy    .. old initial SPM deformation 
%    Ylesion  .. template files
% ______________________________________________________________________
%
%  The TR is expected to be between 1.0 and 1.5 mm. Higher resolution are 
%  possible, but not will not be become a standard CAT templates.
%
%  RRs higher than the TR allow small improvements, but with improper costs.
%  Tests with 0.5 mm RR lead to an shooting error (to high determinant).
%
%  Lower RR lead to smoother deformations and were the focus of the
%  optimization and allow further deformation cases.
%
%  There are different ways to use lower resolutions: 
%    (LR = lowest resolution - should be better than 3 mm!):
%    1) flexible step size with fixed lower limit: 
%       a) LR : (LR - TR)/5 : TR 
%    2) upper limit
%       a) fixed levels without/without limit by TR oder by user (deformation limit)
%          max( RR , 3.0 : -0.5 : 1.0 )
%          max( RR , 2.5 : -0.5 : 0.5 )
%       b) dynamic levels
%    3) TR limit
%       a)  TR+SS*5 : -SS : TR    with SS = 0.25 - 0.5
%    4) LR limit
%
%  Although, it is possible to use the deformation to estimate finer
%  rigid/affine transformation, we do not do this because rigid/affine 
%  output should run without spatial registration.
%
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id$

%#ok<*ASGLU>
  
  
  if ~nargin, help cat_main_registration; end

  
  
  % if there is a breakpoint in this file set debug=1 and do not clear temporary variables 
  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

  
  % Lesions
  if ~exist('Ylesion','var') || sum(Ylesion(:))==0, Ylesion = []; end
  
  
  % initialize output
  Affine        = dres.Affine;
  trans         = struct();
  reg           = struct();
  
  
  % ########################
  % default settings ############# mostly shooting ... not here 
  def.nits      = 64;             % Shooting iterations
  def.opt.ll1th = 0.001;          % smaller better/slower
  def.opt.ll3th = 0.002;          % 
  def.fast      = 0;              % no report, no output
  def.affreg    = 1;              % additional affine registration
  def.iterlim   = inf;            %
  def.regra     = 1;              % use affine registration: 0 - rigid, 1 - affine 
  def.clsn      = 2;              % use only GM and WM ... may extend to CSF ################### DEPENDING ON TEMPLATE #############
  if ~isfield(job.extopts,'reg'), job.extopts.reg = struct(); end
  if isstruct( job.extopts.reg )
    dreg        = cat_io_checkinopt(job.extopts.reg,def);
  else 
    dreg        = def; 
  end
  dres.Affine0  = dres.Affine; 
  % ########################
  
  
  % limit number of classes
  if dreg.clsn>0 && numel(Ycls)>dreg.clsn, Ycls(dreg.clsn + 1:end) = []; end
      
  
  
  % error in parallel processing - can't find shooting files (2016/12)
  % switch to shooting directory
  if ~exist('spm_shoot_defaults','file') || ~exist('spm_shoot_update','file') 
    olddir = cd; cd(fullfile(spm('dir'),'toolbox','Shoot'));
  end
  if job.extopts.subfolders, mrifolder = 'mri'; else, mrifolder = ''; end

  

  % this is just for me to create different templates and can be removed in a final version
  if numel(job.extopts.vox)>1 || numel(job.extopts.regstr)>1 || (isfield(job,'export') && job.export)
    job.extopts.multigreg = 1;
    export                = job.export;  % write files in sub-directories for tests
  else
    job.extopts.multigreg = 0;
    export                = 0; 
  end

  
  
% #################### THIS SHOULD BE EVALUATED IN FUTURE #################  
% additional affine registration of the GM segment (or the brainmask) ?
% only if not the prior (long-pipeline) is used
% RD202101 - this (GM registration) cause problems (e.g. BUSS_2020)
% - We need to rething if a final registration would be usefull because we 
%   now have a much better (or even final) brainmask as well as a (first)
%   segmentation.
% - GM registration would be better for unusal affine GM analysis and would 
%   reduce the costs of a typical GM focused deformation but it does not
%   fit to previous definitions rising different issues.
%   But an update of with new brain mask may be interesting, but the method
%   has to be choosen well - maybe the affreg of the brainmasks? or with T1
%   or the (smoothed) p0 Segmentation. maffreg will allways try to balance 
%   the different classes (what we probably not want)
%{
  if dreg.affreg && ( ~isfield(job,'useprior') || isempty(job.useprior) || ~exist(char(job.useprior),'file') )
    % Create maffreg obj structure similar to cat_run_job but only for the 
    % GM segment.
    obj.image         = res.image; 
    obj.image.pinfo   = repmat([255;0],1,size(Ycls{1},3));
    obj.image.dt      = [spm_type('UINT8') spm_platform('bigend')];
    obj.image.dat     = cat_vol_ctype(single(Ycls{1})/2 + single(Ycls{2}));  
    if isfield(obj.image,'private'), obj.image = rmfield(obj.image,'private'); end
    obj.samp          = 1.5; 
    obj.fwhm          = 1;
    if isfield(obj,'tpm'), obj = rmfield(obj,'tpm'); end

    % call maffreg
    obj.tpm = spm_load_priors8(res.tpm(1:2));
    
    wo = warning('QUERY','MATLAB:RandStream:ActivatingLegacyGenerators'); wo = strfind( wo.state , 'on');
    if wo, warning('OFF','MATLAB:RandStream:ActivatingLegacyGenerators'); end
    res.Affine = spm_maff8(obj.image,obj.samp,obj.fwhm,obj.tpm,res.Affine,job.opts.affreg,80);
    if wo, warning('ON','MATLAB:RandStream:ActivatingLegacyGenerators'); end
    Affine =  res.Affine;
    spm_progress_bar('Clear');
  end
%}
% #########################################################################
  
  

  % Helping function(s):
  % BB2dim estimate the new size of an image given by the old and new
  % boudnary box (oldbb/newbb), image resolution (oldres/newres), and the
  % size of the old image. MNI space has always a 0-slice resulting in odd
  % dimensions.
  BB2dim = @(oldbb,newbb,olddims,oldres,newres) floor( ( olddims * oldres - sum( ( abs(oldbb) - abs(newbb) ) ) ) / newres / 2 ) * 2 + 1;
         
  
                    
  % boundary box settings:    0 - Template, 1 - TPM, +[1x1] - SPM MNI + #; -[1x1] - (TMP>0.1) + #; [2x3] - own BB       
  % boundary box definition:  to write a correct x=-1 output image, we have to be sure that the x value of the bb is negative
  if dres.bb(1)<dres.bb(2), bbt=dres.bb(1); dres.bb(1)=dres.bb(2); dres.bb(2)=bbt; clear bbt; end
  if isfield( job.extopts , 'bb' ) 
    if numel( job.extopts.bb ) == 1 
      switch job.extopts.bb 
        case 1 % TPM
          resbb = dres.bb;
        case 0 % Template
          resbb = spm_get_bbox( dres.tmp2{1}(1) );
        otherwise
          if job.extopts.bb > 1 
            % dynamic boundary based on the human SPM TPM that requires a
            % similar lower limit (they just extend the BB to cover the
            % full head but do not extend below so strong the cerbellum)
            b     = job.extopts.bb; 
            resbb = [-72 -110 -72; -72 -72 103] + [ -b -b  0 ; b b b*2]; 
          else
            % it is possible to use the template to create an dynamic 
            % automatic mapping
            b     = -job.extopts.bb; 
% ######################## THIS NEED TO BE FINISHED #######################            
            %Ytmp  = 
            %tmpbb = 
% #########################################################################        
            tmpbb = [-72 -110 -72; -72 -72 103]; 
            resbb = tmpbb + [ -b -b -b ; b b b]; 
          end
      end
    elseif numel( job.extopts.bb ) == 6
      resbb = reshape( job.extopts.bb(:) , 2 , 3) ; % own
    else
      error('BB has to be a 2x3 matrix.'); 
    end
    if resbb(1)<resbb(2), bbt=resbb(1); resbb(1)=resbb(2); resbb(2)=bbt; clear bbt; end
  else
    resbb = dres.bb; 
  end
  
  
 
  Vtmp   = spm_vol(job.extopts.templates{1});
  % this is the main loop for different parameter
  for regstri = numel(job.extopts.regstr):-1:1
    for voxi = numel(job.extopts.vox):-1:1
      % set registration parameter
      [reg,res,job] = regsetup(dreg,dres,job,regstri,voxi);
      
      tpmM   = res.tpm(1).mat; 
      tmpM   = Vtmp(1).mat; 
  
      % resolutions:
      tmpres = abs(tmpM(1));                                                                   % template resolution 
     %regres = reg(regstri).opt.rres; if isinf(regres), regres = tmpres; end                   % registration resolution
      newres = job.extopts.vox(voxi); if isinf(newres), newres = tmpres; end                   % output resolution

      % image dimension 
      idim = res.image(1).dim(1:3);                                                            % (interpolated) input image resolution
      tdim = res.tmp2{1}(1).dim;                                                               % registration template image size
      odim = BB2dim(res.bb,resbb,res.tmp2{1}(1).dim,tmpres,newres);                            % output image size
      
      % mat matrices for different spaces
      % M0 for the individual volume
      % M1* for the different normalized spaces
      % Mad for the Dartel template with different BB as the TPM
      M0   = res.image.mat;                                                                    
      imat = spm_imatrix(tmpM);   imat(1:3) = imat(1:3) + imat(7:9).*(tdim - (newres/tmpres*odim))/2; imat(7:9) = imat(7:9) * newres/tmpres; M1   = spm_matrix(imat); 
      imat = spm_imatrix(eye(4)); imat(1:3) = imat(1:3) + imat(7:9).*(tdim - (newres/tmpres*odim))/2; imat(7:9) = imat(7:9) * newres/tmpres; Mad  = spm_matrix(imat);
      
      % affine and rigid parameters for registration 
      [M3,R]          = spm_get_closest_affine( affind(rgrid( idim ) ,M0) , affind(Yy,tpmM) , single(Ycls{1})/255); clear M3; % here the TPM space is used
      Mrigid          = M0\inv(R)*M1;                                                          % transformation from subject to registration space (rigid)
      Maffine         = M0\inv(res.Affine)*M1;                                                 % individual to registration space (affine)
      mat0a           = res.Affine\M1;                                                         % mat0 for affine output
      mat0r           = R\M1;                                                                  % mat0 for rigid ouput
      
      % settings for new boundary box of the output images 
      if isfield(res,'imagesc'); trans.native.Vo = res.imagec(1); else, trans.native.Vo = res.image0(1); end
      trans.native.Vi = res.image(1);
      trans.affine    = struct('odim',odim,'mat',M1,'mat0',mat0a,'M',Maffine,'A',res.Affine);  % structure for cat_io_writenii
      trans.rigid     = struct('odim',odim,'mat',M1,'mat0',mat0r,'M',Mrigid ,'R',R);           % structure for cat_io_writenii
      
      if debug && export
        % template creation
        fn = {'GM','WM','CSF'}; clsi=1; 
        cat_io_writenii(trans.native.Vo,single(Ycls{clsi})/255,fullfile(mrifolder,'debug'),sprintf('p%d',clsi),...
          sprintf('%s tissue map',fn{clsi}),'uint8',[0,1/255],[0 0 0 3],trans);
      end 
      
      % rigid vs. affine input in registration: 
      if dreg.regra % affine
        Maffinerigid = Maffine;
        TAR          = res.Affine;
      else % rigid
        Maffinerigid = Mrigid; 
        TAR          = R;
      end


      
      try 
        %  main functions
        if res.do_dartel
          if job.extopts.regstr(regstri)==0       % Dartel
            [trans,reg] = run_Dartel(Ycls,Ylesion,job,reg,res,trans,Mad,Maffinerigid,TAR,regstri,voxi);
          elseif job.extopts.regstr(regstri)>0    % Geodesic Shooting 
            [trans,reg] = run_Shooting(Ycls,Ylesion,job,reg,res,trans,Maffinerigid,regstri,voxi); 
          end
          
          % Report - display registration profil 
          if job.extopts.verb>1 || export
            report(job,reg,regstri);
          end
          
          % Export
          if export
            write_nii(Ycls,job,trans,reg.testfolder,reg(regstri));
          end
        end
      catch
%% ------------------------------------------------------------------------
%  Catch problems and use the original input (e.g. from the unified
%  segmentation) and prepare the output by just changing BB and vox
%  and prepare the output structure we use to write the nifties.
%  This uses the old deformation that was in general done by the Unified
%  Segmentation. 
%  ------------------------------------------------------------------------
    
        cat_io_addwarning('cat_main_registration:regError',...
          'Registration problem due to given input segmentation. \\nUse previous registration.',2,[1 1]);
        res.Affine = res.Affine0; 
        
        if ~exist('res','var'), res = dres; end
        tpmM   = res.tpm(1).mat; 
        
        % resolutions
        % - in this case tpmres == tmpres
        tpmres = abs(tpmM(1));                                                        % template resolution 
        newres = job.extopts.vox(voxi); if isinf(newres), newres = tpmres; end        % output resolution
 
        % image size/dimension 
        idim   = res.image(1).dim(1:3);                                               % (interpolated) input image size
        tdim   = res.tpm(1).dim;                                                      % tpm/tmp size 
        odim   = BB2dim(res.bb,resbb,res.tmp2{1}(1).dim,tpmres,newres);               % output size 
        
        % mat matrices for different spaces
        % - here M1 == M1t
        imat   = spm_imatrix(tpmM); imat(1:3) = imat(1:3) + imat(7:9).*(tdim - (newres/tpmres*odim))/2; imat(7:9) = imat(7:9) * newres/tpmres; 
        M1     = spm_matrix(imat);                                                    % warped matrix
        M0     = res.image.mat;                                                       % for (interpolated) individual volume
        
        % estimate the inverse transformation 
        yid             = spm_diffeo('invdef', Yy , odim, inv(tpmM\M1), M1\res.Affine*M0); 
        yi              = spm_diffeo('invdef', yid, idim, inv(M1\res.Affine*M0), eye(4)); clear yid; 
        yi2             = spm_diffeo('invdef', yi , odim, eye(4), eye(4)); 
        w               = max( eps , abs(spm_diffeo('def2det', yi2 ) ) ); 
        % Adaption to avoid boundary effects by correcting the voxel close
        % to the image boundary that are effected by the interpolation of
        % the field. 
        msk             = cat_vol_morph( isnan( yi2(:,:,:,1)) ,'d',3); w( msk ) = NaN; 
        wa              = cat_vol_approx(w,'nn',1,4); bg = cat_stat_nanmean( wa(msk(:)) ); 
        msk             = cat_vol_smooth3X(msk,2); w( isnan(w) ) = wa( isnan(w) ); 
        w               = wa .* msk + (1-msk) .* w; clear msk wa; 
        % use half registration resolution to define the amout of smoothing to reduce registration artefacts
        fs              = newres / 2; w = w - bg; spm_smooth(w,w,fs); w = bg + w; % spm smoothing boudary problem where values outside are 0 
       
        % affine and rigid parameters
        [M3,R]          = spm_get_closest_affine( affind(rgrid( idim ) ,M0) , affind(Yy,tpmM) , single(Ycls{1})/255); clear M3; 
        Mrigid          = M0\inv(R)*M1;                                        % transformation from subject to registration space
        Maffine         = M0\inv(res.Affine)*M1;                               % individual to registration space
        mat0r           = R\M1;                                                % mat0 for rigid ouput
        mat0a           = res.Affine\M1;                                       % mat0 for affine output
        M2              = inv(Maffine);                                        % warped - if we use the US than we may have to use rigid
        
        % final structure
        if isfield(res,'imagesc'); trans.native.Vo = res.imagec(1); else, trans.native.Vo = res.image0(1); end
        trans.native.Vi = res.image(1);
        trans.affine    = struct('odim',odim,'mat',M1,'mat0',mat0a,'M',Maffine);    % structure for cat_io_writenii
        trans.rigid     = struct('odim',odim,'mat',M1,'mat0',mat0r,'M',Mrigid);     % structure for cat_io_writenii
       %trans.warped    = struct('y',yid,               'odim',odim,'M0',M0,'M1',M1,'M2',M2,'dartel',res.do_dartel);          % simple version with push artefacts for tests 
        trans.warped    = struct('y',yi ,'yi',yi2,'w',w,'odim',odim,'M0',M0,'M1',M1,'M2',M2,'dartel',res.do_dartel,'fs',fs);  % nicer version with interpolation
        if job.output.jacobian.warped
          trans.jc      = struct('odim',odim,'dt2',w); 
        end
        
        % export for tests
        if export
          write_nii(Ycls,job,trans,sprintf('US_tr%3.1f_or%3.1f',tpmres,job.extopts.vox(1))); ping
        end
      end
    end
  end

  % back to old directory 
  if exist('olddir','var'), cd(olddir); end
end
%=======================================================================
function [trans,reg] = run_Shooting(Ycls,Ylesion,job,reg,res,trans,Maffinerigid,regstri,voxi)
%  ------------------------------------------------------------------------
%  Shooting
%  ------------------------------------------------------------------------

  % if there is a breakpoint in this file set debug=1 and do not clear temporary variables 
  dbs   = dbstatus; debug = 0; for dbsi=1:numel(dbs), if strcmp(dbs(dbsi).name,mfilename); debug = 1; break; end; end

  
% ####################### THIS NEED FURTHER CLEANUP #######################  
  
  % helping functions
  %BB2dim = @(oldbb,newbb,olddims,oldres,newres) floor( ( olddims * oldres - sum( ( abs(oldbb) - abs(newbb) ) ) )/newres/2)*2 + 1;
  getmm  = @(bb) [[bb(1,1) bb(2,1) bb(1,1) bb(2,1) bb(1,1) bb(2,1) bb(1,1) bb(2,1); 
                   bb(1,2) bb(1,2) bb(2,2) bb(2,2) bb(1,2) bb(1,2) bb(2,2) bb(2,2);
                   bb(1,3) bb(1,3) bb(1,3) bb(1,3) bb(2,3) bb(2,3) bb(2,3) bb(2,3)]; 
                   ones(1,8)]; % definition from the SPM registration function 

  if isempty( Ylesion ), clear Ylesion; end
                 
  n1     = reg(regstri).clsn;
  Vtmp   = spm_vol(job.extopts.templates{1});
  tmpM   = Vtmp(1).mat; 
  idim   = res.image(1).dim(1:3);
  odim   = trans.rigid.odim; 
  vx_vol = sqrt(sum(res.image(1).mat(1:3,1:3).^2));
  
  % resolutions:
  tmpres = abs(tmpM(1));                                                                   % template resolution 
  regres = reg(regstri).opt.rres; if isinf(regres), regres = tmpres; end                   % registration resolution
  newres = job.extopts.vox(voxi); if isinf(newres), newres = tmpres; end                   % output resolution

  M0     = res.image.mat; 
  M1     = trans.rigid.mat;
  imat   = spm_imatrix(tmpM); imat(7:9) = imat(7:9) * regres/tmpres; M1r = spm_matrix(imat); 
  
  
  
  % multiresolution parameter
  % this part may require further work 
  tempres2 = reg(regstri).opt.resfac * regres;  % registration resolution
  if numel(reg(regstri).opt.ll3th)~=numel(tempres2), reg(regstri).opt.ll3th = repmat(reg(regstri).opt.ll3th(1),numel(tempres2)); end


  tmpbb   = spm_get_bbox( res.tmp2{1}(1) ); % Template
  if tmpbb(1)<tmpbb(2), bbt=tmpbb(1); tmpbb(1)=tmpbb(2); tmpbb(2)=bbt; clear bbt; end
  mmtmp   = getmm(tmpbb);

  % registration space size
  rdims   = zeros([numel(reg(regstri).opt.resfac),3]); 
  for ri=1:numel(reg(regstri).opt.resfac), rdims(ri,:) = floor(res.tmp2{1}(1).dim * tmpres/regres / reg(regstri).opt.resfac(ri)); end
  
  % redi
  rred    = zeros([numel(reg(regstri).opt.resfac),3]); 
  for ri=1:numel(reg(regstri).opt.resfac), rred(ri,:) = (regres / reg(regstri).opt.resfac(ri)) ./ vx_vol; end

  imat      = spm_imatrix(tmpM); 
%    imat(1:3) = imat(1:3) + imat(7:9).*(odim - (tmpres/regres* rdims(ri,:)))/2; 
  imat(7:9) = imat(7:9) * regres/tmpres; M1rr=spm_matrix(imat); 
  vx2rr     = M1rr\mmtmp; 

  Mys       = cell(size(reg(regstri).opt.resfac));
  Mads      = cell(size(reg(regstri).opt.resfac));
  Mrregs    = cell(size(reg(regstri).opt.resfac)); 
  %%
  for ri=1:numel(reg(regstri).opt.resfac)
    vx3rr = ones(4,8); vx3rr([5,13,21,29]) = rdims(ri,1); vx3rr([10,14,26,30]) = rdims(ri,2); vx3rr([19,23,27,31]) = rdims(ri,3); % registration image
    vxtpm = tmpM\mmtmp; % registration image
    Mads{ri} = (tmpM\mmtmp)/vx3rr;
    if reg.regra
      if ri==1, mat0reg = res.Affine\M1rr * vx2rr/vxtpm; end 
      Mrregs{ri} = M0\inv(res.Affine)*M1rr*vx2rr/vx3rr;  
    else
      if ri==1, mat0reg = R\M1rr * vx2rr/vxtpm; end 
      Mrregs{ri} = M0\inv(R)*M1rr*vx2rr/vx3rr; %;    
    end
    if ri==1, Mys{ri} = eye(4); else,  Mys{ri} =  Mads{ri}/Mads{ri-1}; end
    %{
    if ri==1, 
      Mys{ri} = eye(4); 
    else, 
      Mys{ri}= eye(4); Mys{ri}(1:12) = Mys{ri}(1:12) * reg(regstri).opt.resfac(ri)/reg(regstri).opt.resfac(ri-1);  
      Mys{ri}(13:15) = -reg(regstri).opt.resfac(ri)/reg(regstri).opt.resfac(ri-1) / 2; 
    end;
    %}
  end
  % if saffine, Mrregs{1} = M0\inv(res.Affine)*M1rr*vx2rr/vx3rr; end
  %%
% ######################################################################### 

  


  if job.extopts.verb
    if reg(regstri).opt.stepsize>10^-3 || regres~=tmpres
      if job.extopts.regstr(regstri)>0 && job.extopts.regstr(regstri)<=1
        stime   = cat_io_cmd(sprintf('Optimized Shooting registration with %0.2f:%0.2f:%0.2f mm (regstr=%0.2f)',...
          tempres2(1),diff(tempres2(1:2)),tempres2(end),job.extopts.regstr(regstri))); 
      else
        stime   = cat_io_cmd(sprintf('Optimized Shooting registration with %0.2f:%0.2f:%0.2f mm',...
          tempres2(1),diff(tempres2(1:2)),tempres2(end))); 
      end
    else
      stime   = cat_io_cmd(sprintf('Default Shooting registration with %0.2f mm',reg(regstri).opt.rres)); 
    end
    fprintf('\n  Template: "%s"\n',job.extopts.templates{1})
  end    

  
% shooting parameter ... again ????
% ######################################################################### 
  sd = spm_shoot_defaults;           % load shooting defaults
  if reg.fast, reg(regstri).opt.nits = min(reg(regstri).opt.nits,5*reg.fast); end  % at least 5 iterations to use each tempalte
  if (job.extopts.regstr(regstri)>0 || regres~=tmpres) && job.extopts.regstr(regstri)~=4 
    %% need finer schedule for coarse to fine for ll3 adaptive threshold
    nits       = reg(regstri).opt.nits;        % default was 24  
    onits      = 24;                           % use John's interation number for the form and just interpolate the curve to modify different interation numbers
    lam        = 0.25;                         % Decay of coarse to fine schedule (default 0.5) ... update default to 0.25 to have smoother changes
    inter      = 32;                           % Scaling of parameters at first iteration
    sd.sched   = (inter-1) * exp(-lam*((1:(onits))-1) )+1; %sd.sched = sd.sched/sd.sched(end);
    sd.sched   = ((sd.sched - min(sd.sched)) / (max(sd.sched) - min(sd.sched)) * (inter-1) + 1); 
    sd.sched   = interp1(sd.sched,1:(numel(sd.sched)-1)/(nits-1):numel(sd.sched));  
    maxoil     = 8;                            % Maximum number of time steps for integration
    sd.eul_its = round((0:(nits-1))*(maxoil-0.5001)/(nits-1)+1); % Start with fewer steps
    nits       = numel(sd.sched)-1;            % Shooting iterations 
    tmpl_no    = floor(((1:nits)-1)/(nits-1)*(numel(res.tmp2)-0.51))+1; % Sort out which template for each iteration (default = round with more hr-iter)
  else
    nits       = reg(regstri).opt.nits;         
    tmpl_no    = round(((1:nits)-1)/(nits-1)*(numel(res.tmp2)-0.51))+1; 
  end
% ######################################################################### 





  %% The actual work
  % ---------------------------------------------------------------------
  it = 1; reg(regstri).dtc = zeros(1,5); ll  = zeros(1,2);
  while it<=nits 
    itime = clock;  

    if it==1 || (tmpl_no(it)~=tmpl_no(it-1)) 
      ti  = tmpl_no(it); 
      if debug && it>1, fo=f{1}; end %#ok<NASGU> % just for debugging    

      
      % load rigide/affine data
      f = cell(1,n1+1); f{n1+1} = ones(rdims(ti,1:3),'single'); 
      for k1=1:n1
        Yclsk1 = single(Ycls{k1}); 
        f{k1}  = zeros(rdims(ti,1:3),'single');
        spm_smooth(Yclsk1,Yclsk1,rred(tmpl_no(it)) * 0.1*max(0,size(rred,1) - tmpl_no(it))); % RD202101: Advanced Shooting: smoothing for resolution reduction and template level
        for i=1:rdims(ti,3)
          f{k1}(:,:,i) = single(spm_slice_vol(Yclsk1,Mrregs{ti}*spm_matrix([0 0 i]),rdims(ti,1:2),[1,NaN])/255); 
        end
        msk         = ~isfinite(f{k1});
        f{k1}(msk)  = 0;
        clear msk; 
      end
      if debug, fx = f{1}; end %#ok<NASGU> % just for debugging

      
      % template
      g = cell(1,n1+1); g{n1+1} = ones(rdims(ti,1:3),'single');
      for k1=1:n1
        g{k1} = zeros(rdims(ti,1:3),'single');
        tpm2k1 = res.tmp2{ti}(k1).private.dat(:,:,:,k1); 
        spm_smooth(tpm2k1,tpm2k1,repmat( tempres2(tmpl_no(it))/tmpres ,1,3)); % RD202101: Advanced Shooting: smoothing for resolution reduction but not template level
        for i=1:rdims(ti,3)
          g{k1}(:,:,i) = single(spm_slice_vol(tpm2k1,Mads{ti}*spm_matrix([0 0 i]),rdims(ti,1:2),[1,NaN]));
        end
        g{k1}(isnan(g{k1}(:))) = min(g{k1}(:)); % remove boundary interpolation artefact
        g{n1+1} = g{n1+1} - g{k1};
        if debug && k1==1, gx = g{1}; end %#ok<NASGU> % just for debugging
        g{k1} = spm_bsplinc(log(g{k1}), sd.bs_args);
      end
      g{n1+1} = log(max(g{n1+1},eps)); 

      
      
% -------------------------------------------------------------------------
%  RD202101: Advanced Shooting 
%            Seperate large ventricle 
%  ------------------------------------------------------------------------
%  In hydrocephalus with extremely large ventricles the outer part of the
%  ventriclar CSF aspire torwards the cranial CSF causing severe problems.
%  Hence, we try to seperate the ventricle as a seperate additional layer 
%  in segmentation and template.
%  This is only necessary for very large ventricle and the errors we made
%  are less important then the problems in the original setting.
%  However, a better ventricular segmentation should allow improvements.
%  ------------------------------------------------------------------------
      hydro = 0; 
      if job.extopts.regstr(regstri)~=4 
        if n1>2
          % use the CSF class directly
          Yven = cat_vol_morph( cat_vol_morph( cat_vol_morph(f{3},'l',[10 0.1]), 'o' ,1) ,'d'); 
        else
          % if CSF==BG than find the large central CSF region within the brain
          Ygwm = f{2} + f{1}; 

          % the brainmask requires large close in case of open ventricles that are connected to extra-ventricluar CSF 
          Yb   = cat_vol_morph(Ygwm,'ldc',12/tempres2(ti));
          Ybd  = cat_vol_smooth3X(Yb,6/tempres2(ti)); Ybd = Ybd/max(Ybd(:)); % get some distance to the skull

          % seperate the biggest part
          Yven = (Yb - Ygwm)>0.9 & Ybd>0.9; 
          Yven = cat_vol_morph(cat_vol_morph( cat_vol_morph( Yven ,'lo',1) ,'d',max(1,1/rred(ti))),'l'); 
          Yven = min(1,cat_vol_smooth3X(Yven,4/tempres2(ti))*1.5);  
          %Yven = Yven ./ max(Yven(:)); 
        end

        if sum(Yven(:))/sum(Yb(:))>0.1 % large ventricle
          hydro = 1; 
          if job.extopts.verb && ti==1
            %cat_io_cprintf([0 0 1], '  == Very large ventricles detected. Use seperate class and adapt schedule. == \n');  
            cat_io_addwarning([mfilename ':Hydrocephalus'],sprintf( ...
              'Very large ventricles detected (%0.2f%% of TIV). Use seperate class and adapt schedule. \\\\n', ...
              sum(Yven(:))/sum(Yb(:))>0.1),1,[0 1]);   
          end
        

%  ------------------------------------------------------------------------
%  RD202101: Advanced Shooting: 
%  ------------------------------------------------------------------------
%  Update schedule to have again stronger deformations to improve the 
%  addaption to the new template and resolution. 
%  This seams to help especially in worse cases with large ventricles.
%  ------------------------------------------------------------------------
          if it > 1 
            % at least the number of iterations for this template but maximal 4 
            fac = min( (max(tmpl_no) - tmpl_no(it) + 2) * 2, sum(tmpl_no==tmpl_no(it)) ); 
            fac = log(fac) / log(sd.sched(it));
            sd.sched = max( sd.sched , sd.sched.^fac ); 
          end

          % seperate template ventricle
          % no labopen because the ventricle is to small
          Ygven     = exp(g{3}) .* cat_vol_smooth3X(Yven,2/tempres2(ti)); 
          Ygven     = cat_vol_morph( cat_vol_morph( cat_vol_morph(Ygven>0.5,'l') ,'d',max(1,1/rred(ti))),'l');  
          Ygven     = exp(g{3}) .* cat_vol_smooth3X(Ygven,4/tempres2(ti));
          % moreover ... 
          if 0
            Ygvenoc   = Ygven; 
            Ygven     = min(1,Ygven / mean(Ygven(Ygven(:)>0.2))) * 0.99;
            Ygvenoc   = Ygven - Ygvenoc; 
          end

          % seperate individual ventricle 
          f{n1+2} = ones(rdims(ti,1:3),'single');
          f{n1+1} = max(eps,Yven - Ygwm);

          g{n1+2} = double(log(max(eps,exp(g{3}) - Ygven))); % temove ventricle from background / CSF
          g{n1+1} = double(log(max(eps,Ygven))); % own class
        end
        if ~debug, clear Yven Ygven Ygwm Yb Ybd; end
      end
      
      
      % set lesions area by template tissues to avoid transformation
      if exist('Ylesion','var') && sum(Ylesion(:))>0
        if job.extopts.verb && ti==1
          cat_io_cprintf([0 0 1], '  == Use lesion masking to limit deformations to non-masked areas == \n');  
        end

        Yclsk1 = single(Ylesion); 
        if reg(regstri).opt.resfac(ti)>1, spm_smooth(Yclsk1,Yclsk1,repmat((reg(regstri).opt.resfac(ti)-1) * 2,1,3)); end
        ls = zeros(rdims(ti,1:3),'single');
        for i=1:rdims(ti,3)
          ls(:,:,i) = single(spm_slice_vol(Yclsk1,Mrregs{ti}*spm_matrix([0 0 i]),rdims(ti,1:2),[1,NaN])); 
        end
        ls(isnan(ls(:))) = min(ls(:)); 
        for k1=1:n1
          tpm2k1 = res.tmp2{ti}(k1).private.dat(:,:,:,k1); 
          t = zeros(rdims(ti,1:3),'single');
          for i=1:rdims(ti,3)
            t(:,:,i) = single(spm_slice_vol(tpm2k1,Mads{ti}*spm_matrix([0 0 i]),rdims(ti,1:2),[1,NaN]));
          end
          t(isnan(t(:))) = min(t(:)); 
          f{k1} = f{k1} .* (1-ls) + t .* ls;
        end
      end
      
      
      % set background (and 
      for k1=1:numel(f)-1
        f{end}      = f{end} - f{k1}; 
        msk         = ~isfinite(f{k1});
        f{end}(msk) = 0.00001;
      end
      

      % loading segmentation and creating of images vs. updating these maps
      ll  = zeros(1,2);
      if it==1
        % create shooting maps
        y   = affind( squeeze( reshape( affind( spm_diffeo('Exp',zeros([rdims(ti,:),3],'single'),[0 1]), ...
              mat0reg), [rdims(ti,:),1,3] ) ) , inv(mat0reg)); clear def;                          %#ok<NASGU>  % deformation field
        u   = zeros([rdims(ti,:) 3],'single');                                                     %#ok<PREALL> % flow field
        dt  = ones(rdims(ti,:),'single');                                                          %#ok<PREALL> % jacobian
      elseif any(rdims(ti,:)~=rdims(ti-1,:))
        % updates only for changed resolutions

        if 0 
            %% just test code to verify the mapping/interpolation between resolution levels
            fx = cell(0); 
            for k1=1:numel(fo)
              for i=1:rdims(ti,3)
                fx{k1}(:,:,i) = single(spm_slice_vol(fo{k1},Mys{ti}*spm_matrix([0 0 i]),...
                  rdims(ti,1:2),[1,NaN])) / Mys{ti}(1); 
              end
            end
            fx(~isfinite(fx))=eps;
            ds('d2sm','',1,f{1},(fx{1} - f{1})/2+0.5,30), cat_stat_nansum(abs(fx{1}(:) - f{1}(:)))    
        end
        
        % size update u - flow field
        uo = u; 
        u  = zeros([rdims(ti,:) 3],'single');
        for k1=1:3
          for i=1:rdims(ti,3)
            u(:,:,i,k1) = single(spm_slice_vol(uo(:,:,:,k1),Mys{ti}*spm_matrix([0 0 i]),...
              rdims(ti,1:2),[1,NaN])) / Mys{ti}(1); % (tempres(ti) / tempres(ti-1))^2; % adapt for res 
          end
        end
        u(~isfinite(u))=eps;
        if ~debug, clear uo; end
        
        % update the deformation y and the determinatent by re-estimation 
        % rather than interpolation of low resolution data
        [y,J] = spm_shoot3d(u,prm,int_args); 
        dt    = spm_diffeo('det',J); clear J %#ok<NASGU>
      
      end
    end



    % More regularisation in the early iterations, as well as a less accurate approximation in the integration.
    % No, similar regularisation works in our case better and avoid to trap into local maxima.  
    vxreg    = repmat(reg(regstri).opt.vxreg,1,3);  % repmat(tempres(ti)^3,1,3)
    prm      = [vxreg, sd.rparam * sd.sched(it) * prod(vxreg)]; 
    int_args = [sd.eul_its(it), sd.cyc_its]; 

    % Gauss-Newton iteration to re-estimate deformations for this subject
    if job.extopts.verb
      if reg(regstri).opt.stepsize<=10^-3
        cat_io_cprintf(sprintf('g%d',5+2*(it==1 || (tmpl_no(it)~=tmpl_no(it-1)))),sprintf('% 5d |',it));
      else
        cat_io_cprintf(sprintf('g%d',5+2*(it==1 || (tmpl_no(it)~=tmpl_no(it-1)))),sprintf('% 5d | %0.2f |',it,tempres2(ti)));
      end
    end

    llo=ll; 

    [txt,u,ll(1),ll(2)] = evalc('spm_shoot_update(g,f,u,y,dt,prm,sd.bs_args,sd.scale)'); 
    [y,J] = spm_shoot3d(u,prm,int_args); 
    dt    = spm_diffeo('det',J); clear J
    extdt = max([dt(:)', 1/max(eps,dt(isfinite(dt(:))))]); 
    if job.extopts.verb
      if job.extopts.expertgui>1
        cat_io_cprintf(sprintf('g%d',5+2*(it==1 || (tmpl_no(it)~=tmpl_no(it-1)))),sprintf('%7.4f%8.4f%8.4f | %8.4f | %3.1f', ...
          ll(1)/numel(u), ll(2)/numel(u), (ll(1)+ll(2))/numel(u), sd.sched(it), log10(extdt) ) );
      else
        cat_io_cprintf(sprintf('g%d',5+2*(it==1 || (tmpl_no(it)~=tmpl_no(it-1)))),sprintf('%7.4f%8.4f%8.4f | %8.4f ', ...
          ll(1)/numel(u), ll(2)/numel(u), (ll(1)+ll(2))/numel(u), sd.sched(it) ) );
      end
    end

    wit = 0;
    if job.extopts.regstr(regstri)~=4 && ( extdt>20 || any(~isfinite(dt(:))) ) 
    % ---------------------------------------------------------------------
    % RD202101: Advanced Shooting 
    %           Correction of local hot spots by smoothing
    % ---------------------------------------------------------------------
    % In some cases but in particular hydrocephalus the deformation runs 
    % localy into some minima and get suck at some anatomical features 
    % and the determinant increase stronly, pointing to problematic areas. 
    % However, a very high/low determinat will also occure as wanted in 
    % in the ventricular regions. 
    % The idear is now to detect such regions and filter them lightly.
    % Alhough this will result in smoother deformations and lower overlap
    % it will at least support that we become a result and do not run
    % into a servere determinate error (with inf or nan)
    % It is also better to have a smoother deformation that does not fit 
    % than strong local outliers. 
    % ---------------------------------------------------------------------

      vr   = 3/tempres2(ti); 
      mxu  = dt./cat_vol_smooth3X(dt,vr) > 1.5;                       % we can use local increase of the determinate to detect local outlier
      mxu  = cat_vol_smooth3X(mxu,2*vr);                              % add some neighborhood
      for i=1:3, u(:,:,:,i) = cat_vol_laplace3R(u(:,:,:,i)/10,mxu>0.01, 0.1)*10; end  % filtering of the deformation 
      for i=1:3, ux = u(:,:,:,i);  uxs = cat_vol_smooth3X(ux,max(0.2,min(2,sd.sched(it)/8))); u(:,:,:,i) = uxs.*mxu + ux.*(1-mxu); end % filtering of the deformation 

      % re-estimate the determinat
      [y,J] = spm_shoot3d(u,prm,int_args); 
      dt    = spm_diffeo('det',J); clear J

      extdt = max([dt(:)', 1/max(eps,dt(isfinite(dt(:))))]); 
      if job.extopts.verb && job.extopts.expertgui>1, fprintf(' > %0.1f', log10(extdt) ); end

      dtl = 60;
      while ( any(~isfinite(dt(:))) && wit<20 ) || ... % values that we have to filter (whatever it takes)
            ( extdt>dtl             && wit<4  )        % regions that we should fitler

        % creating a filter mask
        mxu  = dt./cat_vol_smooth3X(dt,1*vr) > 1.5;    % we can use local increase of the determinate to detect local outlier
        mxu  = cat_vol_smooth3X(mxu,2*vr) > 0.01;      % add some neighborhood
        mxu2 = dt>(dtl*2) | dt<1/(dtl*2);           % critical regions with high determinant
        for i=1:3, mxu2 = mxu2 + abs(cat_vol_div((u(:,:,:,1)))); end  % local divergence in u is also a good indicator
     %   for i=1:3, mxu2 = mxu2 + ( cat_vol_div(u(:,:,:,1))) + ( cat_vol_div(-u(:,:,:,1))); end % local divergence in u is also a good indicator
        mxu  = ( mxu | mxu2>(1 - exp(g{3})) ) & exp(g{3})<0.5;           % prefare filtering of regions out of the ventricles were we expect strong defs
        mxu  = mxu | ~isfinite(dt) | isnan(dt)| dt>1000 | dt<1/1000;     % regions we have to filter anyway
   %     mxu  = cat_vol_smooth3X(mxu,4) > 0.1;                            % add some neighborhood

        % filtering of the deformation 
        for i=1:3
          if wit>0, u(:,:,:,i) = cat_vol_smooth3X(u(:,:,:,i),1); end     % global filtering in the worst-case
          u(:,:,:,i) = cat_vol_laplace3R(u(:,:,:,i)/10,mxu , 0.02)*10;   % local filtering 
        end

        % parameter updates 
        sd.sched = sd.sched.^1.1;  % + min(0.3,0.4*(1-it/nits)));      % 1.2
        prm      = [vxreg, sd.rparam * sd.sched(it+1) * prod(vxreg)]; 

        % re-estimate the determinat
        [y,J] = spm_shoot3d(u,prm,int_args); 
        dt    = spm_diffeo('det',J); clear J

        wit   = wit + 1;    
        extdt = max([dt(:)', 1/max(eps,dt(isfinite(dt(:))))]); 
        if job.extopts.verb && job.extopts.expertgui>1, fprintf(' > %0.1f', log10(extdt) ); end
      end
    end
    if job.extopts.verb, fprintf('\n'); end


    % save iteration parameter for later analysis
    if it==1 || (tmpl_no(it)~=tmpl_no(it-1)) 
      reg(regstri).ll(ti,1:4)  = [ll(1)/numel(dt) ll(2)/numel(dt) (ll(1)+ll(2))/numel(dt) ll(2)]; 
      dtx = dt; 
      dtx(dtx>eps & dtx<1)     = 1./dtx(dtx>eps & dtx<1); 
      reg(regstri).dtc(ti)     = mean(abs(dtx(:)-1)); 
      reg(regstri).rmsdtc(ti)  = mean((dtx(:)-1).^2).^0.5;
      dtg = cat_vol_grad(single(dtx)); 
      reg(regstri).rmsgdt(ti)  = mean((dtg(:)).^2).^0.5;
      clear dtx;
      clear dtg; 
    end

    % default Shooting error detection
    if any(~isfinite(dt(:)) | dt(:)>100 | dt(:)<1/100)
      cat_io_cprintf('err',sprintf('Problem with Shooting (dets: %g .. %g)\n', min(dt(:)), max(dt(:)))); %it=nits;
    end

    % avoid unneccessary iteration
    if job.extopts.regstr(regstri)>0 && job.extopts.regstr(regstri)~=4 && wit==0 && ~hydro && ll(1)<0.075 && ...
        ( ti>1 || (ti==1 && ll(1)/numel(u)<1 && ll(1)/max(eps,llo(1))<1 && ll(1)/max(eps,llo(1))>(1-0.01) )) && ...
        ( ll(1)/numel(u)<1 && ll(1)/max(eps,llo(1))<1 && ll(1)/max(eps,llo(1))>(1-reg(regstri).opt.ll1th) )
      it = max(it+1,find([tmpl_no,nits]>tmpl_no(it),1,'first')); 
      reg(regstri).ll(ti,1:4) = [ll(1)/numel(dt) ll(2)/numel(dt) (ll(1)+ll(2))/numel(dt) ll(2)]; 
      reg(regstri).dtc(ti) = mean(abs(dt(:)-1)); 
    else
      it = it+1; 
    end
  end

  % save some parameter for later ..
  dtx = dt; 
  dtx(dtx>eps & dtx<1)       = 1./dtx(dtx>eps & dtx<1); 
  reg(regstri).rmsdt         = mean((dtx(:)-1).^2).^0.5; 
  reg(regstri).dt            = mean(abs(dtx(:)-1));
  reg(regstri).dtc(ti+1)     = mean(abs(dtx(:)-1)); 
  reg(regstri).rmsdtc(ti+1)  = mean((dtx(:)-1).^2).^0.5; 
  dtg = cat_vol_grad(single(dtx)); 
  reg(regstri).rmsgdt(ti+1)  = mean((dtg(:)).^2).^0.5;
  clear dtg; 
  reg(regstri).ll(ti+1,1:4)  = [ll(1)/numel(dt) ll(2)/numel(dt) (ll(1)+ll(2))/numel(dt) ll(2)]; 
  clear dt1; 

  % preparte output
  if job.extopts.verb
    if job.extopts.regstr(regstri)==0
      cat_io_cmd(sprintf('Dartel registration with %0.2f mm takes',tempres(1)));
    elseif reg(regstri).opt.stepsize>10^-3  
      cat_io_cmd(sprintf('Shooting registration with %0.2f:%0.2f:%0.2f mm takes',tempres2(1),diff(tempres2(1:2)),tempres2(end))); 
    else
      cat_io_cmd(sprintf('Shooting registration with %0.2f mm takes',tempres2(1))); 
    end
    itime = cat_io_cmd(sprintf('  Prepare output'),'','',job.extopts.verb,stime);
  end

  
  
  %% Modulation using spm_diffeo and push introduces aliasing artefacts,
  % thus we use the def2det function of the inverted deformations to obtain the old and 
  % in my view a more appropriate jacobian determinant 
  % The 2nd reason to use the old modulation is compatibility with cat_vol_defs.m
  yi  = spm_diffeo('invdef', y  , idim, Maffinerigid * inv(M1r\M1), inv(M1r\M1)); %#ok<MINV>
  yi2 = spm_diffeo('invdef', yi , odim, eye(4), eye(4)); 
  w   = max( eps , abs(spm_diffeo('def2det', yi2 ) ) ); 

  
  % avoid boundary effects that are not good for the global measurements 
  % use half registration resolution to define the amout of smoothing to reduce registration artefacts
  vxd = M1(1)./M1r(1); 
  msk = cat_vol_morph( isnan( yi2(:,:,:,1)) ,'d',2/vxd); w( msk ) = NaN; 
  wa  = cat_vol_approx(w,'nn',1,4); bg = cat_stat_nanmean( wa(msk(:)) ); 
  msk = cat_vol_smooth3X(msk,1/vxd); w( isnan(w) ) = wa( isnan(w) ); 
  w   = wa .* msk + (1-msk) .* w; clear msk wa; 
  fs  = newres / 2; w = w - bg; spm_smooth(w,w,fs); w = bg + w; % spm smoothing boudary problem

  % yi2 for fast high quality output
  trans.warped = struct('y',yi,'yi',yi2,'w',w,'odim',odim,'M0',M0,'M1',M1,'M2',inv(Maffinerigid),'dartel',2,'fs',fs);

  if job.output.jacobian.warped
    % RD202101: not sure if the transfomration is correct ...
    uo  = zeros([odim 3],'single');
    for k1=1:3
      for i=1:odim
        uo(:,:,i,k1)  = single(spm_slice_vol(u(:,:,:,k1),(M1r\M1)*spm_matrix([0 0 i]),odim(1:2),[1,NaN]));
      end
    end
    uo(~isfinite(uo)) = eps;
    trans.jc = struct('u',uo,'odim',odim,'dt2',w); 
  end 

  
  if job.extopts.verb
    cat_io_cmd('','','',job.extopts.verb,itime); 
    cat_io_cmd(' ','',''); cat_io_cmd('','','',job.extopts.verb,stime); 
  end

  
  if debug
  % just for fast internal tests
    write_nii(Ycls,job,trans,'debug');
  end
end
%=======================================================================
function [trans,reg] = run_Dartel(Ycls,Ylesion,job,reg,res,trans,Mad,Maffinerigid,TAR,regstri,voxi)
%  ------------------------------------------------------------------------
%  Dartel spatial normalization to given template
%  ------------------------------------------------------------------------

  if isempty( Ylesion ), clear Ylesion; end

  Vtmp   = spm_vol(job.extopts.templates{1});
  tmpM   = Vtmp(1).mat; 
  
  M0   = res.image.mat; 
  M1   = trans.rigid.mat;
 
  idim = res.image(1).dim(1:3);
  odim = trans.rigid.odim; 
  rdim = trans.rigid.odim;
  
  % resolutions:
  tmpres = abs(tmpM(1));                                                                   % template resolution 
  regres = reg(regstri).opt.rres; if isinf(regres), regres = tmpres; end                   % registration resolution
  newres = job.extopts.vox(voxi); if isinf(newres), newres = tmpres; end                   % output resolution

  stime = cat_io_cmd(sprintf('Dartel registration with %0.2f mm on a %0.2f mm Template',newres,tmpres)); 
  fprintf('\n  Template: "%s"\n',job.extopts.templates{1})

  reg(regstri).opt.rres = newres; % 

  % dartel parameter 1
  rform = 0;          % regularization form: 0 - Linear Elastic Energy
  code  = 2;          % multinomial
  lmreg = 0.01;       % LM regularization
  cyc   = 3;          % cycles
  its   = 3;          % relaxation iterations (inner iteration)
  n1    = reg.clsn;   % use GM/WM for dartel
  if reg.fast, its = min(3,max(1,min(its,reg.fast))); end % subiteration

  % rparam .. regularization parameters: mu, lambda, id
  % K      .. time steps
  param = struct('K',{0 0 1 2 4 6},'its',its, ...
    'rparam',{[4 2 1e-6],[2 1 1e-6],[1 0.5 1e-6],[0.5 0.25 1e-6],[0.25 0.125 1e-6],[0.25 0.125 1e-6]});

  % initialize varibles and load anatomical image in registration space 
  f = zeros([rdim(1:3) 2],'single');
  g = zeros([rdim(1:3) 2],'single');
  u = zeros([rdim(1:3) 3],'single');
  for k1=1:n1
    for i=1:rdim(3)
      f(:,:,i,k1) = single(spm_slice_vol(single(Ycls{k1}),Maffinerigid*spm_matrix([0 0 i]),rdim(1:2),[1,NaN])/255);
    end
  end

  if exist('Ylesion','var') && sum(Ylesion(:))>0
    Ylesion = single(Ylesion); 
    ls = zeros(rdim(1:3),'single');
    for i=1:rdim(3)
      ls(:,:,i) = single(spm_slice_vol(single(Ylesion),Maffinerigid*spm_matrix([0 0 i]),rdim(1:2),[1,NaN]));
    end
  end

  % iterative processing
  % ---------------------------------------------------------------------
  it0 = 0;  % main iteration number for output
  reg(regstri).dtc = zeros(1,6);
  for it = 1:numel(param)
    prm   = [rform, param(it).rparam, lmreg, cyc, its, param(it).K, code];
    % load new template for this iteration
    for k1=1:n1
      for i=1:rdim(3)
        g(:,:,i,k1) = single(spm_slice_vol(res.tmp2{it}(k1),Mad*spm_matrix([0 0 i]),rdim(1:2),[1,NaN]));
      end

      if exist('Ylesion','var') && sum(Ylesion(:))>0
        f(:,:,:,k1) = f(:,:,:,k1) .* (1-ls) + g(:,:,:,k1) .* ls;
      end
    end

    for j = 1:param(it).its
      it0 = it0 + 1;
      [u,ll] = dartel3(u,f,g,prm);
      reg(regstri).lld(it0,:)  = ll ./ [prod(rdim) prod(rdim) newres^3]; 
      reg(regstri).lldf(it0,:) = [reg(regstri).lld(it0,1) / prod(regres),ll(1),ll(2),ll(1)+ll(2),ll(3)];
      cat_io_cprintf(sprintf('g%d',5+2*(mod(it0,its)==1)),...
        sprintf('% 5d | %6.4f | %8.0f %8.0f %8.0f %8.3f \n',it0,reg(regstri).lldf(it0,:)));

      if it0==1 % simplified! use values after first iteration rather than before
        [y0, dt] = spm_dartel_integrate(reshape(u,[odim(1:3) 1 3]),[1 0], 6); clear y0; 
        reg(regstri).ll(1,:)    =  reg(regstri).lldf(it0,:); 
        dtx = dt; 
        dtx(dtx>eps & dtx<1)    = 1./dtx(dtx>eps & dtx<1); 
        reg(regstri).dtc(1)     =  mean(abs(dtx(:)-1)); 
        reg(regstri).rmsdtc(1)  =  mean((dtx(:)-1).^2).^0.5;
        dtg = cat_vol_grad(single(dtx)); 
        reg(regstri).rmsgdt(1)  = mean((dtg(:)).^2).^0.5;
        clear dtg dt dtx; 
      end
    end 

    [y0, dt] = spm_dartel_integrate(reshape(u,[odim(1:3) 1 3]),[1 0], 6); clear y0; 
    reg(regstri).ll(it+1,:)    =  reg(regstri).lldf(it0,:); 
    reg(regstri).dtc(it+1)     =  mean(abs(dt(:)-1)); 
    reg(regstri).rmsdtc(it+1)  =  mean((dt(:)-1).^2).^0.5;
    dtg = cat_vol_grad(single(dt)); 
    reg(regstri).rmsgdt(it+1)  = mean((dtg(:)).^2).^0.5;
    clear dtg; 
  end
  reg(regstri).rmsdt         = mean((dt(:)-1).^2).^0.5; 
  reg(regstri).dt            = mean(abs(dt(:)-1));
  clear f g;


  % jacobian 
  if job.output.jacobian.warped || (isfield(job.extopts,'multigreg') && job.extopts.multigreg)
    if any(odim ~= idim)
      uo  = zeros([odim 3],'single');
      for k1=1:3
        for i=1:odim
          uo(:,:,i,k1)  = single(spm_slice_vol(u(:,:,:,k1)  ,spm_matrix([0 0 i]),odim(1:2),[1,NaN]));
        end
      end
      uo(~isfinite(uo))=eps; %dto(~isfinite(dto))=eps;
    else
      uo = u; 
    end
    if job.extopts.bb
      trans.jc = struct('u',uo,'odim',idim); 
    end
    clear uo
  end

  % deformation
  y0      = spm_dartel_integrate(reshape(u,[rdim(1:3) 1 3]),[0 1], 6); clear u;
  prm     = [3 3 3 0 0 0];
  Coef    = cell(1,3);
  Coef{1} = spm_bsplinc(y0(:,:,:,1),prm);
  Coef{2} = spm_bsplinc(y0(:,:,:,2),prm);
  Coef{3} = spm_bsplinc(y0(:,:,:,3),prm);
  clear y0;
  [t1,t2] = ndgrid(1:idim(1),1:idim(2),1); t3 = 1:idim(3);
  Yy = zeros([idim 3],'single');
  for z=1:idim(3)
    [t11,t22,t33] = defs2(Coef,z,Maffinerigid,prm,t1,t2,t3);
    Yy(:,:,z,1) = t11;
    Yy(:,:,z,2) = t22;
    Yy(:,:,z,3) = t33;
  end
  clear Coef t1 t2 t3 t11 t22 t33 z
  M = eye(4);
  for i=1:size(Yy,3)
    t1          = Yy(:,:,i,1);
    t2          = Yy(:,:,i,2);
    t3          = Yy(:,:,i,3);
    Yy(:,:,i,1) = M(1,1)*t1 + M(1,2)*t2 + M(1,3)*t3 + M(1,4);
    Yy(:,:,i,2) = M(2,1)*t1 + M(2,2)*t2 + M(2,3)*t3 + M(2,4);
    Yy(:,:,i,3) = M(3,1)*t1 + M(3,2)*t2 + M(3,3)*t3 + M(3,4);
  end
  clear t1 t2 t3 M; 

  
  % Modulation using spm_diffeo and push introduces aliasing artefacts,
  % thus we use the def2det function of the inverted deformations to obtain the old and 
  % in my view a more appropriate jacobian determinant 
  % The 2nd reason to use the old modulation is compatibility with cat_vol_defs.m
  yi2 = spm_diffeo('invdef' , Yy, odim,  M1\M1, eye(4)); 

  w   = max( eps , abs(spm_diffeo('def2det', yi2 ) ) ); % .* prod( sqrt(sum( M1(1:3,1:3).^2))); 
  % avoid boundary effects that are not good for the global measurements 
  w(:,:,[1 end]) = NaN; w(:,[1 end],:) = NaN; w([1 end],:,:) = NaN;
  % use half registration resolution to define the amout of
  % smoothing to reduce registration artefacts
  fs = newres / 2;
  spm_smooth(w,w,fs);

  
  trans.warped = struct('y',Yy,'yi',yi2,'w',w,'odim',odim,'M0',M0,'M1',M1,'M2',M1\TAR*M0,'dartel',1,'fs',fs);


  if job.extopts.verb<1, fprintf(sprintf('%s',repmat('\b',1,it0*47 + 2))); fprintf('\n'); end
  fprintf('%s %4.0fs\n',repmat(' ',1,66),etime(clock,stime)); 
end
%=======================================================================
%################################# NEED CLEANUP
function [reg,res,job] = regsetup(reg,res,job,regstri,voxi)
% Definition of the registration parameters reg of the i-th registration 
% regstri.

  
  tpmM   = res.tpm(1).mat; 

  if regstri < numel(job.extopts.regstr) && job.extopts.verb && ...
    (numel(job.extopts.regstr) || numel(job.extopts.vox)) 
    fprintf('\n\n'); 
  end 


  % set dartel/shooting templates
  if job.extopts.regstr(regstri)==0
    job.extopts.templates = job.extopts.darteltpms;
  else
    job.extopts.templates = job.extopts.shootingtpms;
  end
  res.tmp2 = cell(1,numel(job.extopts.templates)); 
  Vtmp = spm_vol(job.extopts.templates{1}); tmpM = Vtmp(1).mat; 
  n1   = max(min(2,numel(Vtmp)),numel(Vtmp)-1);  % use GM and WM for shooting


  % registration main parameter
  lowres                     = 2.5;                   % lowest resolution .. best between 2 and 3 mm 
  tpmres                     = abs(tpmM(1));          % TPM resolution 
  tempres                    = abs(tmpM(1));          % template resolution 
  reg(regstri).opt.nits      = reg.nits;              % registration iteration (shooting default = 24)
  reg(regstri).opt.vxreg     = tpmres;                % regularisation parameter that original depend on the template resolution
  reg(regstri).opt.rres      = tempres;               % final registration resolution 
  reg(regstri).opt.stepsize  = (lowres - reg(regstri).opt.rres)/4;  % stepsize of reduction 
  reg(regstri).opt.resfac    = (lowres : -reg(regstri).opt.stepsize : reg(regstri).opt.rres) / reg(regstri).opt.rres; % reduction factor 
  reg(regstri).opt.ll1th     = reg.opt.ll1th;                 % smaller better/slower
  reg(regstri).opt.ll3th     = reg.opt.ll3th;                 % smaller better/slower 
  reg(regstri).opt.regstr    = job.extopts.regstr;   
  reg.fast                   = reg.iterlim;                  % limit iterations per template level to test if processing work in principle 
  

  if job.extopts.regstr(regstri)==0
  % Dartel
    res.do_dartel            = 1; 
    reg(regstri).opt.rres    = tempres; %job.extopts.vox(voxi);
  elseif job.extopts.regstr(regstri)>0 && job.extopts.regstr(regstri)<=1
  % Optimized Shooting - manual limit 
    reg(regstri).opt.stepsize  = (lowres - tempres)/4;  % stepsize of reduction 
    reg(regstri).opt.resfac    = (lowres : -reg(regstri).opt.stepsize : tempres) / tempres; % reduction factor 

    reg(regstri).opt.ll1th     = 0.0010 + 0.10*(1-job.extopts.regstr(regstri));   % smaller better/slower
    reg(regstri).opt.ll3th     = 0.0001 + 0.10*(1-job.extopts.regstr(regstri));   % smaller better/slower 
  elseif job.extopts.regstr(regstri)==4 
  % Default Shooting  
    reg(regstri).opt.rres        = tempres;           % registration resolution depending on template resolution 
    reg(regstri).opt.stepsize    = 0;                 % stepsize of reduction 
    reg(regstri).opt.nits        = 24;                % Dartel default iteration number
    reg(regstri).opt.resfac      = ones(1,5);         % reduction factor 
    reg(regstri).opt.ll1th       = 0;                 % smaller better/slower
    reg(regstri).opt.ll3th       = 0;                 % smaller better/slower 
  elseif job.extopts.regstr(regstri)==5 % this may not work because you finally need another interpolation! 
  % based on vox  
    tempres                      = job.extopts.vox(voxi);   % template resolution 
    reg(regstri).opt.rres        = job.extopts.vox(voxi);   % registration resolution depending on template resolution 
    reg(regstri).opt.stepsize    = (tempres*2 - tempres)/4;  % stepsize of reduction 
    reg(regstri).opt.resfac      = (tempres*2 : -reg(regstri).opt.stepsize : job.extopts.vox(voxi)) / job.extopts.vox(voxi); % reduction factor 
  elseif job.extopts.regstr(regstri)==2 
  % Optimized Shooting - manual limit 
    reg(regstri).opt.stepsize    = (lowres - tempres)/4;  % stepsize of reduction 
    reg(regstri).opt.resfac      = (lowres : -reg(regstri).opt.stepsize : tempres) / tempres; % reduction factor 
  elseif job.extopts.regstr(regstri)==3
  % Optimized Shooting - dynamic limit (depending on template resolution)  
    reg(regstri).opt.stepsize    = (tempres/2 - tempres)/4;  % stepsize of reduction 
    reg(regstri).opt.resfac      = (tempres/2 : -reg(regstri).opt.stepsize : tempres) / tempres; % reduction factor 
  else
    % futher test cases
    highres = 1.0; % 0.5
    switch job.extopts.regstr(regstri)
    % -----------------------------------------------------------------
    % absolute fixed resolutions and reduction
    % -----------------------------------------------------------------
    % This allows identical smooth iterations for all levels and some
    % kind of frequency/deformation limit. 
    % The resolution levels are 1.0:0.5:3.0 mm, but can be changed to 
    % to 0.5:0.5:2.5 mm.
    % default = 12 | 22, expert = 11:13 | 21:23
    % -----------------------------------------------------------------
      case {10,11,12,13,14,15,16,17} 
        % independent of the TR and therefore can include interpolation
        reg(regstri).opt.rres     = 1.0 + 0.5 * (job.extopts.regstr(regstri) - 11); 
        highres                   = min(highres, reg(regstri).opt.rres);
        reg(regstri).opt.stepsize = 0.5; 
        reg(regstri).opt.ll1th    = 0.005 * reg(regstri).opt.rres;                 % smaller better/slower
        reg(regstri).opt.ll3th    = 0.010 * reg(regstri).opt.rres;                 % smaller better/slower 
        reg(regstri).opt.resfac   = max( highres + 4*reg(regstri).opt.stepsize:-reg(regstri).opt.stepsize : highres , reg(regstri).opt.rres ) / reg(regstri).opt.rres;   
      case {21,22,23,24,25}
        % dependent on the TR without interpolation interpolation
        reg(regstri).opt.rres     = tempres; %max(tempres,1.0 + 0.5 * (job.extopts.regstr(regstri) - 21)); 
        reg(regstri).opt.stepsize = 0.5; 
        reg(regstri).opt.ll1th    = 0.005 * reg(regstri).opt.rres;                 % smaller better/slower
        reg(regstri).opt.ll3th    = 0.010 * reg(regstri).opt.rres;                 % smaller better/slower 
        reg(regstri).opt.resfac   = ones(1,5); %max( highres + 4*reg(regstri).opt.stepsize : -reg(regstri).opt.stepsize : highres , reg(regstri).opt.rres ) / reg(regstri).opt.rres;   
    % -----------------------------------------------------------------
    % There are further cases, but they all have some drawbacks:
    % * Using a fixed stepsize with different final resolutions run 
    %   into the problematic low resolution (<3 mm).
    % * Using a additive or multiplicative template depending reduction
    %   will be bad for the GUI and it is better to have some fixed
    %   levels.
    % -----------------------------------------------------------------
      otherwise
         error('cat_main_registration:incorrectparameter','Incorrect value of "regres".\n');
    end

    % not required from a theoretic point of view but ...
    reg(regstri).opt.ll1th = reg(regstri).opt.ll1th * tempres/tpmres; 
    reg(regstri).opt.ll3th = reg(regstri).opt.ll3th * tempres/tpmres; 
  end


  %% manual setting of shooting parameter  
  if 0
    job.extopts.vox(voxi)          = 1.5;                   %#ok<UNRCH> % output resolution ...
    reg(regstri).opt.stepsize      = max(eps,0.5);          % stepsize of reduction
    reg(regstri).opt.rres          = 1.0;                   % registration resolution expert parameter ...
    reg(regstri).opt.vxreg         = 1.5;                   % regularisation parameter that original depend on the template resolution
    reg(regstri).opt.nits          = 5;                     % registration iteration (shooting default = 24)
    reg(regstri).opt.ll1th         = 0.01;                  % smaller better/slower
    reg(regstri).opt.ll3th         = 0.04;                  % smaller better/slower 
    res.do_dartel                  = 2;                     % method, 1 - Dartel, 2 - Shooting
    reg.fast                       = 10;                    % inner iterations limit
  end


  %% set default dartel/shooting templates in debug mode 
 % ... error for bad template numbers ... use defaults? if res.tmp2
  if 0 %debug
    % only in case of the default templates
    job.extopts.templates = templates; 
    if job.extopts.expertgui==2 && ...
       (res.do_dartel==1 && job.extopts.regstr(regstri)>0) || ...
       (res.do_dartel==2 && job.extopts.regstr(regstri)==0)
      if res.do_dartel==2 && job.extopts.regstr(regstri)==0
        cat_io_cprintf('warn','Switch to default Dartel Template.\n');
        job.extopts.templates      = cat_vol_findfiles(fullfile(spm('dir'),'toolbox','cat12','templates_volumes'),'Template_*_IXI555_MNI152.nii'); 
        job.extopts.templates(end) = []; 
        reg(regstri).opt.rres = job.extopts.vox(voxi); 
      elseif res.do_dartel==1 && job.extopts.regstr(regstri)>0
        cat_io_cprintf('warn','Switch to default Shooting Template.\n');
        job.extopts.templates = cat_vol_findfiles(fullfile(spm('dir'),'toolbox','cat12','templates_volumes'),'Template_*_IXI555_MNI152_GS.nii',struct('depth',1)); 
      end
    end
    res.tmp2 = cell(1,numel(job.extopts.templates)); 
  end
  run2 = struct(); 
  for j=1:numel(res.tmp2)
    for i=1:n1, run2(i).tpm = sprintf('%s,%d',job.extopts.templates{j},i);end
    res.tmp2{j} = spm_vol(char(cat(1,run2(:).tpm)));
  end
  
  % preparte output directory
  [temppp,tempff] = spm_fileparts(job.extopts.templates{1}); clear temppp;  
  tempres2 = reg(regstri).opt.resfac * reg(regstri).opt.rres; 
  newres   = job.extopts.vox(voxi); 
  regres   = reg(regstri).opt.rres;
  if job.extopts.regstr(regstri)==0
    reg.testfolder = sprintf('Dartel_%s_rr%0.1f_default',tempff,newres);
  elseif job.extopts.regstr(regstri)==4
    reg.testfolder = sprintf('Shooting_%s_rr%0.1f_or%0.1f_default',tempff,regres,newres);
  else 
    if reg(regstri).opt.stepsize>10^-3 
      reg.testfolder = sprintf('Shooting_%s_tr%0.1f_rr%0.1f-%0.1f_or%0.1f_regstr%0.1f%s',tempff,tempres,tempres2([1,5]),newres,job.extopts.regstr(regstri));
    else
      reg.testfolder = sprintf('Shooting_%s_tr%0.1f_rr%0.1f_or%0.1f_regstr%0.1f%s',tempff,tempres,regres,newres,job.extopts.regstr(regstri));
    end
  end
end
%=======================================================================
function report(job,reg,regstri)
% report function   

  % display registration power profil 
  if job.extopts.expertgui==2
    mdisplay = [2 2 2 2 2]; 
  elseif job.extopts.expertgui==1
    mdisplay = [0 2 0 2 0]; 
  else
    mdisplay = [0 1 0 1 0]; 
  end  
  if job.extopts.regstr(regstri)==0, dartelfac = 1.5; else, dartelfac = 1.0; end
  if mdisplay(1)
    fprintf('Registration power: \n'); 
    fprintf('%30s','Jacobian determinant: '); 
    QMC   = cat_io_colormaps('marks+',17);
    reg(regstri).reldtc = reg(regstri).dtc / max(reg(regstri).dtc); 
    color = @(QMC,m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
    if mdisplay(1)>1
      for dti=1:numel(reg(regstri).dtc)-1
        cat_io_cprintf( color(QMC,( 1 - reg(regstri).reldtc(dti)/dartelfac/0.25 ) *6),sprintf('%0.3f ',reg(regstri).reldtc(dti))); 
      end
      fprintf('| ');
    end
    cat_io_cprintf( color(QMC,(reg(regstri).dt - 0.05)/dartelfac/0.25 * 6), sprintf(' %0.6f ',reg(regstri).dt));
    fprintf('\n'); 
  end

  if mdisplay(2)
    fprintf('%30s','Jacobian determinant (RMS): '); 
    QMC   = cat_io_colormaps('marks+',17);
    reg(regstri).relrmsdtc = reg(regstri).rmsdtc; %/max(reg(regstri).rmsdtc); 
    color = @(QMC,m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
    if mdisplay(2)>1
      for dti=1:numel(reg(regstri).relrmsdtc)-1
        cat_io_cprintf( color(QMC,reg(regstri).relrmsdtc(dti)/dartelfac/0.5*6),sprintf('%0.3f ',reg(regstri).relrmsdtc(dti))); 
      end
      fprintf('| ');
    end
    cat_io_cprintf( color(QMC,(reg(regstri).rmsdt)/dartelfac/0.5 * 6), sprintf(' %0.6f ',reg(regstri).rmsdt));
    fprintf('\n'); 
  end

  if mdisplay(3)
    fprintf('%30s','Jacobian determinant'' (RMS): '); 
    QMC   = cat_io_colormaps('marks+',17);
    color = @(QMC,m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
    if mdisplay(3)>1
      for dti=1:numel(reg(regstri).rmsgdt)-1
        cat_io_cprintf( color(QMC,reg(regstri).rmsgdt(dti)/dartelfac/0.5*6),sprintf('%0.3f ',reg(regstri).rmsgdt(dti))); 
      end
      fprintf('| ');
    end
    cat_io_cprintf( color(QMC,(reg(regstri).rmsgdt(end))/dartelfac/0.5 * 6), sprintf(' %0.6f ',reg(regstri).rmsgdt(end)));
    fprintf('\n'); 
  end

  if mdisplay(4)
    % this work very well
    fprintf('%30s','Template Matching: '); 
    QMC   = cat_io_colormaps('marks+',17);
    color = @(QMC,m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
    if mdisplay(4)>1
      for dti=1:size(reg(regstri).ll,1)-1 
        cat_io_cprintf( color(QMC, (reg(regstri).ll(dti,1) - 0.05) / 0.15 * 6),sprintf('%0.3f ',reg(regstri).ll(dti,1))); 
      end
      fprintf('| ');
    end
    cat_io_cprintf( color(QMC,(reg(regstri).ll(end,1) - 0.05)/0.15 * 6), sprintf(' %0.6f ',reg(regstri).ll(end,1)));
    fprintf('\n'); 
  end

  reg(regstri).cbr  = diff(reg(regstri).relrmsdtc(1:end-1)) ./ -diff(reg(regstri).ll(1:end-1,1)');
  reg(regstri).scbr = sum(diff(reg(regstri).relrmsdtc(1:end-1)) ./ -diff(reg(regstri).ll(1:end-1,1)')); 
  reg(regstri).mcbr = mean(diff(reg(regstri).relrmsdtc(1:end-1)) ./ -diff(reg(regstri).ll(1:end-1,1)')); 
  if mdisplay(5)
    % this work very well
    fprintf('%30s','Cost Benefit Ratio (CBR): ');
    QMC   = cat_io_colormaps('marks+',17);
    color = @(QMC,m) QMC(max(1,min(size(QMC,1),round(((m-1)*3)+1))),:);
    if mdisplay(5)>1
      for dti=1:size(reg(regstri).cbr,2)
        cat_io_cprintf( color(QMC, (reg(regstri).cbr(dti))/2),sprintf('%0.3f ',reg(regstri).cbr(dti))); 
      end
      fprintf('| ');
    end
    cat_io_cprintf( color(QMC,( reg(regstri).mcbr)/2), sprintf(' %0.6f ', reg(regstri).mcbr));
    fprintf('\n'); 
  end
end
%=======================================================================
function write_nii(Ycls,job,trans,testfolder,reg)
% Fast export of output data to test and run different settings in the same
% preprocessing.
  
  trans.warped.verb = 1; % test output
  
  stime = cat_io_cmd(sprintf('Write Output with %0.2f mm',job.extopts.vox(1)));
  
  if job.extopts.subfolders, mrifolder = 'mri'; else, mrifolder = ''; end
  if ~exist('testfolder','var'); testfolder = sprintf('backup_or%3.1f',job.extopts.vox(1)); end
  [pth,nam] = spm_fileparts(trans.native.Vo.fname); 
  
  % registration information 
  if exist('reg','var')
    cat_io_xml(fullfile(pth,mrifolder,testfolder,['reg_', nam, '.xml']),reg)
  end
  
  % tissue ouptut
  fn = {'GM','WM','CSF'};
  for clsi = 1 %:numel(Ycls)
    if ~isempty(Ycls{clsi})
      cat_io_writenii(trans.native.Vo,single(Ycls{clsi})/255,fullfile(mrifolder,testfolder),sprintf('p%d',clsi),...
        sprintf('%s tissue map',fn{clsi}),'uint8' ,[0,1/255],[0 0 0 3],trans);
      cat_io_writenii(trans.native.Vo,single(Ycls{clsi})/255,fullfile(mrifolder,testfolder),sprintf('p%d',clsi),...
        sprintf('%s tissue map',fn{clsi}),'uint8' ,[0,1/255],[0 1 0 0],trans);
      cat_io_writenii(trans.native.Vo,single(Ycls{clsi})/255,fullfile(mrifolder,testfolder),sprintf('p%d',clsi),...
        sprintf('%s tissue map',fn{clsi}),'uint16',[0,1/255],[0 0 3 0],trans); 
    end
  end

  if job.output.jacobian.warped
  %% write jacobian determinant
    if isfield(trans.jc,'dt2'); % % shooting
      dt2 = trans.jc.dt2; 
      dx = 10; % smaller values are more accurate, but large look better; 
      [D,I] = cat_vbdist(single(~(isnan(dt2) | dt2<0 | dt2>100) )); D=min(1,D/min(dx,max(D(:)))); 
      dt2 = dt2(I); dt2 = dt2 .* ((1-D) + D); 
      dt2(isnan(dt2))=1; 
    else % dartel
      [y0, dt2] = spm_dartel_integrate(reshape(trans.jc.u,[trans.warped.odim(1:3) 1 3]),[1 0], 6); clear y0;  
    end
    
    % create nifti
    N         = nifti;
    N.dat     = file_array(fullfile(pth,mrifolder,testfolder,['wj_', nam, '.nii']),trans.warped.odim(1:3),...
                [spm_type('float32') spm_platform('bigend')],0,10/256^2,0);
    N.mat     = trans.warped.M1;
    N.mat0    = trans.warped.M1;
    N.descrip = ['Jacobian' trans.native.Vo.descrip];
    create(N);
    N.dat(:,:,:) = dt2;
  end
  
  if job.output.warps(1)
  %% deformations y - dartel > subject
    Yy2       = spm_diffeo('invdef',trans.warped.y,trans.warped.odim,eye(4),trans.warped.M0);
    N         = nifti;
    N.dat     = file_array(fullfile(pth,mrifolder,testfolder,['y_', nam, '.nii']),[trans.warped.odim(1:3),1,3],'float32',0,1,0);
    N.mat     = trans.warped.M1;
    N.mat0    = trans.warped.M1;
    N.descrip = 'Deformation';
    create(N);
    N.dat(:,:,:,:,:) = reshape(Yy2,[trans.warped.odim,1,3]);
    clear Yy2; 
  end
  
  
  if job.output.warps(2)
  %% deformation iy - subject > dartel
    if any(trans.native.Vo.dim~=trans.native.Vi.dim)
      vx_voli  = sqrt(sum(trans.native.Vi.mat(1:3,1:3).^2));  
      vx_volo  = sqrt(sum(trans.native.Vo.mat(1:3,1:3).^2));
      eyev     = eye(4); eyev([1 6 11]) = eyev([1 6 11]) .* vx_volo./vx_voli; 
      Yy2      = zeros([trans.native.Vo.dim 1 3],'single');                        
      for k1 = 1:3
        for i = 1:trans.native.Vo.dim(3)
          Yy2(:,:,i,:,k1) = trans.warped.M1(k1,4) + trans.warped.M1(k1,k1) * ...
            single(spm_slice_vol(trans.warped.y(:,:,:,k1),eyev*spm_matrix([0 0 i]), ...
            trans.native.Vo.dim(1:2),[1,NaN])); % adapt for res
        end
      end
    else 
      yn     = numel(trans.warped.y); 
      p      = ones([4,yn/3],'single'); 
      p(1,:) = trans.warped.y(1:yn/3);
      p(2,:) = trans.warped.y(yn/3+1:yn/3*2);
      p(3,:) = trans.warped.y(yn/3*2+1:yn);
      p      = trans.warped.M1(1:3,:) * p;

      Yy2                = zeros([trans.native.Vo.dim(1:3),1,3],'single'); 
      Yy2(1:yn/3)        = p(1,:);
      Yy2(yn/3+1:yn/3*2) = p(2,:);
      Yy2(yn/3*2+1:yn)   = p(3,:);
    end
    clear p; 

    % f2 = spm_diffeo('resize', f1, dim)
    % write new output
    Ndef      = nifti;
    Ndef.dat  = file_array(fullfile(pth,mrifolder,testfolder,['iy_', nam, '.nii']),[trans.native.Vo.dim,1,3],...
                [spm_type('float32') spm_platform('bigend')],0,1,0);
    Ndef.mat  = res.image0(1).mat;
    Ndef.mat0 = res.image0(1).mat;
    Ndef.descrip = 'Inverse Deformation';
    create(Ndef);
    Ndef.dat(:,:,:,:,:) = Yy2;
    clear Yy2;
  end
  
  cat_io_cmd('','',''); cat_io_cmd('','','',job.extopts.verb,stime);

end
%=======================================================================
function x = rgrid(d)
  x = zeros([d(1:3) 3],'single');
  [x1,x2] = ndgrid(single(1:d(1)),single(1:d(2)));
  for i=1:d(3)
      x(:,:,i,1) = x1;
      x(:,:,i,2) = x2;
      x(:,:,i,3) = single(i);
  end
end
%=======================================================================
function y1 = affind(y0,M)
  y1 = zeros(size(y0),'single');
  for d=1:3
      y1(:,:,:,d) = y0(:,:,:,1)*M(d,1);
      y1(:,:,:,d) = y1(:,:,:,d) + y0(:,:,:,2)*M(d,2);
      y1(:,:,:,d) = y1(:,:,:,d) + y0(:,:,:,3)*M(d,3) + M(d,4);
  end
end
%=======================================================================
function [x1,y1,z1] = defs2(sol,z,M,prm,x0,y0,z0)
  iM = inv(M);
  z01 = z0(z)*ones(size(x0));

  x1a  = iM(1,1)*x0 + iM(1,2)*y0 + iM(1,3)*z01 + iM(1,4);
  y1a  = iM(2,1)*x0 + iM(2,2)*y0 + iM(2,3)*z01 + iM(2,4);
  z1a  = iM(3,1)*x0 + iM(3,2)*y0 + iM(3,3)*z01 + iM(3,4);

  x1 = spm_bsplins(sol{1},x1a,y1a,z1a,prm);
  y1 = spm_bsplins(sol{2},x1a,y1a,z1a,prm);
  z1 = spm_bsplins(sol{3},x1a,y1a,z1a,prm);
end
%=======================================================================
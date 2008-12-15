function varargout = cg_vbm8_run(job,arg)
% Segment a bunch of images
% FORMAT cg_vbm8_run(job)
% job.channel(n).vols{m}
% job.channel(n).biasreg
% job.channel(n).biasfwhm
% job.channel(n).write
% job.tissue(k).tpm
% job.tissue(k).ngaus
% job.tissue(k).native
% job.tissue(k).warped
% job.warp.affreg
% job.warp.reg
% job.warp.samp
% job.warp.write
% job.warp.bb
% job.warp.vox
%
% See the user interface for a description of the fields.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% based on John Ashburners version of
% spm_preproc8_run.m 2281 2008-10-01 12:52:50Z john $
%
% Christian Gaser
% $Id$

rev = '$Rev$';

if nargin==1,
    varargout{:} = run_job(job);
elseif strcmpi(arg,'check'),
    varargout{:} = check_job(job);
elseif strcmpi(arg,'vfiles'),
    varargout{:} = vfiles_job(job);
elseif strcmpi(arg,'vout'),
    varargout{:} = vout_job(job);
else
    error('Unknown argument ("%s").', arg);
end
return
%_______________________________________________________________________

%_______________________________________________________________________
function vout = run_job(job)

vout   = vout_job(job);
tpm    = strvcat(cat(1,job.tissue(:).tpm));
tpm    = spm_load_priors8(tpm);

% check whether estimation & write should be done
estwrite = isfield(job.warp,'reg');

nit = 1;

for iter=1:nit,
    if nit>1,
        % Sufficient statistics for possible generation of group-specific
        % template data.
        SS = zeros([size(tpm.dat{1}),numel(tpm.dat)],'single');
    end
    for subj=1:numel(job.channel(1).vols),
        if estwrite % estimate and write segmentations
            images = '';
            for n=1:numel(job.channel),
                images = strvcat(images,job.channel(n).vols{subj});
            end
            obj.image    = spm_vol(images);
            spm_check_orientations(obj.image);

            obj.fudge    = 5;
            obj.biasreg  = cat(1,job.channel(:).biasreg);
            obj.biasfwhm = cat(1,job.channel(:).biasfwhm);
            obj.tpm      = tpm;
            obj.lkp      = [];
            if all(isfinite(cat(1,job.tissue.ngaus))),
                for k=1:numel(job.tissue),
                    obj.lkp = [obj.lkp ones(1,job.tissue(k).ngaus)*k];
                end;
            end
            obj.reg      = job.warp.reg;
            obj.samp     = job.warp.samp;

            if iter==1,
                % Initial affine registration.
                Affine  = eye(4);
                if ~isempty(job.warp.affreg),
                    if job.warp.affmethod == 0
                        Affine  = spm_maff8(obj.image(1),job.warp.samp,obj.fudge*8,tpm,Affine,job.warp.affreg); % Close to rigid
                        Affine  = spm_maff8(obj.image(1),job.warp.samp,obj.fudge,  tpm,Affine,job.warp.affreg)
                    else
                        VG = spm_vol(fullfile(spm('Dir'),'templates','T1.nii'));
                        VF = spm_vol(obj.image(1));
                    
                        % smooth source with 8mm
                        VF1 = spm_smoothto8bit(VF,8);

                        % Rescale images so that globals are better conditioned
                        VF1.pinfo(1:2,:) = VF1.pinfo(1:2,:)/spm_global(VF1);
                        VG.pinfo(1:2,:)  = VG.pinfo(1:2,:)/spm_global(VG);

                        fprintf('Coarse Affine Registration..\n');
                        aflags    = struct('sep',8, 'regtype',job.warp.affreg,...
                            'WG',[],'WF',[],'globnorm',0);
                        aflags.sep = max(aflags.sep,max(sqrt(sum(VG(1).mat(1:3,1:3).^2))));
                        aflags.sep = max(aflags.sep,max(sqrt(sum(VF(1).mat(1:3,1:3).^2))));

                        M = eye(4);
                        spm_chi2_plot('Init','Affine Registration','Mean squared difference','Iteration');
                        [M,scal]  = spm_affreg(VG, VF1, aflags, M);
 
                        fprintf('Fine Affine Registration..\n');
                        aflags.WG  = spm_vol(fullfile(spm('Dir'),'apriori','brainmask.nii'));
                        aflags.sep = aflags.sep/2;
                        [Affine,scal]   = spm_affreg(VG, VF1, aflags, M,scal);
                    end
                end;
                obj.Affine = Affine;
            else
                % Load results from previous iteration for use with next round of
                % iterations, with the new group-specific tissue probability map.
                [pth,nam] = fileparts(job.channel(1).vols{subj});
                res       = load(fullfile(pth,[nam '_seg8.mat']));
                obj.Affine = res.Affine;
                obj.Twarp  = res.Twarp;
                obj.Tbias  = res.Tbias;
                if ~isempty(obj.lkp),
                    obj.mg     = res.mg;
                    obj.mn     = res.mn;
                    obj.vr     = res.vr;
                end
            end

            res = spm_preproc8(obj);

            try
                [pth,nam] = fileparts(job.channel(1).vols{subj});
                savefields(fullfile(pth,[nam '_seg8.mat']),res);
            catch
            end

        else % only write segmentations
            [pth,nam] = fileparts(job.channel(1).vols{subj});
            seg8_name = fullfile(pth,[nam '_seg8.mat']);
            if exist(seg8_name)
                res = load(seg8_name);
                % use path of mat-file in case that image was moved
                for i=1:numel(res.image)
					[image_pth,image_nam,image_ext]=fileparts(res.image(i).fname);
                	res.image(i).fname = fullfile(pth,[image_nam,image_ext]);
                end
            else
                error(['Can''t load file ' seg8_name]);  
                return
            end
        end
          
        if iter==nit,
            % Final iteration, so write out the required data.
            tmp1 = [cat(1,job.tissue(:).native) cat(1,job.tissue(:).warped)];
            tmp2 =  cat(1,job.channel(:).write);
            tmp3 = job.warp.write;
            cg_vbm8_write(res,tmp1, tmp2, tmp3, job.warp);
        else
            % Not the final iteration, so compute sufficient statistics for
            % re-estimating the template data.
            N    = numel(job.channel);
            K    = numel(job.tissue);
            cls  = spm_preproc_write8(res,zeros(K,4),zeros(N,2),[0 0]);
            for k=1:K,
                SS(:,:,:,k) = SS(:,:,:,k) + cls{k};
            end
        end

    end
    if iter<nit && nit>1,
         % Treat the tissue probability maps as Dirichlet priors, and compute the 
         % MAP estimate of group tissue probability map using the sufficient
         % statistics.
         alpha = 1.0;
         for k=1:K,
             SS(:,:,:,k) = SS(:,:,:,k) + spm_bsplinc(tpm.V(k),[0 0 0  0 0 0])*alpha + eps;
         end

         s = sum(SS,4);
         for k=1:K,
             tmp        = SS(:,:,:,k)./s;
             tpm.bg(k)  = mean(mean(tmp(:,:,1)));
             tpm.dat{k} = spm_bsplinc(log(tmp+tpm.tiny),[ones(1,3)*(tpm.deg-1)  0 0 0]);
         end
    end
end
return
%_______________________________________________________________________

%_______________________________________________________________________
function msg = check_job(job)
msg = {};
if numel(job.channel) >1,
    k = numel(job.channel(1).vols);
    for i=2:numel(job.channel),
        if numel(job.channel(i).vols)~=k,
            msg = {['Incompatible number of images in channel ' num2str(i)]};
            break
        end
    end
elseif numel(job.channel)==0,
    msg = {'No data'};
end
return
%_______________________________________________________________________

%_______________________________________________________________________
function savefields(fnam,p)
if length(p)>1, error('Can''t save fields.'); end;
fn = fieldnames(p);
if numel(fn)==0, return; end;
for i=1:length(fn),
    eval([fn{i} '= p.' fn{i} ';']);
end;
if spm_matlab_version_chk('7') >= 0
    save(fnam,'-V6',fn{:});
else
    save(fnam,fn{:});
end;

return;
%_______________________________________________________________________

%_______________________________________________________________________
function vout = vout_job(job)

n     = numel(job.channel(1).vols);
parts = cell(n,4);

channel = struct('biasfield',{},'biascorr',{});
for i=1:numel(job.channel),
    for j=1:n,
        [parts{j,:}] = spm_fileparts(job.channel(i).vols{j});
    end
    if job.channel(i).write(1),
        channel(i).biasfield = cell(n,1);
        for j=1:n
            channel(i).biasfield{j} = fullfile(parts{j,1},['BiasField_',parts{j,2},'.nii']);
        end
    end
    if job.channel(i).write(2),
        channel(i).biascorr = cell(n,1);
        for j=1:n
            channel(i).biascorr{j} = fullfile(parts{j,1},['m',parts{j,2},'.nii']);
        end
    end
end

for j=1:n,
    [parts{j,:}] = spm_fileparts(job.channel(1).vols{j});
end
param = cell(n,1);
for j=1:n
    param{j} = fullfile(parts{j,1},[parts{j,2},'_seg8.mat']);
end

tiss = struct('c',{},'rc',{},'wc',{},'mwc',{});
for i=1:numel(job.tissue),
    if job.tissue(i).native(1),
        tiss(i).c = cell(n,1);
        for j=1:n
            tiss(i).c{j} = fullfile(parts{j,1},['p',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).native(2),
        tiss(i).rc = cell(n,1);
        for j=1:n
            tiss(i).rc{j} = fullfile(parts{j,1},['rp',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).warped(1),
        tiss(i).wc = cell(n,1);
        for j=1:n
            tiss(i).wc{j} = fullfile(parts{j,1},['wp',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).warped(2),
        tiss(i).mwc = cell(n,1);
        for j=1:n
            tiss(i).mwc{j} = fullfile(parts{j,1},['mwp',num2str(i),parts{j,2},'.nii']);
        end
    end
end

if job.warp.write(1),
    invdef = cell(n,1);
    for j=1:n
        invdef{j} = fullfile(parts{j,1},['iy_',parts{j,2},'.nii']);
    end
else
    invdef = {};
end

if job.warp.write(2),
    fordef = cell(n,1);
    for j=1:n
        fordef{j} = fullfile(parts{j,1},['y_',parts{j,2},'.nii']);
    end
else
    fordef = {};
end

vout  = struct('channel',channel,'tiss',tiss,'param',{param},'invdef',{invdef},'fordef',{fordef});
%_______________________________________________________________________

%_______________________________________________________________________
function vf = vfiles_job(job)
vout = vout_job(job);
vf   = vout.param;
if ~isempty(vout.invdef), vf = {vf{:}, vout.invdef{:}}; end
if ~isempty(vout.fordef), vf = {vf{:}, vout.fordef{:}}; end
for i=1:numel(vout.channel),
    if ~isempty(vout.channel(i).biasfield), vf = {vf{:}, vout.channel(i).biasfield{:}}; end
    if ~isempty(vout.channel(i).biascorr),  vf = {vf{:}, vout.channel(i).biascorr{:}};  end
end

for i=1:numel(vout.tiss)
    if ~isempty(vout.tiss(i).c),   vf = {vf{:}, vout.tiss(i).c{:}};   end
    if ~isempty(vout.tiss(i).rc),  vf = {vf{:}, vout.tiss(i).rc{:}};  end
    if ~isempty(vout.tiss(i).wc),  vf = {vf{:}, vout.tiss(i).wc{:}};  end
    if ~isempty(vout.tiss(i).mwc), vf = {vf{:}, vout.tiss(i).mwc{:}}; end
end
vf = reshape(vf,numel(vf),1);
%_______________________________________________________________________

%_______________________________________________________________________

%=======================================================================
function t = transf(B1,B2,B3,T)
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t = zeros(size(B1,1),size(B2,1),size(B3,1));
end;
return;
%=======================================================================

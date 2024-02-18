function cat_tst_histBWPT(S,qa,msave,ff)
% cat_test_BWPT_hist. Print figure of histogram for thickness phantom.

  s = 0; 

  if ~exist('msave','var'), msave = 0; end
  if ~exist('ff','var')
    ff = ''; 
  else
    [~,ff] = spm_fileparts(ff); 
  end
  
  % smoothing for robust peaks
  if s > 0
    M  = spm_mesh_smooth(S.rh);
    T  = single(spm_mesh_smooth(M,double(S.rh.th1),s));
  else
    T  = S.rh.th1;
  end

  
  if msave 
    fh = figure(4949); 
  else
    fh = figure; 
  end
  
  ss = 0.01;
  hx = ss/2:ss:5; 
  hy = hist(T,hx); 
  hy(2:end-1) = hy(2:end-1)*0.5 + hy(1:end-2)*0.25 + hy(3:end)*0.25; %smooth
  gh = histogram(T,hx); 
 
  fh.Position = [1151 1 560 420];   
  gh.EdgeColor = 'none'; gh.FaceColor = [0 0.2 0.4]; gh.FaceAlpha = .9;
  grid on; xlim([0 4]); ylim([1 1.2] .* ylim);
  
  [hy15,hx15] = max(hy .* (hx>1.25 & hx<=1.75)); hx15x = hx(hx15) - 1.5; 
  [hy20,hx20] = max(hy .* (hx>1.75 & hx<=2.25)); hx20x = hx(hx20) - 2.0; 
  [hy25,hx25] = max(hy .* (hx>2.25 & hx<=2.75)); hx25x = hx(hx25) - 2.5; 

  % print title
  justnow = datestr(datetime,'YYYYmmdd-HHMM'); 
  if exist('qa','var')
    title(sprintf('histogram thickness BWPT'), ...
      sprintf('%s [th=%0.0f, RMSE=%0.3f mm, Int: %0.3f, Pos: %0.3f, IS: %0.03f%%]', ...
      justnow, sum([hy15,hy20,hy25]), sum([hx15x,hx20x,hx25x].^2).^0.5, ...
      mean([qa.rh.createCS_final.RMSE_Ym_white, ...
            qa.rh.createCS_final.RMSE_Ym_layer4, ...
            qa.rh.createCS_final.RMSE_Ym_pial]), ...
      mean([qa.rh.createCS_final.RMSE_Ypp_white, ...
            qa.rh.createCS_final.RMSE_Ypp_central, ...
            qa.rh.createCS_final.RMSE_Ypp_pial]), ...
      mean([qa.rh.createCS_final.white_self_interections, ...
            qa.rh.createCS_final.pial_self_interections]) ));
  else
    title(sprintf('histogram thickness BWPT'), ...
      sprintf('%s [th=%d, RMSE=%0.3f mm]', ...
      justnow, sum([hy15,hy20,hy25]), sum([hx15x,hx20x,hx25x].^2).^0.5));
  end

  % add datatips for 3 peaks
  dp = datatip(gh,'dataindex',hx15); dp.Location = "northwest";
  dp = datatip(gh,'dataindex',hx20);
  dp = datatip(gh,'dataindex',hx25);

  if msave
    %%
    sdir = fullfile(spm('dir'),'toolbox','cat12','internal','CSx_pbt_test',justnow);
    mkdir(sdir); 
    copyfile(which('cat_vol_pbtsimple'),sdir);
    copyfile(which('cat_surf_createCS3'),sdir);
    print( fh, fullfile(fileparts(sdir),sprintf('BWPT_result_%s%s.png',ff,justnow)) ,'-dpng')
  end
  if msave > 1
    close( fh )
  end
end
% Parameters
pixel_um   = 0.88;          % µm
pixel_m    = pixel_um * 1e-6;
tau_fs     = 40;            % fs
tau_s      = tau_fs * 1e-15;

% Energies
E1_mJ = 8.3;      E1_J = E1_mJ * 1e-3;       % for input.tiff
E2_mJ = 7.636;    E2_J = E2_mJ * 1e-3;       % for output.tiff

% -------- Case 1: input.tiff --------
[in1, bg] = LoadLaserImages('input.tiff', 'background.tiff');
corr1 = SubtractBackground(in1, bg);
[~, I1_Wcm2] = ComputePeakIntensity(corr1, pixel_m, E1_J, tau_s);

% -------- Case 2: output.tiff (same background) --------
[in2, bg] = LoadLaserImages('output.tiff', 'background.tiff');
corr2 = SubtractBackground(in2, bg);
[~, I2_Wcm2] = ComputePeakIntensity(corr2, pixel_m, E2_J, tau_s);

% ---- Show both intensity maps side-by-side (W/cm^2) ----
figure; tiledlayout(1,2, 'Padding','compact', 'TileSpacing','compact');

nexttile;
imagesc(I1_Wcm2); axis image; colormap hot; cb1=colorbar;
cb1.Label.String = 'Peak Intensity (W/cm^2)';
title(sprintf('input.tiff  (%.3f mJ)', E1_mJ));
xlabel('X pixel'); ylabel('Y pixel');

nexttile;
imagesc(I2_Wcm2); axis image; colormap hot; cb2=colorbar;
cb2.Label.String = 'Peak Intensity (W/cm^2)';
title(sprintf('output.tiff (%.3f mJ)', E2_mJ));
xlabel('X pixel'); ylabel('Y pixel');

% ---- Print peak intensities to Command Window ----
Ipk1 = max(I1_Wcm2(:));
Ipk2 = max(I2_Wcm2(:));
fprintf('Peak intensity (input.tiff, %.3f mJ):  %.3g W/cm^2\n',  E1_mJ, Ipk1);
fprintf('Peak intensity (output.tiff, %.3f mJ): %.3g W/cm^2\n',  E2_mJ, Ipk2);

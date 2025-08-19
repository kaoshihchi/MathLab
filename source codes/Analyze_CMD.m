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

% -------- Dipole phase from intensity --------
global q_e m_e hbar
q_e  = 1.602176634e-19;      % elementary charge (C)
m_e  = 9.1093837015e-31;     % electron mass (kg)
hbar = 1.054571817e-34;      % reduced Planck constant (J*s)
mu_0 = 4*pi*1e-7;            % vacuum permeability (H/m)
c    = 2.99792458e8;         % speed of light (m/s)

q          = 25;             % harmonic order
lambda_nm  = 800;            % driving wavelength (nm)
omega_d    = 2*pi*c/(lambda_nm*1e-9);
I_p_eV     = 15.7596;        % ionization potential of Ar (eV)
I_p        = I_p_eV * q_e;   % ionization potential (J)

I1_Wm2 = I1_Wcm2 * 1e4;      % convert W/cm^2 -> W/m^2
I2_Wm2 = I2_Wcm2 * 1e4;
E1 = sqrt(2*mu_0*c*I1_Wm2);  % electric field amplitude (V/m)
E2 = sqrt(2*mu_0*c*I2_Wm2);

Phi1_l = zeros(size(E1));    % long-trajectory phase (rad)
Phi1_s = zeros(size(E1));    % short-trajectory phase (rad)
for idx = 1:numel(E1)
    temp = HHG_DipolePhase(q, E1(idx), I_p, omega_d);
    Phi1_l(idx) = temp(1);
    Phi1_s(idx) = temp(2);
end

Phi2_l = zeros(size(E2));
Phi2_s = zeros(size(E2));
for idx = 1:numel(E2)
    temp = HHG_DipolePhase(q, E2(idx), I_p, omega_d);
    Phi2_l(idx) = temp(1);
    Phi2_s(idx) = temp(2);
end

% ---- Show dipole phase maps ----
figure; tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile;
imagesc(Phi1_l); axis image; colorbar;
title('Long \phi (input)');

nexttile;
imagesc(Phi1_s); axis image; colorbar;
title('Short \phi (input)');

nexttile;
imagesc(Phi2_l); axis image; colorbar;
title('Long \phi (output)');

nexttile;
imagesc(Phi2_s); axis image; colorbar;
title('Short \phi (output)');

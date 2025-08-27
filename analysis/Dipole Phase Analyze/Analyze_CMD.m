close all; clear;

%% -------- Constants (SI) --------
global q_e m_e hbar
q_e  = 1.602176634e-19;      % C
m_e  = 9.1093837015e-31;     % kg
hbar = 1.054571817e-34;      % J*s
mu_0 = 4*pi*1e-7;            % H/m
c    = 2.99792458e8;         % m/s

%% -------- User parameters --------
pixel_um = 0.88;                       % µm
pixel_m  = pixel_um * 1e-6;
tau_fs   = 40;                         % fs
tau_s    = tau_fs * 1e-15;

E1_mJ = 8.3;      E1_J = E1_mJ * 1e-3; % input.tiff
E2_mJ = 7.6;    E2_J = E2_mJ * 1e-3; % output.tiff

q         = 47;                        % harmonic order
lambda_nm = 808;                       % nm
omega_d   = 2*pi*c/(lambda_nm*1e-9);
I_p_eV    = 15.7596;                   % eV (Ar)
I_p       = I_p_eV * q_e;              % J

%% -------- Load images once --------
bg  = im2double(imread('background.tiff'));
in1 = im2double(imread('input.tiff'));
in2 = im2double(imread('output_0.3psi.tiff'));

% Background subtraction (clip at 0)
corr1 = max(in1 - bg, 0);
corr2 = max(in2 - bg, 0);

%% -------- Peak intensity maps (W/cm^2) --------
% If you already have these helpers, keep them. Otherwise:
% [~, I_Wcm2] = ComputePeakIntensity(corr, pixel_m, E_J, tau_s);
[~, I1_Wcm2] = ComputePeakIntensity(corr1, pixel_m, E1_J, tau_s);
[~, I2_Wcm2] = ComputePeakIntensity(corr2, pixel_m, E2_J, tau_s);

% Display intensity maps
figure; tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
nexttile; imagesc(I1_Wcm2); axis image; colormap hot; cb=colorbar;
cb.Label.String='Peak Intensity (W/cm^2)';
title(sprintf('input.tiff  (%.3f mJ)', E1_mJ)); xlabel('X pixel'); ylabel('Y pixel');

nexttile; imagesc(I2_Wcm2); axis image; colormap hot; cb=colorbar;
cb.Label.String='Peak Intensity (W/cm^2)';
title(sprintf('output.tiff (%.3f mJ)', E2_mJ)); xlabel('X pixel'); ylabel('Y pixel');

Ipk1 = max(I1_Wcm2(:)); Ipk2 = max(I2_Wcm2(:));
fprintf('Peak intensity (input.tiff,  %.3f mJ): %.3g W/cm^2\n', E1_mJ, Ipk1);
fprintf('Peak intensity (output.tiff, %.3f mJ): %.3g W/cm^2\n', E2_mJ, Ipk2);

%% -------- Use I_low to skip pixels before dipole phase --------
[I_low_Wm2, ~, Up_eV] = Find_I_low(q, omega_d);
fprintf('I_low = %.3e W/m^2  (U_p = %.2f eV)\n', I_low_Wm2, Up_eV);

% Convert intensity maps once to W/m^2
I1_Wm2 = I1_Wcm2 * 1e4;
I2_Wm2 = I2_Wcm2 * 1e4;

% Masks of pixels to compute
mask1 = I1_Wm2 >= I_low_Wm2;
mask2 = I2_Wm2 >= I_low_Wm2;

% ---- E-field from intensity (your requested formula) ----
% I = E^2 / (2*mu0*c)  =>  E = sqrt(2*mu0*c*I)
E1 = zeros(size(I1_Wm2)); E1(mask1) = sqrt(2*mu_0*c*I1_Wm2(mask1));
E2 = zeros(size(I2_Wm2)); E2(mask2) = sqrt(2*mu_0*c*I2_Wm2(mask2));

%% -------- Dipole phase (skip low pixels) --------
Phi1_l = NaN(size(I1_Wm2)); Phi1_s = NaN(size(I1_Wm2));
Phi2_l = NaN(size(I2_Wm2)); Phi2_s = NaN(size(I2_Wm2));

% Compute only where needed
Phi1_l(mask1) = compute_phase_map(E1(mask1), q, I_p, omega_d, 'long');
Phi1_s(mask1) = compute_phase_map(E1(mask1), q, I_p, omega_d, 'short');
Phi2_l(mask2) = compute_phase_map(E2(mask2), q, I_p, omega_d, 'long');
Phi2_s(mask2) = compute_phase_map(E2(mask2), q, I_p, omega_d, 'short');

% ---- Show dipole phase maps ----
figure; tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile;
imagesc(Phi1_l); axis image; 
cb = colorbar; cb.Label.String = 'Phase (rad)';
title('\phi_{dipole\_long} (input)');

nexttile;
imagesc(Phi1_s); axis image; 
cb = colorbar; cb.Label.String = 'Phase (rad)';
title('\phi_{dipole\_short} (input)');

nexttile;
imagesc(Phi2_l); axis image; 
cb = colorbar; cb.Label.String = 'Phase (rad)';
title('\phi_{dipole\_long} (output)');

nexttile;
imagesc(Phi2_s); axis image; 
cb = colorbar; cb.Label.String = 'Phase (rad)';
title('\phi_{dipole\_short} (output)');


%% -------- Local helper (vector-friendly) --------
function phi_out = compute_phase_map(E_flat, q, I_p, omega_d, whichBranch)
    % E_flat: vector of E-field values at masked pixels (V/m)
    % whichBranch: 'long' or 'short'
    phi_out = NaN(size(E_flat));
    % Optionally: use parfor if you have PCT
    for i = 1:numel(E_flat)
        tmp = HHG_DipolePhase(q, E_flat(i), I_p, omega_d);
        if strcmp(whichBranch,'long'),  phi_out(i) = tmp(1);
        else,                          phi_out(i) = tmp(2);
        end
    end
end

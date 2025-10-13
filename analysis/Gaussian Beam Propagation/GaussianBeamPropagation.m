% GaussianBeamPropagation_R2023a.m
% Simulate Gaussian beam propagation and peak intensity along z.
% Assumptions:
% - z = 0 at focus
% - Given "focal spot FWHM 30 um" is the *intensity FWHM diameter* at focus
% - Temporal and spatial profiles are Gaussian (separable)

clear; clc;

%% ----- Parameters -----
lambda   = 810e-9;         % wavelength [m]
fwhm_d0  = 30e-6;          % focal FWHM diameter at z=0 [m]
E_pulse  = 2e-3;           % pulse energy [J]
tau_fwhm = 35e-15;         % pulse duration (intensity FWHM) [s]

z_min = 0e-3;              % [m]
z_max = 5e-3;              % [m]
Nz    = 2001;              % number of z samples

%% ----- Conversions and core formulas -----
% Intensity Gaussian: I(r) = I0 * exp(-2 r^2 / w^2)
% FWHM diameter relation at a given z: d_FWHM(z) = sqrt(2*ln(2)) * 2 * w(z) / 2?  -> d = w * sqrt(2 ln2)
% Derivation: I(r_FWHM/2) = I0/2 -> r_half = w*sqrt(ln2/2) -> diameter d = 2*r_half = w*sqrt(2*ln2)
% => w0 from FWHM diameter d0:
w0 = fwhm_d0 / sqrt(2*log(2));     % 1/e^2 radius at focus [m]

% Rayleigh length
zR = pi * w0^2 / lambda;           % [m]

% z grid
z = linspace(z_min, z_max, Nz);

% Beam radius vs z
wz = w0 .* sqrt(1 + (z./zR).^2);   % [m]

% FWHM diameter vs z
fwhm_d = wz .* sqrt(2*log(2));     % [m]

% Peak intensity vs z (spatiotemporal Gaussian):
% Total energy E = I0_peak * (∫_A e^{-2r^2/w^2} dA) * (∫_t e^{-4 ln2 t^2 / tau_fwhm^2} dt)
% Spatial integral = (pi * w^2) / 2
% Temporal integral = tau_fwhm * sqrt(pi) / (2*sqrt(ln 2))
% => I0_peak(z) = E / [ (pi*w(z)^2/2) * (tau_fwhm * sqrt(pi)/(2*sqrt(ln2))) ]
I0_peak = E_pulse ./ ( (pi.*wz.^2./2) .* (tau_fwhm.*sqrt(pi)./(2*sqrt(log(2)))) ); % [W/m^2]

% Also provide in W/cm^2 for convenience
I0_peak_Wcm2 = I0_peak / 1e4;

%% ----- Plot 1: FWHM vs z -----
fig1 = figure('Color','w','Position',[100 100 720 480]);
plot(z*1e3, fwhm_d*1e6, 'LineWidth', 2);
grid on;
xlabel('z (mm)','Interpreter','none');
ylabel('FWHM diameter (\mum)','Interpreter','tex');
title(sprintf('Gaussian Beam FWHM vs z (\\lambda = %.0f nm, FWHM@focus = %.0f \\mum)', ...
      lambda*1e9, fwhm_d0*1e6));
% Save as SVG
print(fig1, 'gaussian_fwhm_vs_z.svg', '-dsvg');

%% ----- Plot 2: Peak intensity vs z -----
fig2 = figure('Color','w','Position',[120 120 720 480]);
plot(z*1e3, I0_peak_Wcm2, 'LineWidth', 2);
grid on;
xlabel('z (mm)','Interpreter','none');
ylabel('Peak Intensity (W/cm^2)','Interpreter','none');
title(sprintf('Peak Intensity vs z (E=%.1f mJ, \\tau_{FWHM}=%.0f fs)', ...
      E_pulse*1e3, tau_fwhm*1e15));
% Optionally show on log scale (uncomment if desired)
% set(gca,'YScale','log');

% Save as SVG
print(fig2, 'gaussian_peak_intensity_vs_z.svg', '-dsvg');

%% ----- Console summary -----
fprintf('--- Summary ---\n');
fprintf('lambda       = %.0f nm\n', lambda*1e9);
fprintf('FWHM@focus   = %.1f um (diameter)\n', fwhm_d0*1e6);
fprintf('w0 (1/e^2 r) = %.2f um\n', w0*1e6);
fprintf('zR           = %.3f mm\n', zR*1e3);
fprintf('I0@focus     = %.3e W/cm^2\n', I0_peak_Wcm2(1));
fprintf('Saved: gaussian_fwhm_vs_z.svg, gaussian_peak_intensity_vs_z.svg\n');

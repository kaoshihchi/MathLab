%% Make Super-Gaussian Spectrum TXT (run this)
% Generates a two-column text file:
%   col 1: wavelength [nm]
%   col 2: intensity  [arb. units]
% Compatible with: data = load(filename);  % numeric-only file (no header)

clear; clc; close all;

% ---------- User Parameters ----------
out_file   = '131900_average spectrum.txt'; % output path
center_nm  = 810;     % center wavelength [nm]
fwhm_nm    = 60;    % FWHM [nm]
order_n    = 4;       % 1 = Gaussian, >1 = super-Gaussian (e.g., 2, 4, 8)
amplitude  = 1.0;     % peak amplitude
offset     = 0.0;     % baseline offset
noise_std  = 0.00;    % Gaussian noise std (0 = no noise)
seed       = 42;      % RNG seed for reproducibility

lam_min    = 700;     % lower bound [nm]
lam_max    = 900;     % upper bound [nm]
step_nm    = 1.0;     % sampling step [nm]

% ---------- Generate ----------
if ~isempty(seed)
    rng(seed);  %#ok<RAND> % reproducible noise
end

gen = SpectrumGenerator(center_nm, fwhm_nm, order_n, amplitude, offset, noise_std);

lambda_nm = (lam_min:step_nm:lam_max).';
intensity = gen.evaluate(lambda_nm);
intensity = max(intensity, 0);  % keep nonnegative

% ---------- Save (space-separated, no header) ----------
data = [lambda_nm, intensity];
writematrix(data, out_file, 'Delimiter', ' ');

fprintf('Wrote %d rows to: %s\n', size(data,1), out_file);
fprintf('Center = %.2f nm,  FWHM = %.2f nm,  Order n = %g\n', center_nm, fwhm_nm, order_n);

% ---------- Plot ----------
figure('Color','w');
plot(lambda_nm, intensity, 'LineWidth', 1.8);
grid on;
xlabel('Wavelength (nm)');
ylabel('Intensity (arb. units)');
title(sprintf('Super-Gaussian Spectrum  (\\lambda_0=%.1f nm, FWHM=%.1f nm, n=%g)', center_nm, fwhm_nm, order_n));
xlim([lam_min lam_max]);

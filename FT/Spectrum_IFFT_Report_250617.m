% Initialization
close all;
clear;
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, '..', 'theory'));
warning('off', 'all');

% Set default font properties for all figures
set(groot, 'defaultAxesFontName', 'Helvetica');
set(groot, 'defaultTextFontName', 'Helvetica');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextFontSize', 12);

% Constants
c = 299792458; % Speed of light in vacuum (m/s)

% Shared optimization options
optimOptions = struct('MaxIter', 100, 'MaxFunEvals', 500, 'MaxTime', 120);

% Load Data
filename = fullfile(scriptDir, '131900_average spectrum.txt');

% Use fileparts to separate out the parts of the file path
[~, name, ~] = fileparts(filename);
data = load(filename);
if false
    wavelength = data(:, 1) .* 1e9;
else
    wavelength = data(:, 1);
end

intensity = data(:, 2);

figure
plot(wavelength, intensity)
xlabel('Wavelength (nm)');
ylabel('Intensity (counts)');
title('Original Input Spectrum');
grid on;
%xlim([700, 900]);
% Calculate background noise
bk_indices = wavelength >= 700 & wavelength <= 750;
bk_mean = mean(intensity(bk_indices));
bk_std = std(intensity(bk_indices)); 
clear bk_indices; 

% Make intensity below few times bk_std as zero to filter DC signal
% avoiding intensity < 0
bk_indices = intensity <= (bk_mean + 2*bk_std); 
intensity_corr = intensity - bk_mean; 
intensity_corr(bk_indices) = 0; 
clear bk_indices; 

figure
plot(wavelength, intensity, wavelength, intensity_corr);
legend('Input spectrum', 'BK corrected spectrum');
xlabel('Wavelength (nm)');
ylabel('Normalized Intensity');
title('Background corrected spectrum');

% Convert to frequency domain
frequency = c ./ (wavelength * 1e-9);
% normalized_intensity = (intensity-bk_mean) ./ max(intensity-bk_mean);
normalized_intensity = (intensity_corr) ./ max(intensity_corr);

% Electric field amplitude vs. Frequency
electric_field_amplitude = abs(sqrt(normalized_intensity));

% High-resolution frequency domain
N_points = 1000 * length(frequency);
f_min = min(frequency);
f_max = max(frequency);
even_frequencies = linspace(-f_max, f_max, N_points);
extended_electric_field_amplitude = zeros(size(even_frequencies));

% Create a logical mask for the range where interpolation is valid
interp_mask = even_frequencies >= f_min & even_frequencies <= f_max;


extended_electric_field_amplitude(even_frequencies >= f_min & even_frequencies <= f_max) = interp1(frequency, electric_field_amplitude, even_frequencies(even_frequencies >= f_min & even_frequencies <= f_max), 'linear'); 
% even_electric_field_amplitude = interp1(frequency, electric_field_amplitude, even_frequencies, 'linear');

% Time-domain analysis
E_t = ifft(extended_electric_field_amplitude);

% Plot: Electric Field Amplitude vs. Frequency
figure('Name', 'E(f) - Electric Field Amplitude vs Frequency', 'NumberTitle', 'off');
plot(even_frequencies, abs(extended_electric_field_amplitude), 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('|E(f)| (arb. units)');
title('Electric Field Amplitude vs. Frequency');
grid on;
xlim([0, 1.5e15]);

delta_f = even_frequencies(2) - even_frequencies(1);
f_s = delta_f * N_points; % Sampling rate
dt = 1 / f_s; 
T_total = 1 / delta_f;
t = linspace(0, N_points, N_points) .* dt;
I_t = abs(E_t(1:end)).^2;
t = t(1:end);

% Normalize E_t before plotting if needed
E_t_norm = E_t / max(abs(E_t));

% Plot: Electric Field Amplitude vs. Time
figure('Name', 'E(t) - Electric Field Amplitude vs Time', 'NumberTitle', 'off');
plot(t * 1e15, abs(E_t_norm), 'LineWidth', 2); % time in fs
xlabel('Time (fs)');
ylabel('|E(t)| (arb. units)');
title('Electric Field Amplitude vs. Time');
grid on;
xlim([min(t * 1e15), max(t * 1e15)]);


% Shift intensity in time domain
[max_val, max_idx] = max(I_t);
shift_amount = floor(length(I_t) / 2) - max_idx;
I_t_shifted = circshift(I_t, [0, shift_amount]);
I_t_shifted = I_t_shifted ./ max(I_t_shifted); 
t_shifted = linspace(-T_total / 2, T_total / 2, length(I_t));

% There is large noise at t_shifted around 0 fs, so we need to do some
% filter. 
[I_t_shifted_filtered, t_shifted_filtered] = filterData(I_t_shifted, t_shifted, 0e-15);
[gaussianFit, gaussianGof] = createFit(t_shifted_filtered, I_t_shifted_filtered, optimOptions);

figure('Name', 'IFFT Check', 'NumberTitle', 'off');
subplot(2, 1, 1)
plot(even_frequencies, extended_electric_field_amplitude, '-', 'LineWidth',2);
xlabel('Frequency (Hz)');
ylabel('Intensity (arb. units)');
title('Extended Intensity in Frequency Domain I(f)');
%xlim([-100E-15 100E-15]);
grid on;

subplot(2, 1, 2)
plot(t_shifted_filtered .* 1e15, I_t_shifted_filtered, 'o', 'LineWidth',2);
xlabel('Time (fs)');
ylabel('Intensity (arb. units)');
title('IFFT Intensity in Time Domain I(t)');
xlim([-100 100]);
grid on;
% Frequency domain analysis of shifted signal
E_f = fft(E_t);
T = max(t) - min(t); 
N_FFT = length(t); 
df = N_FFT / T; 
frequencies = (df / N_FFT) .* (-1/2 * N_FFT:N_FFT/2 - 1);
E_f_magnitude = abs(E_f);

% Convert to wavelength domain and plot I(lambda)
positive_frequencies = frequencies(frequencies > 0);
wavelengths_lambda = c ./ positive_frequencies;
E_lambda = E_f(frequencies > 0);
I_lambda = abs(E_lambda).^2;

figure
subplot(3, 1, 1)
plot(t, abs(E_t))
xlabel('Time (s)');
ylabel('Intensity (arb. units)');
title('Intensity in Time Domain I(t)');
grid on; 
subplot(3, 1, 2)
plot(even_frequencies, extended_electric_field_amplitude);
xlim([0 1.5e15])
xlabel('Frequency (Hz)');
ylabel('Intensity (arb. units)');
title('Input Intensity in Frequency Domain I(f)');
grid on;
subplot(3, 1, 3)
plot(frequencies, E_f_magnitude)
xlim([0 1.5e15])
xlabel('Frequency (Hz)');
ylabel('Intensity (arb. units)');
title('Output Intensity in Frequency Domain I(f), After IFFT & FFT');
grid on;
figure
subplot(2, 1, 1)
plot(wavelength, intensity_corr)
xlabel('Wavelength (nm)');
ylabel('Intensity (arb. units)');
title('Background Corrected Intensity in Wavelength I(\lambda)');
grid on;
xlim([750 850]);
subplot(2, 1, 2)
plot(wavelengths_lambda .* 1e9, I_lambda)
xlabel('Wavelength (nm)');
ylabel('Intensity (arb. units)');
title('Output Intensity in Wavelength I(\lambda), After IFFT & FFT');
grid on;
xlim([750 850]);

% Define the Super-Gaussian fit type
ftype = fittype(@(A, B, C, n, D, x) superGaussian(x, A, B, C, n, D));

% Set fitting options
opts = fitoptions('Method', 'NonlinearLeastSquares');
opts.StartPoint = [1, 810, 60, 2, 0];  % Initial guess
opts.MaxIter = optimOptions.MaxIter;                   % Limit number of iterations
opts.MaxFunEvals = optimOptions.MaxFunEvals;
opts.MaxTime = optimOptions.MaxTime;               % Optional: limit function evaluations
opts.Display = 'Off';                 % Suppress fitting display

% Super-Gaussian fitting on normalized intensity
[IinFitResult, ~] = fit(wavelength, normalized_intensity, ftype, opts);

% Super-Gaussian fitting after IFFT & FFT process
[Iout1FitResult, ~] = fit(wavelengths_lambda' * 1e9, I_lambda', ftype, opts);


% Combined figure for Intensity Spectrum vs Wavelength
figure('Name', 'Combined Intensity Spectrum vs Wavelength', 'NumberTitle', 'off');
subplot(2, 1, 1);
plot(wavelength, normalized_intensity, 'r-', wavelength, superGaussian(wavelength, IinFitResult.A, IinFitResult.B, IinFitResult.C, IinFitResult.n, IinFitResult.D), 'b--', 'LineWidth', 2);
xlabel('Wavelength (nm)');
ylabel('Normalized Intensity');
title('Input Intensity Spectrum vs Wavelength');
grid on;
xlim([700, 900]);
subplot(2, 1, 2);
plot(wavelengths_lambda * 1e9, I_lambda, 'r-', wavelengths_lambda * 1e9, superGaussian(wavelengths_lambda * 1e9, Iout1FitResult.A, Iout1FitResult.B, Iout1FitResult.C, Iout1FitResult.n, Iout1FitResult.D), 'b--', 'LineWidth', 2);
xlabel('Wavelength (nm)');
ylabel('Intensity (arb. units)');
title('Reconstructed Intensity Spectrum (After IFFT & FFT) vs Wavelength I(\lambda)');
grid on;
xlim([700, 900]);
% Save spectrum data
saveSpectrumData(wavelength, normalized_intensity, [char(name), '_BKCorrectedSpectrumData.txt']);
disp(IinFitResult.B)
disp(IinFitResult.C * 1.825)
disp(Iout1FitResult.B)
disp(Iout1FitResult.C * 1.825)
% Define the fitting range
t_lim = 0.1e-12; % Limit parameter for fitting

% Filter the data based on the fitting range
fittingRangeIndices = t_shifted_filtered >= -t_lim & t_shifted_filtered <= t_lim;
t_shifted_fit = t_shifted_filtered(fittingRangeIndices) * 1e15;     % (fs)
I_t_shifted_fit = I_t_shifted_filtered(fittingRangeIndices);
I_t_shifted_fit = I_t_shifted_fit ./ max(I_t_shifted_fit); 

[xData, yData] = prepareCurveData( t_shifted_fit, I_t_shifted_fit );

% Set up fittype and options.
ft = fittype( 'A*exp(-2*(abs((x)/C).^(2)))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 10];
opts.Robust = 'Bisquare';
opts.StartPoint = [1.5 50];
opts.Upper = [2 100];
opts.MaxIter = optimOptions.MaxIter;
opts.MaxFunEvals = optimOptions.MaxFunEvals;
opts.MaxTime = optimOptions.MaxTime;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'SSA' );
h = plot( fitresult, xData, yData );
legend( h, 'I_t_shifted_fit vs. t_shifted_fit', 'SSA', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 't_shifted_fit', 'Interpreter', 'none' );
ylabel( 'I_t_shifted_fit', 'Interpreter', 'none' );
grid on
% Combined figure for Shifted Intensity in Time Domain
figure('Name', 'Shifted Intensity in Time Domain', 'NumberTitle', 'off');
subplot(2, 1, 1);
plot(t_shifted_fit .* 1e-15, I_t_shifted_fit, 'ro', t_shifted, GaussianAmp(t_shifted, fitresult.A, 0, fitresult.C * 1e-15, 0), 'b--', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Intensity (arb. units)');
title('IFFT Intensity in Time Domain I(t)');
xlim([-100E-15 100E-15]);
grid on;
subplot(2, 1, 2);
semilogy(t_shifted_filtered, I_t_shifted_filtered);
%plot(t_shifted_filtered, I_t_shifted_filtered);
xlabel('Time (s)');
ylabel('Intensity (arb. units, log scale)');
title('IFFT Intensity in Time Domain I(t) on Logarithmic Scale');
grid on;
xlim([0 1e-12]);
ylim([1e-5 1]);
saveIFFTPulse(t_shifted_fit', I_t_shifted_fit', [char(name), '_IFFT_TimeProfile.txt']); 
disp(num2str(fitresult.C * sqrt(2*log(2)), '%.2f'))
% Save data, and down size data size to 1/1000
% Assuming t_shifted_filtered and I_t_shifted_filtered are your original arrays
N = 10; % Downsampling factor

% Downsampling by selecting every N-th element
t_shifted_filtered_downsampled = t_shifted_filtered(1:N:end)';
I_t_shifted_filtered_downsampled = I_t_shifted_filtered(1:N:end)';

[t_selected, I_selected] = selectDataInRange(t_shifted_filtered_downsampled, I_t_shifted_filtered_downsampled, 0, 5e-12); 

saveIFFTData(t_selected .*1e12, I_selected, [char(name), '_IFFTProfile.txt']);

figure
semilogy(t_selected, I_selected, '-o');
%plot(t_shifted_filtered, I_t_shifted_filtered);
xlabel('Time (s)');
ylabel('Intensity (arb. units, log scale)');
title('IFFT Intensity in Time Domain I(t) on Logarithmic Scale');
grid on;
xlim([0 5e-12]);
ylim([1e-6 1]);
function y = superGaussian(x, A, B, C, n, D)
    y = A * exp(-(abs((x - B) / C).^(2*n))) + D;
end

function y = GaussianAmp(x, A, B, C, D)
    y = A * exp(-2 * (abs((x - B) / C).^(2))) + D;
end

function [I_t_shifted_filtered, t_shifted_filtered] = filterData(I_t_shifted, t_shifted, t_bound)
    % Convert fs to seconds for comparison
    time_lower_bound = -t_bound; % -5 fs in seconds
    time_upper_bound = t_bound;  % 5 fs in seconds

    % Create a mask for values outside the range -5 fs to 5 fs
    mask = t_shifted < time_lower_bound | t_shifted > time_upper_bound;

    % Apply the mask to filter the data
    I_t_shifted_filtered = I_t_shifted(mask);
    t_shifted_filtered = t_shifted(mask);
end

function saveSpectrumData(array1, array2, filename)
    % Combine the arrays into a matrix for easy writing
    combinedData = [array1, array2];
    
    % Define main headers
    mainHeaders = {'Wavelength', 'Intensity'};
    
    % Define units (to be placed in the second row)
    units = {'nm', 'arb. units'};
    
    % Check if the filename is provided without '.txt' extension and append if necessary
    if ~endsWith(filename, '.txt')
        filename = [filename, '.txt'];
    end
    
    % Open the file for writing
    fid = fopen(filename, 'w');
    
    % Write main headers and units to the file
    fprintf(fid, '%s,%s\n', mainHeaders{:});
    fprintf(fid, '%s,%s\n', units{:});
    
    % Close the file
    fclose(fid);
    
    % Append data
    writematrix(combinedData, filename, 'WriteMode', 'append');
end

function saveIFFTData(array1, array2, filename)
    % saveIFFTData(wavelength, normalized_intensity, [char(name), '_Y-axis Line Profile.txt']);
    % Combine the arrays into a matrix for easy writing
    combinedData = [array1, array2];
    
    % Define main headers
    mainHeaders = {'Time', 'Intensity'};
    
    % Define units (to be placed in the second row)
    units = {'ps', 'arb. units'};
    
    % Check if the filename is provided without '.txt' extension and append if necessary
    if ~endsWith(filename, '.txt')
        filename = [filename, '.txt'];
    end
    
    % Open the file for writing
    fid = fopen(filename, 'w');
    
    % Write main headers and units to the file
    fprintf(fid, '%s,%s\n', mainHeaders{:});
    fprintf(fid, '%s,%s\n', units{:});
    
    % Close the file
    fclose(fid);
    
    % Append data
    writematrix(combinedData, filename, 'WriteMode', 'append');
end

function [t_selected, I_selected] = selectDataInRange(t_array, I_array, lowerBound, upperBound)
    % Find indices where t_array is within the specified range
    inRangeIndices = t_array >= lowerBound & t_array <= upperBound;

    % Select data from both arrays using the found indices
    t_selected = t_array(inRangeIndices);
    I_selected = I_array(inRangeIndices);
end

function saveIFFTPulse(array1, array2, filename)
    % saveIFFTData(wavelength, normalized_intensity, [char(name), '_Y-axis Line Profile.txt']);
    % Combine the arrays into a matrix for easy writing
    combinedData = [array1, array2];
    
    % Define main headers
    mainHeaders = {'Time', 'Intensity'};
    
    % Define units (to be placed in the second row)
    units = {'fs', 'arb. units'};
    
    % Check if the filename is provided without '.txt' extension and append if necessary
    if ~endsWith(filename, '.txt')
        filename = [filename, '.txt'];
    end
    
    % Open the file for writing
    fid = fopen(filename, 'w');
    
    % Write main headers and units to the file
    fprintf(fid, '%s,%s\n', mainHeaders{:});
    fprintf(fid, '%s,%s\n', units{:});
    
    % Close the file
    fclose(fid);
    
    % Append data
    writematrix(combinedData, filename, 'WriteMode', 'append');
end

clear; clc;
set(0,'defaultfigurevisible','off');

% Add theory functions to path
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, '..', 'theory'));

% Global constants
global eV fs c q_e m_e mu_0 epsilon_0 h hbar

eV   = 1.602e-19;
fs   = 1e-15;
c    = 2.99792458e8;
q_e  = 1.602e-19;
m_e  = 9.109e-31;
mu_0 = 4*pi*1e-7;
epsilon_0 = 8.854e-12;
h    = 6.626e-34;
hbar = h/(2*pi);

% Laser parameters
lambda  = 800e-9;                 % wavelength (m)
omega_d = 2*pi*c/lambda;         % angular frequency (rad/s)
tau     = 40/1.2*fs;                  % pulse duration (1/e) (s)

% Time vector for field calculation
t = (-200:0.2:200)*fs;

% Ionization potentials (converted to Joules)
Eion_Ar = [15.76, 27.63, 40.74, 59.81, 75.02]*eV;
Eion_Ne = [21.56, 40.96, 63.45, 97.11, 126.21]*eV;
Eion_He = [24.59, 54.42, 1000, 1000, 1000]*eV;  % large values suppress >2+

% Intensity range in W/cm^2
I_Wcm2 = linspace(5e14, 5e15, 50);

Z_Ar = zeros(size(I_Wcm2));
Z_Ne = zeros(size(I_Wcm2));
Z_He = zeros(size(I_Wcm2));

for k = 1:length(I_Wcm2)
    I = I_Wcm2(k)*1e4;          % convert to W/m^2
    E_peak = sqrt(2*mu_0*c*I);  % peak electric field (V/m)
    E_t = E_peak * exp(-(t/tau).^2) .* exp(1i*omega_d*t);

    res = TunnelingIonizationRate_Linear(E_t,omega_d,t,Eion_Ar(1),Eion_Ar(2),Eion_Ar(3),Eion_Ar(4),Eion_Ar(5),0);
    Z_Ar(k) = res(1);

    res = TunnelingIonizationRate_Linear(E_t,omega_d,t,Eion_Ne(1),Eion_Ne(2),Eion_Ne(3),Eion_Ne(4),Eion_Ne(5),0);
    Z_Ne(k) = res(1);

    res = TunnelingIonizationRate_Linear(E_t,omega_d,t,Eion_He(1),Eion_He(2),Eion_He(3),Eion_He(4),Eion_He(5),0);
    Z_He(k) = res(1);
end

% Plot ionization ratios
figure;
plot(I_Wcm2, Z_Ar, 'r', I_Wcm2, Z_Ne, 'g', I_Wcm2, Z_He, 'b', 'LineWidth', 1.5);
xlabel('Laser intensity (W/cm^2)');
ylabel('Average ionization state');
legend({'Ar','Ne','He'}, 'Location', 'northwest');
title('Ionization ratio vs laser intensity');

% Save plot to .fig file
outputFile = fullfile(scriptDir, 'ionization_ratio.fig');
savefig(gcf, outputFile);
disp(['Ionization ratio figure saved to ', outputFile]);


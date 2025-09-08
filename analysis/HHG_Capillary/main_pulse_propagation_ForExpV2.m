clear all
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'theory'));
% MainPulsePropagation_HHG.m
%  Part I:
%    Calculate the main beam loss propagated in a capillary waveguide.
%    1. Consider ionization loss, ATI heating loss, inverse Bremsstrahlung
%       heating loss, Thomson scattering loss, and capillary guiding loss.
%    2. Consider the plasma dispersion and the capillary dispersion.
%    3. Use a transverse-average electric field E(z,t) in the capillary.
%    4. SI units.
%
%  Part II:
%    Calculate the HHG driven by the main pulse.
%
%    Syntax:
%       MainPulsePropagation_HHG('ParameterFilename',FigureSwitch)
%
%       FigureSwitch: figure switch of ionization
%                     = 1: figure ON, = 0: figure OFF.
%
%    Returned value:
%       colume 1: position 'z' (m)
%       colume 2: laser pulse energy 'LaserEnergy' at z (J)
%       colume 3: average ionization state 'Z_ion' (arb. units)
%       colume 4: electron density 'N_e' (1/m^3)
%       colume 5: plasma group-velocity dispersion 'GVD' at z (sec^2/m)
%       colume 6: accumulated group-delay dispersion (GDD) 'D' at z (sec^2)
%       colume 7: laser pulse duration 'tau' at z (sec)
%       colume 8: laser peak intensity 'I_peak' at z (W/m^2)
%       colume 9: laser peak electric field 'E_peak' at z (V/m)
%       colume 10: electron temperature 'kT' (k_B times T) (J)
%       colume 11: total ionization energy 'Energy_Ionization' at (z,z+dz) (J)
%       colume 12: total ATI heating energy 'Energy_ATI' at (z,z+dz) (J)
%       colume 13: total IB heating energy 'Energy_IB' at (z,z+dz) (J)
%       colume 14: total Thomson scattering energy 'Energy_TS' at (z,z+dz) (J)
%       colume 15: capillary attenuated energy 'Energy_capillary' at (z,z+dz) (J)
%
%    Exported file:
%       'MainPulsePropagation_HHG_results_SI_units.txt': resuls with SI units
%       'MainPulsePropagation_HHG_results_usual_units.txt': resuls with usual units
%
%    Example:
%       a = MainPulsePropagation_HHG('MainPulsePropagation_HHG_Parameters.txt',0);
%
%     Author: Hsu-hsin Chu (2021/7/13)
% 
%    Fixed: wavenumber mismatch of capillary (20230224)

%function result = MainPulsePropagation_HHG(ParameterFilename,FigureSwitch)

%%
close all

q_HHG = linspace(0,0,11);
for qq = 1 : 1
% Unit conversion
   global eV
   global cm
   global mm
   global um
   global nm
   global fs
   global ps
   global mJ
   eV   = 1.602 * 10^(-19);          % 1 (eV) = 1.6*10^(-19) (J)
   cm   = 10^(-2);                   % 1 (cm) = 10^(-2) (m)
   mm   = 10^(-3);                   % 1 (mm) = 10^(-3) (m)
   um   = 10^(-6);                   % 1 (um) = 10^(-6) (m)
   nm   = 10^(-9);                   % 1 (nm) = 10^(-9) (m)
   fs   = 10^(-15);                  % 1 (fs) = 10^(-15)(sec)
   ps   = 10^(-12);                  % 1 (ps) = 10^(-12)(sec)
   mJ   = 10^(-3);                   % 1 (mJ) = 10^(-3) (J)

% Physical constants
   global c
   global q_e
   global m_e
   global mu_0
   global epsilon_0
   global h
   global hbar
   c    = 2.99792458 * 10^8;             % speed of light (m/sec)
   q_e  = 1.602 * 10^(-19);          % electron charge (Coul)
   m_e  = 9.101 * 10^(-31);          % electron mass (kg)
   mu_0 = 4 * pi * 10^(-7);          % permeability of vacuum (Nt/Amp^2)
   epsilon_0 = 8.854*10^-12;         % permittivity of vacuum (F/m)
   h    = 6.626 * 10^(-34);          % Planck's constant (J-sec)
   hbar = h/(2*pi);
    % Load simulation parameters from CSV file
    params = readtable(fullfile(fileparts(mfilename('fullpath')), 'Parameters.csv'));
    if iscell(params.OnOff)
        active = params(strcmpi(params.OnOff,'TRUE'), :);
    else
        active = params(params.OnOff == 1, :);
    end
    if isempty(active)
        error('No active parameter set found in Parameters.csv');
    end
    lambda      = active.lambda(1)      * nm;
    tau_0       = active.tau_0(1)       * fs; % FWHM
    PulseEnergy = active.PulseEnergy(1) * mJ;
    R_capillary = active.R_capillary(1) * um;
    L_capillary = active.L_capillary(1) * mm;
    n_capillary = active.n_capillary(1);
    a_capillary = active.a_capillary(1);
    n_gas       = active.n_gas(1)       / cm^3;
    dz          = active.dz(1)          * mm;
    time_initial = active.time_initial(1) * fs;
    time_final   = active.time_final(1)   * fs;
    DeltaTime    = active.DeltaTime(1)    * fs;
    Gas          = active.Gas{1};
    q            = active.q(1);
    dz2_token    = active.dz2{1};
    Ne_exp_avg   = active.Ne_exp_avg(1) * 1E6;    % Experiment (cm^-3)
    t_shift      = active.t_shift(1) * fs;             % (fs)

    % Ionization potentials for the selected gas
    switch Gas
        case 'Ar'
            E_ion_0 = 15.7596 * eV;
            E_ion_1 = 27.6297 * eV;
            E_ion_2 = 40.74   * eV;
            E_ion_3 = 59.81   * eV;
            E_ion_4 = 75.02   * eV;
        case 'He'
           E_ion_0    = 24.587 * eV;          % ionization potential (J) ( 0  -> 1+ )
           E_ion_1    = 54.418 * eV;          % ionization potential (J) ( 1+ -> 2+ )
           E_ion_2    = 1000   * eV;          % ionization potential (J) ( 2+ -> 3+ )
           E_ion_3    = 1200   * eV;          % ionization potential (J) ( 3+ -> 4+ )
           E_ion_4    = 2000   * eV;          % ionization potential (J) ( 4+ -> 5+ )
        case 'Ne'
           E_ion_0    = 21.5645 * eV;          % ionization potential (J) ( 0  -> 1+ )
           E_ion_1    = 40.9630 * eV;          % ionization potential (J) ( 1+ -> 2+ )
           E_ion_2    = 63.45   * eV;          % ionization potential (J) ( 2+ -> 3+ )
           E_ion_3    = 97.12   * eV;          % ionization potential (J) ( 3+ -> 4+ )
           E_ion_4    = 126.21  * eV;          % ionization potential (J) ( 4+ -> 5+ )
        otherwise
            error('Unsupported gas type: %s', Gas);
    end
    I_p = eval(active.I_p{1});
    if ischar(dz2_token) || isstring(dz2_token)
        dz2 = eval(dz2_token);
    else
        dz2 = dz2_token * mm;
    end

    FigureSwitch = 1;

% Laser parameters
   omega_d = 2 * pi * c / lambda;    % laser angular frequency (rad/sec)
   k_0     = 2 * pi / lambda;        % laser wavenumber in vacuum (1/m)
   I_peak_ini = 2 * PulseEnergy / (pi^1.5 * (0.64*R_capillary)^2 * tau_0); % 2P_0 / pi*r^2
                                           % laser peak intensity (W/m^2)
   E_peak_ini = sqrt(2*mu_0*c*I_peak_ini); % laser peak electric field (V/m)
   
% Set the coordinate and variables
   z = [0:dz:L_capillary];    % Generate longitudinal coordinate z (unit: m)
   N = size(z,2);             % number of data points
   Z_ion       = zeros(1,N);  % average ionization state
   N_e         = zeros(1,N);  % electron density N_e(z) (1/m^3)
   Energy_Ionization  = zeros(1,N);   % total ionization energy at (z,z+dz) (J)
   Energy_ATI  = zeros(1,N);  % total ATI energy at (z,z+dz) (J)
   LaserEnergy = zeros(1,N);  % laser pulse energy (J)
   tau         = zeros(1,N);  % laser pulse duration at z (sec)
   I_peak      = zeros(1,N);  % laser peak intensity at z (W/m^2)
   E_peak      = zeros(1,N);  % laser peak electric field at z (V/m)
   kT          = zeros(1,N);  % electron temperature (k_B times T) (J)
   a_IB        = zeros(1,N);  % inverse Bremsstrahlung absorption coefficient (1/m)
   Energy_IB   = zeros(1,N);  % total IB absorbed energy at (z,z+dz) (J)
   a_TS        = zeros(1,N);  % Thomson scattering coefficient (1/m)
   Energy_TS   = zeros(1,N);  % total Thomson scattered energy at (z,z+dz) (J)
   Energy_capillary = zeros(1,N);   % capillary attenuated energy at (z,z+dz) (J)
   GVD         = zeros(1,N);  % group-velocity dispersion at z (sec^2/m)
   D           = zeros(1,N);  % accumulated group-delay dispersion at z (sec^2)
   n_plasma_d  = zeros(1,N);  % plasma refractive index of the driving pulse at z
   C           = zeros(1,N);  % accumulated group delay at z (sec)
   k_d_plasma  = zeros(1,N);  % wavenumber of the driving pulse in the plasma (1/m)
   k_d_total   = zeros(1,N);  % total wavenumber of the driving pulse (1/m)
   phi_d_prop  = zeros(1,N);  % accumulated phase shift of the driving pulse due to propagation (rad)
   v_d         = zeros(1,N);  % phase velocity of the driving pulse (m/sec)
   HHG_gating  = zeros(5,300); % 1 is traced phase, 2 is HHG for long, 3 is HHG for short
% Calculate the initial condition of the main pulse propagation

   LaserEnergy(1) = PulseEnergy;
   tau(1)    = tau_0;
   I_peak(1) = I_peak_ini;
   E_peak(1) = E_peak_ini;

   % Generate initial waveform E(t) at z=0. (complex form, unit: V/m)
   time = [time_initial:DeltaTime:time_final];   % time sequence (sec)
   E_t = ElectricField_d(time,E_peak(1),omega_d,tau_0,0,0,0);
   
   % Calculate the above-threshold ionization (ATI) heating
   ATI_result = TunnelingIonizationRate_Linear...
     (E_t,omega_d,time,E_ion_0,E_ion_1,E_ion_2,E_ion_3,E_ion_4,FigureSwitch);
   if FigureSwitch
      sgtitle('ionization at z = 0');
   end
              
   % average ionizatin state (relative electron density) (arb. units)
   Z_ion(1) = ATI_result(1);
 
   % electron density (1/m^3) at position z = 0.
   N_e(1) = n_gas * Z_ion(1);
   
   % ionization energy per electron (unit: J) at position z = 0.
   Energy_Ionization_e = ATI_result(2);
   
   % ATI energy per electron (unit: J) at position z = 0.
   Energy_ATI_e = ATI_result(3);
   
   % total Ionization energy in a capillary section z = 0~dz (J)
   Energy_Ionization(1) = Energy_Ionization_e * N_e(1) * pi * R_capillary^2 * dz;

   % total ATI energy in a capillary section z = 0~dz (J)
   Energy_ATI(1) = Energy_ATI_e * N_e(1) * pi * R_capillary^2 * dz;
   
   % electron temperature (k_B times T) (J)
   kT(1) = Energy_ATI_e * (2/3);
   
   % inverse Bremsstrahlung absorption coefficient (1/m)
   a_IB(1) = InverseBremsstrahlungCoefficient(N_e(1),Z_ion(1),kT(1),omega_d);
   
   % total IB absorbed energy in a capillary section z = 0~dz (J)
   Energy_IB(1) = LaserEnergy(1) * a_IB(1) * dz;
   
   % Thomson scattering coefficient (1/m)
   a_TS(1) = ThomsonScatteringCoefficient(N_e(1));

   % total Thomson scattered energy in a capillary section z = 0~dz (J)
   Energy_TS(1) = LaserEnergy(1) * a_TS(1) * dz;

   % capillary attenuation coefficient (1/m)
   
   % capillary attenuated energy in a capillary section z = 0~dz (J)
   Energy_capillary(1) = LaserEnergy(1) * a_capillary * dz;
   
   % plasma group-velocity dispersion (GVD) at z (sec^2/m)
   GVD(1) = GroupVelocityDispersion(omega_d,N_e(1));
   
   % waveguide group-velocity dispersion (sec^2/m)
   u_11 = 2.405;
   GVD_capillary = -(u_11^2 * c)/(R_capillary^2 * omega_d^3);
   
   % accumulated group-delay dispersion (GDD) at z (sec^2)
   D(1) = (GVD(1) + GVD_capillary) * dz;

   % plasma refractive index of the driving pulse at z
   n_plasma_d(1) = sqrt(1-(q_e^2 * N_e(1))/(epsilon_0 * m_e * omega_d^2));
   
   % accumulated group delay at z (sec)
   C(1) = dz/(c*n_plasma_d(1)) + (u_11^2*c)/(2*R_capillary^2*omega_d^2)*dz;
   
   % wavenumber of the driving pulse in plasma (1/m)
   k_d_plasma(1) = k_0 * n_plasma_d(1);
   
   % wavenumber variation of the driving pulse in the capillary (1/m)
   k_d_capillary = -(u_11^2*lambda)/(4*pi*R_capillary^2) + 0*k_d_plasma;
   
   % total wavenumber of the driving pulse (1/m)
   k_d_total(1) = k_d_plasma(1) + k_d_capillary(1);
   
   % accumulated phase shift of the driving pulse due to propagation (rad)
   phi_d_prop(1) = k_d_total(1)*dz;
   
   % phase velocity of the driving pulse (m/sec)
   v_d(1)        = omega_d/k_d_total(1);

 
% Calculate the main pulse propagation step-by-step
   for j = 1:N-1
      disp([num2str(j),'/',num2str(N-1)]);
      LaserEnergy(j+1) = LaserEnergy(j) - Energy_Ionization(j) - ...
                         Energy_ATI(j) - Energy_IB(j) - Energy_TS(j) - ...
                         Energy_capillary(j);
      tau(j+1)    = sqrt(tau(j)^2 + D(j)^2/tau(j)^2);
      I_peak(j+1) = 2 * LaserEnergy(j+1)/(pi^1.5 * (0.64*R_capillary)^2 * tau(j+1));
      E_peak(j+1) = sqrt(2*mu_0*c*I_peak(j+1));      % peak electric field (V/m)
      E_t = ElectricField_d(time,E_peak(j+1),omega_d,tau_0,D(j),0,0);
      ATI_result  = TunnelingIonizationRate_Linear...
         (E_t,omega_d,time,E_ion_0,E_ion_1,E_ion_2,E_ion_3,E_ion_4,and(FigureSwitch,j==N-1));
      Z_ion(j+1)  = ATI_result(1);   % average ionization state
      N_e(j+1) = n_gas * Z_ion(j+1);  % electron density (1/m^3)
      Energy_Ionization_e = ATI_result(2);  % ionization energy per electron (J)
      Energy_ATI_e = ATI_result(3);         % ATI energy per electron (J)
      Energy_Ionization(j+1) = Energy_Ionization_e * N_e(j+1) * pi * R_capillary^2 * dz;
      Energy_ATI(j+1) = Energy_ATI_e * N_e(j+1) * pi * R_capillary^2 * dz;
      kT(j+1) = Energy_ATI_e * (2/3);       % electron temperature (k_B times T) (J)
      a_IB(j+1) = InverseBremsstrahlungCoefficient(N_e(j+1), Z_ion(j+1), kT(j+1),omega_d);
      Energy_IB(j+1) = LaserEnergy(j+1) * a_IB(j+1) * dz;
      a_TS(j+1) = ThomsonScatteringCoefficient(N_e(j+1));
      Energy_TS(j+1) = LaserEnergy(j+1) * a_TS(j+1) * dz;
      Energy_capillary(j+1) = LaserEnergy(j+1) * a_capillary * dz;
      GVD(j+1) = GroupVelocityDispersion(omega_d,N_e(j+1));
      D(j+1) = D(j) + (GVD(j+1)+GVD_capillary)*dz;
      n_plasma_d(j+1) = sqrt(1-(q_e^2 * N_e(j+1))/(epsilon_0 * m_e * omega_d^2));
      C(j+1) = C(j) + dz/(c*n_plasma_d(j+1)) + (u_11^2*c)/(2*R_capillary^2*omega_d^2)*dz;
      k_d_plasma(j+1) = k_0 * n_plasma_d(j+1);
      k_d_total(j+1)  = k_d_plasma(j+1) + k_d_capillary(j+1);
      phi_d_prop(j+1) = phi_d_prop(j) + k_d_total(j+1)*dz;
      v_d(j+1)        = omega_d/k_d_total(j+1);
   end
   if FigureSwitch
      sgtitle('ionization at z = L');
   end
   
% --------------------------------------------------------------------
% Preserve raw ATI result before applying experimental corrections
Z_ion_raw = Z_ion;             % relative electron density from ATI theory
N_e_raw   = n_gas .* Z_ion;    % uncorrected electron density [m^{-3}]

% Calibrate theoretical plasma density with experimental average
%---------------------------------------------------------------------
   % f_avg = Ne_exp_avg / <N_e_raw(z)> ; Ne(z) = f_avg * N_e_raw(z)
   Ne_th_avg = trapz(z, N_e_raw) / L_capillary;      % model average [m^{-3}]
   f_avg = max(0, Ne_exp_avg / max(Ne_th_avg, eps)); % experimental correction
   Z_ion = max(0, f_avg .* Z_ion_raw);               % corrected charge state
   N_e   = n_gas .* Z_ion;                           % corrected electron density

   % Recompute plasma dependent quantities with calibrated density
   GVD = arrayfun(@(ne) GroupVelocityDispersion(omega_d, ne), N_e);
   n_plasma_d = sqrt(1 - (q_e^2 .* N_e) ./ (epsilon_0 * m_e * omega_d^2));
   k_d_plasma = k_0 .* n_plasma_d;
   k_d_total  = k_d_plasma + k_d_capillary;
   phi_d_prop = cumsum(k_d_total * dz);
   v_d        = omega_d ./ k_d_total;
   C = zeros(1, N);
   for j = 1:N-1
       C(j+1) = C(j) + dz/(c*n_plasma_d(j+1)) + (u_11^2*c)/(2*R_capillary^2*omega_d^2)*dz;
   end
   D = cumsum((GVD + GVD_capillary) * dz);
   
%---------------------------------------------------------------------

% Export results (SI units)
   result_SI_temp(1,:) = z;            % position (m)
   result_SI_temp(2,:) = LaserEnergy;  % laser pulse energy (J)
   result_SI_temp(3,:) = Z_ion;        % average ionization state
   result_SI_temp(4,:) = N_e;          % electron density N_e(z) (1/m^3, f_{avg} corrected)
   result_SI_temp(5,:) = GVD;          % group-velocity dispersion at z (sec^2/m)
   result_SI_temp(6,:) = D;            % accumulated group-delay dispersion at z (sec^2)
   result_SI_temp(7,:) = tau;          % laser pulse duration at z (sec)
   result_SI_temp(8,:) = I_peak;       % laser peak intensity at z (W/m^2)
   result_SI_temp(9,:) = E_peak;       % laser peak electric field at z (V/m)
   result_SI_temp(10,:) = kT;          % electron temperature (k_B times T) (J)
   result_SI_temp(11,:) = Energy_Ionization; % total ionization energy at (z,z+dz) (J)
   result_SI_temp(12,:) = Energy_ATI;        % total ATI energy at (z,z+dz) (J)
   result_SI_temp(13,:) = Energy_IB;         % total IB absorbed energy at (z,z+dz) (J)
   result_SI_temp(14,:) = Energy_TS;         % total Thomson scattered energy at (z,z+dz) (J)
   result_SI_temp(15,:) = Energy_capillary;  % capillary attenuated energy at (z,z+dz) (J)

   result_SI = result_SI_temp;
   
   fid_output = fopen('MainPulsePropagation_HHG_results_SI_units.txt','w');
   fprintf(fid_output,'position   laser_energy  ion_state_Z  e_density   GVD         GDD        pulse_duration  peak_intensity  peak_E_field  e_temp_kT  ionization  ATI_energy  IB_energy  TS_energy  capillary\r\n');
   fprintf(fid_output,'(m)        (J)           (arb.units)  (1/m^3)     (s^2/m)     (s^2)      (sec)           (W/m^2)         (V/m)         (J)        energy(J)   (J)         (J)        (J)        loss(J)\r\n');
   fprintf(fid_output,'%5.3E  %5.3E     %5.3f        %5.3E  % 5.3E  % 5.3E  %5.3E       %5.3E       %5.3E     %5.3E  %5.3E   %5.3E   %5.3E  %5.3E  %5.3E\r\n',result_SI);
   fclose(fid_output);

% Export results (usual units)
   result_temp(1,:) = z/mm;           % position (mm)
   result_temp(2,:) = LaserEnergy/mJ; % laser pulse energy (mJ)
   result_temp(3,:) = Z_ion;          % average ionization state
   result_temp(4,:) = N_e*cm^3;       % electron density N_e(z) (1/cm^3, f_{avg} corrected)
   result_temp(5,:) = GVD/fs^2*mm;    % group-velocity dispersion at z (fs^2/mm)
   result_temp(6,:) = D/fs^2;         % accumulated group-delay dispersion at z (fs^2)
   result_temp(7,:) = tau/fs;         % laser pulse duration at z (fs)
   result_temp(8,:) = I_peak*cm^2;    % laser peak intensity at z (W/cm^2)
   result_temp(9,:) = E_peak;         % laser peak electric field at z (V/m)
   result_temp(10,:) = kT/eV;         % electron temperature (k_B times T) (eV)
   result_temp(11,:) = Energy_Ionization/mJ; % total ionization energy at (z,z+dz) (mJ)
   result_temp(12,:) = Energy_ATI/mJ;        % total ATI energy at (z,z+dz) (mJ)
   result_temp(13,:) = Energy_IB/mJ;         % total IB absorbed energy at (z,z+dz) (mJ)
   result_temp(14,:) = Energy_TS/mJ;         % total Thomson scattered energy at (z,z+dz) (mJ)
   result_temp(15,:) = Energy_capillary/mJ;  % capillary attenuated energy at (z,z+dz) (mJ)

   fid_output = fopen('MainPulsePropagation_HHG_results_usual_units.txt','w');
   fprintf(fid_output,'position   laser_energy  ion_state_Z  e_density   GVD         GDD        pulse_duration  peak_intensity  peak_E_field  e_temp_kT  ionization  ATI_energy  IB_energy  TS_energy  capillary\r\n');
   fprintf(fid_output,'(mm)       (mJ)          (arb.units)  (1/cm^3)    (fs^2/mm)   (fs^2)     (fs)            (W/cm^2)        (V/m)         (eV)       energy(mJ)  (mJ)        (mJ)       (mJ)       loss(mJ)\r\n');
   fprintf(fid_output,'%6.3f     %6.3f        %5.3f        %5.3E  % 5.3E  % 5.3E  %6.3f          %5.3E       %5.3E     %5.3E  %5.3E   %5.3E   %5.3E  %5.3E  %5.3E\r\n',result_temp);
   fclose(fid_output);
   %%
% ================= Plot Driving Pulse Propagation =================
if FigureSwitch
    fig3 = figure('Color','w');
    tiledlayout(fig3,4,3,'Padding','compact','TileSpacing','compact');

    % 1) Laser energy vs z
    nexttile; 
    plot(z/mm, LaserEnergy/mJ, 'LineWidth',1.2);
    xlabel('z (mm)'); ylabel('Laser energy (mJ)');
    title('Pulse energy remaining');

    % 2) Electric field amplitude vs z
    nexttile; 
    plot(z/mm, E_peak, 'LineWidth',1.2);
    xlabel('z (mm)'); ylabel('E_{peak} (V/m)');
    title('Peak electric field amplitude');

    % 3) Intensity vs z
    nexttile; 
    plot(z/mm, I_peak*cm^2, 'LineWidth',1.2);
    xlabel('z (mm)'); ylabel('I_{peak} (W/cm^2)');
    title('Peak intensity');

    % 4) Electron density vs z
    nexttile; 
    plot(z/mm, N_e/cm^3, 'LineWidth',1.2);
    xlabel('z (mm)'); ylabel('N_e (cm^{-3})');
    legend('Electron density','Location','northeast');
    title('Plasma density from ionization');

    % 5) Group-velocity dispersion vs z
    nexttile; 
    plot(z/mm, GVD*mm/fs^2, 'LineWidth',1.2); hold on;
    plot(0, GVD_capillary*mm/fs^2, 'r*', 'MarkerSize',8);
    hold off;
    xlabel('z (mm)'); ylabel('GVD (fs^2/mm)');
    legend('Plasma','Waveguide','Location','northeast');
    title('Group-velocity dispersion contributions');

    % 6) Accumulated GDD vs z
    nexttile; 
    plot(z/mm, D/fs^2, 'LineWidth',1.2);
    xlabel('z (mm)'); ylabel('GDD (fs^2)');
    title('Accumulated group-delay dispersion');

    % 7) Pulse duration vs z
    nexttile; 
    plot(z/mm, tau/fs, 'LineWidth',1.2);
    xlabel('z (mm)'); ylabel('\tau (fs)');
    title('Temporal broadening of pulse');

    % 8) Ionization energy loss per length
    nexttile; 
    plot(z/mm, Energy_Ionization/dz, 'LineWidth',1.2);
    xlabel('z (mm)'); ylabel('Loss (mJ/mm)');
    title('Ionization energy loss');

    % 9) ATI energy loss per length
    nexttile; 
    plot(z/mm, Energy_ATI/dz, 'LineWidth',1.2);
    xlabel('z (mm)'); ylabel('Loss (mJ/mm)');
    title('Above-threshold ionization loss');

    % 10) Inverse Bremsstrahlung loss per length
    nexttile; 
    plot(z/mm, Energy_IB/dz, 'LineWidth',1.2);
    xlabel('z (mm)'); ylabel('Loss (mJ/mm)');
    title('Inverse Bremsstrahlung loss');

    % 11) Thomson scattering loss per length
    nexttile; 
    plot(z/mm, Energy_TS/dz, 'LineWidth',1.2);
    xlabel('z (mm)'); ylabel('Loss (mJ/mm)');
    title('Thomson scattering loss');

    % 12) Capillary attenuation loss per length
    nexttile; 
    plot(z/mm, Energy_capillary/dz, 'LineWidth',1.2);
    xlabel('z (mm)'); ylabel('Loss (mJ/mm)');
    title('Capillary guiding loss');

    % ---- Global figure title ----
    sgtitle({'Driving pulse propagation in capillary waveguide', ...
             'Part I, raw evolution'});
    set(gcf,'WindowState','maximized');  % R2018a+ for regular figures
    drawnow;

end

   %%
                    

% Curve fitting
   % E_peak(z) = E_peak_cfit(z) [V/m]
   E_peak_cfit = fit(z',E_peak','exp1');

   % I_peak(z) = I_peak_cfit(z) [W/m^2]
   I_peak_cfit = fit(z',I_peak','exp1');

   % tau(z) = tau_cfit(z) [s]
   tau_cfit = fit(z',(tau/fs)','poly3');
   tau_fun = @(z) tau_cfit(z)*fs;

   % N_e(z) = N_e_cfit(z) [m^{-3}]
   N_e_cfit = fit(z',N_e','poly3');

   % D(z) = D_cfit(z) [s^2]
   D_cfit = fit(z',(D/fs^2)','poly3');
   D_fun = @(z) D_cfit(z)*fs^2;

   % n_plasma_d(z) = n_plasma_d_cfit(z)
   n_plasma_d_cfit = fit(z',n_plasma_d','poly3');

   % C(z) = C_cfit(z) [s]
   C_cfit = fit(z',(C/fs)','poly3');       % unit: fs
   C_fun = @(z) C_cfit(z)*fs;              % unit: sec
   
   % total wavenumber of the driving pulse k_d_total_cfit(z)
   k_d_total_cfit = fit(z',k_d_total','poly3');    % unit: 1/m
   
   % accumulated phase shift of the driving pulse due to propagation (rad)
   phi_d_prop_cfit = fit(z',phi_d_prop','poly3');    % unit: rad

   % phase velocity of the driving pulse v_d_cfit(z)
   v_d_cfit = fit(z',v_d','poly3');    % unit: m/sec
%%
if FigureSwitch
    fig_fit = figure('Color','w');
    tiledlayout(fig_fit,4,3,'Padding','compact','TileSpacing','compact');

    % Convenience for dynamic title
    try
        gas_str = Gas;
    catch, gas_str = 'Gas'; end
    try
        sgt = sprintf(['Driving pulse propagation — fits vs data\n' ...
                       '\\lambda = %.0f nm, R_{cap} = %.0f \\mu m, L_{cap} = %.1f mm, gas = %s'], ...
                       lambda*1e9, R_capillary*1e6, L_capillary*1e3, gas_str);
    catch
        sgt = 'Driving pulse propagation — fits vs data';
    end

    lw_fit  = 1.6; % line width for fits
    lw_data = 1.2; % line width for raw data

    % (1,1) Peak electric field
    nexttile; 
    plot(z, E_peak_cfit(z), '--', 'LineWidth', lw_fit); hold on;
    plot(z, E_peak,          '-',  'LineWidth', lw_data); hold off; grid on;
    ylabel('E_{peak} (V/m)');
    title('Peak electric field');
    legend('fit','data','Location','best');

    % (2,1) Peak intensity (W/cm^2)
    nexttile;
    plot(z, I_peak_cfit(z)*cm^2, '--', 'LineWidth', lw_fit); hold on;
    plot(z, I_peak*cm^2,          '-',  'LineWidth', lw_data); hold off; grid on;
    ylabel('I_{peak} (W/cm^2)');
    title('Peak intensity');
    legend('fit','data','Location','best');

    % (3,1) Electron density (cm^{-3})
    nexttile;
    plot(z, N_e_cfit(z)*cm^3, '--', 'LineWidth', lw_fit); hold on;
    plot(z, N_e*cm^3,          '-',  'LineWidth', lw_data); hold off; grid on;
    ylabel('N_e (cm^{-3})');
    title('Electron density');
    legend('fit','data','Location','best');

    % (1,2) Plasma refractive index
    nexttile;
    plot(z, n_plasma_d_cfit(z), '--', 'LineWidth', lw_fit); hold on;
    plot(z, n_plasma_d,          '-',  'LineWidth', lw_data); hold off; grid on;
    ylabel('n_{plasma}');
    title('Plasma refractive index');
    legend('fit','data','Location','best');

    % (2,2) Group-delay dispersion D (fs^2)
    % D_cfit fitted on D/fs^2, so evaluate then multiply
    nexttile;
    plot(z, D_cfit(z), '--', 'LineWidth', lw_fit); hold on;        % already fs^2
    plot(z, D/fs^2,    '-',  'LineWidth', lw_data); hold off; grid on;
    ylabel('D (fs^2)');
    title('Group-delay dispersion (GDD)');
    legend('fit','data','Location','best');

    % (3,2) Pulse duration (fs)
    nexttile;
    plot(z, tau_cfit(z), '--', 'LineWidth', lw_fit); hold on;       % already fs
    plot(z, tau/fs,      '-',  'LineWidth', lw_data); hold off; grid on;
    ylabel('\tau (fs)');
    title('Pulse duration');
    legend('fit','data','Location','best');

    % (4,2) Group delay C (fs)
    nexttile;
    plot(z, C_cfit(z), '--', 'LineWidth', lw_fit); hold on;         % already fs
    plot(z, C/fs,       '-',  'LineWidth', lw_data); hold off; grid on;
    xlabel('position z (m)'); ylabel('C (fs)');
    title('Group delay');
    legend('fit','data','Location','best');

    % (1,3) Total wavenumber k_d_total (1/m)
    nexttile;
    plot(z, k_d_total_cfit(z), '--', 'LineWidth', lw_fit); hold on;
    plot(z, k_d_total,          '-',  'LineWidth', lw_data); hold off; grid on;
    ylabel('k_{d,total} (1/m)');
    title('Total wavenumber');
    legend('fit','data','Location','best');

    % (2,3) Accumulated phase from propagation \phi_{d,prop} (rad)
    nexttile;
    plot(z, phi_d_prop_cfit(z), '--', 'LineWidth', lw_fit); hold on;
    plot(z, phi_d_prop,          '-',  'LineWidth', lw_data); hold off; grid on;
    ylabel('\phi_{d,prop} (rad)');
    title('Accumulated propagation phase');
    legend('fit','data','Location','best');

    % (3,3) Phase velocity v_d (m/s)
    nexttile;
    plot(z, v_d_cfit(z), '--', 'LineWidth', lw_fit); hold on;
    plot(z, v_d,          '-',  'LineWidth', lw_data); hold off; grid on;
    ylabel('v_d (m/s)');
    title('Phase velocity');
    legend('fit','data','Location','best');

    % Bottom-right tile left blank -> place a summary note if you like:
    nexttile; axis off;
    text(0,0.9,'Notes:', 'FontWeight','bold');
    text(0,0.7,'• Fits: dashed; Data: solid');
    text(0,0.5,'• Units shown on each axis');
    text(0,0.3,'• Same z-range (meters) across panels');

    sgtitle(sgt, 'Interpreter','tex');
    set(gcf,'WindowState','maximized');  % R2018a+ for regular figures
    drawnow;
end

   %%
% Part II: HHG calculation ------------------------------------------------

% Set the coordinate and variables
   lambda_q = lambda/q;        % HHG wavelength (m)
   omega_q = q * omega_d;      % HHG angular frequency (rad/sec)
   z2 = [0:dz2:L_capillary];   % Generate coordinate z2 (unit: m)
   N2 = size(z2,2);            % number of data points

% Dipole phase calculation
   % Find the minimum intensity of the driving laser:
   % U_q = q hbar omega_d = I_p + 3.17 Up(I_min)
   U_q = q*hbar*omega_d;    % HHG photon energy (J) 
   I_min = (U_q-I_p)/3.17*(2*m_e*omega_d^2)/(q_e^2*mu_0*c);  % unit: W/m^2
   
   % dipole phase Phi_dipole(I_d2) as a function of laser intensity I
   N_dipole = 1000;    % number of data points
   % laser intensity range for dipole phase calculation (W/m^2)
   %I_dipole = linspace(I_min,I_min*10,N_dipole);
   I_dipole = linspace(min(I_peak)*0.9,max(I_peak)*1.1,N_dipole);  
   E_dipole = sqrt(2*mu_0*c*I_dipole); % laser amplitude (V/m)
   Phi_dipole_l = zeros(1,N_dipole);  % long-trajectory dipole phase (rad)
   Phi_dipole_s = zeros(1,N_dipole);  % short-trajectory dipole phase (rad)
   for j=1:N_dipole
       temp = HHG_DipolePhase(q,E_dipole(j),I_p,omega_d);
       Phi_dipole_l(j) = temp(1);   % long-trajectory dipole phase (rad)
       Phi_dipole_s(j) = temp(2);   % short-trajectory dipole phase (rad)
   end
   
   % alpha coefficient
   alpha_l = zeros(1,N_dipole);  % unit: m^2/W
   alpha_s = zeros(1,N_dipole);  % unit: m^2/W
   dI = I_dipole(2)-I_dipole(1);
   alpha_l(1:N_dipole-1) = (Phi_dipole_l(2:N_dipole)-Phi_dipole_l(1:N_dipole-1))/dI;
   alpha_l(N_dipole) = 2*alpha_l(N_dipole-1) - alpha_l(N_dipole-2);
   alpha_s(1:N_dipole-1) = (Phi_dipole_s(2:N_dipole)-Phi_dipole_s(1:N_dipole-1))/dI;
   alpha_s(N_dipole) = 2*alpha_s(N_dipole-1) - alpha_s(N_dipole-2);
   
% curve fitting of the dipole phase and the alpha coefficient
   %I_dipole_l = I_dipole; 
   Phi_dipole_l_temp = Phi_dipole_l; 
   [I_dipole_l, Phi_dipole_l] = cleanNAN(I_dipole, Phi_dipole_l_temp); 

   %I_dipole_s = I_dipole; 
   Phi_dipole_s_temp = Phi_dipole_s; 
   [I_dipole_s, Phi_dipole_s] = cleanNAN(I_dipole, Phi_dipole_s_temp); 

   Phi_dipole_l_cfit = fit(I_dipole_l',Phi_dipole_l','fourier3');  % unit: rad
   Phi_dipole_s_cfit = fit(I_dipole_s',Phi_dipole_s','power2');    % unit: rad
   
   alpha_l_temp = alpha_l; 
   [I_dipole_l, alpha_l] = cleanNAN(I_dipole, alpha_l_temp); 

   alpha_s_temp = alpha_s; 
   [I_dipole_s, alpha_s] = cleanNAN(I_dipole, alpha_s_temp); 
   
   
   alpha_l_cfit = fit(I_dipole_s',alpha_l','poly5');       % unit: m^2/W
   alpha_s_cfit = fit(I_dipole_l',alpha_s','poly5');       % unit: m^2/W
   %%
   if FigureSwitch
    fig_dip = figure('Color','w');
    tiledlayout(fig_dip,1,4,'Padding','compact','TileSpacing','compact');

    % intensity in W/cm^2
    I_dipole_cm2_l = I_dipole_l * 1e-4;
    I_dipole_cm2_s = I_dipole_s * 1e-4;

    % (1) Long-trajectory dipole phase
    nexttile;
    plot(I_dipole_cm2_l, Phi_dipole_l, 'o', 'LineWidth',1.2); hold on;
    plot(I_dipole_cm2_l, Phi_dipole_l_cfit(I_dipole_l), 'r--', 'LineWidth',1.6); hold off; grid on;
    xlabel('Intensity (W/cm^2)'); ylabel('\Phi_{dipole,long} (rad)');
    title('Long-trajectory dipole phase');
    legend('data','fit','Location','best');

    % (2) Short-trajectory dipole phase
    nexttile;
    plot(I_dipole_cm2_s, Phi_dipole_s, 'o', 'LineWidth',1.2); hold on;
    plot(I_dipole_cm2_s, Phi_dipole_s_cfit(I_dipole_s), 'r--', 'LineWidth',1.6); hold off; grid on;
    xlabel('Intensity (W/cm^2)'); ylabel('\Phi_{dipole,short} (rad)');
    title('Short-trajectory dipole phase');
    legend('data','fit','Location','best');

    % (3) Long-trajectory alpha coefficient
    nexttile;
    plot(I_dipole_cm2_l, alpha_l, 'o', 'LineWidth',1.2); hold on;
    plot(I_dipole_cm2_l, alpha_l_cfit(I_dipole_s), 'r--', 'LineWidth',1.6); hold off; grid on;
    xlabel('Intensity (W/cm^2)'); ylabel('\alpha_{long} (m^2/W)');
    title('Long-trajectory \alpha');
    legend('data','fit','Location','best');

    % (4) Short-trajectory alpha coefficient
    nexttile;
    plot(I_dipole_cm2_s, alpha_s, 'o', 'LineWidth',1.2); hold on;
    plot(I_dipole_cm2_s, alpha_s_cfit(I_dipole_l), 'r--', 'LineWidth',1.6); hold off; grid on;
    xlabel('Intensity (W/cm^2)'); ylabel('\alpha_{short} (m^2/W)');
    title('Short-trajectory \alpha');
    legend('data','fit','Location','best');

    % ---- Super title ----
    sgtitle('Dipole phase and \alpha coefficients vs intensity');
    set(gcf,'WindowState','maximized');  % R2018a+ for regular figures
    drawnow;
end

   %%
% wavenumber of the driving pulse   (1-by-N2 array)
   % wavenumber of the driving pulse in plasma (1/m)
   k_d_plasma2 = k_0 * n_plasma_d_cfit(z2)';
   
   % total wavenumber of the driving pulse (1/m)
   k_d_total2 = k_d_plasma2 + k_d_capillary;
   
   % phase velocity of the driving pulse (m/sec)
   v_d2 = omega_d./k_d_total2;

% wavenumber of the HHG   (1-by-N2 array)
   % plasma refractive index of the HHG
   n_plasma_q = sqrt(1-(q_e^2 * N_e_cfit(z2)')/(epsilon_0*m_e*omega_q^2));
   
   % wavenumber of the driving pulse in plasma (1/m)
   k_q_plasma = q * k_0 * n_plasma_q;
   
   % wavenumber variation of the HHG in the capillary (1/m)
   k_q_capillary = -(u_11^2*lambda_q)/(4*pi*R_capillary);
   
   % total wavenumber of the HHG (1/m)
   k_q_total = k_q_plasma + k_q_capillary;
   
   % phase velocity of the HHG (m/sec)
   v_q = omega_q./k_q_total;
   
% wavenumber mismatch due to plasma dispersion (1-by-N2 array)
   Delta_k_plasma = q*k_d_plasma2 - k_q_plasma;
   
% wavenumber mismatch due to wavguide dispersion (1-by-N2 array)
   Delta_k_capillary = q*k_d_capillary - k_q_capillary + 0.0*Delta_k_plasma;
   
% wavenumber mismatch due to dipole phase (1-by-N2 array)
   % peak intensity of the driving pulse as a function of z (1-by-N2 array)
   I_d2 = I_peak_cfit(z2)';
   % peak intensity variation dI/dz as a function of z (1-by-N2 array)
   dI_dz = (I_peak_cfit(z2+dz2)' - I_peak_cfit(z2)') / dz2;
   % wavenumber mismatch of the long-trajectory dipole phase as a function of z (1-by-N2 array)
   Delta_k_dipole_l = alpha_l_cfit(I_d2)' .* dI_dz;
   % wavenumber mismatch of the short-trajectory dipole phase as a function of z (1-by-N2 array)
   Delta_k_dipole_s = alpha_s_cfit(I_d2)' .* dI_dz;

% total wavenumber mismatch (1-by-N2 array)
   Delta_k_total_l = Delta_k_plasma + Delta_k_capillary + Delta_k_dipole_l;
   Delta_k_total_s = Delta_k_plasma + Delta_k_capillary + Delta_k_dipole_s;

   
   
% dephasing length (1-by-N2 array)
   L_dephasing_l = pi ./ abs(Delta_k_total_l);  % unit: m
   L_dephasing_s = pi ./ abs(Delta_k_total_s);  % unit: m
   
% accumulated phase mismatch (rad)
   Delta_Phi_plasma    = cumsum(Delta_k_plasma)*dz2;     
   Delta_Phi_capillary = cumsum(Delta_k_capillary)*dz2;
   Delta_Phi_dipole_l  = cumsum(Delta_k_dipole_l)*dz2;
   Delta_Phi_dipole_s  = cumsum(Delta_k_dipole_s)*dz2;
   Delta_Phi_total_l   = cumsum(Delta_k_total_l)*dz2;
   Delta_Phi_total_s   = cumsum(Delta_k_total_s)*dz2;

   figure; sgtitle('The wavefront is at t=0')
   subplot(5,3,1), plot(z2/mm,Delta_k_plasma);
      xlabel('z (mm)');
      ylabel('\Delta k_{plasma} (1/m)');
      hold on, plot(0,Delta_k_capillary,'r*'), hold off;
      legend('plasma','waveguide','Location','NorthEast');
      title('wavenumber mismatch');
   subplot(5,3,4), plot(z2/mm,Delta_k_dipole_l);
      xlabel('z (mm)');
      ylabel('\Delta k_{dipole\_long} (1/m)');
      legend('dipole\_long');
   subplot(5,3,7), plot(z2/mm,Delta_k_dipole_s);
      xlabel('z (mm)');
      ylabel('\Delta k_{dipole\_short} (1/m)');
      legend('dipole\_short');
   subplot(5,3,10), plot(z2/mm,Delta_k_total_l);
      xlabel('z (mm)');
      ylabel('\Delta k_{total\_long} (1/m)');
      legend('total\_long');
   subplot(5,3,13), plot(z2/mm,Delta_k_total_s);
      xlabel('z (mm)');
      ylabel('\Delta k_{total\_short} (1/m)');
      legend('total\_short');
   subplot(5,3,2), plot(z2/mm,Delta_Phi_plasma/pi,z2/mm,Delta_Phi_capillary/pi,'r*');
      xlabel('z (mm)');
      ylabel('\Delta \Phi (\pi)');
      legend('plasma','waveguide','Location','NorthEast');
      title('accumulated phase mismatch');
   subplot(5,3,5), plot(z2/mm,Delta_Phi_dipole_l/pi);
      xlabel('z (mm)');
      ylabel('\Delta \Phi_{dipole\_long} (\pi)');
      legend('dipole\_long');
   subplot(5,3,8), plot(z2/mm,Delta_Phi_dipole_s/pi);
      xlabel('z (mm)');
      ylabel('\Delta \Phi_{dipole\_short} (\pi)');
      legend('dipole\_short');
   subplot(5,3,11), plot(z2/mm,Delta_Phi_total_l/pi);
      xlabel('z (mm)');
      ylabel('\Delta \Phi_{total\_long} (\pi)');
      legend('total\_long');
   subplot(5,3,14), plot(z2/mm,Delta_Phi_total_s/pi);
      xlabel('z (mm)');
      ylabel('\Delta \Phi_{total\_short} (\pi)');
      legend('total\_short');
   subplot(5,3,3), plot(z2/mm,L_dephasing_l/mm);
      xlabel('z (mm)');
      ylabel('L_{dephasing\_long} (mm)');
      legend('long');
      title('dephasing length');
   subplot(5,3,6), plot(z2/mm,L_dephasing_s/mm);
      xlabel('z (mm)');
      ylabel('L_{dephasing\_short} (mm)');
      legend('short');
   set(gcf,'WindowState','maximized');  % R2018a+ for regular figures
   drawnow;

  % Export results (usual units)
   result_Dispersion(1,:) = z2/mm;                      % position (mm)
   result_Dispersion(2,:) = Delta_Phi_plasma/pi;       % plasma dispersion (rad/pi)
   result_Dispersion(3,:) = Delta_Phi_capillary/pi;    % capillary dispersion (rad/pi)
   result_Dispersion(4,:) = Delta_Phi_dipole_l/pi;     % dipole phase - long (rad/pi)
   result_Dispersion(5,:) = Delta_Phi_dipole_s/pi;     % dipole phase - long (rad/pi)
   result_Dispersion(6,:) = Delta_Phi_total_l/pi;      % total phase - long (rad/pi)
   result_Dispersion(7,:) = Delta_Phi_total_s/pi;      % total phase - short (rad/pi))
   result_Dispersion(8,:) = L_dephasing_l/mm;          % dephasing length (mm)
   result_Dispersion(9,:) = L_dephasing_s/mm;          % dephasing length (mm)
 

   fid_output = fopen('DispersionSource.txt','w');
   fprintf(fid_output,'position   Delta_Phi_plasma  Delta_Phi_capillary  Delta_Phi_dipole_l   Delta_Phi_dipole_s         Delta_Phi_total_l        Delta_Phi_total_s     L_dephasing_l   L_dephasing_s  \r\n');
   fprintf(fid_output,'(mm)       (rad/pi)          (rad/pi)             (rad/pi)             (rad/pi)                   (rad/pi)                 (rad/pi)              (mm)            (mm)           \r\n');
   fprintf(fid_output,'%5.9e       %5.9e            %5.9e                %5.9e                %5.9e                      %5.9e                    %5.9e                 %5.9e           %5.9e           \r\n',result_Dispersion);
   fclose(fid_output);     
% HHG calculation


% index i for traing different wavefront 
% ---- Safe timing shift used ONLY in Part II (HHG observer) ----
% Negative t_shift => waveform earlier. Adds to C so (t-C) is unchanged.
% ElectricField_shifted = @(t,Epk,om,tau0,D,C,phi) ...
%     ElectricField_d(t, Epk, om, tau0, D, C + t_shift, phi);


for ii = 150:150  
% Trace a fixed wavefront of the driving pulse
   % propagation time interval in each section dz2
   dt2_dwf = dz2 ./ v_d2;
   % arrival time t_d(z2) to position z2
   t2_dwf  = cumsum(dt2_dwf)-dt2_dwf(1) + (ii-150)*0.1*fs;  % + is backward, - is forward 
   % the driving field E_d2 at (z2,t=t_d2(z2))
   E_d2_dwf = zeros(1,N2);
   
% Trace a fixed wavefront of the HHG
   % propagation time interval in each section dz2
   dt2_qwf = dz2 ./ v_q;
   % arrival time t_q(z2) to position z2
   t2_qwf  = cumsum(dt2_qwf)-dt2_qwf(1) + (ii-150)*0.1*fs; %(pi/omega_q);
   % the driving field E_d2 at (z2,t=t_q(z2))
   E_d2_qwf = zeros(1,N2);
   z2_qwf = t2_qwf.*v_q; 
   
   for j = 1:N2
      % Generate the driving field E_d2 at (z2,t=t_d2(z2)) E_d2(z)
      E_d2_dwf(j) = ElectricField_d(t2_dwf(j),E_peak_cfit(z2(j)),omega_d,...
                                tau_0,D_fun(z2(j))-D_fun(z2(1)),...
                                C_fun(z2(j))-C_fun(z2(1)),...
                                phi_d_prop_cfit(z2(j))-phi_d_prop_cfit(z2(1)));
      % Generate the driving field E_d2 at (z2,t=t_q(z2))
      E_d2_qwf(j) = ElectricField_d(t2_qwf(j),E_peak_cfit(z2(j)),omega_d,...
                                tau_0,D_fun(z2(j))-D_fun(z2(1)),...
                                C_fun(z2(j))-C_fun(z2(1)),...
                                phi_d_prop_cfit(z2(j))-phi_d_prop_cfit(z2(1)));

   end
   I_d2_dwf = abs(E_d2_dwf).^2/(2*mu_0*c);    % driving intensity with fixed driving wavefront (W/m^2)
   I_d2_qwf = abs(E_d2_qwf).^2/(2*mu_0*c);    % driving intensity with fixed harmonic wavefront (W/m^2)
   Phi_d2_dwf = angle(E_d2_dwf);   % driving field phase with fixed driving wavefront
   Phi_d2_qwf = angle(E_d2_qwf);   % driving field phase with fixed harmonic wavefront

 
% % --------------------------------------------------------------------------

% % ===== Intensity-ratio wavefront density model =====
% % Calibrated peak density profile on z2 (already includes f_peak and f_avg
% % because N_e was fitted after calibration):
% N_e_peak_z2 = N_e_cfit(z2)';         % [1/m^3] calibrated peak density
% 
% % Intensities:
% I_peak_z2   = I_peak_cfit(z2)';      % [W/m^2] peak intensity along z2
% I_qwf       = abs(E_d2_qwf).^2 / (2*mu_0*c);   % [W/m^2] already computed
% 
% % Ratio r \in [0,1]
% r = I_qwf ./ max(I_peak_z2, eps);
% r = max(0, min(1, r));
% 
% % Nonlinearity (tunable). Start with 1.5; use 1.0 if you want linear scaling.
% gamma = 1;
% 
% % Wavefront electron density (heuristic):
% N_e_qwf = N_e_peak_z2 .* (r.^gamma);   % [1/m^3]

% If (for any reason) you prefer to explicitly carry the same calibration
% constants here instead of relying on N_e_cfit being calibrated, you can use:
% Ne_qwf = (f_avg * f_peak) .* (n_gas .* Z_ion_cfit(z2)') .* (r.^gamma);
% (only if you have a Z_ion fit; otherwise the first line is preferred)
% ================================================

N_e_qwf = zeros(N2, 1);


% phase of the local harmonic field Phi_HHG(z) = q*Phi_drive(z) + dipole
% Trace a fixed harmonic wavefront
   % long-trajectory emission
   Phi_LH_l = q*Phi_d2_qwf + Phi_dipole_l_cfit(I_d2_qwf)';
   % short-trajectory emission
   Phi_LH_s = q*Phi_d2_qwf + Phi_dipole_s_cfit(I_d2_qwf)';
   
% Calculate the local harmonic field
   % long-trajectory emission
   %Phi_LH_l = Phi_LH_l - Phi_LH_l(1);
   % Ionization at t_shift wavefront
   E_t2 = zeros(N,size(time,2));
   ionizationResult = zeros(8,size(time,2));
   n_0 = zeros(N,size(time,2));
   n_1 = zeros(N,size(time,2));
   n_2 = zeros(N,size(time,2));
   n_3 = zeros(N,size(time,2));  
   n_e = zeros(N,size(time,2));   
   delta_tw = zeros(N,1);
   n_0_qwf = zeros(N,1); 
   n_1_qwf = zeros(N,1);    
   n_2_qwf = zeros(N,1);   
   n_3_qwf = zeros(N,1);    
   E_t2_qwf = zeros(N,1);
   for jj = 1:N     % propagation along z
        % E_t2(z, t)
        E_t2(jj,:) = ElectricField_d(time+C_fun(z2(jj))-C_fun(z2(1)),E_peak(jj),omega_d,tau_0,D(jj)-D(1),C_fun(z2(jj))-C_fun(z2(1)),0); 
        % E_t      = ElectricField_d(time,                           E_peak(1) ,omega_d,tau_0,0         ,0                         ,0);
        ionizationResult = TunnelingIonizationRate_Linear2((E_t2(jj,:)),omega_d,time+C_fun(z2(jj)),E_ion_0,E_ion_1,E_ion_2,E_ion_3,E_ion_4,FigureSwitch);
        ionizationResult = ionizationResult';
        n_e(jj,:) = ionizationResult(1,:);        
        n_0(jj,:) = ionizationResult(2,:);
        n_1(jj,:) = ionizationResult(3,:);
        n_2(jj,:) = ionizationResult(4,:);
        n_3(jj,:) = ionizationResult(5,:);  
        % the observer time at each z for the harmonic wavefront
        % delta_tw(jj) = C_fun(z2(jj))-C_fun(z2(1)) - t2_qwf(jj);
        % C_fun: relative delay of the pulse envelope center. Subtracting t2_qwf(jj) aligns your laser field’s “local clock” with the harmonic wavefront clock.
        delta_tw(jj) = (C_fun(z2(jj)) - C_fun(z2(1)) - t_shift) - t2_qwf(jj); 
        % If delta_tw(jj) = 0, then idx = center index → the pulse peak.
        idx = round(0.5*size(time,2)) - round(delta_tw(jj)/DeltaTime);
        n_0_qwf(jj) = n_0(jj, idx);
        n_1_qwf(jj) = n_1(jj, idx);
        n_2_qwf(jj) = n_2(jj, idx);
        n_3_qwf(jj) = n_3(jj, idx);
        E_t2_qwf(jj) =  E_t2(jj, idx);
   end
   % The plasma density at t_shift
   Z_ion_qwf = n_1_qwf + 2*n_2_qwf + 3*n_3_qwf; 
   N_e_qwf = n_gas .* Z_ion_qwf' .* f_avg; 
  
 %%
% neutral ~ 1+
   switch I_p
    case E_ion_0
        W_n_1_qwf = StaticIonizationRate(E_ion_0, abs(E_d2_qwf));
        n_source  = n_gas .* n_0_qwf' .* W_n_1_qwf;
    case E_ion_1
        W_n_1_qwf = StaticIonizationRate(E_ion_1, abs(E_d2_qwf));
        n_source  = n_gas .* n_1_qwf' .* W_n_1_qwf;
    case E_ion_2
        W_n_1_qwf = StaticIonizationRate(E_ion_2, abs(E_d2_qwf));
        n_source  = n_gas .* n_2_qwf' .* W_n_1_qwf;
    otherwise
        error('Unsupported ionization potential I_p = %g', I_p);
    end

    E_LH_l = n_source .* abs(E_d2_qwf).^5 .* exp(1i*Phi_LH_l);    % arb. units
    % short-trajectory emission
    E_LH_s = n_source .* abs(E_d2_qwf).^5 .* exp(1i*Phi_LH_s);    % arb. units

% Calculate the accumulated harmonic field (Bug: n_source is at t_shift, but driving field is at t0)
   E_HHG_l = cumsum(E_LH_l);
   E_HHG_s = cumsum(E_LH_s);
   E_q_PhaseMatched = n_source.*abs(E_d2_qwf).^5;% phase matched HHG field (as a function of z, column vector)(Z_ion - 1)
   E_q_PhaseMatched_final = cumsum(E_q_PhaseMatched);
   E_q_PhaseMatched_final_max = max(E_q_PhaseMatched_final); %max(E_q_PhaseMatched_final);
%%   
   figure; % driving field check
   subplot(5,1,1), plot(z2/mm,abs(E_d2_dwf),z2/mm,abs(E_d2_qwf),z/mm,E_peak);
      xlabel('z (mm)');
      ylabel('|E| (V/m)');
      legend('fixed d-wf','fixed q-wf','E_{peak}','Location','SouthWest');
      title('driving field amplitude')
   subplot(5,1,2)
        plot(z2/mm, (t2_dwf - (C_fun(z2)'-C_fun(0)' - t_shift))/fs, ...
             z2/mm, (t2_qwf - (C_fun(z2)'-C_fun(0)' - t_shift))/fs, 'LineWidth',1.2);
        xlabel('z (mm)');
        ylabel('\Delta t (fs)');
        legend('t2_{dwf} - C(z)','t2_{qwf} - C(z)','Location','SouthWest');
        title('time difference');
   subplot(5,1,3), plot(z2/mm,Phi_d2_dwf,z2/mm,Phi_d2_qwf);
      xlabel('z (mm)');
      ylabel('\Phi_{d2} (rad)');
      legend('fixed d-wf','fixed q-wf','Location','SouthWest');
      title('phase of the driving pulse')
   subplot(5,1,4), plot(z2/mm,n_source*1E-6);
      xlabel('z (mm)');
      ylabel('Source Rate(Hz cm^{-3})');
      title('Source');
   subplot(5,1,5), plot(z2/mm,W_n_1_qwf);
      xlabel('z (mm)');
      ylabel('Ionization rate (Hz)');
      title('Ionization rate');
   set(gcf,'WindowState','maximized');  % R2018a+ for regular figures
    drawnow;
 %%  
  %% ================= q-th HHG wavefront diagnostics (baseline: t_shift = 0) =================
% Purpose: Same panels as Figure 10, but evaluated strictly at the q-wavefront with t_shift = 0.
% This lets you compare "t = 0" vs your shifted observer (Figure 10).

% --- 0) Build the q-wf driving field with NO observer shift ---
E_d2_qwf0 = zeros(1, N2);
for j = 1:N2
    % note: NO -t_shift here
    E_d2_qwf0(j) = ElectricField_d(t2_qwf(j), ...
        E_peak_cfit(z2(j)), omega_d, ...
        tau_0, D_fun(z2(j)) - D_fun(z2(1)), ...
        C_fun(z2(j)) - C_fun(z2(1)), ...
        phi_d_prop_cfit(z2(j)) - phi_d_prop_cfit(z2(1)));
end
I_qwf0  = abs(E_d2_qwf0).^2 / (2*mu_0*c);    % [W/m^2]
Phi_d2_qwf0 = angle(E_d2_qwf0);              % driving phase at q-wf, t_shift = 0

% --- 1) Wavefront-aligned populations at t_shift = 0 (indexing via delta_tw0) ---
% Reuse your ionizationResult engine, sample the same way you did but WITHOUT t_shift.
n_0_qwf0 = zeros(N2,1); n_1_qwf0 = zeros(N2,1);
n_2_qwf0 = zeros(N2,1); n_3_qwf0 = zeros(N2,1);
E_t2_qwf0 = zeros(N2,1);

for jj = 1:N2
    % local field trace at z2(jj), same as before
    Et_local = ElectricField_d(time + C_fun(z2(jj)) - C_fun(z2(1)), ...
        E_peak_cfit(z2(jj)), omega_d, ...
        tau_0, D_fun(z2(jj)) - D_fun(z2(1)), ...
        C_fun(z2(jj)) - C_fun(z2(1)), 0);
    ionRes = TunnelingIonizationRate_Linear2(Et_local, omega_d, ...
        time + C_fun(z2(jj)), E_ion_0, E_ion_1, E_ion_2, E_ion_3, E_ion_4, FigureSwitch).';
    % Align observer to q-wavefront WITHOUT t_shift:
    delta_tw0 = (C_fun(z2(jj)) - C_fun(z2(1))) - t2_qwf(jj);
    idx0 = round(0.5*size(time,2)) - round(delta_tw0/DeltaTime);
    idx0 = max(1, min(length(time), idx0)); % clamp for safety

    n_0_qwf0(jj) = ionRes(2, idx0);
    n_1_qwf0(jj) = ionRes(3, idx0);
    n_2_qwf0(jj) = ionRes(4, idx0);
    n_3_qwf0(jj) = ionRes(5, idx0);
    E_t2_qwf0(jj) = Et_local(idx0);
end

Z_ion_qwf0 = n_1_qwf0 + 2*n_2_qwf0 + 3*n_3_qwf0;
Ne_qwf0    = n_gas .* Z_ion_qwf0' .* f_avg;     % [1/m^3] calibrated wavefront density (t_shift = 0)

% --- 2) Plasma/total phases along q-wf (t_shift = 0) ---
n_plasma_d_qwf0 = sqrt(1 - (q_e^2 .* Ne_qwf0) ./ (epsilon_0 * m_e * omega_d^2));
n_plasma_q_qwf0 = sqrt(1 - (q_e^2 .* Ne_qwf0) ./ (epsilon_0 * m_e * (q*omega_d)^2));

k_d_plasma_qwf0 = k_0     .* n_plasma_d_qwf0;
k_q_plasma_qwf0 = (q*k_0) .* n_plasma_q_qwf0;

k_d_total_qwf0 = k_d_plasma_qwf0 + k_d_capillary;  % same capillary term
k_q_total_qwf0 = k_q_plasma_qwf0 + k_q_capillary;

Delta_k_plasma_qwf0 = q*k_d_plasma_qwf0 - k_q_plasma_qwf0;   % [1/m]
Phi_plasma_qwf0     = cumsum(Delta_k_plasma_qwf0) * dz2;     % [rad]

Phi_dip_long0     = Phi_dipole_l_cfit(I_qwf0)';              % [rad]
Phi_dip_long0_rel = Phi_dip_long0 - Phi_dip_long0(1);

Phi_total_long_qwf0 = Phi_plasma_qwf0 + Delta_Phi_capillary + Phi_dip_long0;
Phi_total_long_qwf0_rel = Phi_total_long_qwf0 - Phi_total_long_qwf0(1);

% --- 3) Source & ionization rate (q-wf, t_shift = 0) ---
% pick the correct ground/state for I_p:
switch I_p
    case E_ion_0
        W_qwf0 = StaticIonizationRate(E_ion_0, abs(E_d2_qwf0));
        n_source0 = n_gas .* n_0_qwf0' .* W_qwf0;
    case E_ion_1
        W_qwf0 = StaticIonizationRate(E_ion_1, abs(E_d2_qwf0));
        n_source0 = n_gas .* n_1_qwf0' .* W_qwf0;
    case E_ion_2
        W_qwf0 = StaticIonizationRate(E_ion_2, abs(E_d2_qwf0));
        n_source0 = n_gas .* n_2_qwf0' .* W_qwf0;
    otherwise
        error('Unsupported ionization potential I_p = %g', I_p);
end

% --- 4) Accumulated HHG (long/short) at q-wf, t_shift = 0 ---
Phi_LH_l0 = q*Phi_d2_qwf0 + Phi_dipole_l_cfit(I_qwf0)';   % long
Phi_LH_s0 = q*Phi_d2_qwf0 + Phi_dipole_s_cfit(I_qwf0)';   % short

E_LH_l0 = n_source0 .* abs(E_d2_qwf0).^5 .* exp(1i*Phi_LH_l0);
E_LH_s0 = n_source0 .* abs(E_d2_qwf0).^5 .* exp(1i*Phi_LH_s0);

E_HHG_l0 = cumsum(E_LH_l0);
E_HHG_s0 = cumsum(E_LH_s0);

E_q_PM0      = n_source0 .* abs(E_d2_qwf0).^5;
E_q_PM0_final = cumsum(E_q_PM0);
E_q_PM0_max   = max(E_q_PM0_final) + eps;

% --- 5) Plot (same panel order/style as Figure 10) ---
z2_mm = z2 / mm;
fig9b = figure('Name','Wavefront diagnostics (baseline t_shift=0)','Color','w');
tiledlayout(fig9b,4,2,'Padding','compact','TileSpacing','compact');

% (1) Intensity at wavefront
nexttile; plot(z2_mm, I_qwf0*1e-4, 'LineWidth',1.2); grid on;  % -> W/cm^2
xlabel('z (mm)'); ylabel('I_{qwf} (W/cm^2)');
title('1) Intensity at wavefront (t_{shift}=0)');

% (2) Plasma density
nexttile; plot(z2_mm, Ne_qwf0*1e-6, 'LineWidth',1.2); grid on; % -> cm^{-3}
xlabel('z (mm)'); ylabel('N_e (cm^{-3})');
title('2) Plasma density (t_{shift}=0)');

% (2.5) Plasma phase
nexttile; plot(z2_mm, Phi_plasma_qwf0/pi, 'LineWidth',1.2); grid on;
xlabel('z (mm)'); ylabel('\Delta\Phi_{plasma} (\pi)');
title('2.5) Plasma phase (t_{shift}=0)');

% (3) Dipole phase (long)
nexttile; plot(z2_mm, Phi_dip_long0_rel/pi, 'LineWidth',1.2); grid on;
xlabel('z (mm)'); ylabel('\Delta\Phi_{dip,long} (\pi)');
title('3) Dipole phase (long, t_{shift}=0)');

% (3.5) Total phase (long)
nexttile; plot(z2_mm, Phi_total_long_qwf0_rel/pi, 'LineWidth',1.2); grid on;
xlabel('z (mm)'); ylabel('\Delta\Phi_{total,long} (\pi)');
title('3.5) Total phase (long, t_{shift}=0)');

% (4) HHG source
nexttile; plot(z2_mm, n_source0, 'LineWidth',1.2); grid on;
xlabel('z (mm)'); ylabel('Source (arb./s)');
title('4) HHG source (t_{shift}=0)');

% (5) Ionization rate
nexttile; plot(z2_mm, W_qwf0, 'LineWidth',1.2); grid on;
xlabel('z (mm)'); ylabel('W (1/s)');
title('5) Ionization rate (t_{shift}=0)');

% (6) Accumulated HHG field (long)
nexttile; plot(z2_mm, abs(E_HHG_l0), 'LineWidth',1.2); grid on;
xlabel('z (mm)'); ylabel('|E_{HHG,long}| (arb. units)');
title('6) Accumulated HHG field (long, t_{shift}=0)');

% ---- sgtitle ----
if exist('q','var')
    sgtitle(sprintf('q-th HHG wavefront diagnostics (baseline t_{shift}=0, q=%d)', q));
else
    sgtitle('q-th HHG wavefront diagnostics (baseline t_{shift}=0)');
end
set(gcf,'WindowState','maximized');  % R2018a+ for regular figures
    drawnow;
%% 
%    % Export results (usual units)
%    result_Dispersion(1,:) = z2/mm;                      % position (mm)
%    result_Dispersion(2,:) = Delta_Phi_plasma/pi;       % plasma dispersion (rad/pi)
%    result_Dispersion(3,:) = Delta_Phi_capillary/pi;    % capillary dispersion (rad/pi)
%    result_Dispersion(4,:) = Phi_dip_long0_rel/pi;     % dipole phase - long (rad/pi)
%    result_Dispersion(5,:) = Phi_/pi;     % dipole phase - long (rad/pi)
%    result_Dispersion(6,:) = Delta_Phi_total_L_qwf/pi;      % total phase - long (rad/pi)
%    result_Dispersion(7,:) = Delta_Phi_total_S_qwf/pi;      % total phase - short (rad/pi))
%    result_Dispersion(8,:) = I_peak*1e-4;                   % peak intensity
%    result_Dispersion(9,:) = I_d2_qwf*1e-4;                 % qfw intensity
% 
%    fid_output = fopen('DispersionSource_qwf.txt','w');
%    fprintf(fid_output,'position   Delta_Phi_plasma  Delta_Phi_capillary  Delta_Phi_dipole_l   Delta_Phi_dipole_s         Delta_Phi_total_l        Delta_Phi_total_s    I_peak    I_d2_qwf     \r\n');
%    fprintf(fid_output,'(mm)       (rad/pi)          (rad/pi)             (rad/pi)             (rad/pi)                   (rad/pi)                 (rad/pi)             (W/cm^2)  (W/cm^2)     \r\n');
%    fprintf(fid_output,'%5.9e       %5.9e            %5.9e                %5.9e                %5.9e                      %5.9e                    %5.9e                %5.9e     %5.9e        \r\n',result_Dispersion);
%    fclose(fid_output);          
%       
 %%     
 % ==== Figure: q-th HHG wavefront diagnostics (q-wf) ====
fig9 = figure('Color','w'); 
tiledlayout(fig9,2,2,'Padding','compact','TileSpacing','compact');

% 1) Long-trajectory LH phase (relative) along q-wf
nexttile; 
plot(z2/mm, (Phi_LH_l - Phi_LH_l(1))/pi, 'LineWidth', 1.2); grid on;
xlabel('z (mm)');
ylabel('\Delta \Phi_{LH} (rad/\pi)');
title('Long-trajectory LH phase (q-wf)');

% 2) Short-trajectory LH phase (relative) along q-wf
nexttile; 
plot(z2/mm, (Phi_LH_s - Phi_LH_s(1))/pi, 'LineWidth', 1.2); grid on;
xlabel('z (mm)');
ylabel('\Delta \Phi_{LH} (rad/\pi)');
title('Short-trajectory LH phase (q-wf)');

% 3) Source density rate along q-wf
%    n_source is already built from q-wf populations * W_n_1_qwf (q-wf rate)
nexttile; 
plot(z2/mm, 1e-6*n_source, 'LineWidth', 1.2); hold on;
plot(z2/mm, W_n_1_qwf, 'LineWidth', 1.2); grid on; hold off;
xlabel('z (mm)');
ylabel('source density rate (1/s/cm^{-3})');
legend('n_{source} (×10^{-6})','W(|E_{qwf}|)','Location','best');
title('Source @ q-wf');

% 4) Accumulated HHG field magnitude (normalized) along q-wf
long_end  = abs(E_HHG_l(end))/E_q_PhaseMatched_final_max;
short_end = abs(E_HHG_s(end))/E_q_PhaseMatched_final_max;
nexttile; 
plot(z2/mm, abs(E_HHG_l)/E_q_PhaseMatched_final_max, 'LineWidth', 1.2); hold on;
plot(z2/mm, abs(E_HHG_s)/E_q_PhaseMatched_final_max, 'LineWidth', 1.2); grid on; hold off;
xlabel('z (mm)');
ylabel('|E_{HHG}| (arb. units)');
legend(sprintf('long, end=%.3g', long_end), sprintf('short, end=%.3g', short_end), ...
       'Location','southeast');
title('Accumulated HHG field (q-wf)');

% ---- Title (shows q and t_shift if available) ----
if exist('q','var') && exist('t_shift','var')
    sgtitle(sprintf('q-th HHG wavefront diagnostics (q = %d, t_{shift} = %.1f fs)', q, t_shift/1e-15));
elseif exist('q','var')
    sgtitle(sprintf('q-th HHG wavefront diagnostics (q = %d, no t_{shift})', q));
else
    sgtitle('q-th HHG wavefront diagnostics (q-wf)');
end
set(gcf,'WindowState','maximized');  % R2018a+ for regular figures
    drawnow;

% Export results (usual units)
   result_HHG(1,:) = z2/mm;           % position (mm)
   result_HHG(2,:) = E_q_PhaseMatched_final/E_q_PhaseMatched_final_max;
   result_HHG(3,:) = abs(E_HHG_l)/E_q_PhaseMatched_final_max;
   result_HHG(4,:) = abs(E_HHG_s)/E_q_PhaseMatched_final_max;              % ELH_WF_S
   result_HHG(5,:) = Phi_LH_l;
   result_HHG(6,:) = Phi_LH_s;
   result_HHG(7,:) = n_source; 
   fid_output = fopen(['HHG_output_q=',num2str(q),'.txt'],'w');
   fprintf(fid_output,'position     E_HHG_pm        E_HHG_l           E_HHG_s        Phi_LH_l        Phi_LH_s     n_source \r\n');
   fprintf(fid_output,'(mm)         (arb.units)     (arb.units)      (arb.units)    (rad)           (rad)         (1/s)    \r\n');
   fprintf(fid_output,'%6.9f        %6.9f           %6.9f            %6.9f         %6.9f           %6.9f          %6.9f    \r\n', result_HHG);
   fclose(fid_output);        

   q_HHG(qq) = q;
   HHG_final_l(qq) = abs(E_HHG_l(end))/E_q_PhaseMatched_final_max ;
   HHG_final_s(qq) = abs(E_HHG_s(end))/E_q_PhaseMatched_final_max ;
   
  
   
   HHG_gating(1,ii) = (ii-150)*0.1*fs;
   HHG_gating(2,ii) = abs(E_HHG_l(end))/E_q_PhaseMatched_final_max;
   HHG_gating(3,ii) = abs(E_HHG_s(end))/E_q_PhaseMatched_final_max;
   HHG_gating(4,ii) = (Phi_LH_l(end) - Phi_LH_l(1))/pi;
   HHG_gating(5,ii) = (Phi_LH_s(end) - Phi_LH_s(1))/pi;
   
%   close all
end
 

end
% Export results (usual units)
   result_HHG_gating(1,:) = HHG_gating(1,:) ;         
   result_HHG_gating(2,:) = HHG_gating(2,:) ;      
   result_HHG_gating(3,:) = HHG_gating(3,:) ;        
   result_HHG_gating(4,:) = HHG_gating(4,:) ;      
   result_HHG_gating(5,:) = HHG_gating(5,:) ;   
   
   fid_output = fopen(['HHG_gating.txt'],'w');
   fprintf(fid_output,'t     E_HHG_l        E_HHG_s       Phi_LH_l       Phi_LH_s      \r\n');
   fprintf(fid_output,'(fs)     (arb.units)    (arb.units)   (pi)          (pi)           \r\n');
   fprintf(fid_output,'%6.9f     %6.9f          %6.9f        %6.9f         %6.9f          \r\n', result_HHG_gating);
   fclose(fid_output);   
%-------------------------------------------------------------------------------------------

% Export results (usual units)
   result_HHG_final(1,:) = q_HHG ;         % position (mm)
   result_HHG_final(2,:) = HHG_final_l;           % position (mm)
   result_HHG_final(3,:) = HHG_final_s;           % position (mm)
   
   fid_output = fopen(['HHG_final.txt'],'w');
   fprintf(fid_output,'q_HHG     HHG_final_l        HHG_final_s           \r\n');
   fprintf(fid_output,' (arb.units)  (arb.units)  (arb.units)  \r\n');
   fprintf(fid_output,'%d     %6.9f           %6.9f        \r\n', result_HHG_final);
   fclose(fid_output);   
%%
% ================= Wavefront diagnostics vs z (using shifted observer) =================
z2_mm = z2 / mm;

% 1) Intensity at wavefront (harmonic-frame)
I_qwf = abs(E_t2_qwf).^2 / (2*mu_0*c);          % [W/m^2]

% 2) Plasma density along z AT THE WAVEFRONT (already t_shift-aligned above)
% Ne_z2 = N_e_cfit(z2)';                        % [1/m^3] (fit, not used below)
Ne_qwf = N_e_qwf;                                % [1/m^3] (from your wavefront sampling)

% --- Rebuild plasma refractive indices using the wavefront density ---
n_plasma_d_qwf = sqrt(1 - (q_e^2 .* Ne_qwf) ./ (epsilon_0 * m_e * omega_d^2));
n_plasma_q_qwf = sqrt(1 - (q_e^2 .* Ne_qwf) ./ (epsilon_0 * m_e * (q*omega_d)^2));

% Wavenumbers in plasma
k_d_plasma_qwf = k_0     .* n_plasma_d_qwf;     % [1/m]
k_q_plasma_qwf = (q*k_0) .* n_plasma_q_qwf;     % [1/m]

% Waveguide terms (reuse your definitions; broadcast if scalar)
if isscalar(k_d_capillary), k_d_cap_qwf = repmat(k_d_capillary, size(z2));
else,                       k_d_cap_qwf = k_d_capillary; end
if isscalar(k_q_capillary), k_q_cap_qwf = repmat(k_q_capillary, size(z2));
else,                       k_q_cap_qwf = k_q_capillary; end

k_d_total_qwf = k_d_plasma_qwf + k_d_cap_qwf;   % [1/m]
k_q_total_qwf = k_q_plasma_qwf + k_q_cap_qwf;   % [1/m]

% Plasma-only mismatch and accumulated plasma phase at the wavefront
Delta_k_plasma_qwf = q*k_d_plasma_qwf - k_q_plasma_qwf;     % [1/m]
Phi_plasma_qwf     = cumsum(Delta_k_plasma_qwf) * dz2;       % [rad]

% ================= Dipole contribution on q-wf: Δk_dip = α(I) * dI/dz =================
% α_l(I), α_s(I) evaluated at wavefront intensity
alpha_l_qwf = alpha_l_cfit(I_qwf)';   % [m^2/W]
alpha_s_qwf = alpha_s_cfit(I_qwf)';   % [m^2/W]

% dI/dz on the q-wavefront (central difference; ends use one-sided)
dI_dz_qwf = zeros(size(I_qwf));
if isscalar(dz2)
    if numel(I_qwf) >= 3
        dI_dz_qwf(1)        = (I_qwf(2) - I_qwf(1)) / dz2;
        dI_dz_qwf(2:end-1)  = (I_qwf(3:end) - I_qwf(1:end-2)) / (2*dz2);
        dI_dz_qwf(end)      = (I_qwf(end) - I_qwf(end-1)) / dz2;
    else
        dI_dz_qwf(:) = 0;  % degenerate grid
    end
else
    % nonuniform grid support
    dI_dz_qwf = gradient(I_qwf(:), z2(:)); 
    dI_dz_qwf = dI_dz_qwf(:).';
end

Delta_k_dipole_l_qwf = alpha_l_qwf .* dI_dz_qwf';   % [1/m]
Delta_k_dipole_s_qwf = alpha_s_qwf .* dI_dz_qwf';   % [1/m]

% ================= Total mismatch at q-wf (t_shift observer) =================
% Capillary mismatch term (geometry-only) reused
Delta_k_capillary_qwf = q*k_d_cap_qwf - k_q_cap_qwf;         % [1/m]

Delta_k_total_l_qwf = Delta_k_plasma_qwf + Delta_k_capillary_qwf + Delta_k_dipole_l_qwf;  % [1/m]
Delta_k_total_s_qwf = Delta_k_plasma_qwf + Delta_k_capillary_qwf + Delta_k_dipole_s_qwf;  % [1/m]

% Accumulated total phase on q-wf (optional diagnostics)
Delta_Phi_total_l_qwf = cumsum(Delta_k_total_l_qwf) * dz2;   % [rad]
Delta_Phi_total_s_qwf = cumsum(Delta_k_total_s_qwf) * dz2;   % [rad]

% Dipole phase (long) evaluated at the wavefront intensity
Phi_dip_long     = Phi_dipole_l_cfit(I_qwf)';                 % [rad]
Phi_dip_long_rel = Phi_dip_long - Phi_dip_long(1);

% Total phase (long) using q-wf plasma phase + same capillary + dipole
Phi_total_long_qwf     = Phi_plasma_qwf + Delta_Phi_capillary + Phi_dip_long;  % keep your plotting ref
Phi_total_long_qwf_rel = Phi_total_long_qwf - Phi_total_long_qwf(1);

% 4) HHG source and 5) Ionization rate at the wavefront (already qwf-based)
Source_HHG = n_source;               % [arb./s]
W_qwf      = W_n_1_qwf;              % [1/s]

% 6) Accumulated HHG (long) magnitude (normalized)
Eacc_long      = abs(E_HHG_l);
Eacc_long_norm = Eacc_long / (max(Eacc_long) + eps);

% (Optional) coherence lengths at q-wf
L_dephasing_l_qwf = pi ./ max(abs(Delta_k_total_l_qwf), realmin);  % [m]
L_dephasing_s_qwf = pi ./ max(abs(Delta_k_total_s_qwf), realmin);  % [m]

%%
% ---------- Plots ----------
figure('Name','Wavefront variation vs z','Color','w');
tiledlayout(4,3, 'Padding','compact','TileSpacing','compact');  % <- was (4,2)

% 1. Intensity (q-wf)
nexttile; plot(z2_mm, I_qwf*1e-4, 'LineWidth', 1.2); % W/cm^2
xlabel('z (mm)'); ylabel('I_{qwf} (W/cm^2)'); title('1) Intensity at wavefront');

% 2. Plasma density (q-wf)
nexttile; plot(z2_mm, Ne_qwf*1e-6, 'LineWidth', 1.2); % cm^{-3}
xlabel('z (mm)'); ylabel('N_e (cm^{-3})'); title('2) Plasma density');

% 2.5 Plasma phase (q-wf)
nexttile; plot(z2_mm, Phi_plasma_qwf/pi, 'LineWidth', 1.2);
xlabel('z (mm)'); ylabel('\Delta\Phi_{plasma} (\pi)'); title('2.5) Plasma phase (qwf)');

% 3. Dipole phase (long, q-wf)
nexttile; plot(z2_mm, Phi_dip_long_rel/pi, 'LineWidth', 1.2);
xlabel('z (mm)'); ylabel('\Delta\Phi_{dip,long} (\pi)'); title('3) Dipole phase (long)');

% 3.5 Total phase (long, q-wf)
nexttile; plot(z2_mm, Phi_total_long_qwf_rel/pi, 'LineWidth', 1.2);
xlabel('z (mm)'); ylabel('\Delta\Phi_{total,long} (\pi)'); title('3.5) Total phase (long)');

% 4. HHG source (q-wf)
nexttile; plot(z2_mm, Source_HHG, 'LineWidth', 1.2);
xlabel('z (mm)'); ylabel('Source (arb./s)'); title('4) HHG source');

% 5. Ionization rate (q-wf)
nexttile; plot(z2_mm, W_qwf, 'LineWidth', 1.2);
xlabel('z (mm)'); ylabel('W (1/s)'); title('5) Ionization rate');

% 6. Accumulated HHG field (long)
nexttile; plot(z2_mm, Eacc_long, '-o', 'LineWidth', 1.2);
xlabel('z (mm)'); ylabel('|E_{HHG,long}| (arb. units)'); title('6) Accumulated HHG field');

% 7. Total wavenumber mismatch (q-wf) — long & short on the same axes
nexttile;
plot(z2_mm, pi ./ Delta_k_total_l_qwf .*1E3, 'LineWidth', 1.2); hold on;
plot(z2_mm, pi ./ Delta_k_total_s_qwf .*1E3, 'LineWidth', 1.2);
yline(0,'k:'); hold off; grid on;
xlabel('z (mm)'); ylabel('Dephasing length (mm)');
legend('long','short','0','Location','best');
title('7) Dephasing length (q-wf)');

% leave the last tile empty (or use it for something else)
nexttile; axis off;

sgtitle(sprintf('Wavefront diagnostics (t\\_shift = %.1f fs, q = %d)', t_shift/fs, q));
set(gcf,'WindowState','maximized');  % R2018a+ for regular figures
drawnow;


%%
% outDir = save_all_figs_with_sgtitle('D:\Google_BackUp\Capillary HHG');
% Lib=======================================================================================
%%
function outDir = save_all_figs_with_sgtitle(outRoot)
% save_all_figs_with_sgtitle  Save all open figures using sgtitle/tiledlayout title as the filename.
%
%   outDir = save_all_figs_with_sgtitle(outRoot)
%
%   - Creates a timestamped subfolder under outRoot.
%   - For each open figure:
%       * Tries to read the tiledlayout Title (if present)
%       * Else tries sgtitle (suptitle) text
%       * Else falls back to Figure Name
%       * Else 'Figure_<N>'
%     Then saves: <title>.png (300 dpi) and <title>.fig
%
%   Example:
%       outDir = save_all_figs_with_sgtitle('D:\HHG_runs');

    if nargin < 1 || isempty(outRoot)
        error('Please provide an output root folder, e.g., D:\HHG_runs');
    end
    if ~exist(outRoot, 'dir')
        mkdir(outRoot);
    end

    % Timestamped folder
    ts = datestr(now, 'yyyymmdd_HHMMSS');
    outDir = fullfile(outRoot, ts);
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    % Get figures in creation order
    figs = findall(0, 'Type', 'figure');
    if isempty(figs)
        warning('No open figures to save.');
        return;
    end
    figs = flipud(figs); % so Figure 1..N order

    for k = 1:numel(figs)
        fig = figs(k);

        % -------- Get a global title string --------
        titleStr = '';

        % 1) tiledlayout Title (preferred when using tiledlayout)
        tl = findall(fig, 'Type', 'tiledlayout');
        if ~isempty(tl)
            try
                tStr = tl(1).Title.String;
                if iscell(tStr), tStr = strjoin(tStr, ' - '); end
                if ~isempty(strtrim(tStr))
                    titleStr = tStr;
                end
            end
        end

        % 2) sgtitle / suptitle (if used)
        if isempty(titleStr)
            % sgtitle usually creates a special axes tagged 'suptitle'
            supAx = findall(fig, 'Type', 'axes', 'Tag', 'suptitle');
            if ~isempty(supAx)
                % Text object is its child
                txt = findall(supAx(1), 'Type', 'text');
                if ~isempty(txt)
                    tStr = txt(1).String;
                    if iscell(tStr), tStr = strjoin(tStr, ' - '); end
                    if ~isempty(strtrim(tStr))
                        titleStr = tStr;
                    end
                end
            else
                % Some MATLAB versions tag the text directly
                supTxt = findall(fig, 'Type', 'text', 'Tag', 'suptitle');
                if ~isempty(supTxt)
                    tStr = supTxt(1).String;
                    if iscell(tStr), tStr = strjoin(tStr, ' - '); end
                    if ~isempty(strtrim(tStr))
                        titleStr = tStr;
                    end
                end
            end
        end

        % 3) Fallback: Figure Name
        if isempty(titleStr)
            fName = get(fig, 'Name');
            if ~isempty(strtrim(fName))
                titleStr = fName;
            end
        end

        % 4) Final fallback
        if isempty(titleStr)
            titleStr = sprintf('Figure_%d', k);
        end

        % -------- Sanitize to a safe filename --------
        safeName = titleStr;
        % Replace newline/whitespace runs with single space
        safeName = regexprep(safeName, '\s+', ' ');
        % Remove illegal filename characters
        safeName = regexprep(safeName, '[/\\:*?"<>|]', '');
        % Trim
        safeName = strtrim(safeName);
        % Collapse spaces to underscores
        safeName = regexprep(safeName, '\s', '_');
        % Limit length if needed (Windows path limits)
        if numel(safeName) > 180
            safeName = safeName(1:180);
        end

        % -------- Save --------
        pngPath = fullfile(outDir, [safeName '.png']);
        figPath = fullfile(outDir, [safeName '.fig']);

        try
            % Use exportgraphics for high-quality output (R2020a+)
            exportgraphics(fig, pngPath, 'Resolution', 300);
        catch
            % Fallback
            saveas(fig, pngPath);
        end

        savefig(fig, figPath);

        fprintf('Saved: %s\n', pngPath);
    end

    fprintf('\nAll figures saved to: %s\n', outDir);
end



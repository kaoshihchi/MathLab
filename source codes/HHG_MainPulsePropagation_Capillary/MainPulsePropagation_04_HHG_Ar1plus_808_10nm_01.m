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
clear all

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
     
%Load the parameter file
%    fid = fopen(ParameterFilename,'r');
%    while feof(fid) == 0
%       line = fgetl(fid);
%       eval(line);
%    end
%    fclose(fid);
   
%%   
% % %%
    FigureSwitch = 1;

% Laser parameters
   lambda = 808 * nm;         % laser wavelength (m)
   tau_0  = 35 * fs;          % initial pulse duration (1/e intensity) (unit: sec)
   PulseEnergy = 16 * 10^-3;  % pulse energy (unit: J)

% Capillary parameters
   R_capillary = 75 * um;     % capillary radius (unit: m)
   L_capillary = 10  * mm;    % capillary length (unit: m)
   n_capillary = 1.4531;      % refractive index of the capillary
                              % fused silica at 810 nm: 1.4531
                              % bk7 at 810 nm: 1.51
   a_capillary = 28;        % guiding loss coefficient (unit: 1/m)

% gas target parameters
   %n_gas = 2.37 * 10^16 / cm^3;  % gas density (unit: m^-3) 
   %n_gas = 5.87 * 10^16 / cm^3;  % gas density (unit: m^-3) 
   n_gas = 5 * 10^16 / cm^3;  % gas density (unit: m^-3) % 9.33 * 10^15  91-0.683

% Calculation parameters
   dz = 0.05 * mm;        % Spatial resolution of the ionization calculation (m)

%---------------------------------------------------------------------------

% Parameters for ATI heating calculation
   time_initial    = -100 * fs;      % initial time (sec)
   time_final      =  100 * fs;      % final time (sec)
   DeltaTime       = 0.005 * fs;      % time interval (sec)

% Ar ionization potential
   E_ion_0    = 15.7596 * eV;          % ionization potential (J) ( 0  -> 1+ )
   E_ion_1    = 27.6297 * eV;          % ionization potential (J) ( 1+ -> 2+ )
   E_ion_2    = 40.74   * eV;          % ionization potential (J) ( 2+ -> 3+ )
   E_ion_3    = 59.81   * eV;          % ionization potential (J) ( 3+ -> 4+ )
   E_ion_4    = 75.02   * eV;          % ionization potential (J) ( 4+ -> 5+ )
%---------------------------------------------------------------------------

% Parameters for HHG
   I_p = E_ion_1;   % ionizatiob potential of 
   q   = 81 + 0*(2*qq - 50);           % harmonic order
   dz2 = dz; %0.001 * mm;     % Spatial resolution of the HHG calculation (m)

%---------------------------------------------------------------------------


% Note: ionization potential
% He ionization potential
   % E_ion_0    = 24.587 * eV;          % ionization potential (J) ( 0  -> 1+ )
   % E_ion_1    = 54.418 * eV;          % ionization potential (J) ( 1+ -> 2+ )
   % E_ion_2    = 1000   * eV;          % ionization potential (J) ( 2+ -> 3+ )
   % E_ion_3    = 1200   * eV;          % ionization potential (J) ( 3+ -> 4+ )
   % E_ion_4    = 2000   * eV;          % ionization potential (J) ( 4+ -> 5+ )

% Ar ionization potential
   % E_ion_0    = 15.7596 * eV;          % ionization potential (J) ( 0  -> 1+ )
   % E_ion_1    = 27.6297 * eV;          % ionization potential (J) ( 1+ -> 2+ )
   % E_ion_2    = 40.74   * eV;          % ionization potential (J) ( 2+ -> 3+ )
   % E_ion_3    = 59.81   * eV;          % ionization potential (J) ( 3+ -> 4+ )
   % E_ion_4    = 75.02   * eV;          % ionization potential (J) ( 4+ -> 5+ )

% Ne ionization potential
   % E_ion_0    = 21.5645 * eV;          % ionization potential (J) ( 0  -> 1+ )
   % E_ion_1    = 40.9630 * eV;          % ionization potential (J) ( 1+ -> 2+ )
   % E_ion_2    = 63.45   * eV;          % ionization potential (J) ( 2+ -> 3+ )
   % E_ion_3    = 97.12   * eV;          % ionization potential (J) ( 3+ -> 4+ )
   % E_ion_4    = 126.21  * eV;          % ionization potential (J) ( 4+ -> 5+ )
%-------------------------------------------------------------------------------------------

% Laser parameters
   omega_d = 2 * pi * c / lambda;    % laser angular frequency (rad/sec)
   k_0     = 2 * pi / lambda;        % laser wavenumber in vacuum (1/m)
   I_peak_ini = PulseEnergy / (pi^1.5 * R_capillary^2 * tau_0);
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
   % a_capillary = CapillaryAttenuationCoefficient(lambda,n_capillary,R_capillary);
   
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
      I_peak(j+1) = LaserEnergy(j+1)/(pi^1.5 * R_capillary^2 * tau(j+1));
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
   
% Export results (SI units)
   result_SI_temp(1,:) = z;            % position (m)
   result_SI_temp(2,:) = LaserEnergy;  % laser pulse energy (J)
   result_SI_temp(3,:) = Z_ion;        % average ionization state
   result_SI_temp(4,:) = N_e;          % electron density N_e(z) (1/m^3)
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
   result_temp(4,:) = N_e*cm^3;       % electron density N_e(z) (1/cm^3)
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
   
% Plot results
   if FigureSwitch
   figure;
   subplot(4,3,1),  plot(z/mm,LaserEnergy/mJ),
                    ylabel('laser energy (mJ)');
   subplot(4,3,4),  plot(z/mm,E_peak),
                    ylabel('E_{peak} (V/m)');
   subplot(4,3,7),  plot(z/mm,I_peak),
                    ylabel('I_{peak} (W/m^2)');
   subplot(4,3,10), plot(z/mm,N_e/cm^(-3)),
                    ylabel('N_e (cm^{-3})');
                    xlabel('position z (mm)');
                    legend('electron density','Location','NorthEast');
   subplot(4,3,2),  plot(z/mm,GVD*mm/fs^2),
                    ylabel('GVD (fs^2/mm)');
                    hold on, plot(0,GVD_capillary*mm/fs^2,'r*'), hold off;
                    legend('plasma','waveguide','Location','NorthEast');
   subplot(4,3,5),  plot(z/mm,D/fs^2),
                    ylabel('accumulated GDD (fs^2)');
   subplot(4,3,8),  plot(z/mm,tau/fs),
                    ylabel('pulse duration (fs)');
   subplot(4,3,11), plot(z/mm,Energy_Ionization/dz),
                    ylabel('ionization loss (mJ/mm)');
                    xlabel('position z (mm)');
   subplot(4,3,3),  plot(z/mm,Energy_ATI/dz),
                    ylabel('ATI loss (mJ/mm)'),
   subplot(4,3,6),  plot(z/mm,Energy_IB/dz),
                    ylabel('IB loss (mJ/mm)'),
   subplot(4,3,9),  plot(z/mm,Energy_TS/dz),
                    ylabel('TS loss (mJ/mm)'),
   subplot(4,3,12), plot(z/mm,Energy_capillary/dz),
                    ylabel('capillary loss (mJ/mm)'),
                    xlabel('position z (mm)');
   sgtitle('Driving pulse propagation');
   end
                    

% Curve fitting
   % peak electric field (cfit object) E_peak_cfit(z)
   % This is a cfit object, which can be used as a function:
   %    E_peak = E_peak_cfit(z),
   E_peak_cfit = fit(z',E_peak','exp1');   % unit: V/m 
   
   % peak laser intensity (cfit object) I_peak_cfit
   % This is a cfit object, which can be used as a function:
   %    I_peak = I_peak_cfit(z),
   I_peak_cfit = fit(z',I_peak','exp1');   % unit: W/m^2
   
   % pulse duration (cfit object) tau_cfit
   tau_cfit = fit(z',(tau/fs)','poly3');   % unit: fs
   tau_fun = @(z) tau_cfit(z)*fs;         % unit: sec 

   % electron density (cfit object) N_e_cfit(z)
   N_e_cfit = fit(z',N_e','poly3');        % unit: m^-3
   
   % accumulated GDD (cfit object) D_cfit(z)
   D_cfit = fit(z',(D/fs^2)','poly3');     % unit: fs^2
   D_fun = @(z) D_cfit(z)*fs^2;            % unit: s^2
   
   % plasma refractive index of the driving pulse (cfit object) n_plasma_cfit(z)
   n_plasma_d_cfit = fit(z',n_plasma_d','poly3');
   
   % accumulated group delay (cfit object) C_cfit(z)
   C_cfit = fit(z',(C/fs)','poly3');       % unit: fs
   C_fun = @(z) C_cfit(z)*fs;              % unit: sec
   
   % total wavenumber of the driving pulse k_d_total_cfit(z)
   k_d_total_cfit = fit(z',k_d_total','poly3');    % unit: 1/m
   
   % accumulated phase shift of the driving pulse due to propagation (rad)
   phi_d_prop_cfit = fit(z',phi_d_prop','poly3');    % unit: rad

   % phase velocity of the driving pulse v_d_cfit(z)
   v_d_cfit = fit(z',v_d','poly3');    % unit: m/sec

   if FigureSwitch
      figure;
      subplot(4,3,1),  plot(E_peak_cfit,z,E_peak),
                       ylabel('E_{peak} (V/m)');
                       xlabel('position z (m)');
                       title('peak electric field');
      subplot(4,3,4),  plot(I_peak_cfit,z,I_peak),
                       ylabel('I_{peak} (W/m^2)');
                       xlabel('position z (m)');
                       title('peak intensity');
      subplot(4,3,7),  plot(N_e_cfit,z,N_e),
                       ylabel('N_e (m^{-3})');
                       xlabel('position z (m)');
                       title('electron density');
      subplot(4,3,2),  plot(n_plasma_d_cfit,z,n_plasma_d),
                       ylabel('n_{plasma}');
                       xlabel('position z (m)');
                       title('plasma refractive index');
      subplot(4,3,5),  plot(D_cfit,z,D/fs^2),
                       ylabel('D (fs^2)');
                       xlabel('position z (m)');
                       title('group-delay-dispersion (GDD)');
      subplot(4,3,8),  plot(tau_cfit,z,tau/fs),
                       ylabel('\tau (fs)');
                       xlabel('position z (m)');
                       title('pulse duration');
      subplot(4,3,11),  plot(C_cfit,z,C/fs),
                       ylabel('C (fs)');
                       xlabel('position z (m)');
                       title('group delay');
      subplot(4,3,3),  plot(k_d_total_cfit,z,k_d_total),
                       ylabel('k_{d\_total}(z) (1/m)');
                       xlabel('position z (m)');
                       title('total wavenumber');
      subplot(4,3,6),  plot(phi_d_prop_cfit,z,phi_d_prop),
                       ylabel('\phi_{d\_prop}(z) (rad)');
                       xlabel('position z (m)');
                       title('accumulated phase due to propagation');
      subplot(4,3,9),  plot(v_d_cfit,z,v_d),
                       ylabel('v_d (m/sec)');
                       xlabel('position z (m)');
                       title('phase velocity');
      sgtitle('Driving pulse propagation');
   end
   
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
   Phi_dipole_l_cfit = fit(I_dipole',Phi_dipole_l','fourier3');  % unit: rad
   Phi_dipole_s_cfit = fit(I_dipole',Phi_dipole_s','power2');    % unit: rad
   alpha_l_cfit = fit(I_dipole',alpha_l','poly5');       % unit: m^2/W
   alpha_s_cfit = fit(I_dipole',alpha_s','poly5');       % unit: m^2/W
   
   if FigureSwitch
   figure;
   subplot(1,4,1), plot(Phi_dipole_l_cfit,I_dipole,Phi_dipole_l);
      xlabel('intensity (W/m^2)'), ylabel('\Phi_{dipole\_long} (rad)');
      title('long-trajectory dipole phase');
   subplot(1,4,2), plot(Phi_dipole_s_cfit,I_dipole,Phi_dipole_s);
      xlabel('intensity (W/m^2)'), ylabel('\Phi_{dipole\_short} (rad)');
      title('short-trajectory dipole phase');
   subplot(1,4,3), plot(alpha_l_cfit,I_dipole,alpha_l);
      xlabel('intensity (W/cm^2)'), ylabel('\alpha_{long} (m^2/W)');
      title('long-trajectory \alpha');
   subplot(1,4,4), plot(alpha_s_cfit,I_dipole,alpha_s);
      xlabel('intensity (W/cm^2)'), ylabel('\alpha_{short} (m^2/W)');
      title('short-trajectory \alpha');
   sgtitle('Dipole phase calculation');
   end
   
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

   figure;
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
   subplot(5,3,2), plot(z2/mm,Delta_Phi_plasma,z2/mm,Delta_Phi_capillary,'r*');
      xlabel('z (mm)');
      ylabel('\Delta \Phi (rad)');
      legend('plasma','waveguide','Location','NorthEast');
      title('accumulated phase mismatch');
   subplot(5,3,5), plot(z2/mm,Delta_Phi_dipole_l);
      xlabel('z (mm)');
      ylabel('\Delta \Phi_{dipole\_long} (rad)');
      legend('dipole\_long');
   subplot(5,3,8), plot(z2/mm,Delta_Phi_dipole_s);
      xlabel('z (mm)');
      ylabel('\Delta \Phi_{dipole\_short} (rad)');
      legend('dipole\_short');
   subplot(5,3,11), plot(z2/mm,Delta_Phi_total_l);
      xlabel('z (mm)');
      ylabel('\Delta \Phi_{total\_long} (rad)');
      legend('total\_long');
   subplot(5,3,14), plot(z2/mm,Delta_Phi_total_s);
      xlabel('z (mm)');
      ylabel('\Delta \Phi_{total\_short} (rad)');
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
      % Generate the driving field E_d2 at (z2,t=t_d2(z2))
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

   
% phase of the local harmonic field
% Trace a fixed harmonic wavefront
   % long-trajectory emission
   Phi_LH_l = q*Phi_d2_qwf + Phi_dipole_l_cfit(I_d2_qwf)';
   % short-trajectory emission
   Phi_LH_s = q*Phi_d2_qwf + Phi_dipole_s_cfit(I_d2_qwf)';
   
% Calculate the local harmonic field
   % long-trajectory emission
   %Phi_LH_l = Phi_LH_l - Phi_LH_l(1);
   E_t2 = zeros(N,size(time,2));
   ionizationResult = zeros(8,size(time,2));
   n_0 = zeros(N,size(time,2));
   n_1 = zeros(N,size(time,2));
   n_2 = zeros(N,size(time,2));
   n_3 = zeros(N,size(time,2));  
   n_e = zeros(N,size(time,2));   
   delta_tw = zeros(N,1);
   n_1_qwf = zeros(N,1);    
   n_2_qwf = zeros(N,1);   
   n_3_qwf = zeros(N,1);    
   E_t2_qwf = zeros(N,1);
   for jj = 1:N
        E_t2(jj,:) = ElectricField_d(time+C_fun(z2(jj))-C_fun(z2(1)),E_peak(jj),omega_d,tau_0,D(jj)-D(1),C_fun(z2(jj))-C_fun(z2(1)),0); 
        % E_t      = ElectricField_d(time,                           E_peak(1) ,omega_d,tau_0,0         ,0                         ,0);
        ionizationResult = TunnelingIonizationRate_Linear2((E_t2(jj,:)),omega_d,time+C_fun(z2(jj)),E_ion_0,E_ion_1,E_ion_2,E_ion_3,E_ion_4,FigureSwitch);
        ionizationResult = ionizationResult';
        n_e(jj,:) = ionizationResult(1,:);        
        n_0(jj,:) = ionizationResult(2,:);
        n_1(jj,:) = ionizationResult(3,:);
        n_2(jj,:) = ionizationResult(4,:);
        n_3(jj,:) = ionizationResult(5,:);        
        delta_tw(jj) = C_fun(z2(jj))-C_fun(z2(1)) - t2_qwf(jj);
        n_1_qwf(jj) = n_1(jj, round(0.5*size(time,2))- 1*round(delta_tw(jj)/DeltaTime));
        n_2_qwf(jj) = n_2(jj, round(0.5*size(time,2))- 1*round(delta_tw(jj)/DeltaTime));
        n_3_qwf(jj) = n_3(jj, round(0.5*size(time,2))- 1*round(delta_tw(jj)/DeltaTime));
        E_t2_qwf(jj) =  E_t2(jj, round(0.5*size(time,2)) - 1*round(delta_tw(jj)/DeltaTime));
   end
  
 %%
  % 1+ ~ 2+
    W_n_1_qwf = StaticIonizationRate(E_ion_1, abs(E_d2_qwf));
    n_source = n_gas.*n_1_qwf'.*W_n_1_qwf;
    E_LH_l = n_source.*abs(E_d2_qwf).^5 .* exp(1i*Phi_LH_l);    % arb. units
    % short-trajectory emission
    E_LH_s = n_source.*abs(E_d2_qwf).^5 .* exp(1i*Phi_LH_s);    % arb. units
   
%    % 2+ ~ 3+
%    W_n_1_qwf = StaticIonizationRate(E_ion_2, abs(E_d2_qwf));
%    n_source = n_gas.*n_2_qwf'.*W_n_1_qwf;
%    E_LH_l = n_source.*abs(E_d2_qwf).^5 .* exp(1i*Phi_LH_l);    % arb. units
%    % short-trajectory emission
%    E_LH_s = n_source.*abs(E_d2_qwf).^5 .* exp(1i*Phi_LH_s);    % arb. units
   
% Calculate the accumulated harmonic field
   E_HHG_l = cumsum(E_LH_l);
   E_HHG_s = cumsum(E_LH_s);
   E_q_PhaseMatched = n_source.*abs(E_d2_qwf).^5;% phase matched HHG field (as a function of z, column vector)(Z_ion - 1)
   E_q_PhaseMatched_final = cumsum(E_q_PhaseMatched);
   E_q_PhaseMatched_final_max = max(E_q_PhaseMatched_final); %max(E_q_PhaseMatched_final);
%%   
   figure; % driving field check
   subplot(3,1,1), plot(z2/mm,abs(E_d2_dwf),z2/mm,abs(E_d2_qwf),z/mm,E_peak);
      xlabel('z (mm)');
      ylabel('|E| (V/m)');
      legend('fixed d-wf','fixed q-wf','E_{peak}','Location','SouthWest');
      title('driving field amplitude')
   subplot(3,1,2), plot(z2/mm,(t2_dwf-(C_fun(z2)-C_fun(0))')/fs,z2/mm,(t2_qwf-(C_fun(z2)-C_fun(0))')/fs);
      xlabel('z (mm)');
      ylabel('\Delta t (fs)');
      legend('t2_{dwf} - C(z) (fixed d-wf)','t2_{qwf} - C(z) (fixed q-wf)','Location','SouthWest');
      title('time difference')
   subplot(3,1,3), plot(z2/mm,Phi_d2_dwf,z2/mm,Phi_d2_qwf);
      xlabel('z (mm)');
      ylabel('\Phi_{d2} (rad)');
      legend('fixed d-wf','fixed q-wf','Location','SouthWest');
      title('phase of the driving pulse')
 %%  
   figure
   subplot(2,2,1)
   plot(z2/mm,I_d2_qwf,z2/mm,I_peak);
      xlabel('z (mm)');
      ylabel('I (W/m^2)');
      legend('fixed q-wf','peak','Location','NorthWest');
      title('fixed q-wf driving field intensity')
   subplot(2,2,2)
   
   Delta_Phi_dipole_Peak = Phi_dipole_l_cfit(I_peak)-(Phi_dipole_l_cfit(I_peak(1)));  
   Delta_Phi_dipole_L_qwf = Phi_dipole_l_cfit(I_d2_qwf)-(Phi_dipole_l_cfit(I_d2_qwf(1)));
   Delta_Phi_dipole_S_qwf = Phi_dipole_s_cfit(I_d2_qwf)-(Phi_dipole_s_cfit(I_d2_qwf(1)));   

   plot(z2/mm,Delta_Phi_dipole_L_qwf/pi,z2/mm,Delta_Phi_dipole_Peak/pi);
      xlabel('z (mm)');
      ylabel('\Delta\Phi(rad/\pi)');
      legend('fixed q-wf','peak','Location','NorthWest');
      title('fixed q-wf driving field dipole')
   subplot(2,2,3)
   plot(z2/mm,Delta_Phi_plasma/pi);
      xlabel('z (mm)');
      ylabel('\Delta\Phi_p(rad/\pi)');
      title('plasma phase')
   subplot(2,2,4)
   Delta_Phi_total_L_qwf = Delta_Phi_plasma + Delta_Phi_capillary + Delta_Phi_dipole_L_qwf' ;
   Delta_Phi_total_S_qwf = Delta_Phi_plasma + Delta_Phi_capillary + Delta_Phi_dipole_S_qwf' ;   
   plot(z2/mm, (Delta_Phi_total_L_qwf - Delta_Phi_total_L_qwf(1))/pi );
      xlabel('z (mm)');
      ylabel('\Delta\Phi_{total}(rad/\pi)');
      title('total phase') 
 
   % Export results (usual units)
   result_Dispersion(1,:) = z2/mm;                      % position (mm)
   result_Dispersion(2,:) = Delta_Phi_plasma/pi;       % plasma dispersion (rad/pi)
   result_Dispersion(3,:) = Delta_Phi_capillary/pi;    % capillary dispersion (rad/pi)
   result_Dispersion(4,:) = Delta_Phi_dipole_L_qwf/pi;     % dipole phase - long (rad/pi)
   result_Dispersion(5,:) = Delta_Phi_dipole_S_qwf/pi;     % dipole phase - long (rad/pi)
   result_Dispersion(6,:) = Delta_Phi_total_L_qwf/pi;      % total phase - long (rad/pi)
   result_Dispersion(7,:) = Delta_Phi_total_S_qwf/pi;      % total phase - short (rad/pi))
   result_Dispersion(8,:) = I_peak*1e-4;                   % peak intensity
   result_Dispersion(9,:) = I_d2_qwf*1e-4;                 % qfw intensity

   fid_output = fopen('DispersionSource_qwf.txt','w');
   fprintf(fid_output,'position   Delta_Phi_plasma  Delta_Phi_capillary  Delta_Phi_dipole_l   Delta_Phi_dipole_s         Delta_Phi_total_l        Delta_Phi_total_s    I_peak    I_d2_qwf     \r\n');
   fprintf(fid_output,'(mm)       (rad/pi)          (rad/pi)             (rad/pi)             (rad/pi)                   (rad/pi)                 (rad/pi)             (W/cm^2)  (W/cm^2)     \r\n');
   fprintf(fid_output,'%5.9e       %5.9e            %5.9e                %5.9e                %5.9e                      %5.9e                    %5.9e                %5.9e     %5.9e        \r\n',result_Dispersion);
   fclose(fid_output);          
      
 %%     
   figure    
   subplot(2,2,1), plot(z2/mm,(Phi_LH_l-Phi_LH_l(1))/pi);grid on
      xlabel('z (mm)');
      ylabel('\Delta \Phi_{LH} (rad/\pi)');
      title('phase of the long-trajectory LH field (fixed q-wf)');
   subplot(2,2,2), plot(z2/mm,(Phi_LH_s-Phi_LH_s(1))/pi);grid on
      xlabel('z (mm)');
      ylabel('\Delta \Phi_{LH} (rad/\pi)');
      title('phase of the short-trajectory LH field (fixed q-wf)');
   subplot(2,2,3), plot(z2/mm,1e-6*n_source,z2/mm,W_n_1_qwf);grid on
      xlabel('z (mm)');
      ylabel('source density rate (1/s/cm^{-3})');
      title('He^{1+} \times W(|E_{qwf}|)');           
   subplot(2,2,4), plot(z2/mm,abs(E_HHG_l)/E_q_PhaseMatched_final_max,z2/mm,abs(E_HHG_s)/E_q_PhaseMatched_final_max);grid on
   % subplot(2,2,4), plot(z2/mm,abs(E_HHG_l)/E_q_PhaseMatched_final,z2/mm,abs(E_HHG_s)/E_q_PhaseMatched_final);grid on
      xlabel('z (mm)');
      ylabel('|E_{HHG}| (arb. units)');
      legend(['long, end=',num2str(abs(E_HHG_l(end))/E_q_PhaseMatched_final_max)],['short, end=',num2str(abs(E_HHG_s(end))/E_q_PhaseMatched_final_max)],'Location','southeast');
      % legend(['long, end=',num2str(abs(E_HHG_l(end))/E_q_PhaseMatched_final_max)],['short, end=',num2str(abs(E_HHG_s(end))/E_q_PhaseMatched_final_max)],'Location','southeast');
      title('accumulated harmonic field (fixed q-wf)');

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
%-------------------------------------------------------------------------------------------


%-------------------------------------------------------------------------
% TunnelingIonizationRate_Linear:
%     evaluation of the optical-field ionization rate, electron density,
%     energy distribution, averaged energy, and etc.
%
%     1. Assume linear-polarized Gaussian pulse.
%     2. MKS unit
%-------------------------------------------------------------------------
function result = TunnelingIonizationRate_Linear(E_t,omega_d,time,...
                      E_ion_0,E_ion_1,E_ion_2,E_ion_3,E_ion_4,FigureSwitch)

 % Unit conversion
   global eV   
   global fs   

 % Physical constants
   global c           %#ok<NUSED>
   global q_e
   global m_e
   global mu_0        %#ok<NUSED>
   global epsilon_0   %#ok<NUSED>
   global h           %#ok<NUSED>
   global hbar        %#ok<NUSED>
     
 % Physical Parameters
   global E_H
   global omega_a
   global EField_0
   E_H      = 13.6 * eV;            % ionization potential of hydrogen (J)
   omega_a  = 4.134 * 10^16;        % atomic frequency unit (Hz)
   EField_0 = 5.142 * 10^11;        % atomic field strength (V/m)
   
   N = size(E_t,2);                 % number of the time points
   DeltaTime = time(2) - time(1);   % delta t (unit: sec)   
   
 % Calculation of the electric field of the laser pulse as a function of time (V/m)
   % real form of the electric field (V/m)
   EField_t     = real(E_t);
   % absolute value of the electric field (V/m)   
   EField_t_abs = abs(EField_t);
   % envelope of the electric field (V/m)
   EField_t_envelope = abs(E_t);
   % phase of the electric field (rad)
   EPhase_t = angle(E_t);
   
 % Ionization rate as a function of time (row vector) (Hz)
   W_0_t = StaticIonizationRate(E_ion_0, EField_t_abs);  % 0  -> 1+
   W_1_t = StaticIonizationRate(E_ion_1, EField_t_abs);  % 1+ -> 2+
   W_2_t = StaticIonizationRate(E_ion_2, EField_t_abs);  % 2+ -> 3+
   W_3_t = StaticIonizationRate(E_ion_3, EField_t_abs);  % 3+ -> 4+
   W_4_t = StaticIonizationRate(E_ion_4, EField_t_abs);  % 4+ -> 5+
   
 % Calculation of the relative ion density as a function of time (row vector)
   %col = ones(N,1);
   
   % atom density (relative)
   %disp('atom');
   n_0 = ones(1,N);
   for i = 2:N
       n_0(i) = n_0(i-1) - n_0(i-1) * W_0_t(i-1) * DeltaTime;
       if n_0(i)<0
           n_0(i) = 0;
       end
   end
   
   % 1+ ion density (relative)
   %disp('ion+1');
   n_1 = zeros(1,N);
   for i = 2:N
       n_1(i) = n_1(i-1) - n_1(i-1) * W_1_t(i-1) * DeltaTime  + n_0(i-1) * W_0_t(i-1) * DeltaTime;
       if n_1(i)<0
           n_1(i) = 0;
       end
   end
           
   % 2+ ion density (relative)
   %disp('ion+2');
   n_2 = zeros(1,N);
   for i = 2:N
       n_2(i) = n_2(i-1) - n_2(i-1) * W_2_t(i-1) * DeltaTime + n_1(i-1) * W_1_t(i-1) * DeltaTime;
       if n_2(i)<0
           n_2(i) = 0;
       end
   end
   
   % 3+ ion density (relative)
   %disp('ion+3');
   n_3 = zeros(1,N);
   for i = 2:N
       n_3(i) = n_3(i-1) - n_3(i-1) * W_3_t(i-1) * DeltaTime + n_2(i-1) * W_2_t(i-1) * DeltaTime;
       if n_3(i)<0
           n_3(i) = 0;
       end
   end
   
   % 4+ ion density (relative)
   %disp('ion+4');
   n_4 = zeros(1,N);
   for i = 2:N
       n_4(i) = n_4(i-1) - n_4(i-1) * W_4_t(i-1) * DeltaTime + n_3(i-1) * W_3_t(i-1) * DeltaTime;
       if n_4(i)<0
           n_4(i) = 0;
       end
   end
   
   % 5+ ion density (relative)
   %disp('ion+5');
   n_5 = zeros(1,N);
   for i = 2:N
       n_5(i) = n_5(i-1) + n_4(i-1) * W_4_t(i-1) * DeltaTime;
   end
      
   % error correction
   n_total = n_0 + n_1 + n_2 + n_3 + n_4 + n_5;
   % During the evaluation of n_i, the error is accumulated.
   % This makes the relative total ion density to be unconserved.
   % (i.e. 'n_total' is no longer being 1.)
   % In order to reduce this error, we rescale n_i (i = 1, 2... 8, 9) by n_total.
   n_0 = n_0 ./ n_total;
   n_1 = n_1 ./ n_total;
   n_2 = n_2 ./ n_total;
   n_3 = n_3 ./ n_total;
   n_4 = n_4 ./ n_total;
   n_5 = n_5 ./ n_total;
   %n_total = n_0 + n_1 + n_2 + n_3 + n_4 + n_5;
   
 % electron absorbed energy from ATI heating
   Energy_ATI_t = ( q_e^2 .* EField_t_envelope .^2 ) / ( 2 * m_e * omega_d^2 ) .* sin(EPhase_t).^2;
   % electron ATI energy as a function of time (J) (row vector)
   
   n_re = n_0*0 + n_1*1 + n_2*2 + n_3*3 + n_4*4 + n_5*5;
   % relative electron density as a function of time n_e(t)
   
   dne_dt(2:N) = (n_re(2:N)-n_re(1:N-1))/DeltaTime;
   % variation of relative electron density as a function of time (dn_e/dt)
   
   Energy_ATI_e = sum(dne_dt .* Energy_ATI_t) / sum(dne_dt);  % (unit: J)
   % average ATI energy per electron
   
 % electron absorbed energy for ionization
   Energy_Ionization_e = n_1(N)*E_ion_0 + ...
                         n_2(N)*(E_ion_0 + E_ion_1) + ...
                         n_3(N)*(E_ion_0 + E_ion_1 + E_ion_2) + ...
                         n_4(N)*(E_ion_0 + E_ion_1 + E_ion_2 + E_ion_3) + ...
                         n_5(N)*(E_ion_0 + E_ion_1 + E_ion_2 + E_ion_3 + E_ion_4);
   % average ionization energy per electron
   
 % Return results
   result = [n_re(N), Energy_Ionization_e, Energy_ATI_e];         % unit: [arb. units, J]
 % Data plotting
   if FigureSwitch
      t = time ./ fs;
      figure;
      subplot(4,1,1), plot(t,EField_t),
                      ylabel('laser electric field (V/m)');
      subplot(4,1,2), plot(t,n_0,t,n_1,t,n_2,t,n_3,t,n_4,t,n_5),
                      ylabel('relative ion population');
      subplot(4,1,3), plot(t, n_re), ylabel('relative electron density n_e');
      subplot(4,1,4), plot(t, Energy_ATI_t / eV), ylabel('Energy_{ATI} (eV)'),
      xlabel('time (fs)');
 
      % Export results (usual units)
      result_population(1,:) = t ;         % position (mm)
      result_population(2,:) = n_0;        % position (mm)
      result_population(3,:) = n_1;        % position (mm)
      result_population(4,:) = n_2;        % position (mm)
      result_population(5,:) = n_3;        % position (mm)
      result_population(6,:) = n_4;        % position (mm)
      result_population(7,:) = n_5;        % position (mm)
      result_population(8,:) = n_re;        % position (mm)
      result_population(9,:) = Energy_ATI_t / eV;        % position (mm)
      result_population(10,:) = EField_t;        % position (mm)
      fid_output = fopen(['Population_',num2str(1),'.txt'],'w');       
      fprintf(fid_output,'t     n_0          n_1        n_2         n_3         n_4         n_5         n_re        Energy_ATI_t  EField_t \r\n');
      fprintf(fid_output,'(fs)  (arb.units) (arb.units) (arb.units) (arb.units) (arb.units) (arb.units) (arb.units) (eV)          (V/m)    \r\n');
      fprintf(fid_output,'%f    %f          %f          %f          %f          %f          %f          %f          %f            %f       \r\n', result_population);     
      fclose(fid_output);    
      %pause
   end
end

function result = TunnelingIonizationRate_Linear2(E_t,omega_d,time,...
                      E_ion_0,E_ion_1,E_ion_2,E_ion_3,E_ion_4,FigureSwitch)

 % Unit conversion
   global eV   
   global fs   

 % Physical constants
   global c           %#ok<NUSED>
   global q_e
   global m_e
   global mu_0        %#ok<NUSED>
   global epsilon_0   %#ok<NUSED>
   global h           %#ok<NUSED>
   global hbar        %#ok<NUSED>
     
 % Physical Parameters
   global E_H
   global omega_a
   global EField_0
   E_H      = 13.6 * eV;            % ionization potential of hydrogen (J)
   omega_a  = 4.134 * 10^16;        % atomic frequency unit (Hz)
   EField_0 = 5.142 * 10^11;        % atomic field strength (V/m)
   
   N = size(E_t,2);                 % number of the time points
   DeltaTime = time(2) - time(1);   % delta t (unit: sec)   
   
 % Calculation of the electric field of the laser pulse as a function of time (V/m)
   % real form of the electric field (V/m)
   EField_t     = real(E_t);
   % absolute value of the electric field (V/m)   
   EField_t_abs = abs(EField_t);
   % envelope of the electric field (V/m)
   EField_t_envelope = abs(E_t);
   % phase of the electric field (rad)
   EPhase_t = angle(E_t);
   
 % Ionization rate as a function of time (row vector) (Hz)
   W_0_t = StaticIonizationRate(E_ion_0, EField_t_abs);  % 0  -> 1+
   W_1_t = StaticIonizationRate(E_ion_1, EField_t_abs);  % 1+ -> 2+
   W_2_t = StaticIonizationRate(E_ion_2, EField_t_abs);  % 2+ -> 3+
   W_3_t = StaticIonizationRate(E_ion_3, EField_t_abs);  % 3+ -> 4+
   W_4_t = StaticIonizationRate(E_ion_4, EField_t_abs);  % 4+ -> 5+
   
 % Calculation of the relative ion density as a function of time (row vector)
   %col = ones(N,1);
   
   % atom density (relative)
   %disp('atom');
   n_0 = ones(1,N);
   for i = 2:N
       n_0(i) = n_0(i-1) - n_0(i-1) * W_0_t(i-1) * DeltaTime;
       if n_0(i)<0
           n_0(i) = 0;
       end
   end
   
   % 1+ ion density (relative)
   %disp('ion+1');
   n_1 = zeros(1,N);
   for i = 2:N
       n_1(i) = n_1(i-1) - n_1(i-1) * W_1_t(i-1) * DeltaTime  + n_0(i-1) * W_0_t(i-1) * DeltaTime;
       if n_1(i)<0
           n_1(i) = 0;
       end
   end
           
   % 2+ ion density (relative)
   %disp('ion+2');
   n_2 = zeros(1,N);
   for i = 2:N
       n_2(i) = n_2(i-1) - n_2(i-1) * W_2_t(i-1) * DeltaTime + n_1(i-1) * W_1_t(i-1) * DeltaTime;
       if n_2(i)<0
           n_2(i) = 0;
       end
   end
   
   % 3+ ion density (relative)
   %disp('ion+3');
   n_3 = zeros(1,N);
   for i = 2:N
       n_3(i) = n_3(i-1) - n_3(i-1) * W_3_t(i-1) * DeltaTime + n_2(i-1) * W_2_t(i-1) * DeltaTime;
       if n_3(i)<0
           n_3(i) = 0;
       end
   end
   
   % 4+ ion density (relative)
   %disp('ion+4');
   n_4 = zeros(1,N);
   for i = 2:N
       n_4(i) = n_4(i-1) - n_4(i-1) * W_4_t(i-1) * DeltaTime + n_3(i-1) * W_3_t(i-1) * DeltaTime;
       if n_4(i)<0
           n_4(i) = 0;
       end
   end
   
   % 5+ ion density (relative)
   %disp('ion+5');
   n_5 = zeros(1,N);
   for i = 2:N
       n_5(i) = n_5(i-1) + n_4(i-1) * W_4_t(i-1) * DeltaTime;
   end
      
   % error correction
   n_total = n_0 + n_1 + n_2 + n_3 + n_4 + n_5;
   % During the evaluation of n_i, the error is accumulated.
   % This makes the relative total ion density to be unconserved.
   % (i.e. 'n_total' is no longer being 1.)
   % In order to reduce this error, we rescale n_i (i = 1, 2... 8, 9) by n_total.
   n_0 = n_0 ./ n_total;
   n_1 = n_1 ./ n_total;
   n_2 = n_2 ./ n_total;
   n_3 = n_3 ./ n_total;
   n_4 = n_4 ./ n_total;
   n_5 = n_5 ./ n_total;
   %n_total = n_0 + n_1 + n_2 + n_3 + n_4 + n_5;
   
 % electron absorbed energy from ATI heating
   Energy_ATI_t = ( q_e^2 .* EField_t_envelope .^2 ) / ( 2 * m_e * omega_d^2 ) .* sin(EPhase_t).^2;
   % electron ATI energy as a function of time (J) (row vector)
   
   n_re = n_0*0 + n_1*1 + n_2*2 + n_3*3 + n_4*4 + n_5*5;
   % relative electron density as a function of time n_e(t)
   
   dne_dt(2:N) = (n_re(2:N)-n_re(1:N-1))/DeltaTime;
   % variation of relative electron density as a function of time (dn_e/dt)
   
   Energy_ATI_e = sum(dne_dt .* Energy_ATI_t) / sum(dne_dt);  % (unit: J)
   % average ATI energy per electron
   
 % electron absorbed energy for ionization
   Energy_Ionization_e = n_1(N)*E_ion_0 + ...
                         n_2(N)*(E_ion_0 + E_ion_1) + ...
                         n_3(N)*(E_ion_0 + E_ion_1 + E_ion_2) + ...
                         n_4(N)*(E_ion_0 + E_ion_1 + E_ion_2 + E_ion_3) + ...
                         n_5(N)*(E_ion_0 + E_ion_1 + E_ion_2 + E_ion_3 + E_ion_4);
   % average ionization energy per electron
   
 % Return results
   result = [n_re', n_0', n_1', n_2', n_3', n_4'];         % unit: [arb. units, J]
%  % Data plotting
%    if FigureSwitch
%       t = time ./ fs;
%       figure;
%       subplot(4,1,1), plot(t,EField_t),
%                       ylabel('laser electric field (V/m)');
%       subplot(4,1,2), plot(t,n_0,t,n_1,t,n_2,t,n_3,t,n_4,t,n_5),
%                       ylabel('relative ion population');
%       subplot(4,1,3), plot(t, n_re), ylabel('relative electron density n_e');
%       subplot(4,1,4), plot(t, Energy_ATI_t / eV), ylabel('Energy_{ATI} (eV)'),
%       xlabel('time (fs)');
%  
%       % Export results (usual units)
%       result_population(1,:) = t ;         % position (mm)
%       result_population(2,:) = n_0;        % position (mm)
%       result_population(3,:) = n_1;        % position (mm)
%       result_population(4,:) = n_2;        % position (mm)
%       result_population(5,:) = n_3;        % position (mm)
%       result_population(6,:) = n_4;        % position (mm)
%       result_population(7,:) = n_5;        % position (mm)
%       result_population(8,:) = n_re;        % position (mm)
%       result_population(9,:) = Energy_ATI_t / eV;        % position (mm)
%       result_population(10,:) = EField_t;        % position (mm)
%       fid_output = fopen(['Population_',num2str(1),'.txt'],'w');       
%       fprintf(fid_output,'t     n_0          n_1        n_2         n_3         n_4         n_5         n_re        Energy_ATI_t  EField_t \r\n');
%       fprintf(fid_output,'(fs)  (arb.units) (arb.units) (arb.units) (arb.units) (arb.units) (arb.units) (arb.units) (eV)          (V/m)    \r\n');
%       fprintf(fid_output,'%f    %f          %f          %f          %f          %f          %f          %f          %f            %f       \r\n', result_population);     
%       fclose(fid_output);     
%    end

end

%-------------------------------------------------------------------------
% StaticIonizationRate:
%    Evaluate the ionization rate as a function of ionization potential
%    'E_ion' and static electric field 'EField'.
%-------------------------------------------------------------------------
function result = StaticIonizationRate(E_ion, EField)

   global E_H
   global omega_a
   global EField_0

   a = E_ion / E_H;
   EField = EField + (EField==0)*10^-10;
   b = EField_0 ./ EField;
   
   result = 4 * omega_a * a^(5/2) * b .* exp( (-2/3) * a^(3/2) .* b);

end


%-------------------------------------------------------------------------
% InverseBremsstrahlungCoefficient:
%    Evaluate the inverse Bremsstrahlung absorption coefficient as a
%    function of electron density 'N_e', average ionization state 'Z_ion',
%    and plasma temperature 'kT'.
%-------------------------------------------------------------------------
function result = InverseBremsstrahlungCoefficient(N_e, Z_ion, kT, omega_d)

   global c
   global q_e
   global m_e
   global epsilon_0
   
   % Coulomb logarithm
   Lambda = 6*pi*((epsilon_0)*kT)^1.5 / (q_e^3 * N_e^0.5 * Z_ion);
   
   % plasma frequency
   omega_p = sqrt((q_e^2 * N_e)/(epsilon_0 * m_e));
   
   % refractive index
   n_plasma = sqrt(1-(omega_p^2/omega_d^2));
   
   % IB absorption coefficient (1/m)
   a_IB = (1/(3*c*omega_d^2*n_plasma)) * (q_e^6*Z_ion*N_e^2*log(Lambda))/...
          (2*pi*epsilon_0^2*m_e*kT)^1.5;
   
   result = a_IB;

end


%------------------------------------------------------------------------------
% ThomsonScatteringCoefficient:
%    Evaluate the Thomson scattering coefficient as a function of electron
%    density 'N_e'.
%------------------------------------------------------------------------------
function result = ThomsonScatteringCoefficient(N_e)

   global c
   global q_e
   global m_e
   global epsilon_0
   
   % Thomson scattering cross-section
   sigma_TS = (8*pi/3) * (q_e^4)/(4*pi*epsilon_0*m_e*c^2)^2;
   
   % Thomsone scattering coefficient
   a_TS = N_e * sigma_TS;
      
   result = a_TS;

end


%------------------------------------------------------------------------------
% CapillaryAttenuationCoefficient:
%    Assume EH_11 mode.
%    Evaluate the capillary attenuation coefficient as a function of laser
%    wavelength 'lambda', capillary refractive index 'n_capillary', and
%    capillary radius 'R_capillary'.
%------------------------------------------------------------------------------
function result = CapillaryAttenuationCoefficient(lambda,n_capillary,R_capillary)

   u_11 = 2.405;     % the 1st root of J_0(u_01) = 0.
   
   nu_1 = (1/2) * (n_capillary^2 + 1) / sqrt(n_capillary^2 - 1);
   
   % capillary attenuation coefficient
   a_capillary = (u_11/2/pi)^2 * (lambda^2)/(R_capillary^3) * nu_1;
      
   result = a_capillary;

end

%------------------------------------------------------------------------------
% GroupVelocityDispersion:
%    Evaluate the plasma group-velocity dispersion as a function of laser
%    angular frequency 'omega_d', and electron density 'N_e'.
%------------------------------------------------------------------------------
function result = GroupVelocityDispersion(omega_d,N_e)

   global c
   global q_e
   global m_e
   global epsilon_0
   
   % plasma frequency 
   omega_p = sqrt((N_e * q_e^2)/(epsilon_0 * m_e));
   
   % plasma group-velocity dispersion
   GVD = -(omega_p^2)/c/(omega_d^2 - omega_p^2)^1.5;
      
   result = GVD;

end


%--------------------------------------------------------------------------
% ElectricField_d:
%    Generate the complete electric field for ionization and HHG.
%    E_d(z,t) = 
%       ElectricField_d(t,E_peak(z),omega_d,tau_0,D(z),C(z),phi_d_prop(z))
%
%       t: time (sec)
%       E_peak(z): peak electric field (V/m)
%       omega_d: laser angular frequency (rad/sec)
%       tau_0: initial pulse duration (sec)
%       D(z): accumulated group-delay dispersion (GDD) (sec^2)
%       C(z): accumulated group delay (sec)
%       phi_d_prop(z): accumulated phase shift due to propagation
%                      (k(omega_d)z) (rad)
%
%    The output is a complex number. (unit: V/m)
%    Note: For the ionization calculation, set C=0 and phi_d_prop=0.
%--------------------------------------------------------------------------
function E_d = ElectricField_d(t,E_peak,omega_d,tau_0,D,C,phi_d_prop)

   % phase of the electric field (rad)
   phi_d = (1/2)*atan(D/tau_0^2) - D/2/(tau_0^4+D^2)*(t-C).^2 + ...
              phi_d_prop - omega_d*t;
   
   % Electric field (V/m)
   E_d = E_peak * exp(-(t-C).^2/2/(tau_0^2+D^2/tau_0^2)) .* exp(1i*phi_d);
end



%-------------------------------------------------------------------------
% HHG_DipolePhase:
%    Calculate the intrinsic dipole phase of the HHG.
%
%    HHG_DipolePhase(q,E_d,I_p,omega_d);
%
%       q: harmonic order
%       E_d: driving electric field amplitude (V/m)
%       I_p: material ionization potential (J)
%       omega_d: laser angular frequency (rad/sec)
%
%    output: [Phi_dipole_l  Phi_dipole_s]
%       Phi_dipole_l: long-trajectory dipole phase (rad)
%       Phi_dipole_s: short-trajectory dipole phase (rad)
%
% Author: Hsu-hsin Chu (2021/7/11)
%-------------------------------------------------------------------------
function result = HHG_DipolePhase(q,E_d,I_p,omega_d)

% Physical constants
   global q_e
   global m_e
   global hbar

% Find the ionization time t_0_q for a given harmonic order q and laser
% amplitude E_d: t_0_q(q,E_d)
%    long trajectory:  t_0 = 0 ~ 0.05T
%    short trajectory: t_0 = 0.05T ~ 0.25T
   t_0_q_l = IonizationTime_long(q,E_d,I_p,omega_d);   % unit: sec
   t_0_q_s = IonizationTime_short(q,E_d,I_p,omega_d);  % unit: sec

% Find the recombination time t_r_q for a given harmonic order q and laser
% amplitude E_d: t_r_q(q,E_d)
   t_r_q_l = RecombinationTime(omega_d,t_0_q_l);  % unit: sec
   t_r_q_s = RecombinationTime(omega_d,t_0_q_s);  % unit: sec

% Find the dipole phase for a given harmonic order q as a function of laser
% amplitude E_d: Phi_dipole(E_d)
   Phi_dipole_l = q*omega_d*t_r_q_l - (I_p/hbar)*(t_r_q_l-t_0_q_l) - ...
          (q_e^2*E_d^2)/(2*hbar*m_e*omega_d^2) * (...
            (-1/2/omega_d)*sin(omega_d*t_r_q_l)*cos(omega_d*t_r_q_l) + ...
            (2/omega_d)*cos(omega_d*t_r_q_l)*sin(omega_d*t_0_q_l) - ...
            (3/2/omega_d)*sin(omega_d*t_0_q_l)*cos(omega_d*t_0_q_l) + ...
            (1/2+sin(omega_d*t_0_q_l)^2)*(t_r_q_l-t_0_q_l) );
   Phi_dipole_s = q*omega_d*t_r_q_s - (I_p/hbar)*(t_r_q_s-t_0_q_s) - ...
          (q_e^2*E_d^2)/(2*hbar*m_e*omega_d^2) * (...
            (-1/2/omega_d)*sin(omega_d*t_r_q_s)*cos(omega_d*t_r_q_s) + ...
            (2/omega_d)*cos(omega_d*t_r_q_s)*sin(omega_d*t_0_q_s) - ...
            (3/2/omega_d)*sin(omega_d*t_0_q_s)*cos(omega_d*t_0_q_s) + ...
            (1/2+sin(omega_d*t_0_q_s)^2)*(t_r_q_s-t_0_q_s) );
   
   result = [Phi_dipole_l, Phi_dipole_s];
end
   
  
% -----------------------------------------------------------------------
% Find the recombination time t_r(omega_d,t_0) (unit: sec)
%------------------------------------------------------------------------
function result = RecombinationTime(omega_d,t_0)
   global fs
   
   % unit conversion
   t_0_fs     = t_0/fs;          % unit: fs
   omega_d_fs = omega_d*fs;      % unit: red/fs
   T_fs = 2*pi/omega_d_fs;       % period, unit: fs
   
   % electron position x(t_0,t)
   x_fun = @(t_fs) cos(omega_d_fs*t_fs) - cos(omega_d_fs*t_0_fs) + ...
                   omega_d_fs * sin(omega_d_fs*t_0_fs)*(t_fs - t_0_fs);
   
   % recombination time, determined from x(t_0,t_r) = 0.
   t_r_fs = fzero(x_fun,[t_0_fs+0.001 T_fs]);  % unit: fs
   result = t_r_fs * fs;                       % unit: sec

  % Ref: https://www.mathworks.com/help/matlab/ref/fzero.html
end


% -----------------------------------------------------------------------
% Find the ionization time of the "LONG" trajectory as a function of the
% harmonic order q and the laser amplitude E_d:
%    t_0_q_long(q,E_d)  unit: sec
% long trajectory:  t_0 = 0 ~ 0.05T
% -----------------------------------------------------------------------
function result = IonizationTime_long(q,E_d,I_p,omega_d)
% Unit conversion
   global fs
   global eV

% Physical constants
   global q_e
   global m_e
   global hbar

% angular frequency and period
   omega_d_fs = omega_d*fs;      % unit: rad/fs
   T_fs = 2*pi/omega_d_fs;       % period, unit: fs
   
% photon energy U(t_0,I)    unit: J
   U_fun = @(t_0_fs) I_p + q_e^2*E_d^2/2/m_e/omega_d^2 * ...
      (sin(omega_d*RecombinationTime(omega_d,t_0_fs*fs)) - sin(omega_d_fs*t_0_fs))^2;

% Error function (unit: eV)
   Err = @(t_0_fs) (U_fun(t_0_fs) - q*hbar*omega_d)/eV;
  
% Find the ionization time t_0_q_fs for q-th harmonic
   % long trajectory:  t_0 = 0 ~ 0.05T
   t_0_q_fs = fzero(Err,[0 0.05*T_fs]);    % unit: fs
   result = t_0_q_fs * fs;                 % unit: sec
end


% -----------------------------------------------------------------------
% Find the ionization time of the "SHORT" trajectory as a function of the
% harmonic order q and the laser amplitude E_d:
%    t_0_q_short(q,E_d)  unit: sec
% short trajectory:  t_0 = 0.05T ~ 0.25T
% -----------------------------------------------------------------------
function result = IonizationTime_short(q,E_d,I_p,omega_d)
% Unit conversion
   global fs
   global eV

% Physical constants
   global q_e
   global m_e
   global hbar

% angular frequency and period
   omega_d_fs = omega_d*fs;      % unit: rad/fs
   T_fs = 2*pi/omega_d_fs;       % period, unit: fs
   
% photon energy U(t_0,I)    unit: J
   U_fun = @(t_0_fs) I_p + q_e^2*E_d^2/2/m_e/omega_d^2 * ...
      (sin(omega_d*RecombinationTime(omega_d,t_0_fs*fs)) - sin(omega_d_fs*t_0_fs))^2;

% Error function (unit: eV)
   Err = @(t_0_fs) (U_fun(t_0_fs) - q*hbar*omega_d)/eV;
  
% Find the ionization time t_0_q_fs for q-th harmonic
   % short trajectory:  t_0 = 0.05T ~ 0.25T
   t_0_q_fs = fzero(Err,[0.05*T_fs 0.249*T_fs]);    % unit: fs
   result = t_0_q_fs * fs;                          % unit: sec
end





function propagation = runHHGDrivingPulsePropagation(params, figureSwitch)
%RUNHHGDRIVINGPULSEPROPAGATION Propagate the driving laser pulse.
%   This helper encapsulates Part I of HHG_Gouy_PM_v4. It initialises all
%   arrays, evaluates the plasma response section by section along the
%   propagation axis, stores intermediate quantities, and saves the
%   resulting diagnostics to text files. The implementation is intentionally
%   faithful to the original code in order to guarantee reproducibility.
%
%   INPUTS:
%     params       - structure produced by LOADHHGPM_PARAMETERS.
%     figureSwitch - logical flag controlling the generation of figures.
%
%   OUTPUT:
%     propagation  - structure containing all state variables required by
%                    the subsequent HHG analysis stages. The structure keeps
%                    the naming of the original arrays for clarity.

   arguments
       params struct
       figureSwitch (1,1) double {mustBeMember(figureSwitch,[0 1])}
   end

   % Unpack frequently used parameters from the structure for readability.
   lambda_d    = params.lambda_d;
   PulseEnergy = params.PulseEnergy;
   tau_0       = params.tau_0;
   z_ini       = params.z_ini;
   L           = params.L;
   dz          = params.dz;
   time_initial = params.time_initial;
   DeltaTime    = params.DeltaTime;
   time_final   = params.time_final;
   E_ion_0 = params.E_ion_0;
   E_ion_1 = params.E_ion_1;
   E_ion_2 = params.E_ion_2;
   E_ion_3 = params.E_ion_3;
   E_ion_4 = params.E_ion_4;
   N_gas   = params.N_gas;

   % Global constants -------------------------------------------------------
   global c mu_0 q_e m_e epsilon_0 cm mm um fs mJ

   omega_d = 2 * pi * c / lambda_d;  % laser angular frequency (rad/sec)
   k_0     = 2 * pi / lambda_d;      % laser wavenumber in vacuum (1/m)

   % laser peak intensity at the focal spot (W/m^2)
   I_0 = 2 * PulseEnergy / (pi^1.5 * tau_0 * params.w0^2);
   % laser peak electric field at the focal spot (V/m)
   E_0 = sqrt(2*mu_0*c*I_0);

   % Set the coordinate and variables --------------------------------------
   z_fin = z_ini + L;           % final position (m)
   z     = z_ini:dz:z_fin;      % longitudinal coordinate (unit: m)
   N_z   = numel(z);            % number of data points along z

   tau         = zeros(1,N_z);
   q           = zeros(1,N_z);
   b           = zeros(1,N_z);
   zz          = zeros(1,N_z);
   w           = zeros(1,N_z);
   R           = zeros(1,N_z);
   LaserEnergy = zeros(1,N_z);
   I_peak      = zeros(1,N_z);
   E_peak      = zeros(1,N_z);

   Z_ion       = zeros(3,N_z);
   N_e         = zeros(3,N_z);
   N_e_ave     = zeros(3,N_z);
   N_atom      = zeros(3,N_z);
   N_ion1      = zeros(3,N_z);
   N_ion2      = zeros(3,N_z);

   Energy_Ionization  = zeros(1,N_z);
   Energy_ATI  = zeros(1,N_z);
   kT          = zeros(1,N_z);
   a_IB        = zeros(1,N_z);
   Energy_IB   = zeros(1,N_z);
   a_TS        = zeros(1,N_z);
   Energy_TS   = zeros(1,N_z);
   GVD_plasma  = zeros(1,N_z);
   GVD_Gouy    = zeros(1,N_z);
   D           = zeros(1,N_z);
   n_plasma_d  = zeros(1,N_z);
   C_Gouy      = zeros(1,N_z);
   C           = zeros(1,N_z);
   k_d_plasma  = zeros(1,N_z);
   k_d_total   = zeros(1,N_z);
   Delta_phi_d_Gouy = zeros(1,N_z);
   phi_d_Gouy  = zeros(1,N_z);
   phi_d_prop  = zeros(1,N_z);
   v_d         = zeros(1,N_z);
   ID          = zeros(3,N_z);
   f_plasma    = ones(1,N_z)*1e9;

   time = time_initial:DeltaTime:time_final;
   N_time = numel(time);

   t_z  = zeros(N_time,N_z);
   t_z(:,1) = time';
   N_atom_r0 = zeros(N_time,N_z);
   N_ion1_r0 = zeros(N_time,N_z);
   N_ion2_r0 = zeros(N_time,N_z);
   N_e_r0    = zeros(N_time,N_z);

   % Initial condition of the main pulse propagation -----------------------
   b(1)  = pi * params.w0^2 / lambda_d;
   q(1)  = z(1) - 1i * b(1);
   zz(1) = real(q(1));
   w(1)  = sqrt( (2/k_0) * ((zz(1)^2+b(1)^2)/b(1)) );
   R(1)  = (zz(1)^2 + b(1)^2)/zz(1);
   LaserEnergy(1) = PulseEnergy;
   tau(1)    = tau_0;
   I_peak(1) = I_0 * (params.w0/w(1))^2;
   E_peak(1) = E_0 * (params.w0/w(1));

   E_t = zeros(3,N_time);
   E_t(1,:) = hhgElectricFieldInitial(time,0       ,LaserEnergy(1),q(1),omega_d,tau_0);
   E_t(2,:) = hhgElectricFieldInitial(time,w(1)/2 ,LaserEnergy(1),q(1),omega_d,tau_0);
   E_t(3,:) = hhgElectricFieldInitial(time,w(1)   ,LaserEnergy(1),q(1),omega_d,tau_0);
   Efield_ini = E_t(1,:);

   ATI_result_1 = TunnelingIonizationRate_Linear(...
     E_t(1,:),omega_d,time,E_ion_0,E_ion_1,E_ion_2,E_ion_3,E_ion_4,figureSwitch);
   if figureSwitch
      sgtitle('ionization at z = z_{ini} and r = 0.');
      savefig(gcf,'Fig_1_Ionization_at_z_ini','compact');
   end
   ATI_result_1_0 = ATI_result_1;

   ATI_result_2 = TunnelingIonizationRate_Linear(...
     E_t(2,:),omega_d,time,E_ion_0,E_ion_1,E_ion_2,E_ion_3,E_ion_4,0);
   ATI_result_2_0 = ATI_result_2;

   ATI_result_3 = TunnelingIonizationRate_Linear(...
     E_t(3,:),omega_d,time,E_ion_0,E_ion_1,E_ion_2,E_ion_3,E_ion_4,0);
   ATI_result_3_0 = ATI_result_3;

   Z_ion(:,1) = [ATI_result_1_0{1}; ATI_result_2_0{1}; ATI_result_3_0{1}];
   N_e(:,1)   = N_gas * Z_ion(:,1);

   N_atom(:,1) = N_gas * [ATI_result_1_0{10}; ATI_result_2_0{10}; ATI_result_3_0{10}];
   N_ion1(:,1) = N_gas * [ATI_result_1_0{11}; ATI_result_2_0{11}; ATI_result_3_0{11}];
   N_ion2(:,1) = N_gas * [ATI_result_1_0{12}; ATI_result_2_0{12}; ATI_result_3_0{12}];
   N_e_ave(:,1)= N_gas * [ATI_result_1_0{13}; ATI_result_2_0{13}; ATI_result_3_0{13}];

   Energy_Ionization_e = ATI_result_1_0{2};
   Energy_ATI_e        = ATI_result_1_0{3};

   N_atom_r0(:,1) = ATI_result_1_0{4}' * N_gas;
   N_ion1_r0(:,1) = ATI_result_1_0{5}' * N_gas;
   N_ion2_r0(:,1) = ATI_result_1_0{6}' * N_gas;
   N_e_r0(:,1)    = ATI_result_1_0{9}' * N_gas;

   rr = (0:1*um:3*params.w0)';
   N_rr = numel(rr);
   N_e_fit     = zeros(N_rr,N_z);
   N_e_ave_fit = zeros(N_rr,N_z);
   N_ion1_fit  = zeros(N_rr,N_z);
   N_ion2_fit  = zeros(N_rr,N_z);
   R_N_ion2    = zeros(1,N_z);

   DensityFitResult = hhgDensityDistribution(N_e_ave(:,1),w(1),rr,figureSwitch,1,N_z);
   N_e_ave_fit(:,1) = DensityFitResult{2};
   if figureSwitch
      sgtitle('fitted time-averaged electron density at z_{ini}');
   end

   DensityFitResult = hhgDensityDistribution(N_ion1(:,1),w(1),rr,figureSwitch,1,N_z);
   N_ion1_fit(:,1)  = DensityFitResult{2};
   if figureSwitch
      sgtitle('fitted final 1+ ion density at z_{ini}');
   end

   DensityFitResult = hhgDensityDistribution(N_ion2(:,1),w(1),rr,figureSwitch,1,N_z);
   N_ion2_fit(:,1)  = DensityFitResult{2};
   R_N_ion2(1)      = DensityFitResult{1};
   if figureSwitch
      sgtitle('fitted final 2+ ion density at z_{ini}');
   end

   Energy_Ionization(1) = Energy_Ionization_e * N_e(1,1) * pi * w(1)^2 * dz;
   Energy_ATI(1)        = Energy_ATI_e        * N_e(1,1) * pi * w(1)^2 * dz;
   kT(1)                = Energy_ATI_e * (2/3);
   a_IB(1)              = InverseBremsstrahlungCoefficient(N_e(1,1),Z_ion(1,1),kT(1),omega_d);
   Energy_IB(1)         = LaserEnergy(1) * a_IB(1) * dz;
   a_TS(1)              = ThomsonScatteringCoefficient(N_e(1,1));
   Energy_TS(1)         = LaserEnergy(1) * a_TS(1) * dz;

   GVD_plasma(1) = hhgGroupVelocityDispersionPlasma(omega_d,N_e_ave(1,1));
   GVD_Gouy(1)   = hhgGroupVelocityDispersionGouy(omega_d,q(1));
   D(1)          = (GVD_plasma(1) + GVD_Gouy(1)) * dz;
   C_Gouy(1)     = hhgGroupDelayGouy(omega_d,q(1),dz);
   n_plasma_d(1) = sqrt(1-(q_e^2 * N_e_ave(1,1))/(epsilon_0 * m_e * omega_d^2));
   C(1)          = dz/(c*n_plasma_d(1)) + C_Gouy(1);
   k_d_plasma(1) = k_0 * n_plasma_d(1);
   k_d_Gouy(1)   = -b(1)/(zz(1)^2 + b(1)^2);
   k_d_total(1)  = k_d_plasma(1) + k_d_Gouy(1);
   Delta_phi_d_Gouy(1) = - atan((zz(1) + dz/n_plasma_d(1)) / b(1)) + atan(zz(1)/b(1));
   phi_d_Gouy_ini = -atan(zz(1)/b(1));
   phi_d_Gouy(1)  = phi_d_Gouy_ini;
   phi_d_prop(1)  = phi_d_Gouy(1) + k_d_plasma(1) * dz + Delta_phi_d_Gouy(1);
   v_d(1)         = omega_d/k_d_total(1);
   ID(:,1)        = hhgIonizationDefocusing(N_e_ave(:,1),w(1),dz,omega_d,figureSwitch,1,N_z);
   f_plasma(1)    = ID(3,1);

   for j = 1:N_z-1
      disp([num2str(j),'/',num2str(N_z-1)]);
      t_z(:,j+1) = time' + C(j);

      q(j+1)  = ((q(j)+dz)^(-1) - f_plasma(j)^(-1))^(-1);
      zz(j+1) = real(q(j+1));
      b(j+1)  = -imag(q(j+1));
      w(j+1)  = sqrt( (2/k_0) * ((zz(j+1)^2+b(j+1)^2)/b(j+1)) );
      R(j+1)  = (zz(j+1)^2 + b(j+1)^2)/zz(j+1);

      LaserEnergy(j+1) = LaserEnergy(j) - Energy_Ionization(j) - ...
                         Energy_ATI(j) - Energy_IB(j) - Energy_TS(j);
      tau(j+1)    = sqrt(tau(j)^2 + D(j)^2/tau(j)^2);
      I_peak(j+1) = 2*LaserEnergy(j+1)/(pi^1.5 * tau(j+1) * w(j+1)^2);
      E_peak(j+1) = sqrt(2*mu_0*c*I_peak(j+1));

      E_t(1,:) = hhgElectricField(t_z(:,j+1)',0       ,LaserEnergy(j+1),q(j+1),omega_d,tau_0,C(j),D(j),phi_d_prop(j));
      E_t(2,:) = hhgElectricField(t_z(:,j+1)',w(j+1)/2,LaserEnergy(j+1),q(j+1),omega_d,tau_0,C(j),D(j),phi_d_prop(j));
      E_t(3,:) = hhgElectricField(t_z(:,j+1)',w(j+1)  ,LaserEnergy(j+1),q(j+1),omega_d,tau_0,C(j),D(j),phi_d_prop(j));

      ATI_result_1 = TunnelingIonizationRate_Linear(...
                     E_t(1,:),omega_d,t_z(:,j+1)',E_ion_0,E_ion_1,...
                     E_ion_2,E_ion_3,E_ion_4,and(figureSwitch,j==N_z-1));
      if and(figureSwitch,j==N_z-1)
          sgtitle('ionization at z_{final} and r = 0.');
          savefig(gcf,'Fig_1_Ionization_at_z_fin','compact');
      end

      ATI_result_2 = TunnelingIonizationRate_Linear(...
                     E_t(2,:),omega_d,t_z(:,j+1)',E_ion_0,E_ion_1,...
                     E_ion_2,E_ion_3,E_ion_4,0);

      ATI_result_3 = TunnelingIonizationRate_Linear(...
                     E_t(3,:),omega_d,t_z(:,j+1)',E_ion_0,E_ion_1,...
                     E_ion_2,E_ion_3,E_ion_4,0);

      Z_ion(1,j+1) = ATI_result_1{1};
      Z_ion(2,j+1) = ATI_result_2{1};
      Z_ion(3,j+1) = ATI_result_3{1};

      N_e(:,j+1) = N_gas * Z_ion(:,j+1);
      N_e_ave(1,j+1) = ATI_result_1{13} * N_gas;
      N_e_ave(2,j+1) = ATI_result_2{13} * N_gas;
      N_e_ave(3,j+1) = ATI_result_3{13} * N_gas;
      DensityFitResult = hhgDensityDistribution(N_e_ave(:,j+1),w(j+1),rr,figureSwitch,j+1,N_z);
      N_e_ave_fit(:,j+1) = DensityFitResult{2};
      if and(figureSwitch,j+1==N_z)
         sgtitle('fitted time-averaged electron density at z_{fin}');
      end

      N_atom(1,j+1) = ATI_result_1{10} * N_gas;
      N_atom(2,j+1) = ATI_result_2{10} * N_gas;
      N_atom(3,j+1) = ATI_result_3{10} * N_gas;

      N_ion1(1,j+1) = ATI_result_1{11} * N_gas;
      N_ion1(2,j+1) = ATI_result_2{11} * N_gas;
      N_ion1(3,j+1) = ATI_result_3{11} * N_gas;
      DensityFitResult  = hhgDensityDistribution(N_ion1(:,j+1),w(j+1),rr,figureSwitch,j+1,N_z);
      N_ion1_fit(:,j+1) = DensityFitResult{2};
      if and(figureSwitch,j+1==N_z)
         sgtitle('fitted final 1+ ion density at z_{fin}');
      end

      N_ion2(1,j+1) = ATI_result_1{12} * N_gas;
      N_ion2(2,j+1) = ATI_result_2{12} * N_gas;
      N_ion2(3,j+1) = ATI_result_3{12} * N_gas;
      DensityFitResult  = hhgDensityDistribution(N_ion2(:,j+1),w(j+1),rr,figureSwitch,j+1,N_z);
      N_ion2_fit(:,j+1) = DensityFitResult{2};
      R_N_ion2(j+1)   = DensityFitResult{1};
      if and(figureSwitch,j+1==N_z)
         sgtitle('fitted final 2+ ion density at z_{fin}');
      end

      Energy_Ionization_e = ATI_result_1{2};
      Energy_ATI_e = ATI_result_1{3};
      N_atom_r0(:,j+1) = ATI_result_1{4}' * N_gas;
      N_ion1_r0(:,j+1) = ATI_result_1{5}' * N_gas;
      N_ion2_r0(:,j+1) = ATI_result_1{6}' * N_gas;
      N_e_r0(:,j+1)    = ATI_result_1{9}' * N_gas;

      Energy_Ionization(j+1) = Energy_Ionization_e * N_e(1,j+1) * pi * w(j+1)^2 * dz;
      Energy_ATI(j+1)        = Energy_ATI_e        * N_e(1,j+1) * pi * w(j+1)^2 * dz;
      kT(j+1) = Energy_ATI_e * (2/3);
      a_IB(j+1) = InverseBremsstrahlungCoefficient(N_e(1,j+1),Z_ion(1,j+1),kT(j+1),omega_d);
      Energy_IB(j+1) = LaserEnergy(j+1) * a_IB(j+1) * dz;
      a_TS(j+1) = ThomsonScatteringCoefficient(N_e(1,j+1));
      Energy_TS(j+1) = LaserEnergy(j+1) * a_TS(j+1) * dz;

      GVD_plasma(j+1) = hhgGroupVelocityDispersionPlasma(omega_d,N_e_ave(1,j+1));
      GVD_Gouy(j+1)   = hhgGroupVelocityDispersionGouy(omega_d,q(j+1));
      D(j+1) = D(j) + (GVD_plasma(j+1) + GVD_Gouy(j+1)) * dz;
      C_Gouy(j+1) = hhgGroupDelayGouy(omega_d,q(j+1),dz);
      n_plasma_d(j+1) = sqrt(1-(q_e^2 * N_e_ave(1,j+1))/(epsilon_0 * m_e * omega_d^2));
      C(j+1) = C(j) + dz/(c*n_plasma_d(j+1)) + C_Gouy(j+1);
      k_d_plasma(j+1) = k_0 * n_plasma_d(j+1);
      k_d_Gouy(j+1) = -b(j+1)/(zz(j+1)^2 + b(j+1)^2);
      k_d_total(j+1) = k_d_plasma(j+1) + k_d_Gouy(j+1);
      Delta_phi_d_Gouy(j+1) = - atan((zz(j+1) + dz/n_plasma_d(j+1)) / b(j+1)) + atan(zz(j+1)/b(j+1));
      phi_d_Gouy(j+1) = phi_d_Gouy(j) + Delta_phi_d_Gouy(j+1);
      phi_d_prop(j+1) = phi_d_Gouy(j+1) + k_d_plasma(j+1) * dz + Delta_phi_d_Gouy(j+1);
      v_d(j+1) = omega_d/k_d_total(j+1);
      ID(:,j+1) = hhgIonizationDefocusing(N_e_ave(:,j+1),w(j+1),dz,omega_d,figureSwitch,j+1,N_z);
      f_plasma(j+1) = ID(3,j+1);
   end

   % Export propagation results (SI units) ---------------------------------
   result_SI = [z; LaserEnergy; Z_ion; N_e(1,:); C; GVD_plasma; GVD_Gouy; D; ...
                tau; w; I_peak; E_peak; kT; Energy_Ionization; Energy_ATI; ...
                Energy_IB; Energy_TS; N_e_ave(1,:)];
   fid_output = fopen('Results_1_DrivingPulsePropagation_SI_units.txt','w');
   fprintf(fid_output,'position    laser_energy  ion_state_Z1  ion_state_Z2  ion_state_Z3  e_density   GD         GVD_plasma   GVD_Gouy    GDD      pulse_duration  beam_size_w  peak_intensity  peak_E_field  e_temp_kT   ionization  ATI_energy  IB_energy  TS_energy  N_e_ave\r\n');
   fprintf(fid_output,'(m)         (J)           (arb.units)   (arb.units)   (arb.units)   (1/m^3)     (sec)       (s^2/m)     (s^2/m)    (s^2)         (sec)          (m)           (W/m^2)         (V/m)         (J)      energy(J)     (J)         (J) (J)     (1/m^3)\r\n');
   fprintf(fid_output,'%6.4E  %7.5E   %6.4f        %6.4f        %6.4f        %5.3E  % 5.3E  % 5.3E  % 5.3E  % 5.3E  %8.6E   %5.3E      %5.3E       %5.3E   %5.3E   %5.3E   %5.3E   %5.3E  %5.3E  %5.3E\r\n',result_SI);
   fclose(fid_output);

   result = zeros(20,N_z);
   result(1,:)  = z/mm;
   result(2,:)  = LaserEnergy/mJ;
   result(3,:)  = Z_ion(1,:);
   result(4,:)  = Z_ion(2,:);
   result(5,:)  = Z_ion(3,:);
   result(6,:)  = N_e(1,:)*cm^3;
   result(7,:)  = C/fs;
   result(8,:)  = GVD_plasma/fs^2*mm;
   result(9,:)  = GVD_Gouy/fs^2*mm;
   result(10,:) = D/fs^2;
   result(11,:) = tau/fs;
   result(12,:) = w/mm;
   result(13,:) = I_peak*cm^2;
   result(14,:) = E_peak;
   result(15,:) = kT/eV;
   result(16,:) = Energy_Ionization/mJ;
   result(17,:) = Energy_ATI/mJ;
   result(18,:) = Energy_IB/mJ;
   result(19,:) = Energy_TS/mJ;
   result(20,:) = N_e_ave(1,:)*cm^3;

   fid_output = fopen('Results_2_DrivingPulsePropagation_usual_units.txt','w');
   fprintf(fid_output,'position    laser_energy  ion_state_Z1  ion_state_Z2  ion_state_Z3  e_density   GD         GVD_plasma   GVD_Gouy    GDD      pulse_duration  beam_size_w  peak_intensity  peak_E_field  e_temp_kT   ionization  ATI_energy  IB_energy  TS_energy  N_e_ave\r\n');
   fprintf(fid_output,'(mm)        (mJ)          (arb.units)   (arb.units)   (arb.units)   (1/cm^3)    (fs)       (fs^2/mm)    (fs^2/mm)  (fs^2)         (fs)          (mm)         (W/cm^2)         (V/m)        (eV)     energy(mJ)     (mJ)        (mJ) (mJ)     (1/m^3)\r\n');
   fprintf(fid_output,'%6.4E  %7.5E   %6.4f        %6.4f        %6.4f        %5.3E  % 5.3E  % 5.3E  % 5.3E  % 5.3E  %8.6E   %5.3E      %5.3E       %5.3E   %5.3E   %5.3E   %5.3E   %5.3E  %5.3E  %5.3E\r\n',result);
   fclose(fid_output);

   if figureSwitch
      figure;
      subplot(5,3,1),  plot(z/mm,LaserEnergy/mJ), ylabel('laser energy (mJ)'), xlabel('position z (mm)');
      subplot(5,3,4),  plot(z/mm,E_peak),        ylabel('E_{peak} (V/m)'),    xlabel('position z (mm)');
      subplot(5,3,7),  plot(z/mm,I_peak),        ylabel('I_{peak} (W/m^2)'),  xlabel('position z (mm)');
      subplot(5,3,10), plot(z/mm,w/mm),          ylabel('beam size w(z) (mm)'), xlabel('position z (mm)');
      subplot(5,3,13), plot(z/mm,phi_d_Gouy),    ylabel('Gouy phase (rad)'),  xlabel('position z (mm)');
      subplot(5,3,2), plot(z/mm,N_e(1,:)/cm^(-3),z/mm,N_e_ave(1,:)/cm^(-3)), ylabel('N_e(r=0) (cm^{-3})');
      xlabel('position z (mm)'); legend('final','time-average','Location','NorthEast');
      subplot(5,3,5),  plot(z/mm,GVD_plasma*mm/fs^2,z/mm,GVD_Gouy*mm/fs^2), ylabel('GVD (fs^2/mm)');
      xlabel('position z (mm)'); legend('plasma','Gouy','Location','NorthEast');
      subplot(5,3,8),  plot(z/mm,D/fs^2), ylabel('accumulated GDD (fs^2)'); xlabel('position z (mm)');
      subplot(5,3,11), plot(z/mm,tau/fs), ylabel('pulse duration (fs)'); xlabel('position z (mm)');
      subplot(5,3,3),  plot(z/mm,Energy_Ionization/dz), ylabel('ionization loss (mJ/mm)'); xlabel('position z (mm)');
      subplot(5,3,6),  plot(z/mm,Energy_ATI/dz), ylabel('ATI loss (mJ/mm)'), xlabel('position z (mm)');
      subplot(5,3,9),  plot(z/mm,Energy_IB/dz), ylabel('IB loss (mJ/mm)'),  xlabel('position z (mm)');
      subplot(5,3,12), plot(z/mm,Energy_TS/dz), ylabel('TS loss (mJ/mm)'),  xlabel('position z (mm)');
      sgtitle('Driving pulse propagation');
      savefig(gcf,'Fig_2_DrivingPulsePropagation','compact');
   end

   if figureSwitch
      figure;
      subplot(2,1,1), imagesc([z_ini z_fin]/mm,[0 max(rr)]/um,N_e_ave_fit*cm^3);
      title('fitted time-averaged electron density'); xlabel('z (mm)'); ylabel('r (\mum)');
      aa = colorbar; aa.Label.String = '(cm^{-3})';
      subplot(2,1,2), imagesc([z_ini z_fin]/mm,[0 max(rr)]/um,N_ion2_fit*cm^3);
      title('fitted final 2+ ion density'); xlabel('z (mm)'); ylabel('r (\mum)');
      aa = colorbar; aa.Label.String = '(cm^{-3})';
      sgtitle('fitted electron/ion density distributions');
      savefig(gcf,'Fig_3_FittedhhgDensityDistributions','compact');
   end

   propagation = struct('z', z, 'z_ini', z_ini, 'z_fin', z_fin, 'N_z', N_z, ...
                        'tau', tau, 'q', q, 'b', b, 'zz', zz, 'w', w, 'R', R, ...
                        'LaserEnergy', LaserEnergy, 'I_peak', I_peak, ...
                        'E_peak', E_peak, 'Z_ion', Z_ion, 'N_e', N_e, ...
                        'N_e_ave', N_e_ave, 'N_atom', N_atom, 'N_ion1', N_ion1, ...
                        'N_ion2', N_ion2, 'Energy_Ionization', Energy_Ionization, ...
                        'Energy_ATI', Energy_ATI, 'kT', kT, 'a_IB', a_IB, ...
                        'Energy_IB', Energy_IB, 'a_TS', a_TS, 'Energy_TS', Energy_TS, ...
                        'GVD_plasma', GVD_plasma, 'GVD_Gouy', GVD_Gouy, 'D', D, ...
                        'n_plasma_d', n_plasma_d, 'C_Gouy', C_Gouy, 'C', C, ...
                        'k_d_plasma', k_d_plasma, 'k_d_Gouy', k_d_Gouy, 'k_d_total', k_d_total, ...
                        'Delta_phi_d_Gouy', Delta_phi_d_Gouy, 'phi_d_Gouy', phi_d_Gouy, ...
                        'phi_d_prop', phi_d_prop, 'v_d', v_d, 'ID', ID, ...
                        'f_plasma', f_plasma, 'time', time, 'N_time', N_time, ...
                        't_z', t_z, 'N_atom_r0', N_atom_r0, 'N_ion1_r0', N_ion1_r0, ...
                        'N_ion2_r0', N_ion2_r0, 'N_e_r0', N_e_r0, 'rr', rr, ...
                        'N_e_fit', N_e_fit, 'N_e_ave_fit', N_e_ave_fit, ...
                        'N_ion1_fit', N_ion1_fit, 'N_ion2_fit', N_ion2_fit, ...
                        'R_N_ion2', R_N_ion2, 'I_0', I_0, 'E_0', E_0, ...
                        'omega_d', omega_d, 'k_0', k_0, 'Efield_ini', Efield_ini);
end

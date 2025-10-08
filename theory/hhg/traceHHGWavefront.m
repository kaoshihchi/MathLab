function [hhg, result_scalar] = traceHHGWavefront(params, propagation, dipole, figureSwitch)
%TRACEHHGWAVEFRONT Follow a fixed HHG wavefront and evaluate phase matching.
%   This helper encapsulates Part II of HHG_Gouy_PM_v4.

   arguments
       params struct
       propagation struct
       dipole struct
       figureSwitch (1,1) double {mustBeMember(figureSwitch,[0 1])}
   end

   global mu_0 c q_e m_e epsilon_0 hbar fs mm cm

   lambda_d = params.lambda_d;
   m        = params.m;
   I_p      = params.I_p;
   t0       = params.t0;
   dz       = params.dz;
   DeltaTime = params.DeltaTime;

   z     = propagation.z;
   N_z   = propagation.N_z;
   R     = propagation.R;
   R_N_ion2 = propagation.R_N_ion2;
   C     = propagation.C;
   D     = propagation.D;
   q     = propagation.q;
   w     = propagation.w;
   LaserEnergy = propagation.LaserEnergy;
   phi_d_prop  = propagation.phi_d_prop;
   E_peak      = propagation.E_peak;
   k_d_plasma  = propagation.k_d_plasma;
   k_d_total   = propagation.k_d_total;
   k_d_Gouy    = k_d_total - k_d_plasma;
   phi_d_Gouy  = propagation.phi_d_Gouy;
   N_e_ave     = propagation.N_e_ave;
   N_ion1_r0   = propagation.N_ion1_r0;
   N_e_r0      = propagation.N_e_r0;
   N_time      = propagation.N_time;
   tau_0       = params.tau_0;
   omega_d     = propagation.omega_d;
   k_0         = propagation.k_0;

   lambda_m = lambda_d/m;
   omega_m  = m * omega_d;

   U_m = m*hbar*omega_d;
   I_min = (U_m-I_p)/3.17*(2*m_e*omega_d^2)/(q_e^2*mu_0*c);
   if min(propagation.I_peak) < I_min
       error('traceHHGWavefront:InsufficientIntensity', ...
             'Laser intensity is too low for the requested harmonic order.');
   end

   I_dipole = dipole.I_dipole;
   Phi_dipole_l = dipole.Phi_dipole_l;
   Phi_dipole_s = dipole.Phi_dipole_s;

   q_m         = zeros(1,N_z);
   b_m         = zeros(1,N_z);
   zz_m        = zeros(1,N_z);
   w_m         = zeros(1,N_z);
   R_m         = zeros(1,N_z);
   n_plasma_m  = zeros(1,N_z);
   k_m_plasma  = zeros(1,N_z);
   k_m_Gouy    = zeros(1,N_z);
   k_m_total   = zeros(1,N_z);
   Delta_phi_m_Gouy = zeros(1,N_z);
   phi_m_Gouy  = zeros(1,N_z);
   ID_m        = zeros(3,N_z);
   f_m_plasma  = zeros(1,N_z);
   v_m         = zeros(1,N_z);

   t_HWF       = zeros(1,N_z);
   Delta_t_HWF = zeros(1,N_z);
   i_t_HWF     = zeros(1,N_z);
   N_ion1_r0_HWF = zeros(1,N_z);
   N_e_r0_HWF    = zeros(1,N_z);
   E_d_HWF       = zeros(1,N_z);
   I_d_HWF       = zeros(1,N_z);
   Phi_d_HWF     = zeros(1,N_z);
   W_HWF         = zeros(1,N_z);

   R_m(1)  = R(1);
   w_m(1)  = R_N_ion2(1);
   q_m(1)  = (1/R_m(1) + (1i*lambda_m)/(pi*w_m(1)^2))^(-1);
   zz_m(1) = real(q_m(1));
   b_m(1)  = -imag(q_m(1));
   w0_m    = sqrt(lambda_m*b_m(1)/pi); %#ok<NASGU>
   phi_m_Gouy(1) = -atan(zz_m(1)/b_m(1));

   t_HWF(1) = t0;

   for j = 1:N_z
      Delta_t_HWF(j) = t_HWF(j) - (C(j) - C(1));
      i_t_HWF(j) = (N_time + 1)/2 + Delta_t_HWF(j)/DeltaTime;

      N_ion1_r0_HWF(j) = interp1((1:N_time)',N_ion1_r0(:,j),i_t_HWF(j),'linear');
      N_e_r0_HWF(j)    = interp1((1:N_time)',N_e_r0(:,j),i_t_HWF(j),'linear');

      n_plasma_m(j) = sqrt(1-(q_e^2 * N_e_r0_HWF(j))/(epsilon_0 * m_e * omega_m^2));
      k_m_plasma(j) = m * k_0 * n_plasma_m(j);

      Delta_phi_m_Gouy(j) = - atan((zz_m(j) + dz/n_plasma_m(j))/b_m(j)) ...
                            + atan(zz_m(j) / b_m(j));
      k_m_Gouy(j) = -b_m(j) ./ (zz_m(j).^2 + b_m(j).^2);
      k_m_total(j) = k_m_plasma(j) + k_m_Gouy(j);
      v_m(j) = omega_m ./ k_m_total(j);

      ID_m(:,j) = hhgIonizationDefocusing(N_e_ave(:,j),w(j),dz,omega_m,0,j,N_z);
      f_m_plasma(j) = ID_m(3,j);

      if j < N_z
         q_m(j+1)  = ((q_m(j)+dz)^(-1) - f_m_plasma(j)^(-1))^(-1);
         zz_m(j+1) = real(q_m(j+1));
         b_m(j+1)  = -imag(q_m(j+1));
         w_m(j+1)  = sqrt( (2/k_0) * ((zz_m(j+1)^2+b_m(j+1)^2)/b_m(j+1)) );
         R_m(j+1)  = (zz_m(j+1)^2 + b_m(j+1)^2)/zz_m(j+1);

         phi_m_Gouy(j+1) = phi_m_Gouy(j) + Delta_phi_m_Gouy(j);
         t_HWF(j+1) = t_HWF(j) + dz/v_m(j);
      end

      E_d_HWF(j) = hhgElectricField(t_HWF(j),0,LaserEnergy(j),q(j), ...
                                   omega_d,tau_0,C(j)-C(1),D(j)-D(1), ...
                                   phi_d_prop(j));
      I_d_HWF(j) = abs(E_d_HWF(j)).^2 / (2*mu_0*c);
      Phi_d_HWF(j) = angle(E_d_HWF(j));
      if j > 1
         Delta = Phi_d_HWF(j) - Phi_d_HWF(j-1);
         if Delta > 5
            Phi_d_HWF(j) = Phi_d_HWF(j) - 2*pi();
         end
         if Delta < -5
            Phi_d_HWF(j) = Phi_d_HWF(j) + 2*pi();
         end
      end
      W_HWF(j) = StaticIonizationRate(I_p, abs(E_d_HWF(j)));
   end

   Delta_k_plasma = m * k_d_plasma - k_m_plasma;
   Delta_Phi_plasma   = cumsum(Delta_k_plasma) * dz;
   Delta_Phi_plasma   = Delta_Phi_plasma - Delta_Phi_plasma(1);

   Delta_k_Gouy = m * k_d_Gouy - k_m_Gouy;
   Delta_Phi_Gouy = m * (phi_d_Gouy-phi_d_Gouy(1)) - ...
                    (phi_m_Gouy-phi_m_Gouy(1));

   Phi_dipole_l_HWF = interp1(I_dipole,Phi_dipole_l,I_d_HWF,'spline');
   Delta_Phi_dipole_l_HWF = Phi_dipole_l_HWF - Phi_dipole_l_HWF(1);
   Phi_dipole_s_HWF = interp1(I_dipole,Phi_dipole_s,I_d_HWF,'spline');
   Delta_Phi_dipole_s_HWF = Phi_dipole_s_HWF - Phi_dipole_s_HWF(1);

   Delta_Phi_total_l_HWF = Delta_Phi_plasma + Delta_Phi_Gouy + Delta_Phi_dipole_l_HWF;
   Delta_Phi_total_s_HWF = Delta_Phi_plasma + Delta_Phi_Gouy + Delta_Phi_dipole_s_HWF;

   Phi_LH_l = m * Phi_d_HWF + Phi_dipole_l_HWF;
   Phi_LH_s = m * Phi_d_HWF + Phi_dipole_s_HWF;

   E_LH_l = N_ion1_r0_HWF .* W_HWF .* abs(E_d_HWF).^5 .* exp(1i*Phi_LH_l);
   E_LH_s = N_ion1_r0_HWF .* W_HWF .* abs(E_d_HWF).^5 .* exp(1i*Phi_LH_s);
   E_LH_perfect = N_ion1_r0_HWF .* W_HWF .* abs(E_d_HWF).^5;

   E_HHG_l       = cumsum(E_LH_l);
   E_HHG_s       = cumsum(E_LH_s);
   E_HHG_perfect = cumsum(E_LH_perfect);

   if figureSwitch
      figure;
      subplot(5,3,1), plot(z/mm,Delta_t_HWF/fs); xlabel('z (mm)'); ylabel('\Delta t (fs)'); title('time difference t_{HWF}(z) - C(z)');
      subplot(5,3,4), plot(z/mm,abs(E_d_HWF),z/mm,E_peak); xlabel('z (mm)'); ylabel('|E| (V/m)'); legend('E_{HWF}','E_{peak}','Location','SouthWest'); title('driving field amplitude');
      subplot(5,3,7), plot(z/mm,Phi_d_HWF); xlabel('z (mm)'); ylabel('\Phi_{d\_HWF} (rad)'); title('driving field phase');
      subplot(5,3,10), plot(z/mm,Delta_Phi_plasma); xlabel('z (mm)'); ylabel('\Delta \Phi_{plasma} (rad)'); legend('plasma','Location','NorthEast'); title('accumulated plasma phase mismatch');
      subplot(5,3,13), plot(z/mm,Delta_Phi_Gouy); xlabel('z (mm)'); ylabel('\Delta \Phi_{Gouy} (rad)'); legend('Gouy','Location','NorthEast'); title('accumulated Gouy phase mismatch');
      subplot(5,3,2), plot(z/mm,Delta_Phi_dipole_l_HWF); xlabel('z (mm)'); ylabel('\Delta \Phi_{dipole\_l\_HWF} (rad)'); legend('long','Location','NorthWest'); title('accumulated dipole phase mismatch');
      subplot(5,3,5), plot(z/mm,Delta_Phi_dipole_s_HWF); xlabel('z (mm)'); ylabel('\Delta \Phi_{dipole\_s\_HWF} (rad)'); legend('short','Location','NorthWest'); title('accumulated dipole phase mismatch');
      subplot(5,3,8), plot(z/mm,Phi_LH_l-Phi_LH_l(1),z/mm,Delta_Phi_total_l_HWF); xlabel('z (mm)'); ylabel('\Delta \Phi_{LH\_long} (rad)'); legend('\Phi_{LH\_l}(z)-\Phi_{LH\_l}(0)','total phase mismatch','Location','NorthEast'); title('phase of the long-trajectory LH field');
      subplot(5,3,11), plot(z/mm,Phi_LH_s-Phi_LH_s(1),z/mm,Delta_Phi_total_s_HWF); xlabel('z (mm)'); ylabel('\Delta \Phi_{LH\_short} (rad)'); legend('\Phi_{LH\_s}(z)-\Phi_{LH\_s}(0)','total phase mismatch','Location','NorthEast'); title('phase of the short-trajectory LH field');
      subplot(5,3,14), plot(z/mm,abs(E_HHG_perfect),z/mm,abs(E_HHG_l),z/mm,abs(E_HHG_s)); xlabel('z (mm)'); ylabel('|E_{HHG}| (arb. units)'); legend('perfect','long','short','Location','NorthWest'); title('accumulated harmonic field');
      subplot(5,3,3), plot(z/mm,N_ion1_r0_HWF*cm^3); xlabel('z (mm)'); ylabel('N_1_{HWF} (cm^{-3})'); title('1+ ion density met by the HHG wavefront');
      subplot(5,3,6), plot(z/mm,W_HWF*fs); xlabel('z (mm)'); ylabel('W_{HWF} (1/fs)'); title('ionization rate met by the HHG wavefront');
      subplot(5,3,9), plot(z/mm,N_ion1_r0_HWF .* W_HWF); xlabel('z (mm)'); ylabel('N_{source} (arb. units)'); title('source density');
      sgtitle('On the Fixed HHG Wavefront');
      savefig(gcf,'Fig_5_PM_HHG','compact');
   end

   result_PM_HHG = z/mm;
   result_PM_HHG(2,:)  = abs(E_d_HWF);
   result_PM_HHG(3,:)  = Phi_d_HWF;
   result_PM_HHG(4,:)  = Delta_Phi_plasma;
   result_PM_HHG(5,:)  = Delta_Phi_Gouy;
   result_PM_HHG(6,:)  = Delta_Phi_dipole_l_HWF;
   result_PM_HHG(7,:)  = Delta_Phi_dipole_s_HWF;
   result_PM_HHG(8,:)  = Phi_LH_l - Phi_LH_l(1);
   result_PM_HHG(9,:)  = Phi_LH_s - Phi_LH_s(1);
   result_PM_HHG(10,:) = abs(E_HHG_perfect);
   result_PM_HHG(11,:) = abs(E_HHG_l);
   result_PM_HHG(12,:) = abs(E_HHG_s);
   result_PM_HHG(13,:) = N_ion1_r0_HWF * cm^3;
   result_PM_HHG(14,:) = W_HWF * fs;

   fid_output = fopen('Results_4_PM_HHG.txt','w');
   fprintf(fid_output,'position    |E_d_HWF|   Phi_d_HWF  D_Phi_plasma D_Phi_Gouy  D_Phi_dipole_l D_Phi_dipole_s  D_Phi_LH_l  D_Phi_LH_s |E_HHG_perfect|  |E_HHG_l|   |E_HHG_s|   N_ion1_HWF   W_HWF\r\n');
   fprintf(fid_output,'  (mm)        (V/m)       (rad)        (rad)       (rad)        (rad)         (rad)          (rad)        (rad)        (arb.u)       (arb.u)     (arb.u)     (cm^-3)     (1/fs)\r\n');
   fprintf(fid_output,'%6.4E  %5.3E  % 5.3E   % 5.3E  % 5.3E   % 5.3E    % 5.3E     % 5.3E   % 5.3E   % 5.3E    % 5.3E  % 5.3E % 5.3E  % 5.3E\r\n',result_PM_HHG);
   fclose(fid_output);

   result_scalar = abs(E_HHG_l(N_z)/E_HHG_perfect(N_z));

   hhg = struct('Delta_t_HWF', Delta_t_HWF, 'E_d_HWF', E_d_HWF, ...
                'Phi_d_HWF', Phi_d_HWF, 'Delta_Phi_plasma', Delta_Phi_plasma, ...
                'Delta_Phi_Gouy', Delta_Phi_Gouy, ...
                'Delta_Phi_dipole_l_HWF', Delta_Phi_dipole_l_HWF, ...
                'Delta_Phi_dipole_s_HWF', Delta_Phi_dipole_s_HWF, ...
                'Phi_LH_l', Phi_LH_l, 'Phi_LH_s', Phi_LH_s, ...
                'E_HHG_l', E_HHG_l, 'E_HHG_s', E_HHG_s, ...
                'E_HHG_perfect', E_HHG_perfect, ...
                'N_ion1_r0_HWF', N_ion1_r0_HWF, 'W_HWF', W_HWF);
end

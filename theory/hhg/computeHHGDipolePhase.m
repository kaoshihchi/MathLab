function dipole = computeHHGDipolePhase(params, propagation, figureSwitch)
%COMPUTEHHGDIPOLEPHASE Evaluate long/short dipole phase curves.

   arguments
       params struct
       propagation struct
       figureSwitch (1,1) double {mustBeMember(figureSwitch,[0 1])}
   end

   global mu_0 c cm

   m   = params.m;
   I_p = params.I_p;
   omega_d = propagation.omega_d;

   I_peak = propagation.I_peak;
   N_dipole = 100;
   I_dipole = linspace(min(I_peak)*0.8,max(I_peak)*1.2,N_dipole);
   E_dipole = sqrt(2*mu_0*c*I_dipole);
   Phi_dipole_l = zeros(1,N_dipole);
   Phi_dipole_s = zeros(1,N_dipole);
   for j = 1:N_dipole
       temp = HHG_DipolePhase(m,E_dipole(j),I_p,omega_d);
       Phi_dipole_l(j) = temp(1);
       Phi_dipole_s(j) = temp(2);
   end

   dI = I_dipole(2)-I_dipole(1);
   alpha_l = zeros(1,N_dipole);
   alpha_s = zeros(1,N_dipole);
   alpha_l(1:N_dipole-1) = (Phi_dipole_l(2:N_dipole)-Phi_dipole_l(1:N_dipole-1))/dI;
   alpha_l(N_dipole) = 2*alpha_l(N_dipole-1) - alpha_l(N_dipole-2);
   alpha_s(1:N_dipole-1) = (Phi_dipole_s(2:N_dipole)-Phi_dipole_s(1:N_dipole-1))/dI;
   alpha_s(N_dipole) = 2*alpha_s(N_dipole-1) - alpha_s(N_dipole-2);

   if figureSwitch
      figure;
      subplot(1,4,1), plot(I_dipole*cm^2,Phi_dipole_l,'-o');
         xlabel('intensity (W/cm^2)'), ylabel('\Phi_{dipole\_long} (rad)');
         title('long-trajectory dipole phase');
      subplot(1,4,2), plot(I_dipole*cm^2,Phi_dipole_s,'-o');
         xlabel('intensity (W/cm^2)'), ylabel('\Phi_{dipole\_short} (rad)');
         title('short-trajectory dipole phase');
      subplot(1,4,3), plot(I_dipole,alpha_l,'-o');
         xlabel('intensity (W/m^2)'), ylabel('\alpha_{long} (m^2/W)');
         title('long-trajectory \alpha');
      subplot(1,4,4), plot(I_dipole,alpha_s,'-o');
         xlabel('intensity (W/m^2)'), ylabel('\alpha_{short} (m^2/W)');
         title('short-trajectory \alpha');
      sgtitle('Dipole phase calculation');
      savefig(gcf,'Fig_4_DipolePhaseCalculation','compact');
   end

   result_DipolePhase      = I_dipole;
   result_DipolePhase(2,:) = Phi_dipole_l;
   result_DipolePhase(3,:) = Phi_dipole_s;
   result_DipolePhase(4,:) = alpha_l;
   result_DipolePhase(5,:) = alpha_s;

   fid_output = fopen('Results_3_DipolePhase.txt','w');
   fprintf(fid_output,'intensity  Phi_dipole_l  Phi_dipole_s   alpha_l     alpha_s\r\n');
   fprintf(fid_output,' (W/m^2)      (rad)         (rad)       (m^2/W)     (m^2/W)\r\n');
   fprintf(fid_output,'%5.3E   % 5.3E    % 5.3E  % 5.3E  % 5.3E\r\n',result_DipolePhase);
   fclose(fid_output);

   dipole = struct('I_dipole', I_dipole, 'Phi_dipole_l', Phi_dipole_l, ...
                   'Phi_dipole_s', Phi_dipole_s, 'alpha_l', alpha_l, ...
                   'alpha_s', alpha_s);
end

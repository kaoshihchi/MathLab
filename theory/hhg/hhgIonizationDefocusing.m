function result = hhgIonizationDefocusing(N_e,w,dz,omega,figureSwitch,j,N_z)
%HHGIONIZATIONDEFOCUSING Effective focal length from ionization defocusing.
%   Direct transcription of the original IonizationDefocusing helper.

   global c q_e m_e epsilon_0 um

   k = omega/c;
   r = [0; w/2; w]/um;
   n_plasma = sqrt(1-(N_e * q_e^2)/(epsilon_0 * m_e * omega^2));
   Psi_plasma = k * (n_plasma-1) * dz;

   ft = fittype('A-B*x^2');
   Psi_plasma_cfit = fit(r,Psi_plasma,ft,'StartPoint',[0,1]);

   if and(figureSwitch,j==1)
       figure;
       plot(Psi_plasma_cfit,r,Psi_plasma);
       xlabel('r (\mum)'), ylabel('phase shift (rad)');
       legend('Location','NorthWest');
       sgtitle('Ionization defocusing calculation at z_{ini}');
       savefig(gcf,'Fig_IonizationDefocusingCalculation_ini','compact');
   end

   if and(figureSwitch,j==N_z)
       figure;
       plot(Psi_plasma_cfit,r,Psi_plasma);
       xlabel('r (\mum)'), ylabel('phase shift (rad)');
       legend('Location','NorthWest');
       sgtitle('Ionization defocusing calculation at z_{final}');
       savefig(gcf,'Fig_IonizationDefocusingCalculation_final','compact');
   end

   A = Psi_plasma_cfit.A;
   B = Psi_plasma_cfit.B / um^2;
   f_plasma = k/(2*B);

   result = [A;B;f_plasma];
end

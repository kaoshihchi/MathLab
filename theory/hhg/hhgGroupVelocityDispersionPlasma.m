function GVD_plasma = hhgGroupVelocityDispersionPlasma(omega_d,N_e)
%HHGGROUPVELOCITYDISPERSIONPLASMA Plasma contribution to GVD.

   global c q_e m_e epsilon_0

   omega_p = sqrt((N_e * q_e^2)/(epsilon_0 * m_e));
   GVD_plasma = -(omega_p^2)/c/(omega_d^2 - omega_p^2)^1.5;
end

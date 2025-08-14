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

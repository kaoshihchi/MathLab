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

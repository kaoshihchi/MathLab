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

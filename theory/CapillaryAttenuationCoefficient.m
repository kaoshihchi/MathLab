function result = CapillaryAttenuationCoefficient(lambda,n_capillary,R_capillary)

   u_11 = 2.405;     % the 1st root of J_0(u_01) = 0.
   
   nu_1 = (1/2) * (n_capillary^2 + 1) / sqrt(n_capillary^2 - 1);
   
   % capillary attenuation coefficient
   a_capillary = (u_11/2/pi)^2 * (lambda^2)/(R_capillary^3) * nu_1;
      
   result = a_capillary;

end

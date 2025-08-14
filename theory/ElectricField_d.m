function E_d = ElectricField_d(t,E_peak,omega_d,tau_0,D,C,phi_d_prop)

   % phase of the electric field (rad)
   phi_d = (1/2)*atan(D/tau_0^2) - D/2/(tau_0^4+D^2)*(t-C).^2 + ...
              phi_d_prop - omega_d*t;
   
   % Electric field (V/m)
   E_d = E_peak * exp(-(t-C).^2/2/(tau_0^2+D^2/tau_0^2)) .* exp(1i*phi_d);
end

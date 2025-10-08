function E_d = hhgElectricField(t,r,LaserEnergy,q,omega_d,tau_0,C,D,phi_d_prop)
%HHGELECTRICFIELD Full complex field used along propagation.
%   Mirrors the ElectricField_d helper embedded in HHG_Gouy_PM_v4.

   global c mu_0

   tau_z = sqrt(tau_0^2 + D^2/tau_0^2);
   zz = real(q);
   b = -imag(q);
   k_0 = omega_d/c;
   w  = sqrt( (2/k_0) * ((zz^2+b^2)/b) );
   w0 = sqrt( (2/k_0) * b );

   E_peak = sqrt((4*mu_0*c*LaserEnergy)/(pi^(1.5)*tau_z*w^2));

   phi_d = (1/2)*atan(D/tau_0^2) - D/(2*(tau_0^4+D^2))*(t-C).^2 ...
           + phi_d_prop - omega_d*t;

   E_d = E_peak * exp(1i*k_0*r^2/(2*q)) * exp(-(t-C).^2/(2*tau_z^2)) ...
         .* exp(1i*phi_d);
end

function E_d_ini = hhgElectricFieldInitial(t,r,LaserEnergy,q,omega_d,tau_0)
%HHGELECTRICFIELDINITIAL Initial complex Gaussian field at z = z_ini.
%   This is a direct extraction of the ElectricField_d_ini local function
%   from HHG_Gouy_PM_v4. The notation is preserved for transparency.

   global c mu_0

   zz = real(q);
   b = -imag(q);
   k_0 = omega_d/c;
   w  = sqrt( (2/k_0) * ((zz^2+b^2)/b) );
   w0 = sqrt( (2/k_0) * b );

   E_peak = sqrt((4*mu_0*c*LaserEnergy)/(pi^(1.5)*tau_0*w^2));
   E_0    = (w/w0) * E_peak;

   E_d_ini = E_0 * (-1i*b/q) * exp(1i*k_0*r^2/(2*q)) * ...
             exp(-(t.^2)/(2*tau_0^2)) .* exp(-1i * omega_d * t);
end

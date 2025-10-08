function GVD_Gouy = hhgGroupVelocityDispersionGouy(omega_d,q)
%HHGGROUPVELOCITYDISPERSIONGOUY Gouy phase contribution to GVD.

   global c

   zz = real(q);
   b  = -imag(q);
   w0 = sqrt(2*c*b/omega_d);

   GVD_Gouy = 4 * omega_d * (12 * c^3 * w0^6 * zz^2 - c * w0^10 * omega_d^2)...
         / (4 * c^2 * zz^2 + w0^4 * omega_d^2)^3;
end

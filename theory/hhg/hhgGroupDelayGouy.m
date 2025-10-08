function C_Gouy = hhgGroupDelayGouy(omega_d,q,dz)
%HHGGROUPDELAYGOUY Group delay accumulated across section dz.

   global c

   z1 = real(q);
   z2 = z1 + dz;
   b  = -imag(q);
   w0 = sqrt(2*c*b/omega_d);

   C_Gouy = (2 * c * w0^2 * z2) / (4 * c^2 * z2^2 + w0^4 * omega_d^2) - ...
            (2 * c * w0^2 * z1) / (4 * c^2 * z1^2 + w0^4 * omega_d^2);
end

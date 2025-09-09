% ElectricField_d
% -------------------------------------------------------------------------
% Generate the complex electric field of a Gaussian laser pulse 
% propagating through a dispersive medium (plasma + capillary), 
% including envelope delay, group-delay dispersion, and accumulated 
% propagation phase.
%
% SYNTAX:
%   E_d = ElectricField_d(t, E_peak, omega_d, tau_0, D, C, phi_d_prop)
%
% INPUTS:
%   t          : Time array [s]
%   E_peak     : Peak electric field amplitude at position z [V/m]
%   omega_d    : Driving laser angular frequency [rad/s]
%   tau_0      : Initial pulse duration (FWHM/1.665) at z=0 [s]
%   D          : Accumulated group-delay dispersion (GDD) at position z [s^2]
%   C          : Group delay (arrival time shift of the pulse envelope) [s]
%   phi_d_prop : Accumulated propagation phase of the driving field [rad]
%
% OUTPUT:
%   E_d        : Complex electric field E(t) at the given position [V/m]
%
% MODEL:
%   - The Gaussian pulse envelope is broadened by GDD D, resulting in 
%     a stretched temporal width:
%         τ_eff = sqrt( τ_0^2 + D^2 / τ_0^2 ).
%   - The instantaneous phase contains:
%       • quadratic chirp term from GDD
%       • group delay shift C
%       • accumulated propagation phase phi_d_prop
%       • oscillation term at frequency ω_d
%
% EQUATIONS:
%   φ_d(t) = 0.5*atan(D/τ_0^2) - (D / (2(τ_0^4+D^2))) * (t-C)^2 
%            + φ_d_prop - ω_d t
%
%   E_d(t) = E_peak * exp(-(t-C)^2 / (2τ_eff^2)) * exp(i φ_d(t))
%
% NOTES:
%   - This field is expressed in complex form. The physical, real-valued 
%     electric field is Re[E_d(t)].
%   - Designed for use in HHG propagation models, where plasma dispersion, 
%     capillary dispersion, and accumulated phase are provided externally.
%
% -------------------------------------------------------------------------


function E_d = ElectricField_d(t,E_peak,omega_d,tau_0,D,C,phi_d_prop)

   % phase of the electric field (rad)
   phi_d = (1/2)*atan(D/tau_0^2) - D/2/(tau_0^4+D^2)*(t-C).^2 + ...
              phi_d_prop - omega_d*t;
   
   % Electric field (V/m)
   E_d = E_peak * exp(-(t-C).^2/2/(tau_0^2+D^2/tau_0^2)) .* exp(1i*phi_d);
end

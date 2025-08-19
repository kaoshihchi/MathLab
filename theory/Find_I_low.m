function [I_low_Wm2, I_low_Wcm2, Up_eV] = Find_I_low(q, omega_d)
% Find_I_low  Closed-form solution for the lower intensity threshold I_low.
% Inputs:
%   q        : harmonic order (integer)
%   omega_d  : driving angular frequency (rad/s)
% Outputs:
%   I_low_Wm2  : I_low in W/m^2
%   I_low_Wcm2 : I_low in W/cm^2
%   Up_eV      : ponderomotive energy at I_low (eV)

    % ---- Physical constants (SI) ----
    qe   = 1.602176634e-19;       % C
    me   = 9.1093837015e-31;      % kg
    hbar = 1.054571817e-34;       % J*s
    c    = 299792458;             % m/s
    eps0 = 8.8541878128e-12;      % F/m
    eV   = 1.602176634e-19;       % J

    % ---- From (q*hbar*omega / (2*Up))^(1/2) = 1.2596  ----
    K = 1.2596;
    Up_J  = (q*hbar*omega_d) / (2*K^2);                        % J
    I_low_Wm2 = Up_J * (2*me*omega_d^2*c*eps0) / (qe^2);       % W/m^2

    % Convenience outputs
    I_low_Wcm2 = I_low_Wm2 / 1e4;                              % W/cm^2
    Up_eV      = Up_J / eV;                                    % eV
end

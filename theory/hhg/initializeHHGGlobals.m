function constants = initializeHHGGlobals()
%INITIALIZEHHGGLOBALS Set unit conversion factors and physical constants.
%   The original HHG_Gouy_PM_v4 implementation relies on a collection of
%   global variables for unit conversion and core physical constants. This
%   helper gathers those assignments in a single place so that the
%   refactored workflow can initialise them once and share the same values
%   across the newly extracted functions.
%
%   OUTPUT:
%     constants : structure containing the numerical values for quick
%                 reference inside the analysis script if desired.
%
%   The function keeps the "global" assignments to preserve compatibility
%   with legacy routines that still expect them.

   % Unit conversion factors ------------------------------------------------
   global eV cm mm um nm fs mJ

   eV   = 1.602e-19;    % 1 (eV) = 1.6e-19 (J)
   cm   = 1e-2;         % 1 (cm) = 1e-2 (m)
   mm   = 1e-3;         % 1 (mm) = 1e-3 (m)
   um   = 1e-6;         % 1 (um) = 1e-6 (m)
   nm   = 1e-9;         % 1 (nm) = 1e-9 (m)
   fs   = 1e-15;        % 1 (fs) = 1e-15 (sec)
   mJ   = 1e-3;         % 1 (mJ) = 1e-3 (J)

   % Physical constants -----------------------------------------------------
   global c q_e m_e mu_0 epsilon_0 h hbar

   c    = 2.99792458e8;    % speed of light (m/sec)
   q_e  = 1.602e-19;       % electron charge (Coul)
   m_e  = 9.101e-31;       % electron mass (kg)
   mu_0 = 4*pi*1e-7;       % permeability of vacuum (Nt/Amp^2)
   epsilon_0 = 8.854e-12;  % permittivity of vacuum (F/m)
   h    = 6.626e-34;       % Planck's constant (J-sec)
   hbar = h/(2*pi);

   constants = struct('eV', eV, 'cm', cm, 'mm', mm, 'um', um, 'nm', nm, ...
                      'fs', fs, 'mJ', mJ, 'c', c, 'q_e', q_e, 'm_e', m_e, ...
                      'mu_0', mu_0, 'epsilon_0', epsilon_0, 'h', h, ...
                      'hbar', hbar);
end

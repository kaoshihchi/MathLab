function result = IonizationTime_short(q,E_d,I_p,omega_d)
% IonizationTime_short
% Returns ionization time t0 (in seconds) for the q-th harmonic (short trajectory).
% If no root is found in [0.05*T, 0.249*T], returns NaN.
%
% Inputs:
%   q        : harmonic order (integer)
%   E_d      : driving-field amplitude (V/m)
%   I_p      : ionization potential (J)
%   omega_d  : angular frequency (rad/s)
%
% Output:
%   result   : ionization time t0 in seconds (NaN if not found)

   % ---- Globals (same as your code) ----
   global fs         % 1 fs in seconds
   global eV         % 1 eV in Joules
   global q_e        % electron charge (C)
   global m_e        % electron mass (kg)
   global hbar       % reduced Planck constant (J*s)

   % ---- Angular frequency and period in fs ----
   omega_d_fs = omega_d * fs;      % rad/fs
   T_fs       = 2*pi / omega_d_fs; % fs

   % ---- Photon energy U(t0) [J] ----
   % NOTE: RecombinationTime(omega_d, t0_sec) must be defined elsewhere.
   U_fun = @(t0_fs) I_p + (q_e^2 * E_d^2) / (2 * m_e * omega_d^2) .* ...
       ( sin( omega_d * RecombinationTime(omega_d, t0_fs * fs) ) ...
       - sin( omega_d_fs * t0_fs ) ).^2;

   % ---- Error function in eV (root when U = q*hbar*omega_d) ----
   Err = @(t0_fs) ( U_fun(t0_fs) - q * hbar * omega_d ) / eV;

   % ---- Bracket for short trajectory ----
   a = 0.05  * T_fs;
   b = 0.249 * T_fs;   % keep your original upper bound

   % ---- Safe evaluation at endpoints ----
   fa = safeEval(Err, a);
   fb = safeEval(Err, b);

   % Default: no solution
   t0_fs = NaN;

   % Proceed only if endpoints are finite and exhibit a sign change (or zero at an endpoint)
   if isfinite(fa) && isfinite(fb) && (fa * fb <= 0)
       try
           t0_fs = fzero(Err, [a b]);  % fs
       catch
           % leave t0_fs = NaN on failure
       end
   else
       % Optional: enable this warning for debugging
       % warning('IonizationTime_short:NoRoot', ...
       %   'No sign change in [0.05*T, 0.249*T]. Returning NaN.');
   end

   % ---- Convert to seconds ----
   if isnan(t0_fs)
       result = NaN;
   else
       result = t0_fs * fs;  % s
   end
end

% -------- Local helper: safe evaluation that returns NaN on error --------
function val = safeEval(fun, x)
    try
        val = fun(x);
        if ~isfinite(val), val = NaN; end
    catch
        val = NaN;
    end
end

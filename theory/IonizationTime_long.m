function result = IonizationTime_long(q,E_d,I_p,omega_d)
% IonizationTime_long
% Returns ionization time t0 (in seconds) for the q-th harmonic (long trajectory).
% If no root is found in [0, 0.05*T], returns NaN.
%
% Inputs:
%   q        : harmonic order (integer)
%   E_d      : driving-field amplitude (V/m)
%   I_p      : ionization potential (J)
%   omega_d  : angular frequency (rad/s)
%
% Output:
%   result   : ionization time t0 in seconds (NaN if not found)

   % ---- Globals (kept as in your code) ----
   global fs         % 1 fs in seconds
   global eV         % 1 eV in Joules
   global q_e        % electron charge (C)
   global m_e        % electron mass (kg)
   global hbar       % reduced Planck (J*s)

   % ---- Angular frequency and period in fs ----
   omega_d_fs = omega_d * fs;     % rad/fs
   T_fs       = 2*pi / omega_d_fs; % fs

   % ---- Photon energy function U(t0) [J] ----
   % NOTE: RecombinationTime(omega_d, t0_sec) must be defined elsewhere.
   U_fun = @(t0_fs) I_p + (q_e^2 * E_d^2) / (2 * m_e * omega_d^2) .* ...
       ( sin( omega_d * RecombinationTime(omega_d, t0_fs * fs) ) ...
       - sin( omega_d_fs * t0_fs ) ).^2;

   % ---- Error function in eV (root when U = q*hbar*omega_d) ----
   Err = @(t0_fs) ( U_fun(t0_fs) - q * hbar * omega_d ) / eV;

   % ---- Bracket for long trajectory ----
   a = 0;
   b = 0.05 * T_fs;

   % ---- Safe evaluation at endpoints ----
   fa = safeEval(Err, a);
   fb = safeEval(Err, b);

   % Default: no solution
   t0_fs = NaN;

   % Proceed only if endpoints are finite and exhibit a sign change
   if isfinite(fa) && isfinite(fb) && (fa * fb <= 0)
       % Try fzero; if it fails, keep NaN
       try
           t0_fs = fzero(Err, [a b]);  % fs
       catch
           % leave t0_fs = NaN
       end
   else
       % No sign change or invalid endpoints -> NaN
       % (Optional) Uncomment to see a warning:
       % warning('IonizationTime_long:NoRoot',...
       %     'No sign change in [0, 0.05*T]. Returning NaN.');
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

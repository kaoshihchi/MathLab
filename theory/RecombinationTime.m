function result = RecombinationTime(omega_d, t_0)
% RecombinationTime
% Returns recombination time t_r (in seconds) for a given ionization time t_0 (seconds).
% If no root of x(t0, t_r)=0 exists in (t0, T], returns NaN.
%
% Inputs:
%   omega_d : angular frequency (rad/s)
%   t_0     : ionization time (s)
%
% Output:
%   result  : recombination time t_r (s), or NaN if not found.

    global fs   % 1 fs in seconds

    % ---- Unit conversions ----
    t_0_fs     = t_0 / fs;          % fs
    omega_d_fs = omega_d * fs;      % rad/fs
    T_fs       = 2*pi / omega_d_fs; % fs
    eps_fs     = 1e-3;              % small offset (fs) to avoid singular bracket

    % ---- Electron displacement x(t0, t) in fs-domain (unchanged math) ----
    x_fun = @(t_fs) ...
        cos(omega_d_fs * t_fs) - cos(omega_d_fs * t_0_fs) + ...
        omega_d_fs * sin(omega_d_fs * t_0_fs) .* (t_fs - t_0_fs);

    % ---- Try the original bracket first: [t0+eps, T] ----
    a = t_0_fs + eps_fs;
    b = T_fs;

    % Safe evaluation at endpoints
    fa = safeEval(x_fun, a);
    fb = safeEval(x_fun, b);

    t_r_fs = NaN; % default

    if isfinite(fa) && isfinite(fb) && (fa * fb <= 0)
        % Endpoints bracket a root (or zero at an endpoint)
        t_r_fs = tryFzero(x_fun, [a b]);
    else
        % No sign change at endpoints: scan for a bracketing sub-interval
        % (keeps math unchanged; only improves robustness)
        Nscan = 200;  % number of scan points
        ts = linspace(a, b, Nscan+1);
        vals = arrayfun(@(t) safeEval(x_fun, t), ts);

        % Find first adjacent pair with sign change (or zero)
        idx = find(isfinite(vals(1:end-1)) & isfinite(vals(2:end)) & (vals(1:end-1).*vals(2:end) <= 0), 1, 'first');

        if ~isempty(idx)
            aa = ts(idx);
            bb = ts(idx+1);
            t_r_fs = tryFzero(x_fun, [aa bb]);
        else
            % No bracket found -> NaN
            % (Optional) warning:
            % warning('RecombinationTime:NoRoot', 'No sign change found; returning NaN.');
        end
    end

    % ---- Convert to seconds ----
    if isnan(t_r_fs)
        result = NaN;
    else
        result = t_r_fs * fs;
    end
end

% ---------- helpers ----------
function val = safeEval(fun, x)
    try
        val = fun(x);
        if ~isfinite(val), val = NaN; end
    catch
        val = NaN;
    end
end

function root = tryFzero(fun, bracket)
    try
        root = fzero(fun, bracket);
    catch
        root = NaN;
    end
end

function result = RecombinationTime(omega_d,t_0)
   global fs
   
   % unit conversion
   t_0_fs     = t_0/fs;          % unit: fs
   omega_d_fs = omega_d*fs;      % unit: red/fs
   T_fs = 2*pi/omega_d_fs;       % period, unit: fs
   
   % electron position x(t_0,t)
   x_fun = @(t_fs) cos(omega_d_fs*t_fs) - cos(omega_d_fs*t_0_fs) + ...
                   omega_d_fs * sin(omega_d_fs*t_0_fs)*(t_fs - t_0_fs);
   
   % recombination time, determined from x(t_0,t_r) = 0.
   t_r_fs = fzero(x_fun,[t_0_fs+0.001 T_fs]);  % unit: fs
   result = t_r_fs * fs;                       % unit: sec

  % Ref: https://www.mathworks.com/help/matlab/ref/fzero.html
end

function result = IonizationTime_long(q,E_d,I_p,omega_d)
% Unit conversion
   global fs
   global eV

% Physical constants
   global q_e
   global m_e
   global hbar

% angular frequency and period
   omega_d_fs = omega_d*fs;      % unit: rad/fs
   T_fs = 2*pi/omega_d_fs;       % period, unit: fs
   
% photon energy U(t_0,I)    unit: J
   U_fun = @(t_0_fs) I_p + q_e^2*E_d^2/2/m_e/omega_d^2 * ...
      (sin(omega_d*RecombinationTime(omega_d,t_0_fs*fs)) - sin(omega_d_fs*t_0_fs))^2;

% Error function (unit: eV)
   Err = @(t_0_fs) (U_fun(t_0_fs) - q*hbar*omega_d)/eV;
  
% Find the ionization time t_0_q_fs for q-th harmonic
   % long trajectory:  t_0 = 0 ~ 0.05T
   t_0_q_fs = fzero(Err,[0 0.05*T_fs]);    % unit: fs
   result = t_0_q_fs * fs;                 % unit: sec
end

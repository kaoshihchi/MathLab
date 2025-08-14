function result = HHG_DipolePhase(q,E_d,I_p,omega_d)

% Physical constants
   global q_e
   global m_e
   global hbar

% Find the ionization time t_0_q for a given harmonic order q and laser
% amplitude E_d: t_0_q(q,E_d)
%    long trajectory:  t_0 = 0 ~ 0.05T
%    short trajectory: t_0 = 0.05T ~ 0.25T
   t_0_q_l = IonizationTime_long(q,E_d,I_p,omega_d);   % unit: sec
   t_0_q_s = IonizationTime_short(q,E_d,I_p,omega_d);  % unit: sec

% Find the recombination time t_r_q for a given harmonic order q and laser
% amplitude E_d: t_r_q(q,E_d)
   t_r_q_l = RecombinationTime(omega_d,t_0_q_l);  % unit: sec
   t_r_q_s = RecombinationTime(omega_d,t_0_q_s);  % unit: sec

% Find the dipole phase for a given harmonic order q as a function of laser
% amplitude E_d: Phi_dipole(E_d)
   Phi_dipole_l = q*omega_d*t_r_q_l - (I_p/hbar)*(t_r_q_l-t_0_q_l) - ...
          (q_e^2*E_d^2)/(2*hbar*m_e*omega_d^2) * (...
            (-1/2/omega_d)*sin(omega_d*t_r_q_l)*cos(omega_d*t_r_q_l) + ...
            (2/omega_d)*cos(omega_d*t_r_q_l)*sin(omega_d*t_0_q_l) - ...
            (3/2/omega_d)*sin(omega_d*t_0_q_l)*cos(omega_d*t_0_q_l) + ...
            (1/2+sin(omega_d*t_0_q_l)^2)*(t_r_q_l-t_0_q_l) );
   Phi_dipole_s = q*omega_d*t_r_q_s - (I_p/hbar)*(t_r_q_s-t_0_q_s) - ...
          (q_e^2*E_d^2)/(2*hbar*m_e*omega_d^2) * (...
            (-1/2/omega_d)*sin(omega_d*t_r_q_s)*cos(omega_d*t_r_q_s) + ...
            (2/omega_d)*cos(omega_d*t_r_q_s)*sin(omega_d*t_0_q_s) - ...
            (3/2/omega_d)*sin(omega_d*t_0_q_s)*cos(omega_d*t_0_q_s) + ...
            (1/2+sin(omega_d*t_0_q_s)^2)*(t_r_q_s-t_0_q_s) );
   
   result = [Phi_dipole_l, Phi_dipole_s];
end

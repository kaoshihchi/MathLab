function result = TunnelingIonizationRate_Linear(E_t,omega_d,time,...
                      E_ion_0,E_ion_1,E_ion_2,E_ion_3,E_ion_4,FigureSwitch)

 % Unit conversion
   global eV   
   global fs   

 % Physical constants
   global c           %#ok<NUSED>
   global q_e
   global m_e
   global mu_0        %#ok<NUSED>
   global epsilon_0   %#ok<NUSED>
   global h           %#ok<NUSED>
   global hbar        %#ok<NUSED>
     
 % Physical Parameters
   global E_H
   global omega_a
   global EField_0
   E_H      = 13.6 * eV;            % ionization potential of hydrogen (J)
   omega_a  = 4.134 * 10^16;        % atomic frequency unit (Hz)
   EField_0 = 5.142 * 10^11;        % atomic field strength (V/m)
   
   N = size(E_t,2);                 % number of the time points
   DeltaTime = time(2) - time(1);   % delta t (unit: sec)   
   
 % Calculation of the electric field of the laser pulse as a function of time (V/m)
   % real form of the electric field (V/m)
   EField_t     = real(E_t);
   % absolute value of the electric field (V/m)   
   EField_t_abs = abs(EField_t);
   % envelope of the electric field (V/m)
   EField_t_envelope = abs(E_t);
   % phase of the electric field (rad)
   EPhase_t = angle(E_t);

   I_t = abs(E_t).^2 / (2*mu_0*c);
   
 % Ionization rate as a function of time (row vector) (Hz)
   W_0_t = StaticIonizationRate(E_ion_0, EField_t_abs);  % 0  -> 1+
   W_1_t = StaticIonizationRate(E_ion_1, EField_t_abs);  % 1+ -> 2+
   W_2_t = StaticIonizationRate(E_ion_2, EField_t_abs);  % 2+ -> 3+
   W_3_t = StaticIonizationRate(E_ion_3, EField_t_abs);  % 3+ -> 4+
   W_4_t = StaticIonizationRate(E_ion_4, EField_t_abs);  % 4+ -> 5+
   
 % Calculation of the relative ion density as a function of time (row vector)
   %col = ones(N,1);
   
   % atom density (relative)
   %disp('atom');
   n_0 = ones(1,N);
   for i = 2:N
       n_0(i) = n_0(i-1) - n_0(i-1) * W_0_t(i-1) * DeltaTime;
       if n_0(i)<0
           n_0(i) = 0;
       end
   end
   
   % 1+ ion density (relative)
   %disp('ion+1');
   n_1 = zeros(1,N);
   for i = 2:N
       n_1(i) = n_1(i-1) - n_1(i-1) * W_1_t(i-1) * DeltaTime  + n_0(i-1) * W_0_t(i-1) * DeltaTime;
       if n_1(i)<0
           n_1(i) = 0;
       end
   end
           
   % 2+ ion density (relative)
   %disp('ion+2');
   n_2 = zeros(1,N);
   for i = 2:N
       n_2(i) = n_2(i-1) - n_2(i-1) * W_2_t(i-1) * DeltaTime + n_1(i-1) * W_1_t(i-1) * DeltaTime;
       if n_2(i)<0
           n_2(i) = 0;
       end
   end
   
   % 3+ ion density (relative)
   %disp('ion+3');
   n_3 = zeros(1,N);
   for i = 2:N
       n_3(i) = n_3(i-1) - n_3(i-1) * W_3_t(i-1) * DeltaTime + n_2(i-1) * W_2_t(i-1) * DeltaTime;
       if n_3(i)<0
           n_3(i) = 0;
       end
   end
   
   % 4+ ion density (relative)
   %disp('ion+4');
   n_4 = zeros(1,N);
   for i = 2:N
       n_4(i) = n_4(i-1) - n_4(i-1) * W_4_t(i-1) * DeltaTime + n_3(i-1) * W_3_t(i-1) * DeltaTime;
       if n_4(i)<0
           n_4(i) = 0;
       end
   end
   
   % 5+ ion density (relative)
   %disp('ion+5');
   n_5 = zeros(1,N);
   for i = 2:N
       n_5(i) = n_5(i-1) + n_4(i-1) * W_4_t(i-1) * DeltaTime;
   end
      
   % error correction
   n_total = n_0 + n_1 + n_2 + n_3 + n_4 + n_5;
   % During the evaluation of n_i, the error is accumulated.
   % This makes the relative total ion density to be unconserved.
   % (i.e. 'n_total' is no longer being 1.)
   % In order to reduce this error, we rescale n_i (i = 1, 2... 8, 9) by n_total.
   n_0 = n_0 ./ n_total;
   n_1 = n_1 ./ n_total;
   n_2 = n_2 ./ n_total;
   n_3 = n_3 ./ n_total;
   n_4 = n_4 ./ n_total;
   n_5 = n_5 ./ n_total;
   %n_total = n_0 + n_1 + n_2 + n_3 + n_4 + n_5;
   
 % electron absorbed energy from ATI heating
   Energy_ATI_t = ( q_e^2 .* EField_t_envelope .^2 ) / ( 2 * m_e * omega_d^2 ) .* sin(EPhase_t).^2;
   % electron ATI energy as a function of time (J) (row vector)
   
   n_re = n_0*0 + n_1*1 + n_2*2 + n_3*3 + n_4*4 + n_5*5;
   % relative electron density as a function of time n_e(t)
   
   dne_dt(2:N) = (n_re(2:N)-n_re(1:N-1))/DeltaTime;
   % variation of relative electron density as a function of time (dn_e/dt)
   
   Energy_ATI_e = sum(dne_dt .* Energy_ATI_t) / sum(dne_dt);  % (unit: J)
   % average ATI energy per electron
   
 % electron absorbed energy for ionization
   Energy_Ionization_e = n_1(N)*E_ion_0 + ...
                         n_2(N)*(E_ion_0 + E_ion_1) + ...
                         n_3(N)*(E_ion_0 + E_ion_1 + E_ion_2) + ...
                         n_4(N)*(E_ion_0 + E_ion_1 + E_ion_2 + E_ion_3) + ...
                         n_5(N)*(E_ion_0 + E_ion_1 + E_ion_2 + E_ion_3 + E_ion_4);
   % average ionization energy per electron
   
 % Return results
   result = [n_re(N), Energy_Ionization_e, Energy_ATI_e];         % unit: [arb. units, J]
 % Data plotting
   if FigureSwitch
      t = time ./ fs;
      figure;
      subplot(5,1,1), plot(t,EField_t),
                      ylabel('laser electric field (V/m)');
      subplot(5,1,2), plot(t,I_t .* 1E-4),
                      ylabel('laser intensity (W*cm^{-2})');
      subplot(5,1,3), plot(t,n_0,t,n_1,t,n_2,t,n_3,t,n_4,t,n_5),
                      ylabel('relative ion population');
      subplot(5,1,4), plot(t, n_re), ylabel('relative electron density n_e');
      subplot(5,1,5), plot(t, Energy_ATI_t / eV), ylabel('Energy_{ATI} (eV)'),
      xlabel('time (fs)');
 
      % Export results (usual units)
      result_population(1,:) = t ;         % position (mm)
      result_population(2,:) = n_0;        % position (mm)
      result_population(3,:) = n_1;        % position (mm)
      result_population(4,:) = n_2;        % position (mm)
      result_population(5,:) = n_3;        % position (mm)
      result_population(6,:) = n_4;        % position (mm)
      result_population(7,:) = n_5;        % position (mm)
      result_population(8,:) = n_re;        % position (mm)
      result_population(9,:) = Energy_ATI_t / eV;        % position (mm)
      result_population(10,:) = EField_t;        % position (mm)
      fid_output = fopen(['Population_',num2str(1),'.txt'],'w');       
      fprintf(fid_output,'t     n_0          n_1        n_2         n_3         n_4         n_5         n_re        Energy_ATI_t  EField_t \r\n');
      fprintf(fid_output,'(fs)  (arb.units) (arb.units) (arb.units) (arb.units) (arb.units) (arb.units) (arb.units) (eV)          (V/m)    \r\n');
      fprintf(fid_output,'%f    %f          %f          %f          %f          %f          %f          %f          %f            %f       \r\n', result_population);     
      fclose(fid_output);    
      %pause
   end
end


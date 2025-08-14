function result = StaticIonizationRate(E_ion, EField)

   global E_H
   global omega_a
   global EField_0

   a = E_ion / E_H;
   EField = EField + (EField==0)*10^-10;
   b = EField_0 ./ EField;
   
   result = 4 * omega_a * a^(5/2) * b .* exp( (-2/3) * a^(3/2) .* b);

end

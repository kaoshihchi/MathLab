function result = hhgDensityDistribution(N_density,w,rr,figureSwitch,j,N_z)
%HHGDENSITYDISTRIBUTION Gaussian fit of radial density profile.

   ft = fittype('A*exp(-(x.^2)/(B^2))','coefficients',{'A','B'},...
                'independent',{'x'});

   r = [0;w/2;w];
   N_density_cfit = fit(r,N_density,ft,'StartPoint',[N_density(1,1),w]);
   R_density = N_density_cfit.B;

   if and(figureSwitch,j==1)
       figure;
       plot(N_density_cfit,r,N_density);
       xlabel('r (m)'), ylabel('ion density (m^{-3})');
       legend('Location','NorthEast');
       savefig(gcf,'Fig_Fitted2+IonDensity_ini','compact');
   end

   if and(figureSwitch,j==N_z)
       figure;
       plot(N_density_cfit,r,N_density);
       xlabel('r (m)'), ylabel('ion density (m^{-3})');
       legend('Location','NorthEast');
       savefig(gcf,'Fig_Fitted2+IonDensity_final','compact');
   end

   result = {R_density,N_density_cfit(rr)};
end

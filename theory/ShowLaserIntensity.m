function ShowLaserIntensity(intensity_Wcm2)
%SHOWLASERINTENSITY Display peak intensity map in W/cm^2.
%
%   SHOWLASERINTENSITY(intensity_Wcm2) displays the peak intensity map and
%   annotates the peak value.

    if isempty(intensity_Wcm2)
        error('Intensity map not provided.');
    end
    figure;
    imagesc(intensity_Wcm2);
    axis image; colormap hot;
    cb = colorbar; cb.Label.String = 'Peak Intensity (W/cm^2)';
    title('Laser Peak Intensity (Gaussian temporal, background-subtracted)');
    xlabel('X pixel'); ylabel('Y pixel');

    Ipk = max(intensity_Wcm2(:));
    text(0.02, 0.98, sprintf('Peak: %.3g W/cm^2', Ipk), ...
         'Units','normalized','Color','w','FontWeight','bold', ...
         'HorizontalAlignment','left','VerticalAlignment','top');
end

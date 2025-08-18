function [intensity_Wm2, intensity_Wcm2] = ComputePeakIntensity(correctedImage, pixelSize_m, pulseEnergy_J, tauFWHM_s)
%COMPUTEPEAKINTENSITY Convert corrected counts to peak intensity maps.
%   [intensity_Wm2, intensity_Wcm2] = COMPUTEPEAKINTENSITY(correctedImage, pixelSize_m, pulseEnergy_J, tauFWHM_s)
%   returns the peak intensity in W/m^2 and W/cm^2 for each pixel assuming a
%   Gaussian temporal profile.

    if isempty(correctedImage)
        error('No corrected image provided.');
    end

    totalCounts = sum(correctedImage(:));
    if totalCounts <= 0
        error('Total counts are zero after background subtraction.');
    end
    energyPerCount = pulseEnergy_J / totalCounts;      % J / count
    E_pixel = correctedImage * energyPerCount;         % J per pixel

    A_pixel = pixelSize_m^2;                           % m^2
    gaussianFactor = sqrt(4*log(2)/pi);
    intensity_Wm2 = (E_pixel ./ A_pixel) .* (gaussianFactor / tauFWHM_s);
    intensity_Wcm2 = intensity_Wm2 / 1e4;              % W/cm^2
end

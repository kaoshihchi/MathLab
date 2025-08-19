function ShowLaserImage(correctedImage)
%SHOWLASERIMAGE Display background-subtracted counts.
%
%   SHOWLASERIMAGE(correctedImage) displays the corrected image using a hot
%   colormap.

    if isempty(correctedImage)
        error('Corrected image not available.');
    end
    figure;
    imagesc(correctedImage);
    axis image; colormap hot; colorbar;
    title('Laser Profile (Background Subtracted Counts)');
    xlabel('X pixel'); ylabel('Y pixel');
end

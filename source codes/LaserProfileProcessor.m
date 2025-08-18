classdef LaserProfileProcessor
    properties
        inputFile           % Filename of input tiff
        backgroundFile      % Filename of background tiff
        inputImage          % raw (double)
        backgroundImage     % raw (double)
        correctedImage      % background-subtracted counts (>=0)
        intensity_Wm2       % peak intensity map (W/m^2)
        intensity_Wcm2      % peak intensity map (W/cm^2)
    end
    
    methods
        % --- Constructor ---
        function obj = LaserProfileProcessor(inputFile, backgroundFile)
            obj.inputFile = inputFile;
            obj.backgroundFile = backgroundFile;
        end
        
        % --- Load Images ---
        function obj = loadImages(obj)
            assert(isfile(obj.inputFile),    "File not found: %s", obj.inputFile);
            assert(isfile(obj.backgroundFile),"File not found: %s", obj.backgroundFile);
            obj.inputImage      = double(imread(obj.inputFile));
            obj.backgroundImage = double(imread(obj.backgroundFile));
        end
        
        % --- Subtract Background ---
        function obj = subtractBackground(obj)
            if isempty(obj.inputImage) || isempty(obj.backgroundImage)
                error('Images not loaded. Run loadImages() first.');
            end
            obj.correctedImage = obj.inputImage - obj.backgroundImage;
            obj.correctedImage(obj.correctedImage < 0) = 0; % clip negatives
        end
        
        % --- Compute Peak Intensity from energy & pulse width ---
        % pixelSize_m:   pixel pitch in meters (e.g., 0.88e-6)
        % pulseEnergy_J: pulse energy in Joules (e.g., 8.3e-3)
        % tauFWHM_s:     FWHM pulse duration in seconds (e.g., 40e-15)
        function obj = computeIntensity(obj, pixelSize_m, pulseEnergy_J, tauFWHM_s)
            if isempty(obj.correctedImage)
                error('No corrected image. Run subtractBackground() first.');
            end
            
            % --- Energy normalization (counts -> Joules) ---
            totalCounts = sum(obj.correctedImage(:));
            if totalCounts <= 0
                error('Total counts are zero after background subtraction.');
            end
            energyPerCount = pulseEnergy_J / totalCounts;      % J / count
            E_pixel = obj.correctedImage * energyPerCount;     % J per pixel
            
            % --- Convert to peak intensity for Gaussian temporal shape ---
            % I0 = (E_pixel / A_pixel) * sqrt(4 ln 2 / pi) / tauFWHM
            A_pixel = pixelSize_m^2;                           % m^2
            gaussianFactor = sqrt(4*log(2)/pi);
            obj.intensity_Wm2 = (E_pixel ./ A_pixel) .* (gaussianFactor / tauFWHM_s);
            obj.intensity_Wcm2 = obj.intensity_Wm2 / 1e4;      % W/cm^2
        end
        
        % --- Show background-subtracted counts ---
        function showImage(obj)
            if isempty(obj.correctedImage)
                error('Corrected image not available. Run subtractBackground() first.');
            end
            figure;
            imagesc(obj.correctedImage);
            axis image; colormap hot; colorbar;
            title('Laser Profile (Background Subtracted Counts)');
            xlabel('X pixel'); ylabel('Y pixel');
        end
        
        % --- Show intensity (W/cm^2) ---
        function showIntensity(obj)
            if isempty(obj.intensity_Wcm2)
                error('Intensity not computed. Run computeIntensity() first.');
            end
            figure;
            imagesc(obj.intensity_Wcm2);
            axis image; colormap hot; 
            cb = colorbar; cb.Label.String = 'Peak Intensity (W/cm^2)';
            title('Laser Peak Intensity (Gaussian temporal, background-subtracted)');
            xlabel('X pixel'); ylabel('Y pixel');
            
            % annotate peak value
            Ipk = max(obj.intensity_Wcm2(:));
            text(0.02, 0.98, sprintf('Peak: %.3g W/cm^2', Ipk), ...
                'Units','normalized','Color','w','FontWeight','bold', ...
                'HorizontalAlignment','left','VerticalAlignment','top');
        end
    end
end

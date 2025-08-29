classdef SpectrumGenerator
    % SpectrumGenerator
    %   Super-Gaussian family:
    %   I(λ) = A * exp( - |(λ - λ0)/C|^(2n) ) + D
    %   where C = FWHM / ( 2 * (ln 2)^(1/(2n)) )
    %
    %   Set order_n = 1 for Gaussian; larger n → flatter top.

    properties
        center_nm (1,1) double = 815
        fwhm_nm   (1,1) double = 35.8
        order_n   (1,1) double = 1
        amplitude (1,1) double = 1.0
        offset    (1,1) double = 0.0
        noise_std (1,1) double = 0.0
    end

    methods
        function obj = SpectrumGenerator(center_nm, fwhm_nm, order_n, amplitude, offset, noise_std)
            if nargin >= 1, obj.center_nm = center_nm; end
            if nargin >= 2, obj.fwhm_nm   = fwhm_nm;   end
            if nargin >= 3, obj.order_n   = order_n;   end
            if nargin >= 4, obj.amplitude = amplitude; end
            if nargin >= 5, obj.offset    = offset;    end
            if nargin >= 6, obj.noise_std = noise_std; end
        end

        function C = widthParam(obj)
            % Scale parameter from FWHM for super-Gaussian order n
            C = obj.fwhm_nm / ( 2 * (log(2))^(1/(2*obj.order_n)) );
        end

        function I = evaluate(obj, lambda_nm)
            % Evaluate intensity at given wavelengths [nm]
            C = obj.widthParam();
            I = obj.amplitude .* exp( -abs((lambda_nm - obj.center_nm)./C).^(2*obj.order_n) ) + obj.offset;

            if obj.noise_std > 0
                I = I + obj.noise_std .* randn(size(lambda_nm));
            end
        end
    end
end

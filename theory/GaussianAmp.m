function y = GaussianAmp(x, A, B, C, D)
%GAUSSIANAMP Gaussian amplitude function.
%   y = GAUSSIANAMP(x, A, B, C, D) computes the value of a Gaussian
%   amplitude profile defined by
%       y = A * exp(-2 * (abs((x - B) / C).^2)) + D.
%
%   Inputs:
%       x - Evaluation points.
%       A - Amplitude of the peak.
%       B - Center position.
%       C - Width parameter.
%       D - Baseline offset.
%
%   Output:
%       y - Calculated intensity.

    y = A * exp(-2 * (abs((x - B) / C).^2)) + D;
end

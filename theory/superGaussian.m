function y = superGaussian(x, A, B, C, n, D)
%SUPERGAUSSIAN Super-Gaussian function.
%   y = SUPERGAUSSIAN(x, A, B, C, n, D) computes the value of a
%   super-Gaussian profile defined by
%       y = A * exp(-(abs((x - B) / C).^(2*n))) + D.
%
%   Inputs:
%       x - Evaluation points.
%       A - Amplitude of the peak.
%       B - Center position.
%       C - Width parameter.
%       n - Super-Gaussian order.
%       D - Baseline offset.
%
%   Output:
%       y - Calculated intensity.

    y = A * exp(-(abs((x - B) / C).^(2*n))) + D;
end

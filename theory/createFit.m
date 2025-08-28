function [fitresult, gof] = createFit(t_shifted_filtered, I_t_shifted_filtered, optimOptions)
%CREATEFIT Fit Gaussian amplitude to shifted intensity data.
%   [FITRESULT, GOF] = CREATEFIT(T_SHIFTED_FILTERED, I_T_SHIFTED_FILTERED, OPTIONS)
%   fits the data in T_SHIFTED_FILTERED and I_T_SHIFTED_FILTERED using a
%   Gaussian amplitude model. OPTIONS is a struct containing optimization
%   limits such as MaxIter, MaxFunEvals and MaxTime. The implementation
%   uses streamlined data prep and relaxed tolerances for faster
%   convergence.
%
%   Example:
%       opts = struct('MaxIter',100,'MaxFunEvals',500,'MaxTime',120);
%       [fitresult, gof] = createFit(t, y, opts);

% Prepare data for fitting without helper utilities
xData = t_shifted_filtered(:);
yData = I_t_shifted_filtered(:);

% Set up fittype and default options
ft = fittype('A*exp(-2*((x-B)/C).^2)+D', 'independent', 'x', 'dependent', 'y');
opts = fitoptions('Method', 'NonlinearLeastSquares');
opts.Display = 'Off';
opts.Lower = [0 -1e-15 10e-15 -0.1];
opts.Upper = [2 1e-15 50e-15 0.1];
opts.StartPoint = [1.2 0 20e-15 0];
opts.TolFun = 1e-6;
opts.TolX = 1e-6;

% Apply shared optimization limits
opts.MaxIter = optimOptions.MaxIter;
opts.MaxFunEvals = optimOptions.MaxFunEvals;
if isfield(optimOptions, 'MaxTime')
    opts.MaxTime = optimOptions.MaxTime;
end

% Fit model to data
[fitresult, gof] = fit(xData, yData, ft, opts);
end

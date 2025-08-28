function fitFcn = timed_lsqcurvefit(x, y, fitType, varargin)
%TIMED_LSQCURVEFIT Fit data with runtime limits.
%   FITFCN = TIMED_LSQCURVEFIT(X, Y, FITTYPE) fits the data (X,Y) using
%   LSQCURVEFIT with a model specified by FITTYPE. The function returns a
%   function handle FITFCN that evaluates the fitted model at new query
%   points. Additional name-value pairs control the solver:
%       'MaxIterations'          - maximum number of iterations
%       'MaxFunctionEvaluations' - maximum number of function evaluations
%       'MaxTime'                - maximum runtime in seconds
%       'Downsample'             - use every N-th point from the data
%
%   Supported FITTYPE strings: 'exp1', 'poly3', 'poly5', 'power2',
%   and 'fourier3'.
%
%   Example:
%       f = timed_lsqcurvefit(x, y, 'poly3', 'MaxTime', 60);
%       yfit = f(x);
%
%   This utility requires the Optimization Toolbox.

p = inputParser;
addParameter(p, 'MaxIterations', 400);
addParameter(p, 'MaxFunctionEvaluations', 800);
addParameter(p, 'MaxTime', 120);
addParameter(p, 'Downsample', 1);
parse(p, varargin{:});

x = x(:);  % ensure column
y = y(:);

% Downsample data if requested
if p.Results.Downsample > 1
    x = x(1:p.Results.Downsample:end);
    y = y(1:p.Results.Downsample:end);
end

switch fitType
    case 'exp1'
        model = @(c, xdata) c(1) .* exp(c(2) .* xdata);
        c0 = [y(1); 0];
    case 'poly3'
        model = @(c, xdata) polyval(c, xdata);
        c0 = zeros(4,1);
    case 'poly5'
        model = @(c, xdata) polyval(c, xdata);
        c0 = zeros(6,1);
    case 'power2'
        model = @(c, xdata) c(1) .* xdata.^c(2) + c(3);
        c0 = [1; 1; 0];
    case 'fourier3'
        model = @(c, xdata) c(1) + ...
            c(2).*cos(xdata.*c(8)) + c(3).*sin(xdata.*c(8)) + ...
            c(4).*cos(2*xdata.*c(8)) + c(5).*sin(2*xdata.*c(8)) + ...
            c(6).*cos(3*xdata.*c(8)) + c(7).*sin(3*xdata.*c(8));
        c0 = zeros(8,1); c0(8) = 1;
    otherwise
        error('timed_lsqcurvefit:UnsupportedType', 'Unsupported fit type %s', fitType);
end

options = optimoptions('lsqcurvefit', ...
    'MaxIterations', p.Results.MaxIterations, ...
    'MaxFunctionEvaluations', p.Results.MaxFunctionEvaluations, ...
    'MaxTime', p.Results.MaxTime);

coeff = lsqcurvefit(model, c0, x, y, [], [], options);
fitFcn = @(xquery) model(coeff, xquery);
end


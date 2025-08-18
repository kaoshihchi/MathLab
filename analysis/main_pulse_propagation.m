function result = main_pulse_propagation()
% Load parameters from CSV file
% Table format: Parameter,Value,OnOff
UseParametersFile = true;
if UseParametersFile
    params = readtable('analysis/Parameters.csv');
    for idx = 1:height(params)
        if params.OnOff(idx)
            assignin('base', params.Parameter{idx}, params.Value(idx));
        end
    end
end

% Placeholder for analysis code
result = [];
end

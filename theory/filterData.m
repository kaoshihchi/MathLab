function [I_t_shifted_filtered, t_shifted_filtered] = filterData(I_t_shifted, t_shifted, t_bound)
%FILTERDATA Remove data within a symmetric time bound.
%   [I_F, T_F] = FILTERDATA(I_T_SHIFTED, T_SHIFTED, T_BOUND) removes
%   elements of the arrays whose time values fall within [-T_BOUND, T_BOUND].
%
%   Inputs:
%       I_t_shifted - Intensity array.
%       t_shifted   - Time array (seconds).
%       t_bound     - Bound for exclusion (seconds).
%
%   Outputs:
%       I_t_shifted_filtered - Filtered intensity array.
%       t_shifted_filtered   - Corresponding time array.

    time_lower_bound = -t_bound;
    time_upper_bound = t_bound;

    mask = t_shifted < time_lower_bound | t_shifted > time_upper_bound;

    I_t_shifted_filtered = I_t_shifted(mask);
    t_shifted_filtered = t_shifted(mask);
end

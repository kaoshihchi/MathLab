function [t_selected, I_selected] = selectDataInRange(t_array, I_array, lowerBound, upperBound)
%SELECTDATAINRANGE Select data within specified bounds.
%   [T_SELECTED, I_SELECTED] = SELECTDATAINRANGE(T_ARRAY, I_ARRAY,
%   LOWERBOUND, UPPERBOUND) returns elements of the arrays where T_ARRAY
%   is between LOWERBOUND and UPPERBOUND.

    inRangeIndices = t_array >= lowerBound & t_array <= upperBound;
    t_selected = t_array(inRangeIndices);
    I_selected = I_array(inRangeIndices);
end

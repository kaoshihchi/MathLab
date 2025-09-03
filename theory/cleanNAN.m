function [x_clean, y_clean] = cleanNAN(x, y)
    validMask = isfinite(x) & isfinite(y);
    x_clean = x(validMask);
    y_clean = y(validMask);
end
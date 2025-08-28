function saveIFFTData(array1, array2, filename)
%SAVEIFFTDATA Save time-domain IFFT data to a text file.
%   SAVEIFFTDATA(array1, array2, filename) writes two-column data with
%   headers 'Time' (ps) and 'Intensity' (arb. units).
%
%   Inputs:
%       array1   - Time array (ps).
%       array2   - Intensity array (arb. units).
%       filename - Output file name. '.txt' is appended if missing.

    combinedData = [array1, array2];
    mainHeaders = {'Time', 'Intensity'};
    units = {'ps', 'arb. units'};

    if ~endsWith(filename, '.txt')
        filename = [filename, '.txt'];
    end

    fid = fopen(filename, 'w');
    fprintf(fid, '%s,%s\n', mainHeaders{:});
    fprintf(fid, '%s,%s\n', units{:});
    fclose(fid);

    writematrix(combinedData, filename, 'WriteMode', 'append');
end

function saveIFFTPulse(array1, array2, filename)
%SAVEIFFTPULSE Save IFFT pulse data to a text file.
%   SAVEIFFTPULSE(array1, array2, filename) writes two-column data with
%   headers 'Time' (fs) and 'Intensity' (arb. units).
%
%   Inputs:
%       array1   - Time array (fs).
%       array2   - Intensity array (arb. units).
%       filename - Output file name. '.txt' is appended if missing.

    combinedData = [array1, array2];
    mainHeaders = {'Time', 'Intensity'};
    units = {'fs', 'arb. units'};

    if ~endsWith(filename, '.txt')
        filename = [filename, '.txt'];
    end

    fid = fopen(filename, 'w');
    fprintf(fid, '%s,%s\n', mainHeaders{:});
    fprintf(fid, '%s,%s\n', units{:});
    fclose(fid);

    writematrix(combinedData, filename, 'WriteMode', 'append');
end

function saveSpectrumData(array1, array2, filename)
%SAVESPECTRUMDATA Save wavelength and intensity arrays to a text file.
%   SAVESPECTRUMDATA(array1, array2, filename) writes two-column data
%   with headers 'Wavelength' (nm) and 'Intensity' (arb. units).
%
%   Inputs:
%       array1   - Wavelength array (nm).
%       array2   - Intensity array (arb. units).
%       filename - Output file name. '.txt' is appended if missing.

    combinedData = [array1, array2];
    mainHeaders = {'Wavelength', 'Intensity'};
    units = {'nm', 'arb. units'};

    if ~endsWith(filename, '.txt')
        filename = [filename, '.txt'];
    end

    fid = fopen(filename, 'w');
    fprintf(fid, '%s,%s\n', mainHeaders{:});
    fprintf(fid, '%s,%s\n', units{:});
    fclose(fid);

    writematrix(combinedData, filename, 'WriteMode', 'append');
end

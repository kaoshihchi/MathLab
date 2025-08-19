function [inputImage, backgroundImage] = LoadLaserImages(inputFile, backgroundFile)
%LOADLASERIMAGES Read input and background TIFF images as doubles.
%   [inputImage, backgroundImage] = LOADLASERIMAGES(inputFile, backgroundFile)
%   reads the specified TIFF files and returns them as double precision
%   matrices.

    assert(isfile(inputFile), "File not found: %s", inputFile);
    assert(isfile(backgroundFile), "File not found: %s", backgroundFile);
    inputImage      = double(imread(inputFile));
    backgroundImage = double(imread(backgroundFile));
end

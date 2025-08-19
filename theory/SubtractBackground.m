function correctedImage = SubtractBackground(inputImage, backgroundImage)
%SUBTRACTBACKGROUND Subtract background from input image and clip negatives.
%   correctedImage = SUBTRACTBACKGROUND(inputImage, backgroundImage)
%   subtracts the backgroundImage from inputImage and sets negative values to
%   zero.

    if isempty(inputImage) || isempty(backgroundImage)
        error('Input and background images must be provided.');
    end
    correctedImage = inputImage - backgroundImage;
    correctedImage(correctedImage < 0) = 0;
end

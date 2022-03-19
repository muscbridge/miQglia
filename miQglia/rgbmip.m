function outputImg = rgbmip(inputImg)
%% Maximum Intensity Proection (MIP) of RGB Image Stack
% This function open an RGB image, splits it into red, green and blue
% channels, then computes the MIP. The resulting output is a 2D image.
% Parameters
% ----------
%   inputImg : 3D array containing input image
%
% Returns
% -------
%   outputImg : 2D array of returned MIP
% Author: Siddhartha Dhiman

    [rows, columns, numSlices] = size(inputImg);
    if numSlices <= 1 | numSlices > 3
        error('Error. \nInput image must be a 3-dimensional array, input image contains %.0f dimensions.', numSlices);
    end
    outputImg = zeros(rows, columns, class(inputImg));
    for col = 1 : columns
        for row = 1 : rows
            zSlice = inputImg(row, col, :);
            maxValue = max(zSlice);
            outputImg(row, col) = maxValue;
        end
    end
end


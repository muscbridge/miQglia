function outputImg = bwclearborder(inputImg)
%% Remove Objects at Image Borders
% This function removes objects touching the border of a logical image
% without cropping.
% ----------
%   inputImg : 2D logical array containing input image
%
% Returns
% -------
%   outputImg : 2D logical array same size an input
% Author: Siddhartha Dhiman

    [rows, columns, numSlices] = size(inputImg);
    if numSlices ~= 1
        error('Error. \nInput image must be a 2-D image, input image contains %.0f dimensions.', numSlices);
    end
    if ~isa(inputImg, 'logical')
        error('Input image needs to be a "logical" array, user entered "%s" array', class(inputImg));
    end
    bw = imclearborder(inputImg);
    [x, y] = size(bw);
    marginX = (rows - x) / 2;
    marginY = (columns - y) / 2;
    if (marginX > 0 & marginY > 0)
        outputImg = logical(padarray(bw, [marginX, marginY], 'direction', 'both', 'padval', 0));
    else
        outputImg = bw;
    end
    


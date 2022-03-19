function outputImg = fft_bandpass(inputImg, l , u, N)
%% FFT Bandpass Filter for Grayscale (2D) Images
% This function applies an FFT bandpass filter to 2D grayscale images to
% remove noise (small features) while preserving structures in the image
% Parameters
% ----------
%   inputImg : 2D array containing grayscale input image
%
% Returns
% -------
%   outputImg : 2D array of filtered image
% Author: Siddhartha Dhiman

    [rows, columns, numSlices] = size(inputImg);
    if numSlices ~= 1
        error('Error. \nInput image must be a 1-dimensional array, input image contains %.0f dimensions.', numSlices);
    end
    imgType = class(inputImg);
    % 2-D fast Fourier transform 
    nx = rows;
    ny = columns;
    fft_img = fft2(inputImg, 2*nx-1, 2*ny-1);
    fft_img = fftshift(fft_img);
    % Initialize filter.
    filter1 = ones(2*nx-1,2*ny-1);
    filter2 = ones(2*nx-1,2*ny-1);
    filter3 = ones(2*nx-1,2*ny-1);
    for i = 1:2*nx-1
        for j =1:2*ny-1
            dist = ((i-(nx+1))^2 + (j-(ny+1))^2)^.5;
            % Create Butterworth filter.
            filter1(i,j)= 1/(1 + (dist/u)^(2*N));
            filter2(i,j) = 1/(1 + (dist/l)^(2*N));
            filter3(i,j)= 1.0 - filter2(i,j);
            filter3(i,j) = filter1(i,j).*filter3(i,j);
        end
    end
    % Apply filter to image
    outputImg = fft_img + filter3.*fft_img;
    % Inverse Fourier tranform to return image
    outputImg = ifftshift(outputImg);
    outputImg = ifft2(outputImg, 2*nx-1, 2*ny-1);
    outputImg = real(outputImg(1:nx,1:ny));
    outputImg = cast(outputImg, imgType);
end

# miQglia
A set of MATLAB scripts for quantification of microglia images acquired through
brightfield microscopy.

![Alt Text](docs/readme_animation.gif)

## Image Processing
Input images are preprocessed to enhace SNR for better cell segmentation. The
current preprocessing steps are:

1. Conversion to maximum intensity projection (MIP)
2. Denoising with 2D Weiner filter
3. Background feature removal with a top-hat filter
4. FFT bandpass filter
5. Image sharpening

The culmination of image processing in this order yielded the best segmentation
results.

## Usage
Call the function `ramificationStats()` to compute process and quantify an image
containing microglia.
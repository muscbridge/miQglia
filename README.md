# miQglia
A set of MATLAB scripts for quantification of microglia images acquired through
brightfield microscopy.

>Image processing and enhancement

<img src="docs/readme_animation.gif" alt="Image processing" width=512/>

> Marker-controller watershed segmentation of microglia

<img src="docs/I_min_Contour.jpg" alt="Watershed contour" width=512/>
<img src="docs/Topoography_Cell_1.jpg" alt="Watershed topography" width=512/>

> Fractal analysis

<img src="docs/soma_10_boxes.gif" alt="Fractal dimension" width=512/>

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

## Watershed Image Segmentation
> Steps in marker-controlled watershed segmentation

<img src="docs/Watershed_Segmentation.png" alt="Watershed Segmentation" width=512/>

The processed image is subtly blurred with a Gaussian filter, using the MATLAB
function `imgaussfilt`. A microglia mask (MM) is then created by thresholding
the Gaussian smoothed image, which is also subsequently used to create a soma
seed mask (SSM) by applying a higher threshold to demarcate bright regions
likely to be soma. Regions in SMM were shrunk to point-pixel and regions outside
of MM were removed to identify and index soma locations (SL). Then, a composite
image is created for marker-controlled watershed segmentation by imposing SL
onto inverted gaussian-smoothed as regional minima. The composite image
undergoes watershed segmentation to identify microglia. Microglia touching the
border or not overlapping with SL are first removed, then eroded by two pixels
to clean up watershed segmentation results. Microglia are then quantified by
metrics in the table below and any microglia with exclusion criteria are
rejected and removed.

## Microglia Quantification
### Skeleton Analysis
Metrics are based on skeletonized microglia to capture their complicated network
process

<dl>
  <dt><strong>1. BranchPoints</strong></dt>
  <dd>Number of branch points in microglia</dd>
  <dt><strong>2. EndPoints</strong></dt>
  <dd>Length of largest process on microglia</dd>
  <dt><strong>3. Branches</strong></dt>
  <dd>Number of processes in microglia</dd>
  <dt><strong>4. MaxBranchLength</strong></dt>
  <dd>Length of largest process on microglia</dd>
  <dt><strong>5. AverageBranchLength</strong></dt>
  <dd>Average process length on microglia</dd>
  <dt><strong>6. TotalBranchLength</strong></dt>
  <dd>Summation of all processes lengths on microglia</dd>
</dl>

### Fractal Analysis
A set of measures to quantify complexity of microglia structure.

<dl>
  <dt><strong>1. FractalDimension</strong></dt>
  <dd>Ratio quantifying cell complexity</dd>
  <dt><strong>2. Circularity</strong></dt>
  <dd>Roundness of microglia. Ranges from 0 (perfect circle) to 1 (ellipsoid)</dd>
  <dt><strong>3. Span Ratio</strong></dt>
  <dd>Measures microglia elongation</dd>
  <dt><strong>4. Density</strong></dt>
  <dd>Quantifies of microglia size</dd>
</dl>

## Usage
Change the variable `pixel_micron` in L31 of _ramificationStats.m_ to your pixel
to micron ratio. Then, call the function `ramificationStats('Path to input image', 'Path to output directory')` to compute, process, and quantify an image containing microglia.

Example:
```
>> cell_stats = ramificationStats("D:\Datasets\Morphology\data\Benchmark\Input\w1t5m18_s20_40x_Sub2.TIF",'D:\SystemFiles\siddh\\Projects\\Mouse_Histology\Output');
```

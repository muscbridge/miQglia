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

function props = ramificationStatsOld(impath, outpath)
%% Skeleton and Fractal Analysis of Micrographs
% A function designed to extract skeleton and fractal parameters from RGB
% mictographs of ramified cells
%
% Author: Siddhartha Dhiman
%
% Parameters
% ----------
%   impath : str
%       Path to input image file (.tif)
%   outpath : str
%       Path to output directory write out files
% Returns
% -------
%   props : struct
%       Struct containing quantified parameters
%--------------------------------------------------------------------------

    %% Begin
    if ~isfile(impath)
        error('Input file %s does not exist', impath);
    end
    if ~isdir(outpath)
        error('Output path %s is not a directory. Please specify a valid directory.', outpath);
    end
    if ~exist(outpath, 'dir')
        error('Output directory %s does not exist.', outpath)
    end

    %% Main
    [fpath,fname,fext] = fileparts(impath);
    img = imread(impath);
    % Maximum intensity projection to convert to grayscale
    mip = rgbmip(img);
    % CLAHE
    eqI = adapthisteq(mip);
    % FFT bandpass filter to remove noise
    fftI = fft_bandpass(eqI, 10, 100, 1);
    % Sharpen image
    shI = imsharpen(fftI, 'Radius', 3, 'Amount', 0.6);
    % Denoise image with 2D Weiner filter
    dnI = wiener2(shI, [10,10]);
    % Median filter
    mdI = medfilt2(dnI, [5,5]);

    % Binarize image using Otsu's thresholding
    bwI = ~imbinarize(mdI, 'global');
    % Remove objects at border
    bwI_bod = imclearborder(bwI);
    % Median filter image again to remove specks
    bwI_med = medfilt2(bwI_bod, [3,3]);
    % Morphologically close image
    bwI_close = imclose(bwI_med, strel('disk', 5));
    % Remove small objects
    bwI_open = bwareaopen(bwI_close, 1000);
    % Skeletonize image
    skel = bwskel(bwI_open);
    % Find objects in image
    props = regionprops(skel, 'Centroid', 'Area', 'BoundingBox', 'Image');
    % Populate `props` with skeleton quantification
    for i = 1:size(props, 1)
        props(i).Index = i;
        skelI = props(i).Image;
        props(i).BranchPoints = nnz(bwmorph(skelI, 'branchpoints'));
        props(i).EndPoints = nnz(bwmorph(skelI, 'endpoints'));
        skel_props = regionprops(logical(skelI - bwmorph(skelI, 'branchpoints')), 'Area', 'Image');
        props(i).Branches = size(skel_props, 1);
        props(i).MaxBranchPoint = max(extractfield(skel_props, 'Area'));
        props(i).AverageBranchPoint = mean(extractfield(skel_props, 'Area'));
        props(i).TotalBranchPoint = sum(extractfield(skel_props, 'Area'));
        % Save each soma
        x = fix(props(i).BoundingBox(1)- 25);
        y = fix(props(i).BoundingBox(2) - 25);
        bx = fix(props(i).BoundingBox(3) + 50);
        by = fix(props(i).BoundingBox(4) + 50);
        props(i).RawImage = imcrop(mip, [x, y, bx, by]);
        img_crop = imcrop(bwI_open, [x, y, bx, by]);
        img_crop = imclearborder(img_crop);
    %     img_crop = imfill(img_crop, 'holes'); % fill holes
        img_crop = bwareafilt(img_crop, 1); % keep the biggest object in bounding box
    %     img_crop = bwperim(img_crop);
        props(i).PerimeterImage = img_crop;
        img_props = regionprops(img_crop, 'Area', 'Circularity', 'ConvexArea', 'ConvexImage', 'Image', 'MajorAxisLength', 'MinorAxisLength');
        props(i).FractalDimension = hausDim(img_crop);
        props(i).Circularity = img_props(1).Circularity;
        props(i).SpanRatio = img_props(1).MajorAxisLength /  img_props(1).MinorAxisLength;
        props(i).Density = nnz(bwperim(img_crop)) / img_props(1).ConvexArea;
        fname_crop = fullfile(outpath, strcat('S', num2str(i), '_', fname, '.tif'));
%         imwrite(img_crop, fname_crop);
    end
    % Filter results based on thresholds
    rmIdx = extractfield(props, 'EndPoints') == 2 & extractfield(props, 'MaxBranchPoint') < 50;
    props = props(~rmIdx);
    totalEndPoints = sum(extractfield(props, 'EndPoints'));
    rmIdx = extractfield(props, 'TotalBranchPoint') < 50;
    props = props(~rmIdx);
    totalBranchPoint = sum(extractfield(props, 'TotalBranchPoint'));
    
    props = rmfield(props, {'Image', 'PerimeterImage', 'RawImage'});
    % Write properties
    fname_table = fullfile(outpath, strcat('Stats_', fname, '.csv'));
    table_header = {'Index', 'Centroid', 'BoundingBox', 'Area', 'BranchPoints', 'EndPoints', 'Branches', 'MaxBranchPoint', 'AverageBranchPoint', 'TotalBranchPoint', 'FractalDimension', 'Circularity', 'SpanRatio', 'Density'};
    props = orderfields(props, table_header);
    writetable(struct2table(props), fname_table);
 
    % Make nice figure for QC
%     fout = fullfile(outpath, strcat('Processing_Steps_', fname, '.png'));
%     fig1 = figure('visible', 'off', 'WindowState','maximized');
%     t = tiledlayout(3, 9, 'TileSpacing', 'none', 'Padding', 'tight');
%     nexttile([3,3]); imshow(img); title(sprintf('INPUT: %s', replace(fname, '_', ' ')));
%     nexttile; imshow(mip); title('1. Maximum Intensity Projection');
%     nexttile; imshow(eqI); title('2. CLAHE');
%     nexttile; imshow(fftI); title('3. FFT Bandpass Filter');
%     nexttile([3,3]); imshow(labeloverlay(labeloverlay(img, bwI_open, 'Colormap', [0 1 0], 'Transparency', 0.8), skel, 'Colormap', [1 0 0], 'Transparency', 0));
%     title('Detected Soma and Skeleton Overlay');
%     nexttile; imshow(shI); title('4. Unsharpen Mask');
%     nexttile; imshow(dnI); title('5. Denoising');
%     nexttile; imshow(mdI); title('6. Median Filter');
%     nexttile; imshow(bwI_bod); title('7. Remove Border Objects');
%     nexttile; imshow(bwI_med); title('8. Median Filter');
%     nexttile; imshow(bwI_close); title('9. Morphological Closing');
%     set(gcf, 'Position', get(0, 'Screensize'));
%     saveas(fig1, fout, 'png');
    
    fout = fullfile(outpath, strcat('Labels', '_', fname, '.png'));
    fig2 = figure('visible', 'off', 'WindowState','maximized');
    imshow(labeloverlay(labeloverlay(img, bwI_open, 'Colormap', [0 1 0], 'Transparency', 0.8), skel, 'Colormap', [1 0 0], 'Transparency', 0));
    hold on;
    for i = 1:size(props, 1)
        rectangle('Position', props(i).BoundingBox, 'EdgeColor', 'r',...
            'LineWidth', 0.5, 'LineStyle', ':');
        text(props(i).BoundingBox(1)+5, props(i).BoundingBox(2)+10,...
            num2str(props(i).Index), 'FontSize', 12);
    end
    hold off;
    title(strcat(fname, ': Detected Somas and Labels'), 'Interpreter', 'none');
    saveas(fig2, fout, 'png');
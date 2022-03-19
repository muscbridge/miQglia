function props = ramificationStats(impath, outpath)
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
    %% Define Fixed Variables
    pixel_micron = 0.25;    % Number of microns per pixel in XY. I measured 40 pixels for 10 um (microns)
    tophat_filter = 2;      % Size of top-hat filter. Ideally the spatial size of dendrites (in microns)
    obj_remove = 50;        % Remove features same size or under (in microns)
    min_branch = 70;        % Minimum length of total branch distance (microns)

    %% Internal Pixel <--> Pixel Conversion (DO NOT ALTER!)
    obj_remove_pixel = round(obj_remove / pixel_micron, 0);
    pixel_micron_diag = sqrt(2* pixel_micron^2); % Diagonal distance of voxel
    pixel_micron_area =  pixel_micron^2; % Area of one pixel
    nhood1 = 10; % Neighborhood size of disk-shaped structuring element for initial erosion
    nhood2 = 3;  % Neighborhood size of disk-shaped structuring element closing
    nhood_th = round(tophat_filter / pixel_micron);
    %% Main
    [fpath,fname,fext] = fileparts(impath);
    img = imread(impath);
    % Maximum intensity projection to convert to grayscale
    mip = rgbmip(img);
    mip_ = uint8(max(mip(:))) - mip;     % inverted image;
    % Denoise image with 2D Weiner filter
    dnI = wiener2(mip_, [10,10]);
    % Remove redundant background
    thI = imtophat(dnI, strel('disk', nhood_th));
    % FFT bandpass filter to remove noise
    fftI = fft_bandpass(dnI, 10, 100, 1);       % GOOD FOR SOMA DETECTION
    % Sharpen image
    shI = imsharpen(fftI, 'Radius', 3, 'Amount', 0.6);

    % Median filter
%     mdI = medfilt2(dnI, [5,5]);
    % Invert final image so cells are bright and background is dark
%     I = uint8(255) - mdI;
    % Compute gradient magnitude as segmentation function. The gradient is
    % high at the borders of the objects and low (mostly) inside the 
    % objects. This is used as a function to direct watershed segmentation.
%     gmag = imgradient(I);
    gmag = imgradient(shI);
    % Morphologically reconstruct image through a series of erosions and
    % dilationsso that cellular features are flat i.e. minimal intensity
    % variations
    se = strel('disk', nhood1);
    Io = imopen(shI,se);  % Open image
    Ie = imerode(shI,se); % Erode opened image
    Iobr = imreconstruct(Ie,shI); % Morphologically reconstruct prior opening-closing onto original iamge
    
    Iobrd = imdilate(Iobr,se);
    Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
    Iobrcbr = imcomplement(Iobrcbr);
    
    % Compute regions of local minima usually where the soma is located
%     fgm = imextendedmax(Iobrcbr, 1, 4);
    fgm = imregionalmax(fftI);
    se2 = strel(ones(nhood2,nhood2));
    fgm2 = imclose(fgm, se2);
    fgm3 = imerode(fgm2, se2);
    fgm4 = bwareaopen(fgm3, obj_remove_pixel);
    
    % Compute background lines that separate image features
%     bw = imbinarize(Iobrcbr);
    bw = imbinarize(Iobrcbr);
    D = bwdist(bw);
    DL = watershed(D);
    bgm = DL == 0;
 
    gmag2 = imimposemin(gmag, bgm | fgm4);
    L = watershed(gmag2);
    labels = imdilate(L==0,ones(3,3)) + 2*bgm + 3*fgm4;
    Lrgb = label2rgb(L,'jet','w','shuffle');

    img_props = zeros(size(mip), 'logical');
    img_skel = zeros(size(mip), 'logical');
    props = struct([]);
    minI = uint8(255) - mip;
    binWatershed = L > 0;
    binWatershed = imclearborder(binWatershed);
    L = L .* uint8(binWatershed);
    iterL = unique(L); iterL(1) = [];
    cnt = 0;
    for i = 1:length(iterL)
        mask = L == iterL(i);
        maskI = minI .* uint8(mask);
        % Binarize image using Otsu's thresholding
        bwI = imbinarize(maskI, 'global');
        % Morphologically close image
        bwI_close = imclose(bwI, strel('disk', 5));
        bwI_holes = imfill(bwI_close, 'holes');
        % Remove small features
        bwI_open = bwareaopen(bwI_holes, 1000);
        % Go through IF loop only when voxels are detected AND there's one
        % continuous structure
        if nnz(bwI_open) > 0 & max(bwlabel(bwI_open), [], 'all') == 1
            cnt = cnt + 1;
            img_props = img_props + bwI_open;
            % Skeletonize image
            skel = bwskel(bwI_open);
            img_skel = img_skel + skel;
            props(i).Index = cnt;
            props(i).BranchPoints = nnz(bwmorph(skel, 'branchpoints'));
            props(i).EndPoints = nnz(bwmorph(skel, 'endpoints'));
            skel_props_ = regionprops(logical(skel), 'BoundingBox');
            skel_props = regionprops(logical(skel - bwmorph(skel, 'branchpoints')), 'BoundingBox', 'Area', 'Perimeter', 'Image');
            props(i).Branches = size(skel_props, 1);
            props(i).MaxBranchLength = max(extractfield(skel_props, 'Perimeter')/2);
            props(i).AverageBranchLength = mean(extractfield(skel_props, 'Perimeter')/2);
            props(i).TotalBranchLength = sum(extractfield(skel_props, 'Perimeter')/2);
            maskI_props = regionprops(bwI_open, 'BoundingBox', 'Area', 'Circularity', 'ConvexArea', 'ConvexImage', 'Image', 'MajorAxisLength', 'MinorAxisLength');
            x = fix(maskI_props.BoundingBox(1));
            y = fix(maskI_props.BoundingBox(2));
            bx = fix(maskI_props.BoundingBox(3));
            by = fix(maskI_props.BoundingBox(4));
            soma = imcrop(bwI_open, [x, y, bx, by]);
            props(i).BoundingBox = maskI_props(1).BoundingBox;
            props(i).Area = maskI_props(1).Area;
            props(i).FractalDimension = hausDim(soma);
            props(i).Circularity = maskI_props(1).Circularity;
            props(i).SpanRatio = maskI_props(1).MajorAxisLength /  maskI_props(1).MinorAxisLength;
            props(i).Density = nnz(bwperim(bwI_open)) / maskI_props(1).ConvexArea;
        else
            continue;
        end
    end
    img_props = logical(img_props);
    img_skel = logical(img_skel);
    % Compress props by removing empty fields
    props = props(all(~cellfun(@isempty,struct2cell(props))));
    % Convert to microns
    for i = 1:length(props)
        props(i).MaxBranchLength = props(i).MaxBranchLength * pixel_micron;
        props(i).AverageBranchLength = props(i).AverageBranchLength * pixel_micron;
        props(i).TotalBranchLength = props(i).TotalBranchLength  * pixel_micron;
        props(i).Area = props(i).Area * pixel_micron_area;
    end
        
    % Filter results based on thresholds
    rmIdx = extractfield(props, 'EndPoints') == 2 & extractfield(props, 'MaxBranchLength') < min_branch;
    props = props(~rmIdx);
    totalEndPoints = sum(extractfield(props, 'EndPoints'));
    rmIdx = extractfield(props, 'TotalBranchLength') < min_branch;
    props = props(~rmIdx);
    totalBranchLength = sum(extractfield(props, 'TotalBranchLength'));
    
    % Reindex
    for i = 1:length(props)
        props(i).Index = i;
    end

    fout = fullfile(outpath, strcat('Labels', '_', fname, '.png'));
    fig2 = figure('visible', 'off', 'WindowState','maximized');
    imshow(labeloverlay(labeloverlay(img, img_props, 'Colormap', [0 1 0], 'Transparency', 0.85), img_skel, 'Colormap', [1 0 0], 'Transparency', 0.10));
    hold on;
    for i = 1:size(props, 2)
        rectangle('Position', props(i).BoundingBox, 'EdgeColor', 'r',...
            'LineWidth', 0.5, 'LineStyle', ':');
        text(props(i).BoundingBox(1)+5, props(i).BoundingBox(2)+10,...
            num2str(props(i).Index), 'FontSize', 12);
    end
    hold off;
    title(strcat(fname, ': Detected Somas and Labels'), 'Interpreter', 'none');
    saveas(fig2, fout, 'png');
    
    % Write properties
    props = rmfield(props, 'BoundingBox');
    fname_table = fullfile(outpath, strcat('Stats_', fname, '.csv'));
    writetable(struct2table(props), fname_table);
    
%     % TEST
%     minI = uint8(255) - mdI;
%     for i = 1:max(L);
%         idx = L == i;
%         A = minI .* uint8(idx); figure; imshow(A);
%         
%         
%         
%     
%     
%     % Binarize image using Otsu's thresholding
%     bwI = ~imbinarize(mdI, 'global');
%     % Remove objects at border
%     bwI_bod = imclearborder(bwI);
%     % Median filter image again to remove specks
%     bwI_med = medfilt2(bwI_bod, [3,3]);
%     % Morphologically close image
%     bwI_close = imclose(bwI_med, strel('disk', 5));
%     % Remove small objects
%     bwI_open = bwareaopen(bwI_close, 1000);
%     % Skeletonize image
%     skel = bwskel(bwI_open);
%     % Find objects in image
%     props = regionprops(skel, 'Centroid', 'Area', 'BoundingBox', 'Image');
%     % Populate `props` with skeleton quantification
%     for i = 1:size(props, 1)
%         props(i).Index = i;
%         skelI = props(i).Image;
%         props(i).BranchPoints = nnz(bwmorph(skelI, 'branchpoints'));
%         props(i).EndPoints = nnz(bwmorph(skelI, 'endpoints'));
%         skel_props = regionprops(logical(skelI - bwmorph(skelI, 'branchpoints')), 'Area', 'Image');
%         props(i).Branches = size(skel_props, 1);
%         props(i).MaxBranchPoint = max(extractfield(skel_props, 'Area'));
%         props(i).AverageBranchPoint = mean(extractfield(skel_props, 'Area'));
%         props(i).TotalBranchPoint = sum(extractfield(skel_props, 'Area'));
%         % Save each soma
%         x = fix(props(i).BoundingBox(1)- 25);
%         y = fix(props(i).BoundingBox(2) - 25);
%         bx = fix(props(i).BoundingBox(3) + 50);
%         by = fix(props(i).BoundingBox(4) + 50);
%         props(i).RawImage = imcrop(mip, [x, y, bx, by]);
%         img_crop = imcrop(bwI_open, [x, y, bx, by]);
%         img_crop = imclearborder(img_crop);
%     %     img_crop = imfill(img_crop, 'holes'); % fill holes
%         img_crop = bwareafilt(img_crop, 1); % keep the biggest object in bounding box
%     %     img_crop = bwperim(img_crop);
%         props(i).PerimeterImage = img_crop;
%         img_props = regionprops(img_crop, 'Area', 'Circularity', 'ConvexArea', 'ConvexImage', 'Image', 'MajorAxisLength', 'MinorAxisLength');
%         props(i).FractalDimension = hausDim(img_crop);
%         props(i).Circularity = img_props(1).Circularity;
%         props(i).SpanRatio = img_props(1).MajorAxisLength /  img_props(1).MinorAxisLength;
%         props(i).Density = nnz(bwperim(img_crop)) / img_props(1).ConvexArea;
%         fname_crop = fullfile(outpath, strcat('S', num2str(i), '_', fname, '.tif'));
%         imwrite(img_crop, fname_crop);
%     end
%     % Filter results based on thresholds
%     rmIdx = extractfield(props, 'EndPoints') == 2 & extractfield(props, 'MaxBranchPoint') < 50;
%     props = props(~rmIdx);
%     totalEndPoints = sum(extractfield(props, 'EndPoints'));
%     rmIdx = extractfield(props, 'TotalBranchPoint') < 50;
%     props = props(~rmIdx);
%     totalBranchPoint = sum(extractfield(props, 'TotalBranchPoint'));
%     
%     props = rmfield(props, {'ConvexImage', 'Image', 'PerimeterImage', 'RawImage'});
%     % Write properties
%     fname_table = fullfile(outpath, strcat('Stats_', fname, '.csv'));
%     table_header = {'Index', 'Centroid', 'BoundingBox', 'Area', 'BranchPoints', 'EndPoints', 'Branches', 'MaxBranchPoint', 'AverageBranchPoint', 'TotalBranchPoint', 'FractalDimension', 'Circularity', 'SpanRatio', 'Density'};
%     props = orderfields(props, table_header);
%     writetable(struct2table(props), fname_table);
%  
%     % Make nice figure for QC
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
%     
%     fout = fullfile(outpath, strcat('Labels', '_', fname, '.png'));
%     fig2 = figure('visible', 'off', 'WindowState','maximized');
%     imshow(labeloverlay(labeloverlay(img, bwI_open, 'Colormap', [0 1 0], 'Transparency', 0.8), skel, 'Colormap', [1 0 0], 'Transparency', 0));
%     hold on;
%     for i = 1:size(props, 1)
%         rectangle('Position', props(i).BoundingBox, 'EdgeColor', 'r',...
%             'LineWidth', 0.5, 'LineStyle', ':');
%         text(props(i).BoundingBox(1)+5, props(i).BoundingBox(2)+10,...
%             num2str(props(i).Index), 'FontSize', 12);
%     end
%     hold off;
%     title(strcat(fname, ': Detected Somas and Labels'), 'Interpreter', 'none');
%     saveas(fig2, fout, 'png');
% end

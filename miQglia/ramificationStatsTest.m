function props = ramificationStatsTest(impath, outpath)
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
    pixel_micron = 0.25;   % Number of microns per pixel in XY. I measured 40 pixels for 10 um (microns)
    tophat_filter = 10;    % Size of top-hat filter. Ideally the spatial size of dendrites (in microns)
    obj_remove = 50;       % Remove features same size or under (in microns)
    min_branch = 5;        % Minimum length of total branch distance (microns)
    min_area = 100;        % Minimum area of microglia to keep in final stats (microns)
    Nr = 0;                % Number of microglia to select randomly in image. Set to zero disable random selection.
    Rseed = 1;             % Random number generation seed. A positive integer generates the same seed each
                           % run (to test multiple runs). Set to 0 to truly randomize blindly.

    %% Internal Pixel <--> Pixel Conversion (DO NOT ALTER!)
    obj_remove_pixel = round(obj_remove / pixel_micron, 0);
    pixel_micron_diag = sqrt(2* pixel_micron^2); % Diagonal distance of voxel
    pixel_micron_area =  pixel_micron^2; % Area of one pixel
    nhood1 = 10; % Neighborhood size of disk-shaped structuring element for initial erosion
    nhood2 = 3;  % Neighborhood size of disk-shaped structuring element closing
    nhood_th = round(tophat_filter / pixel_micron);
    S = rng(Rseed, 'twister');  % Seed random number generator
    
    %% Main
    [fpath,fname,fext] = fileparts(impath);
    img = imread(impath);
    % Maximum intensity projection to convert to grayscale
    mip = rgbmip(img);
    mip_ = uint8(max(mip(:))) - mip;     % inverted image;
    % Denoise image with 2D Weiner filter
    dnI = wiener2(mip_, [5,5]);
    % Remove redundant background
    thI = imtophat(dnI, strel('disk', nhood_th));
    % FFT bandpass filter to remove noise
    fftI = fft_bandpass(thI, 50, 100, 1);       % GOOD FOR SOMA DETECTION
    % Sharpen image
    shI = imsharpen(fftI, 'Radius', 3, 'Amount', 0.6);
    
    % Remove extreme pixel values and normalize from 0 to 1
    bot = double(prctile(shI(:),1));
    top = double(prctile(shI(:),99));
    I = (double(shI)-bot) / (top - bot);
    I(I>1) = 1; % do the removing
    I(I<0) = 0; % do the removing
    
    %% Create cellular mask
    mask = I > .10;
    mask = bwareaopen(mask, 500);
    
    %% Find Seeds
    I_smooth = imgaussfilt(I, 2);
%     seeds = imregionalmax(I_smooth);
%     seeds = imextendedmax(I_smooth, double(prctile(I_smooth(:),85)), 4);
    seeds = I_smooth >= 0.75;
    % Remove small objects
    seeds = bwareaopen(seeds, 200);
    % ADD METHOD TO REMOVE NON-ELLIPTICAL OBJECTS
    % Connect seeds close to each other
    seeds = bwmorph(seeds, 'bridge');
    seeds = bwmorph(imfill(seeds, 'holes'), 'shrink', Inf);
    seeds(mask==0)=0; % remove seeds outside of our cyto mask
    [X, Y] = find(seeds);
    seedsInd = sub2ind(size(mask), X, Y);
    
    %% Watershed
    I_smooth = imgaussfilt(I, 0.1); % don't smooth too much for watersheding
    I_min = imimposemin(max(I_smooth(:)) - I_smooth, seeds); % set locations of seeds to be -Inf (cause matlab watershed)
    L = watershed(I_min);
    L(mask==0)=0; % remove areas that aren't in our cellular mask
    % Remove segmentations that don't overlap with seeds
    L_mask = unique(L(seeds));
    L_labels = unique(L(:));
    for i = length(L_labels)
        if ismember(L_labels(i), L_mask)
            L(L==L_labels(i)) = L_labels(i);
        else
            L(L==L_labels(i)) = 0;
        end
    end
    binWatershed = L > 0;
    binWatershed = imclearborder(binWatershed);
%     binWatershed = imerode(binWatershed, strel('disk', 3));
    binWatershed = bwmorph(binWatershed, 'shrink', 2);
    L = L .* uint8(binWatershed);
    
%     figure;       % for debugging
%     imshow(labeloverlay(labeloverlay(img, imdilate(bwperim(L), ones(2,2)), 'Colormap', [0 0 1], 'Transparency', 0.5), mask, 'Colormap', [0 1 0], 'Transparency', 1));
%     hold on;
%     plot(Y,X,'x','markersize',10,'markeredgecolor','g','markeredgecolor','r');
%     plot(Y,X,'o','markersize',10,'markeredgecolor','g','markeredgecolor','r')
%     hold off;
    
    % Restore image to unint8
    I = uint8(I*255);
    
    img_props = zeros(size(mip), 'logical');
    img_skel = zeros(size(mip), 'logical');
    props = struct([]);
    iterL = unique(L); iterL(1) = [];
    cnt = 0;
    for i = 1:length(iterL)
        cell = L == iterL(i);
        bwI_open = bwareafilt(cell, 1); % Keep the largest object
        tmp = regionprops(bwI_open, 'PixelIdxList');
        % Go through IF loop only when voxels are detected AND there's one
        % continuous structure
        if nnz(bwI_open) > 0  & any(ismember(seedsInd, tmp.PixelIdxList)) %max(bwlabel(bwI_open), [], 'all') == 1
            cnt = cnt + 1;
            img_props = img_props + bwI_open;
            % Skeletonize image
            skel = bwskel(bwI_open);
            img_skel = img_skel + skel;
            props(i).Index = cnt;
            props(i).File = fname;
            props(i).BranchPoints = nnz(bwmorph(skel, 'branchpoints'));
            props(i).EndPoints = nnz(bwmorph(skel, 'endpoints'));
%             skel_props_ = regionprops(logical(skel), 'BoundingBox');
            skel_props = regionprops(logical(skel - bwmorph(skel, 'branchpoints')), 'BoundingBox', 'Area', 'Perimeter', 'Image');
            props(i).Branches = size(skel_props, 1);
            props(i).MaxBranchLength = max(extractfield(skel_props, 'Perimeter')/2);
            props(i).AverageBranchLength = mean(extractfield(skel_props, 'Perimeter')/2);
            props(i).TotalBranchLength = sum(extractfield(skel_props, 'Perimeter')/2);
            maskI_props = regionprops(bwI_open, 'BoundingBox','PixelIdxList', 'Area', 'Circularity', 'ConvexArea', 'ConvexImage', 'Image', 'MajorAxisLength', 'MinorAxisLength');
            x = fix(maskI_props.BoundingBox(1));
            y = fix(maskI_props.BoundingBox(2));
            bx = fix(maskI_props.BoundingBox(3));
            by = fix(maskI_props.BoundingBox(4));
            soma = imcrop(bwI_open, [x, y, bx, by]);
            props(i).PixelIdxList = maskI_props.PixelIdxList;
            props(i).BoundingBox = maskI_props(1).BoundingBox;
            props(i).Area = maskI_props(1).Area;
            props(i).FractalDimension = hausDim(soma);
            props(i).Circularity = maskI_props(1).Circularity;
            props(i).SpanRatio =  maskI_props(1).MinorAxisLength / maskI_props(1).MajorAxisLength;
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
    rmIdx = extractfield(props, 'Area') <= min_area;
    props = props(~rmIdx);
    
    % Randomly select microglia (if specified)
    if Nr > 0
        if Nr > length(props)
            Nr = length(props);
        end
        rng(S);
        Ridx = randperm(length(props), Nr);
        kpIdx = zeros(length(props), 1, 'logical');
        for k = 1:length(Ridx)
            kpIdx(Ridx(k)) = true;
        end
        props = props(kpIdx);
    end
       
    % Remove seeds not present in image properties (props)
    seedsRm = logical(ismember(seedsInd, vertcat(props.PixelIdxList)));
    X = X(seedsRm);
    Y = Y(seedsRm);
    % Remove connected component image (img_props) objects not present in
    % image properties (props)
    img_props = bwselect(img_props, Y, X);
    % Remove skeleton (img_skel) objects not present in image properties
    % (props)
    [rrm, crm] = ind2sub(size(mask), find(img_props));
    img_skel = bwselect(img_skel, crm, rrm);
    
    % Reindex
    for i = 1:length(props)
        props(i).Index = i;
    end
    
    fout = fullfile(outpath, strcat('Labels', '_', fname, '.png'));
    fig2 = figure('visible', 'off', 'WindowState','maximized');
    imshow(labeloverlay(labeloverlay(img, imdilate(bwperim(img_props), ones(2,2)), 'Colormap', [0 0 1], 'Transparency', 0.5), img_skel, 'Colormap', [1 0 0], 'Transparency', 0.5));
%     imshow(labeloverlay(labeloverlay(img, img_props, 'Colormap', [0 1 0], 'Transparency', 0.85), img_skel, 'Colormap', [1 0 0], 'Transparency', 0.10));
    hold on;
    plot(Y,X,'.','markersize',30,'markeredgecolor','r');
%     plot(Y,X,'o','markersize',10,'markeredgecolor','g','markeredgecolor','r')
    for i = 1:length(props)
        text(Y(i)+5, X(i)-5,...
            num2str(props(i).Index), 'FontSize', 20, 'Color', '#00FF00');
    end
    hold off;
    title(strcat(fname, ': Detected Somas and Labels'), 'Interpreter', 'none');
    saveas(fig2, fout, 'png');
    
    % Write properties
    props = rmfield(props, {'BoundingBox', 'PixelIdxList'});
    fname_table = fullfile(outpath, strcat('Stats_', fname, '.csv'));
    writetable(struct2table(props), fname_table);
    
function selectedVoxel = interactiveVoxelEditor(voxelData, backgroundImg)
    % voxelData: 3D binary volume
    % backgroundImg: 2D grayscale image (same XY size as voxelData)

    % Initialize figure
    hFig = figure('Name', 'Interactive Voxel Editor', 'NumberTitle', 'off',...
                  'KeyPressFcn', @keyPress);

    % Calculate initial labeling and projection
    labeledVolume = bwlabeln(voxelData);
    projection = squeeze(max(labeledVolume, [], 3));
    overlayImg = im_merge(cat(3,mat2gray(backgroundImg)*5,mat2gray(projection)),[0 0.6 1;1 0 1]);

    % Initialize selectedVoxel and label tracking
    selectedVoxel = false(size(voxelData));
    selectedLabels = [];
    labelSelectMode = false;

    % Display initial projection
    hAx = axes('Parent', hFig);
    hImg = imshow(renderLabelOverlay(projection, selectedLabels, backgroundImg, false), [], 'Parent', hAx);
    title(hAx, 'Max Projection');

    % Set mouse click callback
    set(hImg, 'ButtonDownFcn', @clickToToggleLabel);

    % UI buttons
    uicontrol('Style', 'pushbutton', 'String', 'Add ROI',...
        'Units', 'normalized', 'Position', [0.05, 0.02, 0.2, 0.05],...
        'Callback', @(~,~)updateROI(true));

    uicontrol('Style', 'pushbutton', 'String', 'Remove ROI',...
        'Units', 'normalized', 'Position', [0.28, 0.02, 0.2, 0.05],...
        'Callback', @(~,~)updateROI(false));

    uicontrol('Style', 'togglebutton', 'String', 'Select Labels',...
        'Units', 'normalized', 'Position', [0.51, 0.02, 0.2, 0.05],...
        'Callback', @(btn,~)toggleLabelSelect(btn));

    % Wait for figure to close before returning
    uiwait(hFig);

    function updateROI(isAdd)
        [~, ROImask] = get_ROI(overlayImg);
        if isempty(ROImask)
            return;
        end

        xyMask = max(ROImask, [], 3);
        mask3D = repmat(xyMask, [1, 1, size(voxelData, 3)]);

        if isAdd
            voxelData(mask3D > 0) = 1;
        else
            voxelData(mask3D > 0) = 0;
        end

        labeledVolume = bwlabeln(voxelData);
        projection = squeeze(max(labeledVolume, [], 3));
        hImg.CData = renderLabelOverlay(projection, selectedLabels, backgroundImg, true); % show boundary
    end

    function toggleLabelSelect(button)
        labelSelectMode = get(button, 'Value');
        % Don't change current image rendering here
        hImg.CData = renderLabelOverlay(projection, selectedLabels, backgroundImg, true);
    end

    function clickToToggleLabel(~, event)
        if ~labelSelectMode
            return;
        end
        clickPos = round(event.IntersectionPoint(1:2));
        x = min(max(clickPos(1), 1), size(projection, 2));
        y = min(max(clickPos(2), 1), size(projection, 1));
        label = projection(y, x);
        if label == 0
            return;
        end
        if ismember(label, selectedLabels)
            selectedLabels(selectedLabels == label) = [];
        else
            selectedLabels(end+1) = label;
        end
        selectedVoxel = ismember(labeledVolume, selectedLabels);
        hImg.CData = renderLabelOverlay(projection, selectedLabels, backgroundImg, false); % no boundary
    end

    function rgbImg = renderLabelOverlay(labelProj, selected, bgImg, showBoundary)
    bgNorm = uint8(mat2gray(bgImg) * 255);  % Normalize grayscale to 0-255
    rgbImg = repmat(bgNorm, 1, 1, 3);       % RGB version of background

    if showBoundary
        boundary = bwperim(labelProj > 0);  % Draw boundaries of any non-zero labels
        rgbImg(:,:,1) = rgbImg(:,:,1);      % Red stays unchanged
        rgbImg(:,:,2) = max(rgbImg(:,:,2), uint8(255 * boundary));  % Green lines
        rgbImg(:,:,3) = rgbImg(:,:,3);      % Blue unchanged
    else
        mask = ismember(labelProj, selected);
        rgbImg(:,:,1) = rgbImg(:,:,1) + uint8(127 * mask);
        rgbImg(:,:,2) = rgbImg(:,:,2) .* uint8(~mask);
        rgbImg(:,:,3) = rgbImg(:,:,3) .* uint8(~mask);
    end
end

    function keyPress(~, event)
        if strcmp(event.Key, 'q')
            uiresume(hFig);
            close(hFig);
        end
    end
end

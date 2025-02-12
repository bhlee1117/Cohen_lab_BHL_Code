function Result = removeROI(Result, roiIndex)
    % Function to remove a specific ROI from the Result structure
    % Inputs:
    %   Result - Original structure containing 28 ROIs
    %   roiIndex - Index of the ROI to remove (1 to 28)
    % Outputs:
    %   Result - Updated structure with the specified ROI removed

    % Validate the ROI index
    if roiIndex < 1 || roiIndex > numel(Result.ROIpoly)
        error('Invalid ROI index. Must be between 1 and %d.', numel(Result.ROIpoly));
    end

    % Remove the ROI from each field
    Result.ROIpoly(roiIndex) = [];                      % Remove from ROIpoly
    Result.ftprnt(:, :, roiIndex) = [];                 % Remove from ftprnt
    Result.traces(roiIndex, :) = [];                    % Remove from traces
    remove_dist=Result.dist_order(roiIndex);
    Result.dist_order(roiIndex) = [];                   % Remove from dist_order
    Result.dist_order(Result.dist_order>remove_dist)=Result.dist_order(Result.dist_order>remove_dist)-1;
    Result.traces_bvMask(roiIndex, :) = [];             % Remove from traces_bvMask
    Result.normTraces(roiIndex, :) = [];                % Remove from normTraces
    Result.spike(roiIndex, :) = [];                     % Remove from spike
    Result.interDendDist(roiIndex) = [];                % Remove from interDendDist
    if isfield(Result,'F0_PCA')
    Result.F0_PCA(roiIndex) = [];                       % Remove from F0_PCA
    end

    % Display a message
    fprintf('ROI %d has been successfully removed.\n', roiIndex);
end
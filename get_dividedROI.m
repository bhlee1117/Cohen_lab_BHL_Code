function [roiOut, mapping] = get_dividedROI(roiIn, maxSize)
% Splits large ROI masks in roiIn into smaller masks, each ≤ maxSize pixels.
% Preferentially splits along the longer dimension of the ROI.
% 
% INPUTS:
%   roiIn   - X x Y x N matrix. Each slice is an ROI with non-zero values.
%   maxSize - Maximum allowed number of pixels for a single ROI.
%
% OUTPUTS:
%   roiOut  - X x Y x M matrix. Each slice is a split or original ROI with preserved values.
%   mapping - M x 1 vector. mapping(i) tells which original ROI (1 to N) this slice came from.

    [X, Y, N] = size(roiIn);
    roiList = {};
    roiMap = [];

    for k = 1:N
        roiSlice = roiIn(:, :, k);
        bw = roiSlice > 0;
        cc = bwconncomp(bw);  % Identify connected components

        for j = 1:cc.NumObjects
            pixIdx = cc.PixelIdxList{j};
            nPix = numel(pixIdx);

            if nPix <= maxSize
                temp = zeros(X, Y);
                temp(pixIdx) = roiSlice(pixIdx);
                roiList{end+1} = temp;
                roiMap(end+1, 1) = k;
            else
                % Split large component into smaller pieces
                [rows, cols] = ind2sub([X, Y], pixIdx);
                numChunks = ceil(nPix / maxSize);

                % Determine longer direction
                rowRange = range(rows);
                colRange = range(cols);
                if rowRange >= colRange
                    featureVec = rows;
                else
                    featureVec = cols;
                end

                % Sort and split based on longer axis
                [~, sortIdx] = sort(featureVec);
                rows = rows(sortIdx);
                cols = cols(sortIdx);
                pixIdx = pixIdx(sortIdx);

                % Equal-size chunks based on sorted order
                edges = round(linspace(1, nPix+1, numChunks+1));
                for m = 1:numChunks
                    idxRange = edges(m):(edges(m+1)-1);
                    selectedIdx = pixIdx(idxRange);
                    temp = zeros(X, Y);
                    temp(selectedIdx) = roiSlice(selectedIdx);
                    roiList{end+1} = temp;
                    roiMap(end+1, 1) = k;
                end
            end
        end
    end

    % Convert list to 3D matrix
    roiOut = cat(3, roiList{:});
    mapping = roiMap;
end

function [kymo kymoROI]=polyLineKymo_lineInput(mov,dL,dP,lineCoord)
% POLYLINEKYMO_LINEINPUT  Generate kymograph ROIs along a polyline path.
%
%   [kymo, kymoROI] = polyLineKymo_lineInput(mov, dL, dP, lineCoord)
%   samples rectangular regions of interest (ROIs) along a user-defined
%   polyline and applies them to an image stack or movie to create a
%   kymograph.
%
%   INPUTS:
%       mov        - 3D array (y × x × frames) representing the image stack.
%       dL         - Step size along the polyline in pixels (spacing of ROIs).
%       dP         - Perpendicular width of each rectangular ROI in pixels.
%       lineCoord  - N×2 matrix of (x,y) coordinates defining the polyline.
%
%   OUTPUTS:
%       kymo       - Kymograph extracted by sampling all ROIs across frames.
%       kymoROI    - Cell array of polygon coordinates (Xs,Ys) for each ROI.
%
%   The function:
%       • Divides each polyline segment into steps spaced by dL.
%       • Builds rectangular ROIs of width dP perpendicular to the line.
%       • Clamps ROI boundaries to the movie dimensions to avoid out-of-bounds.
%       • Calls APPLY_CLICKY to extract intensity profiles and generate
%         the kymograph.
%
%   Example:
%       [kymo, roi] = polyLineKymo_lineInput(movieStack, 2, 5, linePts);
%       imagesc(kymo); colormap(gray);
    X = lineCoord(:,1);
    Y = lineCoord(:,2);
    dX = diff(X);
    dY = diff(Y);
    L  = hypot(dX, dY);          % Segment lengths
    nSeg = numel(L);

    % Unit vectors along and perpendicular to segments
    dX0 = dX ./ L;
    dY0 = dY ./ L;
    dXp = -dY0;
    dYp =  dX0;

    % Determine image bounds only once
    [ySize, xSize, ~] = size(mov);

    % Estimate total number of ROIs for preallocation
    totalSteps = sum(floor(L./dL)+1);
    roi  = cell(totalSteps,1);
    trail = zeros(totalSteps,2);

    c = 1;
    for k = 1:nSeg
        % Positions along this segment
        lSteps = 0:dL:L(k);
        nSteps = numel(lSteps);

        % Precompute cos/sin-like offsets for efficiency
        for j = 1:nSteps
            Xc = X(k) + lSteps(j)*dX0(k);
            Yc = Y(k) + lSteps(j)*dY0(k);

            % Save trail
            trail(c,:) = [Xc Yc];

            % Compute rectangle corners
            dL2 = dL/2; dP2 = dP/2;
            XoffL = dX0(k)*dL2;  YoffL = dY0(k)*dL2;
            XoffP = dXp(k)*dP2;  YoffP = dYp(k)*dP2;

            Xs = [Xc - XoffL - XoffP;
                  Xc - XoffL + XoffP;
                  Xc + XoffL + XoffP;
                  Xc + XoffL - XoffP;
                  Xc - XoffL - XoffP];

            Ys = [Yc - YoffL - YoffP;
                  Yc - YoffL + YoffP;
                  Yc + YoffL + YoffP;
                  Yc + YoffL - YoffP;
                  Yc - YoffL - YoffP];

            % Clamp to image bounds
            Xs = min(max(Xs,1), xSize);
            Ys = min(max(Ys,1), ySize);

            roi{c} = [Xs, Ys];
            c = c + 1;
        end
    end

    % Compute kymograph
    kymo = apply_clicky(roi', mov,'no');  % 'no' flag optional
    kymoROI = roi';
end

function fibStruct = findFib(refImg, rawImg, fibW, minL, maxL, dW, varargin);

% fibStruct = findFib(refImg, rawImg, fibW, minL, maxL, dW);
% find all the fibers in an image.
% INPUTS
% refImg: reference monochrome image used for finding the fibers.  If
% fibers have multiple colors, construct this image by combining suitably
% scaled color channels.  
% rawImg: raw data image from which fiber profiles are to be extracted.
% Should have same xSize and ySize as refImg, but can have up to three
% color planes.
% fibW: estimate of fiber width, in pixels
% minL: minimum fiber length, in pixels
% maxL: maximum fiber length, in pixels
% dW: 2*dW+1 is the width of the returned fiber images.
% varargin: if included, should be an option (OPT) structure with 2 fields:
% (1) OPT.DEBUG_PLOTS (true = show plots/images [slower]; false = don't show); (2)
% OPT.USE_IMPROFILE (true = use improfile [slower]; false = direct call to
% interp2). If not included, defaults are OPT.DEBUG_PLOTS = false and
% OPT.USE_IMPROFILE = false
%
% OUTPUTS
% fibStruct: struct array with one element per fiber, containing subfields
% X: 2-element array with x-coords of ends
% Y: 2-element array with y-coords of ends
% L: length of fiber
% Img: 2*dW + 1 by L by 3 image of the fiber
% edges: 2-element array on whether the ends of the fiber are within the
% field of view (0 = in-field, 1 = out-of-field)
%
% Adam Cohen 24 Feb. 2022
% Modified Eric M. Moult 07 March 2022

if isempty(varargin)
    OPT.DEBUG_PLOTS = false;
    OPT.USE_IMPROFILE = false;
else
    TMP_OPT = varargin{1};
    try 
        OPT.DEBUG_PLOTS = TMP_OPT.DEBUG_PLOTS;
    catch
        OPT.DEBUG_PLOTS = false;
    end

    try 
        OPT.USE_IMPROFILE = TMP_OPT.USE_IMPROFILE;
    catch
        OPT.USE_IMPROFILE = false;
    end
end

theta = [0:179];
[ySize, xSize, nChan] = size(rawImg);

if nChan == 3;
    colorImg = zeros(ySize, xSize, 3);
    for j = 1:3;
        colorImg(:,:,j) = mat2gray(rawImg(:,:,j));
    end;
else
    colorImg = mat2gray(rawImg);
end;

if OPT.DEBUG_PLOTS
    figure(100)
    subplot(1,2,1);
    imshow2(colorImg);
end

thresh = graythresh(refImg);  % Automatically set a grayscale threshold.
mask = refImg > thresh;

if OPT.DEBUG_PLOTS
    subplot(1,2,2);
    imshow2(mask)
end

stats = regionprops(mask, 'Area', 'BoundingBox');
a = [stats(:).Area];
stats(a < 80) = [];  % Get rid of regions with < 80 pixels
nROI = length(stats);  

fibC = 0;  % fiber count.
if OPT.DEBUG_PLOTS
    figure(101); clf
    imshow2(colorImg); hold all
end
'Finding the fibers'
% tic
for j = 1:nROI;  % Go through each ROI.
    coords = round(stats(j).BoundingBox);  % Set the dimensions of the box around each fiber
    xMin = max(coords(1)-dW,1); xMax = min(coords(1) + coords(3) + dW, xSize);  % Make the box a little bigger
    yMin = max(coords(2)-dW,1); yMax = min(coords(2) + coords(4) + dW, ySize);

    subImg = refImg(yMin:yMax, xMin:xMax);
    [ySizeSub, xSizeSub] = size(subImg);
    L2 = (ySizeSub^2 + xSizeSub^2)^.5;  % diagonal length of image
    yCent = ySizeSub/2;  % center of the sub-image
    xCent = xSizeSub/2;
    [radonImg, xp] = radon(subImg, theta);

    if OPT.DEBUG_PLOTS
        figure(102);
        subplot(1,3,1); hold off% Show the radon image
        imagesc(theta, xp, radonImg); colormap('hot'); hold all                    
        xlabel('Angle (deg)')
        ylabel('Position (pix)')
        pbaspect([1 1 1])
        title('Radon transform')

        subplot(1,3,2); hold off% Show the calculated line
        imagesc(subImg); hold all
        axis('off')
        daspect([1 1 1])
    end

    % Find all the maxima in the Radon image which could be fibers
    filt = fspecial('gaussian', [8*fibW, 1], fibW);  % filter on fiber width
    filt = filt - mean(filt);
    radonImg2 = imfilter(radonImg, filt, 'replicate');
    radonMask = radonImg2 > (mean(radonImg2(:)) + 4*std(radonImg2(:)));  % adjust this threshold
    
    if OPT.DEBUG_PLOTS
        figure(104); 
        subplot(1,2,1); imshow2(radonImg2, []);
        subplot(1,2,2); imshow2(radonMask, []);
    end
    b = regionprops(radonMask, radonImg2, 'MaxIntensity', 'Area');
    b([b(:).Area] < 20) = []; % adjust this minimum area threshold

    nReg = length(b);  % number of possible fibers within the ROI
    for k = 1:nReg;
        [xM, thetaM] = find(radonImg2 == b(k).MaxIntensity);  % location of the maximum
        if ~((thetaM == 1) | (thetaM == 180));  % don't allow fibers at the boundaries
            rM = xp(xM);  % radial distance of closest approach of the line to the center of the image
            qM = theta(thetaM);  % angle of the line
            
            if OPT.DEBUG_PLOTS
                figure(102); subplot(1,3,1);
                plot(qM, rM, 'cx');
            end

            % Extract an image of the fiber.
            %prof = zeros(2*dW+1,maxL);  % fiber profile image
            % Assume the fiber has maximum length
            x02 = xCent + rM*cosd(qM) + xMin;  % coordinates of line on whole refImg
            y02 = yCent - rM*sind(qM) + yMin;
            X2 = [x02 - maxL*cosd(qM + 90)/2, x02 + maxL*cosd(qM + 90)/2];
            Y2 = [y02 + maxL*sind(qM + 90)/2, y02 - maxL*sind(qM + 90)/2];
            % Get the maximal length fiber profile
            prof = getOneFib(X2, Y2, dW, refImg, OPT);
            
            % Find the ends of the fiber
            fibLine = prof(dW+1,:) - (prof(dW+1-2*fibW,:) + prof(dW+1+2*fibW,:))/2;  % background-subtracted line-profile
            fibLine2 = fibLine;
            fibLine2(isnan(fibLine)) = 0;  % clean up NaNs
            fibLine2 = smooth(fibLine2, 11); % smooth out small fluctuations
            maxI = nanmean(fibLine((maxL/2 - 9):(maxL/2+11)));  % average intensity near the middle of the line
            rightEdge = find(fibLine2(maxL/2+1:end) < 0.1*maxI, 1)+maxL/2;
            rightEdge = min(rightEdge+dW, length(fibLine));  % make the edges of the profile a little beyond the end of the fiber
            leftEdge = maxL/2+1 - find(fibLine2(maxL/2:-1:1) < 0.1*maxI, 1);
            leftEdge = max(leftEdge-dW, 1);
            if (rightEdge-leftEdge > minL) & (rightEdge-leftEdge < maxL)
                fibC = fibC + 1;  % fiber index count
                fibStruct(fibC).theta = qM;
                allQ(fibC) = qM;
                
                % plot the fiber on the small ROI
                x0 = xCent + rM*cosd(qM);  % coordinates of closest approach to center
                y0 = yCent - rM*sind(qM);

                X = [x0 - L2*cosd(qM + 90), x0 + L2*cosd(qM + 90)];
                Y = [y0 + L2*sind(qM + 90), y0 - L2*sind(qM + 90)];
                
                if OPT.DEBUG_PLOTS
                    figure(102); subplot(1,3,2);
                    plot(X, Y, 'g-');
                end
                
                % plot the fiber on the whole FOV
                if OPT.DEBUG_PLOTS
                    figure(101); hold all  % plot the location of the fiber on the main image
                end
                
                x02 = xCent + rM*cosd(qM) + xMin;  % coordinates of closest approach to center
                y02 = yCent - rM*sind(qM) + yMin;
                X2 = [x02 - (maxL/2-leftEdge)*cosd(qM + 90), x02 + (rightEdge-maxL/2)*cosd(qM + 90)];
                Y2 = [y02 + (maxL/2-leftEdge)*sind(qM + 90), y02 - (rightEdge-maxL/2)*sind(qM + 90)];
                
                if OPT.DEBUG_PLOTS
                    plot(X2, Y2, 'g-', 'LineWidth', 0.5);
                end
                
                fibStruct(fibC).X = X2;
                fibStruct(fibC).Y = Y2;
                allX0(fibC,:) = X2;  % store the locations of the line ends
                allY0(fibC,:) = Y2;
                
                prof = getOneFib(X2, Y2, dW, colorImg, OPT);
                
                if isnan(fibLine(leftEdge))
                    leftOff(fibC) = 1;
                else
                    leftOff(fibC) = 0;
                end;
                if isnan(fibLine(rightEdge))
                    rightOff(fibC) = 1;
                else
                    rightOff(fibC) = 0;
                end;
%                 outImgs{fibC} = prof;
                fibStruct(fibC).img = prof;
                fibStruct(fibC).leftOff = leftOff(fibC);
                fibStruct(fibC).rightOff = rightOff(fibC);
                
                if OPT.DEBUG_PLOTS
                    figure(102); 
                    subplot(1,3,3); hold off
                    imshow2(mat2gray(prof));
                    title(['Fiber ' num2str(fibC) ' Left Off: ' num2str(leftOff(fibC)) ' Right Off: ' num2str(rightOff(fibC))])
                    drawnow
                end
            end;
        end;
    end;
end;
% toc
% remove duplicates
'Remove the duplicates'
% tic
distMat = zeros(fibC, fibC, 2);  % Distance of each end of fiber n to the line-segment m
if OPT.DEBUG_PLOTS
    figure(101); hold all
end
for m = 1:fibC;
    Lm(m) = (diff(allX0(m,:))^2 + diff(allY0(m,:))^2).^.5;  % length of line-segment m
    fibStruct(m).L = Lm;
    for n = 1:fibC;
        for p = 1:2;  % the two ends of fiber n
            vecm = [diff(allX0(m,:)), diff(allY0(m,:))];  % Vector of segment m
            vecp1 = [allX0(n,p) - allX0(m,1), allY0(n,p) - allY0(m,1)]; % point on n to first point on m
            vecp2 = [allX0(n,p) - allX0(m,2), allY0(n,p) - allY0(m,2)]; % point on n to second point on m
            dot1 = vecm*vecp1';  
            dot2 = -vecm*vecp2';
            if (dot1 < 0) | (dot2 < 0); % The dot products are negative if the closest approach of n to the m-line lies outside the edges of m.
                dp = min([vecp1*vecp1', vecp2*vecp2'])^.5;  % closest approach of point on n to either end of m
            else;  % orthogonal distance of n to the m-segment
                dp = abs((allX0(m,2) - allX0(m,1))*(allY0(m,1) - allY0(n,p)) - ...
                    (allX0(m,1) - allX0(n,p))*(allY0(m,2) - allY0(m,1)))/Lm(m);  % closest approach of point on n to m-segment
            end;
            distMat(m,n,p) = dp;
        end;
    end;
end;
goodFib = ones(fibC,1);  % 
for m = 1:fibC;
    for n = setdiff(1:fibC, m);
        if ((distMat(m,n,1) < 10) | (distMat(m,n,2) < 10))& (abs(sind(allQ(m) - allQ(n))) < .1);
            if OPT.DEBUG_PLOTS
                plot(allX0(m,:), allY0(m,:), 'r-');
                plot(allX0(n,:), allY0(n,:), 'c-');
                title([num2str(m) ' ' num2str(n)])
            end
            if Lm(m) >= Lm(n);
                goodFib(n) = 0;
            else
                goodFib(m) = 0;
            end;
            if OPT.DEBUG_PLOTS
                drawnow
            end
        end;
    end;
end;
goodFib = logical(goodFib);
fibStruct = fibStruct(goodFib);
allQ = allQ(goodFib);
allX0 = allX0(goodFib,:);
allY0 = allY0(goodFib,:);
fibC = sum(goodFib);

% toc

if OPT.DEBUG_PLOTS
    figure(101); clf
    imshow2(colorImg, []); hold all
    for m = 1:fibC;
        plot(allX0(m,:), allY0(m,:), 'g-');
    end;
    hold off
end

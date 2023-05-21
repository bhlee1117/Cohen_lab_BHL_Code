dList = dir('int_*');
nDir = length(dList);
ySize = 2048; xSize = 2048;

theta = [0:179];
pixSize = 1;  % need to update with real scalebar


for d = 1:nDir;  % Go through all the directories
    cd(dList(d).name);
    fList = dir('*Unmixing.czi');
    nFiles = length(fList);
    for f = 1:nFiles;  % go through all the unmixed files in each directory
        fName = fList(f).name;
        dat = bfopen(fName);  % open each Unmixed image.
        [ySize, xSize] = size(dat{1}{1,1});  % image dimensions
        L = sqrt(xSize^2 + ySize^2);  % image diagonal length
        
        colorImg = zeros(ySize, xSize, 3);
        for j = 1:3;
            colorImg(:,:,j) = mat2gray(double(dat{1}{j,1}));  % make a 3-color image
        end;
        figure(f)
        imshow2(colorImg)
        
        refImg = max(colorImg,[], 3);  % define a one-color reference image for finding fibers.
        figure(f + nFiles);
        subplot(1,2,1);
        imshow2(refImg);
        thresh = graythresh(refImg);  % Automatically set a grayscale threshold.
        mask = refImg > thresh;
        subplot(1,2,2);
        imshow2(mask)
        stats = regionprops(mask, 'Area', 'BoundingBox')
        a = [stats(:).Area];
        stats(a < 10) = [];  % Get rid of regions with < 10 pixels
        nROI = length(stats);  
        for j = 1:nROI;  % Go through each ROI.
            coords = round(stats(j).BoundingBox);  % Set the dimensions of the box around each fiber
            xMin = max(coords(1)-10,1); xMax = min(coords(1) + coords(3) + 10, xSize);
            yMin = max(coords(2)-10,1); yMax = min(coords(2) + coords(4) + 10, ySize);
            
            subImg = refImg(yMin:yMax, xMin:xMax);
            [ySizeSub, xSizeSub] = size(subImg);
            yCent = ySizeSub/2;  % Find the center of the sub-movie
            xCent = xSizeSub/2;
            [radonImg, xp] = radon(subImg, theta);
            [maxI, idx] = max(radonImg, [], 'all', 'linear');  % Find the linear index of max value in Radon xform
%         allIntens(j,c) = maxI;
            [xM, thetaM] = ind2sub(size(radonImg), idx);  % convert to subscripts
            rM = xp(xM);  % closest approach of the line to the center of the image
            qM = theta(thetaM);  % angle of the line
            subplot(1,3,1); hold off % Show the maximum on the radon image
            imagesc(theta, xp, radonImg); colormap('hot'); hold all
            plot(qM, rM, 'cx'); hold off
            xlabel('Angle (deg)')
            ylabel('Position (pix)')
            pbaspect([1 1 1])
            title('Radon transform')
            subplot(1,3,2); hold off;  % Show the calculated line
            imagesc(subImg); hold all
            axis('off')
            daspect([1 1 1])
            x0 = xCent + rM*cosd(qM);  % coordinates of closest approach to center
            y0 = yCent - rM*sind(qM);
%             allX0(j,c) = x0;
%             allY0(j,c) = y0;
%             allQ(j,c) = qM + 90;
%             allRM(j,c) = rM;
    %     plot(xCent, yCent, 'co');
    %     plot([xCent, x0], [yCent, y0], 'b.-')
            X = [x0 - L*cosd(qM+90), x0 + L*cosd(qM+90)];
            Y = [y0 + L*sind(qM + 90), y0 - L*sind(qM + 90)];
            plot(X, Y, 'g-', 'LineWidth', 1);
            subplot(1,3,3)
            imagesc(colorImg); hold all
            axis('off')
            daspect([1 1 1])
            x02 = x0 + coords(1);
            y02 = y0 + coords(2);
            X2 = [x02 - L*cosd(qM+90), x02 + L*cosd(qM+90)];
            Y2 = [y02 + L*sind(qM + 90), y02 - L*sind(qM + 90)];
            plot(X2, Y2, 'g-', 'LineWidth', 1);

            % Extract an image of the fiber.
            subMov = colorImg(yMin:yMax,xMin:xMax,:);
            
    % geometry to rotate the image
    dX = diff(X);
    dY = diff(Y);
    theta = atan2d(dY, dX);
    img2 = imrotate(inImgs, theta);
    [ySize2, xSize2, ~] = size(img2);
    
    % geometry to rotate the line
    X0 = xSize/2; Y0 = ySize/2;
    xRel = X - X0;
    yRel = Y - Y0;
    rotMat = [cosd(-theta) -sind(-theta); sind(-theta) cosd(-theta)];
    coordsF = rotMat*([xRel'; yRel']);
%     XF = coordsF(1,:) + X0;
%     YF = coordsF(2,:) + Y0;
    XF = coordsF(1,:) + xSize2/2;
    YF = coordsF(2,:) + ySize2/2;
    XF = max(XF, 1); YF = max(YF, 1);
    XF = min(XF, xSize2); YF = min(YF, ySize2);

    yBot = max(YF(1)-dW, 1);
    yTop = min(YF(1)+dW, ySize2);
    outImgs{j} = img2(round(yBot:yTop),round(XF(1):XF(2)),:);
    for k = 1:3
        outImgs{j}(:,:,k) = mat2gray(outImgs{j}(:,:,k));
    end;
    subplot(nROI,1,j)
    imshow2(outImgs{j});
            
            
            
            
            

        % Extract the profile of the line, averaging over a stripe of width
        % 2*dW+1
            dW = 2; % number of pixels to average on either side of the line;
            tmp = [];
            for k = -dW:dW;
                x1 = xCent + (rM+k)*cosd(qM);  % coordinates of closest approach to center
                y1 = yCent - (rM+k)*sind(qM);  
                plot([x1 - L*cosd(qM+90), x1 + L*cosd(qM+90)], [y1 + L*sind(qM + 90), y1 - L*sind(qM + 90)], 'b-')
                if k == -dW;
                    tmp = (1/(2*dW+1))*improfile(subImg, [x1 - L*cosd(qM+90), x1 + L*cosd(qM+90)], [y1 + L*sind(qM + 90), y1 - L*sind(qM + 90)]);
                else
                    tmp = tmp + (1/(2*dW+1))*improfile(subImg, [x1 - L*cosd(qM+90), x1 + L*cosd(qM+90)], [y1 + L*sind(qM + 90), y1 - L*sind(qM + 90)]);
                end;
            end;
        % background subtraction
        for k = [-dW-3 -dW-2 dW+2 dW+3];
            x1 = xCent + (rM+k)*cosd(qM);  % coordinates of closest approach to center
            y1 = yCent - (rM+k)*sind(qM); 
            plot([x1 - L*cosd(qM+90), x1 + L*cosd(qM+90)], [y1 + L*sind(qM + 90), y1 - L*sind(qM + 90)], 'm-')
            tmp = tmp - (1/4)*improfile(subImg, [x1 - L*cosd(qM+90), x1 + L*cosd(qM+90)], [y1 + L*sind(qM + 90), y1 - L*sind(qM + 90)]);
        end;
        hold off
        
        
        
        profile{j} = tmp;
        subplot(1,3,3);
        plot((1:length(profile{j}))*pixSize, profile{j});
        xlabel('Position (\mum)')
        title('Line profile')
        pbaspect([1 1 1])
        
        
    end;
        
end;

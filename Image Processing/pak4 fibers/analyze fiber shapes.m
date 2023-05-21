% Make a pretty picture of one fiber

cd 'X:\Lab\Labmembers\Dingchang Lin\Data\HCBI\20210906_2h'
addpath(pwd);
cd '.\503_1&2'
ySize = 2048; xSize = 2048;
fName = '9_Unmixing.czi';
dat = bfopen(fName);  % open each Unmixed image.
colorImg = zeros(ySize, xSize, 3);
for j = 1:3;
    colorImg(:,:,j) = mat2gray(double(dat{1}{j,1}));  % make a 3-color image
end;
figure(1)
imshow2(colorImg)
dW = 45; % half-width in pixels of a fiber
[h, fiberImgs] = getFibers(colorImg, dW, 1); 
img = fiberImgs{1};
colorImg2 = zeros(size(img));
colorImg2(:,:,:) = ...
    reshape([0 .8 1], [1 1 3]).*img(:,:,1) + ...
    reshape([.7 .7 0], [1 1 3]).*img(:,:,2) + ...
    reshape([1 0 0], [1 1 3]).*img(:,:,3);
imshow2(colorImg2)
omeMeta = dat{1,4};
voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER); % in µm
voxelSizeXdouble = voxelSizeX.doubleValue();
hold all
dX = 5/voxelSizeXdouble;  % 5 um scalebar
plot([20 20+dX], [85 85], 'LineWidth', 3, 'Color', 'white'); hold off
saveas(gca, [fName(1:3) 'pretty pic.fig'])
saveas(gca, [fName(1:3) 'pretty pic.png'])



%%
% quantify the fluorescence traces

[~, ySize, ~] = size(img);
r = drawrectangle(gca);
yMin = r.Position(2);
yMax = yMin + r.Position(4);
yMin = round(max(yMin, 1));
yMax = round(min(yMax, ySize));
iPlot = squeeze(median(img(yMin:yMax,:,:), 1));
yLims(:,j) = [yMin; yMax];
save('9_allTraces.mat', 'iPlot', 'voxelSizeXdouble');

xaxis = (1:length(iPlot))*voxelSizeXdouble;
figure(2); clf
plot(xaxis, iPlot)
xlabel('Position (\mum)')
saveas(gca, '9_fiber_traces.fig')
saveas(gca, '9_fiber_traces.png')



for n = 1:nDir  % Combine all the fiber images for each condition into a single file
    n
    cd(dList(n).name)
    fList = dir('*fibers.mat');
    nFiles = length(fList);
    allFibers = [];
    for p = 1:nFiles;
        fName = fList(p).name;
        load(fName)
        allFibers = [allFibers, fiberImgs];
    end;
    nFib = length(allFibers);
    save('allTraces.mat', 'allFibers', '-append');
    cd ..
end;

% Decide which fibers are good
for n = 1:nDir
    n
    cd(dList(n).name)
    load('allTraces.mat');
    goodDat = zeros(nFib,1);
    for j = 1:nFib;
        subplot(2,1,1);
        img = allFibers{j};
        [ySize(j), xSize(j), ~] = size(img);
        imshow2(img); hold all
        plot([1, xSize(j)], [yLims(1,j), yLims(1,j)], 'w-.')
        plot([1, xSize(j)], [yLims(2,j), yLims(2,j)], 'w-.'); hold off
        subplot(2,1,2);
        plot(iPlot{j})
        colororder(diag([1 1 1]));
        title(num2str(j))
        prompt = 'Is the image ok? Y/N [Y]: ';  % Get user input on whether the image is ok
        str = input(prompt,'s');
        if isempty(str) || strcmp(str, 'Y') || strcmp(str, 'y')
            str = 'Y';
            goodDat(j) = 1;
        else
            goodDat(j) = 0;
        end;
    end;
    save('allTraces.mat','xSize', 'ySize', 'goodDat', '-append');
    cd ..
end;

colors = {'r', 'g', 'b'};

filt = diff(diff(fspecial('gaussian', [100 1], 4))); 
filt = filt/sum(filt.^2).^.5;
figure(80); plot(filt)
dyeIdx = [3 2 1 3 3];  % order of the dye transitions (blue, green, red, blue, blue) 
winW = 20; % window width for finding local maxima in second derivative.

for n = 1:nDir
    n
    cd(dList(n).name)
%     load('allImgs.mat');
    load('allTraces.mat');
    xSwitch = zeros(5, 2, nFib);  % coordinates of the switches in the three color channels
    for j = 1:nFib;
        if goodDat(j);
            clf
            subplot(2,1,1);
            img = allFibers{j};
            [ySize(j), xSize(j), ~] = size(img);
            imshow2(img);
            subplot(2,1,2);
            dat = iPlot{j};
            plot(mat2gray(dat)); hold all
            colororder(diag([1 1 1]));
            datSpp = imfilter(dat, filt, 'replicate');
            plot(1+mat2gray(datSpp))
            title(['Dir ' num2str(n) ' Fib ' num2str(j) ' of ' num2str(nFib)])
            for k = 1:5;
                roi = drawpoint;
                xMin = round(max(roi.Position(1)-winW, 1));
                xMax = round(min(roi.Position(1)+winW,xSize(j)));
                [~, idx] = max(datSpp(xMin:xMax,dyeIdx(k)));
                xSwitch(k,1,j) = idx+xMin-1;
                plot([xSwitch(k,1,j), xSwitch(k,1,j)], [0, 2], 'k-');
            end;
            prompt = 'Is there a second side (y/n) [y]: ';  % Get user input on whether the image is ok
            str = input(prompt,'s');
            if isempty(str) || strcmp(str, 'Y') || strcmp(str, 'y')
                for k = 1:5;
                    roi = drawpoint;
                    xMin = round(max(roi.Position(1)-winW, 1));
                    xMax = round(min(roi.Position(1)+winW,xSize(j)));
                    [~, idx] = max(datSpp(xMin:xMax,dyeIdx(k)));
                    xSwitch(k,2,j) = idx+xMin-1;
                    plot([xSwitch(k,2,j), xSwitch(k,2,j)], [0, 2], 'k-');
                end;
            end;
            hold off
            pause(2)
        end;
    end;
    save('allTraces.mat', 'xSwitch', 'goodDat', '-append');
    cd ..
end


%% Do some quality control on the data

% load('allTraces.mat');    
for j = 1:nFib;
    if goodDat(j);
        clf
        subplot(2,1,1);
        img = allFibers{j};
        [ySize(j), xSize(j), ~] = size(img);
        imshow2(img);
        subplot(2,1,2);
        dat = iPlot{j};
        plot(mat2gray(dat)); hold all
        colororder(diag([1 1 1]));
        title([' Fib ' num2str(j) ' of ' num2str(nFib)])
        for k = 1:5;
            plot([xSwitch(k,1,j), xSwitch(k,1,j)], [0, 1], 'k-');
        end;
        if max(xSwitch(:,2,j) > 0)
            for k = 1:5;
                plot([xSwitch(k,2,j), xSwitch(k,2,j)], [0, 1], 'k-');
            end;
        end;
        hold off
        pause
    end
end


%% Look at all the double-ended fibers
figure(16); clf
for j = 1:nFib;
    if (max(xSwitch(:,1,j)) > 0) && (max(xSwitch(:,2,j)) > 0);
        subplot(nFib, 1, j); 
        img = allFibers{j};
        [ySize, xSize] = size(img);
        imshow2(img(:,[1:xSwitch(1,2,j)+20 xSwitch(1,1,j)-20:end],:));
        text(xSize/4, ySize, num2str(j))
    end;
end;


%% normalize each fiber end by the first and last points
xNorm = zeros(size(xSwitch)); % fibers with normalized lengths
xNorm2 = zeros(size(xSwitch)); % each fiber has the longer end in k = 1, shorter end in k = 2;
longEnd = zeros(nFib, 1);  % 0 if only one end, 1 if k = 1 side is longer, 2 if k = 2 side is longer
for j = 1:nFib
    for k = 1:2;
        if max(xSwitch(:,k,j)) > 0;
            xNorm(:,k,j) = (xSwitch(:,k,j) - xSwitch(1,k,j))/(xSwitch(5,k,j) - xSwitch(1,k,j));
        end;
    end;
end;
for j = 1:nFib
    if abs((xSwitch(5,1,j) - xSwitch(1,1,j))) >= abs((xSwitch(5,2,j) - xSwitch(1,2,j)))
        xNorm2(:,1,j) = (xSwitch(:,1,j) - xSwitch(1,1,j))/(xSwitch(5,1,j) - xSwitch(1,1,j));
        xNorm2(:,2,j) = (xSwitch(:,2,j) - xSwitch(1,2,j))/(xSwitch(5,2,j) - xSwitch(1,2,j));
        longEnd(j) = 1;
    elseif abs((xSwitch(5,1,j) - xSwitch(1,1,j))) < abs((xSwitch(5,2,j) - xSwitch(1,2,j)))
        xNorm2(:,1,j) = (xSwitch(:,2,j) - xSwitch(1,2,j))/(xSwitch(5,2,j) - xSwitch(1,2,j));
        xNorm2(:,2,j) = (xSwitch(:,1,j) - xSwitch(1,1,j))/(xSwitch(5,1,j) - xSwitch(1,1,j));
        longEnd(j) = 2;
    end;
    if (xSwitch(1,1,j) == 0) || (xSwitch(1,2,j) == 0)
        longEnd(j) = 0;
    end;
end;
xNorm = xNorm2;  % use this switch if plotting long end vs short end below

% The actual order of dye addition was:
% JF503, JF669, JFX608, JF503, JF669
% The corresponding color palette is (RGB):
% (0 204 205), (255 0 0), (255 255 0)
% Plot the trajectories for the fibers with both ends labeled
figure(14); clf
c = 1;
for j = 1:nFib
    plot(0, 0, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); hold all
    plot(1, 1, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); 
    if (max(xNorm(:,1,j)) > 0) && (max(xNorm(:,2,j)) > 0);
        plot(xNorm(:,1,j), xNorm(:,2,j), '-', 'Color', [.5 .5 .5])
        plot(xNorm(2,1,j), xNorm(2,2,j), 'o', 'MarkerFaceColor', [.7 .7 0], 'MarkerEdgeColor', [.7 .7 0])
        plot(xNorm(3,1,j), xNorm(3,2,j), 'o', 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'c'); 
        plot(xNorm(4,1,j), xNorm(4,2,j), 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
        c = c + 1;
    end;
end;
hold off
xlabel('Faster end position')
ylabel('Slower end position')
daspect([1 1 1])
saveas(gca, 'Transitions on the two ends by speed.fig')
saveas(gca, 'Transitions on the two ends by speed.png')

%% Analyze the mean and std and covariance of the transition measurements

% Pull the data from the figure
open('X:\Lab\Labmembers\Dingchang Lin\Data\HCBI\20210906_2h\503_1&2\Transitions on the two ends.fig');
h = gcf;
xDat = zeros(5,17);
yDat = zeros(5,17);
c = 1;
for j = 1:length(h.Children.Children);
    d = h.Children.Children(j).XData;
    d2 = h.Children.Children(j).YData;
    if length(d) == 5;
        xDat(:,c) = d';
        yDat(:,c) = d2';
        c = c + 1;
    end;
end;
figure(2); clf
plot(xDat, yDat, 'o-')
mean(xDat,2)  % 0 .233 .475 .737 1
std(xDat, [], 2)  % 0 .0334 .043 .0382 0
mean(std(xDat(2:4,:), [], 2))*60*8    % 18.33 min variability
mean(yDat,2) % 0 .217 .429 .678 1
std(yDat, [], 2) % 0 .0496 .0547 .0506 0
mean(std(yDat(2:4,:), [], 2))*60*8  % 24.77 min variability
mean((xDat+yDat)/2,2)  % 0 .225 .452 .708 1
std((xDat+yDat)/2, [], 2) % 0 .0363 .0369 .0369

% calculate the standard deviation in timing two different ways
sX1 = std([xDat, yDat], [], 2)*8*60   % 0 20.3 25.8 25.6
mean(sX1(2:4)) % 23.9 minutes  (don't include first and last data points in the avg)
sX2 = [std(xDat, [], 2); std(yDat, [], 2)]*8*60
mean(sX2([2:4 7:9]))  % 21.55 minutes

% calculate the covariance in the timing errors
dxDat = xDat - mean(xDat,2);
dyDat = yDat - mean(yDat,2);

mean(dxDat.*dyDat,2)./(std(dxDat,[], 2).*std(dyDat, [], 2))  %NaN 0.484 123 0.348 NaN
nanmean(mean(dxDat.*dyDat,2)./(std(dxDat,[], 2).*std(dyDat, [], 2)))  % 0.318



%% Look at the non-normalized growth rates on the two ends
c = 1; L1 = []; L2 = [];
for j = 1:nFib
    if (max(xSwitch(:,1,j)) > 0) && (max(xSwitch(:,2,j)) > 0);
        L1(c) = abs(xSwitch(5,1,j)-xSwitch(1,1,j))*allFibPixDim{1}(j);
        L2(c) = abs(xSwitch(5,2,j)-xSwitch(1,2,j))*allFibPixDim{1}(j);
        c = c + 1;
    end;
end;
figure(15); clf
plot(max(L1, L2), min(L1, L2), 'bo') 
xlabel('L1 (\mum)')
ylabel('L2 (\mum)')
xlim([0 22])
ylim([0 22])
saveas(gca, 'Lengths of the two ends.fig')
saveas(gca, 'Lengths of the two ends.png')


figure(16); clf
tau = 8; % hours between start and finish of the experiment
plot(max(L1, L2)/tau, min(L1, L2)/tau, 'bo') 
xlabel('End 1 growth rate (\mum/h)')
ylabel('End 2 growth rate (\mum/h)')
fastMean = mean(max(L1, L2)/tau)  % 1.88 um/h
fastStd = std(max(L1, L2)/tau)  % 0.59 um/h
slowMean = mean(min(L1, L2)/tau)  % 1.15 um/h
slowStd = std(min(L1, L2)/tau)  % 0.51 um/h
overallMean = mean([L1 L2]/tau) % 1.52 um/h
overallStd = std([L1 L2]/tau)  % 0.66 um/h
xlim([0 3])
ylim([0 3])
saveas(gca, 'Growth rates of the two ends.fig')
saveas(gca, 'Growth rates of the two ends.png')



%% Interpolate all the fibers to the same x axis between transitions 1 and 5, so they can be compared and
% averaged.
nInterp = 400;  % number of points to interpolate between first and third dye transitions
xInterp = -200:600; % Capture most of the filaments

% cd(dList(n).name)
% load('allTraces.mat');
goodDat = logical(goodDat);
iPlot2 = zeros(length(xInterp), 3, 2, nFib);  % 3 colors, 2 ends per fiber
for j = 1:nFib;
    if goodDat(j);
        for k = 1:2;
            if max(xSwitch(:,k,j) > 0)
                x1 = xSwitch(1,k,j);
                x5 = xSwitch(5,k,j);
                dx = nInterp/(x5-x1);
                xMeas = ((1:length(iPlot{j})) - x1)*dx;
                iPlot2(:,:,k,j) = mat2gray(interp1(xMeas, iPlot{j}, xInterp, 'linear', 0));
            else
                iPlot2(:,:,k,j) = NaN;
            end;
        end;
    else;
        iPlot2(:,:,:,j) = NaN;
    end;
end;
% cd ..
save('allTraces.mat', 'iPlot2', 'xInterp', 'nInterp', 'xNorm', 'xNorm2', '-append');

for k = 1:2;
    figure(18+k); clf
    subplot(3,1,1)
    plot(xInterp/nInterp, squeeze(iPlot2(:,1,k,:)), 'Color', [.5 1 1]); hold all
    plot(xInterp/nInterp, nanmean(iPlot2(:,1,k,:),4),  'LineWidth', 2, 'Color', [0 .6 .6]); hold off
    xlim([-.1 1.1])
    subplot(3,1,2)
    plot(xInterp/nInterp, squeeze(iPlot2(:,3,k,:)), 'Color', [1 .5 .5]); hold all
    plot(xInterp/nInterp, nanmean(iPlot2(:,3,k,:),4),  'LineWidth', 2, 'Color', [.6 0 0]); hold off
    xlim([-.1 1.1])
    subplot(3,1,3)
    plot(xInterp/nInterp, squeeze(iPlot2(:,2,k,:)), 'Color', [.85 .85 0.35]); hold all
    plot(xInterp/nInterp, nanmean(iPlot2(:,2,k,:),4),  'LineWidth', 2, 'Color', [.6 .6 0]); hold off
    xlim([-.1 1.1])
%     saveas(gca, [num2str(k) 'NormalizedProfiles.fig'])
%     saveas(gca, [num2str(k) 'NormalizedProfiles.png'])
end;

% Make a plot showing the average profiles for the long and short ends

longAvg = zeros(length(xInterp),3);
shortAvg = zeros(length(xInterp),3);
for j = 1:nFib;
    if longEnd(j) == 1;
        longAvg = longAvg + iPlot2(:,:,1,j);
        shortAvg = shortAvg + iPlot2(:,:,2,j);
    elseif longEnd(j) == 2;
        longAvg = longAvg + iPlot2(:,:,2,j);
        shortAvg = shortAvg + iPlot2(:,:,1,j);
    end;
end;
longAvg = mat2gray(longAvg);
shortAvg = mat2gray(shortAvg);
% colOrder = {'b', 'y', 'r'};
colMat = [0 .6 .6; .6 .6 0; .6 0 0];
subplot(2,1,1)
for j = 1:3;
    plot(xInterp(201:600)/nInterp, mat2gray(longAvg(201:600,j)), 'Color', colMat(j,:)); hold all
end
title('Long end'); hold off
subplot(2,1,2)
for j = 1:3;
    plot(xInterp(201:600)/nInterp, mat2gray(shortAvg(201:600,j)), 'Color', colMat(j,:)); hold all
end
title('Short end')
hold off
saveas(gca, 'NormalizedProfiles by length.fig')
saveas(gca, 'NormalizedProfiles by length.png')



for j = 1:nFib;
    subplot(4,6,j)
    colororder(diag([1 1 1]));
    plot(xInterp/nInterp, squeeze(iPlot2(:,:,1,j))); hold on
    plot(xInterp/nInterp, squeeze(iPlot2(:,:,2,j))); hold off
    title(num2str(j))
end;
saveas(gca, 'All fiber end pairs.fig')
saveas(gca, 'All fiber end pairs.png')
    
%% Show a particular fiber comparing the traces on the two ends
j = 10;
% img1 = mat2gray(allFibers{j}(:,end-430:end-30,:));
% img2 = mat2gray(allFibers{j}(:,411:-1:11,:));

img = fiberImgs{1};  % this is from the file 18_Unmixing.czi
img1 = mat2gray(img(:,end-430:end-30,:)); 
img2 = mat2gray(img(:,401:-1:1,:));

img12 = [img1; img2];
colorImg = zeros(size(img12));
colorImg(:,:,:) = ...
    reshape([0 .8 1], [1 1 3]).*img12(:,:,1) + ...
    reshape([.7 .7 0], [1 1 3]).*img12(:,:,2) + ...
    reshape([1 0 0], [1 1 3]).*img12(:,:,3);

figure(20); clf
subplot(3,1,1)
imshow2(colorImg(1:91,:,:)); 
subplot(3,1,2)
imshow2(colorImg(92:end,:,:)); 
hold all
dX = 5/allFibPixDim{1}(j);  % 5 um scalebar
plot([20 20+dX], [85 85], 'LineWidth', 3, 'Color', 'white'); hold off
subplot(3,1,3)
colororder([0 .4 .8;  .7 .7 0; 1 0 0]);
plot(xInterp/nInterp, mat2gray(squeeze(iPlot2(:,:,1,j))), 'LineWidth', 1.5); hold on
plot(1.05*xInterp/nInterp, mat2gray(squeeze(iPlot2(:,:,2,j))), '-.', 'LineWidth', 1.5); hold off
xlim([-.2 1.25])
saveas(gca, 'CompareEndsExample.fig')
saveas(gca, 'CompareEndsExample.png')

%% Show some beautiful images of fibers.
% files 11, 12, 15 are nice

%%
% Quantify the steepness of the transitions by ligning up the fibers at
% each of the key transition
% Fiber 15,1 is backwards
plot(squeeze(xSwitch(:,1,:)))
iPlot{15} = flipud(iPlot{15});
plot(iPlot{15})
xSwitch(:,1,15) = length(iPlot{15}) + 1 - xSwitch(:,1,15);
plot(squeeze(xSwitch(:,1,:)))

% Quantify the steepness of each dye transition
figure(16); clf
L = length(iPlot2);
subplot(2,1,1)
plot(nanmean(iPlot2(:,:,1,:),4))
subplot(2,1,2)
plot(nanmean(iPlot2(:,:,2,:),4))
figure(17); clf
alignedPlots = zeros(4*nInterp, 4);
colororder(diag([1 1 1]));
xSwitch2 = (xSwitch - xSwitch(1,:,:))./(xSwitch(5,:,:) - xSwitch(1,:,:));
fibOrder = [3 2 1 3];
colOrder = {'r-', 'y-', 'b-', 'r-'};
for k = setdiff(1:nFib, [2, 3, 7, 9, 13, 14, 18]);
    for m = 1:2;
        if xSwitch(1,m,k) > 0;
            for j = 1:4;
    %           xAx = (1:L)/nInterp - xSwitch2(j,m,k);
    %           plot(xAx, j + mat2gray(iPlot2(:,fibOrder(j),m,k)), colOrder{fibOrder(j)}); hold all
                idx = (0:L-1) - round(nInterp*xSwitch2(j,m,k) - nInterp*nanmean(xSwitch2(j,m,:),3)) + nInterp;
                alignedPlots(idx,j) = alignedPlots(idx,j) + mat2gray(iPlot2(:,fibOrder(j),m,k));
                plot(idx, j + mat2gray(iPlot2(:,fibOrder(j),m,k)), colOrder{fibOrder(j)}); hold all
            end;
%             title(['Fiber ' num2str(k) ' end ' num2str(m)])
%             pause
        end;
    end;
end;
alignedPlots = alignedPlots(451:1200,:);

figure(18); clf
% xAx = (1:length(alignedPlots))/nInterp -.3;
xAx = (1:length(alignedPlots));
subplot(2,1,1)
plot(xAx, mat2gray(alignedPlots(:,2)), 'Color', [.5 .5 0]); hold all
plot(xAx, mat2gray(alignedPlots(:,3)), 'c-');
plot(xAx, mat2gray(alignedPlots(:,4)), 'r-'); hold off;

subplot(2,1,2)
d2 = diff(diff(alignedPlots));
d2S = imfilter(d2, fspecial('Gaussian', [100, 1], 4), 'replicate');
plot(xAx(200:300), d2S(200:300,2)/max(d2S(200:300,2)), 'Color', [.5 .5 0]); hold all
plot(xAx(275:375), d2S(275:375,3)/max(d2S(275:375,3)), 'c-');
plot(xAx(375:475), d2S(375:475,4)/max(d2S(375:475,4)), 'r-'); hold off
xlim([0 800])
%%


for j = 1:nDir;  % Get the number of fibers in each directory.
    nFib(j) = size(iPlot2{j}, 3);
    nGoodFib(j) = sum(allSwitch{j} > 0);
end;
allFibs = zeros(length(xInterp),3,sum(nFib));

c = 1;  % put all the fiber profiles into a single array.
allSwitchMat = [];
allFibs = [];
for j = 1:nDir;
    allFibs(:,:,c:c+nFib(j)-1) = iPlot2{j};
    allSwitchMat = [allSwitchMat; allSwitch{j}];
    c = c+nFib(j);
end;
% Get rid of the NaN lines (bad fibers)
badDat = squeeze(isnan(max(allFibs(:,1,:), [], 1)));
allFibs(:,:,badDat) = [];
allFibs = permute(allFibs,[1,3,2]);

%%
% The actual order of dye addition was:
% JF525, JF669, JFX608, JF525
% The corresponding color palette is (RGB):
% (0 204 205), (255 0 0), (255 255 0)
colorImg = zeros(size(allFibs));
colorImg(:,:,2) = allFibs(:,:,1)*204/255;
colorImg(:,:,3) = allFibs(:,:,1)*205/255;
colorImg(:,:,1) = allFibs(:,:,3);
colorImg(:,:,1) = colorImg(:,:,1) + allFibs(:,:,2);
colorImg(:,:,2) = colorImg(:,:,2) + allFibs(:,:,2);

figure(1); clf
subplot(5,1,1:4)
imshow2(flipud(colorImg))
daspect([1 10 1])
hold all;
plot(1:length(allSwitchMat), 3001-(1000+nInterp*allSwitchMat), 'w-', 'LineWidth', 2)
plot(1:length(allSwitchMat), 3001-(1000), 'k.-', 'LineWidth',1)
plot(1:length(allSwitchMat), 3001-(1000+nInterp), 'k.-', 'LineWidth',1)
hold off
ylim([1200 2200])
subplot(5,1,5)
nVals = cumsum([1 nGoodFib]);
for j = 1:nDir;
    plot([nVals(j) nVals(j+1)-1], [tau(j), tau(j)], 'k-'); hold all
end;
ylim([0 10])
xlim([0 sum(nGoodFib)])
daspect([2 1 1])
hold off
xlabel('Fiber number')
ylabel('Time (h)')
saveas(gca, 'AllFiberImages.fig')
saveas(gca, 'AllFiberImages.png')



plot(nanmean(allFibs(:,1,:),3))

avgI = zeros(length(xInterp),3,nDir);
for j = 1:nDir;
    avgI(:,:,j) = nanmean(iPlot2{j},3);
end;

colors = generateColorSpecLocal(3*nDir);
figure(1); clf
subplot(3,1,1);
set(gca, 'ColorOrder', colors,  'NextPlot', 'replacechildren');
plot(xInterp/nInterp, mat2gray(squeeze(avgI(:,1,:))), 'LineWidth', 2)
xlim([-200, 800]/nInterp)
title('JF525');
subplot(3,1,2);
set(gca, 'ColorOrder', colors(2*nDir+1:end,:),  'NextPlot', 'replacechildren');
plot(xInterp/nInterp, mat2gray(squeeze(avgI(:,3,:))), 'LineWidth', 2)
xlim([-200, 800]/nInterp)
title('JF669');
subplot(3,1,3);
set(gca, 'ColorOrder', colors(nDir+1:end,:),  'NextPlot', 'replacechildren');
plot(xInterp/nInterp, mat2gray(squeeze(avgI(:,2,:))), 'LineWidth', 2)
xlim([-200, 800]/nInterp)
title('JFX608');
xlabel('Fractional position along filament')
legend('2 h', '4 h', '6 h', '7 h', '8 h', '9 h')
saveas(gca, 'Avg color transitions.fig')
saveas(gca, 'Avg color transitions.png')

for j = 1:nDir;
    subplot(nDir,1,j);
    imshow(shiftdim(iPlot2{j},2))
    daspect([1 .1 1])
end;
    


for j = 1:nFib
    subplot(2,1,1)
    imshow2(allFibers{j})
    subplot(2,1,2)
    plot(iPlot{j})
    colororder(diag([1 1 1]));
    title(num2str(j))
    pause
end;


for j = 1:nFib;
    yMin = yLims(1,j);
    yMax = yLims(2,j);
    iPlot{j} = squeeze(median(allFibers{j}(yMin:yMax,:,:), 1));
end;
        
        



dat = bfopen('1_Unmixing.czi');
[ySize, xSize] = size(dat{1}{1});
nChan = 4; % number of color channels
mov = zeros(ySize, xSize, nChan);
for j = 1:nChan;
    mov(:,:,j) = double(dat{1}{j,1});
end

for j = 1:nChan;
    subplot(2,2,j);
    imshow2(mov(:,:,j), [])
end;

colorImg = zeros(ySize, xSize, 3);
for j = 1:3;
    colorImg(:,:,j) = mat2gray(mov(:,:,j));
end;

h = getFibers(colorImg);


for j = 1:nROI;
    traces{j} = squeeze(median(subImg{j}, 1));
end;

figure(5); clf
for j = 1:nROI;
    subplot(3,3,j)
    plot(traces{j})
end;

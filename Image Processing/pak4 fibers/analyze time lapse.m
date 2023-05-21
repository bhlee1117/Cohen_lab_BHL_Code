%% Load the video
vid = VideoReader('13-26hr_10min_interval_6fps.avi');
nFrames = vid.FrameRate*vid.Duration;
mov = vid.read;

% Pull out the scalebar from the red channel
scalebar = mov(1150:end, 1100:end, 1, 1);
imshow2(scalebar)
pixSize = 50/(143-16); % microns

mov = squeeze(mov(:,:,2,:)); % All the movie data is in the green channel
maxImg = max(mov, [], 3);  % Take a maximum projection to see where each fiber goes\
figure(1); clf
moviefixsc(mov);


c = 1;
allROI = [];  % cell array, hand-drawn boundary
allbBox = [];  % rectangular bounding box
allProfile = [];  % cell array, intensity profile along fiber
allIntens = [];  % integrated intensity along fiber
allX0 = []; % coordinates of fiber closest to image center
allY0 = [];
allQ = []; % fiber angle
allRM = []; % projected distance of closest approach of fiber to center
theta = [0:179];   % Angles for Radon transform

%% uncomment to save an example of the analysis
% vid = VideoWriter('TrackingDemo.avi', 'Uncompressed AVI');
% vid.FrameRate = 10;
% open(vid);

%%
% Run this block once per fiber
[subMov, bBox, roi] = getROI(mov, maxImg);  % draw an ROI around a fiber.
for c = 1:12
    roi = allROI{c};
    [subMov, bBox] = applyROI(mov, maxImg, roi);
    [ySize, xSize] = size(subMov(:,:,1));
    yCent = ySize/2;  % Find the center of the sub-movie
    xCent = xSize/2;

    L = sqrt(xSize^2 + ySize^2)/2;  % Longest-possible half-line length
    profile = [];
    subMov2 = zeros(size(subMov));
    for j = 1:nFrames
        subMov2(:,:,j) = subMov(:,:,j) - .5*medfilt2(subMov(:,:,j), [10 10], 'symmetric');

        [radonMov, xp] = radon(subMov2(:,:,j), theta);
        [maxI, idx] = max(radonMov, [], 'all', 'linear');  % Find the linear index of max value in Radon xform
        allIntens(j,c) = maxI;
        [xM, thetaM] = ind2sub(size(radonMov), idx);  % convert to subscripts
        rM = xp(xM);  % closest approach of the line to the center of the image
        qM = theta(thetaM);  % angle of the line
        subplot(1,3,1);  % Show the maximum on the radon image
        imagesc(theta, xp, radonMov); colormap('hot'); hold all
        plot(qM, rM, 'cx'); hold off
        xlabel('Angle (deg)')
        ylabel('Position (pix)')
        pbaspect([1 1 1])
        title('Radon transform')
        subplot(1,3,2);  % Show the calculated line
        imagesc(subMov(:,:,j)); hold all
        axis('off')
        daspect([1 1 1])
        x0 = xCent + rM*cosd(qM);  % coordinates of closest approach to center
        y0 = yCent - rM*sind(qM);
        allX0(j,c) = x0;
        allY0(j,c) = y0;
        allQ(j,c) = qM + 90;
        allRM(j,c) = rM;
    %     plot(xCent, yCent, 'co');
    %     plot([xCent, x0], [yCent, y0], 'b.-')
        plot([x0 - L*cosd(qM+90), x0 + L*cosd(qM+90)], [y0 + L*sind(qM + 90), y0 - L*sind(qM + 90)], 'g-', 'LineWidth', 1);
        title(['Frame = ' num2str(j)])

        % Extract the profile of the line, averaging over a stripe of width
        % 2*dW+1
        dW = 2; % number of pixels to average on either side of the line;
        tmp = [];
        for k = -dW:dW;
            x1 = xCent + (rM+k)*cosd(qM);  % coordinates of closest approach to center
            y1 = yCent - (rM+k)*sind(qM);  
            plot([x1 - L*cosd(qM+90), x1 + L*cosd(qM+90)], [y1 + L*sind(qM + 90), y1 - L*sind(qM + 90)], 'b-')
            if k == -dW;
                tmp = (1/(2*dW+1))*improfile(subMov(:,:,j), [x1 - L*cosd(qM+90), x1 + L*cosd(qM+90)], [y1 + L*sind(qM + 90), y1 - L*sind(qM + 90)]);
            else
                tmp = tmp + (1/(2*dW+1))*improfile(subMov(:,:,j), [x1 - L*cosd(qM+90), x1 + L*cosd(qM+90)], [y1 + L*sind(qM + 90), y1 - L*sind(qM + 90)]);
            end;
        end;
        % background subtraction
        for k = [-dW-3 -dW-2 dW+2 dW+3];
            x1 = xCent + (rM+k)*cosd(qM);  % coordinates of closest approach to center
            y1 = yCent - (rM+k)*sind(qM); 
            plot([x1 - L*cosd(qM+90), x1 + L*cosd(qM+90)], [y1 + L*sind(qM + 90), y1 - L*sind(qM + 90)], 'm-')
            tmp = tmp - (1/4)*improfile(subMov(:,:,j), [x1 - L*cosd(qM+90), x1 + L*cosd(qM+90)], [y1 + L*sind(qM + 90), y1 - L*sind(qM + 90)]);
        end;
        hold off
        profile{j} = tmp;

        subplot(1,3,3);
        plot((1:length(profile{j}))*pixSize, profile{j});
        xlabel('Position (\mum)')
        title('Line profile')
        pbaspect([1 1 1])
        pause(.1)
    %     frame = getframe(gcf);
    %     writeVideo(vid, frame);
    end;
    allbBox(:,c) = bBox;
    allROI{c} = roi;
    allProfile{c} = profile;
    % close(vid)
end;
c = c+1;
%%
save('TimeLapseDat3.mat', 'allbBox', 'allIntens', 'allProfile', 'allROI', 'allX0', 'allY0', 'allQ', 'allRM', 'nFrames', 'pixSize', 'scalebar', 'theta');

c = c - 1;

% Look at all the selected ROIs
figure(12); clf
imshow2(maxImg, []); hold all
for j = 1:12;
    plot([allROI{j}(end,1); allROI{j}(:,1)], [allROI{j}(end,2); allROI{j}(:,2)])
    text(mean(allROI{j}(:,1)), mean(allROI{j}(:,2)), num2str(j), 'FontSize', 15, 'Color', 'blue')
end;
title('Analyzed fibers')
hold off
saveas(gca, 'Fiber map.fig')
saveas(gca, 'Fiber map.png')

% Plot the kymographs along the fiber axes
for j = 1:c;
    for k = 1:nFrames;
        fibL(k) = length(allProfile{j}{k});
    end;
    maxL = max(fibL);
    kymo = zeros(maxL, nFrames);
    for k = 1:nFrames;
        kymo(1:fibL(k),k) = allProfile{j}{k};
    end;
    figure(1)
    subplot(3,4,j)
    imshow2(kymo, [])
    title(num2str(j))
end;
saveas(gca, 'allKymos.fig')
saveas(gca, 'allKymos.png')

% Angle as a function of time
figure(2); clf
for j = 1:c;
    subplot(3,4,j)
    plot(allQ(:,j));
end;
    

% Now analyze the trajectories one at a time.
len = [];
c = 1;
figure(5); clf
subplot(2,2,1)
plot(allIntens(:,c));
title('Intensity')
subplot(2,2,2)
plot(allQ(:,c))
title('Angle')
thresh = 15;
subplot(2,2,3:4)
for j = 1:nFrames;
    plot(allProfile{c}{j})
    x1(j) = find(allProfile{c}{j} > thresh, 1, 'first');
    x2(j) = find(allProfile{c}{j} > thresh, 1, 'last');
    len(j,c) = x2(j) - x1(j);
    hold all
    plot(x1(j):x2(j), allProfile{c}{j}(x1(j):x2(j)), 'r-', 'LineWidth', 2);
    hold off
    title(num2str(j));
    ylim([-20 200])
    pause(.1)
end;
figure(6); clf
tau = (0:(nFrames-1))/6; % time in h, frame every 10 min
plot(tau, pixSize*len)
xlabel('Time (h)')
ylabel('Length (\mum)')

kymo = zeros(300,nFrames);
for j = 1:nFrames
    trace = allProfile{c}{j};
    ctr = round(mean([x2(j) x1(j)]));
    halfL = min([length(trace)-ctr ctr])-1;
    kymo((151-halfL):(151+halfL),j) = trace(ctr-halfL:ctr+halfL);
end;
figure(11); clf
kymo(isnan(kymo(:))) = 0;
kymo2 = zeros(size(kymo,1), size(kymo,2), 3);
kymo2(:,:,2) = mat2gray(kymo);
hold all
imshow2(kymo2, [])
% colormap('winter')
% caxis([-5 175])
daspect([1 3 1])
barX0 = 5;
barY0 = 170;
barDX = 12; % 2 h
barDY = 10/pixSize;  % 10 microns
plot([barX0 barX0 barX0+barDX], [barY0-barDY barY0 barY0], 'LineWidth', 2, 'Color', 'w')
hold off
title(['Fiber ' num2str(c) ', tScale = 2 h yScale = 10 \mum'])
saveas(gca, ['Fiber_' num2str(c) '_kymo.fig'])
saveas(gca, ['Fiber_' num2str(c) '_kymo.png'])

%%
c = 2;
figure(4); clf
bBox = allbBox(:,c);
subMov = mov(bBox(2):(bBox(2)+bBox(4)),bBox(1):(bBox(1)+bBox(3)),:);
moviefixsc(subMov);
goodF = 28:72;

figure(5); clf
subplot(2,3,1);
imshow2(maxImg, []); hold all
plot(allROI{c}(:,1), allROI{c}(:,2))
hold off
subplot(2,3,2)
plot(allIntens(:,c));
title('Intensity')
subplot(2,3,3)
plot(allQ(:,c))
title('Angle')
thresh = 8;
subplot(2,3,4:6)
for j = 1:nFrames;
    plot(allProfile{c}{j})
    ylim([-20 30])
    if ismember(j, goodF);
        x1 = find(allProfile{c}{j} > thresh, 1, 'first');
        x2 = find(allProfile{c}{j} > thresh, 1, 'last');
        if ~(isempty(x2) || isempty(x1));
            len(j,c) = x2 - x1;
        else
            len(j,c) = NaN;
        end;
        hold all
        plot(x1:x2, allProfile{c}{j}(x1:x2), 'r-', 'LineWidth', 2);
        hold off
    else
        len(j,c) = NaN;
    end;
    title(num2str(j));
    pause(.1)
end;
figure(6); clf
tau = (0:(nFrames-1))/6; % time in h, frame every 10 min
plot(tau, pixSize*len(:,1:c))
xlabel('Time (h)')
ylabel('Length (\mum)')

%%
c = 3;
figure(4); clf
bBox = allbBox(:,c);
subMov = mov(bBox(2):(bBox(2)+bBox(4)),bBox(1):(bBox(1)+bBox(3)),:);
moviefixsc(subMov);
goodF = 1:78;

figure(5); clf
subplot(2,3,1);
imshow2(maxImg, []); hold all
plot(allROI{c}(:,1), allROI{c}(:,2))
hold off
subplot(2,3,2)
plot(allIntens(:,c));
title('Intensity')
subplot(2,3,3)
plot(allQ(:,c))
title('Angle')
thresh = 12;
subplot(2,3,4:6)
for j = 1:nFrames;
    plot(allProfile{c}{j})
    ylim([-20 100])
    if ismember(j, goodF);
        x1 = find(allProfile{c}{j} > thresh, 1, 'first');
        x2 = find(allProfile{c}{j} > thresh, 1, 'last');
        if ~(isempty(x2) || isempty(x1));
            len(j,c) = x2 - x1;
        else
            len(j,c) = NaN;
        end;
        hold all
        plot(x1:x2, allProfile{c}{j}(x1:x2), 'r-', 'LineWidth', 2);
        hold off
    else
        len(j,c) = NaN;
    end;
    title(num2str(j));
    pause(.1)
end;
figure(6); clf
tau = (0:(nFrames-1))/6; % time in h, frame every 10 min
plot(pixSize*len(:,1:c))
xlabel('Time (h)')
ylabel('Length (\mum)')

%%
c = 4;
figure(4); clf
bBox = allbBox(:,c);
subMov = mov(bBox(2):(bBox(2)+bBox(4)),bBox(1):(bBox(1)+bBox(3)),:);
moviefixsc(subMov);
% goodF = [9:25 35:64];
goodF = 1:78;

figure(5); clf
subplot(2,3,1);
imshow2(maxImg, []); hold all
plot(allROI{c}(:,1), allROI{c}(:,2))
hold off
subplot(2,3,2)
plot(allIntens(:,c));
title('Intensity')
subplot(2,3,3)
plot(allQ(:,c))
title('Angle')
thresh1 = 15;
thresh2 = 20;
subplot(2,3,4:6)
for j = 1:nFrames;
    trace = allProfile{c}{j};
    trace = trace - prctile(trace, 10);
    traceS = medfilt1(trace, 10);
    plot(trace); hold all
    plot(traceS);
    ylim([-20 75])
    if ismember(j, goodF);
        l1 = find(traceS > thresh1, 1, 'first');
        l2 = find(traceS > thresh2, 1, 'last');
        if ~(isempty(l1) || isempty(l2));
            x1(j) = l1; 
            x2(j) = l2;
            len(j,c) = x2(j) - x1(j);
            hold all
            plot(x1(j):x2(j), trace(x1(j):x2(j)), 'r-', 'LineWidth', 2);
            hold off
        else
            len(j,c) = NaN;
            x1(j) = NaN; x2(j) = NaN;
        end;
    else
        len(j,c) = NaN;
    end;
    title(num2str(j));
    pause(.1)
end;

kymo = zeros(400,nFrames);
for j = 1:nFrames
    trace = allProfile{c}{j};
    if ~isnan(x2(j));
        ctr(j) = round(mean([x2(j) x1(j)]));
    else
        ctr(j) = round(length(trace)/2);
    end;
end;

ctr2 = ctr;
ctr2(27) = ctr(27)+30;
ctr2(28:29) = ctr(28:29)+30;
for j = 1:nFrames;
    trace = allProfile{c}{j};
    halfL = min([length(trace)-ctr2(j) ctr2(j)])-1;
    kymo((201-halfL):(201+halfL),j) = trace(ctr2(j)-halfL:ctr2(j)+halfL);
end;

figure(11); clf
kymo(isnan(kymo(:))) = 0;
kymo2 = zeros(size(kymo,1), size(kymo,2), 3);
kymo2(:,:,2) = 1.2*kymo/max(kymo(:)) + 0.1;
hold all
imshow2(kymo2, [])
daspect([1 5 1])
barX0 = 10;
barY0 = 300;
barDX = 6; % 1 h
barDY = 10/pixSize;  % 10 microns
plot([barX0 barX0 barX0+barDX], [barY0-barDY barY0 barY0], 'LineWidth', 2, 'Color', 'w')
hold off
title(['Fiber ' num2str(c) ', tScale = 1 h yScale = 10 \mum'])
xlim([0 68])
saveas(gca, ['Fiber_' num2str(c) '_kymo.fig'])
saveas(gca, ['Fiber_' num2str(c) '_kymo.png'])

figure(12); clf
dW = 25;
trace = allProfile{c}{1};
L = length(trace)/2;
profiles = zeros(2*L+1,2*dW+1, nFrames);
for j = 1:nFrames;
    x0 = allX0(j,c);
    y0 = allY0(j,c);
    qM = allQ(j,c);
    rM = ((x0-xCent)^2 + (y0 - yCent)^2).^.5;
    trace = allProfile{c}{j};
    imshow2(subMov(:,:,j)); hold all
    plot([x0 - L*cosd(qM), x0 + L*cosd(qM)], [y0 + L*sind(qM), y0 - L*sind(qM)], 'g-', 'LineWidth', 1);   
    for k = (-dW):dW;
        x1 = xCent + (rM+k)*cosd(qM+90);  % coordinates of closest approach to center
        y1 = yCent - (rM+k)*sind(qM+90);  
        tmp = improfile(subMov(:,:,j), [x1 - L*cosd(qM), x1 + L*cosd(qM)], [y1 + L*sind(qM), y1 - L*sind(qM)]);
        profiles(1:length(tmp),k+dW+1,j) = tmp;
    end;
    pause(.01)
end;
profiles(isnan(profiles(:))) = 0;

profMax = max(profiles, [], 3);
moviefixsc(profiles);

[~, intens] = clicky(profiles, profMax);
figure(13); clf
plot(intens(:,1) + intens(:,2));
clf
for j = 1:nFrames;
    imshow2(profiles(:,:,j), []);
    title(num2str(j))
    pause
end;


clicky(profiles(:,:,1:40));

%% 
c = 7;
bBox = allbBox(:,c);
subMov = mov(bBox(2):(bBox(2)+bBox(4)),bBox(1):(bBox(1)+bBox(3)),:);
figure(13); clf
[ySize, xSize, ~] = size(subMov);
yCent = ySize/2; xCent = xSize/2;
dW = 25;
trace = allProfile{c}{1};
L = length(trace)/2;
profiles = zeros(2*L+1,2*dW+1, nFrames);  % Calculate the profiles along many lines parallel to the fiber
for j = 1:nFrames;
    x0 = allX0(j,c);
    y0 = allY0(j,c);
    qM = allQ(j,c);
    rM = allRM(j,c);
    trace = allProfile{c}{j};
    imshow2(subMov(:,:,j)); hold all
    plot([x0 - L*cosd(qM), x0 + L*cosd(qM)], [y0 + L*sind(qM), y0 - L*sind(qM)], 'g-', 'LineWidth', 1);   
    for k = (-dW):dW;
        x1 = xCent + (rM+k)*cosd(qM-90);  % coordinates of closest approach to center
        y1 = yCent - (rM+k)*sind(qM-90);
%         plot([xCent x1], [yCent y1], 'c-')
        plot([x1 - L*cosd(qM), x1 + L*cosd(qM)], [y1 + L*sind(qM), y1 - L*sind(qM)], 'b-');
        tmp = improfile(subMov(:,:,j), [x1 - L*cosd(qM), x1 + L*cosd(qM)], [y1 + L*sind(qM), y1 - L*sind(qM)]);
        profiles(1:length(tmp),k+dW+1,j) = tmp;
    end;
    hold off
    pause(.1)
end;
profiles(isnan(profiles(:))) = 0;
figure(17); clf
imshow2(squeeze(profiles(:,dW+1,:)), [])

profMax = max(profiles, [], 3);
moviefixsc(profiles);

% Select the two regions of the cell on either side of the fiber
[~, intens] = clicky(profiles, profMax);
figure(13); clf
plot(intens(:,1) + intens(:,2));

figure(5); clf
subplot(2,3,1);
imshow2(maxImg, []); hold all
plot(allROI{c}(:,1), allROI{c}(:,2))
hold off
subplot(2,3,2)
plot(allIntens(:,c));
title('Intensity')
subplot(2,3,3)
plot(allQ(:,c))
title('Angle')
thresh1 = 35;
thresh2 = 35;
subplot(2,3,4:6)
for j = 1:nFrames;
    trace = allProfile{c}{j};
    trace = trace - prctile(trace, 10);
    traceS = medfilt1(trace, 10);
    plot(trace); hold all
    plot(traceS);
    ylim([-20 75])
    if ismember(j, goodF);
        l1 = find(traceS > thresh1, 1, 'first');
        l2 = find(traceS > thresh2, 1, 'last');
        if ~(isempty(l1) || isempty(l2));
            x1(j) = l1; 
            x2(j) = l2;
            len(j,c) = x2(j) - x1(j);
            hold all
            plot(x1(j):x2(j), trace(x1(j):x2(j)), 'r-', 'LineWidth', 2);
            hold off
        else
            len(j,c) = NaN;
            x1(j) = NaN; x2(j) = NaN;
        end;
    else
        len(j,c) = NaN;
    end;
    title(num2str(j));
    pause(.1)
end;

kymo = zeros(400,nFrames);
for j = 1:nFrames
    trace = allProfile{c}{j};
    if ~isnan(x2(j));
        ctr(j) = round(mean([x2(j) x1(j)]));
    else
        ctr(j) = round(length(trace)/2);
    end;
end;

ctr2 = ctr;
for j = 1:nFrames;
    trace = allProfile{c}{j};
    halfL = min([length(trace)-ctr2(j) ctr2(j)])-1;
    kymo((201-halfL):(201+halfL),j) = trace(ctr2(j)-halfL:ctr2(j)+halfL);
end;

figure(11); clf
kymo(isnan(kymo(:))) = 0;
hold all
kymo2 = zeros(size(kymo,1), size(kymo,2),3);
kymo2(:,:,2) = kymo/max(kymo(:));
imshow2(kymo2, [])
daspect([1 5 1])
barX0 = 50;
barY0 = 250;
barDX = 6; % 1 h
barDY = 10/pixSize;  % 10 microns
plot([barX0 barX0 barX0+barDX], [barY0-barDY barY0 barY0], 'LineWidth', 2, 'Color', 'w')
hold off
title(['Fiber ' num2str(c) ', tScale = 1 h yScale = 10 \mum'])
saveas(gca, 'Fiber 7 kymo.fig')
saveas(gca, 'Fiber 7 kymo.png')


figure(6); clf
tau = (0:(nFrames-1))/6; % time in h, frame every 10 min
ax = plotyy(tau, intens(:,1)+intens(:,2), tau(52:end), pixSize*len(52:end,c))
xlabel('Time (h)')
ylabel(ax(1), 'Cytoplasmic intensity (A.U.)')
ylabel(ax(2), 'Length (\mum)')
saveas(gca, 'fiber 7 length and cell brightness.fig')
saveas(gca, 'fiber 7 length and cell brightness.png')

figure(7); clf
m = 1;
[ySize, xSize, ~] = size(profiles);
cmin = min(profiles(:));
cmax = max(profiles(:));
for j = 1:6:nFrames;
    colorimg = zeros(ySize, xSize, 3);
    colorimg(:,:,2) = profiles(:,:,j)/255;
    subplot(3,6,m);
    imshow2(colorimg);
    m = m+1;
end;

%%
c = 10;
figure(4); clf
bBox = allbBox(:,c);
subMov = mov(bBox(2):(bBox(2)+bBox(4)),bBox(1):(bBox(1)+bBox(3)),:);
moviefixsc(subMov);
% goodF = [9:25 35:64];
goodF = 1:78;

figure(5); clf
subplot(2,3,1);
imshow2(maxImg, []); hold all
plot(allROI{c}(:,1), allROI{c}(:,2))
hold off
subplot(2,3,2)
plot(allIntens(:,c));
title('Intensity')
subplot(2,3,3)
plot(allQ(:,c))
title('Angle')
thresh1 = 15;
thresh2 = 20;
subplot(2,3,4:6)
for j = 1:nFrames;
    trace = allProfile{c}{j};
    trace = trace - prctile(trace, 10);
    traceS = medfilt1(trace, 10);
    plot(trace); hold all
    plot(traceS);
    ylim([-20 75])
    if ismember(j, goodF);
        l1 = find(traceS > thresh1, 1, 'first');
        l2 = find(traceS > thresh2, 1, 'last');
        if ~(isempty(l1) || isempty(l2));
            x1(j) = l1; 
            x2(j) = l2;
            len(j,c) = x2(j) - x1(j);
            hold all
            plot(x1(j):x2(j), trace(x1(j):x2(j)), 'r-', 'LineWidth', 2);
            hold off
        else
            len(j,c) = NaN;
            x1(j) = NaN; x2(j) = NaN;
        end;
    else
        len(j,c) = NaN;
    end;
    title(num2str(j));
    pause(.1)
end;

kymo = zeros(400,nFrames);
for j = 1:nFrames
    trace = allProfile{c}{j};
    if ~isnan(x2(j));
        ctr(j) = round(mean([x2(j) x1(j)]));
%     else
%         ctr(j) = round(length(trace)/2);
    end;
end;

ctr2 = ctr;
for j = 1:nFrames;
    trace = allProfile{c}{j};
    halfL = min([length(trace)-ctr2(j) ctr2(j)])-1;
    kymo((201-halfL):(201+halfL),j) = trace(ctr2(j)-halfL:ctr2(j)+halfL);
end;

figure(11); clf
kymo(isnan(kymo(:))) = 0;
kymo2 = zeros(size(kymo,1), size(kymo,2), 3);
kymo2(:,:,2) = 1.2*kymo/max(kymo(:)) + 0.1;
hold all
imshow2(kymo2, [])
daspect([1 5 1])
barX0 = 10;
barY0 = 250;
barDX = 6; % 1 h
barDY = 10/pixSize;  % 10 microns
plot([barX0 barX0 barX0+barDX], [barY0-barDY barY0 barY0], 'LineWidth', 2, 'Color', 'w')
hold off
title(['Fiber ' num2str(c) ', tScale = 1 h yScale = 10 \mum'])
xlim([0 65])
ylim([125 275])
saveas(gca, ['Fiber_' num2str(c) '_kymo.fig'])
saveas(gca, ['Fiber_' num2str(c) '_kymo.png'])

%%



figure(6); clf
tau = (0:(nFrames-1))/6; % time in h, frame every 10 min
plot(tau, pixSize*len(:,1:c))
plot( pixSize*len(:,1:c))
xlabel('Time (h)')
ylabel('Length (\mum)')


c = 6;
figure(4); clf
bBox = allbBox(:,c);
subMov = mov(bBox(2):(bBox(2)+bBox(4)),bBox(1):(bBox(1)+bBox(3)),:);
moviefixsc(subMov);
goodF = [1:78];
figure(5); clf
subplot(2,3,1);
imshow2(maxImg, []); hold all
plot(allROI{c}(:,1), allROI{c}(:,2))
hold off
subplot(2,3,2)
plot(allIntens(:,c));
title('Intensity')
subplot(2,3,3)
plot(allQ(:,c))
title('Angle')
thresh1 = 15;
thresh2 = 20;
subplot(2,3,4:6)
for j = 1:nFrames;
    trace = allProfile{c}{j};
    trace = trace - prctile(trace, 10);
    traceS = medfilt1(trace, 10);
    plot(trace); hold all
    plot(traceS);
    ylim([-20 75])
    if ismember(j, goodF);
        l1 = find(traceS > thresh1, 1, 'first');
        l2 = find(traceS > thresh2, 1, 'last');
        if ~(isempty(l1) || isempty(l2));
            x1(j) = l1; 
            x2(j) = l2;
            len(j,c) = x2(j) - x1(j);
            hold all
            plot(x1(j):x2(j), trace(x1(j):x2(j)), 'r-', 'LineWidth', 2);
            hold off
        else
            len(j,c) = NaN;
            x1(j) = NaN; x2(j) = NaN;
        end;
    else
        len(j,c) = NaN;
    end;
    title(num2str(j));
    pause(.1)
end;

kymo = zeros(400,nFrames);
for j = 1:nFrames
    trace = allProfile{c}{j};
    if ~isnan(x2(j));
        ctr = round(mean([x2(j) x1(j)]));
        halfL = min([length(trace)-ctr ctr])-1;
        kymo((201-halfL):(201+halfL),j) = trace(ctr-halfL:ctr+halfL);
    end;
end;
figure(11); clf
kymo = kymo(101:300,63:end);
imshow2(kymo, [])
colormap('hot')
caxis([-5 100])
daspect([1 6 1])



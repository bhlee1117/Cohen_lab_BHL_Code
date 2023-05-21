% AnalyzeTrajecsBeadsinMucin.m
% Accompanies file AnalyzeDataBeadsinMucinByFrame
% Loads data saved from there to calculate diffusion and other interesting
% statistics
% JHH 22 March 2009

clear all;
pack;

load trajecsBeadsMucinBig;
% trajecsBeadsMucin is a MATLAB contains the following variables:
% centers, cell array with position of centers of particles for each file
    %  centers{iFile}(:,1) is the x-coordinates
    %  centers{iFile}(:,2) is the y-coordinates
    %  centers{iFile}(:,3) is the brightnesses
    %  centers{iFile}(:,4) is the square of the radius of gyration
    %  centers{iFile}(:,5) is the frame in which the particle was found
% tracks, the linked trajectory of the particles.  
    % It has the same structure as centers, but with a sixth column
    % containing a unique particle ID for each particle that persists
    % between frames. tracks is sorted by particle ID.
% path, file directory path of raw data files
% fileNames, list of file names
% nFiles, number of files analyzed in centers, tracks
% sizeX,
% sizeY,
% nFrames, array of number of frames analyzed for each file
% secsPerFrame, number of seconds per frame; an array of size nFiles.

% nTotFrames = nFrames*ones(nFiles, 1);
nTotFrames = nFrames;
clear nFrames;

colors = ['r', 'y', 'g', 'c', 'b', 'm', 'k', 'r'];

% For plotting purposes, make a third row in groupNamesByFile with full
% name. Rows 1 and 2 correspond to groupings 1 and 2 found in rows of
% groups
groupNamesByFile = cell(length(groupNames)+1, nFiles);
for iFile = 1:nFiles
    for iGroupDim = 1:length(groupNames)
        groupNamesByFile{iGroupDim,iFile} = groupNames{iGroupDim}{groups(iGroupDim,iFile)};
    end
    groupNamesByFile{end,iFile} = [groupNamesByFile{1,iFile}, ', ', groupNamesByFile{2,iFile}];
end

% Build array of structures trajec, with iFile elements.
% Fields include:
% trajec(iFile).nParticles, number of particles tracked in iFile
% trajec(iFile).centers, 5-column array centers 
% trajec(iFile).nFrames, array of size trajec(iFile).nParticles containing
% number of frames each particle is tracked for

% Stores trajectory data into trajec(iFile).centers, a cell array of length
% nParticles, with each cell element containing the information in tracks

for iFile = 1:nFiles;
    
    trajec(iFile).nParticles = max(tracks{iFile}(:,6));
    trajec(iFile).centers = cell(trajec(iFile).nParticles,1);
    trajec(iFile).nFrames = zeros(trajec(iFile).nParticles,1);

    for iParticle = 1:trajec(iFile).nParticles;
      
        trajec(iFile).centers{iParticle} = tracks{iFile}(tracks{iFile}(:,6) == iParticle,1:5);
        trajec(iFile).nFrames(iParticle) = sum(tracks{iFile}(:,6) == iParticle);
        
    end; %iParticle loop
end;

% Plot trajectories in space
for iFile = 1:nFiles
    figure(iFile)
    subplot(2,2,1)
    for iParticle = 1:trajec(iFile).nParticles;
        plot(trajec(iFile).centers{iParticle}(:,1), trajec(iFile).centers{iParticle}(:,2), colors(mod(iParticle,length(colors))+1));
        hold on;
    end % iParticle loop
    title(plotNames{iFile})
    xlim([0,sizeX]);
    ylim([0,sizeY]);
    axis equal;
    hold off;
end % iFile loop

% Remove mean drift per file
% Calculate mean displacement between each frame and subtract off

for iFile = 1:nFiles
    
    driftX = zeros(nTotFrames(iFile), 1);
    driftY = zeros(nTotFrames(iFile), 1);
    countFrames = zeros(nTotFrames(iFile), 1);
    countFrames(1) = 1;
    trajec(iFile).centersNoDrift = cell(trajec(iFile).nParticles,1);
   
    for iParticle = 1:trajec(iFile).nParticles
        nFramesTmp = trajec(iFile).nFrames(iParticle);
        framesTmp = trajec(iFile).centers{iParticle}(:,5);
        driftX(framesTmp(2:end)) = driftX(framesTmp(2:end)) + diff(trajec(iFile).centers{iParticle}(:,1));
        driftY(framesTmp(2:end)) = driftY(framesTmp(2:end)) + diff(trajec(iFile).centers{iParticle}(:,2));
        countFrames(framesTmp(2:end)) = countFrames(framesTmp(2:end)) + ones(size(framesTmp(2:end)));
    end % iParticle loop
    
    driftX = driftX./countFrames;
    driftY = driftY./countFrames;
    
    driftXCum = cumsum(driftX);
    driftYCum = cumsum(driftY);
    
    figure(iFile)
    subplot(2,2,2)
    for iParticle = 1:trajec(iFile).nParticles
        trajec(iFile).centersNoDrift{iParticle}(:,1) = trajec(iFile).centers{iParticle}(:,1) - driftXCum(1:trajec(iFile).nFrames(iParticle));
        trajec(iFile).centersNoDrift{iParticle}(:,2) = trajec(iFile).centers{iParticle}(:,2) - driftYCum(1:trajec(iFile).nFrames(iParticle));
        plot(trajec(iFile).centersNoDrift{iParticle}(:,1), trajec(iFile).centersNoDrift{iParticle}(:,2), colors(mod(iParticle,length(colors))+1));
        hold on;
    end % iParticle loop

    title(['Corrected trajectories'])
    xlim([0,sizeX])
    ylim([0,sizeY])
    axis equal;
    hold off;
    
    subplot(2,2,3)
    plot(driftX, 'b*')
    hold on;
    plot(driftY, 'g*')
    hold on;
    plot(driftXCum, 'b')
    hold on;
    plot(driftYCum, 'g')
    hold on;
    
% %     plot(trajec(iFile).driftx, 'r*')
% %     hold on;
% %     plot(trajec(iFile).drifty, 'k*')
% %     hold on;
% %     plot(cumsum(trajec(iFile).driftx), 'r')
% %     hold on;
% %     plot(cumsum(trajec(iFile).drifty), 'k')
% %     hold on;
% %     
    title('Mean drift')
    xlabel('Frame number')
    ylabel('Mean displacement (pixels)')
    legend('Mean x drift per frame', 'Mean y drift per frame', 'Cumulative x drift', 'Cumulative y drift')
    hold off;
    
end % iFile loop

clear driftX driftY driftXCum driftYCum countFrames


% Take a look at statistics of trajectories
% Number of particles tracked per file
figure(nFiles + 1);
subplot(2,1,1)
for iFile = 1:nFiles
    plot(groups(1,iFile), trajec(iFile).nParticles, [colors(groups(2,iFile)), '*'])
    hold on;
end
xlim([0,nGroups(1)+1]);
title('Number of particles tracked')
set(gca,'XTick', 0:nGroups(1)+1)
set(gca,'XTickLabel', [' ', groupNames{1}, ' '])
legend(groupNamesByFile{3,:});


%Histogram of number of frames each particle is tracked, per file
figure(nFiles + 1);
subplot(2,1,2)
for iFile = 1:nFiles
    [count, binIndex] = histc(trajec(iFile).nFrames, 1:nTotFrames);
    plot(1:nTotFrames, count, colors(groups(1,iFile)))
    hold on;
end
title('Length of trajectory');
xlabel('Number of frames tracked')
ylabel('Number of particles')
legend(groupNamesByFile{3, :});


%First difference diffusion coefficients
for iFile = 1:nFiles;
    trajec(iFile).dx = zeros([trajec(iFile).nParticles, 1]);
    trajec(iFile).dy = zeros([trajec(iFile).nParticles, 1]);

    for iParticle = 1:trajec(iFile).nParticles;
        deltaX = diff(trajec(iFile).centersNoDrift{iParticle}(:,1));
        deltaY = diff(trajec(iFile).centersNoDrift{iParticle}(:,2));
        deltaT = diff(trajec(iFile).centers{iParticle}(:,5)); %in units of frames
        
        trajec(iFile).dx(iParticle) = mean((deltaX.^2)./(2*deltaT));
        trajec(iFile).dy(iParticle) = mean((deltaY.^2)./(2*deltaT));
    end 
    
    trajec(iFile).dx = trajec(iFile).dx*(micronsPerPixel(iFile)^2)/secsPerFrame(iFile);
    trajec(iFile).dy = trajec(iFile).dy*(micronsPerPixel(iFile)^2)/secsPerFrame(iFile); 
    dxMean(iFile) = sum(trajec(iFile).dx.*trajec(iFile).nFrames(:))/sum(trajec(iFile).nFrames);
    dyMean(iFile) = sum(trajec(iFile).dy.*trajec(iFile).nFrames(:))/sum(trajec(iFile).nFrames);
end

figure(nFiles + 2);
subplot(2,1,1)
for iFile = 1:nFiles
    plot(groups(1, iFile), dxMean(iFile) + dyMean(iFile), [colors(groups(2, iFile)), '*'])
    hold on;
end
xlim([0, nGroups(1)+1]);
set(gca,'XTick', 0:nGroups(1)+1)
set(gca,'XTickLabel', [' ', groupNames{1}, ' '])
legend(groupNamesByFile{2,1:max(groups(2,:))});
ylabel('Mean Diffusion Coefficient (\mum^2/s)');

binSize = .01;
binEdges = 0:binSize:0.1;

for iFile = 1:nFiles
    subplot(2,2,3)
    [count, binIndex] = histc(trajec(iFile).dx, binEdges);
    plot(binEdges+binSize/2, count, colors(iFile));
    hold on;
    xlim([0, binEdges(end)])
    title('Histogram of diffusion coefficients in X')
    xlabel('Mean D_x (\mum^2/s)')
    ylabel('Number of particles')
    
    subplot(2,2,4)
    [count, binIndex] = histc(trajec(iFile).dy, binEdges);
    plot(binEdges+binSize/2, count, colors(iFile));
    hold on;
    xlim([0, binEdges(end)])
    title('Histogram of diffusion coefficients in Y')
    legend(groupNamesByFile{3,:});
    xlabel('Mean D_y (\mum^2/s)')
    ylabel('Number of particles')
    
end

%Plot first difference diff coeff vs. brightnesses and area of particle
figure(nFiles + 3)
for iFile = 1:nFiles
    for iParticle = 1:trajec(iFile).nParticles
        subplot(2,1,1)
        plot(trajec(iFile).centers{iParticle}(:, 3), trajec(iFile).dx(iParticle), [colors(iFile), '*'])
        hold on;
        subplot(2,1,2)
        plot(trajec(iFile).centers{iParticle}(:, 4), trajec(iFile).dx(iParticle), [colors(iFile), '*'])
        hold on;
    end % iParticle loop
end % iFile loop
subplot(2,1,1)
xlabel('Brightness of particle')
ylabel('Diffusion coefficient')
subplot(2,1,2)
xlabel('Radius of particle')
ylabel('Diffusion coefficient')
legend(groupNamesByFile{3,:})


% Extract diffusion coefficients:
for iFile = 1:nFiles;
    trajec(iFile).deltaR2 = zeros(nTotFrames(iFile),1);
    tally = zeros(nTotFrames(iFile),1);
    for iParticle = 1:trajec(iFile).nParticles
        nFramesTmp = trajec(iFile).nFrames(iParticle);
        for lag = 1:(nFramesTmp-1);
            
            tau = 1:lag:nFramesTmp;
            tau1 = tau(1:(end - 1));
            tau2 = tau(2:end);
            deltaX = trajec(iFile).centersNoDrift{iParticle}(tau2,1) - trajec(iFile).centersNoDrift{iParticle}(tau1,1);
            deltaY = trajec(iFile).centersNoDrift{iParticle}(tau2,2) - trajec(iFile).centersNoDrift{iParticle}(tau1,2);
            tally(lag) = tally(lag) + length(tau1);
            trajec(iFile).deltaR2(lag) = trajec(iFile).deltaR2(lag) + deltaX'*deltaX + deltaY'*deltaY;
        end;
    end;
    trajec(iFile).deltaR2 = trajec(iFile).deltaR2./tally;
end;

for iFile = 1:nFiles;
    if groups(1, iFile) == 1
        figure(nFiles+8);
    else
        figure(nFiles+9);
    end
    
    timeAxis = (1:nTotFrames(iFile))*secsPerFrame(iFile);
    plot(log(timeAxis), log(trajec(iFile).deltaR2), [colors(iFile), '*']);
    hold all;
    trajec(iFile).linFits = polyfit(timeAxis, trajec(iFile).deltaR2', 1);
%     plot(timeAxis, timeAxis*trajec(iFile).linFits(1) + trajec(iFile).linFits(2), colors(iFile))
    
%     figure(nFiles+7)
%     plot(trajec(iFile).linFits(1), trajec(iFile).linFits(2), [colors(iFile), '*'])
%     hold on;
end;

figure(nFiles+8)
hold off;
title('<r^2> vs t')
legend(groupNamesByFile{3,groups(1,:) == 1})
xlabel('Time (seconds)')
ylabel('<r^2> (pixels^2)')

figure(nFiles+9)
hold off;
title('<r^2> vs t')
legend(groupNamesByFile{3,groups(1,:) == 2})
xlabel('Time (seconds)')
ylabel('<r^2> (pixels^2)')

figure(nFiles+7)
hold off;
title('Linear Coefficients')
legend(groupNamesByFile{3,:)
xlabel('')
ylabel('<r^2> (pixels^2)')


figure(nFiles+10)
for iFile = 1:nFiles
    plot(trajec(iFile).nFrames, trajec(iFile).dx + trajec(iFile).dy, [colors(iFile), '*'])
    hold on
end
title('Diffusion Coefficient vs. Num Frames Tracked');
xlabel('Number of frames tracked')
ylabel('Diffusion coefficient (\mum^2/s)')












% figure(8);


% for k = 1:nfiles;
%     [n, b] = histc(trajec(k).Dx + trajec(k).Dy);
%     plot(b, n, clist{k});
%     hold all;
% end;
% hold off;
% title('Histogram of diffusion coefficients');

% % %Length of trajectory vs. diffusion coefficient
% % figure(9);
% % for k = 1:length(noPK);
% %     plot(trajec(noPK(k)).nframes, trajec(noPK(k)).Dx + trajec(noPK(k)).Dy, 'b*');
% %     hold on;
% % end;
% % 
% % for k = 1:length(PK0);
% %     plot(trajec(PK0(k)).nframes, trajec(PK0(k)).Dx + trajec(PK0(k)).Dy, 'g*');
% %     hold on;
% % end;
% % 
% % for k = 1:length(PK1);
% %     plot(trajec(PK1(k)).nframes, trajec(PK1(k)).Dx + trajec(PK1(k)).Dy, 'r*');
% %     hold on;
% % end;
% % 
% % for k = 1:length(PK6);
% %     plot(trajec(PK6(k)).nframes, trajec(PK6(k)).Dx + trajec(PK6(k)).Dy, 'k*');
% %     hold off;
% % end;
% % 
% % title('Diffusion Coefficient vs. Length particle tracked');
% % xlabel('Number of frames tracked')
% % ylabel('Diffusion coefficient (\mum^2/s)')
% % 
% % 
% % % Extract the initial intensities of the particles that are tracked in the
% % % first frame of the movies
% % boxsize = 2;
% % for k = 1:nfiles;
% %     figure(4 + k)
% %     m = mean(mean(firstimg(:,:,k)));
% %     s = std(reshape(firstimg(:,:,k), xsize*ysize,1));
% %     tmp = imfilter(firstimg(:,:,k), fspecial('Gaussian', 5,2), 'replicate');
% %     imshow(tmp, [m - 2*s, m + 2*s]);
% %     line([30,30+38], [20,20], 'Color', 'k', 'lineWidth', 4)
% %     text(13, 40, '10 \mum', 'fontSize', 20, 'Color', 'k')
% % 
% %     hold all;
% %     trajec(k).initI = zeros(trajec(k).nparticles,1);
% %     for j = 1:trajec(k).nparticles;
% %         if trajec(k).centersNoDrift{j}(1,5) == 1;
% %             left = floor(trajec(k).centersNoDrift{j}(1,1) - boxsize);
% %             right = ceil(trajec(k).centersNoDrift{j}(1,1) + boxsize);
% %             top = ceil(trajec(k).centersNoDrift{j}(1,2) + boxsize);
% %             bottom = floor(trajec(k).centersNoDrift{j}(1,2) - boxsize);
% %             spot = squeeze(firstimg(bottom:top, left:right, k));
% %             plot(trajec(k).centersNoDrift{j}(1,1), trajec(k).centersNoDrift{j}(1,2), 'o');
% %             trajec(k).initI(j) = intens(k,1)*mean(mean(spot))/laserpwr(k);
% %         end;
% %     end;
% %     hold off;
% %     pause
% % end;
% % 
% % 
% % 
% % for iFile = 1:nFiles;
% %     
% %     trajec(iFile).dxCov = zeros(nTotFrames(iFile),1);
% %     trajec(iFile).dyCov = zeros(nTotFrames(iFile),1);
% %     trajec(iFile).dxdyCov = zeros(nTotFrames(iFile),1);
% %     
% %     for iParticle = 1:trajec(iFile).nParticles;
% %         
% %         dx = diff(trajec(iFile).centersNoDrift{iParticle}(:,1));
% %         dy = diff(trajec(iFile).centersNoDrift{iParticle}(:,2));
% %         
% %         trajecLength = trajec(iFile).nFrames(iParticle) - 1;
% %         iCov = (trajecLength):(2*trajecLength - 1); %only need second half, because array is symmetric
% %         iTrajec = 1:trajecLength;
% %         
% %         trajecCov = xcov(dx, 'unbiased');
% %         trajec(iFile).dxCov(iTrajec) = trajec(iFile).dxCov(iTrajec) + trajec(iFile).nFrames(iParticle)*trajecCov(iCov); % weight trajectory covariance by length of trajectory
% %         
% %         trajecCov = xcov(dy, 'unbiased');
% %         trajec(iFile).dyCov(iTrajec) = trajec(iFile).dyCov(iTrajec) + trajec(iFile).nFrames(iParticle)*trajecCov(iCov); % weight trajectory covariance by length of trajectory
% %         
% %         trajecCov = xcov(dx, dy, 'unbiased');
% %         trajec(iFile).dxdyCov(iTrajec) = trajec(iFile).dxdyCov(iTrajec) + trajec(iFile).nFrames(iParticle)*trajecCov(iCov); % weight trajectory covariance by length of trajectory
% %         
% %     end; % iParticle loop
% %     
% %     trajec(iFile).dxCov = trajec(iFile).dxCov/sum(trajec(iFile).nFrames); % normalize by sum of all trajectory lengths
% %     trajec(iFile).dyCov = trajec(iFile).dyCov/sum(trajec(iFile).nFrames); 
% %     trajec(iFile).dxdyCov = trajec(iFile).dxdyCov/sum(trajec(iFile).nFrames); 
% % 
% % end; % iFile loop
% % 
% % figure(nFiles+4);
% % for iFile = 1:7;
% %     subplot(2,1,1)
% %     plot(trajec(iFile).dxCov, colors(iFile));
% %     hold on;
% %     plot(trajec(iFile).dyCov, colors(iFile));
% %     
% %     subplot(2,1,2)
% %     plot(trajec(iFile).dxdyCov, colors(iFile));
% %     hold on;
% % end;
% % % xlim([0,20]);
% % hold off;
% % 
% % 

% 
% figure(6);
% for k = 1:nfiles;
%     tmp = trajec(k).initI;
%     tmp(tmp == 0) = [];
%     [n, b] = hist(tmp);
%     plot(b, n/sum(n), clist{k});
%     hold all;
% end;
% hold off;
% title('Histogram of initial intensities');    
% 
% % Plot the initial intensity of particles that appeared in the first frame
% % as a function of their x diffusion coefficient.
% figure(7);
% for k = 1:nfiles;
%     plot(trajec(k).Dx(trajec(k).initI ~= 0), trajec(k).initI(trajec(k).initI ~=0), [clist{k} 'x']);
%     hold all;
% end;
% hold off;
% 
% for k = 1:nfiles;
%     trajec(k).dxcov = zeros(ntotframes(k),1);
%     trajec(k).dycov = zeros(ntotframes(k),1);
%     trajec(k).dxdycov = zeros(ntotframes(k),1);
%     for j = 1:trajec(k).nparticles;
%         dx = diff(trajec(k).centers{j}(:,1));
%         dy = diff(trajec(k).centers{j}(:,2));
%         L = trajec(k).nframes(j) - 1;
%         tmp = xcov(dx, 'unbiased');
%         tmpindx = (L):(2*L - 1);
%         trajindx = 1:L;
%         trajec(k).dxcov(trajindx) = trajec(k).dxcov(trajindx) + trajec(k).nframes(j)*tmp(tmpindx);
%         tmp = xcov(dy, 'unbiased');
%         trajec(k).dycov(trajindx) = trajec(k).dycov(trajindx) + trajec(k).nframes(j)*tmp(tmpindx);
%         tmp = xcov(dx, dy, 'unbiased');
%         trajec(k).dxdycov(trajindx) = trajec(k).dxdycov(trajindx) + trajec(k).nframes(j)*tmp(tmpindx);
%     end;
%     trajec(k).dxcov = trajec(k).dxcov/sum(trajec(k).nframes);
%     trajec(k).dycov = trajec(k).dycov/sum(trajec(k).nframes);
%     trajec(k).dxdycov = trajec(k).dxdycov/sum(trajec(k).nframes);
% end;
% 
% figure(8);
% for k = 1:nfiles;
%     plot(trajec(k).dxcov, clist{k});
%     plot(trajec(k).dycov, [clist{k} '*']);
%     hold all;
% end;
% hold off;
% 
% figure(9);
% for k = 1:nfiles;
%     plot(trajec(k).dxdycov, clist{k});
%     hold all;
% end;
% hold off;
% 
% for k = 1:nfiles;
%     k
%     trajec(k).deltar2 = zeros(ntotframes(k),1);
%     tally = zeros(ntotframes(k),1);
%     for j = 1:trajec(k).nparticles
%         n = trajec(k).nframes(j);
%         for lag = 1:(n-1);
%             tau = 1:lag:n;
%             tau1 = tau(1:(end - 1));
%             tau2 = tau(2:end);
%             deltax = trajec(k).centers{j}(tau2,1) - trajec(k).centers{j}(tau1,1);
%             deltay = trajec(k).centers{j}(tau2,2) - trajec(k).centers{j}(tau1,2);
%             tally(lag) = tally(lag) + length(tau1);
%             trajec(k).deltar2(lag) = trajec(k).deltar2(lag) + deltax'*deltax + deltay'*deltay;
%         end;
%     end;
%     trajec(k).deltar2 = trajec(k).deltar2./tally;
% end;
% 
% figure(10);
% for k = 1:nfiles;
%     plot(trajec(k).deltar2, clist{k});
%     hold all;
% end;
% hold off;
% title('<r^2> vs t')

% % % This algorithm only works assuming that tracks is sorted by particle ID.
% % % for iFile = 1:nFiles;
% % %     trajec(iFile).nParticles = max(tracks{iFile}(:,6));
% % %     trajec(iFile).centers = cell(trajec(iFile).nParticles,1);
% % %     iEnd = 1;
% % %     iStart = 1;
% % %     for iParticle = 1:trajec(iFile).nParticles;
% % %         while tracks{iFile}(iEnd,6) == iParticle;
% % %             iEnd = iEnd + 1;
% % %             if iEnd > size(tracks{iFile},1);
% % %                 break;
% % %             end;
% % %         end;
% % %         trajec(iFile).centers{iParticle} = tracks{iFile}(iStart:(iEnd - 1),:);
% % %         trajec(iFile).nFrames(iParticle) = length(trajec().centers{j}(:,1));
% % %         iStart = iEnd;
% % %     end;
% % % end;
% % % 

% % nParticlesBar = zeros(nGroups(1), sum(groups(1,:) == mode(groups(1,:))))
% % count = zeros(nGroups(1), 1)
% % barGroup = zeros(sum(groups(1,:) == mode(groups(1,:))), 1)
% % 
% % for iFile = 1:nFiles
% %     tmpGroup = groups(1,iFile)
% %     nParticlesBar(tmpGroup, count(tmpGroup)+1) = trajec(iFile).nParticles
% %     barGroup() =  groups(2,iFile);
% %     count(tmpGroup) = count(tmpGroup)+1
% % end
% figure(5);
% plot(1:length(noPK), meandx(noPK) + meandy(noPK), 'b*',...
%     1:length(PK0), meandx(PK0) + meandy(PK0), 'g*',...
%     1:length(PK1), meandx(PK1) + meandy(PK1), 'r*',...
%     1:length(PK6), meandx(PK6) + meandy(PK6), 'k*');
% legend('no PK', 'PK, t = 0h', 'PK, t = 1h', 'PK, t = 5.5h')
% title('Diffusion coefficients')
% figure;
% % Calculate the drift per frame.
% for k = 1:nFiles;
%       trajec(k).nperframe = zeros(nTotFrames(k),1);
%       trajec(k).driftx = zeros(nTotFrames(k)-1,1);
%       trajec(k).drifty = zeros(nTotFrames(k)-1,1);
%       for j = 1:nTotFrames(k)-1;
%           trajec(k).nperframe(j) = sum(tracks{k}(:,5) == j);
%           trajec(k).driftx(j) = sum( diff(tracks{k}(:,1)) .* (tracks{k}(1:end-1,5) == j)...
%                                 .* (tracks{k}(1:end-1,6) == tracks{k}(2:end,6)))/...
%                                 sum((tracks{k}(1:end-1,5) == j)...
%                                 .* (tracks{k}(1:end-1,6) == tracks{k}(2:end,6)));
%           trajec(k).drifty(j) = sum( diff(tracks{k}(:,2)) .* (tracks{k}(1:end-1,5) == j)...
%                                 .* (tracks{k}(1:end-1,6) == tracks{k}(2:end,6)))/...
%                                 sum((tracks{k}(1:end-1,5) == j)...
%                                 .* (tracks{k}(1:end-1,6) == tracks{k}(2:end,6)));                            
%       end;
%       
%       subplot(2,1,1)
% %       plot(trajec(k).driftx, colors(k))
% %       hold on;
%       plot(cumsum(trajec(k).driftx), colors(k))
%       hold on; 
%       
%       subplot(2,1,2)
% %       plot(trajec(k).drifty, colors(k))
% %       hold on;
%       
%       plot(cumsum(trajec(k).drifty), colors(k))
%       hold on;
%       
% end;
% legend(groupNamesByFile{3,:})
% 
% 
% 
% 

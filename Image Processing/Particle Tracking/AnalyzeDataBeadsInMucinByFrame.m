% Using the particle tracking software to analyze diffusing bead data
% Reading in one frame at a time to reduce memory usage
% JHH 22 March 2009

% Variables to change before running the script:
% path, file directory path of raw data files
% bgFiles, list of background file names
% fileNames, list of file names
% plotNames, list of file names, for plotting
% bgPairs, array containg the index of background file for each data file
% groups, 2D array whose rows represent different groupings (can only handle
% two ways to group)
% groupNames, cell arrya for titles for the groups
% nFiles, number of files analyzed in centers, tracks
% sizeX, width of image in pixels (640 for Genwac, 512 for Andor)
% sizeY, height of image in pixels (480 for Genwac, 512 for Andor)
% nFrames, array of size nFiles, number of frames analyzed for each file
% secsPerFrame, array of size nFiles, number of seconds per frame for each file 

clear all; 
close all;

path = 'C:\Users\Jennifer\Desktop\temp_matlab\With_no_protease\';

bgFiles = {'1p5MpH7_200nm_bg', '1p5MpH7_200nm_bg_2', '1p5MpH5_200nm_bg'}
fileNames = {'1p5MpH7_200nm_30Hz_60x', ...
    '1p5MpH7_200nm_3Hz_60x', ...
    '1p5MpH7_200nm_0p3Hz_60x', ...
    '1p5MpH7_200nm_0p3Hz_60x_2', ...
    '1p5MpH5_200nm_30Hz_60x', ...
    '1p5MpH5_200nm_3Hz_60x', ...
    '1p5MpH5_200nm_0p3Hz_60x', ...
    '1p5MpH5_200nm_0p3Hz_60x_2'}

plotNames = {'1.5% pH7 Mucin, 200nm beads, 30Hz, 60x', ...
    '1.5% pH7 Mucin, 200nm beads, 3Hz, 60x', ...
    '1.5% pH7 Mucin, 200nm beads, 0.3Hz, 60x', ...
    '1.5% pH7 Mucin, 200nm beads, 0.3Hz, 60x, 2', ...
    '1.5% pH5 Mucin, 200nm beads, 30Hz, 60x', ...
    '1.5% pH5 Mucin, 200nm beads, 3Hz, 60x', ...
    '1.5% pH5 Mucin, 200nm beads, 0.3Hz, 60x', ...
    '1.5% pH5 Mucin, 200nm beads, 0.3Hz, 60x, 2'}

bgPairs = [1, 1, 1, 2, 3, 3, 3, 3]; %index of background file for each data file

groups = [1, 1, 1, 1, 2, 2, 2, 2; %grouped by pH
          1, 2, 3, 3, 1, 2, 3, 3]; %grouped by acquisition speed
nGroups = max(groups');
groupNames{1} = {'pH7', 'pH5'};
groupNames{2} = {'30 Hz', '3 Hz', '0.3 Hz'};

deltaT = 1/30; %default frame rate of Genwac
nFiles = size(fileNames,2);
nFilesBG = size(bgFiles,2);

sizeX = 640;
sizeY = 480;

nFrames = 500*ones(nFiles,1);
nFramesBG = 10*ones(nFilesBG,1);
baseline = 0;
colors = ['r', 'y', 'g', 'c', 'b', 'm', 'k', 'r'];

tracks = cell(nFiles, 1);
% intens = zeros(nFiles, nFrames); % mean intensity of each frame for each file
deltaI = [2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5]; % number of standard deviations for thresholding
% histImg = zeros(sizeY, sizeX, nFiles); % max intensity bin of each pixel
% firstImg = zeros(sizeY, sizeX, nFiles);
centers = cell(nFiles, 1);

secsPerFrame = [1/30, 1/3, 1/0.3, 1/0.3, 1/30, 1/3, 1/0.3, 1/0.3];
micronsPerPixel = 6/76*ones(nFiles, 1); %if 60x objective is use with the genwac.

smoothRad = 5;
pkSize = 13;
winSize = 19; %must be odd for cntrd to work without warnings
maxDisp = 10;
intensMean = zeros(nFiles, max(nFrames));

for iFile = 1:nFiles;
% for iFile = 6:8;
    disp('Loading data from')
    fileNames{iFile}
    bgDat = loadTIFFStackQuiet([path bgFiles{bgPairs(iFile)}], sizeY, sizeX, nFramesBG(bgPairs(iFile)));
    bgMean = mean(bgDat, 3);
    
    for iFrame = 1:nFrames(iFile)
        pack
        iFrame
        tmpFrame = imread([path fileNames{iFile}], 'tif', iFrame);
        tmpFrameBGSub = double(tmpFrame)-bgMean; % subtract off background from separate file
        tmpFrameSmooth = imfilter(tmpFrameBGSub, fspecial('Gaussian', smoothRad, smoothRad/2.5), 'replicate'); % smooth images with gaussian filter

        % smooth images with gaussian filter
        
        frameIntensMean = mean(tmpFrameSmooth(:)); % mean intensity of frame after gaussian smoothing
%         tmpFrameNorm = tmpFrameSmooth; %Do not normalize by /frameIntensMean; % normalize data by mean intensity of the frame
        frameIntensSig = std(tmpFrameSmooth(:)); % std dev of intensity of frame after smoothing
        thresh = frameIntensMean + deltaI(iFile)*frameIntensSig; % set threshold to be some number of std dev above mean intensity (default is 2.5)

        intensMean(iFile, iFrame) = frameIntensMean;
%         disp('Finding Centers')
        %  the tmpCenters variable has 5 columns, with the format:
        %  tmpCenters(:,1) is the x-coordinates
        %  tmpCenters(:,2) is the y-coordinates
        %  tmpCenters(:,3) is the brightnesses
        %  tmpCenters(:,4) is the square of the radius of gyration
        %  tmpCenters(:,5) is the frame in which the particle was found

        %    clear M

        mask = tmpFrameSmooth > thresh; % binary image, 1 = above threshold, 0 = below threshold
        tmpFrameT = double(bwmorph(mask, 'clean')).*tmpFrameSmooth; % remove isolated pixels with bwmorph, then only retain pixel intensities above threshold
        tmpPeaks = pkfnd(tmpFrameT, thresh, pkSize); % tmppklst contains x and y coordinates of peaks
        tmpCenters = cntrd(tmpFrameT, tmpPeaks, winSize); % find centroids to subpixel accuracy, tmpcenters contains x, y, brightness, square of radius
        tmpPeaks = [tmpPeaks iFrame*ones(size(tmpPeaks,1),1)]; % add third column to tmplklst containing frame number of peak
        tmpCenters = [tmpCenters iFrame*ones(size(tmpCenters,1),1)]; %add fifth column to tmpcenters containing frame number of center
        centers{iFile} = [centers{iFile}; tmpCenters]; %concatenate centers from newest frame onto centers{iFiles}

        %Plot if first frame
        if iFrame == 1
            
            figure(iFile)
            subplot(2,2,1)
            imshow(tmpFrame, []);
            title(['Raw data: ' fileNames{iFile}])

            subplot(2,2,2)
            imshow(bgMean, []);
            title('Background correction')

            subplot(2,2,3)
            imshow(tmpFrameBGSub(:,:,1), []);
            title('After background correction')

            subplot(2,2,4)
            imshow(tmpFrameSmooth(:,:,1), []);
            title('After smoothing')
            hold on
            plot(tmpPeaks(:,1), tmpPeaks(:,2), 'r*')
            hold on;
            plot(centers{iFile}(:,1), centers{iFile}(:,2), 'bo')
        end %if iFrame == 1
    end %iFrame

    disp('Linking particles')

    param.mem = 0; %number of frames particle can be lost
    param.dim = 2; %dimensionality ot position data, default is 2
    param.good = 5; %minimum number of frames particle is tracked to be considered "good"
    param.quiet = 0; %don't suppress text

    % tracks is the linked trajectory of the particles.  It has the same
    % structure as centers, but with a sixth column containing a unique
    % particle ID for each particle that persists between frames.
    % tracks is sorted by particle ID, trackSort is sorted by frames.
    
    tracks{iFile} = track(centers{iFile}, maxDisp, param);
    
end %iFile loop

save trajecsBeadsMucinBig, centers, tracks, path, fileNames, plotNames, groups, nGroups, groupNames, nFiles, sizeX, sizeY, nFrames, secsPerFrame, micronsPerPixel;

%Make movies of tracking
for iFile = 1:nFiles
    
    [junk, ind] = sort(tracks{iFile}(:,5), 1, 'ascend');
    tracksSort = tracks{iFile}(ind,:);
    clear junk ind;

    figure(nFiles+iFile);
    clear M;
    frameCount = 1;

    for iFrame = 1:nFrames(iFile)
        tmpFrame = imread([path fileNames{iFile}], 'tif', iFrame);
        imshow(tmpFrame, []);
        hold all
        while tracksSort(frameCount,5) == iFrame
            plot(tracksSort(frameCount,1), tracksSort(frameCount,2), [colors(mod(tracksSort(frameCount,6),length(colors))+1) 'x']);
            frameCount = frameCount + 1;
            if frameCount > length(tracksSort);
                break
            end
        end
        M(iFrame) = getframe;
        hold off
    end %iFrame movie loop

    disp('Writing AVI');
    movie2avi(M, [fileNames{iFile}(1:end) '.avi'], 'FPS', 20, 'compression', 'None');
end %iFile loop



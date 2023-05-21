% Using the particle tracking software to analyze diffusing bead data
% Reading in one frame at a time to reduce memory usage
% JHH 07 April 2009

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
% sizeX, array for width of image in pixels (640 for Genwac, 512 for Andor)
% sizeY, array for height of image in pixels (480 for Genwac, 512 for Andor)
% nFrames, array of size nFiles, number of frames analyzed for each file
% secsPerFrame, array of size nFiles, number of seconds per frame for each file 

clear all; 
close all;

path = 'C:\Users\Jennifer\Desktop\EPB\';

% bgFiles = {'1p5MpH7_200nm_bg', '1p5MpH7_200nm_bg_2', '1p5MpH5_200nm_bg'}
fileNames = {'movie_200nm_1to50_1mgPmlTrypsin_in_DGel_60x_30Hz', ...
    'movie_200nm_1to50_1mgPmlTrypsin_in_DGel_60x_3Hz', ...
    'movie_100nm_1mgPmlTrypsin_in_DGel_60x_3Hz_2halves', ...
    '36nm_beads_low_conc_1mgpermilTryps_in_digested_gel_1s_frames.tif'}

plotNames = {'5% gelatin, 200nm beads, 1 mg/mL trypsin in digested gel', ...
    '5% gelatin, 200nm beads, 1 mg/mL trypsin in digested gel', ...
    '5% gelatin, 100nm beads, 1 mg/mL trypsin in digested gel', ...
    '5% gelatin, 36nm beads, 1 mg/mL trypsin in digested gel'}

% bgPairs = [1, 1, 1, 2, 3, 3, 3, 3]; %index of background file for each data file

groups = [1, 1, 2, 3; %grouped by bead size
          1, 2, 2, 3]; %grouped by acquisition speed
      
nGroups = max(groups');
groupNames{1} = {'200nm', '100nm', '36nm'};
groupNames{2} = {'30 Hz', '3 Hz', '1 Hz'};

deltaT = 1/30; %default frame rate of Genwac
nFiles = size(fileNames,2);
% nFilesBG = size(bgFiles,2);

sizeX = [640, 640, 640, 512];
sizeY = [480, 480, 480, 512];

nFrames = [50, 500, 20, 200];
% nFrames = [500, 500, 400, 200];


% nFramesBG = 10*ones(nFilesBG,1);
baseline = 0;
colors = ['r', 'y', 'g', 'c', 'b', 'm', 'k', 'r'];

tracks = cell(nFiles, 1);
% intens = zeros(nFiles, nFrames); % mean intensity of each frame for each file
deltaI = [2.0, 2.5, 1.5, 2.5, 2.5, 2.5, 2.5, 2.5]; % number of standard deviations for thresholding
% histImg = zeros(sizeY, sizeX, nFiles); % max intensity bin of each pixel
% firstImg = zeros(sizeY, sizeX, nFiles);
centers = cell(nFiles, 1);

secsPerFrame = [1/30, 1/3, 1/3, 1];
micronsPerPixel = [6/76, 6/76, 6/76, 512*16/60]; %if 60x objective is use with the genwac.

smoothRad = [5, 5, 5, 5];
pkSize = [9, 9, 9, 9];
winSize = [15, 15, 15, 15]; %must be odd for cntrd to work without warnings
maxDisp = [15, 15, 15, 5];
intensMean = zeros(nFiles, max(nFrames));

% for iFile = 1:nFiles;
for iFile = 4:4;
    disp('Loading data from')
    fileNames{iFile}
%     bgDat = loadTIFFStackQuiet([path bgFiles{bgPairs(iFile)}], sizeY, sizeX, nFramesBG(bgPairs(iFile)));
%     bgMean = mean(bgDat, 3);
%     
    for iFrame = 1:nFrames(iFile)
        pack
        iFrame
        tmpFrame = imread([path fileNames{iFile}], 'tif', iFrame);
        tmpFrameBGSub = double(tmpFrame);%-bgMean; % subtract off background from separate file
        tmpFrameSmooth = imfilter(tmpFrameBGSub, fspecial('Gaussian', smoothRad(iFile), smoothRad(iFile)/2.5), 'replicate'); % smooth images with gaussian filter

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
        tmpPeaks = pkfnd(tmpFrameT, thresh, pkSize(iFile)); % tmppklst contains x and y coordinates of peaks
        tmpCenters = cntrd(tmpFrameT, tmpPeaks, winSize(iFile)); % find centroids to subpixel accuracy, tmpcenters contains x, y, brightness, square of radius
        tmpPeaks = [tmpPeaks iFrame*ones(size(tmpPeaks,1),1)]; % add third column to tmplklst containing frame number of peak
        tmpCenters = [tmpCenters iFrame*ones(size(tmpCenters,1),1)]; %add fifth column to tmpcenters containing frame number of center
        centers{iFile} = [centers{iFile}; tmpCenters]; %concatenate centers from newest frame onto centers{iFiles}

        %Plot if first frame
        if iFrame == 1
            
            figure(iFile)
            subplot(2,2,1)
            imshow(tmpFrame, []);
            title(['Raw data: ' plotNames{iFile}])

%             subplot(2,2,2)
%             imshow(bgMean, []);
%             title('Background correction')
% 
%             subplot(2,2,3)
%             imshow(tmpFrameBGSub(:,:,1), []);
%             title('After background correction')

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
    param.dim = 2; %dimensionality of position data, default is 2
    param.good = 10; %minimum number of frames particle is tracked to be considered "good"
    param.quiet = 0; %don't suppress text

    % tracks is the linked trajectory of the particles.  It has the same
    % structure as centers, but with a sixth column containing a unique
    % particle ID for each particle that persists between frames.
    % tracks is sorted by particle ID, trackSort is sorted by frames.
    
    tracks{iFile} = track(centers{iFile}, maxDisp(iFile), param);
    
end %iFile loop

save trajecsBeadsGel, centers, tracks, path, fileNames, plotNames, groups, nGroups, groupNames, nFiles, sizeX, sizeY, nFrames, secsPerFrame, micronsPerPixel;

%Make movies of tracking
for iFile = 4:4
    
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



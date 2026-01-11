  aviFile='/Volumes/cohen_lab/Lab/Papers/2025 Voltron Optopatch prism dendrites in vivo/Submission 1 Science/MovieSs/MovieS1_CA3EC_axonProjection.avi';
  mp4File='/Volumes/cohen_lab/Lab/Papers/2025 Voltron Optopatch prism dendrites in vivo/Submission 1 Science/MovieSs/MovieS1_CA3EC_axonProjection.mp4';
    if ~isfile(aviFile)
        error('Input AVI file not found.');
    end

    % Read input video
    vIn = VideoReader(aviFile);

    % Create MP4 writer (H.264)
    vOut = VideoWriter(mp4File, 'MPEG-4');
    vOut.FrameRate = vIn.FrameRate;
    vOut.Quality   = 95;   % 0–100, higher = better quality

    open(vOut);

    fprintf('Converting %s → %s\n', aviFile, mp4File);

    while hasFrame(vIn)
        frame = readFrame(vIn);
        writeVideo(vOut, frame);
    end

    close(vOut);

    fprintf('Done.\n');.


  %%

    
        maxMB = 50;
    inputMP4='/Volumes/cohen_lab/Lab/Papers/2025 Voltron Optopatch prism dendrites in vivo/Submission 1 Science/MovieSs/MovieS3_Raw_voltage_video.mp4';
    outputMP4='/Volumes/cohen_lab/Lab/Papers/2025 Voltron Optopatch prism dendrites in vivo/Submission 1 Science/MovieSs/SizeReduced/MovieS3_Raw_voltage_video_compressed.mp4';

    if ~isfile(inputMP4)
        error('Input file not found.');
    end

    % Target size (bytes)
    maxBytes = maxMB * 1024^2;

    % Read input
    vIn = VideoReader(inputMP4);
    duration = vIn.Duration;   % seconds

    % Initial bitrate estimate (bits/s)
    targetBitrate = (maxBytes * 8) / duration;

    fprintf('Target bitrate: %.2f Mbps\n', targetBitrate / 1e6);

    % Initial encoding parameters
    qualityList = [75 60 55 50 45 40 35 28 20];   % fallback qualities
    frameRateOut = min(vIn.FrameRate, 30);  % cap FPS for compression

    for q = qualityList
        fprintf('Trying Quality = %d\n', q);

        vOut = VideoWriter(outputMP4, 'MPEG-4');
        vOut.FrameRate = frameRateOut;
        vOut.Quality   = q;

        open(vOut);

        vIn.CurrentTime = 0;
        while hasFrame(vIn)
            frame = readFrame(vIn);
            writeVideo(vOut, frame);
        end

        close(vOut);

        fileInfo = dir(outputMP4);
        outMB = fileInfo.bytes / 1024^2;

        fprintf('Output size: %.2f MB\n', outMB);

        if fileInfo.bytes <= maxBytes
            fprintf('Success: file is under %d MB\n', maxMB);
            return;
        end
    end

    warning('Could not reach target size (%d MB). Final size = %.2f MB.', ...
            maxMB, outMB);
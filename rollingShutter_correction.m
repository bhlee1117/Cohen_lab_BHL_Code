function spikemovCorr=rollingShutter_correction(spikeMov,sampleRate,cameraType)

spikeMov = double(spikeMov);
[ysize, xsize, ~] = size(spikeMov);
fprintf(['Sampling Rate: %.20f\n'],sampleRate);
    %sampleRate = 1/exposuretime; %Sampling rate or 1/Device_Data{1, 4}.exposuretime of the output file
    xAxImg = 1/sampleRate:1/sampleRate:size(spikeMov,3)/sampleRate; %TimeAxis
    center = floor(ysize/2);
spikemovCorr = zeros(size(spikeMov));

switch cameraType
    case 'fusion' % Correction for rolling shutter (Byung Hun), Orca Fusion

    for n = 1:ysize
        pixelMat = spikeMov(n,:,:);
        pixelMat = reshape(pixelMat, size(pixelMat,2), size(pixelMat,3))';
        xq = xAxImg - 4.88e-6 * (n-1);
        pixelInterp = interp1(xAxImg, pixelMat, xq, 'cubic');
        pixelInterp = reshape(pixelInterp', 1, size(spikeMov, 2), size(spikeMov, 3));
        spikemovCorr(n,:,:) = pixelInterp;
    end

    case 'flash' % Correction for rolling shutter (Urs)
    
    for n = 1:center
        
        pixelMat = spikeMov([center + n-1, center - (n-1)],:,:);
        pixelMat = reshape(pixelMat, 2*size(pixelMat,2), size(pixelMat,3))';
        xq = xAxImg - 9.74436e-6 * (n-1);
        pixelInterp = interp1(xAxImg, pixelMat, xq, 'cubic');
        pixelInterp = reshape(pixelInterp', 2, size(spikeMov, 2), size(spikeMov, 3));
        spikemovCorr([center + n-1, center - (n-1)],:,:) = pixelInterp;
    end
end
end
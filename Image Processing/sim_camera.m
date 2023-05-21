function [mov, trueImg, trueIntens] = sim_camera(cellCoords, cellSize, cellOffset, cellAmp, spikeExponent, bkg, scale, addNoise, ySize, xSize, nFrames)
% function [mov, trueImg, trueIntens] = sim_camera(cellCoords, cellSize, cellOffset, cellAmp, spikeExponent, bkg, scale, addNoise, ySize, xSize, nFrames)
% Simulate camera data from spiking cells, including subthreshold dynamics.
%
% mov: simulated movie output
% trueImg: matrix containing the image of each cell, dimensions: [ySize
% xSize, nCells]
% trueIntens: intensity trace for each cell, dimensions: [nCells, nFrames]
%
% cellCoords: nCells x 2 vector giving [row column] integer coordinates of
% cell centers.  Must lie inside [1 ySize], [1 xSize].
% cellSize: nCells x 1 vector giving radius, in pixels, of each cell.  If a
% scalar, then all cells have the same radius.
% cellOffset: Baseline brightness of each cell.  If a scalar, then all
% cells have the same baseline brightness
% cellAmp: Maximum spike amplitude.  If a scalar, then all cells have the same max spike amplitude.
% spikeExponent: measure of how spiky the time-trace for each cell is.  If a
% scalar, then all cells have the same spike exponent.
% bkg: background counts, where spike intensities range from 0 to 1.
% scale: scale factor for overall signal, used to determine poisson counts
% for noise
% addNoise: Boolean. 0 --> no noise, 1 --> add Poisson noise.
% ySize, xSize, nFrames: dimensions of out.
% AEC 10/16/2014

nCells = size(cellCoords, 1);
mov = zeros(ySize, xSize, nFrames);
trueImg = zeros(ySize, xSize, nCells);
trueIntens = zeros(nCells, nFrames);

if length(cellSize) == 1;   % If all the cells are the same size, apply to all of them
    cellSize = cellSize*ones(nCells, 1);
end;
if length(spikeExponent) == 1;  % If all spike exponents are the same, apply to all of them
    spikeExponent = spikeExponent*ones(nCells, 1);
end;
if length(cellOffset) == 1;  % If all cells are the same brightness, apply to all of them
    cellOffset = cellOffset*ones(nCells, 1);
end;
if length(cellAmp) == 1;  % If all cells are the same amplitude spikes, apply to all of them
    cellAmp = cellAmp*ones(nCells, 1);
end;


for j = 1:nCells
    timeTrace = cellOffset(j) + cellAmp(j)*rand(1,1,nFrames).^spikeExponent(j);
    
    % make a static image of the cell
    cellImg = zeros(ySize, xSize);
    
    cellImg(cellCoords(j,1), cellCoords(j,2)) = 1;
    cellImg = imdilate(cellImg, strel('ball', cellSize(j), 1, 0))-1;  % This makes a ball-shaped cell, with max 1 and zero background
    
    % modulate the image in time and add to the movie
    mov = mov + repmat(cellImg, [1, 1, nFrames]).*repmat(timeTrace, [ySize, xSize, 1]);
    trueImg(:,:,j) = cellImg;
    trueIntens(j,:) = squeeze(timeTrace)';
end;
mov = mov + bkg;
mov = mov*scale;

if addNoise;
    mov = poissrnd(mov);  % .^(addNoise/2)??
end;
    
    

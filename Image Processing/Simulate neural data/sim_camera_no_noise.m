function [mov, trueImg, trueIntens] = sim_camera_no_noise(cellCoords, cellSize, cellBrightness, spikeExponent, ySize, xSize, nFrames)
% function [mov, trueImg, trueIntens] = sim_camera_no_noise(cellCoords, cellSize, cellBrightness, spikeExponent, ySize, xSize, nFrames)
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
% cellBrightness: Baseline brightness of each cell.  Fluctuations in
% intensity range from 0 to 1, so baseline brightness should be of same
% order.  If a scalar, then all cells have the same baseline brightness.
% spikeExponent: measure of how spiky the time-trace for each cell is.  If a
% scalar, then all cells have the same spike exponent.
% ySize, xSize, nFrames: dimensions of out.
% AEC 10/16/2014

nCells = size(cellCoords, 1);
mov = zeros(ySize, xSize, nFrames);
trueImg = zeros(ySize, xSize, nCells);
trueIntens = zeros(nCells, nFrames);

if length(cellSize) == 1;
    cellSize = cellSize*ones(nCells, 1);
end;
if length(spikeExponent) == 1;
    spikeExponent = spikeExponent*ones(nCells, 1);
end;
if length(cellBrightness) == 1;
    cellBrightness = cellBrightness*ones(nCells, 1);
end;


for j = 1:nCells
    timeTrace = cellBrightness(j) + rand(1,1,nFrames).^spikeExponent(j);
    cellImg = zeros(ySize, xSize);
    cellImg(cellCoords(j,1), cellCoords(j,2)) = 1;
    cellImg = imdilate(cellImg, strel('ball', cellSize(j), 1, 0));
    mov = mov + repmat(cellImg, [1, 1. nFrames]).*repmat(timeTrace, [ySize, xSize, 1]);
    trueImg(:,:,j) = cellImg;
    trueIntens(j,:) = squeeze(timeTrace)';
end;

    
    

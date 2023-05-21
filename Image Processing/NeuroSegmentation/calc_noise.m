function [spaceNoise, timeNoise, spaceErrMap, timeErrMap] = calc_noise(calcCellImg, trueCellImg, calcTimeTrace, trueTimeTrace)
% function [spaceNoise, timeNoise, spaceErrMap, timeErrMap] = calc_noise(calcCellImg, trueCellImg, calcTimeTrace, trueTimeTrace)
% Calculates the accuracy with which a segmentation algorithm extracts the
% spatial and temporal components of a neuron firing trace.  
% Can omit either the space or time components of the calculation by giving
% [] as input argument.
%
% Allows for arbitrary scaling and offset between calculated and true
% signals.  Corrects for these with a linear least-squares fit.  
%
% Inputs:
% calcCellImg: calculated spatial map of the cell
% trueCellImg: underlying true spatial map of the cell
% calcTimeTrace: calcuated intensity trace
% trueTimeTrace: true time trace.
%
% Outputs:
% Noise parameters are: Var(true - calc)/Var(true).
% Error parameters are: true - calc
% In both cases, these are calculated after the linear mapping of calc onto
% true.
%
% AEC 11/16/2014

% make both time traces into column vectors
if ~(isempty(calcTimeTrace) | isempty(trueTimeTrace));  % If time variables given

    if isrow(calcTimeTrace)
        calcTimeTrace = calcTimeTrace';
    end;
    if isrow(trueTimeTrace)
        trueTimeTrace = trueTimeTrace';
    end;

    pTime = polyfit(calcTimeTrace, trueTimeTrace, 1);

    % make the scaled versions
    calcTimeScaled = pTime(1)*calcTimeTrace + pTime(2);

    timeNoise = var(calcTimeScaled - trueTimeTrace)/var(trueTimeTrace);
    timeErrMap = calcTimeScaled - trueTimeTrace;
else
    timeNoise = [];
    timeErrMap = [];
end;

if ~(isempty(calcCellImg) | isempty(trueCellImg));  % If space variables given    
    pSpace = polyfit(calcCellImg(:), trueCellImg(:), 1);
    calcCellScaled = pSpace(1)*calcCellImg + pSpace(2);
    spaceNoise = var(calcCellScaled(:) - trueCellImg(:))/var(trueCellImg(:));
    spaceErrMap = calcCellScaled - trueCellImg;
else
    spaceNoise = [];
    spaceErrMap = [];
end;

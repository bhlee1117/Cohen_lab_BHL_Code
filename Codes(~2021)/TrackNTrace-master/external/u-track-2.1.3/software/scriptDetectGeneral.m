
%% movie information
%
% Copyright (C) 2014 LCCB 
%
% This file is part of u-track.
% 
% u-track is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% u-track is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with u-track.  If not, see <http://www.gnu.org/licenses/>.
% 
% 
movieParam.imageDir = 'C:\Users\user\Desktop\�����\U track Practice from 2014_12_29'; %directory where images are
movieParam.filenameBase = 'old_gfp_004_proj_1'; %image file name base
movieParam.firstImageNum = 0; %number of first image in movie
movieParam.lastImageNum = 300; %number of last image in movie
movieParam.digits4Enum = 3; %number of digits used for frame enumeration (1-4).

%% detection parameters
detectionParam.psfSigma = 1; %point spread function sigma (in pixels)
detectionParam.testAlpha = struct('alphaR',0.05,'alphaA',0.2,'alphaD',0.2,'alphaF',0); %alpha-values for detection statistical tests
detectionParam.visual = 1; %1 to see image with detected features, 0 otherwise
detectionParam.doMMF = 1; %1 if mixture-model fitting, 0 otherwise
detectionParam.bitDepth = 16; %Camera bit depth
detectionParam.alphaLocMax = 0.05; %alpha-value for initial detection of local maxima
detectionParam.numSigmaIter = 0; %maximum number of iterations for PSF sigma estimation
detectionParam.integWindow = 0; %number of frames before and after a frame for time integration

detectionParam.calcMethod = 'g';

% %absolute background info and parameters...
% background.imageDir = 'C:\kjData\Galbraiths\data\alphaVY773AandCellEdge\140109_Cs1C4_Y773A\bgAlphaVY773A\';
% background.filenameBase = 'crop_140109_Cs1C4_mEos2AvBeta3Y773A_';
% background.alphaLocMaxAbs = 0.01;
% detectionParam.background = background;

% detectionParam.maskLoc = 'C:\kjData\Javitch\140115_data\40ms1-5mWPd80G500\maskTest.tif';

%% additional input

%saveResults
saveResults.dir = 'C:\kjData\test\140825_Sungsoo\'; %directory where to save input and output
saveResults.filename = 'detectionMaxProject.mat'; %name of file where input and output are saved
% saveResults = 0;

%verbose state
verbose = 1;

%% run the detection function
[movieInfo,exceptions,localMaxima,background,psfSigma] = ...
    detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults,verbose);

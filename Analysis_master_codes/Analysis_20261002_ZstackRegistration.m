%% Section 0 : Setup
clear
clc;

% Remove any shadow of the built-in 'graph' (old cvx/sdpt3 defines graph.m,
% which breaks the centroid clustering in extractVoltronST / extractGluSNFR3).
gshadow = which('graph','-all');
for gi = 1:numel(gshadow)
    if ~contains(gshadow{gi}, fullfile('toolbox','matlab'))   % keep MATLAB's built-in
        rmpath(fileparts(gshadow{gi}));
        fprintf('Removed shadowing graph.m folder from path: %s\n', fileparts(gshadow{gi}));
    end
end

[~, ~, raw] = xlsread('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Prism_V2+Glu_Data_Arrangement_V2.xlsx', 'Sheet1', 'C5:AA45');
load('/Volumes/BHL18TB_D2/20260203_SD_V2+iGluSNFR4/25X_transformationMatrix.mat');   % tformReg
fpath           = raw(:,1)';
V2moviemaxTime  = 15000;   % voltage movie chunk length (frames)
GlumoviemaxTime = 5000;    % glutamate movie chunk length (frames)
o_laser         = 200001;  % 607/488 modulation onset (DAQ sample index, ~2 s after acq start)
mod488          = 200001;
foi             = [14:19];   % files of interest
set(0,'DefaultFigureWindowStyle','docked')
StructureFilePath=raw(:,6);

%% Section 1 : Options
% Structure z-stacks: one xy position, frames1.tiff (blue) and, when present, frames1_Or.tiff (orange)
opt.blueFile     = 'frames1.tiff';
opt.orangeFile   = 'frames1_Or.tiff';
opt.zReg         = 'translation';   % inter-z registration (StackReg-style): 'translation' or 'rigid'
opt.zRegMaxShift = 150;             % px; larger slice-to-slice steps are treated as failures
opt.zRegMinPeak  = 0.03;            % min phase-correlation peak to accept a slice-to-slice step
opt.residualSearch_px = 40;         % search window for the orange->blue offset left after tformReg
opt.bgRadius_px  = 50;              % Fiji 'Subtract Background' rolling-ball radius (px), orange only
opt.fijiExe      = '/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx';
opt.pixel_um     = [];              % um/px written to the TIFFs ([] = uncalibrated; not in output_data.mat)
opt.zStep_um     = [];              % z step written to the TIFFs ([] = uncalibrated)

%% Section 2 : Register each structure z-stack
% Structure paths of the files of interest; several recordings can share one z-stack (f16/f17)
sDir = cell(numel(foi),1);
for i = 1:numel(foi)
    p = StructureFilePath{foi(i)};
    if ~(ischar(p) || isstring(p)), continue; end
    p = char(p);
    if ~isfolder(p), p = fileparts(p); end           % entry may point at a file inside the folder
    sDir{i} = p;
end
has = ~cellfun(@isempty, sDir);
[uDir, ~, grp] = unique(sDir(has), 'stable');
fUse = foi(has);

for u = 1:numel(uDir)
    tic;
    d = uDir{u};
    [~, name] = fileparts(d);
    outDir = d;                                      % outputs go into the structure folder itself
    fprintf('%s (files %s)\n', name, mat2str(fUse(grp==u)));

    % Blue: inter-z registration from the most structured slice outwards
    B = readTiffStack(fullfile(d, opt.blueFile));
    [B, zShiftB, anchor] = stackReg(B, opt.zReg, opt.zRegMaxShift, opt.zRegMinPeak);
    fprintf('  blue: %d slices, anchor %d, max |shift| %.1f px\n', size(B,3), anchor, max(vecnorm(zShiftB(:,1:2),2,2)));
    writeStackTiff(B, fullfile(outDir,[name '_blue_reg.tif']), opt.pixel_um, opt.zStep_um);
    Reg = struct('structureFolder',d, 'files',fUse(grp==u), 'zShiftBlue',zShiftB, 'anchor',anchor, 'opt',opt);

    % Orange: camera transform into blue coordinates, own inter-z registration anchored at blue's
    % anchor slice (same moment), residual offset correction, Fiji background subtraction
    if exist(fullfile(d, opt.orangeFile), 'file')
        O = readTiffStack(fullfile(d, opt.orangeFile));
        roi = cameraROI(d);
        O = orangeToBlue(O, tformReg, roi, roi);
        [O, zShiftO] = stackReg(O, opt.zReg, opt.zRegMaxShift, opt.zRegMinPeak, anchor);
        [res, okRes, pkRes] = channelOffset({max(B,[],3)}, {max(O,[],3)}, [0 0], opt.residualSearch_px);
        if okRes && any(res)
            O = imtranslate(O, res, 'FillValues', 0);
        end
        fprintf('  orange: max |shift| %.1f px, residual offset after tformReg [dcol drow] = %s (peak %.3f%s)\n', ...
            max(vecnorm(zShiftO(:,1:2),2,2)), mat2str(res), pkRes, repmat(', not applied',1,~okRes));
        tmp = [tempname '.tif'];
        writeStackTiff(O, tmp, opt.pixel_um, opt.zStep_um);
        fijiSubtractBackground(tmp, fullfile(outDir,[name '_orange_reg_bgsub.tif']), opt.bgRadius_px, opt.fijiExe);
        delete(tmp);
        Reg.zShiftOrange = zShiftO; Reg.cameraROI = roi; Reg.tformReg = tformReg;
        Reg.residualOffset = res; Reg.residualApplied = okRes && any(res); Reg.residualPeak = pkRes;
    else
        fprintf('  no %s: blue only\n', opt.orangeFile);
    end
    save(fullfile(outDir,[name '_RegInfo.mat']), 'Reg');
    fprintf('  saved %s (%.0f s)\n', outDir, toc);
end

%% Local functions
function S = readTiffStack(f)
info = imfinfo(f);
S = zeros(info(1).Height, info(1).Width, numel(info), 'uint16');
for z = 1:numel(info), S(:,:,z) = imread(f, z, 'Info', info); end
end

function roi = cameraROI(d)
% Blue camera ROI [x0 w y0 h] from output_data.mat; the orange camera uses the same ROI
ws = warning('off','all'); S = load(fullfile(d,'output_data.mat'),'Device_Data'); warning(ws);
roi = [];
for k = 1:numel(S.Device_Data)
    v = S.Device_Data{k};
    if isstruct(v) && isfield(v,'deviceType') && strcmp(v.deviceType,'Camera'), roi = double(v.ROI); return; end
end
error('No camera ROI in %s', fullfile(d,'output_data.mat'));
end

function Ob = orangeToBlue(O, tformReg, roiO, roiB)
% Warp orange slices into blue camera coordinates with the sensor-level camera registration
% (as transformCamera_O2B: ROIs place both crops on the 2048 px virtual sensor)
Ro = refFromROI(size(O(:,:,1)), roiO);
Rb = refFromROI(size(O(:,:,1)), roiB);
tInv = invert(tformReg);
Ob = zeros(size(O), 'like', O);
for z = 1:size(O,3)
    Ob(:,:,z) = imwarp(O(:,:,z), Ro, tInv, 'OutputView', Rb, 'FillValues', 0);
end
end

function R = refFromROI(imgSize, ROI)
% ROI = [x0 w y0 h] on the virtual sensor (as in transformCamera.m)
R = imref2d(imgSize, [ROI(1)+0.5, ROI(1)+ROI(2)+0.5], [ROI(3)+0.5, ROI(3)+ROI(4)+0.5]);
end

function [S,shift,anchor] = stackReg(S,type,maxShift,minPeak,anchor)
% StackReg-style alignment (Thevenaz et al. 1998): each slice is aligned to its neighbour towards
% the anchor slice, propagating up and down the stack. Steps are estimated between the raw
% neighbouring slices and chained (A_z = A_neighbour * T_step); registering to the already-warped
% neighbour instead fails at large jumps, where its zero-filled border dominates the correlation.
% Each slice is also matched to the slice two steps towards the anchor and the better peak wins:
% a frame exposed during a fast movement is smeared between two positions, and chaining through
% it can lose most of the jump (155809 tile 1, slice 66: 42 of 88 px).
% shift(z,:) = [dx dy theta(deg)] of slice z relative to the anchor frame.
nZ = size(S,3); shift = zeros(nZ,3);
F = zeros(size(S),'single');
for z = 1:nZ, F(:,:,z) = hpf(S(:,:,z)); end
if nargin < 5 || isempty(anchor)
    sharp = squeeze(std(F,0,[1 2]));
    [~,anchor] = max(movmean(sharp,5));           % most structured part of the stack
end
R = imref2d([size(S,1) size(S,2)]);
for z = [anchor+1:nZ, anchor-1:-1:1]
    d = sign(anchor - z);                         % direction towards the anchor
    best = -Inf; A = prevA(shift(z+d,:));         % fallback: no motion relative to the neighbour
    for nb = [z+d, z+2*d]
        if (nb-anchor)*d > 0, continue; end       % past the anchor: not registered yet
        [tf,pk] = imregcorr(F(:,:,z), F(:,:,nb), type, 'Window', true);
        T = tf.A;
        if pk < minPeak || norm(T(1:2,3)) > maxShift || pk <= best, continue; end
        best = pk; A = prevA(shift(nb,:)) * T;
    end
    shift(z,:) = [A(1,3) A(2,3) atan2d(A(2,1),A(1,1))];
end
for z = 1:nZ
    if all(shift(z,:) == 0), continue; end
    S(:,:,z) = imwarp(S(:,:,z), affinetform2d(prevA(shift(z,:))), 'linear', 'OutputView', R, 'FillValues', 0);
end
end

function A = prevA(sh)
% Rebuild the transform of the neighbouring slice (registration failed for this one)
A = [cosd(sh(3)) -sind(sh(3)) sh(1); sind(sh(3)) cosd(sh(3)) sh(2); 0 0 1];
end

function F = hpf(I)
% Band-pass for registration; the median filter removes camera hot pixels, which are
% identical in every tile and otherwise pull the correlation peak to zero shift
I = medfilt2(single(I),[3 3],'symmetric');
F = imgaussfilt(I,1.5) - imgaussfilt(I,15);
end

function [off, ok, pk] = channelOffset(Mb, Mo, guess, L)
% Camera offset [dcol drow] of flipped orange relative to blue (orange pixel p -> blue p + off) from
% the per-tile MIPs, searched within +-L px of guess. Correlation maps of all tiles are summed so the
% fixed camera offset accumulates while false matches (partly different labelled structures) cancel.
[h,w] = size(Mb{1}); nT = numel(Mb);
m = L + max(abs(guess));
C = zeros(2*L+1);
for t = 1:nT
    X = normxcorr2(hpf(Mo{t}(m+1:h-m, m+1:w-m)), hpf(Mb{t}));   % template top-left at (m+1,m+1) = lag 0
    rr = h-m + guess(2) + (-L:L); cc = w-m + guess(1) + (-L:L);
    C = C + X(rr, cc);
end
C = C / nT;
[pk,ix] = max(C(:)); [a,b] = ind2sub(size(C),ix);
off = guess + [b a] - L - 1;
ok = pk > 0.05 && a > 1 && a < 2*L+1 && b > 1 && b < 2*L+1;   % peak inside the window
end

function fijiSubtractBackground(inFile, outFile, radius, fijiExe)
% Runs Fiji headless: Process > Subtract Background (rolling ball, whole stack) on inFile -> outFile
if ~exist(fijiExe,'file')
    warning('Fiji not found at %s; saving without background subtraction', fijiExe);
    copyfile(inFile, strrep(outFile,'_bgsub','')); return
end
mac = [tempname '.ijm'];
fid = fopen(mac,'w');
fprintf(fid, 'setBatchMode(true);\nopen("%s");\n', inFile);
fprintf(fid, 'run("Subtract Background...", "rolling=%g stack");\n', radius);
fprintf(fid, 'saveAs("Tiff", "%s");\nclose("*");\n', outFile);
fclose(fid);
if exist(outFile,'file'), delete(outFile); end
[st,out] = system(sprintf('"%s" --headless -batch "%s"', fijiExe, mac));
delete(mac);
if st ~= 0 || ~exist(outFile,'file')
    warning('Fiji background subtraction failed for %s:\n%s', inFile, out);
else
    fprintf('  background subtracted (rolling ball %g px)\n', radius);
end
end

function writeStackTiff(S, fname, pix_um, dz_um)
% 16-bit multi-page TIFF in ImageJ's layout: header, first IFD + ImageJ description, then all
% pixel data contiguous, then the remaining IFDs. ImageJ/Fiji reads stacks whose description
% says images=N as one contiguous block, so a writer that interleaves IFDs between slices
% (MATLAB's Tiff) makes every slice appear progressively shifted in Fiji.
[H,W,nZ] = size(S);
nb = H*W*2;
if exist(fname,'file'), delete(fname); end
if nb*nZ > 3.9e9                                   % classic TIFF offsets are 32-bit
    writeBigTiff(S, fname, pix_um); return
end
% pix_um / dz_um may be empty (uncalibrated): no unit or spacing then
cal = '';
if ~isempty(pix_um), cal = 'unit=micron\n'; end
if ~isempty(dz_um), cal = [cal sprintf('spacing=%g\\n', dz_um)]; end
desc = [uint8(sprintf(['ImageJ=1.54f\nimages=%d\nslices=%d\n' cal 'loop=false\n'], nZ, nZ)) 0];
if isempty(pix_um), res = uint32([1 1]); else, res = uint32([round(1e6/pix_um) 1e6]); end   % px per um
n0 = 14; n1 = 13;                                  % IFD entries: first page / other pages
ifd0 = 8;
descOff = ifd0 + 2 + 12*n0 + 4;
resOff  = descOff + numel(desc); resOff = resOff + mod(resOff,2);
dataOff = resOff + 16;
ifdOff  = dataOff + nZ*nb + (0:nZ-2)*(2 + 12*n1 + 4);   % IFDs of pages 2..nZ after the data
fid = fopen(fname,'w','l');                        % written strictly in file order
fwrite(fid, uint8('II'), 'uint8'); fwrite(fid, 42, 'uint16'); fwrite(fid, ifd0, 'uint32');
for z = 1:nZ
    nxt = 0; if z < nZ, nxt = ifdOff(z); end
    E = [254 4 1 0; 256 4 1 W; 257 4 1 H; 258 3 1 16; 259 3 1 1; 262 3 1 1];
    if z == 1, E = [E; 270 2 numel(desc) descOff]; end %#ok<AGROW>
    E = [E; 273 4 1 dataOff+(z-1)*nb; 277 3 1 1; 278 4 1 H; 279 4 1 nb; ...
         282 5 1 resOff; 283 5 1 resOff+8; 296 3 1 1]; %#ok<AGROW>
    fwrite(fid, size(E,1), 'uint16');
    for k = 1:size(E,1)
        fwrite(fid, E(k,1), 'uint16'); fwrite(fid, E(k,2), 'uint16'); fwrite(fid, E(k,3), 'uint32');
        if E(k,2) == 3, fwrite(fid, [E(k,4) 0], 'uint16'); else, fwrite(fid, E(k,4), 'uint32'); end
    end
    fwrite(fid, nxt, 'uint32');
    if z == 1
        fwrite(fid, desc, 'uint8');
        fwrite(fid, zeros(1, resOff-descOff-numel(desc)), 'uint8');   % word alignment
        fwrite(fid, [res res], 'uint32');
        for zz = 1:nZ, fwrite(fid, S(:,:,zz).', 'uint16'); end   % row-major pixel data
    end
end
pos = ftell(fid);
fclose(fid);
if nZ > 1 && pos ~= ifdOff(end) + 2 + 12*n1 + 4, error('writeStackTiff: layout mismatch in %s', fname); end
end

function writeBigTiff(S, fname, pix_um)
% > 4 GB: BigTIFF via MATLAB's Tiff (open in Fiji with Bio-Formats); pixel size in the resolution tags
t = Tiff(fname,'w8');
tag.ImageLength = size(S,1); tag.ImageWidth = size(S,2);
tag.Photometric = Tiff.Photometric.MinIsBlack; tag.BitsPerSample = 16;
tag.SamplesPerPixel = 1; tag.SampleFormat = Tiff.SampleFormat.UInt;
tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tag.Compression = Tiff.Compression.None;
if ~isempty(pix_um), tag.XResolution = 1/pix_um; tag.YResolution = 1/pix_um; end
tag.ResolutionUnit = Tiff.ResolutionUnit.None;
for z = 1:size(S,3)
    t.setTag(tag); t.write(S(:,:,z));
    if z < size(S,3), t.writeDirectory(); end
end
t.close();
end


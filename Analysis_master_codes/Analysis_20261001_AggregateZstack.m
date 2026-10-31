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

Snapfolder='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20260902_V2ST_Glu_VR/Snaps';
Tilefolder='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20260902_V2ST_Glu_VR/ZStackTiles';
Savefolder='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20260902_V2ST_Glu_VR/StitchedZstack';
Orangefolder='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20260902_V2ST_Glu_VR';   % <HHMMSS>..._Snap/frames1.bin
DateStr='20260902';                 % snap filenames only carry HHMMSS

opt.sessionKey   = 'Zstack|Ztile';  % regexp on session name to keep (tiling sessions only)
opt.maxGap_s     = 30;              % a pause longer than this splits a session
opt.minImages    = 20;              % ignore blocks with fewer snaps (stray test snaps)
opt.sessions     = [];              % indices of sessions to process (after Section 2); [] = all
% Stage -> camera mapping. Camera columns follow stage +X and rows follow stage +Y.
% The pixel size stored in scan_coordinates.mat (0.7312) does not match the data;
% registration of neighbouring tiles on 20260902 gives ~0.887 um/px.
opt.calibrate    = false;           % true: fit um->px per axis from neighbouring tiles of each session.
                                    % Off: per-session fits use 1-3 pairs and scatter 0.865-0.905 um/px,
                                    % while the optics are the same; one value keeps all outputs on one scale
opt.pixel_um     = 0.887;           % camera pixel size (um/px) used for tile placement and saved metadata
opt.minPeak      = 0.2;             % min normxcorr2 peak to trust a pair
opt.refineShift  = true;            % refine tile offsets by cross-correlation of overlaps
opt.maxShift_px  = 50;              % max correction allowed per tile pair
opt.flatfield    = true;            % BaSiC flat-field correction (Peng et al. 2017; port of BaSiCPy)
opt.darkLevel    = 100;             % camera offset (Orca), subtracted before dividing by the flat-field
opt.basicSize    = 128;             % working resolution for the BaSiC fit
opt.flatEps      = 0.3;             % regularised division gain = f(1+e^2)/(f^2+e^2): 1 at the peak,
                                    % max (1+e^2)/2e ~1.8x at f = e, rolling off in unlit corners
opt.zReg         = 'translation';   % inter-z registration (StackReg-style): 'translation', 'rigid' or 'none'
opt.zRegMaxShift = 150;             % px; larger slice-to-slice jumps are treated as failures. Awake-mouse
                                    % jumps reach ~90 px (155809 tile 1), so keep this well above that
opt.zRegMinPeak  = 0.03;            % min phase-correlation peak to accept a slice-to-slice estimate
opt.blend        = 'feather';       % 'feather' (linear ramp) or 'max'
opt.saveStack    = true;            % outputs: <tag>_blue.tif and <tag>_orange_bgsub.tif (common area only)
% Orange channel: one frame per snap position in frames1.bin, left-right flipped relative to blue.
% Own flat-field and own slice-to-slice z-registration (the cameras do not expose at the same instant,
% so blue's shifts do not fit it), but anchored at blue's anchor slice of each tile and placed at
% blue's tile positions: the cameras are rigidly coupled, so with the same anchor (same moment) the
% tile positions must be the same. Independent tile refinement put tiles up to ~60 px apart.
opt.orange       = true;
opt.orangeKey    = '_Snap';         % folder name postfix
opt.cameraSearch_px = 80;           % search window around opt.cameraOffset for the per-session estimate
opt.orangeBgSub  = true;            % Fiji 'Subtract Background' (rolling ball) on the orange stack only
opt.bgRadius_px  = 50;              % rolling-ball radius in px (Fiji default; 50 px = 44 um)
opt.fijiExe      = '/Applications/Fiji.app/Contents/MacOS/ImageJ-macosx';
opt.cameraOffset = [0 23];          % flipped orange -> blue camera offset [dcol drow]: fallback and search centre
                                    % (all sessions on 20260902 gave [-1..0, 22..24])

if ~exist(Savefolder,'dir'), mkdir(Savefolder); end

%% Section 1 : Group snaps into sessions (name after the timestamp + time gaps)
d = dir(fullfile(Snapfolder,'*.tiff'));
tok = regexp({d.name}, '^(\d{6})(.+)\.tiff$', 'tokens', 'once');
ok = ~cellfun(@isempty,tok);
d = d(ok); tok = vertcat(tok{ok});
tSnap = datetime(strcat(DateStr,tok(:,1)), 'InputFormat','yyyyMMddHHmmss');
sName = tok(:,2);
[tSnap, si] = sort(tSnap); sName = sName(si); d = d(si);

newBlk = [true; ~strcmp(sName(2:end),sName(1:end-1)) | seconds(diff(tSnap)) > opt.maxGap_s];
blkID  = cumsum(newBlk);

Sess = struct('name',{},'files',{},'t',{});
for b = 1:max(blkID)
    id = find(blkID==b);
    if numel(id) < opt.minImages || isempty(regexp(sName{id(1)}, opt.sessionKey, 'once')), continue; end
    Sess(end+1).name = sName{id(1)};
    Sess(end).files  = fullfile(Snapfolder, {d(id).name}');
    Sess(end).t      = tSnap(id);
end
fprintf('Found %d tiling sessions\n', numel(Sess));

%% Section 2 : Match each session to the closest preceding coordinate list
c = dir(fullfile(Tilefolder,'2026*'));
c = c([c.isdir] & arrayfun(@(x) exist(fullfile(Tilefolder,x.name,'scan_coordinates.mat'),'file')==2, c)');
tList = datetime(cellfun(@(s) s(1:15), {c.name}', 'UniformOutput', false), 'InputFormat','yyyyMMdd_HHmmss');
% folders with the same second (e.g. _02) are ordered by suffix so the last export wins
[~,ci] = sortrows([datenum(tList) (1:numel(c))']); c = c(ci); tList = tList(ci);

for s = 1:numel(Sess)
    dt  = seconds(Sess(s).t(1) - tList);
    pre = find(dt >= 0);
    if isempty(pre), error('No coordinate list before session %s', Sess(s).name); end
    nRow = arrayfun(@(k) height(loadScanList(fullfile(Tilefolder,c(k).name))), pre);
    nImg = numel(Sess(s).files);
    good = pre(nRow == nImg);
    if ~isempty(good)
        k = good(find(dt(good)==min(dt(good)),1,'last'));
    else
        k = pre(find(dt(pre)==min(dt(pre)),1,'last'));
        warning('%s: %d snaps but closest list %s has %d positions; using the first min(N).', ...
            Sess(s).name, nImg, c(k).name, nRow(pre==k));
    end
    Sess(s).listFolder = c(k).name;
    [Sess(s).plan, Sess(s).scanInfo] = loadScanList(fullfile(Tilefolder,c(k).name));
    fprintf('%s %-22s (%4d snaps) <- %s (%4d positions, %.0f s before)\n', datestr(Sess(s).t(1),'HH:MM:SS'), ...
        Sess(s).name, nImg, c(k).name, height(Sess(s).plan), dt(k));
end

%% Section 3 : Stitch each session
if ~isempty(opt.sessions), Sess = Sess(opt.sessions); end
for s = 1:numel(Sess)
    tic;
    plan = Sess(s).plan;
    n    = min(height(plan), numel(Sess(s).files));
    plan = plan(1:n,:);
    files = Sess(s).files(1:n);

    % Load all snaps (snap order == plan row order)
    I0 = imread(files{1});
    [h,w] = size(I0);
    Raw = zeros(h,w,n,'like',I0);
    for i = 1:n
        Raw(:,:,i) = imread(files{i});
    end

    tag = [datestr(Sess(s).t(1),'HHMMSS') Sess(s).name];
    outDir = fullfile(Savefolder, tag);
    if ~exist(outDir,'dir'), mkdir(outDir); end

    % Blue channel
    fprintf('  [blue]\n');
    [Stack, St, MIPb] = stitchChannel(Raw, plan, opt, []);
    clear Raw
    Info = struct('session',Sess(s).name, 'snapFiles',{files}, 'listFolder',Sess(s).listFolder, ...
        'scanInfo',Sess(s).scanInfo, 'opt',opt);
    Info = mergeStruct(Info, St);
    save(fullfile(outDir,[tag '_StitchInfo.mat']), 'Info');
    dz = median(abs(diff(St.z_um))); pix = mean(St.pixel_um);
    fprintf('  blue stitched (%.0f s)\n', toc);

    % Orange channel: own flat-field and z-registration, blue's anchors and tile positions
    oDir = '';
    if opt.orange
        [oDir, oSize] = findOrangeFolder(Orangefolder, opt.orangeKey, DateStr, Sess(s).t(1), n);
        if isempty(oDir), warning('%s: no orange %s folder with %d frames found', tag, opt.orangeKey, n); end
    end
    if isempty(oDir)
        if opt.saveStack, writeStackTiff(Stack, fullfile(outDir,[tag '_blue.tif']), pix, dz); end
        clear Stack; continue
    end
    Ro = readBinFrames(fullfile(oDir,'frames1.bin'), oSize, n);
    Ro = flip(Ro,2);                                  % orange camera is mirrored left-right
    fprintf('  [orange] %s\n', oDir);
    [StackO, StO, MIPo] = stitchChannel(Ro, plan, opt, St);
    clear Ro

    % Camera offset (flipped orange -> blue): consensus over tiles, within +-cameraSearch_px of the prior
    [off, okOff, pkOff] = channelOffset(MIPb, MIPo, opt.cameraOffset, opt.cameraSearch_px);
    if ~okOff
        fprintf('  camera offset estimate unclear (peak %.3f); using %s\n', pkOff, mat2str(opt.cameraOffset));
        off = opt.cameraOffset;
    end
    fprintf('  camera offset [dcol drow] = %s (peak %.3f)\n', mat2str(off), pkOff);

    % Common area of both channels: same size and pixel grid, ready to merge in Fiji
    [Bc, Oc, ovRows, ovCols] = cropToOverlap(Stack, StackO, off);
    clear Stack StackO
    InfoOrange = mergeStruct(struct('folder',oDir, 'flip','left-right', 'cameraOffset_px',off, ...
        'cameraOffsetPeak',pkOff, 'cameraOffsetEstimated',okOff), StO);
    Overlap = struct('blueRows',ovRows,'blueCols',ovCols,'orangeRows',ovRows-off(2),'orangeCols',ovCols-off(1),'size',size(Bc));
    save(fullfile(outDir,[tag '_StitchInfo.mat']), 'InfoOrange', 'Overlap', '-append');
    if opt.saveStack
        writeStackTiff(Bc, fullfile(outDir,[tag '_blue.tif']), pix, dz);
        if opt.orangeBgSub
            tmp = [tempname '.tif'];
            writeStackTiff(Oc, tmp, pix, dz);
            fijiSubtractBackground(tmp, fullfile(outDir,[tag '_orange_bgsub.tif']), opt.bgRadius_px, opt.fijiExe);
            delete(tmp);
        else
            writeStackTiff(Oc, fullfile(outDir,[tag '_orange.tif']), pix, dz);
        end
    end
    clear Bc Oc
    fprintf('  saved %s  [%d x %d x %d]  (%.0f s)\n', outDir, numel(ovCols), numel(ovRows), numel(St.z_um), toc);
end

%% Local functions
function flatfield = basicFlatfield(Raw, dark, sz)
% BaSiC flat-field estimation (Peng et al. 2017, Nat Commun; BaSiCPy 'approximate' fit).
% Images (dark level removed) are modelled as I_n = b_n*S + E_n: S is a smooth flat-field,
% sparse in the DCT domain (L1), b_n a per-image baseline and E_n a sparse residual (L1)
% holding the actual structure. Solved by inexact ALM with reweighting of E.
[h,w,n] = size(Raw);
D = zeros(sz,sz,n,'single');
for i = 1:n
    D(:,:,i) = imresize(single(Raw(:,:,i)) - dark, [sz sz], 'bilinear');
end
D = max(D,0);
D = sort(D,3);                                    % per-pixel sorting, as in BaSiC
meanD = mean(D,3); meanD = meanD/mean(meanD(:));
lambda = sum(abs(dct2(meanD)),'all')/400*0.5;     % BaSiC default smoothness weight
Dm = reshape(D,[],n);
weight = ones(size(Dm),'single');
flatLast = ones(sz,'single');
for it = 1:5                                      % reweighting iterations
    [S,~,E] = basicALM(Dm, sz, lambda, weight);
    A = S(:)*mean(Dm,1);                          % low-rank part for normalisation
    weight = 1 ./ (abs(E ./ mean(A(:))) + 0.1);
    weight = weight * numel(weight) / sum(weight(:));
    flat = S / mean(S(:));
    mad = sum(abs(flat(:)-flatLast(:))) / sum(abs(flatLast(:)));
    flatLast = flat;
    if mad < 1e-3, break; end
end
flatfield = imresize(flat, [h w], 'bilinear');
% Normalised to the well-lit peak (BaSiC uses the mean): the centre keeps raw counts and the
% divisor never exceeds 1, so dim edges are brightened, not the centre darkened
flatfield = single(flatfield / prctile(flatfield(:),99.9));
end

function [S,b,E] = basicALM(D, sz, lambda, weight)
% Inexact ALM for D = S(:)*b + E with ||DCT(S)||_1 and ||weight.*E||_1 (inexact_alm_rspca_l1)
n = size(D,2);
normTwo = svds(double(D),1);
mu = 12.5/normTwo; muBar = mu*1e7; rho = 1.5;
dNorm = norm(D,'fro');
W = zeros(sz,'single');                           % DCT coefficients of S
b = ones(1,n,'single');
E = zeros(size(D),'single'); Y = zeros(size(D),'single');
for k = 1:500
    S = idct2(W);
    A = S(:)*b;
    T = reshape(D - A - E + Y/mu, sz, sz, n);
    W = W + dct2(mean(T,3));
    W = max(W - lambda/mu, 0) + min(W + lambda/mu, 0);
    S = idct2(W);
    A = S(:)*b;
    E = E + (D - A - E + Y/mu);
    E = max(E - weight/mu, 0) + min(E + weight/mu, 0);
    R = D - E;
    b = mean(R,1) / mean(R(:)); b(b<0) = 0;     % baseline normalised to mean 1, S carries the scale
    A = S(:)*b;
    Z = D - A - E;
    Y = Y + mu*Z;
    mu = min(mu*rho, muBar);
    if norm(Z,'fro')/dNorm < 1e-6, break; end
end
S = idct2(W);
end

function [plan,info] = loadScanList(folder)
% Older lists only have 'plan' and 'meta'; newer ones also have 'scanInfo'
S = load(fullfile(folder,'scan_coordinates.mat'));
plan = S.plan;
if isfield(S,'scanInfo'), info = S.scanInfo; else, info = S.meta; end
end

function [pairs,ov] = tilePairs(pos,h,w)
% Tile pairs whose footprints overlap, with the overlap box for each
nT = size(pos,1); pairs = []; ov = [];
for a = 1:nT
    for b = a+1:nT
        c0 = max(pos(a,1),pos(b,1)); c1 = min(pos(a,1),pos(b,1)) + w;
        r0 = max(pos(a,2),pos(b,2)); r1 = min(pos(a,2),pos(b,2)) + h;
        if c1-c0 > 30 && r1-r0 > 30
            pairs(end+1,:) = [a b]; %#ok<AGROW>
            ov(end+1,:) = [r0 r1 c0 c1]; %#ok<AGROW>
        end
    end
end
end

function [A,B] = overlapCrop(Ia,Ib,pa,pb,box)
% Crop the shared region of two tiles placed at pa,pb ([col row], px)
pa = round(pa); pb = round(pb); box = round(box);
ra = (box(1):box(2)-1) - pa(2) + 1; ca = (box(3):box(4)-1) - pa(1) + 1;
rb = (box(1):box(2)-1) - pb(2) + 1; cb = (box(3):box(4)-1) - pb(1) + 1;
A = Ia(ra,ca); B = Ib(rb,cb);
end

function F = hpf(I)
% Band-pass for registration; the median filter removes camera hot pixels, which are
% identical in every tile and otherwise pull the correlation peak to zero shift
I = medfilt2(single(I),[3 3],'symmetric');
F = imgaussfilt(I,1.5) - imgaussfilt(I,15);
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

function pxPerUm = calibrateScale(M,pos_um,grid,pxPerUm,minPeak)
% Unconstrained registration of grid neighbours (2x downsampled) -> px per um for each axis.
% An axis without a reliable pair keeps the fallback value.
[h,w] = size(M{1}); f = 0.5;
Ms = cellfun(@(x) imresize(hpf(x),f), M, 'UniformOutput', false);
[hs,ws] = size(Ms{1});
[cc,rr] = meshgrid(1:2*ws-1, 1:2*hs-1);
ovl = (ws-abs(cc-ws)) .* (hs-abs(rr-hs));        % overlap area for each lag
est = {[],[]};
for a = 1:numel(M)
    for b = a+1:numel(M)
        dg = grid(b,:) - grid(a,:);
        if sum(abs(dg)) ~= 1, continue; end      % grid neighbours only
        ax = find(dg);                            % 1: X move, 2: Y move
        X = normxcorr2(Ms{a},Ms{b});
        X(ovl < 0.15*hs*ws) = -1;
        [pk,ix] = max(X(:));
        if pk < minPeak, continue; end
        [yy,xx] = ind2sub(size(X),ix);
        shift = [ws-xx hs-yy]/f;                  % origin of b relative to a, [col row] px
        dStage = pos_um(b,ax) - pos_um(a,ax);
        sz = [w h]; other = 3-ax;
        if abs(shift(other)) > 0.05*sz(other), continue; end   % should move along one axis only
        est{ax}(end+1) = shift(ax)/dStage;
    end
end
for ax = 1:2
    e = est{ax};
    e = e(abs(e/pxPerUm(ax) - 1) < 0.15);          % within 10% of the prior (also fixes the sign)
    if isempty(e), continue; end
    m = median(e);
    e = e(abs(e-m) < 0.1*m);
    pxPerUm(ax) = mean(e);
    fprintf('  axis %d: %d pairs -> %.4f um/px\n', ax, numel(e), 1/pxPerUm(ax));
end
end

function pos = refinePositions(M,pos0,maxShift,minPeak)
% Pairwise shift by normxcorr2 on the overlaps, then least-squares tile positions (tile 1 fixed)
[h,w] = size(M{1}); nT = numel(M);
[pairs,ov] = tilePairs(pos0,h,w);
rows = []; d = []; wt = [];
for p = 1:size(pairs,1)
    a = pairs(p,1); b = pairs(p,2);
    % pad b's crop by maxShift so the template (a's crop) can slide inside it
    box = round(ov(p,:));
    [A,~] = overlapCrop(M{a},M{b},pos0(a,:),pos0(b,:),box);
    boxB = box + [-1 1 -1 1]*maxShift;
    pb = round(pos0(b,:));
    boxB([1 2]) = min(max(boxB([1 2]), pb(2)), pb(2)+h);
    boxB([3 4]) = min(max(boxB([3 4]), pb(1)), pb(1)+w);
    Bp = M{b}((boxB(1):boxB(2)-1)-pb(2)+1, (boxB(3):boxB(4)-1)-pb(1)+1);
    A = hpf(A); Bp = hpf(Bp);
    if std(A(:))==0 || any(size(Bp)<size(A)), continue; end
    X = normxcorr2(A,Bp);
    [pk,ix] = max(X(:)); [yy,xx] = ind2sub(size(X),ix);
    % location of A's top-left inside Bp, relative to where the nominal positions put it
    dy = (yy - size(A,1)) - (round(box(1)) - boxB(1));
    dx = (xx - size(A,2)) - (round(box(3)) - boxB(3));
    if abs(dx) >= maxShift || abs(dy) >= maxShift || pk < minPeak, continue; end
    % content at A-pixel lies dy,dx further into B than nominal => move b by -[dx dy]
    rows(end+1,:) = [a b]; d(end+1,:) = [-dx -dy]; wt(end+1,1) = pk; %#ok<AGROW>
end
pos = pos0;
if isempty(rows), return; end
L = zeros(size(rows,1)+1, nT);
for k = 1:size(rows,1), L(k,rows(k,1)) = -1; L(k,rows(k,2)) = 1; end
L(end,1) = 1;                                % anchor tile 1
Wd = diag([wt; 1]);
for dim = 1:2
    nominal = pos0(rows(:,2),dim) - pos0(rows(:,1),dim);
    rhs = [nominal + d(:,dim); pos0(1,dim)];
    pos(:,dim) = (L'*Wd*L + 1e-6*eye(nT)) \ (L'*Wd*rhs + 1e-6*pos0(:,dim));
end
end

function Stack = assembleStack(Raw, zIdx, pos, wTile, nZ, blend)
% Blend tiles plane by plane at integer positions pos ([col row], 0-based); 0 = no data
[ht,wt] = size(wTile); nT = numel(zIdx);
H = max(pos(:,2)) + ht; W = max(pos(:,1)) + wt;
Stack = zeros(H,W,nZ,'uint16');
for z = 1:nZ
    P = zeros(H,W,'single'); Wsum = zeros(H,W,'single');
    for t = 1:nT
        if z > numel(zIdx{t}), continue; end
        rr = pos(t,2) + (1:ht); cc = pos(t,1) + (1:wt);
        T = single(Raw(:,:,zIdx{t}(z)));
        if strcmpi(blend,'max')
            P(rr,cc) = max(P(rr,cc), T);
        else
            wv = wTile .* imerode(T>0, ones(3));
            P(rr,cc) = P(rr,cc) + T.*wv;
            Wsum(rr,cc) = Wsum(rr,cc) + wv;
        end
    end
    if ~strcmpi(blend,'max'), P = P ./ max(Wsum, eps('single')); end
    Stack(:,:,z) = uint16(P);
end
end

function [Stack, St, MIPr] = stitchChannel(Raw, plan, opt, ref)
% One channel: flat-field, inter-z registration, um->px calibration, tile refinement, assembly.
% Raw is h x w x n in plan row order (snap order). With a reference channel ref (St of blue), the
% z-registration of each tile is anchored at ref.zAnchor and ref's pixel size and tile positions are
% used (no calibration or refinement of its own). MIPr: per-tile MIPs used for registration.
[h,w,n] = size(Raw);
flatfield = ones(h,w,'single');
if opt.flatfield
    flatfield = basicFlatfield(Raw, opt.darkLevel, opt.basicSize);
    gain = max(flatfield,0) * (1+opt.flatEps^2) ./ (flatfield.^2 + opt.flatEps^2);
    for i = 1:n
        Raw(:,:,i) = uint16((single(Raw(:,:,i)) - opt.darkLevel) .* gain + opt.darkLevel);
    end
    fprintf('  flat-field range %.2f-%.2f\n', min(flatfield(:)), max(flatfield(:)));
end

tiles = unique(plan.Tile);
nT = numel(tiles);
nZ = max(plan.ZIndex);
zIdx = cell(nT,1); tilePos_um = zeros(nT,2);      % [X Y] stage position of each tile
tileGrid = zeros(nT,2);                           % [TileX TileY]
MIP = cell(nT,1); zShift = cell(nT,1); zAnchor = zeros(nT,1);
for t = 1:nT
    r = find(plan.Tile == tiles(t));
    [~,o] = sort(plan.ZIndex(r)); r = r(o);
    zIdx{t} = r;
    tilePos_um(t,:) = [plan.X_um(r(1)) plan.Y_um(r(1))];
    tileGrid(t,:)   = [plan.TileX(r(1)) plan.TileY(r(1))];
    % Brain motion between z-planes (awake mouse): align slices within the tile stack
    if ~strcmpi(opt.zReg,'none')
        if isempty(ref), a0 = []; else, a0 = ref.zAnchor(t); end
        [Raw(:,:,r), zShift{t}, zAnchor(t)] = stackReg(Raw(:,:,r), opt.zReg, opt.zRegMaxShift, opt.zRegMinPeak, a0);
        fprintf('  tile %d: z-reg anchor %d, max |shift| %.1f px\n', tiles(t), zAnchor(t), max(vecnorm(zShift{t}(:,1:2),2,2)));
    end
    MIP{t} = max(Raw(:,:,r),[],3);
end
% Registration images: unlit corners (amplified noise after flat-field) set to a constant
MIPr = MIP;
for t = 1:nT
    MIPr{t}(flatfield < opt.flatEps) = median(MIP{t}(flatfield >= opt.flatEps));
end
z_um = accumarray(plan.ZIndex, plan.Z_um, [nZ 1], @mean);

% Stage um -> pixel scale per axis ([col row] = stage [X Y] .* pxPerUm)
pxPerUm = [1 1]./opt.pixel_um;                    % pixel_um: scalar or [X Y]
if ~isempty(ref), pxPerUm = 1./ref.pixel_um; end
if opt.calibrate && nT > 1 && isempty(ref)
    pxPerUm = calibrateScale(MIPr, tilePos_um, tileGrid, pxPerUm, opt.minPeak);
end
pix = 1./pxPerUm;
fprintf('  pixel size: X %.4f, Y %.4f um/px\n', pix);
pos0 = (tilePos_um - min(tilePos_um,[],1)) .* pxPerUm;   % [col row], px

% Refine tile offsets from pairwise overlap registration (global least squares)
pos = pos0;
if ~isempty(ref)
    pos = ref.tilePos_px;                         % reference channel's tile positions
else
    if opt.refineShift && nT > 1
        pos = refinePositions(MIPr, pos0, opt.maxShift_px, opt.minPeak);
        fprintf('  refinement moved tiles by up to %.1f px\n', max(vecnorm(pos-pos0,2,2)));
    end
    pos = round(pos - min(pos,[],1));
end

% Assemble plane by plane; 0 marks pixels without data (borders after z-registration)
[ht,wt] = size(MIP{1});
wTile = blendWeight(ht,wt) .* flatfield;          % overlaps favour the better-lit tile
Stack = assembleStack(Raw, zIdx, pos, wTile, nZ, opt.blend);
St = struct('pixel_um',pix, 'z_um',z_um, 'zShift',{zShift}, 'zAnchor',zAnchor, 'flatfield',flatfield, ...
    'tiles',tiles, 'tileSize',[ht wt], 'tilePos_px',pos, 'tilePos_um',tilePos_um);
end

function [Bc, Oc, rB, cB] = cropToOverlap(B, O, off)
% Crop B and O to their common area; O pixel (r,c) sits at B pixel (r+drow, c+dcol), off = [dcol drow].
% rB, cB: row/column ranges of the overlap in B (in O: rB-drow, cB-dcol).
rB = max(1, 1+off(2)) : min(size(B,1), size(O,1)+off(2));
cB = max(1, 1+off(1)) : min(size(B,2), size(O,2)+off(1));
Bc = B(rB, cB, :);
Oc = O(rB-off(2), cB-off(1), :);
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

function S = mergeStruct(S, T)
f = fieldnames(T);
for k = 1:numel(f), S.(f{k}) = T.(f{k}); end
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

function [oDir, sz] = findOrangeFolder(root, key, DateStr, t0, n)
% Orange acquisition folder (<HHMMSS><name><key>) created closest before the first blue snap
oDir = ''; sz = [];
d = dir(root); d = d([d.isdir]);
tok = regexp({d.name}, ['^(\d{6}).*' key '$'], 'tokens', 'once');
d = d(~cellfun(@isempty,tok)); tok = vertcat(tok{~cellfun(@isempty,tok)});
if isempty(d), return; end
tO = datetime(strcat(DateStr,tok(:,1)), 'InputFormat','yyyyMMddHHmmss');
dt = seconds(t0 - tO);
cand = find(dt >= -5 & arrayfun(@(x) exist(fullfile(root,x.name,'frames1.bin'),'file')==2, d));
[~,o] = sort(dt(cand)); cand = cand(o);
for k = cand'
    f = fullfile(root, d(k).name);
    szk = camFrameSize(f);
    b = dir(fullfile(f,'frames1.bin'));
    if b.bytes == 2*prod(szk)*n, oDir = f; sz = szk; return; end
end
end

function sz = camFrameSize(f)
% [rows cols] of the orange camera from output_data.mat (ROI = [x0 w y0 h]); 1000x1000 if unreadable
sz = [1000 1000];
try
    ws = warning('off','all'); S = load(fullfile(f,'output_data.mat'),'Device_Data'); warning(ws);
    for k = 1:numel(S.Device_Data)
        v = S.Device_Data{k};
        if isstruct(v) && isfield(v,'deviceType') && strcmp(v.deviceType,'Camera')
            sz = double(v.ROI([4 2])); return;
        end
    end
catch
end
end

function M = readBinFrames(fname, sz, n)
% frames1.bin: little-endian uint16, row-major frames (as readBinMov_times)
fid = fopen(fname,'r','l');
M = fread(fid, [prod(sz) n], '*uint16');
fclose(fid);
M = permute(reshape(M, sz(2), sz(1), []), [2 1 3]);
end

function W = blendWeight(h,w)
% Linear ramp to the tile edges so overlaps fade smoothly
[cc,rr] = meshgrid(1:w,1:h);
W = single(min(min(cc,w-cc+1), min(rr,h-rr+1)));
W = W / max(W(:));
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
desc = [uint8(sprintf('ImageJ=1.54f\nimages=%d\nslices=%d\nunit=micron\nspacing=%g\nloop=false\n', nZ, nZ, dz_um)) 0];
res  = uint32([round(1e6/pix_um) 1e6]);            % pixels per micron as a rational
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
tag.XResolution = 1/pix_um; tag.YResolution = 1/pix_um;
tag.ResolutionUnit = Tiff.ResolutionUnit.None;
for z = 1:size(S,3)
    t.setTag(tag); t.write(S(:,:,z));
    if z < size(S,3), t.writeDirectory(); end
end
t.close();
end

cd 'X:\Lab\Labmembers\Dingchang Lin\Data\HCBI\20220109\';
% cd 'G:\Cohen lab\HCBI imaging\20220109';

dList = dir('int_*');
dList = dList([dList(:).isdir]);
nDir = length(dList);
ySize = 2048; xSize = 2048;

theta = [0:179];
fibW = 2; % fiber width, pixels
dW = 5*fibW; % Width of segment to take around each fiber
minL = 50; % minimum fiber length, pixels
maxL = 600; % maximum fiber length, pixels.  Use an even number.

colMat = [0 1 1; 0 1 0; 1 0.7 0; 1 0 0];

% Get one picture from each directory just to check that it looks ok.
for d = 1:nDir;  % Go through all the directories
    cd(dList(d).name);
    fList = dir('*low_Unmixing.czi');
    nFiles = length(fList);
    for f = 1:nFiles;  % go through all the unmixed files in each directory
        fName = fList(f).name;
        dat = bfopen(fName);  % open one Unmixed image.
        nChan = size(dat, 2)
        ['number of channels ' num2str(nChan)]
        omeMeta = dat{1,4};
        voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER); % in µm
        pixSize = voxelSizeX.doubleValue();
        [d pixSize]
        [ySize, xSize] = size(dat{1}{1,1});  % image dimensions
        colorImg = zeros(ySize, xSize, nChan);
        for j = 1:nChan;
            colorImg(:,:,j) = mat2gray(double(dat{1}{j,1}));  % make a 3-color image
        end;
        figure(d); clf
        multicolor(colorImg, colMat)  % show a multicolor image.  Optional to specify the color of each channel
    end;
cd ..
end;

%extract all the fibers at lowmag
for d = 1:nDir;  % Go through all the directories
    cd(dList(d).name);
    fList = dir('*low_Unmixing.czi');
    nFiles = length(fList);
    clear fIdx
    for f = 1:nFiles;
        fIdx(f) = str2num(fList(f).name(1:end-16));
    end;
    [~, idx] = sort(fIdx);
    fList = fList(idx);
    clear('allFibs', 'ySize', 'xSize', 'pixSize', 'nChan')
    for f = 1:nFiles;  % go through all the unmixed files in each directory
        fName = fList(f).name;
        dat = bfopen(fName);  % open each Unmixed image.
        omeMeta = dat{1,4};
        voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER); % in µm
        pixSize(f) = voxelSizeX.doubleValue();
        for j = 1:4;
            subplot(2,2,j);
            imshow2(double(dat{1}{j,1}), []);
        end;
        [ySize(f), xSize(f)] = size(dat{1}{1,1});  % image dimensions
        nChan(f) = size(dat, 2);
        colorImg = zeros(ySize(f), xSize(f), nChan(f));
        for j = 1:nChan(f);
            colorImg(:,:,j) = mat2gray(double(dat{1}{j,1}));  % make a color image
        end;
        h = figure(f); clf
        h.WindowState = 'maximized';
        multicolor(colorImg, colMat); hold all
        refImg = max(colorImg(:,:,2:4),[], 3);  % define a one-color reference image for finding fibers.
        clear fibs
        tic
        fibs = findFib(refImg, colorImg, fibW, minL, maxL, dW);
        toc
        figure(f); hold all
        nFib = length(fibs);
        for j = 1:nFib;
            plot(fibs(j).X, fibs(j).Y, 'm-');
        end;
        hold off
%         saveas(gca, [fName(1:end-16) '_fibs.fig'])
        'saving...'
        tic
        saveas(gca, [fName(1:end-16) '_fibs.png'])        
        toc
        allFibs{f} = fibs;
    end;
    cd ..
    save([dList(d).name '_allImgFibs.mat'], 'pixSize', 'ySize', 'xSize', 'allFibs', 'fibW', 'minL', 'maxL', 'dW')
    close all
end;

%extract all the fibers at medium mag
theta = [0:179];
fibW = 10; % fiber width, pixels
dW = 5*fibW; % Width of segment to take around each fiber
minL = 50; % minimum fiber length, pixels
maxL = 1200; % maximum fiber length, pixels.  Use an even number.
for d = 1:nDir;  % Go through all the directories
    cd(dList(d).name);
    fList = dir('*med_Unmixing.czi');
    nFiles = length(fList);
    clear fIdx
    for f = 1:nFiles;
        fIdx(f) = str2num(fList(f).name(1:end-16));
    end;
    [~, idx] = sort(fIdx);
    fList = fList(idx);
    clear('allFibs', 'ySize', 'xSize', 'pixSize', 'nChan')
    for f = 1:nFiles;  % go through all the unmixed files in each directory
        fName = fList(f).name;
        dat = bfopen(fName);  % open each Unmixed image.
        omeMeta = dat{1,4};
        voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER); % in µm
        pixSize(f) = voxelSizeX.doubleValue();
        for j = 1:4;
            subplot(2,2,j);
            imshow2(double(dat{1}{j,1}), []);
        end;
        [ySize(f), xSize(f)] = size(dat{1}{1,1});  % image dimensions
        nChan(f) = size(dat, 2);
        colorImg = zeros(ySize(f), xSize(f), nChan(f));
        for j = 1:nChan(f);
            colorImg(:,:,j) = mat2gray(double(dat{1}{j,1}));  % make a color image
        end;
        h = figure(f); clf
        h.WindowState = 'maximized';
        multicolor(colorImg, colMat); hold all
        refImg = max(colorImg(:,:,2:4),[], 3);  % define a one-color reference image for finding fibers.
        clear fibs
        tic
        fibs = findFib(refImg, colorImg, fibW, minL, maxL, dW);
        toc
        figure(f); hold all
        nFib = length(fibs);
        for j = 1:nFib;
            plot(fibs(j).X, fibs(j).Y, 'm-');
        end;
        hold off
%         saveas(gca, [fName(1:end-16) '_fibs.fig'])
        'saving...'
        tic
        saveas(gca, [fName(1:end-16) '_fibs.png'])        
        toc
        allFibs{f} = fibs;
    end;
    cd ..
    save([dList(d).name '_allImgFibs.mat'], 'pixSize', 'ySize', 'xSize', 'allFibs', 'fibW', 'minL', 'maxL', 'dW')
    close all
end;


load 'int_3h_no_mitomycin_allImgFibs.mat';
fibS = [];
for j = 1:length(allFibs);
    fibS = [fibS, allFibs{j}];
end;

nFib = length(fibS);
for j = 1:10:nFib;
    multicolor(fibS(j).img, colMat);
    pause
end;
            
allL = [fibS(:).L]*pixSize(1);
hist(allL, 30)


        
        
        profile{j} = tmp;
        subplot(1,3,3);
        plot((1:length(profile{j}))*pixSize, profile{j});
        xlabel('Position (\mum)')
        title('Line profile')
        pbaspect([1 1 1])
        
        
    end;
        
end;

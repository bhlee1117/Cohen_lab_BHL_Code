function [avgImgRed, GFP, traces, cellimgs, filters, rawTrace] = ...
    unmix_bin_optopatch_movies(baseDir,nanz,ncells,nr,nc,cycleLen)


cd(baseDir);
dt = 2;
denoise = 1;

number_of_cells = ncells;

addpath('X:\Lab\Computer Code\Image Processing\Hamamatsu DCIMG reader\')
addpath('X:\Lab\Computer Code\Image Processing\FastICA_25\')

dname = dir;
danz = dname(nanz+2).name;
cd(danz);
prefix = [baseDir filesep 'Results' filesep danz(1:2)];

binName = dir('*.bin');
[movtmp nframe] = readBinMov(binName.name, nr, nc);
movtmp = double(movtmp);
ntmp = length(movtmp);
tifName = dir('*.tif');
GFP = imread(tifName.name);

tmpIntensR = squeeze(mean(mean(movtmp,1),2));
badFrames = find(tmpIntensR > mean(tmpIntensR) + 2*std(tmpIntensR));
ims = movtmp;
ims(:,:,[badFrames badFrames+1 badFrames+2]) =[];
nt = length(ims);
avgImgRed = mean(ims,3);
rawTrace = squeeze(mean(mean(ims,2),1));

if denoise % spatial 2D median filtering of source images
    for it = 1:size(ims,3)
        ims(:,:,it) = medfilt2(ims(:,:,it),[3 3], 'symmetric');
    end
end
tovec = @(m)reshape(m,nr*nc,[]);
toimg = @(m)reshape(m,nr,nc,[]);
ims = tovec(ims);
ims = double(ims);

% high-pass filter the movie so the ICA is done only on spikes
f = fspecial('gaussian', [1 20/dt], 3/dt);  % The parameters of this filter are set for a 500 Hz frame rate
f = f - mean(f);
dims = imfilter(ims, f, 'replicate');   %high-pass filtered version of ims

%Perform time-domain pca
number_of_pcs = 20;
tic
[pcs,eigvals] = eigs(dims'*dims,number_of_pcs);
disp(['pca took ' num2str(toc) ' s']);

% Perform time-domain ICA
number_of_ics = 10;
tic
[ics0,mixmat0,sepmat0] = fastica(pcs',... % ica
    'numofic',number_of_ics,...
    'approach','symm',...
    'verbose','off');
ics0 = ics0';
icskews = skewness(ics0);
[~,iidx] = sort(abs(icskews),'descend'); % sort ic skews, get indexing values
ics = bsxfun(@times,ics0(:,iidx),sign(icskews(iidx))); % sort ics and invert negative skews
mixmat = bsxfun(@times,mixmat0(:,iidx),sign(icskews(iidx))); % mixmat
sepmat = bsxfun(@times,sepmat0(iidx,:),sign(icskews(iidx)')); % sepmat
disp(['ica took ' num2str(toc) ' s']);
clear ics0 mixmat0 sepmat0

%
cellimgs = toimg(dims*ics);
filters = toimg((dims*pcs)*diag(1./diag(eigvals))*sepmat'); % no offset in signals

if denoise  % median filter the filters
    for j = 1:number_of_ics;
        cellimgs(:,:,j) = medfilt2(cellimgs(:,:,j), [3 3], 'symmetric'); 
        filters(:,:,j) = medfilt2(filters(:,:,j), [3 3], 'symmetric');  
    end;
end;
    
traces = ims'*tovec(filters(:,:,1:number_of_cells));

% %% Sort by a combination of skewness in time and skewness in space
% [~,cidx] = sort(skewness(cellimgs),'descend'); % sort cells by skew, get indexing values
% scores = 1./(1:size(ics,2)); % decreasing score for previous ics order
% [~,idx] = sort(scores(cidx)./(1:size(ics,2)),'descend'); % sort the product of score in space and score in time
% 
% % sort ic dimension in ic variables to be used further on
% % ics      = ics     (:,cidx(idx));
% % mixmat   = mixmat  (:,cidx(idx));
% cellimgs = cellimgs(:,cidx(idx));
% filters  = filters (:,cidx(idx));
% % sepmat   = sepmat  (cidx(idx),:);

%
avgimg = mean(ims,2);
stdimg = std(ims, [], 2);


% format output for nr,nc agnostic environment outside
avgimg   = toimg(avgimg  );
stdimg   = toimg(stdimg  );




nics = number_of_cells; % number of independent components shown
g = 1/nics/40; % gap width
w = .08; % width of each image subplot
traces = bsxfun(@minus, traces, mean(traces));
scaling = 1.2/(max(traces(:)) - min(traces(:)));

% standardize window level scale
sc = sort(cellimgs(:));
sclims = sc(ceil(end*[.001 .9999]));
sf = sort(filters(:));
sflims = sf(ceil(end*[.001 .9999]));
dscale = mean(sclims./sflims);

%%
if ~exist(prefix(1:end-2),'file')
    mkdir(prefix(1:end-2))
end
if ~exist('traces','var')
    return
end
save([prefix danz(3:end) '_results.mat'],...
    'traces','cellimgs','filters','avgimg','stdimg','GFP','rawTrace')
if ~size(traces,1)
    return
end

figure(1)
clf
set(gcf,'position',[100 200 1366 768])
for it = 1:nics
    subplot('position',[g (nics-it)/(nics+1)+g w-g 1/(nics+1)-g])
    imshow(cellimgs(:,:,it)/dscale,sflims)
    subplot('position',[w+g (nics-it)/(nics+1)+g w-g 1/(nics+1)-g])
    imshow(filters(:,:,it),sflims)
end
subplot('position',[g nics/(nics+1)+g w-g 1/(nics+1)-g])
imshow(avgimg,[])
subplot('position',[w+g nics/(nics+1)+g w-g 1/(nics+1)-g])
imshow(stdimg,[])
subplot('position',[2*w+g g 1-2*w-2*g 1-2*g])
plot(bsxfun(@plus,nics-1:-1:0,traces(:,1:nics)*scaling));
tl = mat2cell(char(33*ones(1,nics)),1,ones(1,nics));
for itl=1:numel(tl), tl{itl} = ''; end
set(gca,'YTick',[0:nics-1],'XTick',[],'Yticklabel',tl)
axis tight
ylim([min(-.5,min(traces(:,1)*scaling)) max(nics+.5,max(traces(:,1)*scaling)+nics-1)])
drawnow
saveas(gca, [prefix danz(3:end) '_traces.fig'])
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 1366 768]);
print('-dpng','-r1',[prefix danz(3:end) '_traces.png'])
close all


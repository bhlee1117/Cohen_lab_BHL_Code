clear
close all


%% Select directories to analyze
rootPath = 'C:\Joel\Data\';
expdirs = {'2014-07-15 3d Round 1 (a)',...
    };

%% Locate all spikes in all traces:: ONLY RUN ONCE
c = 1; % cell counter
onespike = fspecial('gaussian', [10, 1], .6); % low pass filter
onespike = onespike - mean(onespike);
clear data;
figure(1)
for ii = 1:numel(expdirs)
    exppath = [expdirs{ii} filesep];
    filedirs = dir([rootPath exppath 'Take*']);
    for jj = 1:numel(filedirs)
        filepath = [filedirs(jj).name filesep];
        fnames = dir([rootPath exppath filepath 'Results\*.mat']);
        for kk = 1:numel(fnames)
            load([rootPath exppath filepath 'Results\' fnames(kk).name]);
            [nframes ntraces] = size(traces);
            tracesf = imfilter(traces, onespike, 'replicate')/sum(onespike.^2);
            trtmp = sort(tracesf,1);
            t75 = trtmp(round(.75*end),:);
            t50 = median(tracesf);
            thresh = t50 + 7*(t75 - t50);
            xitions = 1:525:nframes;  % these are the intervals between recordings
            traceN = bsxfun(@plus, mat2gray(bsxfun(@minus, traces, mean(traces))), 1:ntraces);
            plot(traceN);
            for m = 1:ntraces;
                spikeT = spikefind(tracesf(:,m), thresh(m));
%                 plot(tracesf(:,m)); hold all; plot([t(1) t(end)], [thresh(m) thresh(m)]); hold off
                spikeT = setdiff(spikeT, [xitions, xitions + 1, xitions - 1]);
                 if (~isempty(spikeT) || m == 1)
                    hold all;
                    plot(spikeT, traceN(spikeT,m), 'r*')
                    hold off
                    drawnow
                    data(c).trace = traces(:,m);
                    data(c).filter = filters(:,:,m);
                    data(c).cellimg = cellimgs(:,:,m);
                    data(c).avgimg = avgimg;
                    data(c).round = expdirs{ii};
                    data(c).take = filedirs(jj).name;
                    data(c).fov = fnames(kk).name(1:2);
                    data(c).gfp = GFP;
                    data(c).spikeT = spikeT;
                    c = c + 1
                    if isempty(spikeT);  % Break if there are no spikes in the first component
                        break
                    end;
                 else;  % as soon as one trace has no spikes, don't look to higher order traces
                     break;
                 end;
            end;
                        
        end
    end
end

save allcells.mat data



%% Remove bad traces:: ONLY RUN ONCE
for j = 1:length(data)
    data(j).totalspikes = length(data(j).spikeT);
    totalspikes(j) = length(data(j).spikeT);
end
ncells = length(data);
figure(1)
plot(totalspikes)

% Remove if there's too many or too few spikes
data(totalspikes < 5 | totalspikes > 160) = [];

clear totalspikes
for j = 1:length(data)
    data(j).totalspikes = length(data(j).spikeT);
    totalspikes(j) = length(data(j).spikeT);
end
ncells = length(data);
figure(1)
plot(totalspikes)


% Go through by hand to select real traces
figure(1)
for j = 1:ncells
    subplot(1,4,1)
    imshow(data(j).cellimg,[])
    subplot(1,4,2:4)
    plot(data(j).trace)
    title([num2str(j) ' of ' num2str(ncells)])
    output = input('Is this a real trace: ','s');
    if output == 'y'
        badindx(j) = 0;
    else
        badindx(j) = 1;
    end
end

clear totalspikes
badtraces = find(badindx);
data(badtraces) = [];

ncells = length(data);
for j = 1:length(data)
    data(j).totalspikes = length(data(j).spikeT);
    totalspikes(j) = length(data(j).spikeT);
end
figure(1)
plot(totalspikes)

save allcells
%% Start to analyze traces

load(['allcells.mat'])

% Count the frames.  The first 46 cells are shorter.
for j = 1:ncells;
    nframes(j) = length(data(j).trace);
end;
plot(nframes)
t = (1:max(nframes))*2;

% Associate the data to a directory
clear dlist
for j = 1:ncells;
    dlist{j} = data(j).take;
end;
[dlist, dindx, cindx] = unique(dlist, 'stable');  % Find the unique directory entries.
plot(cindx)
dindx2 = [(dindx(2:end)-1); ncells];  % Indices of the end of each take
nruns = length(dindx);

AA = logical([1 0 0 0]); 
allAA = AA(cindx);
BB = logical([0 1 0 0]);
allBB = BB(cindx);
CC = logical([0 0 1 0]);
allCC = CC(cindx);
DD = logical([0 0 0 1]);
allDD = DD(cindx);

% Add a field which keeps track of cell # within each fov;
roinum = zeros(ncells, 1);
roinum(1) = 1;
for j = 2:ncells;
    if strcmp(data(j).fov, data(j-1).fov)
        roinum(j) = roinum(j-1) + 1;
    else
        roinum(j) = 1;
    end;
end;
mean(roinum(allAA)) % =  1.265
mean(roinum(allBB))  % = 1.819
mean(roinum(allCC))
mean(roinum(allDD))

%% Plot out cell images
close all
for j = 1:length(AA)
    figure('units','normalized','position',[.1 .1 .25 .4])
    q = 1;
    for k = dindx(j):dindx2(j);
        if roinum(k) == 1;
            [r, c] = ind2sub([4, 4], q);
            subplot('position',[(r-1)/4 (c-1)/4 .24 .24])
            imshow(double(data(k).gfp).*data(k).cellimg, [])
            q = q + 1;
            if q > 16; break; end;
        end;
    end;
    if AA(j); text(10, 20, '39bCorr', 'color', 'white'); 
    elseif BB(j); text(10, 20, '39b', 'color', 'white');  
    elseif CC(j); text(10, 20, 'RB9d', 'color', 'white');
    elseif DD(j); text(10, 20, 'RB9dCorr', 'color', 'white'); end;
    tmp = char(dlist(j));
    saveas(gca,[tmp(1:6) '_cell_images.png'])
    saveas(gca,[tmp(1:6) '_cell_images.fig'])
end;


%% Plot out individual traces
nplot = 20;

for d = 1:length(AA)
    figure
    if dindx2(d)-dindx(d) < nplot
        nplot2 = dindx2(d)-dindx(d)-1;
    else
        nplot2 = nplot;
    end
    for j = 1:nplot2
        trace = data(dindx(d)+j).trace;
        time = (1/500)*(1:length(trace));
        plot(time,mat2gray(trace)+j)
        hold all
    end
    xlabel('Time (s)')
    ylabel('Normalized Fluorescence')
    if AA(d); title('39bCorr');
    elseif BB(d); title('39b');  
    elseif CC(d); title('RB9d')
    elseif DD(d); title('RB9dCorr')
    end;
    tmp = char(dlist(d));
    saveas(gca,[tmp(1:6) '_raw_traces.png'])
    saveas(gca,[tmp(1:6) '_raw_traces.fig'])
end


%% Remove photobleaching of each trace
rawtraces = zeros(5220, ncells);
for j = 1:ncells;
    rawtraces(1:nframes(j),j) = data(j).trace;
end;
% remove photobleaching by a linear fit to each illumination cycle.
alltraces = zeros(size(rawtraces));
xitions = 1:522:5220;
modelfun = [ones(50,1), [1:25  (-24:0)+522]'];  % Fit to a constant plus a linear slope.
basefun = [ones(522,1) (1:522)'];
for j = 1:8;  % Do the fit for all traces, in chunks corresponding to each illumination cycle
    subdat = rawtraces([xitions(j):(xitions(j) + 24) (xitions(j+1)-25):(xitions(j+1)-1)],:);
    fitfuns = modelfun\subdat;
    alltraces(xitions(j):(xitions(j+1)-1),:) = rawtraces(xitions(j):(xitions(j+1)-1),:) - basefun*fitfuns;
end;
    
% pcolor([alltraces(:,allHeads) alltraces(:,~allHeads)]);  shading 'interp'; caxis([0, 45]);
% colormap('jet')
% plot(alltraces(:,20))

plot(t, mean(alltraces(:,allAA), 2), 'r-', t, mean(alltraces(:,~allAA), 2), 'k-')



%% Find the spikes
clear allspikes
for j = 1:ncells;
    allspikes{j} = data(j).spikeT;
end;

% Fill a matrix with all the spikes
spikemat = zeros(size(rawtraces));
dspikemat = zeros(size(rawtraces));
for j = 1:ncells;
    spikemat(allspikes{j},j) = 1;
    if length(allspikes{j}) >= 2;
        dspikemat(allspikes{j}(2:end),j) = 500./diff(allspikes{j});
    end;
end;
pcolor(spikemat); shading 'interp'
pcolor([spikemat(:,allAA) spikemat(:,~allAA)]);  
set(gca,'ZLim',[0 .1])

sum(allAA)   
sum(~allAA)  
sum(spikemat(:))  

% Calculate total firing density
firingdensity = mean(spikemat, 2);
plot(firingdensity)

%%
% Firing density of heads and tails cells
dirstart = 1; dirstop = 4;
AAdensity = mean(spikemat(:,allAA & (cindx >= dirstart)' & (cindx <= dirstop)'), 2);  %  cindx >= 11 is rounds 2 and 3.
nAA = sum(allAA & (cindx >=dirstart)' & (cindx <= dirstop)')  % = 97
BBdensity = mean(spikemat(:,allBB & (cindx >=dirstart)' & (cindx <= dirstop)'), 2);
nBB = sum(allBB & (cindx >=dirstart)' & (cindx <= dirstop)')  % = 493
CCdensity = mean(spikemat(:,allCC & (cindx >=dirstart)' & (cindx <= dirstop)'), 2);
nCC = sum(allCC & (cindx >=dirstart)' & (cindx <= dirstop)')  % = 493
DDdensity = mean(spikemat(:,allDD & (cindx >=dirstart)' & (cindx <= dirstop)'), 2);
nDD = sum(allDD & (cindx >=dirstart)' & (cindx <= dirstop)')  % = 493

nsmooth = 30;
figure
plot(t, smooth(AAdensity, nsmooth), 'b-', t, smooth(BBdensity, nsmooth), 'g-',...
    t, smooth(CCdensity, nsmooth), 'k-', t, smooth(DDdensity, nsmooth), 'r-');

xlim([0 max(t)])
legend({['n = ' num2str(nAA) ' 39bCorr'],['n = ' num2str(nBB) ' 39b'],...
    ['n = ' num2str(nCC) ' RB9d'],['n = ' num2str(nDD) ' RB9dCorr']}, 'location', 'northwest')
xlabel('Time (ms)')
ylabel('Mean probability of spike')

saveas(gca, 'Biogen lines spiking patterns.fig')
saveas(gca, 'Biogen lines spiking patterns.png')

figure
plot(t, smooth(BBdensity, nsmooth), 'r-', t, smooth(CCdensity, nsmooth), 'g-');
xlim([0 max(t)])
legend({['n = ' num2str(nBB) ' 39b ALS'],['n = ' num2str(nCC) ' 39b retig']}, 'location', 'northwest')
xlabel('Time (ms)')
ylabel('Mean probability of spike')

saveas(gca, 'ALS_v_retig spiking patterns.fig')
saveas(gca, 'ALS_v_retig spiking patterns.png')

figure
plot(t, smooth(AAdensity, nsmooth), 'k-', t, smooth(CCdensity, nsmooth), 'g-');
xlim([0 max(t)])
legend({['n = ' num2str(nHeads) ' 39b WT'],['n = ' num2str(nHeads) ' 39b retig']}, 'location', 'northwest')
xlabel('Time (ms)')
ylabel('Mean probability of spike')

saveas(gca, 'ALS_v_WT spiking patterns.fig')
saveas(gca, 'ALS_v_WT spiking patterns.png')

% Spontaneous firing
spontAAsprob = sum(AAdensity(1:4950))/9.9;
spontBBsprob = sum(BBdensity(1:4950))/9.9;
spontCCprob = sum(CCdensity(1:4950))/9.9;
spontDDprob = sum(DDdensity(1:4950))/9.9;
figure
bar([1;2;3;4],[spontAAsprob;spontBBsprob;spontCCprob;spontDDprob])
ylabel('Firing rate per cell (Hz)')
set(gca,'XTickLabel',[' 10 '; '1872';' 73 ';' 59 '])
saveas(gca,'Heads_v_Tails_spont_firing_rate.png')
saveas(gca,'Heads_v_Tails_spont_firing_rate.fig')

save allcells
%%
% Firing density resolved dish by dish
firingdensities = zeros(max(nframes), nruns);
c = zeros(nruns, 1);
for j = 1:nruns;
    for k = dindx(j):dindx2(j);
        firingdensities(:,j) = firingdensities(:,j) + spikemat(:,k);
        c(j) = c(j) + 1;
    end
    firingdensities(:,j) = firingdensities(:,j)/c(j);
end;

for j = 1:nruns;
    plot(smooth(firingdensities(:,j), 30))
    hold all;
end;
hold off

clf
figure
colororder = {'r-'; 'k-'};
for j = [3:4];
    plot(t, smooth(firingdensities(:,j), 30), colororder{AA(j) + 1});
    hold all
end;
hold off;
legend({'Heads', 'Tails'})
xlabel('Time (ms)')
ylabel('Mean probability of spike')
ylim([0 0.055]);
saveas(gca, 'Dish_by_dish_spike_probabilities.fig')
saveas(gca, 'Dish_by_dish_spike_probabilities.png')


%% Plot out spikes for each excitation
blueOnL = 255; blueOffL = 271;
for j = 1:8
    blueOn(1:blueOnL,j) = (j-1)*522+5021:(j-1)*522 + 5275;
    blueOff(1:blueOffL,j) = (j-1)*523 + blueOnL+1:(j-1)*523 + blueOnL + blueOffL;
end

nback = 10;
nfront = 18;
nmat = nback+nfront+1;

% Average of all the spikes in an epoch
indxAA(1:8) = 0;
indxBB(1:8) = 0;
indxCC(1:8) = 0;
indxDD(1:8) = 0;
for cel = 1:ncells
    spikeT = data(cel).spikeT;
    trace = mat2gray(data(cel).trace);
    for j = 1:8
        tmpspikes = find(spikeT > blueOn(1,j) & spikeT < blueOn(end,j));
        if ~isempty(tmpspikes)
            for s = 1:length(tmpspikes);
                if allAA(cel)
                    indxAA(j) = indxAA(j)+1;
                    spikematAA(1:nmat,j,indxAA(j)) = trace(spikeT(tmpspikes(s))-nback:spikeT(tmpspikes(s))+nfront);
                elseif allBB(cel)
                    indxBB(j) = indxBB(j)+1;
                    spikematBB(1:nmat,j,indxBB(j)) = trace(spikeT(tmpspikes(s))-nback:spikeT(tmpspikes(s))+nfront);
                elseif allCC(cel)
                    indxAA(j) = indxAA(j)+1;
                    spikematCC(1:nmat,j,indxAA(j)) = trace(spikeT(tmpspikes(s))-nback:spikeT(tmpspikes(s))+nfront);
                elseif allDD(cel)
                    indxBB(j) = indxBB(j)+1;
                    spikematDD(1:nmat,j,indxBB(j)) = trace(spikeT(tmpspikes(s))-nback:spikeT(tmpspikes(s))+nfront);
                end
            end
        end
    end
end

kernelsAA = mean(spikematAA,3);
kernelsBB = mean(spikematBB,3);
kernelsCC = mean(spikematCC,3);
kernelsDD = mean(spikematDD,3);


figure
for j = 1:8
    subplot(2,4,j)
    plot((1/500)*(1:nmat),mat2gray(kernelsAA(:,j)),'b',...
        (1/500)*(1:nmat),mat2gray(kernelsBB(:,j)),'g',...
        (1/500)*(1:nmat),mat2gray(kernelsCC(:,j)),'k',...
        (1/500)*(1:nmat),mat2gray(kernelsDD(:,j)),'r')
    title([num2str(j)])
end

saveas(gca,'Average_Excited_Heads_vs_Tails_AP.png')
saveas(gca,'Average_Excited_Heads_vs_Tails_AP.fig')



% Average of first the spike in an epoch
spikeOfInterest = 1;
indxH1(1:8) = 0;
indxT1(1:8) = 0;
for cel = 1:ncells
    spikeT = data(cel).spikeT;
    trace = mat2gray(data(cel).trace);
    for j = 1:8
        tmpspikes = find(spikeT > blueOn(1,j) & spikeT < blueOn(end,j),spikeOfInterest,'last');
        if ~isempty(tmpspikes)
            for s = 1:length(tmpspikes);
                if allAA(cel)
                    indxH1(j) = indxH1(j)+1;
                    spikematH1(1:nmat,j,indxH1(j)) = trace(spikeT(tmpspikes(s))-nback:spikeT(tmpspikes(s))+nfront);
                else
                    indxT1(j) = indxT1(j)+1;
                    spikematT1(1:nmat,j,indxT1(j)) = trace(spikeT(tmpspikes(s))-nback:spikeT(tmpspikes(s))+nfront);
                end
            end
        end
    end
end

kernelsHeads1 = mean(spikematH1,3);
kernelsTails1 = mean(spikematT1,3);


figure
for j = 1:8
    subplot(2,5,j)
    plot((1/500)*(1:nmat),mat2gray(kernelsHeads1(:,j)),'k',...
        (1/500)*(1:nmat),mat2gray(kernelsTails1(:,j)),'r')
    title([num2str(j)])
end

saveas(gca,'First_Excited_Heads_vs_Tails_AP.png')
saveas(gca,'First_Excited_Heads_vs_Tails_AP.fig')



%% Look at average rise of each epoch
nmat2 = 60;
for cel = 1:ncells
    trace = mat2gray(data(cel).trace);
    for j = 1:8
        if allAA(cel)
            spikeriseAA(1:nmat2,j,cel) = trace(blueOn(1:nmat2,j)-15);
        elseif allBB(cel)
            spikeriseBB(1:nmat2,j,cel) = trace(blueOn(1:nmat2,j)-15);
        elseif allCC(cel)
            spikeriseCC(1:nmat2,j,cel) = trace(blueOn(1:nmat2,j)-15);
        elseif allDD(cel)
            spikeriseDD(1:nmat2,j,cel) = trace(blueOn(1:nmat2,j)-15);
        end
    end
end

spikeriseAAavg = mean(spikeriseAA,3);
spikeriseBBavg = mean(spikeriseBB,3);
spikeriseCCavg = mean(spikeriseCC,3);
spikeriseDDavg = mean(spikeriseDD,3);

figure
for j = 1:8
    subplot(2,4,j)
    plot((1/500)*(1:nmat2),mat2gray(spikeriseAAavg(:,j)),'b',...
        (1/500)*(1:nmat2),mat2gray(spikeriseBBavg(:,j)),'g',...
        (1/500)*(1:nmat2),mat2gray(spikeriseCCavg(:,j)),'k',...
        (1/500)*(1:nmat2),mat2gray(spikeriseDDavg(:,j)),'r')
    title([num2str(j)])
end

saveas(gca,'Blue_light_rise_heads_v_tails.png')
saveas(gca,'Blue_light_rise_heads_v_tails.fig')
















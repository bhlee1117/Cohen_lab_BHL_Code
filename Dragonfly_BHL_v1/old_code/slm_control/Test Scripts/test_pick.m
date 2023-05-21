m = vm('R:\S1\FOV2\stack 1620');
mb = m.pblc.blur(1);
moviefixsc(m)
[rois, intens] = nested_clicky_rects(mb,imtophat(max(mb-min(mb)),ones(11)));
apply_clicky_faster(rois,mb)
%%
[~, cellDepthIdx] = max(intens);
cellBrightnessDiff = max(intens)-min(intens);
cellDepthFit = zeros(1,numel(rois));
cellFWHMFit = zeros(1,numel(rois));
cellBrightnessFit = zeros(1,numel(rois));
% zax = (1:m.frames)';
% zax = (linspace(1.82252,1.685,23)'-1.82252)*1000*1.6;
zax = linspace(1.65238,1.53665,20)';
gaus = @(zax,z0)exp(-4*log(2)*(zax-z0(1,:)).^2./z0(2,:).^2).*z0(3,:)+z0(4,:);
for icell = 1:numel(rois)
    vinit = [zax(cellDepthIdx(icell)) .1 cellBrightnessDiff(icell) min(intens(:,icell))]';
    objfun = @(z0)sum((intens(:,icell) - gaus(zax,z0)).^2);
    vmin = fminunc(objfun,vinit,optimset('display','off'));
    cellDepthFit(icell) = vmin(1);
    cellFWHMFit(icell) = vmin(2);
    cellBrightnessFit(icell) = vmin(3);
    clf
    plot(zax,gaus(zax,vmin))
    hold on
    plot(zax,intens(:,icell))
    drawnow
end

%% given n depths, which cells should be measured in each plane?
figure windowstyle docked
nDepths = 4;
testDepths = linspace(zax(1),zax(end),nDepths+2)';
testDepths = testDepths(2:end-1);
% %%
% testDepths = [115 120]'/10;
%%
disp ' '
maxPower = 170;
unitPower = 10;
practicalBrightness = gaus(testDepths,[
    cellDepthFit
    cellFWHMFit % *0+mean(cellFWHMFit)
    cellBrightnessFit
    cellBrightnessFit*0
]);
relativeIllumNeed = unitPower*max(cellBrightnessFit)./practicalBrightness;
clf
plot(cellDepthFit',relativeIllumNeed','o')
ylim([1 5]*unitPower)
ylabel 'Power (mW)'
xlabel 'Depth (frame)'
hold on
idxAtEachDepth = cell(numel(testDepths),1);
isCellAvailable = relativeIllumNeed < 3*unitPower;
for it = 1:numel(cellDepthFit)
    idxAvailableCells = find(isCellAvailable(:));
    if isempty(idxAvailableCells)
        break
    end
    [~, ix] = min(relativeIllumNeed(isCellAvailable(:)));
    [ii, jj] = ind2sub(size(relativeIllumNeed),idxAvailableCells(ix));
    idxAtEachDepth{ii} = [idxAtEachDepth{ii} jj];
    isCellAvailable(:,jj) = false;
    if sum(relativeIllumNeed(ii,idxAtEachDepth{ii})) > maxPower
        idxAtEachDepth{ii} = idxAtEachDepth{ii}(1:end-1);
        isCellAvailable(ii,:) = false;
    end
end
validDepthsMask = ~cellfun(@(el)isempty(el),idxAtEachDepth);
testDepths = testDepths(validDepthsMask);
relativeIllumNeed = relativeIllumNeed(validDepthsMask,:);
idxAtEachDepth = idxAtEachDepth(validDepthsMask);

% isCellAvailable
% disp(it)
% idxAtEachDepth{:}
for ii = 1:numel(testDepths)
% disp(numel(idxAtEachDepth{ii}))
    plot(cellDepthFit(idxAtEachDepth{ii}),relativeIllumNeed(ii,idxAtEachDepth{ii}),'+')
end
%
newDepths = testDepths;
for ii = 1:numel(testDepths)
    fun = @(z) sum(max(cellBrightnessFit)./gaus(z,[
        cellDepthFit(idxAtEachDepth{ii})
        cellFWHMFit(idxAtEachDepth{ii}) % *0+mean(cellFWHMFit(idxAtEachDepth{ii}))
        cellBrightnessFit(idxAtEachDepth{ii})
        cellBrightnessFit(idxAtEachDepth{ii})*0
    ]));
    newDepths(ii) = fminunc(fun,testDepths(ii),optimset('display','off'));

%     testDepths(ii) = mean(cellDepthFit(idxAtEachDepth{ii}));
%     disp(mean(cellDepthFit(idxAtEachDepth{ii})))
end

for ii = 1:numel(testDepths)
    fprintf('test depth = %6.5g, new depth = %6.5g, nc = %2d, p = %5.4g\n',testDepths(ii),newDepths(ii),numel(idxAtEachDepth{ii}),sum(relativeIllumNeed(ii,idxAtEachDepth{ii})))
end
testDepths = newDepths;
%%
figure windowstyle docked
for ii = 1:numel(testDepths)
    subplot(numel(testDepths),1,ii)
    imshow(interp3(double(m.data),1:m.cols,(1:m.rows)',interp1(linspace(zax(1),zax(end),m.frames),1:m.frames,testDepths(ii))),[])
%     hold on
%     xy = cell2mat(cellfun(@(el)mean(el(1:4,:))',rois(idxAtEachDepth{ii}),'uni',false))';
%     plot(xy(:,1),xy(:,2),'o')
    title(sprintf('depth = %6.5g, nc = %2d, p = %5.4g\n',testDepths(ii),numel(idxAtEachDepth{ii}),sum(relativeIllumNeed(ii,idxAtEachDepth{ii}))))
end



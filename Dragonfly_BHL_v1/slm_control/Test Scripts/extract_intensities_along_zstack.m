mov = vm('R:\S1\FOV1\1633 stack');
nested_clicky_rects(mov,max(vm(imtophat(mov.data,ones(11)))));
rois = ans;
intens = apply_clicky_faster(rois,mov);
% select planes
% manually go to offset Arch imaging
% auto click neurons
%%
nclicked = size(intens,2);
ndepths = 4;
ncperdepth = 10;
correctedIntens = intens - linspace(0,1,mov.frames)'*(intens(end,:)-intens(1,:))+intens(1,:);
normCorrectedIntens = (correctedIntens - mean(correctedIntens))./std(correctedIntens);
varOverMin = sum((correctedIntens - min(correctedIntens)).^2)./sum(intens.^2);
[~, sIdxVOM] = sort(varOverMin);
[~, idxFrameMax] = max(normCorrectedIntens);
[~, sIdxFrameMax] = sort(idxFrameMax);
figure windowstyle docked
plot(normCorrectedIntens(:,sIdxFrameMax) + (1:size(intens,2))/2)
%%
rng(142857)
depthKMIdx = kmeans(idxFrameMax',ndepths,'Replicates',10);
fprintf('\n')
for idepth = 1:ndepths
    maskCellsAtDepth = depthKMIdx == idepth;
    fprintf('depth %d, %d cells, z(frame)=%g, z(GFP)=%g, z(Arch)=%g\n',...
        idepth,...
        nnz(maskCellsAtDepth),...
        mean(idxFrameMax(maskCellsAtDepth)),...
        interp1([1 21],[1.82248 1.71083],mean(idxFrameMax(maskCellsAtDepth))),...
        interp1([1 21],[1.82248 1.71083],mean(idxFrameMax(maskCellsAtDepth)))-.0248)
end
%%
figure windowstyle docked
for idepth = 1:ndepths
    maskCellsAtDepth = depthKMIdx == idepth;
    idxCellsAtDepth = find(maskCellsAtDepth);
    subplot(ndepths,1,idepth)
    imshow(interp3(double(mov.data),1:mov.cols,(1:mov.rows)',mean(idxFrameMax(depthKMIdx == idepth))),[])
    title(idepth)
    hold on
    for icell = 1:nclicked
        xv = rois{icell}(:,1);
        yv = rois{icell}(:,2);
        plot(xv, yv,'k'); % , 'Linewidth', 1,'Color',currcolor);
    end
    for icell = 1:numel(idxCellsAtDepth)
        xv = rois{idxCellsAtDepth(icell)}(:,1);
        yv = rois{idxCellsAtDepth(icell)}(:,2);
        plot(xv, yv,'r'); % , 'Linewidth', 1,'Color',currcolor);
    end      
end
figure windowstyle docked
for idepth = 1:ndepths
    maskCellsAtDepth = depthKMIdx == idepth;
    idxCellsAtDepth = find(maskCellsAtDepth);
    subplot(ndepths,1,idepth)
    imshow(interp3(double(mov.data),1:mov.cols,(1:mov.rows)',mean(idxFrameMax(depthKMIdx == idepth))),[])
    title(idepth)
end
%%
histogram(ix)

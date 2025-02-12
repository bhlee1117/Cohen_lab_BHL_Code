function spatial_info = SpatialInfo(SpikeTr, VRbinTr)
pos_bin=max(VRbinTr);
MeanF=mean(SpikeTr,'omitnan');
FR=zeros(1,pos_bin);
PR=zeros(1,pos_bin);

VRbinTr(isnan(SpikeTr))=NaN;

for k=1:pos_bin
    frameInK=find(VRbinTr==k);
PR(k)=length(frameInK)./sum(~isnan(VRbinTr));
FR(k)=mean(SpikeTr(frameInK),'omitnan');
end

spatial_info = sum(PR.*FR.*log2(FR/MeanF),'omitnan')/MeanF;
% bits per spike

end
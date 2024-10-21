
function [tr_sub binary_dilation]=get_subthreshold(tr,sp,dilate,avgwnd)
tr_sub=[];
for n=1:size(tr,1)
    n
tr_sub(n,:)=tr(n,:);
se = strel('square', dilate); % 0 degree means horizontal
seg_wind=[-30:30];
binary_dilation=ones(1,length(tr));
    for s=find(sp)
        try
    tmp=zeros(1,length(tr));
    tmp(s)=1;
    sp_di = imdilate(tmp, se);
    binary_dilation(find(sp_di))=0;
    %binary_dilate=sp_di;
    tr_sub(n,find(sp_di))=NaN;
    tr_sub_seg=tr_sub(n,s+seg_wind);
    valid_point=find(~isnan(tr_sub_seg));
    tr_sub_seg=interp1(valid_point,tr_sub_seg(valid_point),[1:length(seg_wind)],'linear');
    tr_sub(n,s+seg_wind)=tr_sub_seg;
        end
    end

    tr_sub(n,:)=movmean(tr_sub(n,:),avgwnd,'omitnan');

end
end
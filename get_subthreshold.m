function tr_sub=get_subthreshold(tr,sp,dilate,avgwnd)
tr_sub=tr;
se = strel('square', dilate); % 0 degree means horizontal
seg_wind=[-30:30];
    for s=find(sp)
        try
    tmp=zeros(1,length(tr));
    tmp(s)=1;
    sp_di = imdilate(tmp, se);
    tr_sub(find(sp_di))=NaN;
    tr_sub_seg=tr_sub(s+seg_wind);
    valid_point=find(~isnan(tr_sub_seg));
    tr_sub_seg=interp1(valid_point,tr_sub_seg(valid_point),[1:length(seg_wind)],'linear');
    tr_sub(s+seg_wind)=tr_sub_seg;
        end
    end

    tr_sub=movmean(tr_sub,avgwnd);

end
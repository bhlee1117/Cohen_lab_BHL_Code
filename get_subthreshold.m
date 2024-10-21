function [tr_sub binary_dilation]=get_subthreshold(tr,sp,dilate,avgwnd)
<<<<<<< HEAD

tr_sub=tr;
=======
tr_sub=[];
for n=1:size(tr,1)
    n
tr_sub(n,:)=tr(n,:);
>>>>>>> 19582be0c2e64e1213569743106aa681a92c454e
se = strel('square', dilate); % 0 degree means horizontal
seg_wind=[-30:30];
binary_dilation=zeros(1,size(tr,2));
binary_dilation(find(sp))=1;
binary_dilation=imdilate(binary_dilation,se);

tr_sub(:,binary_dilation==1)=NaN;
valid_point=find(binary_dilation==0);
for n=1:size(tr_sub,1)
tr_sub(n,:)=interp1(valid_point,tr_sub(n,valid_point)',[1:size(tr,2)],'linear');
end
    % for s=find(sp)
    %     try
    % tmp=zeros(1,length(tr));
    % tmp(s)=1;
    % sp_di = imdilate(tmp, se);
    % binary_dilation(find(sp_di))=0;
    % %binary_dilate=sp_di;
    % tr_sub(n,find(sp_di))=NaN;
    % tr_sub_seg=tr_sub(n,s+seg_wind);
    % valid_point=find(~isnan(tr_sub_seg));
    % tr_sub_seg=interp1(valid_point,tr_sub_seg(valid_point),[1:length(seg_wind)],'linear');
    % tr_sub(n,s+seg_wind)=tr_sub_seg;
    %     end
    % end

    tr_sub=movmean(tr_sub,avgwnd,2,'omitnan');

end
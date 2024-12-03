function [tr_sub tr_sub2]=get_subthreshold(tr,sp,dilate,avgwnd)
%input: tr = N (neuron) X T (time) voltage trace,
%       sp = N (neuron) X T (time) binary spike trace
%       dilate = number of peri-spike frames to omit, ex) 7: -3 frame to 3 frame
%       avgwnd = window size to smooth after interpolation

tr_sub=NaN(size(tr));
se = strel('square', dilate);
seg_wind=[-30:30];
binary_dilation=zeros(1,size(tr,2));
binary_dilation(find(sp))=1;
binary_dilation=imdilate(binary_dilation,se);
 
valid_point=find(binary_dilation==0);
for n=1:size(tr_sub,1)
tr_sub(n,:)=interp1(valid_point,tr(n,valid_point)',[1:size(tr,2)],'linear');
%tr_sub(n,:)=interp1(valid_point,tr(n,valid_point)',[1:size(tr,2)],'spline');
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
    tr_sub2=tr_sub;

    tr_sub=movmean(tr_sub,avgwnd,2,'omitnan');

end
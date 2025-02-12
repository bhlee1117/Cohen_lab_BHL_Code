function dirtMov_dilate = tracking_dirt(mov_res,detection_threshold)
load(fullfile('/Volumes/BHL18TB_D2/20240807/201405BHLm141_N1_VR_LowStim','disk_ref'),'z2');
bound=6;
sz=size(mov_res);
mov_int=imresize(mov_res(bound:end-bound,bound:end-bound,:),1/2);
mov_int_filt=imgaussfilt3(mov_int,[8 8 5])-imgaussfilt3(mov_int,[3 3 5]);

match_img=[]; match_data=[];
for k=1:size(mov_int_filt,3)
    [ match_data{k}, match_img(:,:,k)] = matchPattern(mov_int_filt(:,:,k), z2, detection_threshold, 2);
end
match_data=cellfun(@(x) 2*(x+[21 21 0]),match_data,'UniformOutput',false);
detectedPtl=trackParticle(match_data,5,20);
traveldist=cell2mat(cellfun(@(x) sum(sqrt(sum((x(end,1:2)-x(1,1:2)).^2,2)),1),detectedPtl,'UniformOutput',false));
traveltime=cell2mat(cellfun(@(x) size(x,1),detectedPtl,'UniformOutput',false));
valid_ptl=traveldist>50 & traveltime>70;
filteredPtl=detectedPtl(valid_ptl);
dirtMov=zeros(size(mov_res));

ref_im=abs(sqrt(mean(mov_res(:,:,2:round(sz/3)).*mov_res(:,:,1:round(sz/3)-1),3)));

figure(222); clf; filteredPtl_exp=[];
imshow2(ref_im,[]); hold all
for p=1:length(filteredPtl)
    t_exp=[round(filteredPtl{p}(1,3) - range(filteredPtl{p}(:,3))*0.4):1:round(filteredPtl{p}(end,3) + range(filteredPtl{p}(:,3))*0.4)];

    p1 = polyfit(filteredPtl{p}(:,3), filteredPtl{p}(:,1), 1);
    p2 = polyfit(filteredPtl{p}(:,3), filteredPtl{p}(:,2), 1);

    filteredPtl_exp{p}(:,1)=round(p1(1)*t_exp+p1(2));
    filteredPtl_exp{p}(:,2)=round(p2(1)*t_exp+p2(2));
    filteredPtl_exp{p}(:,3)=t_exp;
    outofrngInd=find(filteredPtl_exp{p}(:,1)+bound-1<1 | filteredPtl_exp{p}(:,1)+bound-1>sz(1) | filteredPtl_exp{p}(:,2)+bound-1<1 | filteredPtl_exp{p}(:,2)+bound-1>sz(2) | t_exp' < 1 | t_exp' > size(mov_res,3));
    filteredPtl_exp{p}(outofrngInd,:)=[];
    plot(filteredPtl_exp{p}(:,2)+bound-1,filteredPtl_exp{p}(:,1)+bound-1,'r');
    plot(filteredPtl{p}(:,2)+bound-1,filteredPtl{p}(:,1)+bound-1,'wo');
    for dd=1:size(filteredPtl_exp{p},1)
        dirtMov(filteredPtl_exp{p}(dd,1)+bound-1,filteredPtl_exp{p}(dd,2)+bound-1,filteredPtl_exp{p}(dd,3))=1;
    end
end
drawnow;
se = strel('disk',30);
dirtMov_dilate = imdilate(dirtMov, se);
end
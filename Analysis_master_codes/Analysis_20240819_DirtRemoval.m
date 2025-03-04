clear
clc;
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:Z31');

ref_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,9),'UniformOutput',false);

oblique_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,10),'UniformOutput',false);
PeriSoma_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,11),'UniformOutput',false);
basal_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false);
apical_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,13),'UniformOutput',false);
distal_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,14),'UniformOutput',false);

fpath=raw(:,1)';
StructureData=raw(:,8);
BadROI=cellfun(@(x) (str2num(num2str(x))),raw(:,17),'UniformOutput',false);
EndFrame=cell2mat(raw(:,15));
ifmotionReject=cell2mat(raw(:,16));
ifdirtRemov=cell2mat(raw(:,18));
Pixelsize=cell2mat(raw(:,6));
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;

%%
load(fullfile(fpath{20},'disk_ref'),'z2');

for f=[27]%length(fpath)
    load(fullfile(fpath{f},'PC_Result.mat'),'Result');
    load([fpath{f} '/output_data.mat'])
    sz=double(Device_Data{1, 3}.ROI([2 4])); blueDMDcontour=[];
    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

    bound=10;
    frm_end=EndFrame(f);
    f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
    take_window=repmat([1 time_segment],length(f_seg)-1,1);
    take_window(2:end,1)=take_window(2:end,1)+overlap; take_window(1:end-1,2)=take_window(1:end-1,2)+overlap;
    take_window(end)=mod(f_seg(end),time_segment);
    take_window(take_window==0)=time_segment;
    Result.dirtTrace=[];

    for j=1:length(f_seg)-1
        j
        mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
        load([fpath{f} '/mcTrace' num2str(j,'%02d') '.mat']);

        mov_mc=mov_mc(:,:,[take_window(j,1):take_window(j,2)]);
        mc=mcTrace.xymean([take_window(j,1):take_window(j,2)],:);
        Blue_on_Seg=unique(ceil(find(Result.Blue)/time_segment));

        if ismember(j,Blue_on_Seg)
            isfirst_seg=double(j==1);
            [~, blueomitTr]=get_blueoffTrace(squeeze(mean(mov_mc,[1 2])), ...
                Result.Blue(f_seg(j)+take_window(j,1)-isfirst_seg:f_seg(j)+take_window(j,2)-isfirst_seg),30);
            bkg = zeros(1, size(mov_mc,3));
            bkg(1,:)=movmedian(blueomitTr,4000,'omitnan');
        else
            bkg = zeros(2, size(mov_mc,3));
            bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
            bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
        end

        mov_res=mov_mc;
        mov_res= mov_mc-median(mov_mc,3);
        mov_res = SeeResiduals(mov_res,mc);
        mov_res = SeeResiduals(mov_res,mc.^2);
        mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,end));
        mov_res= SeeResiduals(mov_res,bkg,1);
        if sum(isnan(mov_res(:)))>0
            mov_res= mov_mc-median(mov_mc,3);
        end
        mov_int=imresize(mov_res(bound:end-bound,bound:end-bound,:),1/2);
        mov_int_filt=imgaussfilt3(mov_int,[8 8 5])-imgaussfilt3(mov_int,[3 3 5]);

        match_img=[]; match_data=[];
        for k=1:size(mov_int_filt,3)
            [ match_data{k}, match_img(:,:,k)] = matchPattern(mov_int_filt(:,:,k), z2, 0.3, 2);
        end
        match_data=cellfun(@(x) 2*(x+[21 21 0]),match_data,'UniformOutput',false);
        detectedPtl=trackParticle(match_data,5,20);
        traveldist=cell2mat(cellfun(@(x) sum(sqrt(sum((x(end,1:2)-x(1,1:2)).^2,2)),1),detectedPtl,'UniformOutput',false));
        traveltime=cell2mat(cellfun(@(x) size(x,1),detectedPtl,'UniformOutput',false));
        valid_ptl=traveldist>50 & traveltime>70;
        filteredPtl=detectedPtl(valid_ptl); 
        dirtMov=zeros(size(mov_res));

        figure(2); clf; filteredPtl_exp=[];
        imshow2(Result.ref_im(bound:end-bound,bound:end-bound),[]); hold all
        for p=1:length(filteredPtl)
            t_exp=[round(filteredPtl{p}(1,3) - range(filteredPtl{p}(:,3))*0.4):1:round(filteredPtl{p}(end,3) + range(filteredPtl{p}(:,3))*0.4)];

            p1 = polyfit(filteredPtl{p}(:,3), filteredPtl{p}(:,1), 1);
            p2 = polyfit(filteredPtl{p}(:,3), filteredPtl{p}(:,2), 1);

            filteredPtl_exp{p}(:,1)=round(p1(1)*t_exp+p1(2));
            filteredPtl_exp{p}(:,2)=round(p2(1)*t_exp+p2(2));
            filteredPtl_exp{p}(:,3)=t_exp;
            outofrngInd=find(filteredPtl_exp{p}(:,1)+bound-1<1 | filteredPtl_exp{p}(:,1)+bound-1>sz(2) | filteredPtl_exp{p}(:,2)+bound-1<1 | filteredPtl_exp{p}(:,2)+bound-1>sz(1) | t_exp' < 1 | t_exp' > size(mov_res,3));
            filteredPtl_exp{p}(outofrngInd,:)=[];
            plot(filteredPtl_exp{p}(:,2)+bound-1,filteredPtl_exp{p}(:,1)+bound-1,'r'); 
            plot(filteredPtl{p}(:,2)+bound-1,filteredPtl{p}(:,1)+bound-1,'wo'); 
            for dd=1:size(filteredPtl_exp{p},1)
            dirtMov(filteredPtl_exp{p}(dd,1)+bound-1,filteredPtl_exp{p}(dd,2)+bound-1,filteredPtl_exp{p}(dd,3))=1;
            end
        end
        drawnow
        se = strel('disk',30);
        dirtMov_dilate = imdilate(dirtMov, se);
        Result.dirtTrace=[Result.dirtTrace (tovec(dirtMov_dilate)'*tovec(Result.ftprnt))'];
    end
    save(fullfile(fpath{f},'PC_Result.mat'),'Result','fpath','-v7.3')
    fprintf('Variable "Result" is saved at: %s \n', fpath{f});
end
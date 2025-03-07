
clear
clc;
cd '/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'C5:M164');

save_to='/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,5));
CamType=raw(:,3);
StructureData=raw(:,10);
BadROI=cellfun(@(x) (str2num(num2str(x))),raw(:,13),'UniformOutput',false);
time_segment=25000;
%% Stim

i=102; bound=6;
    load([fpath{i} '/OP_Result.mat'])
    cd(fpath{i});
    load(fullfile(fpath{i},"output_data.mat"))
    sz=double(Device_Data{1, 3}.ROI([2 4]));

    Result.Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
    DMDtrigger=Device_Data{1, 2}.buffered_tasks(1, 2).channels(1, 3).data;
    DMDtrigger=DMDtrigger(CamTrigger);
    Result.DMDtrigger=[0 (DMDtrigger(2:end)-DMDtrigger(1:end-1))>0];
    Result.Blue=Result.Blue(CamTrigger);
    [blueDMDimg bluePatt]=get_blueDMDPatt(Device_Data,'stack');

    DMDbluetrace=(Result.Blue>0).*cumsum(Result.DMDtrigger)+1;
    DMDbluetrace(Result.Blue==0)=0;

    mov_mc=readBinMov_BHL(fpath{i},3);
    load(fullfile(fpath{i},'mcTrace01.mat'))
    mean_F=squeeze(mean(mov_mc(bound:end-bound,bound:end-bound,:),[1 2]));
    [~, blueOff]=get_blueoffTrace(mean_F,[Result.Blue],70);
    [y_fit]=expfitDM_2(find(blueOff)',mean_F(find(blueOff)),[1:size(mov_mc,3)]',1000);

mov_res= mov_mc-mean(mov_mc,3);
mov_res = SeeResiduals(mov_res,y_fit);
mov_res = SeeResiduals(mov_res,mcTrace.xymean(:,:));
mov_res = SeeResiduals(mov_res,mcTrace.xymean(:,:).^2);
mov_res = SeeResiduals(mov_res,mcTrace.xymean(:,1).*mcTrace.xymean(:,end));
dirtMov_dilate = tracking_dirt(mov_res,0.3);
% %%
% [u,s,v] = svds(tovec(mov_res(:,:,1:5000)-mean(mov_res(:,:,1:5000),3)),20);
% reshape_u=reshape(u,sz(2),sz(1),[]);
% bvMask=[];
% [~, bvMask]=get_ROI(max(abs(reshape_u),[],3),bvMask);
% 
% Result.bvMask=bvMask;
% Result.traces_bvMask=(tovec(mov_res.*double(max(Result.bvMask,[],3)==0))'*tovec(Result.ftprnt))';
 Result.dirtTrace=(tovec(dirtMov_dilate)'*tovec(Result.ftprnt))';
%%
VoltageTr=Result.traces_bvMask;
VoltageTr(Result.dirtTrace>0)=NaN;
Result.F0_PCA=get_F0PCA(VoltageTr,3);

for b=1:max(DMDbluetrace)
plot(DMDbluetrace)



%%

for d=1:max(DMDbluetrace)
    BlueOntime=(DMDbluetrace==d);
    Onset_vec=((BlueOntime(2:end)-BlueOntime(1:end-1))==1);

end

%%
ref_imSTA=max(cat(3,SpTAmovie{1},SpTAmovie{2}),[],3);
    ref_imSTA_gauss=imgaussfilt(ref_imSTA,6);
    ref_imSTA=mat2gray(ref_imSTA-ref_imSTA_gauss);
    level = graythresh(ref_imSTA);
    ref_imSTA_bin=ref_imSTA>level*0.8;

    se = strel('sphere', 1);
    ref_imSTA_bin = imdilate(ref_imSTA_bin, se);
    bwSeg=bwlabeln(ref_imSTA_bin);
    segments = regionprops3(bwSeg,'Volume','EquivDiameter');
    segments = table2array(segments);
    %bwlist=find(arrayfun(@(x) x.Volume>4000, segments) & arrayfun(@(x) x.EquivDiameter>10, segments));
    bwlist = find(segments(:,1)>20 & segments(:,2)>1);

    se = strel('sphere',1);
    dendrite_bin=double(ismember(bwSeg,bwlist));
    dendrite_bin= imdilate(dendrite_bin,se);
    dendrite_bin= imgaussfilt3(dendrite_bin,1)>0.1;

    imshow2(imfuse(ref_imSTA,dendrite_bin),[])

    STAmovie=cellfun(@(x) x.*dendrite_bin,STAmovie,'UniformOutput',false);    
%%

F_refImg_filt=imgaussfilt(F_refImg,4);
AvgImg_filt=imgaussfilt(mean(mov_mc,3),4);
[~, kymoROI]=polyLineKymo3(mean(mov_mc,3),5,10);
PCvoltage=[]; STAvoltage=[];

kymoCoord=cell2mat(cellfun(@(x) mean(x([1:end 1],:),1),kymoROI,'UniformOutput',false)');
PattCoord=cell2mat(cellfun(@(x) mean(x([1:end 1],:),1),bluePatt,'UniformOutput',false)');
[~, l]=min(distance_mat(PattCoord(:,[2 1]),kymoCoord),[],2);
kymolength=[0; cumsum(sqrt(sum((kymoCoord(2:end,:)-kymoCoord(1:end-1,:)).^2,2)))]*1.17;

figure(4); clf; 
tiledlayout(max(DMDbluetrace)+1,6);
for p=1:max(DMDbluetrace)
    nexttile([1 3])
    STAimg=max(STAmovie{p},[],3)./AvgImg_filt;
    imshow2(imgaussfilt(STAimg,2),[2 10]*0.01); hold all
    plot(bluePatt{p}(:,2),bluePatt{p}(:,1),'color',[1 1 1],'LineWidth',1.5)
    plot(kymoCoord(:,1),kymoCoord(:,2),'r')
    title('Max projection of Stim. triggered average movie')

    [V D eigTrace]=get_eigvector(tovec(STAmovie{p}));
    [icsTrace, ~, sepmat]=sorted_ica(V(:,1:5),5);
    icsImg=toimg(tovec(STAmovie{p})*icsTrace,size(STAmovie{p},1),size(STAmovie{p},2));

    nexttile([1 3])
    % imshow2(max(pcafilt(STAmovie{p},3),[],3)./F_refImg_filt,[1 3]); hold all
    % plot(bluePatt{p}(:,2),bluePatt{p}(:,1),'r')
    PCimg=toimg(tovec(STAmovie{p}./AvgImg_filt)*movmean(V(:,1),5),size(STAmovie{p},1),size(STAmovie{p},2));
    imshow2(imgaussfilt(PCimg,2),[2 30]*0.01); hold all; colormap(turbo)
    plot(bluePatt{p}(:,2),bluePatt{p}(:,1),'color',[1 1 1],'LineWidth',1.5)
    plot(kymoCoord(:,1),kymoCoord(:,2),'r')
    title('PC#1 of Stim. triggered average movie')

    PCvoltage=[PCvoltage; apply_clicky(kymoROI,imgaussfilt(PCimg,2),'no')];
    STAvoltage=[STAvoltage; apply_clicky(kymoROI,imgaussfilt(STAimg,2),'no')];

end

for p=1:2
nexttile([1 2])
plot(kymolength,STAvoltage(p,:)); hold all; yyaxis right
plot(kymolength,PCvoltage(p,:));
plot([kymolength(l(p)) kymolength(l(p))],[0 0.1],'r')
legend({'Max Proj.','PC #1'})

expDecayModel = @(b, x) b(1) * exp(-x/b(2)) + b(3); % b(1): amplitude, b(2): growth rate, b(3): offset
lb = [0, 30, 0]; % Lower limit: amplitude, decay rate, and offset must be non-negative
ub = [Inf, Inf, Inf]; % No upper limits
options = optimoptions('lsqcurvefit', 'Display', 'off');
b_estimated = lsqcurvefit(expDecayModel, [0.1, 100, 0.2], kymolength(l(p))-kymolength([l(p):-1:1]), STAvoltage(p,[l(p):-1:1])', lb, ub, options);
Ld_bw(p,1)=b_estimated(2);
b_estimated = lsqcurvefit(expDecayModel, [0.1, 100, 0.2], kymolength([l(p):end])-kymolength(l(p)), STAvoltage(p,[l(p):end])', lb, ub, options);
Ld_fw(p,1)=b_estimated(2);
b_estimated = lsqcurvefit(expDecayModel, [0.1, 100, 0.2], kymolength(l(p))-kymolength([l(p):-1:1]), PCvoltage(p,[l(p):-1:1])', lb, ub, options);
Ld_bw(p,2)=b_estimated(2);
b_estimated = lsqcurvefit(expDecayModel, [0.1, 100, 0.2], kymolength([l(p):end])-kymolength(l(p)), PCvoltage(p,[l(p):end])', lb, ub, options);
Ld_fw(p,2)=b_estimated(2);
end

nexttile([1 2])
scatter(tovec(repmat([1:4],2,1)'),tovec([Ld_fw Ld_bw]'),40,repmat(distinguishable_colors(4),2,1),'filled')
xlim([0.5 4.5])
set(gca,'XTick',[1:4],'XTickLabel',{'Outbound max.','Outbound PC#1','Inbound max.','Inbound PC#1'})
ylabel('Length constant (\mum)')

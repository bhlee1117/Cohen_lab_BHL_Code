clear
clc;
cd '/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'B5:K154');

save_to='/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,5));
CamType=raw(:,3);
StructureData=raw(:,10);

%% Low stim
for i=139
    load([fpath{i} '/Result.mat'])
    load(fullfile(fpath{i},"output_data.mat"))
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    mov_mc=double(readBinMov([fpath{i} '/mc_ShutterReg01.bin'],sz(2),sz(1)));
end
[blueDMDimg bluePatt]=get_blueDMDPatt(Device_Data);

%%
bound=5;
mean_F=squeeze(mean(mov_mc(bound:end-bound,bound:end-bound,:),[1 2]));
[~, ~, blueOff]=get_blueoffTrace(mean_F,[Result.Blue],700,200);
[y_fit]=expfitDM_2(find(blueOff)',mean_F(find(blueOff)),[1:size(mov_mc,3)]',[1000 100]);
bkg(1,:)=y_fit;
figure(6); clf;
plot(bkg)
hold all
plot(mean_F)
mov_res=mov_mc-mean(mov_mc,3);
mov_res = SeeResiduals(mov_res,Result.mc);
mov_res = SeeResiduals(mov_res,Result.mc.^2);
mov_res = SeeResiduals(mov_res,Result.mc(:,1).*Result.mc(:,end));
mov_res= SeeResiduals(mov_res,bkg,1);


%% Regress out
regout_ind=[3 4 5 6 7 8];
mov_res=SeeResiduals(mov_res,V_ics(:,regout_ind));

mov_res_mask=mov_res(:,:,:).*double(max(Result.bvMask,[],3)==0);
subMov=tovec(imresize(imgaussfilt3(mov_res_mask(bound:end-bound,bound:end-bound,:),[1 1 0.1]),1/2));
subMov=subMov-mean(subMov,2);
covMat=subMov'*subMov;

[V, D] = eig(covMat);
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V = V.*vSign;

nPCs=28;
eigImg=toimg(tovec(mov_res(:,:,:))*V(:,1:nPCs),size(mov_res,1),size(mov_res,2));
figure(4); clf; 
for n=1:nPCs
    nexttile([1 1])
    imshow2(eigImg(:,:,n),[])
    title(['PC #', num2str(n), ' Fraction : ' num2str(D(n)/sum(D),2)])
end

[V_ics, mixmat, sepmat]=sorted_ica(V(:,1:nPCs),10);
icsImg=toimg(tovec(mov_res(:,:,:))*V_ics,size(mov_res,1),size(mov_res,2));
figure(5); clf; 
for n=1:size(V_ics,2)
    nexttile([1 1])
    %show_footprnt_contour(Result.bvMask,icsImg(:,:,n))
    imshow2(icsImg(:,:,n),[])
    title(['ICS #', num2str(n)])
end
colormap('gray')
figure(7); clf;
plot(rescale2(V_ics,1)+[1:size(V_ics,2)])

%%
movVec=tovec(mov_res);
ROI_ind=[1 2]; %neuron #1 and neuron #2;
ics_ind={[2],1}; %neuron #1 and neuron #2;
ftprnt = zeros(sz(2)*sz(1),2);
for r=1:2
    ROI_soma=ROI_ind(r);
    mask = poly2mask(Result.ROIpoly{ROI_soma}(:,1), Result.ROIpoly{ROI_soma}(:,2), sz(2), sz(1));
    pixelList=find(tovec(squeeze(mask)));
    subMov = movVec(pixelList,:);
    icsTrace=subMov*V_ics;
    coeff=subMov*mean(icsTrace(:,ics_ind{r})*V_ics(:,ics_ind{r})',1)';
    ftprnt(pixelList,r)=coeff;
    %V_trace(:,r)=mean(icsTrace(:,ics_ind{r})*V_ics(:,ics_ind{r})',1)';
    Vics_trpace(:,r)=mean(V_ics(:,ics_ind{r}),2);
end
ftprnt=imgaussfilt(toimg(ftprnt,sz(2),sz(1)),2);
V_trace=-movVec'*tovec(ftprnt);

figure(13); clf; tiledlayout(2,2)
show_footprnt(ftprnt,Result.ref_im)
nexttile([1 2]);
plot(rescale2(V_trace,1)+[1 2])

V_trace=V_trace'./get_threshold(V_trace',1);
tr_norm= Vics_trace'-movprc(Vics_trace',100,30,2);
tr_norm= tr_norm./get_threshold(tr_norm,1);
spike= find_spike_bh(tr_norm,5,2);

show_traces_spikes(V_trace,spike,Result.Blue)

%% putative co-inhibition
toi=[3730 3801;2459 2500;6405 6444];
t_seg=zeros(1,size(mov_res,3));
for t=1:size(toi,1);
t_seg(1,toi(t,1):toi(t,2))=1;
end
mov_seg=imgaussfilt(mov_res(:,:,find(t_seg>0)),2);

mov_res_mask=mov_seg.*double(max(Result.bvMask,[],3)==0);
subMov=tovec(imresize(imgaussfilt3(mov_res_mask(bound:end-bound,bound:end-bound,:),[1 1 0.1]),1/2));
subMov=subMov-mean(subMov,2);
covMat=subMov'*subMov;

[V, D] = eig(covMat);
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V = V.*vSign;

nPCs=28;
eigImg=toimg(tovec(mov_seg(:,:,:))*V(:,1:nPCs),size(mov_res,1),size(mov_res,2));
figure(4); clf; 
for n=1:nPCs
    nexttile([1 1])
    imshow2(eigImg(:,:,n),[])
    title(['PC #', num2str(n), ' Fraction : ' num2str(D(n)/sum(D),2)])
end

[V_ics, mixmat, sepmat]=sorted_ica(V(:,1:nPCs),10);
icsImg=toimg(tovec(mov_seg(:,:,:))*V_ics,size(mov_res,1),size(mov_res,2));
figure(5); clf; 
for n=1:size(V_ics,2)
    nexttile([1 1])
    %show_footprnt_contour(Result.bvMask,icsImg(:,:,n))
    imshow2(icsImg(:,:,n),[])
    title(['ICS #', num2str(n)])
end
colormap('gray')
figure(7); clf;
plot(rescale2(V_ics,1)+[1:size(V_ics,2)])


%% Low Stim plot
nROI=size(Result.ftprnt,3);
SkelDend = Skeletonize_dendrite(Result.ref_im,4,0.01,25);
interDendDist=[];
for i=1:nROI
    i
    for j=1:nROI
        [interDendDist(i,j), ~]=geodesic_distance(SkelDend,get_coord(Result.ftprnt(:,:,i)),get_coord(Result.ftprnt(:,:,j)));
    end
end
coord_1d=dim_reduce(get_coord(Result.ftprnt));
[~, dist_order]=sort(coord_1d,'descend');
som_roi=find(dist_order==1);
geodist=interDendDist(1,:)'.*sign(coord_1d-coord_1d(1));
show_footprnt_contour(Result.ftprnt(:,:,dist_order),Result.ref_im)

tr_norm= Result.traces_bvMask-movprc(Result.traces_bvMask,100,30,2);
tr_norm= tr_norm./get_threshold(tr_norm,1);
Result.spike= find_spike_bh(tr_norm,4,1.5);
F0=tovec(Result.ref_im)'*tovec(Result.ftprnt.*double(max(Result.bvMask,[],3)==0));
NormTrace=Result.traces_bvMask./F0';
subth_trace=get_subthreshold(NormTrace,Result.spike(1,:),5,10);

figure(6); clf;
tiledlayout(5,1)
nexttile([1 1])
imshow2(Result.ref_im,[]); hold all
plot(bluePatt{1}([1:end 1],2),bluePatt{1}([1:end 1],1),'color',[0 0.6 1],'LineWidth',2)
ax1=[];
title_str={'No Stim.','Low Stim.','No Stim.','Low Stim.'};
for n=1:4
ax1=[ax1 nexttile([1 1])];
t=[1:2000]+2000*(n-1);
imagesc(NormTrace(dist_order,t),[-2 4]*0.01); hold all
sp=find(Result.spike(1,t));
plot(sp,ones(length(sp),1)+nROI-1,'color','k','marker','^','MarkerFaceColor',[1 0 0],'LineStyle','none')
colormap(ax1(n),turbo);
title(title_str{n})
end
%%
spike_time=tovec(find(Result.spike(1,:))'+[-5:50]);
spike_time(spike_time<1 | spike_time>size(Result.traces,2))=[];
spike_time=unique(spike_time);
spike_erode_trace=zeros(1,size(Result.traces,2));
spike_erode_trace(spike_time)=1;
[~, BlueOffTime]=get_blueoffTrace(zeros(1,length(Result.Blue)),Result.Blue,100);
BlueTime= {BlueOffTime,Result.Blue>0}; % Blue off: b=1; Blue on: b=2; Excitation: inh=1; inhibition: inh=2;

bound=13;
mov_blueOff=mov_res(:,:,BlueTime{1} & spike_erode_trace==0); %off
mov_blueOff_msk=mov_blueOff.*(max(Result.bvMask,[],3)==0);
mov_blueOn=mov_res(:,:,BlueTime{2} & spike_erode_trace==0); %on
mov_blueOn_msk=mov_blueOn.*(max(Result.bvMask,[],3)==0);

subMov=mov_blueOff_msk(bound:end-bound,bound:end-bound,:);
subMov2=mov_blueOn_msk(bound:end-bound,bound:end-bound,:);
subMov=tovec(imgaussfilt3(subMov,[1.5 1.5 0.1]));
subMov2=tovec(imgaussfilt3(subMov2,[1.5 1.5 0.1]));

subMov=subMov-mean(subMov,2);
covMat=subMov'*subMov;
[V, D] = eig(covMat);
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V = V.*vSign;

nPCs=12;
eigImg=toimg(tovec(mov_blueOff)*V(:,1:nPCs),size(mov_res,1),size(mov_res,2));
figure(9); clf; 
for n=1:nPCs
    nexttile([1 1])
    imshow2(imgaussfilt(eigImg(:,:,n),1),[])
    title(['Blue Off PC #', num2str(n), ' Fraction : ' num2str(D(n)/sum(D),2)])
end
%colormap(turbo)

subMov2=subMov2-mean(subMov2,2);
covMat2=subMov2'*subMov2;
[V2, D] = eig(covMat2);
D = diag(D);
D = D(end:-1:1);
V2 = V2(:,end:-1:1);
vSign = sign(max(V2) - max(-V2));  % make the largest value always positive
V2 = V2.*vSign;

nPCs=12;
eigImg=toimg(tovec(mov_blueOn)*V2(:,1:nPCs),size(mov_res,1),size(mov_res,2));
figure(10); clf; 
for n=1:nPCs
    nexttile([1 1])
    imshow2(imgaussfilt(eigImg(:,:,n),1),[])
    title(['Blue On PC #', num2str(n), ' Fraction : ' num2str(D(n)/sum(D),2)])
end
%colormap(turbo)

%%
putative_inh_time=[1930 2784 3820 6320]; nTau=[-200:200];
Rfixed = imref2d([size(Result.ref_im,1) size(Result.ref_im,2)]);
inverseTform = invert(Result.tform);
revertedStruct = mat2gray(imwarp(Result.Structure, inverseTform,'OutputView',Rfixed));
revertedStruct_bin = imwarp(Result.Structure_bin, inverseTform,'OutputView',Rfixed);
revertedStruct_bin=revertedStruct_bin>0.05;
for j=1:length(putative_inh_time)
mov_func=imgaussfilt3(-mov_res(:,:,putative_inh_time(j)+nTau),[7 7 0.1])./imgaussfilt(Result.ref_im,3);
  colorMov=grs2rgb(tovec(mov_func.*revertedStruct_bin),jet(300),-0.04,0.07);
  colorMov=reshape(colorMov,size(mov_res,1),size(mov_res,2),length(nTau),[]);
  colorMov=permute(colorMov,[1 2 4 3]);
 mov_show=colorMov.*revertedStruct*5;
 mov_show=mov_show(bound:end-bound,bound:end-bound,:,:);
 writeMov4d([fpath{100} '/PutativeInhibitMov_' num2str(j)],mov_show,nTau,10,1,[])
end
%%
[~, BlueOffTime]=get_blueoffTrace(zeros(1,length(Result.Blue)),Result.Blue,100);
BlueTime= {BlueOffTime,Result.Blue>0}; % Blue off: b=1; Blue on: b=2; Excitation: inh=1; inhibition: inh=2;
thres=[0.005 0.002];
basalSub=movmean(mean(subth_trace(BA_Rois{1},:),1),5,'omitnan'); basalSub(isnan(basalSub))=0;
apicalSub=movmean(mean(subth_trace(BA_Rois{2},:),1),5,'omitnan'); apicalSub(isnan(apicalSub))=0;

basalTr=[]; apicalTr=[];
transBasal=[]; transApical=[];
for b=1:2
% basalTr{b}=basalSub(BlueTime{b})-median(basalSub(BlueTime{b}));
% apicalTr{b}=apicalSub(BlueTime{b})-median(apicalSub(BlueTime{b}));
basalTr{b}=basalSub(BlueTime{b})-movmedian(basalSub(BlueTime{b}),500);
apicalTr{b}=apicalSub(BlueTime{b})-movmedian(apicalSub(BlueTime{b}),500);
    for inh=1:2
        signfunc=(1-(inh-1)*2);
transBasal{b,inh}=detect_transient(basalTr{b}*signfunc,thres,ones(1,length(basalTr{b})));
transApical{b,inh}=detect_transient(apicalTr{b}*signfunc,thres,ones(1,length(apicalTr{b})));

bAPtrans=unique(transBasal{b,inh}.interval.*spike_erode_trace(BlueTime{b})); bAPtrans=bAPtrans(2:end);
transBasal{b,inh}.amp(bAPtrans)=[]; transBasal{b,inh}.length(bAPtrans)=[]; transBasal{b,inh}.int(bAPtrans)=[];
transBasal{b,inh}.interval(ismember(transBasal{b,inh}.interval,bAPtrans))=0;
bAPtrans=unique(transApical{b,inh}.interval.*spike_erode_trace(BlueTime{b}));  bAPtrans=bAPtrans(2:end);
transApical{b,inh}.amp(bAPtrans)=[]; transApical{b,inh}.length(bAPtrans)=[]; transApical{b,inh}.int(bAPtrans)=[];
transApical{b,inh}.interval(ismember(transApical{b,inh}.interval,bAPtrans))=0;
    end
end

figure(8); clf;
nexttile([1 1])
scatter(transBasal{1,1}.int,transBasal{1,1}.amp,'o','filled'); hold all
scatter(transBasal{2,1}.int,transBasal{2,1}.amp,'o','filled')
xlabel('Area under transient')
ylabel('Amplitude of transient (\DeltaF/F)')
title('Basal excitation')
legend({'Blue off','Blue on'})

nexttile([1 1])
scatter(transApical{1,1}.int,transApical{1,1}.amp,'o','filled'); hold all
scatter(transApical{2,1}.int,transApical{2,1}.amp,'o','filled')
xlabel('Area under transient')
ylabel('Amplitude of transient (\DeltaF/F)')
title('Apical excitation')
legend({'Blue off','Blue on'})

nexttile([1 1])
scatter(transBasal{1,2}.int,-transBasal{1,2}.amp,'o','filled'); hold all
scatter(transBasal{2,2}.int,-transBasal{2,2}.amp,'o','filled')
xlabel('Area under transient')
ylabel('Amplitude of transient (\DeltaF/F)')
title('Basal inhibition')
legend({'Blue off','Blue on'})

nexttile([1 1])
scatter(transApical{1,2}.int,-transApical{1,2}.amp,'o','filled'); hold all
scatter(transApical{2,2}.int,-transApical{2,2}.amp,'o','filled')
xlabel('Area under transient')
ylabel('Amplitude of transient (\DeltaF/F)')
title('Apical inhibition')
legend({'Blue off','Blue on'})

%%
%% Skewness
mov_blueOff=mov_res(:,:,find(BlueTime{1}==1 & spike_erode_trace==0)); %off
mov_blueOff=mov_blueOff.*(max(Result.bvMask,[],3)==0);
mov_blueOff=imgaussfilt3(mov_blueOff,[2 2 0.1]);
mov_blueOn=mov_res(:,:,find(BlueTime{2}==1 & spike_erode_trace==0)); %on
mov_blueOn=mov_blueOn.*(max(Result.bvMask,[],3)==0);
mov_blueOn=imgaussfilt3(mov_blueOn,[2 2 0.1]);

skew_blueOff=skewness(mov_blueOff,1,3);
skew_blueOff(isnan(skew_blueOff))=0;
skew_blueOn=skewness(mov_blueOn,1,3);
skew_blueOn(isnan(skew_blueOn))=0;

figure(11); clf;
nexttile([1 1])
imshow2(skew_blueOff,[])
nexttile([1 1])
imshow2(skew_blueOn,[])
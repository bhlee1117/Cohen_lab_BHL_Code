clear
clc;
cd '/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'B5:K128');

save_to='/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,5));
CamType=raw(:,3);
StructureData=raw(:,10);
%% blood vessel masking
[u,s,v] = svds(tovec(mov_res),20);
reshape_u=reshape(u,sz(2),sz(1),[]);
Result.bvMask=[];
[~, Result.bvMask]=get_ROI(max(abs(reshape_u),[],3),Result.bvMask);
%% Low stim
for i=100
    load([fpath{i} '/Result.mat'])
    load(fullfile(fpath{i},"output_data.mat"))
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    mov_mc=double(readBinMov([fpath{i} '/mc_ShutterReg01.bin'],sz(2),sz(1)));
end
[blueDMDimg bluePatt]=get_blueDMDPatt(Device_Data);

%%
bound=5;
mean_F=squeeze(mean(mov_mc(bound:end-bound,bound:end-bound,:),[1 2]));
[~, blueOff]=get_blueoffTrace(mean_F,[Result.Blue],70);
[y_fit]=expfitDM_2(find(blueOff)',mean_F(find(blueOff)),[1:size(mov_mc,3)]',1000);
bkg(1,:)=y_fit;
mov_res=mov_mc-mean(mov_mc,3);
mov_res = SeeResiduals(mov_res,Result.mc);
mov_res = SeeResiduals(mov_res,Result.mc.^2);
mov_res = SeeResiduals(mov_res,Result.mc(:,1).*Result.mc(:,end));
mov_res= SeeResiduals(mov_res,bkg,1);

%%

Struct_valid=find(1-cell2mat(cellfun(@(x) sum(isnan(x)), StructureData, 'UniformOutput', false)));

for i=Struct_valid(39)'
    load(fullfile(fpath{i},'Result.mat'))
    StructureStack=mat2gray(double(tiffreadVolume(StructureData{i})));
    StructureStack(StructureStack==0)=median(StructureStack(:));
    StructureStack=StructureStack(:,:,27:190);
    %StructureStack_med=medfilt2_mov(StructureStack,[15 15]);
    illumination_field=imgaussfilt(max(StructureStack,[],3),50);
    StructureStack=StructureStack./illumination_field;
    StructureStack_Gauss=imgaussfilt3(StructureStack,[7 7 0.1]);
    %StructureStack_med(StructureStack_med==0)=median(StructureStack_med(:));
    %StructureStack=(StructureStack-StructureStack_med)./StructureStack_med;
    StructureStack_filt=(StructureStack-StructureStack_Gauss);
    StructureStack_filt=mat2gray(StructureStack_filt);
    StructureStack_bin=[]; level=[];
    level = graythresh(StructureStack_filt);
    StructureStack_bin=StructureStack_filt>level*1.25;
    moviefixsc(StructureStack_bin)

    se = strel('sphere', 1);
    StructureStack_bin = imdilate(StructureStack_bin, se);
    bwSeg=bwlabeln(StructureStack_bin);
    segments = regionprops3(bwSeg,'Volume','EquivDiameter');
    segments = table2array(segments);
    %bwlist=find(arrayfun(@(x) x.Volume>4000, segments) & arrayfun(@(x) x.EquivDiameter>10, segments));
    bwlist = find(segments(:,1)>7000 & segments(:,2)>10);

    se = strel('sphere',1);
    dendrite_bin=double(ismember(bwSeg,bwlist));
    dendrite_bin= imdilate(dendrite_bin,se);
    dendrite_bin= imgaussfilt3(dendrite_bin,2);

    figure(3); clf;
    imshow2(max(dendrite_bin,[],3),[])
    g=1; ROIrmv=[];
    while g
        h = drawpolygon('Color','r');
        if size(h.Position,1)==1 %no more ROI
            g=0;
        else
            ROIrmv=[ROIrmv; {h.Position}];
            hold all
            plot(h.Position(:,1),h.Position(:,2))
        end
    end
    ROIrmvmask=roi2mask(ROIrmv,size(dendrite_bin,1),size(dendrite_bin,2));
    close(figure(3));
    dendrite_bin(repmat(ROIrmvmask,1,1,size(dendrite_bin,3)))=0;
    figure(3); clf;
    imshow2(max(dendrite_bin,[],3),[])

    StructureStack_final = double(StructureStack).* dendrite_bin;
    figure(2); clf;
    imshow2(imfuse(mat2gray(max(StructureStack_final,[],3)),mat2gray(max(StructureStack,[],3))),[])

    rot_ang=90;
    Structure_ref=(imrotate(StructureStack_final,rot_ang));
    ref_img=Result.ref_im; ref_img(ref_img<prctile(ref_img(:),20))=median(ref_img(:)); ref_img=ref_img-prctile(ref_img(:),20);
    [RegImg,tformReg]=imReg(ref_img,max(Structure_ref,[],3));
    saveastiff(uint16(mat2gray(Structure_ref)*255), [fpath{i} 'Structure.tiff']);

    Result.Structure=max(Structure_ref,[],3);
    Result.Structure_bin=max(imrotate(dendrite_bin,rot_ang),[],3);
    Result.tform=tformReg;
    save(fullfile(fpath{i},'Result.mat'),'Result','fpath','-v7.3')
end
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

tr_norm= Result.traces-movprc(Result.traces,100,30,2);
tr_norm= tr_norm./get_threshold(tr_norm,1);
Result.spike= find_spike_bh(tr_norm,5,2);
F0=tovec(Result.ref_im)'*tovec(Result.ftprnt);
NormTrace=Result.traces./F0';
subth_trace=get_subthreshold(NormTrace,Result.spike(1,:),5,10);

ROIinDMD=[4 5 6 8 12];

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

%% pca analysis
bound=7;
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
%%
BA_Rois={20,[5 8]}; 

spike_time=tovec(find(Result.spike(1,:))'+[-5:50]);
spike_time(spike_time<1 | spike_time>size(Result.traces,2))=[];
spike_time=unique(spike_time);
spike_erode_trace=zeros(1,size(Result.traces,2));
spike_erode_trace(spike_time)=1;

figure(6); clf;
tiledlayout(4,4)
nexttile([1 1])

basalTr=mean(NormTrace(BA_Rois{1},:),1); apicalTr=mean(NormTrace(BA_Rois{2},:),1);
basalSubdVdt=get_slope(mean(subth_trace(BA_Rois{1},:),1),10);
apicalSubdVdt=get_slope(mean(subth_trace(BA_Rois{2},:),1),10);

[BAcorr_bAP]=corr(basalTr(spike_time)',apicalTr(spike_time)');
[BAcorr_nonbAP]=corr(basalTr(~spike_erode_trace)',apicalTr(~spike_erode_trace)');
plot(basalTr,apicalTr,'k.'); hold all
ll=plot(basalTr(spike_time),apicalTr(spike_time),'ro');
xlabel('Basal branch V')
ylabel('Apical branch V')
legend(ll,'bAPs')

nexttile([1 1])
BlueOffdat=[basalTr(find(Result.Blue==0 & spike_erode_trace==0)); apicalTr(find(Result.Blue==0 & spike_erode_trace==0))];
BlueOffdat=BlueOffdat-median(BlueOffdat,2);
BlueOndat=[basalTr(find(Result.Blue>0 & spike_erode_trace==0)); apicalTr(find(Result.Blue>0 & spike_erode_trace==0))];
BlueOndat=BlueOndat-median(BlueOndat,2);

% BlueOffSubdat=[basalSubdVdt(find(Result.Blue==0 & spike_erode_trace==0)); apicalSubdVdt(find(Result.Blue==0 & spike_erode_trace==0))];
% BlueOnSubdat=[basalSubdVdt(find(Result.Blue>0 & spike_erode_trace==0)); apicalSubdVdt(find(Result.Blue>0 & spike_erode_trace==0))];

h=histogram(BlueOffdat(1,:),150,'normalization','probability'); hold all
histogram(BlueOndat(1,:),h.BinEdges,'normalization','probability')
legend({'Blue off','Blue on'})
xlabel('\DeltaF/F-\DeltaF/F_{median}')
ylabel('Probability')
title('Basal branch')

nexttile([1 1])
histogram(BlueOffdat(2,:),h.BinEdges,'normalization','probability'); hold all
histogram(BlueOndat(2,:),h.BinEdges,'normalization','probability')
legend({'Blue off','Blue on'})
xlabel('\DeltaF/F-\DeltaF/F_{median}')
ylabel('Probability')
title('Apical branch')

nexttile([1 1])
h2=histogram(V(find(Result.Blue==0 & spike_erode_trace==0),2),50,'normalization','cdf'); hold all
histogram(V(find(Result.Blue>0 & spike_erode_trace==0),2),h2.BinEdges,'normalization','cdf')
legend({'Blue off','Blue on'})
xlabel('PC #2')
ylabel('Probability')
title('PC#2, Positive when apical is high')

ax1=nexttile([1 4]);
imagesc(NormTrace(dist_order,:),[-2 4]*0.01); hold all
scatter(find(Result.spike(1,:)),ones(sum(Result.spike(1,:)),1)*nROI,'r','^','filled')
title('Voltage kymograph')

ax2=nexttile([2 4]);
% imagesc([moving_corr_size/2+1:9999-moving_corr_size/2],t_lag,BAcorr_mov',[-0.001 0.002])

l=plot(subth_trace(dist_order,:)');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2)); hold all
ylabel('V')
linkaxes([ax1 ax2],'x')
lll=plot(V(:,2),'k');
yyaxis right
plot(Result.Blue,'color',[0 0.6 1])
legend(lll,'PC #2, Positive when apical is high')
xlabel('Time (ms)')
ylabel('Blue')
title('Subthreshold')

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
bound=13;
mov_blueOff=mov_res(:,:,BlueTime{1}); %off
mov_blueOff=mov_blueOff.*(max(Result.bvMask,[],3)==0);
mov_blueOn=mov_res(:,:,BlueTime{2}); %on
mov_blueOn=mov_blueOn.*(max(Result.bvMask,[],3)==0);

inhibit_time=unique([find(ismember(bwlabel(transApical{2,2}.interval),find(transApical{2,2}.amp>0.012))) find(ismember(bwlabel(transBasal{2,2}.interval),find(transBasal{2,2}.amp>0.012)))]);
excite_time=unique([find(ismember(bwlabel(transApical{1,1}.interval),find(transApical{1,1}.amp>0.012))) find(ismember(bwlabel(transBasal{1,1}.interval),find(transBasal{1,1}.amp>0.012)))]);

mov_blueOn_inh=mov_blueOn(:,:,inhibit_time);
mov_blueOff_exc=mov_blueOff(:,:,excite_time);

subMov=mov_blueOn_exc(bound:end-bound,bound:end-bound,:);
subMov2=mov_blueOn_inh(bound:end-bound,bound:end-bound,:);
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
eigImg=toimg(tovec(mov_blueOn_exc)*V(:,1:nPCs),size(mov_res,1),size(mov_res,2));
figure(9); clf; 
for n=1:nPCs
    nexttile([1 1])
    imshow2(imgaussfilt(eigImg(:,:,n),2),[])
    title(['Exc. PC #', num2str(n), ' Fraction : ' num2str(D(n)/sum(D),2)])
end
colormap(turbo)

subMov2=subMov2-mean(subMov2,2);
covMat2=subMov2'*subMov2;
[V2, D] = eig(covMat2);
D = diag(D);
D = D(end:-1:1);
V2 = V2(:,end:-1:1);
vSign = sign(max(V2) - max(-V2));  % make the largest value always positive
V2 = V2.*vSign;

nPCs=12;
eigImg=toimg(tovec(mov_blueOn_inh)*V2(:,1:nPCs),size(mov_res,1),size(mov_res,2));
figure(10); clf; 
for n=1:nPCs
    nexttile([1 1])
    imshow2(imgaussfilt(eigImg(:,:,n),2),[])
    title(['inh. PC #', num2str(n), ' Fraction : ' num2str(D(n)/sum(D),2)])
end
colormap(turbo)

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
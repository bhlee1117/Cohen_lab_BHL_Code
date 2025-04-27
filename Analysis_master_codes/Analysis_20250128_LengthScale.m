%% Set the path
clear
clc;
cd '/Users/bhlee1117/Documents/GitHub/Cohen_lab_BHL_Code/Analysis_master_codes';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'C5:P164');

save_to='/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,5));
CamType=raw(:,3);
StructureData=raw(:,10);
time_segment=25000;
PixelSize=cell2mat(cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false));
[~, unqInd] = unique([Mouse NeuronInd] ,'row');

%% Get F0PCA image
bound=7; i=72;
SameCellInd=find(Mouse==Mouse(i) & NeuronInd==NeuronInd(i));
mov_res_small=[];
for j=SameCellInd'
j
    load([fpath{j} '/OP_Result.mat'])

    mov_mc=readBinMov_BHL(fpath{j},3);
    mean_F=squeeze(mean(mov_mc(bound:end-bound,bound:end-bound,:),[1 2]));
    [~, blueOff]=get_blueoffTrace(mean_F,[Result.Blue],70);
    [y_fit]=expfitDM_2(find(blueOff)',mean_F(find(blueOff)),[1:size(mov_mc,3)]',1000);

    mov_res= mov_mc-mean(mov_mc,3);
    mov_res = SeeResiduals(mov_res,y_fit);
    mov_res = SeeResiduals(mov_res,Result.mc(:,:));
    mov_res = SeeResiduals(mov_res,Result.mc(:,:).^2);
    mov_res = SeeResiduals(mov_res,Result.mc(:,1).*Result.mc(:,end));

    if isfield(Result,'dirtTrace')
    validFrm=find(sum(isnan(Result.dirtTrace),1)==0);
    else
    validFrm=[1:size(mov_res,3)];
    end
    mov_res=movmean(mov_res,15,3);
    mov_res=mov_res(:,:,validFrm);
    mov_res_small=cat(3,mov_res_small,imresize(mov_res,1/4));
end

[V D]=get_eigvector(tovec(mov_res_small),10);
F0img=sqrt(sum((V.^2).*D(1:10)',2));
F0img=toimg(F0img,size(mov_res_small,1),size(mov_res_small,2));
F0img=imresize(F0img,4);   

for j=SameCellInd'
    load([fpath{j} '/OP_Result.mat'])
    Result.F0PCAimg=F0img;
    save([fpath{j} '/OP_Result.mat'],'Result','-v7.3')
end


%% Load the data

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
    [Result.blueDMDimg Result.bluePatt]=get_blueDMDPatt(Device_Data,'stack');

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

%dirtMov_dilate = tracking_dirt(mov_res,0.3);
% %%
% [u,s,v] = svds(tovec(mov_res(:,:,1:5000)-mean(mov_res(:,:,1:5000),3)),20);
% reshape_u=reshape(u,sz(2),sz(1),[]);
% bvMask=[];
% [~, bvMask]=get_ROI(max(abs(reshape_u),[],3),bvMask);
% 
% Result.bvMask=bvMask;
 %Result.traces_bvMask=(tovec(mov_res.*double(max(Result.bvMask,[],3)==0))'*tovec(Result.ftprnt))';
% Result.dirtTrace=(tovec(dirtMov_dilate)'*tovec(Result.ftprnt))';
 %% Get STAs
 i=102; bound=6; nTau=[-30:50];
 load([fpath{i} '/OP_Result.mat'])

%SomTr=Result.traces_bvMask(1,:);
%SomTr=SomTr/get_threshold(SomTr,1);
%Result.spike(1,:)=find_spike_bh(SomTr,5,3);
VoltageTr=Result.normTraces;
VoltageTr(Result.dirtTrace>0)=NaN;
%Result.F0_PCA=get_F0PCA(VoltageTr,3);
VoltageTr=movmean(VoltageTr./Result.F0_PCA,3,2,'omitnan');
%VoltageTr=pcafilterTrace(VoltageTr,[1:15]);

NaNFrm=find(sum(Result.dirtTrace,1)>0);
% mov_res_filt=-mov_res.*dou;
% mov_res_filt(:,:,NaNFrm)=NaN;
% mov_vec_sub=imresize(mov_res_filt(:,:,setdiff([1:size(mov_res,3)],NaNFrm)),1/4);
% [V eigVal eigTrace]=get_eigvector(tovec(mov_vec_sub),10);
% 
% mov2dfilt = eigTrace*V';
% mov2dfilt = reshape(mov2dfilt', size(mov_vec_sub,1), size(mov_vec_sub,2), []);
% mov2dfilt = mov2dfilt + mean(mov_vec_sub,3);
% mov2dfilt=imresize(mov2dfilt,4);
% mov_vec=zeros(size(mov2dfilt,1).*size(mov2dfilt,2),size(mov_res,3));
% mov_vec(:,setdiff([1:size(mov_res,3)],NaNFrm))=tovec(mov2dfilt);

mov_vec=tovec(-mov_res);
mov_vec(:,NaNFrm)=NaN;

STAmov=[]; STAtr=[]; N_avg=[];
for b=1:max(DMDbluetrace)
PattBlue=double(DMDbluetrace==b);
onsetTime=find((PattBlue(2:end)-PattBlue(1:end-1))==1)+1;
PattBlueOnset=ind2vec(length(DMDbluetrace),onsetTime,1);
[~, spMat]=get_STA(Result.spike(1,:),PattBlueOnset,-nTau(1),nTau(end));
Nsp=squeeze(sum(spMat,3));
[~, m]=get_STA(mov_vec,PattBlueOnset,-nTau(1),nTau(end));
STAmov(:,:,:,b)=toimg(squeeze(mean(m(:,Nsp==0,:),2,'omitnan')),sz(2),sz(1));
[~, trMat]=get_STA(VoltageTr,PattBlueOnset,-nTau(1),nTau(end));
STAtr(:,:,b)=squeeze(mean(trMat(:,Nsp==0,:),2,'omitnan'));
N_avg(b)=length(find(Nsp==0));
end

STAmov=STAmov./Result.F0PCAimg;
STAmov=STAmov-median(STAmov(:,:,1:15,:),3,'omitnan');
STAtr=STAtr-median(STAtr(:,1:15,:),2,'omitnan');
%%
b_show=[1 2 5 7 8];
avg_frame=[40:50]; cax_sub=[-0.001 0.007]; bin_tick=40;
avg_frame_trace=[35:45];

STAmov_filt=zeros(size(STAmov));
for b=b_show
STAmov_filt(:,:,:,b)=pcafilt(imaveragefilt(STAmov(:,:,:,b),30),3);
end

STAimg_sub=squeeze(mean(STAmov_filt(:,:,avg_frame,:),3,'omitnan'));
mask=max(Result.ftprnt>0,[],3);

colorSTA=grs2rgb(tovec(STAimg_sub),colormap('turbo'),cax_sub(1),cax_sub(2));
colorSTA=reshape(colorSTA,size(STAmov,1),size(STAmov,2),10,[]);
colorSTA=permute(colorSTA,[1 2 4 3]);
STAimg_sub_Struc=colorSTA.*mat2gray(Result.ref_im-100)*4;

blueCoord=get_coord(Result.blueDMDimg);
ftprntCoord=get_coord(Result.ftprnt);
D=distance_BH(blueCoord,ftprntCoord);
D2=distance_BH(blueCoord,ftprntCoord(Result.dist_order,:));

figure(13); clf; tiledlayout(length(b_show),2);
Dsign=ones(1,size(Result.interDendDist,2));
Dsign(Result.dist_order(1:find(Result.dist_order==1)-1))=-1;
dendaxis=Result.interDendDist(1,:).*Dsign;
dendaxis=dendaxis*PixelSize(i);
dendaxis=dendaxis(Result.dist_order);

for b=b_show 
nexttile([1 1])
Show_kymo=pcafilterTrace(STAtr(Result.dist_order,:,b),[1:10]);
imagesc(nTau,[1:length(dendaxis)],Show_kymo,cax_sub); hold all
%[binnedZ binX binY]=show3Dbinning({[[NaN; dendaxis'] [nTau; Show_kymo]]}, nTau(1:3:end), [min(dendaxis):39:max((dendaxis))], 'image'); hold all
%pcolor(binX',binY',binnedZ');
cb=colorbar;
xlabel('Peri-stimulation time (ms)')
cb.Label.String='\DeltaF/F';
shading interp
set(gca,'YDir','reverse')

[~, roi_stim]=min(D2(b,:));
%plot(2,dendaxis(roi_stim),'>','MarkerFaceColor',[1 0 0]);
plot(2,(roi_stim),'>','MarkerFaceColor',[1 0 0]);
set(gca,'YTick',[1 find(Result.dist_order==1) length(dendaxis)],'YTickLabel',num2str([min(dendaxis) 0 max((dendaxis))]',3))

nexttile([1 1])
imshow2(STAimg_sub_Struc(:,:,:,b),[]); hold all

hold all
plot(Result.bluePatt{b}{1}(:,2),Result.bluePatt{b}{1}(:,1),'color',[0 0.6 1],'LineWidth',1.5)
title(['Number of stimulation averaged: ', num2str(N_avg(b))])
end
colormap(turbo);

figure(16); clf;
[~, arg]=min(D(b_show,:),[],2);
dff_dist=squeeze(mean(STAtr(:,avg_frame_trace,b_show),2,'omitnan'));
l=plot(Result.interDendDist(arg,:)'*PixelSize(i),dff_dist,'.','markersize',25);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(gen_colormap([1 0.5 0; 0 0.5 1],length(b_show)),2));

g=1; dat=[];
for b=b_show
dat{g}=[Result.interDendDist(:,arg(g))*PixelSize(i) dff_dist(:,g)];
g=g+1;
end
[M S B N]=binning_data(dat,[-bin_tick:bin_tick*2:400+bin_tick]);
N=sum(cell2mat(cellfun(@double,N,'UniformOutput',false)'),1);
hold all; 
errorbar(B',M',S'/sqrt(N'),'color',[0 0 0],'LineWidth',2);

exp_model = @(params, x) params(3) + params(1) * exp(params(2) * x);
initial_guess = [max(M) - min(M), -1 / (max(B) - min(B))  , min(M)];

options = optimoptions('lsqcurvefit', 'FunctionTolerance', 1e-6, 'OptimalityTolerance', 1e-6, 'MaxIterations', 10000);
params_fit = lsqcurvefit(exp_model, initial_guess, B, M, [], [], options);

% Compute fitted values
y_fit = exp_model(params_fit, B);

% Compute R^2
SS_res = sum((M - y_fit).^2);  % Residual sum of squares
SS_tot = sum((M - mean(M)).^2);  % Total sum of squares
R_squared = 1 - (SS_res / SS_tot);

% Compute Length Scale (1/|b|)
length_scale = 1 / abs(params_fit(2));

title(['Fit func. : a*exp(bx)+c , 1/|b| = ' num2str(length_scale,3) ' \mum, R^2 = ' num2str(R_squared,3)])
colormap(gen_colormap([1 0.5 0; 0 0.5 1]));
cb=colorbar;
cb.Label.String='Basal to apical';
cb.Ticks=[];
xlabel('Pairwise distance from stimulated ROI (\mum)')
ylabel('\DeltaF/F')
%% Get F0PCA image % Cell145
bound=7; i=145;
SameCellInd=find(Mouse==Mouse(i) & NeuronInd==NeuronInd(i));
mov_res_small=[];
for j=[145]%SameCellInd'
j
    load([fpath{j} '/OP_Result.mat'])

    mov_mc=readBinMov_BHL(fpath{j},3);
    mean_F=squeeze(mean(mov_mc(bound:end-bound,bound:end-bound,:),[1 2]));
    [~, blueOff]=get_blueoffTrace(mean_F,[Result.Blue],70);
    [y_fit]=expfitDM_2(find(blueOff)',mean_F(find(blueOff)),[1:size(mov_mc,3)]',1000);

    mov_res= mov_mc-mean(mov_mc,3);
    mov_res = SeeResiduals(mov_res,y_fit);
    mov_res = SeeResiduals(mov_res,Result.mc(:,:));
    mov_res = SeeResiduals(mov_res,Result.mc(:,:).^2);
    mov_res = SeeResiduals(mov_res,Result.mc(:,1).*Result.mc(:,end));

    dirtMov_dilate = tracking_dirt(mov_res,0.3);
    Result.dirtTrace=(tovec(dirtMov_dilate)'*tovec(Result.ftprnt))';
    save([fpath{j} '/OP_Result.mat'],'Result','-v7.3')

    validFrm=find(sum(isnan(Result.dirtTrace),1)==0);
    mov_res=movmean(mov_res,15,3);
    mov_res=mov_res(:,:,validFrm);
    mov_res_small=cat(3,mov_res_small,imresize(mov_res,1/4));   
end

[V D]=get_eigvector(tovec(mov_res_small),10);
F0img=sqrt(sum((V.^2).*D(1:10)',2));
F0img=toimg(F0img,size(mov_res_small,1),size(mov_res_small,2));
F0img=imresize(F0img,4);   

Result.F0PCAimg=F0img;

for j=SameCellInd'
    load([fpath{j} '/OP_Result.mat'])
    Result.F0PCAimg=F0img;
    save([fpath{j} '/OP_Result.mat'],'Result','-v7.3')
end


%% Load the data

i=145; bound=6;
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
    [Result.blueDMDimg Result.bluePatt]=get_blueDMDPatt(Device_Data,'stack');

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
Result.normTraces=-Result.traces_bvMask./get_threshold(-Result.traces_bvMask,1);
Result.spike=find_spike_bh(Result.normTraces,5,3);
Result.traces_bvMask=(tovec(mov_res.*double(max(Result.bvMask,[],3)==0))'*tovec(Result.ftprnt))';
Result.dirtTrace=(tovec(dirtMov_dilate)'*tovec(Result.ftprnt))';

save([fpath{i} '/OP_Result.mat'],'Result','-v7.3')
 %% Get STAs
 i=145; bound=6; nTau=[-30:50];
 load([fpath{i} '/OP_Result.mat'])

%SomTr=Result.traces_bvMask(1,:);
%SomTr=SomTr/get_threshold(SomTr,1);
%Result.spike(1,:)=find_spike_bh(SomTr,5,3);
VoltageTr=Result.normTraces;
VoltageTr(Result.dirtTrace>0)=NaN;
%Result.F0_PCA=get_F0PCA(VoltageTr,3);
VoltageTr=VoltageTr./Result.F0_PCA;

NaNFrm=find(sum(Result.dirtTrace,1)>0);
% mov_res_filt=-mov_res.*dou;
% mov_res_filt(:,:,NaNFrm)=NaN;
% mov_vec_sub=imresize(mov_res_filt(:,:,setdiff([1:size(mov_res,3)],NaNFrm)),1/4);
% [V eigVal eigTrace]=get_eigvector(tovec(mov_vec_sub),10);
% 
% mov2dfilt = eigTrace*V';
% mov2dfilt = reshape(mov2dfilt', size(mov_vec_sub,1), size(mov_vec_sub,2), []);
% mov2dfilt = mov2dfilt + mean(mov_vec_sub,3);
% mov2dfilt=imresize(mov2dfilt,4);
% mov_vec=zeros(size(mov2dfilt,1).*size(mov2dfilt,2),size(mov_res,3));
% mov_vec(:,setdiff([1:size(mov_res,3)],NaNFrm))=tovec(mov2dfilt);

mov_vec=tovec(-mov_res);
mov_vec(:,NaNFrm)=NaN;

STAmov=[]; STAtr=[]; N_avg=[];
for b=1:max(DMDbluetrace)
PattBlue=double(DMDbluetrace==b);
onsetTime=find((PattBlue(2:end)-PattBlue(1:end-1))==1)+1;
PattBlueOnset=ind2vec(length(DMDbluetrace),onsetTime,1);
[~, spMat]=get_STA(Result.spike(1,:),PattBlueOnset,-nTau(1),nTau(end));
Nsp=squeeze(sum(spMat,3));
[~, m]=get_STA(mov_vec,PattBlueOnset,-nTau(1),nTau(end));
STAmov(:,:,:,b)=toimg(squeeze(mean(m(:,Nsp==0,:),2,'omitnan')),sz(2),sz(1));
[~, trMat]=get_STA(VoltageTr,PattBlueOnset,-nTau(1),nTau(end));
STAtr(:,:,b)=squeeze(mean(trMat(:,Nsp==0,:),2,'omitnan'));
N_avg(b)=length(find(Nsp==0));
end

STAmov=STAmov./Result.F0PCAimg;
STAmov=STAmov-median(STAmov(:,:,10:30,:),3,'omitnan');
STAtr=STAtr-median(STAtr(:,10:30,:),2,'omitnan');
%%
b_show=[1:8];
avg_frame=[40:45]; cax_sub=[-0.001 0.007];
figure(11); clf;

STAmov_filt=zeros(size(STAmov));
for b=b_show
STAmov_filt(:,:,:,b)=pcafilt(imaveragefilt(STAmov(:,:,:,b),10),3);
end

STAimg_sub=squeeze(mean(STAmov_filt(:,:,avg_frame,:),3,'omitnan'));
mask=max(Result.ftprnt>0,[],3);

colorSTA=grs2rgb(tovec(STAimg_sub),colormap('hot'),cax_sub(1),cax_sub(2));
colorSTA=reshape(colorSTA,size(STAmov,1),size(STAmov,2),max(DMDbluetrace),[]);
colorSTA=permute(colorSTA,[1 2 4 3]);
STAimg_sub_Struc=colorSTA.*mat2gray(Result.ref_im-100).*double(max(Result.ftprnt,[],3)>0).*double(max(Result.bvMask,[],3)==0)*3;

blueCoord=get_coord(Result.blueDMDimg);
ftprntCoord=get_coord(Result.ftprnt(:,:,Result.dist_order));
D=distance_BH(blueCoord,ftprntCoord);

figure(13); clf; tiledlayout(length(b_show),2);
Dsign=ones(1,size(Result.interDendDist,2));
Dsign(Result.dist_order(1:find(Result.dist_order==1)-1))=-1;
dendaxis=Result.interDendDist(1,:).*Dsign;
dendaxis=dendaxis*PixelSize(i);
dendaxis=dendaxis(Result.dist_order);

STAtr=STAtr(Result.dist_order,:,:);
l=plot(nTau,movmean(squeeze(mean(STAtr([15 21],:,:),1,'omitnan')),10,1,'omitnan'));
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(max(DMDbluetrace)),2));

figure(14); clf; tiledlayout(4,2);
bvBoundary=bwboundaries(max(Result.bvMask,[],3)>0);
for b=b_show
    nexttile([1 1]);
    imshow2(STAimg_sub_Struc(:,:,:,b),cax_sub); hold all
    plot(Result.bluePatt{b}{1}(:,2),Result.bluePatt{b}{1}(:,1),'color',[0 0.6 1],'linewidth',1.5)
    for bv=1:length(bvBoundary)
    plot(bvBoundary{bv}(:,2),bvBoundary{bv}(:,1),'r')
    end
end
colormap(hot);
colorbar;

%% Get F0PCA image % Cell145
bound=7; i=145;
SameCellInd=find(Mouse==Mouse(i) & NeuronInd==NeuronInd(i));
mov_res_small=[];
for j=[145]%SameCellInd'
j
    load([fpath{j} '/OP_Result.mat'])

    mov_mc=readBinMov_BHL(fpath{j},3);
    mean_F=squeeze(mean(mov_mc(bound:end-bound,bound:end-bound,:),[1 2]));
    [~, blueOff]=get_blueoffTrace(mean_F,[Result.Blue],70);
    [y_fit]=expfitDM_2(find(blueOff)',mean_F(find(blueOff)),[1:size(mov_mc,3)]',1000);

    mov_res= mov_mc-mean(mov_mc,3);
    mov_res = SeeResiduals(mov_res,y_fit);
    mov_res = SeeResiduals(mov_res,Result.mc(:,:));
    mov_res = SeeResiduals(mov_res,Result.mc(:,:).^2);
    mov_res = SeeResiduals(mov_res,Result.mc(:,1).*Result.mc(:,end));

    dirtMov_dilate = tracking_dirt(mov_res,0.3);
    Result.dirtTrace=(tovec(dirtMov_dilate)'*tovec(Result.ftprnt))';
    save([fpath{j} '/OP_Result.mat'],'Result','-v7.3')

    validFrm=find(sum(isnan(Result.dirtTrace),1)==0);
    mov_res=movmean(mov_res,15,3);
    mov_res=mov_res(:,:,validFrm);
    mov_res_small=cat(3,mov_res_small,imresize(mov_res,1/4));   
end

[V D]=get_eigvector(tovec(mov_res_small),10);
F0img=sqrt(sum((V.^2).*D(1:10)',2));
F0img=toimg(F0img,size(mov_res_small,1),size(mov_res_small,2));
F0img=imresize(F0img,4);   

Result.F0PCAimg=F0img;

for j=SameCellInd'
    load([fpath{j} '/OP_Result.mat'])
    Result.F0PCAimg=F0img;
    save([fpath{j} '/OP_Result.mat'],'Result','-v7.3')
end


%% Load the data

i=145; bound=6;
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
    [Result.blueDMDimg Result.bluePatt]=get_blueDMDPatt(Device_Data,'stack');

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
Result.normTraces=-Result.traces_bvMask./get_threshold(-Result.traces_bvMask,1);
Result.spike=find_spike_bh(Result.normTraces,5,3);
Result.traces_bvMask=(tovec(mov_res.*double(max(Result.bvMask,[],3)==0))'*tovec(Result.ftprnt))';
Result.dirtTrace=(tovec(dirtMov_dilate)'*tovec(Result.ftprnt))';

save([fpath{i} '/OP_Result.mat'],'Result','-v7.3')
 %% Get STAs
 i=145; bound=6; nTau=[-30:50];
 load([fpath{i} '/OP_Result.mat'])

%SomTr=Result.traces_bvMask(1,:);
%SomTr=SomTr/get_threshold(SomTr,1);
%Result.spike(1,:)=find_spike_bh(SomTr,5,3);
VoltageTr=Result.normTraces;
VoltageTr(Result.dirtTrace>0)=NaN;
%Result.F0_PCA=get_F0PCA(VoltageTr,3);
VoltageTr=VoltageTr./Result.F0_PCA;

NaNFrm=find(sum(Result.dirtTrace,1)>0);
% mov_res_filt=-mov_res.*dou;
% mov_res_filt(:,:,NaNFrm)=NaN;
% mov_vec_sub=imresize(mov_res_filt(:,:,setdiff([1:size(mov_res,3)],NaNFrm)),1/4);
% [V eigVal eigTrace]=get_eigvector(tovec(mov_vec_sub),10);
% 
% mov2dfilt = eigTrace*V';
% mov2dfilt = reshape(mov2dfilt', size(mov_vec_sub,1), size(mov_vec_sub,2), []);
% mov2dfilt = mov2dfilt + mean(mov_vec_sub,3);
% mov2dfilt=imresize(mov2dfilt,4);
% mov_vec=zeros(size(mov2dfilt,1).*size(mov2dfilt,2),size(mov_res,3));
% mov_vec(:,setdiff([1:size(mov_res,3)],NaNFrm))=tovec(mov2dfilt);

mov_vec=tovec(-mov_res);
mov_vec(:,NaNFrm)=NaN;

STAmov=[]; STAtr=[]; N_avg=[];
for b=1:max(DMDbluetrace)
PattBlue=double(DMDbluetrace==b);
onsetTime=find((PattBlue(2:end)-PattBlue(1:end-1))==1)+1;
PattBlueOnset=ind2vec(length(DMDbluetrace),onsetTime,1);
[~, spMat]=get_STA(Result.spike(1,:),PattBlueOnset,-nTau(1),nTau(end));
Nsp=squeeze(sum(spMat,3));
[~, m]=get_STA(mov_vec,PattBlueOnset,-nTau(1),nTau(end));
STAmov(:,:,:,b)=toimg(squeeze(mean(m(:,Nsp==0,:),2,'omitnan')),sz(2),sz(1));
[~, trMat]=get_STA(VoltageTr,PattBlueOnset,-nTau(1),nTau(end));
STAtr(:,:,b)=squeeze(mean(trMat(:,Nsp==0,:),2,'omitnan'));
N_avg(b)=length(find(Nsp==0));
end

STAmov=STAmov./Result.F0PCAimg;
STAmov=STAmov-median(STAmov(:,:,10:30,:),3,'omitnan');
STAtr=STAtr-median(STAtr(:,10:30,:),2,'omitnan');
%%
b_show=[1:8];
avg_frame=[40:45]; cax_sub=[-0.001 0.007];
figure(11); clf;

STAmov_filt=zeros(size(STAmov));
for b=b_show
STAmov_filt(:,:,:,b)=pcafilt(imaveragefilt(STAmov(:,:,:,b),10),3);
end

STAimg_sub=squeeze(mean(STAmov_filt(:,:,avg_frame,:),3,'omitnan'));
mask=max(Result.ftprnt>0,[],3);

colorSTA=grs2rgb(tovec(STAimg_sub),colormap('hot'),cax_sub(1),cax_sub(2));
colorSTA=reshape(colorSTA,size(STAmov,1),size(STAmov,2),max(DMDbluetrace),[]);
colorSTA=permute(colorSTA,[1 2 4 3]);
STAimg_sub_Struc=colorSTA.*mat2gray(Result.ref_im-100).*double(max(Result.ftprnt,[],3)>0).*double(max(Result.bvMask,[],3)==0)*3;

blueCoord=get_coord(Result.blueDMDimg);
ftprntCoord=get_coord(Result.ftprnt(:,:,Result.dist_order));
D=distance_BH(blueCoord,ftprntCoord);

figure(13); clf; tiledlayout(length(b_show),2);
Dsign=ones(1,size(Result.interDendDist,2));
Dsign(Result.dist_order(1:find(Result.dist_order==1)-1))=-1;
dendaxis=Result.interDendDist(1,:).*Dsign;
dendaxis=dendaxis*PixelSize(i);
dendaxis=dendaxis(Result.dist_order);

STAtr=STAtr(Result.dist_order,:,:);
l=plot(nTau,movmean(squeeze(mean(STAtr([15 21],:,:),1,'omitnan')),10,1,'omitnan'));
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(max(DMDbluetrace)),2));

figure(14); clf; tiledlayout(4,2);
bvBoundary=bwboundaries(max(Result.bvMask,[],3)>0);
for b=b_show
    nexttile([1 1]);
    imshow2(STAimg_sub_Struc(:,:,:,b),cax_sub); hold all
    plot(Result.bluePatt{b}{1}(:,2),Result.bluePatt{b}{1}(:,1),'color',[0 0.6 1],'linewidth',1.5)
    for bv=1:length(bvBoundary)
    plot(bvBoundary{bv}(:,2),bvBoundary{bv}(:,1),'r')
    end
end
colormap(hot);
colorbar;


%% Load the data

i=152; bound=6;
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
    [Result.blueDMDimg Result.bluePatt]=get_blueDMDPatt(Device_Data,'stack');

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
Result.normTraces=-Result.traces_bvMask./get_threshold(-Result.traces_bvMask,1);
Result.spike=find_spike_bh(Result.normTraces,5,3);
Result.traces_bvMask=(tovec(mov_res.*double(max(Result.bvMask,[],3)==0))'*tovec(Result.ftprnt))';
Result.dirtTrace=(tovec(dirtMov_dilate)'*tovec(Result.ftprnt))';

save([fpath{i} '/OP_Result.mat'],'Result','-v7.3')
 %% Get STAs
 i=145; bound=6; nTau=[-30:50];
 load([fpath{i} '/OP_Result.mat'])

%SomTr=Result.traces_bvMask(1,:);
%SomTr=SomTr/get_threshold(SomTr,1);
%Result.spike(1,:)=find_spike_bh(SomTr,5,3);
VoltageTr=Result.normTraces;
VoltageTr(Result.dirtTrace>0)=NaN;
%Result.F0_PCA=get_F0PCA(VoltageTr,3);
VoltageTr=VoltageTr./Result.F0_PCA;

NaNFrm=find(sum(Result.dirtTrace,1)>0);
% mov_res_filt=-mov_res.*dou;
% mov_res_filt(:,:,NaNFrm)=NaN;
% mov_vec_sub=imresize(mov_res_filt(:,:,setdiff([1:size(mov_res,3)],NaNFrm)),1/4);
% [V eigVal eigTrace]=get_eigvector(tovec(mov_vec_sub),10);
% 
% mov2dfilt = eigTrace*V';
% mov2dfilt = reshape(mov2dfilt', size(mov_vec_sub,1), size(mov_vec_sub,2), []);
% mov2dfilt = mov2dfilt + mean(mov_vec_sub,3);
% mov2dfilt=imresize(mov2dfilt,4);
% mov_vec=zeros(size(mov2dfilt,1).*size(mov2dfilt,2),size(mov_res,3));
% mov_vec(:,setdiff([1:size(mov_res,3)],NaNFrm))=tovec(mov2dfilt);

mov_vec=tovec(-mov_res);
mov_vec(:,NaNFrm)=NaN;

STAmov=[]; STAtr=[]; N_avg=[];
for b=1:max(DMDbluetrace)
PattBlue=double(DMDbluetrace==b);
onsetTime=find((PattBlue(2:end)-PattBlue(1:end-1))==1)+1;
PattBlueOnset=ind2vec(length(DMDbluetrace),onsetTime,1);
[~, spMat]=get_STA(Result.spike(1,:),PattBlueOnset,-nTau(1),nTau(end));
Nsp=squeeze(sum(spMat,3));
[~, m]=get_STA(mov_vec,PattBlueOnset,-nTau(1),nTau(end));
STAmov(:,:,:,b)=toimg(squeeze(mean(m(:,Nsp==0,:),2,'omitnan')),sz(2),sz(1));
[~, trMat]=get_STA(VoltageTr,PattBlueOnset,-nTau(1),nTau(end));
STAtr(:,:,b)=squeeze(mean(trMat(:,Nsp==0,:),2,'omitnan'));
N_avg(b)=length(find(Nsp==0));
end

STAmov=STAmov./Result.F0PCAimg;
STAmov=STAmov-median(STAmov(:,:,10:30,:),3,'omitnan');
STAtr=STAtr-median(STAtr(:,10:30,:),2,'omitnan');
%%
b_show=[1:8];
avg_frame=[40:45]; cax_sub=[-0.001 0.007];
figure(11); clf;

STAmov_filt=zeros(size(STAmov));
for b=b_show
STAmov_filt(:,:,:,b)=pcafilt(imaveragefilt(STAmov(:,:,:,b),10),3);
end

STAimg_sub=squeeze(mean(STAmov_filt(:,:,avg_frame,:),3,'omitnan'));
mask=max(Result.ftprnt>0,[],3);

colorSTA=grs2rgb(tovec(STAimg_sub),colormap('hot'),cax_sub(1),cax_sub(2));
colorSTA=reshape(colorSTA,size(STAmov,1),size(STAmov,2),max(DMDbluetrace),[]);
colorSTA=permute(colorSTA,[1 2 4 3]);
STAimg_sub_Struc=colorSTA.*mat2gray(Result.ref_im-100).*double(max(Result.ftprnt,[],3)>0).*double(max(Result.bvMask,[],3)==0)*3;

blueCoord=get_coord(Result.blueDMDimg);
ftprntCoord=get_coord(Result.ftprnt(:,:,Result.dist_order));
D=distance_BH(blueCoord,ftprntCoord);

figure(13); clf; tiledlayout(length(b_show),2);
Dsign=ones(1,size(Result.interDendDist,2));
Dsign(Result.dist_order(1:find(Result.dist_order==1)-1))=-1;
dendaxis=Result.interDendDist(1,:).*Dsign;
dendaxis=dendaxis*PixelSize(i);
dendaxis=dendaxis(Result.dist_order);

STAtr=STAtr(Result.dist_order,:,:);
l=plot(nTau,movmean(squeeze(mean(STAtr([15 21],:,:),1,'omitnan')),10,1,'omitnan'));
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(max(DMDbluetrace)),2));

figure(14); clf; tiledlayout(4,2);
bvBoundary=bwboundaries(max(Result.bvMask,[],3)>0);
for b=b_show
    nexttile([1 1]);
    imshow2(STAimg_sub_Struc(:,:,:,b),cax_sub); hold all
    plot(Result.bluePatt{b}{1}(:,2),Result.bluePatt{b}{1}(:,1),'color',[0 0.6 1],'linewidth',1.5)
    for bv=1:length(bvBoundary)
    plot(bvBoundary{bv}(:,2),bvBoundary{bv}(:,1),'r')
    end
end
colormap(hot);
colorbar;

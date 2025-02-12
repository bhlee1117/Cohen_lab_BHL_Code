
clear
clc;
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_Optopatch_Prism';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'B5:P175');

fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,5));
CamType=raw(:,3);
StructureData=raw(:,10);
StimROI=raw(:,6);
StimWfn=raw(:,7);
isGoodCell=cell2mat(raw(:,11));
PixelSize=cell2mat(cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false));
refROI=cellfun(@(x) (str2num(num2str(x))),raw(:,14),'UniformOutput',false);
maintrunkROI=cellfun(@(x) (str2num(num2str(x))),raw(:,15),'UniformOutput',false);
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
title_str={'Basal','Apical','Peri-Soma'};
[~, unqInd] = unique([Mouse NeuronInd] ,'row');
set(0,'DefaultFigureWindowStyle','docked')


%%
f=35;
cd(fpath{f});
load(fullfile(fpath{f},'OP_Result.mat'),'Result')

figure(7); clf;
nexttile([1 1])
imshow2(Result.Structure,[])
hold all
plot(Result.bluePatt{1}(:,2),Result.bluePatt{1}(:,1),'color',[0 0.3 1])
figure(10); clf;
ftprnt=Result.ftprnt(:,:,Result.dist_order);
ftprnt_show=cat(3,ftprnt(:,:,3),max(ftprnt(:,:,[15 16 17]),[],3));
show_footprnt_contour(ftprnt_show,Result.Structure)

normTr=Result.normTraces./Result.F0_PCA;
normTr=normTr(Result.dist_order,:);

tr=[normTr(3,:); mean(normTr([15 16 17],:),1)];
tr=tr./get_threshold(tr,1);
sp=[];
sp(1,:)=find_spike_bh(tr(1,:),15,5);
[~, maxarg]=max(reshape(tr(2,find(sp(1,:))'+[-4:4]),[],9),[],2);
sp(2,:)=ind2vec(size(tr,2),find(sp(1,:))+maxarg-5,1);


figure(8); clf; %tiledlayout(2,1);
t=[1:size(normTr,2)]/1000;
nexttile([1 1])
plot(t,normTr(3,:),'color',[0.3 0.1 0.5]); hold all
plot(t,mean(normTr([15 16 17],:),1)+0.04,'r'); hold all
plot(t,Result.Blue/100-0.02,'color',[0 0.3 1]);
sp_time=find(sp(1,:));
plot(sp_time(maxarg<5)/1000,0.08,'.','color',[0. 0.3 1],'markersize',24); hold all
plot(sp_time(maxarg>=5)/1000,0.08,'.','color',[1 0.3 0],'markersize',24); hold all
legend({'Soma','Dendrite'})
axis off


figure(9); clf; tiledlayout(3,4);
mainaxis=[2 3 4 8 11 13 15 16 18 19 20];
t=[1:size(normTr,2)]/1000;
nTau=[-5:50];
g=1;
for b=[1 6 53 64]
ax1=nexttile(g,[2 1]);
bwBlue=bwlabel(Result.Blue>0);
t_show=(find(bwBlue==b,1)+nTau);
normTrfilt=pcafilterTrace(normTr,3);
plot(nTau,normTr(3,t_show),'color',[0.3 0.1 0.5]); hold all
plot(nTau,mean(normTr([15 16 17],t_show),1),'r'); hold all
plot(nTau,Result.Blue(t_show)/100-0.02,'color',[0 0.3 1]);
ylim([-0.02 0.035])
axis off

Dsign=ones(1,size(Result.interDendDist,2));
Dsign(Result.dist_order(1:find(Result.dist_order==1)-1))=-1;
dendaxis=Result.interDendDist.*Dsign;
dendaxis=dendaxis(Result.dist_order)*PixelSize(f);

ax2=nexttile(g+8,[1 1]);
imagesc(nTau,dendaxis(mainaxis),normTrfilt(mainaxis,t_show),[0 0.04])
xlabel('Time (ms)')
ylabel(['Distance from' newline 'soma (\mum)'])
linkaxes([ax1 ax2],'x')
xlim([nTau(1) nTau(end)])
g=g+1;
end
colormap(turbo)
colorbar
ax1=nexttile(4,[2 1]);
legend({'Soma','Dendrite'})


figure(11); clf; tiledlayout(3,2);
ax1=nexttile(1,[2 1]);
[dSTA]=get_STA(normTr,ind2vec(size(tr,2),sp_time(maxarg<5),1),10,20);
l=plot([-10:20],dSTA');
ylabel('\DeltaF/F')
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(gen_colormap([0 0 0; 0.5 0 0; 1 0 0],size(normTr,1)),2));
ax2=nexttile(5,[1 1]);
imagesc([-10:20],dendaxis,dSTA)
xlabel('Time (ms)')
ylabel(['Distance from' newline 'soma (\mum)'])

ax1=nexttile(2,[2 1]);
[STA]=get_STA(normTr,ind2vec(size(tr,2),sp_time(maxarg>=5),1),10,20);
l=plot([-10:20],STA');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(gen_colormap([0 0 0; 0.5 0 0; 1 0 0],size(normTr,1)),2));
ylabel('\DeltaF/F')
ax2=nexttile(6,[1 1]);
imagesc([-10:20],dendaxis,STA)
xlabel('Time (ms)')
ylabel(['Distance from' newline 'soma (\mum)'])
colormap(turbo)

%%

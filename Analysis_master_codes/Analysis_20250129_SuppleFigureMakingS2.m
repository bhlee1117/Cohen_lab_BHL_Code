clear; clc;

%load('/Volumes/BHL18TB_D2/YQ601_PlaceCellResults/PF_Result_20250130.mat')

[~, fpath, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'PlaceCellData_Arrangement.xlsx'], 'Sheet1', 'C8:F42');

[~, ~, PF_list]=xlsread(['/Volumes/BHL18TB_D2/Arranged_Data/NaV_InactivationResult/' ...
    'PlaceFieldData_Arrangement_20210112.xlsx'], 'Sheet1', 'A5:E24');

[~, ~, rawPrism] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:AA31');
fpath_prism=rawPrism(:,2);

for i=1:length(fpath)
    Result_tmp=load(fullfile(fpath{i},'PC_Result.mat'));
    PC_Result{i}=Result_tmp.Result;
end

InductionNeuron=cell2mat(raw(:,4));

PF_list= cell2mat(PF_list);
PC_list=[1 1;3 1;3 3;4 2;4 3;5 1;6 1;9 1;12 3;18 1;20 1;23 1;27 1;30 1;32 1;32 2;35 1;35 2];
cmap=[0.5 0.05 0.15;0.1 0.1 0.1]/5;

place_bin=150; 
PlaceFieldList=cellfun(@(x) (str2num(num2str(x))),rawPrism(:,22),'UniformOutput',false);
PlaceFieldBin=cellfun(@(x) (str2num(num2str(x))),rawPrism(:,23),'UniformOutput',false);
%%
% % Show all the place cell with complex spiking map
% 
% 
% f1=figure(1); clf;
% coarse_bin=3;
% for i=1:length(PC_Result)
%     [sz_fprntY, sz_fprntX, nNeurons]=size(PC_Result{i}.c_ftprnt);
%     StimN=find(mean(PC_Result{i}.spike(:,find(PC_Result{i}.Blue>0)),2)>mean(PC_Result{i}.spike,2)*5)';
%     UnstimN=setdiff([1:nNeurons],StimN);
%     for n=[StimN UnstimN]
%         nexttile([1 1])
%         pos_bin=size(PC_Result{i}.Lap_CS,2);
% 
%         Lap_tmp_CS=repmat(PC_Result{i}.Lap_CS(:,:,n),1,3);
%         Lap_tmp_CS=movmean(Lap_tmp_CS,coarse_bin,2);
%         Lap_tmp_CS=(Lap_tmp_CS(:,pos_bin+1:2*pos_bin));
% 
%         Lap_tmp_SS=repmat(PC_Result{i}.Lap_FR(:,:,n),1,3);
%         Lap_tmp_SS=movmean(Lap_tmp_SS,coarse_bin,2);
%         Lap_tmp_SS=(Lap_tmp_SS(:,pos_bin+1:2*pos_bin));
% 
%         composite_LR=(squeeze(sum(cat(3,Lap_tmp_CS,Lap_tmp_SS).*reshape(cmap,1,1,[],3),3)));
%         %imagesc(reshape(rescale(composite_LR(:)),[],pos_bin,3))
%         imagesc(composite_LR)
% 
%         if ismember(n,UnstimN)
%             title(['Unstimulated ' num2str(i) '-' num2str(n)])
%         else
%             title([num2str(i) '-' num2str(n)])
%         end
%     end
% end

foi2=[1 4 5 6 8 10 11 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27];
FR_cannula=[]; FR_prism=[]; SI_cannula=[]; SI_prism=[]; FramN=[]; Brst_prism=[]; Brst_cannula=[];
BrstSz_prism=[]; BrstSz_cannula=[];
for p=1:length(PC_Result)
    FR_cannula=[FR_cannula; mean(PC_Result{p}.spike(:,PC_Result{p}.Blue==0),2,'omitnan')*1000];

    for n=1:size(PC_Result{p}.spike,1)
        sp_tr=PC_Result{p}.spike(n,:);
        sp_tr(find(PC_Result{p}.Blue>0))=NaN;
        ss_time=find(sp_tr); % BS is subclass of SS
        brst=bwlabel((ss_time(2:end)-ss_time(1:end-1))<=20); % SSs that have an ISI shorter than 20 ms are BS.
        brstprop = regionprops("table",brst,"Area");
        Brst_cannula=[Brst_cannula; max(brst)/sum(PC_Result{p}.Blue==0)*1000];
        BrstSz_cannula=[BrstSz_cannula; histcounts(brstprop.Area+1,[0.5:10.5],'Normalization','probability')];
    end

    FR_cannula=[FR_cannula; mean(PC_Result{p}.spike(:,PC_Result{p}.Blue==0),2,'omitnan')*1000];

     binTrack=(ceil(PC_Result{p}.VR(5,:)/((115)/150)));
    for n=1:size(PC_Result{p}.normTraces,1)
        if n==InductionNeuron(p)
validFrame=(PC_Result{p}.VR(2,:)>=2);
validFrame(find(PC_Result{p}.Blue>0,1):end)=0;
        else
validFrame=(PC_Result{p}.VR(2,:)>=2);
        end
        FramN=[FramN; sum(validFrame)];
        SI_cannula = [SI_cannula; SpatialInfo(PC_Result{p}.spike(n,validFrame), binTrack(validFrame))];
    end
end

for f=foi2
    load(fullfile(fpath_prism{f},'PC_Result.mat'),'Result');
    FR_prism=[FR_prism; mean(Result.spike(1,Result.Blue==0 & Result.motionReject==0),2,'omitnan')*1000];

        sp_tr=double(Result.spike(n,:));
        sp_tr(find(Result.Blue>0))=NaN;
        ss_time=find(sp_tr); % BS is subclass of SS
        brst=bwlabel((ss_time(2:end)-ss_time(1:end-1))<=20); % SSs that have an ISI shorter than 20 ms are BS.
        brstprop = regionprops("table",brst,"Area");
        Brst_prism=[Brst_prism; max(brst)/sum(Result.Blue==0)*1000];
        BrstSz_prism=[BrstSz_prism; histcounts(brstprop.Area+1,[0.5:10.5],'Normalization','probability')];

    binTrack=(ceil(Result.VR(5,:)/((115)/150)));
    validFrame=(Result.VR(2,:)>=2 & Result.Blue==0);
    SI_prism = [SI_prism; SpatialInfo(Result.spike(1,validFrame), binTrack(validFrame))];
end

figure(6); clf; tiledlayout(1,5);
cmap=[0 0.3 1; 1 0.3 0];
nexttile([1 1]);
Boxplot_wPoints(FR_prism, FR_cannula, cmap);
ylabel('Firing rate (Hz)')
set(gca,'XTick',[1 2],'XTickLabel',{'Prism','Cannula'})

nexttile([1 1]);
Boxplot_wPoints(Brst_prism, Brst_cannula, cmap)
ylabel('Burst rate (Hz)')
set(gca,'XTick',[1 2],'XTickLabel',{'Prism','Cannula'})

nexttile([1 2]);
M_prism=mean(BrstSz_prism,1,'omitnan');
M_cannula=mean(BrstSz_cannula,1,'omitnan');
S_prism=std(BrstSz_prism,0,1,'omitnan')./sqrt(size(BrstSz_prism,1));
S_cannula=std(BrstSz_cannula,0,1,'omitnan')./sqrt(size(BrstSz_cannula,1));
b=bar([1:10],[M_prism; M_cannula],'grouped','FaceAlpha',0.5,'BarWidth',1); hold all
b(1).FaceColor=cmap(1,:); b(2).FaceColor=cmap(2,:);
errorbar([1:10]-0.15,M_prism', S_prism','color',cmap(1,:),'LineStyle','none')
errorbar([1:10]+0.15,M_cannula', S_prism','color', cmap(2,:),'LineStyle','none')
ylabel('Probability')
xlabel('Burst size (# spike)')
legend(b,{'Prism','Cannula'})
box off; ylim([0 0.6]);

nexttile([1 1]);
Boxplot_wPoints(SI_prism, SI_cannula(FramN>100000), cmap)
ylabel('Spatial information (bits/spike)')
set(gca,'XTick',[1 2],'XTickLabel',{'Prism','Cannula'})

%% Show good PCs

show_PC=10; % index list from PFs
pos_bin=150;
PC_ind=PC_list(PF_list(show_PC,1),:);
nVRbin=[-50:50];

Result_ST=importdata(fullfile(fpath{PC_ind(1)},'PC_Result.mat'),'Result');
Result_Pr=importdata(fullfile(fpath_prism{4},'PC_Result.mat'),'Result');

%show_frame=find(PC_Result{PC_ind(1)}.Blue==0);
show_frame=find(Result_ST.VR(8,:)~=25 & Result_ST.VR(8,:)>=3);
StimLap=unique(Result_ST.VR(8,:).*(Result_ST.Blue>0));
Endlap=max(Result_ST.VR(8,:));
%LOI=setdiff([PF_list(show_PC,4):PF_list(show_PC,5)],StimLap);
LOI=setdiff([3:Endlap],[25]);
repLap=ringmovMean(Result_ST.Lap_FR(:,:,PC_ind(2))*1000,7);
Show_Lap=repmat(repLap(LOI,:),1,3);

if PF_list(show_PC,2)<PF_list(show_PC,3)
PFcenter=round(mean(PF_list(show_PC,2:3)));
else %end is after teleport
PFcenter=mean([PF_list(show_PC,2) PF_list(show_PC,3)+pos_bin]);
PFcenter=round(mod(PFcenter,pos_bin));
end
PFcenter_pr=30;
%show_frame_pr=[find(Result_Pr.VR(8,:)>=3,1)-10000:size(Result_Pr.normTraces,2)];
show_frame_pr=[1:size(Result_Pr.normTraces,2)];
LapFR_pr=PlaceTrigger_average(Result_Pr.spike(1,show_frame_pr),150,Result_Pr.VR(:,show_frame_pr),0,115)*1000;
repLap_pr=ringmovMean(LapFR_pr,7);
Show_Lap_pr=repmat(repLap_pr,1,3);

figure(3); clf; tiledlayout(8,2);
nexttile(1,[5 1])
tr=mean(Result_Pr.normTraces([1 2],show_frame_pr),1,'omitnan');
VRtr=Result_Pr.VR(:,show_frame_pr); VRtr(8,:)=VRtr(8,:)-2;
show_traces_align_Position(zscore(tr),Result_Pr.spike(1,show_frame_pr),PFcenter_pr/150*115,[-2500:2500],VRtr,1)
%show_traces_align_Position(zscore(PC_Result{PC_ind(1)}.normTraces(PC_ind(2),:)),PC_Result{PC_ind(1)}.spike(PC_ind(2),:),PFcenter/150*115,[-2000:2000],PC_Result{PC_ind(1)}.VR(:,:),1)

nexttile(2,[5 1])
show_traces_align_Position(zscore(Result_ST.normTraces(PC_ind(2),show_frame)),Result_ST.spike(PC_ind(2),show_frame),PFcenter/150*115,[-2500:2500],Result_ST.VR(:,show_frame),1)
%show_traces_align_Position(zscore(PC_Result{PC_ind(1)}.normTraces(PC_ind(2),:)),PC_Result{PC_ind(1)}.spike(PC_ind(2),:),PFcenter/150*115,[-2000:2000],PC_Result{PC_ind(1)}.VR(:,:),1)

ax1=[];
ax1=[ax1 nexttile(11,[2 1])];
imagesc(nVRbin/150*200,[1:size(LapFR_pr,1)],Show_Lap_pr(:,PFcenter_pr+pos_bin+nVRbin)); 
xlabel('Position (cm)');
ylabel('Laps')
cb=colorbar;
cb.Label.String='Firing rate (Hz)';
colormap('hot')

ax1=[ax1 nexttile(12,[2 1])];
imagesc(nVRbin/150*200,[1:length(LOI)],Show_Lap(:,PFcenter+pos_bin+nVRbin)); 
xlabel('Position (cm)');
ylabel('Laps')
cb=colorbar;
cb.Label.String='Firing rate (Hz)';
colormap('hot')

ax2=[];
ax2=[ax2 nexttile(15,[1 1])];
plot(nVRbin/150*200,mean(Show_Lap_pr(17:28,PFcenter_pr+pos_bin+nVRbin),1,'omitnan'))
xlabel('Position (cm)')
ylabel('Firing rate (Hz)')

ax2=[ax2 nexttile(16,[1 1])];
plot(nVRbin/150*200,mean(Show_Lap(12:24,PFcenter+pos_bin+nVRbin),1,'omitnan'))
xlabel('Position (cm)')
ylabel('Firing rate (Hz)')

linkaxes([ax1 ax2],'x')
axis tight

%% Place cell data
% Peak firing rate, mean firing rate, place field width, spatial
% information

nVRbin=[-50:50];
PF_MeanFiringRate=[];
PF_peakFiringRate=[];
PlaceField=[];
spatial_info=[];
% Prism
foi=find(cell2mat(cellfun(@(x) sum(isnan(x)),PlaceFieldList,'UniformOutput',false))==0);

spatial_info{1}=NaN(max(foi),1);
PF_peakFiringRate{1}=NaN(max(foi),1);
PF_MeanFiringRate{1}=NaN(max(foi),1);
for f=foi'
    load(fullfile(fpath_prism{f},'PC_Result.mat'),'Result');
    binTrack=(ceil(Result.VR(5,:)/((115)/150)));
    if ~isfield(Result,'LapFR')
        Result.LapFR=PlaceTrigger_average(Result.spike(1,:),150,Result.VR,-0.002,115)*1000;
    end
    for p=1:length(PlaceFieldList{f})/2
        Lapvec=(Result.VR(8,:)>PlaceFieldList{f}(2*(p-1)+1) & Result.VR(8,:)<PlaceFieldList{f}(2*(p-1)+2));
        if PlaceFieldBin{f}(2*(p-1)+1)>PlaceFieldBin{f}(2*(p-1)+2) % front is bigger than end
            Pvec=~(binTrack<(PlaceFieldBin{f}(2*(p-1)+1)) & binTrack>(PlaceFieldBin{f}(2*(p-1)+2)));
            PFcenter=mod(mean([PlaceFieldBin{f}(2*(p-1)+1) PlaceFieldBin{f}(2*(p-1)+2)+place_bin]),place_bin);
        else
            Pvec=(binTrack>(PlaceFieldBin{f}(2*(p-1)+1)) & binTrack<(PlaceFieldBin{f}(2*(p-1)+2)));
            PFcenter=mean([PlaceFieldBin{f}(2*(p-1)+1) PlaceFieldBin{f}(2*(p-1)+2)]);
        end
        PFvec=double(Lapvec & Pvec);
        PFvec(Result.motionReject)=NaN;
        SpikeTr=double(Result.spike(1,:));
        SpikeTr(Result.motionReject)=NaN;
        PF_MeanFiringRate{1}(f)=mean(SpikeTr(PFvec>0),'omitnan')*1000;

        LapFR_rep=repmat(Result.LapFR(PlaceFieldList{f}(2*(p-1)+1):PlaceFieldList{f}(2*(p-1)+2),:),1,3);
        PF_FR=movmean(mean(LapFR_rep(:,round(PFcenter)+nVRbin+place_bin),1,'omitnan'),5);
        PlaceField{1}{f}=[nVRbin' PF_FR'];
        
        PF_peakFiringRate{1}(f)=max(PF_FR);
        spatial_info{1}(f) = SpatialInfo(Result.spike(1,Lapvec), binTrack(Lapvec));
    end
end

% Cannula
            
spatial_info{2}=NaN(size(PF_list,1),1);
for p=1:size(PF_list,1)
    PC_ind=PC_list(PF_list(p,1),:); % file ,N th neuron

    binTrack=(ceil(PC_Result{PC_ind(1)}.VR(5,:)/((115)/150)));
    if ~isfield(PC_Result{PC_ind(1)},'Lap_FR')
        PC_Result{PC_ind(1)}.Lap_FR=PlaceTrigger_average(PC_Result{PC_ind(1)}.spike,150,PC_Result{PC_ind(1)}.VR,-0.002,115)*1000;
    end

      Lapvec=(PC_Result{PC_ind(1)}.VR(8,:)>PF_list(p,4) & PC_Result{PC_ind(1)}.VR(8,:)<PF_list(p,5));
        if PF_list(p,2)>PF_list(p,3) % front is bigger than end
            Pvec=~(binTrack<PF_list(p,2) & binTrack>PF_list(p,3));
            PFcenter=mod(mean([PF_list(p,2) PF_list(p,3)+place_bin]),place_bin);
        else
            Pvec=(binTrack>PF_list(p,2) & binTrack<PF_list(p,3));
            PFcenter=mean([PF_list(p,2) PF_list(p,3)]);
        end
        PFvec=double(Lapvec & Pvec);
        SpikeTr=double(PC_Result{PC_ind(1)}.spike(PC_ind(2),:));
        PF_MeanFiringRate{2}(p)=mean(SpikeTr(PFvec>0),'omitnan')*1000;

        BlueLap=unique(PC_Result{PC_ind(1)}.VR(8,:).*(PC_Result{PC_ind(1)}.Blue>0));
        Lap_list=setdiff([PF_list(p,4):PF_list(p,5)],BlueLap);
        LapFR_rep=repmat(PC_Result{PC_ind(1)}.Lap_FR(Lap_list,:,PC_ind(2)),1,3);
        PF_FR=movmean(mean(LapFR_rep(:,round(PFcenter)+nVRbin+place_bin),1,'omitnan'),5);
        PlaceField{2}{p}=[nVRbin' PF_FR'];
        
        PF_peakFiringRate{2}(p)=max(PF_FR);
        spatial_info{2}(p) = SpatialInfo(PC_Result{PC_ind(1)}.spike(PC_ind(2),Lapvec), binTrack(Lapvec));
end

figure(4); clf;
cmap=[0. 0.3 1; 1 0.3 0]; 
cmap_line=[0.7 0.8 1; 1 0.8 0.7]; h=[];
for c=1:2
    if c==1
PFs=cell2mat(cellfun(@(x) x(:,2),PlaceField{c}(foi),'UniformOutput',false));
    else
PFs=cell2mat(cellfun(@(x) x(:,2),PlaceField{c},'UniformOutput',false));        
    end
M=mean(PFs,2,'omitnan'); S=std(PFs,0,2,'omitnan')./sqrt(size(PFs,2));
plot(nVRbin/150*200,PFs,'color',cmap_line(c,:));
hold all
h(c)=errorbar_shade(nVRbin/150*200,M',S',cmap(c,:));
end
legend(h,{'Prism','Cannula'})
xlabel('VR position (cm)')
ylabel('Firing rate (Hz)');
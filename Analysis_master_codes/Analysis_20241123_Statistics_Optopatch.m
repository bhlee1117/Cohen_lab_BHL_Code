
clear
clc;
cd '/Users/bhlee1117/Documents/BHL/Matlab_project/20241020_OptogeneticAnalysis';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'B5:O175');

fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,5));
CamType=raw(:,3);
StructureData=raw(:,10);
StimROI=raw(:,6);
StimWfn=raw(:,7);
PixelSize=cell2mat(cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false));
refROI=cellfun(@(x) (str2num(num2str(x))),raw(:,14),'UniformOutput',false);

place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
title_str={'Basal','Apical','Peri-Soma'};
set(0,'DefaultFigureWindowStyle','docked')

%%
foi=89;
for f=foi
    if contains(fpath{f},'BHL18TB_D2')
        targetFile = strrep(fpath{f}, 'BHL18TB_D2','cohen_lab/Lab/Labmembers/Byung Hun Lee/Data');
        load(fullfile(targetFile,'OP_Result.mat'))
    else
        load(fullfile(fpath{f},'OP_Result.mat'))
    end

    figure(10); clf; tiledlayout(2,1)
    nexttile([1 1])
    show_footprnt_contour(Result.ftprnt,Result.ref_im)
    nexttile([1 1])
    nROI=size(Result.traces,1);
    plot(rescale2(Result.normTraces,2)'+[1:nROI]);
end

%%
for f=46:49;
load(fullfile(fpath{f},'OP_Result.mat'),'Result');
    interDendDist=[];
    SkelDend = Skeletonize_dendrite(Result.ref_im,6,0.02,10);
    nROI=size(Result.ftprnt,3);
    for i=1%:nROI
        i
        for j=1:nROI
            [interDendDist(i,j), path]=geodesic_distance(SkelDend,get_coord(Result.ftprnt(:,:,i)),get_coord(Result.ftprnt(:,:,j)));
        end
    end
    Result.interDendDist=interDendDist;
    save(fullfile(fpath{f},'OP_Result.mat'),'Result','-v7.3');
end
%%
[~, unqInd] = unique([Mouse NeuronInd] ,'row');
StimROI_Ind={'Soma','Distal Dend'};
StimWfn_Ind={'Ramp Stim','Short Pulse'};
CatResult=[];
for i=unqInd(25);
    SameCellInd=find(Mouse==Mouse(i) & NeuronInd==NeuronInd(i));
    isSoma = ~cellfun(@isempty,strfind(StimROI(SameCellInd), StimROI_Ind{1}));
    isRamp = ~cellfun(@isempty,strfind(StimWfn(SameCellInd), StimWfn_Ind{1}));

    ROIwvf_ind=(double(isSoma)*1+double(isRamp)*2)+1; % 1 = Dend, SP; 2 = Soma, SP; 3 = Dend, RP; 4 = Soma, RP;
    g=ones(1,4);
    for j=1:length(SameCellInd)
        f2read=SameCellInd(j);
        load(fullfile(fpath{f2read},'OP_Result.mat'),'Result');
        CatResult{ROIwvf_ind(j),g(ROIwvf_ind(j))}=Result;
        g(ROIwvf_ind(j))=g(ROIwvf_ind(j))+1;
    end

    nROI=size(Result.ftprnt,3);

CS_STAmat=cellfun(@(x) reshape(x.normTraces(:,find(x.SpClass(2,:)==1)'+[-100:100]),nROI,[],201),CatResult(find(cellfun(@isempty,CatResult)==0)),'UniformOutput',false);
CS_STAmat=cat(2,CS_STAmat{:});
CS_STAmat=squeeze(mean(CS_STAmat,2));
F_ref=mean(CS_STAmat(:,100+[20:30]),2,'omitnan')';

end
%%
i=46; bound=6;
load(fullfile(fpath{i},'OP_Result.mat'),'Result');
load([fpath{i} '/mcTrace' num2str(1,'%02d') '.mat']);
mov_res=readBinMov_BHL(fpath{i});
mean_F=squeeze(mean(mov_res(bound:end-bound,bound:end-bound,:),[1 2]));

mov_res= mov_res-mean(mov_res,3);
bkg = zeros(1, size(mov_res,3));

[~, blueOff]=get_blueoffTrace(mean_F,[Result.Blue],70);
[y_fit]=expfitDM_2(find(blueOff)',mean_F(find(blueOff)),[1:size(mov_res,3)]',1000);
bkg(1,:)=y_fit;
mov_res = SeeResiduals(mov_res,Result.mc);
mov_res = SeeResiduals(mov_res,Result.mc.^2);
mov_res = SeeResiduals(mov_res,Result.mc(:,1).*Result.mc(:,end));
mov_res= SeeResiduals(mov_res,bkg,1);
CS_s=find(Result.SpClass(2,:)==1);
CS_STAmov=reshape(mov_res(:,:,CS_s'+[-100:100]),size(mov_res,1),size(mov_res,2),[],201);
CS_STAmov=squeeze(mean(CS_STAmov,3));
CS_STAmov=-(CS_STAmov-mean(CS_STAmov(:,:,1:30),3));
CS_STAmov=CS_STAmov.*double(max(Result.bvMask,[],3)==0);
STAmovVec=tovec(CS_STAmov);

%%

% CS_avg=Result.normTraces(:,CS_s'+[-100:100]);
% CS_avg=squeeze(mean(reshape(Result.normTraces(:,CS_s'+[-100:100]),size(Result.ftprnt,3),[],201),2));
% CS_avg=CS_avg-mean(CS_avg(:,1:30),2);

coord_1d=dim_reduce(get_coord(CatResult{2}.ftprnt));
coord_1d=coord_1d-coord_1d(1);
[~, dsort]=sort(coord_1d,'ascend');

% nROI=size(Result.ftprnt,3);
% maxfrm=[]; F0slope=[]; Fslope=[]; Fslope_weight=[];
% F0=tovec(Result.ref_im)-100;
% F0_filter=tovec(imgaussfilt(Result.ref_im,2))-100;
%F_ref=mean(STAmovVec(:,100+[20:30]),2,'omitnan')'*tovec(Result.ftprnt);
%F_ref=mean(CS_avg(:,100+[20:30]),2,'omitnan')';


figure(13); clf;
for n=1:nROI
    
px=find(tovec(Result.ftprnt(:,:,n)>0) & tovec(max(Result.bvMask,[],3)==0));
[~, maxfrm(n)]=max(mean(STAmovVec(px,1:105),1,'omitnan'));

dF=STAmovVec(:,maxfrm(n));
%dF=max(STAmovVec,[],2);
F0_weight=F0.*rescale(tovec(Result.ftprnt(:,:,n)));
dF_weight=dF.*rescale(tovec(Result.ftprnt(:,:,n)));

px2=px(find(F0(px)>prctile(F0(px),30) & F0(px)<prctile(F0(px),98)));
px2_weight=px(find(F0_weight(px)>prctile(F0_weight(px),30) & F0_weight(px)<prctile(F0_weight(px),98)));

nexttile([1 1])
plot(F0(px),dF(px),'.'); hold all
%plot(F0_weight(px),dF_weight(px),'.'); hold all
[p]=polyfit(F0(px2), dF(px2), 1); hold all
[p_weight]=polyfit(F0_weight(px2_weight), dF_weight(px2_weight), 1); hold all
% Get fitted values
y_fit = polyval(p, F0(px2));
plot(F0(px2),y_fit,'r')
SS_res = sum((dF(px2) - y_fit).^2);
SS_tot = sum((dF(px2) - mean(dF(px2))).^2);
R_squared = 1 - (SS_res / SS_tot);

y_fit_weight = polyval(p, F0_weight(px2_weight));
SS_res = sum((dF_weight(px2_weight) - y_fit_weight).^2);
SS_tot = sum((dF_weight(px2_weight) - mean(dF(px2_weight))).^2);
R_squared_weight = 1 - (SS_res / SS_tot);

title(['d(\DeltaF)/dF0 : ' num2str(p(1),2) ', R^2 : ' num2str(R_squared,2)])
Fslope(n)=p(1);
Fslope_weight(n,1)=p_weight(1);
Rsq(n)=R_squared;
F0slope(n,1)=CS_avg(n,maxfrm(n))/Fslope(n);
F0slope_weight(n,1)=CS_avg(n,maxfrm(n))/Fslope_weight(n);
end

figure(14); clf;
nexttile([1 1])
l=plot(CS_avg(dsort,:)'./F0slope_weight(dsort)');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2))
xlabel('Time (ms)'); ylabel('\DeltaF/F (robust)');

nexttile([1 1])
imagesc(CS_avg(dsort,:)./F0slope_weight(dsort));
colormap(turbo);
xlabel('Time (ms)'); ylabel('Basal to apical ROI');

nexttile([1 1])
Dsign=ones(1,size(Result.interDendDist,2));
Dsign(dsort(1:find(dsort==1)-1))=-1;
dendaxis=Result.interDendDist(1,:).*Dsign;

plot(dendaxis(dsort),max(CS_avg(dsort,100:104)./F0slope_weight(dsort),[],2),'.-'); hold all
plot(dendaxis([19 1:8]),max(CS_avg([19 1:8],100:104)./F0slope_weight([19 1:8]),[],2),'ro')
xlabel('Distance from soma (\mum)'); ylabel('\DeltaF/F (robust)');

nexttile([1 1])
l=plot(CS_avg(dsort,:)'./F_ref(dsort));
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2))
xlabel('Time (ms)'); ylabel('Spike Amplitude (\DeltaF/F_{ref})');

nexttile([1 1])
imagesc(CS_avg(dsort,:)./F_ref(dsort)');
colormap(turbo);
xlabel('Time (ms)'); ylabel('Basal to apical ROI');

nexttile([1 1])
Dsign=ones(1,size(Result.interDendDist,2));
Dsign(dsort(1:find(dsort==1)-1))=-1;
dendaxis=Result.interDendDist(1,:).*Dsign;

plot(dendaxis(dsort),max(CS_avg(dsort,100:104)./F_ref(dsort)',[],2),'.-'); hold all
plot(dendaxis([19 1:8]),max(CS_avg([19 1:8],100:104)./F_ref([19 1:8])',[],2),'ro')
xlabel('Distance from soma (\mum)'); ylabel('Spike Amplitude (\DeltaF/F_{ref})');
%%

figure(15); clf; tiledlayout(5,1); ax1=[];
ax1=[ax1 nexttile([1 1])];
NormTrace=CatResult{4,1}.normTraces./F_ref';
NormTrace=NormTrace-prctile(NormTrace,10,2);
plot(NormTrace([1],:)')
title('Soma Stimulation')
ylabel('\DeltaF/F_{ref}')

ax1=[ax1 nexttile([1 1])];
imagesc(NormTrace(dsort,:),[-0.5 2.5]);
if find(dsort==1)~=1
    set(gca,'YTick',[1 find(dsort==1) nROI],'YTickLabel',num2str([min(interDendDist(1,:).*Dsign) 0 max((interDendDist(1,:).*Dsign))]',3))
else
    set(gca,'YTick',[1 nROI],'YTickLabel',num2str([0 max((interDendDist(1,:)))]',3))
end

ax1=[ax1 nexttile([1 1])];
NormTrace=CatResult{3,2}.normTraces./F_ref';
NormTrace=NormTrace-prctile(NormTrace,10,2);
plot(NormTrace([1],:)')
title('Distal Dendrite Stimulation')
ylabel('\DeltaF/F_{ref}')

ax1=[ax1 nexttile([1 1])];
imagesc(NormTrace(dsort,:),[-0.5 2.5]); colormap(turbo);
if find(dsort==1)~=1
    set(gca,'YTick',[1 find(dsort==1) nROI],'YTickLabel',num2str([min(interDendDist(1,:).*Dsign) 0 max((interDendDist(1,:).*Dsign))]',3))
else
    set(gca,'YTick',[1 nROI],'YTickLabel',num2str([0 max((interDendDist(1,:)))]',3))
end

ax1=[ax1 nexttile([1 1])];
plot(CatResult{4,1}.Blue)
xlabel('Time (ms)')
ylabel('Time (ms)')
linkaxes(ax1,'x')
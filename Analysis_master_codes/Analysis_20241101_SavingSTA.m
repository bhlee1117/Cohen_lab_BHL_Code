clear
clc;
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:S31');

ref_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,10),'UniformOutput',false);
basal_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,11),'UniformOutput',false);
apical_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false);
fpath=raw(:,1)';
StructureData=raw(:,10);
BadROI=cellfun(@(x) (str2num(num2str(x))),raw(:,15),'UniformOutput',false);
EndFrame=cell2mat(raw(:,13));
ifmotionReject=cell2mat(raw(:,14));
ifdirtRemov=cell2mat(raw(:,16));
Pixelsize=cell2mat(raw(:,6));
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
title_str={'Basal','Apical','Peri-Soma'};
set(0,'DefaultFigureWindowStyle','docked')
foi=[1 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27];

%% Save the STA movies
nTau={[-30:20],[-70:100],[-30:20]}; %SS, CS, dSP
bound=6;
%f_tmp='/Volumes/BHL18TB_D1/20240218/134705BHLm117_FOV2_VR2';
for f=[27]
    f
    load(fullfile(fpath{f},'PC_Result.mat'),'Result')
    StackedSpike=[];
    alignmovlist=dir(fullfile(fpath{f},'STA_Mat*.bin'));
    alignmov_ind=ones(1,3);
    for l=1:length(alignmovlist)
        delete(fullfile(fpath{f},alignmovlist(l).name));
    end

    for c=1:3 % Get frames of spikes of each classes
        fileInd{c}=ceil(find(Result.SpClass(c,:))/time_segment);
        frameInd{c}=mod(find(Result.SpClass(c,:)),time_segment);
        STA_movie_align{c}=[];
    end
    %line up movie segments at simple spike

    load([fpath{f} '/output_data.mat'])
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    frm_end=EndFrame(f);
    f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
    g=zeros(1,3);
    ss=[0 0 0];
    filetoOpen=unique(cell2mat(fileInd));
    for j=filetoOpen
        j
        mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
        load([fpath{f} '/mcTrace' num2str(j,'%02d') '.mat']);

        if j==length(f_seg)-1
            mc=mcTrace.xymean;
        else
            mov_mc=mov_mc(:,:,1:time_segment);
            mc=mcTrace.xymean(1:time_segment,:);
        end

        mov_res= mov_mc-mean(mov_mc,3);
        bkg = zeros(1, size(mov_mc,3));
        bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
        mov_res = SeeResiduals(mov_res,mc);
        mov_res = SeeResiduals(mov_res,mc.^2);
        mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,end));
        mov_res= SeeResiduals(mov_res,bkg,1);
        
        if sum(isnan(mov_res(:)))>0
            mov_res= mov_mc-mean(mov_mc,3);
        end

        for c=1:3
            for k=find(fileInd{c}==j)
                if frameInd{c}(k)+nTau{c}(1)>1 && frameInd{c}(k)+nTau{c}(end)<size(mov_res,3)
                    mov_seg=mov_res(:,:,frameInd{c}(k)+nTau{c});
                    %mov_seg=mat2gray(mov_seg);
                    STA_movie_align{c} = cat(3,STA_movie_align{c}, mov_seg);
                    g(c)=g(c)+1;
                    StackedSpike{c}(1,g(c))=1;
                    StackedSpike{c}(2,g(c))=time_segment*(j-1)+frameInd{c}(k);
                else
                    STA_movie_align{c} = cat(3,STA_movie_align{c}, zeros(sz(2),sz(1),length(nTau{c})));
                    g(c)=g(c)+1;
                    StackedSpike{c}(1,g(c))=0;
                    StackedSpike{c}(2,g(c))=NaN;
                end
            end

            if ~isempty(STA_movie_align{c})

                % write_tif_stack(STA_movie_align{c},fullfile(f_tmp,[alignedMovFN{c} '_' num2str(alignmov_ind(c)) '.tiff']))
                % alignmovlist=dir(fullfile(f_tmp,[alignedMovFN{c} '*.tiff']));
                %write_tif_stack(STA_movie_align{c},fullfile(fpath{f},[alignedMovFN{c} '_' num2str(alignmov_ind(c)) '.tiff']))
                MovtoWrite=vm(STA_movie_align{c}-min(STA_movie_align{c}(:)));
                Movinfo=whos('MovtoWrite');
                if Movinfo.bytes > 6.0*10^9 
                    MovtoWrite.transpose.savebin(fullfile(fpath{f},[alignedMovFN{c} '_' num2str(alignmov_ind(c)) '.bin']))
                    disp('Move on to the next bin')
                    alignmov_ind(c) = alignmov_ind(c) + 1;
                    STA_movie_align{c}=[];
                else if j==filetoOpen(end)
                    MovtoWrite.transpose.savebin(fullfile(fpath{f},[alignedMovFN{c} '_' num2str(alignmov_ind(c)) '.bin']))
                    disp(['Saving ' alignedMovFN{c} ' of this recording is finished.'])    
                end
                end
            end
            pause(10);
        end
        g
    end

    Result.StackedSpike=StackedSpike;
    save(fullfile(fpath{f},'PC_Result.mat'),'Result','-v7.3')
end

%%
nTau_bAP=[-25:15];
nTau=[-70:50];
for f=[27]
f
    load(fullfile(fpath{f},'PC_Result.mat'),'Result')
 load([fpath{f} '/output_data.mat'])
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    
alignmovlist=dir(fullfile(fpath{f},[alignedMovFN{1} '*.bin']));
AlignMov=[];
nReadmov=min([length(alignmovlist) 4]);
for l=1:nReadmov
AlignMov=cat(3,AlignMov,readBinMov(fullfile(fpath{f},alignmovlist(l).name),sz(2),sz(1)));
%AlignMov=cat(3,AlignMov,readtiff(fullfile(fpath{f},alignmovlist(l).name)));
end

if length(alignmovlist)>4
AlignMovrsh=double(reshape(AlignMov,size(AlignMov,1),size(AlignMov,2),51,[]));
else
    if isfield(Result,'StackedSpike')
AlignMovrsh=double(reshape(AlignMov,size(AlignMov,1),size(AlignMov,2),[],size(Result.StackedSpike{1},2)));
    else
AlignMovrsh=double(reshape(AlignMov,size(AlignMov,1),size(AlignMov,2),[],sum(Result.SpClass(1,:))));        
    end
end
%AlignMovrsh=AlignMovrsh-mean(AlignMovrsh(:,:,1:15,:),3);

clear AlignMov

alignmovlist=dir(fullfile(fpath{f},[alignedMovFN{2} '*.bin']));
AlignMovCS=[];
nReadmov=min([length(alignmovlist) 2]);
for l=1:nReadmov
AlignMovCS=cat(3,AlignMovCS,readBinMov(fullfile(fpath{f},alignmovlist(l).name),sz(2),sz(1)));
%AlignMovCS=cat(3,AlignMovCS,readtiff(fullfile(fpath{f},alignmovlist(l).name)));
end

if length(alignmovlist)>4
AlignMovrshCS=double(reshape(AlignMovCS,size(AlignMovCS,1),size(AlignMovCS,2),171,[]));
else
    if isfield(Result,'StackedSpike')
AlignMovrshCS=double(reshape(AlignMovCS,size(AlignMovCS,1),size(AlignMovCS,2),[],size(Result.StackedSpike{2},2)));
    else
AlignMovrshCS=double(reshape(AlignMovCS,size(AlignMovCS,1),size(AlignMovCS,2),[],sum(Result.SpClass(2,:))));        
    end
end

som_spike=find(Result.spike(1,:)); 
bAP_ref=[];
for s=som_spike
    isnearby=sum(ismember(s+nTau_bAP,som_spike))>1;
    isnearbyCS=sum(ismember(s+nTau_bAP,find(Result.CStrace)))>1;
    if ~isnearby & ~isnan(s) & ~isnearbyCS
        bAP_ref=[bAP_ref s];
    end
end
sp_na=sum((bAP_ref'+nTau_bAP)<0 | (bAP_ref'+nTau_bAP)>size(Result.traces,2),2)==0;
bAP_ref=bAP_ref(sp_na);

% Isolated Complex spike
    CS_s=[];
    C_list=find(Result.SpClass(2,:));
    CS_label=max(Result.spike,[],1).*bwlabel(Result.CStrace);
    CS_list=[];
    for g=1:length(C_list)
        s=C_list(g);
        if sum((s+nTau)<0 | (s+nTau)>size(Result.traces,2),2)==0
            s_tmp=Result.spike(1,:);
            s_tmp(find(CS_label==g))=0;
            isnearbyCS=max(bwlabel(Result.CStrace(s+nTau)))>1;
            isnearbyS=sum(ismember(s+nTau,find(s_tmp)))>0;
            if ~isnearbyCS & ~isnearbyS
                CS_s=[CS_s s];
                CS_list=[CS_list g];
            end
        end
    end
    if ~isempty(CS_s)
        sp_na=sum((CS_s'+nTau)<0 | (CS_s'+nTau)>size(Result.traces,2),2)==0;
        CS_s=CS_s(sp_na);
    end


if isfield(Result,'StackedSpike')
Refindex=ismember(Result.StackedSpike{1,1}(2,1:size(AlignMovrsh,4)),bAP_ref);
RefindexCS=ismember(Result.StackedSpike{1,2}(2,1:size(AlignMovrshCS,4)),CS_s);
else
Refindex=ismember(find(Result.SpClass(1,:)) ,bAP_ref);    
Refindex=Refindex(1:size(AlignMovrsh,4));

RefindexCS=ismember(find(Result.SpClass(2,:)) ,CS_s);    
RefindexCS=RefindexCS(1:size(AlignMovrshCS,4));
end

STAmov=mean(-AlignMovrsh(:,:,:,find(Refindex)),4);
STAmovvm=vm(STAmov); 
STAmovvm.transpose.savebin(fullfile(fpath{f},['STAmovieSS' '.bin']))
STAmov=STAmov-mean(STAmov(:,:,1:15),3);

STAmovCS=mean(-AlignMovrshCS(:,:,:,find(RefindexCS)),4);
STAmovCSvm=vm(STAmovCS); 
STAmovCSvm.transpose.savebin(fullfile(fpath{f},['STAmovieCS' '.bin']))
STAmovCS=STAmovCS-mean(STAmovCS(:,:,1:15),3);
%[~, maxfrm]=max(squeeze(mean(STAmov,[1 2])));
%STAmovVec=tovec(mean(STAmov(:,:,maxfrm+[-1:1]),3,'omitnan'));
STAmovVec=tovec(STAmov);

Fslope=[]; Rsq=[]; bound=6;
figure(f); clf;
F0=tovec(Result.ref_im);
for n=1:size(Result.ftprnt,3)
    
px=find(tovec(Result.ftprnt(:,:,n)>0));
[~, maxfrm]=max(mean(STAmovVec(px,:),1,'omitnan'));

dF=STAmovVec(:,maxfrm);
F0_weight=F0.*rescale(tovec(Result.ftprnt(:,:,n)));
dF_weight=dF.*rescale(tovec(Result.ftprnt(:,:,n)));

px2=px(find(F0_weight(px)>prctile(F0_weight(px),30) & F0_weight(px)<prctile(F0_weight(px),98)));

nexttile([1 1])
plot(F0(px),dF(px),'.'); hold all
plot(F0_weight(px),dF_weight(px),'.'); hold all
[p]=polyfit(F0_weight(px2), dF_weight(px2), 1); hold all
% Get fitted values
y_fit = polyval(p, F0_weight(px2));
plot(F0_weight(px2),y_fit,'r')
SS_res = sum((dF_weight(px2) - y_fit).^2);
SS_tot = sum((dF_weight(px2) - mean(dF_weight(px2))).^2);
R_squared = 1 - (SS_res / SS_tot);

title(['d(\DeltaF)/dF0 : ' num2str(p(1),2) ', R^2 : ' num2str(R_squared,2)])
Fslope(n)=p(1);
Rsq(n)=R_squared;
end
Result.dFF_slope=Fslope;
Result.dFF_slopefit=Rsq;
save(fullfile(fpath{f},'PC_Result.mat'),'Result','-v7.3');
fprintf('Variable "Result" is saved at: %s \n', fpath{f});
pause(1);
end

%%
clf;
plot(F_ref,Fslope,'.')
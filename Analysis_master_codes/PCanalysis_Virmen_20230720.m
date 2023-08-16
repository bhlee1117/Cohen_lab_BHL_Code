%% Load data
clear
cd /Volumes/BHL_WD18TB/20230717
load('pcResult_20230716.mat')
%%
for i=1:length(fpath)
    Result{i}.spike=zeros(size(Result{i}.traces));
%     tmp=squeeze(Result{i}.traces_res) - movmedian(squeeze(Result{i}.traces_res),10,2); tmp=tmp./get_threshold(tmp,1);
%     Result{i}.spike=find_spike_bh(tmp,4,3);
    tmp=squeeze(Result{i}.traces_res) - movmedian(squeeze(Result{i}.traces_res),250,2); tmp=tmp./get_threshold(tmp,1);
    Result{i}.spike=(Result{i}.spike+find_spike_bh(tmp,5.5,2))>0;
end

%%
place_bin=40;
grps={[1 2 3 4]};
im_corr_th=[3 3 3 3];
cmap = jet(256);
cmap(1:find(cmap(:,1)==0.5),1) = 0.5;
clear ind_list Virmen

for i=1%:length(grps)
  
    goi=grps{i}(1);

    [sz1 sz2]=size(Result{goi}.ref_im);
    match_list=zeros(size(Result{goi}.centers,1),2,length(grps{i})-1);
    Virmen{1}=Result{goi}.Virmen(:,1:end-2);
    dist_thres=5;
    for g=2:length(grps{i})
        goi=grps{i}(g);
        dist=distance_mat(Result{grps{i}(1)}.centers,Result{goi}.centers);
        [d arg]=min((dist),[],2);
        match_list(find((d)<dist_thres),1,g-1)=[1:sum(d<dist_thres)]';
        match_list(arg(find((d)<dist_thres)),2,g-1)=[1:sum(d<dist_thres)]';
        Virmen{g}=Result{goi}.Virmen(:,1:end-2);
    end

    ind_list(:,1)=find(sum(match_list(:,1,:)==0,3)==0);

    for c=1:length(ind_list(:,1))
        for g=2:length(grps{i})
            c_tmp=find(match_list(:,2,g-1)==match_list(ind_list(c,1),1,g-1));
            if isempty(c_tmp)
                ind_list(c,g)=NaN;
            else
                ind_list(c,g)=c_tmp;
            end
        end
    end

    traces_sp=[]; traces_F=[];
    for g=1:length(grps{i})
        goi=grps{i}(g);
        trace_length=size(Result{goi}.traces_res,2);
        sp=double(Result{goi}.spike(ind_list(:,g),:)); sp(:,end:size(Result{goi}.Virmen,1))=0;
        F=Result{goi}.traces(ind_list(:,g),:);  F(:,end:size(Result{goi}.Virmen,1))=0;
        %F=F-movmedian(F,400,2); 
        F=F./get_threshold(F,1);

        clear im_corr

        for n=1:size(Result{goi}.centers,1)
        im_corr(n,:)=movmean(Result{goi}.im_corr{n}(1:trace_length),20); 
        end
        im_corr=abs(zscore(im_corr,0,2));
        F(im_corr>im_corr_th(g))=NaN; sp(im_corr>im_corr_th(g))=NaN;
        
        traces_sp=[traces_sp sp(:,1:end-2)]; traces_F=[traces_F F(:,1:end-2)];
    end
%traces_F=traces_F-movmedian(traces_F,200,2);    
Virmen_align=cell2mat(Virmen)';
rmv_frame=find(sum(isnan(Virmen_align),2)>0);
traces_sp(:,rmv_frame)=[]; traces_F(:,rmv_frame)=[];
Virmen_align(rmv_frame,:)=[];
Lick_trace=zeros(size(Virmen_align,1),1);
Lick_trace(find(rescale(Virmen_align(:,9))>0.4 & rescale(Virmen_align(:,9)) <0.5))=1;
figure(1); clf;
    show_traces_spikes(traces_F,traces_sp,Virmen_align(:,5))
    figure(2); clf;
    imshow2(Result{1}.FOV,[])
    hold all
    plot(Result{1}.centers(ind_list(:,1),1),Result{1}.centers(ind_list(:,1),2),'ro')
    text(Result{1}.centers(ind_list(:,1),1)+2,Result{1}.centers(ind_list(:,1),2) ...
        ,num2str([1:size(ind_list,1)]'),'color','w','FontSize',10)

%     figure(2); clf;
%     show_footprnt(Result{goi}.c_ftprnt,Result{goi}.ref_im);
    %show_traces_spikes(traces_F([18 19 22],:),traces_sp([18 19 22],:),Virmen_align(:,5))


    % for r=1:size(Result{i}.rm_ind,1)
    % traces_F(:,[Result{i}.rm_ind(r,1):Result{i}.rm_ind(r,2)])=NaN;
    % traces_sp(:,[Result{i}.rm_ind(r,1):Result{i}.rm_ind(r,2)])=NaN;
    % end
% 
traces_sub=get_subthreshold(traces_F,traces_sp,2);
traces_sub=movmedian(traces_F,50,2);


lap_length=115;
 Lap_FR=get_LapFR_Virmen(movsum(traces_sp,100,2,'omitnan'),place_bin,lap_length,Virmen_align,0.005);
 Lap_sub=get_LapFR_Virmen(movmean(traces_sub,100,2,'omitnan'),place_bin,lap_length,Virmen_align,0.005);
 Lap_blue=get_LapFR_Virmen(Virmen_align(:,12)'>0,place_bin,lap_length,Virmen_align,0);
% [trace_place SI boot_SI]=traces2place_Virmen(traces_sp,place_bin,Virmen_align,0.005,100);
% Pct=prctile(boot_SI',99.9)'; put_pc=find((SI-Pct)>0);
% trace_place=movmean(repmat(trace_place,1,3),6,2);
% trace_place=trace_place(:,place_bin+1:2*place_bin);
% trace_place_norm=(trace_place-min(trace_place,[],2));
% trace_place_norm=trace_place_norm./(max(trace_place,[],2)-min(trace_place,[],2));
% [~, arg]=max(trace_place_norm,[],2);
% [~, order]=sort(arg,'ascend');
% figure; tiledlayout(1,2); nexttile([1 1]);
% imagesc(trace_place_norm(order,:))
% set(gca,'ytick',[1:size(traces_F,1)],'yticklabel',num2str([order]))
% colormap(turbo)
% 
% trace_place_PC=movmean(repmat(trace_place(put_pc,:),1,3),6,2);
% trace_place_PC=trace_place_PC(:,place_bin+1:2*place_bin);
% trace_place_PC_norm=(trace_place_PC-min(trace_place_PC,[],2));
% trace_place_PC_norm=trace_place_PC_norm./(max(trace_place_PC,[],2)-min(trace_place_PC,[],2));
% [~, arg]=max(trace_place_PC_norm,[],2);
% %arg=sum(trace_place_PC_norm.*[1:place_bin],2,'omitnan')./sum(trace_place_PC_norm,2,'omitnan');
% [~, order_PC]=sort(arg,'ascend');
% nexttile([1 1]);
% imagesc(trace_place_PC_norm(order_PC,:))
% set(gca,'ytick',[1:size(trace_place_PC_norm,1)],'yticklabel',num2str([put_pc(order_PC)]))
% colormap(turbo)

end
figure;
moviesc(Lap_FR)
%%
imaged(im_merge(cat(3,Lap_FR(:,:,9),Lap_sub(:,:,9),Lap_blue),[1 0 0.5; 0 0.2 1]),[])
show_traces_align_Position_Virmen(traces_F,traces_sp,115*0.7,[-1000:3000],Virmen_align,2,Virmen_align(:,12)'>0)
show_traces_align_Position_Virmen(traces_sub,traces_sp,115*0.3,[-1000:3000],Virmen_align,9,Virmen_align(:,12)'>0)
show_traces_place_Virmen(traces_F,traces_sp,4,Virmen_align,[1:4],22)

%Lap_FR=show_traces_place(traces_F,traces_sp,place_bin,Ard_data,[30:36],[n]);

%%
for n=[2]
    [~, pf]=max(trace_place_norm(n,:));
    Lap_FR=show_traces_place(traces_F,traces_sp,place_bin,Ard_data,[30:36],[n]);
    nexttile([1 1])
    imagesc(trace_place_norm(n,:))
    colormap('jet')

    figure;
    imagesc((Lap_FR),[10 50])
end

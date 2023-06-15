%% Load data
clear
cd /Volumes/ByungHun4TB/Cohen_lab_Data/20230603
load('20230607_Result.mat')

%%
place_bin=50;
grps={[3 4 5]};
im_corr_th=[0.992 0.992 0.992];
cmap = jet(256);
cmap(1:find(cmap(:,1)==0.5),1) = 0.5;
clear ind_list Virmen

for i=1%:length(grps)
    %Result{i}.rm_ind=manual_trace_deletion(Result{i});
    goi=grps{i}(1);
    %Result{i}.rm_ind=[];
    %Result{i}.rm_ind=[Result{i}.rm_ind; manual_trace_deletion(Result{i})];
    [sz1 sz2]=size(Result{goi}.ref_im);
    match_list=zeros(size(Result{goi}.coord,1),2,length(grps{i})-1);
    Virmen{1}=Result{goi}.Virmen(:,1:end-2)';

    for g=2:length(grps{i})
        goi=grps{i}(g);
        xcorrRefIm=normxcorr2(Result{grps{i}(1)}.ref_im,Result{goi}.ref_im);
        [~, ind]=max(xcorrRefIm(:));
        [shifty shiftx]=ind2sub(size(xcorrRefIm,1),ind);
        shifty=sz1-shifty; shiftx=sz2-shiftx;
        dist=distance_mat(Result{grps{i}(1)}.coord,Result{goi}.coord+[shiftx shifty]);
        [d arg]=min((dist),[],2);
        match_list(find((d)<1.5),1,g-1)=[1:sum(d<1.5)]';
        match_list(arg(find((d)<1.5)),2,g-1)=[1:sum(d<1.5)]';
        Virmen{g}=Result{goi}.Virmen(:,1:end-2)';
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
        sp=double(Result{goi}.spike(ind_list(:,g),:)); sp(:,end:size(Result{goi}.Virmen,1))=0;
        F=Result{goi}.traces_res(ind_list(:,g),:);  F(:,end:size(Result{goi}.Virmen,1))=0;
        %F=F-movmedian(F,200,2); 
        F=F./get_threshold(F,1);
        im_corr=movmean(Result{goi}.im_corr,20); 
        F(:,im_corr<im_corr_th(g))=NaN; sp(:,im_corr<im_corr_th(g))=NaN;
        traces_sp=[traces_sp sp(:,1:end-2)]; traces_F=[traces_F F(:,1:end-2)];
    end
Virmen_align=cell2mat(Virmen');
rmv_frame=find(sum(isnan(Virmen_align),2)>0);
traces_sp(:,rmv_frame)=[]; traces_F(:,rmv_frame)=[];
Virmen_align(rmv_frame,:)=[];

    show_traces_spikes(traces_F,traces_sp,[Virmen_align(:,3)])

    % for r=1:size(Result{i}.rm_ind,1)
    % traces_F(:,[Result{i}.rm_ind(r,1):Result{i}.rm_ind(r,2)])=NaN;
    % traces_sp(:,[Result{i}.rm_ind(r,1):Result{i}.rm_ind(r,2)])=NaN;
    % end

Lap_FR=get_LapFR_Virmen(traces_sp,place_bin,Virmen_align,-1,[1:size(traces_sp)]);
[trace_place SI boot_SI]=traces2place(traces_sp,place_bin,Virmen_align,0,100);
Pct=prctile(boot_SI',99)'; put_pc=find((SI-Pct)>0);
trace_place=movmean(repmat(trace_place,1,3),6,2);
trace_place=trace_place(:,place_bin+1:2*place_bin);
trace_place_norm=(trace_place-min(trace_place,[],2));
trace_place_norm=trace_place_norm./(max(trace_place,[],2)-min(trace_place,[],2));
[~, arg]=max(trace_place_norm,[],2);
[~, order]=sort(arg,'ascend');
figure; tiledlayout(1,2); nexttile([1 1]);
imagesc(trace_place_norm(order,:))
set(gca,'ytick',[1:size(traces_F,1)],'yticklabel',num2str([order]))
colormap(turbo)

trace_place_PC=movmean(repmat(trace_place(put_pc,:),1,3),6,2);
trace_place_PC=trace_place_PC(:,place_bin+1:2*place_bin);
trace_place_PC_norm=(trace_place_PC-min(trace_place_PC,[],2));
trace_place_PC_norm=trace_place_PC_norm./(max(trace_place_PC,[],2)-min(trace_place_PC,[],2));
[~, arg]=max(trace_place_PC_norm,[],2);
%arg=sum(trace_place_PC_norm.*[1:place_bin],2,'omitnan')./sum(trace_place_PC_norm,2,'omitnan');
[~, order_PC]=sort(arg,'ascend');
nexttile([1 1]);
imagesc(trace_place_PC_norm(order_PC,:))
set(gca,'ytick',[1:size(trace_place_PC_norm,1)],'yticklabel',num2str([put_pc(order_PC)]))
colormap(turbo)


end
%%
show_traces_align_Position_Virmen(traces_F,traces_sp,110*0.3,[-500:3000],Virmen_align,42)
show_traces_place_Virmen(traces_F,traces_sp,10,Virmen_align,[1:10],42)

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

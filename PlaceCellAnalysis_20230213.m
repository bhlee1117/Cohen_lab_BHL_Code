%% Load data
clear
cd /Volumes/ByungHun4TB/Cohen_lab_Data/20221228_BHLm008_Treadmill
load('Treadmill_BHLm008_20230209_result.mat')

%%
place_bin=50;
grps={[1 2],[3 4]};
im_corr_th=[0.994 0.994];
cmap = jet(256); vel_thresh=40;
cmap(1:find(cmap(:,1)==0.5),1) = 0.5;
clear ind_list Arduinos

for i=1%:length(grps)
    %Result{i}.rm_ind=manual_trace_deletion(Result{i});
    %Result{i}.rm_ind=[];
    %Result{i}.rm_ind=[Result{i}.rm_ind; manual_trace_deletion(Result{i})];

    goi=grps{i}(1);
    [sz1 sz2]=size(Result{goi}.ref_im);
    match_list=zeros(size(Result{goi}.coord,1),2,length(grps{i})-1);
    Arduinos{1}=Result{goi}.Arduino;

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
        Arduinos{g}=Result{goi}.Arduino;
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

    traces_sp=[]; traces_F=[]; trace_imcorr=[];
    for g=1:length(grps{i})
        goi=grps{i}(g);
        sp=double(Result{goi}.spike(ind_list(:,g),:)); sp(:,end:size(Result{goi}.Arduino,1))=0;
        F=Result{goi}.traces_res_hi(ind_list(:,g),:);  F(:,end:size(Result{goi}.Arduino,1))=0;
        im_corr=movmean(Result{goi}.im_corr,20);
        trace_imcorr=[trace_imcorr im_corr];
        F(:,im_corr<im_corr_th(g))=NaN; sp(:,im_corr<im_corr_th(g))=NaN;
        traces_sp=[traces_sp sp]; traces_F=[traces_F F];

    end
    Arduino_align=realign_arduino(Arduinos);

    Arduino_align=cell2mat(Arduino_align');
    show_traces_spikes(traces_F,traces_sp,[Arduino_align(:,2) Arduino_align(:,3)*3000])

    % for r=1:size(Result{i}.rm_ind,1)
    % traces_F(:,[Result{i}.rm_ind(r,1):Result{i}.rm_ind(r,2)])=NaN;
    % traces_sp(:,[Result{i}.rm_ind(r,1):Result{i}.rm_ind(r,2)])=NaN;
    % end
Lap_imcorr=get_LapFR(trace_imcorr,place_bin,Arduino_align,vel_thresh,1);
Lap_FR=get_LapFR(traces_sp,place_bin,Arduino_align,vel_thresh,[1:size(traces_sp)]);
[trace_place SI put_pc]=traces2place(traces_sp,place_bin,Arduino_align,75);
%trace_place=trace_place(find(max(trace_place,[],2)>13),:);
trace_place=movmean(repmat(trace_place,1,3),5,2,'omitnan');
trace_place=trace_place(:,place_bin+1:2*place_bin);
trace_place_norm=(trace_place-min(trace_place,[],2));
trace_place_norm=trace_place_norm./(max(trace_place,[],2)-min(trace_place,[],2));

%[~, arg]=max(trace_place_norm(put_pc,:),[],2);
[~, arg]=max(trace_place_norm,[],2);
[~, order]=sort(arg,'ascend');
figure;
%imagesc(trace_place_norm(put_pc(order),:))
imagesc(trace_place_norm(order,:))
%set(gca,'ytick',[1:size(traces_F,1)],'yticklabel',num2str([put_pc(order)]))
set(gca,'ytick',[1:size(traces_F,1)],'yticklabel',num2str([order]))
colormap(turbo)


end
%%
for i=[21 41 42] 
    show_traces_place(traces_F,traces_sp,place_bin,Arduino_align,[1:20],order(i))
end
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

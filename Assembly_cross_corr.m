function [Assembles]=Assembly_cross_corr(Result,wind_size,corr_th,N_patt)
%wind_size=100;
interval=wind_size/4;
bin=[interval+1:interval:size(Result.spike,2)-wind_size]; wind=[0:wind_size];
%gaussFilter = gausswin(15,3);
template=smoothdata(movmean(Result.spike,30,2),2,'gaussian',10);
target=movmean(Result.spike,5,2);

coord=Result.coord;
for b=1:length(bin)

    patt=target;
    patt(:,bin(b)+[0:wind_size*2]-interval)=0;
    patt=patt(:,1:size(target,2));
Temp_neuron(b)=sum(sum(template(:,bin(b)+wind),2)>0);
CC = normxcorr2(template(:,bin(b)+wind),patt); % normxcorr2 requires a smaller template than the image!
xcorr_img(b,:)=CC(size(target,1),:);
end
xcorr_img=xcorr_img(:,wind_size/2+1:end-wind_size/2);
%%
xcorr_img_th=zeros(size(xcorr_img,1),size(xcorr_img,2));
for i=1:size(xcorr_img,1)
    [pks s_tmp width prom]=findpeaks(xcorr_img(i,:),'MinPeakDistance',interval);
    s_tmp=s_tmp(find(pks>corr_th));
    xcorr_img_th(i,s_tmp)=xcorr_img(i,s_tmp);
end
[corr_list, corr_pair]=sort(xcorr_img_th(:),'descend');

% [tmp_count, tmp_list]=sort(sum(xcorr_img_th>0,2),'descend');

[~, tmp_list]=sort(Temp_neuron,'descend'); 
N_matchPatt=sum(xcorr_img_th>0,2); tmp_count=N_matchPatt(tmp_list);

for i=1:N_patt
   % figure;

    cmap=jet(tmp_count(i));
    
    tmp=template(:,bin(tmp_list(i))+wind);
    tmp_ext=template(:,bin(tmp_list(i))+[0:wind_size*1.5]-wind_size/4);

    %centroid=tmp_ext*([0:wind_size*1.5]-wind_size/4)'./sum(tmp_ext,2); MaxF=sum(tmp,2);
    
    [MaxF, centroid]=max(tmp_ext,[],2); centroid(find(sum(tmp,2)==0))=NaN;
    [~, order]=sort(centroid,'ascend'); ind_ass=find(~isnan(centroid));
%     clear xcc
%     for j=2:length(ind_ass)
%         xcc(j-1,1)=max(xcorr(tmp_ext(order(j-1),:),tmp_ext(order(j),:)));
%     end

    
    c_ftprnt=Result.c_ftprnt;
    c_ftprnt=c_ftprnt(:,:,order);
    frm_ind=find(xcorr_img_th(tmp_list(i),:)>0);
    coord_m=mean(coord(ind_ass,:),1); coord_list=coord(order(1:length(ind_ass)),:);
    %W=[1:length(ind_ass)]-mean([1:length(ind_ass)]);
    C_time=centroid(order(1:length(ind_ass)))-min(centroid(order(1:length(ind_ass))));
    cmap_frnt=zeros(size(centroid,1),3)+0.2;
    jet_time=jet(round(max(C_time))+1);
    cmap_frnt(1:sum(~isnan(centroid)),:)=jet_time(round(C_time)+1,:);

    %W=1./(C_time(2:end)-C_time(1:end-1)); W(W==Inf)=1;
    %W=xcc./range(xcc);
    %W=MaxF(order(2:length(ind_ass))).*MaxF(order(1:length(ind_ass)-1)); W=W./range(W);
    %W=C_time(2:end)-C_time(1:end-1);
    %Seq_vec=-sum((coord_m-coord(order(1:length(ind_ass)),:)).*W')./std(coord(order(1:length(ind_ass)),:),0,1);
    %R=coord_list(2:end,:)-coord_list(1:end-1,:); %Rij
    clear R T
    for c1=1:size(coord_list,1)
        R(c1,:,1)=coord_list(:,1)-coord_list(c1,1); %X
        R(c1,:,2)=coord_list(:,2)-coord_list(c1,2); %Y
        T(c1,:)=C_time-C_time(c1);
    end
    %Seq_vec=sum(Seq_vec.*W./sqrt(sum(Seq_vec.^2,2))); Seq_vec=Seq_vec./sqrt(sum(Seq_vec.^2,2));
    %Seq_vec=sum(Seq_vec); Seq_vec=Seq_vec./sqrt(sum(Seq_vec.^2,2));
    %RT=R.*W; %R*T;
    RT=squeeze(sum(R.*T,2));
    Seq_vec=sum(RT); %Seq_vec=(Seq_vec./sqrt(sum(Seq_vec.^2,2)));
    [theta,rho]=cart2pol(RT(:,1),RT(:,2)); [Seq_vec_theta,Seq_vec_rho]=cart2pol(Seq_vec(:,1),Seq_vec(:,2));
    

    nexttile([1 1])
    ftax=imshow2(squeeze(sum(c_ftprnt.*reshape(cmap_frnt,1,1,[],3),3)),[]); hold all;
    cb=colorbar; cb.Ticks=[0 0.5 1]; cb.TickLabels=num2str([0 max(C_time)/Result.frm_rate/2 max(C_time)/Result.frm_rate]','%2.2f');
    cb.Label.String ='Time delay (s)';
    colormap(cb,'jet')
    q=quiver(30,30,Seq_vec(:,1),Seq_vec(:,2),0.3,'w');
    q.ShowArrowHead = 'on'; q.LineWidth = 1.5; q.MarkerSize=12; q.MaxHeadSize=2;
    
    nexttile([1 1])
    polarplot(-theta,rho,'.','markersize',12); hold all;
    polarplot(-[0 Seq_vec_theta],[0 Seq_vec_rho],'r','markersize',12,'marker','*')
    
    nexttile([1 1])
    imagesc(tmp_ext(order,:))
    hold all
    line([wind_size/4 wind_size/4],[0 size(target,1)],'color','y','linestyle','--','linewidth',2)
    line([size(tmp,2)+wind_size/4 size(tmp,2)+wind_size/4],[0 size(target,1)],'color','y','linestyle','--','linewidth',2)
    set(gca,'Xtick',[wind_size/4 size(tmp,2)+wind_size/4],'Xticklabel',[bin(tmp_list(i)) bin(tmp_list(i))+wind_size-1])

    nexttile([1 3])
    imagesc(target(order,:))
    hold all
    line([bin(tmp_list(i)) bin(tmp_list(i))],[0 size(target,1)],'color','y','linestyle','--','linewidth',2)
    line([bin(tmp_list(i))+wind_size bin(tmp_list(i))+wind_size],[0 size(target,1)],'color','y','linestyle','--','linewidth',2)
    colormap('gray')

    for j=1:tmp_count(i)
        line([frm_ind(j)-wind_size/2 frm_ind(j)-wind_size/2],[0 size(target,1)],'color',cmap(j,:),'linewidth',2)
        line([frm_ind(j)+wind_size/2 frm_ind(j)+wind_size/2],[0 size(target,1)],'color',cmap(j,:),'linewidth',2)
    end

    for j=1:tmp_count(i)

        ax1=nexttile([1 1]);
        try
            imagesc(target(order,[frm_ind(j)-wind_size/2:frm_ind(j)+wind_size/2]))
            set(gca,'Xtick',[1 wind_size],'Xticklabel',[frm_ind(j)-wind_size/2 frm_ind(j)+wind_size/2])
        catch
            if frm_ind(j)-wind_size/2<1
                imagesc(target(order,[1:frm_ind(j)+wind_size/2]))
                set(gca,'Xtick',[1 frm_ind(j)+wind_size/2-1],'Xticklabel',[1 frm_ind(j)+wind_size/2])
            else
                imagesc(target(order,[frm_ind(j)-wind_size/2:size(target,2)]))
                set(gca,'Xtick',[1 size(target,2)-frm_ind(j)+wind_size/2+1],'Xticklabel',[frm_ind(j)-wind_size/2 size(target,2)])
            end
        end
        ax1.XColor=cmap(j,:); ax1.YColor=cmap(j,:); ax1.LineWidth=1.5;
        title(['Corr Coeff = ' num2str(xcorr_img_th(tmp_list(i),frm_ind(j)))])

    end

    Assembles{i}.order=order;
    Assembles{i}.template=tmp_ext;
    Assembles{i}.maxtime=C_time;
    Assembles{i}.Corr_frm=[xcorr_img_th(tmp_list(i),frm_ind)' frm_ind'];
end
end



%imagesc(xcorr_img_th)
function [cluster]=Find_cluster_online(traces,spike)
 tr=traces-median(traces,2); 
 tr=tr./get_threshold(traces,1);

% ax1 = subplot(8,2,1:12);
% cmap=distinguishable_colors(3);
% % tr=Result{f}.traces-median(Result{f}.traces,2); fprnt=Result{f}.c_ftprnt;
% % tr=tr./get_threshold(Result{f}.traces,1);
% t=[1:size(tr,2)]/800;
% %noi=[1:size(tr,1)];
%  s=4;
%  noi=(cluster_order{s})';
% %noi=list(o)';
% 
% line_color=zeros(size(noi,2),3); line_color2=[];
% order=[];
% 
% % for j=2:size(Result{f}.DMDPatt,2)
% % line_color(Result{f}.clist_p{j},:)=repmat(cmap(j,:),length(Result{f}.clist_p{j}),1);
% % end
% 
% lines=plot(t(1:end-2),tr(noi,1:end-2)'+[1:size(noi,2)]*scale);
% arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(line_color,2))
% hold all
% S=tr; S(~spike)=NaN;
% plot(t(1:end-2),S(noi,1:end-2)'+[1:size(noi,2)]*scale,'r.')
% set(gca,'ytick',[1:size(noi,2)]*scale,'yticklabel',noi)
% axis tight
% 
% ax2 = subplot(8,2,14:16);
% plot(t(1:end-2),sum(tr(:,1:end-2)))
% hold all
% yyaxis right
% plot(t,movsum(sum(Result{f}.spike),30,2))
% axis tight
% 
% linkaxes([ax1 ax2],'x')
%%
clear synch_trigg synch_trigg_spike
cluster.thres=size(traces,1)*0.6;
synch=bwlabel(movsum(sum(spike),30,2)>cluster.thres);
for i=1:max(synch)
    s=find(synch==i);
ind=[s(1)-200:s(1)+100];
cluster.synch_time(i)=s(1);
try
synch_trigg(:,:,i)=reshape(tr(:,ind),size(tr,1),length(ind));
synch_trigg_spike(:,:,i)=reshape(spike(:,ind),size(tr,1),length(ind));
end
end

%%
clear corrMat
for s=1:size(synch_trigg,3)
for i=1:size(synch_trigg,1)
    for j=1:size(synch_trigg,1)
        x1=zscore(synch_trigg(i,:,s));
        x2=zscore(synch_trigg(j,:,s));
[cluster.corrMat(i,j,:,s),lags] = xcorr(x1,x2,100);
    end
end

[r c]=find(max(cluster.corrMat(:,:,[80:120],s),[],3)>200);
x=find(r>=c);
r(x)=[]; c(x)=[]; 
cluster.cluster_N{s}=unique([r;c]);
cluster.cluster_N{s}(sum(synch_trigg_spike(cluster.cluster_N{s},[150:230],s),2)==0)=[];

[m arg]=max(cluster.corrMat(cluster.cluster_N{s}(1),cluster.cluster_N{s},:,s),[],3);
[~, o]=sort(arg,'ascend');
cluster.cluster_order{s}=cluster.cluster_N{s}(o);

clear c2
for n=1:length(cluster.cluster_N{s})
[r2]=find(synch_trigg_spike(cluster.cluster_N{s}(n),[150:230],s));
c2(n)=r2(1); end
[~, o]=sort(c2,'ascend');
cluster.cluster_order{s}=cluster.cluster_N{s}(o);

end
%%
cluster.list=unique(cell2mat(cluster.cluster_N'));
cluster.counts=histcounts(cell2mat(cluster.cluster_N'),cluster.list);
[~, o]=sort(cluster.counts,'ascend');
cluster.list_sorted=cluster.list(o);
end
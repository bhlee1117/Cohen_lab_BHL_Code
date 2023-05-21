
load('E:\004158_BHLm009_ConstStim_Bp05_FOV3\result.mat')
 tr=traces-median(traces,2); 
 tr=tr./get_threshold(traces,1);
%%
figure;
f=7; scale=10;

ax1 = subplot(8,2,1:12);
cmap=distinguishable_colors(3);
 tr=Result{f}.traces-median(Result{f}.traces,2); fprnt=Result{f}.c_ftprnt;
 tr=tr./get_threshold(Result{f}.traces,1);
t=[1:size(tr,2)]/800;
noi=[1:size(tr,1)];
s=4;
noi=(cluster_order{s})';

tiledlayout(10,4)
% 

line_color=zeros(size(noi,2),3); line_color2=[];
order=[];
for j=2:size(Result{f}.DMDPatt,2)
order=[order; noi(find(ismember(noi,Result{f}.clist_p{j})))'];
line_color(find(ismember(noi,Result{f}.clist_p{j})),:)=repmat(cmap(j,:),sum(ismember(Result{f}.clist_p{j},noi)),1);
line_color2=[line_color2; repmat(cmap(j,:),sum(ismember(Result{f}.clist_p{j},noi)),1)];
end
order=[order; setdiff([noi],order)'];
line_color2(end+1:size(noi,2),:)=0;


ax1 = nexttile([7 4]);
lines=plot(t,tr(order,:)'+[1:size(noi,2)]*scale);
%arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(line_color,2))
try arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(line_color2,2))
catch arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(jet(length(noi)),2))
end
axis tight
hold all
S=tr; S(~Result{f}.spike)=NaN;
plot(t,S(order,:)'+[1:size(noi,2)]*scale,'r.')
set(gca,'ytick',[1:size(noi,2)]*scale,'yticklabel',order)
title(fpath_mac2window(fpath{f}),'Interpreter', 'none')
hold all
[a b]=unique(bwlabel(Result{f}.Blue(1,:)>0));
line([t(b(2:end))' t(b(2:end))']',repmat([0 size(noi,2)*scale],size(b,1)-1,1)','color','r')


ax2 = nexttile([2 4]);
plot(t(1:end-2),sum(tr(:,1:end-2)))
hold all
yyaxis right
plot(t,movsum(sum(Result{f}.spike),30,2))
axis tight

linkaxes([ax1 ax2],'x')
%%
clear synch_trigg synch_trigg_spike
synch=bwlabel(movsum(sum(Result{f}.spike),30,2)>80);
for i=1:max(synch)
    s=find(synch==i);
ind=[s(1)-200:s(1)+100];
try
synch_trigg(:,:,i)=reshape(tr(:,ind),size(tr,1),length(ind));
synch_trigg_spike(:,:,i)=reshape(Result{f}.spike(:,ind),size(tr,1),length(ind));
end
end

%%
clear corrMat cluster_N
for s=1:size(synch_trigg,3)
for i=1:size(synch_trigg,1)
    for j=1:size(synch_trigg,1)
        x1=zscore(synch_trigg(i,:,s));
        x2=zscore(synch_trigg(j,:,s));
[corrMat(i,j,:,s),lags] = xcorr(x1,x2,100,'normalized');
    end
end

[r c]=find(max(corrMat(:,:,[80:120],s),[],3)>0.5);
x=find(r>=c);
r(x)=[]; c(x)=[]; 
cluster_N{s}=unique([r;c]);
cluster_N{s}(sum(synch_trigg_spike(cluster_N{s},[170:230],s),2)==0)=[];

[m arg]=max(corrMat(cluster_N{s}(1),cluster_N{s},:,s),[],3);
[~, o]=sort(arg,'ascend');
cluster_order{s}=cluster_N{s}(o);

clear c2
for n=1:length(cluster_N{s})
[r2]=find(synch_trigg_spike(cluster_N{s}(n),[170:230],s));
c2(n)=r2(1); end
[~, o]=sort(c2,'ascend');
cluster_order{s}=cluster_N{s}(o);

end
%%
list=unique(cell2mat(cluster_N'));
h=histcounts(cell2mat(cluster_N'),list);
[~, o]=sort(h,'ascend');
noi=list(o);
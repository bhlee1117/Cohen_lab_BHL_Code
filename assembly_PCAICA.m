function [CAs ic_sub]=assembly_PCAICA(Result,bin_time,theshold,n_comp,figure_on)
% Analysis to find cell assemblies
% 2023/01/15, Byung Hun Lee

%%
i=1;
frm_rate=1/Result.frm_rate;
%bin_time=20*1e-3; 
bin=bin_time/frm_rate;
S_trace=movsum(Result.spike,bin,2);
%S_trace=movsum(Result.traces_hi,bin,2);
%S_trace=zscore(S_trace(:,[1:bin:end]),0,2);
S_trace=zscore(S_trace,0,2);
figure; 
subplot(2,1,1)
imagesc(S_trace)
colormap('jet')
%covMat=corr(S_trace');
covMat=S_trace*S_trace';
[V, D, W] = eig(covMat);

    D = diag(D); 
    D = D(end:-1:1);
    V = V(:,end:-1:1);
    vSign = sign(max(V) - max(-V));  % make the largest value always positive
    V = V.*vSign;

    PCA_trace=zscore(V(:,1:10)'*S_trace,0,2);
subplot(2,1,2)
plot(D/sum(D)); hold all;
D_th=min(find(cumsum(D/sum(D))>0.95));
plot(D_th,D(D_th)/sum(D),'ro')
plot(n_comp,D(n_comp)/sum(D),'ko')


    idx_pc = [1:n_comp];
[icasig, A,W] = fastica((V(:,idx_pc)'*S_trace));
n_ic = size(W,1);
V_ic = (W * V(:,idx_pc)')';
ic_sub = (V_ic(:,1:n_ic)'*S_trace)';

%%
if figure_on
    figure;
    i=1; scale=10;
CAs=zeros(size(V_ic,1),n_ic); order=[];
cmap=jet(n_ic);
%cmap=cmap(randperm(n_ic,n_ic),:);

%cmap=distinguishable_colors(n_ic); %n_add=length(find(sum(cmap,2)<0.5));
%cmap=distinguishable_colors(n_ic+n_add);
%cmap(find(sum(cmap,2)<0.5),:)=1-cmap(find(sum(cmap,2)<0.5),:);

line_color2=[];
t=[1:size(Result.traces_hi,2)];%/Result.frm_rate;

tr=Result.traces-movmedian(Result.traces,300,2); fprnt=Result.c_ftprnt;
tr=tr./get_threshold(Result.traces,1);
%V_ic=zscore(abs(V_ic),0,1);
%V_ic=V_ic-median(V_ic,1);
V_ic_norm=V_ic./std(V_ic,0,1);
for j=1:n_ic
    cell=find(V_ic_norm(:,j)>theshold);
    CAs(cell,j)=1;

end
for c=find(sum(CAs,2)>0)'
    [~ , n]=max(V_ic_norm(c,:));
    CAs(c,:)=0; CAs(c,n)=1;
end
for j=1:n_ic
    cell=find(CAs(:,j));
    order=[order; cell];
    line_color2=[line_color2; repmat(cmap(j,:),length(cell),1)];
end
noi=find(sum(CAs,2)>0)';
line_color2(end+1:size(noi,2),:)=0;


ax2 = subplot(8,2,1:12);
S=tr; S(~Result.spike)=NaN;
lines=plot(t,tr(order,:)'+[1:size(order,1)]*scale);
%lines=plot(t,S(order,:)'+[1:size(order,1)]*scale,'.');
arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(line_color2,2))

axis tight
hold all
plot(t,S(order,:)'+[1:size(order,1)]*scale,'r.')
set(gca,'ytick',[1:size(order,1)]*scale,'yticklabel',order)
hold all
[a b]=unique(bwlabel(Result.Blue(1,:)>0));
line([t(b(2:end))' t(b(2:end))']',repmat([0 size(noi,2)*scale],size(b,1)-1,1)','color','r')

ax3 =  subplot(8,2,13:15);
ic_trace=abs(ic_sub-median(ic_sub,1))./range(ic_sub,1);
lines_ic=plot([1:size(ic_sub,1)],ic_trace+[1:n_ic]);
arrayfun(@(l,c) set(l,'Color',c{:}),lines_ic,num2cell(cmap,2))
axis tight

linkaxes([ax2 ax3],'x')

c_ftprnt=Result.c_ftprnt;
figure;
colr=zeros(size(V_ic,1),3)+1;
for i=1:n_ic
colr(find(CAs(:,i)),:)=repmat(cmap(i,:),sum(CAs(:,i)),1);
end
imshow2(squeeze(sum(c_ftprnt.*reshape(colr,1,1,[],3),3)),[]);
for i=1:n_ic
    Cell=find(CAs(:,i));
text(Result.coord(Cell,1)',Result.coord(Cell,2)',num2str(repmat(i,length(Cell),1)),'color','w')
end
end
end
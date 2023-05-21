% Analysis to find cell assemblies
% 2023/01/15, Byung Hun Lee
clear
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20221227_BHLm008_Cluster'
load('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20221227_BHLm008_Cluster/162338_BHLm008_ContStim_Bp15_Or089_FOV1result.mat')

%%
fp='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20221227_BHLm008_Cluster/162338_BHLm008_ContStim_Bp15_Or089_FOV1';
Sz = importdata([fp '/experimental_parameters.txt']);
sz1=Sz.data(1); sz2=Sz.data(2);
tic;
mov_mc=double(readBinMov([fp '/mc' '.bin'],sz2,sz1));
toc;
%%
mov_hi=mov_mc-movmean(mov_mc,70,3);
mov_res=mov_mc-movmean(mov_mc,200,3);
%%

n_comp=15;
[CAs ic_sub]=assembly_PCAICA(Result{1},30*1e-3,2,n_comp,1);

order=[];
for i=1:n_comp; order=[order; find(CAs(:,i))];
end
S_trace=movsum(Result{1}.spike(order,:),30/1.25,2);
imshow2(corrcoef(S_trace'),[])
hold all
for i=1:n_comp
    line([sum(CAs(:,1:i),[1 2])+0.5 sum(CAs(:,1:i),[1 2])+0.5],[0.5 112.5],'color','r')
    line([0.5 112.5],[sum(CAs(:,1:i),[1 2])+0.5 sum(CAs(:,1:i),[1 2])+0.5],'color','r')
end
%%
%tr=Result{1}.traces_hi./get_threshold(Result{1}.traces_hi,1);
tr=movsum(Result{1}.spike,24,2);
CS_vec=cal_cross_corr(tr,Result{1}.coord,0.7);
%imshow2(mean(mov_mc,3),[])
hold all
q=quiver(CS_vec(:,1),CS_vec(:,2),CS_vec(:,3),CS_vec(:,4),0,'r');
q.ShowArrowHead = 'on';
q.Marker = '.';
%%
PC=15; wind=[-150:300];
movSTA=zeros(sz2,sz1,length(wind));
g=0; noi=find(sum(CAs(:,PC),2)>0);
sp_trace=movsum(Result{1}.spike(noi,:),10/1.25,2)>0;
sp_t=find(sum(sp_trace,1)>length(noi)*0.8);
%sp_t=locmax;
n_sp=length(sp_t);

for i=1:n_sp
    i
    try
        movSTA = movSTA + mov_res(:,:,sp_t(i)+wind);
        g=g+1;
    end
end
movSTA=movSTA./g;
figure;
moviefixsc(-movSTA,[0 100])
hold all
plot(Result{1}.coord(noi,1),Result{1}.coord(noi,2),'ro','markersize',12)

%%

i=1; f_ref=3515;
tr_tmp=Result{i}.traces_hi;
tr_sp=Result{i}.spike;
clear CS_spike
for j=1:size(tr_tmp,1)
sp_l=find(tr_sp(j,:))-f_ref; sp_l(sp_l<0)=[]; 
delt(j)=min(sp_l);
[CS_list CS_spike(j,:)]=find_CS(tr_tmp(j,:),tr_sp(j,:),12,4);
end
[~, order]=sort(delt,'ascend');
show_traces_spikes(tr_tmp(order,:),tr_sp(order,:),Result{i}.Blue)

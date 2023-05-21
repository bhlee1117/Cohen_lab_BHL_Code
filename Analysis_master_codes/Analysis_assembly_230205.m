% Analysis to find cell assemblies
% 2023/02/05, Byung Hun Lee
clear
addpath(genpath('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Code'))

cd '/Volumes/ByungHun4TB/Cohen_lab_Data/20230206_BHLm017_Treadmill_Assembly'
load('223426_BHLm017_Contst_Bp05_Awakeresult.mat')

%%
fp='/Volumes/ByungHun4TB/Cohen_lab_Data/20230206_BHLm017_Treadmill_Assembly/223426_BHLm017_Contst_Bp05_Awake';
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
[Patt]=Assembly_cross_corr(Result{1},100,0.25,2);
%STA=make_STA(mov_res,Patt{14}.Corr_frm(:,2)'-50,[-100:150]);
%%
for i=1:length(Result_perm)
xcorr_img_th(:,:,i)=apply_CCtemplate(Result_perm{i},Patt,0.25);
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
clear;
tmp=importdata('/Volumes/BHL_WD18TB/20230901_PP68_Adaptation_IsoKetCpp/Result_V2CheRiffST_20240110.mat');
Ad_Result=tmp.Result;
%%
figure(1); 
for n=[89]%:length(Ad_Result)
pass=0; 
tr=Ad_Result{n}.trace; tr=tr-prctile(tr,20);
Ad_Result{n}.normtrace=tr./get_threshold(tr,1);
tr=tr-movmedian(tr,200); tr=tr./get_threshold(tr,1);
Ad_Result{n}.spike=find_spike_bh(tr,5,3);
Ad_Result{n}.Subth=get_subthreshold(Ad_Result{n}.normtrace,Ad_Result{n}.spike,3,5);


CS_par=[8 1];
while pass~=1    
clf;
Sub_hi=Ad_Result{n}.Subth-movprc(Ad_Result{n}.Subth,200,30);
plot(Ad_Result{n}.normtrace); hold all
plot(find(Ad_Result{n}.spike),Ad_Result{n}.normtrace(find(Ad_Result{n}.spike)),'ro')
plot(Sub_hi,'c')
xlim([3000 12000])    

[trans tr_trace]=detect_transient(Sub_hi,CS_par,Ad_Result{n}.spike);
CS_ind=find(trans.spike_number>2 & trans.mean_ISI<15);
Ad_Result{n}.CS_trace=ismember(tr_trace,CS_ind);
tmp_tr=Ad_Result{n}.normtrace; tmp_tr(Ad_Result{n}.CS_trace==0)=NaN;
plot(tmp_tr,'r')

pass=input('Pass? ');
if pass==0
par1=input('New CS par top');
par2=input('New CS par down');
CS_par=[par1 par2];
end

end

if isempty(pass)
    break;
end

end

%%
Result=Ad_Result; %fpath=tmp.fpath;
save('/Volumes/BHL_WD18TB/20230901_PP68_Adaptation_IsoKetCpp/Result_V2CheRiffST_20240110.mat','Result','-v7.3')

run_line('/Users/bhlee1117/Documents/GitHub/Cohen_lab_BHL_Code/Analysis_master_codes/Analysis_20250331_FigureMakingSeeSaw.m',[1 2])

%% Simulate PSP
dt=0.001; T=15; tax=[0:dt:T-dt]; BlueOnTime=[5 10];
E_spk=[rand(1,length(tax))]< 10*dt; I_spk=[rand(1,length(tax))]< 10*dt;
Seesaw_spk=[rand(1,length(tax))]< 10*dt;
O_stim=double(tax>BlueOnTime(1) & tax<BlueOnTime(2));

V_opto=simulate_PSP([zeros(1,length(tax))], [zeros(1,length(tax))], O_stim, -65, dt, T, 0, 0.2);
V_exp=simulate_PSP(E_spk, [zeros(1,length(tax))], O_stim, -65, dt, T, 0, 0.2);
V_inh=simulate_PSP([zeros(1,length(tax))], I_spk , O_stim, -65, dt, T, 0, 0.2);
figure(302); clf;
V_total=simulate_PSP(E_spk, I_spk, O_stim, -65, dt, T, 0, 0.2);

V_BasalEI=simulate_PSP(Seesaw_spk, [zeros(1,length(tax))], O_stim, -65, dt, T, 0, 0.2);
V_ApicalEI=simulate_PSP([zeros(1,length(tax))], Seesaw_spk , O_stim, -65, dt, T, 0, 0.2);

V_BasalE=simulate_PSP(Seesaw_spk, [zeros(1,length(tax))], O_stim, -65, dt, T, 0, 0.2);
V_ApicalE=simulate_PSP(Seesaw_spk, [zeros(1,length(tax))] , O_stim, -65, dt, T, 0, 0.2);

V_total=simulate_PSP(E_spk, I_spk, O_stim, -65, dt, T, 0, 0.2);

V_rest_opto=median(V_opto(BlueOnTime(1)/dt:BlueOnTime(2)/dt));
[~, ESPAmp texp]=get_STA(V_exp,E_spk,5,5);
[~, ISPAmp tinh]=get_STA(V_inh,I_spk,5,5);
ESPAmpOff=squeeze(max(ESPAmp(:,texp<BlueOnTime(1)/dt | texp>BlueOnTime(2)/dt,:),[],3))+65;
ESPAmpOn=squeeze(max(ESPAmp(:,texp>BlueOnTime(1)/dt & texp<BlueOnTime(2)/dt,:),[],3))-V_rest_opto;

ISPAmpOff=(squeeze(min(ISPAmp(:,tinh<BlueOnTime(1)/dt | tinh>BlueOnTime(2)/dt,:),[],3))+65);
ISPAmpOn=(squeeze(min(ISPAmp(:,tinh>BlueOnTime(1)/dt & tinh<BlueOnTime(2)/dt,:),[],3))-V_rest_opto);

figure(303); clf;
Boxplot_wPoints2({ESPAmpOff,ISPAmpOff,ESPAmpOn,ISPAmpOn},[0.7 0 0; 0 0 0.7;1 0 0; 0 0 1]);
set(gca,'XTick',[1:4],'XTickLabel',{'EPSP, Blue off','IPSP, Blue off','EPSP, Blue on','IPSP, Blue on'});
ylabel('V_{PSP} - V_{rest} (mV)')

Case_sp=[]; V_case=[];
Case_sp{1}=zeros(4,4000); Case_sp{1}([1 4],2001:2200)=repmat([rand(1,200)]< 30*dt,2,1); %Basal E/I, Apical E/I
Case_sp{2}=zeros(4,4000); Case_sp{2}(1:4,:)=repmat([rand(1,4000)]< 15*dt,4,1); 
Case_sp{2}(1,2001:2200)=Case_sp{2}(1,2001:2200)+[rand(1,200)]< 90*dt; Case_sp{2}(3,2001:2200)=0;
Case_sp{3}=zeros(4,4000); Case_sp{3}(1:4,:)=repmat([rand(1,4000)]< 10*dt,4,1); 
Case_sp{3}(4,2001:2200)=Case_sp{3}(4,2001:2200)+[rand(1,200)]< 90*dt; Case_sp{3}(2,2001:2200)=0;

for c=1:3
V_case{c,1,1}=simulate_PSP(Case_sp{c}(1,:), Case_sp{c}(2,:) , zeros(1,4000), -65, dt, 4, 0, 0.2);
V_case{c,2,1}=simulate_PSP(Case_sp{c}(3,:), Case_sp{c}(4,:) , zeros(1,4000), -65, dt, 4, 0, 0.2);
V_case{c,1,2}=simulate_PSP(Case_sp{c}(1,:), Case_sp{c}(2,:) , ones(1,4000), -65, dt, 4, 0, 0.2);
V_case{c,2,2}=simulate_PSP(Case_sp{c}(3,:), Case_sp{c}(4,:) , ones(1,4000), -65, dt, 4, 0, 0.2);
end
V_case_mean=cellfun(@(x) mean(x(1,2100:2300)),V_case);

figure(304); clf; tiledlayout(1,3); ax1=[];
for c=1:3
ax1=[ax1 nexttile([1 1])];
plot([1 2],V_case_mean(c,:,1),'color',[1 0 1]); hold all
plot([1 2],V_case_mean(c,:,2),'color',[0 1 1]);
set(gca,'XTick',[1:2],'XTickLabel',{'Basal','Apical'});
ylabel('Vm (mV)');
box off
end
linkaxes(ax1);
ylim([-70 -40])
legend({'Blue off','Blue on'})
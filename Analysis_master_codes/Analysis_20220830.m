%%
clear
load('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20220627_Treadmill/013646_20220626_M2_Treadmill_Stim_p4/20220830_voltage_trace.mat')
load('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20220627_Treadmill/013646_20220626_M2_Treadmill_Stim_p4/settings.mat')
fpth='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20220627_Treadmill/';
fnm='20220626_M2_last.csv';

[Arduino_data reward_pos lap_dist]=match_treadmill_DAQ([fpth fnm],1e-5,DAQ_data,1.27*1e-3,2); % time, treadmill, Reward, Run
Blue=DAQ_waves.amplitude(4,round([1:size(volt_hi,2)]*1.27*1e-3/1e-5));
%%
figure;
scale=0.4;
tiledlayout(10,4)
ax1 = nexttile([8 4]);
plot((volt_hi+scale*[1:size(volt_hi,1)]')')
hold all
line(repmat([0 size(volt_hi,2)],[size(volt_hi,1) 1])',repmat([1:size(volt_hi,1)]*scale,[2 1]),'color',[0.7 0.7 0.7])
axis tight

% ax2 = nexttile([4 4]);
% plot((volt+0.4*[1:size(volt,1)]')')
% hold all
% line(repmat([0 size(volt,2)],[size(volt,1) 1])',repmat([1:size(volt,1)]*0.4,[2 1]),'color',[0.7 0.7 0.7])
% axis tight

ax3 = nexttile([1 4]);
plot(Blue)
ylabel('Blue LED')
axis tight

ax4 = nexttile([1 4]);
plot(Arduino_data(:,2))
ylabel('Position')
axis tight

linkaxes([ax1 ax3 ax4],'x')
%linkaxes([ax1 ax2],'xy')

%%

figure;
noi=67;
tiledlayout(4,1)
ax1 = nexttile([1 1]);
plot((volt_hi(noi,:)))
axis tight

ax2=nexttile([2 1]);
[wwt fs]=cwt(volt_hi(noi,:),800);
imagesc(abs(wwt));
set(gca,'ytick',[1:round(size(wt,1)/10):size(wt,1)],'yticklabel',num2str(fs([1:round(size(wt,1)/10):size(wt,1)]),'%1.2f'))
ylabel('Frequency (Hz)')

ax3 = nexttile([1 1]);
% plot(Arduino_data(:,4),'-')
% ylabel('Run')
% ylim([-0.2 1.5])
plot(Blue)
ylabel('Blue LED')

linkaxes([ax1 ax2 ax3],'x')
%%
for i=1:size(volt,1)
    [wt(:,:,i) fs]=cwt(volt_hi(i,:),800);
    burst_pow(i,:)=sum(abs(wt(6:15,:,i)),1);
    theta_pow(i,:)=sum(abs(wt(50:62,:,i)),1);

    thres=get_threshold(volt_hi(i,:));
    spike(i,:)=zeros(1,size(volt_hi,2));
    spike_tr(i,:)=zeros(1,size(volt_hi,2));
    [pks s_tmp width prom]=findpeaks(volt_hi(i,:));
    s_tr=s_tmp;
    s_tmp=s_tmp(find(pks>thres));
    spike_tr(i,s_tr)=1;
    spike(i,s_tmp)=1;
end
burst_pow_norm=burst_pow-median(burst_pow,2); burst_pow_norm=burst_pow_norm./range(burst_pow_norm,2);
theta_pow_norm=theta_pow-median(theta_pow,2); theta_pow_norm=theta_pow_norm./range(theta_pow_norm,2);

[trans trans_traces]=detect_transient(volt_hi_g,burst_pow_norm);
%%
nice_ind=[4 5 7 8 9 26 29 32 35 36 37 38 45 48 52 53 57 58 59 60 61 62 63 64 65 66 67];
volt_hi_g=volt_hi(nice_ind,:);

figure;
tiledlayout(10,4)
ax1 = nexttile([8 4]);

plot((volt_hi_g+scale*[1:size(volt_hi_g,1)]')','color','k')
hold all
%plot((burst_pow_norm(nice_ind,:)*0.5+scale*[1:size(volt_hi_g,1)]')','color','r')
%plot((theta_pow_norm(nice_ind,:)+scale*[1:size(volt_hi_g,1)]')','color',[0 0.5 1])
axis tight

ax4 = nexttile([1 4]);
plot(Blue)
ylabel('Blue LED')
axis tight

linkaxes([ax1 ax4],'x')

%% detect transient
tic;
trans_traces2=trans_traces;
trans_spike=trans_traces.*spike(nice_ind,:);
for i=1:length(trans)
    i
[n ind]=groupcounts(trans_spike(i,:)');
%t=intersect(find(trans(i).int<0.7),ind(find(n<2)))';
t=unique([find(trans(i).int<0.7) setdiff([1:length(trans(i).length)],ind(find(n>1))')]);
for j=t
trans_traces2(i,trans_traces(i,:)==j)=0;
end
end
toc;

figure;tiledlayout(10,4)
ax1 = nexttile([8 4]);
[m, order]=sort(sum(spike(nice_ind,5500:6500),2),'ascend');
%order=[1:size(volt_hi_g,1)];
plot((volt_hi_g(order,:)+scale*[1:size(volt_hi_g,1)]')','color',[0.8 0.8 0.8])
hold all

tmp=(volt_hi_g(order,:)+scale*[1:size(volt_hi_g,1)]')';
tmp(trans_traces2(order,:)'==0)=NaN;
plot(tmp,'color',[1 .4 0])
text(zeros(size(volt_hi_g,1),1),[1:size(volt_hi_g,1)]'*scale,num2str(nice_ind(order)'))

tmp2=(volt_hi_g(order,:)+scale*[1:size(volt_hi_g,1)]')';
tmp2((spike(nice_ind(order),:)==0)')=NaN;
plot(tmp2,'r.')
axis tight
ax4 = nexttile([1 4]);
plot(Arduino_data(:,2))
ylabel('Position')
axis tight

linkaxes([ax1 ax4],'x')

%% Place triggered average trace; Voltage trace, Spike
place_bin=20;
[place_volt order] = show_place_field(volt_hi_g,Arduino_data(:,2),Arduino_data(:,4),place_bin,lap_dist,reward_pos);
[place_spike order2] = show_place_field(spike(nice_ind,:),Arduino_data(:,2),Arduino_data(:,4),place_bin,lap_dist,reward_pos);

figure;
tiledlayout(10,4)
ax1 = nexttile([4 4]);
imagesc(place_volt(order,:)); hold all;
line(repmat(ceil(reward_pos/(lap_dist/place_bin)),1,2),[0 length(order)+0.5],'color','r');
ax2 = nexttile([4 4]);
imagesc(place_spike(order,:)); hold all;
line(repmat(ceil(reward_pos/(lap_dist/place_bin)),1,2),[0 length(order)+0.5],'color','r');
%%
% now plot by every lap
place_bin=25;
volt_cs=(volt_hi_g); volt_cs(trans_traces2'==0)=0;
volt_ss=(spike(nice_ind,:)); volt_ss(trans_traces2'==1)=0;

[lap_place_volt temp]=show_lap_field(spike,Arduino_data(:,2),Arduino_data(:,4),nice_ind,place_bin,lap_dist,reward_pos,1);
[lap_place_volt temp]=show_lap_field(volt_ss,Arduino_data(:,2),Arduino_data(:,4),[1:27],place_bin,lap_dist,reward_pos,1);
[lap_place_volt temp]=show_lap_field(volt_cs,Arduino_data(:,2),Arduino_data(:,4),[1:27],place_bin,lap_dist,reward_pos,1);
%%
target_bin=12; noi=[62 66];
show_lap_rawtrace(volt_hi,spike,Arduino_data(:,2),target_bin,noi,place_bin,lap_dist)
%%
spike_number=movsum(spike(nice_ind(order),3500:7500),400,2)';
lines = plot(spike_number);
colr=jet(size(spike_number,2));
% lines = plot(idx_t*dt,traces_plot+cumsum(range(traces_plot)));
arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(colr,2))

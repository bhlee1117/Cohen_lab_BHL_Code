function [output reward_pos lap_each_dist]=match_treadmill_DAQ(filename,dt,DAQ_data,frm_rate,thres)
% Match DAQ and arduino signal
%

% INPUTS 
% Arduino data of treadmill and DAQ data of reward signal

% OUTPUTS
% Interpolated track, time, reward signal


% MODIFICATION HISTORY : 
%     2022.07.03.
%     Byung Hun Lee, Created
%     2022.08.03.
%     Corrected interpolation, Byung Hun Lee


T = readtable(filename,'NumHeaderLines',3);
%%
l=[3:3:size(T,1)]; ll=find(~isnan(T.Var2)); l_ind=3-ll(1);
t_series=(T.Var1(l)); 
try
t_series=milliseconds(t_series-t_series(1))*10^-3;
catch
t_series=datetime(T.Var1(l),'InputFormat','mm:ss.SSS');
t_series=milliseconds(t_series-t_series(1))*10^-3;
end

Reward=T.Var3(l(2:end)-l_ind+1)-T.Var3(l(1:end-1)-l_ind+1); R_time=t_series(find(Reward==1));

Tread=T.Var2(l-l_ind); 
dif=Tread(2:end)-Tread(1:end-1); 
dif(find(abs(dif)>5000 | isnan(dif)))=0;
Tread=[0;cumsum(dif)];
%converto 100 Hz
[a b]=unique(t_series(1:end-1));
remap_t=[0.01:0.01:max(t_series)];
Tread=interp1(a,Tread(b),remap_t); %dif=Tread(2:end)-Tread(1:end-1); 

Tread_run=movmean(Tread,4);
R_t=zeros(size([0.01:0.01:max(t_series)],2),1);
rr=round(R_time/0.01); rr(rr==0)=[];
R_t(rr)=1; 
%R_t(17743)=1;
%R_t([13062 14254 16953 20067])=1;
ind_Reward=find(R_t);

difrun=Tread_run(2:end)-Tread_run(1:end-1);
run=[0 abs(difrun>thres)];

Rx=Tread(find(R_t)); lap_each_dist=(Rx(2:end)-Rx(1:end-1)); lap_dist=mean(lap_each_dist);

% exp_lap_dist=6774; 
% skip_R=find(abs(lap_each_dist-exp_lap_dist)>exp_lap_dist*0.2);
% if ~isempty(skip_R)
%     for sk=1:length(skip_R)
%         [~, ind_inc]=min(abs(Tread-(Tread(ind_Reward(skip_R(sk)))+exp_lap_dist)));
%         R_t(ind_inc)=1;
%     end
% end
ind_Reward=find(R_t);
Rx=Tread(find(R_t)); lap_each_dist=(Rx(2:end)-Rx(1:end-1));

reward_pos=Tread(ind_Reward(1));
Tread_reg=Tread;
for i=1:sum(R_t)-1 % register rotary encoder error
    Tread_reg(ind_Reward(i):ind_Reward(i+1))=Tread_reg(ind_Reward(i)-1)+...
        (Tread(ind_Reward(i):ind_Reward(i+1))-Tread(ind_Reward(i)))/(Tread(ind_Reward(i+1))-Tread(ind_Reward(i)))*lap_dist;
end
Tread_reg(ind_Reward(sum(R_t)):end)=Tread_reg(ind_Reward(sum(R_t))-1)+Tread_reg(ind_Reward(sum(R_t)):end)-Tread_reg(ind_Reward(sum(R_t)));

%%
t_DAQ=[dt:dt:dt*(size(DAQ_data,1))]; %*(frm_rate)/(frm_rate+0.00002);
s=find(DAQ_data(:,1)==1); 
DAQ_Arduino_t_offset=t_DAQ(s(1))-R_time(1);
[a start_point]=min(abs(remap_t+DAQ_Arduino_t_offset));
[a end_point]=min(abs(remap_t+DAQ_Arduino_t_offset-t_DAQ(end)));

t_output=[frm_rate:frm_rate:t_DAQ(end)];
[a b]=unique(remap_t(start_point:end_point)+DAQ_Arduino_t_offset);
c=Tread_reg(start_point:end_point);
tread_output=interp1(a,c(b),t_output,'linear');
while sum(tread_output>lap_dist)>0
    ll=find(tread_output>lap_dist);
    tread_output(ll(1):end)=tread_output(ll(1):end)-lap_dist;
end

%tread_output=interp1(Tread(start_point:end_point),t_output);

d=run(start_point:end_point);
run_output=interp1(a,double(d(b)),t_output,'linear');

reward_output=zeros(size(t_output,2),1);
reward_output(round((R_time+DAQ_Arduino_t_offset)/frm_rate))=1;
reward_output=reward_output(1:length(t_output));

output=[t_output' tread_output' reward_output ceil(run_output)'];
output(isnan(output))=0;

gitter=0.07;

DAQ_data_ds=sum(reshape(DAQ_data(:,1),round(frm_rate/(dt)),[]))>0;
N_match=sum(DAQ_data_ds(ceil(gitter/frm_rate):end).*output(1:end-ceil(gitter/frm_rate)+1,3)');

if N_match~=sum(output(:,3))
 plot(DAQ_data_ds(ceil(gitter/frm_rate):end))
 hold all
 plot(output(1:end-ceil(gitter/frm_rate)+1,3))
 error('Check the reward')
end
disp(['Matched reward ' num2str(N_match) '/' num2str(sum(output(:,3)))])


end
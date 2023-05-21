[fname,fpath] = uigetfile('*.*'); if fname == 0, return;end
% cd(fpath)
mov = vm(fpath);

[roi,intens]=clicky_faster(mov);



dt=1; %ms


n_frame_ignore = 100;

filt_order = 4;
low_f_lim = 100;
high_f_lim = inf;

t_total = length(intens)*dt*1e-3;
t_bin = 50e-3;
n_bin = t_total/t_bin;
t_edge = 0:t_bin:t_total;
x_t = t_edge(1:end-1)+t_bin/2;

fr = zeros(length(x_t),length(roi));
for i=1:length(roi)
intens_filt = [butterworth_filt(fft_clean(intens(n_frame_ignore+1:end,i)),...
    filt_order,[low_f_lim high_f_lim],dt/1e-3)];
spk_t = spikefindhyst2(-intens_filt,4.5*median(abs(squeeze(intens_filt))),3*median(abs(intens_filt)),1);

% figure(50)
% clf
% plot(fft_clean(intens));hold on;plot(spk_t+100,ones(1,length(spk_t))*(min(intens(101:end))-10),'*')

[cnt,edge]=histcounts(spk_t*dt*1e-3,t_edge);

fr(:,i) = sgolayfilt(cnt/t_bin,3,floor(length(cnt)/5/2)*2+1);
end

figure(21);clf
plot(x_t,fr)
xlabel('Time (s)');
ylabel('Firing Rate (Hz)')
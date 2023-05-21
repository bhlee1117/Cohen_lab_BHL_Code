%% multiple file mean
% fpath = uigetdir();
% cd(fpath)
%%
fpath = 'R:\';
folder_prefix = {
    'M-YQ0201-7_FOV7_d1_c1_.2s_.5c_1.2V_t150ms_wt60s_c_rep'
    'M-YQ0201-7_FOV7_d1_c1_.2s_.5c_1.2V_t150ms_wt60s_cw_rep'
    'M-YQ0201-7_FOV7_d1_c1_.2s_.5c_1.2V_t150ms_wt60s_sw_rep'
    };
    
    
figure(2423);clf
intens_all = {};
subplot(3,3,1)
for jj = 1:length(folder_prefix)
flist = dir([fpath '\*' folder_prefix{jj} '*']);

intens_all = [];
spk_t={};
for ii=1:length(flist)
    mov = vm([fpath '\' flist(ii).name]);
    if ii==1&&jj==1
        [roi,~] = clicky_faster(mov);close
    end
    intens = apply_clicky_faster(roi,mov,'no');
    
    n_frame_ignore = 100;
    dt = 1e-3;
    filt_order = 4;
    low_f_lim = 50;
    high_f_lim = inf;
    intens_filt = butterworth_filt(intens(n_frame_ignore+1:end),...
    filt_order,[low_f_lim high_f_lim],1/dt);
    spk_t{ii} = spikefindhyst2(-intens_filt,4.5*median(abs(squeeze(intens_filt))),3*median(abs(intens_filt)),1)+n_frame_ignore;
    intens(1:50)=nan;
    intens_all = [intens_all,butterworth_filt(mat2gray(intens),4,[50,inf],1/dt)];
end

subplot(6,3,jj)
% % pat_intens = waveform_readout(fullfile(fpath,flist(ii).name),3);
% plot((1:length(pat_intens))*dt,pat_intens)
% xlim([0 1.5])
% ylim([0 10])
% y_lim = get(gca,'ylim');
% hold on;
% plot([.525 .525],ylim,'k--')
% xlim([0 1.5])
% ylabel('Blue LED Voltage')

subplot(6,3,[3 6]+jj)
intens_all(1:50,:)=nan;
plot((1:length(intens_all))*dt,intens_all+(0:size(intens_all,2)-1))
y_lim = get(gca,'ylim');
hold on;
plot([.55 .55],ylim,'k--')
xlim([0 1.5])

% subplot(6,3,9+jj)
% raster_plot(spk_t,1e-3,'k')
% xlim([0 1.5])
% y_lim = get(gca,'ylim');
% hold on;
% plot([.525 .525],ylim,'k--')
% ylabel('Trial No.')

subplot(6,3,[12 15]+jj)
t_total = length(intens)*dt;
t_bin = 50e-3;
n_bin = t_total/t_bin;
t_edge = 0:t_bin:t_total;
x_t = t_edge(1:end-1)+t_bin/2;
[cnt,edge]=histcounts(cell2mat(spk_t)*dt,t_edge);
histogram('binedges',t_edge,'bincounts',cnt/t_bin/length(flist))
xlim([0 1.5])
ylim([0 80])
y_lim = get(gca,'ylim');
hold on;
plot([.55 .55],ylim,'k--')
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')
end

figure(2423)
set(findobj(gcf,'type','axes'),'fontsize',14)

% saveas(gcf,[folder_prefix{1}(1:end-5) '.fig'])
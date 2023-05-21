%% multiple file mean
fpath = uigetdir();
% cd(fpath)
%%
folder_prefix = {
    'M-YQ0201-7_FOV7_d1_c1_wt60s_c_only_V_inc*_rep2'
    'M-YQ0201-7_FOV7_d1_c1_wt60s_w_only_V_inc*_rep2'
    };
    
    
figure(2423);clf
% intens_all = {};
% subplot(3,2,1)
for jj = 1:length(folder_prefix)
flist = dir([fpath '\*' folder_prefix{jj} '*']);

intens_all = [];
pat_intens_all = [];
spk_t={};
color_order = colormap(winter(length(flist)));
h_LED_trace = line(nan,nan);
h_fl_trace = line(nan,nan);

for ii=1:length(flist)
    mov = vm([fpath '\' flist(ii).name]);
    if ii==1&&jj==1
        [roi,~] = clicky_faster(mov);close
    end
    intens = apply_clicky_faster(roi,mov,'no');
    
    n_frame_ignore = 100;
    dt = 1e-3;
    filt_order = 4;
    low_f_lim = 25;
    high_f_lim = inf;
    intens_filt = butterworth_filt(fft_clean(intens(n_frame_ignore+1:end)),...
    filt_order,[low_f_lim high_f_lim],1/dt);
    spk_t{ii} = spikefindhyst2(-intens_filt,4.5*median(abs(squeeze(intens_filt))),3*median(abs(intens_filt)),1)+n_frame_ignore;
    intens(1:50)=nan;
    intens_all = [intens_all,mat2gray(intens)];
%     pat_intens = waveform_readout(fullfile(fpath,flist(ii).name),3);
%     pat_intens_all = [pat_intens_all pat_intens(:,jj)];
    

end

% subplot(6,2,jj)    
% h_LED_trace=line((1:length(pat_intens))*dt,pat_intens_all);
% xlim([0 1.5])
% ylim([0 max(pat_intens_all,[],'all')])
% xlim([0 1.5])
% ylabel('Blue LED Voltage')
% set(h_LED_trace,{'Color'},num2cell(color_order,2))

subplot(1,2,jj)
intens_all(1:50,:)=nan;
h_fl_trace = plot((1:length(intens_all))*dt,intens_all+(0:size(intens_all,2)-1));
xlim([0 1.5])
set(h_fl_trace,{'Color'},num2cell(color_order,2))

% subplot(6,2,6+jj)
% raster_plot(spk_t,dt,color_order)
% xlim([0 1.5])
% ylabel('Trial No.')
% xlabel('Time (s)')


% subplot(6,2,[8 10]+jj)
% t_on = find(pat_intens(:,jj))*dt;
% t_total = length(intens)*dt;
% t_edge = [0 t_on(1)-dt t_on(end) t_total];
% [cnt,edge]= cellfun(@(spk_t) histcounts(spk_t*dt,t_edge),spk_t,'uniformoutput',false);
% F_t_on = cellfun(@(cnt,edge) cnt(2)/(edge(3)-edge(2)),cnt,edge);
% plot(max(pat_intens_all,[],1),F_t_on)
% % xlim([0 1.5])
% ylim([0 50])
% xlabel('LED Voltage')
% ylabel('Firing Rate (Hz)')
end

figure(2423)
set(findobj(gcf,'type','axes'),'fontsize',14)

% saveas(gcf,[folder_prefix{1}(1:end-5) '.fig'])
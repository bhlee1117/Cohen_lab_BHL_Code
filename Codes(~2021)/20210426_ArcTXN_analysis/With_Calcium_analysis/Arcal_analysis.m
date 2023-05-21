%% This code is written for analyzing Calcium imaging + Arc-TXN population
% MODIFICATION HISTORY :
%           Written by Byung Hun Lee, Deptartment of Physics and Astronomy,
%           Seoul National University, 2020/08/14
%           Modified C_df, 2020/09/28
%           Add baseline correction from both starting part and end part,
%           2020/10/1
%% Load data
clear
%load('E:\BACKUP\대학원\연구실\MY_Projects\In_vivo_imaging\VR\Data\20200928_Data.mat')
for i=1:6
    [fnm{1,i} pth{1,i}]=uigetfile('Multiselect','On');
end
for i=1:6
    Full_result{i}=combine_var(fullfile(pth{1,i},fnm{1,i}));
end
%% Generate the identified cell matrix
%dat=show_cdf_arc(Full_result,[1:6]);
dat=deconvolution_batch(dat);
dat=fourier_batch(dat);
%% show cal trace examples
m=2;
%plot_contours_multipleday(Full_result{m}.SP_identified,Full_result{m}.options,...
%                        [Full_result{m}.cell_list(:,1:2)/2 Full_result{m}.cell_list(:,3)],Full_result{m}.list_identified  ,Full_result{m}.regist.argz,[1:3])
%plot_contours_multipleday_im(imR,imG,identified_SpC,options,cell_list,identified_list,argz,doi)
m=6; day=1;
rng=find(dat{m}.Arc_post_ref(:,2)==2); % show the trace of ~pre%~post neurons.
% %rng=[1:100]';
 show_transient_peaks(dat,m,day,rng(29)',1/30,[0 0 0])
%% Place field Spatial correlation
clear P
PV_corr=[]; sw=0; %PV correlation between Arc TXN - No TXN : sw = 1;
post_ref=1; PC_day=[]; t_mouse=0; groups={[1 1],[1 2]}; image_on=1;
xtick={'A1-Arc^-','A1-Arc^+'};
if sw
    P{1}=[]; P{2}=[]; 
    for i=[1:6]
        PV_corr=plot_place_field(dat,i,1,sw,image_on);
        %PV_corr=plot_place_field2(dat,i,1,sw,image_on);
        if t_mouse
            P{i,1}=[PV_corr{1,1}]; P{i,2}=[PV_corr{1,2}]; else% No TXN vs TXN
            P{1}=[P{1};PV_corr{1,1}]; P{2}=[P{2};PV_corr{1,2}]; end
    end
    p=plot_PV_corr(P,{'No TXN','TXN'},sw,t_mouse);
else
    for i=[1:6]
        PV_corr=[PV_corr; plot_place_field(dat,i,1,sw,image_on)]; end   % Day 1, 2, 3
    plot_PV_corr(PV_corr(:,1:2),{'A - A','A - B'},sw)
end
%plot_joyplot_place(place_cdf,ll)
%[SI lambda_p]=plot_spatial_information(dat,20,Full_result,groups,post_ref,PC_day,t_mouse,sw,5,xtick);
%% Fraction TXN
clear compare compare_ch
xtick={'Day 1 Pre VR','Day 1 Post VR','Day 2 Pre VR','Day 2 Post VR','Day 3 Pre VR','Day 3 Post VR'};
cmap=[repmat([0.1 0.2 0.2;0.2 0.8 0.8],2,1);0.3 0.2 0.1 ;0.8 0.5 0.3]; cmap_pie=[0.5 0.5 0.5;distinguishable_colors(3)];
[cond_mat chan_mat Arcclass]=cal_cond_mat(Full_result,[2 4]);
[Fraction TXN_frac]=plot_ArcTxn(Arcclass,0,[repmat([0.1 0.2 0.2;0 0.8 0.8],2,1);0.3 0.2 0;0.9 0.4 0],50,xtick);
group=1; class=1; ch_class=1; %1: Overlap ch, 2: Reactivation Ch, Reactivation Neg.
color=[0.2 0.8 0.8; 0.5 0.3 0.5];
[p_mat_ch p_value_ch compare_ch]=sign_matrix_chance(cond_mat,chan_mat,group,class,ch_class,[6]); %group,class,ch_class
plot_cmpr_int(compare_ch,[1 2;1 3;2 3],[0.2 0.8 0.8;0.1 0.2 0.2],class,0,30,{'A1{\bf\cap}A2','A1{\bf\cap}B','A2{\bf\cap}B'},1)
% 1:sw, 10: y_lim, Chance
Venn_data=plot_venn_Arcal(Arcclass,1,[2 4]);
%% Plot firing rate graph, No TXN vs TXN
cat=5; sw=0; post_ref=1; t_mouse=1;% sw: 1=different figures for each group., 0 : diff. fig. for each days.
tick={'Neg-Neg','Neg-Pos','Pos-Neg','Pos-Pos'};
%acf=Plot_acf(dat,xtick,1,0.5); %just to see Arc- vs Arc+
% [FR SP]=plot_FR_bar(dat, {'Arc^-','Arc^+'},[0.1 0.1 0.1; 1 0 0],8,cat,t_mouse);
% plot_percent_TXN(FR,6,cat,t_mouse)
%[HR]=plot_histo_bin_trace_nogroup(dat,tick,1,t_mouse);
%F=Plot_four(dat,tick,0.2,15,post_ref);
[IBI IBI_pooled Pow_spec_pooled Brst_pooled theta_Brst_pooled]=burst_detection(dat,0,2.5,post_ref);
%% Plot firing rate graph, in groups
groups={[1 1 2 2 3 2],[1 2 2 2 3 2]}; %[Day cond Day cond]. Rows for number of groups. Columns/2 : number of conditions.
cat=5; sw=1; post_ref=1; t_mouse=0;% sw: 1=different figures for each group., 0 : diff. fig. for each days.
% post_ref : only consider post VR (ex, ~pre & post , pre & post -> TXN)
% t_mouse  : calculate p-value by mouse(=1) or by pooled neurons(=0)
if sw==1
    xtick={'Day 1','Day 2','Day 3'}; else xtick={'+++','-++'}; end
%fourier=Plot_fft_cond(dat,groups,xtick,0.1,15,sw,post_ref,t_mouse);
%[HR]=plot_histo_bin_trace(dat,groups,xtick,sw,1,1,t_mouse);
%plot_freq_power(HR,[3:8],[1 2],'ranksum')
%FR2=plot_FR_Cond(dat,groups,xtick,4.5,sw,cat,post_ref,t_mouse);
%[MFC fr_mat]=plot_MFC_Cond(dat,groups,xtick,1,sw);
%[IBI_cond thet_cond p_value]=burst_detection_cond(IBI,dat,groups,sw,1,0.1,t_mouse,xtick,post_ref);
%acf=Plot_acf_cond(dat,groups,xtick,1,0.5,sw,[],post_ref); % show groups day by day
%[p M_inform lambda]=plot_inform(dat,post_ref,1,1) % day_sw = 1 : A1,A2 ; A3 ,  0 = A2; A3.
%[p M_inform lambda]=plot_inform_theta(dat,IBI,post_ref,0.5,1) % day_sw = 1 : A1,A2 ; A3 ,  0 = A2; A3.
%[norm_d polar_d]=plot_3d_FR(dat,cat,post_ref,[1 2 3],pi/8);
 %[norm_d polar_d]=plot_3d_FR_theta(dat,IBI,post_ref,[1 2 3],pi/8);
%% Synchronization
tic;
groups={[1 1 2 1],[1 2 2 1],[1 1 2 2],[1 2 2 2]}; sprt=[]; % no binning
groups_pc={[1 0],[1 1]};
xtick={'A2 Arc^-','A2 Arc^+'}; post_ref=1;
g_show=[1 1;2 2]; t_mouse=0;
[corr_mat corr_mat_pooled group_div pc p_val]=cal_plot_synchro_cond(dat,groups,sprt,post_ref);
%[corr_mat_pc corr_mat_pooled_pc group_div_pc p_val_pc]=cal_plot_synchro_cond_pc(dat,groups_pc,sprt);

d=plot_syncro(corr_mat,corr_mat_pooled,group_div,sprt,'x.tif',1,2,1,xtick,g_show,t_mouse);
% fnm=['\\Neurobiophysics\Byunghun_Lee\Reports\20210105_Report자료\m6d3'];
% export_to_Gephi(corr_mat,p_val,6,3,0.1,group_div,pc,[2 1],fnm) 
toc;
%% plot graph properties

repd=[1 2 3]; t_mouse=0; sw=1; dot_sw=0; group_show=[1];
if sw
ticks={'Day 1','Day 2','Day 3'}; 
else ticks={'--','+-','++','-+'};
end
for corr_th=[0.1]
    [centrality cc_dist estrada_index Mo]=cal_parms_graphs(corr_mat,group_div,corr_th,p_val);
   [centrality_pc cc_dist_pc estrada_index_pc Mo_pc]=cal_parms_graphs(corr_mat_pc,group_div_pc,corr_th,p_val_pc);
    p_value=plot_degrees(centrality,repd,t_mouse,sw,dot_sw,group_show,ticks);
     plot_cluster_coeff(cc_dist,repd,t_mouse,sw,dot_sw,group_show,ticks)
    plot_errorbar([1 2],{Mo.modular(:,2) Mo_pc.modular(:,2)},'ttest2',0.5,'Modularity',{'A2 Arc^{+/-}','A1 PC/nonPC'},1)
    %plot_errorbar([1 2],[Mo.modular_norm(:,2) Mo_pc.modular_norm(:,2)],'ttest',0.5,'Normalized Modularity',{'A2 Arc^{+/-}','A1 PC/nonPC'},1)
    xtickangle(45)
end
% sw rep_mouse rep_day xtick
%% Event interval
t_mouse=1;
groups={[1 1],[1 2]}; bins=[0:2/30:1]; ticks={'D1 Arc^-','D1 Arc^+'};
I=cal_spike_interval(dat,groups,0);
%plot_hist_interval(I,[1 2],bins,ticks,7,'Probability')
acf=Plot_acf_cond(dat,groups,ticks,1,0.5,0,1/30,1);
%plot_population_vector(dat,groups,5,5,1)
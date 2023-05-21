clear;
[fnmLR pthLR]=uigetfile('*.mat','[Low Rank] Select analyzed calcium data','Multiselect','on');
[fnmHR pthHR]=uigetfile('*.mat','[High Rank] Select analyzed calcium data','Multiselect','on');
timescale=1.089; %Sec
meet_time=180; %Sec
%%
for i=1:size(fnmLR,2)
    load([pthLR fnmLR{1,i}])
    Astro_cal.DFF{1,i}=full(C_df);
    [Coor,json_file] = plot_contours(A_or,Cn,options,1); 
    colormap('gray')
    Astro_cal.F{1,i}=getframe;
end

for i=1:size(fnmHR,2)
    load([pthHR fnmHR{1,i}])
    Astro_cal.DFF{2,i}=full(C_df);
    [Coor,json_file] = plot_contours(A_or,Cn,options,1); 
    colormap('gray')
    Astro_cal.F{2,i}=getframe;
    
end
close all
%% plot statistics 
stim_time=180; %sec
plot_dff_astro(Astro_cal.DFF(:,1),timescale,stim_time,{'cage1-mouse5','cage1-mouse3'},[0.8 0.2 0.8;0.2 0.8 0.8;0.2 0.8 0.8;0.4 0.7 0.7])
plot_dff_astro(Astro_cal.DFF(:,2),timescale,stim_time,{'cage2-mouse2','cage2-mouse1'},[0.8 0.2 0.8;0.2 0.8 0.8;0.2 0.8 0.8;0.4 0.7 0.7])
p_value=plot_dff_bar(Astro_cal.DFF,stim_time,timescale,{'High','Low','High','Low'},[0.4 0.2 0.4;0.2 0.4 0.4;1 0.4 1;0.4 1 1]);
%%
cal_im_show(Astro_cal.DFF,2,1,timescale)
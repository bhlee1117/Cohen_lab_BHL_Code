clear; clc;
Ca_path='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Updates/2024/20240408_Movs_Figs/BHLm125_FOV1_Ca_100ms98-2.tif';
V_path='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Updates/2024/20240408_Movs_Figs/BHLm125_FOV1_V2_100ms_30mW-1.tif';
Ca_mov=readtiff(Ca_path);
V_img=readtiff(V_path);
%%
Ca_avgImg=max(Ca_mov,[],3);
[RegImg,tformReg]=imReg(Ca_avgImg,V_img);
%%
figure; clf;
cmap=[0 0.5 0.25;1 0.4 0];
imshow2(im_merge(cat(3,RegImg,V_img),cmap),[])

%%
[tr1 poly1]=polyLineKymo3(Ca_mov,10,10);
[tr2 poly2]=polyLineKymo3(Ca_mov,10,10);
%%
figure; clf;
ax1=nexttile([1 1]);
imagesc(rescale2(tr1',2))
colormap(turbo)
ax2=nexttile([1 1]);
imagesc(rescale2(tr2',2))
colormap(turbo)
linkaxes([ax1 ax2],'x')
plot(tr2)
l=plot(tr2-median(tr2,1));
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(26),2))
axis off
axis on
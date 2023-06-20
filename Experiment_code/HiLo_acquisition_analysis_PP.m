
% HiLo acquisition
%Configure_Controller(Stage);
Stage=xx.GetDevice("Scientifica_Stage_Controller");
  
Stage.Update_Current_Position_Microns();
last_z=Stage.z
  
Stage.Update_Current_Position_Microns();
first_z=Stage.z

diff_z=last_z-first_z

Acquire_1P_HiLo_ZStack(xx,600,60,'HiLo_Cell2',1)



%% HiLo data loading
cameradata=load(fullfile('output.mat'));
nrow=cameradata.output{1, 2}(1);ncol=cameradata.output{1,2}(2);


% Read the data
dat = readBinMov('frames.bin',nrow,ncol);

[ySize, xSize, nFrames] = size(dat);

% matching the data to the camera view
for j=1:size(dat,3);
    dat(:,:,j) = imrotate(flip(dat(:,:,j),2),90);
end
% 
% avgImg = mean(dat,3);
% imshow(avgImg,[])
% 
% % Max Project of WF
% figure;imshow(mat2gray(max(dat,[],3)));%caxis([0 0.5]);
% saveas(gcf,'max_project.png');
% saveas(gcf,'max_project.fig');
% 
% 
% figure;imshow(mat2gray(mean(dat,3)));%caxis([0 0.5]);
% saveas(gcf,'mean_project.png');
% saveas(gcf,'mean_project.fig');

clear HiloMov
close all
% HiLo movie process
counter = 1;
for ii = 1:2:size(dat,3)
    H=hilospeckle(dat(:,:,ii), dat(:,:,ii+1));
    HiloMov(:,:,counter) = H.HiloFinal;
    counter = counter+1;
end

% save('HiLo_im.mat','HiloMov');
% moviefixsc(HiloMov, [0, 6000])

% Max Project of HiLo processed mov
figure;imshow(mat2gray(max(HiloMov,[],3)));caxis([0 0.5]);
save HiloMov HiloMov
saveas(gcf,'max_project_Hilo.png');
saveas(gcf,'max_project_Hilo.fig');
% % ylim([150 600])
% % xlim([50 500])
% 
figure;imshow(mat2gray(mean(HiloMov,3)));caxis([0 0.5]);
saveas(gcf,'mean_project_Hilo.png');
saveas(gcf,'mean_project_Hilo.fig');
% 
% %%Hilo movie filtered
% HiloMov_filtered = zeros(nrow,ncol,size(dat,3)/2);  %define the y- and x-size
% 
% for ii = 1:size(HiloMov,3)
%     HiloMov_filtered (:,:,ii) = medfilt2(mat2gray(HiloMov(:,:,ii)));
% end
% 
% figure;imshow(mat2gray(max(HiloMov_filtered(1:end,:,:),[],3)))%,caxis([0 0.75])
% saveas(gca,'max_Hilo_medfilt.fig')
% saveas(gca,'max_Hilo_medfilt.png')
% 
% figure;imshow(mat2gray(mean(HiloMov_filtered(1:end,:,:),3)))%,caxis([0 0.75])
% saveas(gca,'mean_Hilo_medfilt.fig')
% saveas(gca,'mean_Hilo_medfilt.png')
% 
% close all
% delete(Stage.serial_com) %%deleting stage for 2P imaging

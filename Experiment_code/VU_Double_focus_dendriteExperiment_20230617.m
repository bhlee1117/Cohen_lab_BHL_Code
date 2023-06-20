% Experiment code
% z-stack
clear Zstack1 Zstack2
addpath(genpath('C:\Users\Labmember\Data\ByungHun\GitHub_code'))
Stage=xx.GetDevice("Scientifica_Stage_Controller");
cam=xx.GetDevice('Cam_Controller');
Shutters=xx.GetDevice('Shutter_Device');
dmd = xx.GetDevice('DMD_Device'); %dmd(1,1) is the ALP for the Behavior
Saveto=xx.datafolder;
   cam1 = cam(1,1);
   cam2 = cam(1,2);
%%
thickness= 100; % um
step= 1; %um
scale=100;
numslices= thickness/step;
Stage.Update_Current_Position_Microns();
last_z=Stage.z;
first_z=Stage.z + thickness*scale;


    %Checks if enable594 is ON
    if  Shutters(1, 1).State
    else
        Shutters(1, 1).State=1;
    end
tic;
for z = 1:numslices
    Stage.Update_Current_Position_Microns();
    CurrentZ=Stage.z;
    Zstack1(:,:,z) = cam1.Snap();
    Zstack2(:,:,z) = cam2.Snap();

    pause(cam1.exposuretime) % Pause to not saturate the COM
    Stage.Step_Fixed(3,step*scale);
    disp([num2str(z) 'th slice, z = ' num2str(CurrentZ/100) 'um'])
end
toc;
tt=clock; t = datetime(tt);
t.Format='ddHHmm';
save([Saveto '\' char(t) '_Zstack.mat'],'Zstack1','Zstack2','-v7.3')
%%
[fname fpath]=uigetfile();
load(fullfile(fpath,fname))
im=snap.img;
%%
DMDmask = GenerateDMD_dendrite(im,10,'skeleton',0.1);
boundaries = bwboundaries(DMDmask); MaskEdges=[];
for i=1:length(boundaries)
MaskEdges=[MaskEdges; [boundaries{i}; NaN(1,2)]];
end
transform=dmd(1,1).tform;
%transform=snap.tform;
Rfixed = imref2d([dmd(1,1).Dimensions(2) dmd(1,1).Dimensions(1)]);

newSize = [2304 2304];
rmask = imresize(DMDmask, newSize);
% apply transformation from image space to patterning device space
tmask=imwarp(rmask,transform,'OutputView',Rfixed);

dmd(1,1).Target=tmask;
dmd(1,1).Write_Static();
if  Shutters(1, 1).State
    else
    Shutters(1, 1).State=1;
end
pause(cam1.exposuretime)
CurrentIm = cam1.Snap();
Shutters(1, 1).State=0;

figure(2); clf;
ax1=nexttile([1 1]);
imshow2(im,[]); hold all;
plot(MaskEdges(:,2),MaskEdges(:,1),'r')
title(['Before masking'])
ax2=nexttile([1 1]);
imshow2(CurrentIm,[]); hold all
plot(MaskEdges(:,2),MaskEdges(:,1),'r')
title(['After masking'])
linkaxes([ax1 ax2],'xy')
tt=clock; t = datetime(tt);
t.Format='ddHHmm';
saveas(gca,[Saveto '\' char(t) 'Masking.fig'])
saveas(gca,[Saveto '\' char(t) 'Masking.png'])

%% Clicky

fpath=uigetdir;
load([fpath '/output_data.mat'])
sz=double(Device_Data{1, 3}.ROI([2 4]));
mov=double(readBinMov([fpath '/frames1.bin'],sz(2),sz(1)));
%Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
DAQ_rate=Device_Data{1, 2}.buffered_tasks(1, 1).rate;
CamDAQ_rate=Device_Data{1, 2}.Counter_Inputs.rate;
CamTrig=Device_Data{1, 2}.Counter_Inputs.data;
CamTrig=find(CamTrig(2:end)-CamTrig(1:end-1)>0);
Frm_rate=(CamTrig(2)-CamTrig(1))/CamDAQ_rate;

% t=[1:size(mov,3)]/CamDAQ_rate;
% plot(t,mat2gray(squeeze(mean(mov,[1 2]))),[1:length(Blue)]/DAQ_rate,mat2gray(Blue),'--')
% saveas(gca,[fpath 'DAQvsF.fig'])

clicky(mov)
saveas(gca,fullfile(fpath,'Clicky_pre.fig'))
saveas(gca,fullfile(fpath,'Clicky_pre.png'))
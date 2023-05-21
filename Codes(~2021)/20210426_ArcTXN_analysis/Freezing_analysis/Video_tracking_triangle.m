%%

%     Website : http://neurobiophysics.snu.ac.kr/
%
% This code traces mouse trajectory in tri-angle shaped context.

% INPUTS 

% Image stack
% 

% OUTPUTS
% Filtered image stack

% MODIFICATION HISTORY : 
%           Written by Byung Hun Lee, Deptartment of Physics and Astronomy, Seoul National University, 6/11/2018.

clear, clc, close all
pth='H:\Image_data\Mouse_behavior\20190220_OE222325\';
fnm='OE25_Ret.MOV';
v=VideoReader([pth fnm]);
fps=4;
maxtime=180;
tail=200;
width=76;
freez_time=0.5;%sec
%%
figure(1)
vim=readFrame(v);
image(vim);
poly=roipoly;
[row col]=find(poly==1);
minr=min(row); maxr=max(row);
minc=min(col); maxc=max(col);
pix_scale=width/(maxr-minr);
close(figure(1))
%%
detector=vision.ForegroundDetector(...
    'NumGaussians',2,'NumTrainingFrames',500,'MinimumBackgroundRatio',0.9);
blobAnalyser=vision.BlobAnalysis('BoundingBoxOutputPort',true,'AreaOutputPort',true,...
    'CentroidOutputPort',true,'MinimumBlobArea',400);
%%
clear track ave_disp disp
frameNumber=0;
aviobj = VideoWriter([pth,fnm,'_track.avi']);
open(aviobj);
g=1;
freez_rate=0;
for i=1:maxtime*fps
     v.CurrentTime = i/fps;
     im=readFrame(v);
     im(~poly)=255;
    frame=imcrop(im,[round(minc) round(minr) round(maxc-minc) round(maxr-minr)]);
    frame2=rgb2gray(frame);
    mask=frame2<50;
    [area,centroid,bboxes]=blobAnalyser.step(mask);
    track(i,1:size(centroid,2))=centroid(1,1:end);
    frame_markers=insertMarker(frame,centroid);
    if i>1
    disp(i,1)=pix_scale*sqrt((track(i,1)-track(i-1,1))^2+(track(i,2)-track(i-1,2))^2);
    end
    figure(1)
    subplot(4,1,1:3)
    axis equal
    imagesc(frame_markers)
    xlim([1 size(frame_markers,2)])
    ylim([1 size(frame_markers,1)])
    hold all
    text(size(frame2,2)-100,20,[num2str(round(i/fps)) 'sec'],'color','black','fontweight','bold')
    
    if i>tail
    plot(track(i-tail:i,1),track(i-tail:i,2),'g','linewidth',1.5)
    else
    plot(track(1:i,1),track(1:i,2),'g','linewidth',1.5)
    end
    
    


if mod(i,freez_time*fps)==0
    
    ave_disp(g,1)=mean(disp(i-freez_time*fps+1:i,1));
    subplot(4,1,4)
    plot([freez_time*fps:freez_time*fps:g*freez_time*fps],ave_disp,'r','linewidth',2)
    xlim([0,maxtime*fps])
    ylim([0,10])
    hold all
%    freez_rate=sum(ave_disp<freez_thr)/size(ave_disp,1);
    
    
        g=g+1;
end
    F=figure(1);
writeVideo(aviobj,getframe(F));
if mod(i,100)==0
    close(figure(1))
    figure(1)    
end
end
 close(aviobj);
 close(figure(1))
save([pth,fnm,'.mat'],'track','disp','ave_disp');
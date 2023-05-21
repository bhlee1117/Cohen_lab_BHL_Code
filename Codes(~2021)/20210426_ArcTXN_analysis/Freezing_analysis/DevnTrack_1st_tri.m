%%
clear, clc, close all
pth='C:\Users\Administrator\Documents\카카오톡 받은 파일\';
fnm='KakaoTalk_20200819_145932883.mp4';
v=VideoReader([pth fnm]);
fps=30;
maxtime=180;
tail=200;
width=76;
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
imtool(rgb2gray(imcrop(vim,[minc minr maxc-minc maxr-minr])))
thr=input('input the threshold value\n');
%%

% aviobj = VideoWriter([pth,fnm,'_devtrack.avi']);
% open(aviobj);
clear fz

for i=1:maxtime*fps
 v.CurrentTime = i/fps;
 im=readFrame(vf);
 im(~poly)=255;
frame=rgb2gray(imcrop(im,[round(minc) round(minr) round(maxc-minc) round(maxr-minr)]));
binary_frame=frame<thr;
 v.CurrentTime = (i+1)/fps;
 im=readFrame(v);
 im(~poly)=255;
frame2=rgb2gray(imcrop(im,[round(minc) round(minr) round(maxc-minc) round(maxr-minr)]));
binary_frame2=frame2<thr;
clear dev
dev=abs(binary_frame2-binary_frame);
sumdev(i,1)=pix_scale^2*sum(sum(dev));

    mask=frame2<thr;
    [area,centroid,bboxes]=blobAnalyser.step(mask);
    if isempty(centroid)
        centroid=cent_temp;
    else
        cent_temp=centroid;
    end
    track(i,1:size(centroid,2))=centroid(1,1:end);
    frame_markers=insertMarker(frame,centroid);
    if i>1
    disp(i,1)=pix_scale*sqrt((track(i,1)-track(i-1,1))^2+(track(i,2)-track(i-1,2))^2);
    end
%     subplot(4,1,1:3)
%     imagesc(dev)
%     colormap('pink')
%     axis equal
%     hold all
%     plot(centroid(1,1),centroid(1,2),'g+')
%     
%     text(size(frame2,2)-50,20,[num2str(round(i/fps)) 'sec'],'color','w','fontweight','bold')
%     
%     if i>tail
%     plot(track(i-tail:i,1),track(i-tail:i,2),'g','linewidth',1.5)
%     else
%     plot(track(1:i,1),track(1:i,2),'g','linewidth',1.5)
%     end
%     
% subplot(3,1,3)
% plot([1/fps:1/fps:i/fps],sumdev(1:i,1)','r')
% hold all
% ylim([0 max(sumdev)])
% xlim([0 maxtime])

% F=figure(1);
% writeVideo(aviobj,getframe(F));
% if mod(i,50)==0
%     close(figure(1))
% end
end
%  close(figure(1));
 save([pth,fnm,'.mat'],'sumdev');
%  close(aviobj);
 

 
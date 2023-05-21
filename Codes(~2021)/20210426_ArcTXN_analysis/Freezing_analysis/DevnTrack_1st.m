%%
clear, clc, close all
% pth='H:\Image_data\Mouse_behavior\20200305_POE12_ROE4546\';
% fnm='POE10_Ret4.MOV';
[fnm pth]=uigetfile('*.MOV','Multiselect','On');
if ischar(fnm)
    fnm={fnm};
end
%%
sw_video=false;
fps=4;
maxtime=180;
tail=200;
width=50;
thr=50;


figure(1)
for p=1:length(fnm)
v=VideoReader([pth fnm{1,p}]);
vim=readFrame(v);
image(vim);
poly{p}=roipoly;
[row col]=find(poly{p}==1);
minr(p,1)=min(row); maxr(p,1)=max(row);
minc(p,1)=min(col); maxc(p,1)=max(col);
pix_scale(p,1)=width/(maxr(p,1)-minr(p,1));
close(figure(1))
end
%%
detector=vision.ForegroundDetector(...
    'NumGaussians',2,'NumTrainingFrames',500,'MinimumBackgroundRatio',0.9);
blobAnalyser=vision.BlobAnalysis('BoundingBoxOutputPort',true,'AreaOutputPort',true,...
    'CentroidOutputPort',true,'MinimumBlobArea',400);

%%%
% imtool(rgb2gray(imcrop(vim,[minc minr maxc-minc maxr-minr])))
% thr=input('input the threshold value\n');
%%
for p=1:length(fnm)
    if sw_video
aviobj = VideoWriter([pth,fnm{1,p},'_devtrack.avi']);
open(aviobj);
    end
clear fz
v=VideoReader([pth fnm{1,p}]);
for i=1:maxtime*fps
 v.CurrentTime = i/fps;
frame=rgb2gray(imcrop(readFrame(v),[round(minc(p,1)) round(minr(p,1)) round(maxc(p,1)-minc(p,1)) round(maxr(p,1)-minr(p,1))]));
binary_frame=frame<thr;
 v.CurrentTime = (i+1)/fps;
frame2=rgb2gray(imcrop(readFrame(v),[round(minc(p,1)) round(minr(p,1)) round(maxc(p,1)-minc(p,1)) round(maxr(p,1)-minr(p,1))]));
binary_frame2=frame2<thr;
clear dev
dev=abs(binary_frame2-binary_frame);
sumdev(i,1)=pix_scale(p,1)^2*sum(sum(dev));

    mask=frame2<50;
    [area,centroid,bboxes]=blobAnalyser.step(mask);
    if isempty(centroid)
        centroid=cent_temp;
    else
        cent_temp=centroid;
    end
    track(i,1:size(centroid,2))=centroid(1,1:end);
    frame_markers=insertMarker(frame,centroid);
    if i>1
    disp(i,1)=pix_scale(p,1)*sqrt((track(i,1)-track(i-1,1))^2+(track(i,2)-track(i-1,2))^2);
    end
    
    if sw_video
   subplot(4,1,1:3)
   imagesc(dev)
   colormap('pink')
   axis equal tight
   hold all
   plot(centroid(1,1),centroid(1,2),'g+')
    
   text(size(frame2,2)-50,20,[num2str(round(i/fps)) 'sec'],'color','w','fontweight','bold')
    
    if i>tail
    plot(track(i-tail:i,1),track(i-tail:i,2),'g','linewidth',1.5)
    else
    plot(track(1:i,1),track(1:i,2),'g','linewidth',1.5)
    end
    
subplot(3,1,3)
plot([1/fps:1/fps:i/fps],sumdev(1:i,1)','r')
hold all
ylim([0 max(sumdev)])
xlim([0 maxtime])

F=figure(1);
writeVideo(aviobj,getframe(F));
if mod(i,50)==0
    close(figure(1))
end
    end
end
 close(figure(1));
 %save([pth,fnm{1,p},'.mat'],'sumdev');
 save(['H:\Image_data\OlympusTPM\Invivo_images_CFC\track_behavior\',fnm{1,p},'.mat'],'sumdev','track');
 if sw_video
 close(aviobj);
 end
 end

 
clear
Path_n='H:\Image_data\Mouse_behavior\20180607-08\';
File_n='OE1_CFC.MOV';
File_R='OE1_Ret.MOV';
v=VideoReader([Path_n File_n]);
vR=VideoReader([Path_n File_R]);
v.CurrentTime = 1;
vR.CurrentTime = 1;
vim=readFrame(v);
vimR=readFrame(vR);
%%
figure(1)
image(vim);
poly=roipoly;
[row col]=find(poly==1);
minr=min(row); maxr=max(row);
minc=min(col); maxc=max(col);
close(figure(1))

figure(1)
image(vimR);
polyR=roipoly;
[rowR colR]=find(polyR==1);
minrR=min(rowR); maxrR=max(rowR);
mincR=min(colR); maxcR=max(colR);
close(figure(1))
%%
thr=input('input the threshold value\n');
aviobj = VideoWriter([Path_n,File_n,'_dev.avi']);
open(aviobj);
clear fz
maxtime=240;
for i=1:maxtime
 v.CurrentTime = i+1;
vim=readFrame(v);
grim=heaviside(double(rgb2gray(vim(minr:maxr,minc:maxc,:)))-thr);
 v.CurrentTime = i;
 vim=readFrame(v);
 grim_b=heaviside(double(rgb2gray(vim(minr:maxr,minc:maxc,:)))-thr);
dev_grim(:,:,i)=abs(grim-grim_b);
fz(i)=sum(sum(dev_grim(:,:,i)));
subplot(3,1,1:2)
figure(1)
imagesc(dev_grim(:,:,i))
axis equal
colormap('gray')
hold all
text(size(dev_grim,2)-100,20,[num2str(i) 'sec'],'color','white')
subplot(3,1,3)
plot([1:1:i],fz(:)','r')
hold all
ylim([0 max(fz)])
xlim([0 maxtime])

F=figure(1);
writeVideo(aviobj,getframe(F));
if mod(i,50)==0
    close(figure(1))
end
end
 close(figure(1));
 save([Path_n,File_n,'.mat'],'dev_grim');
 close(aviobj);
 
 %%
 aviobj = VideoWriter([Path_n,File_R,'_dev.avi']);
open(aviobj);
clear fz
maxtime=180;
for i=1:maxtime
 vR.CurrentTime = i+1;
vimR=readFrame(vR);
grimR=heaviside(double(rgb2gray(vimR(minrR:maxrR,mincR:maxcR,:)))-thr);
 vR.CurrentTime = i;
 vimR=readFrame(vR);
 grim_bR=heaviside(double(rgb2gray(vimR(minrR:maxrR,mincR:maxcR,:)))-thr);
dev_grimR(:,:,i)=abs(grimR-grim_bR);
fzR(i)=sum(sum(dev_grimR(:,:,i)));
subplot(3,1,1:2)
figure(1)
imagesc(dev_grimR(:,:,i))
axis equal
colormap('gray')
hold all
text(size(dev_grimR,2)-100,20,[num2str(i) 'sec'],'color','white')
subplot(3,1,3)
plot([1:1:i],fzR(:)','r')
hold all
ylim([0 max(fzR)])
xlim([0 maxtime])

F=figure(1);
writeVideo(aviobj,getframe(F));
if mod(i,50)==0
    close(figure(1))
end
end
 close(figure(1));
 save([Path_n,File_R,'.mat'],'dev_grimR');
 close(aviobj);
 
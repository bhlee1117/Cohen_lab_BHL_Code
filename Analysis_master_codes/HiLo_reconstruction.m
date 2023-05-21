%%
% 2022/08/17 Hilo reconstruction


%% Load images
%clear
[fpth]=uigetfile_n_dir();
tic;
clear Hilo_im im Hilo_im_histmat
for i=2:length(fpth)
    mov=vm(fpth{i});
    [inten order]=sort(sum(tovec(double(mov)),1),'descend');
    order(inten<10^4)=[];
    s=[1:floor(length(order)/3):length(order)];
    %im(:,:,1:3,i)=double(mov(:,:,order(s(1:3))));
    im(:,:,1:3,i)=double(mov(:,:,[1:3]));
end
toc;
%%

for i=1:size(im,4)
H=hilospeckle(squeeze(im(:,:,1,i)),squeeze(im(:,:,2,i)));
H2=hilospeckle(squeeze(im(:,:,1,i)),squeeze(im(:,:,3,i)));
H.targetThickness=3;    
H.runHiLo
H2.targetThickness=3;    
H2.runHiLo
Hilo_im(:,:,i)=mean(cat(3,H.HiloFinal,H2.HiloFinal),3);
end
Hilo_im=uint16(Hilo_im);
close all
for i=2:size(Hilo_im,3)
Hilo_im_histmat(:,:,i)=imhistmatch(Hilo_im(:,:,i),Hilo_im(:,:,i-1));
end
figure;
moviefixsc(Hilo_im_histmat)
figure;
subplot(4,1,1:3)
imshow2(max(Hilo_im_histmat,[],3),[])
subplot(4,1,4)
imshow2(max(Hilo_im_histmat,[],3),[])
imshow2(imresize(max(imrotate3(Hilo_im_histmat,90,[1 0 0]),[],3),[150 1024]),[])
%%
writeMov('Hilo_10x_oldV',Hilo_im(:,:,:),[],20,200/70,[100 30000],[],[])
save('HiLo_im','Hilo_im','Hilo_im_histmat')

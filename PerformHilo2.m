fname='frames.bin'; 
path= pwd;
% addpath('X:\Lab\Labmembers\David Wong-Campos\Sandbox\HiLo_SIM\@hilospeckle')
load('output.mat')


nRow = output{2}(1);
nCol = output{2}(2);
[mov1,nframe] = readBinMov(fullfile(path,fname),nCol,nRow);

if nRow ~= nCol
    mov1 = mov1(1:min(nCol,nRow),1:min(nCol,nRow),:);
end

mov2 = zeros(nRow,nCol,size(mov1,3));
% %Clean signal
for ii = 1:size(mov1,3)
    mov2(:,:,ii) = medfilt2(mov1(:,:,ii),[8,8]);
end
%rotate and flip to get right direction
HiloMov = zeros(nRow,nCol,size(mov1,3)/2);
% Perform HiLo

    H=hilospeckle(mov1,mov1)%, mov2(:,:,1:2:end));
    HiloMov= H.HiloFinal;
    HiloMov = double(HiloMov);
    %% Test

        %H=hilospeckle(uniform_image,structured_speckle_image);
        A = medfilt2(mov1,[16,16]);
        figure(1)
        imshow(A(1000:1600,1000:1600),[])
        figure(2)
        imshow(uniform_image,[])
        %%
H=hilospeckle(A(1000:1600,1000:1600),A(1000:1600,1000:1600));



%% Discard last pixels
HiloMov_clean = HiloMov(1:end-20,1:end-20,:);
% moviefixsc(imrotate(fliplr(HiloMov_clean),90),[0 5000])

HiloMov_clean = permute(HiloMov_clean,[2 1 3]);
% HiloMov_clean = HiloMov_clean(:,:,2:end);
moviefixsc(HiloMov_clean,[0 7000])
moviefixsc(mov2,[0 10000])

%% Comparison Hilo and WF
WF = permute(double(mov1(10:end-10,10:end-10,1:120)),[2 1 3]);
    H=hilospeckle(WF,WF)%, mov2(:,:,1:2:end));
    HiloMov= H.HiloFinal;
    HiloMov = double(HiloMov);

figure;
moviefixsc(WF)

forvid = -WF(200:1600,500:1900,:);
v = VideoWriter('HiLoMov.avi');
open(v);

for k = 1:80
   imshow(forvid(:,:,k),[])
   set(gcf,'Color','w')
   frame = getframe(gcf);
   
   writeVideo(v,frame);
end

close(v);

close(v);
%%
Z = HiloMov;
Z_index = uint8((Z - min(Z(:))) * (255 / (max(Z(:)) - min(Z(:)))));
%Z_index = Z_index(200:2000,200:2000,:);

disp('16 bit, grayscale image');
clear options;
saveastiff(uint16(single(Z_index)/255*(2^16-1)), 'Neuron_Zstack.tif');

%%
[nRow,nCol,~] = size(HiloMov_clean);

HiloMov_filtered_Side = zeros(nRow,nCol,size(mov1,3)/2);
HiloMov_filtered_FromTop = zeros(nRow,nCol,size(mov1,3)/2);

filt3 =  fspecial('gaussian', [10 10], 4);

for ii = 1:size(HiloMov_clean,3)
     HiloMov_filtered_Side (:,:,ii) =mat2gray(medfilt2(HiloMov_clean(:,:,ii)).^0.8); %square parts of the signal to enhance contrast
     HiloMov_filtered_FromTop (:,:,ii) =medfilt2(HiloMov_clean(:,:,ii));
%     HiloMov_filtered_Side (:,:,ii) = imfilter(mat2gray(HiloMov_clean(:,:,ii)), filt3, 'replicate');
end
%%
%moviefixsc(permute(HiloMov_filtered_FromTop, [2,1,3]),[0 1000])
moviefixsc(HiloMov_filtered_FromTop)
%%
%% % Make a movies WF

saveavi(vm(mat2gray(HiloMov_filtered_FromTop)),'Filtered Hilo.avi')
saveavi(vm(mat2gray(HiloMov_clean)),'HiLo.avi')

%%
v = VideoWriter('HiLoClean.avi');
v.FrameRate = 25;

Mov_rep = mat2gray(HiloMov_clean);
open(v);
counter = 1;
for k = 1:size(Mov_rep,3) 
    if counter == 1
    h=imshow(Mov_rep(:,:,k),[0 0.2])
    set(gcf,'Color','w');
    else
        set(h,'CData',Mov_rep(:,:,k))
        
    end
%     set(gca,'Position',[0 0.0 0.8 0.8]); %Gets rid of the white margins
    counter = counter+1;
    axis square
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);
%%
imshow(mat2gray(max(movi,[],3)))

%%
Side1 = -imrotate(reshape(max(movi,[],1),[1701,120,1]),90);
imagesc(Side1)

%% bare HiLo
figure;imshow(mat2gray(max(HiloMov_clean,[],3))),caxis([0 0.2])
saveas(gca,'Hilo.fig')
saveas(gca,'Hilo.png')
%% Filtered HiLo
figure;imshow(mat2gray(max(HiloMov_filtered_FromTop,[],3))),caxis([0 0.2])
saveas(gca,'Hilo_medfilt.fig')
saveas(gca,'Hilo_medfilt.png')
%% Raw video
figure;imshow(mat2gray(max(mov1(:,:,1:2:end),[],3)')),caxis([0,1])
saveas(gca,'WF_maxIntens.fig')
saveas(gca,'WF_maxIntens.png')
%%
close all
Stage_Zstep =  output{1}(3,:);
steps = mean(diff(Stage_Zstep)/10);

MovieToSlice = mat2gray(HiloMov_clean);
MovieToSlice = MovieToSlice(10:end,10:end,:);

[nRow,nCol,nframes] = size(MovieToSlice);
filt3 =  fspecial('gaussian', [3 3], 2);

Depth = 1:steps:nframe;
Transverse = (1:nRow)*6.5/25;
Side1 = -imrotate(reshape(max(MovieToSlice,[],1),[nRow,nframes,1]),90);
Side1 = imfilter(Side1, filt3, 'replicate');

figure;imagesc(Transverse,Depth,Side1),colormap gray

% caxis([-2.25 0])
ylabel('Depth ({\mu}m)')
xlabel('Transverse axis ({\mu}m)')
saveas(gca,'HiLo_SideProj1.fig')
saveas(gca,'HiLo_SideProj1.png')
% Side projection 2

Side2 = -imrotate(reshape(max(MovieToSlice,[],2),[nRow,nframes,1]),90);
Side2 = imfilter(Side2, filt3, 'replicate');

figure;imagesc(Transverse,Depth,Side2),colormap gray
% caxis([-1.25 0])
ylabel('Depth ({\mu}m)')

%
xlabel('Transverse axis ({\mu}m)')
saveas(gca,'HiLo_SideProj2.fig')
saveas(gca,'HiLo_SideProj2.png')
%%
Stage_Zstep =  output{1}(3,:);
steps = mean(diff(Stage_Zstep)/10);

MovieToSlice = mat2gray(mov1(:,:,1:2:end));
[nRow,nCol] = size(MovieToSlice);
Depth = 1:steps:nframe;
Transverse = (1:nRow)*6.5/25;
figure;imagesc(Transverse,Depth,-imrotate(reshape(max(MovieToSlice,[],1),[nRow,nframes,1]),90)),colormap gray
ylabel('Depth ({\mu}m)')
xlabel('Transverse axis ({\mu}m)')
saveas(gca,'WF_SideProj1.fig')
saveas(gca,'WF_SideProj1.png')
% Other side

figure;imagesc(Transverse,Depth,-imrotate(reshape(max(MovieToSlice,[],2),[nRow,nframes,1]),90)),colormap gray
ylabel('Depth ({\mu}m)')
xlabel('Transverse axis ({\mu}m)')
saveas(gca,'WF_SideProj2.fig')
saveas(gca,'WF_SideProj2.png')
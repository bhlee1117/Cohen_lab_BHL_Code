fname='frames.bin'; 
path= pwd;
addpath('X:\Lab\Labmembers\David Wong-Campos\Sandbox\HiLo_SIM\@hilospeckle')
load('output.mat')


nRow = output{2}(1);
nCol = output{2}(2);
[mov1,nframe] = readBinMov(fullfile(path,fname),nCol,nRow);

if nRow ~= nCol
    mov1 = mov1(1:min(nCol,nRow),1:min(nCol,nRow),:);
end


%rotate and flip to get right direction
HiloMov = zeros(nRow,nCol,size(mov1,3)/2);

% Perform HiLo

    H=hilospeckle(mov1(:,:,1:2:end), mov1(:,:,2:2:end));
    HiloMov= H.HiloFinal;
    HiloMov = double(HiloMov);
%% Discard last pixels
HiloMov_clean = double(HiloMov(1:end-20,1:end-20,:));
% moviefixsc(imrotate(fliplr(HiloMov_clean),90),[0 5000])

HiloMov_clean = permute(HiloMov_clean,[2 1 3]);
% HiloMov_clean = HiloMov_clean(:,:,2:end);
moviefixsc(HiloMov_clean,[0 2000])
%% Comparison Hilo and WF
WF = permute(double(mov1(:,:,1:2:end)),[2 1 3]);
moviefixsc(WF)
figure;moviefixsc(WF)


%%
[nRow,nCol,~] = size(HiloMov_clean);

HiloMov_filtered_Side = zeros(nRow,nCol,size(mov1,3)/2);
HiloMov_filtered_FromTop = zeros(nRow,nCol,size(mov1,3)/2);

for ii = 1:size(HiloMov_clean,3)
    HiloMov_filtered_Side (:,:,ii) =mat2gray(medfilt2(HiloMov_clean(:,:,ii)).^0.8); %square parts of the signal to enhance contrast
    HiloMov_filtered_FromTop (:,:,ii) =medfilt2(HiloMov_clean(:,:,ii));

end
%%
moviefixsc(permute(HiloMov_filtered_FromTop, [2,1,3]),[0 1000])
%%
%% % Make a movies WF

saveavi(vm(mat2gray(HiloMov_filtered_FromTop)),'Filtered Hilo.avi')
saveavi(vm(mat2gray(HiloMov_clean)),'HiLo.avi')

%%
v = VideoWriter('WF.avi');
v.FrameRate = 25;

Mov_rep = mat2gray(mov1(:,:,1:2:end));
open(v);
counter = 1;
for k = 1:size(Mov_rep,3) 
    if counter == 1
    h=imshow(Mov_rep(:,:,k),[0 0.8])
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

MovieToSlice = mat2gray(HiloMov_filtered_Side);
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
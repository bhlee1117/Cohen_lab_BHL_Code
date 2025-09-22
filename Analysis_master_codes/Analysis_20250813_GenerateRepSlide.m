clear; clc;
structurepath=['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20240730_BHLm149_lowStim/Structure/Fused_162941BHLm149_N3_05um_step.tif'];
swcpath='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20240730_BHLm149_lowStim/Structure/Fused_162941BHLm149_N3_05um_step.swc';
treepoint=load_tree(swcpath);
StructureStack=mat2gray(double(tiffreadVolume(structurepath)));
moviename=['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Group Meeting Presentation/20250815_Movs_Figs/SWCmovie'];
%%
swcpoint=[];
[Nx Ny Nstack]=size(StructureStack);
[x, y, z] = meshgrid(1:Nx, 1:Ny, 1:Nstack); % 50x50x50 grid
swcpoint(:,1) = treepoint.X; % Extract X coordinates from treepoint
swcpoint(:,2) = treepoint.Y; % Extract Y coordinates from treepoint
swcpoint(:,3) = treepoint.Z; % Extract Z coordinates from treepoint
swcpoint(:,4) = treepoint.R*1.5; % Extract Z coordinates from treepoint
swcpoint(1,4)=30;
pixelsize=0.468;
f=figure(1);
set(f, 'Color', 'k');
myVideo = VideoWriter([moviename],"MPEG-4"); %open video file
myVideo.FrameRate = 50;  %can adjust this, 5 - 10 works well for me
myVideo.Quality= 100;
open(myVideo)
showDR=[prctile(StructureStack(:),0.5) prctile(StructureStack(:),99.9)];

for i=1:Nstack
    clf;
    imshow2(StructureStack(:,:,i),showDR);
    axis tight off equal
    hold on
    text(Ny*0.63,Nx*0.78,['100 \mum'],'color',[1 1 1],'FontSize',12);
    plot(Ny*0.6+[0 100/pixelsize],[Nx Nx]*0.8,'color','w','LineWidth',2);
    text(Ny*0.18,Nx*0.23,[num2str(i*pixelsize,3) ' \mum'],'color',[1 1 1],'FontSize',12);
    xlim([220 1110]);
    ylim([360 1600])
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    pause(0.05);
end
pointshow=zeros(size(swcpoint,1),1);
for i=[Nstack:-1:1]
    clf;
    imshow2(StructureStack(:,:,i),showDR);
    axis tight off equal
    hold on
    text(Ny*0.63,Nx*0.78,['100 \mum'],'color',[1 1 1],'FontSize',12);
    plot(Ny*0.6+[0 100/pixelsize],[Nx Nx]*0.8,'color','w','LineWidth',2);
    text(Ny*0.18,Nx*0.23,[num2str(i*pixelsize,3) ' \mum'],'color',[1 1 1],'FontSize',12);
    hold all
    pointshow=((round(swcpoint(:,3)/pixelsize)+5)==i) | pointshow;
    scatter(swcpoint(find(pointshow),1)/pixelsize,swcpoint(find(pointshow),2)/pixelsize,swcpoint(find(pointshow),4)/pixelsize,[1 0.8 0],'filled');
    xlim([220 1110]);
    ylim([360 1580])
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    pause(0.05);
end

% Overlay 3D points
clf;
pointPlot=scatter3(swcpoint(:,1)/pixelsize, swcpoint(:,2)/pixelsize, swcpoint(:,3)/pixelsize, swcpoint(:,4)/pixelsize, [1 0.8 0], 'filled'); % Red filled circles
hold all;
plot(Ny*0.6+[0 100/pixelsize],[Nx Nx]*0.8,'color','r','LineWidth',2);
plot([Ny Ny]*0.6+[100/pixelsize],[Nx]*0.8+[0 -100/pixelsize],'color','g','LineWidth',2);
plot3([Ny Ny]*0.6+[100/pixelsize],[Nx Nx]*0.8+[0],[0 0]+[0 100/pixelsize],'color','b','LineWidth',2);
view([0 -90]); axis equal off;
xlim([220 1110]);
ylim([360 1580]);
zlim([0 300]);

% Animate rotation around x-axis
pointsCentered=swcpoint(:,1:3)/pixelsize;
pointsCentered=pointsCentered-median(pointsCentered,1);
% Animate rotation around x-axis using view
for az = [0:-3:-90]
    view(az, -90); % Rotate azimuth, keep elevation at 90 for x-axis rotation
    drawnow;
    pause(0.1); % Slow rotation for visibility
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end

for az = [-90:1:90 90:-1:-90]
    view(-90, az); % Rotate azimuth, keep elevation at 90 for x-axis rotation
    xlim([220 1110]);
    ylim([360 1580]);
    zlim([0 300]);
    drawnow;
    pause(0.1); % Slow rotation for visibility
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end
close(myVideo);
%close(figure(1));
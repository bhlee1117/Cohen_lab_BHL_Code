%% Compile the mex function
mex wooRaytrace.cpp

%% Test in 2D

gridSize = [8 8 1];
gridBounds = [-4 -4 -0.5 4 4 0.5];
gridImage = zeros(8,8);
Ntrials = 100;

figure;
set(gca,'XLim',[-8 8],'YLim',[-8 8],'XTick',-8:8,'YTick',-8:8);
box on;
grid on;
axis square;
xlabel('X');
ylabel('Y');
hold all;
rng(2);
for nTrial = 1:Ntrials
    lineCoord = [rand(1)*16-8 rand(1)*16-8 0 rand(1)*16-8 rand(1)*16-8 0];
    fprintf('Line coordinates:%s\n',num2str(lineCoord));
    indexes = wooRaytrace(gridSize,gridBounds,lineCoord);
    gridImage(indexes) = 1+mod(1:numel(indexes),2);
    im = imagesc([-3.5 3.5],[-3.5 3.5],gridImage);
    im.AlphaData = 0.5;
    plot([lineCoord(1) lineCoord(4)]',[lineCoord(2) lineCoord(5)]','-r*');
    waitforbuttonpress
    gridImage(indexes) = 0;
    cla;
end

%% Test in 3D

gridSize = [8 8 8];
gridBounds = [-4 -4 -4 4 4 4];
gridImage = zeros(8,8,8);
Ntrials = 100;

figure;
ha1 = subplot(1,3,1);
xlabel('X');
ylabel('Y');
ha2 = subplot(1,3,2);
xlabel('X');
ylabel('Z');
for ha = [ha1 ha2]
    axes(ha);
    hold all;
    set(ha,'XLim',[-8 8],'YLim',[-8 8],'XTick',-8:8,'YTick',-8:8);
    box on;
    grid on;
    axis square;
end
ha3 = subplot(1,3,3);
axes(ha3);
hold all;
xlabel('X');
ylabel('Y');
zlabel('Z');
view(45,45);
box on;
grid on;
axis equal;


rng(1);
for nTrial = 1:Ntrials
    lineCoord = rand(1,6)*16-8;
    fprintf('Line coordinates:%s\n',num2str(lineCoord));
    indexes = wooRaytrace(gridSize,gridBounds,lineCoord);
    gridImage(indexes) = 1+mod(1:numel(indexes),2);
    % draw two projections
    axes(ha1);
    im = imagesc([-3.5 3.5],[-3.5 3.5],max(gridImage,[],3));
    im.AlphaData = 0.5;
    plot([lineCoord(1) lineCoord(4)]',[lineCoord(2) lineCoord(5)]','-r*');
    
    axes(ha2);
    im = imagesc([-3.5 3.5],[-3.5 3.5],squeeze(max(gridImage,[],1))');
    im.AlphaData = 0.5;
    plot([lineCoord(1) lineCoord(4)]',[lineCoord(3) lineCoord(6)]','-r*');
    % 3d plot
    axes(ha3);
    for nVoxel = 1:size(indexes,1)
        [dy, dx, dz] = ind2sub([8 8 8],indexes(nVoxel,1));
        vx = -4 + [(dx-1)*1 (dx-1)*1 + 1];
        vy = -4 + [(dy-1)*1 (dy-1)*1 + 1];
        vz = -4 + [(dz-1)*1 (dz-1)*1 + 1];
        fv.vertices = [[vx(1) vy(1) vz(1)];[vx(2) vy(1) vz(1)];[vx(2) vy(2) vz(1)];[vx(1) vy(2) vz(1)]; ...
            [vx(1) vy(1) vz(2)];[vx(2) vy(1) vz(2)];[vx(2) vy(2) vz(2)];[vx(1) vy(2) vz(2)]];
        fv.faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
        h = patch(fv, 'FaceColor', 'blue', 'EdgeColor', 'white');
        set(ha3,'XLim',[-4 4],'YLim',[-4 4],'ZLim',[-4 4],'XTick',-4:4,'YTick',-4:4,'ZTick',-4:4);
        h.FaceAlpha = 0.2;
    end
    plot3([lineCoord(1) lineCoord(4)]',[lineCoord(2) lineCoord(5)]',[lineCoord(3) lineCoord(6)]','-r*');
    waitforbuttonpress
    %pause(2);
    gridImage(indexes) = 0;
    cla(ha1);
    cla(ha2);
    cla(ha3);
end

%% Test performance

Ntrials = 1000000;
gridSize = [128 128 128];
gridBounds = [-4 -4 -4 4 4 4];
lineCoords = rand(Ntrials,6) * 16 - 8;
tic
for nTrial = 1:Ntrials
    indexes = wooRaytrace(gridSize,gridBounds,lineCoords(nTrial,:));
end
toc

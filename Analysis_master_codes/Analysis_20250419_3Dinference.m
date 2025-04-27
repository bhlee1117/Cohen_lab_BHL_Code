%% Load data
clear; clc;
load('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20250421_3Dinference/bhlm113_n2_results_summary.mat')
load('/Volumes/BHL18TB_D2/20240207/150121BHLm113_N2_SomStim_SP/OP_Result.mat')
%%
figure;
imagesc(B)
graph_adj = full(graph_adj);
    [row, col] = find(graph_adj);
    edges = [row, col];

    node_positions=double(node_positions);
    figure; clf;
    for i=1:size(edges,1)
plot3([node_positions(edges(i,1), 1), node_positions(edges(i,2), 1)], ...
             [node_positions(edges(i,1), 2), node_positions(edges(i,2), 2)], ...
             [node_positions(edges(i,1), 3), node_positions(edges(i,2), 3)], ...
             'k-', 'LineWidth', 1.5); hold all
    end


    %%
TwDTr=Result.normTraces(:,11:18399+10);
bwSp=bwlabel(Result.Blue(11:18399+10)).*Result.spike(1,11:18399+10);
sp_list=[];
for b=1:max(bwSp)
sp_list=[sp_list find(bwSp==b,1)];
end
TwDSp=ind2vec(18399,sp_list,1);
denoised_voltages_STA=get_STA(denoised_voltages(:,1:18399),TwDSp,50,50);
%% Preprocessing (once)
scale_factors = [4, 1.87, 1.87];
voxel_data_um = imresize3(voxel_data, scale_factors .* size(voxel_data));

% Reorient: voxel_data (Z Y X) → (Y X Z)
[Y, X, Z] = ndgrid(1:size(voxel_data_um,2), 1:size(voxel_data_um,3), 1:size(voxel_data_um,1));
fv = isosurface(X, Y, Z, permute(voxel_data_um, [2 3 1]), 50);

% Reorder node positions: [Y X Z] → [X Y Z]
nodeXYZ = node_positions(:, [1 2 3]);

% Movie setup
nFrames = size(denoised_voltages_STA, 2);
v = VideoWriter('voltage_dendrite_rotation.mp4', 'MPEG-4');
v.FrameRate = 10;  % Adjust as needed
open(v);

% Create figure (only once)
fig = figure(99); clf;
axis equal off tight;
colormap turbo;
p = patch(fv);
set(p, 'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceAlpha', 1);
colorbar;
camlight headlight;
lighting gouraud;

% Set static limits and view
caxis([min(denoised_voltages_STA(:)), max(denoised_voltages_STA(:))]);
% Initialize view from XY plane
initialAz = 2;
initialelev = 90;

for t = 1:nFrames
    % Update color
    F = scatteredInterpolant(nodeXYZ(:,1), nodeXYZ(:,2), nodeXYZ(:,3), ...
                             denoised_voltages_STA(:, t), 'natural', 'none');
    vertexColors = F(fv.vertices(:,1), fv.vertices(:,2), fv.vertices(:,3));
    set(p, 'FaceVertexCData', vertexColors);

    % Rotate around Y-axis (update azimuth)
    az = initialAz;  % Rotate 2° per frame
    elev = initialelev + (t)*1.6;
    view(az, elev);           % azimuth, elevation

    title([num2str(t-50) ' ms'])
    drawnow;
    frame = getframe(fig);
    writeVideo(v, frame);
end

close(v);


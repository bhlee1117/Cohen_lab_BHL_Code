function boundary = semiauto_roi_test(im0)

params.equalization_cliplim.value = 0.01;
params.equalization_cliplim.on = true;

params.background_size.value = 50;
params.background_size.on = true;

params.median_size.value = 15;
params.median_size.on = true;

params.gaussian_sigma.value = 1;
params.gaussian_sigma.on = true;

params.minimum_area.value = 7^2;
params.minimum_area.on = true;

params.maximum_area.value = 20^2;
params.maximum_area.on = true;

params.minimum_signal.value = .1;
params.minimum_signal.on = true;

[label_matrix, ~] = find_regions(im0./max(im0,[],'all'), params);

boundary = bwboundaries(label_matrix);

figure
for k = 1:length(boundary)
   boundary_k = boundary{k};
   plot(boundary_k(:,2), boundary_k(:,1), 'r', 'LineWidth', 2)
   hold on
end
axis equal
end
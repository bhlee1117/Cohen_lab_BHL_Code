params.equalization_cliplim.value = 0.01;
params.equalization_cliplim.on = false;

params.background_size.value = 40;
params.background_size.on = true;

params.median_size.value = 21;
params.median_size.on = false;

params.gaussian_sigma.value = 4;
params.gaussian_sigma.on = true;

params.minimum_area.value = 7^2;
params.minimum_area.on = true;

params.maximum_area.value = 30^2;
params.maximum_area.on = true;

params.minimum_signal.value = .01;
params.minimum_signal.on = false;

[label_matrix, ~] = find_regions(im0./max(im0,[],'all'), params);

boundary = bwboundaries(label_matrix);

figure;imshow(im0,[],'initialmagnification','fit');hold on
for k = 1:length(boundary)
   boundary_k = boundary{k};
   plot(boundary_k(:,2), boundary_k(:,1), 'r', 'LineWidth', 2)

end
axis equal
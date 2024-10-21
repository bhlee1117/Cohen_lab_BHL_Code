function circular_mean=mean_angle(angles,dim)

% Step 1: Convert angles to Cartesian coordinates
x = cos(angles);
y = sin(angles);

% Step 2: Compute the mean of Cartesian components
mean_x = mean(x,dim,'omitnan');
mean_y = mean(y,dim,'omitnan');

% Step 3: Compute the circular mean
circular_mean = atan2(mean_y, mean_x);

end
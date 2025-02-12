function [mean_amplitudes std_amplitudes x_bin_centers indicies binnedData]=binning_data(Data,x_bin_edges)

% 2024.10.27 Byung Hun Lee
% Binning data, input: Data(cell variable, N x 2, matrix, 1 st column: x
% axis (axis to bin), 2nd column: measure))

x_bin_centers = mean([x_bin_edges(1:end-1); x_bin_edges(2:end)],1);

mean_amplitudes = zeros(1, length(x_bin_centers));
std_amplitudes = zeros(1, length(x_bin_centers));


for j = 1:length(x_bin_centers)
    bin_values = [];  % To collect all amplitude values in this bin
    
    for f = 1:length(Data)
        if ~isempty(Data{f})
        x = Data{f}(:,1);  % x-values
        amplitude = Data{f}(:,2);  % amplitude values
        bin_indices = (x >= x_bin_edges(j)) & (x < x_bin_edges(j + 1));
        indicies{j,f}=bin_indices;

        bin_values = [bin_values, amplitude(bin_indices)'];  % Append amplitudes in this bin
        end
    end

    if ~isempty(bin_values) 
        mean_amplitudes(j) = mean(bin_values,'omitnan'); % Calculate mean and standard deviation if there are values in this bin
        std_amplitudes(j) = std(bin_values,'omitnan');
        binnedData{j}=bin_values;
    else
        mean_amplitudes(j) = NaN;
        std_amplitudes(j) = NaN;
    end
end
end
function p_values = get_pValue(data, isPaired)
% computePairwisePValues - Calculates pairwise p-values for paired or unpaired data.
%
%   p_values = computePairwisePValues(data, isPaired)
%     - data: numeric matrix (paired) or cell array (unpaired)
%     - isPaired: logical flag
%     - p_values: symmetric matrix of pairwise p-values

if isPaired
    numGroups = size(data, 2);
else
    numGroups = numel(data);
end

p_values = nan(numGroups);

for i = 1:numGroups
    for j = i+1:numGroups
        if isPaired
            pairedData = [data(:,i), data(:,j)];
            mask = all(~isnan(pairedData), 2);
            [~, p] = ttest(pairedData(mask,1), pairedData(mask,2));
        else
            xi = data{i}(~isnan(data{i}));
            xj = data{j}(~isnan(data{j}));
            [p] = ranksum(xi, xj);
        end
        p_values(i,j) = p;
        p_values(j,i) = p;
    end
end
end

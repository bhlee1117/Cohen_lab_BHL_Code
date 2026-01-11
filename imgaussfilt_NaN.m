function filteredA2=imgaussfilt_NaN(A2,sigma)

%sigma = 1; % Standard deviation
filterSize = sigma*2+2; % Size of the filter
G = fspecial('gaussian', filterSize, sigma);

% Replace NaN values with zeros and keep track of NaN positions
filteredA2=[];
for n=1:size(A2,3)
    A=A2(:,:,n);
nanMask = isnan(A);
A(nanMask) = 0;

% Apply Gaussian filter
filteredA = imfilter(A, G, 'replicate');

% Create a normalization matrix that considers the NaN values
normalizationMatrix = imfilter(double(~nanMask), G, 'replicate');

% Adjust the filtered result to account for the NaN positions
filteredA = filteredA./ normalizationMatrix;
filteredA(nanMask) = NaN;
filteredA(filteredA==inf)=NaN;
filteredA2(:,:,n)=filteredA;
end
end
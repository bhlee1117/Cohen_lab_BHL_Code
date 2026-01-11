function MovM=ringmovGaussian(M,sigma)

MovM=repmat(M,1,3);

% Build 1D Gaussian kernel
filterSize = ceil(6*sigma);
x = -filterSize:filterSize;
g = exp(-(x.^2)/(2*sigma^2));
g = g / sum(g); % normalize

% Convolve ignoring NaNs
A_filtered = nan(size(MovM));
for r = 1:size(MovM,1)
    row = MovM(r,:);
    valid = ~isnan(row);

    if any(valid)
        % Convolution only on valid values
        num = conv(row .* valid, g, 'same');
        den = conv(valid, g, 'same');
        A_filtered(r,:) = num ./ den;

        % Where all neighbors were NaN, stay NaN
        A_filtered(r,den == 0) = NaN;
    end
end




% gF = fspecial('gaussian', 7, floor(g/2));
% MovM = imfilter(MovM, gF(4,:), 'same', 'replicate');
MovM=A_filtered(:,size(M,2)+1:2*size(M,2),:);
%MovM(isnan(M))=NaN;
end
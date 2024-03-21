function result = movprc(data, windowSize, percentile ,dim)
if nargin<4
    dim=1;
end

if dim==2
    n = size(data,2);
    result = zeros(size(data,1), n);
    
    for j=1:size(data,1)
    for i = 1:n
        % Define window boundaries
        startIdx = max(1, i - floor(windowSize/2));
        endIdx = min(n, i + floor(windowSize/2));
        
        % Extract the window of data
        windowData = data(j,startIdx:endIdx);
        
        % Calculate percentile and store it
        result(j,i) = prctile(windowData, percentile);
    end
    end
end

if dim==1;
    n = size(data,1);
    result = zeros(n,size(data,2));
    
    for j=1:size(data,2)
    for i = 1:n
        % Define window boundaries
        startIdx = max(1, i - floor(windowSize/2));
        endIdx = min(n, i + floor(windowSize/2));
        
        % Extract the window of data
        windowData = data(startIdx:endIdx,j);
        
        % Calculate percentile and store it
        result(i,j) = prctile(windowData, percentile);
    end
    end
end
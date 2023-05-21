function locs = peaksAround(data, locsBefore, neighborhood)

locs = zeros(1,length(locsBefore));
for i = 1:length(locsBefore)
    locArea = (locsBefore(i)-neighborhood):(locsBefore(i)+neighborhood);
    [~,locMax] = max(data(locArea));
    locs(i) = locArea(locMax);
end

function ind=bwlabel_firstind(trace)
[~, ind]=unique(bwlabel(trace));
ind=ind(2:end);
end

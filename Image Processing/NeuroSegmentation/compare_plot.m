function compareplot(m)
plot(1.3*m/(max(m(:))-min(m(:))));
hold on
plot(bsxfun(@minus,1.3*m/(max(m(:))-min(m(:))),1:size(m,2)));
hold off


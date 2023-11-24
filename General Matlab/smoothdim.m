function out = smoothdim(A,smoothing,dim)
shiftvect = zeros(1,length(size(A)));
shiftvect(dim) = 1;
B = A;
for i = 1:(smoothing-1)
    B = B + circshift(A,i*shiftvect);
end
B = B/smoothing;

out = B;
end

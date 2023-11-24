function out = fillmat(A,Bsize)
% does repmat on A to give it dimensions Bsize (only works if all
% dimensions of A are 1 or equal to correpsonding B dimension)
if length(Bsize) > length(size(A)) x = ones(1,length(Bsize) - length(size(A)) );
else x = [];
end
out = repmat(A,Bsize-[size(A), x]+1);
end

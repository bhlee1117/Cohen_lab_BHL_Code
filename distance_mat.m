function dist=distance_mat(X,Y)
%dist = distance_mat(X, Y) returns a matrix dist where dist(i,j) is the
%   Euclidean distance between the i-th row of X and the j-th row of Y.
%
%   INPUTS:
%       X - An N×D matrix, where each row is a D-dimensional point.
%       Y - An M×D matrix, where each row is a D-dimensional point.
%
%   OUTPUT:
%       dist - An N×M matrix of Euclidean distances.
%
%   Example:
%       X = [0 0; 1 1];
%       Y = [1 0; 2 2; -1 1];
%       d = distance_mat(X, Y);
%       % d(i,j) is the distance between X(i,:) and Y(j,:)
for i=1:size(X,1)
    for j=1:size(Y,1)
        dist(i,j)=sqrt(sum((X(i,:)-Y(j,:)).^2));
    end
end
end
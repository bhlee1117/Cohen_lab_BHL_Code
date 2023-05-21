function [moviefilt, eigvecs, eigvals] = pcafilt(movie, npcs);
% function [moviefilt, eigvecs, eigvals] = pcafilt(movie, npcs);
% uses time-domain PCA to de-noise a movie of a spike.
% npcs is the number of principal components to keep.  npcs must be less
% than or equal to the number of frames in the movie.  If npcs = nframes,
% then the movie is returned unaltered.
% returns the filtered movie, and all eigenvectors and eigenvalues

[ysize, xsize, nframes] = size(movie);

if npcs > nframes;
    'Error: npcs must be <= number of frames'
    return
end;
avgimg = mean(movie, 3);
mov2d = reshape(movie, [ysize*xsize, nframes]);
avgimg2d = reshape(avgimg, [ysize*xsize, 1]);
dmov2d = mov2d - repmat(avgimg2d, [1 nframes]);
covmat = dmov2d'*dmov2d;
[V, D] = eig(covmat);
eigvals = flipud(diag(D));
V = fliplr(V);
% correct the sign of all the components
for j = 1:nframes
    V(:,j) = V(:,j)*sign(V(find(abs(V(:,j)) == max(abs(V(:,j)))),j));
end;
eigvecs = V;

coeffs = mov2d*V(:,1:npcs);
mov2dfilt = coeffs*V(:,1:npcs)';
moviefilt = reshape(mov2dfilt, [ysize, xsize, nframes]);
moviefilt = moviefilt + repmat(avgimg, [1 1 nframes]);

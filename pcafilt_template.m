function [mov_resfilt] = pcafilt_template(mov_res, templates)
% [mov_resfilt, eigvecs, eigvals] = pcafilt_template(mov_res, templates)
% uses time-domain PCA to de-noise a mov_res of a spike.
% npcs is the number of principal components to keep.  npcs must be less
% than or equal to the number of frames in the mov_res.  If npcs = nframes,
% then the mov_res is returned unaltered.
% returns the filtered mov_res, and all eigenvectors and eigenvalues

[ysize, xsize, nframes] = size(mov_res);

avgimg = mean(mov_res, 3);
mov2d = reshape(mov_res, [ysize*xsize, nframes]);
avgimg2d = reshape(avgimg, [ysize*xsize, 1]);
dmov2d = mov2d - repmat(avgimg2d, [1 nframes]);
% covmat = dmov2d'*dmov2d;
% [V, D] = eig(covmat);
% eigvals = flipud(diag(D));
% V = V(:,end:-1:1);
% vSign = sign(max(V) - max(-V));  % make the largest value always positive
% V = V.*vSign;
% eigvecs = V;
EigV=tovec(templates);

[Q, ~] = qr(EigV, 0); % Orthonormal basis in columns of Q
coeffs = Q'*mov2d;

mov2dfilt = (Q*coeffs);
mov_resfilt = reshape(mov2dfilt, [ysize, xsize, nframes]);
mov_resfilt = mov_resfilt + repmat(avgimg, [1 1 nframes]);

for n=1:size(coeffs,1)
[fraction(n)] = get_variance(dmov2d', Q(:,n)');
end
disp(['Variance explained by template is ' num2str(sum(fraction))]);
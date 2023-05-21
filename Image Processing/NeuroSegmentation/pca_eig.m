function [u, s, v] = pca_eig(a,n)
%PCA_EIG Computes N Principal Components from matrix A.
%   [u, s, v] = pca_eig(a,n) computes a truncated, rank-N approximation to
%   the singular value decomposition of input matrix A. Matrices U, S, and
%   V are returned such that U*S*V' approximates A, U and V have
%   orthonormal columns, and S is a diagonal matrix with singular values
%   associated with each component.
%
%   2014 Vicente Parot
%   Cohen Lab - Harvard University
%

%% check input and output format
    thistic = tic;
    disp 'computing pca ...'
    swap = size(a,1) < size(a,2);
    if swap
        a = a';
    end
    [v, eigvals] = eigs(a'*a,n);
    s = sqrt(eigvals);
    u = a/(s*v');
    if swap
        aux = u;
        u = v;
        v = aux;
    end
    disp(['pca took ' num2str(toc(thistic)) ' s']);
end

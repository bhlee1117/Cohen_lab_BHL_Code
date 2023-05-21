function [u, s, v, a] = empca(a, ncomps, itertol, maxiters)
%EMPCA	Expectation-Maximization Principal Component Analysis
%   [U, S, V] = EMPCA(A,N) calculates N principal components of matrix A,
%   and returns U, S, V that approximate the N-rank truncation of the
%   singular value decomposition of A. S is a diagonal matrix with singular
%   values corresponding to the contribution of each principal component to
%   matrix A. U and V have orthonormal columns. 
% 
%   [U, S, V, E] = EMPCA(A,N) also returns the residual error matrix E
%   resulting from the PCA decomposition, such that A == U*S*V' + E.
% 
%   [...] = EMPCA(A,N,TOL) keeps principal components when the angular
%   change of U in the last EM iteration is smaller than TOL, instead of
%   the default value of 1e-6. 
% 
%   [...] = EMPCA(A,N,TOL,MAXITER) keeps principal components after MAXITER
%   EM iterations even if it does not convergence with angular change below
%   TOL. If omitted, a maximum of 100 EM iterations are computed.
% 
%   Matrix A must be a 2D matrix, and it can be of gpuArray class. N and
%   MAXITER must be a positive integers, TOL must be a real positive
%   scalar.
% 
%   This function implements the expectation maximization principal
%   component analysis algorithm by Stephen Bailey, available in 
%   http://arxiv.org/pdf/1208.4122v2.pdf
% 
%   Bailey, Stephen. "Principal Component Analysis with Noisy and/or
%   Missing Data." Publications of the Astronomical Society of the Pacific
%   124.919 (2012): 1015-1023.  
% 
%   2014 Vicente Parot
%   Cohen Lab
%   Harvard University
% 

%% notes
% keep notation from a = u*s*v'.
% regarding empca paper notation:
%   x   : a
%   phi : u
%   c   : sv
%   
    
%% parameters
if ~exist('itertol','var')
    itertol = 1e-6; % an eigenvector is found when 1-cos(angle) between EM iterations is below itertol
end
if ~exist('maxiters','var')
    maxiters = 100; % or when maxiters is reached, whatever first
end

gpu = isa(a,'gpuArray');
if gpu
    u = gpuArray.zeros(size(a,1),ncomps); % allocate memory for results
    sv = gpuArray.zeros(size(a,2),ncomps);
else
    u = zeros(size(a,1),ncomps); % allocate memory for results
    sv = zeros(size(a,2),ncomps);
end

normc = @(m)bsxfun(@rdivide,m,sqrt(sum(m.^2))); % returns normalized columns

%% empca
for comp = 1:ncomps % for each component
    if gpu
        u(:,comp) = normc(gpuArray.randn([size(a,1) 1]));
    else
        u(:,comp) = normc(randn([size(a,1) 1]));
    end
    for iter = 1:maxiters % repeat until u does not change or too many iterations
        u0 = u(:,comp); % store last iteration's u for comparison
        sv(:,comp) = a'*u(:,comp);       % E-step
        u(:,comp) = normc(a*sv(:,comp)); % M-step
        if 1-(u0'*u(:,comp)) < itertol % check convergence
            break % iter
        end
    end
    disp(['eigenvector ' num2str(comp) ' kept after ' num2str(iter) ' iterations'])
    a = a - u(:,comp)*sv(:,comp)'; % update a removing converged component and leaving residual
end
s = diag(sqrt(sum(sv.^2)));
v = normc(sv);

function out = hill(beta, L);
% function out = hill(L, beta);
% Returns the hill function with ligand concentrations L in a vector.
% beta is a 2-element vector containing the k_d and the cooperativity
% exponent, n
% AEC 1 April, 2011

kd = beta(1);
n = beta(2);
tmp = (L.^n)./(kd^n + L.^n);
out = tmp;


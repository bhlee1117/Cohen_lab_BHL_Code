function [q q_norm] = modularity(g, s)
% modularity    - modularity of the graph for given partitioning vector
%
%   q = modularity (g, s) modularity of the graph g for given partitioning
%   vector s. 
%
%         1
%   q = ----- s' * B * s
%        4*m
%
% B: modularity matrix.
% s: partitioning vector, s_i = 1 if vertex i belongs to group 1 and s_i =
% -1 if vertex i belong to group 2.
% m: the total number of edges.
%
% Example:
%  n = size(g,1);
%  s = (rand(1,n) > rand)*2-1; % generate random partitioning
%  q = modularity(g, s);
%
% Ref: Newman M. E. Proc Natl Acad Sci U S A 2006 
%
% See also MODULARITY2, MODMAT.

error(nargchk(2,2,nargin));

adj = adjacency(g);
deg = sum(adj);
m = sum(deg)/2;
s = s(:);
q = s' * modmat(g) * s / (4*m);

%add normalization
        n = length(adj);
        gnorm = zeros(size(adj));
        for k1 = 1:n
            for k2 = 1:n
                gnorm(k1, k2) = deg(k1) * deg(k2)*(s(k1)==s(k2)) / (2 * m)^2;
            end
        end
        norm_const=sum(sum(gnorm));
q_norm=q/(1-norm_const);
end

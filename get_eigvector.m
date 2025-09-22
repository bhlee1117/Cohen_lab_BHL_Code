function [V D eigTrace validInd]=get_eigvector(subMov,n)
% [V D eigTrace validInd]=get_eigvector(subMov,n)
% 2024.09.01, Byung Hun Lee

sz=size(subMov);
% if sz(1)>sz(2)
%     subMov=subMov';
% end
validInd=sum(isnan(subMov),1)==0;
subMov=subMov(:,validInd);
subMov=subMov-mean(subMov,2,'omitnan');
covMat=subMov*subMov';

[V, D] = eig(covMat);
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V = V.*vSign;

if nargin<2
n=size(V,1);
end
V=V(:,1:n);
eigTrace=subMov'*V;
end
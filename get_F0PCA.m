function [F0_PCA N]=get_F0PCA(VoltageTrace,N)
%Voltage trace = n X T, matrix,
%N= number of components you want to keep

subMov=VoltageTrace(:,sum(isnan(VoltageTrace))==0);
subMov=subMov-mean(subMov,2,'omitnan');
nTime=size(subMov,2);
%covMat=subMov(:,2:end)*subMov(:,1:end-1)';
covMat=subMov*subMov';
[V, D] = eig(covMat);
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V = V.*vSign;

if nargin<2
    N=find(cumsum(D)/sum(D)>0.95,1);
    N=[1:N];
    fprintf('%d component is used\n',N);
end
F0_PCA=sqrt(sum(((V(:,[N]).^2).*D([N])'),2)/nTime);
end
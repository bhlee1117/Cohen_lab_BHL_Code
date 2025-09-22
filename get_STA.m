function [STA newMatrix truePoints] = get_STA(nXT, binaryMatrix, tau1, tau2)
[n, T] = size(nXT);
% Ensure the binary matrix is a row vector

if size(binaryMatrix,1)>size(binaryMatrix,2)
    binaryMatrix=binaryMatrix';
end
if length(binaryMatrix)==T
    truePoints = find(binaryMatrix == 1);
else
    truePoints = binaryMatrix;
    disp('Treat binaryMatrix as index vector');
end

omitPoints=sum((truePoints'+[-tau1:tau2])<=0 | (truePoints'+[-tau1:tau2])>T,2)>0;
truePoints(omitPoints)=[];

newMatrix=reshape(nXT(:,truePoints'+[-tau1:tau2]),n,length(truePoints),[tau1+tau2+1]);
STA=squeeze(mean(newMatrix,2,'omitnan'));
end
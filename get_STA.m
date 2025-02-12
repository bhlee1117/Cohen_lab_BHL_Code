function [STA newMatrix truePoints] = get_STA(nXT, binaryMatrix, tau1, tau2)
    % Ensure the binary matrix is a row vector
    if size(binaryMatrix, 1) ~= 1
        error('The binary matrix must be a 1xT row vector.');
    end


    % Get the size of the input matrix
    [n, T] = size(nXT);

    truePoints = find(binaryMatrix == 1);
    omitPoints=sum((truePoints'+[-tau1:tau2])<=0 | (truePoints'+[-tau1:tau2])>T,2)>0;
    truePoints(omitPoints)=[];

    newMatrix=reshape(nXT(:,truePoints'+[-tau1:tau2]),n,length(truePoints),[tau1+tau2+1]);
    STA=squeeze(mean(newMatrix,2,'omitnan'));
end
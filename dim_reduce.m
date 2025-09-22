function red_x=dim_reduce(x)
if size(x,1)<size(x,2)
x=x';
end
x=x-mean(x,1);
covMat = x'*x;  % PCA within each region
    [V, D] = eig(covMat);
    D = diag(D); 
    D = D(end:-1:1);
    V = V(:,end:-1:1);
    vSign = sign(max(V) - max(-V));  % make the largest value always positive
    V = V.*vSign;
    red_x=x*V(1,:)';
end



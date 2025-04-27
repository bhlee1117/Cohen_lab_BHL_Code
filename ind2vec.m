function vec=ind2vec(N,ind,C,diffC)
% Length of N, 
% change ind to value C, others = diffC
if nargin<4
    diffC=0;
end
ind(ind<=0)=[];
vec=zeros(1,N);
vec(ind)=C;
vec(setdiff([1:N],ind))=diffC;
vec=vec(1:N);

end
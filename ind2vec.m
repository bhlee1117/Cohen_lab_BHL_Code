function vec=ind2vec(N,ind,C,diffC)

if nargin<4
    diffC=0;
end
vec=zeros(1,N);
vec(ind)=C;
vec(setdiff([1:N],ind))=diffC;

end
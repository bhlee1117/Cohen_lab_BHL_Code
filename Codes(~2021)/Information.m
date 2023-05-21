clear W
for i=1:100
    for j=1:100
        N=300000*i;
        a=0.02*j;
   W(i,j)=N*log(N)-a*N*log(a*N)-(N-a*N)*log(N-a*N);
    end
end
%%
imagesc(W)

%%
[X,Y] = meshgrid(0.02:0.02:1,20000:20000:1000000);
mesh(X,Y,W)
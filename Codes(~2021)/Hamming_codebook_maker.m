%%
clear 
N=100;
a=0.3;
b=[0.3 0.4];
HD=round(2*N*a*(1-b));
N_TXN=round(N*a); N_overlap=round(N*a*b);

%%
S=[zeros(1,N_TXN)+1 zeros(1,N-N_TXN)];
S_perm{1}=S;
S_perm{1}(1,randperm(N_TXN,HD(1,1)))=0;
S_perm{1}(1,N_TXN+randperm(N-N_TXN,HD(1,1)))=1;
i=1;


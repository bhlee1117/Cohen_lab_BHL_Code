function A_new=set_edge(A,bound,d)
A_new=A;
A_new(1:bound,:,:)=d;
A_new(:,1:bound,:)=d;
A_new(end-bound+1:end,:,:)=d;
A_new(:,end-bound+1:end,:)=d;
end
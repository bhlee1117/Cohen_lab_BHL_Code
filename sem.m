function S=sem(M,dim)

S=std(M,0,dim,'omitnan');
N=sum(~isnan(M),dim);
S=S./sqrt(N);
end
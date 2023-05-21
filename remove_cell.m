function Result=remove_cell(Result,rmv)

noi=[1:size(Result.centers,1)]';
noi(rmv)=[];
Result.centers(rmv,:)=[];
Result.c_ftprnt(:,:,rmv)=[];
Result.coord(rmv,:)=[];
try
for j=1:length(Result.clist_p)
clist_rmv=[];    
for k=1:length(Result.clist_p{j})
    if sum(ismember(rmv,Result.clist_p{j}(k)))>0
clist_rmv=[clist_rmv k];
    else
Result.clist_p{j}(k)=find(noi==Result.clist_p{j}(k));
end
end
Result.clist_p{j}(clist_rmv)=[];
end
end

Result.traces(rmv,:)=[];
Result.traces_bin(rmv,:)=[];
Result.traces_hi(rmv,:)=[];
Result.traces_bin_hi(rmv,:)=[];
Result.noise(rmv,:)=[];
Result.noise_res(rmv,:)=[];
Result.spike(rmv,:)=[];
Result.traces_res(rmv,:)=[];
Result.traces_res_hi(rmv,:)=[];

end
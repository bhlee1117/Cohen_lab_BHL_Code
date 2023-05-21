function rearr=rearrange_actmap(M,ref,array)
% %%
% M=cred_data;
% ref=5;
% array=[2 1 3];
for i=1:size(array,2)
numb(1,i)=sum(M(:,ref)==array(1,i));
end
%%
g=1;
gg=1;
ggg=1;
for i=1:size(M,1)
    if M(i,ref)==array(1,1)
        rearr(g,:)=M(i,:);
        g=g+1;
    end
    
    if M(i,ref)==array(1,2)
    rearr(numb(1,1)+gg,:)=M(i,:);
    gg=gg+1;
    end
    
    if M(i,ref)==array(1,3)
    rearr(numb(1,1)+numb(1,2)+ggg,:)=M(i,:);
    ggg=ggg+1;
    end
    
end
function refined_data=manual_masking(Y)
cmap=distinguishable_colors(40);
j=1;
imagesc(max(Y,[],3))
axis equal tight off
while menu('more?','Yes','No')==1
    ROI{j}=roipoly;
    refined_data.A_or(j,:)=reshape(ROI{j},size(Y,1)*size(Y,2),1);
    for t=1:size(Y,3)
        Yr(:,t)=reshape(Y(:,:,t),size(Y,1)*size(Y,2),1);
        refined_data.C(j,t)=sum(Yr(find(refined_data.A_or(j,:)==1),t));
    end
    
    refined_data.C_df(j,:)=(refined_data.C(j,:)-mean(refined_data.C(j,:)))/mean(refined_data.C(j,:));
    j=j+1;
end

RGB=zeros(512,512,3);
for i=1:size(refined_data.C_df,1)
    for k=1:3
            RGB(:,:,k)=RGB(:,:,k)+reshape(refined_data.A_or(i,:),512,512)*cmap(i,k)*20;
    end
end
imagesc(RGB)
axis equal tight off
end


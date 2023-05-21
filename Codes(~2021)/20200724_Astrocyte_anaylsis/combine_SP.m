function refined_data=combine_SP(A_or,Coor,Y,C_df)
%% Display
Cent=cell2mat(cellfun(@mean_BH,Coor,'UniformOutput',false));
Cent=reshape(Cent,2,size(Cent,1)/2);
C_df(find(cellfun(@isempty,Coor)),:)=[];
A_or(:,find(cellfun(@isempty,Coor)))=[];

%% Remove strange cells 
sw=1;
while sw 
    imagesc(mean(Y,3))
axis equal tight off
hold all
plot(Cent(1,:),Cent(2,:),'r.','markersize',15)
    [x y button]=ginput(1);
    close all
    if isempty(button)
        sw=0;
    end
    if button
        clear r
        for i=1:size(Cent,2) % Cells
           r(i,1)=(Cent(1,i)-x)^2+(Cent(2,i)-y)^2;
        end
        [m argm]=min(r);
        Cent(:,argm)=[];
        C_df(argm,:)=[];
        A_or(:,argm)=[];
    end
end
%% list Cells
cmap=distinguishable_colors(50);
    imagesc(mean(Y,3))
axis equal tight off
hold all
plot(Cent(1,:),Cent(2,:),'r.','markersize',15)
clear ROI list

j=1;
while menu('more?','Yes','No')==1
    ROI{j}=roipoly;
    list{j}=[];
    for i=1:size(Cent,2)
        if ROI{j}(round(Cent(2,i)),round(Cent(1,i)))
         list{j}(size(list{j},1)+1,1)=i;
        end
    end
    plot(Cent(1,list{j}),Cent(2,list{j}),'color',cmap(j,:),'marker','o','linestyle','none','linewidth',2)
    j=j+1;
end

%% Combine data
for j=1:size(list,2)
    refined_data.C_df(j,:)=mean(full(C_df(list{j},:)),1);
    refined_data.A_or(:,j)=sum(full(A_or(:,list{j})),2);
end
%%
RGB=zeros(512,512,3);
for i=1:size(list,2)
    for k=1:3
            RGB(:,:,k)=RGB(:,:,k)+reshape(refined_data.A_or(:,i),512,512)*cmap(i,k)*20;
    end
end
imagesc(RGB)
axis equal tight off
end


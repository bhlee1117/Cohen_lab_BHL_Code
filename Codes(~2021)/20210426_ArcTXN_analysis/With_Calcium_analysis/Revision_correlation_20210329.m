%% Calculate correlation of images across days.
%% Load images
clear
N=[12 16];
for g=1:2
    for m=1:N(g)
        [fnm{g,m} pth{g,m}]=uigetfile('*.tif','Multiselect','On');
    end
end

%%
for m=1:N(1)
    for i=1:size(fnm{1,m},2)
        for z=1:length(imfinfo([pth{1,m} fnm{1,m}{i}]))
            CA1{m,i}(:,:,z)=imread([pth{1,m} fnm{1,m}{i}],z);
        end
    end
end
for m=1:N(2)
    for i=1:size(fnm{2,m},2)
        for z=1:length(imfinfo([pth{2,m} fnm{2,m}{i}]))
            RSC{m,i}(:,:,z)=imread([pth{2,m} fnm{2,m}{i}],z);
        end
    end
end

%% Image crop or Gaussian Filtering(optional)
for m=1:N(1)
    for i=1:size(fnm{1,m},2)
            gCA1{m,i}=double(CA1{m,i});
            %gCA1{m,i}(gCA1{m,i}==0)=NaN;
            gCA1{m,i}=imgaussfilt3(gCA1{m,i},[2 2 2]);
    end
end
for m=1:N(2)
    for i=1:size(fnm{2,m},2)
           gRSC{m,i}=double(RSC{m,i});
            %gRSC{m,i}(gRSC{m,i}==0)=NaN;
            gRSC{m,i}=imgaussfilt3(gRSC{m,i},[2 2 2]);
    end
end

%% Calculate 3d correlations
for i=1:4
    for n=2:5
    %C{1}(i,n-1)=corr(gCA1{i,1}(:),gCA1{i,n}(:));
    C{1}(i,n-1)=corr2(max(gCA1{i,1},[],3),max(gCA1{i,n},[],3));
    end
end
for i=5:8
    for n=2:4
        %C{1}(i,n-1)=corr(gCA1{i,1}(:),gCA1{i,n}(:));
    C{1}(i,n-1)=corr2(max(gCA1{i,1},[],3),max(gCA1{i,n},[],3));
    end
end
for i=9:12
    for n=2:6
        %C{1}(i,n-1)=corr(gCA1{i,1}(:),gCA1{i,n}(:));
    C{1}(i,n-1)=corr2(max(gCA1{i,1},[],3),max(gCA1{i,n},[],3));
    end
end

for i=1:12
    for n=2:5
   % C{2}(i,n-1)=corr(gRSC{i,1}(:),gRSC{i,n}(:));
        C{2}(i,n-1)=corr2(max(gRSC{i,1},[],3),max(gRSC{i,n},[],3));
    end
end
for i=13:16
    for n=2:6
    %C{2}(i,n-1)=corr(gRSC{i,1}(:),gRSC{i,n}(:));
    C{2}(i,n-1)=corr2(max(gRSC{i,1},[],3),max(gRSC{i,n},[],3));
    end
end

%% Load the Data, CA1 & RSC and plot
cmap=[0.2 0.8 0.8;0.1 0.5 0.5;0.8 0.2 0.8; 0.5 0.1 0.5];
   figure
for i=1:2
    C{i}(C{i}==0)=NaN;
 
    errorbar([1:5],mean(C{i}(1:end-4,:),1,'omitnan'),std(C{i}(1:end-4,:),0,1,'omitnan')./sqrt(sum(~isnan(C{i}(1:end-4,:)),1)),'color',cmap(2*i-1,:))
    hold all
end
ylim([0 1])
xlim([0.5 4.5])
set(gca,'Xtick',[1 2 3 4])
xlabel('Days')
ylabel('Correlation coeffcient')
figure
for i=1:2
    C{i}(C{i}==0)=NaN;
 
    errorbar([1 8 15 22 29],mean(C{i}(end-3:end,:),1,'omitnan'),std(C{i}(end-3:end,:),0,1,'omitnan')./sqrt(sum(~isnan(C{i}(end-3:end,:)),1)),'color',cmap(2*i-1,:))
    hold all
end
set(gca,'Xtick',[1 8 15 22 29])
ylim([0 1])
xlim([0.5 30.5])
xlabel('Days')
ylabel('Correlation coeffcient')

%% OP data
[~, unqInd] = unique([Mouse NeuronInd] ,'row');
f=unqInd(53);

load([fpath{f} '/OP_Result.mat'])
SkelDend = Skeletonize_dendrite(Result.ref_im,8,0.01,25);
%SkelDend=SkelDend.*double(max(Result.ftprnt>0,[],3));
SkelDend(mean(Result.ftprnt,3)==0)=0;
SkelDend(abs(mean(Result.ftprnt,3))>0)=1;
[~, mask]=get_ROI(SkelDend+mat2gray(Result.ref_im)*13);
SkelDend(max(mask,[],3)>0)=0;

intDD=[];
figure(f); clf;
nROI=size(Result.ftprnt,3);
pth=[];
g=1;
for i=1:nROI
    i
    for j=i:nROI
        [intDD(i,j), pth(:,:,g)]=geodesic_distance(SkelDend,get_coord(Result.ftprnt(:,:,i)),get_coord(Result.ftprnt(:,:,j)));
        pth(:,:,g)=pth(:,:,g)*i*j;
        %nexttile([1 1]);
        %imshow2(im_merge(cat(3,Result.ref_im,pth),[1 1 1; 1 0 0]),[])
        g=g+1;
    end
end
intDD=max(cat(3,intDD,intDD'),[],3);
Result.interDendDist=intDD;
save(fullfile(fpath{f},'OP_Result.mat'),'Result','-v7.3')
imagesc(intDD);

SameCellInd=find(Mouse==Mouse(f) & NeuronInd==NeuronInd(f));
for k=SameCellInd'
    load([fpath{k} '/OP_Result.mat'])
    Result.interDendDist=intDD;
    save(fullfile(fpath{k},'OP_Result.mat'),'Result','-v7.3')
    disp(['File #' num2str(k) ' is saved.'])
end

%% PC data

f=11;
load([fpath{f} '/PC_Result.mat'])
SkelDend = Skeletonize_dendrite(Result.ref_im,8,0.01,25);
%SkelDend=SkelDend.*double(max(Result.ftprnt>0,[],3));
SkelDend(mean(Result.ftprnt,3)==0)=0;
SkelDend(abs(mean(Result.ftprnt,3))>0)=1;
[~, mask]=get_ROI(SkelDend+mat2gray(Result.ref_im)*13);
SkelDend(max(mask,[],3)>0)=0;

intDD=[];
figure(199); clf;
nROI=size(Result.ftprnt,3);
for i=1:nROI
    i
    for j=i:nROI
        [intDD(i,j), pth]=geodesic_distance(SkelDend,get_coord(Result.ftprnt(:,:,i)),get_coord(Result.ftprnt(:,:,j)));
        nexttile([1 1]);
        imshow2(im_merge(cat(3,Result.ref_im,pth),[1 1 1; 1 0 0]),[])
    end
end
intDD=max(cat(3,intDD,intDD'),[],3);
Result.interDendDist=intDD;

imagesc(intDD)

save(fullfile(fpath{f},'PC_Result.mat'),'Result','-v7.3')
load([fpath{f} '/PC_Result.mat'])
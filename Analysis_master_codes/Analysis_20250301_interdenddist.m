f=55;
load([fpath{f} '/OP_Result.mat'])
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
            %nexttile([1 1]);
            %imshow2(im_merge(cat(3,Result.ref_im,pth),[1 1 1; 1 0 0]),[])
        end
    end
    intDD=max(cat(3,intDD,intDD'),[],3);
    Result.interDendDist=intDD;

    imagesc(intDD)

    save(fullfile(fpath{f},'OP_Result.mat'),'Result','-v7.3')
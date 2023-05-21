[FileName,PathName]=uigetfile('*.tif','Select the tif file');
PathFileName=[PathName,FileName];

ImageInfo=imfinfo(PathFileName);
ImageNum=size(ImageInfo,1);
SliceAndFrame=strsplit(ImageInfo(1).ImageDescription,{'slices=','frames=','hyperstack='});
SliceNum=str2num(SliceAndFrame{1,2});
FrameNum=str2num(SliceAndFrame{1,3});
PosNum=ImageNum/SliceNum/FrameNum;

for i=1:PosNum
    for j=1:FrameNum
        Stack=[];
        for k=1:SliceNum
            Stack=[Stack;imread(PathFileName,'tif',k+SliceNum*(i-1)+SliceNum*PosNum*(j-1))];
        end
        Stack=permute(reshape(Stack,512,SliceNum,512), [1 3 2]);
        MaxProjection=max(Stack,[],3);
        imwrite(MaxProjection,[PathFileName,'MaxProjectionPos',num2str(i),'.tif'],'WriteMode','append');
    end
end
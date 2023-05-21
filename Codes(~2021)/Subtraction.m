clear
bigkernel=5;
smallkernel=1;
currentPath=uigetdir('','Select the Root Folder');
currentDir=dir(currentPath);
lengthList=[length(currentDir)];
indexList=[3];
splf=split(currentDir(1).folder,'\');
outputdir=[];
for i=1:size(splf)-1
outputdir=strcat(outputdir,[char(splf(i,1)) '\']);
end

mkdir([char(outputdir)],'Subtraction\')
outputdir=[char(outputdir) 'Subtraction\'];
%%
while length(lengthList)
    
    %If the pointer went too far;
    if indexList(end)>lengthList(end)
            %Go back one step move the pointer
            indexList=indexList(1:end-1);
            lengthList=lengthList(1:end-1);
            %Break the loop if the work is done
            if length(lengthList)>0
                indexList(end)=indexList(end)+1;
                currentDir=dir([currentPath,'\..']);
                currentPath=currentDir(1).folder;
            end
            continue;
    end
    
    %If the pointer is on the directory, go one step deeper
    while currentDir(indexList(end)).isdir
        currentPath=[currentPath,'\',currentDir(indexList(end)).name];
        currentDir=dir(currentPath);
        lengthList=[lengthList;length(currentDir)];
        indexList=[indexList;3];
    end
    
    %Find the first tiff file starting with 'img' 
    while currentDir(indexList(end)).name(1:2)~='2x'
        indexList(end)=indexList(end)+1;
        if indexList(end)>lengthList(end)
            break;
        end
    end
    
    if length(lengthList)==0
        break;
    end
    
    %If there is none,
    %Go back one step and move the pointer
    if indexList(end)>lengthList(end)
            indexList=indexList(1:end-1);
            lengthList=lengthList(1:end-1);
            if length(lengthList)>0
                indexList(end)=indexList(end)+1;
                currentDir=dir([currentPath,'\..']);
                currentPath=currentDir(1).folder;
            end
            continue;
    end
    
    %numStack=0;
  %  stack=[];
    %If the right tiff file is found, do the Maximum Intensity Projection
for i=1:size(currentDir,1)-2
        stack(:,:,1)=imread([currentPath,'\',currentDir(indexList(end)).name],1);
        stack(:,:,2)=imread([currentPath,'\',currentDir(indexList(end)).name],2);
        stack=double(stack);
        subtraction=abs(stack(:,:,2)/max(reshape(stack(:,:,2),size(stack,1)*size(stack,2),1))-stack(:,:,1)/max(reshape(stack(:,:,1),size(stack,1)*size(stack,2),1)))*max(reshape(stack(:,:,2),size(stack,1)*size(stack,2),1));
%         Bgimfilt=imgaussfilt(subtraction, bigkernel);
%         subim=subtraction-Bgimfilt;
%         Filt_subtraction=imgaussfilt(subim,smallkernel);

        splf2=split(currentPath,'\');
        mkdir(outputdir,[char(splf2(size(splf2,1),1)) '\'])
        imwrite(uint16(subtraction*10),[outputdir,[char(splf2(size(splf2,1),1)) '\'] 'Sub_' currentDir(indexList(end)).name '.tif'],'tif');
        indexList(end)=indexList(end)+1;
       % numStack=numStack+1;
    end
%     
%     stack=permute(reshape(stack,512,numStack,512),[1 3 2]);
%     projection=max(stack,[],3);
    

    indexList(end)=indexList(end)+1;
    
    %Go back one step and move the pointer
    indexList=indexList(1:end-1);
    lengthList=lengthList(1:end-1);
    
    if length(indexList)>0
        indexList(end)=indexList(end)+1;
    end
    
    currentDir=dir([currentPath,'\..']);
    currentPath=currentDir(1).folder;
    
end

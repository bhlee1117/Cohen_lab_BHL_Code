currentPath=uigetdir('','Select the Root Folder');
currentDir=dir(currentPath);
lengthList=[length(currentDir)];
indexList=[3];
mkdir=currentDir(1).folder
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
    while currentDir(indexList(end)).name(1:4)~='img_'
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
    
    numStack=0;
  %  stack=[];
    %If the right tiff file is found, do the Maximum Intensity Projection
    %while [currentDir(indexList(end)).name(1:2)]=='2x'
        stack=imread([currentPath,'\',currentDir(indexList(end)).name],1);
        stack=[stack;imread([currentPath,'\',currentDir(indexList(end)).name],2)];
        imwrite(projection,[outputPath,'\MAX.tif'],'tif');
        indexList(end)=indexList(end)+1;
        numStack=numStack+1;
    %end
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

[fileName,pathName] = uigetfile('*.dcimg','Select a DCImg file');
if isequal(fileName,0)
   disp('User selected Cancel')
else
   disp(['User selected ', fullfile(pathName, fileName)])
   % set default dialog open directory to the present working directory
   lastDir = pathName;
   dcimgfile = fullfile(pathName, fileName);
   dcimgfile = strrep(dcimgfile, '\', '\\');
end

[framedata,totalframes]= dcimgmatlab(1, dcimgfile);
framedatatrans = transpose (framedata);
[ysize, xsize] = size(framedatatrans);

prompt = {'Enter Start Frame Index:','Enter End Frame Index:'};
dlg_title = 'Input';
num_lines = 1;
def = {'1',num2str(totalframes)};
answer = inputdlg(prompt,dlg_title,num_lines,def);

startframestr = (answer{1});
endframestr = (answer{2});

startframe = str2num(['int32(' startframestr ')']);
endframe = str2num(['int32(' endframestr ')']);

if (startframe == - 1)
    startframe = 0;
end



if (endframe == - 1)
    endframe = int32(totalframes(1,1));
end

numFrames = endframe - startframe + 1;

% Preallocate the array
seq1 = uint16(zeros(ysize,xsize,numFrames)); 
seq1(:,:,1) = framedatatrans;

for frame=startframe:endframe
    frame
   % Read each frame into the appropriate frame in memory.
   [framedata,totalframes]= dcimgmatlab(frame, dcimgfile);
   framedatatrans = transpose (framedata);
   seq1(:,:,frame)  = framedatatrans;
   pause
end

figure(1)
for j = startframe:endframe
    imshow(seq1(:,:,j), [])
    title(num2str(j))
    pause
end;

implay(seq1,numFrames);

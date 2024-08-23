clear
clc;
cd '/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'B5:K154');

save_to='/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,5));
CamType=raw(:,3);
StructureData=raw(:,10);
%%

[a, unqInd] = unique([Mouse NeuronInd] ,'row');
figure;
for i=unqInd([41])'
    load(fullfile(fpath{i},'Result.mat'))
    avgImg=Result.ref_im;
    nexttile([1 1])
    load(fullfile(fpath{i},"output_data.mat"))
    disp(fpath{i})

    switch char(CamType(i))
        case 'flash'
            sz=double(Device_Data{1, 4}.ROI([2 4]));
        case 'fusion'
            sz=double(Device_Data{1, 3}.ROI([2 4]));
    end

    ref_time=[2000:4000];
    mov_test=double(readBinMov_times([fpath{i} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1),ref_time));

    [u,s,v] = svds(tovec(mov_test-mean(mov_test,3)),20);
    reshape_u=reshape(u,sz(2),sz(1),[]);

    [~, bvMask]=get_ROI(max(abs(reshape_u),[],3),Result.bvMask);

    SameCellInd=find(Mouse==Mouse(i) & NeuronInd==NeuronInd(i));

    for j=SameCellInd'
        load(fullfile(fpath{j},'Result.mat'))
        [offsetY offsetX]=calculate_shift(avgImg,Result.ref_im)
        [x, y] = meshgrid(1:sz(1), 1:sz(2));
        Result.bvMask=[];
        for n=1:size(bvMask,3)
            Result.bvMask(:,:,n) = interp2(x, y, bvMask(:,:,n), x+offsetX, y+offsetY, 'linear', 0);
        end
        save([fpath{j} '/Result.mat'],'Result','-v7.3')
    end


end
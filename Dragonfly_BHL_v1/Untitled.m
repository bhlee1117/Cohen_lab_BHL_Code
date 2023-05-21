[fname,fpath] = uigetfile('*.*'); if fname == 0, return;end
mov = vm(fpath);


seg_size = 100;
    FOV_size = [size(mov,1) size(mov,2)];
    id11 = [1:seg_size:FOV_size(1)];
    id12 = id11+seg_size-1; id12(id12>FOV_size(1))=FOV_size(1);
    id12 = unique(id12);id11 = id11(1:length(id12));
    id21 = [1:seg_size:FOV_size(2)];
    id22 = id21+seg_size-1; id22(id22>FOV_size(2))=FOV_size(2);
    id22 = unique(id22);id21 = id21(1:length(id22));
    
%%
    mov_new = vm(zeros(size(mov,1),size(mov,2),40));

    for j=1:length(id11)
        for i=1:length(id21)
            if i*j==64;continue;end
            mov_new(id11(i):id12(i),id21(j):id22(j),:)=mat2gray(mov(id11(i):id12(i),id21(j):id22(j),...
            [(j-1)*40*length(id21)+(i-1)*40+1:(j-1)*40*length(id21)+i*40]+40).data);
%         fprintf('[%g,%g] [%g,%g]\n',[id11(i) id12(i),id21(j) id22(j)])
%         sprintf('%g,%g\n',[(j-1)*40*length(id21)+(i-1)*40+1:(j-1)*40*length(id21)+i*40]+40)
%         pause
        end
    end
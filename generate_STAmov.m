function STAmov=generate_STAmov(spike,mov_mc,range)


s=find(spike); g=1;
STAmov=zeros(size(mov_mc,1),size(mov_mc,2),sum(range)+1);

for i=s
    try
        STAmov=STAmov+mov_mc(:,:,i-range(1):i+range(2));
        g=g+1;
    end
end
STAmov=STAmov/g;
disp(['N = ' num2str(g) ' spikes']);
end

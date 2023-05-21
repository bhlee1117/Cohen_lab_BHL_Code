function [STA ST]=generate_STA(spike,volt,range)


s=find(spike); g=1;
for i=s
    try
        ST(:,:,g)=volt(:,i-range(1):i+range(2));
        g=g+1;
    end
end
STA=mean(ST,3);
end

function STA=make_STA(sourceMV,ts,wind)
STA=zeros(size(sourceMV,1),size(sourceMV,2),length(wind));
g=0;
for t=ts
    try
    STA=STA+sourceMV(:,:,t+wind);
    g=g+1;
    end
end
STA=STA./g;
end
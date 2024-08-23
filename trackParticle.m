function filteredParticles=trackParticle(particlesPos,D,minFrames)

pp=[];
for i=1:size(particlesPos,2)
pp=[pp; [particlesPos{i}(:,1:2) repmat(i,size(particlesPos{i},1),1)]];
end
parma=struct('mem',0,'good',0,'dim',2,'quiet',0);
linked=track(pp,D,parma);

filteredParticles=[]; g=1;
for p=unique(linked(:,4))' 
    ind=find(linked(:,4)==p);
    if length(ind)>minFrames
        filteredParticles{g}=linked(ind,1:3);
        g=g+1;
    end
end

%isvalid=cell2mat(cellfun(@(x) sum(x==0,[1 2])==0,filteredParticles,'UniformOutput',false));
%filteredParticles=filteredParticles(isvalid);
end
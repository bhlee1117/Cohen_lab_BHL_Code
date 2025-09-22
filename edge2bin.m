function binn=edge2bin(edge)
binn=mean([edge(2:end); edge(1:end-1)]);
end
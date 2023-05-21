function fpath_window=fpath_mac2window(fpath)
sp=split(fpath,'/');
i=1;
try
while isempty(strfind(sp{i},'Lab'))
i=i+1;
end
ref=i;
fpath_window=['X:\Lab'];
for j=ref+1:length(sp)
fpath_window=[fpath_window '\' sp{j}];
end
catch
fpath_window=[];    
for j=1:length(sp)
fpath_window=[fpath_window '\' sp{j}];
end    
end
end

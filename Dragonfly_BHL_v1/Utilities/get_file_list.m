function file_list = get_file_list(file_path,file_extension)

list = dir([file_path]);
file_list = {};
for ii=1:length(list)
    if contains(list(ii).name,file_extension)
        [~,name,~]=fileparts(list(ii).name);
        file_list = [file_list name];
    end
end

end
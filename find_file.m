function N=find_file(fpath,stringToBeFound)

P = fpath;
F = sprintf(['*' stringToBeFound '*'],9129);
S = dir(fullfile(fpath,F));
try 
    N = S.name;
catch
    disp('No matching files with that string')
    N=[];
end
end
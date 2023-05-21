function movie_to_bin(varargin)

if size(varargin,2) ~= 2
    error('movie_to_bin:inArgs', 'Usage: movie_to_bin(movie, path)')
end

movie = varargin{1};
path = varargin{2}; 

fid = fopen(path, 'w');
fwrite(fid, permute(movie,[2 1 3]), 'ubit1','b');     % DMD reads by rows, Matlab saves by columns, big endian ordering
fclose(fid);
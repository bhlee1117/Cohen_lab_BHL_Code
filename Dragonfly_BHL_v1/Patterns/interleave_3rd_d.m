function inter = interleave_3rd_d(varargin)
%interleave_3rd_d Interleaves input matrices along their third dimension.
%
%   2018-2019 Vicente Parot
%   Cohen Lab - Harvard university
%
%   interleaves series of patterns of equal length, for concurrent
%   acquisition.
    inter = reshape(permute(cell2mat(permute(varargin,[1 4 3 2])),[1 2 4 3]),size(varargin{1},1),size(varargin{1},2),[]);
end

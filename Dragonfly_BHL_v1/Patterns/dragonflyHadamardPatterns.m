function alp_patterns = dragonflyHadamardPatterns(dmd,m_q,pat) % hadamard_mode)
%
%   2018-2019 Vicente Parot
%   Cohen Lab - Harvard university
%
%             roirows = 19*10+0:57*10+3; % extended central quad 2017-03-17
%             roicols = 32*10+0:69*10+5;
%             hadamard_mode = 'voltage';
%             device = [];
if exist('dmd','var') && isvalid(dmd)
    ncols = dmd.device.height;
    nrows = dmd.device.width;
else
    ncols = 1024;
    nrows = 768;
end
% sequence_mode = 'hadamard_23_5'; % structural
sequence_mode = 'hadamard_59_8'; % structural
super_mask = ones(ncols,nrows,'uint8');
switch sequence_mode
    case 'hadamard_59_8'
        clc
        blocksize = [63 14];
        elementsize = 1;
        alp_patterns = hadamard_patterns_scramble_nopermutation(blocksize,elementsize,super_mask);

        alp_patterns = cat(3,alp_patterns(:,:,1:10),alp_patterns);
        alp_patterns = reshape([alp_patterns alp_patterns*0],size(alp_patterns).*[1 1 2]);
    case 'hadamard_23_5'
        clc
        blocksize = [23 5];
        elementsize = 1;
        alp_patterns = hadamard_patterns_scramble_nopermutation(blocksize,elementsize,super_mask);

        alp_patterns = alp_btd_to_logical(alp_patterns);
        npats = ~alp_patterns&any(alp_patterns,3);
        alp_patterns = alp_logical_to_btd(alp_patterns);
        npats = alp_logical_to_btd(npats);
        alp_patterns = interleave_3rd_d(alp_patterns,npats);
        clear npats

        alp_patterns = cat(3,alp_patterns(:,:,1:10),alp_patterns);
        alp_patterns = reshape([alp_patterns alp_patterns*0],size(alp_patterns).*[1 1 2]);
end

end

function inter = interleave_3rd_d(varargin)
%   interleaves series of patterns of equal length, for concurrent
%   acquisition.
    inter = reshape(permute(cell2mat(permute(varargin,[1 4 3 2])),[1 2 4 3]),size(varargin{1},1),size(varargin{1},2),[]);
end

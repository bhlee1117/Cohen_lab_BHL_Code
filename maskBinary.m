function masked = maskBinary(matrix, BinaryMat, Value2put)
% masked = maskBinary(matrix, A, C)
% Input: matrix = 2D or 3D array (rows x cols x T) to mask edges,
%        width = edge boundary width to zero out (default: 5)
% Output: masked = matrix with edge boundaries set to 0
sz=size(matrix);
sz_A=size(BinaryMat);
if sz_A(1)==sz(1) & sz_A(2)==sz(2)
else
    error('Input and binary matrix are not matched')
end
Avec=reshape(BinaryMat,sz(1)*sz(2),[]);
matrixVec=reshape(matrix,sz(1)*sz(2),[]);
matrixVec(Avec>0,:)=Value2put;
masked=reshape(matrixVec,sz(1),sz(2),[]);
end
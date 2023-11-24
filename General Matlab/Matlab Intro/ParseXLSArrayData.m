function [ S ] = ParseXLSArrayData( num,txt )
%ParseXLSArrayData converts imported array data into a useful structure
%   Creates a structure S
%This structure has three field
%The first field is the element in the header of the first text column.
%The second field is the rest of the header and is called header.
%The third field is called data and is an array of the microarray data.

S.header=txt(1,[2,3]);
S.(txt{1,1})=txt(2:end,1);
S.data=num;

end


function zoneAreas=generate_venn(binMat)
% zoneAreas=generate_venn(binMat)
%   Given a 3 x N binary matrix, compute overlaps between the rows
%   and draw a 3-circle Venn diagram using the custom venn() function.
%   Each row represents a class, each column a data item.
% Zone definitions:
% z1   = A only
% z2   = B only
% z3   = C only
% z12  = A & B only
% z13  = A & C only
% z23  = B & C only
% z123 = A & B & C

assert(size(binMat,1) == 3, 'Input must be a 3xN binary matrix');

A = binMat(1,:);
B = binMat(2,:);
C = binMat(3,:);

onlyA   = A & ~B & ~C;
onlyB   = B & ~A & ~C;
onlyC   = C & ~A & ~B;
AB      = A & B & ~C;
AC      = A & C & ~B;
BC      = B & C & ~A;
ABC     = A & B & C;

zoneAreas = [sum(onlyA), sum(onlyB), sum(onlyC), ...
             sum(AB), sum(AC), sum(BC), sum(ABC)];

% Call venn()
figure;
venn(zoneAreas, 'FaceAlpha', {0.6, 0.6, 0.6}, 'EdgeColor', 'black');
title('3-Class Venn Diagram');
end

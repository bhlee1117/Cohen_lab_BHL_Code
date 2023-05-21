function out = mats2RGB(r, g, b);
% function out = mats2RGB(r, b, b);
% Convert three arrays (must be the same size) to an RGB image.
% Each array is scaled so min -> 0, max -> 1.
% AEC 2/13/08

rmin = min(min(r));
rmax = max(max(r));
gmin = min(min(g));
gmax = max(max(g));
bmin = min(min(b));
bmax = max(max(b));

if (rmax - rmin);
    rs = (r - rmin)/(rmax - rmin);
else
    rs = r;
end;
if (gmax - gmin);
    gs = (g - gmin)/(gmax - gmin);
else;
    gs = g;
end;
if (bmax - bmin);
    bs = (b - bmin)/(bmax - bmin);
else;
    bs = b;
end;

sz = size(rs);
out = zeros(sz(1), sz(2), 3);
out(:,:,1) = rs;
out(:,:,2) = gs;
out(:,:,3) = bs;

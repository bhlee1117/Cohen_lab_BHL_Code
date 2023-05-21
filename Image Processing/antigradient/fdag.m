function ghsb = fdag(cat3yx,off)
%FAG    Fourier Domain Anti-Gradient
%   F = FDAG(G) computes the function f which has a gradient
%   as close as possible to g. The vector field g must be an
%   MxNx2 array.
%
%   F = FDAG(G, MU) does the same thing but sets the unrecoverable
%   mean of f to mu. When omitted, the mean of f is set to zero.
%
%   Vicente Parot
%   M+Vision - Colonoscopy
%   Massachusetts Institute of Technology
%   September 2012
%
    if nargin < 2
        off = 0;
    end
    rows = size(cat3yx,1);
    cols = size(cat3yx,2);
    fxsb = ifftshift(fft2(fftshift(cat3yx(:,:,1))));
    fysb = ifftshift(fft2(fftshift(cat3yx(:,:,2))));
    [n m] = meshgrid(1:cols,1:rows);
    m = m - (rows/2+1);
    n = n - (cols/2+1);
    dx = 1/rows;
    dy = 1/cols;
    oscx = sin(2*pi*m*dx);
    oscy = sin(2*pi*n*dy);
    h = oscx + 1i*oscy;
    g = fxsb + 1i*fysb;
    div = g./h;
    div(rows/2+1,cols/2+1) = 1i*off*numel(n);
    div(       1,cols/2+1) = 0;
    div(rows/2+1,       1) = 0;
    div(       1,       1) = 0;
    f = ifftshift(ifft2(fftshift(div)));
    ghsb = imag(f);
end
function [t, uDat, vDat] = aliev2D(t0, tfinal, u0, v0, a, k, c, nx, ny, phi);
%     [t, uDat, vDat] = aliev(t0, tfinal, u0, v0, a, k, c, nx, ny, phi);
%
%     Simulate the Aliev-Panfilov model of the cardiac AP
%     du/dt = F(u, v) = ku(1 - u)(u - a) - uv + c + phi*del^2(u),  [The last term accounts for spatial coupling]
%     dv/dt = H(u, v) = ?(u)(ku - v), 
%     ?(u < 0.05) = 1 and ?(u > 0.05) = 0.1
%     yp = [du/dt   dv/dt];
%
%     simulates on a grid of dimensions [ny, nx]
%     Scalars: phi
%     Dimension [ny, nx]: u0, v0, a, k, c (driving current)
%
%   AEC 9 July 2015

y0 = [reshape(u0, [nx*ny, 1]); reshape(v0, [nx*ny, 1])];
filt = [0 1 0;
        1 -4 1;
        0 1 0];
% Simulate the differential equation.
% tfinal = tfinal*(1+eps);
[t,y] = ode45(@dAliev2D,[t0 tfinal],y0);

nT = length(t);
uDat = reshape(y(:,1:nx*ny), [nT, ny, nx]);
uDat = shiftdim(uDat, 1);
vDat = reshape(y(:,(nx*ny+1):end), [nT, ny, nx]);
vDat = shiftdim(vDat, 1);

function yp = dAliev2D(t, y);

u = y(1:nx*ny);
v = y((nx*ny + 1):end);

% calculate the currents from nearest neighbors
uMat = reshape(u, [ny nx]);
uCurr = imfilter(uMat, filt, 'replicate');  % excess current into each cell.



up = k(:).*u.*(1 - u).*(u - a(:)) - u.*v + c(:) + phi.*uCurr(:);

epsilon = ones(size(u));
epsilon(u >= 0.05) = 0.1;

vp = epsilon.*(k(:).*u - v);

yp = [up; vp];
end

end

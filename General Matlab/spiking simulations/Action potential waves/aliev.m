function [t, y] = aliev(t0, tfinal, u0, v0, a, k, c);
%     [t, y] = aliev(t0, tfinal, u0, v0, a, k, c);
%
%     Simulate the Aliev-Panfilov model of the cardiac AP
%     du/dt = F(u, v) = ku(1 - u)(u - a) - uv + c,  [Adam added the +c term
%     to simulate a drive current.
%     dv/dt = H(u, v) = ?(u)(ku - v), 
%     ?(u < 0.05) = 1 and ?(u > 0.05) = 0.1
%     yp = [du/dt   dv/dt];
%
%   AEC 9 July 2015

y0 = [u0; v0];
% Simulate the differential equation.
% tfinal = tfinal*(1+eps);  % This was in the example code.  Not sure why
% we need it.
[t,y] = ode45(@dAliev,[t0 tfinal],y0);

function yp = dAliev(t, y);

u = y(1);
v = y(2);

yp(1) = k.*u.*(1 - u).*(u - a) - u.*v + c;
if u < 0.05;
    epsilon = 1;
else
    epsilon = 0.1;
end;
yp(2) = epsilon.*(k.*u - v);

yp = yp';
end

end

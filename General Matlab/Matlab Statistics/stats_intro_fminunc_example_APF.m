% Simple least-squares regression as a demo of fminunc
% Generate data
X = (1:1:10)';
B0 = 3;
B1 = -2;
sigmayx = 2;
Y = B0 + B1*X + sigmayx*randn(size(X));
figure(1)
clf
plot(X,Y,'x')
legend('data')

% Do the regression the normal way
X1 = [ones(size(X)),X];
b = X1\Y;
hold all
plot(X,B0+B1*X,'--')
plot(X,b(1)+b(2)*X)
legend({'data','theoretical trend','regression'})

% Create a function to sum squared residuals for a given slope and offset
sumsqr = @(x,y,b) sum((y-(b(1)+x*b(2))).^2);
% Configure options for minimization
fitopts = optimset('Display','iter'); % Output status after each iteration
% Use fminunc to minimize the output of that function, starting from 0 offset and 0 slope
fitres = fminunc(@(b) sumsqr(X,Y,b),[0;0],fitopts);
% Observe fitres == b
% In this case fminunc is massive overkill, but some models are too complex for built-in functions

% Use fmincon() in place of fminunc() to put constraints on the search
% Use optimset() to select various options for either function, such as selecting the minimization algorithm
%   or specifying to use parallel processing (especially for multi-dimensional searches; read help on "matlabpool" first)

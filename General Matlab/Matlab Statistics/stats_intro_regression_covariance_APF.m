% EXAMPLE 1: simple linear regression
% Generate and plot some data
X = (1:1:10)';
B0 = 3;
B1 = -2;
sigmayx = 2;
Y = B0 + B1*X + sigmayx*randn(size(X));
figure(1)
clf
plot(X,Y,'x')
legend('data')

% Least-squares regression, simplest version:
X1 = [ones(size(X)),X];
b = X1\Y;
% If you want some statistics about the results:
[b,bint,r,rint,stats] = regress(Y,X1);
% For user-set confidence intervals:
alpha = 0.05;
[b,bint,r,rint,stats] = regress(Y,X1,alpha);
% b is same as before
% bint is 95% confidence interval for b (or 1-alpha level)
% r is residuals
% rint is 95% confidence interval for r (or 1-alpha level)
% stats gives R^2, F statistic, F-stat p-value, estimate of error variance
hold all
plot(X,B0+B1*X,'--')
plot(X,b(1)+b(2)*X)
legend({'data','theoretical trend','regression'})

% Outlier rejection:
[b_robust,stats_robust] = robustfit(X,Y);
plot(X,b_robust(1)+b_robust(2)*X);
legend({'data','theoretical trend','regression','robust regression'})
% Not that different in this case, but...if we add an outlier:
Y(end) = 10;
figure(2)
clf
plot(X,Y,'x')

[b,bint,r,rint,stats] = regress(Y,X1);
[b_robust,stats_robust] = robustfit(X,Y);
hold all
plot(X,B0+B1*X,'--')
plot(X,b(1)+b(2)*X)
plot(X,b_robust(1)+b_robust(2)*X);
legend({'data','theoretical trend','regression','robust regression'})

% EXAMPLE 2: Covariance/bivariate normal distribution
% Generate and plot some data
clear
mu1 = 10;
mu2 = 20;
MU = [mu1,mu2];
sigma1 = 3;
sigma2 = 1;
rho = 0.7;
SIGMA = [sigma1^2, rho*sigma1*sigma2; ...
         rho*sigma1*sigma2, sigma2^2];
% X1 = mu1+sigma1*randn(100000,1);
% X2 = mu2+rho*sigma2*(X1-mu1)/sigma1+sqrt(1-rho^2)*sigma2*randn(size(X1));
X = mvnrnd(MU,SIGMA,1000);

figure(10)
clf
plot(X(:,1),X(:,2),'.')
axis equal
legend('data')
% Pearson's correlation coefficient (matrix of values):
[R,P,RLO,RUP] = corrcoef(X);
% For user-set confidence intervals:
[R,P,RLO,RUP] = corrcoef(X,'alpha',0.05);
% R gives pairwise correlation coefficients; P is p-value for each; RLO and RUP are confidence intervals for R

% Equivalently, covariance matrix of the samples (differs from R by factors of std(X1) and std(X2))
covmat = cov(X);

% Try linear regressions in each direction
[b12,bint12,~,~,stats12] = regress(X(:,1),[ones(size(X(:,2))),X(:,2)]);
[b21,bint21,~,~,stats21] = regress(X(:,2),[ones(size(X(:,1))),X(:,1)]);
% Note that p-values are all the same (from corrcoef() and both calls of regress()) but lines don't match up
hold all
% Plot the two regression lines
minmax1 = [min(X(:,1)),max(X(:,1))];
minmax2 = [min(X(:,2)),max(X(:,2))];
plot(minmax1,b21(1)+minmax1*b21(2),'--')
plot(b12(1)+minmax2*b12(2),minmax2,'--')
legend({'data','regress 2on1','regress 1on2'})

% Eigenvectors of the covariance matrix (i.e. principal components) yield a more "representative" slope
[V,D] = eig(covmat);
% Care about eigenvector with highest eigenvalue (first principal component), which has slope:
covarslope = V(2,2)/V(1,2);
minmax1 = [min(X(:,1)),max(X(:,1))];
plot(minmax1,(minmax1-mu1)*covarslope+mu2)
legend({'data','regress 2on1','regress 1on2','principal axis'})
% Could also use explicit principal components function
[COEFF,SCORE,latent,tsquare] = princomp(X);
% COEFF contains the same (column) eigenvectors as V, but ordered differently (and perhaps sign-flipped)
% Can also use deming regression (not built-in, programmed by APF) with equal weighting to get the same slope
[b0,b1] = deming(X(:,1),X(:,2),1);
% b1 == covarslope

% EXAMPLE 3: Errors in both variables (not multivariate normal distribution)
% Generate and plot some data
clear
A = (1:0.1:10)';
sigma1 = 3;
sigma2 = 4;
X1 = A+sigma1*randn(size(A));
X2 = 2*A+sigma2*randn(size(A));
figure(20)
clf
plot(X1,X2,'.')
legend('data')
% Try standard linear regressions
[b12,bint12,~,~,stats12] = regress(X1,[ones(size(X2)),X2]);
[b21,bint21,~,~,stats21] = regress(X2,[ones(size(X1)),X1]);
hold all
minmax1 = [min(X1),max(X1)];
minmax2 = [min(X2),max(X2)];
plot(minmax1,b21(1)+minmax1*b21(2),'-')
plot(b12(1)+minmax2*b12(2),minmax2,'--')
legend({'data','regress 2on1','regress 1on2'})
% Neither is quite correct. Proper approach is deming regression (again, not built-in, programmed by APF)
% Third parameter indicates ratio of variances in the two variables - often unknown, but perhaps can be approximated
[b0,b1] = deming(X1,X2,sigma2^2/sigma1^2);
plot(minmax1,b0+b1*minmax1,':')
legend({'data','regress 2on1','regress 1on2','deming'})
% Can get same result a roundabout and non-intuitive way using principal components, when scaling is done appropriately
% Advantage: uses only built-in functions
[V,D] = eig(cov([X1/sigma1,X2/sigma2]));
covarslope = sigma2*V(2,2)/(sigma1*V(1,2));
% b1 == covarslope

% Other functions that may be of use:
%   lscov() - when the variance is non-uniform or the measurements covary
%   fmincon(), fminunc() - for general minimization problems

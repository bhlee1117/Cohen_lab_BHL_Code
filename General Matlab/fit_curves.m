function coeff = fit_curves(y,curves)

%  4-19-2012
%  y ~ [samples, N1] is a set of data curves
%  curves ~ [samples, N2] is a set of component curves for the fit
%  this returns the coefficients of the least squares fit:
%  coeffY = [N2 N1]
%  fits can be recovered as y_fit = curves*coeff

N2 = size(curves,2);

% create orthonormal basis spanning subspace spanned by the curves
w = curves./fillmat2(sqrt(sum(curves.^2,1)),size(curves));

for i = 1:N2
   for j = 1:(i-1)
        w(:,i) = w(:,i) - w(:,j)*sum((w(:,j).*w(:,i)));
        w(:,i) = w(:,i)./sqrt(sum(w(:,i).^2));
   end
end

curvesTOw = w'*curves; % [N2 N2] matrix to convert back from orthogonalized

wTOcurves = curvesTOw^(-1);
coeff = wTOcurves*w'*y;

end
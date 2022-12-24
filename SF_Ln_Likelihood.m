function Ln_Likelihood = SF_Ln_Likelihood(X, Y, Theta)

Beta = Theta(1:3); gamma = Theta(4); sigma2 = Theta(5);
n = size(X, 1);

residual = Y-X*Beta;
Ln_Likelihood = n*log(2/pi)/2 - n*log(sigma2)/2 ...
    +sum(log(normcdf(-residual*sqrt(gamma/(1-gamma))/sqrt(sigma2), 0, 1))) ...
    - residual'*residual/(2*sigma2);
Ln_Likelihood = -Ln_Likelihood/n;
end
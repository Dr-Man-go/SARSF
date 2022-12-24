function Ln_Likelihood = SF_gamma_Ln_Likelihood(gamma, X, Y, Beta, sigma2)
n = size(X, 1);

residual = Y-X*Beta;

Ln_Likelihood = n*log(2/pi)/2 - n*log(sigma2)/2 ...
    +sum(log(normcdf(-residual*sqrt(gamma/(1-gamma))/sqrt(sigma2), 0, 1))) ...
    - residual'*residual/(2*sigma2);
Ln_Likelihood = -Ln_Likelihood/n;
end
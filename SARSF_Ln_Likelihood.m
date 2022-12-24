function Ln_Li = SARSF_Ln_Likelihood(W, X, Y, Theta)
Lambda = Theta(1);           Beta    = Theta(2:size(X,2)+1); 
Gamma  = Theta(size(X,2)+2); Sigma2  = Theta(size(X,2)+3);


n = size(X, 1); S_Lambda = eye(n) - Lambda*W;

residual = Y - Lambda*W*Y - X*Beta;

Ln_Like = n*log(2/pi)/2 - n*log(Sigma2)/2 + log(det(S_Lambda))...
    +sum(log(normcdf(-residual*sqrt(Gamma/(1-Gamma))/sqrt(Sigma2), 0, 1))) ...
    - residual'*residual/(2*Sigma2);
Ln_Li = -Ln_Like/n;
end
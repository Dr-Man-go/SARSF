function Ln_SARSF_Likelihood = MLE_SARSF_concentrate(Lambda, W, X, Y, Precision)
Theta_MLE = MLE_SARSF(W, Lambda, X, Y, Precision);
Ln_SARSF_Likelihood = SARSF_Ln_Likelihood(W, X, Y, Theta_MLE);
end

function Theta_MLE_gamma = MLE_SARSF(W, Lambda, X, Y, Precision)
Theta_MLE_gamma = [Lambda;zeros(5,1)]; Y_L_W_Y = Y - Lambda*W*Y;
Theta_MLE_gamma(2:6) = Corrected_OLS_SF(X, Y_L_W_Y);
rmin=0; rmax=1; options.Display='off'; options.MaxFunEvals=1000;
options.MaxIter=1000;  options.TolX=0.001; options.TolFun=0.001;
Theta_MLE_gamma(5) = fminbnd(@SF_gamma_Ln_Likelihood,rmin,rmax,options, X, Y_L_W_Y,Theta_MLE_gamma(2:4),Theta_MLE_gamma(6));
Theta_MLE_gamma(2:6) = Estimate_SF_Newton(X, Y_L_W_Y, Theta_MLE_gamma(2:6));

Theta_MLE_gamma(5) = fminbnd(@Gamma_SARSF_Ln_Likelihood,rmin,rmax,options,W, X, Y,Theta_MLE_gamma);
Theta_MLE_gamma = SARSF_Gradient_Mango(W, X, Y, Theta_MLE_gamma, Precision);
Theta_MLE_gamma = Grad_SARSF_Mango('SARSF_Ln_Likelihood','SARSF_D_Ln_Likelihood',Theta_MLE_gamma,W, X, Y, Precision);
end


function Ln_Li = Gamma_SARSF_Ln_Likelihood(Gamma, W, X, Y, Theta)
Lambda = Theta(1); Beta = Theta(2:size(X,2)+1); Sigma2 = Theta(size(X,2)+3);
n = size(X, 1); S_Lambda = eye(n) - Lambda*W;

residual = Y - Lambda*W*Y - X*Beta;

Ln_Like = n*log(2/pi)/2 - n*log(Sigma2)/2 + log(det(S_Lambda))...
    +sum(log(normcdf(-residual*sqrt(Gamma/(1-Gamma))/sqrt(Sigma2), 0, 1))) ...
    - residual'*residual/(2*Sigma2);
Ln_Li = -Ln_Like/n;
end
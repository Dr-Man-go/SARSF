function Theta_MLE = MLE_SARSF_Algor_1(W, Theta, X, Y, Precision)
Lambda = Theta(1);
lower  = max(-0.8, Lambda-0.3);upper = min(0.8, Lambda+0.3); options.Display='off'; 
options.MaxFunEvals=1000; options.MaxIter=1000;  options.TolX=0.001; options.TolFun=0.001;
Lambda_opt =  fminbnd(@MLE_SARSF_concentrate,lower,upper,options, W, X, Y, Precision);
Theta_MLE = T_MLE_SARSF_Algor_1(W, [Lambda_opt;zeros(5,1)], X, Y, Precision);
end

function Theta_MLE_gamma = T_MLE_SARSF_Algor_1(W, C2SLSE, X, Y, Precision)
Theta_MLE_gamma = C2SLSE;
lambda = C2SLSE(1); Y_L_W_Y = Y - lambda*W*Y;
Theta_MLE_gamma(2:6) = Corrected_OLS_SF(X, Y_L_W_Y);
rmin=0; rmax=1; options.Display='off'; options.MaxFunEvals=1000;
options.MaxIter=1000;  options.TolX=0.001; options.TolFun=0.001;
Theta_MLE_gamma(5) = fminbnd(@SF_gamma_Ln_Likelihood,rmin,rmax,options, X, Y_L_W_Y,Theta_MLE_gamma(2:4),Theta_MLE_gamma(6));
Theta_MLE_gamma(2:6) = Estimate_SF_Newton(X, Y_L_W_Y, Theta_MLE_gamma(2:6));
Theta_MLE_gamma(5)   = fminbnd(@Find_Gamma_SARSF_Ln_Likelihood,rmin,rmax,options,W, X, Y,Theta_MLE_gamma);
Theta_MLE_gamma = SARSF_Gradient_Mango(W, X, Y, Theta_MLE_gamma, Precision);
end


function Ln_Li = Find_Gamma_SARSF_Ln_Likelihood(Gamma, W, X, Y, Theta)
Lambda = Theta(1); Beta = Theta(2:size(X,2)+1); Sigma2 = Theta(size(X,2)+3);

n = size(X, 1); S_Lambda = eye(n) - Lambda*W;

residual = Y - Lambda*W*Y - X*Beta;

Ln_Like = n*log(2/pi)/2 - n*log(Sigma2)/2 + log(det(S_Lambda))...
    +sum(log(normcdf(-residual*sqrt(Gamma/(1-Gamma))/sqrt(Sigma2), 0, 1))) ...
    - residual'*residual/(2*Sigma2);
Ln_Li = -Ln_Like/n;
end
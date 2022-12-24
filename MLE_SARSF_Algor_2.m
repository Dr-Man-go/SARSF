function Theta_MLE_gamma = MLE_SARSF_Algor_2(W, C2SLSE, X, Y, Precision)
 %temp = Corrected_BO_C2SOLS_SARSF(C2SLSE, W, X, Y);
 %lambda = temp(1);
Theta_MLE_gamma = C2SLSE;
 
lambda = C2SLSE(1); Y_L_W_Y = Y - lambda*W*Y;
Theta_MLE_gamma(2:6) = Corrected_OLS_SF(X, Y_L_W_Y);
rmin=0; rmax=1; options.Display='off'; options.MaxFunEvals=1000;
options.MaxIter=1000;  options.TolX=0.001; options.TolFun=0.001;
Theta_MLE_gamma(5) = fminbnd(@SF_gamma_Ln_Likelihood,rmin,rmax,options, X, Y_L_W_Y,Theta_MLE_gamma(2:4),Theta_MLE_gamma(6));
%Theta_MLE_gamma(5) = fminbnd(@Gamma_SARSF_Ln_Likelihood,rmin,rmax,options,W, X, Y,Theta_MLE_gamma);

Theta_MLE_gamma(2:6) = Estimate_SF_Newton(X, Y_L_W_Y, Theta_MLE_gamma(2:6));

Theta_MLE_gamma(5) = fminbnd(@Gamma_SARSF_Ln_Likelihood,rmin,rmax,options,W, X, Y,Theta_MLE_gamma);
%Theta_MLE_gamma(1) = fminbnd(@Lambda_SARSF_Ln_Likelihood,-1,1,options,W, X, Y,Theta_MLE_gamma);

Theta_MLE_gamma = SARSF_Gradient_Mango(W, X, Y, Theta_MLE_gamma, Precision);
Theta_MLE_gamma = Grad_SARSF_Mango('SARSF_Ln_Likelihood','SARSF_D_Ln_Likelihood',Theta_MLE_gamma,W, X, Y, Precision);
end


function Ln_Li = Gamma_SARSF_Ln_Likelihood(Gamma, W, X, Y, Theta)
Lambda = Theta(1);           Beta    = Theta(2:size(X,2)+1); 
Sigma2  = Theta(size(X,2)+3);


n = size(X, 1); S_Lambda = eye(n) - Lambda*W;

residual = Y - Lambda*W*Y - X*Beta;

Ln_Like = n*log(2/pi)/2 - n*log(Sigma2)/2 + log(det(S_Lambda))...
    +sum(log(normcdf(-residual*sqrt(Gamma/(1-Gamma))/sqrt(Sigma2), 0, 1))) ...
    - residual'*residual/(2*Sigma2);
Ln_Li = -Ln_Like/n;
end
function D_Ln =SF_D_Ln_Likelihood(X, Y, Theta)
n = length(Y);
D_Ln = zeros(size(X,2)+2, 1);


Beta  = Theta(1:3); Residual = Y - X*Beta;
Gamma = Theta(4);     Sigma2 = Theta(5); 

temp = sqrt(Gamma/(1-Gamma))/sqrt(Sigma2);
f = Cal_phi_divide_Phi(-Residual*temp);

D_Ln(1:size(X,2)) =   X'*Residual/Sigma2 + X'*f*temp;
D_Ln(size(X,2)+1) = - f'*Residual/(2*sqrt(Sigma2*Gamma*(1-Gamma))*(1-Gamma));
D_Ln(size(X,2)+2) = - n/(2*Sigma2) + Residual'*Residual/(2*Sigma2^2) + Residual'*f*temp/(2*Sigma2);
D_Ln = -D_Ln/n;
end
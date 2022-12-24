function D_Ln = SARSF_D_Ln_Likelihood(W, X, Y, Theta)
D_Ln = zeros(size(X,2)+3, 1);


Lambda = Theta(1);           Beta    = Theta(2:size(X,2)+1); 
Gamma  = Theta(size(X,2)+2); Sigma2  = Theta(size(X,2)+3);

n = size(X, 1); S_Lambda = eye(n) - Lambda*W; G_Lambda = W/S_Lambda;
WY = W*Y; Residual = Y - Lambda*WY - X*Beta;

temp = sqrt(Gamma/(1-Gamma))/sqrt(Sigma2);
f = Cal_phi_divide_Phi(-Residual*temp);

D_Ln(1) = -trace(G_Lambda) + WY'*Residual/Sigma2 + temp*WY'*f;
D_Ln(2:size(X,2)+1) =   X'*Residual/Sigma2 + X'*f*temp;
D_Ln(size(X,2)+2) = - f'*Residual/(2*sqrt(Sigma2*Gamma*(1-Gamma)^3));
D_Ln(size(X,2)+3) = - n/(2*Sigma2) + Residual'*Residual/(2*Sigma2^2) + Residual'*f*temp/(2*Sigma2);
D_Ln = -D_Ln/n;
end


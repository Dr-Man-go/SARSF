function D2_Ln = SF_D2_Ln_Likelihood(X, Y, Theta)
n = length(Y);
D2_Ln = zeros(size(X,2)+2, size(X,2)+2);


Beta  = Theta(1:3); Residual = Y - X*Beta;
Gamma = Theta(4);     Sigma2 = Theta(5); 

temp = sqrt(Gamma/(1-Gamma))/sqrt(Sigma2);
[f, F] = Cal_D_phi_divide_Phi(-Residual*temp);

% (beta, beta)
D2_Ln(1:size(X,2), 1:size(X,2)) =   X'*diag(F)*X*temp^2 - X'*X/Sigma2;

% (beta, gamma)
D2_Ln(1:size(X,2), 1+size(X,2)) =  X'*f/(2*sqrt(Sigma2*Gamma*(1-Gamma))*(1-Gamma)) - ...
    X'*(Residual.*F)/(2*Sigma2*(1-Gamma)^2);
D2_Ln(1+size(X,2), 1:size(X,2)) =  D2_Ln(1:size(X,2), 1+size(X,2))';

% (beta, sigma2)
D2_Ln(1:size(X,2), 2+size(X,2)) = - X'*Residual/Sigma2^2 - X'*f*temp/(2*Sigma2) +temp^2*X'*(Residual.*F) /(2*Sigma2);
D2_Ln(2+size(X,2), 1:size(X,2)) =  D2_Ln(1:size(X,2), 2+size(X,2))';


% (gamma, gamma)
temp1 = (Gamma*(1-Gamma))^(3/2); temp2 = sqrt(Gamma)*(1-Gamma)^(5/2);
D2_Ln(1+size(X,2), 1+size(X,2)) = (1/temp1 - 3/temp2)*Residual'*f/(4*sqrt(Sigma2))...
    + (Residual.*Residual)'*F/(4*Sigma2*Gamma*(1-Gamma)^3);

% (sigma2, sigma2)
D2_Ln(2+size(X,2), 2+size(X,2)) = n/(2*Sigma2^2) - Residual'*Residual/Sigma2^3 - ...
    3*temp*Residual'*f/(4*Sigma2^2) + (Residual.*Residual)'*F*temp^2/(4*Sigma2^2);

% (sigma2, gamma)
D2_Ln(2+size(X,2), 1+size(X,2)) =  Residual'*f/(4*sqrt(Sigma2)^3*sqrt(Gamma*(1-Gamma)^3))...
    -(Residual.*Residual)'*F/(4*Sigma2^2*(1-Gamma)^2);
D2_Ln(2+size(X,2), 1+size(X,2)) =  D2_Ln(1+size(X,2), 2+size(X,2));

D2_Ln = -D2_Ln/n;
end
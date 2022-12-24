function Theta_COLS = Corrected_OLS_SF(X, Y)
Beta_OLS = (X'*X)\X'*Y; residual = Y-X*Beta_OLS;
m2   =  sum(residual.^2)/length(Y);
m3   =  sum(residual.^3)/length(Y); m3 = m3*(m3<0);

Temp_COLS   = (sqrt(pi/2)*pi*m3/(pi-4))^(2/3);
sigma2_COLS = m2 + 2*Temp_COLS/pi;
Gamma_COLS  = Temp_COLS/sigma2_COLS; Gamma_COLS = Gamma_COLS*(Gamma_COLS<=1) + (Gamma_COLS>1);
Beta_OLS(1) = Beta_OLS(1) + sqrt(2*Gamma_COLS*sigma2_COLS/pi);
Theta_COLS = [Beta_OLS; Gamma_COLS; sigma2_COLS];

end



function Theta_MLE = Estimate_SF_Newton(X, Y, Theta_COLS)
Theta = Theta_COLS;
D_Ln = SF_D_Ln_Likelihood(X, Y, Theta);

for i = 2:2000
    D2_Ln     = SF_D2_Ln_Likelihood(X, Y, Theta(:,i-1));
    if (sum(sum(isnan(D2_Ln)))>0) || (sum(sum(isinf(D2_Ln)))>0)
        Theta_MLE = Theta(:,i-1);
        return;
    end
        
    direction = pinv(D2_Ln)*D_Ln;
    Theta(:,i)  = Theta(:,i-1) - direction;
    if Theta(4,i)>1 || Theta(4,i)<0 || Theta(5,i)<0
        %fprintf('B')
        Theta_MLE = Theta(:,i-1);
        return;
    else    
        D_Ln = SF_D_Ln_Likelihood(X, Y, Theta(:, i));
        if norm(D_Ln)<1e-6
            
            %fprintf('O');disp(norm(D_Ln));
            Theta_MLE = Theta(:,i);
            return;
        end
        
    end
end
Theta_MLE = Theta(:,i);
end






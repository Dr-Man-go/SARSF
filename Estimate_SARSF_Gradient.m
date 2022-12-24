function Theta_MLE = Estimate_SARSF_Gradient(W, X, Y, Theta_MLE)
Theta = Theta_MLE;
new_direction = D_SARSF_Ln_Likelihood(W, X, Y, Theta);

for i = 2:2000
    old = SARSF_Ln_Likelihood(W, X, Y, Theta(:,i-1));
    old_norm = norm(D_SARSF_Ln_Likelihood(W, X, Y, Theta(:,i-1)));
    Theta(:,i)  = Theta(:,i-1) - 0.1*new_direction;
    if Theta(size(X,2)+2,i)>1 || Theta(size(X,2)+2,i)<0 || Theta(size(X,2)+3,i)<=0
        %fprintf('B')
        Theta_MLE = Theta(:,i-1);
        break
    end
    
    new_direction = D_SARSF_Ln_Likelihood(W, X, Y, Theta(:, i));
    new_norm = norm(new_direction);
    new = SARSF_Ln_Likelihood(W, X, Y, Theta(:,i));
    
    if  new-old>=0 || new_norm-old_norm>=0
        %fprintf('B')
        Theta_MLE = Theta(:,i-1);
        break
    end
    
    if norm(new_direction)<1e-4
        %fprintf('O');disp(norm(new_direction));
        Theta_MLE = Theta(:,i);
        break
    end
    
end
end


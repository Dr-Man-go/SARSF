function Theta_MLE = SARSF_Gradient_Mango(W, X, Y, Theta_MLE, precise)
Theta = Theta_MLE;
new_direction = SARSF_D_Ln_Likelihood(W, X, Y, Theta);

for i = 2:1000
    old = SARSF_Ln_Likelihood(W, X, Y, Theta);
    old_norm = norm(new_direction);
    if old_norm<precise
        Theta_MLE = Theta;
        break
    end
    
    for m = 10 : (-1) : 0
        temp = Theta - 0.1*m*new_direction;
        if (temp(5)>0) && (temp(5)<1) && (temp(6)>0) && (temp(6)<11) && (SARSF_Ln_Likelihood(W, X, Y, temp)<old)
            
            new_direction = SARSF_D_Ln_Likelihood(W, X, Y, temp);
            new_norm = norm(new_direction);
            %new_norm-old_norm
            if new_norm<old_norm
                Theta  = temp;
                break;
            end
        end
    end
    if m==0
        Theta_MLE = Theta;
        return;
    end
    
    
    %new_direction = D_SARSF_Ln_Likelihood(W, X, Y, Theta);
    %     new_norm = norm(new_direction);
    %
    %     if  new_norm>old_norm
    %         new_norm - old_norm
    %         Theta_MLE = Theta + 0.1*m*new_direction
    %         break;
    %     end
    %
    
end
end


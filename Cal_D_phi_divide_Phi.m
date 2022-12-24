function [f, F] = Cal_D_phi_divide_Phi(t)
f = Cal_phi_divide_Phi(t);
F = -t.*f - f.^2;
end
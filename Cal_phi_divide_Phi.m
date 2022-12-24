function f = Cal_phi_divide_Phi(t)
f = normpdf(t, 0, 1)./normcdf(t, 0, 1);
end
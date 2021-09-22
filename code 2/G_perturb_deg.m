function G = G_perturb_deg(rho,eps,Kraus)
%This function is used to calculte the hermitian conjugate of the perturbed G 
%We omit the tau part of G^daggar since it is offset in the expression
G = (1-eps)*(Kraus'*rho*Kraus);
end

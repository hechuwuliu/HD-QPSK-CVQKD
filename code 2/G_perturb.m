function G = G_perturb(rho,eps,Kraus,N_cutoff)
%This function is used to calculte the perturbed G 
tau =1/(2*4*(N_cutoff+1))*eye(2*4*(N_cutoff+1));
G = (1-eps)*(Kraus*rho*Kraus')+eps*tau;
end
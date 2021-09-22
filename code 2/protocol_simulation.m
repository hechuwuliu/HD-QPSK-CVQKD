clear all
%some necessary parameters
%Excess noise
xi = 0.002;
%photon number cutoff
N_cutoff = 12;
%p is actually pA/2 in the PRX paper
p = 0.25;
%parameter related to postselection
delta_c = 0;
%transmission distance
L = 350;
%transmission efficiency
eta = 10^(-0.02*L);
%amplitude of states
alpha = 0.65;
%'eps2' is used in construction of the gradient of keyrate function
eps2 = 10^(-12);
%dim_G is dimension of space after G map operating on the original space
dim_G = 2*4*(N_cutoff+1);

%Kraus operator
[Kraus_p,Kraus_q] = K_pq_new(N_cutoff,delta_c);
%Register part of Kraus operator
z0 = [1;0];
z1 = [0;1];
%pinching channel operator
Z_0 = kron(z0*z0',eye(4*(N_cutoff+1)));
Z_1 = kron(z1*z1',eye(4*(N_cutoff+1)));
%Observables
Gamma = Gamma(N_cutoff);
%Expected value
[gamma,~] = gamma_exp(alpha,L,xi,p);

%step I
rho0 = search_rho0(N_cutoff,Gamma,gamma);
[rho,Tr] = algorithm1(N_cutoff,rho0,Kraus_p,Kraus_q,Z_0,Z_1,gamma,Gamma);

%step II
eps = abstractep(rho,Gamma,gamma);
zeta = 4*eps2*(dim_G-1)*log2(dim_G/(eps2*(dim_G-1)));
grad_f_p = (G_perturb_deg(logm(G_perturb(rho,eps2,Kraus_p,N_cutoff))/log(2),eps2,Kraus_p)-...
        G_perturb_deg(logm(Z_0*G_perturb(rho,eps2,Kraus_p,N_cutoff)*Z_0+Z_1*G_perturb(rho,eps2,Kraus_p,N_cutoff)*Z_1)/log(2),eps2,Kraus_p)).';
grad_f_q = (G_perturb_deg(logm(G_perturb(rho,eps2,Kraus_q,N_cutoff))/log(2),eps2,Kraus_q)-...
        G_perturb_deg(logm(Z_0*G_perturb(rho,eps2,Kraus_q,N_cutoff)*Z_0+Z_1*G_perturb(rho,eps2,Kraus_q,N_cutoff)*Z_1)/log(2),eps2,Kraus_q)).';
grad_f = grad_f_p + grad_f_q;
fval = max_gamma_y(N_cutoff,Gamma,grad_f,gamma,eps);
beta = rel_ent(G_perturb(rho,eps2,Kraus_p,N_cutoff)...
        ,Z_0*G_perturb(rho,eps2,Kraus_p,N_cutoff)*Z_0+Z_1*G_perturb...
        (rho,eps2,Kraus_p,N_cutoff)*Z_1)+...
        rel_ent(G_perturb(rho,eps2,Kraus_q,N_cutoff)...
        ,Z_0*G_perturb(rho,eps2,Kraus_q,N_cutoff)*Z_0+Z_1*G_perturb...
        (rho,eps2,Kraus_q,N_cutoff)*Z_1)-trace(rho.'*grad_f) + fval;

%error correction
delta_EC = error_correction(eta,xi,delta_c,alpha);

%tight lower bound of key rate
keyrate = real(0.5*(beta - zeta - delta_EC));

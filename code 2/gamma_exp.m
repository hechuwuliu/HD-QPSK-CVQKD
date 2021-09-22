function [expc,expt] = gamma_exp(alpha,L,ksi,p)
%Calculating expected value
alpha1 = exp(1i*0.25*pi)*alpha;
alpha2 = exp(1i*0.75*pi)*alpha;
alpha3 = exp(1i*1.25*pi)*alpha; 
alpha4 = exp(1i*1.75*pi)*alpha;
eta = 10^(-0.02*L);
exp_q_1 = sqrt(2*eta)*real(alpha1);
exp_q_2 = sqrt(2*eta)*real(alpha2);
exp_q_3 = sqrt(2*eta)*real(alpha3);
exp_q_4 = sqrt(2*eta)*real(alpha4);
exp_p_1 = sqrt(2*eta)*imag(alpha1);
exp_p_2 = sqrt(2*eta)*imag(alpha2);
exp_p_3 = sqrt(2*eta)*imag(alpha3);
exp_p_4 = sqrt(2*eta)*imag(alpha4);
exp_n_1 = eta*(abs(alpha1))^2+eta*ksi*0.5;
exp_n_2 = eta*(abs(alpha2))^2+eta*ksi*0.5;
exp_n_3 = eta*(abs(alpha3))^2+eta*ksi*0.5;
exp_n_4 = eta*(abs(alpha4))^2+eta*ksi*0.5;
exp_d_1 = eta*(alpha1^2+(conj(alpha1))^2);
exp_d_2 = eta*(alpha2^2+(conj(alpha2))^2);
exp_d_3 = eta*(alpha3^2+(conj(alpha3))^2);
exp_d_4 = eta*(alpha4^2+(conj(alpha4))^2);
innp_21 = exp(-0.5*(abs(alpha2)^2+abs(alpha1)^2)+conj(alpha2)*alpha1);
innp_31 = exp(-0.5*(abs(alpha3)^2+abs(alpha1)^2)+conj(alpha3)*alpha1);
innp_41 = exp(-0.5*(abs(alpha4)^2+abs(alpha1)^2)+conj(alpha4)*alpha1);
innp_32 = exp(-0.5*(abs(alpha3)^2+abs(alpha2)^2)+conj(alpha3)*alpha2);
innp_42 = exp(-0.5*(abs(alpha4)^2+abs(alpha2)^2)+conj(alpha4)*alpha2);
innp_43 = exp(-0.5*(abs(alpha4)^2+abs(alpha3)^2)+conj(alpha4)*alpha3);
innp_12 = conj(innp_21);
innp_13 = conj(innp_31); 
innp_14 = conj(innp_41);
innp_23 = conj(innp_32);
innp_24 = conj(innp_42);
innp_34 = conj(innp_43);
expc = [1 p*exp_q_1 p*exp_q_2 (0.5-p)*exp_q_3 (0.5-p)*exp_q_4 p*exp_p_1 ...
    p*exp_p_2 (0.5-p)*exp_p_3 (0.5-p)*exp_p_4 p*exp_n_1 p*exp_n_2 (0.5-p)*exp_n_3 ...
   (0.5-p)*exp_n_4 p*exp_d_1 p*exp_d_2 (0.5-p)*exp_d_3 (0.5-p)*exp_d_4 ...
    p p*(innp_21+innp_12) 1i*p*(innp_21-innp_12) sqrt(p*(0.5-p))*(innp_31+innp_13) ...
    1i*sqrt(p*(0.5-p))*(innp_31-innp_13) sqrt(p*(0.5-p))*(innp_41+innp_14)...
    1i*sqrt(p*(0.5-p))*(innp_41-innp_14) p sqrt(p*(0.5-p))*(innp_32+innp_23) 1i*sqrt(p*(0.5-p))*(innp_32-innp_23)...
    sqrt(p*(0.5-p))*(innp_42+innp_24) 1i*sqrt(p*(0.5-p))*(innp_42-innp_24) 0.5-p...
    (0.5-p)*(innp_43+innp_34) 1i*(0.5-p)*(innp_43-innp_34)...
    0.5-p];
expt = [p*exp_q_1;p*exp_q_2;(0.5-p)*exp_q_3;(0.5-p)*exp_q_4;p*exp_p_1; ...
    p*exp_p_2;(0.5-p)*exp_p_3;(0.5-p)*exp_p_4;p*exp_n_1;p*exp_n_2;(0.5-p)*exp_n_3; ...
   (0.5-p)*exp_n_4;p*exp_d_1;p*exp_d_2;(0.5-p)*exp_d_3;(0.5-p)*exp_d_4; ...
    p*(innp_21+innp_12);1i*p*(innp_21-innp_12);sqrt(p*(0.5-p))*(innp_31+innp_13); ...
    1i*sqrt(p*(0.5-p))*(innp_31-innp_13);sqrt(p*(0.5-p))*(innp_41+innp_14);...
    1i*sqrt(p*(0.5-p))*(innp_41-innp_14);sqrt(p*(0.5-p))*(innp_32+innp_23);1i*sqrt(p*(0.5-p))*(innp_32-innp_23);...
    sqrt(p*(0.5-p))*(innp_42+innp_24);1i*sqrt(p*(0.5-p))*(innp_42-innp_24);...
    (0.5-p)*(innp_43+innp_34);1i*(0.5-p)*(innp_43-innp_34)].';
end


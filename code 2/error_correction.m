function ec = error_correction(eta,ksi,delta_c,alpha)
%this function calculates the error correction
beta = 0.95;

d0 = @(x) (1/sqrt(pi*(eta*ksi+1)))*exp(-(x-sqrt(eta)*alpha).^2/(eta*ksi+1));
d1 = @(x) (1/sqrt(pi*(eta*ksi+1)))*exp(-(x+sqrt(eta)*alpha).^2/(eta*ksi+1));

i0_non = integral(@(x)d0(x),-delta_c,delta_c);
i0_0 = integral(@(x)d0(x),delta_c,Inf);
i0_1 = integral(@(x)d0(x),-Inf,-delta_c);

i1_bot = integral(@(x)d1(x),-delta_c,delta_c);
i1_0 = integral(@(x)d1(x),delta_c,Inf);
i1_1 = integral(@(x)d1(x),-Inf,-delta_c);

pass_ri = 1-i0_non;
pass_le = 1-i1_bot;
pass_tot = 0.5*pass_ri+0.5*pass_le;

h = @(y,z) -y*log2(y) - (z)*log2(z);
HZX = 0.5*h(i0_0/pass_ri,i0_1/pass_ri)+0.5*h(i1_0/pass_le,i1_1/pass_le);

Pz0 = 0.5*(i1_0/pass_le+i0_0/pass_ri);
Pz1 = 0.5*(i0_1/pass_ri+i1_1/pass_le);
HZ = h(Pz0,Pz1);
d = (1-beta)*HZ+beta*HZX;
ec = 2*pass_tot*d;
end


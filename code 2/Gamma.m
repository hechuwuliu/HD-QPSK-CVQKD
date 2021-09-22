function [A] = Gamma(N_cutoff)
d1 = diag([1,0,0,0]);
d2 = diag([0,1,0,0]);
d3 = diag([0,0,1,0]);
d4 = diag([0,0,0,1]);
%a and a_deg are annihilation and creation operators
a = zeros(N_cutoff+1);
for i = 1:N_cutoff
    a(i,i+1) = sqrt(i);
end
a_deg = zeros(N_cutoff+1);
for i = 1:N_cutoff
    a_deg(i+1,i) = sqrt(i);
end
q_op = (a + a_deg)/sqrt(2);
p_op = (a_deg - a)*1i/sqrt(2);
n_op = a_deg*a;
d_op = a_deg*a_deg + a*a;
%construct matrices for partial trace condition
id = eye(N_cutoff+1);
Z = zeros(N_cutoff+1);
I_11 = [id Z Z Z;Z Z Z Z;Z Z Z Z;Z Z Z Z];
I_12 = [Z id Z Z;Z Z Z Z;Z Z Z Z;Z Z Z Z];
I_13 = [Z Z id Z;Z Z Z Z;Z Z Z Z;Z Z Z Z];
I_14 = [Z Z Z id;Z Z Z Z;Z Z Z Z;Z Z Z Z];
I_21 = [Z Z Z Z;id Z Z Z;Z Z Z Z;Z Z Z Z];
I_22 = [Z Z Z Z;Z id Z Z;Z Z Z Z;Z Z Z Z];
I_23 = [Z Z Z Z;Z Z id Z;Z Z Z Z;Z Z Z Z];
I_24 = [Z Z Z Z;Z Z Z id;Z Z Z Z;Z Z Z Z];
I_31 = [Z Z Z Z;Z Z Z Z;id Z Z Z;Z Z Z Z];
I_32 = [Z Z Z Z;Z Z Z Z;Z id Z Z;Z Z Z Z];
I_33 = [Z Z Z Z;Z Z Z Z;Z Z id Z;Z Z Z Z];
I_34 = [Z Z Z Z;Z Z Z Z;Z Z Z id;Z Z Z Z];
I_41 = [Z Z Z Z;Z Z Z Z;Z Z Z Z;id Z Z Z];
I_42 = [Z Z Z Z;Z Z Z Z;Z Z Z Z;Z id Z Z];
I_43 = [Z Z Z Z;Z Z Z Z;Z Z Z Z;Z Z id Z];
I_44 = [Z Z Z Z;Z Z Z Z;Z Z Z Z;Z Z Z id];
A = {eye(4*(N_cutoff+1)).', kron(d1,q_op).', kron(d2,q_op).', kron(d3,q_op).', kron(d4,q_op).',...
    kron(d1,p_op).', kron(d2,p_op).', kron(d3,p_op).', kron(d4,p_op).',...
    kron(d1,n_op).', kron(d2,n_op).', kron(d3,n_op).', kron(d4,n_op).',...
    kron(d1,d_op).', kron(d2,d_op).', kron(d3,d_op).', kron(d4,d_op).',...
    I_11.', (I_21+I_12).', 1i*(I_21-I_12).', (I_31+I_13).', 1i*(I_31-I_13).',...
    (I_41+I_14).', 1i*(I_41-I_14).', I_22.', (I_32+I_23).', 1i*(I_32-I_23).',...
    (I_42+I_24).', 1i*(I_42-I_24).', I_33.', (I_43+I_34).', 1i*(I_43-I_34).',...
    I_44.'};
end


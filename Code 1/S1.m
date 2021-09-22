function [Cons] = S1(Nc)
%%% the constraints are hermitian but not orthogonal
%
% Nc is the photon number cut-off
%
%%% operator for Alice
bra0 = [1,0,0,0];
bra1 = [0,1,0,0];
bra2 = [0,0,1,0];
bra3 = [0,0,0,1];
ket0 = bra0';
ket1 = bra1';
ket2 = bra2';
ket3 = bra3';
A00 = ket0*bra0;
A01 = ket1*bra0 + ket0*bra1;  %01+10
A02 = ket2*bra0 + ket0*bra2;  %02+20
A03 = ket3*bra0 + ket0*bra3;  %03+30
A10 = 1i*ket1*bra0 - 1i*ket0*bra1;  %i01-i10 (ij of matrix A will extract ji block of rho)
A11 = ket1*bra1;
A12 = ket2*bra1 + ket1*bra2;  %12+21
A13 = ket3*bra1 + ket1*bra3;  %13+31
A20 = 1i*ket2*bra0 - 1i*ket0*bra2;  %i02-i20
A21 = 1i*ket2*bra1 - 1i*ket1*bra2;  %i12-i21
A22 = ket2*bra2;
A23 = ket3*bra2 + ket2*bra3;  %23+32
A30 = 1i*ket3*bra0 - 1i*ket0*bra3;  %i03-i30
A31 = 1i*ket3*bra1 - 1i*ket1*bra3;  %i13-i31
A32 = 1i*ket3*bra2 - 1i*ket2*bra3;  %i23-i32
A33 = ket3*bra3;

%%% operator for Bob
% creation and annihilation
hat_a = zeros(Nc+1,Nc+1);
for j = 1:Nc
    hat_a(j,j+1) = sqrt(j);
end
hat_a_dagger = hat_a';
% 4 measurement outcomes operators
hat_q = (hat_a_dagger + hat_a)/sqrt(2);
hat_p = 1i*(hat_a_dagger - hat_a)/sqrt(2);
hat_n = hat_a_dagger * hat_a;
hat_d = hat_a_dagger * hat_a_dagger + hat_a * hat_a;

%%% operator for AB
Cons = cell(32,1);
% constraints related to the experimental values / expected values
Cons{1,1}=kron(A00,hat_q);  %hat_q0: Alice sends state alpha*exp(1i*pi/4) and Bob measures operator q
Cons{2,1}=kron(A11,hat_q);  %hat_q1: Alice sends state alpha*exp(3i*pi/4) and Bob measures operator q
Cons{3,1}=kron(A22,hat_q);  %hat_q2: Alice sends state alpha*exp(5i*pi/4) and Bob measures operator q
Cons{4,1}=kron(A33,hat_q);  %hat_q3: Alice sends state alpha*exp(7i*pi/4) and Bob measures operator q
Cons{5,1}=kron(A00,hat_p);  %hat_p0: Alice sends state alpha*exp(1i*pi/4) and Bob measures operator p
Cons{6,1}=kron(A11,hat_p);  %hat_p1
Cons{7,1}=kron(A22,hat_p);  %hat_p2
Cons{8,1}=kron(A33,hat_p);  %hat_p3
Cons{9,1}=kron(A00,hat_d);  %hat_d0
Cons{10,1}=kron(A11,hat_d);  %hat_d1
Cons{11,1}=kron(A22,hat_d);  %hat_d2
Cons{12,1}=kron(A33,hat_d);  %hat_d3
Cons{13,1}=kron(A00,hat_n);  %hat_n0
Cons{14,1}=kron(A11,hat_n);  %hat_n1
Cons{15,1}=kron(A22,hat_n);  %hat_n2
Cons{16,1}=kron(A33,hat_n);  %hat_n3
% constraints related to the completely positive and trace-preserving mapping
Cons{17,1}=kron(A01,eye(Nc+1));  %01 block + 10 block
Cons{18,1}=kron(A10,eye(Nc+1));  %i01 block -i10 block 
Cons{19,1}=kron(A02,eye(Nc+1));  %02 block +20 block 
Cons{20,1}=kron(A20,eye(Nc+1));  %i02 block -i20 block 
Cons{21,1}=kron(A03,eye(Nc+1));  %03 block +30 block 
Cons{22,1}=kron(A30,eye(Nc+1));  %i03 block -i30 block 
Cons{23,1}=kron(A12,eye(Nc+1));  %12 block +21 block 
Cons{24,1}=kron(A21,eye(Nc+1));  %i12 block -i21 block 
Cons{25,1}=kron(A13,eye(Nc+1));  %13 block +31 block 
Cons{26,1}=kron(A31,eye(Nc+1));  %i13 block -i31 block 
Cons{27,1}=kron(A23,eye(Nc+1));  %23 block +32 block 
Cons{28,1}=kron(A32,eye(Nc+1));  %i23 block -i32 block 
Cons{29,1}=kron(A00,eye(Nc+1));  %00trace
Cons{30,1}=kron(A11,eye(Nc+1));  %11trace
Cons{31,1}=kron(A22,eye(Nc+1));  %22trace
Cons{32,1}=kron(A33,eye(Nc+1));  %33trace
end


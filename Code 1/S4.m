function [Z0Kq,Z1Kq,Z0Kp,Z1Kp] = S4(Nc,delta)
%
% calculate the density matrices of interval operator I0, I1 of measuring q and p under photon number representation
% Z0K means Z0*K
% K = Z0K+Z1K
%
% Nc is the photon number cut-off
% delta is the post-selection parameter
%
syms x
y = sym(zeros(1,Nc+1));
y(1) = exp(-0.5*x^2); %<0|q> photon number 0
for i = 2:Nc+1
    y(i) = x*y(i-1) - diff(y(i-1),x,1);
    y(i) = collect(y(i),x);
end
I0q = zeros(Nc+1,Nc+1);
I1q = zeros(Nc+1,Nc+1);
I0p = zeros(Nc+1,Nc+1);
I1p = zeros(Nc+1,Nc+1);
for i = 1:Nc+1
    for j = i:Nc+1
        effq = 1/sqrt(2^(i+j-2)*factorial(i-1)*factorial(j-1)*pi);
        effp = (-1)^j*(1i)^(i+j)/sqrt(2^(i+j-2)*factorial(i-1)*...
            factorial(j-1)*pi); % i is row£¬j is column£¬|i-1><j-1| element of density matrix under photon number representation
        f = collect(y(i)*y(j),x);
        I0q(i,j)=int(f,x,delta,inf)*effq;
        I0q(j,i)=I0q(i,j);
        I1q(i,j)=(-1)^(i+j)*I0q(i,j);
        I1q(j,i)=I1q(i,j);
        I0p(i,j)=int(f,x,delta,inf)*effp;
        I0p(j,i)=I0p(i,j)*(-1)^(i-j);
        I1p(i,j)=(-1)^(i+j)*I0p(i,j);
        I1p(j,i)=I1p(i,j)*(-1)^(i-j);
    end
end
sqI0q = sqrtm(I0q);
sqI1q = sqrtm(I1q);
sqI0p = sqrtm(I0p);
sqI1p = sqrtm(I1p);
R0 = [1;0];
R1 = [0;1];
A = eye(4);
Z0Kq = kron(R0,kron(A,sqI0q));
Z1Kq = kron(R1,kron(A,sqI1q));
Z0Kp = kron(R0,kron(A,sqI0p));
Z1Kp = kron(R1,kron(A,sqI1p));
end


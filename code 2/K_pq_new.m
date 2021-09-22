function [Kp,Kq] = K_pq_new(Nc,delta)
%计算当Bob测量正则坐标基矢时的Kraus算子
%interval operator
syms x
y = sym(zeros(1,Nc+1));
y(1) = exp(-0.5*x^2);
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
            factorial(j-1)*pi); % i-1为横坐标，j-1为纵坐标，|i-1><j-1|
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
%'A' is Alice's part of the Kraus operator
A = eye(4);
%'z1' and 'z0' is the register part of Kraus operator
z0 = [1;0];
z1 = [0;1];
Kp = kron(kron(z0,A),sqrtm(I0p))+kron(kron(z1,A),sqrtm(I1p));
Kq = kron(kron(z0,A),sqrtm(I0q))+kron(kron(z1,A),sqrtm(I1q));
end
function [D2] = S7(rho0,Z0Kq,Z1Kq,Z0Kp,Z1Kp,Con,Num,Nc,ep1,ep2)
%
% dual problem
%
Num2 = [Num;1];

N = 4*(Nc+1);
entg = eye(2*N)/(2*N);
% for measuring q
Kq = Z0Kq+Z1Kq;
Gepq = (1-ep1)*Kq*rho0*Kq' + ep1*entg;
ZGepq = (1-ep1)*(Z0Kq*rho0*Z0Kq'+Z1Kq*rho0*Z1Kq') + ep1*entg;
fepq = trace(Gepq*logm(Gepq)/log(2))-trace(Gepq*logm(ZGepq)/log(2));
% for measuring p
Kp = Z0Kp+Z1Kp;
Gepp = (1-ep1)*Kp*rho0*Kp' + ep1*entg;
ZGepp = (1-ep1)*(Z0Kp*rho0*Z0Kp'+Z1Kp*rho0*Z1Kp') + ep1*entg;
fepp = trace(Gepp*logm(Gepp)/log(2))-trace(Gepp*logm(ZGepp)/log(2));
% merge
fep = real(fepq+fepp);
zeta = 4*ep1*(2*N-1)*log2(2*N/(ep1*(2*N-1)));

gradTq = (1-ep1)*(Kq'*(logm(Gepq)/log(2))*Kq - Kq'*(logm(ZGepq)/log(2))*Kq);
gradsyq = 0.5*(gradTq + gradTq');
gradTp = (1-ep1)*(Kp'*(logm(Gepp)/log(2))*Kp - Kp'*(logm(ZGepp)/log(2))*Kp);
gradsyp = 0.5*(gradTp + gradTp');
Trrho = real(trace((gradsyq+gradsyp)*rho0));

cvx_begin quiet
variables y(33) z(33)
expressions T(N,N)
T = y(33)*eye(N);
for i = 1:32
    T = T + y(i)*Con{i,1};
end
maximize(Num2.'*y-ep2*sum(z))
subject to
lambda_min((gradsyq+gradsyp)-T) >=0
y<=z
-z<=y
cvx_end

sigma = Num2.'*y-ep2*sum(z);

D2 = fep + sigma - Trrho - zeta;
end


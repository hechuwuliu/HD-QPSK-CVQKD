function [delrho,Tr,ep2] = S5(rho0,Z0Kq,Z1Kq,Z0Kp,Z1Kp,Nc,ep1,Con,Num)
%
% step 1, primal problem
%
N = 4*(Nc+1);
entg = eye(2*N)/(2*N);
% measuring q
Kq = Z0Kq+Z1Kq;
G1q = (1-ep1)*Kq*rho0*Kq'+ep1*entg;
G2q = (1-ep1)*(Z0Kq*rho0*Z0Kq'+Z1Kq*rho0*Z1Kq')+ep1*entg;
L1q = logm(G1q)/log(2);
L2q = logm(G2q)/log(2);
gradTq = (1-ep1)*(Kq'*L1q*Kq - Kq'*L2q*Kq);
gradsyq = 0.5*(gradTq + gradTq');  %It should be symmetrical, but it didn't, so it was forced to do it.
% measuring p
Kp = Z0Kp+Z1Kp;
G1p = (1-ep1)*Kp*rho0*Kp'+ep1*entg;
G2p = (1-ep1)*(Z0Kp*rho0*Z0Kp'+Z1Kp*rho0*Z1Kp')+ep1*entg;
L1p = logm(G1p)/log(2);
L2p = logm(G2p)/log(2);
gradTp = (1-ep1)*(Kp'*L1p*Kp - Kp'*L2p*Kp);
gradsyp = 0.5*(gradTp + gradTp'); 

cvx_begin quiet
variable Q(N,N) hermitian   %del_rho
minimize(real(trace((gradsyq+gradsyp)*Q)))    % Q.'*gradsy.'
subject to
trace(Q*Con{1,1}) == 0
trace(Q*Con{2,1}) == 0
trace(Q*Con{3,1}) == 0
trace(Q*Con{4,1}) == 0
trace(Q*Con{5,1}) == 0
trace(Q*Con{6,1}) == 0
trace(Q*Con{7,1}) == 0
trace(Q*Con{8,1}) == 0
trace(Q*Con{9,1}) == 0
trace(Q*Con{10,1}) == 0
trace(Q*Con{11,1}) == 0
trace(Q*Con{12,1}) == 0
trace(Q*Con{13,1}) == 0
trace(Q*Con{14,1}) == 0
trace(Q*Con{15,1}) == 0
trace(Q*Con{16,1}) == 0
trace(Q*Con{17,1}) == 0
trace(Q*Con{18,1}) == 0
trace(Q*Con{19,1}) == 0
trace(Q*Con{20,1}) == 0
trace(Q*Con{21,1}) == 0
trace(Q*Con{22,1}) == 0
trace(Q*Con{23,1}) == 0
trace(Q*Con{24,1}) == 0
trace(Q*Con{25,1}) == 0
trace(Q*Con{26,1}) == 0
trace(Q*Con{27,1}) == 0
trace(Q*Con{28,1}) == 0
trace(Q*Con{29,1}) == 0
trace(Q*Con{30,1}) == 0
trace(Q*Con{31,1}) == 0
trace(Q*Con{32,1}) == 0
trace(Q*eye(N)) == 0
lambda_min(Q+rho0) >= 0
cvx_end

if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
    if lambda_min(Q+rho0) < 0
        Q = Q - lambda_min(Q+rho0)*eye(N);
    end
    ep2 = abstractep(rho0,Con,Num);
    delrho = Q;
    Tr = real(trace((gradsyq+gradsyp)*Q));
else
    delrho = zeros(N);
    Tr = 0;
    ep2 = abstractep(rho0,Con,Num);
end
end
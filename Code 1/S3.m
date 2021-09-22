function [rho0,ep2] = S3(Nc, Con, Num)
%
% generate a primal density matrix rho0
%
% minimize(0) is the feasibility problem
% maximize(lambda_min(Q)) is to maximize the minimal eigenvalue of Q
% 
N = 4*(Nc+1);
cvx_precision high
cvx_begin quiet
variable Q(N,N) hermitian semidefinite % Q is the variable that can be used as rho0
maximize(lambda_min(Q)) % or use minimize(0)
subject to
trace(Q*Con{1,1}) == Num(1)
trace(Q*Con{2,1}) == Num(2)
trace(Q*Con{3,1}) == Num(3)
trace(Q*Con{4,1}) == Num(4)
trace(Q*Con{5,1}) == Num(5)
trace(Q*Con{6,1}) == Num(6)
trace(Q*Con{7,1}) == Num(7)
trace(Q*Con{8,1}) == Num(8)
trace(Q*Con{9,1}) == Num(9)
trace(Q*Con{10,1}) == Num(10)
trace(Q*Con{11,1}) == Num(11)
trace(Q*Con{12,1}) == Num(12)
trace(Q*Con{13,1}) == Num(13)
trace(Q*Con{14,1}) == Num(14)
trace(Q*Con{15,1}) == Num(15)
trace(Q*Con{16,1}) == Num(16)
trace(Q*Con{17,1}) == Num(17)
trace(Q*Con{18,1}) == Num(18)
trace(Q*Con{19,1}) == Num(19)
trace(Q*Con{20,1}) == Num(20)
trace(Q*Con{21,1}) == Num(21)
trace(Q*Con{22,1}) == Num(22)
trace(Q*Con{23,1}) == Num(23)
trace(Q*Con{24,1}) == Num(24)
trace(Q*Con{25,1}) == Num(25)
trace(Q*Con{26,1}) == Num(26)
trace(Q*Con{27,1}) == Num(27)
trace(Q*Con{28,1}) == Num(28)
trace(Q*Con{29,1}) == Num(29)
trace(Q*Con{30,1}) == Num(30)
trace(Q*Con{31,1}) == Num(31)
trace(Q*Con{32,1}) == Num(32)
trace(Q*eye(N)) == 1
cvx_end

if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
    if lambda_min(Q) < 0 % this is a method mentioned in [Quantum 2, 77 (2018)]
        rho0 = Q - lambda_min(Q)*eye(N); % we make the matrix positive semidefinite
    else
        rho0 = Q;
    end
    ep2 = abstractep(rho0,Con,Num);
else % no feasible density matrix, therefore no key rate
    rho0 = zeros(N);
    ep2 = 0;
end
end


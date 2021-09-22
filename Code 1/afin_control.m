function [keyrate,Tr,ep2] = afin_control(alpha,Nc,L,p0,xi,ep1,delta_c)

eta = 10^(-0.02*L); % transmission efficiency
% constraints
Con = S1(Nc); % Hermitian matrix
Num = S2(alpha,eta,p0,xi); % expected value

% primal density matrix
[rho0,ep2] = S3(Nc,Con,Num);

% interval operators
% Computing matrices of interval operators is time consuming.
% We advise users to run S4 alone and save the results for further use,
% as long as the parameters Nc,delta_c are not changed.
[Z0Kq,Z1Kq,Z0Kp,Z1Kp] = S4(Nc,delta_c);


% After generating the primal density matrix, we can begin the step 1:
% ->1 find a del_rho that makes the linearization term minimum --> S_5
% ->2 search the optimal lambda(belongs to (0,1)) to get the minimal key
% rate with density matrix rho0+lambda*del_rho --> S_6
% ->3 iteration until the linearization term Tr is small enough
CheckTr = -1;
for j = 1:300 %iteration
    [del_rho,Tr,ep2] = S5(rho0,Z0Kq,Z1Kq,Z0Kp,Z1Kp,Nc,ep1,Con,Num);
    if Tr > -10^-6
        break;
    else
        if Tr > CheckTr  % when iteration, we save the one makes Tr minimum.
            CheckTr = Tr;
            Savedoc = rho0;
            Saveep = ep2;
        end
        rho0 = S6(rho0, del_rho, Z0Kq,Z1Kq,Z0Kp,Z1Kp, ep1, Nc);
    end
end
% If the output of the loop is not the best, we will extract the saved optimal matrix.
if Tr < CheckTr || Tr>0
    Tr = CheckTr;
    rho0 = Savedoc;
    ep2 = Saveep;
end


% Step 2. We do the dual problem by S_7
D2 = S7(rho0,Z0Kq,Z1Kq,Z0Kp,Z1Kp,Con,Num,Nc,ep1,ep2);
% Error correction
[EC,pass] = S8(eta,xi,alpha,delta_c);
keyrate = 0.5*(D2 - 2*pass*EC);
end


function [rhoout] = S6(rho0,Q,Z0Kq,Z1Kq,Z0Kp,Z1Kp,ep,Nc)
%
% step 1 and search the optimal density matrix to find the minimal key rate Dq
% Here we use a method inspired by BianrySearch. The complete introduction is in our preprint
%
N = 4*(Nc+1);
entg = eye(2*N)/(2*N);
Kq = Z0Kq+Z1Kq;
Dq=@(x) real(trace(((1-ep)*Kq*x*Kq'+ep*entg)*logm((1-ep)*Kq*x*Kq'+ep*entg)/log(2))-...
    trace(((1-ep)*Kq*x*Kq'+ep*entg)*logm((1-ep)*(Z0Kq*x*Z0Kq'+Z1Kq*x*Z1Kq')+ep*entg)/log(2)));
Kp = Z0Kp+Z1Kp;
Dp=@(x) real(trace(((1-ep)*Kp*x*Kp'+ep*entg)*logm((1-ep)*Kp*x*Kp'+ep*entg)/log(2))-...
    trace(((1-ep)*Kp*x*Kp'+ep*entg)*logm((1-ep)*(Z0Kp*x*Z0Kp'+Z1Kp*x*Z1Kp')+ep*entg)/log(2)));
A = rho0;
B = A + 0.5*Q;
C = A + Q;
a = Dq(A)+Dp(A);
b = Dq(B)+Dp(B);
c = Dq(C)+Dp(C);
while(1)
    if b > 0.5*(a+c)  % due to the numerical instability, we may find unexpected concave.
        break; %That means the difference between a b c reaches the calculating accuracy of computer, so we break it.
    end
    Q = 0.5*Q;
    if a-b>=10^-10 && b-c>=10^-10  % Usually, the accuracy of calculation is 10^-13, the choice 10^-10 here is safe.
        A = B;
        a = b;
        B = A + 0.5*Q;
        b = Dq(B)+Dp(B);
    elseif b-a>=10^-10 && c-b>=10^-10
        C = B;
        c = b;
        B = A + 0.5*Q;
        b = Dq(B)+Dp(B);
    elseif a-b>=10^-10 && c-b>=10^-10
        te1 = Dq(A + 0.5*Q)+Dp(A + 0.5*Q);
        if b - te1 > 10^-10
            C = B;
            c = b;
            B = A + 0.5*Q;
            b = te1;
        else
            te2 = Dq(C - 0.5*Q)+Dp(C - 0.5*Q);
            if b - te2 > 10^-10
                A = B;
                a = b;
                B = C - 0.5*Q;
                b = te2;
            else
                A = A + 0.5*Q;
                a = te1;
                C = C - 0.5*Q;
                c = te2;
            end
        end
    else
        break;
    end
end
% lambda belongs to£¨0,1£©, so we need to avoid choosing the case that lambda = 0(A) or 1(C)
rhoout = B;

end
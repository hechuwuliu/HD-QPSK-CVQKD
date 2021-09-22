function f = max_gamma_y(N_cutoff,Gamma,grad_f,gamma,eps)
%Search for the maximum of the product of gamma and y
gradsy = 0.5*(grad_f+grad_f');%symmetrization
cvx_begin quiet
    variable y(33);
    variable z(33);
    expressions lhs(4*(N_cutoff+1),4*(N_cutoff+1));
    lhs = zeros(4*(N_cutoff+1));
    for i = 1:33
        lhs = lhs + y(i)*Gamma{1,i};
    end
    maximize (gamma*y-eps*sum(z))
    subject to
         lambda_min(gradsy-lhs)>=0;
         abs(y) <= z; 
cvx_end
f = gamma*y-eps*sum(z);
end


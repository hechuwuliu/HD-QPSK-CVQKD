function [EC, pass] = S8(eta,xi,alpha,delta)
%
% Due to the symmetry between operator q and p in phase space, we only calculate the error correction of q.
% The results of q and p are the same.
%
% measure q
syms q
% states to the right of the p axis in phase space == state 0 (pi/4) and state 3 (7pi/4)
Pq0 = 1/sqrt(pi*(xi*eta+1)) * exp((-(q-sqrt(eta)*alpha)^2)/(xi*eta+1)); %sqrt(eta)*alpha = sqrt(2*eta)*(alpha/sqrt(2))
% states to the left of the p axis in phase space == state 1 (3pi/4) and state 2 (5pi/4)
Pq1 = 1/sqrt(pi*(xi*eta+1)) * exp((-(q+sqrt(eta)*alpha)^2)/(xi*eta+1));

ri_bot = int(Pq0,q,-delta,delta);
ri_0 = int(Pq0,q,delta,inf);
ri_1 = int(Pq0,q,-inf,-delta);

le_bot = int(Pq1,q,-delta,delta);
le_0 = int(Pq1,q,delta,inf);
le_1 = int(Pq1,q,-inf,-delta);

pass = 1 - 0.5*ri_bot - 0.5*le_bot; % actually, due to the symmetry of the states about p axis, ri_bot=le_bot , ri_0 = le_1 , ri_1 = le_0
err = 1/pass * ri_0; % 1-err = 1/pass * ri_1   pass = ri_0+ri_1 = le_0+le_1 = 1-ri_bot = 1-le_bot
her = -err*log2(err) - (1-err)*log2(1-err);

Pz0 = 0.5*(le_0+ri_0)/pass; % the probability of measurement result is 0
Pz1 = 0.5*(ri_1+le_1)/pass; % the probability of measurement result is 1
HZ = -Pz0*log2(Pz0)-Pz1*log2(Pz1);
EC = 0.05*HZ + 0.95*her;
pass = double(pass);
EC = double(EC);
end


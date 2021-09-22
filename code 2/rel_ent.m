function f = rel_ent(sigma,tau)
%this function calculates relative entropy
f = trace(sigma*logm(sigma)/log(2))-trace(sigma*logm(tau)/log(2));
end
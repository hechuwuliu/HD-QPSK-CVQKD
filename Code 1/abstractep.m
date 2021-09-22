function [output] = abstractep(A,Con,Num)
%
% abstract the parameter ep2 of imprecision of numerical method used in dual problem
%
B = zeros(33,1);
for i = 1:32
    B(i) = trace(A*Con{i,1}) - Num(i);
end
B(33) = trace(A)-1;
output = max(abs(B));
end


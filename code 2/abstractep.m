function [output] = abstractep(A,Con,Num)
%abstract imprecision ep
B = zeros(33,1);
for i = 1:33
    B(i) = trace(A*Con{i}.') - Num(i);
end
output = max(abs(B));
end


function [Num] = S2(a, eta, p0, xi)
%%% The values of constraints
%
% a is the amplitude of the optical pulses / coherent states
% eta is the transmission efficiency
% p0 = 0.25 is the probability of sending each kind of coherent states
% xi is the excess noise
%
p1 = 0.5-p0;
p2 = p0;
p3 = p1;

Num = zeros(32,1);

% The values to the right of the constraint equations
Num(1) = p0*sqrt(eta)*a;
Num(2) = -p1*sqrt(eta)*a;
Num(3) = -p2*sqrt(eta)*a;
Num(4) = p3*sqrt(eta)*a;
Num(5) = p0*sqrt(eta)*a;
Num(6) = p1*sqrt(eta)*a;
Num(7) = -p2*sqrt(eta)*a;
Num(8) = -p3*sqrt(eta)*a;
Num(9) = 0;
Num(10) = 0;
Num(11) = 0;
Num(12) = 0;
Num(13) = p0*(eta*a^2+eta*xi/2);
Num(14) = p1*(eta*a^2+eta*xi/2);
Num(15) = p2*(eta*a^2+eta*xi/2);
Num(16) = p3*(eta*a^2+eta*xi/2);
Num(17) = 2*sqrt(p0*p1)*exp(-a^2)*cos(a^2);
Num(18) = 2*sqrt(p0*p1)*exp(-a^2)*sin(a^2);
Num(19) = 2*sqrt(p0*p2)*exp(-2*a^2);
Num(20) = 0;
Num(21) = 2*sqrt(p0*p3)*exp(-a^2)*cos(a^2);
Num(22) = -2*sqrt(p0*p3)*exp(-a^2)*sin(a^2);
Num(23) = 2*sqrt(p1*p2)*exp(-a^2)*cos(a^2);
Num(24) = 2*sqrt(p1*p2)*exp(-a^2)*sin(a^2);
Num(25) = 2*sqrt(p1*p3)*exp(-2*a^2);
Num(26) = 0;
Num(27) = 2*sqrt(p2*p3)*exp(-a^2)*cos(a^2);
Num(28) = 2*sqrt(p2*p3)*exp(-a^2)*sin(a^2);
Num(29) = p0;
Num(30) = p1;
Num(31) = p2;
Num(32) = p3;
end


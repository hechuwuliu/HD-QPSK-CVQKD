function rho0 = search_rho0(N_cutoff,Gamma,gamma)
%this function searches rho0 
N = 4*(N_cutoff+1);
%search for the start point of rho0 using cvx package
cvx_precision high
cvx_begin quiet
    variable rho(N,N)  complex semidefinite
    minimize (0)
    rho == hermitian_semidefinite(N);
    subject to
        trace(rho) == 1;
        trace(rho*((Gamma{2}).')) == gamma(2);
        trace(rho*((Gamma{3}).')) == gamma(3);
        trace(rho*((Gamma{4}).')) == gamma(4);
        trace(rho*((Gamma{5}).')) == gamma(5);
        trace(rho*((Gamma{6}).')) == gamma(6);
        trace(rho*((Gamma{7}).')) == gamma(7);
        trace(rho*((Gamma{8}).')) == gamma(8);
        trace(rho*((Gamma{9}).')) == gamma(9);
        trace(rho*((Gamma{10}).')) == gamma(10);
        trace(rho*((Gamma{11}).')) == gamma(11);
        trace(rho*((Gamma{12}).')) == gamma(12);
        trace(rho*((Gamma{13}).')) == gamma(13);
        trace(rho*((Gamma{14}).')) == gamma(14);
        trace(rho*((Gamma{15}).')) == gamma(15);
        trace(rho*((Gamma{16}).')) == gamma(16);
        trace(rho*((Gamma{17}).')) == gamma(17);
        trace(rho*((Gamma{18}).')) == gamma(18);
        trace(rho*((Gamma{19}).')) == gamma(19);
        trace(rho*((Gamma{20}).')) == gamma(20);
        trace(rho*((Gamma{21}).')) == gamma(21);
        trace(rho*((Gamma{22}).')) == gamma(22);
        trace(rho*((Gamma{23}).')) == gamma(23);
        trace(rho*((Gamma{24}).')) == gamma(24);
        trace(rho*((Gamma{25}).')) == gamma(25);
        trace(rho*((Gamma{26}).')) == gamma(26);
        trace(rho*((Gamma{27}).')) == gamma(27);
        trace(rho*((Gamma{28}).')) == gamma(28);
        trace(rho*((Gamma{29}).')) == gamma(29);
        trace(rho*((Gamma{30}).')) == gamma(30);
        trace(rho*((Gamma{31}).')) == gamma(31);
        trace(rho*((Gamma{32}).')) == gamma(32);
        trace(rho*((Gamma{33}).')) == gamma(33);
cvx_end

if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
    rho0 = rho;
else
    rho0 = zeros(N);
end
end
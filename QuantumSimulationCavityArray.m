function [OP ,P1m]= QuantumSimulationCavityArray(wc,D,g,Nphoton ,L,A,Sp,N_ex ,Hhopp ,t, Model)
HJC = zeros(( Nphoton +1)^L*2^L,( Nphoton +1)^L*2^L);
for i=1:L
    switch Model
        case 'Jaynes -Cummings '
            HJC = HJC + wc*A{i}'*A{i} + (D+wc)*Sp{i}*Sp{i}' + g*(Sp{i}*A{i}+Sp{i}'*A{i }');
        case 'Rabi'
            HJC = HJC + wc*A{i}'*A{i} + (D+wc)*Sp{i}*Sp{i}' + g*(Sp{i}+Sp{i}')*(A{i}+A{ i}');
    end
end
H = HJC + Hhopp;
dt = t(2)-t(1);
U = expm(-1i*H*dt);
up = [1 0]';
down = [0 1]';
Fock = eye(Nphoton +1);
theta1 = 0.5* atan (2*g*sqrt (1)/D);
phi_1m = cos(theta1)*kron(down ,Fock (:,2))... 
-sin(theta1)*kron(up ,Fock (:,1));
PSI_0 = phi_1m;
for k=1:L-1
    PSI_0 = kron(PSI_0 ,phi_1m);
end
dn_T = zeros(size(t));
P1m = zeros(size(t));
for n=1: length(t)
    if n==1
        PSI = PSI_0;
    else
        PSI = U*PSI;
    end
    dn = 0;
    for i=1:L
        dn = dn+PSI'*N_ex{i}^2*PSI-(PSI'*N_ex{i}*PSI)^2;
    end
    P1m(n) = abs(PSI_0'*PSI)^2;
    dn_T(n) = dn;
end
OP = mean(dn_T);
end



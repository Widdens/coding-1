Nsim = 25;
 L = 2;
 wc = 1;
 g = 1e-2*wc;
 J = 1e-4*wc;
 Nph = 2;
 dimFock = Nph+1;
dimT = 2* dimFock;
Deltai = 10^( -2)*g;
Deltaf = 10^(+2)*g;
xi = log10(Deltai/g);
xf = log10(Deltaf/g);
dx = (xf -xi)/(Nsim -1);
x = xi:dx:xf;
OP_JC = zeros(size(x));
OP_R = zeros(size(x));
 A = cell(1,L);
 Sp = cell(1,L);
 N_ex = cell(1,L);
 Iatom = eye(2);
 Icav = eye(dimFock );
 Is = eye(2* dimFock);
 for i=1:L 
A{i} = acav(i,L,Nph ,Is,Iatom);
Sp{i} = sigmap(i,L,Is,Icav);
N_ex{i} = A{i}'*A{i}+Sp{i}*Sp{i}';
 end
Ad = ones(L);
Ad = triu(Ad)-eye(L); 
Hhopp = zeros(dimT^L,dimT^L); 
for  i=1:L
    for j=1:L
       Hhopp = Hhopp - J*Ad(i,j)*A{i}'*A{j} - J*Ad(i,j)*A{i}*A{j}';
    end
end
Nt = 10000;
ti = 0.01/J;
tf = 1/J;
dt = (tf -ti)/(Nt -1);
 t = ti:dt:tf;
 Lambda_R = zeros(Nsim ,Nt); 
 Lambda_JC = zeros(Nsim ,Nt); 
 parfor n=1: Nsim 44
 D = g*10^(x(n));
 Model = 'Rabi';
[OP_R(n), P1m_R] = QuantumSimulationCavityArray(wc,D,g,Nph,L,A,Sp ,N_ex,Hhopp,t,Model);
Model = 'Jaynes -Cummings ';
[OP_JC(n),P1m_JC] = QuantumSimulationCavityArray(wc,D,g,Nph ,L,A,Sp ,N_ex ,Hhopp,t,Model);
Lambda_R(n,:) = -1/L*log2(P1m_R); % Rate function for the Rabi model
Lambda_JC(n,:) = -1/L*log2(P1m_JC); % Rate function for the Jaynes -Cummings model
 end 
 figure () 
 box on 
 hold on 
 plot(J*t,Lambda_R(end ,:),'b-','Linewidth',3)
 plot(J*t,Lambda_JC(end ,:),'r-','Linewidth',3) 
 hold off 
 xlabel('$Jt$','Interpreter','LaTex','Fontsize', 30) 
 ylabel('$\Lambda(t)$','Interpreter','LaTex','Fontsize', 30) 
 set(gca ,'fontsize' ,21) 
 legend ({'$\mbox{RH}$','$\mbox{JCH}$'},'Interpreter','latex','Fontsize', 21,'Location',' best')

 figure () 
 box on 
 hold on 
 plot(x,real(OP_R),'.b','Markersize' ,30) 
 plot(x,real(OP_JC),'.r','Markersize' ,30) 
 hold off 
 xlabel('$\mbox{Log}_{10}(\Delta/g)$','Interpreter','LaTex','Fontsize', 30) 
 ylabel('$\mbox{OP}$','Interpreter','LaTex','Fontsize', 30) 
 set(gca ,'fontsize' ,21) 
 legend ({'$\mbox{RH}$','$\mbox{JCH}$'},'Interpreter','latex','Fontsize', 21,'Location',' best')

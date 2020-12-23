
clear all
Nspins=6;
H0=zeros(2^Nspins,2^Nspins);
H1=zeros(2^Nspins,2^Nspins);
J=0;
alpha=0.2;
for j=1:Nspins
    for i=j+1:Nspins
        Jij=abs(i-j)^(-alpha);
        J=J+(Nspins-1)^(-1)*Jij;
    end
end
 B=J/0.42;
 Sx=[0 1;1 0];
 Sz=[1 0;0 -1];
 for i=1:Nspins
     Szi=getSci(Sz,i,Nspins);
     H0=H0-B*Szi;
     for j=1:Nspins
         if i~=j
             Sxi=getSci(Sx,i,Nspins);
             Sxj=getSci(Sx,j,Nspins);
             Vij=abs(i-j)^(-alpha)/J;
             H1=H1-Vij*Sxi*Sxj;
         end
     end
 end
 H=H0+H1;
 xr=[1 1]'/sqrt(2);
 x1=[-1 1]'/sqrt(2);
 Xr=xr;
 X1=x1;
 for n=1:Nspins-1
     Xr=kron(Xr,xr);
     X1=kron(X1,x1);
 end
 PSI_0=Xr;
 ti=0;
 tf=22;
 Nt=10000;
 dt=(tf-ti)/(Nt-1);
 t=ti:dt:tf;
 U=expm(-1i*H*dt);
 Mx=zeros(size(t));
 Lambda=zeros(size(t));
 SSx=0;
 for i=1:Nspins
     Sxi=getSci(Sx,i,Nspins);
     SSx=SSx+Sxi/Nspins;
 end
 for n=1:length(t)
     if n==1
         PSI=PSI_0;
     else
         PSI=U*PSI;
     end
     Pr=abs(Xr'*PSI)^2;
     P1=abs(X1'*PSI)^2;
     Lambda(n)=min(-Nspins^(-1)*log(Pr),-Nspins^(-1)*log(P1));
     Mx(n)=PSI'*SSx*PSI;
 end
 
 figure()
 plot(B*t,Lambda,'b-','Linewidth',3)
 xlabel('$B t$','Interpreter','LaTex','Fontsize',30)
 set(gca,'fontsize',21)
 xlim([0 5])
 
 figure()
 box on
 plot(t*B,real(Mx),'r-','Linewidth',3)
 xlabel('$B t$','Interpreter','LaTex','Fontsize',30)
 ylabel('$\Lambda(t)$','Intetpreter','LaTex','Fontsize',30)
 set(gca,'fontsize',21)
 xlim([0 5]) 

  figure()
  box on 
  plot(t*B,real(Mx),'r-','Linewidth',2)z
  xlabel('$B t$','Interpreter','LaTex','Fontsize',30)
  ylabel('$\langle M_x\rangle','Interpreter','LaTex','Fontsize',30)
  set(gca,'fontsize',21)
  xlim([0 100])
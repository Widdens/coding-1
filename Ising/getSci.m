function Sci=getSci(sc,i,Nspins)
Is=eye(2);
Op_total=zeros(2,2,Nspins);%�����������cell�ɴ�ͬ��Ч��
for site=1:Nspins
    Op_total(:,:,site)=Is+double(eq(i,site))*(sc-Is);
end
Sci=Op_total(:,:,1);
for site=2:Nspins
    Sci=kron(Sci,Op_total(:,:,site));
end
end
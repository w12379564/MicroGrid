clear;

T=24;
baseMVA=100;
%微网参数
dsmax_base=[0.36]/baseMVA;
Loadcoe=[0.6,0.61,0.62,0.63,0.64,0.65,0.7,0.8,0.9,1.1,1.5,1.6,1.2,1.1,1,0.9,0.8,0.85,1,1.5,1.55,1.1,0.8,0.7];
%Loadcoe=[1];
nds=1;
dsmax=zeros(nds,T);
for t=1:T
    dsmax(:,t)=dsmax_base*Loadcoe(t);
end

%场景数
Wcoe=[1.4,1.4,1.3,1.2,1.1,1,0.8,0.7,0.6,0.55,0.5,0.45,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.35,1.4];
%Wcoe=[1];

K=5;
Wact_base=[0.1448,0.1327,0.18,0.1673,0.1852]/baseMVA;
phi=[0.1;0.2;0.4;0.2;0.1];

Wact=zeros(T,K);
for t=1:T
    Wact(t,:)=Wact_base*Wcoe(t);
end


PDGmax=[0.4]/baseMVA;
Rdgumax=0.1*PDGmax;
Rdgdmax=Rdgumax;
Wmax=0.3/baseMVA;
Us=[26];
ODG=[13];

%储能
PESmax=[0.05]/baseMVA;
CESmax=[0.2]/baseMVA;
CES0=[0]/baseMVA;


%variables
%UL problem
Bds=sdpvar(nds,1,T,'full');
PDG=sdpvar(nds,1,T,'full');
ds=sdpvar(nds,1,T,'full');
Pes=sdpvar(nds,1,T,'full');

rdsk=sdpvar(nds,K,T,'full');
rdgk=sdpvar(nds,K,T,'full');
w=sdpvar(nds,1,T,'full');
wpk=sdpvar(nds,K,T,'full');

dsn=0.002*ones(1,1,T);
lamda=[10.46,10.46,10.46,10.46,10.46,10.46,10.46,13.58,13.58,13.58,18.57,18.57,13.58,13.58,13.58,13.58,13.58,13.58,13.58,18.57,18.57,13.58,13.58,10.46];
Constraints=[];
for t=1:T
%UL Constraints

Constraints=[Constraints,0<=PDG(:,:,t)<=PDGmax];
Constraints=[Constraints,0<=ds(:,:,t)<=dsmax(:,t)];
Constraints=[Constraints,0<=w(:,:,t)<=Wmax];

%储能约束
Constraints=[Constraints,-PESmax<=Pes(:,1,t)<=PESmax];
CES=CES0;
for tt=1:t
    CES=CES-Pes(:,1,tt);
end
Constraints=[Constraints,0<=CES<=CESmax];

%微网功率平衡
Constraints=[Constraints,(dsn(:,:,t)+PDG(:,:,t)+w(:,:,t)+ Pes(:,:,t)-ds(:,:,t))==0];
for k=1:K
Constraints=[Constraints,0<=wpk(:,k,t)<=Wact(t,k)];
Constraints=[Constraints,(rdgk(:,k,t)+rdsk(:,k,t)+(Wact(t,k)-w(:,:,t)-wpk(:,k,t)))==0];
Constraints=[Constraints,-(dsmax(:,t)-ds(:,:,t))<=rdsk(:,k,t)<=ds(:,:,t)];
Constraints=[Constraints,-PDG(:,:,t)<=rdgk(:,k,t)<=(PDGmax-PDG(:,:,t))];
Constraints=[Constraints,-Rdgdmax<=rdgk(:,k,t)<=Rdgumax];
end

end

Objective=0;

for t=1:T
    Objective=Objective -(Us*ds(:,:,t)-lamda(t)*dsn(:,:,t)-ODG*PDG(:,:,t)-Us*rdsk(:,:,t)*phi - ODG*rdgk(:,:,t)*phi);
end

ops = sdpsettings('solver','gurobi','verbose',2);
%ops.gurobi.MIPGap=0.01;
optimize(Constraints,Objective,ops);

ds_data=reshape(double(ds),nds,T);
w_data=reshape(double(w),nds,T);
PDG_data=reshape(double(PDG),nds,T);
Pes_data=reshape(double(Pes),nds,T);

rdsk_e=zeros(1,T);
for t=1:T
    rdsk_e(t)=double(rdsk(:,:,t)*phi);
end

w_e=zeros(1,T);
for t=1:T
    w_e(t)=double((Wact(t,:)-wpk(:,:,t))*phi);
end

wpk_data=zeros(1,T);
for t=1:T
    wpk_data(t)=double(wpk(:,:,t)*phi);
end

rdgk_e=zeros(1,T);
for t=1:T
    rdgk_e(t)=double(rdgk(:,:,t)*phi);
end

profit=zeros(1,T);
for t=1:T
    profit(t) = double(Us*ds(:,:,t)-lamda(t)*dsn(:,:,t)-ODG*PDG(:,:,t)-Us*rdsk(:,:,t)*phi - ODG*rdgk(:,:,t)*phi);
end
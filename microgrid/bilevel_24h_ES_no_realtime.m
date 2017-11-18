clear;

% define named indices into bus, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = indexBus;

[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = indexBrch;

[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = indexGen;

%network
data=loadcase('case24_ieee_rts');

baseMVA=data.baseMVA;
bus=data.bus;
branch=data.branch;
gen=data.gen;

nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of lines

[ Bbus, Bf,Bft, Pbusinj, Pfinj ] = calBdc( baseMVA, bus, branch );

% generator info
% on = find(gen(:, GEN_STATUS) > 0);      %% which generators are on?
%此处修改发电机数量
% on = [1;2;3;4];
% gbus = gen(on, GEN_BUS);                %% what buses are they at?
gbus=[1;1;2;2;7;13;15;18;21;22];

% form net complex bus power injection vector
nb = size(bus, 1);
ngon = size(gbus, 1);
Cg = sparse(gbus, (1:ngon)', ones(ngon, 1), nb, ngon);  %% connection matrix of generation

dsbus=[18];
nds=size(dsbus,1);
Cds=sparse(dsbus, (1:nds)', ones(nds, 1), nb, nds);  %% connection matrix of ds

dcbus=[1;2;8;14;19;20];
ndc=size(dcbus,1);
Cdc=sparse(dcbus, (1:ndc)', ones(ndc, 1), nb, ndc);  %% connection matrix of ds

T=24;

%微网参数
dsmax_base=[0.36]/baseMVA;
Loadcoe=[0.6,0.61,0.62,0.63,0.64,0.65,0.7,0.8,0.9,1.1,1.5,1.6,1.2,1.1,1,0.9,0.8,0.85,1,1.5,1.55,1.1,0.8,0.7];
%Loadcoe=[1];

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



%联络线
PTLmax=0.2/baseMVA;


Pgmax=[40;152;40;152;300;591;60;400;400;300]/baseMVA;
dc_base=[110.16;98.94;174.42;197.88;184.62;130.56]/baseMVA;

dc=zeros(ndc,T);

for t=1:T
    dc(:,t)=dc_base*Loadcoe(t);
end

BaseG=[1;2;3;4;5;6;8;9];
PearkG=[7;10];

Rgumax=zeros(ngon,1);
Rgumax(BaseG)=0*Pgmax(BaseG);
Rgumax(PearkG)=0.2*Pgmax(PearkG);
Rgdmax=Rgumax;

coe=1;
Rdsumax=0/baseMVA;
Rdsdmax=Rdsumax;


%cost
Og=[21;18.57;21;13.58;22;18.57;24.57;13.58;10.46;10.46];

Uc=[23.5;20;25;21;22.3;22];

Vc=1000*ones(ndc,1);

pi_=pi*ones(nb,1);

Bdsmax=30;


M=1e3;

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

%LL Problem
%variables
Pg=sdpvar(ngon,1,T,'full');
rgk=sdpvar(ngon,K,T,'full');

dsn=sdpvar(nds,1,T,'full');
dsnk=sdpvar(nds,K,T,'full');

theta=sdpvar(nb,1,T,'full');
thetak=sdpvar(nb,K,T,'full');

%dual variables
%h1
lamda=sdpvar(nb,1,T,'full');
%h2
lamdak=sdpvar(nb,K,T,'full');

%g1
Ug_down=sdpvar(ngon,1,T,'full');
%g2
Ug_up=sdpvar(ngon,1,T,'full');

%g3
Urgk1_down=sdpvar(ngon,K,T,'full');
%g4
Urgk1_up=sdpvar(ngon,K,T,'full');

%g5
Urgk2_down=sdpvar(ngon,K,T,'full');
%g6
Urgk2_up=sdpvar(ngon,K,T,'full');

%g7
Udsn_down=sdpvar(nds,1,T,'full');
%g8
Udsn_up=sdpvar(nds,1,T,'full');

%g9
Udsnk1_down=sdpvar(nds,K,T,'full');
%g10
Udsnk1_up=sdpvar(nds,K,T,'full');

%g11
Udsnk2_down=sdpvar(nds,K,T,'full');
%g12
Udsnk2_up=sdpvar(nds,K,T,'full');


%g29
Utheta_down=sdpvar(nb,1,T,'full');
%g30
Utheta_up=sdpvar(nb,1,T,'full');
%g31
Uthetak_down=sdpvar(nb,K,T,'full');
%g32
Uthetak_up=sdpvar(nb,K,T,'full');


%bin variables
%g1
binUg_down=binvar(ngon,1,T,'full');
%g2
binUg_up=binvar(ngon,1,T,'full');

%g3
binUrgk1_down=binvar(ngon,K,T,'full');
%g4
binUrgk1_up=binvar(ngon,K,T,'full');

%g5
binUrgk2_down=binvar(ngon,K,T,'full');
%g6
binUrgk2_up=binvar(ngon,K,T,'full');

%g7
binUdsn_down=binvar(nds,1,T,'full');
%g8
binUdsn_up=binvar(nds,1,T,'full');

%g9
binUdsnk1_down=binvar(nds,K,T,'full');
%g10
binUdsnk1_up=binvar(nds,K,T,'full');

%g11
binUdsnk2_down=binvar(nds,K,T,'full');
%g12
binUdsnk2_up=binvar(nds,K,T,'full');


%g29
binUtheta_down=binvar(nb,1,T,'full');
%g30
binUtheta_up=binvar(nb,1,T,'full');
%g31
binUthetak_down=binvar(nb,K,T,'full');
%g32
binUthetak_up=binvar(nb,K,T,'full');

Constraints=[];
for t=1:T
%UL Constraints
Constraints=[Constraints,0<=Bds(:,:,t)<=Bdsmax];
% Constraints=[Constraints,Bds(:,:,t)==Us];
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
Constraints=[Constraints,(rdgk(:,k,t)+rdsk(:,k,t)-dsnk(:,k,t)+(Wact(t,k)-w(:,:,t)-wpk(:,k,t)))==0];
Constraints=[Constraints,-(dsmax(:,t)-ds(:,:,t))<=rdsk(:,k,t)<=ds(:,:,t)];
Constraints=[Constraints,-PDG(:,:,t)<=rdgk(:,k,t)<=(PDGmax-PDG(:,:,t))];
Constraints=[Constraints,-Rdgdmax<=rdgk(:,k,t)<=Rdgumax];
end

%LL Constraints
%KKT Pg
Constraints=[Constraints,(Og - Cg'*lamda(:,:,t) - Ug_down(:,:,t) + Ug_up(:,:,t) - sum(Urgk1_down(:,:,t),2) +  sum(Urgk1_up(:,:,t),2)) == 0];
%KKT rgk
for k=1:K
    Constraints=[Constraints,(phi(k)*Og - Cg'*lamdak(:,k,t) - Urgk1_down(:,k,t) + Urgk1_up(:,k,t) - Urgk2_down(:,k,t) + Urgk2_up(:,k,t)) == 0];
end

%KKT dsn
Constraints=[Constraints,(-Bds(:,:,t) + Cds'*lamda(:,:,t) - Udsn_down(:,:,t) + Udsn_up(:,:,t) + sum(Udsnk1_down(:,:,t),2) - sum(Udsnk1_up(:,:,t),2)) == 0];
%KKT dsnk
for k=1:K
    Constraints=[Constraints,(phi(k)*Bds(:,:,t) - Cds'*lamdak(:,k,t) - Udsnk1_down(:,k,t) + Udsnk1_up(:,k,t) -Udsnk2_down(:,k,t) + Udsnk2_up(:,k,t)) == 0];
end


%KKT theta
Constraints=[Constraints,(Bbus'*lamda(:,:,t) - sum(Bbus'*lamdak(:,:,t),2) - Utheta_down(:,:,t) + Utheta_up(:,:,t))==0];
%KKT thetak
for k=1:K
    Constraints=[Constraints,(Bbus'*lamdak(:,k,t) - Uthetak_down(:,k,t) + Uthetak_up(:,k,t)) == 0];
end

%h1
Constraints=[Constraints,(Cg*Pg(:,:,t) - Cds*dsn(:,:,t) - Cdc*dc(:,t) - Bbus*theta(:,:,t)) ==0];

for k=1:K
    %h2
    Constraints=[Constraints,(Cg*rgk(:,k,t) + Cds*dsnk(:,k,t) - Bbus*(thetak(:,k,t)-theta(:,:,t)))==0];
end
%g1
Constraints=[Constraints,Pg(:,:,t)>=0, Ug_down(:,:,t)>=0, Pg(:,:,t)<=binUg_down(:,:,t)*M , Ug_down(:,:,t)<=(1-binUg_down(:,:,t))*M ];
%g2
Constraints=[Constraints,(Pgmax-Pg(:,:,t))>=0, Ug_up(:,:,t)>=0, (Pgmax-Pg(:,:,t))<=binUg_up(:,:,t)*M , Ug_up(:,:,t)<=(1-binUg_up(:,:,t))*M ];

for k=1:K
    %g3
    Constraints=[Constraints,(rgk(:,k,t)+Pg(:,:,t))>=0,Urgk1_down(:,k,t)>=0,(rgk(:,k,t)+Pg(:,:,t))<=binUrgk1_down(:,k,t)*M ,Urgk1_down(:,k,t)<=(1-binUrgk1_down(:,k,t))*M ];
    %g4
    Constraints=[Constraints,(Pgmax-Pg(:,:,t)-rgk(:,k,t))>=0,Urgk1_up(:,k,t)>=0,(Pgmax-Pg(:,:,t)-rgk(:,k,t))<=binUrgk1_up(:,k,t)*M ,Urgk1_up(:,k,t)<=(1-binUrgk1_up(:,k,t))*M ];
    %g5
    Constraints=[Constraints,(rgk(:,k,t)+Rgdmax)>=0,Urgk2_down(:,k,t)>=0,(rgk(:,k,t)+Rgdmax)<=binUrgk2_down(:,k,t)*M ,Urgk2_down(:,k,t)<=(1-binUrgk2_down(:,k,t))*M ];
    %g6
    Constraints=[Constraints,(Rgumax-rgk(:,k,t))>=0,Urgk2_up(:,k,t)>=0,(Rgumax-rgk(:,k,t))<=binUrgk2_up(:,k,t)*M ,Urgk2_up(:,k,t)<=(1-binUrgk2_up(:,k,t))*M ];
end

%g7
Constraints=[Constraints,(dsn(:,:,t)+PTLmax)>=0,Udsn_down(:,:,t)>=0,(dsn(:,:,t)+PTLmax)<=binUdsn_down(:,:,t)*M ,Udsn_down(:,:,t)<=(1-binUdsn_down(:,:,t))*M ];
%g8
Constraints=[Constraints,(PTLmax-dsn(:,:,t))>=0,Udsn_up(:,:,t)>=0,(PTLmax-dsn(:,:,t))<=binUdsn_up(:,:,t)*M ,Udsn_up(:,:,t)<=(1-binUdsn_up(:,:,t))*M ];

for k=1:K
    %g9
    Constraints=[Constraints,(dsnk(:,k,t)+PTLmax-dsn(:,:,t))>=0,Udsnk1_down(:,k,t)>=0,(dsnk(:,k,t)+PTLmax-dsn(:,:,t))<=binUdsnk1_down(:,k,t)*M ,Udsnk1_down(:,k,t)<=(1-binUdsnk1_down(:,k,t))*M ];
    %g10
    Constraints=[Constraints,(dsn(:,:,t)+PTLmax-dsnk(:,k,t))>=0,Udsnk1_up(:,k,t)>=0,(dsn(:,:,t)+PTLmax-dsnk(:,k,t))<=binUdsnk1_up(:,k,t)*M ,Udsnk1_up(:,k,t)<=(1-binUdsnk1_up(:,k,t))*M ];
    %g11
    Constraints=[Constraints,(dsnk(:,k,t)+Rdsdmax)>=0,Udsnk2_down(:,k,t)>=0,(dsnk(:,k,t)+Rdsdmax)<=binUdsnk2_down(:,k,t)*M ,Udsnk2_down(:,k,t)<=(1-binUdsnk2_down(:,k,t))*M ];
    %g12
    Constraints=[Constraints,(Rdsumax-dsnk(:,k,t))>=0,Udsnk2_up(:,k,t)>=0,(Rdsumax-dsnk(:,k,t))<=binUdsnk2_up(:,k,t)*M ,Udsnk2_up(:,k,t)<=(1-binUdsnk2_up(:,k,t))*M ];
end

%g29
Constraints=[Constraints,(theta(:,:,t)+pi)>=0,Utheta_down(:,:,t)>=0,(theta(:,:,t)+pi)<=binUtheta_down(:,:,t)*M ,Utheta_down(:,:,t)<=(1-binUtheta_down(:,:,t))*M ];
%g30
Constraints=[Constraints,(pi-theta(:,:,t))>=0,Utheta_up(:,:,t)>=0,(pi-theta(:,:,t))<=binUtheta_up(:,:,t)*M ,Utheta_up(:,:,t)<=(1-binUtheta_up(:,:,t))*M ];
for k=1:K
    %g31
    Constraints=[Constraints,(thetak(:,k,t)+pi)>=0,Uthetak_down(:,k,t)>=0,(thetak(:,k,t)+pi)<=binUthetak_down(:,k,t)*M ,Uthetak_down(:,k,t)<=(1-binUthetak_down(:,k,t))*M ];
    %g32
    Constraints=[Constraints,(pi-thetak(:,k,t))>=0,Uthetak_up(:,k,t)>=0,(pi-thetak(:,k,t))<=binUthetak_up(:,k,t)*M ,Uthetak_up(:,k,t)<=(1-binUthetak_up(:,k,t))*M ];
end


%ref bus
Constraints=[Constraints,theta(1,:,t)==0];
for k=1:K
    Constraints=[Constraints,thetak(1,k,t)==0];
end

end


Objective=0;

for t=1:T
    Objective = Objective - 100*(lamda(:,:,t)'*Cdc*dc(:,t)...
        -Ug_up(:,:,t)'*Pgmax - sum(Urgk1_up(:,:,t)'*Pgmax) - sum(Urgk2_down(:,:,t)'*Rgdmax) - sum(Urgk2_up(:,:,t)'*Rgumax)...
        -Utheta_down(:,:,t)'*pi_ - Utheta_up(:,:,t)'*pi_ - sum(Uthetak_down(:,:,t)'*pi_) - sum(Uthetak_up(:,:,t)'*pi_)...
        - Og'*Pg(:,:,t) -( Og'*rgk(:,:,t))*phi...
        + Us*ds(:,:,t) - ODG*PDG(:,:,t) - (Us*rdsk(:,:,t) + ODG*rdgk(:,:,t))*phi);
end

ops = sdpsettings('solver','gurobi','verbose',2);
%ops.gurobi.MIPGap=0.01;
optimize(Constraints,Objective,ops);

profit=zeros(1,T);

for t=1:T
    profit(t) = double(Us*ds(:,:,t)-lamda(:,:,t)'*Cds*dsn(:,:,t)-ODG*PDG(:,:,t)-Us*rdsk(:,:,t)*phi+sum(diag(lamdak(:,:,t)'*Cds*dsnk(:,:,t))) - ODG*rdgk(:,:,t)*phi);
end

% t=19;
% p19=double(lamda(:,:,t)'*Cdc*dc(:,t)...
%         -Ug_up(:,:,t)'*Pgmax - sum(Urgk1_up(:,:,t)'*Pgmax) - sum(Urgk2_down(:,:,t)'*Rgdmax) - sum(Urgk2_up(:,:,t)'*Rgumax)...
%         -Utheta_down(:,:,t)'*pi_ - Utheta_up(:,:,t)'*pi_ - sum(Uthetak_down(:,:,t)'*pi_) - sum(Uthetak_up(:,:,t)'*pi_)...
%         - Og'*Pg(:,:,t) -( Og'*rgk(:,:,t))*phi...
%         + Us*ds(:,:,t) - ODG*PDG(:,:,t) - (Us*rdsk(:,:,t) + ODG*rdgk(:,:,t))*phi)

Bds_data=reshape(double(Bds),nds,T);
priceDA=reshape(double(lamda),nb,T);
Pg_data=reshape(double(Pg),ngon,T);
dsn_data=reshape(double(dsn),nds,T);
ds_data=reshape(double(ds),nds,T);
w_data=reshape(double(w),nds,T);
PDG_data=reshape(double(PDG),nds,T);
Pes_data=reshape(double(Pes),nds,T);
dsnk_data=double(dsnk);

rdsk_e=zeros(1,T);
for t=1:T
    rdsk_e(t)=double(rdsk(:,:,t)*phi);
end

dsnk_e=zeros(1,T);
for t=1:T
    dsnk_e(t)=double(dsnk(:,:,t)*phi);
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

Wact_e=zeros(1,T);
for t=1:T
    Wact_e(t)=double(Wact(t,:)*phi);
end

W_delta=zeros(T,K);
for t=1:T
W_delta(t,:)=double(Wact(t,:)-w(:,:,t));
end

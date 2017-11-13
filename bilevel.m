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

nw=1;
wbus=3;
Cw=sparse(wbus, (1:nw)', ones(nw, 1), nb, nw);  %% connection matrix of wind farm

T=24;


%上限
dsmax_base=[808.86]/baseMVA;
dc_base=[110.16;98.94;174.42;197.88;184.62;130.56]/baseMVA;

Loadcoe=[0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.8,0.9,1.1,1.3,1.5,1.2,1.1,1,0.9,0.8,0.85,1,1.1,1.2,1.1,0.8,0.7];

dsmax=zeros(nds,T);
dc=zeros(ndc,T);

for t=1:T
    dsmax(:,t)=dsmax_base*Loadcoe(t);
    dc(:,t)=dc_base*Loadcoe(t);
end

Pgmax=[40;152;40;152;300;591;60;400;400;300]/baseMVA;
PDGmax=[700]/baseMVA;
Wmax=550/baseMVA;
%储能
PESmax=[25]/baseMVA;
CESmax=[100]/baseMVA;
CES0=[50]/baseMVA;

BaseG=[1;2;3;4;5;6;8;9];
PearkG=[7;10];

Rgumax=zeros(ngon,1);
Rgumax(BaseG)=0*Pgmax(BaseG);
Rgumax(PearkG)=0.2*Pgmax(PearkG);
Rgdmax=Rgumax;


coe=1;
Rdsumax=0.1*dsmax*coe;
Rdsdmax=Rdsumax;



%cost
Og=[15;12.46;15;12.46;16;13.58;18.57;12;12;0];

Us=[26];

ODG=[13];

Uc=[23.5;20;25;21;22.3;22];

Vc=1000*ones(ndc,1);

pi_=pi*ones(nb,1);

Bdsmax=30;

%场景数
%sigma=100MW
% K=5;
% 
% phi=[0.1;0.2;0.4;0.2;0.1];
% 
% Wact_base=[124.2,213.65,300,386.35,475.8]/baseMVA;

% K=1;
% phi=[1];
% Wact_base=[300]/baseMVA;

%sigma=20MW
K=5;
Wact_base=[264.8,282.75,300,317.25,335.2]/baseMVA;
phi=[0.1;0.2;0.4;0.2;0.1];

Wcoe=[1.5,1.4,1.3,1.2,1.1,1,0.8,0.7,0.6,0.55,0.5,0.45,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.35,1.4];

Wact=zeros(T,K);
for t=1:T
    Wact(t,:)=Wact_base*Wcoe(t);
end


M=1e3;

%variables
%UL problem
Bds=sdpvar(nds,1,T,'full');
PDG=sdpvar(nds,1,T,'full');
ds=sdpvar(nds,1,T,'full');
Pes=sdpvar(nds,1,T,'full');

%LL Problem
Pg=sdpvar(ngon,1,T,'full');
rgk=sdpvar(ngon,K,T,'full');

dsn=sdpvar(nds,1,T,'full');
rdsk=sdpvar(nds,K,T,'full');

w=sdpvar(nw,1,T,'full');
wpk=sdpvar(nw,K,T,'full');

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
Urdsk1_down=sdpvar(nds,K,T,'full');
%g10
Urdsk1_up=sdpvar(nds,K,T,'full');

%g11
Urdsk2_down=sdpvar(nds,K,T,'full');
%g12
Urdsk2_up=sdpvar(nds,K,T,'full');


%g23
Uw_down=sdpvar(nw,1,T,'full');
%g24
Uw_up=sdpvar(nw,1,T,'full');
%g25
Uwpk_down=sdpvar(nw,K,T,'full');
%g26
Uwpk_up=sdpvar(nw,K,T,'full');

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
binUrdsk1_down=binvar(nds,K,T,'full');
%g10
binUrdsk1_up=binvar(nds,K,T,'full');

%g11
binUrdsk2_down=binvar(nds,K,T,'full');
%g12
binUrdsk2_up=binvar(nds,K,T,'full');


%g23
binUw_down=binvar(nw,1,T,'full');
%g24
binUw_up=binvar(nw,1,T,'full');
%g25
binUwpk_down=binvar(nw,K,T,'full');
%g26
binUwpk_up=binvar(nw,K,T,'full');


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
Constraints=[Constraints,0<=Bds(:,1,t)<=Bdsmax];
%Constraints=[Constraints,Bds==Us];
Constraints=[Constraints,0<=PDG(:,1,t)<=PDGmax];
Constraints=[Constraints,0<=ds(:,1,t)<=dsmax(:,t)];
%储能约束
Constraints=[Constraints,-PESmax<=Pes(:,1,t)<=PESmax];
CES=CES0;
for tt=1:t
    CES=CES-Pes(:,1,tt);
end
Constraints=[Constraints,0<=CES<=CESmax];
%上层功率平衡约束
Constraints=[Constraints,ds(:,1,t) == dsn(:,1,t) + PDG(:,1,t) + Pes(:,1,t)];


%LL Constraints
%KKT Pg
Constraints=[Constraints,(Og - Cg'*lamda(:,1,t) - Ug_down(:,1,t) + Ug_up(:,1,t) - sum(Urgk1_down(:,:,t),2) +  sum(Urgk1_up(:,:,t),2)) == 0];
%KKT rgk
for k=1:K
    Constraints=[Constraints,(phi(k)*Og - Cg'*lamdak(:,k,t) - Urgk1_down(:,k,t) + Urgk1_up(:,k,t) - Urgk2_down(:,k,t) + Urgk2_up(:,k,t)) == 0];
end

%KKT dsn
Constraints=[Constraints,(-Bds(:,1,t) + Cds'*lamda(:,1,t) - Udsn_down(:,1,t) + Udsn_up(:,1,t) + sum(Urdsk1_down(:,:,t),2) - sum(Urdsk1_up(:,:,t),2)) == 0];
%KKT rdsk
for k=1:K
    Constraints=[Constraints,(phi(k)*Us - Cds'*lamdak(:,k,t) - Urdsk1_down(:,k,t) + Urdsk1_up(:,k,t) -Urdsk2_down(:,k,t) + Urdsk2_up(:,k,t)) == 0];
end

%KKT w
Constraints=[Constraints,(-Cw'*lamda(:,1,t) + sum(Cw'*lamdak(:,:,t),2) - Uw_down(:,1,t) + Uw_up(:,1,t)) == 0];
%KKT wpk
for k=1:K
   Constraints=[Constraints,(Cw'*lamdak(:,k,t) - Uwpk_down(:,k,t) + Uwpk_up(:,k,t))==0]; 
end

%KKT theta
Constraints=[Constraints,(Bbus'*lamda(:,1,t) - sum(Bbus'*lamdak(:,:,t),2) - Utheta_down(:,1,t) + Utheta_up(:,1,t))==0];
%KKT thetak
for k=1:K
    Constraints=[Constraints,(Bbus'*lamdak(:,k,t) - Uthetak_down(:,k,t) + Uthetak_up(:,k,t)) == 0];
end

%h1
Constraints=[Constraints,(Cg*Pg(:,1,t) + Cw*w(:,1,t) - Cds*dsn(:,1,t) - Cdc*dc(:,t) - Bbus*theta(:,1,t)) ==0];

for k=1:K
    %h2
    Constraints=[Constraints,(Cg*rgk(:,k,t) + Cds*rdsk(:,k,t) ...
        + Cw*(Wact(t,k)-w(:,1,t)-wpk(:,k,t)) - Bbus*(thetak(:,k,t)-theta(:,1,t)))==0];
end
%g1
Constraints=[Constraints,Pg(:,1,t)>=0, Ug_down(:,1,t)>=0, Pg(:,1,t)<=binUg_down(:,1,t)*M , Ug_down(:,1,t)<=(1-binUg_down(:,1,t))*M ];
%g2
Constraints=[Constraints,(Pgmax-Pg(:,1,t))>=0, Ug_up(:,1,t)>=0, (Pgmax-Pg(:,1,t))<=binUg_up(:,1,t)*M , Ug_up(:,1,t)<=(1-binUg_up(:,1,t))*M ];

for k=1:K
    %g3
    Constraints=[Constraints,(rgk(:,k,t)+Pg(:,1,t))>=0,Urgk1_down(:,k,t)>=0,(rgk(:,k,t)+Pg(:,1,t))<=binUrgk1_down(:,k,t)*M ,Urgk1_down(:,k,t)<=(1-binUrgk1_down(:,k,t))*M ];
    %g4
    Constraints=[Constraints,(Pgmax-Pg(:,1,t)-rgk(:,k,t))>=0,Urgk1_up(:,k,t)>=0,(Pgmax-Pg(:,1,t)-rgk(:,k,t))<=binUrgk1_up(:,k,t)*M ,Urgk1_up(:,k,t)<=(1-binUrgk1_up(:,k,t))*M ];
    %g5
    Constraints=[Constraints,(rgk(:,k,t)+Rgdmax)>=0,Urgk2_down(:,k,t)>=0,(rgk(:,k,t)+Rgdmax)<=binUrgk2_down(:,k,t)*M ,Urgk2_down(:,k,t)<=(1-binUrgk2_down(:,k,t))*M ];
    %g6
    Constraints=[Constraints,(Rgumax-rgk(:,k,t))>=0,Urgk2_up(:,k,t)>=0,(Rgumax-rgk(:,k,t))<=binUrgk2_up(:,k,t)*M ,Urgk2_up(:,k,t)<=(1-binUrgk2_up(:,k,t))*M ];
end

%g7
Constraints=[Constraints,(dsn(:,1,t)+PDGmax)>=0,Udsn_down(:,1,t)>=0,(dsn(:,1,t)+PDGmax)<=binUdsn_down(:,1,t)*M ,Udsn_down(:,1,t)<=(1-binUdsn_down(:,1,t))*M ];
%g8
Constraints=[Constraints,(dsmax(:,t)-dsn(:,1,t))>=0,Udsn_up(:,1,t)>=0,(dsmax(:,t)-dsn(:,1,t))<=binUdsn_up(:,1,t)*M ,Udsn_up(:,1,t)<=(1-binUdsn_up(:,1,t))*M ];

for k=1:K
    %g9
    Constraints=[Constraints,(rdsk(:,k,t)+dsmax(:,t)-dsn(:,1,t))>=0,Urdsk1_down(:,k,t)>=0,(rdsk(:,k,t)+dsmax(:,t)-dsn(:,1,t))<=binUrdsk1_down(:,k,t)*M ,Urdsk1_down(:,k,t)<=(1-binUrdsk1_down(:,k,t))*M ];
    %g10
    Constraints=[Constraints,(dsn(:,1,t)-rdsk(:,k,t))>=0,Urdsk1_up(:,k,t)>=0,(dsn(:,1,t)-rdsk(:,k,t))<=binUrdsk1_up(:,k,t)*M ,Urdsk1_up(:,k,t)<=(1-binUrdsk1_up(:,k,t))*M ];
    %g11
    Constraints=[Constraints,(rdsk(:,k,t)+Rdsdmax(:,t))>=0,Urdsk2_down(:,k,t)>=0,(rdsk(:,k,t)+Rdsdmax(:,t))<=binUrdsk2_down(:,k,t)*M ,Urdsk2_down(:,k,t)<=(1-binUrdsk2_down(:,k,t))*M ];
    %g12
    Constraints=[Constraints,(Rdsumax(:,t)-rdsk(:,k,t))>=0,Urdsk2_up(:,k,t)>=0,(Rdsumax(:,t)-rdsk(:,k,t))<=binUrdsk2_up(:,k,t)*M ,Urdsk2_up(:,k,t)<=(1-binUrdsk2_up(:,k,t))*M ];
end

%g23
Constraints=[Constraints,w(:,1,t)>=0,Uw_down(:,1,t)>=0,w(:,1,t)<=binUw_down(:,1,t)*M ,Uw_down(:,1,t)<=(1-binUw_down(:,1,t))*M ];
%g24
Constraints=[Constraints,(Wmax-w(:,1,t))>=0,Uw_up(:,1,t)>=0,(Wmax-w(:,1,t))<=binUw_up(:,1,t)*M ,Uw_up(:,1,t)<=(1-binUw_up(:,1,t))*M ];
for k=1:K
    %g25
    Constraints=[Constraints,wpk(:,k,t)>=0,Uwpk_down(:,k,t)>=0,wpk(:,k,t)<=binUwpk_down(:,k,t)*M ,Uwpk_down(:,k,t)<=(1-binUwpk_down(:,k,t))*M ];
    %g26
    Constraints=[Constraints,(Wact(t,k)-wpk(:,k,t))>=0,Uwpk_up(:,k,t)>=0,(Wact(t,k)-wpk(:,k,t))<=binUwpk_up(:,k,t)*M ,Uwpk_up(:,k,t)<=(1-binUwpk_up(:,k,t))*M ];
end
%g29
Constraints=[Constraints,(theta(:,1,t)+pi)>=0,Utheta_down(:,1,t)>=0,(theta(:,1,t)+pi)<=binUtheta_down(:,1,t)*M ,Utheta_down(:,1,t)<=(1-binUtheta_down(:,1,t))*M ];
%g30
Constraints=[Constraints,(pi-theta(:,1,t))>=0,Utheta_up(:,1,t)>=0,(pi-theta(:,1,t))<=binUtheta_up(:,1,t)*M ,Utheta_up(:,1,t)<=(1-binUtheta_up(:,1,t))*M ];
for k=1:K
    %g31
    Constraints=[Constraints,(thetak(:,k,t)+pi)>=0,Uthetak_down(:,k,t)>=0,(thetak(:,k,t)+pi)<=binUthetak_down(:,k,t)*M ,Uthetak_down(:,k,t)<=(1-binUthetak_down(:,k,t))*M ];
    %g32
    Constraints=[Constraints,(pi-thetak(:,k,t))>=0,Uthetak_up(:,k,t)>=0,(pi-thetak(:,k,t))<=binUthetak_up(:,k,t)*M ,Uthetak_up(:,k,t)<=(1-binUthetak_up(:,k,t))*M ];
end

%ref bus
Constraints=[Constraints,theta(1,1,t)==0];
for k=1:K
    Constraints=[Constraints,thetak(1,k,t)==0];
end

end

Objective=0;

for t=1:T
    Objective = Objective -(    sum( diag( (-lamdak(:,:,t)')*Cw*Wact(t,:) ) ) + lamda(:,:,t)'*Cdc*dc(:,t)...
    - Ug_up(:,:,t)'*Pgmax - sum(Urgk1_up(:,:,t)'*Pgmax) - sum(Urgk2_down(:,:,t)'*Rgdmax) - sum(Urgk2_up(:,:,t)'*Rgumax)...
    - sum(Urdsk1_down(:,:,t)'*dsmax(:,t))...
    -Uw_up(:,:,t)'*Wmax - sum( diag(Uwpk_up(:,:,t)'*Wact(t,:)) )...
    -Utheta_down(:,:,t)'*pi_ - Utheta_up(:,:,t)'*pi_ - sum(Uthetak_down(:,:,t)'*pi_) - sum(Uthetak_up(:,:,t)'*pi_) -Uc'*dc(:,t)...
    - Og'*Pg(:,:,t) + Uc'*dc(:,t) -( Og'*rgk(:,:,t))*phi...
    + Us'*(Cds'*Cds)*ds(:,:,t) - Us'*(Cds'*Cds)*rdsk(:,:,t)*phi - ODG'*PDG(:,:,t));
end


%ops = sdpsettings('solver','cplex','verbose',2);
ops = sdpsettings('solver','gurobi','verbose',2);
ops.gurobi.MIPGap=0.015;
optimize(Constraints,Objective,ops);


double(Objective)
% double(Bds)
% double(lamda)
% double(Pg)
% double(dsn)
% double(ds)
% double(w)

profit=zeros(1,T);

for t=1:T
    profit(t) = double((Cds*Us)'*Cds*ds(:,:,t)-lamda(:,:,t)'*Cds*dsn(:,:,t)-(Cds*Us)'*Cds*rdsk(:,:,t)*phi+sum(diag(lamdak(:,:,t)'*Cds*rdsk(:,:,t)))-ODG'*PDG(:,:,t));
end

Bds=reshape(double(Bds),nds,T);
priceDA=reshape(double(lamda),nb,T);
Pg=reshape(double(Pg),ngon,T);
dsn=reshape(double(dsn),nds,T);
ds=reshape(double(ds),nds,T);
w=reshape(double(w),nw,T);
Pes=reshape(double(Pes),nds,T);
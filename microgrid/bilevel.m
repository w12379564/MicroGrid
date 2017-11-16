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


%微网参数
dsmax=[0.221]/baseMVA;
PDGmax=[0.2]/baseMVA;
Rdgumax=0.1*PDGmax;
Rdgdmax=Rdgumax;
Wmax=0.23/baseMVA;
Us=[26];
ODG=[13];

%联络线
PTLmax=0.2/baseMVA;


Pgmax=[40;152;40;152;300;591;60;400;400;300]/baseMVA;
dc=[110.16;98.94;174.42;197.88;184.62;130.56]/baseMVA;

BaseG=[1;2;3;4;5;6;8;9];
PearkG=[7;10];

Rgumax=zeros(ngon,1);
Rgumax(BaseG)=0*Pgmax(BaseG);
Rgumax(PearkG)=0.2*Pgmax(PearkG);
Rgdmax=Rgumax;

coe=1;
Rds=[1];
Rdsumax=zeros(nds,1);
Rdsumax(Rds)=0.1*dsmax(Rds)*coe;
Rdsdmax=Rdsumax;


%cost
Og=[15;12.46;15;12.46;16;13.58;18.57;12;12;0];

Uc=[23.5;20;25;21;22.3;22];

Vc=1000*ones(ndc,1);

pi_=pi*ones(nb,1);

Bdsmax=30;

%场景数
%sigma=100MW
K=5;
Wact=[0.02,0.04,0.09,0.12,0.18]/baseMVA;
phi=[0.1;0.2;0.4;0.2;0.1];

M=1e3;

%variables
%UL problem
Bds=sdpvar(nds,1,'full');
PDG=sdpvar(nds,1,'full');
ds=sdpvar(nds,1,'full');

rdsk=sdpvar(nds,K,'full');
rdgk=sdpvar(nds,K,'full');
w=sdpvar(nds,1,'full');
wpk=sdpvar(nds,K,'full');

%LL Problem
%variables
Pg=sdpvar(ngon,1,'full');
rgk=sdpvar(ngon,K,'full');

dsn=sdpvar(nds,1,'full');
dsnk=sdpvar(nds,K,'full');

theta=sdpvar(nb,1,'full');
thetak=sdpvar(nb,K,'full');

%dual variables
%h1
lamda=sdpvar(nb,1,'full');
%h2
lamdak=sdpvar(nb,K,'full');

%g1
Ug_down=sdpvar(ngon,1,'full');
%g2
Ug_up=sdpvar(ngon,1,'full');

%g3
Urgk1_down=sdpvar(ngon,K,'full');
%g4
Urgk1_up=sdpvar(ngon,K,'full');

%g5
Urgk2_down=sdpvar(ngon,K,'full');
%g6
Urgk2_up=sdpvar(ngon,K,'full');

%g7
Udsn_down=sdpvar(nds,1,'full');
%g8
Udsn_up=sdpvar(nds,1,'full');

%g9
Udsnk1_down=sdpvar(nds,K,'full');
%g10
Udsnk1_up=sdpvar(nds,K,'full');

%g11
Udsnk2_down=sdpvar(nds,K,'full');
%g12
Udsnk2_up=sdpvar(nds,K,'full');


%g29
Utheta_down=sdpvar(nb,1,'full');
%g30
Utheta_up=sdpvar(nb,1,'full');
%g31
Uthetak_down=sdpvar(nb,K,'full');
%g32
Uthetak_up=sdpvar(nb,K,'full');


%bin variables
%g1
binUg_down=binvar(ngon,1,'full');
%g2
binUg_up=binvar(ngon,1,'full');

%g3
binUrgk1_down=binvar(ngon,K,'full');
%g4
binUrgk1_up=binvar(ngon,K,'full');

%g5
binUrgk2_down=binvar(ngon,K,'full');
%g6
binUrgk2_up=binvar(ngon,K,'full');

%g7
binUdsn_down=binvar(nds,1,'full');
%g8
binUdsn_up=binvar(nds,1,'full');

%g9
binUdsnk1_down=binvar(nds,K,'full');
%g10
binUdsnk1_up=binvar(nds,K,'full');

%g11
binUdsnk2_down=binvar(nds,K,'full');
%g12
binUdsnk2_up=binvar(nds,K,'full');


%g29
binUtheta_down=binvar(nb,1,'full');
%g30
binUtheta_up=binvar(nb,1,'full');
%g31
binUthetak_down=binvar(nb,K,'full');
%g32
binUthetak_up=binvar(nb,K,'full');

Constraints=[];
%UL Constraints
Constraints=[Constraints,0<=Bds<=Bdsmax];
% Constraints=[Constraints,Bds==Us];
Constraints=[Constraints,0<=PDG<=PDGmax];
Constraints=[Constraints,0<=ds<=dsmax];
Constraints=[Constraints,0<=w<=Wmax];
Constraints=[Constraints,0<=wpk<=w];

%微网功率平衡
Constraints=[Constraints,(dsn+PDG+w-ds)==0];
for k=1:K
Constraints=[Constraints,(rdgk(:,k)+rdsk(:,k)-dsnk(:,k)+(Wact(k)-w-wpk(:,k)))==0];
Constraints=[Constraints,-(dsmax-ds)<=rdsk(:,k)<=ds];
Constraints=[Constraints,-PDG<=rdgk(:,k)<=(PDGmax-PDG)];
Constraints=[Constraints,-Rdgdmax<=rdgk(:,k)<=Rdgumax];
end

%LL Constraints
%KKT Pg
Constraints=[Constraints,(Og - Cg'*lamda - Ug_down + Ug_up - sum(Urgk1_down,2) +  sum(Urgk1_up,2)) == 0];
%KKT rgk
for k=1:K
    Constraints=[Constraints,(phi(k)*Og - Cg'*lamdak(:,k) - Urgk1_down(:,k) + Urgk1_up(:,k) - Urgk2_down(:,k) + Urgk2_up(:,k)) == 0];
end

%KKT dsn
Constraints=[Constraints,(-Bds + Cds'*lamda - Udsn_down + Udsn_up + sum(Udsnk1_down,2) - sum(Udsnk1_up,2)) == 0];
%KKT dsnk
for k=1:K
    Constraints=[Constraints,(phi(k)*Us - Cds'*lamdak(:,k) - Udsnk1_down(:,k) + Udsnk1_up(:,k) -Udsnk2_down(:,k) + Udsnk2_up(:,k)) == 0];
end


%KKT theta
Constraints=[Constraints,(Bbus'*lamda - sum(Bbus'*lamdak,2) - Utheta_down + Utheta_up)==0];
%KKT thetak
for k=1:K
    Constraints=[Constraints,(Bbus'*lamdak(:,k) - Uthetak_down(:,k) + Uthetak_up(:,k)) == 0];
end

%h1
Constraints=[Constraints,(Cg*Pg - Cds*dsn - Cdc*dc - Bbus*theta) ==0];

for k=1:K
    %h2
    Constraints=[Constraints,(Cg*rgk(:,k) + Cds*dsnk(:,k) - Bbus*(thetak(:,k)-theta))==0];
end
%g1
Constraints=[Constraints,Pg>=0, Ug_down>=0, Pg<=binUg_down*M , Ug_down<=(1-binUg_down)*M ];
%g2
Constraints=[Constraints,(Pgmax-Pg)>=0, Ug_up>=0, (Pgmax-Pg)<=binUg_up*M , Ug_up<=(1-binUg_up)*M ];

for k=1:K
    %g3
    Constraints=[Constraints,(rgk(:,k)+Pg)>=0,Urgk1_down(:,k)>=0,(rgk(:,k)+Pg)<=binUrgk1_down(:,k)*M ,Urgk1_down(:,k)<=(1-binUrgk1_down(:,k))*M ];
    %g4
    Constraints=[Constraints,(Pgmax-Pg-rgk(:,k))>=0,Urgk1_up(:,k)>=0,(Pgmax-Pg-rgk(:,k))<=binUrgk1_up(:,k)*M ,Urgk1_up(:,k)<=(1-binUrgk1_up(:,k))*M ];
    %g5
    Constraints=[Constraints,(rgk(:,k)+Rgdmax)>=0,Urgk2_down(:,k)>=0,(rgk(:,k)+Rgdmax)<=binUrgk2_down(:,k)*M ,Urgk2_down(:,k)<=(1-binUrgk2_down(:,k))*M ];
    %g6
    Constraints=[Constraints,(Rgumax-rgk(:,k))>=0,Urgk2_up(:,k)>=0,(Rgumax-rgk(:,k))<=binUrgk2_up(:,k)*M ,Urgk2_up(:,k)<=(1-binUrgk2_up(:,k))*M ];
end

%g7
Constraints=[Constraints,(dsn+PTLmax)>=0,Udsn_down>=0,(dsn+PTLmax)<=binUdsn_down*M ,Udsn_down<=(1-binUdsn_down)*M ];
%g8
Constraints=[Constraints,(PTLmax-dsn)>=0,Udsn_up>=0,(PTLmax-dsn)<=binUdsn_up*M ,Udsn_up<=(1-binUdsn_up)*M ];

for k=1:K
    %g9
    Constraints=[Constraints,(dsnk(:,k)+PTLmax-dsn)>=0,Udsnk1_down(:,k)>=0,(dsnk(:,k)+PTLmax-dsn)<=binUdsnk1_down(:,k)*M ,Udsnk1_down(:,k)<=(1-binUdsnk1_down(:,k))*M ];
    %g10
    Constraints=[Constraints,(dsn+PTLmax-dsnk(:,k))>=0,Udsnk1_up(:,k)>=0,(dsn+PTLmax-dsnk(:,k))<=binUdsnk1_up(:,k)*M ,Udsnk1_up(:,k)<=(1-binUdsnk1_up(:,k))*M ];
    %g11
    Constraints=[Constraints,(dsnk(:,k)+Rdsdmax)>=0,Udsnk2_down(:,k)>=0,(dsnk(:,k)+Rdsdmax)<=binUdsnk2_down(:,k)*M ,Udsnk2_down(:,k)<=(1-binUdsnk2_down(:,k))*M ];
    %g12
    Constraints=[Constraints,(Rdsumax-dsnk(:,k))>=0,Udsnk2_up(:,k)>=0,(Rdsumax-dsnk(:,k))<=binUdsnk2_up(:,k)*M ,Udsnk2_up(:,k)<=(1-binUdsnk2_up(:,k))*M ];
end

%g29
Constraints=[Constraints,(theta+pi)>=0,Utheta_down>=0,(theta+pi)<=binUtheta_down*M ,Utheta_down<=(1-binUtheta_down)*M ];
%g30
Constraints=[Constraints,(pi-theta)>=0,Utheta_up>=0,(pi-theta)<=binUtheta_up*M ,Utheta_up<=(1-binUtheta_up)*M ];
for k=1:K
    %g31
    Constraints=[Constraints,(thetak(:,k)+pi)>=0,Uthetak_down(:,k)>=0,(thetak(:,k)+pi)<=binUthetak_down(:,k)*M ,Uthetak_down(:,k)<=(1-binUthetak_down(:,k))*M ];
    %g32
    Constraints=[Constraints,(pi-thetak(:,k))>=0,Uthetak_up(:,k)>=0,(pi-thetak(:,k))<=binUthetak_up(:,k)*M ,Uthetak_up(:,k)<=(1-binUthetak_up(:,k))*M ];
end


%ref bus
Constraints=[Constraints,theta(1)==0];
for k=1:K
    Constraints=[Constraints,thetak(1,k)==0];
end

dualObj= lamda'*Cdc*dc...
    - Ug_up'*Pgmax - sum(Urgk1_up'*Pgmax) - sum(Urgk2_down'*Rgdmax) - sum(Urgk2_up'*Rgumax)...
    -Udsn_down'*PTLmax - Udsn_up'*PTLmax - sum(Udsnk1_down'*PTLmax) - sum(Udsnk1_up'*PTLmax) - sum(Udsnk2_down'*Rdsdmax) - sum(Udsnk2_up'*Rdsumax)...
    -Utheta_down'*pi_ - Utheta_up'*pi_ - sum(Uthetak_down'*pi_) - sum(Uthetak_up'*pi_) -Uc'*dc;

T=lamda'*Cdc*dc...
    - Ug_up'*Pgmax - sum(Urgk1_up'*Pgmax) - sum(Urgk2_down'*Rgdmax) - sum(Urgk2_up'*Rgumax)...
    -Utheta_down'*pi_ - Utheta_up'*pi_ - sum(Uthetak_down'*pi_) - sum(Uthetak_up'*pi_)...
    - Og'*Pg -( Og'*rgk)*phi;

UL = T + Us*ds - ODG*PDG - (Us*rdsk + ODG*rdgk)*phi;

Objective = -UL;

ops = sdpsettings('solver','gurobi','verbose',2);
optimize(Constraints,Objective,ops);

double(Bds)
double(lamda)
double(Pg)
double(dsn)
ds=double(ds)
double(w)
double(PDG)
double((Wact-wpk)*phi);
double(Objective)

double(Og'*Pg-Bds'*dsn-Uc'*dc...
    +(Og'*rgk+Us'*dsnk)*phi)

double(dualObj)
double(Us*ds-lamda'*Cds*dsn-ODG*PDG-Us*rdsk*phi+sum(diag(lamdak'*Cds*dsnk)) - ODG*rdgk*phi)

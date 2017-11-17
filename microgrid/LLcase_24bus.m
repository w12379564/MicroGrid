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

%上限
Pgmax=[40;152;40;152;300;591;60;400;400;300]/baseMVA;
dsmax=[0.221]/baseMVA;
PDGmax=[0.2]/baseMVA;
dc=[110.16;98.94;174.42;197.88;184.62;130.56]/baseMVA;
Wmax=0.13/baseMVA;

%联络线
PTLmax=0.2/baseMVA;

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

Us=[26];

ODG=[13];

Uc=[23.5;20;25;21;22.3;22];

Vc=1000*ones(ndc,1);

pi_=pi*ones(nb,1);

Bdsmax=30;

%场景数
%sigma=100MW
K=5;
Wact=[124.2,213.65,300,386.35,475.8]/baseMVA;
phi=[0.1;0.2;0.4;0.2;0.1];

% %sigma=20MW
% K=5;
% Wact=[264.8,282.75,300,317.25,335.2]/baseMVA;
% phi=[0.1;0.2;0.4;0.2;0.1];

Bds=12*ones(nds,1);

%variables
Pg=sdpvar(ngon,1,'full');
rgk=sdpvar(ngon,K,'full');

dsn=sdpvar(nds,1,'full');
dsnk=sdpvar(nds,K,'full');

theta=sdpvar(nb,1,'full');
thetak=sdpvar(nb,K,'full');

Constraints=[];
%h1
Constraints=[Constraints,(Cg*Pg - Cds*dsn - Cdc*dc - Bbus*theta) ==0];

for k=1:K
    %h2
    Constraints=[Constraints,(Cg*rgk(:,k) + Cds*dsnk(:,k) - Bbus*(thetak(:,k)-theta)) == 0];
end
%g1
Constraints=[Constraints,Pg>=0];
%g2
Constraints=[Constraints,Pgmax-Pg>=0];

for k=1:K
    %g3
    Constraints=[Constraints,rgk(:,k)+Pg>=0];
    %g4
    Constraints=[Constraints,Pgmax-Pg-rgk(:,k)>=0];
    %g5
    Constraints=[Constraints,rgk(:,k)+Rgdmax>=0];
    %g6
    Constraints=[Constraints,Rgumax-rgk(:,k)>=0];
end

%g7
Constraints=[Constraints,dsn+PTLmax>=0];
%g8
Constraints=[Constraints,PTLmax-dsn>=0];

for k=1:K
    %g9
    Constraints=[Constraints,dsnk(:,k)+PTLmax-dsn>=0];
    %g10
    Constraints=[Constraints,dsn+PTLmax-dsnk(:,k)>=0];
    %g11
    Constraints=[Constraints,dsnk(:,k)+Rdsdmax>=0];
    %g12
    Constraints=[Constraints,Rdsumax-dsnk(:,k)>=0];
end


%g29
Constraints=[Constraints,theta+pi>=0,theta(1)==0];
%g30
Constraints=[Constraints,pi-theta>=0];
for k=1:K
    %g31
    Constraints=[Constraints,thetak(:,k)+pi>=0,thetak(1,k)==0];
    %g32
    Constraints=[Constraints,pi-thetak(:,k)>=0];
end


Objective=Og'*Pg-Bds'*dsn-Uc'*dc...
    +(Og'*rgk+Bds'*dsnk)*phi;



% Objective=Og'*Pg-Bsd'*ds-Uc'*dc;

optimize(Constraints,Objective);

double(Pg)
double(dsn)
double(dc)
double(Objective)




function [ Bbus, Bf,Bft, Pbusinj, Pfinj ] = calBdc( baseMVA, bus, branch )
%copy from matpower makeBdc.m
%   此处显示详细说明
%% define named indices into bus, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = indexBus;

[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = indexBrch;
%%
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of lines

%% check that bus numbers are equal to indices to bus (one set of bus numbers)
if any(bus(:, BUS_I) ~= (1:nb)')
    error('makeBdc: buses must be numbered consecutively in bus matrix')
end

%% for each branch, compute the elements of the branch B matrix and the phase
%% shift "quiescent" injections, where
%%
%%      | Pf |   | Bff  Bft |   | Vaf |   | Pfinj |
%%      |    | = |          | * |     | + |       |
%%      | Pt |   | Btf  Btt |   | Vat |   | Ptinj |
%%
stat = branch(:, BR_STATUS);                    %% ones at in-service branches
b = stat ./ branch(:, BR_X);                    %% series susceptance
tap = ones(nl, 1);                              %% default tap ratio = 1
i = find(branch(:, TAP));                       %% indices of non-zero tap ratios
tap(i) = branch(i, TAP);                        %% assign non-zero tap ratios
b = b ./ tap;


%% build connection matrix Cft = Cf - Ct for line and from - to buses
f = branch(:, F_BUS);                           %% list of "from" buses
t = branch(:, T_BUS);                           %% list of "to" buses
i = [(1:nl)'; (1:nl)'];                         %% double set of row indices
Cft = sparse(i, [f;t], [ones(nl, 1); -ones(nl, 1)], nl, nb);    %% connection matrix

Bft=sparse([f;t],[t;f],[b;b]);  %%

%% build Bf such that Bf * Va is the vector of real branch powers injected
%% at each branch's "from" bus
Bf = sparse(i, [f; t], [b; -b]);    % = spdiags(b, 0, nl, nl) * Cft;

%% build Bbus
Bbus = Cft' * Bf;

%% build phase shift injection vectors
Pfinj = b .* (-branch(:, SHIFT) * pi/180);      %% injected at the from bus ...
    % Ptinj = -Pfinj;                           %% ... and extracted at the to bus
Pbusinj = Cft' * Pfinj;                         %% Pbusinj = Cf * Pfinj + Ct * Ptinj;


end


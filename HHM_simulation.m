function solution = HH_simul(Data, dic1, dic2, Iapp, gpe_gpe, gpe_stn, stn_gpe, gg, gs, sg)
% HH_simul implements Hodgkin-Huxley dynamics for a GPe-STN network
%
% INPUTS:
%   Data      - State vector including voltages, gating vars, calcium, synaptic vars
%   dic1/dic2 - Parameter dictionaries for GPe and STN
%   Iapp      - External applied current
%   gpe_gpe, gpe_stn, stn_gpe - Synaptic connectivity matrices
%   gg, gs, sg - Synaptic conductances for GPe→GPe, STN→GPe, GPe→STN
%
% OUTPUT:
%   solution  - Time derivative of the state vector

%=============================%
%       Data Unpacking       %
%=============================%
N = length(gpe_gpe); % Number of neurons

% GPe variables
GPe   = Data(1:N);                  % Membrane potential
ngpe  = Data(N+1:2*N);              % n gating
hgpe  = Data(2*N+1:3*N);            % h gating
rgpe  = Data(3*N+1:4*N);            % r gating
CAgpe = Data(4*N+1:5*N);            % Intracellular calcium
sgpe  = Data(5*N+1:6*N);            % Synaptic gating

% STN variables
STN   = Data(6*N+1:7*N);
nstn  = Data(7*N+1:8*N);
hstn  = Data(8*N+1:9*N);
rstn  = Data(9*N+1:10*N);
CAstn = Data(10*N+1:11*N);
sstn  = Data(11*N+1:12*N);

%=============================%
%     Initialize Derivatives %
%=============================%
dx0 = zeros(N,1); dx1 = zeros(N,1); dx2 = zeros(N,1);
dx3 = zeros(N,1); dx4 = zeros(N,1); dx5 = zeros(N,1);

dy0 = zeros(N,1); dy1 = zeros(N,1); dy2 = zeros(N,1);
dy3 = zeros(N,1); dy4 = zeros(N,1); dy5 = zeros(N,1);

%=============================%
%     Main Update Loop       %
%=============================%
for i = 1:N
    %% ---------- GPe Neuron Dynamics ----------
    % Equilibrium gating variables
    ninf = equilibx(GPe(i), dic1("tetn"), dic1("sign"));
    hinf = equilibx(GPe(i), dic1("teth"), dic1("sigh"));
    rinf = equilibx(GPe(i), dic1("tetr"), dic1("sigr"));
    minf = equilibx(GPe(i), dic1("tetm"), dic1("sigm"));
    sinf = equilibx(GPe(i), dic1("tets"), dic1("sigs"));
    ainf = equilibx(GPe(i), dic1("teta"), dic1("siga"));

    % Ionic Currents
    IL  = dic1("gL")  * (GPe(i) - dic1("vL"));   % Leak
    IK  = coef(dic1("gK"),  ngpe(i), 1, 4) * (GPe(i) - dic1("vK"));   % K+
    INa = dic1("gNa") * coef(minf, hgpe(i), 3, 1) * (GPe(i) - dic1("vNa")); % Na+
    IT  = dic1("gT")  * coef(ainf, rgpe(i), 3, 1) * (GPe(i) - dic1("vCa")); % T-type Ca
    ICa = coef(dic1("gCa"), sinf, 1, 2) * (GPe(i) - dic1("vCa"));           % Ca
    Iahp = dic1("gAHP") * (CAgpe(i) / (dic1("k1") + CAgpe(i))) * (GPe(i) - dic1("vK")); % After-hyperpolarization

    % Synaptic Inputs to GPe
    Is_g = sg * (GPe(i) - dic1("vSG")) * dot(stn_gpe(i,:), sstn);
    Ig_g = gg * (GPe(i) - dic1("vGG")) * dot(gpe_gpe(:,i), sgpe);

    % Gating Variable Dynamics
    dx1(i) = dic1("fi_n") * (ninf - ngpe(i)) / taux(GPe(i), dic1("tau_n0"), dic1("tau_n1"), dic1("tet_tn"), dic1("sig_tn"));
    dx2(i) = dic1("fi_h") * (hinf - hgpe(i)) / taux(GPe(i), dic1("tau_h0"), dic1("tau_h1"), dic1("tet_th"), dic1("sig_th"));
    dx3(i) = dic1("fi_r") * (rinf - rgpe(i)) / dic1("tau_r");
    dx4(i) = dic1("eps")   * (-ICa - IT - dic1("kCa") * CAgpe(i));
    dx5(i) = dic1("alpha") * equilibx(GPe(i)-dic1("tetg"), dic1("tetHg"), dic1("sigHg")) * (1 - sgpe(i)) - dic1("beta") * sgpe(i);

    % GPe Membrane Potential Equation
    dx0(i) = -IL - IK - INa - ICa - IT - Iahp - Is_g - Ig_g + Iapp;

    %% ---------- STN Neuron Dynamics ----------
    ninfs  = equilibx(STN(i), dic2("tetn"), dic2("sign"));
    hinfs  = equilibx(STN(i), dic2("teth"), dic2("sigh"));
    rinfs  = equilibx(STN(i), dic2("tetr"), dic2("sigr"));
    minfs  = equilibx(STN(i), dic2("tetm"), dic2("sigm"));
    sinfs  = equilibx(STN(i), dic2("tets"), dic2("sigs"));
    ainfs  = equilibx(STN(i), dic2("teta"), dic2("siga"));
    binfs  = equilibb(rstn(i), dic2("tetb"), dic2("sigb"));

    ILs  = dic2("gL")  * (STN(i) - dic2("vL"));
    IKs  = coef(dic2("gK"),  nstn(i), 1, 4) * (STN(i) - dic2("vK"));
    INas = dic2("gNa") * coef(minfs, hstn(i), 3, 1) * (STN(i) - dic2("vNa"));
    ITs  = dic2("gT")  * coef(ainfs, binfs, 3, 2) * rstn(i) * (STN(i) - dic2("vCa"));
    ICas = coef(dic2("gCa"), sinfs, 1, 2) * (STN(i) - dic2("vCa"));
    Iahps = dic2("gAHP") * (CAstn(i) / (dic2("k1") + CAstn(i))) * (STN(i) - dic2("vK"));

    % Synaptic Input to STN
    Ig_sS = gs * (STN(i) - dic2("vGS")) * dot(gpe_stn(:,i), sgpe);

    % STN Gating Dynamics
    dy1(i) = dic2("fi_n") * (ninfs - nstn(i)) / taux(STN(i), dic2("tau_n0"), dic2("tau_n1"), dic2("tet_tn"), dic2("sig_tn"));
    dy2(i) = dic2("fi_h") * (hinfs - hstn(i)) / taux(STN(i), dic2("tau_h0"), dic2("tau_h1"), dic2("tet_th"), dic2("sig_th"));
    dy3(i) = dic2("fi_r") * (rinfs - rstn(i)) / taux(STN(i), dic2("tau_r0"), dic2("tau_r1"), dic2("tet_tr"), dic2("sig_tr"));
    dy4(i) = dic2("eps")  * (-ICas - ITs - dic2("kCa") * CAstn(i));
    dy5(i) = dic2("alpha") * equilibx(STN(i)-dic2("tetg"), dic2("tetHg"), dic2("sigHg")) * (1 - sstn(i)) - dic2("beta") * sstn(i);

    % STN Membrane Potential Equation
    dy0(i) = -ILs - IKs - INas - ICas - ITs - Iahps - Ig_sS + 10.7;
end

%=============================%
%     Pack the Derivatives   %
%=============================%
solution = [dx0; dx1; dx2; dx3; dx4; dx5; dy0; dy1; dy2; dy3; dy4; dy5];
end

% Script Sim_Graph
% WRITTEN BY: YIDA WANG
% November 26, 2019
clc
help Sim_Graph
% Initialize constant values
num_inv = 100e6;
num_dev = num_inv*2;
T = 300;            % K
k = 1.38e-23;       % J/K
q = 1.6e-19;        % C
eps0 = 8.854e-14;   % F/cm
epsSi = 11.68*eps0; % F/cm
epsOx = 3.9*eps0;   % F/cm
ni = 1.45e10;       % 1/cm^3
Nc = 10^25;         % 1/m^3
Nv = Nc;            % 1/m^3
Nsource = 0.1*Nc;   % 1/m^3
Ndrain = 0.1*Nv;    % 1/m^3
bni = 10^16;        % 1/m^3
tox = 2e-7;         % cm
Cox = epsOx/tox;    % F
Wn = 110e-7;        % cm
Ln = 22e-7;         % cm
Lp = Ln;            % cm
un = 1000;          % cm^2/V-s 
up = 0.5*un;        % cm^2/V-s
Vgs = 0;            % V
I0 = 1e-6;          % A
T_CLK = 1e-9;       % s
td = 10e-12;

FT1 = 0.0001;       % Initial failure target
FT2 = 0.1;          % Final failure target
pts = 1000;         % 1000 points
Pcir = linspace(FT1,FT2,pts);
Pdev = 1 - (1 - Pcir).^(1/num_dev);
Eb = -log(Pdev)*k*T;
Nchannel = (bni^2)./(Nsource*Pdev);
psiB = ((k*T)/q)*log(Nchannel/bni);
Vth = 2*psiB + sqrt(2*epsSi*q.*Nchannel*2.*psiB)/Cox;
Wp = (un/up)*Wn;
beta = un*Cox*Wn/Ln;
CL = Cox*(Wp + Wn)*Ln;
VDD = (CL*log(2)/td)/beta + Vth;
tsi = sqrt((2*epsSi*psiB*10^-6)./(q*Nchannel));
n = 1 + (epsSi/epsOx)*(tox./tsi);
Ids = I0*exp((q*(Vgs - Vth))./(n*k*T));
m = 0.1*num_dev;
Etotal = m*CL*VDD.^2 + num_dev.*VDD.*Ids*T_CLK;
figure(1)
plot(Pcir,VDD,'r');
title('Probability of Failure vs. Supply Voltage');
xlabel('Probability of Circuit Failure');
ylabel('Supply Voltage');
legend('VDD');
figure(2)
plot(Pcir,Etotal,'b');
title('Probability of Failure vs. Total Energy Dissipation');
xlabel('Probability of Circuit Failure');
ylabel('Total Energy Dissipation');
legend('Etotal');

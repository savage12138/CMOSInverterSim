% Script Sim_Delay_Graph
% WRITTEN BY: YIDA WANG
% November 26, 2019
clc
help Sim_Delay_Graph
% Initialize constant values
num_inv = 100e6;
num_dev = num_inv*2;
T = 300;            % K
k = 1.38e-23;       % J/K
q = 1.6e-19;        % C
eps0 = 8.854e-12;   % F/m
epsSi = 11.68*eps0; % F/m
epsOx = 3.9*eps0;   % F/m
ni = 1.45e10;       % 1/cm^3
Nc = 10^25;         % 1/m^3
Nv = Nc;            % 1/m^3
Nsource = 0.1*Nc;   % 1/m^3
Ndrain = 0.1*Nv;    % 1/m^3
bni = 10^16;        % 1/m^3
tox = 2e-9;         % m
Cox = epsOx/tox;    % F/m^2
Wn = 110e-9;        % m
Ln = 22e-9;         % m
Lp = Ln;            % m
un = 0.1;           % m^2/V-s 
up = 0.5*un;        % m^2/V-s
Vgs = 0;            % V
I0 = 1e-6;          % A
T_CLK = 1e-9;       % s
Pcir = 0.001;

DT1 = 1e-12;       % Initial failure target
DT2 = 1000e-12;          % Final failure target
pts = 1000;         % 1000 points
td = linspace(DT1,DT2,pts);
Pdev = 1 - (1 - Pcir)^(1/num_dev);
Eb = -log(Pdev)*k*T;
Nchannel = (bni^2)/(Nsource*Pdev);
psiB = ((k*T)/q)*log(Nchannel/bni);
Vth = 2*psiB + sqrt(2*epsSi*q*Nchannel*2*psiB)/Cox;
Wp = (un/up)*Wn;
beta = un*Cox*Wn/Ln;
CL = Cox*(Wp + Wn)*Ln;
VDD = (CL*log(2)./td)/beta + Vth;
tsi = sqrt((2*epsSi*psiB)/(q*Nchannel));
n = 1 + (epsSi/epsOx)*(tox/tsi);
Ids = I0*exp((q*(Vgs - Vth))/(n*k*T));
m = 0.1*num_dev;
Etotal = m*CL*VDD.^2 + num_dev*VDD*Ids*T_CLK;
subplot(1,2,1);
plot(td,VDD,'r');
title({'Delay', 'vs', 'Supply Voltage'});
xlabel('Delay (s)');
ylabel('Supply Voltage (V)');
legend('VDD');
subplot(1,2,2);
plot(td,Etotal,'b');
title({'Delay', 'vs', 'Total Energy Dissipation'});
xlabel('Delay (s)');
ylabel('Total Energy Dissipation (J)');
legend('Etotal');

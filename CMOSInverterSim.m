% Script CMOS Inverter Simulation
% WRITTEN BY: YIDA WANG
% November 26, 2019

clc
help CMOSInverterSim

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

% Request User Input Probability Failure of the Circuit
Pcir = input('Enter circuit failure probability: ');
disp('---------------------------------------------------------------')

disp('a) Find out the height of the energy barrier.')
disp('First, we need to calculate the probability failure')
disp('of a single device.')
Pdev = 1 - (1 - Pcir)^(1/num_dev);
fprintf('Probability failure of a single device is %.4e.\n',Pdev)
disp('Then use Pdev to calculate the height of the energy barrier.')
Eb = -log(Pdev)*k*T;
fprintf('The height of the energy barrier per device is %.4e Joules\n',Eb)
disp('---------------------------------------------------------------')

disp('b) Compute the channel doping(Nchannel).')
disp('Eb = -ln(Pdev)*k*T, Ebi = -ln((ni^2)/(NA*ND))*k*T')
disp('Pdev = (ni^2)/(NA*ND), Nchannel = ND = (ni^2)/(NA*Pdev)')
Nchannel = (bni^2)/(Nsource*Pdev);
fprintf('The channel doping is %.2e 1/m^3.\n',Nchannel)
disp('---------------------------------------------------------------')

disp('c) Find the threshold voltage.')
disp('First calculate psiB,')
psiB = ((k*T)/q)*log(Nchannel/bni);
disp('Then calculate Cox.')
disp('Then plug everything in the formula for threshold voltage.');
Vth = 2*psiB + sqrt(2*epsSi*q*Nchannel*2*psiB)/(Cox);
fprintf('The threshold voltage is %.5f V.\n',Vth)
disp('---------------------------------------------------------------')

disp('d) Find the size of the PFETs')
disp('Wp/Wn = un/up');
Wp = (un/up)*Wn;
fprintf('The size of the PFETs is %.2e m or %d nm.\n',Wp, Wp*10^9)
disp('---------------------------------------------------------------')

disp('e) Compute the supply voltage (VDD)')
disp('First, calculate beta.')
beta = un*Cox*Wn/Ln;
disp('Calculate CL.')
CL = Cox*(Wp + Wn)*Ln;
disp('td = CL*ln(2)/(beta*(VDD-Vth))')
% Request User Input delay
td = input('Enter target delay: ');
VDD = (CL*log(2)/td)/beta + Vth;
fprintf('The supply voltage VDD is %.5f V.\n',VDD)
disp('---------------------------------------------------------------')

disp('f) Compute the subthreshold leakage(Ids) at Vgs = 0.')
disp('First, calculate tsi.')
tsi = sqrt((2*epsSi*psiB)/(q*Nchannel));
n = 1 + (epsSi/epsOx)*(tox/tsi);
fprintf('Then, Calculate the subthreshold slop factor n = %.3f.\n',n)
disp('Plug everything in the Ids formula.')
Ids = I0*exp((q*(Vgs - Vth))/(n*k*T));
fprintf('Subthreshold leakage is %.3e A or %.2f fA.\n',Ids,Ids*10^15);
disp('---------------------------------------------------------------')

disp('g) Compute the total energy dissipation.')
disp('Total Energy = Dynamic Energy + Static Energy')
m = 0.1*num_inv;
Edy = m*CL*VDD^2;
fprintf('Dynamic Energy is %.3e J or %.2f pJ.\n',Edy, Edy*10^12)
Est = num_inv*VDD*Ids*T_CLK;
fprintf('Static Energy is %.3e J or %.2f pJ.\n',Est, Est*10^12)
Etotal = Edy + Est;
fprintf('Total Energy is %.3e J or %.2f pJ.\n',Etotal, Etotal*10^12)
disp('---------------------------------------------------------------')

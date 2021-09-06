clear all; close all; clc;
[c1, c2, c3, c4, c5,c6,c7] = MatlabDefaultPlotColors();
addpath('../LPM');
%% definitions
TsL = 1/100;
TsH = TsL/2;
n = 20;
R=3;
tperiod = 10;
Per = 5;
nT = Per-1;
%% systems analyzed
s = tf('s');
P = 1/(s^2+2*0.1*(2*pi*10)*s+(2*pi*10)^2);
Pss = ss(c2d(P,TsL));
%% input signal multisine
Np = tperiod/TsH; % amount of samples in period
Ntot = Per*Np; % total amount of samples
NnH = floor(Np/2); % nyquist freq
NnL = floor(NnH)/(TsL/TsH);
fH = linspace(0, 1 - 1/NnH, NnH) * (1/TsH)/2; % available frequencies for high sampling rate

fL = linspace(0,1-1/NnL,NnL)*(1/TsL)/2;
A = ones(NnL,1); % amplitude distribution


tp = (0:TsH:(Np-1)*TsH)';
t = (0:TsH:(Ntot-1)*TsH)';

% custom multisine
up = zeros(Np,1); % periodic signal
for k = 1:NnL
   up = up+A(k)*sin(2*pi*fL(k)*tp+rand*2*pi);
end
u = repmat(up,Per,1);

%% simulate
simoutput = sim('simFile');
yH = simoutput.output;
uL = simoutput.outputLow;

UL = fft(uL);
UH=fft(u);
figure
stem(fH,abs(UH(1:Per:Per*NnH))); hold on;
stem(fL,abs(UL(1:Per:Per*NnL)));
legend('input uH','input uL')
%% LPM
[P_LPM,THz] = LPMOpenLoopPeriodicFastBLA(u,yH,n,R,Per,nT);
%% plotting
figure
bodemag(Pss); hold on;
semilogx(fH*2*pi,20*log10(abs(squeeze(P_LPM))))
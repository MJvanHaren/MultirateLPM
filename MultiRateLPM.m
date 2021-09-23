clear all; close all; clc;
[c1, c2, c3, c4, c5,c6,c7] = MatlabDefaultPlotColors();
addpath('../LPM');
%% definitions
TsL = 1/100;
TsH = TsL/2;
n = 8;
degLPM=3;
tperiod = 10;
Per = 10;
nT = Per-1;
reference = 'simple';
%% systems analyzed
s = tf('s');
% GctH = (1/s^2+1/(s^2+2*0.1*(2*pi*10)*s+(2*pi*10)^2));
% C = 50*(1/(2/3*2*pi)*s+1)/(1/(6*2*pi)*s+1);
% P = P*1/(s^2+2*0.1*(2*pi*75)*s+(2*pi*75)^2)

G = 5.6e3*(0.9*s^2+2*s+5.4e4)/(1.1*0.9*s^4+(1.1+0.9)*2*s^3+(1.1+0.9)*5.4e4*s^2);
C = 0.05*(1/(1*2*pi)*s+1)/(1/(5*2*pi)*s+1);

Pdt = ss(c2d(G,TsH,'tustin'));
Cdt = ss(c2d(C,TsL,'tustin'));
%% input signal multisine
Np = tperiod/TsH;   % amount of samples in period

Ntot = Per*Np;      % total amount of samples
NtotL = Per*tperiod/TsL;

NnH = floor(Np/2);          % nyquist freq high sampling rate
NnL = floor(NnH)/(TsL/TsH); % nyquist freq low sampling rate

fH = linspace(0, 1 - 1/NnH, NnH) * (1/TsH)/2;   % available frequencies for high sampling rate
fL = linspace(0,1-1/NnL,NnL)*(1/TsL)/2;         % available frequencies for low sampling rate

tp = (0:TsH:(Np-1)*TsH)'; % periodic time signal
rp = zeros(Np,1); % periodic signal

if strcmp(reference,'high')
    A = ones(NnH,1); % amplitude distribution
    for k = 1:NnH
        rp = rp+A(k)*sin(2*pi*fH(k)*tp+rand*2*pi);
    end
elseif strcmp(reference,'simple')
%     A = ones(5,1);
%     fs = [1 5 10 25 60];
    A = ones(3,1);
    fs = [5 10 20];
    for k = 1:length(A)
        rp = rp+A(k)*sin(2*pi*fs(k)*tp+rand*2*pi);
    end
else
    A = ones(NnL,1); % amplitude distribution
%     A = linspace(1,0.001,NnL);
    for k = 1:NnL
        rp = rp+A(k)*sin(2*pi*fL(k)*tp+rand*2*pi);
    end
end

t = (0:TsH:(Ntot-1)*TsH)'; % time signal


% custom multisine
rH = repmat(rp,Per,1);

%% simulate
simoutput = sim('simFile2');
y = simoutput.high(:,1);
uH = simoutput.high(:,2);
eH = simoutput.high(:,3);
uL = simoutput.low(:,2);
eL = simoutput.low(:,1);
%% LPM and ETFE
% [P_LPM,THz] = LPMOpenLoopPeriodicFastBLA(r,y,n,degLPM,Per,nT);
[P_LPM] = LPMClosedLoopPeriodicFastBLA(uH,y,rH,n,degLPM,Per,nT);
%% plotting
figure
Ptrue = bode(Pdt,fH*2*pi);
semilogx(fH,20*log10(squeeze(Ptrue))); hold on;
semilogx(fH,20*log10(abs(squeeze(P_LPM))));
legend('True system DT','LPM')
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
xline(fL(end),'-','Nyquist Low')
xline(fH(end),'-','Nyquist High')

figure
subplot(2,2,1)
RH = fft(rH)/sqrt(Ntot);
stem(fH,abs(RH(1:Per:Per*NnH))); hold on;
set(gca,'xscale','log','yscale','log');
legend('input rH')
xlabel('Frequency [Hz]');
ylabel('abs [-]');
xline(fL(end),'-','Nyquist Low')
xline(fH(end),'-','Nyquist High')

subplot(2,2,2)
EH = fft(eH)/sqrt(Ntot);
EL = fft(eL)/sqrt(NtotL);
stem(fH,abs(EH(1:Per:Per*NnH))); hold on;
stem(fL,abs(EL(1:Per:Per*NnL)));
legend('eH','downsampled error eL')
set(gca,'xscale','log','yscale','log');
xlabel('Frequency [Hz]');
ylabel('abs [-]');
xline(fL(end),'-','Nyquist Low')
xline(fH(end),'-','Nyquist High')

subplot(2,2,3)
UH = fft(uH)/sqrt(Ntot);
UL = fft(uL)/sqrt(NtotL);
stem(fH,abs(UH(1:Per:Per*NnH))); hold on;
stem(fL,abs(UL(1:Per:Per*NnL)));
legend('upsampled input uH','input uL')
set(gca,'xscale','log','yscale','log');
xlabel('Frequency [Hz]');
ylabel('abs [-]');
xline(fL(end),'-','Nyquist Low')
xline(fH(end),'-','Nyquist High')

subplot(2,2,4)
Y = fft(y)/sqrt(Ntot);
stem(fH,abs(Y(1:Per:Per*NnH))); hold on;
legend('output yH')
set(gca,'xscale','log','yscale','log');
xlabel('Frequency [Hz]');
ylabel('abs [-]');
xline(fL(end),'-','Nyquist Low')
xline(fH(end),'-','Nyquist High')
clear all; close all; clc;
[c1, c2, c3, c4, c5,c6,c7] = MatlabDefaultPlotColors();
addpath('../LPM');
%% definitions
TsL = 1/100;
TsH = TsL/4;
tperiod = 10;
Per = 3;
reference = 'low';
%% systems analyzed
s = tf('s');
G = 5.6e3*(0.9*s^2+2*s+5.4e4)/(1.1*0.9*s^4+(1.1+0.9)*2*s^3+(1.1+0.9)*5.4e4*s^2);
C = 0.05*(1/(1*2*pi)*s+1)/(1/(5*2*pi)*s+1);

Pdt = ss(c2d(G,TsH));
Cdt = ss(c2d(C,TsL));
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

farray = logspace(-1,log10(4*fH(end)),40);
farray = [10 20 40];
Nf = length(farray);
A = ones(Nf,1); % amplitude distribution
for k = 1:Nf
    rp(:,k) = A(k)*sin(2*pi*farray(k)*tp+rand*2*pi);
end

t = (0:TsH:(Ntot-1)*TsH)'; % time signal


% custom multisine
ra = repmat(rp,Per,1);

%% simulate
y = zeros(Ntot,Nf);
Pcal = zeros(Nf,1);
Rcal = Pcal;
for k = 1:Nf
    rH = ra(:,k);
    simoutput = sim('simFile2');
    
    y(:,k) = simoutput.output(:,1);
    uH(:,k) = simoutput.output(:,2);
    uL(:,k) = simoutput.output(:,3);
    rL(:,k) = simoutput.output(:,4);
    
    Pcal(k) = powerNorm(y(:,k))/powerNorm(rH);
    Rcal(k) = norm(y(:,k),2)/norm(rH,2);
end

%% plotting
figure(1);clf;
Ptrue = bode(G*C/(1+G*C),fH*2*pi); 
% semilogx(fH,20*log10(squeeze(Ptrue)));hold on;
semilogx(farray,20*log10(Pcal)); hold on;
semilogx(farray,20*log10(Rcal));

legend('True system CT','PFG','RFG')
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
xline(fL(end),'-','Nyquist Low')
xline(fH(end),'-','Nyquist High')

figure
subplot(2,2,1)
RH = fft(rH)/sqrt(Ntot);
RL = fft(rH)/sqrt(NtotL);
stem(fH,abs(RH(1:Per:Per*NnH))); hold on;
stem(fL,abs(RL(1:Per:Per*NnL)));
set(gca,'xscale','log','yscale','log');
legend('input rH','Downsampled error eL')
xlabel('Frequency [Hz]');
ylabel('abs [-]');

subplot(2,2,2)
UH = fft(uH)/sqrt(Ntot);
UL = fft(uH)/sqrt(NtotL);
stem(fH,abs(UH(1:Per:Per*NnH))); hold on;
stem(fL,abs(UL(1:Per:Per*NnL)));
legend('input uH','Downsampled input uL')
set(gca,'xscale','log','yscale','log');
xlabel('Frequency [Hz]');
ylabel('abs [-]');

subplot(2,2,3)
Y = fft(y)/sqrt(Ntot);
stem(fH,abs(Y(1:Per:Per*NnH))); hold on;
legend('output yH')
set(gca,'xscale','log','yscale','log');
xlabel('Frequency [Hz]');
ylabel('abs [-]');

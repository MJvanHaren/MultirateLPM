clear all; close all; clc;
addpath('../LPM/')
%% complex sinusoid
n=3; % window size
degLPM = 2;
per = 1;
nT = per-1;
j = sqrt(-1);
TsH = 1/144; fs = 1/TsH;
Tper = 100;
Tend = per*Tper;
tper = (0:TsH:Tper)';
t = (0:TsH:Tend)';
Nper = length(tper);
N = length(t);
omegaSin = [5*2*pi]; %f0 should be integer for now! Assumption w0<w_{n,L}
% omegaSin = [4*2*pi 10*2*pi 20*2*pi];
% omegaSin = [10*2*pi (10-fs)*2*pi (10-2*fs)*2*pi]; %W_MR
Nsin = length(omegaSin);
omega = linspace(-pi,pi,N)/TsH;
omegaRes = omega(end)-omega(end-1);
f=omega/(2*pi);

rH = zeros(Nper,1);
for i = 1:Nsin
%     rH = rH+exp(j*omegaSin(i)*tper);
    rH = rH+sin(omegaSin(i)*tper);
end
rH = [repmat(rH(1:end-1),per,1);rH(end)];
%% system
s = tf('s');
m = 1;
c = 0.5;
k = 25;
P = 1/(m*s^2+c*s+k);
Pd = c2d(P,TsH);
J = [0 0 1;1 0 0;-1 1 0];
Gd = lft(Pd,J);
%% downsampling
F = 10;
fsL = fs/F;
if max(omegaSin)/2/pi > fsL/2
    error('different frequency(s) plz, this is above low nyquist');
else
    If0 = find(~omega);
    If0 = 0;
    IexcitedFrequenciesBase = [If0+omegaSin/omegaRes (If0+omegaSin/omegaRes)+(N-1)/F:(N-1)/F:0.5*(N-1)];
    IexcitedFrequenciesNeg = [If0+(fsL*2*pi-omegaSin)/omegaRes  (If0+(fsL*2*pi-omegaSin)/omegaRes)+(N-1)/F:(N-1)/F:0.5*(N-1)]; 
    IexFreqs = round(sort([IexcitedFrequenciesBase IexcitedFrequenciesNeg]')); % TODO: fix/remove round
end
rL = rH(1:F:end);
NL = length(rL);
omegaL = linspace(-pi,pi,NL)/(TsH*F);
fL=omegaL/(2*pi);
K = 2;
% uL = C*rL;
rLH = repelem(rL,F);
rLH = rLH(1:end-(F-1));
% yH = lsim(Gd,uH);
%% simulate
simoutput = sim('MRStandardSimulation');
yH = simoutput.zeta;


%% DFT and plotting
% time domain
% figure(1);clf;
% plot(t,real(rH));

figure(2);clf;
subplot(131)
RH = fftshift(fft(rH))/sqrt(N);
% RLH = fftshift(fft(rLH))/sqrt(N);
RL = fftshift(fft(rL))/sqrt(NL);
NUH = fftshift(fft(simoutput.nuph))/sqrt(N);
YH = fftshift(fft(yH))/sqrt(N);
stem(f,abs(RH)); hold on;
stem(fL,abs(RL));
% stem(f,abs(RLH));

legend('sinusoid high freq','Downsampled sinusoid','Upsampled sinusoid')
xline(fL(end))
xline(-fL(end))


subplot(132)
stem(f,abs(YH)); 
xline(fL(end))
xline(-fL(end))
legend('output of system with natural frequency at 0.8Hz')
set(gca,'yscale','log')
set(gca,'xscale','log')

subplot(133)
stem(f,abs(NUH));
xline(fL(end))
xline(-fL(end))
legend('NUP')
set(gca,'yscale','log')
set(gca,'xscale','log')

figure
stem(f,abs(YH));
xline(fL(end))
xline(-fL(end))
legend('output of system with natural frequency at 0.8Hz')
set(gca,'yscale','log')
set(gca,'xscale','log')
%% LPM
% [P_LPM,THz] = LPMOpenLoopPeriodicFastBLA(rH,yH,n,degLPM,(N-1)/F,(N-1)/F-1);
P_LPM = MRLPMOpenLoopFastBLA(rH,yH,n,degLPM,IexFreqs);


fh = linspace(0, 1 - 1/(floor(Nper/2)), (floor(Nper/2))) * (1/TsH)/2;   % available frequencies for high sampling rate
figure
Ptrue = bode(Pd*K/(1+Pd*K),fh*2*pi);
semilogx(fh,20*log10(squeeze(Ptrue))); hold on;
semilogx(f(find(~omega)+IexFreqs),20*log10(abs(squeeze(P_LPM))),'o');
legend('True system DT','LPM')
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
xline(fsL/2,'-','Nyquist Low')
xline(fs/2,'-','Nyquist High')

clear all; close all; clc;
addpath('../LPM/')
%% complex sinusoid
n=5;            % window size
degLPM = 2;     % degree of polynomial estimator
j = sqrt(-1);
fs = 144; TsH = 1/fs;
t = (0:TsH:100-TsH)';
N = length(t);
Nnyquist = floor(N/2);
f = linspace(0, 1 - 1/Nnyquist, Nnyquist) * (1/TsH)/2;
omegaRes = (f(end)-f(end-1))*2*pi; % resolution of freq. grid

F = 10; % down- and upsampling factor
fsL = fs/F; TsL = 1/fsL; % low sampling frequency and corresponding sampling time.
fL = linspace(0, 1 - 1/(Nnyquist/F), Nnyquist/F) * (1/TsL)/2;

omegaSin = f(2:3:720)*2*pi; % input design TODO: check if f=0 can be incorporated (and f=fNyquistLow ?)
Nsin = length(omegaSin);    
rH = zeros(N,1);
for k = 1:Nsin
%     rH = rH+exp(j*omegaSin(i)*tper); %complex sinus
    rH = rH+sin(omegaSin(k)*t+rand*2*pi); % real sinus
end
% rH = [repmat(rH(1:end-1),per,1);rH(end)];
%% system
s = tf('s');
m = 1;
c = 0.5;
k = 25;
P = 1/(m*s^2+c*s+k);
Pd = c2d(P,TsH);
J = [0 0 1;1 0 0;-1 1 0];
Gd = lft(Pd,J);
K = 2;%controller 
%% downsampling

if max(omegaSin)/2/pi > fsL/2
    error('different frequency(s) plz, this is above low nyquist');
else
    IexFreqsPos=[];
     IexFreqsNeg=[];
    for k = 1:length(omegaSin)
        IexFreqsPos = [IexFreqsPos omegaSin(k)/omegaRes (omegaSin(k)/omegaRes)+(N-1)/F:(N-1)/F:0.5*(N-1)];
        IexFreqsNeg = [IexFreqsNeg (fsL*2*pi-omegaSin(k))/omegaRes  ((fsL*2*pi-omegaSin(k))/omegaRes)+(N-1)/F:(N-1)/F:0.5*(N-1)]; 
    end
    IexFreqs = round(sort([IexFreqsPos IexFreqsNeg]')); % TODO: fix/remove round
end
if length(unique(IexFreqs))<length(IexFreqs)
    error('duplicate frequency entered');
end

rL = rH(1:F:end);
NL = length(rL);
rLH = repelem(rL,F);
% rLH = rLH(1:end-(F-1));
%% simulate
simoutput = sim('MRStandardSimulation');
yH = simoutput.zeta;


%% DFT and plotting
% time domain

figure(2);clf;
subplot(131)
RH = fft(rH)/sqrt(N);
RLH = fft(rLH)/sqrt(N);
RL = fft(rL)/sqrt(NL);
NUH = fft(simoutput.nuph)/sqrt(N);
YH = fft(yH)/sqrt(N);
stem(f,abs(RH(1:Nnyquist))); hold on;
stem(fL,abs(RL(1:Nnyquist/F)));
stem(f,abs(RLH(1:Nnyquist)));

legend('sinusoid high freq','Downsampled sinusoid','Down+Upsampled sinusoid')
xline(fL(end))


subplot(132)
stem(f,abs(YH(1:Nnyquist))); 
xline(fL(end))
legend('output of system with natural frequency at 0.8Hz')
set(gca,'yscale','log')
set(gca,'xscale','log')

subplot(133)
stem(f,abs(NUH(1:Nnyquist)));
xline(fL(end))
legend('NUP')
set(gca,'yscale','log')
set(gca,'xscale','log')

figure
stem(f,abs(YH(1:Nnyquist)));
xline(fL(end))
legend('output of system with natural frequency at 0.8Hz')
set(gca,'yscale','log')
set(gca,'xscale','log')
%% LPM
P_LPM = MRLPMOpenLoopFastBLA(rLH,yH,n,degLPM,IexFreqs);
P_LPM_arbitrary = LPMOpenLoopArbitrary(rLH,yH,n,degLPM);
freqSpacing = 128;
P_ETFE_strL = etfe([yH rLH],[],freqSpacing);
P_ETFE_str = etfe([yH rH],[],freqSpacing);
P_ETFE = P_ETFE_str.ResponseData;
P_ETFEL = P_ETFE_strL.ResponseData;
f_ETFE = P_ETFE_str.Frequency/pi*fs/2;

%% plotting
figure(3);clf;
Ptrue = bode(Pd*K/(1+Pd*K),f*2*pi);
semilogx(f,20*log10(squeeze(Ptrue))); hold on;
semilogx(f(IexFreqs),20*log10(abs(squeeze(P_LPM))),'o');
semilogx(f,20*log10(abs(squeeze(P_LPM_arbitrary))));
semilogx(f_ETFE,20*log10(abs(squeeze(P_ETFE))));
semilogx(f_ETFE,20*log10(abs(squeeze(P_ETFEL))));
legend('True system DT','MRLPM','LPM for arbitrary excitations','ETFE with rH','ETFE with down- and upsampled rLH');
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
xline(fsL/2,':','Nyquist Low','HandleVisibility','off');
xline(fs/2,':','Nyquist High','HandleVisibility','off');

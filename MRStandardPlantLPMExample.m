clear all; close all; clc;
addpath('../LPM/')
addpath('../local methods/LPM-LRM/')
[c1,c2,c3,c4,c5,c6,c7] = MatlabDefaultPlotColors();
%% inputs
n = 3;                  % window size
degLPM = 1;             % degree of polynomial estimator
F = 5;                  % down- and upsampling factor
per = 6;                % amunt of periods in signal r(k)
perSkip = 2;            % amount of (first) periods to remove from data 
fs = 250; TsH = 1/fs;   % high/base sampling frequency
Tend = 20-TsH;          % time of simulation
%% definitions
tp = (0:TsH:Tend)';
Np = length(tp);
Nnyquist = floor(Np/2);
f = linspace(0, 1 - 1/Nnyquist, Nnyquist) * (1/TsH)/2;  % TODO what about nyquist?
fsL = fs/F; TsL = 1/fsL;                                % low sampling frequency and corresponding sampling time.
fL = linspace(0, 1 - 1/(Nnyquist/F), Nnyquist/F) * (1/TsL)/2;
NnL = length(fL);
t = (0:TsH:per*(tp(end)+TsH)-TsH)';
tL = (0:TsL:Tend)';
N = length(t);

fSin = f(1:(length(f)/F+1)); % TODO: check if need to add -1? (otherwise duplicates due to aliasing/imaging)
fSinMulti = reshape(fSin(1:end-1),F,[])'; % uncorrelated inputs
% fSinMulti = repmat(fSin',1,F); % correlated inputs
Nsin = size(fSinMulti,1);   

rLift = zeros(length(tL),F);
for k = 1:Nsin
%     rH = rH+exp(j*omegaSin(i)*tper); %complex sinus
%     rH = rH+sin(2*pi*fSin(k)*tp+rand*2*pi); % real sinus
    rLift = rLift+sin(2*pi*tL*fSinMulti(k,:)+rand(1,F)*2*pi);
end
rLift = repmat(rLift,per,1);
rUnLift=[];
for i=1:N/F
    rUnLift = [rUnLift; rLift(i,:)'];
end
rH = rUnLift;
%% system
s = tf('s');
load('newController.mat')
P = 5.6e3*(0.9*s^2+2*s+5.4e4)/(1.1*0.9*s^4+(1.1+0.9)*2*s^3+(1.1+0.9)*5.4e4*s^2);
% K = 0.01*(1/(0.3*2*pi)*s+1)/(1/(3*2*pi)*s+1);
K = shapeit_data.C_tf;

Pd = c2d(P,TsH,'zoh');
Kd = c2d(K,TsL,'zoh');

J = [0 1 1; % r as disturbance at high rate
    1 0 0;
    0 1 1;
    -1 0 0];
Gd = lft(Pd,J);


LiftPd = liftfrd(Pd ,F,f);
%% downsampling
rL = rH(1:F:end);
NL = length(rL);
%% simulate
% simoutput = sim('MRStandardSimulation');
simoutput = sim('simFile2');
yH = simoutput.zeta(:,1);
uH = simoutput.zeta(:,2);
uL = uH(1:F:end);
noise = simoutput.noise;
Noise = fft(noise)/sqrt(N);
%% lifting of signals
yLifted = zeros(N/F-N/F/per*perSkip,F);
uLifted = zeros(N/F-N/F/per*perSkip,F);
rLifted = zeros(N/F-N/F/per*perSkip,F);

for i=N/F/per*perSkip+1:N/F
    yLifted(i-N/F/per*perSkip,:) = yH((i-1)*F+1:i*F);
    uLifted(i-N/F/per*perSkip,:) = uH((i-1)*F+1:i*F);
    rLifted(i-N/F/per*perSkip,:) = rH((i-1)*F+1:i*F);
end
%% lifted LPM
% 1: regular solution: ETFE estimate
PETFE = etfe([yH rH],256,Np/2+1)/etfe([uH rH],256,Np/2+1);
PETFE.FrequencyUnit = 'Hz';
PETFE.Frequency = PETFE.Frequency/pi*fs/2;

% 2: standard closed loop identification using uH,rH and yH, i.e. ignoring LPTV
PLPM = LPMClosedLoopPeriodicFastBLA(uH(1+N/per*perSkip:end),yH(1+N/per*perSkip:end),rH(1+N/per*perSkip:end),n,degLPM,per-perSkip,per-1-perSkip);

% 3: double unlift on rLift-> yLift and rLift-> uLift
PSLifted3 = LPMOpenLoopPeriodicFastBLA(rLifted,yLifted,n,degLPM,per-perSkip,per-1-perSkip);
SLifted3 = LPMOpenLoopPeriodicFastBLA(rLifted,uLifted,n,degLPM,per-perSkip,per-1-perSkip);
for i = 1:F
    for ii=1:F
        PSLifted(i,ii) = frd(squeeze(PSLifted3(i,ii,:))',[fL fs/2/F],TsL,'FrequencyUnit','Hz');
        SLifted(i,ii) = frd(squeeze(SLifted3(i,ii,:))',[fL fs/2/F],TsL,'FrequencyUnit','Hz');
    end
end
%% unlift using Brittani2009 (6.10)
Gori_LPM = frd(squeeze(PLPM)',[f fs/2],TsH,'FrequencyUnit','Hz'); % does not need to unlift

Gori_Lifted = frd(zeros(1,1,length(f)+1),[f fs/2],'FrequencyUnit','Hz');
for k = 1:F
    unLiftedSystems(k) = unliftfrd(PSLifted(:,k),F,[f fs/2],[fL fs/2/F])/unliftfrd(SLifted(:,k),F,[f fs/2],[fL fs/2/F]);
end
table = [1 F:-1:2];
for kk = 1:F % loop over parts [0-f/2/F] [f/2/F-f/F] etc.
    for k = 1:F % loop over columns of lifted system
        if mod(kk,2) % normal, do not mirror indices
            indices = ((kk-1)*NnL)+(k:F:(NnL+1));
        else % 'even' parts, i.e. need to mirror the shifting of samples
            indices = ((kk-1)*NnL)+(table(k):F:(NnL+1));
        end
        Gori_Lifted.ResponseData(1,1,indices) = unLiftedSystems(k).ResponseData(1,1,indices);
    end
end
%% plotting of estimated P_H
[md, pd]=bode(Pd,[f fs/2]*2*pi); 
[mETFE, pETFE]=bode(PETFE);
[mLPM, pLPM]=bode(Gori_LPM);
[mLifted, pLifted]=bode(Gori_Lifted);

figure(1);clf % comparison plot
semilogx([f fs/2],20*log10(abs(squeeze(md))),':','Color',c4); hold on;
plETFE = semilogx([f fs/2],20*log10(abs(squeeze(mETFE))),'.','Color',c1);
plLPM = semilogx([f fs/2],20*log10(abs(squeeze(mLPM))),'o','Color',c2);
plLifted=semilogx([f fs/2],20*log10(abs(squeeze(mLifted))),'^','Color',c3);
xline([fs/F fs/F/2 fs/F*1.5 fs/F*2 fs/F*2.5 fs/F*3],':'); 
legend([plETFE,plLPM,plLifted],{'ETFE solution','Normal LPM on rH, uH, yH','Lifted-MR LPM solution'});
xlabel('Frequency [Hz]')
ylabel('Magntiude [dB]');

figure(2); clf; % difference plot
plETFE = semilogx([f fs/2], 20*log10(abs(squeeze(md)-squeeze(mETFE))),'.'); hold on;
plLPM = semilogx([f fs/2], 20*log10(abs(squeeze(md)-squeeze(mLPM))),'o');
plLifted=semilogx([f fs/2], 20*log10(abs(squeeze(md)-squeeze(mLifted))),'^');
plNoise = semilogx(f, 20*log10(abs(Noise(1:per:per*Np/2))),':');
xline([fs/F fs/F/2 fs/F*1.5 fs/F*2 fs/F*2.5 fs/F*3],':');
legend([plETFE,plLPM,plLifted, plNoise],{'ETFE solution','Normal LPM on rH, uH, yH','Lifted-MR LPM solution','Output noise term'});
xlabel('Frequency [Hz]')
ylabel('Difference between true system and non-parametric estimate [dB]');
%% PFG calucation
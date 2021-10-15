clear all; close all; clc;
addpath('../LPM/')
addpath('../local methods/LPM-LRM/')
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
% P = 1/(s^2+0.05*s+100);
% K = tf(3,1);
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
PLiftedLPM2 = LPMClosedLoopPeriodicFastBLA(uH(1+N/per*perSkip:end),yH(1+N/per*perSkip:end),rH(1+N/per*perSkip:end),n,degLPM,per-perSkip,per-1-perSkip);



% 5: double unlift on rLift-> yLift and rLift-> uLift
PLiftedLPMN1 = LPMOpenLoopPeriodicFastBLA(rLifted,yLifted,n,degLPM,per-perSkip,per-1-perSkip);
PLiftedLPMN2 = LPMOpenLoopPeriodicFastBLA(rLifted,uLifted,n,degLPM,per-perSkip,per-1-perSkip);
for i = 1:F
    for ii=1:F
        frd41(i,ii) = frd(squeeze(PLiftedLPMN1(i,ii,:))',[fL fs/2/F],TsL,'FrequencyUnit','Hz');
        frd42(i,ii) = frd(squeeze(PLiftedLPMN2(i,ii,:))',[fL fs/2/F],TsL,'FrequencyUnit','Hz');
    end
end
%% unlift using Brittani2009 (6.10)
% In = 10;
% LiftPd.ResponseData(:,:,1:In) = LiftPd.ResponseData(:,:,1:In) - repmat(reshape(logspace(4,0,In),1,1,[]),F,F,1).*randn(F,F,In); % disturb original Lifted system
Gori_frd0 = unliftfrd(LiftPd(:,1),F,[f fs/2],[fL fs/2/F]); % system lifted and then unlifted
Gori_frd2 = frd(squeeze(PLiftedLPM2)',[f fs/2],TsH,'FrequencyUnit','Hz'); % does not need to unlift
Gori_frd4 = unliftfrd(frd41(:,4),F,[f fs/2],[fL fs/2/F])/unliftfrd(frd42(:,5),F,[f fs/2],[fL fs/2/F]); 
Gori_frd5 = frd(zeros(1,1,length(f)+1),[f fs/2],'FrequencyUnit','Hz');

for k = 1:F
    tempFRD = unliftfrd(frd41(:,k),F,[f fs/2],[fL fs/2/F])/unliftfrd(frd42(:,k),F,[f fs/2],[fL fs/2/F]); 
    Gori_frd5.ResponseData(1,1,k:F:length(f)+1) = tempFRD.ResponseData(1,1,k:F:end);
end
%% plotting of estimated P_H
opts = bodeoptions;
opts.FreqUnits = 'Hz';

figure(1);clf % comparison plot
xline([fs/F fs/F/2 fs/F*1.5 fs/F*2 fs/F*2.5 fs/F*3],':'); hold on
bode(Pd,[f fs/2]*2*pi,opts); 
bode(Gori_frd0,'o'); 
bode(PETFE,'rx');
bode(Gori_frd2,'y^');
bode(Gori_frd4,'k+');
bode(Gori_frd5,'md');

figure(2); clf; % difference plot
[md, pd]=bode(Pd,[f fs/2]*2*pi,opts); 
[m0, p0]=bode(Gori_frd0); 
[mETFE, pETFE, wETFE]=bode(PETFE);
[m2, p2]=bode(Gori_frd2);
[m4, p4]=bode(Gori_frd4);
[m5, p4]=bode(Gori_frd5);
xline([fs/F fs/F/2 fs/F*1.5 fs/F*2 fs/F*2.5 fs/F*3],':'); hold on
% semilogx([f fs/2], 20*log10(abs(squeeze(md)-squeeze(m1))),'o');
% semilogx([f fs/2], 20*log10(abs(squeeze(md)-squeeze(mETFE))),'.');
semilogx([f fs/2], 20*log10(abs(squeeze(md)-squeeze(m4))),'s');
semilogx([f fs/2], 20*log10(abs(squeeze(md)-squeeze(m5))),'^');
%% PFG calucation

%% OLD: DFT and plotting
% figure(1);clf;
% subplot(121)
% RH = fft(rH)/sqrt(N);
% RLHZOH = fft(rLHZOH)/sqrt(N);
% RLH = fft(rLH)/sqrt(N);
% RL = fft(rL)/sqrt(NL);
% YH = fft(yH)/sqrt(N);
% stem(f,abs(RH(1:per:per*Nnyquist))); hold on;
% stem(fL,abs(RL(1:per:NnL*per)));
% stem(f,abs(RLHZOH(1:per:per*Nnyquist)));
% stem(f,abs(RLH(1:per:per*Nnyquist)));
% legend('sinusoid high freq','Downsampled sinusoid','Down+Upsampled+ZOH sinusoid','Down+Upsampled sinusoid')
% xline(fL(end))
% 
% subplot(122)
% stem(f,abs(YH(1:Nnyquist))); 
% xline(fL(end))
% legend('output of system with natural frequency at 0.8Hz')
% set(gca,'yscale','log')
% set(gca,'xscale','log')


% nu_p,h = u_h = either nu_h or nu_h+omega_h
% zeta_h = [y_h; u_h] = [psi_p,h; either nu_h or nu_h+omega_h]
% psi_h = either r_h-y_h or -y_h = either omega_h-psi_p,h or -psi_p,h


% J = [0  0 1; % r as reference
%      1  0 0;
%      0  0 1;
%      -1 1 0];
% 

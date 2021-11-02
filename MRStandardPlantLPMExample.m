clear all; close all; clc;
addpath('../LPM/')
addpath('../local methods/LPM-LRM/')
[c1,c2,c3,c4,c5,c6,c7] = MatlabDefaultPlotColors();
%% inputs
F = 5;                 % down- and upsampling factor
n = F-1;                % window size
degLPM = 2;             % degree of polynomial estimator
per = 6;                % amunt of periods in signal r(k)
perSkip = 2;            % amount of (first) periods to remove from data (>2 for robust frm!)
fs = 200; TsH = 1/fs;   % high/base sampling frequency
Tend = 30-TsH;          % time of simulation
M = 3;                  % number of realizations of random phase experiments (TODO)
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

%% reference design
fSin = f(1:(length(f)/F+1)); % only exciting low frequency band, wil alias/image to others
Nsin = length(fSin);

j=sqrt(-1);
for p = 1:F
    for q = 1:F
        T(p,q) = F^(-0.5)*exp(j*2*pi*(p-1)*(q-1)/F); % orthogonal multisines (Dobrowiecki et al., 2006) TODO: fix? or change tLow etc?
%         T(p,q) = F^(-0.5)*exp(j*2*pi*(p)*(q)/F); % such that frist column and row are not equal numbers?
    end
end

RSISO = ones(2*NnL,1).*exp(j*rand(2*NnL,1)*2*pi); % random phase (SISO) multisines (Definition 3.1 pintelon2012)

DR = zeros(F,F,2*NnL);
Dphi = zeros(F,F,2*NnL);
for k = 1:2*NnL
    for p = 1:F
        for q = 1:F
            if p==q
                DR(p,q,k) = RSISO(k);
                Dphi(p,q,k) = exp(j*rand*2*pi);
            end
        end
    end
    Rfat(:,:,k) = DR(:,:,k)*T*Dphi(:,:,k); % Pintelon2012 (3-31) full random orthogonal multisines
end

rOrthogonal = real(ifft(Rfat,[],3));
rOrthogonalPeriodic = repmat(rOrthogonal,1,1,per);

for i = 1:F
    rH(:,i) = liftsig(squeeze(rOrthogonalPeriodic(:,i,:))',-F); % inverse lifting for simulation/identification
end

%% system
s = tf('s');
P = 5.6e3*(0.9*s^2+2*s+5.4e4)/(1.1*0.9*s^4+(1.1+0.9)*2*s^3+(1.1+0.9)*5.4e4*s^2);
K = 0.002*(1/(0.3*2*pi)*s+1)/(1/(3*2*pi)*s+1);


Pd = c2d(P,TsH,'zoh');
Kd = c2d(K,TsL,'zoh');

J = [0 1 1; % r as disturbance at high rate
    1 0 0;
    0 1 1;
    -1 0 0];
Gd = lft(Pd,J);
LiftPd = liftfrd(frd(Pd,f,'FrequencyUnit','Hz') ,F,f);
LiftS = liftfrd(frd(1/(1+Pd*d2d(Kd,TsH)),f,'FrequencyUnit','Hz') ,F,f);
LiftPS = liftfrd(frd(Pd/(1+Pd*d2d(Kd,TsH)),f,'FrequencyUnit','Hz') ,F,f);
%% simulate
for i = 1:F
    rHexp = rH(:,i);
    simoutput = sim('simFile2');
    yH(:,i) = simoutput.zeta(:,1);
    uH(:,i) = simoutput.zeta(:,2);
    % lifting of signals
    yLifted(:,:,i) = liftsig(yH(Np*perSkip+1:end,i),F);
    uLifted(:,:,i) = liftsig(uH(Np*perSkip+1:end,i),F);
    rLifted(:,:,i) = liftsig(rH(Np*perSkip+1:end,i),F);
end
noise = simoutput.noise;
Noise = fft(noise)/sqrt(N);

%% Non-parametric estimates
% 1: regular solution: ETFE estimate
PETFE = etfe([yH(1+Np*perSkip:end,1) rH(1+Np*perSkip:end,1)],120,Np/2+1)/etfe([uH(1+Np*perSkip:end,1) rH(1+Np*perSkip:end,1)],120,Np/2+1);
PETFE.FrequencyUnit = 'Hz';
PETFE.Frequency = PETFE.Frequency/pi*fs/2;

% 2: standard closed loop identification using uH,rH and yH, i.e. ignoring LPTV
PLPM = LPMClosedLoopPeriodicFastBLA(uH(1+Np*perSkip:end,1),yH(1+Np*perSkip:end,1),rH(1+Np*perSkip:end,1),n,degLPM,per-perSkip,per-1-perSkip);
Gori_LPM = frd(squeeze(PLPM)',[f fs/2],TsH,'FrequencyUnit','Hz');

% 3: Lift two LPTV open loop transfer functions to LTI transfer functions, LPM, divide and inverse lift
[PSLifted,~,~] = LPMOpenLoopPeriodicRobustFRM(rLifted,yLifted,n,degLPM,per-perSkip,per-1-perSkip,T,Dphi);
[SLifted,~,~] = LPMOpenLoopPeriodicRobustFRM(rLifted,uLifted,n,degLPM,per-perSkip,per-1-perSkip,T,Dphi);

[PSLiftedRep,~,~] = LPMOpenLoopPeriodicRobustFRMRepF(rLifted,yLifted,n,degLPM,per-perSkip,per-1-perSkip,T,Dphi,F);
[SLiftedRep,~,~] = LPMOpenLoopPeriodicRobustFRMRepF(rLifted,uLifted,n,degLPM,per-perSkip,per-1-perSkip,T,Dphi,F); % mirrored for rows?

for i = 1:F % transpose because f*ck matlab
        for ii=1:F
            PSLiftedFRD(i,ii) = frd(squeeze(PSLifted(i,ii,:))',[fL fs/2/F],TsL,'FrequencyUnit','Hz');
            SLiftedFRD(i,ii) = frd(squeeze(SLifted(i,ii,:))',[fL fs/2/F],TsL,'FrequencyUnit','Hz');
            PSLiftedRepFRD(i,ii,:) = frd(squeeze(PSLiftedRep(i,ii,:))',[f fs/2],TsH,'FrequencyUnit','Hz');
            SLiftedRepFRD(i,ii,:) = frd(squeeze(SLiftedRep(i,ii,:))',[f fs/2],TsH,'FrequencyUnit','Hz');
        end
end
PLiftedFRD = PSLiftedFRD/SLiftedFRD;
PMRLiftedLPM = unliftfrd(PLiftedFRD(:,1),F,[f fs/2],[fL fs/2/F]); % unlift using Brittani2009 (6.10)

PLiftedRepFRD = PSLiftedRepFRD/SLiftedRepFRD;
PMRLiftedLPMRep = unliftfrdWithoutRep(PLiftedRepFRD(:,1),F,[f fs/2]);

% 4: same as 3 but with ETFE instead of LPM and using only one experiment (the first)
ryLiftedID = iddata(yLifted(:,:,1),rLifted(:,:,1),TsL);
ruLiftedID = iddata(uLifted(:,:,1),rLifted(:,:,1),TsL);
PSLiftedETFE = etfe(ryLiftedID,120,Np/2/F+1);
SLiftedETFE = etfe(ruLiftedID,120,Np/2/F+1);
PLiftedETFE = PSLiftedETFE/SLiftedETFE;
PETFELifted = unliftfrd(PLiftedETFE(:,1),F,[f fs/2],[fL fs/2/F]); % unlift using Brittani2009 (6.10)
%% plotting of estimated P_H
[md, pd]=bode(Pd,[f fs/2]*2*pi); 
[mETFE, pETFE]=bode(PETFE);
[mLPM, pLPM]=bode(Gori_LPM);
[mLifted, pLifted]=bode(PMRLiftedLPM);
[mLiftedETFE, pLiftedETFE]=bode(PETFELifted);

figure(1);clf % comparison plot
semilogx([f fs/2],20*log10(abs(squeeze(md))),':','Color',c4); hold on;
plETFE = semilogx([f fs/2],20*log10(abs(squeeze(mETFE))),'.','Color',c1);
plLPM = semilogx([f fs/2],20*log10(abs(squeeze(mLPM))),'o','Color',c2);
plLifted=semilogx([f fs/2],20*log10(abs(squeeze(mLifted))),'^','Color',c3);
plLiftedETFE=semilogx([f fs/2],20*log10(abs(squeeze(mLiftedETFE))),'x','Color',c5);
xline([fs/F fs/F/2 fs/F*1.5 fs/F*2 fs/F*2.5 fs/F*3],':'); 
legend([plETFE,plLPM,plLifted,plLiftedETFE],{'ETFE solution','Normal LPM on rH, uH, yH','Lifted-MR LPM solution','Lifted-MR ETFE solution'});
xlabel('Frequency [Hz]')
ylabel('Magntiude [dB]');

figure(2); clf; % difference plot
plETFE = semilogx([f fs/2], 20*log10(abs(squeeze(md)-squeeze(mETFE))),'.'); hold on;
plLPM = semilogx([f fs/2], 20*log10(abs(squeeze(md)-squeeze(mLPM))),'o');
plLifted=semilogx([f fs/2], 20*log10(abs(squeeze(md)-squeeze(mLifted))),'^'); hold on;
plLiftedETFE=semilogx([f fs/2], 20*log10(abs(squeeze(md)-squeeze(mLiftedETFE))),'x');
plNoise = semilogx(f, 20*log10(abs(Noise(1:per:per*Np/2))),':');
xline([fs/F fs/F/2 fs/F*1.5 fs/F*2 fs/F*2.5 fs/F*3],':');
legend([plETFE,plLPM,plLifted,plLiftedETFE, plNoise],{'ETFE solution','Normal LPM on rH, uH, yH','Lifted-MR LPM solution','Lifted-MR ETFE solution','Output noise term'});
xlabel('Frequency [Hz]')
ylabel('Difference between true system and non-parametric estimate [dB]');
%% PFG calucation
% ZOH filter response for cf calculation
z = tf('z',TsH);
Izoh=0;
for fc=0:F-1
    Izoh = Izoh+z^(-fc);
end
IzohResp = freqresp(Izoh,[f fs/2]*2*pi);
PdResp = cat(3,nan,freqresp(Pd,[f(2:end) fs/2]*2*pi));
PdResp = RepSignalFreqDomain(PdResp,F);

% Qd (feedback connection) calculation
Kresp = squeeze(freqresp(Kd,[fL fs/2/F]*2*pi));
Plow = DownSampleFrequencyDomain(squeeze(PMRLiftedLPM.ResponseData),F);
Qd = (1./(1+Kresp.*Plow)).*Kresp; % TODO: fix Plow see paper oomen2007 % kresp validated
QdRep = RepSignalFreqDomain(Qd,F);

% c_f(omega0) calculation
c = zeros(1,F,length(f));
ctrue = zeros(1,F,length(f));
GLiftedRep = RepSignalFreqDomain(squeeze(PMRLiftedLPM.ResponseData),F);
IzohRespRep = RepSignalFreqDomain(squeeze(IzohResp),F);

for k = 1:length(f)
    for fc = 0:F-1
        if fc==0
            c(:,fc+1,k) = GLiftedRep(k)-1/F*GLiftedRep(k)*IzohRespRep(k)*QdRep(k)*GLiftedRep(k); 
        else
            c(:,fc+1,k) = -1/F*GLiftedRep(k+fc/F*length(f))*IzohRespRep(k+fc/F*length(f))*QdRep(k)*GLiftedRep(k); % TODO: change fc/F*length(f)?
        end
    end
end

% Aomega0 & PFG calculation
Aomega0 = zeros(1,1,length(f));
for k = 1:length(f)
    for fc = 1:F
        Aomega0(:,:,k) = Aomega0(:,:,k)+ norm(c(:,fc,k),2).^2;
    end
    PFG(k) = sqrt(max(eig(Aomega0(:,:,k))));
end

%% plotting of alleged PFG
lowComp = inv(1+Kd*d2d(Pd,TsL))*Kd;

figure(3); clf;
plot([fL fs/2/F]*2*pi,20*log10(abs(Qd)),'Color',c2);
hold on
bodemag(d2d(lowComp,TsL));

figure(4);clf;
semilogx(f*2*pi,20*log10(abs(squeeze(c(:,1,:)))),'color',c2); hold on
bodemag(Pd*(1-1/F*Izoh*d2d(lowComp,TsH)*Pd)); % c_0(omega0)

figure(5);clf;
semilogx(f*2*pi,20*log10(abs(PFG)),'color',c2); hold on

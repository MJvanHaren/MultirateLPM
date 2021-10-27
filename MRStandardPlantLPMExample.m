clear all; close all; clc;
addpath('../LPM/')
addpath('../local methods/LPM-LRM/')
[c1,c2,c3,c4,c5,c6,c7] = MatlabDefaultPlotColors();
%% inputs
F = 5;                 % down- and upsampling factor
n = F-1;                % window size
degLPM = 2;             % degree of polynomial estimator
per = 8;                % amunt of periods in signal r(k)
perSkip = 4;            % amount of (first) periods to remove from data 
fs = 200; TsH = 1/fs;   % high/base sampling frequency
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

%% reference design
fSin = f(1:(length(f)/F+1)); % only exciting low frequency band, wil alias/image to others
Nsin = length(fSin);

j=sqrt(-1);
for p = 1:F
    for q = 1:F
%         T(p,q) = F^(-0.5)*exp(j*2*pi*(p-1)*(q-1)/F); % orthogonal multisines (Dobrowiecki et al., 2006) TODO: fix? or change tLow etc?
        T(p,q) = F^(-0.5)*exp(j*2*pi*(p)*(q)/F);
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
%% simulate
% simoutput = sim('MRStandardSimulation');
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
for i = 1:F
    % 1: regular solution: ETFE estimate
    PETFE(i,1) = etfe([yH(1+Np*perSkip:end,i) rH(1+Np*perSkip:end,i)],120,Np/2+1)/etfe([uH(1+Np*perSkip:end,i) rH(1+Np*perSkip:end,i)],120,Np/2+1);

    % 2: standard closed loop identification using uH,rH and yH, i.e. ignoring LPTV
    PLPM(i,:,:) = LPMClosedLoopPeriodicFastBLA(uH(1+Np*perSkip:end,i),yH(1+Np*perSkip:end,i),rH(1+Np*perSkip:end,i),n,degLPM,per-perSkip,per-1-perSkip);
end
PETFE.FrequencyUnit = 'Hz';
PETFE.Frequency = PETFE.Frequency/pi*fs/2;

% 3: double unlift on rLift-> yLift and rLift-> uLift
ZkhPS = zeros(2*F,F,NnL+1);
ZkhS = zeros(2*F,F,NnL+1);
for i = 1:F
    [~,~,~,ZkhPS(:,i,:)] = LPMOpenLoopPeriodicFastBLA(rLifted(:,:,i),yLifted(:,:,i),n,degLPM,per-perSkip,per-1-perSkip);
%     PSLifted3(:,i,:)
    [~,~,~,ZkhS(:,i,:)]= LPMOpenLoopPeriodicFastBLA(rLifted(:,:,i),uLifted(:,:,i),n,degLPM,per-perSkip,per-1-perSkip);
%      SLifted3(:,i,:)
end
for k = 1:NnL+1
   temp(:,:,k) =  ZkhPS(1:F,:,k)/ZkhPS(F+1:end,:,k);
end

for i = 1:F % transpose because f*ck matlab
        for ii=1:F
            PSLifted(i,ii) = frd(squeeze(PSLifted3(i,ii,:))',[fL fs/2/F],TsL,'FrequencyUnit','Hz');
            SLifted(i,ii) = frd(squeeze(SLifted3(i,ii,:))',[fL fs/2/F],TsL,'FrequencyUnit','Hz');
        end
end

% 4: same as 3 but with ETFE instead of LPM
for i = 1:F
    ryLiftedID = iddata(yLifted(:,:,i),rLifted(:,i,i),TsL);
    ruLiftedID = iddata(uLifted(:,:,i),rLifted(:,i,i),TsL);
    PSLiftedETFE(:,i) = etfe(ryLiftedID,120,Np/2/F+1);
    SLiftedETFE(:,i) = etfe(ruLiftedID,120,Np/2/F+1);
end

%% unlift using Brittani2009 (6.10)
Gori_LPM = frd(squeeze(PLPM)',[f fs/2],TsH,'FrequencyUnit','Hz'); % does not need to unlift

Gori_Lifted = frd(zeros(1,1,length(f)+1),[f fs/2],'FrequencyUnit','Hz');
Gori_Lifted_ETFE = frd(zeros(1,1,length(f)+1),[f fs/2],'FrequencyUnit','Hz');
oriSystTest = PSLifted/SLifted;
for k = 1:F
    unLiftedSystems(k) = unliftfrd(PSLifted(:,k),F,[f fs/2],[fL fs/2/F])/unliftfrd(SLifted(:,k),F,[f fs/2],[fL fs/2/F]);
    unLiftedSystemsETFE(k) = unliftfrd(PSLiftedETFE(:,k),F,[f fs/2],[fL fs/2/F])/unliftfrd(SLiftedETFE(:,k),F,[f fs/2],[fL fs/2/F]);
    unLiftedSystemsTest(k) = unliftfrd(oriSystTest(:,k),F,[f fs/2],[fL fs/2/F]);
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
        Gori_Lifted_ETFE.ResponseData(1,1,indices) = unLiftedSystemsETFE(k).ResponseData(1,1,indices);
    end
end
Gori_Lifted = unLiftedSystems(3);
%% plotting of estimated P_H
[md, pd]=bode(Pd,[f fs/2]*2*pi); 
[mETFE, pETFE]=bode(PETFE);
% [mLPM, pLPM]=bode(Gori_LPM);
[mLifted, pLifted]=bode(Gori_Lifted);
[mLiftedETFE, pLiftedETFE]=bode(Gori_Lifted_ETFE);
[mtest, ~]=bode(unLiftedSystemsTest(3));

figure(1);clf % comparison plot
semilogx([f fs/2],20*log10(abs(squeeze(md))),':','Color',c4); hold on;
plETFE = semilogx([f fs/2],20*log10(abs(squeeze(mETFE))),'.','Color',c1);
% plLPM = semilogx([f fs/2],20*log10(abs(squeeze(mLPM))),'o','Color',c2);
plLifted=semilogx([f fs/2],20*log10(abs(squeeze(mLifted))),'^','Color',c3);
plLiftedETFE=semilogx([f fs/2],20*log10(abs(squeeze(mLiftedETFE))),'x','Color',c5);
xline([fs/F fs/F/2 fs/F*1.5 fs/F*2 fs/F*2.5 fs/F*3],':'); 
% legend([plETFE,plLPM,plLifted,plLiftedETFE],{'ETFE solution','Normal LPM on rH, uH, yH','Lifted-MR LPM solution','Lifted-MR ETFE solution'});
xlabel('Frequency [Hz]')
ylabel('Magntiude [dB]');

figure(2); clf; % difference plot
% plETFE = semilogx([f fs/2], 20*log10(abs(squeeze(md)-squeeze(mETFE))),'.'); hold on;
% plLPM = semilogx([f fs/2], 20*log10(abs(squeeze(md)-squeeze(mLPM))),'o');
plLifted=semilogx([f fs/2], 20*log10(abs(squeeze(md)-squeeze(mLifted))),'^'); hold on;
plLiftedETFE=semilogx([f fs/2], 20*log10(abs(squeeze(md)-squeeze(mLiftedETFE))),'x');
plLiftedtest=semilogx([f fs/2], 20*log10(abs(squeeze(md)-squeeze(mtest))),'s');
plNoise = semilogx(f, 20*log10(abs(Noise(1:per:per*Np/2))),':');
xline([fs/F fs/F/2 fs/F*1.5 fs/F*2 fs/F*2.5 fs/F*3],':');
% legend([plETFE,plLPM,plLifted,plLiftedETFE, plNoise],{'ETFE solution','Normal LPM on rH, uH, yH','Lifted-MR LPM solution','Lifted-MR ETFE solution','Output noise term'});
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
Plow = DownSampleFrequencyDomain(squeeze(Gori_Lifted.ResponseData),F);
Qd = (1./(1+Kresp.*Plow)).*Kresp; % TODO: fix Plow see paper oomen2007 % kresp validated
QdRep = RepSignalFreqDomain(Qd,F);

% c_f(omega0) calculation
c = zeros(1,F,length(f));
ctrue = zeros(1,F,length(f));
GLiftedRep = RepSignalFreqDomain(squeeze(Gori_Lifted.ResponseData),F);
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

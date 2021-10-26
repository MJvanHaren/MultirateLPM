clear all; close all; clc;
addpath('../LPM/')
addpath('../local methods/LPM-LRM/')
[c1,c2,c3,c4,c5,c6,c7] = MatlabDefaultPlotColors();
%% inputs
F = 10;                 % down- and upsampling factor
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

fSin = f(1:(length(f)/F+1));
fSinMulti = reshape(fSin(1:end-1),F,[])'; % uncorrelated inputs
% fSinMulti = repmat(fSin',1,F);
Nsin = size(fSinMulti,1);   

rLift = zeros(length(tL),F);
for k = 1:Nsin
    rLift = rLift+sin(2*pi*tL*fSinMulti(k,:)+rand(1,F)*2*pi);
end
rLift = repmat(rLift,per,1);
rH = liftsig(rLift,-F); % inverse lifting for simulation/identification
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
simoutput = sim('simFile2');
yH = simoutput.zeta(:,1);
uH = simoutput.zeta(:,2);
rH = simoutput.omega;
noise = simoutput.noise;
Noise = fft(noise)/sqrt(N);
%% lifting of signals
yLifted = liftsig(yH(Np*perSkip+1:end),F);
uLifted = liftsig(uH(Np*perSkip+1:end),F);
rLifted = liftsig(rH(Np*perSkip+1:end),F);
%% lifted LPM
% 1: regular solution: ETFE estimate
PETFE = etfe([yH(1+Np*perSkip:end) rH(1+Np*perSkip:end)],120,Np/2+1)/etfe([uH(1+Np*perSkip:end) rH(1+Np*perSkip:end)],120,Np/2+1);
PETFE.FrequencyUnit = 'Hz';
PETFE.Frequency = PETFE.Frequency/pi*fs/2;

% 2: standard closed loop identification using uH,rH and yH, i.e. ignoring LPTV
PLPM = LPMClosedLoopPeriodicFastBLA(uH(1+Np*perSkip:end),yH(1+Np*perSkip:end),rH(1+Np*perSkip:end),n,degLPM,per-perSkip,per-1-perSkip);

% 3: double unlift on rLift-> yLift and rLift-> uLift
PSLifted3 = LPMOpenLoopPeriodicFastBLA(rLifted,yLifted,n,degLPM,per-perSkip,per-1-perSkip);
SLifted3 = LPMOpenLoopPeriodicFastBLA(rLifted,uLifted,n,degLPM,per-perSkip,per-1-perSkip);
for i = 1:F
    for ii=1:F
        PSLifted(i,ii) = frd(squeeze(PSLifted3(i,ii,:))',[fL fs/2/F],TsL,'FrequencyUnit','Hz');
        SLifted(i,ii) = frd(squeeze(SLifted3(i,ii,:))',[fL fs/2/F],TsL,'FrequencyUnit','Hz');
    end
end

% 4: same as 3 but with ETFE instead of LPM
ryLiftedID = iddata(yLifted,rLifted,TsL);
ruLiftedID = iddata(uLifted,rLifted,TsL);
PSLiftedETFE = etfe(ryLiftedID,120,Np/2/F+1);
SLiftedETFE = etfe(ruLiftedID,120,Np/2/F+1);

% 5: temp trials
PSLifted51 = MRSIMOLPMOpenLoop(rLifted,yLifted,20,degLPM,per-perSkip,per-1-perSkip);
SLifted51 = MRSIMOLPMOpenLoop(rLifted,uLifted,20,degLPM,per-perSkip,per-1-perSkip);
for i = 1:F
        PSLifted52(i,1) = frd(squeeze(PSLifted51(i,1,:))',[fL fs/2/F],TsL,'FrequencyUnit','Hz');
        SLifted52(i,1) = frd(squeeze(SLifted51(i,1,:))',[fL fs/2/F],TsL,'FrequencyUnit','Hz');
        PLifted5(i,1) = frd(squeeze(PSLifted51(i,1,:)./SLifted51(i,1,:))',[fL fs/2/F],TsL,'FrequencyUnit','Hz');
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
%% plotting of estimated P_H
[md, pd]=bode(Pd,[f fs/2]*2*pi); 
[mETFE, pETFE]=bode(PETFE);
[mLPM, pLPM]=bode(Gori_LPM);
[mLifted, pLifted]=bode(Gori_Lifted);
[mLiftedETFE, pLiftedETFE]=bode(Gori_Lifted_ETFE);

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
plLifted=semilogx([f fs/2], 20*log10(abs(squeeze(md)-squeeze(mLifted))),'^');
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

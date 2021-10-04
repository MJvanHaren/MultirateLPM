clear all; close all; clc;
addpath('../LPM/')
%% complex sinusoid
n=3;            % window size
degLPM = 1;     % degree of polynomial estimator
j = sqrt(-1);
fs = 250; TsH = 1/fs;
tp = (0:TsH:10-TsH)';
Np = length(tp);
per = 4;
t = (0:TsH:per*(tp(end)+TsH)-TsH)';
N = length(t);

Nnyquist = floor(Np/2);
f = linspace(0, 1 - 1/Nnyquist, Nnyquist) * (1/TsH)/2; % TODO what about nyquist?

fRes = (f(end)-f(end-1)); % resolution of freq. grid

F = 2; % down- and upsampling factor
fsL = fs/F; TsL = 1/fsL; % low sampling frequency and corresponding sampling time.
fL = linspace(0, 1 - 1/(Nnyquist/F), Nnyquist/F) * (1/TsL)/2;
NnL = length(fL);

exfIL = 1:1:CommonElemTol(f, fsL/2, fRes/10); % input frequencies on low input spectrum
fSin = f(exfIL); % TODO: check if need to add -1? (otherwise duplicates due to aliasing/imaging)

Nsin = length(fSin);    
rH = zeros(Np,1);
for k = 1:Nsin
%     rH = rH+exp(j*omegaSin(i)*tper); %complex sinus
    rH = rH+sin(2*pi*fSin(k)*tp+rand*2*pi); % real sinus
end
rH = repmat(rH,per,1);

ETFEfreqSpacing =750;
%% system
s = tf('s');

P = 5.6e3*(0.9*s^2+2*s+5.4e4)/(1.1*0.9*s^4+(1.1+0.9)*2*s^3+(1.1+0.9)*5.4e4*s^2);
K = 0.01*(1/(0.3*2*pi)*s+1)/(1/(3*2*pi)*s+1);

% P = 1/(s^2+10*s+100);

Pd = c2d(P,TsH,'zoh');
Kd = c2d(K,TsL);


J = [   0  0 0 1;
        1  0 1 0;
        0  0 0 1;
        -1 1 -1 0];
Gd = lft(Pd,J);

LiftSys = liftfrd(Pd,F,f);
%% downsampling
idx = find(fSin>=fsL/2);

rL = rH(1:F:end);
NL = length(rL);
rLHZOH = repelem(rL,F);
rLH = reshape([rL'; zeros(F-1,length(rL))],1,[])';

%% simulate
simoutput = sim('MRStandardSimulation');
yH = simoutput.zeta(:,1);
uH = simoutput.zeta(:,2);
uL = uH(1:F:end);
nH = simoutput.omegah(:,2);
%% lifting of signals
yLifted = zeros(N/F,F);
rLifted = zeros(N/F,F);
uLifted = zeros(N/F,F);
for i=1:N/F
    yLifted(i,:) = yH((i-1)*F+1:i*F);
    uLifted(i,:) = uH((i-1)*F+1:i*F);
    rLifted(i,:) = rH((i-1)*F+1:i*F);
end
%% DFT and plotting
figure(1);clf;
subplot(121)
RH = fft(rH)/sqrt(N);
RLHZOH = fft(rLHZOH)/sqrt(N);
RLH = fft(rLH)/sqrt(N);
RL = fft(rL)/sqrt(NL);
NUH = fft(simoutput.nuph)/sqrt(N);
YH = fft(yH)/sqrt(N);
stem(f,abs(RH(1:per:per*Nnyquist))); hold on;
stem(fL,abs(RL(1:per:NnL*per)));
stem(f,abs(RLHZOH(1:per:per*Nnyquist)));
stem(f,abs(RLH(1:per:per*Nnyquist)));
legend('sinusoid high freq','Downsampled sinusoid','Down+Upsampled+ZOH sinusoid','Down+Upsampled sinusoid')
xline(fL(end))


subplot(122)
stem(f,abs(YH(1:Nnyquist))); 
xline(fL(end))
legend('output of system with natural frequency at 0.8Hz')
set(gca,'yscale','log')
set(gca,'xscale','log')

%% lifted LPM
[PLiftedLPM] = LPMClosedLoopPeriodicFastBLA(uL,yLifted,rL,n,degLPM,per,per-1); % TODO: take rH or rLifted instead of rL?
% yrdata = iddata(yLifted,rLifted,TsL);
% urdata = iddata(uL,rL,TsL);
% PLiftedETFE = etfe(yrdata,[],250)/etfe(urdata,[],250);
% PLiftedETFE = PLiftedETFE.ResponseData;
Ptrue = bode(Pd,f*2*pi);
PLiftedLPM = squeeze(PLiftedLPM)';
figure(3);clf;
for i = 1:F
        subplot(F,1,i)
        semilogx(fL,20*log10(abs(PLiftedLPM(1:end-1,i)))); hold on;
%         semilogx(fL,20*log10(abs(squeeze(PLiftedETFE(i,:,:)))));
        semilogx(fL,20*log10(squeeze(abs(LiftSys(i,1).ResponseData(1:end-1)))));
end
legend('LPM Lifted system','Wataru lifted system')

temp = frd(LiftSys(:,1).ResponseData(:,:,1:end),[fL fs/2/F],TsL,'FrequencyUnit','Hz');
Gori_frd = unliftfrd(temp,F,[f fs/2],[fL fs/2/F]);

opts = bodeoptions;
opts.FreqUnits = 'Hz';
figure(4);clf
bode(Gori_frd,'o'); hold on;
bode(Pd,f*2*pi,opts);
xline([fs/F fs/F/2 fs/F*1.5 fs/F*2 fs/F*2.5 fs/F*3],':')
%% LPM
% exfIH = find(abs(RLH(1:per*Nnyquist))>1e-10); % TODO: change 1e-6 to variable
% exffIH = find(abs(RLH(1:per:per*Nnyquist))>1e-10); % TODO: change 1e-6 to variable
% Ptrue = bode(Pd,f*2*pi);
% Izoh = 0;
% z = tf('z',TsH);
% for fc=0:F-1
%     Izoh = Izoh+z^(-fc);
% end
% ZOHrespW = abs(squeeze(freqresp(Izoh,(f(1):(f(2)-f(1))/per:fs-(f(2)-f(1))/per)*2*pi)));
% ZOHresp = abs(squeeze(freqresp(Izoh,f*2*pi)));
% ZOHThreshold = db2mag(-10);
% exfIZOH = find(ZOHresp>ZOHThreshold);
% [P_MRLPM_rLHZOH] = MRLPMClosedLoopFastBLA2(uH,yH,rLHZOH,n,degLPM,per,per-1,(exfIZOH-1)*per+1);
% % [P_MRLPM_weighted] = MRLPMClosedLoopFastBLAWeighted(uH,yH,rLHZOH,n,degLPM,per,per-1,ZOHrespW);
% figure(2); clf;
% semilogx(f,20*log10(squeeze(Ptrue))); hold on
% semilogx(f(exfIZOH),20*log10(squeeze(abs(P_MRLPM_rLHZOH))),'o'); 
% % semilogx(f(exffIH),20*log10(squeeze(abs(P_MRLPM_rLHZOH))),'o'); 
% % semilogx(f,20*log10(squeeze(abs(P_MRLPM_weighted))),'o'); 
% semilogx(f,20*log10(squeeze(abs(ZOHresp))),':'); 
% xline(fL(end),':','Nyquist Low')
% yline(mag2db(ZOHThreshold),'--')
% axis([2e-1 f(end) -60 60])

%% functions
function [AI, BI] = CommonElemTol(A, B, Tol)
   A  = A(:);
   B  = B(:);
   nA = numel(A);
   M  = zeros(1, nA);
   
   % Collect the index of the first occurrence in B for every A:
   for iA = 1:nA
      dist = abs(A(iA) - B);             % EDITED: Of course abs() is needed
      Ind  = find(dist < Tol, 1);        % Absolute tolerance
      % Ind = find(dist ./ A(iA) < Tol, 1);  % Relative tolerance
      
      if ~isempty(Ind)
         M(iA) = Ind;
      end
   end
   AI = find(M);        % If any occurrence was found, this A exists
   if isempty(AI)       % prevent: Empty matrix: 0-by-1
      AI = [];
   end
   BI = M(AI);          % at this index in B
end

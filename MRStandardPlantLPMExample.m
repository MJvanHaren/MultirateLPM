clear all; close all; clc;
addpath('../LPM/')
%% complex sinusoid
n=8;            % window size
degLPM = 2;     % degree of polynomial estimator
j = sqrt(-1);
fs = 250; TsH = 1/fs;
tp = (0:TsH:10-TsH)';
Np = length(tp);
per = 6;
t = (0:TsH:per*(tp(end)+TsH)-TsH)';
N = length(t);

Nnyquist = floor(Np/2);
f = linspace(0, 1 - 1/Nnyquist, Nnyquist) * (1/TsH)/2; % TODO what about zhe negative frequencies? (is in some sense taken into account, see line 46-47)
fRes = (f(end)-f(end-1)); % resolution of freq. grid

F = 10; % down- and upsampling factor
fsL = fs/F; TsL = 1/fsL; % low sampling frequency and corresponding sampling time.
fL = linspace(0, 1 - 1/(Nnyquist/F), Nnyquist/F) * (1/TsL)/2;

% fSin = f(2:3:720); % input design TODO: check if f=0 can be incorporated (and f=fNyquistLow ?)
exfIL = 25:1:CommonElemTol(f, fsL/2, fRes/10); % input frequencies on low input spectrum
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
% m = 1;
% c = 0.5;
% k = 250;
% P = 1/(m*s^2+c*s+k);
% Pd = c2d(P,TsH);

P = 5.6e3*(0.9*s^2+2*s+5.4e4)/(1.1*0.9*s^4+(1.1+0.9)*2*s^3+(1.1+0.9)*5.4e4*s^2);
K = 0.01*(1/(0.3*2*pi)*s+1)/(1/(3*2*pi)*s+1);
Pd = c2d(P,TsH);
Kd = c2d(K,TsL);

% J = [0 0 1;1 0 0;-1 1 0];
% J = [0 0 1;1 0 0;0 0 1;-1 1 0];
J = [   0  0 0 1;
        1  0 1 0;
        0  0 0 1;
        -1 1 -1 0];
Gd = lft(Pd,J);
%% downsampling
idx = find(fSin>=fsL/2);

% if length(unique(f(exfI)))<length(f(exfI))
%     error('duplicate frequency entered'); % TODO: FIX, not doing anything now?
% end

rL = rH(1:F:end);
NL = length(rL);
rLHZOH = repelem(rL,F);
rLH = reshape([rL'; zeros(F-1,length(rL))],1,[])';
%% simulate
simoutput = sim('MRStandardSimulation');
yH = simoutput.zeta(:,1);
uH = simoutput.zeta(:,2);
nH = simoutput.omegah(:,2);
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
stem(fL,abs(RL(1:per:per*Nnyquist/F)));
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
%% LPM
exfIH = find(abs(RLH(1:per*Nnyquist))>1e-10); % TODO: change 1e-6 to variable
exffIH = find(abs(RLH(1:per:per*Nnyquist))>1e-10); % TODO: change 1e-6 to variable
Ptrue = bode(Pd,f*2*pi);
% [P_MRLPM_rH] = MRLPMClosedLoopFastBLA2(uH,yH,rH,n,degLPM,per,per-1,exfIH,F);
% [P_MRLPM_rLH] = MRLPMClosedLoopFastBLA2(uH,yH,rLH,n,degLPM,per,per-1,exfIH,F);
[P_MRLPM_rLHZOH] = MRLPMClosedLoopFastBLA2(uH,yH,rLHZOH,n,degLPM,per,per-1,exfIH,F);
figure(2); clf;
semilogx(f,20*log10(squeeze(Ptrue))); hold on
% semilogx(f(exffIH),20*log10(squeeze(abs(P_MRLPM_rH)))); 
% semilogx(f(exffIH),20*log10(squeeze(abs(P_MRLPM_rLH)))); 
semilogx(f(exffIH),20*log10(squeeze(abs(P_MRLPM_rLHZOH))),'o'); 


% P_LPMrH = LPMClosedLoopPeriodicFastBLA(uH,yH,rH,n,degLPM,per,per-1);
% P_LPMrLH = LPMClosedLoopPeriodicFastBLA(uH,yH,rLH,n,degLPM,per,per-1);
% P_LPMrLHZOH = LPMClosedLoopPeriodicFastBLA(uH,yH,rLHZOH,n,degLPM,per,per-1);
% P_ETFE = abs(squeeze(etfe([yH rH],[],ETFEfreqSpacing).ResponseData))./abs(squeeze(etfe([uH rH],[],ETFEfreqSpacing).ResponseData));
% f_ETFE = etfe([yH rH],[],ETFEfreqSpacing).Frequency/pi*fs/2;

% figure(3);clf;
% semilogx(f,20*log10(squeeze(Ptrue))); hold on
% semilogx(f,20*log10(abs(squeeze(P_LPMrH))));
% semilogx(f,20*log10(abs(squeeze(P_LPMrLH))));
% semilogx(f,20*log10(abs(squeeze(P_LPMrLHZOH))));
% semilogx(f_ETFE,20*log10(abs(squeeze(P_ETFE))));
% legend('True plant','LPM with rH','LPM with rLH','LPM with rLHZOH','ETFE')

% figure(4);clf
% semilogx(f,20*log10(abs(abs(squeeze(P_LPMrH))-abs(squeeze(Ptrue)))));hold on;
% semilogx(f,20*log10(abs(abs(squeeze(P_LPMrLH))-abs(squeeze(Ptrue)))));
% semilogx(f,20*log10(abs(abs(squeeze(P_LPMrLHZOH))-abs(squeeze(Ptrue)))));
% legend('LPM with rH','LPM with rLH','LPM with rLHZOH')
%%
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

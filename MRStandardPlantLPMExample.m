clear all; close all; clc;
addpath('../LPM/')
addpath('../local methods/LPM-LRM/')
%% complex sinusoid
n = 3;            % window size
degLPM = 2;       % degree of polynomial estimator
F = 5;            % down- and upsampling factor
j = sqrt(-1);
fs = 250; TsH = 1/fs;
Tend = 20-TsH;
tp = (0:TsH:Tend)';

Np = length(tp);
Nnyquist = floor(Np/2);
f = linspace(0, 1 - 1/Nnyquist, Nnyquist) * (1/TsH)/2;  % TODO what about nyquist?
fsL = fs/F; TsL = 1/fsL;                                % low sampling frequency and corresponding sampling time.
fL = linspace(0, 1 - 1/(Nnyquist/F), Nnyquist/F) * (1/TsL)/2;
NnL = length(fL);
per = 4;
t = (0:TsH:per*(tp(end)+TsH)-TsH)';
tLift = repmat((0:TsL:Tend)',1,F);

N = length(t);
fRes = (f(end)-f(end-1)); % resolution of freq. grid

exfIL = 1:1:CommonElemTol(f, fsL/2, fRes/10); % input frequencies on low input spectrum
fSin = f(exfIL); % TODO: check if need to add -1? (otherwise duplicates due to aliasing/imaging)
fSinMulti = reshape(fSin(1:end-1),F,[])';
% fSinMulti = repmat(fSin',1,F);

Nsin = length(fSinMulti);    
rH = zeros(Np,1);
rLift = zeros(size(tLift));
for k = 1:Nsin
%     rH = rH+exp(j*omegaSin(i)*tper); %complex sinus
%     rH = rH+sin(2*pi*fSin(k)*tp+rand*2*pi); % real sinus
    rLift = rLift+sin(2*pi*tLift(:,1)*fSinMulti(k,:)+rand(1,F)*2*pi);
end
rH = repmat(rH,per,1);
rLift = repmat(rLift,per,1);
rUnLift=[];
for i=1:N/F
    rUnLift = [rUnLift; rLift(i,:)'];
end
rH = rUnLift;
ETFEfreqSpacing =750;
%% system
s = tf('s');
P = 5.6e3*(0.9*s^2+2*s+5.4e4)/(1.1*0.9*s^4+(1.1+0.9)*2*s^3+(1.1+0.9)*5.4e4*s^2);
K = 0.01*(1/(0.3*2*pi)*s+1)/(1/(3*2*pi)*s+1);

Pd = c2d(P,TsH,'zoh');
Kd = c2d(K,TsL,'zoh');

LiftSys = liftfrd(Pd ,F,f);
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
yLifted = zeros(N/F,F);
uLifted = zeros(N/F,F);
rLifted = zeros(N/F,F);

for i=1:N/F
    yLifted(i,:) = yH((i-1)*F+1:i*F);
    uLifted(i,:) = uH((i-1)*F+1:i*F);
    rLifted(i,:) = rH((i-1)*F+1:i*F);
end
%% lifted LPM
% 1: standard closed loop identification using uH,rH and yH, i.e. ignoring LPTV
PLiftedLPM1 = LPMClosedLoopPeriodicFastBLA(uH,yH,rH,n,degLPM,per,per-1);

% 2: closed loop lifted standard. Results good approximation of lifted
% system, but not so good for unlifted
PLiftedLPM2 = LPMClosedLoopPeriodicFastBLA(uLifted,yLifted,rLifted,n,degLPM,per,per-1);
for i = 1:F
    for ii=1:F
        frd2(i,ii) = frd(squeeze(PLiftedLPM2(i,ii,:))',[fL fs/2/F],TsL,'FrequencyUnit','Hz'); % create frd model of response data
    end
end


% 3: ETFE estimate
iddataPS = iddata(yLifted,rLifted,TsL);
iddataS = iddata(uLifted,rLifted,TsL);
PETFE = etfe(iddataPS,10,Np)/etfe(iddataS,10,Np);
%% unlift using Brittani2009 (6.10)
In = 501;
j = sqrt(-1);
% LiftSys.ResponseData(:,:,1:In) = LiftSys.ResponseData(:,:,1:In) - repmat(reshape(logspace(2,-4,In),1,1,[]),F,F,1).*randn(F,F,In)-repmat(reshape(logspace(2,-4,In),1,1,[]),F,F,1).*randn(F,F,In)*j ; % disturb original Lifted system
Gori_frd0 = unliftfrd(LiftSys(:,1),F,[f fs/2],[fL fs/2/F]); % system lifted and then unlifted
Gori_frd1 = frd(PLiftedLPM1,[f fs/2],TsH,'FrequencyUnit','Hz'); % does not need to unlift
Gori_frd2 = unliftfrd(frd2(:,1),F,[f fs/2],[fL fs/2/F]); % equal to solution 3
%% plotting
opts = bodeoptions;
opts.FreqUnits = 'Hz';

figure(2);clf
bode(Pd,[f fs/2]*2*pi,opts); hold on;
bode(Gori_frd0,'o'); hold on;
bode(Gori_frd1,'y^'); hold on;
bode(Gori_frd2,'md'); hold on;
xline([fs/F fs/F/2 fs/F*1.5 fs/F*2 fs/F*2.5 fs/F*3],':')

%% check simulink vs matlab
% check plant difference
difOutput = yH-lsim(Pd,uH);
if (any(abs(difOutput) > 1e-5)), error('Difference in simulink and lsim'); end
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

%% DFT and plotting
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

% J = [0 1 1; % r as disturbance at high rate
%     1 0 0;
%     0 1 1;
%     -1 0 0];
% J = [0  0 1; % r as reference
%      1  0 0;
%      0  0 1;
%      -1 1 0];
% Gd = lft(Pd,J);

clc; clear all; close all;
addpath('../LPM/')
%% complex sinusoid
n = 3;            % window size
degLPM = 2;     % degree of polynomial estimator
fs = 250; Ts = 1/fs;
Tend = 10-Ts;
tp = (0:Ts:Tend)';
Np = length(tp);
Nnyquist = floor(Np/2);
f = linspace(0, 1 - 1/Nnyquist, Nnyquist) * (1/Ts)/2; % TODO what about nyquist?
per = 4;
t = (0:Ts:per*(tp(end)+Ts)-Ts)';

r = zeros(size(tp));
for k = 1:Nnyquist
    r = r+sin(2*pi*f(k)*tp+rand*2*pi);
end
r = repmat(r,per,1);
%%
s = tf('s');

P =[5.6e3*(0.9*s^2+2*s+5.4e4)/(1.1*0.9*s^4+(1.1+0.9)*2*s^3+(1.1+0.9)*5.4e4*s^2);1/(1.3*s^2+0.09*s+98);1/(0.3*s^2+0.2*s+35)];
K = [0.01*(1/(0.3*2*pi)*s+1)/(1/(3*2*pi)*s+1) 0.5 2/3];

Pd = c2d(P,Ts,'zoh');
Kd = c2d(ss(K),Ts);

%% sim
simoutput = sim('simFile');
y = simoutput.data(:,1:3);
u = simoutput.data(:,end);
%% LPM
[PLiftedLPM,Syr,Sur,Psi] = LPMClosedLoopPeriodicFastBLA(u,y,r,n,degLPM,per,per-1); % TODO: take rH or rLifted instead of rL?

for i = 1:3
    Plift(i,:) = frd(squeeze(PLiftedLPM(i,:,:))',[f fs/2],Ts,'FrequencyUnit','Hz');
end

%%
figure
opts = bodeoptions;
opts.frequnits = 'Hz';

bode(Plift,'ro'); hold on
bode(Pd,f*2*pi,opts)

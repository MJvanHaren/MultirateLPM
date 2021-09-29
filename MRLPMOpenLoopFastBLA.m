function [G_LPM] = MRLPMOpenLoopFastBLA(u,y,n,R,Ifreqs)
% This script will calculate a local polynomial model for the given reference, in- and output of a MR system.
% The system is assumed to be in open loop
% Inputs:
%     u : Input to (open loop) plant in discrete time domain, high sampled
%     y : Output of plant to given input signal u, high sampled
%     F : Up- and downsampling factor
%     n : Window size (left and right) for frequency bin k for model
%     R : Degree of polynomial, e.g. G(omega_k+r) = G(omega)+sum_s=1^R g_s(k)*r^s
%     Ifreqs : indices of excited frequencies
% Outputs:
%     G_LPM : Estimated plant dynamics
%% errors
dof = 2*n+1-(R+1)*2;
if dof<1
    error(['Not high enough DOF = ',num2str(dof),' < 1']);
end
%% define variables.
N = length(u);      % Total amount of samples

Nu = size(u,2); % number of inputs
Ny = size(y,2); % number of outputs

if Nu>1 || Ny>1
    error('Not yet defined for MIMO systems');
end

K1 = @(r) (r*ones(R+1,1)).^((0:R)'); % basis for LPM

Uf=fft(u)/sqrt(N); % correct Pintelon2012 (7-66)
Yf=fft(y)/sqrt(N); % correct Pintelon2012 (7-66)
Zf = [Yf';Uf'];

Yk = Yf(Ifreqs,:)'; % up to nyquist frequency (!!!!???) ESSENTIAL
Uk = Uf(Ifreqs,:)';
Zk = [Yk;Uk];       % Pintelon 2012 (7-48)
%% suppress noise transient contribution
Nk = length(Ifreqs);
thetaHat = zeros(Ny+Nu,R+1,Nk);
for k = 1:Nk
    if k ==1
        r = [1:Ifreqs(k)-1 Ifreqs(k)+1:Ifreqs(k+1)-1]-Ifreqs(k);
    elseif k==Nk
        r = [Ifreqs(k-1)+1:Ifreqs(k)-1 Ifreqs(k)+1:floor(length(Yf)/2)-1]-Ifreqs(k); % TODO -1?
    else
        r = [Ifreqs(k-1)+1:Ifreqs(k)-1 Ifreqs(k)+1:Ifreqs(k+1)-1]-Ifreqs(k);
    end
    
    nT =length(r);                  % total window size for noise estimation (left+right size)
    Kn = zeros((R+1),nT);           % reset Kn for every iteration k (+1??)
    for i = 1:nT                    % freq bin leave out middle frequnecy Pintelon2012 (7-72)
        Kn(:,i) = K1(r(i));         % Pintelon2012 between (7-71) and (7-72)
    end
    
    % scaling, see Pintelon2012 (7-25)
    Dscale = zeros(R+1);
    for i = 1:(R+1)
        Dscale(i,i) = norm(Kn(i,:),2);
    end
    
    Kn = Dscale\Kn;
    
    [U_k,S_k,V_k] = svd(Kn'); % better computational feasability Pintelon 2012 (7-24)
    thetaHat(:,:,k) = Zf(:,Ifreqs(k)+r)*U_k/S_k'*V_k'; % use Zf, this is similar to arbitary excitation. 
    thetaHat(:,:,k) = thetaHat(:,:,k)/Dscale;
end
THz = squeeze(thetaHat(:,1,:));
Zkh = Zk-THz;

%% loop over frequency bins
Grz_LPM = zeros(2*Ny,Nu,Nk);
G_LPM = zeros(Ny,Nu,Nk);
Psi = zeros(Ny+Nu,R+1,Nk);
for k = 1:Nk
    if k<n+1 % left border Pintelon2012 (7-29)
        p = n-k+1;
        r=-n+p:n+p;
    elseif k>Nk-n % right border Pintelon2012 (7-29)
        p=-n+Nk-k;
        r=-n+p:n+p;
    else % everything else
        r = -n:n;
    end
    L = zeros(Nu*(R+1),2*n+1); % reset Kn for every iteration k
    for i = 1:2*n+1
        L(:,i) = kron(K1(r(i)),Uk(:,k+r(i))); % yes?
    end
    
    % scaling, see Pintelon2012 (7-25)
    Dscale = zeros(Nu*(R+1));
    for i = 1:Nu*(R+1)
        Dscale(i,i) = norm(L(i,:),2);
    end
    
    L = Dscale\L;
    
    [U_k,S_k,V_k] = svd(L'); % better computational feasability Pintelon 2012 (7-24)
    Psi(:,:,k) = Zkh(:,k+r)*U_k/S_k'*V_k';
    Psi(:,:,k) = Psi(:,:,k)/Dscale;
    Grz_LPM(:,:,k) = Psi(:,1:Nu,k);% calculate LPM estimate of system
    G_LPM(:,:,k) = Grz_LPM(1:Ny,:,k)/(Grz_LPM(end-Nu+1:end,:,k));
end
end


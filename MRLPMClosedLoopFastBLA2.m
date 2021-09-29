function [G_LPM] = MRLPMClosedLoopFastBLA2(u,y,r,n,R,P,nT,exfI,F)
% This script will calculate a local polynomial model for the given reference, in- and output of a system.
% The system is assumed to be in closed loop, see figure 7-4 Pintelon2012.
% Inputs:
%     u : Input to (open loop) plant in time domain
%     y : Output of plant to given input signal u
%     r : Reference to the plant (not disturbance! see Pintelon 2012 Fig. 7-4)
%     n : Window size (left and right) for frequency bin k
%     R : Degree of polynomial, e.g. G(omega_k+r) = G(omega)+sum_s=1^R g_s(k)*r^s
%     P : number of periods in signal
%     nT : window size for noise transient estimation
%     exfI : indices of excited frequencies
%     F : up- and downsampling factor
% Outputs:
%     G_LPM : Estimated plant dynamics
%% errors
dof = 2*n+1-(R+1)*2;
if dof<1
    error(['Not high enough DOF = ',num2str(dof),' < 1']);
end
%% define variables.
N = length(u);      % Total amount of samples
Np = N/P;           % Samples in one period
Nn = length(exfI);  % amount of frequencies excited and below high nyquist

Nu = size(u,2); % number of inputs
Ny = size(y,2); % number of outputs

thetaHat = zeros(2*Ny,(Nu+1)*(R+1),Nn); % Pintelon2012 (7-6)
K1 = @(r) (r*ones(R+1,1)).^((0:R)'); % basis for LPM

Uf=fft(u)/sqrt(N); % correct Pintelon2012 (7-66)
Yf=fft(y)/sqrt(N); % correct Pintelon2012 (7-66)
Rf=fft(r)/sqrt(N); % correct Pintelon2012 (7-66)
Zf = [Yf';Uf'];

superBin = sort(unique(exfI+(-nT:nT)));

% Yk = Yf(superBin,:)'; % up to nyquist frequency, ESSENTIAL
% Uk = Uf(superBin,:)';
% Zk = [Yk;Uk];       % Pintelon 2012 (7-48)
Rk = Rf(exfI,:)';
%% suppress noise transient contribution
thetaHat = zeros(Ny+Nu,R+1,Nn);

for k = 1:Nn
    if k<nT+1                           % left border Pintelon2012 (7-29)
        p = nT-k+1;
        r=[-nT+p:-1+p 1+p:nT+p];        % freq bin leave out middle frequnecy Pintelon2012 (7-72)
    elseif k>Nn-nT                      % right border Pintelon2012 (7-29)
        p=-nT+Nn-k;
        r=[-nT+p:-1+p 1+p:nT+p];        % freq bin leave out middle frequnecy Pintelon2012 (7-72)
    else                                % everything else
        r = [-nT:-1 1:nT];              % freq bin leave out middle frequnecy Pintelon2012 (7-72)
    end
    
    Kn = zeros((R+1),2*nT);          % reset Kn for every iteration k (+1??)
    for i = 1:2*nT                   % freq bin leave out middle frequnecy Pintelon2012 (7-72)
        Kn(:,i) = K1(r(i));          % Pintelon2012 between (7-71) and (7-72)
    end
    
    % scaling, see Pintelon2012 (7-25)
    Dscale = zeros(R+1);
    for i = 1:(R+1)
        Dscale(i,i) = norm(Kn(i,:),2);
    end
    
    Kn = Dscale\Kn;
    
    [U_k,S_k,V_k] = svd(Kn'); % better computational feasability Pintelon 2012 (7-24)
    thetaHat(:,:,k) = Zf(:,exfI(k)+r)*U_k/S_k'*V_k';
    thetaHat(:,:,k) = thetaHat(:,:,k)/Dscale;
end
THz = squeeze(thetaHat(:,1,:));
Zkh = Zf(:,exfI)-THz;


%% loop over frequency bins
Grz_LPM = zeros(2*Ny,Nu,Nn);
G_LPM = zeros(Ny,Nu,Nn);
Psi = zeros(Ny+Nu,R+1,Nn);
subBinIndices = [0; find(diff(exfI)>P)];

for m = 1:length(subBinIndices)
    if m == length(subBinIndices)
        binLength = length(Rk)-subBinIndices(m);
    else
        binLength = subBinIndices(m+1)-subBinIndices(m);
    end
    for k = 1:binLength
        if k<n+1 % left border Pintelon2012 (7-29)
            p = n-k+1;
            r=-n+p:n+p;
        elseif k>binLength-n % right border Pintelon2012 (7-29)
            p=-n+binLength-k;
            r=-n+p:n+p;
        else % everything else
            r = -n:n;
        end
        L = zeros(Nu*(R+1),2*n+1); % reset Kn for every iteration k
        for i = 1:2*n+1
            L(:,i) = kron(K1(r(i)),Rk(:,subBinIndices(m)+k+r(i))); % yes?
        end
        
        % scaling, see Pintelon2012 (7-25)
        Dscale = zeros(Nu*(R+1));
        for i = 1:Nu*(R+1)
            Dscale(i,i) = norm(L(i,:),2);
        end
        
        L = Dscale\L;
        
        [U_k,S_k,V_k] = svd(L'); % better computational feasability Pintelon 2012 (7-24)
        Psi(:,:,subBinIndices(m)+k) = Zkh(:,subBinIndices(m)+k+r)*U_k/S_k'*V_k';
        Psi(:,:,subBinIndices(m)+k) = Psi(:,:,subBinIndices(m)+k)/Dscale;
        Grz_LPM(:,:,subBinIndices(m)+k) = Psi(:,1:Nu,subBinIndices(m)+k);% calculate LPM estimate of system
        G_LPM(:,:,subBinIndices(m)+k) = Grz_LPM(1:Ny,:,subBinIndices(m)+k)/(Grz_LPM(end-Nu+1:end,:,subBinIndices(m)+k));
    end
end
end


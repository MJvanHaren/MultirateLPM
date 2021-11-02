function [G_LPM,THz,Zkh] = LPMOpenLoopPeriodicRobustFRMRepF(u,y,n,R,P,nT,T,Dphi,F)
% This script will calculate a local polynomial model for the given reference, in- and output of a system.
% The system is assumed to be in open loop
% Inputs:
%     u : Input to (open loop) plant in time domain (N*Nu*Nu) (Last Nu here are the amount of experiments, which is equal to Nu)
%     y : Output of plant to given input signal u (N*Ny*Nu)
%     n : Window size (left and right) for frequency bin k for model
%     R : Degree of polynomial, e.g. G(omega_k+r) = G(omega)+sum_s=1^R g_s(k)*r^s
%     T : Orthogonal matrix used for shaping the multiple experiments (Nu*Nu) (Pintelon2012 (2-78))
%     Dphi : frequency dependent diagonal matrix (Nu*Nu(*Nn)) used for shaping the input
% Outputs:
%     G_LPM : Estimated plant dynamics
%     T_LPM : Estimated transient dynamics
%% errors
dof = 2*n+1-(R+1)*2;
if dof<1
    error(['Not high enough DOF = ',num2str(dof),' < 1']);
end
%% define variables.
N = length(u);      % Total amount of samples
Nnori = floor((N/P)/2)+1; % amount of samples per period up to nyquist

Nu = size(u,2); % number of inputs
Ny = size(y,2); % number of outputs

K1 = @(r) (r*ones(R+1,1)).^((0:R)'); % basis for LPM

Nn = (Nnori-1)*F+1; 
% Nn = ((Nnori*P-1)*F+1)/P;

Zkh = zeros(2*Nu,Nu,Nn);

Dphi = repmat(Dphi,1,1,F);

for expNr = 1:Nu
    Uf=fft(u(:,:,expNr))/sqrt(N); % correct Pintelon2012 (7-66)
    Yf=fft(y(:,:,expNr))/sqrt(N); % correct Pintelon2012 (7-66)

    Yk = [RepSignalFreqDomain(Yf(1:floor(N/2)+1,:),F)' Yf(floor(N/2)+2,:)']; % up to nyquist frequency, 1:P:end is correct
    Uk = [ RepSignalFreqDomain(Uf(1:floor(N/2)+1,:),F)' Uf(floor(N/2)+2,:)'];
%     Yk = [RepSignalFreqDomain(Yf(1:P*Nnori,:),F)']; % up to nyquist frequency, 1:P:end is not correct due to flipping. BUT this seems the right implementation
%     Uk = [ RepSignalFreqDomain(Uf(1:P*Nnori,:),F)'];
    
    Zk = [Yk;Uk];       % Pintelon 2012 (7-48)
    %% suppress noise transient contribution
    thetaHat = zeros(Ny+Nu,R+1,Nn); % Pintelon2012 (7-6)
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
        thetaHat(:,:,k) = Zk(:,P*(k-1)+r+1)*U_k/S_k'*V_k';
        thetaHat(:,:,k) = thetaHat(:,:,k)/Dscale;
    end
    THz = squeeze(thetaHat(:,1,:));
    Zkh(:,expNr,:) = Zk(:,1:P:Nn*P)-THz; %TODO: performance is better without -THz?
%     Zkh(:,expNr,:) = Zk(:,1:P:end);
end
ZRkh = zeros(2*Nu,Nu,Nn);
for k = 1:Nn
    ZRkh(:,:,k) = Zkh(:,:,k)*(T*Dphi(:,:,k))';
    G_LPM(:,:,k) = ZRkh(1:Nu,:,k)/ZRkh(Nu+1:end,:,k);
end
end


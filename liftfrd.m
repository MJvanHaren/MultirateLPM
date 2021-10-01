function Glift_frd = liftfrd(G,T,freq)
%liftfrd - Apply lifting for frequency responce data
%
% Glift_frd = liftfrd(G,T,freq)
% G        : Discrete-time frequency response data
%            (or state space model for debug)
% T        : Number of samples for lifting
% freq     : Equidistant frequency grid (for debug with state space)
% Glift_frd: Lifted frequency response data
% 
% Reference: S. Bittanti and P. Colaneri, Periodic systems: Filtering and control, 2009.
%            Section 6.2.1, p.174
% Author   : Wataru Ohnishi (UTokyo&TU/e) and Nard Strijbosch (TU/e), 2020
%%%%

if nargin < 3, freq = G.freq; end
[n,m] = size(freq); if n < m, freq = freq.'; end
flag_frd = isa(G,'frd');

% add DC and Nyquist freq on frequency grid
if freq(1) ~= 0, freq = [0; freq]; end
if freq(end) ~= 1/G.Ts/2, freq = [freq; 1/G.Ts/2]; end

% low-rate z with response data on frequency grid
z = frd(exp(1j*2*pi*freq(1:(length(freq)-1)/T+1)*G.Ts*T),freq(1:(length(freq)-1)/T+1),G.Ts*T,'FrequencyUnit','Hz');

for i = 1:T
    for j = 1:T
        if i < j % top right
            Glift_frd(i,j) = calcGs(G,T,T-abs(i-j),flag_frd,freq);
            Glift_frd(i,j) = Glift_frd(i,j)/z;
        else
            Glift_frd(i,j) = calcGs(G,T,abs(i-j),flag_frd,freq);
        end
    end
end

end

function Gs = calcGs(G,T,s,flag_frd,freq)
fn = 1/G.Ts/2; % original nyquist frequency [Hz]
phi = exp(2*pi*1j/T);

if flag_frd % input is frd 
    freq = G.freq;
    % fill nan at DC and nyquist
    if abs(freq(1) - 0) > eps*100
        G = fcat(G,frd(nan,0,'FrequencyUnit','Hz'));
        freq = [0; freq];
    end
    if abs(freq(end) - fn) > eps*100
        G = fcat(G,frd(nan,fn,'FrequencyUnit','Hz'));
        freq = [freq; fn];
    end
    
    % check if the frequency grid is equidistant
    if ~all(diff(diff(freq)) < eps*1000)
        error('frequency grid must be equidistant grid')
    end
    
    % freqresp for over nyquist
    resp = squeeze(G.resp);
    resp2 = [resp;real(flipud(resp(1:end-1)))-1j*imag(flipud(resp(1:end-1)));resp(2:end)];
%     resp3= [resp;conj(flipud(resp(1:end-1)));resp(2:end)]; % equivalent to resp2
    freq2 = [freq;freq(2:end)+freq(end);freq(2:end)+2*freq(end)];
    
    kk_fn = length(freq)-1; % number of samples up to Nyquist freq
    out = zeros(length(freq),1);
    out(1) = resp(1);
    out(end) = resp(end);
    idx = 2:length(freq);
    for k = 0:T-1
        out(idx) = out(idx) + resp2(idx+k*kk_fn*2/T).*phi.^(k*s).*exp(1j*2*pi.*freq(idx)*G.Ts).^s/T;
    end
    
else % input is model
    out = zeros(length(freq),1);
    for k = 0:T-1
        out = out + squeeze(freqresp(G,freq+k*fn*2/T,'Hz')).*phi.^(k*s).*exp(1j*2*pi.*freq*G.Ts).^s/T;
    end
    
end
% extract up to new Nyquist freq
Gs = frd(out(1:(length(freq)-1)/T+1),freq(1:(length(freq)-1)/T+1),G.Ts*T,'FrequencyUnit','Hz');

end
function Gori_frd = unliftfrdWithoutRep(GliftSIMO,F,freqH)
% unliftfrd, calculating original system using lifted representation
%
% Gori_frd = unliftfrd(GliftSIMO,F,freq)
% GliftSIMO     : Lifted SIMO (!) FRD (!) model (TODO: rewrite for already fully lifted (MIMO) system?)
% F             : Number of samples for lifting
% freq          : Equidistant frequency grid (for debug with state space)
% Gori_frd      : Lifted frequency response data
%
% Reference     : S. Bittanti and P. Colaneri, Periodic systems: Filtering and control, 2009.
%                   Section 6.2.1, p.174, eq 6.10 + W. Ohnishi Multirate State Tracking for Improving Intersample Behavior in Iterative Learning Control, 2021
% Author        : Max van Haren 2021 TU/e
%%%%

% reshape frequency vectors
freqH = reshape(freqH,1,[]);
GliftSIMORespRep = squeeze(GliftSIMO.ResponseData); 

% sampling times
TsH = 1/(freqH(end)-freqH(end-1));                    % sampling time high sampling frequency
NfH = length(freqH);            % amount of frequencies on high frequency grid (=F*(NfL-1)+1)
z = exp(1j*2*pi*freqH*TsH);     % high rate z

% GliftSIMORespRep = RepSignalFreqDomain(GliftResponse',F)';      % repeat amount of times

Gori_resp = zeros(1,NfH); % original system response
for f = 0:F-1
    Gori_resp = Gori_resp + GliftSIMORespRep(f+1,:).*z.^(-f); % (6.10) Bittanti2009 
end
Gori_frd = frd(Gori_resp,freqH,TsH,'FrequencyUnit','Hz'); % original system frd
end
function Gori_frd = unliftfrd(GliftSIMO,F,freqH,freqL)
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
if ~isa(GliftSIMO,'frd'),error('please specify an FRD model as the lifted system input'); end % check if Glift is frd TODO: rewrite for other models?
if nargin < 3, error('specify at least 3 arguments');
elseif nargin < 4, freqL = GliftSIMO.freq; end % if a freq is not specified, take from FRD model
if ~all(diff(diff(freqL)) < eps*100),error('frequency grid must be equidistant'); end % check for equidistant frequency grid

% reshape frequency vectors
freqH = reshape(freqH,1,[]);
freqL = reshape(freqL,1,[]);

% sampling times
TsL = GliftSIMO.Ts;             % sampling time low sampling frequency
TsH = TsL/F;                    % sampling time high sampling frequency
NfL = length(freqL);            % amount of frequencies on low frequency grid
NfH = length(freqH);            % amount of frequencies on high frequency grid (=F*(NfL-1)+1)
z = exp(1j*2*pi*freqH*TsH);    % high rate z

GliftResponse = squeeze(GliftSIMO.ResponseData(:,:,1:NfL));     % response data of lifted system
% GliftResponse(:,1) = GliftResponse(:,1)/2;                    % divide by 2 due to duplicate frequency?

GliftSIMORespRep = zeros(F,NfH); % response of lifted system F times beyond nyquist of low sampling frequency
% GliftSIMORespRep(:,1:NfL) = [nan(F,1) GliftResponse(:,2:end)]; % first block
GliftSIMORespRep(:,1:NfL) = [GliftResponse(:,1:end)]; % first block
for f=1:F-1
    if mod(f,2)==0 % third, fifth, .... block
        GliftSIMORespRep(:,f*(NfL-1)+2:(f+1)*(NfL-1)+1) = GliftResponse(:,2:end);
    else % second, fourth, ... block
%         GliftSIMORespRep(:,f*(NfL-1)+2:(f+1)*(NfL-1)+1) = conj(fliplr([nan(F,1) GliftResponse(:,2:end-1)]));
        GliftSIMORespRep(:,f*(NfL-1)+2:(f+1)*(NfL-1)+1) = conj(fliplr([GliftResponse(:,1:end-1)]));
    end
end

Gori_resp = zeros(1,NfH); % original system response
for f = 0:F-1
    Gori_resp = Gori_resp + GliftSIMORespRep(f+1,:).*z.^(-f); % (6.10) Bittanti2009 TRANSPOSE OR NOT MAKES A BIG DIFFERNECE
end
Gori_frd = frd(Gori_resp,freqH,TsH,'FrequencyUnit','Hz'); % original system frd
end
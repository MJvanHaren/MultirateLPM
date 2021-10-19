function signalRep = RepSignalFreqDomain(signal,reps)
% repeats signal in frequency domain. 
% INPUTS:   signal : input signal which should be repeated (\in \mathbb{R}^{NL \times ni}). should contain DC and nyquist frequency
%           reps   : amount of repetitions
% OUTPUTS:  signalRep : signal, repeated reps times using input signal (\in \mathbb{R}^{(NL-1)*reps+1 \times ni})
NL = length(signal);
signal = reshape(signal,NL,[]);
ni = size(signal,2);

signalRep = zeros((NL-1)*reps+1,ni);
signalRep(1:NL,:) = signal;

for i = 1:reps-1
    if mod(i,2)
    	signalRep(i*(NL-1)+2:(i+1)*(NL-1)+1,:) = conj(flipud(signal(1:NL-1,:)));
    else
        signalRep(i*(NL-1)+2:(i+1)*(NL-1)+1,:) = signal(2:NL,:);
    end
end
end


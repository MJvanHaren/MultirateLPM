function Y = RepSignalFreqDomain(U,reps)
% repeats signal in frequency domain. 
% INPUTS:   U : input signal which should be repeated (\in \mathbb{R}^{NL \times ni}). should contain DC and nyquist frequency
%           reps   : amount of repetitions
% OUTPUTS:  Y : signal, repeated reps times using input signal (\in \mathbb{R}^{(NL-1)*reps+1 \times ni})
NL = length(U);
U = reshape(U,NL,[]);
ni = size(U,2);

Y = zeros((NL-1)*reps+1,ni);
Y(1:NL,:) = U;

for i = 1:reps-1
    if mod(i,2) % second, fourth, sixth, ....
    	Y(i*(NL-1)+2:(i+1)*(NL-1)+1,:) = conj(flipud(U(1:NL-1,:)));
    else % third, fifth, seventh, .....
        Y(i*(NL-1)+2:(i+1)*(NL-1)+1,:) = U(2:NL,:);
    end
end
end


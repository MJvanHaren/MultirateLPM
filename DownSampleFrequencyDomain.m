function downResponse = DownSampleFrequencyDomain(response,F)
% downsamples complex response in frequency domain. 
% INPUTS:   response : input response which should be downsampled (\in \mathbb{R}^{NL \times ni}). should contain DC and nyquist frequency
%           F :  downsampling factor
% OUTPUTS:  downResponse : F times down sampled complex response from response

[NH, ni] = size(response);
NL = (NH-1)/F+1;
downResponse = response(1:NL,:);
negResponse = conj((response));
for f = 1:F-1
%         downResponse = downResponse+conj(flipud(response(1+(f-1)/F*(NH-1):NL+(f-1)/F*(NH-1),:)));
%         downResponse = downResponse+negResponse(1+(f-1)/F*(NH-1):NL+(f-1)/F*(NH-1),:);
        downResponse = downResponse+conj(response(1+(f)/F*(NH-1):NL+(f)/F*(NH-1),:));
end
end


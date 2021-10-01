function out = liftsig(in,T)
% out = liftsig(in,T)
% in : colmum vector
% OR
% in.time 
% in.signals.values
% T : lifting sample
%     T > 0 (lifting)
%     T < 0 (inv lifting)
% out : lifted signal
% Author: Wataru Ohnishi, The University of Tokyo & TU Eindhoven, 2020

if isfield(in,'signals'), data = in.signals.values; else, data = in; end

[Ndata,Dim] = size(data);
if T >= 0 % lifting
    if mod(Ndata,T) ~= 0, error('error in dimension!'); end

    Ndata_lift = Ndata/T;
    Dim_lift = Dim*T;
    data = reshape(data.',[Dim_lift,Ndata_lift]).';

    if isfield(in,'time')
        out.signals.values = data;
        out.time = in.time(1:T:end);
    else
        out = data;
    end
else % inverse lifting 
    data = reshape(data.',[],1);
    if isfield(in,'time')
        Ts = in.time(2)-in.time(1);
        out.signals.values = data;
        out.time = linspace(in.time(1),in.time(end)+(abs(T)-1)*Ts/abs(T),length(data));
    else
        out = data;
    end
end
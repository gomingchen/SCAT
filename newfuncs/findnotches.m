function ipL = findnotches(d, wh)
% function to find the notches in the frequency channels
% d, the above-threshold values
% wh, which threshold to use, there are usually 10 thresholds
% Fc, the frequency bands of all frequency channels

% ipL, the index of notch frequency, Fc(ipL) are the frequencies

W = d(:,wh);
%%% find NaN
[gg1, ~] = find(~isnan(W));
gdif = diff(gg1);
dg = find(gdif>1);
if ~isempty(dg)
    ip_s = gg1(dg);
    ip_e = gg1(dg+1);
    ipL = round((ip_s + ip_e)/2);
    
    %%% eliminate phantom notches
    dL = diff(ipL);
    nk = find(dL<0.5*mode(dL));
    ipL(nk) = [];
    %%%
    
else
    ipL = NaN;
    disp('There are no notches detected. Please use different threhold or check your data');
end

end

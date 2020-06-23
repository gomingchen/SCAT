

function time_in_us = amp_lat_trading(cc, ra, ALT_coef)
% Determine amplitude latency trading effect -25us/dB

% cc, the output waveform after the filterbank
% ra, the threshold used to find the peaks(pulse and echo), in decimal.
% ALT_coef, the coefficient of amplitude latency trading effect


[yu, ~] = envelope(cc,80,'peak');
kk = find(abs(yu-ra)<0.01);

thd = max(yu(kk(1):end))*ra;
[amps,locs] = findpeaks(yu(kk(1):end), 'NPeaks', 2, 'MinPeakHeight',thd);

if length(amps) == 2
    % make sure the peaks found are true peak values of pulse and echo
    amps(1) = max(amps(1), max(yu(kk(1):kk(1)+locs(1)+500)));
    
    if (kk(1)+locs(2)-500)<0
        amps(2) = max(amps(2), max(yu(1: kk(1)+locs(2)+500)));
    else
        amps(2) = max(amps(2), max(yu(kk(1)+locs(2)-500: kk(1)+locs(2)+500)));
    end
    %
    delta_db = 20*log10(amps(1)/amps(2)); % (amp_of_pulse)/(amp_of_echo)
    time_in_us = delta_db*ALT_coef; % latency increase due to the amplitude decrease
else
    time_in_us = NaN;
end


end

function time_in_us = amp_lat_trading(cc, ra, ALT_coef)
% Determine amplitude latency trading effect -25us/dB

% cc, the output waveform after the filterbank
% ra, the threshold used to find the peaks(pulse and echo), in decimal.
% ALT_coef, the coefficient of amplitude latency trading effect


[yu, ~] = envelope(cc,80,'peak');
kk = find(abs(yu-ra)<0.01);

thd = max(yu(kk(1):end))*ra;
[amps,~] = findpeaks(yu(kk(1):end), 'NPeaks', 2, 'MinPeakHeight',thd);

if length(amps) == 2
    delta_db = 20*log10(amps(1)/amps(2)); % (amp_of_pulse)/(amp_of_echo)
    time_in_us = delta_db*ALT_coef; % latency increase due to the amplitude decrease
else
    time_in_us = NaN;
end


end
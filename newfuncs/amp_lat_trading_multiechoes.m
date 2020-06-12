

function [time_in_us, locs] = amp_lat_trading_multiechoes(cc, ra, TH, ALT_coef, N)
% Determine amplitude latency trading effect -25us/dB

% cc, the output waveform after the filterbank
% ra, the threshold used to find the peaks(pulse and echo), in decimal.
% ALT_coef, the coefficient of amplitude latency trading effect


[yu, ~] = envelope(cc,80,'peak');
kk = find(abs(yu-ra)<0.01);

thd = max(yu(kk(1):end))*ra;
THD = TH*.9/100*max(yu(kk(1):end));
u = 1;
v = 0; % v counts the number of times the while loop runs
while u
    % Attention! amps includes broadcast + echoes
    [amps,locs] = findpeaks(yu(kk(1):end), 'NPeaks',N+1, 'MinPeakHeight',THD, 'MinPeakProminence', thd, 'MinPeakDistance', 500);
    if length(amps)<N+1
        u = 1;
        if THD <= thd
            thd = .9*thd;
        else
            THD = .9*THD;
        end
    elseif length(amps)>N+1
            u = 1;
            
            THD = 1.1*THD;
        else
        u = 0;
    end
    v = v + 1;
    if v>30
        disp('too many loops for while, please check your settings!');
        break;
        
    end
end


locs = locs(2:end) + kk(1) - 1;
    if length(amps)>1
        delta_db = 20.*log10(amps(1)./amps(2:end)); % (amp_of_pulse)/(amp_of_echo)
        time_in_us = delta_db.*ALT_coef; % latency increase due to the amplitude decrease
    else
        time_in_us = NaN(size(amps));
    end

end
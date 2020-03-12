function [delayNN, ip_L, Fc] = getGlintSpace2(cfg, ts, wavPS)
%%% Input the config file and time series of the broadcast-echo duo and get
%%% the delay in microseconds between two glints
figure
SIML = runBiscatMain2_noplot(cfg, ts);
Fc = SIML.coch.Fc;
wavPS.simStruct = SIML;


for t = 1:10

    wavPS.NoT = t; % number of thresholds that will be used
    [echoL, firstGapL] = linear_separate_window_10thresholds(wavPS);
    d1(:,t) = firstGapL;
    Mat(t,:) = echoL;
end

%%% Find the delay between target and a receiver (could be a ear in
%%% binaural setting or a bat if it's treated as one receiver)

[~, col] = find(isnan(d1));
[~, ia, ~] = unique(col);

%%% rule out the possibility that there is only one NaN in column col(1)
if ia(1) == 1
    colstart = col(1) + 1;
else
    colstart = col(1);
end

ip_L = findnotches(d1, colstart);
freq_notch = Fc(ip_L);
difF = diff(freq_notch);

%%% Plot sample figure of the dechirped image and found notches
figure
scatter(zeros(length(Fc),1), Fc, 'ob');
hold on,
scatter(d1(:,colstart), Fc, 'ob');
if ~isempty(ip_L)
    scatter(d1(ip_L,1), Fc(ip_L), '*r');
end
text(d1(ip_L,1)-500, Fc(ip_L), 'notch', 'Color', 'red', 'FontSize', 14);
title('Dechirped Sound Image');
set(gca, 'XTickLabel', 1:length(ts.data)/ts.fs*1E3, ...
    'YTick', 20E3:10E3:100E3, 'YTickLabel', 20:10:100, ...
    'FontSize', 12);
xlabel('time (ms)');
ylabel('frequency (kHz)');
hold off

%%%%


if ~isempty(difF)
M = mode(difF);
delayNN = 1/M;

else
    delayNN = NaN;
end

end


function [g, st] = measureWidth(sig, Necho, threshold)
% measure the widths of incoming echoes. If any echo is longer than the
% broadcast more than 500 us, exclude those from next step
% Comment: the glint spacing is the time difference of arrival of two
% echoes from the two glints at the receiver. So if the glint spacing is
% 500us, it's derived from 1/delta_f, where delta_f is the frequency
% difference of neighboring notches. If the distance between two glints is
% d, then glint spacing = 2*d/340 in microseconds.

% inputs, sig, ts in main program
%         Necho, number of echoes
%         threshold, glint spacing threshold in microseconds

% outputs, the index of echoes (!!, so it starts from 1) that have shorter
% glint spacing than 500 us, though I used 450 instead because of errors.
%         st, starting positions (samples) of echoes

% Another comment, sig.data includes broadcast!

[yu, ~] = envelope(sig.data,60,'peak');

kk = find(yu > 0.15*max(yu));
difk = diff(kk);
gk = find(difk>1);

st = [kk(1); kk(gk+1)];
ep = [kk(gk); kk(end)];

width = ep - st;
midp = st + round(0.5*width);
w = zeros(Necho+1,1);
for i = 1:length(width)
    if midp(i)<1E3
        [M,~] = max(sig.data(1:midp(i)+1E3));
        wk = find(sig.data(1:midp(i)+1E3)>.15*M);
    else
        [M,~] = max(sig.data(midp(i)-1E3:midp(i)+1E3));
        wk = find(sig.data(midp(i)-1E3:midp(i)+1E3)>.15*M);
%         figure, plot(sig.data(midp(i)-1E3:midp(i)+1E3));
%         pie = NaN(length(sig.data(midp(i)-1E3:midp(i)+1E3)),1);
%         pie(wk) = 1;
%         hold on, plot(pie,'.r');
    end
    
    w(i) = wk(end)-wk(1);
    
end
W0 = w(1);
% compare the width of each echo with that of the broadcast

g = find(w(2:end)<W0+(threshold-50)*1E-6*sig.fs);


end
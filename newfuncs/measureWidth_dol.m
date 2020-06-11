
function [g, st] = measureWidth_dol(sig, threshold)
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
maxi = max(sig.data);
thresh = .1*maxi;
kk = find(sig.data>=thresh);
st = kk(1);
WLB  = 4E-3*sig.fs;

while ~isempty(kk)

     b = kk(1);
     kk = b + WLB - 1 + find(sig.data(b+WLB:end)>=thresh);
     if ~isempty(kk)
          st = [st kk(1)]; % it's almost guaranteed to be not empty
     end
end

ep = zeros(1, length(st));

for i = 1:length(st)
    ke = st(i) - 1 + find(sig.data(st(i):st(i)+500)>=thresh);
    ep(i) = ke(end);
end

w = ep - st; % widths of the dolphin clicks

W0 = w(1);
% compare the width of each echo with that of the broadcast

g = find(w(2:end)<W0+(threshold-50)*1E-6*sig.fs);


end
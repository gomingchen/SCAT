
function f = splitRange(brd)

% Find the frequency range of first harmonic and the second harmonic
% according to the broadcast, brd
% inputs, brd, the broadcast, first row is the 1st harmonic, 2nd row is 2nd
% harmonic

b1 = brd(:,1); % first harmonic
k1 = ~isnan(b1);
g1 = find(k1);
fs_1H = g1(1);
fe_1H = g1(end);

b2 = brd(:,2); % second harmonic
k2 = ~isnan(b2);
g2 = find(k2);
fs_2H = g2(1);
fe_2H = g2(end);

f = [fs_1H fe_1H; fs_2H fe_2H];
end